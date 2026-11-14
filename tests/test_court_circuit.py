"""Puissance de court-circuit a chaque noeud (:mod:`g2egrid.court_circuit`)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from g2egrid.court_circuit import puissance_court_circuit
from g2egrid.tables import Tables

SB = 100e6
Z_IB = complex(0.001, 0.02)  # pu base Sn = 1 GVA
Z_L = complex(0.02, 0.2)  # pu base Sb


def deux_bus(**modifs) -> Tables:
    """Bus infini (1 GVA) au bus 1, ligne Z_L vers le bus 2 (PQ)."""
    champs = dict(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3, 20e3, 1.0, 0.0, 10e6, 0.0, 0.0, 0.0),
        ],
        T_Lines_pu=[(1, 1, 2, 1.0, Z_L.real, Z_L.imag, 0.0, 20e3, SB)],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, Z_IB.real, Z_IB.imag)],
        T_GFL1s=[(1, 2, 20e6, 40e3)],
        Sb=SB,
    )
    return Tables(**(champs | modifs))


def test_thevenin_analytique_et_scr():
    r = puissance_court_circuit(deux_bus())
    z1 = Z_IB * SB / 1e9
    n1, n2 = r.noeud(1), r.noeud(2)
    assert n1.Scc == pytest.approx(SB / abs(z1), rel=1e-6)
    assert n2.Scc == pytest.approx(SB / abs(z1 + Z_L), rel=1e-6)
    assert n2.XR == pytest.approx((z1 + Z_L).imag / (z1 + Z_L).real, rel=1e-6)
    assert n2.Ik == pytest.approx(n2.Scc / (math.sqrt(3) * 20e3))
    # le GFL de 20 MVA au bus 2 : SCR = Scc / 20 MVA ; seul GFL -> SCRIF = SCR
    assert n2.S_gfl == 20e6 and n2.SCR == pytest.approx(n2.Scc / 20e6)
    assert n1.SCR is None
    assert n2.SCRIF == pytest.approx(n2.SCR)
    assert r.SCR_min == pytest.approx(n2.SCR)


def test_gfl_en_source_de_courant():
    """CEI 60909 : le GFL ajoute Imax.Sn a son propre noeud."""
    base = puissance_court_circuit(deux_bus()).noeud(2).Scc
    avec = puissance_court_circuit(deux_bus(), gfl="courant").noeud(2).Scc
    assert avec - base == pytest.approx(1.2 * 20e6, rel=1e-6)


def test_gfm_trois_hypotheses():
    t = deux_bus(T_GFL1s=[], T_GFM1s=[(1, 2, 20e6, 40e3)])
    sans = puissance_court_circuit(t, gfm="aucun").noeud(2).Scc
    tension = puissance_court_circuit(t, gfm="tension").noeud(2).Scc
    courant = puissance_court_circuit(t, gfm="courant").noeud(2).Scc
    z_f = complex(0.01 + 0.003, 0.1 + 0.01)  # filtre du GFM, base Sn
    assert tension - sans == pytest.approx(20e6 / abs(z_f), rel=0.05)
    assert courant - sans == pytest.approx(1.2 * 20e6, rel=1e-6)
    assert puissance_court_circuit(t).noeud(2).S_gfm == 20e6


def test_hypotheses_rendent_les_modes_demandes():
    """Les hypotheses exportees (atelier, rapport) sont les modes de l'appel.

    Non-regression : la liste des GFL du calcul du SCRIF ecrasait le mode
    ``gfl``, et ``hypotheses["gfl"]`` valait une liste au lieu de "aucun".
    """
    for mode in ("aucun", "courant"):
        r = puissance_court_circuit(deux_bus(), gfl=mode)
        assert r.hypotheses["gfl"] == mode
        assert r.hypotheses["gfm"] == "tension"  # le defaut


def test_scrif_de_deux_gfl():
    t = deux_bus(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3, 20e3, 1.0, 0.0, 10e6, 0.0, 0.0, 0.0),
            (3, 3, 20e3, 1.0, 0.0, 10e6, 0.0, 0.0, 0.0),
        ],
        T_Lines_pu=[
            (1, 1, 2, 1.0, Z_L.real, Z_L.imag, 0.0, 20e3, SB),
            (2, 2, 3, 1.0, Z_L.real, Z_L.imag, 0.0, 20e3, SB),
        ],
        T_GFL1s=[(1, 2, 20e6, 40e3), (2, 3, 10e6, 40e3)],
    )
    r = puissance_court_circuit(t)
    a, b = r.noeud(2), r.noeud(3)
    assert b.Scc < a.Scc  # plus loin de la source
    # reseau radial : la tension du bus 3 suit celle du bus 2 (IF_32 = 1) ;
    # celle du bus 2 ne suit celle du bus 3 qu'en partie (IF_23 = Z_22/Z_33)
    z1 = Z_IB * SB / 1e9
    z22, z33 = abs(z1 + Z_L), abs(z1 + 2 * Z_L)
    assert a.SCRIF == pytest.approx(a.Scc / (20e6 + 10e6), rel=1e-6)
    assert b.SCRIF == pytest.approx(b.Scc / (10e6 + z22 / z33 * 20e6), rel=1e-6)
    assert a.SCRIF < a.SCR and b.SCRIF < b.SCR
    assert r.SCRIF_min == pytest.approx(min(a.SCRIF, b.SCRIF))


def test_ile_sans_source_de_tension():
    """Un bus isole : Scc quasi nul, et pas d'erreur."""
    t = deux_bus()
    t.Y_DER = [*t.Y_DER, (3, 3, 20e3, 1.0, 0.0, 0.0, 0.0, 1e6, 0.2e6)]
    t.T_Loads = []
    r = puissance_court_circuit(t)
    assert r.noeud(3).Scc < 1.0  # VA


def test_vue_par_l_impedance_proche_du_thevenin_sur_un_reseau_passif():
    """Sans convertisseur, l'impedance dq a 1 Hz est celle des branches a
    50 +- 1 Hz : la Scc vue par l'impedance colle a la statique."""
    t = deux_bus(
        T_GFL1s=[],
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3, 20e3, 1.0, 0.0, 0.0, 0.0, 5e6, 1e6),
        ],
    )
    r = puissance_court_circuit(t, dynamique=True)
    n2 = r.noeud(2)
    assert n2.Scc_dyn == pytest.approx(n2.Scc, rel=0.05)


def test_api_court_circuit():
    pytest.importorskip("fastapi")
    from fastapi import HTTPException

    from g2egrid.tables import tables_vers_json
    from g2egrid_api import app as routes
    from g2egrid_api.schemas import CorpsCourtCircuit

    d = routes.court_circuit(CorpsCourtCircuit(tables=tables_vers_json(deux_bus())))
    assert [x["bus"] for x in d["noeuds"]] == [1, 2]
    assert d["noeuds"][1]["SCR"] == pytest.approx(d["SCR_min"])
    assert d["seuils"]["faible"] == 3.0
    assert all(np.isfinite(x["Scc"]) for x in d["noeuds"])
    with pytest.raises(HTTPException) as err:
        routes.court_circuit(
            CorpsCourtCircuit(tables=tables_vers_json(deux_bus()), gfm="magie")
        )
    assert err.value.status_code == 422


# ---------------------------------------------------------------------------
# Machines synchrones : x''d + facteurs de correction de la CEI 60909-0, x'd en
# option. References : un Zbus ecrit a la main (independant du module) et
# pandapower.shortcircuit.
# ---------------------------------------------------------------------------

#: jeu de Kundur des machines du catalogue (SM_DEFAULT), recopie ici a dessein
LL, LAD, LFD, L1D, RA = 0.15, 1.66, 0.165, 0.1713, 0.003
XD_TRANSITOIRE = LL + LAD * LFD / (LAD + LFD)  # x'd = 0,3000
XD_SOUS_TRANSITOIRE = LL + 1.0 / (1 / LAD + 1 / LFD + 1 / L1D)  # x''d = 0,2300
SIN_PHI = math.sqrt(1 - 0.9**2)


def zbus_main(t: Tables, x_machine: float, k_machine=None, k_trafo=None):
    """Zbus des branches de ``t`` (lignes et transformateurs en pu base
    systeme) et des machines ``Ra + j.x`` ; ``k_machine`` et ``k_trafo`` :
    facteurs de correction par bus de machine et par identifiant de transfo."""
    k_machine, k_trafo = k_machine or {}, k_trafo or {}
    n = max(max(r[1], r[2]) for r in [*t.T_Lines_pu, *t.T_Transfos])
    Y = np.zeros((n, n), dtype=complex)
    for _id, i, j, _km, r, x, *_ in t.T_Lines_pu:
        y = 1 / complex(r, x)
        Y[i - 1, i - 1] += y
        Y[j - 1, j - 1] += y
        Y[i - 1, j - 1] -= y
        Y[j - 1, i - 1] -= y
    for tid, i, j, Sn, _u1, _u2, r, x in t.T_Transfos:
        assert Sn == t.Sb  # base du transfo = base systeme
        y = 1 / (complex(r, x) * k_trafo.get(tid, 1.0))
        Y[i - 1, i - 1] += y
        Y[j - 1, j - 1] += y
        Y[i - 1, j - 1] -= y
        Y[j - 1, i - 1] -= y
    for _id, b, Sn in t.T_SMs:
        z = complex(RA, x_machine) * k_machine.get(b, 1.0) * t.Sb / Sn
        Y[b - 1, b - 1] += 1 / z
    return np.linalg.inv(Y + 1e-9 * np.eye(n))


def scc_main(Z: np.ndarray, t: Tables, c: float = 1.0) -> np.ndarray:
    return c / np.abs(np.diag(Z)) * t.Sb


def wscc9() -> Tables:
    from g2egrid.catalogue import charger

    return charger("wscc9")


def test_machine_transitoire_est_l_ancien_calcul():
    """``machine="transitoire"`` : x'd, sans facteur de correction. Reference :
    Zbus a la main en x'd = 0,3000 sur WSCC 9 bus ; et la valeur que rendait le
    code avant le passage a x''d (noeud 4 : 952,22 MVA, relevee le 2026-10-25)."""
    t = wscc9()
    r = puissance_court_circuit(t, machine="transitoire")
    ref = scc_main(zbus_main(t, XD_TRANSITOIRE), t)
    got = np.array([n.Scc for n in r.noeuds])
    assert got == pytest.approx(ref, rel=1e-9)
    assert r.noeud(4).Scc == pytest.approx(952.22e6, rel=1e-5)
    assert r.hypotheses["machine"] == "transitoire"
    assert r.hypotheses["correction_cei"] is False
    assert r.hypotheses["facteurs_cei"] == {}


def test_machine_sous_transitoire_sans_correction():
    """x''d seule (``correction_cei=False``) contre le Zbus a la main en
    x''d = 0,2300, et plus de Scc qu'en x'd (la machine est plus raide)."""
    t = wscc9()
    r = puissance_court_circuit(t, correction_cei=False)
    ref = scc_main(zbus_main(t, XD_SOUS_TRANSITOIRE), t)
    assert np.array([n.Scc for n in r.noeuds]) == pytest.approx(ref, rel=1e-9)
    ancien = puissance_court_circuit(t, machine="transitoire")
    assert all(a.Scc > b.Scc for a, b in zip(r.noeuds, ancien.noeuds))


def test_defaut_est_sous_transitoire_avec_correction():
    r = puissance_court_circuit(wscc9())
    assert r.hypotheses["machine"] == "soustransitoire"
    assert r.hypotheses["correction_cei"] is True
    assert r.hypotheses["cos_phi"] == 0.9


@pytest.mark.parametrize("c", [1.0, 1.1])
def test_blocs_centrale_ks_et_defaut_aux_bornes(c):
    """Chaque machine de WSCC 9 bus est un bloc centrale (machine + transfo
    elevateur) : K_S = c / (1 + |x''d - x_T|.sin(phi)) sur la machine ET son
    transfo pour un defaut cote reseau ; aux bornes de la machine, K_G sur elle
    et le transfo nu, les autres blocs gardant K_S (CEI 60909-0, 6.7.1)."""
    t = wscc9()
    x_t = {tid: x for tid, *_r, x in t.T_Transfos}
    ks = {tid: c / (1 + abs(XD_SOUS_TRANSITOIRE - x_t[tid]) * SIN_PHI) for tid in x_t}
    kg = c / (1 + XD_SOUS_TRANSITOIRE * SIN_PHI)
    r = puissance_court_circuit(t, c=c)
    # le transfo d'identifiant k relie la machine du bus k : memes cles
    Z = zbus_main(t, XD_SOUS_TRANSITOIRE, k_machine=ks, k_trafo=ks)
    ref = scc_main(Z, t, c)
    for b in range(4, 10):  # cote reseau
        assert r.noeud(b).Scc == pytest.approx(ref[b - 1], rel=1e-9), b
    for b in (1, 2, 3):  # aux bornes : K_G sur la machine, transfo nu
        Zb = zbus_main(
            t, XD_SOUS_TRANSITOIRE, k_machine=ks | {b: kg}, k_trafo=ks | {b: 1.0}
        )
        assert r.noeud(b).Scc == pytest.approx(scc_main(Zb, t, c)[b - 1], rel=1e-9), b
    f = r.hypotheses["facteurs_cei"]
    assert f["KS"]["1"] == pytest.approx(ks[1], rel=1e-9)
    assert f["KG"] == {} and f["kT"] == {}


def reseau_direct() -> Tables:
    """Bus infini 110 kV (1 GVA, X/R = 20) - transfo de reseau 110/20 kV (40 MVA,
    Ucc 10 %, R 0,5 %) - ligne de 5 km - machine de 50 MVA en direct."""
    return Tables(
        Sb=SB,
        Y_DER=[
            (1, 1, 110e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3, 20e3, 1.0, 0.0, 0.0, 0.0, 5e6, 1e6),
            (3, 2, 20e3, 1.0, 0.0, 30e6, 0.0, 0.0, 0.0),
        ],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.02)],
        T_Transfos=[(1, 1, 2, 40e6, 110e3, 20e3, 0.005, 0.10)],
        T_Lines_SI=[(1, 2, 3, 5.0, 0.2, 0.35, 0.0)],
        T_SMs=[(1, 3, 50e6)],
    )


def test_generateur_direct_kg_et_transformateur_de_reseau_kt():
    """Generateur raccorde directement : K_G ; transformateur de reseau : K_T =
    0,95.c / (1 + 0,6.x_T). Reference : Zbus a la main."""
    t = reseau_direct()
    c = 1.1
    kg = c / (1 + XD_SOUS_TRANSITOIRE * SIN_PHI)
    kt = 0.95 * c / (1 + 0.6 * 0.10)
    z_l = complex(0.2, 0.35) * 5 / (20e3**2 / SB)
    z_t = complex(0.005, 0.10) * (SB / 40e6)
    z_ib = complex(0.001, 0.02) * SB / 1e9
    z_m = complex(RA, XD_SOUS_TRANSITOIRE) * kg * SB / 50e6
    Y = np.zeros((3, 3), dtype=complex)
    for i, j, z in ((0, 1, z_t * kt), (1, 2, z_l)):
        Y[i, i] += 1 / z
        Y[j, j] += 1 / z
        Y[i, j] -= 1 / z
        Y[j, i] -= 1 / z
    Y[0, 0] += 1 / z_ib
    Y[2, 2] += 1 / z_m
    ref = c / np.abs(np.diag(np.linalg.inv(Y))) * SB
    r = puissance_court_circuit(t, c=c)
    assert np.array([n.Scc for n in r.noeuds]) == pytest.approx(ref, rel=1e-9)
    f = r.hypotheses["facteurs_cei"]
    assert f["KG"]["3"] == pytest.approx(kg) and f["KS"] == {}
    assert f["kT"]["T1_1_2"] == pytest.approx(kt)


def test_api_transmet_l_option_machine():
    pytest.importorskip("fastapi")
    from g2egrid.tables import tables_vers_json
    from g2egrid_api import app as routes
    from g2egrid_api.schemas import CorpsCourtCircuit

    tj = tables_vers_json(wscc9())
    d = routes.court_circuit(CorpsCourtCircuit(tables=tj))
    assert d["hypotheses"]["machine"] == "soustransitoire"
    d = routes.court_circuit(CorpsCourtCircuit(tables=tj, machine="transitoire"))
    assert d["hypotheses"]["machine"] == "transitoire"
    assert d["noeuds"][3]["Scc"] == pytest.approx(952.22e6, rel=1e-5)


def test_machine_et_cos_phi_invalides():
    with pytest.raises(ValueError, match="machine"):
        puissance_court_circuit(reseau_direct(), machine="magie")
    with pytest.raises(ValueError, match="cos_phi"):
        puissance_court_circuit(reseau_direct(), cos_phi=0.0)


def test_wscc9_contre_pandapower():
    """K_S, K_G (aux bornes) et x''d sur WSCC 9 bus, contre
    ``pandapower.shortcircuit`` (case max, c = 1,1 ; pandapower 3.5.6 :
    ``uv run --with pandapower pytest ...``). Ecart mesure : < 1e-5 sur chaque
    noeud, cote reseau comme aux bornes des machines ; le meme calcul sans
    facteur de correction s'ecarte de 1 a 2 % (ce que les facteurs corrigent)."""
    pp = pytest.importorskip("pandapower")
    import pandapower.shortcircuit as sc

    t = wscc9()
    un = {1: 16.5, 2: 18.0, 3: 13.8}
    net = pp.create_empty_network(sn_mva=t.Sb / 1e6, f_hz=50)
    bus = {b: pp.create_bus(net, vn_kv=un.get(b, 230.0)) for b in range(1, 10)}
    zb = 230.0**2 / (t.Sb / 1e6)
    for _id, i, j, _km, r, x, *_ in t.T_Lines_pu:
        pp.create_line_from_parameters(
            net, bus[i], bus[j], 1.0, r * zb, x * zb, 0.0, max_i_ka=1.0
        )
    trafo = {}
    for tid, i, j, Sn, u1, u2, r, x in t.T_Transfos:
        trafo[tid] = pp.create_transformer_from_parameters(
            net, bus[j], bus[i], Sn / 1e6, u2 / 1e3, u1 / 1e3,
            vk_percent=100 * math.hypot(r, x), vkr_percent=100 * r,
            pfe_kw=0.0, i0_percent=0.0, oltc=True,
        )  # fmt: skip
    for sid, b, Sn in t.T_SMs:
        pp.create_gen(
            net, bus[b], p_mw=0.0, vm_pu=1.0, sn_mva=Sn / 1e6, vn_kv=un[b],
            cos_phi=0.9, xdss_pu=XD_SOUS_TRANSITOIRE,
            rdss_ohm=RA * un[b] ** 2 / (Sn / 1e6), power_station_trafo=trafo[sid],
        )  # fmt: skip
    net.trafo["power_station_unit"] = True  # pas de K_T en plus du K_S
    sc.calc_sc(net, case="max", fault="3ph", ip=False, ith=False)
    avec = puissance_court_circuit(t, c=1.1)
    sans = puissance_court_circuit(t, c=1.1, correction_cei=False)
    for n in avec.noeuds:
        ik_pp = net.res_bus_sc.ikss_ka.at[bus[n.bus]] * 1e3
        coef = math.sqrt(3) * un.get(n.bus, 230.0) * 1e3
        assert n.Scc / coef == pytest.approx(ik_pp, rel=1e-5), n.bus
        ecart = sans.noeud(n.bus).Scc / coef / ik_pp - 1
        assert 0.001 < ecart < 0.03, (n.bus, ecart)


def test_direct_et_transformateur_contre_pandapower():
    """K_G (generateur en direct) et K_T (transfo de reseau) contre
    pandapower ; le bus infini devient un ext_grid de S''k = c.Sn/|z| (pandapower
    multiplie son impedance par c, le module non : meme impedance au final).
    Ecart mesure : < 1e-9 aux trois noeuds ; sans facteur de correction, 1,1 % au
    noeud 2 (K_T) et 0,2 % au noeud 3 (K_G)."""
    pp = pytest.importorskip("pandapower")
    import pandapower.shortcircuit as sc

    t = reseau_direct()
    c = 1.1
    net = pp.create_empty_network(sn_mva=100.0, f_hz=50)
    b1 = pp.create_bus(net, vn_kv=110.0)
    b2 = pp.create_bus(net, vn_kv=20.0)
    b3 = pp.create_bus(net, vn_kv=20.0)
    pp.create_ext_grid(
        net, b1, s_sc_max_mva=c * 1000.0 / abs(complex(0.001, 0.02)), rx_max=0.05
    )
    pp.create_transformer_from_parameters(
        net, b1, b2, 40.0, 110.0, 20.0, vk_percent=100 * math.hypot(0.005, 0.10),
        vkr_percent=0.5, pfe_kw=0.0, i0_percent=0.0,
    )  # fmt: skip
    pp.create_line_from_parameters(net, b2, b3, 5.0, 0.2, 0.35, 0.0, max_i_ka=1.0)
    pp.create_gen(
        net, b3, p_mw=30.0, vm_pu=1.0, sn_mva=50.0, vn_kv=20.0, cos_phi=0.9,
        xdss_pu=XD_SOUS_TRANSITOIRE, rdss_ohm=RA * 20.0**2 / 50.0,
    )  # fmt: skip
    sc.calc_sc(net, case="max", fault="3ph", ip=False, ith=False)
    r = puissance_court_circuit(t, c=c)
    for bus_pp, n, un in ((b1, 1, 110e3), (b2, 2, 20e3), (b3, 3, 20e3)):
        ik_pp = net.res_bus_sc.ikss_ka.at[bus_pp] * 1e3
        assert r.noeud(n).Scc / (math.sqrt(3) * un) == pytest.approx(ik_pp, rel=1e-6)
