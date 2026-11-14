"""Disjoncteurs, ilotage, repere a part, manoeuvres en simulation."""

from __future__ import annotations

import json

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

import g2egrid.catalogue as cat
from g2egrid.disjoncteurs import basculer, en_service, nom_element, ouvrables
from g2egrid.elements import Repere
from g2egrid.tables import assembler, tables_depuis_json, tables_vers_json
from g2egrid.temporel import Manoeuvre, ModeleNonLineaire, mesures
from g2egrid.zones import aplatir, reseau_depuis_json, reseau_vers_json
from g2egrid_api import analyse

# BDF2 a pas fixe : une manoeuvre excite les resonances du reseau, que Radau
# suivrait a la tolerance pres (cf. g2egrid.integrateur)
RAPIDE = dict(methode="BDF2", pas=1e-4)


@pytest.fixture(scope="module")
def wscc9():
    return cat.charger("wscc9")


# ------------------------------------------------------------ regime permanent
def test_ouvrir_une_ligne(wscc9):
    t = basculer(wscc9, "T_Lines_pu", 3, True)  # L3 : 5-7
    net = assembler(t)
    assert "L3_5_7" not in {e.name for e in net.elements}
    # la ligne emporte sa demi-charge capacitive : il ne reste que celle de 4-5
    assert net.bus_specs[5].C == pytest.approx(0.176 / 2)
    assert net.run_loadflow().converged
    assert wscc9.OUVERTS == []  # l'original n'a pas bouge


def test_refermer_rend_le_reseau_d_origine(wscc9):
    t = basculer(basculer(wscc9, "T_SMs", 2, True), "T_SMs", 2, False)
    assert t.OUVERTS == []
    a = assembler(wscc9).build_state_space().A
    b = assembler(t).build_state_space().A
    assert np.array_equal(a, b)


def test_declencher_le_slack(wscc9):
    t, rap = en_service(basculer(wscc9, "T_SMs", 1, True))
    assert int(t.bus(1)[1]) == 3 and int(t.bus(2)[1]) == 1  # SM2, le plus gros
    assert rap.iles == [{"bus": list(range(1, 10)), "reference": 2}]
    d = analyse.t_eig(basculer(wscc9, "T_SMs", 1, True))
    assert d["stable"]


def test_ile_sans_formeur_hors_tension(wscc9):
    """Le bus 6 ne tient qu'a 4-6 et 6-9 : les ouvrir le laisse sans source."""
    t = basculer(basculer(wscc9, "T_Lines_pu", 2, True), "T_Lines_pu", 4, True)
    eff, rap = en_service(t)
    assert rap.hors_tension == [6]
    assert 6 not in {int(r[0]) for r in eff.Y_DER}
    assert assembler(t).run_loadflow().converged


def test_ilot_avec_formeur(wscc9):
    """Ouvrir le transfo de SM3 : deux iles vivantes, donc deux jauges."""
    t = basculer(wscc9, "T_Transfos", 3, True)
    _, rap = en_service(t)
    assert [i["bus"] for i in rap.iles] == [[1, 2, 4, 5, 6, 7, 8, 9], [3]]
    assert analyse.t_eig(t)["jauge"] == 2
    msgs = analyse.valider(t)
    assert not any(m["gravite"] == "erreur" for m in msgs)


def test_disjoncteurs_en_json_et_en_zones():
    t = basculer(cat.charger("wscc9"), "T_Lines_pu", 1, True)
    assert tables_depuis_json(tables_vers_json(t)).OUVERTS == [("T_Lines_pu", 1)]

    r = cat.charger_reseau("wscc3sm_zones")
    (z,) = r.zones.values()
    z.tables.OUVERTS = [("T_Lines_SI", 1)]
    plat = aplatir(r)
    assert sorted(plat.OUVERTS) == [
        ("T_Lines_SI", 101),
        ("T_Lines_SI", 201),
        ("T_Lines_SI", 301),
    ]
    relu = reseau_depuis_json(reseau_vers_json(r))
    assert list(relu.zones.values())[0].tables.OUVERTS == [("T_Lines_SI", 1)]


def test_noms_des_elements(wscc9):
    noms = {o["nom"] for o in ouvrables(wscc9)}
    assert {"L3_5_7", "T3_3_9", "SM2", "Load5"} <= noms
    assert nom_element(wscc9, "T_Transfos", 1) == "T1_1_4"


# ------------------------------------------------------------ repere a part
def test_repere_a_part_meme_spectre_plus_un_zero():
    """Un Repere qui suit la machine de reference : memes modes, plus un zero
    (l'ecart d'angle entre le repere et la machine n'a pas de rappel)."""
    net0 = cat.construire("wscc9")
    ss0 = net0.build_state_space()
    m0 = ModeleNonLineaire(net0)
    theta_g = m0._angle_repere(m0.c.u0())

    net = cat.construire("wscc9")
    net.solve_equilibrium()
    rep = Repere("Repere", wb=net.wb)
    rep.theta0 = theta_g
    rep.x0, rep.z0 = np.array([theta_g]), np.zeros(0)
    rep.u0 = np.array([1.0, 1.0, theta_g])
    net.elements.append(rep)
    net.reference, net.repere_suit = "Repere", "SM1"
    assert ModeleNonLineaire(net).residu_equilibre()[0] < 1e-6
    ss1 = net.build_state_space()
    assert len(ss1.states) == len(ss0.states) + 1

    l0 = np.append(np.linalg.eigvals(ss0.A), 0.0)
    l1 = np.linalg.eigvals(ss1.A)
    cout = np.abs(l0[:, None] - l1[None, :]) / np.maximum(1.0, np.abs(l0[:, None]))
    i, j = linear_sum_assignment(cout)
    assert cout[i, j].max() < 1e-6


# ------------------------------------------------------------ manoeuvres
def _simuler(man, duree=0.25, n=51, quoi=()):
    m = ModeleNonLineaire(cat.construire("wscc9"))
    u = m.scenario(man)
    tr = m.simuler((0.0, duree), uext_fn=u, t_eval=np.linspace(0, duree, n), **RAPIDE)
    return m, tr, (mesures(m, tr, u, quoi=list(quoi)) if quoi else {})


def test_ouverture_de_ligne_en_simulation():
    m, tr, mes = _simuler(Manoeuvre(0.05, "ouverture", "L3_5_7"), quoi=["V_bus5"])
    # deux RESEAUX (deux cablages) ; les segments peuvent etre plus nombreux,
    # une butee d'excitation qui joue en ajoute (cf. test_limiteurs)
    assert len({id(sg[2].c) for sg in tr.segments}) == 2
    apres = tr.segments[-1][2]
    assert "L3_5_7" not in {e.name for e in apres.c.els}
    # la ligne a repris sa demi-charge capacitive aux deux bouts
    assert apres.net.buses[5].pdict["C"] == pytest.approx(
        m.net.buses[5].pdict["C"] - 0.306 / 2  # la demi-charge de 5-7
    )
    i = tr.etats.index("L3_5_7.i_d")
    assert np.isfinite(tr.x[:10, i]).all() and np.isnan(tr.x[-1, i])
    assert np.isfinite(mes["V_bus5"]).all()


def test_declenchement_d_une_machine_en_simulation():
    """La production perdue est reprise : la frequence baisse, SM3 monte."""
    _m, tr, mes = _simuler(
        Manoeuvre(0.05, "ouverture", "SM2"), duree=0.6, n=61, quoi=["fm_bus5", "P_SM3"]
    )
    assert mes["P_SM3"][-1] > mes["P_SM3"][0] + 0.1
    assert mes["fm_bus5"][-1] < 60.0 - 0.02
    assert np.isnan(tr.x[-1, tr.etats.index("SM2.Dwr")])


def test_declenchement_du_slack_passe_la_main_au_repere():
    _m, tr, _ = _simuler(Manoeuvre(0.05, "ouverture", "SM1"))
    apres = tr.segments[-1][2]
    assert apres.net.reference == "Repere" and apres.net.repere_suit == "SM2"
    # l'angle du repere prolonge celui de la machine de reference, sans saut
    k = tr.segments[0][1]  # premier instant apres la manoeuvre
    th_ref = tr.x[k - 1, tr.etats.index("SM1.theta")]
    th_rep = tr.x[k, tr.etats.index("Repere.theta")]
    assert abs(th_rep - th_ref) < 1e-3
    assert np.isfinite(tr.x[-1, tr.etats.index("SM3.Dwr")])


def test_echelon_de_charge():
    """+10 % : l'impedance de la charge change a la tension d'equilibre."""
    m = ModeleNonLineaire(cat.construire("wscc9"))
    x, z, u = m.c.x0(), m.c.z0(), m.c.u0()
    m2, *_ = m.appliquer(Manoeuvre(0.0, "charge", "Load5", dp=0.1, dq=0.1), x, z, u)
    avant = next(e for e in m.c.els if e.name == "Load5")
    apres = next(e for e in m2.c.els if e.name == "Load5")
    Z0 = complex(avant.pdict["R"], avant.pdict["X"])
    Z1 = complex(apres.pdict["R"], apres.pdict["X"])
    assert Z1 == pytest.approx(Z0 / 1.1)
    assert avant.pdict["R"] == pytest.approx(Z0.real)  # l'original intact


def test_saut_de_phase():
    m = ModeleNonLineaire(cat.construire("wscc9"))
    x = m.c.x0()
    _, x2, *_ = m.appliquer(
        Manoeuvre(0.0, "saut_de_phase", "bus5", angle=10.0), x, m.c.z0(), m.c.u0()
    )
    i, j = m.etats.index("bus5.vd"), m.etats.index("bus5.vq")
    a0, a1 = np.angle(complex(x[i], x[j])), np.angle(complex(x2[i], x2[j]))
    assert np.degrees(a1 - a0) == pytest.approx(10.0)
    assert abs(complex(x2[i], x2[j])) == pytest.approx(abs(complex(x[i], x[j])))


def test_api_temporel_avec_manoeuvre():
    d = analyse.temporel(
        cat.charger("wscc9"),
        "SM1.p_ref",
        duree=0.2,
        n=41,
        evenements=[{"instant": 0.05, "genre": "ouverture", "cible": "SM2"}],
        mesures=["P_SM2", "V_bus5"],
    )
    assert d["voies"] == ["non_lineaire"]
    assert d["mesures"]["P_SM2"][-1] is None and d["mesures"]["P_SM2"][0] is not None
    assert None not in d["mesures"]["V_bus5"]
    json.dumps(d, allow_nan=False)
    s = analyse.signaux(cat.charger("wscc9"))["manoeuvres"]
    assert {"nom": "SM2", "genre": "unite"} in s["elements"] and 5 in s["bus"]


# ------------------------------------------------------------------ le continu
@pytest.fixture(scope="module")
def b4():
    return cat.charger("cigre_b4")


def _spectre_max(t):
    net = assembler(t)
    net.solve_equilibrium()
    return float(net.build_state_space().eig().real.max()), net


def test_cable_continu_est_ouvrable(b4):
    """Un cable continu s'ouvre comme une ligne : meme table d'OUVERTS, nom
    assemble ``CDC<ID>_<de>_<a>``, genre ``cable``."""
    o = {(x["table"], x["id"]): x for x in ouvrables(b4)}
    assert o[("T_LinesDC", 2)]["nom"] == "CDC2_202_211"
    assert o[("T_LinesDC", 2)]["genre"] == "cable"
    t = basculer(b4, "T_LinesDC", 2, True)
    en, rap = en_service(t)
    assert [r[0] for r in en.T_LinesDC] == [r[0] for r in b4.T_LinesDC if r[0] != 2]
    assert len(b4.T_LinesDC) == 13  # l'original n'a pas bouge
    assert len(rap.iles_dc) == 2
    lam, net = _spectre_max(t)
    assert "CDC2_202_211" not in {e.name for e in net.elements}
    assert lam < 1e-6


def test_cable_ouvert_qui_coupe_le_continu_donne_des_iles(b4):
    t = b4.copie(OUVERTS=[("T_LinesDC", 3), ("T_LinesDC", 4), ("T_LinesDC", 5)])
    _, rap = en_service(t)
    iles = sorted(rap.iles_dc, key=lambda i: i["bus"][0])
    assert [i["bus"] for i in iles][-1] == [213, 214]
    assert any("coupe en 3 iles" in m for m in rap.messages)
    assert all(i["references"] for i in iles)
    lam, _ = _spectre_max(t)
    assert lam < 1e-6


def test_ile_continue_sans_reference_refusee():
    """Sur mtdc3t, couper le cable 1 laisse {10, 12} (deux stations a
    puissance imposee) sans rien pour tenir la tension : refus clair, par
    ``en_service`` comme par le load flow."""
    t = cat.charger("mtdc3t").copie(OUVERTS=[("T_LinesDC", 1)])
    with pytest.raises(ValueError, match=r"ile continue \[10, 12\] sans reference"):
        en_service(t)
    with pytest.raises(ValueError, match="aucune reference de tension"):
        assembler(t.copie(OUVERTS=[], T_LinesDC=t.T_LinesDC[1:])).solve_equilibrium()


def test_ile_alternative_morte_retire_les_noeuds_dc_orphelins(b4):
    """Le transformateur 7 isole {21, 22, 122} sans formeur : l'ile meurt, et
    avec elle la station 122, dont le noeud continu 222 et ses deux cables
    (4 : 202-222, 8 : 222-231) n'ont plus de raison d'etre."""
    t = b4.copie(OUVERTS=[("T_Transfos", 7)])
    en, rap = en_service(t)
    assert 222 not in {int(r[0]) for r in en.T_BusDC}
    assert {int(r[0]) for r in en.T_LinesDC} == {1, 2, 3, 5, 6, 7, 9, 10, 11, 12, 13}
    assert not [r for r in en.T_ShuntsDC if int(r[1]) == 222]
    assert any("noeud continu 222 retire" in m for m in rap.messages)
    assert any("station GFL122" in m for m in rap.messages)
    lam, net = _spectre_max(t)
    assert 222 not in net.v_dc and lam < 1e-6


def test_noeud_de_type_un_sans_tenant_devient_passif(b4):
    """Declencher GFL112 (VdcQ) laisse son noeud 212 de type 1 sans tenant :
    il devient passif (type 2) au lieu d'un slack ideal que le modele
    dynamique n'a pas."""
    t = b4.copie(OUVERTS=[("T_GFL2s", 4)])
    en, rap = en_service(t)
    assert {int(r[0]): int(r[1]) for r in en.T_BusDC}[212] == 2
    assert any("noeud continu 212" in m and "type 2" in m for m in rap.messages)
    assert next(r for r in b4.T_BusDC if int(r[0]) == 212)[1] == 1
    lam, _ = _spectre_max(t)
    assert lam < 1e-6


# ----------------------------------- la perte d'une station : spectre instable
def test_perte_d_une_station_est_dite_et_le_spectre_est_physique():
    """hvdc2t sans la station PV (GFM3 n°1) : le spectre devient instable
    (+0,67 1/s vers 214 Hz). Diagnostic : ce n'est pas un artefact
    d'assemblage (pas de mode de jauge, pas de noeud flottant) mais le
    convertisseur restant (GFL, regulant Vdc) a puissance modifiee : le meme
    reseau INTACT, avec la puissance de la station forcee a zero, l'est aussi,
    et il ne l'est pas a la puissance nominale. Le rapport le dit."""
    t = cat.charger("hvdc2t")
    lam0, _ = _spectre_max(t)
    assert lam0 < 1e-6

    _, rap = en_service(t.copie(OUVERTS=[("T_GFM3s", 1)]))
    msg = [m for m in rap.messages if m.startswith("station GFM2")]
    assert msg and "GFL3 (VdcQ)" in msg[0] and "imposait sa puissance" in msg[0]
    assert "instable" in msg[0]
    lam, _ = _spectre_max(t.copie(OUVERTS=[("T_GFM3s", 1)]))
    assert lam > 0.1

    # meme station presente, mais a puissance nulle : meme famille
    y0 = [(*r[:5], 0.0, *r[6:]) if int(r[0]) == 2 else r for r in t.Y_DER]
    lam_p0, _ = _spectre_max(t.copie(Y_DER=y0))
    assert lam_p0 > 0.1


def test_perte_d_une_station_qui_tient_vdc_refusee():
    """Sans GFL2 (VdcQ), hvdc2t n'a plus rien pour tenir sa tension continue."""
    t = cat.charger("hvdc2t").copie(OUVERTS=[("T_GFL2s", 1)])
    with pytest.raises(ValueError, match="sans reference de tension"):
        en_service(t)
