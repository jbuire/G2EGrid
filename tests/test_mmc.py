"""Le modele moyen du MMC (:mod:`g2egrid.elements.mmc`) : equilibre exact,
pertes, stabilite sur le CIGRE B4, blocage, et ce que la regulation
d'energie fait aux statismes raides."""

import warnings

import numpy as np
import pytest
from scipy.optimize import root

from g2egrid.catalogue import charger, construire
from g2egrid.domaines import domaine_etat
from g2egrid.element import InitContext
from g2egrid.elements import MMC, MMCEnergie
from g2egrid.robustesse import _spectre
from g2egrid.tables import assembler
from g2egrid.temporel import Manoeuvre, ModeleNonLineaire, protections_dc


def _en_mmc(t, modele="mmc"):
    return t.copie(GFL_PARAMS=dict(t.GFL_PARAMS) | {"modele": modele})


@pytest.fixture(scope="module")
def dcs1():
    net = assembler(_en_mmc(charger("cigre_b4_dcs1")))
    net.run_loadflow()
    return net, ModeleNonLineaire(net)


def _v(tr, nom, t):
    j = tr.etats.index(nom)
    return float(tr.x[min(np.searchsorted(tr.t, t), len(tr.t) - 1), j])


@pytest.mark.parametrize("classe", [MMC, MMCEnergie])
def test_bilan_de_puissance(classe):
    """A l'equilibre, la station prend au continu ce qu'elle rend a
    l'alternatif plus ses pertes cuivre : ``R1.|i|^2 + Rdc.idc^2``."""
    el = classe("M", 1, 2, Sn=800e6, Sb=100e6, wb=2 * np.pi * 50)
    S = complex(4.0, 1.5)  # pu systeme, injecte au bus alternatif
    ctx = InitContext(v={"t": 1.02 * np.exp(0.3j), "dc": 1.01}, s={"t": S})
    x = el.initialize(ctx)
    kS, p = el.pdict["kS"], el.pdict
    i2 = x["isd"] ** 2 + x["isq"] ** 2
    pertes = p["R1_pu"] * i2 + p["Rdc_pu"] * x["idc"] ** 2
    assert 1.01 * x["idc"] == pytest.approx(S.real / kS + pertes, rel=1e-12)
    assert x["theta_pll"] == pytest.approx(0.3)
    # energie libre : W suit le carre de la tension des bras
    v_z = 1.01 - p["Rdc_pu"] * x["idc"]
    attendu = p["W0"] if classe.regule else p["W0"] * v_z**2
    assert x["W"] == pytest.approx(attendu)


def test_assemble_depuis_le_reglage():
    t = charger("cigre_b4_dcs1")
    noms = {type(e).__name__ for e in assembler(_en_mmc(t)).elements}
    assert "MMC" in noms and "GFL2DC" not in noms
    noms = {type(e).__name__ for e in assembler(_en_mmc(t, "mmc_energie")).elements}
    assert "MMCEnergie" in noms
    with pytest.raises(ValueError, match="modele 'lcc' inconnu"):
        assembler(_en_mmc(t, "lcc"))


def test_mmc_sans_port_continu_refuse():
    t = charger("cigre_b4")
    t = t.copie(PARAMS_PAR_BUS=dict(t.PARAMS_PAR_BUS) | {21: {"modele": "mmc"}})
    with pytest.raises(ValueError, match="forcement un port continu"):
        assembler(t)


def test_dcs1_point_fixe_et_load_flow(dcs1):
    """Le load flow mixte prend ses pertes au modele : le point est un point
    fixe du modele non lineaire, et la tension de Bm-C1 reste celle du
    modele de reference (1,0108 pu)."""
    net, m = dcs1
    dx, dz = m.residu_equilibre()
    assert dx < 1e-6 and dz < 1e-10
    v = dict(zip(sorted(net.bus_specs_dc), net.lf.Vdc))
    assert v[221] == pytest.approx(1.0108, abs=3e-4)


@pytest.mark.parametrize("modele", ["mmc", "mmc_energie"])
def test_dcs1_stable(modele):
    sp = _spectre(_en_mmc(charger("cigre_b4_dcs1"), modele))
    assert float(sp.lam.real.max()) < 0


def test_cigre_complet_en_mmc_stable():
    """La variante du catalogue : stations terrestres en MMC a energie libre."""
    net = construire("cigre_b4_mmc")
    assert sum(type(e).__name__ == "MMC" for e in net.elements) == 7
    sp = _spectre(charger("cigre_b4_mmc"))
    assert float(sp.lam.real.max()) < 0


def test_energie_regulee_destabilise_les_statismes():
    """Regulee en energie, la station cache ses sous-modules au reseau
    continu : les noeuds des statismes raides (Bb-B1, Bb-B2, Bm-B3) n'ont
    plus que la capacite des lignes, et leur boucle part instable."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sp = _spectre(_en_mmc(charger("cigre_b4"), "mmc_energie"))
    instables = np.flatnonzero(sp.lam.real > 0)
    assert len(instables)
    dominants = {sp.etats[int(np.argmax(sp.p[k]))] for k in instables}
    assert dominants & {"bus211.v", "bus213.v", "bus215.v"}, dominants


def test_domaine_des_etats():
    net = construire("cigre_b4_mmc")
    el = next(e for e in net.elements if type(e).__name__ == "MMC")
    assert domaine_etat(el, "idc") == "DC"
    assert domaine_etat(el, "W") == "DC"
    assert domaine_etat(el, "isd") == "AC"


def test_defaut_franc_bloque_et_eteint(dcs1):
    """Defaut franc au noeud C1 : Cm-C1 se bloque, ses bras eteignent le
    courant continu en quelques millisecondes ; l'energie des sous-modules
    de la station bloquee ne bouge plus."""
    net, m = dcs1
    sc = m.scenario(Manoeuvre(0.02, "defaut_dc", "bus221"))
    tr = m.simuler((0, 0.05), uext_fn=sc, methode="BDF2", pas=2.5e-5,
                   protections=protections_dc(net))  # fmt: skip
    assert tr.declenchements[0][1] == "sous_tension_dc:bus221"
    assert tr.declenchements[0][2] == [("blocage", "GFL121")]
    t_b = tr.declenchements[0][0]
    assert abs(_v(tr, "GFL121.idc", t_b + 0.01)) < 1e-3
    w = [_v(tr, "GFL121.W", t_b + d) for d in (0.01, 0.025)]
    assert w[0] == pytest.approx(w[1], rel=1e-3)


def test_defaut_resistif_tenu(dcs1):
    """Defaut de 0,5 pu elimine en 10 ms : aucune station ne se bloque,
    l'energie des sous-modules reprend sa valeur."""
    net, m = dcs1
    sc = m.scenario(
        Manoeuvre(0.02, "defaut_dc", "bus221", R=0.5),
        Manoeuvre(0.03, "elimination_dc", "bus221"),
    )
    tr = m.simuler((0, 0.3), uext_fn=sc, methode="BDF2", pas=5e-5,
                   protections=protections_dc(net))  # fmt: skip
    assert not tr.declenchements
    w0 = _v(tr, "GFL121.W", 0.0)
    assert _v(tr, "GFL121.W", 0.299) == pytest.approx(w0, rel=2e-2)


@pytest.mark.parametrize("modele", ["mmc", "mmc_energie"])
def test_jacobien_contre_differences_finies(modele):
    """Le jacobien symbolique de ``MMC`` et de ``MMCEnergie`` colle aux
    differences finies, au point d'equilibre du CIGRE B4 a un MMC.

    La derivee de ``dx/dt`` est prise a ``z`` reresolu (les algebriques suivent
    l'etat) : c'est la matrice ``A`` que l'analyse modale lit."""
    classe = {"mmc": MMC, "mmc_energie": MMCEnergie}[modele]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        net = assembler(_en_mmc(charger("cigre_b4_dcs1"), modele))
        net.solve_equilibrium()
    stations = [e for e in net.all_elements if type(e) is classe]
    assert stations

    for el in stations[:2]:
        cm = el.compiled()
        A, _, _, _ = cm.linearize(el.x0, el.z0, el.u0, el.p)

        def dxdt(x, _el=el, _cm=cm):
            z = _el.z0
            if _cm.n_z:
                z = root(lambda zz: _cm.g(x, zz, _el.u0, _el.p), _el.z0, tol=1e-13).x
            return _cm.f(x, z, _el.u0, _el.p)

        h = 1e-7
        An = np.zeros((cm.n_x, cm.n_x))
        for j in range(cm.n_x):
            xp, xm = el.x0.copy(), el.x0.copy()
            xp[j] += h
            xm[j] -= h
            An[:, j] = (dxdt(xp) - dxdt(xm)) / (2 * h)
        err = np.max(np.abs(A - An)) / max(np.max(np.abs(A)), 1.0)
        assert err < 1e-6, f"{el.name}: erreur jacobien {err:.2e}"


def test_kvdc_different_de_un_refuse():
    """``kVdc`` ne pese que le courant rendu au noeud, pas la tension des
    boucles : un MMC construit a la main avec ``kVdc != 1`` serait faux. Il est
    refuse (comme l'assemblage le fait deja par ``_verifier_bases_dc``)."""
    with pytest.raises(ValueError, match="kVdc"):
        MMC("M", 1, 2, Sn=800e6, Sb=100e6, wb=2 * np.pi * 50, kVdc=0.5)
    MMC("M", 1, 2, Sn=800e6, Sb=100e6, wb=2 * np.pi * 50, kVdc=1.0)
