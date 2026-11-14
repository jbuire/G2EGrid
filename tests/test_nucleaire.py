"""Groupes nucleaires : lignes d'arbre publiees, catalogue des machines et
reseaux des theses de Boudot et de Mazauric.

References :

* rapport EPRI 1013460 (2006), annexe A : inerties [kg.m2] et raideurs [N.m/rad]
  du groupe nucleaire de 1500 MVA ET du First Benchmark Model. Le second est
  aussi publie en pu (IEEE, 1977) : la conversion physique -> pu se verifie
  donc sur une donnee independante ;
* Mudau et Field (KNS 2018), tableau 3 : la ligne a trois masses de l'AM600 ;
* these de Boudot, tableaux IV.1 et A.1 a A.5 ; these de Mazauric, tableaux
  III-4 et III-9.
"""

import copy

import numpy as np
import pytest
from g2egrid.catalogue import charger
from g2egrid.tables import assembler

from g2egrid.elements import (
    ARBRE_AM600,
    ARBRE_FBM,
    ARBRE_NUCLEAIRE_EPRI,
    LIGNES_ARBRE,
    LigneArbre,
    ligne_a_l_echelle,
    ligne_depuis_physique,
)
from g2egrid.reseau.machines import (
    BUS_NE39,
    MACHINES,
    fiche_machine,
    lister_machines,
    retard_moyen_ieeeg1,
    retard_moyen_ieeeg3,
)
from g2egrid.reseau.reglage import standard_depuis_machine

WB50, WB60 = 2.0 * np.pi * 50.0, 2.0 * np.pi * 60.0


def _resolu(t):
    net = assembler(t)
    net.run_loadflow()
    net.solve_equilibrium()
    return net


def _modes(net):
    return np.linalg.eigvals(net.build_state_space().A)


# ----------------------------------------------------------------------
# lignes d'arbre
# ----------------------------------------------------------------------
def test_conversion_physique_retrouve_le_fbm_en_pu():
    """Les inerties et raideurs PHYSIQUES du FBM (EPRI 1013460, annexe A)
    redonnent ses valeurs en pu de 1977 : la conversion est la bonne."""
    lg = ligne_depuis_physique(
        masses=ARBRE_FBM.masses,
        generateur="GEN",
        J=(1166.62, 1953.92, 10783.35, 11104.15, 10906.73, 429.698),
        K=(45693376.74, 82682689.53, 123182507.3, 167732543.6, 6675403.95),
        F=ARBRE_FBM.F,
        Sn_VA=892.4e6,
        vitesse_tr_min=3600.0,
        f0=60.0,
    )
    assert lg.H == pytest.approx(ARBRE_FBM.H, rel=2e-5)
    assert lg.K[:4] == pytest.approx(ARBRE_FBM.K[:4], rel=2e-5)
    # troncon de l'excitatrice : le rapport part de 2,82 pu, la publication de
    # 1977 donne 2,822 -- un arrondi de la source, pas de la conversion
    assert lg.K[4] == pytest.approx(2.82, rel=1e-6)
    assert lg.K[4] == pytest.approx(ARBRE_FBM.K[4], rel=1e-3)
    assert lg.D == pytest.approx(ARBRE_FBM.D, rel=2e-5)


def test_ligne_nucleaire_epri():
    """1500 MVA, 1800 tr/min : 316 584 kg.m2 au total, soit H = 3,75 s ; huit
    modes de torsion, dont quatre sous 25 Hz (la ligne entiere)."""
    lg = ARBRE_NUCLEAIRE_EPRI
    assert len(lg.masses) == 9 and lg.generateur == "GEN"
    wm = 2.0 * np.pi * 30.0
    assert lg.H_total == pytest.approx(316583.89 * wm**2 / (2 * 1500e6), rel=1e-9)
    assert lg.H_total == pytest.approx(3.7495, abs=1e-4)
    f = lg.frequences_libres(WB60)
    assert f[0] == pytest.approx(0.0, abs=1e-5)
    assert f[1:] == pytest.approx(
        [8.15, 15.83, 22.05, 23.71, 50.64, 52.53, 54.48, 133.93], abs=6e-3
    )
    # puissance repartie egalement entre les huit demi-corps, comme le rapport
    assert lg.F == pytest.approx((0.125,) * 8 + (0.0,))
    # les quatre arbres physiques : entre deux rotors
    assert lg.arbres == ("HPB_LP1A", "LP1B_LP2A", "LP2B_LP3A", "LP3B_GEN")


def test_ligne_am600_et_son_ecart_a_l_article():
    """Le tableau 3 de l'article donne 18,87 et 28,80 Hz. L'article annonce 18,2
    et 23,9 Hz : le premier mode est retrouve a 4 %, le second ne l'est PAS
    (20 %). L'ecart est dans la source ; ce test le tient en vue."""
    f = ARBRE_AM600.frequences_libres(WB50)
    assert f[1:] == pytest.approx([18.87, 28.80], abs=6e-3)
    assert f[1] == pytest.approx(18.2, rel=0.04)
    assert abs(f[2] / 23.9 - 1.0) > 0.15
    assert ARBRE_AM600.H_total == pytest.approx(4.70, abs=5e-3)


@pytest.mark.parametrize("nom", sorted(LIGNES_ARBRE))
def test_mise_a_l_echelle_garde_les_frequences(nom):
    """Meme deformee, memes frequences en hertz, inertie imposee -- y compris
    en changeant la frequence du reseau."""
    ligne, f0 = LIGNES_ARBRE[nom]
    ref = ligne.frequences_libres(2 * np.pi * f0)
    for f_cible in (50.0, 60.0):
        lg = ligne_a_l_echelle(ligne, 6.5, rapport_f0=f0 / f_cible)
        assert lg.H_total == pytest.approx(6.5, rel=1e-12)
        assert lg.frequences_libres(2 * np.pi * f_cible) == pytest.approx(ref, abs=1e-6)
        # amortissement : toujours proportionnel a l'inertie
        assert np.array(lg.D) / np.array(lg.H) == pytest.approx(
            np.array(ligne.D) / np.array(ligne.H)
        )
    with pytest.raises(ValueError):
        ligne_a_l_echelle(ligne, 0.0)


def test_accouplement_inconnu_refuse():
    with pytest.raises(ValueError, match="accouplements inconnus"):
        LigneArbre(
            masses=("T", "G"),
            generateur="G",
            H=(1.0, 1.0),
            K=(10.0,),
            D=(0.0, 0.0),
            Dm=(0.0,),
            F=(1.0, 0.0),
            accouplements=("G_T",),
        )


# ----------------------------------------------------------------------
# catalogue des machines
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nom", sorted(MACHINES))
def test_fiche_convertie_et_relue(nom):
    """Chaque fiche se convertit en parametres de Park, et la conversion
    inverse redonne la fiche."""
    f = fiche_machine(nom)
    f0 = f.f0 or 50.0
    park = f.park(f0)
    assert all(v > 0 for k, v in park.items() if k != "Ra")
    relu = standard_depuis_machine(park, f0=f0)
    for k, v in f.standard.items():
        assert relu[k] == pytest.approx(v, rel=1e-9, abs=1e-12), k


def test_smr_de_la_these_de_boudot():
    """Tableau IV.1 : six turbo-alternateurs de 25 a 245,5 MVA."""
    smr = [f for f in lister_machines("nucleaire") if f.nom.startswith("smr_")]
    assert [f.Sn_VA / 1e6 for f in smr] == [25.0, 35.29, 51.2, 75.0, 100.0, 245.5]
    assert [f.H for f in smr] == [5.016, 4.389, 5.078, 6.187, 4.985, 4.624]
    f = fiche_machine("smr_245")
    assert (f.standard["Xd"], f.standard["Xd_pp"], f.standard["Td0_p"]) == (
        1.71,
        0.21,
        7.1,
    )
    # la source ne fixe pas la frequence : il faut la donner
    with pytest.raises(ValueError, match="preciser f0"):
        f.park()
    # la ligne pretee depend de la frequence du reseau
    assert f.reglages(f0=50.0, arbre=True)["arbre"] == "am600"
    assert f.reglages(f0=60.0, arbre=True)["arbre"] == "nucleaire_epri"


def test_machines_du_39_noeuds():
    """Tableau A.1 et figure A-1 : cinq nucleaires, trois thermiques, une
    hydraulique, l'interconnexion."""
    filieres = {b: fiche_machine(n).filiere for b, n in BUS_NE39.items()}
    assert sorted(b for b, f in filieres.items() if f == "nucleaire") == [
        31,
        32,
        35,
        37,
        38,
    ]
    assert sorted(b for b, f in filieres.items() if f == "thermique") == [33, 34, 36]
    assert filieres[30] == "hydraulique" and filieres[39] == "interconnexion"
    assert fiche_machine("ne39_g1").Sn_VA == 10000e6
    assert fiche_machine("ne39_g9").H == 3.45
    with pytest.raises(KeyError, match="fiches"):
        fiche_machine("inconnue")


def test_reductions_des_regulateurs():
    """IEEEG1 du 39 noeuds : 0,2 - 1 + 0,6 + (0,18 + 0,275 + 0,57 + 0,435) =
    1,26 s. IEEEG3 : (0,05 + 0,2.10)/0,04 + 1,5.0,75 = 52,375 s."""
    tg = retard_moyen_ieeeg1(
        0.2, 1.0, 0.6, ((0.3, 0.6), (0.25, 0.5), (0.3, 0.8), (0.15, 1.0))
    )
    assert tg == pytest.approx(1.26, abs=1e-12)
    assert fiche_machine("ne39_g2").regulation["Tg"] == pytest.approx(1.26)
    assert retard_moyen_ieeeg3(0.05, 0.04, 0.2, 10.0, 0.75, 0.5, 1.0, 1.5, 1.0) == (
        pytest.approx(52.375)
    )


def test_avance_retard_du_st1_par_la_boucle_de_stabilisation():
    """``Ka.(1 + s.Tc)/(1 + s.Tb)`` (ST1 : Tc = 1, Tb = 10) est, terme a terme,
    la fonction de transfert de l'amplificateur boucle par ``s.Kfd/(1 + s.Tfd)``
    avec ``Tfd = Tc`` et ``Kfd = (Tb - Tc)/Ka``."""
    reg = fiche_machine("ne39_g4").regulation
    Ka, Kfd, Tfd = reg["Ka"], reg["Kfd"], reg["Tfd"]
    for w in (0.01, 0.3, 2.0, 40.0):
        s = 1j * w
        boucle = Ka / (1.0 + Ka * s * Kfd / (1.0 + s * Tfd))
        st1 = Ka * (1.0 + s * 1.0) / (1.0 + s * 10.0)
        assert boucle == pytest.approx(st1, rel=1e-12)


# ----------------------------------------------------------------------
# la cle "arbre" des tables
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def smib_am600():
    t = copy.deepcopy(charger("smib"))
    t.SM_PARAMS = dict(t.SM_PARAMS) | {
        "arbre": "am600",
        "H": 6.5,
        "K_PSS": 0.0,
        "Tg": 0.5,
    }
    return _resolu(t)


def test_cle_arbre_donne_une_machine_multi_masses(smib_am600):
    sm = next(e for e in smib_am600.elements if e.name.startswith("SM"))
    assert type(sm).__name__ == "MultiMasses_HP_LP_GEN"
    p = sm.pdict
    # le H de la machine est l'inertie de TOUTE la ligne
    assert p["H"] + p["H_HP"] + p["H_LP"] == pytest.approx(6.5, rel=1e-12)
    # et les frequences de l'AM600 sont conservees
    assert sm.ligne.frequences_libres(WB50)[1:] == pytest.approx(
        [18.87, 28.80], abs=6e-3
    )


def test_cle_arbre_modes_et_couples_d_equilibre(smib_am600):
    """Deux modes de torsion faiblement amortis pres de 18,9 et 28,8 Hz ; a
    l'equilibre chaque troncon transmet le couple des etages en amont."""
    lam = _modes(smib_am600)
    assert lam.real.max() < 1e-6
    tors = sorted(x.imag / 2 / np.pi for x in lam if x.imag > 60 and abs(x.real) < 1)
    assert tors == pytest.approx([18.88, 28.81], abs=0.05)
    sm = next(e for e in smib_am600.elements if e.name.startswith("SM"))
    p, x0 = sm.pdict, dict(zip(sm.states, sm.x0))
    Pm = x0["Pm"]  # a vitesse nominale : couple = puissance
    assert p["K_HP_LP"] * x0["d_HP_LP"] == pytest.approx(0.30 * Pm, rel=1e-9)
    assert p["K_LP_GEN"] * x0["d_LP_GEN"] == pytest.approx(1.00 * Pm, rel=1e-9)


def test_cle_arbre_erreurs():
    t = copy.deepcopy(charger("smib"))
    t.SM_PARAMS = dict(t.SM_PARAMS) | {"arbre": "inconnue"}
    with pytest.raises(ValueError, match="ligne d'arbre 'inconnue' inconnue"):
        assembler(t)
    t.SM_PARAMS["arbre"] = "am600"
    t.T_SMs = [(*r[:3], 0.002, 0.15) for r in t.T_SMs]
    with pytest.raises(ValueError, match="transformateur integre"):
        assembler(t)


# ----------------------------------------------------------------------
# reseaux du catalogue
# ----------------------------------------------------------------------
def _amortissement_min(lam, fmin=0.1, fmax=3.0):
    em = [x for x in lam if 2 * np.pi * fmin < x.imag < 2 * np.pi * fmax]
    return min(-x.real / abs(x) for x in em)


@pytest.fixture(scope="module")
def tables_39():
    return charger("ieee39_boudot")


def test_ieee39_boudot_donnees(tables_39):
    """Puissances nominales du tableau A.1, regulations des tableaux A.2 a A.5."""
    t = tables_39
    assert t.f0 == 60.0
    Sn = {int(b): s for _i, b, s, *_ in t.T_SMs}
    assert Sn == {int(b): fiche_machine(n).Sn_VA for b, n in BUS_NE39.items()}
    for b, nom in BUS_NE39.items():
        reg = t.PARAMS_PAR_BUS[b]
        f = fiche_machine(nom)
        assert reg["H"] == f.H
        assert reg["Ka"] == 200.0 and reg["Kfd"] == pytest.approx(0.045)
        assert reg["K_PSS"] == f.regulation["K_PSS"]


def test_ieee39_boudot_stable_et_amorti(tables_39):
    """Avec ses stabilisateurs, le benchmark est amorti a plus de 10 % sur tous
    ses modes electromecaniques ; c'est aussi le controle que l'avance-retard
    de l'AVR est bien traduite (sans elle, le reseau est instable)."""
    lam = _modes(_resolu(copy.deepcopy(tables_39)))
    assert lam.real.max() < 1e-6
    assert _amortissement_min(lam) > 0.10


def test_un_smr_remplace_une_centrale(tables_39):
    """La "modularisation" de la these : la centrale du noeud 37 (700 MVA)
    devient trois SMR de 245,5 MVA, agreges en une machine."""
    t = copy.deepcopy(tables_39)
    smr = fiche_machine("smr_245")
    t.T_SMs = [(i, b, 3 * smr.Sn_VA if int(b) == 37 else s) for i, b, s in t.T_SMs]
    t.PARAMS_PAR_BUS[37] = t.PARAMS_PAR_BUS[37] | smr.reglages(f0=60.0)
    lam = _modes(_resolu(t))
    assert lam.real.max() < 1e-6
    assert _amortissement_min(lam) > 0.08


def test_kundur_mazauric():
    """50 Hz ; plan de production du tableau III-9 ; groupes nucleaires a 4 % de
    statisme et 10 s ; un seul stabilisateur, sur GEN_B1."""
    t = charger("kundur_mazauric")
    assert t.f0 == 50.0
    Pg = {int(r[0]): r[5] for r in t.Y_DER}
    assert (Pg[2], Pg[3], Pg[4], Pg[7]) == (700e6, 720e6, 700e6, 500e6)
    for b in (2, 4):
        assert (t.PARAMS_PAR_BUS[b]["mp"], t.PARAMS_PAR_BUS[b]["Tg"]) == (4.0, 10.0)
    assert [b for b in (1, 2, 3, 4) if t.PARAMS_PAR_BUS[b]["K_PSS"] > 0] == [3]
    assert (t.PARAMS_PAR_BUS[1]["H"], t.PARAMS_PAR_BUS[3]["H"]) == (6.5, 6.175)
    net = _resolu(copy.deepcopy(t))
    lam = _modes(net)
    assert lam.real.max() < 1e-6
    # la machine de reference fournit le reste : 174 MW dans la these
    sm1 = next(e for e in net.elements if e.name == "SM1")
    P1 = dict(zip(sm1.states, sm1.x0))["Pm"] * 900.0
    assert P1 == pytest.approx(174.0, abs=20.0)


def test_cle_arbre_aller_retour_json(smib_am600):
    """La cle ``arbre`` est une chaine parmi des reglages numeriques : elle
    survit a l'ecriture et a la relecture d'un reseau."""
    from g2egrid.tables import tables_depuis_json, tables_vers_json

    t = copy.deepcopy(charger("smib"))
    t.PARAMS_PAR_BUS = {1: {"arbre": "am600", "H": 6.5}}
    relu = tables_depuis_json(tables_vers_json(t))
    assert relu.PARAMS_PAR_BUS[1] == {"arbre": "am600", "H": 6.5}
    sm = next(e for e in assembler(relu).all_elements if e.name == "SM1")
    assert sm.troncons == ("HP_LP", "LP_GEN")
