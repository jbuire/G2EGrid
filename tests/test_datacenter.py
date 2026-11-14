"""Tests des blocs d'un centre de donnees et du reseau type.

Quatre familles de verifications, chacune contre une reference independante :

* le BILAN DE PUISSANCE du redresseur (une identite algebrique du pont de
  diodes ideal : ce qui entre cote alternatif sort cote continu, aux pertes
  ``Rdc.idc^2`` pres) ;
* les JACOBIENS des nouveaux blocs contre les differences finies ;
* la COHERENCE du modele linearise avec le modele non lineaire, sur un petit
  echelon : l'ecart doit etre du second ordre ;
* la physique du GROUPE : le mode de torsion tombe a ``f_tors`` et la butee de
  carburant plafonne le couple moteur.
"""

import warnings

import numpy as np
import pytest
from scipy.optimize import root

from g2egrid import Network
from g2egrid.datacenter import (
    NOMS,
    courbe_charge_ia,
    reseau_datacenter,
    scenario_perte_reseau,
)
from g2egrid.elements import (
    IB,
    BatterieDCDC,
    ChargeACPref,
    GroupeElectrogene,
    Line,
    OnduleurVF,
    RedresseurDiode,
)
from g2egrid.elements.datacenter import reactif_commutation
from g2egrid.fatigue import evaluer_fatigue
from g2egrid.temporel import (
    Courbe,
    Evenement,
    ModeleNonLineaire,
    simuler_lineaire,
)


@pytest.fixture(scope="module")
def net():
    return reseau_datacenter()


@pytest.fixture(scope="module")
def modele(net):
    return ModeleNonLineaire(net)


def element(net, nom):
    return next(e for e in net.elements if e.name == nom)


def sortie(net, nom, sortie_):
    """Une sortie de service d'un element, evaluee a l'equilibre."""
    el = element(net, nom)
    y = el.compiled().h(el.x0, el.z0, el.u0, el.p)
    return float(y[el.n_port_signaux() + el.svc_outputs.index(sortie_)])


# ----------------------------------------------------------------------
# le reseau type
# ----------------------------------------------------------------------
def test_reseau_type_equilibre_et_stabilite(net, modele):
    """Le load flow mixte converge, l'equilibre est un point fixe du modele non
    lineaire (``|dx/dt|`` ~ 1e-13) et le modele linearise est stable."""
    assert net.lf.converged
    dxdt, ecart = modele.residu_equilibre()
    assert dxdt < 1e-9 and ecart < 1e-9
    ss = net.build_state_space()
    assert ss.is_stable()
    # les deux seuls modes nuls : l'etat de charge et l'angle de l'onduleur
    zeros = ss.eig()[np.abs(ss.eig()) < 1e-6]
    assert len(zeros) == 3  # soc, theta de l'onduleur, theta du bus infini


def test_reseau_type_puissances(net):
    """La charge appelle 0,30 pu ; l'onduleur la fournit (pertes comprises), le
    reseau amont fournit le reste, le groupe flotte a vide."""
    assert sortie(net, "IT", "p") == pytest.approx(0.30, abs=1e-9)
    assert sortie(net, "OND", "P") * 0.4 == pytest.approx(0.30, rel=5e-3)  # Sn/Sb = 0,4
    p_red = sortie(net, "RED", "p")
    assert p_red > 0.30  # pertes de l'onduleur + du pont
    assert p_red < 0.33
    assert sortie(net, "GE", "P") * 0.5 < 0.025  # 0,02 pu systeme en attente
    # la batterie est en veille : ni courant ni derive de l'etat de charge
    assert sortie(net, "BAT", "p") == 0.0
    assert sortie(net, "BAT", "p_bat") == 0.0


# ----------------------------------------------------------------------
# redresseur
# ----------------------------------------------------------------------
def _pont_physique(Xcom, U=400.0, S=2e6, f=50.0):
    """Chute de commutation et inductance continue d'un pont triphase, calculees
    en GRANDEURS PHYSIQUES puis ramenees en pu de la base continue (tension
    redressee a vide, meme puissance) : la reference independante du code.
    ``Xcom`` en pu de la base alternative du pont (``U``, ``S``)."""
    Vd0 = 3.0 * np.sqrt(2.0) / np.pi * U  # V
    Zdc = Vd0**2 / S  # ohm, base continue
    Xc = Xcom * U**2 / S  # ohm, reactance de commutation par phase
    Rcom = 3.0 * Xc / np.pi  # ohm, chute moyenne : 3.w.Lc.Id/pi
    L2 = 2.0 * Xc  # ohm : deux phases en serie cote continu, hors empietement
    return Rcom / Zdc, L2 / Zdc


def _tan_phi_fourier(mu, n=400_000):
    """``tan(phi1)`` du fondamental du courant de phase d'un pont a diodes
    (angle de retard nul, empietement ``mu``), par serie de Fourier NUMERIQUE de
    la forme d'onde : montee ``(1 - cos(th + pi/3))/(1 - cos(mu))`` pendant
    l'empietement, palier, descente symetrique ; alternance negative opposee.
    Tension de phase ``cos(th)``."""
    th = np.linspace(-np.pi, np.pi, n, endpoint=False)

    def positive(x):
        i = np.zeros_like(x)
        a, b = -np.pi / 3, np.pi / 3
        m = (x >= a) & (x < a + mu)
        i[m] = (1 - np.cos(x[m] - a)) / (1 - np.cos(mu))
        i[(x >= a + mu) & (x < b)] = 1.0
        m = (x >= b) & (x < b + mu)
        i[m] = 1 - (1 - np.cos(x[m] - b)) / (1 - np.cos(mu))
        return i

    i = positive(th) - positive(np.mod(th, 2 * np.pi) - np.pi)
    a1, b1 = np.mean(i * np.cos(th)), np.mean(i * np.sin(th))
    return b1 / a1  # i ~ A.cos(th - phi) : retard phi


def test_redresseur_bilan_de_puissance(net):
    """Pont ideal : ``P_ac = (Vr - Rcom.idc).idc`` et ``P_dc = v_dc.idc``, la
    chute ``Rcom`` et l'inductance continue etant celles du pont PHYSIQUE
    (``3.Xc/pi`` ohms, ``2.Lc``) ramenees en pu. La tension continue est
    ``Vr - Rd.idc`` au terme de lissage ``eps^2/(4 Rd idc)`` pres."""
    red = element(net, "RED")
    p = red.pdict
    kS = p["kS"]
    idc = red.x0[0]
    vdc = net.v_dc[10]
    vr = sortie(net, "RED", "vr")
    rcom, l2 = _pont_physique(p["Xcom"])
    rd = p["Rdc"] + rcom
    p_ac = sortie(net, "RED", "p") / kS
    assert p_ac == pytest.approx((vr - rcom * idc) * idc, rel=1e-12)
    biais = p["eps_v"] ** 2 / (4.0 * rd * idc)
    assert vdc == pytest.approx(vr - rd * idc + biais, rel=1e-9)
    assert p_ac - vdc * idc == pytest.approx(p["Rdc"] * idc**2 - biais * idc, rel=1e-6)
    # le biais de lissage est negligeable devant la chute ohmique (< 0,2 %)
    assert biais < 2e-3 * rd * idc
    # constante de temps du courant continu : d(didc/dt)/d(idc) = -wb.Rd/Lt
    A = red.compiled().linearize(red.x0, red.z0, red.u0, red.p)[0]
    assert A[0, 0] == pytest.approx(-p["wb"] * rd / (p["Ld"] + l2), rel=1e-9)


def test_redresseur_calage_kv_tphi(net):
    """``kv`` est cale pour que le pont delivre la tension continue du bus ; le
    reactif du load flow fixe ``tphi``, et le constructeur du reseau impose au
    load flow le reactif PHYSIQUE de commutation, verifie ici par la serie de
    Fourier du courant de phase (empietement ``mu`` tire de ``Vd/Vd0``)."""
    red = element(net, "RED")
    p = red.pdict
    V = abs(net.lf.V[net._bus_index()[4]])
    vr = sortie(net, "RED", "vr")
    assert p["kv"] * V == pytest.approx(vr, rel=1e-12)
    idc = red.x0[0]
    rcom, _ = _pont_physique(p["Xcom"])
    mu = np.arccos(2.0 * (vr - rcom * idc) / vr - 1.0)
    tphi = _tan_phi_fourier(mu)
    assert p["tphi"] == pytest.approx(tphi, rel=1e-6)
    assert reactif_commutation(red) == pytest.approx(sortie(net, "RED", "q"), rel=1e-6)
    # ordre de grandeur : 0,145 pour Xcom = 0,03 pu a pleine charge (mu ~ 12 deg)
    assert 0.10 < p["tphi"] < 0.20


def test_redresseur_rejette_les_parametres_non_physiques():
    with pytest.raises(ValueError, match="Rdc"):
        RedresseurDiode("R", 1, 2, Sn=1e6, Sb=1e6, wb=314.0, Rdc=0.0)
    with pytest.raises(ValueError, match="Ld"):
        RedresseurDiode("R", 1, 2, Sn=1e6, Sb=1e6, wb=314.0, Ld=-1.0)


def test_redresseur_se_bloque_quand_la_tension_continue_depasse_vr(net, modele):
    """Si la batterie releve le bus continu au-dessus de ``Vr``, les diodes se
    bloquent : ``idc`` tombe presque a zero (le terme de lissage laisse passer
    moins de 0,5 % du nominal) et la batterie porte la charge."""
    sc = modele.scenario(Evenement(0.2, "BAT.v_ref", 0.14, "rampe", 0.02))
    t = np.linspace(0.0, 1.5, 151)
    tr = modele.simuler((0, 1.5), uext_fn=sc, t_eval=t, rtol=1e-5, atol=1e-7)
    assert tr.etat("RED.idc")[0] > 0.7
    assert abs(tr.etat("RED.idc")[-1]) < 5e-3
    # la batterie debite la puissance que le pont ne fournit plus : sa
    # puissance DC vaut la charge amont (onduleur + pertes), 0,30 a 0,32 pu
    i_bat = tr.etat("BAT.i")[-1]
    v_dc = tr.etat("bus10.v")[-1]
    assert 0.30 < i_bat * v_dc < 0.33


# ----------------------------------------------------------------------
# jacobiens des nouveaux blocs
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nom", ["RED", "BAT", "OND", "IT", "GE"])
def test_jacobien_vs_differences_finies(net, nom):
    el = element(net, nom)
    cm = el.compiled()
    A, _, _, _ = cm.linearize(el.x0, el.z0, el.u0, el.p)

    def dxdt(x):
        z = el.z0
        if cm.n_z:
            z = root(lambda zz: cm.g(x, zz, el.u0, el.p), el.z0, tol=1e-13).x
        return cm.f(x, z, el.u0, el.p)

    # differences finies centrees d'ordre 4 : la troncature de l'ordre 2
    # (h^2 fois la derivee troisieme) depasse 1e-9 sur la caracteristique
    # lissee du pont de diodes
    h = 1e-5
    An = np.zeros((cm.n_x, cm.n_x))
    for j in range(cm.n_x):
        e = np.zeros(cm.n_x)
        e[j] = h
        x0 = el.x0
        An[:, j] = (
            -dxdt(x0 + 2 * e) + 8 * dxdt(x0 + e) - 8 * dxdt(x0 - e) + dxdt(x0 - 2 * e)
        ) / (12 * h)
    err = np.max(np.abs(A - An)) / max(np.max(np.abs(A)), 1.0)
    assert err < 1e-9, f"{nom}: erreur jacobien {err:.2e}"


# ----------------------------------------------------------------------
# batterie
# ----------------------------------------------------------------------
def test_batterie_unilaterale_inerte_au_dessus_de_vref(net, modele):
    """Un echelon de charge de 0,1 pu ne fait pas debiter la batterie : son
    statisme est unilateral, le bus reste au-dessus de ``v_ref``."""
    sc = modele.scenario(Evenement(0.1, "IT.p_ref", 0.1, "rampe", 0.05))
    t = np.linspace(0.0, 3.0, 301)
    tr = modele.simuler((0, 3.0), uext_fn=sc, t_eval=t, rtol=1e-5, atol=1e-7)
    assert np.abs(tr.etat("BAT.i")).max() < 1e-3
    # courant residu du maximum lisse ~ Kv.eps^2/(4|x|) ~ 6e-4 pu, compense par
    # p_ref a l'equilibre : l'etat de charge perd moins de 1e-5 en 3 s
    assert np.abs(tr.etat("BAT.soc") - 0.9).max() < 1e-5
    # le pont, lui, a pris l'echelon (+0,1 pu cote charge = +0,10x cote continu)
    assert tr.etat("RED.idc")[-1] > tr.etat("RED.idc")[0] + 0.2


def test_batterie_decharge_et_etat_de_charge():
    """Batterie seule sur un bus continu tenu par une source : a ``t = 1 s`` sa
    puissance programmee passe de 0 a 0,1 pu ; elle debite ``P`` et son etat de
    charge baisse de ``P_bat/Eb`` par seconde (pertes internes nulles). Elle ne
    peut pas DEMARRER en decharge : l'etat de charge n'a pas d'equilibre quand la
    batterie debite, le point de depart est la veille."""
    from g2egrid.elements import BusDC, SourceDC
    from g2egrid.loadflow import DC_SLACK

    net = Network(f0=50.0, Sb=1e6)
    net.add(BusDC(10, 640e3, C=1.0, V=1.0))
    net.bus_specs_dc[10].btype = DC_SLACK
    net.add(SourceDC("S", 10, R=0.01, X=0.1, wb=net.wb, E=1.0))
    net.add(BatterieDCDC("B", 10, Eb=100.0, Rint=0.0, vref=0.5))
    net.run_loadflow()
    net.solve_equilibrium()
    m = ModeleNonLineaire(net)
    assert m.residu_equilibre()[0] < 1e-9
    t = np.linspace(0.0, 20.0, 201)
    tr = m.simuler(
        (0, 20.0),
        uext_fn=m.echelon("B.p_ref", 0.1, 1.0),
        t_eval=t,
        rtol=1e-6,
        atol=1e-8,
    )
    # le courant de la batterie vaut P/v, et la source absorbe ce que la batterie debite
    v = tr.etat("bus10.v")[-1]
    assert tr.etat("B.i")[-1] * v == pytest.approx(0.1, rel=1e-3)
    soc = tr.etat("B.soc")
    # dsoc/dt = - P/Eb = -0,1/100 = -1e-3 par seconde, Rint nulle
    k = t >= 5.0
    pente = (soc[k][0] - soc[k][-1]) / (t[k][-1] - t[k][0])
    assert pente == pytest.approx(0.1 / 100.0, rel=2e-3)


def test_batterie_parametres_invalides():
    with pytest.raises(ValueError):
        BatterieDCDC("B", 10, Eb=0.0)
    with pytest.raises(ValueError):
        BatterieDCDC("B", 10, soc0=1.2)
    with pytest.raises(ValueError):
        BatterieDCDC("B", 10, Tc=0.0)


# ----------------------------------------------------------------------
# onduleur et charge
# ----------------------------------------------------------------------
def test_onduleur_tient_v_et_f_quelle_que_soit_la_charge(net, modele):
    """Un echelon de charge de 0,1 pu ne bouge ni la tension ni la frequence de
    l'ilot : l'onduleur n'a pas de statisme (``mp = nq = 0``)."""
    ond = element(net, "OND")
    assert ond.pdict["mp"] == 0.0 and ond.pdict["nq"] == 0.0
    sc = modele.scenario(Evenement(0.1, "IT.p_ref", 0.1, "rampe", 0.05))
    t = np.linspace(0.0, 2.0, 201)
    tr = modele.simuler(
        (0, 2.0),
        uext_fn=sc,
        t_eval=t,
        rtol=1e-6,
        atol=1e-8,
        sorties=["OND.Vt", "OND.w"],
    )
    assert np.abs(tr.sorties["OND.Vt"][-1] - tr.sorties["OND.Vt"][0]) < 2e-3
    assert np.abs(tr.sorties["OND.w"] - 1.0).max() < 1e-9


def test_onduleur_refuse_un_statisme():
    with pytest.raises(ValueError, match="statisme"):
        OnduleurVF("O", 5, 10, Sn=1e6, Sb=1e6, wb=314.0, mp=0.01)


def test_charge_ac_suit_p_ref(net, modele):
    """``p_ref`` est une consigne absolue : la puissance consommee la suit."""
    sc = modele.scenario(Evenement(0.1, "IT.p_ref", 0.05, "rampe", 0.1))
    t = np.linspace(0.0, 1.0, 101)
    tr = modele.simuler(
        (0, 1.0), uext_fn=sc, t_eval=t, rtol=1e-6, atol=1e-8, sorties=["IT.p"]
    )
    assert tr.sorties["IT.p"][-1] == pytest.approx(0.35, abs=1e-4)
    assert isinstance(element(net, "IT"), ChargeACPref)


# ----------------------------------------------------------------------
# linearise contre non lineaire
# ----------------------------------------------------------------------
def test_lineaire_colle_au_non_lineaire_sur_un_petit_echelon(net, modele):
    """Echelon de 0,01 pu (3 % de la charge) : le modele lineaire reproduit le
    non lineaire a 1 % de l'excursion pres, sur le pont, le bus continu et le
    groupe. L'ecart est le terme du second ordre."""
    sc = modele.scenario(Evenement(1.0, "IT.p_ref", 0.01, "rampe", 0.05))
    t = np.linspace(0.0, 6.0, 1201)
    tn = modele.simuler((0, 6.0), uext_fn=sc, t_eval=t, rtol=1e-7, atol=1e-9)
    tl = simuler_lineaire(modele, t, sc)
    for nom in ("RED.idc", "bus10.v", "GE.Dwr", "GE.fuel", "OND.pm"):
        a, b = tn.etat(nom), tl.etat(nom)
        excursion = np.abs(a - a[0]).max()
        assert excursion > 1e-6, nom
        assert np.abs(a - b).max() < 1e-2 * excursion, nom


def test_simuler_lineaire_sans_entree_reste_a_l_equilibre(modele):
    t = np.linspace(0.0, 2.0, 41)
    tl = simuler_lineaire(modele, t, modele.scenario())
    assert np.abs(tl.x - modele.c.x0()[None, :]).max() < 1e-12


def test_simuler_lineaire_demande_un_pas_constant(modele):
    with pytest.raises(ValueError, match="pas constant"):
        simuler_lineaire(modele, np.array([0.0, 0.1, 0.3]), modele.scenario())


# ----------------------------------------------------------------------
# courbes d'entree
# ----------------------------------------------------------------------
def test_courbe_interpolations():
    t, v = [0.0, 1.0, 2.0, 4.0], [1.0, 3.0, 2.0, 2.0]
    lin = Courbe("x", t, v, "lineaire")
    pal = Courbe("x", t, v, "palier")
    pch = Courbe("x", t, v, "pchip")
    assert lin.valeur(0.5) == pytest.approx(2.0)
    assert pal.valeur(0.99) == 1.0 and pal.valeur(1.0) == 3.0
    assert pch.valeur(1.0) == pytest.approx(3.0) and 1.0 < pch.valeur(0.5) < 3.0
    # avant / apres : valeurs de bord tenues
    assert lin.valeur(-5.0) == 1.0 and lin.valeur(9.0) == 2.0
    # pchip : monotone entre deux points, sans depassement
    assert 2.0 <= pch.valeur(3.0) <= 2.0 + 1e-12
    assert pch.discontinuites() == ()
    assert lin.discontinuites() == (0.0, 1.0, 2.0, 4.0)


def test_courbe_valide_ses_entrees():
    with pytest.raises(ValueError):
        Courbe("x", [0.0, 1.0, 1.0], [1.0, 2.0, 3.0])
    with pytest.raises(ValueError):
        Courbe("x", [0.0, 1.0], [1.0, np.nan])
    with pytest.raises(ValueError):
        Courbe("x", [0.0], [1.0])
    with pytest.raises(ValueError, match="interpolation"):
        Courbe("x", [0.0, 1.0], [1.0, 2.0], "cubique")


def test_scenario_courbe_absolue_et_avertissement(modele):
    """Une courbe absolue REMPLACE la consigne ; si elle ne part pas de
    l'equilibre, le scenario le dit."""
    t = np.linspace(0.0, 2.0, 21)
    sc = modele.scenario(Courbe("IT.p_ref", t, 0.30 + 0.01 * t))
    i = modele.c.iu("IT.p_ref")
    assert sc(0.0)[i] == pytest.approx(0.30)
    assert sc(2.0)[i] == pytest.approx(0.32)
    with pytest.warns(UserWarning, match="saut de consigne"):
        modele.scenario(Courbe("IT.p_ref", t, 0.35 + 0 * t))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        # une courbe relative n'a pas a partir de l'equilibre
        modele.scenario(Courbe("IT.p_ref", t, 0.05 + 0 * t, absolue=False))
    with pytest.raises(ValueError, match="deux courbes absolues"):
        modele.scenario(
            Courbe("IT.p_ref", t, 0.30 + 0 * t), Courbe("IT.p_ref", t, 0.30 + 0 * t)
        )


def test_courbe_charge_ia():
    t = np.linspace(0.0, 300.0, 3001)
    P = courbe_charge_ia(t, P_moy=0.30, amplitude=0.25, periode=4.0)
    assert P[0] == pytest.approx(0.30)  # part de l'equilibre
    assert P.min() > 0.0
    assert P.max() == pytest.approx(0.30 * 1.25, rel=0.06)
    # reproductible, et la graine change le bruit
    assert np.array_equal(
        P, courbe_charge_ia(t, P_moy=0.30, amplitude=0.25, periode=4.0)
    )
    assert not np.array_equal(P, courbe_charge_ia(t, graine=1))
    # avant le demarrage de l'entrainement, rien ne bouge
    assert np.allclose(P[t < 19.0], 0.30)


# ----------------------------------------------------------------------
# groupe electrogene
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def groupe_sur_reseau():
    net = Network(f0=50.0, Sb=2.5e6)
    wb = net.wb
    net.add_bus(1, 400.0, btype=2, V=1.0, Pg=0.5, C=1e-3)
    net.add_bus(2, 400.0, btype=1, V=1.0, C=1e-3)
    net.add(GroupeElectrogene("GE", 1, Sn=2.5e6, Sb=2.5e6, wb=wb))
    net.add(Line("L1", 1, 2, R=0.02, X=0.2, wb=wb))
    net.add(IB("GRID", 2, R=1e-4, X=1e-3, wb=wb))
    net.run_loadflow()
    net.solve_equilibrium()
    return net


def test_groupe_dimensions_et_equilibre(groupe_sur_reseau):
    """17 etats de l'alternateur + 5 de l'entrainement ; couple d'arbre =
    couple electrique a l'equilibre ; torsion statique ``C/Ksh``."""
    ge = element(groupe_sur_reseau, "GE")
    assert len(ge.states) == 17 + 5
    assert ge.states[:17] == type(ge).__mro__[1].states
    x = dict(zip(ge.states, ge.x0))
    assert x["Dwe"] == pytest.approx(0.0, abs=1e-12)
    assert x["Dwr"] == pytest.approx(0.0, abs=1e-12)
    c_arbre = sortie(groupe_sur_reseau, "GE", "C_arbre")
    assert c_arbre == pytest.approx(sortie(groupe_sur_reseau, "GE", "Ce"), rel=1e-12)
    assert x["dsh"] * ge.pdict["Ksh"] == pytest.approx(c_arbre, rel=1e-12)
    assert ModeleNonLineaire(groupe_sur_reseau).residu_equilibre()[0] < 1e-9


def test_groupe_mode_de_torsion_a_f_tors(groupe_sur_reseau):
    """La raideur tiree de ``f_tors`` redonne le mode de torsion a cette
    frequence (a l'amortissement et au couplage au reseau pres) et l'amortissement
    relatif demande."""
    ss = groupe_sur_reseau.build_state_space()
    lam = np.linalg.eigvals(ss.A)
    cible = 2 * np.pi * 25.0
    proche = lam[(lam.imag > 0) & (np.abs(lam.imag - cible) < 0.1 * cible)]
    assert len(proche) == 1
    lam0 = proche[0]
    assert abs(lam0.imag) / (2 * np.pi) == pytest.approx(25.0, rel=0.02)
    zeta = -lam0.real / abs(lam0)
    assert zeta == pytest.approx(0.03, rel=0.15)


def test_groupe_butee_de_carburant(groupe_sur_reseau):
    """Une consigne de puissance trop haute amene la cremaillere a ``Fmax`` : le
    mode de butee bascule, le couple moteur est plafonne et l'integrateur ne
    s'emballe pas (le couple d'arbre ne depasse pas ``Fmax`` de plus de 1 %)."""
    ge = element(groupe_sur_reseau, "GE")
    fmax = ge.pdict["Fmax"]
    m = ModeleNonLineaire(groupe_sur_reseau)
    t = np.linspace(0.0, 12.0, 1201)
    tr = m.simuler(
        (0, 12.0),
        uext_fn=m.echelon("GE.p_ref", 0.6, 1.0),
        t_eval=t,
        rtol=1e-6,
        atol=1e-8,
        sorties=["GE.C_arbre"],
    )
    assert any(b[1:3] == ("GE", "mode_fuel") and b[3] == 1.0 for b in tr.basculements)
    assert tr.etat("GE.fuel").max() <= fmax * (1 + 1e-3)
    assert tr.sorties["GE.C_arbre"].max() <= fmax * 1.01
    # puis le retour de consigne relache la butee
    tr2 = m.simuler(
        (0, 25.0),
        uext_fn=m.scenario(
            Evenement(1.0, "GE.p_ref", 0.6), Evenement(10.0, "GE.p_ref", -0.6)
        ),
        t_eval=np.linspace(0.0, 25.0, 2501),
        rtol=1e-6,
        atol=1e-8,
    )
    assert any(b[1:3] == ("GE", "mode_fuel") and b[3] == 0.0 for b in tr2.basculements)
    assert tr2.etat("GE.fuel")[-1] == pytest.approx(tr2.etat("GE.fuel")[0], abs=2e-3)


def test_groupe_parametres_invalides():
    with pytest.raises(ValueError):
        GroupeElectrogene("G", 1, Sn=1e6, Sb=1e6, wb=314.0, He=-1.0)
    with pytest.raises(ValueError, match="Fmax"):
        GroupeElectrogene("G", 1, Sn=1e6, Sb=1e6, wb=314.0, Fmax=0.0, Fmin=0.5)


def test_groupe_couple_d_arbre_excite_par_dT_exc(groupe_sur_reseau):
    """Un couple d'excitation sinusoidal a la frequence de torsion cree une
    resonance : le couple d'arbre depasse nettement l'excitation (facteur de
    surtension ~ 1/(2 zeta) = 16)."""
    m = ModeleNonLineaire(groupe_sur_reseau)
    t = np.linspace(0.0, 4.0, 4001)
    ampl = 0.001

    def uext(tt):
        u = m.uext0()
        u[m.c.iu("GE.dT_exc")] = ampl * np.sin(2 * np.pi * 25.2 * tt)
        return u

    tr = m.simuler(
        (0, 4.0),
        uext_fn=uext,
        t_eval=t,
        rtol=1e-7,
        atol=1e-10,
        sorties=["GE.C_arbre"],
        pas=None,
    )
    dev = tr.sorties["GE.C_arbre"] - tr.sorties["GE.C_arbre"][0]
    assert np.abs(dev[-1000:]).max() > 5 * ampl


# ----------------------------------------------------------------------
# la chaine complete : simulation puis fatigue
# ----------------------------------------------------------------------
def test_chaine_simulation_fatigue(groupe_sur_reseau):
    """Un echelon du groupe puis le comptage de fatigue : le couple d'arbre
    simule donne un dommage fini et positif (la pente de Haibach compte les
    cycles sous le coude), et le modele lineaire (petit signal) donne le meme
    plus grand cycle a 5 % pres."""
    m = ModeleNonLineaire(groupe_sur_reseau)
    t = np.linspace(0.0, 30.0, 6001)
    # consigne de puissance qui alterne : 12 cycles de +-0,05 pu
    p = 0.05 * np.sign(np.sin(2 * np.pi * t / 2.5)) * (t > 2.0)
    sc = m.scenario(Courbe("GE.p_ref", t, m.uext0()[m.c.iu("GE.p_ref")] + p, "palier"))
    tn = m.simuler(
        (0, 30.0), uext_fn=sc, t_eval=t, rtol=1e-6, atol=1e-8, sorties=["GE.C_arbre"]
    )
    tl = simuler_lineaire(m, t, sc, sorties=["GE.C_arbre"])
    rn = evaluer_fatigue(t, tn.sorties["GE.C_arbre"], 2.5e6, seuil_pu=1e-4)
    rl = evaluer_fatigue(t, tl.sorties["GE.C_arbre"], 2.5e6, seuil_pu=1e-4)
    assert rn.n_cycles >= 10
    assert 0.0 < rn.dommage < 1.0
    assert rn.amplitude_max > 0.0
    assert rl.amplitude_max == pytest.approx(rn.amplitude_max, rel=0.05)


def test_ilot_apres_perte_du_reseau_court():
    """Le scenario de reference, en version courte (10 s) : la batterie prend la
    charge avant l'ouverture, le pont se bloque, l'ilot est tenu par le groupe a
    la frequence du statisme, la tension continue reste dans [0,95 ; 1,10]."""
    net = reseau_datacenter()
    m = ModeleNonLineaire(net)
    t = np.linspace(0.0, 10.0, 201)
    P = np.full_like(t, 0.30)
    sc = scenario_perte_reseau(
        m, t, P, t_perte=2.0, t_reprise=4.0, duree_reprise=4.0, anticipation=0.5
    )
    tr = m.simuler(
        (0, 10.0),
        uext_fn=sc.uext,
        t_eval=t,
        rtol=1e-4,
        atol=1e-6,
        sorties=["GE.P", "GE.Omega"],
    )
    # seuls les etats de la liaison ouverte valent NaN apres son retrait
    nan = [tr.etats[j] for j in np.flatnonzero(~np.isfinite(tr.x[-1]))]
    assert nan == ["L_AMONT.i_d", "L_AMONT.i_q"]
    assert np.isfinite(tr.x[t < 1.9]).all()
    # entre la prise de charge par la batterie et la fin du scenario, le pont a ete
    # bloque (idc < 3 %) puis a repris la charge
    idc = tr.etat("RED.idc")
    assert idc[(t > 1.8) & (t < 2.2)].max() < 0.03  # pic de 1,4 % a l'ouverture
    assert idc[-1] > 0.5
    vdc = tr.etat("bus10.v")
    assert 0.95 < vdc.min() and vdc.max() < 1.10
    # fin de la reprise : le groupe porte la charge (0,30 pu systeme + pertes)
    assert tr.sorties["GE.P"][-1] * 0.5 == pytest.approx(0.32, abs=0.02)
    assert NOMS["groupe"] == "GE"
