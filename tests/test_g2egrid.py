"""Tests de validation de g2egrid.

Les references viennent du framework MATLAB G2ELIN (``matlab/core``) : memes
equations, memes valeurs par defaut. Plusieurs tests verifient explicitement
cette equivalence (nombre d'etats par modele, parametres par defaut).
"""

import json
import pathlib

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment, root

from g2egrid import Network
from g2egrid.elements import (
    GFL1,
    GFL2,
    GFM1,
    GFM2,
    GFM3,
    GFM4,
    IB,
    Branch,
    ILoad,
    Line,
    PLoad,
    RLLoad,
    SynchronousMachine,
    Transformer,
    ZLoad,
)
from g2egrid.elements.sync_machine import SM_DEFAULT
from g2egrid.tables import Tables, assembler
from g2egrid.temporel import ModeleNonLineaire
from g2egrid.separation import reponse_indicielle, separer, troncature_equilibree
from g2egrid.frequence import (
    admittance_dq,
    axes_tension,
    gain,
    impedance_dq,
    vecteur_complexe,
)


# ----------------------------------------------------------------------
# fixtures
# ----------------------------------------------------------------------
def build_smib(C=1e-3, Pg=0.8):
    """SMIB : machine (bus 1, PV) -- ligne -- reseau amont (bus 2, slack)."""
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=2, V=1.0, Pg=Pg, C=C)
    net.add_bus(2, 20e3, btype=1, V=1.0, C=C)
    net.add(SynchronousMachine("SM1", 1, Sn=1e6, Sb=1e6, wb=net.wb))
    net.add(Line("L1", 1, 2, R=0.01, X=0.20, wb=net.wb))
    net.add(IB("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
    return net


def build_test08(C=1e-9):
    """Reseau de test_elem_08 : slack 110 kV -- transfo -- charge 20 kV."""
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 110e3, btype=1, V=1.0, C=C)
    net.add_bus(2, 20e3, btype=3, Pl=0.5, Ql=0.1, C=C)
    net.add(Transformer("Tr1", 1, 2, R=0.0002, X=0.024, wb=net.wb))
    net.add(RLLoad("Load2", 2, P=0.5, Q=0.1, wb=net.wb))
    net.add(IB("GRID", 1, R=1e-4, X=1e-3, wb=net.wb))
    return net


def build_converter(cls, X_line=0.20, Pg=0.5, C=1e-3, **kw):
    """Convertisseur (bus 1) -- ligne -- reseau amont (bus 2, slack)."""
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=3, V=1.0, Pg=Pg, C=C)
    net.add_bus(2, 20e3, btype=1, V=1.0, C=C)
    net.add(cls("CONV", 1, Sn=1e6, Sb=1e6, wb=net.wb, **kw))
    net.add(Line("L1", 1, 2, R=X_line / 10, X=X_line, wb=net.wb))
    net.add(IB("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
    return net


def dc_gain(ss, y, u):
    """Gain statique y/u, robuste a un mode integrateur (A singuliere)."""
    s = 1e-4
    G = ss.C @ np.linalg.solve(s * np.eye(len(ss.states)) - ss.A, ss.B) + ss.D
    return G[ss.iy(y), ss.iu(u)]


# ----------------------------------------------------------------------
# equivalence avec les modeles MATLAB
# ----------------------------------------------------------------------
@pytest.mark.parametrize(
    "cls, n_x, n_z",
    [
        (SynchronousMachine, 17, 6),  # SM.m : 17 etats (PSS + AVR 4 etats)
        (IB, 3, 0),  # IB.m : i_d, i_q, theta
        (Branch, 2, 0),  # TransmissionLine.m / Transformer.m
        (RLLoad, 4, 0),  # GenericLoad.m
        (ZLoad, 2, 0),  # AdvancedLoad.m, mode 'Z'
        (ILoad, 2, 0),  # AdvancedLoad.m, mode 'I'
        (PLoad, 2, 0),  # AdvancedLoad.m, mode 'P'
        (GFM1, 13, 0),
        (GFM2, 14, 0),
        (GFM3, 15, 0),
        (GFM4, 16, 0),
        # 3 algebriques : md, mq et w_pll, comme algVec dans symGFL. Les
        # garder algebriques plutot qu'expressions reproduit la substitution
        # de la reference, qui y met ses propres valeurs d'equilibre.
        (GFL1, 12, 3),
        (GFL2, 14, 3),
    ],
)
def test_nombre_etats_conforme_a_matlab(cls, n_x, n_z):
    cm = cls.compiled()
    assert (cm.n_x, cm.n_z) == (n_x, n_z)


def test_defauts_machine_conformes_a_sm_m():
    """Valeurs par defaut reprises litteralement du constructeur de SM.m."""
    attendu = dict(
        Ra=0.003,
        Ll=0.15,
        Lad=1.66,
        Laq=1.61,
        Rg=100.0,
        H=5.0,
        KD=1.0,
        mp=0.5,
        Tg=0.2,
        Tr=0.02,
        Ka=300.0,
        Ta=0.001,
        Ke=1.0,
        Te=0.0001,
        Kfd=0.001,
        Tfd=0.1,
        K_PSS=20.0,
        T_LP=0.03,
        T_HP=2.0,
    )
    for k, v in attendu.items():
        assert SM_DEFAULT[k] == v, k


def test_etats_machine_nommes_comme_matlab():
    assert SynchronousMachine.states == (
        "phi_d",
        "phi_q",
        "phi_fd",
        "phi_1d",
        "phi_1q",
        "phi_2q",
        "Dwr",
        "Pm",
        "Dw1",
        "v1",
        "v2",
        "vpss",
        "e1",
        "e2",
        "efd",
        "e3",
        "theta",
    )


def test_compilation_et_cache():
    cm = SynchronousMachine.compiled()
    assert cm.spec.inputs[-2:] == ("vd_t", "vq_t")
    assert cm.spec.outputs[:2] == ("id_t", "iq_t")
    assert SynchronousMachine.compiled() is cm  # mis en cache


def test_ports_branche():
    cm = Branch.compiled()
    # deux entrees globales avant les tensions de port : la frequence systeme
    # et l'angle du repere de reference (celui du slack), comme theta_g en V1.
    assert cm.spec.inputs == ("w_ref", "theta_g", "vd_f", "vq_f", "vd_t", "vq_t")
    assert cm.spec.outputs[:4] == ("id_f", "iq_f", "id_t", "iq_t")


def test_noeud_sans_capacite_ideal_ou_resistif():
    """Un noeud sans capacite shunt est un noeud ideal par defaut (FR-005) ;
    en mode ``"resistance"``, comme Network.m, il recoit une resistance de
    causalite : les 100 pu de Node.m, ou celle qu'on demande."""
    from g2egrid.elements import Bus, IdealBus, ResistiveBus
    from g2egrid.elements.bus import R_NODE_M

    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=1, V=1.0)  # C = 0 par defaut
    net.add_bus(2, 20e3, btype=3, C=0.05)
    assert isinstance(net.buses[1], IdealBus)
    assert not isinstance(net.buses[2], (IdealBus, ResistiveBus))
    for r_causalite, attendu in ((None, R_NODE_M), (1e3, 1e3)):
        kw = {} if r_causalite is None else {"r_causalite": r_causalite}
        net = Network(f0=50.0, Sb=1e6, noeuds_sans_capacite="resistance", **kw)
        net.add_bus(1, 20e3, btype=1, V=1.0)  # C = 0 par defaut
        net.add_bus(2, 20e3, btype=3, C=0.05)
        assert isinstance(net.buses[1], ResistiveBus)
        assert net.buses[1].pdict["R"] == attendu
        assert isinstance(net.buses[2], Bus)


# ----------------------------------------------------------------------
# load flow
# ----------------------------------------------------------------------
def test_loadflow_vs_matlab():
    """Reference MATLAB : bus 2 |V| = 0.9974 pu, angle = -0.69 deg."""
    lf = build_test08().run_loadflow()
    assert lf.converged and lf.mismatch < 1e-10
    assert abs(lf.V[1]) == pytest.approx(0.9974, abs=5e-5)
    assert np.degrees(np.angle(lf.V[1])) == pytest.approx(-0.69, abs=5e-3)


def test_loadflow_consigne_pv():
    net = build_smib()
    lf = net.run_loadflow()
    assert lf.S[0].real == pytest.approx(0.8, abs=1e-9)  # P impose au bus PV
    assert abs(lf.V[0]) == pytest.approx(1.0, abs=1e-9)  # |V| impose


def test_rg_de_la_machine_hors_du_ybus():
    """Rg est deja dans le courant injecte par la machine : pas au Y-bus."""
    net = build_smib()
    sm = net.elements[0]
    Y = net.build_ybus()
    net.elements.remove(sm)  # meme reseau, sans la machine
    Y0 = net.build_ybus()
    assert np.allclose(Y, Y0, rtol=0, atol=1e-15)


# ----------------------------------------------------------------------
# equilibre
# ----------------------------------------------------------------------
def test_residus_equilibre_nuls():
    net = build_smib()
    net.run_loadflow()
    net.solve_equilibrium()
    for el in net.all_elements:
        cm = el.compiled()
        r = cm.residual(el.x0, el.z0, el.u0, el.p)
        if r.size:
            assert np.max(np.abs(r)) < 1e-8, f"{el.name}: {np.max(np.abs(r))}"


def test_loi_des_noeuds():
    """Somme des courants injectes = courant des shunts du bus.

    Shunts = capacite du noeud (j.w.C.V) + shunts declares par les elements
    (noeud resistif). C'est le test qui prouve que load flow, equilibre par
    element et interconnexion sont coherents entre eux.
    """
    net = build_smib()
    net.run_loadflow()
    net.solve_equilibrium()
    for b, i in net._bus_index().items():
        tot = 0j
        for el in net.elements:
            y = el.compiled().h(el.x0, el.z0, el.u0, el.p)
            for pk, bb in enumerate(el.buses):
                if bb == b:
                    a, c = el.port_output_idx(pk)
                    tot += complex(y[a], y[c])
        ysh = 1j * net.w_sys * net.bus_specs[b].C
        for el in net.all_elements:
            adm = el.shunt_admittance()
            if adm is not None and el.buses[0] == b:
                ysh += adm
        want = ysh * net.lf.V[i]
        assert abs(tot - want) < 1e-9, f"bus {b}: {tot} != {want}"


def test_puissance_machine_conforme_au_loadflow():
    net = build_smib(Pg=0.8)
    net.run_loadflow()
    net.build_state_space()
    sm = net.elements[0]
    y = sm.compiled().h(sm.x0, sm.z0, sm.u0, sm.p)
    assert y[sm.compiled().iy("P")] == pytest.approx(0.8, abs=1e-6)


def test_pss_et_avr_a_lequilibre():
    """PSS a sortie nulle, AVR cale : e1 = Et, e3 = 0, e2 = efd/Ke."""
    net = build_smib()
    net.run_loadflow()
    net.solve_equilibrium()
    sm = net.elements[0]
    x = dict(zip(sm.states, sm.x0))
    for nm in ("Dw1", "v1", "v2", "vpss", "e3", "Dwr"):
        assert abs(x[nm]) < 1e-9, nm
    assert x["e1"] == pytest.approx(abs(net.lf.V[0]), abs=1e-9)
    assert x["e2"] * sm.pdict["Ke"] == pytest.approx(x["efd"], abs=1e-9)


# ----------------------------------------------------------------------
# linearisation
# ----------------------------------------------------------------------
def test_jacobiens_vs_differences_finies():
    """Les jacobiens symboliques doivent coller aux differences finies."""
    net = build_smib()
    net.run_loadflow()
    net.solve_equilibrium()
    for el in net.all_elements:
        cm = el.compiled()
        if cm.n_x == 0:
            continue
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


# ----------------------------------------------------------------------
# modele global
# ----------------------------------------------------------------------
def test_dimensions_modele_global():
    net = build_smib()
    net.run_loadflow()
    ss = net.build_state_space()
    # SM 17 + ligne 2 + IB 3 + 2 bus x 2 = 26
    assert len(ss.states) == 26
    assert "w_ref" in ss.inputs and "SM1.p_ref" in ss.inputs
    assert "SM1.P" in ss.outputs and "bus1.vmag" in ss.outputs
    # les courants de port restent internes
    assert not any(o.endswith(".id_t") for o in ss.outputs)


def test_mode_electromecanique():
    """Un SMIB doit presenter un mode electromecanique vers 1-2 Hz,
    faiblement amorti, et le systeme doit rester stable."""
    net = build_smib()
    net.run_loadflow()
    ss = net.build_state_space()
    ev = ss.eig()
    assert ev.real.max() < 1e-6
    osc = ev[ev.imag > 1e-3]
    f = np.abs(osc.imag) / (2 * np.pi)
    em = osc[(f > 0.5) & (f < 3.0)]
    assert em.size >= 1, "aucun mode electromecanique entre 0.5 et 3 Hz"
    zeta = -em[0].real / abs(em[0])
    assert 0.0 < zeta < 0.40


def test_pss_agit_sur_le_mode_electromecanique():
    """Le PSS doit agir sur le mode local.

    ATTENTION -- resultat a connaitre : avec les constantes de SM.m
    (K_PSS = 20, T_LP = 0.03, T_HP = 2, avance 0.05/0.02 puis retard 3/5.4),
    l'avance de phase ne compense pas le retard GEP de l'excitation rapide :
    sur ce SMIB le PSS *reduit* l'amortissement du mode local (15.3 % -> 9.8 %)
    au lieu de l'augmenter, et c'est un gain NEGATIF qui l'amortit. Le portage
    reproduit fidelement SM.m ; le reglage est a revoir des deux cotes. Le test
    verifie donc l'action du PSS, pas son signe.
    """

    def zeta_local(K_PSS):
        net = Network(f0=50.0, Sb=1e6)
        net.add_bus(1, 20e3, btype=2, V=1.0, Pg=0.8, C=1e-3)
        net.add_bus(2, 20e3, btype=1, V=1.0, C=1e-3)
        net.add(SynchronousMachine("SM1", 1, Sn=1e6, Sb=1e6, wb=net.wb, K_PSS=K_PSS))
        net.add(Line("L1", 1, 2, R=0.01, X=0.20, wb=net.wb))
        net.add(IB("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
        net.run_loadflow()
        ev = net.build_state_space().eig()
        osc = ev[ev.imag > 1e-3]
        f = np.abs(osc.imag) / (2 * np.pi)
        em = osc[(f > 0.5) & (f < 3.0)][0]
        return -em.real / abs(em)

    z0, z20, zneg = zeta_local(0.0), zeta_local(20.0), zeta_local(-20.0)
    assert abs(z20 - z0) > 0.02, "le PSS n'a aucun effet"
    assert zneg > z0 > z20, (zneg, z0, z20)


def test_capacite_bus_ne_pollue_pas_les_modes_lents():
    """La capacite shunt est un element HF : elle ne doit pas deplacer le
    mode electromecanique."""

    def fem(C):
        net = build_smib(C=C)
        net.run_loadflow()
        ev = net.build_state_space().eig()
        osc = ev[ev.imag > 1e-3]
        f = np.abs(osc.imag) / (2 * np.pi)
        return osc[(f > 0.5) & (f < 3.0)][0]

    a, b = fem(1e-2), fem(1e-4)  # C varie d'un facteur 100
    assert abs(a - b) / abs(a) < 2e-2


def test_reponse_indicielle_puissance():
    """Un echelon de p_ref doit augmenter la puissance de la machine."""
    net = build_smib()
    net.run_loadflow()
    ss = net.build_state_space()
    g = dc_gain(ss, "SM1.P", "SM1.p_ref")
    assert g.real > 0.5, f"gain p_ref -> P trop faible ({g:.3f})"


# ----------------------------------------------------------------------
# convertisseurs
# ----------------------------------------------------------------------
@pytest.mark.parametrize(
    "cls, n_states",
    [
        (GFL1, 12 + 2 + 3 + 4),
        (GFM1, 13 + 2 + 3 + 4),
        (GFM2, 14 + 2 + 3 + 4),
        (GFM3, 15 + 2 + 3 + 4),
        (GFM4, 16 + 2 + 3 + 4),
    ],
)
def test_convertisseur_equilibre(cls, n_states):
    net = build_converter(cls)
    net.run_loadflow()
    net.solve_equilibrium()
    for el in net.all_elements:
        cm = el.compiled()
        r = cm.residual(el.x0, el.z0, el.u0, el.p)
        if r.size:
            assert np.max(np.abs(r)) < 1e-8, el.name
    ss = net.build_state_space()
    assert len(ss.states) == n_states
    assert ss.eig().real.max() < 1e-6


def test_gfl_pll_verrouillee():
    """A l'equilibre la PLL est verrouillee sur la tension du condensateur :
    veq = 0, w_pll = 1, et la puissance injectee suit le load flow."""
    net = build_converter(GFL1, Pg=0.5)
    net.run_loadflow()
    net.solve_equilibrium()
    g = net.elements[0]
    cm = g.compiled()
    x = dict(zip(g.states, g.x0))
    y = cm.h(g.x0, g.z0, g.u0, g.p)
    assert x["veq"] == pytest.approx(0.0, abs=1e-9)
    assert y[cm.iy("w_pll_out")] == pytest.approx(1.0, abs=1e-9)
    iD, iQ = y[0], y[1]
    S = net.lf.V[0] * np.conj(complex(iD, iQ))
    assert S.real == pytest.approx(0.5, abs=1e-7)


def test_gfm_ilotage_mode_marginal():
    """Un GFM seul (aucune reference d'angle) presente un mode marginal
    d'angle absolu -- c'est la jauge, pas une instabilite."""
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=1, V=1.0, C=1e-3)
    net.add_bus(2, 20e3, btype=3, Pl=0.5, Ql=0.1, C=1e-3)
    net.add(GFM1("GFM_1", 1, Sn=1e6, Sb=1e6, wb=net.wb))
    net.add(Line("L1", 1, 2, R=0.01, X=0.10, wb=net.wb))
    net.add(RLLoad("Load2", 2, P=0.5, Q=0.1, wb=net.wb))
    net.run_loadflow()
    ev = net.build_state_space().eig()
    assert np.min(np.abs(ev)) < 1e-6, "mode marginal de jauge absent"
    assert ev.real.max() < 1e-6, "systeme instable"


def test_gfl_se_degrade_en_reseau_faible():
    """Comportement de reference du grid-following : la marge de stabilite
    se degrade quand le reseau s'affaiblit (X croissant)."""
    marges = []
    for X in (0.1, 0.4, 1.0):
        net = build_converter(GFL1, X_line=X)
        net.run_loadflow()
        ev = net.build_state_space().eig()
        marges.append(np.max(ev.real[np.abs(ev) > 1e-6]))
    assert marges[0] < marges[1] < marges[2] < 0, marges


def test_pll_trop_rapide_destabilise():
    """Accelerer la PLL (tau_PLL plus court) finit par destabiliser le GFL
    raccorde a un reseau faible."""
    net = build_converter(GFL1, X_line=1.0, tau_PLL=0.002)
    net.run_loadflow()
    ev = net.build_state_space().eig()
    assert np.max(ev.real[np.abs(ev) > 1e-6]) > 0


# ----------------------------------------------------------------------
# charges speciales
# ----------------------------------------------------------------------
def build_charge(cls, P=0.5, Q=0.1, **kw):
    """test_elem_08 avec une charge de la classe donnee."""
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 110e3, btype=1, V=1.0, C=1e-9)
    net.add_bus(2, 20e3, btype=3, Pl=P, Ql=Q, C=1e-9)
    net.add(Transformer("Tr1", 1, 2, R=0.0002, X=0.024, wb=net.wb))
    if cls in (RLLoad, ZLoad):
        kw["wb"] = net.wb
    net.add(cls("Load2", 2, P=P, Q=Q, **kw))
    net.add(IB("GRID", 1, R=1e-4, X=1e-3, wb=net.wb))
    return net


@pytest.mark.parametrize("cls, exposant", [(ZLoad, 2), (ILoad, 1), (PLoad, 0)])
def test_charge_sensibilite_tension(cls, exposant):
    """En regime permanent, P varie comme |V|^n : n = 2 (Z), 1 (I), 0 (P)."""
    net = build_charge(cls)
    net.run_loadflow()
    net.solve_equilibrium()
    el = next(e for e in net.elements if e.name == "Load2")
    cm = el.compiled()
    r = cm.residual(el.x0, el.z0, el.u0, el.p)
    assert np.max(np.abs(r)) < 1e-10
    A, B, C, D = cm.linearize(el.x0, el.z0, el.u0, el.p)
    G = D - C @ np.linalg.solve(A, B)
    vd, vq = el.u0[cm.iu("vd_t")], el.u0[cm.iu("vq_t")]
    V = np.hypot(vd, vq)
    ip = cm.iy("p")
    dP_dV = G[ip, cm.iu("vd_t")] * vd / V + G[ip, cm.iu("vq_t")] * vq / V
    P = vd * el.x0[cm.ix("i_d")] + vq * el.x0[cm.ix("i_q")]
    assert P == pytest.approx(0.5, rel=1e-9)
    assert dP_dV == pytest.approx(exposant * 0.5 / V, abs=1e-9)


#: ILoad etait INSTABLE dans ce montage (C = 1e-9 pu aux deux bus : +61 s^-1 a
#: Q = 0,1, +123 s^-1 a Q = -0,2), ce que masquait l'ancien seuil de
#: ``is_stable`` (1e-9 . max|A| = 314 s^-1). Le modele a ete corrige par le
#: doctorant convertisseurs (revue du 2026-10-25) et le nouveau critere
#: ``is_stable`` le confirme : ILoad est de nouveau teste comme les autres.


@pytest.mark.parametrize(
    "cls, Q",
    [
        (ZLoad, 0.1),
        (PLoad, 0.1),
        (ILoad, 0.1),
        (PLoad, -0.2),
        (ILoad, -0.2),
    ],
)
def test_charge_speciale_stable(cls, Q):
    """Pas d'inductance : un reactif negatif reste stable."""
    net = build_charge(cls, Q=Q)
    net.run_loadflow()
    net.solve_equilibrium()
    assert net.build_state_space().is_stable()


@pytest.mark.parametrize("cls", [ILoad, PLoad])
def test_charge_speciale_accepte_q_negatif(cls):
    """Pas d'inductance : un reactif negatif donne un load flow et un modele."""
    net = build_charge(cls, Q=-0.2)
    net.run_loadflow()
    assert net.build_state_space().A.shape[0] > 0


# ----------------------------------------------------------------------
# reponse frequentielle : nature des sources
# ----------------------------------------------------------------------
def banc_frequentiel(cls, P=0.5, Q=0.0, **kw):
    """Equipement seul face a un reseau infini, a V = 1 pu et angle nul."""
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=3, Pg=P, Qg=Q, C=1e-6)
    net.add_bus(2, 20e3, btype=1, V=1.0, C=1e-6)
    net.add(cls("EQ", 1, Sn=1e6, Sb=1e6, wb=net.wb, **kw))
    net.add(Line("L", 1, 2, R=1e-6, X=1e-5, wb=net.wb))
    net.add(IB("GRID", 2, R=1e-6, X=1e-5, wb=net.wb, Sn=1e10, Sb=1e6))
    net.run_loadflow()
    net.solve_equilibrium()
    return next(e for e in net.elements if e.name == "EQ")


def test_admittance_charge_rl_vaut_son_impedance():
    """Sur une branche R-L, Y(jw) doit rendre exactement 1/(R + jX(w/wb + w0))."""
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 110e3, btype=1, V=1.0, C=1e-9)
    net.add_bus(2, 20e3, btype=3, Pl=0.5, Ql=0.2, C=1e-9)
    net.add(Transformer("Tr1", 1, 2, R=0.0002, X=0.024, wb=net.wb))
    net.add(ZLoad("Load2", 2, P=0.5, Q=0.2, wb=net.wb))
    net.add(IB("GRID", 1, R=1e-4, X=1e-3, wb=net.wb))
    net.run_loadflow()
    net.solve_equilibrium()
    el = next(e for e in net.elements if e.name == "Load2")

    f = np.array([0.0, 0.1, 1.0, 10.0, 50.0])
    Z_direct, Z_inverse = vecteur_complexe(impedance_dq(el, f))
    R, X, wb = el.pdict["R"], el.pdict["X"], el.pdict["wb"]
    attendu = R + 1j * X * (2 * np.pi * f / wb + 1.0)
    assert np.allclose(Z_direct, attendu, rtol=1e-9)
    assert np.max(np.abs(Z_inverse)) < 1e-12  # branche symetrique : pas de couplage


@pytest.mark.parametrize(
    "cls, kw, borne_min, borne_max",
    [
        # GFM : |Y| = 1/(R2 + L2) = 1/0.013, la boucle de tension tient le
        # condensateur et il ne reste que le filtre de raccordement
        (GFM1, {}, 50.0, 200.0),
        (GFM3, {}, 50.0, 200.0),
        # machine synchrone : raideur donnee par le regulateur de tension
        (SynchronousMachine, {}, 20.0, 500.0),
        # GFL a boucles P/Q : |Y| = S/V^2, source de courant
        (GFL1, {}, 0.3, 1.0),
    ],
)
def test_nature_de_source(cls, kw, borne_min, borne_max):
    """|Y| en basse frequence : grand pour une source de tension, petit sinon."""
    el = banc_frequentiel(cls, **kw)
    Y = admittance_dq(el, np.array([0.1]))
    e_mod, _ = axes_tension(el)
    assert borne_min < gain(Y, e_mod)[0] < borne_max


def test_gfm_vaut_l_inverse_de_son_filtre():
    """Sous 1 Hz, l'admittance d'un GFM est celle de sa branche de raccordement.

    A quelques pour cent pres : le condensateur du filtre et le gain fini de la
    boucle de tension ajoutent 0,0006 pu a l'impedance vue du reseau.
    """
    el = banc_frequentiel(GFM1)
    e_mod, _ = axes_tension(el)
    y = gain(admittance_dq(el, np.array([0.01])), e_mod)[0]
    z_filtre = abs(complex(el.pdict["R2_pu"], el.pdict["L2_pu"]))
    assert 1.0 / y == pytest.approx(z_filtre, rel=0.1)


# ----------------------------------------------------------------------
# separation transport / distribution
# ----------------------------------------------------------------------
def reseau_td():
    """Bus infini (1) -- PCC (2) -- noeud de distribution (3) : charge et GFL."""
    net = Network(f0=50.0, Sb=10e6)
    net.add_bus(1, 20e3, btype=1, V=1.0, C=1e-3)
    net.add_bus(2, 20e3, btype=3, C=1e-3)
    net.add_bus(3, 20e3, btype=3, Pg=0.3, Pl=0.6, Ql=0.2, C=1e-3)
    net.add(IB("GRID", 1, R=1e-3, X=1e-2, wb=net.wb))
    net.add(Line("L12", 1, 2, R=0.01, X=0.1, wb=net.wb))
    net.add(Line("L23", 2, 3, R=0.02, X=0.05, wb=net.wb))
    net.add(RLLoad("Load3", 3, P=0.6, Q=0.2, wb=net.wb))
    net.add(GFL1("GFL3", 3, Sn=5e6, Sb=10e6, wb=net.wb))
    net.run_loadflow()
    net.solve_equilibrium()
    ss = net.build_state_space()
    dist = {"L23", "Load3", "GFL3", "bus3"}
    sep = separer(ss, lambda e: e.split(".")[0] in dist, entrees=["GFL3.p_ref"])
    return ss, sep


def test_separation_reconnexion_exacte():
    """Transport + distribution complete reconnectes = modele global."""
    ss, sep = reseau_td()
    A, _, _ = sep.modele_global()
    ordre = np.concatenate([sep.iT, sep.iD])
    assert np.array_equal(A, ss.A[np.ix_(ordre, ordre)])
    # la distribution n'agit que sur la tension du PCC
    assert sep.distribution.sorties == ["d/dt bus2.vd", "d/dt bus2.vq"]


def test_elimination_exacte_redonne_les_valeurs_propres():
    """lambda.I - A_eff(lambda) est singuliere pour chaque valeur propre lente."""
    ss, sep = reseau_td()
    lam = np.linalg.eigvals(ss.A)
    for lk in lam[np.abs(lam) < 50]:
        s = np.linalg.svd(
            lk * np.eye(len(sep.iT)) - sep.a_effectif(lk), compute_uv=False
        )
        assert s[-1] / s[0] < 1e-9


def test_residualisation_ordre_zero_est_le_quasi_statique():
    ss, sep = reseau_td()
    red, _ = troncature_equilibree(sep.distribution, 0, methode="residualisation")
    assert red.ordre == 0
    A0, _, _ = sep.modele_global(red)
    assert np.allclose(
        A0, sep.a_quasi_statique(), rtol=1e-8, atol=1e-8 * np.abs(A0).max()
    )


@pytest.mark.parametrize("ordre", [2, 6, 10])
def test_residualisation_conserve_le_gain_statique(ordre):
    """La residualisation garde le gain statique ; la troncature, non."""
    ss, sep = reseau_td()
    d = sep.distribution
    red, _ = troncature_equilibree(d, ordre, methode="residualisation")
    G0 = d.transfert(np.array([0.0]))[0]
    assert np.allclose(
        red.transfert(np.array([0.0]))[0], G0, rtol=1e-6, atol=1e-9 * np.abs(G0).max()
    )
    assert red.entrees == d.entrees and red.sorties == d.sorties


def test_reponse_indicielle_premier_ordre():
    A, B = np.array([[-2.0]]), np.array([[4.0]])
    t = np.linspace(0.0, 3.0, 301)
    x = reponse_indicielle(A, B, t)[:, 0]
    assert np.allclose(x, 2.0 * (1.0 - np.exp(-2.0 * t)), atol=1e-12)


# ----------------------------------------------------------------------
# repere global : invariance au choix de la reference
# ----------------------------------------------------------------------
def reseau_deux_machines(slack, P1=None, V=None):
    """Deux machines et un GFM ; la machine slack est au choix."""
    net = Network(f0=50.0, Sb=100e6)
    t1, t2 = (1, 2) if slack == 1 else (2, 1)
    net.add_bus(1, 20e3, btype=t1, V=1.02 if V is None else V[0], Pg=P1 or 0.0, C=0.05)
    net.add_bus(2, 20e3, btype=t2, V=1.01 if V is None else V[1], Pg=0.5, C=0.05)
    net.add_bus(3, 20e3, btype=3, Pg=0.3, Pl=1.0, Ql=0.3, C=0.05)
    net.add(SynchronousMachine("SM1", 1, Sn=100e6, Sb=100e6, wb=net.wb))
    net.add(SynchronousMachine("SM2", 2, Sn=100e6, Sb=100e6, wb=net.wb))
    net.add(GFM1("GFM3", 3, Sn=100e6, Sb=100e6, wb=net.wb))
    net.add(Line("L13", 1, 3, R=0.01, X=0.1, wb=net.wb))
    net.add(Line("L23", 2, 3, R=0.01, X=0.1, wb=net.wb))
    net.add(RLLoad("Load3", 3, P=1.0, Q=0.3, wb=net.wb))
    net.run_loadflow()
    net.solve_equilibrium()
    return net


def test_loi_des_noeuds_avec_machine_de_reference():
    """Au bus de chaque machine, la somme des courants injectes alimente la capacite."""
    net = reseau_deux_machines(1)
    for b in (1, 2, 3):
        tot = 0j
        for el in net.elements:
            y = el.compiled().h(el.x0, el.z0, el.u0, el.p)
            for pk, bb in enumerate(el.buses):
                if bb == b:
                    a, c = el.port_output_idx(pk)
                    tot += complex(y[a], y[c])
        bus = net.buses[b]
        assert abs(tot - 1j * bus.pdict["C"] * complex(*bus.x0)) < 1e-9


def test_spectre_independant_de_la_machine_de_reference():
    """Changer de machine slack, au meme point de fonctionnement, ne change pas le spectre.

    C'est la signature d'un repere global coherent : le repere suit le rotor
    de la reference, et tous les elements -- machines comprises -- s'y
    reperent par theta - theta_g.
    """
    a = reseau_deux_machines(1)
    i = a._bus_index()
    V = (abs(a.lf.V[i[1]]), abs(a.lf.V[i[2]]))
    b = reseau_deux_machines(2, P1=a.lf.S[i[1]].real, V=V)
    la = np.linalg.eigvals(a.build_state_space().A)
    lb = np.linalg.eigvals(b.build_state_space().A)
    C = np.abs(la[:, None] - lb[None, :]) / np.maximum(1.0, np.abs(la[:, None]))
    r, c = linear_sum_assignment(C)
    assert C[r, c].max() < 1e-8


# ==========================================================================
# Assemblage a partir des tables (g2egrid.tables)
# ==========================================================================
def _tables_du_notebook(nom):
    """Execute les cellules de donnees d'un notebook de reseau et rend ses globales."""
    chemin = (
        pathlib.Path(__file__).resolve().parents[1] / "reseaux_test" / f"{nom}.ipynb"
    )
    nb = json.loads(chemin.read_text(encoding="utf-8"))
    code = [c for c in nb["cells"] if c["cell_type"] == "code"]
    g = {}
    for c in code[:3]:  # imports, donnees, assemblage
        exec(compile("".join(c["source"]), nom, "exec"), g)
    return g


@pytest.mark.parametrize("nom", ["smib", "wscc9"])
def test_assembleur_identique_au_notebook(nom):
    """``tables.assembler`` construit le meme reseau que la cellule des notebooks.

    Les notebooks gardent leur cellule d'assemblage, qui se lit et s'edite ;
    ``g2egrid.tables`` refait la meme chose comme fonction, pour construire des
    variantes par programme. Les deux chemins doivent coincider, sinon les
    etudes ne parlent pas du meme reseau que les notebooks.
    """
    g = _tables_du_notebook(nom)
    t = Tables(
        **{
            k: g[k]
            for k in (
                "Y_DER",
                "T_Lines_SI",
                "T_Lines_pu",
                "T_Shunts",
                "T_Transfos",
                "T_IBs",
                "T_SMs",
                "T_GFM1s",
                "T_GFM2s",
                "T_GFM3s",
                "T_GFM4s",
                "T_GFL1s",
                "T_GFL2s",
                "T_Loads",
            )
        },
        f0=g["f0"],
        Sb=g["Sb"],
        SM_PARAMS=g.get("SM_PARAMS", {}),
        GFM_PARAMS=g.get("GFM_PARAMS", {}),
        GFL_PARAMS=g.get("GFL_PARAMS", {}),
        # les reseaux de parite MATLAB fixent la resistance de Node.m
        **({"NOEUDS": g["NOEUDS"]} if "NOEUDS" in g else {}),
        **({"R_NOEUDS": float(g["R_NOEUDS"])} if "R_NOEUDS" in g else {}),
    )
    a, b = g["net"], assembler(t)
    assert [el.name for el in a.all_elements] == [el.name for el in b.all_elements]
    for reseau in (a, b):
        reseau.run_loadflow()
        reseau.solve_equilibrium()
    assert np.abs(a.build_state_space().A - b.build_state_space().A).max() < 1e-12


def test_params_par_bus_l_emporte_sur_la_famille():
    """``PARAMS_PAR_BUS`` regle un equipement seul, sans toucher a sa famille."""
    g = _tables_du_notebook("wscc9")
    t = Tables(
        Y_DER=g["Y_DER"],
        T_Lines_pu=g["T_Lines_pu"],
        T_Transfos=g["T_Transfos"],
        T_SMs=g["T_SMs"],
        f0=g["f0"],
        Sb=g["Sb"],
        SM_PARAMS={"H": 5.0},
        PARAMS_PAR_BUS={2: {"H": 9.0}},
    )
    H = {
        el.name: el.pdict["H"]
        for el in assembler(t).elements
        if el.name.startswith("SM")
    }
    assert H == {"SM1": 5.0, "SM2": 9.0, "SM3": 5.0}


# ==========================================================================
# Resolution non lineaire (g2egrid.temporel)
# ==========================================================================
def test_equilibre_est_un_point_fixe_du_reseau_entier():
    """Le systeme couple, resolu au point de linearisation, y reste.

    C'est la verification croisee entre les deux voies : le modele lineaire
    suppose que (x0, z0, u0) annule [f; g] ET la loi de cablage. Si les
    points de fonctionnement par element ne formaient pas un point fixe du
    reseau entier, ce serait invisible en modal (un jacobien est valable
    partout) mais fatal en temporel, qui partirait d'un transitoire parasite.
    """
    m = ModeleNonLineaire(reseau_deux_machines(1))
    dx, dzu = m.residu_equilibre()
    assert dx < 1e-5
    assert dzu < 1e-8


def test_non_lineaire_tend_vers_le_lineaire():
    """L'ecart entre les deux voies decroit comme le CARRE de l'amplitude.

    C'est la signature d'une linearisation correcte : la reponse du modele
    non lineaire a un echelon d'amplitude a vaut celle du modele lineaire
    plus un terme en a^2. Les deux chemins de code -- jacobiens symboliques
    plus complement de Schur d'un cote, equations non lineaires plus Newton
    couple de l'autre -- doivent donc converger l'un vers l'autre.
    """
    net = reseau_deux_machines(1)
    ss = net.build_state_space()
    m = ModeleNonLineaire(net)
    t = np.linspace(0.0, 1.0, 101)
    ecarts = []
    for amp in (0.002, 0.02):
        u = np.zeros(len(ss.inputs))
        u[ss.iu("SM1.p_ref")] = amp
        xl = reponse_indicielle(ss.A, ss.B, t, u)
        tr = m.simuler((0.0, 1.0), uext_fn=m.echelon("SM1.p_ref", amp), t_eval=t)
        ecarts.append(np.abs((tr.x - m.c.x0()) - xl).max())
    # amplitude x 10 -> ecart x 100 a un facteur 3 pres
    assert 30 < ecarts[1] / ecarts[0] < 300


def test_cablage_partage_par_les_deux_voies():
    """``cablage()`` ne depend pas du point de fonctionnement : le modele
    lineaire et le modele non lineaire partent bien du meme."""
    net = reseau_deux_machines(1)
    c = net.cablage()
    ss = net.build_state_space()
    assert c.state_names == ss.states
    assert c.uext_names == ss.inputs
    assert len(c.uext0()) == len(ss.inputs)
    assert c.Ruext.shape == (c.nu, len(c.uext_names))
    assert c.Ry.shape == (c.nu, c.ny)
