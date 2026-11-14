"""Tests de l'alternateur a ligne d'arbre multi-masses.

References calculables a la main ou publiees :

* frequences de torsion a rotor libre du IEEE First Benchmark Model (IEEE SSR
  Working Group, 1977) : 15,71 / 20,21 / 25,55 / 32,28 / 47,46 Hz ;
* frequences des groupes hydrauliques du PNNL (Chatterjee et al., arXiv
  2607.14412, tableau 1) : 17 / 8,14 / 6,26 Hz, a partir des inerties et
  raideurs physiques ;
* contraintes moyennes du tableau 2 du meme papier (arbres creux) ;
* a l'equilibre, chaque troncon transmet la somme des couples moteurs des
  etages en amont (statique des arbres).
"""

from dataclasses import replace

import numpy as np
import pytest
from scipy.optimize import root

from g2egrid import Network
from g2egrid.elements import (
    ARBRE_FBM,
    IB,
    TAU_TORSION,
    Line,
    Transformer,
    ligne_deux_masses,
    turbo_alternateur,
)
from g2egrid.fatigue import (
    couple_base,
    limite_non_resonante_epri,
    securite_goodman,
    tension_cisaillement,
)

WB60 = 2.0 * np.pi * 60.0


# ----------------------------------------------------------------------
# ligne d'arbre seule
# ----------------------------------------------------------------------
def test_fbm_frequences_libres_publiees():
    f = ARBRE_FBM.frequences_libres(WB60)
    assert f[0] == pytest.approx(0.0, abs=1e-5)  # mode rigide
    assert f[1:] == pytest.approx([15.71, 20.21, 25.55, 32.28, 47.46], abs=6e-3)


def test_amortissement_proportionnel_a_l_inertie():
    """``D_r = 4.H_r/tau`` : chaque mode a rotor libre decroit en
    ``exp(-t/tau)``, le mode rigide en ``exp(-2t/tau)``."""
    g = ARBRE_FBM.masses.index("GEN")
    A, _, _ = ARBRE_FBM.modele_mecanique(WB60, D_gen=4.0 * ARBRE_FBM.H[g] / TAU_TORSION)
    lam = np.linalg.eigvals(A)
    osc = lam[np.abs(lam.imag) > 1.0]
    assert len(osc) == 10
    assert osc.real == pytest.approx(-1.0 / TAU_TORSION, rel=1e-9)
    rigide = lam[np.abs(lam.imag) <= 1.0]
    assert rigide.real == pytest.approx([-2.0 / TAU_TORSION], rel=1e-9)
    # et les frequences amorties restent celles de la ligne sans amortissement
    f = np.sort(osc.imag[osc.imag > 0]) / (2 * np.pi)
    w0 = 2 * np.pi * ARBRE_FBM.frequences_libres(WB60)[1:]
    assert f == pytest.approx(np.sqrt(w0**2 - TAU_TORSION**-2) / (2 * np.pi), rel=1e-9)


def test_modele_mecanique_gain_statique():
    """Ligne libre sans amortissement, couple electrique unite : toute la ligne
    decelere d'un bloc, ``dDw/dt = -1/(2.sum(H))``, et chaque troncon transmet
    la part du couple qui freine les masses en amont,
    ``C_i = sum(H_amont)/sum(H) - [troncon en aval de l'alternateur]``."""
    lg = replace(ARBRE_FBM, D=(0.0,) * 6)
    A, B, C = lg.modele_mecanique(WB60)
    H = np.array(lg.H)
    n, g = 6, lg.masses.index("GEN")
    a = -1.0 / (2.0 * H.sum())
    # regime d'acceleration uniforme : les torsions resolvent A.x + B = (a, ..., a)
    d = np.linalg.lstsq(A[:n, n:], np.full(n, a) - B[:n, 0], rcond=None)[0]
    assert A[:n, n:] @ d + B[:n, 0] == pytest.approx(np.full(n, a), rel=1e-12)
    attendu = [H[: i + 1].sum() / H.sum() - (i >= g) for i in range(n - 1)]
    assert C[:, n:] @ d == pytest.approx(attendu, rel=1e-12)


@pytest.mark.parametrize(
    "J_t, J_g, K, Sn, n, f_pub, tol",
    [
        # Unit A : le tableau arrondit a 17 Hz (17,2 par la formule)
        (1.68e6, 33.75e6, 18.66e9, 500e6, 150.0, 17.0, 0.25),
        (90.6e3, 0.653e6, 0.208e9, 60e6, 200.0, 8.14, 0.01),  # Grand Coulee GC-7
        (8428e3, 107e6, 12.09e9, 715e6, 72.0, 6.26, 0.01),  # Grand Coulee GC-19
    ],
)
def test_groupes_hydrauliques_pnnl(J_t, J_g, K, Sn, n, f_pub, tol):
    """La conversion en pu (radians electriques) redonne la frequence de la
    formule physique ``sqrt(K.(1/Jt + 1/Jg))/(2.pi)`` et celle du tableau."""
    ligne, Hg = ligne_deux_masses(J_t, J_g, K, 0.0, Sn, n, f0=60.0)
    f = ligne.frequences_libres(WB60)[1]
    f_phys = np.sqrt(K * (1.0 / J_t + 1.0 / J_g)) / (2.0 * np.pi)
    assert f == pytest.approx(f_phys, rel=1e-12)
    assert f == pytest.approx(f_pub, abs=tol)
    assert Hg == ligne.H[1]


# ----------------------------------------------------------------------
# contraintes (arbres creux du PNNL) et criteres simplifies
# ----------------------------------------------------------------------
@pytest.mark.parametrize(
    "Sn, n, T_pu, D, d, tau_m_pub",
    [
        (60e6, 200.0, 0.692, 0.686, 0.160, 31.4e6),  # GC-7
        (715e6, 72.0, 0.684, 2.540, 2.134, 40.2e6),  # GC-19
    ],
)
def test_contrainte_moyenne_pnnl(Sn, n, T_pu, D, d, tau_m_pub):
    paires = round(60.0 * 60.0 / n)
    tau = tension_cisaillement(T_pu * couple_base(Sn, 60.0, paires), D, d)
    assert tau == pytest.approx(tau_m_pub, abs=0.06e6)


def test_goodman_limites():
    """Sans alterne, ``n = S_su/tau_m`` ; sans moyenne, ``n = S_se/tau_a``."""
    assert securite_goodman(0.0, 40e6, 35e6, 585e6) == pytest.approx(585 / 40)
    assert securite_goodman(7e6, 0.0, 35e6, 585e6) == pytest.approx(5.0)
    assert securite_goodman(3.5e6, 29.25e6, 35e6, 585e6) == pytest.approx(
        1.0 / (0.1 + 0.05)
    )


def test_methode_simplifiee_epri_20_pourcent():
    """Annexe A du rapport EPRI 3002036374 : acier de 80 ksi, Kt = 4, SF = 2,
    tau_m = 50 MPa => environ 20 % de la puissance nominale."""
    assert limite_non_resonante_epri(80e3 * 6894.757) == pytest.approx(0.20, abs=0.005)
    # proportionnelle a l'UTS, inversement a Kt et au coefficient de securite
    assert limite_non_resonante_epri(552e6, Kt=2.0) == pytest.approx(
        2.0 * limite_non_resonante_epri(552e6)
    )


# ----------------------------------------------------------------------
# machine raccordee : reseau machine - transformateur - ligne - reseau infini
# (donnees du FBM, base machine 892,4 MVA)
# ----------------------------------------------------------------------
def _smib(**params):
    Sb = 892.4e6
    net = Network(f0=60.0, Sb=Sb)
    wb = net.wb
    net.add_bus(1, 26e3, btype=2, V=1.0, Pg=0.9, C=1e-2)
    net.add_bus(2, 500e3, btype=3, C=1e-2)
    net.add_bus(3, 500e3, btype=1, V=1.0, C=1e-2)
    net.add(turbo_alternateur("G", 1, Sn=Sb, Sb=Sb, wb=wb, **params))
    net.add(Transformer("TR", 1, 2, R=0.002, X=0.14, wb=wb))
    net.add(Line("L", 2, 3, R=0.02, X=0.50, wb=wb))
    net.add(IB("GRID", 3, R=1e-3, X=0.06, wb=wb))
    net.run_loadflow()
    net.solve_equilibrium()
    return net


@pytest.fixture(scope="module")
def smib():
    """Machine sans PSS, regulateur de vitesse fige (``Tg`` = 1000 s) : la ligne
    d'arbre et le reseau seuls. Le mode electromecanique (1,2 Hz), que le
    regulateur de vitesse raide de ``SM_DEFAULT`` amortissait, est alors
    instable : seuls les modes de torsion sont examines sur ce reseau."""
    return _smib(K_PSS=0.0, Tg=1e3)


def _modes(net):
    lam = np.linalg.eigvals(net.build_state_space().A)
    m = [x for x in lam if 2 * np.pi * 10 < x.imag < 2 * np.pi * 55 and x.real > -5]
    return sorted(m, key=lambda x: x.imag)


def _sortie(net, nom, sortie_):
    """Une sortie de service d'un element, evaluee a l'equilibre."""
    el = next(e for e in net.elements if e.name == nom)
    y = el.compiled().h(el.x0, el.z0, el.u0, el.p)
    return float(y[el.n_port_signaux() + el.svc_outputs.index(sortie_)])


def test_equilibre_couples_de_troncon(smib):
    """Chaque troncon transmet la somme des couples moteurs en amont
    (fractions 0,30 / 0,26 / 0,22 / 0,22), le dernier (excitatrice) rien."""
    G = next(e for e in smib.elements if e.name == "G")
    Pm = G.x0[G.states.index("Pm")]  # a vitesse nominale : couple = puissance
    attendu = np.cumsum([0.30, 0.26, 0.22, 0.22]) * Pm
    obtenu = [_sortie(smib, "G", f"C_{t}") for t in G.troncons]
    assert obtenu[:4] == pytest.approx(attendu, rel=1e-9)
    assert obtenu[4] == pytest.approx(0.0, abs=1e-12)
    # le couple sur l'accouplement turbine - alternateur est le couple
    # electrique d'entrefer, soit la puissance aux bornes (0,9 pu du load flow)
    # augmentee des pertes statoriques (Ra) et du shunt de causalite (Rg)
    phi_d, phi_q = (G.x0[G.states.index(n)] for n in ("phi_d", "phi_q"))
    i_d, i_q = (G.z0[G.algebs.index(n)] for n in ("i_d", "i_q"))
    assert obtenu[3] == pytest.approx(phi_d * i_q - phi_q * i_d, rel=1e-12)
    assert _sortie(smib, "G", "P") == pytest.approx(0.9, rel=1e-9)
    assert 0.0 < obtenu[3] - 0.9 < 0.02


def test_jacobien_multi_masses(smib):
    el = next(e for e in smib.elements if e.name == "G")
    cm = el.compiled()
    A, _, _, _ = cm.linearize(el.x0, el.z0, el.u0, el.p)

    def dxdt(x):
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
    assert err < 1e-9, f"erreur jacobien {err:.2e}"


def test_modes_de_torsion_raccordes(smib):
    """Machine raccordee, sans regulation : cinq modes de torsion a 0,1 Hz pres
    des frequences a rotor libre (la raideur electrique ne pese que sur
    l'alternateur) et de celles que l'EPRI donne pour le FBM (rapport
    3002036374, p. 31 : 15,7 / 20,2 / 25,6 / 32,3 / 47,5 Hz). Sur une ligne sans
    compensation serie, le reseau et les amortisseurs de l'alternateur ne
    changent presque pas l'amortissement : chaque mode garde ``-1/tau`` a 10 %
    pres."""
    torsion = _modes(smib)
    f = [x.imag / (2 * np.pi) for x in torsion]
    assert f == pytest.approx(list(ARBRE_FBM.frequences_libres(smib.wb)[1:]), abs=0.1)
    assert f == pytest.approx([15.7, 20.2, 25.6, 32.3, 47.5], abs=0.1)
    assert [x.real for x in torsion] == pytest.approx([-1.0 / TAU_TORSION] * 5, rel=0.1)


def test_pss_vitesse_desamortit_le_premier_mode_de_torsion():
    """Interaction de torsion connue (Kundur, ch. 15) : un PSS a entree vitesse
    sans filtre de torsion reinjecte la vitesse de l'alternateur, ou les modes de
    torsion sont visibles, et retire de l'amortissement au premier mode
    (15,6 Hz). Avec les regulations de ``SM_DEFAULT`` (PSS K = 20, filtre du
    premier ordre a 5 Hz), il lui en retire plus de la moitie ; la machine
    reste stable."""
    sans_net = _smib(K_PSS=0.0)
    avec_net = _smib()
    assert sans_net.build_state_space().is_stable()
    assert avec_net.build_state_space().is_stable()
    sans, avec = _modes(sans_net)[0], _modes(avec_net)[0]
    assert abs(avec.imag - sans.imag) / (2 * np.pi) < 0.2
    assert avec.real > 0.5 * sans.real
