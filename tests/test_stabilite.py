"""Critere de stabilite ``StateSpace.is_stable`` : un seuil par valeur propre.

L'ancien critere comparait ``max Re(lambda)`` a ``max(tol, 1e-9 . max|A|)``.
Avec un bus de C = 1e-9 pu, ``max|A|`` atteint 3e11 et le seuil 314 s^-1 :
un mode instable de +61 s^-1 (ILoad) etait declare stable. Le critere actuel
(``marge_stabilite``) borne l'erreur de CHAQUE valeur propre par
``facteur . eps . ||A|| . kappa_i`` ; seul le mode de jauge (kappa = 1,
lambda = 0) bruite assez pour en avoir besoin.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
from test_g2egrid import build_charge

from g2egrid import Network
from g2egrid.catalogue import RESEAUX, charger
from g2egrid.elements import Line, PLoad, SynchronousMachine, ZLoad
from g2egrid.network import StateSpace
from g2egrid.tables import assembler


def ss_de(A: np.ndarray) -> StateSpace:
    n = A.shape[0]
    return StateSpace(
        A=A,
        B=np.zeros((n, 0)),
        C=np.zeros((0, n)),
        D=np.zeros((0, 0)),
        states=[f"x{i}" for i in range(n)],
        inputs=[],
        outputs=[],
    )


def matrice_a_jauge_bruitee(mode_extra: float = 0.0, graine: int = 3) -> np.ndarray:
    """Un A de norme 1e12 dont une valeur propre vaut 0 (ou ``mode_extra``).

    Bloc stable fortement non normal (couplage 1e12), mode decouple
    (conditionnement 1), le tout tourne par une orthogonale : ``eig`` rend
    alors le mode a quelques eps.||A|| ~ 1e-4, soit le bruit d'un mode de
    jauge sur un reseau a C = 1e-9 pu.
    """
    petit = np.array([[-1.0, 1e12], [0.0, -2.0]])
    A = np.zeros((3, 3))
    A[:2, :2] = petit
    A[2, 2] = mode_extra
    rng = np.random.default_rng(graine)
    Q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    return Q @ A @ Q.T


def test_jauge_bruitee_ignoree_au_dela_de_tol():
    """Le bruit du mode de jauge depasse ``tol`` (le test serait vide sinon),
    et le critere l'absorbe."""
    A = matrice_a_jauge_bruitee()
    assert np.linalg.eigvals(A).real.max() > 1e-6  # le bruit existe
    assert ss_de(A).is_stable()


def test_mode_instable_de_2_noye_dans_la_norme_est_vu():
    """+2 s^-1 pour ||A|| = 1e12 : l'ancien seuil (1e3) le masquait, celui du
    mode (4e-3 pour kappa = 1) le voit."""
    A = matrice_a_jauge_bruitee(mode_extra=2.0)
    ss = ss_de(A)
    ancien = np.linalg.eigvals(A).real.max() < max(1e-6, 1e-9 * np.abs(A).max())
    assert ancien  # TEMOIN : l'ancien critere declarait stable
    assert not ss.is_stable()
    assert ss.marge_stabilite() > 1.0


def test_mode_defectif_instable_est_vu():
    """Bloc de Jordan en +1 : kappa ~ 1e16 ne doit pas valoir excuse."""
    assert not ss_de(np.array([[1.0, 1.0], [0.0, 1.0]])).is_stable()


def test_marge_est_negative_et_continue_pour_un_modele_stable():
    ss = ss_de(np.diag([-1.0, -2.0, -3.0]))
    assert ss.is_stable() and ss.marge_stabilite() < 0


@pytest.fixture(scope="module")
def charge_c_1e_9():
    """Charge P, C = 1e-9 pu aux deux bus : ||A|| = 3e11, stable."""
    net = build_charge(PLoad)
    net.run_loadflow()
    net.solve_equilibrium()
    return net.build_state_space()


def test_decalage_de_61_declare_instable_malgre_la_norme(charge_c_1e_9):
    """Le cas ILoad / C = 1e-9 pu, SANS dependre d'ILoad (que le doctorant
    convertisseurs corrige) : le spectre du modele P est decale de +61 s^-1,
    l'instabilite mesuree sur ILoad. L'ancien seuil (314 s^-1) la masquait."""
    ss = charge_c_1e_9
    assert ss.is_stable()
    A = ss.A + 61.0 * np.eye(ss.A.shape[0])
    ancien = np.linalg.eigvals(A).real.max() < max(1e-6, 1e-9 * np.abs(A).max())
    assert ancien and np.abs(A).max() > 1e11  # TEMOIN
    instable = dataclasses.replace(ss, A=A)
    assert not instable.is_stable()
    assert instable.marge_stabilite() > 50.0


def reseau_sans_bus_infini(C: float) -> Network:
    """Deux machines, aucun IB : le mode de jauge existe (lambda = 0)."""
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=1, V=1.0, C=C)
    net.add_bus(2, 20e3, btype=2, V=1.0, Pg=0.4, C=C)
    net.add_bus(3, 20e3, btype=3, Pl=0.6, Ql=0.1, C=C)
    net.add(SynchronousMachine("SM1", 1, Sn=1e6, Sb=1e6, wb=net.wb))
    net.add(SynchronousMachine("SM2", 2, Sn=1e6, Sb=1e6, wb=net.wb))
    net.add(Line("L12", 1, 2, R=0.01, X=0.2, wb=net.wb))
    net.add(Line("L23", 2, 3, R=0.01, X=0.2, wb=net.wb))
    net.add(ZLoad("Z3", 3, P=0.6, Q=0.1, wb=net.wb))
    return net


@pytest.mark.parametrize("C", [1e-3, 1e-6, 1e-9])
def test_reseau_sans_bus_infini_stable_malgre_sa_jauge(C):
    net = reseau_sans_bus_infini(C)
    net.run_loadflow()
    net.solve_equilibrium()
    ss = net.build_state_space()
    w = np.linalg.eigvals(ss.A)
    assert np.min(np.abs(w)) < 1e-6  # il y a bien un mode de jauge
    assert ss.is_stable()


def _sans_continu(nom):
    t = charger(nom)
    return not (t.T_BusDC or t.T_VSCs or t.T_LinesDC)


#: ieee39_mixte (3 410 etats) et ieee14_ieee33 coutent des dizaines de secondes.
CATALOGUE_STABLE = [
    n
    for n in RESEAUX
    if _sans_continu(n) and n not in ("ieee39_mixte", "ieee14_ieee33")
]


@pytest.mark.parametrize("nom", CATALOGUE_STABLE)
def test_reseaux_stables_du_catalogue_restent_stables(nom):
    net = assembler(charger(nom))
    net.run_loadflow()
    net.solve_equilibrium()
    assert net.build_state_space().is_stable(), nom
