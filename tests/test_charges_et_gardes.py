"""Garde-fous des charges et calage du bus continu des GFM.

* ``RLLoad``/``ZLoad`` appelees directement : Q <= 0 est refuse par une
  erreur claire (reactance nulle ou negative), comme ``tables.py`` ;
* ``ILoad`` avec Q != 0 : spectre STABLE (le mode haute frequence faussement
  instable venait d'une admittance incrementale non passive) ;
"""

from __future__ import annotations

import numpy as np
import pytest

from g2egrid import Network
from g2egrid.elements import (
    IB,
    ILoad,
    PLoad,
    RLLoad,
    Transformer,
    ZLoad,
)


# ----------------------------------------------------------------------
# RLLoad / ZLoad : Q <= 0 refuse
# ----------------------------------------------------------------------
@pytest.mark.parametrize("cls", [RLLoad, ZLoad])
@pytest.mark.parametrize("Q", [0.0, -0.1])
def test_charge_rl_refuse_un_reactif_nul_ou_negatif(cls, Q):
    with pytest.raises(ValueError, match=r"Q > 0.*PLoad ou ILoad"):
        cls("Load", 2, P=0.5, Q=Q, wb=314.159)


@pytest.mark.parametrize("cls", [RLLoad, ZLoad])
def test_charge_rl_accepte_un_reactif_positif(cls):
    el = cls("Load", 2, P=0.5, Q=0.1, wb=314.159)
    assert el.Q == 0.1


# ----------------------------------------------------------------------
# ILoad : spectre stable avec Q != 0
# ----------------------------------------------------------------------
def _reseau_charge(cls, P, Q, C, X):
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 110e3, btype=1, V=1.0, C=C)
    net.add_bus(2, 20e3, btype=3, Pl=P, Ql=Q, C=C)
    net.add(Transformer("Tr1", 1, 2, R=0.0002, X=X, wb=net.wb))
    net.add(cls("Load2", 2, P=P, Q=Q))
    net.add(IB("GRID", 1, R=1e-4, X=1e-3, wb=net.wb))
    net.run_loadflow()
    net.solve_equilibrium()
    return net


@pytest.mark.parametrize("C", [1e-9, 1e-6, 1e-3, 1e-2])
@pytest.mark.parametrize("X", [0.024, 0.2])
@pytest.mark.parametrize("P,Q", [(0.5, 0.1), (0.0, 0.3), (0.5, -0.2), (0.5, 0.0)])
def test_iload_spectre_stable(C, X, P, Q):
    """Aucune valeur propre a partie reelle positive : comparer a ``0``, pas
    a ``is_stable()`` dont le seuil (>= 1e-9 max|A|, 100 a 300 /s avec un bus
    de 1e-9) masquait l'instabilite de +61 /s. Avant la correction : Re de +2
    a +6000 /s des que Q != 0."""
    ev = _reseau_charge(ILoad, P, Q, C, X).build_state_space().eig()
    assert ev.real.max() < 1e-6


def test_cause_admittance_incrementale_passive():
    """La cause, lue sur le jacobien : l'admittance incrementale ``dI/dV`` de
    la charge (``Tf.B`` sur les entrees de tension) doit avoir une partie
    symetrique semi-definie positive (passive). Celle de l'ancien
    ``courant constant, facteur de puissance constant`` avait pour determinant
    ``-Q^2/4`` (une direction active) ; ``PLoad`` est symetrique."""
    for cls in (ILoad, PLoad):
        net = _reseau_charge(cls, 0.5, 0.1, 1e-6, 0.024)
        el = next(e for e in net.elements if e.name == "Load2")
        cm = el.compiled()
        _, B, _, _ = cm.linearize(el.x0, el.z0, el.u0, el.p)
        Y = el.pdict["Tf"] * B[:, [cm.iu("vd_t"), cm.iu("vq_t")]]
        if cls is PLoad:  # symetrique, de trace nulle : traverse le zero
            assert np.allclose(Y, Y.T, atol=1e-3)
        else:
            sym = (Y + Y.T) / 2
            assert np.linalg.eigvalsh(sym).min() > -1e-9


def test_iload_equilibre_et_sensibilite_statique_inchanges():
    """P(|V|) = P0.|V|/V0 (exposant 1), consigne tenue : ce que l'ILoad
    garantissait avant, y compris a Q != 0."""
    net = _reseau_charge(ILoad, 0.5, 0.1, 1e-9, 0.024)
    el = next(e for e in net.elements if e.name == "Load2")
    cm = el.compiled()
    assert np.abs(cm.residual(el.x0, el.z0, el.u0, el.p)).max() < 1e-10
    vd, vq = el.u0[cm.iu("vd_t")], el.u0[cm.iu("vq_t")]
    V = np.hypot(vd, vq)
    A, B, C, D = cm.linearize(el.x0, el.z0, el.u0, el.p)
    G = D - C @ np.linalg.solve(A, B)
    ip, iq = cm.iy("p"), cm.iy("q")
    # derivee statique le long de |V| (angle fixe)
    dP = G[ip, cm.iu("vd_t")] * vd / V + G[ip, cm.iu("vq_t")] * vq / V
    dQ = G[iq, cm.iu("vd_t")] * vd / V + G[iq, cm.iu("vq_t")] * vq / V
    assert vd * el.x0[0] + vq * el.x0[1] == pytest.approx(0.5, rel=1e-9)
    assert dP == pytest.approx(0.5 / V, abs=1e-9)  # exposant 1
    assert dQ == pytest.approx(0.1 / V, abs=1e-9)


def test_iload_q_nul_identique_a_l_ancien_modele():
    """A Q0 = 0 la correction ne change rien, QUEL QUE SOIT l'angle de repli
    ``phi0`` : les equations sont celles de ``ILoad`` d'avant,
    ``Tf.di/dt = (|v|/V0).conj(S/v) - i``, evaluees en un point quelconque."""
    net = _reseau_charge(ILoad, 0.5, 0.0, 1e-9, 0.024)
    el = next(e for e in net.elements if e.name == "Load2")
    cm = el.compiled()
    p = np.array(el.p, dtype=float)
    p[el.params.index("Q0")] = 0.0
    p[el.params.index("cphi0")], p[el.params.index("sphi0")] = 0.6, 0.8  # autre phi0
    x = np.array([0.31, -0.07])
    u = np.array(el.u0, dtype=float)
    u[cm.iu("dp")] = u[cm.iu("dq")] = 0.0
    u[cm.iu("vd_t")], u[cm.iu("vq_t")] = 0.93, 0.21
    v = complex(u[cm.iu("vd_t")], u[cm.iu("vq_t")])
    i_ref = abs(v) / p[el.params.index("V0")] * np.conj(0.5 / v)
    attendu = (
        np.array([i_ref.real - x[0], i_ref.imag - x[1]]) / p[el.params.index("Tf")]
    )
    f = cm.f(x, el.z0, u, p)
    assert np.allclose(f, attendu, rtol=1e-12, atol=1e-12)
