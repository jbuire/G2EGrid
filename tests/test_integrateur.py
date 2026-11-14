"""Integrateur DAE monolithique (BDF2), evaluation en lot, pseudo-transitoire."""

from __future__ import annotations

import numpy as np
import pytest

import g2egrid.catalogue as cat
from g2egrid.integrateur import (
    equilibre_pseudo_transitoire,
    jacobien_dae,
    residu_dae,
)
from g2egrid.temporel import ModeleNonLineaire


@pytest.fixture(scope="module")
def wscc9():
    return ModeleNonLineaire(cat.construire("wscc9"))


@pytest.mark.parametrize("nom", ["wscc9", "hvdc2t"])
def test_jacobien_dae_contre_differences_finies(nom):
    m = ModeleNonLineaire(cat.construire(nom))
    c = m.c
    rng = np.random.default_rng(0)
    x = c.x0() * (1 + 0.02 * rng.standard_normal(c.nx))
    z = c.z0() * (1 + 0.02 * rng.standard_normal(c.nz))
    u = c.u0() * (1 + 0.02 * rng.standard_normal(c.nu))
    uext, a, hist = m.uext0(), 3000.0, rng.standard_normal(c.nx)

    def F(y):
        return residu_dae(
            m,
            y[: c.nx],
            y[c.nx : c.nx + c.nz],
            y[c.nx + c.nz :],
            uext,
            a * y[: c.nx] + hist,
        )

    y = np.concatenate([x, z, u])
    J = jacobien_dae(m, x, z, u, a).toarray()
    Jn = np.zeros_like(J)
    for j in range(len(y)):
        e = 1e-7 * max(1.0, abs(y[j]))
        yp, ym = y.copy(), y.copy()
        yp[j] += e
        ym[j] -= e
        Jn[:, j] = (F(yp) - F(ym)) / (2 * e)
    assert np.abs(J - Jn).max() / np.abs(Jn).max() < 1e-7


def test_residu_en_lot_egale_la_voie_element_par_element(wscc9):
    """L'evaluation groupee par classe rend exactement ce que rend chaque
    element pris a part."""
    m, c = wscc9, wscc9.c
    rng = np.random.default_rng(1)
    x = c.x0() * (1 + 0.01 * rng.standard_normal(c.nx))
    z, u = c.z0(), c.u0()
    r = residu_dae(m, x, z, u, m.uext0(), np.zeros(c.nx))
    f = np.concatenate(
        [
            cm.f(xi, zi, ui, pi)
            for cm, pi, xi, zi, ui in zip(c.cms, m.p, m._x(x), m._z(z), m._u(u))
            if cm.n_x
        ]
    )
    assert np.allclose(r[: c.nx], -f, rtol=1e-12, atol=1e-12)
    # et a l'equilibre, tout le residu est nul
    r0 = residu_dae(m, c.x0(), z, u, m.uext0(), np.zeros(c.nx))
    assert np.abs(r0).max() < 1e-6


def test_bdf2_rejoint_radau(wscc9):
    """Echelon de consigne sur wscc9 : BDF2 a 200 us colle a Radau serre."""
    m = wscc9
    te = np.linspace(0.0, 0.4, 81)
    u = m.echelon("SM1.p_ref", 0.01, 0.05)
    ref = m.simuler((0, 0.4), uext_fn=u, t_eval=te, rtol=1e-9, atol=1e-11)
    b = m.simuler((0, 0.4), uext_fn=u, t_eval=te, methode="BDF2", pas=2e-4)
    exc = np.abs(ref.x - ref.x[:1]).max(axis=0)
    ecart = np.abs(b.x - ref.x).max(axis=0) / np.maximum(exc, 1e-12)
    assert ecart.max() < 1e-3


def test_pseudo_transitoire_retrouve_l_equilibre(wscc9):
    """D'un point faux de 1 %, la continuation revient sur l'equilibre."""
    m = wscc9
    x_ref = m.c.x0()
    x, _z, _u, r = equilibre_pseudo_transitoire(m, x0=x_ref * 1.01 + 1e-3)
    assert r < 1e-9
    angle = np.array(["theta" in n for n in m.c.state_names])
    assert np.abs(x - x_ref)[~angle].max() < 1e-8


def test_solve_equilibrium_pseudo_transitoire():
    net = cat.construire("wscc9_gfm")
    net.solve_equilibrium(pseudo_transitoire=True)
    dx, dzu = ModeleNonLineaire(net).residu_equilibre()
    assert dx < 1e-8 and dzu < 1e-9
