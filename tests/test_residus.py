"""Residus, grammiens et rapport de modes (``g2egrid.residus``), confrontes a
des references independantes : la fonction de transfert calculee par
resolution directe, et les equations de Lyapunov resolues par scipy."""

from types import SimpleNamespace

import numpy as np
import pytest
from scipy.linalg import solve_continuous_lyapunov

from g2egrid import residus as rs
from g2egrid.catalogue import construire


@pytest.fixture(scope="module")
def smib():
    """Le bus infini garde un angle libre : un mode nul aussi."""
    return construire("smib").build_state_space()


@pytest.fixture(scope="module")
def wscc9():
    """Sans bus infini : un mode de jauge nul a ecarter."""
    return construire("wscc9").build_state_space()


def test_residus_reconstruisent_le_transfert(smib):
    """H(s) - D = sum R_k / (s - lam_k), a une frequence quelconque."""
    s = 2j * np.pi * 1.3 + 0.2
    R = rs.residus(smib)
    lam = np.linalg.eigvals(smib.A)
    somme = np.einsum("kij,k->ij", R, 1 / (s - lam))
    n = smib.A.shape[0]
    direct = smib.C @ np.linalg.solve(s * np.eye(n) - smib.A, smib.B)
    np.testing.assert_allclose(
        somme, direct, rtol=1e-6, atol=1e-8 * np.abs(direct).max()
    )


def _aleatoire(n=12, nu=3, ny=4, graine=0):
    """Un systeme stable quelconque, valeurs propres complexes comprises."""
    rng = np.random.default_rng(graine)
    A = rng.normal(size=(n, n))
    A -= (np.linalg.eigvals(A).real.max() + 0.5) * np.eye(n)
    return SimpleNamespace(
        A=A,
        B=rng.normal(size=(n, nu)),
        C=rng.normal(size=(ny, n)),
        D=np.zeros((ny, nu)),
        states=[f"x{i}" for i in range(n)],
        inputs=[f"u{i}" for i in range(nu)],
        outputs=[f"y{i}" for i in range(ny)],
    )


def test_grammiens_lyapunov():
    ss = _aleatoire()
    g = rs.grammiens(ss)
    assert g.exclus.size == 0
    Wc, Wo = g.physique()
    A, B, C = ss.A, ss.B, ss.C
    Wc_ref = solve_continuous_lyapunov(A, -B @ B.T)
    Wo_ref = solve_continuous_lyapunov(A.T, -C.T @ C)
    np.testing.assert_allclose(Wc, Wc_ref, atol=1e-9 * np.abs(Wc_ref).max())
    np.testing.assert_allclose(Wo, Wo_ref, atol=1e-9 * np.abs(Wo_ref).max())
    h_ref = np.sort(np.sqrt(np.abs(np.linalg.eigvals(Wc_ref @ Wo_ref))))[::-1]
    np.testing.assert_allclose(g.hankel, h_ref, rtol=1e-7)
    assert g.ordre(0.0) == A.shape[0]


@pytest.mark.parametrize("nom", ["smib", "wscc9"])
def test_grammiens_ecartent_la_jauge(nom, smib, wscc9):
    g = rs.grammiens({"smib": smib, "wscc9": wscc9}[nom])
    assert g.exclus.size >= 1
    assert np.all(np.abs(g.exclus) < 1e-3)
    assert np.all(np.isfinite(g.hankel))


def test_rapport_electromecanique(wscc9):
    """Les modes les moins amortis de wscc9 entre 0,5 et 3 Hz sont
    electromecaniques : vitesses rotor, commandables par les consignes des
    machines."""
    r = rs.rapport(wscc9, n_modes=2, f_min=0.5, f_max=3.0)
    assert len(r) == 2
    for m in r:
        assert m.stable and m.valeur.imag > 0
        assert m.elements[0][0].startswith("SM")
        assert m.entrees[0][1] == pytest.approx(100.0)
        assert m.sorties[0][1] == pytest.approx(100.0)
        assert abs(sum(p for _, p in m.elements) - 100) < 1e-6
        assert "Mode" in m.texte()
    assert r[0].zeta <= r[1].zeta


def test_dominance_classement(wscc9):
    r = rs.rapport(wscc9, n_modes=3, critere="dominance", f_min=0.1)
    d = [m.dominance for m in r]
    assert d == sorted(d, reverse=True)
    with pytest.raises(ValueError):
        rs.rapport(wscc9, critere="inconnu")
