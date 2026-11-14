"""Faisceau creux ``(E, J)`` (``g2egrid.descripteur``) contre le modele d'etat
reduit : memes valeurs propres finies, memes participations, et couverture
complete d'un rectangle du plan complexe."""

import numpy as np
import pytest

from g2egrid import descripteur as dd
from g2egrid import modal
from g2egrid.catalogue import construire


@pytest.fixture(scope="module", params=["wscc9", "ieee39"])
def cas(request):
    net = construire(request.param)
    ss = net.build_state_space()
    return dd.descripteur(net), ss


def test_creux(cas):
    d, ss = cas
    assert d.nx == ss.A.shape[0] or d.etats != ss.states  # index 2 : etats elimines
    assert d.densite < 0.05


def test_modes_proches_et_participations(cas):
    d, ss = cas
    lam, _, p = modal.spectre(ss)
    m = dd.modes_proches(d, sigma=-0.5 + 2j * np.pi * 1.0, k=6)
    assert len(m) == 6
    for j, mu in enumerate(m.lam):
        i = int(np.argmin(np.abs(lam - mu)))
        assert abs(lam[i] - mu) < 1e-8 * max(abs(mu), 1.0)
        if d.etats == list(ss.states):
            np.testing.assert_allclose(m.p[j], p[i], atol=1e-6)
    # les six plus proches, dans l'ordre
    ref = np.sort(np.abs(lam - (-0.5 + 2j * np.pi)))[:6]
    np.testing.assert_allclose(np.abs(m.lam - (-0.5 + 2j * np.pi)), ref, rtol=1e-8)


def test_sigma_sur_le_mode_de_jauge(cas):
    """sigma = 0 est exactement une valeur propre (jauge) : pas de bruit."""
    d, ss = cas
    lam = np.linalg.eigvals(ss.A)
    m = dd.modes_proches(d, 0.0, k=8, gauche=False)
    for mu in m.lam:
        assert np.min(np.abs(lam - mu)) < 1e-8 * max(abs(mu), 1.0)


def test_balayage_complet(cas):
    d, ss = cas
    lam = np.linalg.eigvals(ss.A)
    m = dd.balayage(d, f_max=3.0, profondeur=0.5, gauche=False)
    dedans = (
        (lam.imag >= -1e-9) & (lam.imag / (2 * np.pi) <= 3) & (np.abs(lam.real) <= 0.5)
    )
    assert len(m) == dedans.sum()
    for mu in lam[dedans]:
        assert np.min(np.abs(m.lam - mu)) < 1e-6
    assert np.all(np.diff(m.lam.real) <= 0)


def test_k_borne(cas):
    d, _ = cas
    with pytest.raises(ValueError):
        dd.modes_proches(d, 0.0, k=d.nx)
