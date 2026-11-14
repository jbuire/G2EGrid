"""Le critere de la boucle mineure (:mod:`g2egrid.impedance`) : le Nyquist
generalise sur le gain de boucle equipement/reseau doit dire la meme chose
que le modele d'etat complet."""

from __future__ import annotations

import numpy as np
import pytest
from _etudes import dossier_etude
from g2egrid.impedance import bandes, boucle_mineure, impedance_bus, passivite

ICI = dossier_etude("campagne")

import campagne


@pytest.fixture(scope="module")
def instable():
    """wscc9, tirage 28 : GFL2 rend le reseau instable vers 78 Hz."""
    t = campagne.construire(campagne.tirer(28, "wscc9"))
    return campagne.modele_complet(t)[0]


def test_gfl_instable_detecte(instable):
    b = boucle_mineure(instable, "GFL2")
    assert b.forme == "courant"
    assert b.instables == 2
    assert b.applicable
    assert b.nyquist["tours"] == -2
    assert not b.stable


def test_impedance_bus_poles_du_complet(instable):
    """Les poles de Z_t sont les modes du reseau complet : |Z_t| culmine pres
    de la frequence du mode instable."""
    lam = np.linalg.eigvals(instable.build_state_space().A)
    lam = lam[lam.real > 0]
    f0 = abs(lam[0].imag) / (2 * np.pi)
    f = np.array([f0 / 3, f0, 3 * f0])
    bus = next(e for e in instable.all_elements if e.name == "GFL2").buses[0]
    Z = np.abs(impedance_bus(instable, bus, f)).max(axis=(1, 2))
    assert Z[1] > Z[0] and Z[1] > Z[2]


def test_reference_refusee(instable):
    with pytest.raises(ValueError):
        boucle_mineure(instable, instable.reference)


def test_passivite_et_bandes():
    f = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    Y = np.zeros((5, 2, 2))
    Y[:, 0, 0] = Y[:, 1, 1] = [1.0, -1.0, -1.0, 1.0, -1.0]
    p = passivite(Y)
    np.testing.assert_allclose(p, [1, -1, -1, 1, -1])
    assert bandes(f, p < 0) == [(2.0, 3.0), (5.0, 5.0)]


def test_gfl_stable_aucun_tour():
    t = campagne.construire(campagne.tirer(3, "wscc9"))
    b = boucle_mineure(campagne.modele_complet(t)[0], "GFL2")
    assert b.instables == 0 and b.applicable and b.stable
    assert b.nyquist["marge"] > 0


def test_comptage_independant_de_la_grille():
    """smib_gfl est stable : 0 tour, quelle que soit la grille. Un mode du
    reste du reseau a 16 kHz, amorti sur 0,3 Hz, et une resonance au-dela de
    20 kHz faisaient compter -2 tours (P = -2) a la grille par defaut."""
    from g2egrid import catalogue

    net = catalogue.construire("smib_gfl")
    for f in (np.geomspace(0.01, 20000, 2400), np.geomspace(0.1, 5000, 300)):
        b = boucle_mineure(net, "GFL2", f)
        assert b.nyquist["tours"] == 0 and b.coherent and b.stable
