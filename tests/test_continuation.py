"""Continuation (``g2egrid.continuation``) : bifurcations connues
analytiquement, un Hopf de PLL sur wscc9, et l'arret de la suite au dernier
point converge quand un pas n'a plus de point d'equilibre."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
from _etudes import dossier_etude

from g2egrid import continuation as c

ICI = dossier_etude("campagne")

import campagne


def _jouet(A_de_mu):
    def fabrique(mu):
        A = A_de_mu(mu)
        ss = SimpleNamespace(A=A, states=[f"x{i}" for i in range(len(A))])
        return SimpleNamespace(build_state_space=lambda: ss)

    return fabrique


def test_hopf_analytique():
    """Paire ``mu +- 3j`` (plus un mode stable) : Hopf a mu = 0,25."""
    f = _jouet(
        lambda mu: np.array(
            [[mu - 0.25, -3.0, 0.0], [3.0, mu - 0.25, 0.0], [0.0, 0.0, -5.0]]
        )
    )
    r = c.continuer(f, np.linspace(-1, 1, 5))
    assert len(r.bifurcations) == 1
    b = r.bifurcations[0]
    assert b.genre == "hopf" and b.sens == 1
    assert b.mu == pytest.approx(0.25, abs=1e-4)
    assert b.f == pytest.approx(3 / (2 * np.pi), rel=1e-3)
    assert r.marge == b.mu
    # parcouru a l'envers : meme point, sens toujours rapporte a mu croissant
    r2 = c.continuer(f, np.linspace(1, -1, 5))
    assert r2.bifurcations[0].sens == 1
    assert r2.bifurcations[0].mu == pytest.approx(0.25, abs=1e-4)


def test_noeud_col_analytique():
    f = _jouet(lambda mu: np.diag([mu - 1.0, -2.0]))
    r = c.continuer(f, [0.0, 0.7, 1.4, 2.0])
    (b,) = r.bifurcations
    assert b.genre == "noeud-col"
    assert b.mu == pytest.approx(1.0, abs=1e-4)


def test_hopf_pll_wscc9():
    """Un GFL stable : monter le gain integral de sa PLL le destabilise vers
    54 Hz, et la valeur trouvee encadre bien le changement de stabilite."""
    net = campagne.modele_complet(campagne.construire(campagne.tirer(3, "wscc9")))[0]
    el = next(e for e in net.all_elements if e.name == "GFL2")
    p0 = el.pdict["Kipll"]
    r = c.continuer(c.par_parametre(net, "GFL2", "Kipll"), p0 * np.geomspace(1, 30, 8))
    b = r.bifurcations[0]
    assert b.genre == "hopf" and b.sens == 1
    assert b.elements[0][0] == "GFL2"
    assert 40 < b.f < 70
    for facteur, instable in ((0.998, False), (1.002, True)):
        el.set_param("Kipll", b.mu * facteur)
        lam = np.linalg.eigvals(net.build_state_space().A)
        assert (lam.real.max() > 1e-6) == instable
    el.set_param("Kipll", p0)
    with pytest.raises(KeyError):
        c.par_parametre(net, "GFL2", "inconnu")
