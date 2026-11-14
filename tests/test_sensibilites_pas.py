"""Pas des differences finies de ``sensibilites`` : relatif a |p0|, sans changer de signe.

Avant : ``h = pas . max(|p0|, 1)``, soit 1 % de p0 = 1e-4 avec ``pas = 1e-6``
(erreur de troncature en h^2 = 1e-4 relatif, au lieu de 1e-10), et ``p0 - h``
pouvait devenir negatif (resistance, constante d'inertie...), c'est-a-dire
sortir du domaine physique du modele.
"""

from __future__ import annotations

import numpy as np
import pytest
from test_g2egrid import build_smib

from g2egrid.sensibilites import sensibilites


def valeurs_essayees(net, mode, **kw):
    """``{(element, parametre): [valeurs posees par set_param]}`` pendant le calcul."""
    vus: dict[tuple[str, str], list[float]] = {}
    classes = {type(e) for e in net.all_elements}
    origines = {}
    for cls in classes:
        origines[cls] = cls.set_param

        def espion(self, nom, valeur, _o=cls.set_param):
            vus.setdefault((self.name, nom), []).append(float(valeur))
            return _o(self, nom, valeur)

        cls.set_param = espion
    try:
        res = sensibilites(net, mode, methode="differences", **kw)
    finally:
        for cls, o in origines.items():
            cls.set_param = o
    return res, vus


@pytest.fixture
def smib():
    net = build_smib()
    net.run_loadflow()
    net.solve_equilibrium()
    return net


def test_pas_relatif_et_jamais_de_changement_de_signe(smib):
    _, vus = valeurs_essayees(smib, 0)
    p0 = {
        (e.name, n): float(e.p[j])
        for e in smib.all_elements
        if e.x0 is not None
        for j, n in enumerate(e.params)
    }
    petits = 0
    for cle, essais in vus.items():
        # on ne garde que les valeurs differentes de p0 (le retour a p0 est un reglage)
        ecarts = [v - p0[cle] for v in essais if v != p0[cle]]
        if not ecarts:
            continue
        h = max(abs(e) for e in ecarts)
        attendu = 1e-6 * max(abs(p0[cle]), 1e-3)
        if abs(p0[cle]) <= attendu:  # unilateral : jusqu'a p0 + 2h
            attendu *= 2.0
        assert h == pytest.approx(attendu, rel=1e-6), (cle, h, attendu)
        if p0[cle] > 0:
            assert min(essais) >= 0.0, (cle, p0[cle], min(essais))
        if 0 < abs(p0[cle]) < 1:
            petits += 1
    assert petits > 0  # le cas |p0| < 1 est bien couvert


def test_pas_unilateral_si_p0_moins_h_change_de_signe(smib):
    """Une resistance de ligne nulle est derivee vers le haut seulement."""
    el = next(e for e in smib.all_elements if e.name == "L1")
    el.set_param("R", 0.0)
    _, vus = valeurs_essayees(smib, 0, elements=["L1"], parametres=["R"])
    essais = vus[("L1", "R")]
    assert min(essais) >= 0.0 and max(essais) > 0.0


def test_differences_collent_a_l_analytique_sur_les_petits_parametres(smib):
    """Les differences finies (pas relatif) redonnent l'analytique, y compris
    pour les parametres de valeur < 1 ; ecart relatif <= 1e-3 sur les
    sensibilites significatives."""
    ana = {(s.element, s.parametre): s for s in sensibilites(smib, 0)}
    dif = {
        (s.element, s.parametre): s
        for s in sensibilites(smib, 0, methode="differences")
    }
    assert ana.keys() == dif.keys()
    max_norm = max(abs(s.normalisee) for s in ana.values())
    verifies = 0
    for cle, a in ana.items():
        if abs(a.normalisee) < 1e-3 * max_norm:
            continue
        d = dif[cle]
        ref = abs(a.dlambda)
        assert abs(d.dlambda - a.dlambda) <= 1e-3 * ref, (cle, d.dlambda, a.dlambda)
        verifies += 1
    assert verifies >= 5
    assert np.isfinite([s.dlambda for s in dif.values()]).all()
