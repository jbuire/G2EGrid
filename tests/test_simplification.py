"""Hypotheses de simplification du calcul des modes (Tables.SIMPLIFICATIONS,
g2egrid.simplification)."""

from __future__ import annotations

import pytest

from g2egrid.catalogue import charger
from g2egrid.robustesse import _spectre
from g2egrid.simplification import (
    PREREGLES,
    comparer,
    inventaire,
    preregle,
)
from g2egrid.tables import assembler, tables_depuis_json, tables_vers_json


@pytest.fixture(scope="module")
def wscc():
    return charger("wscc9")


@pytest.fixture(scope="module")
def complet(wscc):
    return _spectre(wscc)


def _n_etats(t):
    net = assembler(t)
    return sum(len(el.states) for el in net.all_elements)


def test_par_type_puis_par_element(wscc):
    n = _n_etats(wscc)
    flux = {"SynchronousMachine": {"phi_d": "statique", "phi_q": "statique"}}
    assert _n_etats(wscc.copie(SIMPLIFICATIONS=flux)) == n - 6
    # l'element l'emporte sur son type
    rendu = {**flux, "SM2": {"phi_d": "dynamique", "phi_q": "dynamique"}}
    assert _n_etats(wscc.copie(SIMPLIFICATIONS=rendu)) == n - 4
    net = assembler(wscc.copie(SIMPLIFICATIONS=rendu))
    sm = {el.name: el for el in net.all_elements}
    assert type(sm["SM1"]).dynamiques()["phi_d"] == "statique"
    assert type(sm["SM2"]).dynamiques()["phi_d"] == "dynamique"


def test_nom_inconnu_signale(wscc):
    with pytest.raises(ValueError, match="SM99"):
        assembler(wscc.copie(SIMPLIFICATIONS={"SM99": {"theta": "statique"}}))


def test_aller_retour_json(wscc):
    s = {"Line": {"i_d": "statique", "i_q": "statique"}}
    t2 = tables_depuis_json(tables_vers_json(wscc.copie(SIMPLIFICATIONS=s)))
    assert t2.SIMPLIFICATIONS == s
    assert "SIMPLIFICATIONS" not in tables_vers_json(wscc)


def test_inventaire_et_preregles(wscc):
    inv = inventaire(wscc)
    assert "SynchronousMachine" in inv and "GFL1" not in inv
    assert set(inv["SynchronousMachine"]["elements"]) == {"SM1", "SM2", "SM3"}
    # un preregle ne garde que les types presents
    assert preregle(wscc, "boucles_courant") == {}
    rq = preregle(wscc, "reseau_quasi_statique")
    assert rq["Bus"] == {"vd": "statique", "vq": "statique"}
    assert set(PREREGLES) >= {"reseau_quasi_statique", "flux_stator"}
    with pytest.raises(ValueError):
        preregle(wscc, "inexistant")


def test_reseau_quasi_statique_garde_les_modes_electromecaniques(wscc, complet):
    c = comparer(wscc, preregle(wscc, "reseau_quasi_statique"), complet=complet)
    assert c.n_etats_reduit < c.n_etats
    assert c.retires and not c.perdus
    assert c.stable and c.stable_reduit
    # les modes electromecaniques bougent a peine
    assert c.ecart_max(f_max=5.0) < 0.01
    # les modes retires sont ceux du reseau, portes par les etats retires
    assert all(m.part_retiree >= 0.5 for m in c.retires)


def test_une_simplification_peut_changer_le_verdict():
    """Le retour tachymetrique de l'AVR retire seul (dx/dt = 0) destabilise
    WSCC9 : la comparaison doit le dire."""
    t = charger("wscc9")
    c = comparer(t, {"SynchronousMachine": {"e3": "statique"}})
    assert c.stable and not c.stable_reduit
    assert not c.sure()
