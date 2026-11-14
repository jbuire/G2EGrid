"""Agregation d'unites en equivalents (g2egrid.equivalents)."""

from __future__ import annotations

import numpy as np
import pytest
from _etudes import dossier_etude

dossier_etude("equivalents")

from g2egrid.equivalents import Groupe, agreger, unites
from g2egrid.robustesse import _spectre
from g2egrid.tables import assembler

from reseaux import SITES, bus_unite, kundur_centrales


def _groupes(modele="technologie"):
    return [
        Groupe([bus_unite(s, j) for j in (1, 2, 3)], SITES[s][1], modele) for s in SITES
    ]


def test_regles_d_agregation_des_machines():
    t = kundur_centrales("SM")
    eq = agreger(t, _groupes())
    assert len(eq.T_SMs) == 4 and not eq.T_GFM1s
    us = [
        u for u in unites(t).values() if u.bus in {bus_unite(1, j) for j in (1, 2, 3)}
    ]
    b = max(int(r[0]) for r in t.Y_DER) + 1  # le premier equivalent : site 1
    p = eq.PARAMS_PAR_BUS[b]
    S = sum(u.Sn for u in us)
    assert p["H"] == pytest.approx(sum(u.params["H"] * u.Sn for u in us) / S)
    assert S / p["mp"] == pytest.approx(sum(u.Sn / u.params["mp"] for u in us))


def test_le_point_de_fonctionnement_ht_est_garde():
    t = kundur_centrales("SM+GFL")
    eq = agreger(t, _groupes())
    n1, n2 = assembler(t), assembler(eq)
    l1, l2 = n1.run_loadflow(), n2.run_loadflow()
    i1, i2 = n1._bus_index(), n2._bus_index()
    assert max(abs(l1.V[i1[b]] - l2.V[i2[b]]) for b in range(5, 12)) < 0.01


def test_un_gfm_ne_se_fond_pas_dans_une_machine():
    t = kundur_centrales("SM+GFM")
    zmin = lambda s: min(
        -100 * x.real / abs(x) for x in s.lam if x.imag > 0.3 and x.imag < 2 * np.pi * 5
    )
    z0 = zmin(_spectre(t))
    z_tech = zmin(_spectre(agreger(t, _groupes("technologie"))))
    z_mach = zmin(_spectre(agreger(t, _groupes("machine"))))
    assert abs(z_tech - z0) < 1
    assert z_mach > z0 + 10
