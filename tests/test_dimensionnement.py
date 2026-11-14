"""Les conditions de dimensionnement IB/GFM et IB/GFL
(``etudes/dimensionnement``) : physique et petit signal seulement, le grand
signal est trop long pour la suite de tests."""

import math

import pytest
from _etudes import dossier_etude

ICI = dossier_etude("dimensionnement")

from cas import Point
from contraintes import (
    Specs,
    bande,
    f_resonance,
    petit_signal,
    physiques,
    regime,
)


def _echecs(vs):
    return {v.code for v in vs if not v.ok}


def test_bande_passante_de_l_autoreglage():
    # wn = 3/(zeta.tau) : tau = 0,1 ms, zeta = 0,707 -> 6,75 kHz
    assert bande(0.1e-3, 0.707) == pytest.approx(6753, abs=1)


def test_resonance_lcl():
    p = dict(L1_pu=0.1, L2_pu=0.01, C_pu=0.05)
    assert abs(f_resonance(p) - 50 * math.sqrt(0.11 / (0.1 * 0.01 * 0.05))) < 1e-9


def test_defauts_g2egrid_hors_plages_physiques():
    """Le GFM par defaut a une boucle de courant a 6,8 kHz : irrealiste a
    5 kHz de decoupage ; et un condensateur de 10 %."""
    e = _echecs(physiques(Point("GFM"), Specs()))
    assert {"f_CL_decoupage", "C_filtre"} <= e


def test_reference_gfm_tient_physique_et_petit_signal():
    p = dict(
        C_pu=0.05,
        tau_CL=1.5e-3,
        tau_VL=0.15,
        mp=0.02,
        nq=0.02,
        wf_p=2 * math.pi * 0.8,
        wf_q=2 * math.pi * 0.8,
    )
    pt = Point("GFM", params=p)
    s = Specs()
    assert not _echecs(physiques(pt, s) + regime(pt, s))
    assert not _echecs(petit_signal(pt, s))


def test_separation_usuelle_insuffisante():
    """tau_VL/tau_CL = 10 respecte la regle k = 5 mais reste instable."""
    p = dict(
        C_pu=0.05,
        tau_CL=1.5e-3,
        tau_VL=15e-3,
        mp=0.02,
        nq=0.02,
        wf_p=2 * math.pi * 0.8,
        wf_q=2 * math.pi * 0.8,
    )
    pt = Point("GFM", params=p)
    assert "sep_VL_CL" not in _echecs(physiques(pt, Specs()))
    assert "stabilite" in _echecs(petit_signal(pt, Specs()))
