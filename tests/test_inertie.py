"""Inertie equivalente des boucles de puissance (Kelada, Buire, Hadjsaid
2025) et inertie effective mesuree sur une simulation (:mod:`g2egrid.inertie`)."""

import math

import pytest

import g2egrid.catalogue as cat
from g2egrid.inertie import (
    inertie_equivalente,
    inertie_glissante,
    inertie_reseau,
    reglage_equivalent,
    sorties_inertie,
)
from g2egrid.tables import assembler
from g2egrid.temporel import Manoeuvre, ModeleNonLineaire


@pytest.mark.parametrize(
    "loi, table, extra",
    [
        ("statisme", "T_GFM1s", {}),
        ("statisme2", "T_GFM2s", {}),
        ("dvoc", "T_GFM1s", {}),
        ("appariement", "T_GFM3s", {}),
    ],
)
def test_reglage_puis_inertie(loi, table, extra):
    """Le reglage equivalent d'une loi redonne, par les formules de
    l'article, l'inertie et l'amortissement demandes."""
    H, D = 4.0, 25.0
    t = cat.charger("wscc9_gfm")
    bus = int(t.T_GFM1s[0][1])
    ligne = t.T_GFM1s[0]
    t.T_GFM1s = [r for r in t.T_GFM1s if int(r[1]) != bus]
    setattr(t, table, [*getattr(t, table), ligne])
    wb = 2 * math.pi * t.f0  # WSCC9 est a 60 Hz
    t.PARAMS_PAR_BUS = {bus: reglage_equivalent(loi, H, D, wb=wb, **extra)}
    net = assembler(t)
    el = next(e for e in net.elements if e.name == f"GFM{bus}")
    i = inertie_equivalente(el)
    assert i.H == pytest.approx(H) and i.D == pytest.approx(D)


def test_inertie_du_filtre_de_mesure():
    """Statisme du premier ordre : 2H = 1/(m_p.w_f), D = 1/m_p."""
    net = cat.construire("wscc9_gfm")
    el = next(e for e in net.elements if type(e).__name__ == "GFM1")
    i = inertie_equivalente(el)
    p = el.pdict
    assert i.H == pytest.approx(1 / (2 * p["mp"] * p["wf_p"]))
    assert i.D == pytest.approx(1 / p["mp"])


def test_reglage_inaccessible():
    with pytest.raises(ValueError, match="inaccessible"):
        reglage_equivalent("statisme2", 0.001, 20.0, wf_p=10.0)


def test_inertie_mesuree_des_machines():
    """Sur un echelon de charge, l'inertie mesuree des machines, fenetre par
    fenetre, est leur inertie nominale (puissance d'entrefer)."""
    net = cat.construire("wscc9")
    m = ModeleNonLineaire(net)
    charge = next(e.name for e in net.elements if e.name.startswith("Load"))
    tr = m.simuler(
        (0, 1.2),
        uext_fn=m.scenario(Manoeuvre(0.2, "charge", charge, dp=0.2)),
        methode="BDF2",
        pas=2e-3,
        sorties=sorties_inertie(net),
    )
    g = inertie_glissante(tr, net, 0.2, duree=0.8)
    nominal = inertie_reseau(net)["H_sm"]
    assert g["nominal"]["machines"] == pytest.approx(nominal)
    assert g["cumul"]["machines"][-1] == pytest.approx(nominal, rel=0.03)
