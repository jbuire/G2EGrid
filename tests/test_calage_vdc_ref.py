"""Calage de ``vdc_ref`` des GFM3/GFM4 avec ``Kpdc`` negatif.

Batterie en charge, P < -0,19 pu : la correction CPL de ``Kpdc`` change de
signe et l'ancienne branche (``Kpdc > 1e-12``) ne calait plus ``vdc_ref``.
"""

from __future__ import annotations

import pytest

from g2egrid import Network
from g2egrid.elements import GFM3, GFM4, IB, Line


def _banc_dc(cls, P):
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=3, Pg=P, Qg=0.0, C=1e-6)
    net.add_bus(2, 20e3, btype=1, V=1.0, C=1e-6)
    net.add(cls("EQ", 1, Sn=1e6, Sb=1e6, wb=net.wb))
    net.add(Line("L", 1, 2, R=1e-6, X=1e-5, wb=net.wb))
    net.add(IB("GRID", 2, R=1e-6, X=1e-5, wb=net.wb, Sn=1e10, Sb=1e6))
    net.run_loadflow()
    net.solve_equilibrium(tol=1e-13)
    return next(e for e in net.elements if e.name == "EQ")


@pytest.mark.parametrize("cls", [GFM3, GFM4])
@pytest.mark.parametrize("P", [0.5, -0.1, -0.3, -0.6])
def test_vdc_ref_cale_meme_avec_kpdc_negatif(cls, P):
    """``idc_ref = p_ref/vdc_ref + Kpdc.(vdc_ref - vdc)`` doit valoir ``idc``
    a l'equilibre. Avec ``Kpdc < 0`` (P < -0,19 pu : la correction CPL
    ``Kpdc = 3.Cdc/(tau_DC.wb) - Gdc + P/vdc^2`` change de signe), l'ancienne
    branche ne calait pas ``vdc_ref`` (laisse a vdc = 1) : le Newton local
    devait absorber un residu de l'ordre de 7e-10. Ici l'ecart est a
    l'arrondi."""
    el = _banc_dc(cls, P)
    p = el.pdict
    assert (p["Kpdc"] < 0) == (P < -0.19)
    x = lambda n: el.x0[el.states.index(n)]  # noqa: E731
    u = lambda n: el.u0[el.input_names().index(n)]  # noqa: E731
    idc_ref = u("p_ref") / u("vdc_ref") + p["Kpdc"] * (u("vdc_ref") - x("vdc"))
    assert idc_ref - x("idc") == pytest.approx(0.0, abs=1e-12)
    if p["Kpdc"] < 0:
        assert u("vdc_ref") != pytest.approx(x("vdc"), abs=1e-3)  # bien cale
        # racine physique : la plus proche de vdc
        assert abs(u("vdc_ref") - x("vdc")) < 0.1
