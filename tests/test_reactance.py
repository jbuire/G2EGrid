"""Reactances shunt : une branche R-L vers la terre, pas une capacite negative."""

from __future__ import annotations

import numpy as np
import pytest

import g2egrid.catalogue as cat
from g2egrid.elements import ReactanceShunt
from g2egrid.elements.load import RAPPORT_XR_REACTANCE
from g2egrid.tables import assembler
from g2egrid.temporel import ModeleNonLineaire


def _wscc9_avec(B: float):
    t = cat.charger("wscc9")
    t.T_Shunts = [*t.T_Shunts, (99, 5, B)]
    net = assembler(t)
    net.run_loadflow()
    net.solve_equilibrium()
    return net


@pytest.mark.parametrize("B", [-0.1, -0.3, -1.0])
def test_point_de_fonctionnement_reste_un_equilibre(B):
    """Avant : au-dela de la charge des lignes du noeud (0,24 pu au bus 5),
    la reactance faisait basculer le noeud en resistance de causalite et
    disparaissait du modele -- |dx/dt| valait 2,5e4 a B = -0,3."""
    net = _wscc9_avec(B)
    dx, dzu = ModeleNonLineaire(net).residu_equilibre()
    assert dx < 1e-6 and dzu < 1e-9
    assert type(net.buses[5]).__name__ == "Bus"  # le noeud garde sa capacite
    assert net.build_state_space().is_stable()


def test_la_reactance_consomme_son_admittance():
    """Q consomme = |B|.V^2 (a la resistance pres), lu sur le modele dynamique."""
    B = -0.3
    net = _wscc9_avec(B)
    (x,) = [e for e in net.elements if isinstance(e, ReactanceShunt)]
    V = abs(net.lf.V[net._bus_index()[5]])
    m = ModeleNonLineaire(net)
    tr = m.simuler((0.0, 1e-3), sorties=[f"{x.name}.p", f"{x.name}.q"])
    q = tr.sorties[f"{x.name}.q"][-1]
    p = tr.sorties[f"{x.name}.p"][-1]
    XR = RAPPORT_XR_REACTANCE
    attendu = abs(B) * V**2 / (1 + 1 / XR**2)  # Im(conj(y)).|V|^2
    assert q == pytest.approx(attendu, rel=1e-6)
    assert p == pytest.approx(attendu / XR, rel=1e-6)


def test_un_condensateur_ne_change_rien():
    """B > 0 : de la capacite de noeud, exactement comme avant."""
    t = cat.charger("wscc9")
    n0 = len(assembler(t).all_elements)
    t.T_Shunts = [*t.T_Shunts, (99, 5, 0.2)]
    net = assembler(t)
    assert len(net.all_elements) == n0
    assert net.bus_specs[5].C == pytest.approx(0.088 + 0.153 + 0.2)


def test_la_reactance_est_une_inductance():
    """Son admittance vue du bus decroit avec la frequence : c'est ce qu'une
    capacite negative ne fait pas (elle croit comme w)."""
    x = ReactanceShunt("X", 1, B=-0.5, wb=2 * np.pi * 50)
    R, X = x.pdict["R"], x.pdict["X"]
    y = [abs(1 / complex(R, X * w)) for w in (0.5, 1.0, 2.0)]
    assert y[0] > y[1] > y[2]
    with pytest.raises(ValueError, match="NEGATIVE"):
        ReactanceShunt("C", 1, B=0.5, wb=100.0)


def test_capacite_negative_refusee():
    t = cat.charger("wscc9")
    t.T_Lines_pu = [
        (*r[:6], -1.0, *r[7:]) if int(r[0]) == 1 else r for r in t.T_Lines_pu
    ]
    with pytest.raises(ValueError, match="capacite shunt negative"):
        assembler(t)


def test_simulation_avec_reactance():
    """Un echelon de consigne : la trajectoire reste finie et revient."""
    net = _wscc9_avec(-0.3)
    m = ModeleNonLineaire(net)
    sm = next(n for n in m.entrees if n.endswith("p_ref"))
    tr = m.simuler((0.0, 0.5), uext_fn=m.echelon(sm, 0.01, t0=0.05))
    assert np.isfinite(tr.x).all()
