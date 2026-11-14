"""GFM a oscillateur virtuel (dVOC) et a commande par appariement (matching),
choisis par le reglage ``loi`` : equilibre exact, lois de synchronisation
lues sur le jacobien, reseaux stables, et refus d'une loi sans son bus
continu."""

import numpy as np
import pytest

from g2egrid.catalogue import charger
from g2egrid.elements.gfm import GFMdVOC, GFMMatching
from g2egrid.tables import assembler

TOUTES = ("T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s")


def _avec_loi(nom: str, loi: str, table: str):
    t = charger(nom).copie()
    lignes = [r for k in TOUTES for r in getattr(t, k)]
    for k in TOUTES:
        getattr(t, k).clear()
    getattr(t, table).extend(lignes)
    t.GFM_PARAMS = dict(t.GFM_PARAMS, loi=loi)
    net = assembler(t)
    net.solve_equilibrium()
    return net


CAS = [("dvoc", "T_GFM1s", GFMdVOC), ("matching", "T_GFM3s", GFMMatching)]


@pytest.mark.parametrize("nom", ["smib_gfm", "wscc9_gfm"])
@pytest.mark.parametrize("loi,table,classe", CAS)
def test_equilibre_et_stabilite(nom, loi, table, classe):
    net = _avec_loi(nom, loi, table)
    gfm = [e for e in net.all_elements if e.name.startswith("GFM")]
    assert gfm and all(isinstance(e, classe) for e in gfm)
    for e in gfm:
        f = e.compiled().f(e.x0, e.z0, e.u0, e.p)
        assert np.abs(f).max() < 1e-8
    assert net.build_state_space().is_stable()


def _derivee(e, etat_eq, var):
    """d(f[etat_eq])/d(var) au point d'equilibre, var etat ou entree."""
    A, B = e.compiled().linearize(e.x0, e.z0, e.u0, e.p)[:2]
    i = e.states.index(etat_eq)
    if var in e.states:
        return A[i, e.states.index(var)]
    return B[i, e.inputs.index(var)]


def test_dvoc_lois():
    net = _avec_loi("smib_gfm", "dvoc", "T_GFM1s")
    e = next(e for e in net.all_elements if isinstance(e, GFMdVOC))
    p = e.pdict
    Vr = e.x0[e.states.index("Vr")]
    # w = w_ref + eta.(p_ref/v_ref^2 - pm/Vr^2) ; dtheta/dt = wb.(w - w_ref)
    assert _derivee(e, "theta", "pm") == pytest.approx(-p["wb"] * p["eta"] / Vr**2)
    # dVr/dt contient -wb.eta.qm/Vr
    assert _derivee(e, "Vr", "qm") == pytest.approx(-p["wb"] * p["eta"] / Vr)
    # l'amplitude est rappelee vers v_ref : -2.wb.eta.alpha a v_ref = Vr
    assert _derivee(e, "Vr", "Vr") < 0


def test_matching_frequence_suit_vdc():
    net = _avec_loi("smib_gfm", "matching", "T_GFM3s")
    e = next(e for e in net.all_elements if isinstance(e, GFMMatching))
    vdc = e.x0[e.states.index("vdc")]
    assert _derivee(e, "theta", "vdc") == pytest.approx(e.pdict["wb"] / vdc)
    assert "pm" not in e.states


def test_loi_incompatible():
    with pytest.raises(ValueError, match="matching"):
        _avec_loi("smib_gfm", "dvoc", "T_GFM3s")
    with pytest.raises(ValueError, match="dvoc"):
        _avec_loi("smib_gfm", "matching", "T_GFM1s")
