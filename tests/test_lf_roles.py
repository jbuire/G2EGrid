"""Les roles de commande des convertisseurs au load flow AC/DC : statisme
COURANT continu (``droopI``) et statisme Q-V cote alternatif.

Eprouves sur le CIGRE B4 : la station 111 passe de ``droop`` a ``droopI``, la
station 122 (bus PQ) recoit un statisme Q-V. Le jacobien est confronte aux
differences finies, les lois sont relues sur la solution, et le point de
fonctionnement doit etre un equilibre du modele dynamique.
"""

from __future__ import annotations

import copy

import numpy as np
import pytest

from g2egrid.catalogue import charger
from g2egrid.tables import assembler

from test_loadflow_acdc import _x_solution, jacobien_numerique

NQ = 0.05  # statisme Q-V : pu de tension par pu de reactif, sur Sn


def _b4_roles(droop_i: bool = True, qv: bool = True):
    t = copy.deepcopy(charger("cigre_b4"))
    if droop_i:
        t.T_VSCs = [
            (r[0], r[1], r[2], "droopI") if int(r[1]) == 111 else r for r in t.T_VSCs
        ]
    if qv:
        t.PARAMS_PAR_BUS.setdefault(122, {})["statisme_qv"] = NQ
    return t


@pytest.fixture(scope="module")
def b4():
    net = assembler(_b4_roles())
    net.run_loadflow()
    return net


def test_jacobien_droopI_et_qv_contre_differences_finies():
    net = assembler(_b4_roles())
    pb = net._probleme_acdc()
    assert any(c.mode == "droopI" for c in pb.vscs)
    assert any(c.kqv > 0 for c in pb.vscs)
    rng = np.random.default_rng(7)
    depart = pb.depart()
    for etiquette, x in (
        ("amorce", depart),
        ("perturbe", depart + 0.01 * rng.standard_normal(pb.n_x)),
        ("solution", _x_solution(net, pb)),
    ):
        Ja, Jn = pb.jacobien(x), jacobien_numerique(pb, x)
        err = float(np.abs(Ja - Jn).max()) / max(float(np.abs(Jn).max()), 1e-12)
        assert err < 1e-6, f"{etiquette} : ecart relatif {err:.3e}"


def test_statisme_courant_respecte(b4):
    """``p_ac = v_dc.(i_ref + s.(v_dc - v_ref))`` a la solution."""
    pb = b4._probleme_acdc()
    c = next(c for c in pb.vscs if c.mode == "droopI")
    lf = b4.lf
    vdc = float(lf.Vdc[c.idc])
    i_ref = c.p_ac / c.vdc_ref
    attendu = vdc * (i_ref + c.pente * (vdc - c.vdc_ref))
    P = lf.S.real[c.iac] - c.p_autres
    assert abs(P - attendu) < 1e-9
    # la loi n'est pas le statisme en puissance : il faut l'avoir vraiment
    # employee
    assert abs(attendu - (c.p_ac + c.pente * (vdc - c.vdc_ref))) > 1e-6


def test_statisme_qv_respecte(b4):
    """``q = q_ref + kqv.(v_ref - |V|)`` au bus PQ de la station 122."""
    pb = b4._probleme_acdc()
    c = next(c for c in pb.vscs if c.kqv > 0)
    V = abs(b4.lf.V[c.iac])
    Q = b4.lf.S.imag[c.iac] - c.q_autres
    assert abs(Q - (c.q_ac + c.kqv * (c.v_ref_ac - V))) < 1e-9
    assert abs(V - c.v_ref_ac) > 1e-5  # la tension a bien bouge


def test_point_de_fonctionnement_est_un_equilibre(b4):
    """Le modele dynamique (Kpdc du droopI cale sur la pente du load flow)
    est a l'equilibre au point du load flow, et le reseau reste stable."""
    b4.solve_equilibrium()
    ss = b4.build_state_space()
    lam = np.linalg.eigvals(ss.A)
    lam = lam[np.abs(lam) > 1e-6]
    assert float(lam.real.max()) < 0.0


def test_statisme_qv_refuse_sur_un_bus_pv():
    t = _b4_roles(droop_i=False, qv=False)
    t.PARAMS_PAR_BUS.setdefault(121, {})["statisme_qv"] = NQ  # slack (Vtheta)
    with pytest.raises(ValueError, match="statisme Q-V"):
        assembler(t).run_loadflow()
