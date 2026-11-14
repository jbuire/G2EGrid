"""Les reparations proposees par la validation de l'atelier : chacune fait
disparaitre le probleme qu'elle accompagne."""

from __future__ import annotations

import pytest

pytest.importorskip("fastapi")

from g2egrid.catalogue import charger  # noqa: E402
from g2egrid_api import analyse  # noqa: E402


def _probleme(t, mot, positions=None):
    pbs = analyse.valider(t, positions)
    return next((p for p in pbs if mot in p["message"]), None)


def _appliquer(t, p, positions=None):
    r = p["reparation"]
    return analyse.reparer(t, r["action"], r["args"], positions)


def test_sans_slack_prend_le_plus_gros_formeur():
    t = charger("wscc9").copie()
    t.Y_DER = [(r[0], 2 if r[1] == 1 else r[1], *r[2:]) for r in t.Y_DER]
    p = _probleme(t, "aucun bus slack")
    assert p["reparation"]["action"] == "slack"
    t2 = _appliquer(t, p)
    assert _probleme(t2, "slack") is None
    sn = {int(r[1]): r[2] for r in t.T_SMs}
    slack = next(int(r[0]) for r in t2.Y_DER if r[1] == 1)
    assert sn[slack] == max(sn.values())


def test_deux_slacks():
    t = charger("wscc9").copie()
    t.Y_DER = [(r[0], 1 if int(r[0]) in (1, 2) else r[1], *r[2:]) for r in t.Y_DER]
    t2 = _appliquer(t, _probleme(t, "bus slack dans"))
    assert sum(1 for r in t2.Y_DER if r[1] == 1) == 1
    assert analyse.assembler(t2).run_loadflow().converged


def test_bus_isole_raccorde_au_plus_proche():
    t = charger("wscc9").copie()
    t.Y_DER = [*t.Y_DER, (10, 3, t.Y_DER[4][2], 1.0, 0.0, 0.0, 0.0, 1e6, 5e5)]
    pos = {str(r[0]): (0.1 * int(r[0]), 0.5) for r in t.Y_DER}
    pos["10"] = (0.81, 0.5)  # a cote du bus 8
    p = _probleme(t, "bus 10 isole", pos)
    assert p["reparation"]["args"]["vers"] == 8
    t2 = _appliquer(t, p, pos)
    assert _probleme(t2, "isole") is None
    assert analyse.assembler(t2).run_loadflow().converged


def test_gfl_en_pv_et_pv_sans_equipement():
    t = charger("wscc9").copie()
    t.Y_DER = [(r[0], 2 if int(r[0]) == 5 else r[1], *r[2:]) for r in t.Y_DER]
    t2 = _appliquer(t, _probleme(t, "PV mais sans equipement"))
    assert next(r[1] for r in t2.Y_DER if int(r[0]) == 5) == 3


def test_charge_a_reactif_negatif():
    t = charger("wscc9").copie()
    t.Y_DER = [(*r[:8], -1e6) if int(r[0]) == 5 else r for r in t.Y_DER]
    t2 = _appliquer(t, _probleme(t, "Ql ="))
    assert _probleme(t2, "Ql =") is None
    assert any(int(r[1]) == 5 and r[4] == "P" for r in t2.T_Loads)


def test_noeud_continu_sans_capacite():
    t = charger("hvdc2t").copie()
    t.T_ShuntsDC = []
    t.T_LinesDC = [(*r[:6], 0.0) for r in t.T_LinesDC]
    # les convertisseurs portent leur condensateur : il faut un noeud sans
    b = int(t.T_BusDC[0][0])
    t.T_BusDC = [*t.T_BusDC, (99, 2, t.T_BusDC[0][2], 1.0, 0.0)]
    t.T_LinesDC = [*t.T_LinesDC, (99, b, 99, 1.0, 0.01, 0.1, 0.0)]
    p = _probleme(t, "sans capacite")
    assert p and p["reparation"]["action"] == "cdc"
    t2 = _appliquer(t, p)
    assert _probleme(t2, "sans capacite") is None


def test_reparation_inconnue():
    with pytest.raises(ValueError):
        analyse.reparer(charger("smib"), "inventee", {})


def test_projection_du_mode_sur_le_schema():
    """Le mode inter-machines de wscc9 a 2,2 Hz : SM3 et SM2 en opposition,
    et les participations des bus et des branches font tout le mode."""
    import numpy as np

    t = charger("wscc9")
    lam = analyse.construire(t).etat_stationnaire().eig()
    k = int(
        np.argmin(np.where((lam.imag > 1) & (lam.imag < 20), -lam.real / abs(lam), 9))
    )
    d = analyse.detail_mode(t, k)["projection"]
    total = sum(x["part"] for x in d["bus"].values()) + sum(
        x["part"] for x in d["branches"]
    )
    assert abs(total - 1) < 1e-9
    assert d["bus"]["3"]["etat"] == "SM3.Dwr"
    ecart = abs((d["bus"]["3"]["phase"] - d["bus"]["2"]["phase"] + 180) % 360 - 180)
    assert ecart > 150
