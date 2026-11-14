"""Comparer deux versions d'un reseau : tables, point de fonctionnement, modes."""

from __future__ import annotations

import g2egrid.catalogue as cat
from g2egrid.comparaison import comparer, differences_tables
from g2egrid_api import analyse
from g2egrid_api.atelier import Atelier


def _plus_de_charge(t, bus, facteur):
    b = t.copie()
    b.Y_DER = [
        tuple(r) if r[0] != bus else (*r[:7], r[7] * facteur, r[8]) for r in b.Y_DER
    ]
    return b


def test_un_reseau_face_a_lui_meme():
    a = cat.charger("wscc9")
    c = comparer(a, a.copie())
    assert c.tables == [] and c.disparus == [] and c.apparus == []
    assert c.ecart_tension == 0.0
    assert max(abs(p.lam_a - p.lam_b) for p in c.paires) == 0.0


def test_une_charge_retouchee():
    """Une seule cellule change : elle seule est rapportee, et les modes
    electromecaniques glissent sans disparaitre."""
    a = cat.charger("wscc9")
    b = _plus_de_charge(a, 5, 1.2)
    (d,) = differences_tables(a, b)
    assert d.table == "Y_DER" and d.modifiees == [(5, "Pl_W", 125e6, 150e6)]
    c = comparer(a, b)
    assert c.ecart_tension > 1e-3  # le load flow a bouge
    em = [p for p in c.paires if 0.5 < p.f_a < 3 and p.zeta_a < 30]
    assert len(em) >= 2
    assert all(abs(p.ecart_f) < 5 for p in em)
    assert any(abs(p.ecart_f) > 0.1 for p in em)


def test_lignes_ajoutees_et_reglages():
    a = cat.charger("wscc9")
    b = a.copie()
    b.T_Shunts = [*b.T_Shunts, (99, 5, 0.1)]
    b.SM_PARAMS = {**b.SM_PARAMS, "H": 9.0}
    diffs = {d.table: d for d in differences_tables(a, b)}
    assert diffs["T_Shunts"].ajoutees == [(99, 5, 0.1)]
    assert diffs["SM_PARAMS"].modifiees[0][1] == "H"


def test_api_et_atelier():
    d = analyse.comparer({"catalogue": "wscc9"}, {"catalogue": "wscc9"}, f_max=5)
    assert d["tables"] == [] and d["stable"] == [True, True]
    assert all(max(m["f_a"], m["f_b"]) <= 5 for m in d["paires"])
    a, b = Atelier.catalogue("wscc9"), Atelier.catalogue("wscc9")
    b.charge(5, 20e6, 0.0)
    df = b.comparer(a, f_max=5)
    assert {"f_a", "f_b", "ecart_zeta"} <= set(df.columns)
    assert df.attrs["tables"][0]["table"] == "Y_DER"
