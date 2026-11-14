"""FR-052, routes de l'API des equations modifiables (transport reel,
``fastapi.testclient``) : reseau, element, modifier, ajouter, supprimer,
revenir, verifier, depuis_schema, vers_blocs, et ``detail="liste"`` de
``/api/blocs/depuis_element``. Chaque route a un cas nominal et un cas
d'erreur : 422 avec un message en francais, jamais de trace."""

from __future__ import annotations

import pytest

pytest.importorskip("fastapi")
from fastapi.testclient import TestClient

import g2egrid.reseau.catalogue as cat
from g2egrid.elements import blocs as B
from g2egrid.reseau.catalogue import charger
from g2egrid.reseau.tables import tables_vers_json
from g2egrid_api.app import app


@pytest.fixture(scope="module")
def smib_gfl():
    t = charger("smib_gfl")
    t.T_GFL1s, t.T_GFL2s = t.T_GFL2s, []
    return tables_vers_json(t)


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setattr(cat, "DOSSIER_RESEAUX", tmp_path)
    c = TestClient(app, headers={"Host": "127.0.0.1"})
    if c.get("/api/blocs/types").status_code == 403:
        pytest.skip("atelier protege par jeton dans cet environnement")
    return c


def test_reseau(client, smib_gfl):
    r = client.post("/api/equations/reseau", json={"tables": smib_gfl})
    assert r.status_code == 200
    lignes = r.json()["lignes"]
    champs = {"element", "table", "index", "bus", "classe", "genre", "variable"}
    champs |= {"expression", "modifiable", "modifiee"}
    assert champs <= set(lignes[0])
    assert any(x["element"].startswith("GFL") and x["modifiable"] for x in lignes)
    assert all(x["pourquoi"] for x in lignes if not x["modifiable"])


def test_element_et_erreur(client, smib_gfl):
    r = client.post(
        "/api/equations/element",
        json={"tables": smib_gfl, "table": "T_GFL1s", "index": 0},
    )
    assert r.status_code == 200
    m = r.json()["modele"]
    assert m["genre"] == "equations" and m["etats"][0]["derivee"]
    r = client.post(
        "/api/equations/element",
        json={"tables": smib_gfl, "element": "bus1"},
    )
    assert r.status_code == 422 and "ne se modifie pas" in r.json()["detail"]


def _modifier(client, tables, expr="2*Kipll*veq"):
    return client.post(
        "/api/equations/modifier",
        json={
            "tables": tables,
            "table": "T_GFL1s",
            "index": 0,
            "variable": "M_pll",
            "expression": expr,
        },
    )


def test_modifier_puis_revenir(client, smib_gfl):
    r = _modifier(client, smib_gfl)
    assert r.status_code == 200
    corps = r.json()
    assert corps["tables"]["T_Blocs"] and not corps["tables"].get("T_GFL1s")
    assert any(
        e["nom"] == "M_pll" and e["derivee"] == "2*Kipll*veq"
        for e in corps["modele"]["etats"]
    )
    r = client.post(
        "/api/equations/revenir",
        json={"tables": corps["tables"], "table": "T_Blocs", "index": 0},
    )
    assert r.status_code == 200
    assert r.json()["tables"] == smib_gfl
    # un natif n'a rien a retrouver
    r = client.post(
        "/api/equations/revenir",
        json={"tables": smib_gfl, "table": "T_GFL1s", "index": 0},
    )
    assert r.status_code == 422 and "natif" in r.json()["detail"]


def test_modifier_refus_422(client, smib_gfl):
    r = _modifier(client, smib_gfl, "__import__('os').system('id')")
    assert r.status_code == 422
    assert "detail" in r.json() and "Traceback" not in r.json()["detail"]


def test_ajouter_et_supprimer(client, smib_gfl):
    r = client.post(
        "/api/equations/ajouter",
        json={
            "tables": smib_gfl,
            "table": "T_GFL1s",
            "index": 0,
            "genre": "etat",
            "variable": "xf",
            "expression": "(M_pll - xf)/0.1",
            "depart": 0.0,
        },
    )
    assert r.status_code == 200
    corps = r.json()
    assert any(e["nom"] == "xf" for e in corps["modele"]["etats"])
    r = client.post(
        "/api/equations/supprimer",
        json={
            "tables": corps["tables"],
            "table": "T_Blocs",
            "index": 0,
            "variable": "xf",
        },
    )
    assert r.status_code == 200
    assert not any(e["nom"] == "xf" for e in r.json()["modele"]["etats"])
    r = client.post(
        "/api/equations/supprimer",
        json={
            "tables": corps["tables"],
            "table": "T_Blocs",
            "index": 0,
            "variable": "M_pll",
        },
    )
    assert r.status_code == 422 and "encore employe" in r.json()["detail"]
    r = client.post(
        "/api/equations/ajouter",
        json={
            "tables": smib_gfl,
            "table": "T_GFL1s",
            "index": 0,
            "genre": "etat",
            "variable": "M_pll",
            "expression": "1",
        },
    )
    assert r.status_code == 422 and "deja pris" in r.json()["detail"]


def test_verifier(client):
    m = B.equations_depuis_element("GFL1")
    r = client.post("/api/equations/verifier", json={"modele": m})
    assert r.status_code == 200 and r.json()["ok"] and r.json()["equations"]
    m["etats"][0]["derivee"] = "inconnu*2"
    r = client.post("/api/equations/verifier", json={"modele": m})
    assert r.status_code == 200 and not r.json()["ok"] and r.json()["message"]


def test_depuis_schema_et_vers_blocs(client):
    r = client.post(
        "/api/equations/depuis_schema", json={"schema_bloc": B.schema_gfl()}
    )
    assert r.status_code == 200
    m = r.json()["modele"]
    assert m["genre"] == "equations"
    r = client.post("/api/equations/vers_blocs", json={"modele": m})
    assert r.status_code == 200 and r.json()["schema"]["blocs"]
    r = client.post("/api/equations/depuis_schema", json={"schema_bloc": {"blocs": []}})
    assert r.status_code == 422


def test_blocs_depuis_element_liste_et_compiler(client):
    r = client.post(
        "/api/blocs/depuis_element", json={"classe": "GFL1", "detail": "liste"}
    )
    assert r.status_code == 200
    m = r.json()["modele"]
    assert m["genre"] == "equations" and r.json()["schema"] == m
    r = client.post("/api/blocs/compiler", json={"schema_bloc": m})
    assert r.status_code == 200 and len(r.json()["etats"]) == 12
