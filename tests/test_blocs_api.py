"""FR-051, routes de l'API des blocs (transport reel, ``fastapi.testclient``) :
regrouper, degrouper, aplatir, depuis_element (classe, ou element d'un reseau
avec SES parametres), bibliotheque de modeles persistante. Chaque route a un
cas nominal et un cas d'erreur : 422 avec un message en francais, jamais de
trace ; un modele de la bibliotheque qui n'existe pas est une 404."""

from __future__ import annotations

import copy

import pytest

pytest.importorskip("fastapi")
from fastapi.testclient import TestClient

import g2egrid.reseau.catalogue as cat
from g2egrid.elements.blocs import compiler, exemples
from g2egrid.reseau.catalogue import charger
from g2egrid.reseau.tables import tables_vers_json
from g2egrid_api.app import app


@pytest.fixture(scope="module")
def wscc9_cigre6():
    """Lu AVANT que ``client`` ne detourne le catalogue."""
    return charger("wscc9_cigre6")


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setattr(cat, "DOSSIER_RESEAUX", tmp_path)
    c = TestClient(app, headers={"Host": "127.0.0.1"})
    r = c.get("/api/blocs/types")
    if r.status_code == 403:
        pytest.skip("atelier protege par jeton dans cet environnement")
    return c


def _mini() -> dict:
    return {
        "nom": "mini",
        "parametres": {"k": 2.0, "tau": 0.1},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref"},
            {"id": "qr", "type": "consigne", "nom": "Q_ref"},
            {"id": "g", "type": "gain", "K": "k"},
            {"id": "f", "type": "premier_ordre", "K": 1, "tau": "tau"},
            {"id": "op", "type": "sortie", "signal": "id"},
            {"id": "oq", "type": "sortie", "signal": "iq"},
        ],
        "liens": [
            {"de": "pr", "vers": "g", "entree": 0},
            {"de": "g", "vers": "f", "entree": 0},
            {"de": "f", "vers": "op", "entree": 0},
            {"de": "qr", "vers": "oq", "entree": 0},
        ],
    }


def test_types_categories_et_nouveaux_blocs(client):
    d = client.get("/api/blocs/types").json()
    assert d["types"]["pid"]["categorie"] == "dynamique"
    assert d["types"]["port_entree"]["categorie"] == "structure"
    assert d["types"]["sous_systeme"]["champs"] == ["nom", "schema"]
    assert {c["id"] for c in d["categories"]} >= {"sources", "operations", "structure"}
    assert all(t["categorie"] for t in d["types"].values())
    for k in (
        "pid",
        "avance_retard",
        "second_ordre",
        "fonction_transfert",
        "saturation",
        "implicite",
    ):
        assert k in d["types"]


def test_compiler_un_schema_a_sous_systemes(client):
    r = client.post(
        "/api/blocs/regrouper",
        json={"schema_bloc": _mini(), "ids": ["g", "f"], "nom": "chaine"},
    )
    assert r.status_code == 200, r.text
    s = r.json()["schema"]
    d = client.post("/api/blocs/compiler", json={"schema_bloc": s})
    assert d.status_code == 200, d.text
    assert d.json()["consignes"] == ["P_ref", "Q_ref"]
    assert any("chaine__f" in e["texte"] for e in d.json()["equations"])


def test_regrouper_degrouper_aplatir(client):
    s = _mini()
    r = client.post(
        "/api/blocs/regrouper",
        json={"schema_bloc": s, "ids": ["g", "f"], "nom": "chaine"},
    )
    sch = r.json()["schema"]
    nouveau = ({b["id"] for b in sch["blocs"]} - {b["id"] for b in s["blocs"]}).pop()
    plat = client.post("/api/blocs/aplatir", json={"schema_bloc": sch}).json()["schema"]
    assert {b["id"] for b in plat["blocs"]} >= {f"{nouveau}__g", f"{nouveau}__f"}
    back = client.post(
        "/api/blocs/degrouper", json={"schema_bloc": sch, "id": nouveau}
    ).json()["schema"]
    assert {b["id"] for b in back["blocs"]} == {b["id"] for b in s["blocs"]}
    assert compiler(back).etats == compiler(s).etats


def test_regrouper_a_l_interieur_par_le_niveau_affiche(client):
    """L'interface envoie le niveau affiche ({nom, parametres, blocs, liens},
    ports compris) et remplace blocs/liens par ceux du resultat."""
    sch = client.post(
        "/api/blocs/regrouper",
        json={"schema_bloc": _mini(), "ids": ["g", "f"], "nom": "chaine"},
    ).json()["schema"]
    ss = next(b for b in sch["blocs"] if b["type"] == "sous_systeme")
    niveau = {
        "nom": "mini",
        "parametres": sch["parametres"],
        **copy.deepcopy(ss["schema"]),
    }
    avant = {b["id"] for b in niveau["blocs"]}
    r = client.post(
        "/api/blocs/regrouper",
        json={"schema_bloc": niveau, "ids": ["g"], "nom": "gain"},
    )
    assert r.status_code == 200, r.text
    out = r.json()["schema"]
    assert avant - {b["id"] for b in out["blocs"]} == {"g"}
    assert len({b["id"] for b in out["blocs"]} - avant) == 1
    # un port ne se regroupe pas
    port = next(b["id"] for b in niveau["blocs"] if b["type"] == "port_entree")
    r = client.post(
        "/api/blocs/regrouper", json={"schema_bloc": niveau, "ids": [port], "nom": "p"}
    )
    assert r.status_code == 422 and "port" in r.json()["detail"]
    # chemin explicite depuis le premier niveau
    r = client.post(
        "/api/blocs/regrouper",
        json={"schema_bloc": sch, "ids": ["g"], "nom": "gain", "chemin": [ss["id"]]},
    )
    assert r.status_code == 200, r.text
    assert compiler(r.json()["schema"]).etats == ["x_chaine__f"]


@pytest.mark.parametrize(
    "route, corps, mot",
    [
        ("regrouper", {"ids": [], "nom": "x"}, "aucun bloc"),
        ("regrouper", {"ids": ["zz"], "nom": "x"}, "inconnu"),
        ("regrouper", {"ids": ["g"], "nom": ""}, "nom"),
        ("degrouper", {"id": "g"}, "pas un sous-systeme"),
        ("degrouper", {"id": "absent"}, "pas un sous-systeme"),
    ],
)
def test_erreurs_422_en_francais(client, route, corps, mot):
    r = client.post(f"/api/blocs/{route}", json={"schema_bloc": _mini(), **corps})
    assert r.status_code == 422, r.text
    assert mot in r.json()["detail"] and "Traceback" not in r.text


def test_aplatir_erreur(client):
    s = _mini()
    s["blocs"].append({"id": "p", "type": "port_entree", "index": 0, "nom": "a"})
    r = client.post("/api/blocs/aplatir", json={"schema_bloc": s})
    assert r.status_code == 422 and "sous-systeme" in r.json()["detail"]


def test_depuis_element_classe(client):
    r = client.post(
        "/api/blocs/depuis_element",
        json={"classe": "GFL1", "parametres": {"Kipll": 50.0}},
    )
    assert r.status_code == 200, r.text
    s = r.json()["schema"]
    assert s["parametres"]["Kipll"] == 50.0 and s["amorce"]["classe"] == "GFL1"
    assert compiler(s).consignes == ["p_ref", "q_ref"]
    r = client.post(
        "/api/blocs/depuis_element", json={"classe": "IB", "detail": "equations"}
    )
    assert r.status_code == 200
    assert not any(b["type"] == "sous_systeme" for b in r.json()["schema"]["blocs"])


def test_depuis_element_d_un_reseau_avec_ses_parametres(client, wscc9_cigre6):
    t = copy.deepcopy(wscc9_cigre6)
    t.GFL_PARAMS = {"Kipll": 77.0}
    r = client.post(
        "/api/blocs/depuis_element",
        json={"tables": tables_vers_json(t), "table": "T_GFL1s", "index": 2},
    )
    assert r.status_code == 200, r.text
    s = r.json()["schema"]
    assert s["amorce"]["classe"] == "GFL1" and s["parametres"]["Kipll"] == 77.0
    # le 60 Hz du reseau est celui du schema : wb = 2 pi 60 dans les gains regles
    d0 = client.post("/api/blocs/depuis_element", json={"classe": "GFL1"}).json()[
        "schema"
    ]
    assert d0["parametres"]["KpCL"] != s["parametres"]["KpCL"]


def test_depuis_element_erreurs(client, wscc9_cigre6):
    for corps, mot in [
        ({"classe": "GFM3DC"}, "port continu"),
        ({"classe": "Line"}, "port"),
        ({"classe": "Nexiste"}, "inconnue"),
        ({"classe": "GFL1", "detail": "xx"}, "detail"),
        ({}, "classe"),
    ]:
        r = client.post("/api/blocs/depuis_element", json=corps)
        assert r.status_code == 422, (corps, r.text)
        assert mot in r.json()["detail"] and "Traceback" not in r.text
    t = tables_vers_json(wscc9_cigre6)
    r = client.post(
        "/api/blocs/depuis_element", json={"tables": t, "table": "T_SMs", "index": 9}
    )
    assert r.status_code == 422 and "indice" in r.json()["detail"]
    r = client.post(
        "/api/blocs/depuis_element", json={"tables": t, "table": "T_Nul", "index": 0}
    )
    assert r.status_code == 422 and "table inconnue" in r.json()["detail"]


def test_bibliotheque_aller_retour(client, tmp_path):
    assert client.get("/api/blocs/bibliotheque").json() == []
    s = exemples()["gfl_complet"]["schema"]
    s["titre"] = "Mon GFL"
    r = client.post("/api/blocs/bibliotheque", json={"schema_bloc": s})
    assert r.status_code == 200, r.text
    assert r.json() == {"nom": "gfl_complet"}
    assert (tmp_path / "sauvegardes" / "modeles" / "gfl_complet.json").is_file()
    liste = client.get("/api/blocs/bibliotheque").json()
    assert [(m["nom"], m["titre"]) for m in liste] == [("gfl_complet", "Mon GFL")]
    assert liste[0]["modifie"]
    assert client.get("/api/blocs/bibliotheque/gfl_complet").json() == {"schema": s}
    # remplacer
    s2 = copy.deepcopy(s)
    s2["parametres"]["Kipll"] = 5.0
    client.post("/api/blocs/bibliotheque", json={"schema_bloc": s2})
    assert (
        client.get("/api/blocs/bibliotheque/gfl_complet").json()["schema"][
            "parametres"
        ]["Kipll"]
        == 5.0
    )
    assert len(client.get("/api/blocs/bibliotheque").json()) == 1
    # le modele range ne gene pas la liste des reseaux enregistres
    assert cat.sauvegardes() == ()
    # supprimer
    r = client.delete("/api/blocs/bibliotheque/gfl_complet")
    assert r.status_code == 200 and r.json() == {"nom": "gfl_complet"}
    assert client.get("/api/blocs/bibliotheque").json() == []


def test_bibliotheque_schema_avec_sous_systemes(client):
    s = client.post(
        "/api/blocs/regrouper",
        json={"schema_bloc": _mini(), "ids": ["g", "f"], "nom": "chaine"},
    ).json()["schema"]
    assert client.post("/api/blocs/bibliotheque", json={"schema_bloc": s}).json() == {
        "nom": "mini"
    }
    assert client.get("/api/blocs/bibliotheque/mini").json()["schema"] == s


def test_bibliotheque_erreurs(client):
    for nom in ("../x", "a b", "", "x" * 80, "é"):
        s = _mini()
        s["nom"] = nom
        r = client.post("/api/blocs/bibliotheque", json={"schema_bloc": s})
        assert r.status_code == 422, (nom, r.text)
        assert "nom de modele invalide" in r.json()["detail"]
    r = client.post("/api/blocs/bibliotheque", json={"schema_bloc": {"nom": "x"}})
    assert r.status_code == 422 and "schema-bloc" in r.json()["detail"]
    r = client.get("/api/blocs/bibliotheque/inconnu")
    assert r.status_code == 404 and "absent" in r.json()["detail"]
    r = client.delete("/api/blocs/bibliotheque/inconnu")
    assert r.status_code == 404
    r = client.get("/api/blocs/bibliotheque/a b")
    assert r.status_code == 422
    r = client.delete("/api/blocs/bibliotheque/a%20b")
    assert r.status_code == 422
    # rien n'est sorti du dossier de la bibliotheque
    assert client.get("/api/blocs/bibliotheque").json() == []
