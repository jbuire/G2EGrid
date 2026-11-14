"""FR-051, critere 7 : l'API de l'editeur de blocs et la bibliotheque de modeles.

Tests ecrits A L'AVEUGLE depuis le seul cahier
``.claude/equipe/cahiers/FR-051_editeur_blocs.md`` (tableau des routes). Le
transport est exerce pour de bon, par ``fastapi.testclient.TestClient`` sur
``g2egrid_api.app`` (comme ``tests/test_securite_atelier.py``) : serialisation
des corps comprise.

Pour chaque route, au moins un cas nominal et un cas d'erreur ; les erreurs
sont des 422 avec un message en texte (``detail`` chaine non vide), jamais une
trace (pas de ``Traceback``, pas de 500).

- ``GET /api/blocs/types`` : ``categorie`` parmi les sept du cahier, entrees,
  sorties, champs, exemples ;
- ``POST /api/blocs/compiler`` : un schema a sous-systemes compile, et donne la
  meme reponse que le schema aplati ;
- ``POST /api/blocs/depuis_element`` : par classe, ou par element d'un reseau
  (AVEC SES parametres : un ``tau_PLL`` change dans ``GFL_PARAMS`` change le
  ``Kipll`` du schema) ; refus 422 d'un element a port continu (« continu »
  dans le message) ;
- ``regrouper`` / ``degrouper`` / ``aplatir`` ;
- bibliotheque : enregistrer puis rouvrir rend le MEME schema (egalite
  exacte du JSON), un fichier JSON par modele dans
  ``<dossier_sauvegardes>/modeles/``, lister, supprimer, noms invalides refuses
  (422), rien n'est ecrit hors du dossier ;
- la garde du serveur local (jeton) couvre les nouvelles routes.

Lectures du cahier retenues :

- ``POST /api/blocs/bibliotheque`` nomme le modele d'apres le champ ``nom`` du
  schema ; la reponse est ``{"nom": <ce nom>}`` ;
- un modele absent (lecture ou suppression) est refuse en 404 ou 422 (le
  cahier ne tranche pas) ;
- un nom de plus de 64 caracteres, ou contenant un espace, est invalide
  (regle de ``catalogue._verifier_nom`` : ``[A-Za-z0-9][A-Za-z0-9_-]{0,63}``) ;
- les reponses de ``regrouper``, ``degrouper`` et ``aplatir`` sont
  ``{"schema": ...}``.
"""

from __future__ import annotations

import copy

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")

import g2egrid.catalogue as cat
from fastapi.testclient import TestClient
from g2egrid.tables import assembler, tables_vers_json
from test_fr051_aide import schema_sans_position

from g2egrid.elements import blocs as B
from g2egrid_api.app import COOKIE_JETON, app

CATEGORIES = {
    "sources",
    "operations",
    "dynamique",
    "repere_dq",
    "limiteurs",
    "structure",
    "sorties",
}


@pytest.fixture
def client(monkeypatch, tmp_path):
    """Sans jeton (la garde est testee a part) ; les sauvegardes, donc la
    bibliotheque, vont dans un dossier jetable."""
    monkeypatch.delenv("G2EGRID_JETON", raising=False)
    monkeypatch.setattr(cat, "DOSSIER_RESEAUX", tmp_path)
    return TestClient(app)


def _ss_schema() -> dict:
    """Un schema avec un sous-systeme (id = 3 x (5 - 2) = 9, iq = 3)."""
    interieur = {
        "blocs": [
            {
                "id": "pb",
                "type": "port_entree",
                "index": 1,
                "nom": "b",
                "x": 20,
                "y": 120,
            },
            {
                "id": "pa",
                "type": "port_entree",
                "index": 0,
                "nom": "a",
                "x": 20,
                "y": 20,
            },
            {"id": "s", "type": "somme", "signes": "+-", "x": 170, "y": 70},
            {"id": "g", "type": "gain", "K": 3, "x": 320, "y": 20},
            {
                "id": "p1",
                "type": "port_sortie",
                "index": 1,
                "nom": "d",
                "x": 470,
                "y": 120,
            },
            {
                "id": "p0",
                "type": "port_sortie",
                "index": 0,
                "nom": "r",
                "x": 470,
                "y": 20,
            },
        ],
        "liens": [
            {"de": "pa", "vers": "s", "entree": 0},
            {"de": "pb", "vers": "s", "entree": 1},
            {"de": "s", "vers": "g"},
            {"de": "g", "vers": "p0"},
            {"de": "s", "vers": "p1"},
        ],
    }
    return {
        "nom": "modele_api",
        "parametres": {},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref", "x": 20, "y": 20},
            {"id": "qr", "type": "consigne", "nom": "Q_ref", "x": 20, "y": 120},
            {"id": "a", "type": "constante", "valeur": 5.0, "x": 20, "y": 220},
            {"id": "b", "type": "constante", "valeur": 2.0, "x": 20, "y": 320},
            {
                "id": "ss",
                "type": "sous_systeme",
                "nom": "diff3",
                "x": 250,
                "y": 250,
                "schema": interieur,
            },
            {"id": "oi", "type": "sortie", "signal": "id", "x": 500, "y": 220},
            {"id": "oq", "type": "sortie", "signal": "iq", "x": 500, "y": 320},
        ],
        "liens": [
            {"de": "a", "vers": "ss", "entree": 0},
            {"de": "b", "vers": "ss", "entree": 1},
            {"de": "ss", "vers": "oi", "sortie": 0},
            {"de": "ss", "vers": "oq", "sortie": 1},
        ],
    }


def _plat_schema() -> dict:
    return {
        "nom": "plat_api",
        "parametres": {"tau": 0.05, "k": 2.0},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref", "x": 20, "y": 20},
            {"id": "qr", "type": "consigne", "nom": "Q_ref", "x": 20, "y": 120},
            {"id": "g1", "type": "gain", "K": "k", "x": 150, "y": 20},
            {
                "id": "f",
                "type": "premier_ordre",
                "K": 1,
                "tau": "tau",
                "x": 300,
                "y": 20,
            },
            {"id": "oi", "type": "sortie", "signal": "id", "x": 450, "y": 20},
            {"id": "oq", "type": "sortie", "signal": "iq", "x": 450, "y": 120},
        ],
        "liens": [
            {"de": "pr", "vers": "g1"},
            {"de": "g1", "vers": "f"},
            {"de": "f", "vers": "oi"},
            {"de": "qr", "vers": "oq"},
        ],
    }


def _verifier_erreur_422(r, texte: bool = True) -> str:
    assert r.status_code == 422, (r.status_code, r.text[:300])
    corps = r.json()
    detail = corps["detail"]
    assert "Traceback" not in r.text
    if not texte:  # corps mal forme : la validation pydantic rend une liste
        return str(detail)
    assert isinstance(detail, str) and detail.strip(), corps
    return detail


# =============================================================== types
def test_types_avec_categories_entrees_sorties_champs():
    c = TestClient(app)
    r = c.get("/api/blocs/types")
    assert r.status_code == 200
    d = r.json()
    types = d["types"]
    for nom in (
        "pid",
        "avance_retard",
        "second_ordre",
        "fonction_transfert",
        "saturation",
        "implicite",
        "sous_systeme",
        "port_entree",
        "port_sortie",
        "somme",
        "gain",
    ):
        assert nom in types, nom
    for nom, t in types.items():
        assert t["categorie"] in CATEGORIES, nom
        assert "entrees" in t and "sorties" in t and "champs" in t, nom
    assert {t["categorie"] for t in types.values()} == CATEGORIES
    assert {"Kp", "Ki", "Kd", "Tf"} <= set(types["pid"]["champs"])
    assert {"K", "T1", "T2"} <= set(types["avance_retard"]["champs"])
    assert {"K", "wn", "zeta"} <= set(types["second_ordre"]["champs"])
    assert {"num", "den"} <= set(types["fonction_transfert"]["champs"])
    assert {"min", "max"} <= set(types["saturation"]["champs"])
    assert "expr" in types["implicite"]["champs"]
    assert {"index", "nom"} <= set(types["port_entree"]["champs"])
    assert {"index", "nom"} <= set(types["port_sortie"]["champs"])
    assert {"gfl_complet", "gfm_complet"} <= set(d["exemples"])


# ============================================================= compiler
def test_compiler_un_schema_a_sous_systemes(client):
    s = _ss_schema()
    r = client.post("/api/blocs/compiler", json={"schema_bloc": s})
    assert r.status_code == 200, r.text
    d = r.json()
    assert d["equations"]
    ra = client.post(
        "/api/blocs/compiler",
        json={
            "schema_bloc": client.post(
                "/api/blocs/aplatir", json={"schema_bloc": s}
            ).json()["schema"]
        },
    )
    assert ra.status_code == 200
    assert ra.json() == d  # compiler(s) et compiler(aplatir(s)) : meme reponse


def test_compiler_erreur_422(client):
    temoin = client.post("/api/blocs/compiler", json={"schema_bloc": _ss_schema()})
    assert temoin.status_code == 200, temoin.text  # le schema sain compile
    s = _ss_schema()
    s["liens"].pop(0)  # entree 0 du sous-systeme non reliee
    _verifier_erreur_422(client.post("/api/blocs/compiler", json={"schema_bloc": s}))
    s = _plat_schema()
    s["blocs"].append({"id": "p", "type": "port_entree", "index": 0, "nom": "u"})
    _verifier_erreur_422(client.post("/api/blocs/compiler", json={"schema_bloc": s}))


# ======================================================== depuis_element
def test_depuis_element_par_classe(client):
    r = client.post("/api/blocs/depuis_element", json={"classe": "GFL1"})
    assert r.status_code == 200, r.text
    s = r.json()["schema"]
    assert s["amorce"]["classe"] == "GFL1"
    c = client.post("/api/blocs/compiler", json={"schema_bloc": s})
    assert c.status_code == 200, c.text


def test_depuis_element_detail_et_parametres(client):
    r1 = client.post("/api/blocs/depuis_element", json={"classe": "GFM2"})
    r2 = client.post(
        "/api/blocs/depuis_element",
        json={"classe": "GFM2", "detail": "equations", "parametres": {"mp": 0.0123}},
    )
    assert r1.status_code == r2.status_code == 200
    s1, s2 = r1.json()["schema"], r2.json()["schema"]
    assert s1 != s2
    assert s2["parametres"]["mp"] == pytest.approx(0.0123)
    n1 = sum(b["type"] == "sous_systeme" for b in s1["blocs"])
    n2 = sum(b["type"] == "sous_systeme" for b in s2["blocs"])
    # detail "blocs" : un sous-systeme par equation d'etat ; "equations" les
    # remplace par un bloc fonction chacune
    assert n1 > n2


def test_depuis_element_d_un_reseau_avec_ses_parametres(client):
    t = copy.deepcopy(cat.charger("smib_gfl"))
    kipll = {}
    for tau_pll in (None, 0.1):
        tt = copy.deepcopy(t)
        if tau_pll is not None:
            tt.GFL_PARAMS = dict(tt.GFL_PARAMS, tau_PLL=tau_pll)
        el = next(e for e in assembler(tt).all_elements if e.name == "GFL2")
        r = client.post(
            "/api/blocs/depuis_element",
            json={"tables": tables_vers_json(tt), "table": "T_GFL2s", "index": 0},
        )
        assert r.status_code == 200, r.text
        s = r.json()["schema"]
        assert s["amorce"]["classe"] == "GFL2"
        # le schema porte les parametres de CET element, pas ceux par defaut
        assert s["parametres"]["Kipll"] == pytest.approx(el.pdict["Kipll"], rel=1e-12)
        assert s["parametres"]["R1_pu"] == pytest.approx(el.pdict["R1_pu"], rel=1e-12)
        kipll[tau_pll] = s["parametres"]["Kipll"]
    assert kipll[None] != pytest.approx(kipll[0.1], rel=1e-3)


@pytest.mark.parametrize(
    "corps",
    [
        {"classe": "ClasseQuiNExistePas"},
        {"classe": "GFM3DC"},
        {"classe": "MMC"},
        {"classe": "GFL1", "detail": "n_importe_quoi"},
        {},
        {"table": "T_GFL2s", "index": 0},  # pas de reseau
    ],
    ids=[
        "inconnue",
        "port_continu_gfm3dc",
        "port_continu_mmc",
        "detail",
        "vide",
        "sans_tables",
    ],
)
def test_depuis_element_erreurs_422(client, corps):
    mal_forme = "classe" not in corps and "tables" not in corps
    detail = _verifier_erreur_422(
        client.post("/api/blocs/depuis_element", json=corps), texte=not mal_forme
    )
    if corps.get("classe") in ("GFM3DC", "MMC"):
        assert "continu" in detail


def test_depuis_element_ligne_hors_table_422(client):
    t = tables_vers_json(copy.deepcopy(cat.charger("smib_gfl")))
    for corps in (
        {"tables": t, "table": "T_GFL2s", "index": 7},
        {
            "tables": t,
            "table": "T_Lines_SI",
            "index": 0,
        },  # une ligne n'a pas de modele a blocs
        {"tables": t, "table": "T_Nexiste_Pas", "index": 0},
    ):
        _verifier_erreur_422(client.post("/api/blocs/depuis_element", json=corps))


# ====================================== regrouper, degrouper, aplatir
def test_regrouper_degrouper_aplatir(client):
    s = _plat_schema()
    r = client.post(
        "/api/blocs/regrouper",
        json={"schema_bloc": s, "ids": ["g1", "f"], "nom": "chaine_k_f"},
    )
    assert r.status_code == 200, r.text
    reg = r.json()["schema"]
    ss = [b for b in reg["blocs"] if b["type"] == "sous_systeme"]
    assert len(ss) == 1 and ss[0]["nom"] == "chaine_k_f"
    # le serveur et la fonction Python donnent le meme schema
    assert reg == B.regrouper(copy.deepcopy(s), ["g1", "f"], "chaine_k_f")

    r = client.post(
        "/api/blocs/degrouper", json={"schema_bloc": reg, "id": ss[0]["id"]}
    )
    assert r.status_code == 200, r.text
    assert schema_sans_position(r.json()["schema"]) == schema_sans_position(s)

    r = client.post("/api/blocs/aplatir", json={"schema_bloc": reg})
    assert r.status_code == 200, r.text
    plat = r.json()["schema"]
    assert not {b["type"] for b in plat["blocs"]} & {
        "sous_systeme",
        "port_entree",
        "port_sortie",
    }
    assert f"{ss[0]['id']}__f" in {b["id"] for b in plat["blocs"]}


def test_regrouper_erreurs_422(client):
    s = _plat_schema()
    for corps in (
        {"schema_bloc": s, "ids": ["inconnu"], "nom": "x"},
        {"schema_bloc": s, "ids": [], "nom": "x"},
    ):
        _verifier_erreur_422(client.post("/api/blocs/regrouper", json=corps))


def test_degrouper_erreurs_422(client):
    s = _plat_schema()
    for ident in ("g1", "inconnu"):  # pas un sous-systeme / inexistant
        _verifier_erreur_422(
            client.post("/api/blocs/degrouper", json={"schema_bloc": s, "id": ident})
        )


def test_aplatir_erreur_422(client):
    s = _plat_schema()
    s["blocs"].append({"id": "p", "type": "port_sortie", "index": 0, "nom": "u"})
    _verifier_erreur_422(client.post("/api/blocs/aplatir", json={"schema_bloc": s}))


# ============================================================ bibliotheque
def _temoin(client) -> None:
    """Les routes de la bibliotheque existent : un modele valide s'enregistre et
    se relit (sinon un 404 de route absente passerait pour un refus de nom)."""
    s = _plat_schema()
    assert (
        client.post("/api/blocs/bibliotheque", json={"schema_bloc": s}).status_code
        == 200
    )
    assert client.get("/api/blocs/bibliotheque/plat_api").json()["schema"] == s


def _dossier_modeles():
    return cat.dossier_sauvegardes() / "modeles"


def test_bibliotheque_enregistrer_rouvrir_lister_supprimer(client):
    assert client.get("/api/blocs/bibliotheque").json() == []
    s = _ss_schema()
    r = client.post("/api/blocs/bibliotheque", json={"schema_bloc": s})
    assert r.status_code == 200, r.text
    assert r.json() == {"nom": "modele_api"}
    # un fichier JSON par modele, dans modeles/ du dossier des sauvegardes
    assert [p.name for p in _dossier_modeles().iterdir()] == ["modele_api.json"]

    liste = client.get("/api/blocs/bibliotheque").json()
    assert [m["nom"] for m in liste] == ["modele_api"]
    assert {"nom", "titre", "modifie"} <= set(liste[0])
    assert isinstance(liste[0]["titre"], str)

    r = client.get("/api/blocs/bibliotheque/modele_api")
    assert r.status_code == 200
    assert r.json()["schema"] == s  # enregistrer puis rouvrir : le meme schema

    r = client.delete("/api/blocs/bibliotheque/modele_api")
    assert r.status_code == 200 and r.json() == {"nom": "modele_api"}
    assert client.get("/api/blocs/bibliotheque").json() == []
    assert not list(_dossier_modeles().glob("*.json"))


def test_bibliotheque_plusieurs_modeles_et_ecrasement(client):
    a, b = _ss_schema(), _plat_schema()
    for s in (a, b):
        assert (
            client.post("/api/blocs/bibliotheque", json={"schema_bloc": s}).status_code
            == 200
        )
    assert sorted(p.name for p in _dossier_modeles().glob("*.json")) == [
        "modele_api.json",
        "plat_api.json",
    ]
    assert sorted(m["nom"] for m in client.get("/api/blocs/bibliotheque").json()) == [
        "modele_api",
        "plat_api",
    ]
    b2 = copy.deepcopy(b)
    b2["parametres"]["tau"] = 0.2
    assert (
        client.post("/api/blocs/bibliotheque", json={"schema_bloc": b2}).status_code
        == 200
    )
    assert len(client.get("/api/blocs/bibliotheque").json()) == 2
    assert client.get("/api/blocs/bibliotheque/plat_api").json()["schema"] == b2


def test_bibliotheque_partagee_entre_reseaux(client, monkeypatch, tmp_path):
    """Le modele enregistre reste la quand le reseau change : la bibliotheque
    vit dans le dossier des sauvegardes, pas dans un reseau."""
    s = _plat_schema()
    assert (
        client.post("/api/blocs/bibliotheque", json={"schema_bloc": s}).status_code
        == 200
    )
    neuf = TestClient(app)  # un autre client, sans etat commun
    assert neuf.get("/api/blocs/bibliotheque/plat_api").json()["schema"] == s


@pytest.mark.parametrize(
    "nom",
    ["a" * 65, "a" * 300],
    ids=["65_caracteres", "300_caracteres"],
)
def test_bibliotheque_nom_invalide_refuse_a_l_enregistrement(client, nom):
    s = _plat_schema()
    s["nom"] = nom
    _verifier_erreur_422(
        client.post("/api/blocs/bibliotheque", json={"schema_bloc": s})
    )
    dossier = _dossier_modeles()
    assert not dossier.exists() or not list(dossier.iterdir())


@pytest.mark.parametrize("nom", ["a%20b", "..", "a.b", "%2e%2e"])
def test_bibliotheque_nom_invalide_refuse_en_lecture_et_suppression(client, nom):
    _temoin(client)
    for methode in (client.get, client.delete):
        r = methode(f"/api/blocs/bibliotheque/{nom}")
        assert r.status_code in (404, 422), (methode.__name__, nom, r.status_code)
        assert "Traceback" not in r.text
        assert r.status_code != 200


def test_bibliotheque_nom_traversee_de_dossier_refuse(client, tmp_path):
    """Un nom avec des separateurs de chemin ne sort jamais du dossier."""
    _temoin(client)
    secret = tmp_path / "secret.json"
    secret.write_text('{"schema": {"x": 1}}')
    for chemin in ("..%2Fsecret", "..%2F..%2Fsecret", "%2E%2E%2Fsecret"):
        r = client.get(f"/api/blocs/bibliotheque/{chemin}")
        assert r.status_code in (404, 422), (chemin, r.status_code)
        r = client.delete(f"/api/blocs/bibliotheque/{chemin}")
        assert r.status_code in (404, 422), (chemin, r.status_code)
    assert secret.exists()


def test_bibliotheque_modele_absent(client):
    _temoin(client)
    r = client.get("/api/blocs/bibliotheque/inexistant")
    assert r.status_code in (404, 422) and "Traceback" not in r.text
    r = client.delete("/api/blocs/bibliotheque/inexistant")
    assert r.status_code in (404, 422) and "Traceback" not in r.text


def test_bibliotheque_corps_invalide_422(client):
    _verifier_erreur_422(
        client.post("/api/blocs/bibliotheque", json={"schema_bloc": {}})
    )


# ====================================================== garde du serveur local
ROUTES = [
    ("get", "/api/blocs/types", None),
    ("post", "/api/blocs/compiler", {"schema_bloc": {}}),
    ("post", "/api/blocs/depuis_element", {"classe": "GFL1"}),
    ("post", "/api/blocs/regrouper", {"schema_bloc": {}, "ids": [], "nom": "x"}),
    ("post", "/api/blocs/degrouper", {"schema_bloc": {}, "id": "x"}),
    ("post", "/api/blocs/aplatir", {"schema_bloc": {}}),
    ("get", "/api/blocs/bibliotheque", None),
    ("get", "/api/blocs/bibliotheque/x", None),
    ("post", "/api/blocs/bibliotheque", {"schema_bloc": {}}),
    ("delete", "/api/blocs/bibliotheque/x", None),
]


@pytest.mark.parametrize(
    "methode, url, corps", ROUTES, ids=[f"{m}-{u}" for m, u, _ in ROUTES]
)
def test_les_routes_sont_derriere_le_jeton(monkeypatch, tmp_path, methode, url, corps):
    jeton = "jeton-de-test-0123456789"
    monkeypatch.setenv("G2EGRID_JETON", jeton)
    monkeypatch.setenv("G2EGRID_HOTES", "testserver")
    monkeypatch.setattr(cat, "DOSSIER_RESEAUX", tmp_path)
    c = TestClient(app)
    kw = {} if corps is None else {"json": corps}
    assert getattr(c, methode)(url, **kw).status_code == 403
    c.cookies.set(COOKIE_JETON, jeton)
    r = getattr(c, methode)(url, **kw)
    assert r.status_code != 403
    # la route existe : un 404 n'est admis que pour le modele absent (.../x)
    assert r.status_code != 405
    assert r.status_code != 404 or url.endswith("/x"), "route absente"
