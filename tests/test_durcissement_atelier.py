"""Durcissement de l'atelier (revue 2026-10-25) : un test par point.

Chaque test echoue sans la correction qu'il garde :

* une erreur interne ne renvoie JAMAIS de trace Python (500 generique, trace au
  journal) ;
* les entrees invalides sont des 422 en francais qui disent quoi corriger :
  duree, nombre de points, indice de mode, f_max, ligne de fiche, table
  inconnue, catalogue inconnu, bus inexistant, reseau vide, table mal formee ;
* ``/api/valider`` signale un load flow qui ne converge pas ;
* le corps d'une requete est plafonne (413) et les taches en attente aussi
  (429) ;
* le nom de catalogue de ``/api/code`` (et les chaines de l'export) ne
  peuvent pas ajouter de ligne au script ;
* ``/doc/`` n'est plus une 500 avant la construction de la documentation.
"""

from __future__ import annotations

import ast
import json
import logging
import threading

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")

from fastapi.testclient import TestClient  # noqa: E402

from g2egrid.catalogue import charger  # noqa: E402
from g2egrid_api import analyse, app as module_app, taches  # noqa: E402
from g2egrid_api.app import app  # noqa: E402


@pytest.fixture
def client(monkeypatch):
    monkeypatch.delenv("G2EGRID_JETON", raising=False)
    # une 500 doit etre rendue au client, pas relevee par le client de test
    return TestClient(app, raise_server_exceptions=False)


@pytest.fixture(scope="module")
def wscc9():
    return analyse.tables_vers_json(charger("wscc9"))


def _copie(d):
    return json.loads(json.dumps(d))


def _detail(r):
    d = r.json()["detail"]
    assert isinstance(d, str), d  # une phrase, pas la liste pydantic
    return d


# ------------------------------------------------- 500 : jamais de trace
def test_erreur_interne_sans_trace_dans_la_reponse(client, wscc9, monkeypatch, caplog):
    def casse(*a, **k):
        raise ZeroDivisionError("secret-interne /home/prive/fichier.py")

    monkeypatch.setattr(analyse, "loadflow", casse)
    with caplog.at_level(logging.ERROR, logger="g2egrid_api"):
        r = client.post("/api/loadflow", json={"tables": wscc9})
    assert r.status_code == 500
    assert "Traceback" not in r.text and "secret-interne" not in r.text
    assert "fichier.py" not in r.text and "ZeroDivisionError" not in r.text
    assert "reference" in _detail(r) and "journal" in _detail(r)
    # la trace est au journal, avec la meme reference
    ref = _detail(r).split("reference ")[1].split(")")[0]
    assert any(
        ref in m.getMessage() and m.exc_info and "secret-interne" in str(m.exc_info[1])
        for m in caplog.records
    )


def test_erreur_hors_calcul_aussi_generique(client, monkeypatch):
    """Une route qui ne passe pas par ``_calcul`` : meme 500 generique."""

    def casse():
        raise RuntimeError("secret-catalogue")

    monkeypatch.setattr(analyse, "catalogue", casse)
    r = client.get("/api/reseaux")
    assert r.status_code == 500
    assert "secret-catalogue" not in r.text and "Traceback" not in r.text


def test_tache_en_echec_sans_trace(client, wscc9, monkeypatch):
    def casse(*a, **k):
        raise TypeError("secret-tache /chemin/prive.py")

    monkeypatch.setattr(analyse, "continuation_parametre", casse)
    corps = {
        "tables": wscc9,
        "element": "SM1",
        "parametre": "H",
        "debut": 1,
        "fin": 2,
        "pas": 3,
    }
    d = client.post("/api/taches/continuation", json=corps).json()
    t = taches.lire(d["id"])
    t.arret.wait(0)  # laisse la main au fil
    for _ in range(300):
        if t.terminee:
            break
        threading.Event().wait(0.02)
    assert t.etat == "echec"
    assert "secret-tache" not in t.erreur and "Traceback" not in t.erreur
    assert "reference" in t.erreur


# ---------------------------------------------- 422 : entrees invalides
@pytest.mark.parametrize(
    "corps, mot",
    [
        ({"duree": 0}, "duree"),
        ({"duree": -1}, "duree"),
        ({"n": 1}, " n "),
        ({"n": 0}, " n "),
        ({"pas": 0}, "pas"),
    ],
)
def test_temporel_bornes(client, wscc9, corps, mot):
    r = client.post(
        "/api/temporel", json={"tables": wscc9, "entree": "SM1.p_ref", **corps}
    )
    assert r.status_code == 422
    d = _detail(r)
    assert mot.strip() in d and "Traceback" not in d


@pytest.mark.parametrize("indice", [-1, -5])
def test_mode_indice_negatif_refuse(client, wscc9, indice):
    r = client.post("/api/mode", json={"tables": wscc9, "indice": indice})
    assert r.status_code == 422 and "indice" in _detail(r)


def test_mode_indice_trop_grand(client, wscc9):
    r = client.post("/api/mode", json={"tables": wscc9, "indice": 9999})
    assert r.status_code == 422
    d = _detail(r)
    assert "9999" in d and "de 0 a" in d


@pytest.mark.parametrize("f_max", [0, -1, -0.5])
def test_modal_f_max_strictement_positif(client, wscc9, f_max):
    r = client.post("/api/modal", json={"tables": wscc9, "f_max": f_max})
    assert r.status_code == 422 and "f_max" in _detail(r)


def test_fiche_index_hors_table(client, wscc9):
    r = client.post("/api/fiche", json={"tables": wscc9, "table": "T_SMs", "index": 99})
    assert r.status_code == 422
    assert "T_SMs" in _detail(r) and "99" in _detail(r)
    r = client.post("/api/fiche", json={"tables": wscc9, "table": "T_SMs", "index": -1})
    assert r.status_code == 422


def test_fiche_table_inconnue(client, wscc9):
    r = client.post("/api/fiche", json={"tables": wscc9, "table": "X", "index": 0})
    assert r.status_code == 422
    assert "table inconnue" in _detail(r)


def test_catalogue_inconnu_a_l_import(client, wscc9):
    r = client.post(
        "/api/reseau/importer",
        json={"reseau": {"tables": wscc9}, "source": {"catalogue": "nexistepas"}},
    )
    assert r.status_code == 422
    assert "nexistepas" in _detail(r) and "Traceback" not in r.text


def test_catalogue_inconnu_en_lecture(client):
    r = client.get("/api/reseaux/nexistepas")
    assert r.status_code == 404 and "nexistepas" in _detail(r)


# ------------------------------------- 422 : messages qui disent quoi corriger
def test_bus_inexistant_au_loadflow(client, wscc9):
    t = _copie(wscc9)
    t["T_Lines_pu"][0][1] = 99
    r = client.post("/api/loadflow", json={"tables": t})
    assert r.status_code == 422
    d = _detail(r)
    assert d != "99" and "bus 99" in d and "T_Lines_pu" in d and "Y_DER" in d


def test_reseau_vide(client):
    for route in ("loadflow", "modal"):
        r = client.post(f"/api/{route}", json={"tables": {}})
        assert r.status_code == 422, route
        d = _detail(r)
        assert "reseau vide" in d and "zero-size" not in d


def test_table_mal_formee(client, wscc9):
    t = _copie(wscc9)
    t["T_Lines_pu"][0][0] = "abc"
    r = client.post("/api/loadflow", json={"tables": t})
    assert r.status_code == 422
    d = _detail(r)
    assert "invalid literal" not in d
    assert "T_Lines_pu" in d and "ligne n. 1" in d and "ID" in d and "'abc'" in d


def test_table_trop_courte_n_est_pas_une_500(client, wscc9):
    t = _copie(wscc9)
    t["Y_DER"] = [[1, 2]]
    r = client.post("/api/loadflow", json={"tables": t})
    assert r.status_code == 422 and "au lieu de 9" in _detail(r)


def test_nombre_remplace_par_du_texte(client, wscc9):
    t = _copie(wscc9)
    t["T_Lines_pu"][0][5] = "grand"
    r = client.post("/api/loadflow", json={"tables": t})
    assert r.status_code == 422 and "X_pu" in _detail(r)


def test_corps_de_mauvais_type_en_francais(client, wscc9):
    r = client.post("/api/modal", json={"tables": wscc9, "f_max": "vite"})
    assert r.status_code == 422
    assert "f_max" in _detail(r) and "nombre" in _detail(r)


def _charges_x200(wscc9):
    t = _copie(wscc9)
    t["Y_DER"] = [[*r[:7], r[7] * 200, r[8] * 200, *r[9:]] for r in t["Y_DER"]]
    return t


def test_valider_signale_un_loadflow_qui_diverge(client, wscc9):
    t = _charges_x200(wscc9)
    # le calcul refuse...
    assert client.post("/api/loadflow", json={"tables": t}).status_code == 422
    # ... et la validation le dit au lieu de repondre « valide »
    d = client.post("/api/valider", json={"tables": t}).json()
    assert d["ok"] is False
    pbs = [p for p in d["problemes"] if "converge" in p["message"]]
    assert pbs and pbs[0]["gravite"] == "erreur"
    # un reseau sain reste valide
    assert client.post("/api/valider", json={"tables": wscc9}).json()["ok"] is True


# ----------------------------------------------- 413 et 429 : plafonds
def test_corps_trop_gros_413(client, wscc9, monkeypatch):
    monkeypatch.setenv("G2EGRID_TAILLE_MAX", "20000")
    r = client.post("/api/loadflow", json={"tables": {"nom": "x" * 50000}})
    assert r.status_code == 413
    assert "G2EGRID_TAILLE_MAX" in _detail(r) and "20000" in _detail(r)
    # sous le plafond, la requete passe
    monkeypatch.setenv("G2EGRID_TAILLE_MAX", str(10_000_000))
    assert client.post("/api/loadflow", json={"tables": wscc9}).status_code == 200


def test_corps_trop_gros_sans_longueur_annoncee(client, monkeypatch):
    """Un corps par morceaux n'a pas de Content-Length : les octets sont comptes."""
    monkeypatch.setenv("G2EGRID_TAILLE_MAX", "5000")

    def morceaux():
        for _ in range(20):
            yield b'{"tables": {"nom": "' if _ == 0 else b"x" * 1000

    r = client.post(
        "/api/loadflow",
        content=morceaux(),
        headers={"Content-Type": "application/json"},
    )
    assert r.status_code == 413


def test_le_plafond_laisse_passer_le_plus_gros_reseau_du_catalogue():
    """Marge d'au moins x20 sur le plus gros reseau du catalogue, a plat."""
    plus_gros = 0
    for f in analyse.catalogue():
        r = analyse.depuis_catalogue(f["nom"])
        corps = json.dumps({"tables": r["tables"], "positions": r["positions"]})
        plus_gros = max(plus_gros, len(corps))
    assert plus_gros > 50_000  # ieee39_mixte : de l'ordre de 83 ko
    assert module_app.TAILLE_MAX_DEFAUT >= 20 * plus_gros


def test_trop_de_taches_429(client, wscc9, monkeypatch):
    monkeypatch.setenv("G2EGRID_TACHES_MAX", "2")
    libere = threading.Event()
    lancees = []
    try:
        for _ in range(2):
            lancees.append(taches.lancer("attente", lambda: libere.wait(30)))
        avant = len(taches.liste())
        r = client.post("/api/taches/robustesse", json={"tables": wscc9})
        assert r.status_code == 429
        d = _detail(r)
        assert "trop de calculs" in d and "Arreter" in d
        # le refus ne cree aucune tache
        assert len(taches.liste()) == avant
    finally:
        libere.set()
        for t in lancees:
            for _ in range(300):
                if t.terminee:
                    break
                threading.Event().wait(0.02)
    # une fois les taches finies, on peut a nouveau lancer
    monkeypatch.setenv("G2EGRID_TACHES_MAX", "2")
    t = taches.lancer("court", lambda: 1)
    for _ in range(300):
        if t.terminee:
            break
        threading.Event().wait(0.02)
    assert t.etat == "fini"


# ----------------------------------------------- /api/code : aucune ligne libre
@pytest.mark.parametrize(
    "nom",
    [
        "x\nimport os; os.system('id')",
        "wscc9')\nimport os\n#",
        "a'); import os; ('",
        "../../etc/passwd",
        "",
        "wscc9\n",
    ],
)
def test_code_refuse_un_nom_hostile(client, wscc9, nom):
    r = client.post("/api/code", json={"reseau": {"tables": wscc9}, "nom": nom})
    if nom == "":  # pas de nom : le script embarque ses donnees
        assert r.status_code == 200
        return
    assert r.status_code == 422, r.text
    assert "invalide" in _detail(r)


def test_code_nom_valide_en_litteral(client, wscc9):
    r = client.post("/api/code", json={"reseau": {"tables": wscc9}, "nom": "wscc9"})
    assert r.status_code == 200
    src = r.json()["python"]
    assert "Atelier.catalogue('wscc9')" in src
    ast.parse(src)  # le script reste du Python valide


def test_export_ne_laisse_pas_une_chaine_ajouter_du_code(wscc9):
    """Le nom du reseau et les chaines des tables sont ecrits en litteraux : ni
    guillemet ni saut de ligne ne sort du litteral."""
    t = analyse.tables_depuis_json(_copie(wscc9))
    t.nom = 'x"\nimport os\nos.system("id")\n#'
    t.T_Loads = [(*r[:4], "Z'\nimport sys\n#") for r in t.T_Loads]
    src = analyse.exporter(t)
    arbre = ast.parse(src)
    assert not any(isinstance(n, (ast.Import, ast.ImportFrom)) for n in ast.walk(arbre))
    assert not any(isinstance(n, ast.Call) for n in ast.walk(arbre)) or all(
        getattr(n.func, "id", "") in ("dict", "")
        for n in ast.walk(arbre)
        if isinstance(n, ast.Call)
    )
    # le nom revient intact en relisant le script
    espace: dict = {}
    exec(compile(src, "export", "exec"), {"np": __import__("numpy")}, espace)
    assert espace["NOM_RESEAU"] == t.nom


def test_export_matlab_nom_sur_une_seule_ligne(wscc9):
    from g2egrid.export import vers_matlab

    t = analyse.tables_depuis_json(_copie(wscc9))
    t.nom = "wscc9\nsystem('id')"
    source, _ = vers_matlab(t)
    assert not any(l.startswith("system(") for l in source.splitlines())


# ---------------------------------------------------------------- /doc
def test_doc_avant_construction_est_une_page_claire(client, monkeypatch, tmp_path):
    monkeypatch.setattr(module_app, "DOC_HTML", tmp_path / "pas_construite")
    for chemin in ("/doc/", "/doc/index.html", "/doc/guide/x.html", "/doc"):
        r = client.get(chemin, follow_redirects=True)
        assert r.status_code == 404, chemin
        assert "pas encore construite" in r.text, chemin
        assert "Internal Server Error" not in r.text
