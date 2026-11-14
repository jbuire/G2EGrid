"""Les taches de fond de l'atelier : avancement venu du moteur, flux SSE,
annulation au prochain point de controle."""

from __future__ import annotations

import json
import time

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")

from fastapi.testclient import TestClient  # noqa: E402

from g2egrid import progres  # noqa: E402
from g2egrid.catalogue import charger  # noqa: E402
from g2egrid_api import analyse, taches  # noqa: E402
from g2egrid_api.app import app  # noqa: E402


def _attendre(t, delai=60):
    fin = time.time() + delai
    while not t.terminee and time.time() < fin:
        time.sleep(0.05)
    return t


def test_progres_sans_suivi_ne_fait_rien():
    progres.signaler(0.5, "rien")


def test_tache_suit_et_finit():
    vus = []

    def calcul():
        for i in range(5):
            progres.signaler(i / 5, f"pas {i}")
            vus.append(i)
        return {"ok": True}

    t = _attendre(taches.lancer("essai", calcul))
    assert t.etat == "fini" and t.resultat == {"ok": True} and t.fraction == 1.0
    assert vus == list(range(5))


def test_tache_annulee():
    def calcul():
        for i in range(500):
            progres.signaler(i / 500)
            time.sleep(0.01)
        return "fini trop tot"

    t = taches.lancer("longue", calcul)
    time.sleep(0.1)
    taches.annuler(t.id)
    assert _attendre(t).etat == "annule"


def test_tache_en_echec():
    def calcul():
        raise ValueError("reseau sans slack")

    t = _attendre(taches.lancer("echec", calcul))
    assert t.etat == "echec" and "sans slack" in t.erreur


def test_robustesse_par_flux_sse(monkeypatch):
    monkeypatch.delenv("G2EGRID_JETON", raising=False)
    c = TestClient(app)
    corps = {"tables": analyse.tables_vers_json(charger("smib")), "deconnexion": False}
    d = c.post("/api/taches/robustesse", json=corps).json()
    evenements = []
    with c.stream("GET", f"/api/taches/{d['id']}/flux") as r:
        for ligne in r.iter_lines():
            if ligne.startswith("data: "):
                evenements.append(json.loads(ligne[6:]))
    assert evenements[-1]["etat"] == "fini"
    assert "modes" in evenements[-1]["resultat"]
    assert c.get(f"/api/taches/{d['id']}").json()["etat"] == "fini"
    assert c.post("/api/taches/inconnue", json={}).status_code == 404


def test_route_interrompue_est_annulee_pas_en_echec():
    """Une fonction de route (passee par _calcul) arretee en cours : annulee."""
    from g2egrid_api.app import robustesse
    from g2egrid_api.schemas import CorpsRobustesse

    corps = CorpsRobustesse(tables=analyse.tables_vers_json(charger("ieee39")))
    t = taches.lancer("robustesse", lambda: robustesse(corps))
    fin = time.time() + 60
    while t.fraction == 0 and time.time() < fin:
        time.sleep(0.05)
    taches.annuler(t.id)
    assert _attendre(t).etat == "annule"


def test_genres_limites_et_impedances(monkeypatch):
    """Continuation et boucle mineure, par les taches de l'atelier."""
    monkeypatch.delenv("G2EGRID_JETON", raising=False)
    c = TestClient(app)

    def executer(genre, corps):
        d = c.post(f"/api/taches/{genre}", json=corps).json()
        t = _attendre(taches.lire(d["id"]), 300)
        assert t.etat == "fini", t.erreur
        return t.resultat

    gfl = analyse.tables_vers_json(charger("smib_gfl"))
    els = c.post("/api/elements", json={"tables": gfl}).json()
    assert any(e["nom"] == "GFL2" and "Kipll" in e["parametres"] for e in els)
    co = executer(
        "continuation",
        {"tables": gfl, "element": "GFL2", "parametre": "Kipll", "debut": 10,
         "fin": 1000, "pas": 5, "log": True},
    )  # fmt: skip
    assert len(co["mu"]) == 5 and co["trajectoires"]
    imp = executer("impedance", {"tables": gfl, "element": "GFL2", "n": 200})
    assert imp["applicable"] and imp["stable"] and len(imp["lieux"]) == 2
