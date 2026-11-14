"""La garde du serveur local : jeton echange contre un cookie, controle de
l'hote, et rien de change quand aucun jeton n'est defini."""

from __future__ import annotations

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")

from fastapi.testclient import TestClient  # noqa: E402

from g2egrid_api import __main__ as lanceur  # noqa: E402
from g2egrid_api.app import COOKIE_JETON, app  # noqa: E402

JETON = "jeton-de-test-0123456789"


@pytest.fixture
def client(monkeypatch):
    monkeypatch.setenv("G2EGRID_JETON", JETON)
    monkeypatch.setenv("G2EGRID_HOTES", "testserver")
    return TestClient(app)


def test_sans_jeton_rien_ne_change(monkeypatch):
    monkeypatch.delenv("G2EGRID_JETON", raising=False)
    c = TestClient(app)
    assert c.get("/api/defauts").status_code == 200


def test_refus_sans_cookie(client):
    assert client.get("/api/defauts").status_code == 403
    assert client.get("/").status_code == 403
    assert client.get("/web/app.js").status_code == 403


def test_jeton_echange_contre_cookie(client):
    r = client.get(f"/?jeton={JETON}", follow_redirects=False)
    assert r.status_code == 303 and r.headers["location"] == "/"
    cookie = r.headers["set-cookie"]
    assert "httponly" in cookie.lower() and "samesite=strict" in cookie.lower()
    client.cookies.set(COOKIE_JETON, JETON)
    assert client.get("/api/defauts").status_code == 200
    assert client.get("/").status_code == 200


def test_mauvais_jeton(client):
    assert client.get("/?jeton=faux", follow_redirects=False).status_code == 403
    client.cookies.set(COOKIE_JETON, "faux")
    assert client.get("/api/defauts").status_code == 403


def test_jeton_non_ascii_refuse_sans_erreur(client):
    """Un jeton non ASCII est un mauvais jeton : 403, pas une erreur 500."""
    r = client.get("/?jeton=%C3%A9", follow_redirects=False)
    assert r.status_code == 403


def test_en_tete_pour_les_scripts(client):
    r = client.get("/api/defauts", headers={"X-G2egrid-Jeton": JETON})
    assert r.status_code == 200


def test_hote_etranger_refuse(client):
    """Rebinding DNS : un nom qui pointe sur 127.0.0.1 n'est pas admis."""
    r = client.get(
        "/api/defauts",
        headers={"X-G2egrid-Jeton": JETON, "Host": "attaquant.example:8765"},
    )
    assert r.status_code == 403


def test_fichier_jeton_prive(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    lanceur._ecrire_jeton(9999, "abc")
    f = lanceur.fichier_jeton(9999)
    assert f.read_text() == "abc"
    assert (f.stat().st_mode & 0o077) == 0
    assert lanceur._lire_jeton(9999) == "abc"
