"""L'atelier en une commande : ``uv run --extra api python -m g2egrid_api``.

Demarre le serveur local et ouvre le navigateur sur la page des qu'il
repond -- c'est ce que lancent ``lancer_atelier.bat`` (double-clic sous
Windows) et ``lancer_atelier.sh``, a la racine du depot. Si le port est deja
pris, c'est presque toujours l'atelier lance par un double-clic precedent :
on ouvre alors seulement la page, plutot que d'echouer sur "address already
in use".

SECURITE -- le serveur n'ecoute que la machine locale (``--hote 127.0.0.1``)
et exige un JETON tire au hasard a chaque lancement : l'adresse ouverte le
porte (``/?jeton=...``), et toute autre page du navigateur, qui ne le connait
pas, se voit refuser l'acces (voir ``g2egrid_api.app``). Le jeton est range
dans un fichier lisible du seul utilisateur, pour qu'un second double-clic
rouvre la page du serveur deja lance. ``--sans-jeton`` retire cette garde.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import secrets
import socket
import threading
import time
import webbrowser


def _repond(hote: str, port: int) -> bool:
    with socket.socket() as s:
        s.settimeout(0.3)
        return s.connect_ex((hote, port)) == 0


def _ouvrir_des_que_pret(url: str, hote: str, port: int, delai: float = 60.0) -> None:
    fin = time.monotonic() + delai
    while time.monotonic() < fin:
        if _repond(hote, port):
            webbrowser.open(url)
            return
        time.sleep(0.25)


def fichier_jeton(port: int) -> pathlib.Path:
    base = os.environ.get("XDG_CACHE_HOME") or pathlib.Path.home() / ".cache"
    return pathlib.Path(base) / "g2egrid" / f"atelier-{port}.jeton"


def _ecrire_jeton(port: int, jeton: str) -> None:
    f = fichier_jeton(port)
    f.parent.mkdir(parents=True, exist_ok=True)
    # cree d'emblee en 0600 : jamais lisible par un autre, meme un instant
    fd = os.open(f, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as h:
        h.write(jeton)


def _lire_jeton(port: int) -> str | None:
    try:
        return fichier_jeton(port).read_text().strip() or None
    except OSError:
        return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m g2egrid_api",
        description="Atelier g2egrid : serveur local et navigateur.",
    )
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--hote", default="127.0.0.1")
    ap.add_argument(
        "--sans-navigateur", action="store_true", help="ne pas ouvrir la page"
    )
    ap.add_argument(
        "--sans-jeton",
        action="store_true",
        help="ne pas exiger de jeton (a eviter hors d'un conteneur isole)",
    )
    a = ap.parse_args(argv)
    hote_url = "127.0.0.1" if a.hote in ("0.0.0.0", "::") else a.hote
    base = f"http://{hote_url}:{a.port}/"
    if _repond(hote_url, a.port):
        jeton = _lire_jeton(a.port)
        url = base + (f"?jeton={jeton}" if jeton else "")
        print(f"le port {a.port} repond deja -- l'atelier tourne sans doute : {url}")
        if not a.sans_navigateur:
            webbrowser.open(url)
        return 0
    if a.sans_jeton:
        os.environ.pop("G2EGRID_JETON", None)
        url = base
    else:
        jeton = os.environ.get("G2EGRID_JETON") or secrets.token_urlsafe(24)
        os.environ["G2EGRID_JETON"] = jeton
        _ecrire_jeton(a.port, jeton)
        url = f"{base}?jeton={jeton}"
    if a.hote not in ("127.0.0.1", "localhost", "::1"):
        print(
            f"ATTENTION : l'atelier ecoute sur {a.hote}, pas seulement la machine "
            "locale" + ("" if a.sans_jeton else " (le jeton reste exige)"),
            flush=True,
        )
    if not a.sans_navigateur:
        threading.Thread(
            target=_ouvrir_des_que_pret, args=(url, hote_url, a.port), daemon=True
        ).start()

    import uvicorn

    print(f"atelier g2egrid : {url}   (Ctrl+C pour arreter)", flush=True)
    uvicorn.run("g2egrid_api.app:app", host=a.hote, port=a.port, log_level="warning")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
