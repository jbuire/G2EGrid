"""Recupere MathJax et Mermaid pour une documentation lisible HORS LIGNE.

Sans eux, la documentation charge ces deux bibliotheques depuis un CDN a
l'affichage : c'est le comportement par defaut, et il suffit en ligne. Pour la
consulter sans reseau (ou derriere un pare-feu qui bloque les CDN), lancer une
fois, depuis ``g2egrid/`` ::

    uv run python docs/outils/vendoriser.py

Les fichiers vont dans ``docs/_static/vendor/`` (ignore par git) et
``conf.py`` les prend d'office a la construction suivante. Ils viennent du
registre npm, pas du CDN : c'est le meme contenu, par une autre porte.
"""

from __future__ import annotations

import io
import pathlib
import tarfile
import urllib.request

VENDOR = pathlib.Path(__file__).resolve().parent.parent / "_static" / "vendor"

#: (paquet npm, version, prefixe dans l'archive, dossier local). Mermaid est
#: charge comme module ES par sphinxcontrib-mermaid : il faut sa version
#: ``esm`` ET les fragments qu'elle importe, d'ou un prefixe et non un fichier.
PAQUETS = [
    ("mathjax", "3.2.2", "package/es5/tex-svg-full.js", ""),
    ("mermaid", "11.4.1", "package/dist/", "mermaid"),
]


def main() -> None:
    VENDOR.mkdir(parents=True, exist_ok=True)
    for paquet, version, prefixe, dossier in PAQUETS:
        url = f"https://registry.npmjs.org/{paquet}/-/{paquet}-{version}.tgz"
        with urllib.request.urlopen(url) as r:
            archive = tarfile.open(fileobj=io.BytesIO(r.read()), mode="r:gz")
        n = 0
        for m in archive.getmembers():
            if not m.isfile() or not m.name.startswith(prefixe):
                continue
            if dossier == "mermaid" and (
                "esm.min" not in m.name or m.name.endswith(".map")
            ):
                continue  # seulement la version module minifiee et ses fragments
            rel = pathlib.Path(m.name[len(prefixe) :] or pathlib.Path(m.name).name)
            dst = VENDOR / dossier / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(archive.extractfile(m).read())
            n += 1
        print(f"{paquet} {version} : {n} fichier(s)")


if __name__ == "__main__":
    main()
