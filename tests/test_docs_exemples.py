"""Les exemples de la documentation sont executes tels qu'ecrits.

Chaque page listee est lue, ses blocs ```python executes dans l'ordre, dans
un meme espace de noms, depuis le dossier ``g2egrid/`` -- comme le lecteur le
ferait. Un exemple qui ne tourne plus fait echouer la suite, au lieu de
tromper le lecteur."""

from __future__ import annotations

import os
import pathlib
import re

import pytest

ICI = pathlib.Path(__file__).resolve().parents[1]
DOCS = ICI / "docs"
PAGES = [
    "exemples.md",
    "guide/impedance.md",
    "guide/residus.md",
    "guide/continuation.md",
    "guide/grands_reseaux.md",
    "guide/lois_gfm.md",
    "guide/interoperabilite.md",
    "guide/cigre_b4.md",
    "guide/domaines.md",
    "guide/evenements.md",
    "guide/protections_dc.md",
    "guide/mmc.md",
    "guide/automatiser.md",
    "guide/datacenter.md",
    "guide/schemas_blocs.md",
    "guide/equations_modifiables.md",
]
BLOC = re.compile(r"```python\n(.*?)```", re.S)


@pytest.mark.parametrize("page", PAGES)
def test_exemples_de_la_page(page, monkeypatch, tmp_path):
    pytest.importorskip("fastapi")
    # un dossier jetable qui voit les donnees du depot : les exemples qui
    # ecrivent un fichier ne laissent rien dans l'arbre de travail
    for d in ("tests", "etudes", "reseaux", "reseaux_test"):
        if (ICI / d).exists():
            (tmp_path / d).symlink_to(ICI / d)
    monkeypatch.chdir(tmp_path)
    os.environ.pop("G2EGRID_JETON", None)
    texte = (DOCS / page).read_text(encoding="utf-8")
    blocs = BLOC.findall(texte)
    if not (ICI / "etudes").is_dir() and any("etudes/" in b for b in blocs):
        pytest.skip(f"{page} : ses exemples lisent etudes/, non distribue ici")
    assert blocs, f"{page} : aucun exemple"
    espace: dict = {"__name__": "__doc__"}
    for i, code in enumerate(blocs):
        try:
            exec(compile(code, f"{page}[bloc {i + 1}]", "exec"), espace)
        except Exception as err:
            pytest.fail(
                f"{page}, bloc {i + 1} :\n{code}\n-> {type(err).__name__}: {err}"
            )
