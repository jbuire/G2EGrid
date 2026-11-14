"""Les tests qui s'appuient sur un script d'``etudes/``.

La version publiee de l'outil ne distribue pas les etudes de recherche : ces
tests s'y SAUTENT, avec la raison, au lieu d'echouer a l'import.
"""

from __future__ import annotations

import pathlib
import sys

import pytest

RACINE = pathlib.Path(__file__).resolve().parents[1]


def dossier_etude(nom: str, sur_le_chemin: bool = True) -> pathlib.Path:
    """Le dossier ``etudes/<nom>``, ajoute a ``sys.path`` ; s'il n'existe pas,
    tout le module de test est saute."""
    d = RACINE / "etudes" / nom
    if not d.is_dir():
        pytest.skip(
            f"etudes/{nom} n'est pas distribue avec cette version de l'outil",
            allow_module_level=True,
        )
    if sur_le_chemin:
        sys.path.insert(0, str(d))
    return d
