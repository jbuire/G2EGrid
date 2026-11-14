"""La reorganisation du paquet en sous-paquets (2026-10-28) : les anciens
chemins de modules restent valides, et les couches ne s'importent que vers le
bas.

* chaque ancien nom (``g2egrid.modal``, ``g2egrid.network``...) donne le MEME
  objet module que le nouveau (``g2egrid.analyse.modal``...), par
  ``import``, ``from g2egrid import x`` et ``from g2egrid.x import nom`` ;
* un ``monkeypatch`` sur l'ancien nom agit sur le module reel ;
* le module garde sa vraie identite (``__name__``, ``__spec__``) ;
* ``G2EGRID_AVERTIR_ANCIENS_CHEMINS`` fait avertir a l'import d'un ancien nom ;
* ``g2egrid_api.analyse``, devenu un paquet, expose tous ses anciens noms ;
* aucun module d'une couche basse n'importe, au niveau du module, une couche
  plus haute (les imports tardifs, dans une fonction, sont admis).
"""

from __future__ import annotations

import ast
import importlib
import pathlib
import subprocess
import sys
import warnings

import pytest

import g2egrid
from g2egrid import _compat

RACINE = pathlib.Path(g2egrid.__file__).resolve().parent

#: l'ordre des couches : un module n'importe que sa couche ou une plus basse
COUCHES = [
    "noyau",
    "elements",
    "reseau",
    ("temporel", "analyse"),
    "scenarios",
    "multi_niveau",
    ("formats", "validation"),
]


def _rang(nom: str) -> int | None:
    """Le rang de la couche d'un module ``g2egrid.<couche>...``."""
    parts = nom.split(".")
    if len(parts) < 2 or parts[0] != "g2egrid":
        return None
    for k, c in enumerate(COUCHES):
        if parts[1] in (c if isinstance(c, tuple) else (c,)):
            return k
    return None


@pytest.mark.parametrize("ancien", sorted(_compat.ANCIENS))
def test_ancien_nom_meme_module(ancien):
    nouveau = _compat.ANCIENS[ancien]
    m_ancien = importlib.import_module(ancien)
    m_neuf = importlib.import_module(nouveau)
    assert m_ancien is m_neuf
    assert m_neuf.__name__ == nouveau
    assert m_neuf.__spec__.name == nouveau
    assert sys.modules[ancien] is m_neuf


def test_from_g2egrid_import_ancien_module():
    from g2egrid import catalogue, modal, network
    from g2egrid.analyse import modal as modal2
    from g2egrid.reseau import catalogue as catalogue2
    from g2egrid.reseau import network as network2

    assert modal is modal2 and network is network2 and catalogue is catalogue2
    assert g2egrid.tables is importlib.import_module("g2egrid.reseau.tables")


def test_from_ancien_module_import_nom():
    from g2egrid.reseau.tables import assembler as a2
    from g2egrid.tables import assembler

    assert assembler is a2


def test_monkeypatch_par_l_ancien_nom(monkeypatch):
    import g2egrid.reseau.catalogue as reel
    from g2egrid import catalogue

    monkeypatch.setattr(catalogue, "DOSSIER_RESEAUX", pathlib.Path("/ailleurs"))
    assert reel.DOSSIER_RESEAUX == pathlib.Path("/ailleurs")


def test_avertissement_a_la_demande():
    """Dans un interpreteur neuf : l'avertissement n'existe qu'a l'import."""
    code = (
        "import warnings; warnings.simplefilter('error')\n"
        "try:\n"
        "    import g2egrid.modal\n"
        "except DeprecationWarning as e:\n"
        "    print('AVERTI', e)\n"
    )
    env = {"G2EGRID_AVERTIR_ANCIENS_CHEMINS": "1", "PATH": "/usr/bin:/bin"}
    sortie = subprocess.run(
        [sys.executable, "-c", code],
        cwd=RACINE.parent,
        env=env,
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    assert "AVERTI g2egrid.modal s'appelle maintenant g2egrid.analyse.modal" in sortie


def test_pas_d_avertissement_par_defaut():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        importlib.import_module("g2egrid.residus")


def test_analyse_de_l_atelier_garde_ses_noms():
    from g2egrid_api import analyse

    for nom in analyse.__all__:
        assert hasattr(analyse, nom), nom
    # prives utilises par les tests et l'API
    for nom in ("_tables_de", "_CACHE", "Assemble", "construire", "vider_cache"):
        assert hasattr(analyse, nom), nom


def _imports_de_module(fichier: pathlib.Path) -> list[str]:
    arbre = ast.parse(fichier.read_text(encoding="utf-8"))
    noms = []
    for n in arbre.body:  # niveau module seulement
        if isinstance(n, ast.ImportFrom) and n.module and n.level == 0:
            noms.append(n.module)
        elif isinstance(n, ast.Import):
            noms += [a.name for a in n.names]
    return noms


def test_les_couches_n_importent_que_vers_le_bas():
    montees = []
    for f in sorted(RACINE.rglob("*.py")):
        rel = f.relative_to(RACINE.parent).with_suffix("")
        nom = ".".join(p for p in rel.parts if p != "__init__")
        r = _rang(nom)
        if r is None:
            continue
        for cible in _imports_de_module(f):
            rc = _rang(cible)
            if rc is not None and rc > r:
                montees.append(f"{nom} -> {cible}")
    assert not montees, montees
