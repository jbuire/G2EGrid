"""La version publiee de l'outil (``g2egrid/``, a la racine du depot).

Elle est GENEREE depuis ``g2egrid/`` par ``developpement/publier_g2egrid.py`` :

* le renommage ne touche que ce qui designe l'outil -- le logiciel MATLAB
  d'origine (G2Elin V1/V2, G2ELIN), le depot tiers et le depot de
  developpement gardent leur nom ;
* un fichier enregistre sous l'un des deux noms de l'outil se lit sous
  l'autre ;
* ``g2egrid/`` est A JOUR : apres un changement dans ce qui est publie, il se
  regenere (ce test le dit, avec la commande).

Dans la version publiee elle-meme, le generateur n'est pas distribue : ces
tests s'y sautent.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib

import pytest

PYTHON = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = PYTHON / "developpement" / "publier_g2egrid.py"
ANCIEN = "g2" + "elin"  # ecrit en deux morceaux : ce fichier est lui-meme publie
NEUF = "g2" + "egrid"

if not SCRIPT.is_file():
    pytest.skip(
        "generateur absent : ceci est la version publiee", allow_module_level=True
    )


@pytest.fixture(scope="module")
def pub():
    spec = importlib.util.spec_from_file_location("publier_g2egrid", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.mark.parametrize(
    "avant, apres",
    [
        (f"from {ANCIEN}.reseau import x", f"from {NEUF}.reseau import x"),
        (f"python -m {ANCIEN}_api", f"python -m {NEUF}_api"),
        (f"{ANCIEN.upper()}_JETON", f"{NEUF.upper()}_JETON"),
        (f"reseaux/smib.{ANCIEN}.json", f"reseaux/smib.{NEUF}.json"),
        ("X-G2egrid-Jeton", "X-G2egrid-Jeton"),
        ("exporte par G2Egrid (f0 = 50 Hz)", "exporte par G2Egrid (f0 = 50 Hz)"),
        # ce qui ne designe pas l'outil garde son nom
        ("fidele a symSG de G2Elin V1", "fidele a symSG de G2Elin V1"),
        ("matlab/G2Elin_V2/core", "matlab/G2Elin_V2/core"),
        ("Portage Python de G2ELIN.", "Portage Python de G2ELIN."),
        ("le depot FKELADA/G2ELin", "le depot FKELADA/G2ELin"),
        (
            "git clone https://github.com/jbuire/g2elin_v2.git",
            "git clone https://github.com/jbuire/g2elin_v2.git",
        ),
        ("cd g2elin_v2/g2egrid", "cd g2elin_v2/g2egrid"),
        # le dossier du projet : python/ en developpement, g2egrid/ une fois publie
        ("cd g2egrid\nuv sync", "cd g2egrid\nuv sync"),
        (
            f"dans ``g2egrid/reseaux/<nom>.{ANCIEN}.json``",
            f"dans ``g2egrid/reseaux/<nom>.{NEUF}.json``",
        ),
        (
            "depuis le dossier `g2egrid/` du depot",
            "depuis le dossier `g2egrid/` du depot",
        ),
        ("uv run python -m pytest", "uv run python -m pytest"),
    ],
)
def test_renommage(pub, avant, apres):
    assert pub.renommer(avant) == apres


def test_les_deux_marques_de_format_se_lisent():
    """Un reseau enregistre par la version de developpement s'ouvre dans la
    version publiee, et inversement."""
    from g2egrid.reseau import zones
    from g2egrid.reseau.catalogue import charger_reseau

    d = zones.reseau_vers_json(charger_reseau("smib"))
    for marque in (ANCIEN, NEUF):
        r = zones.reseau_depuis_json(json.loads(json.dumps({**d, "format": marque})))
        assert len(r.transport.Y_DER) == 3
    with pytest.raises(ValueError, match="inconnu"):
        zones.reseau_depuis_json({**d, "format": "autre"})


def test_la_version_publiee_est_a_jour(pub, tmp_path):
    if not pub.SORTIE.is_dir():
        pytest.skip("g2egrid/ n'a pas encore ete genere")
    pub.construire(tmp_path / NEUF)
    ecarts = pub._differences(tmp_path / NEUF, pub.SORTIE)
    assert not ecarts, (
        f"{len(ecarts)} ecart(s) entre python/ et g2egrid/ (ex. {ecarts[:5]}) : "
        "relancer  uv run python developpement/publier_g2egrid.py"
    )


def test_la_version_publiee_ne_contient_que_l_outil(pub, tmp_path):
    racine = tmp_path / NEUF
    bilan = pub.construire(racine)
    presents = {p.name for p in racine.iterdir()}
    assert {NEUF, f"{NEUF}_api", "tests", "reseaux", "docs", "LICENSE", "README.md",
            "CITATION.cff", "pyproject.toml"} <= presents  # fmt: skip
    assert not {"etudes", "developpement", ANCIEN, f"{ANCIEN}_api"} & presents
    assert "GNU GENERAL PUBLIC LICENSE" in (racine / "LICENSE").read_text()
    projet = (racine / "pyproject.toml").read_text()
    assert f'name = "{NEUF}"' in projet and 'license = "GPL-3.0-or-later"' in projet
    # plus aucune mention de l'ancien nom de paquet dans le code publie
    for f in (racine / NEUF).rglob("*.py"):
        assert f"import {ANCIEN}" not in f.read_text(encoding="utf-8"), f
    assert bilan["fichiers"] > 300
