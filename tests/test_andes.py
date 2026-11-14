"""Validation croisee contre ANDES sur le reseau a deux zones de Kundur.

ANDES ne s'installe qu'avec le groupe ``validation`` (``uv sync --group
validation``) : sans lui, ces tests sont sautes.

ANDES n'est PAS importe a la collecte (pas de ``pytest.importorskip``) : son
import modifie ``sympy.Piecewise``, et les modeles compiles par les autres
fichiers de tests en restaient casses. Il ne l'est qu'a l'usage, par
:func:`g2egrid.andes_croise.importer_andes`, qui rend a SymPy sa classe.
"""

from __future__ import annotations

import importlib.util

import pytest

from g2egrid.andes_croise import comparer, kundur_deux_zones

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("andes") is None,
    reason="ANDES absent : uv sync --group validation",
)


@pytest.fixture(scope="module")
def ecarts():
    return comparer(kundur_deux_zones())


def test_load_flow_identique(ecarts):
    assert ecarts["ecart_V"] < 1e-5
    assert ecarts["ecart_angle"] < 1e-5


def test_trois_modes_electromecaniques(ecarts):
    """Deux modes locaux (un par zone) et le mode inter-zones, des deux cotes."""
    assert ecarts["n_andes"] == ecarts["n_g2egrid"] == 3
    f = [p["andes"][0] for p in ecarts["modes"]]
    assert 0.4 < f[0] < 0.7 < 1.0 < f[1] < f[2] < 1.3


def test_memes_frequences_memes_amortissements(ecarts):
    """Frequences a 0,5 % ; amortissement a 0,2 point pour les modes locaux,
    a 1 point pour l'inter-zones (0,35 mesure : le mode le plus sensible a la
    dynamique du reseau et du stator, que seul ce paquet represente)."""
    inter, *locaux = ecarts["modes"]
    for p in ecarts["modes"]:
        assert p["g2egrid"][0] == pytest.approx(p["andes"][0], rel=5e-3)
    for p in locaux:
        assert abs(p["g2egrid"][1] - p["andes"][1]) < 0.2
    assert abs(inter["g2egrid"][1] - inter["andes"][1]) < 1.0


def test_import_d_andes_ne_casse_pas_les_modeles_du_paquet():
    """Temoin du nettoyage de ``sympy.Piecewise`` (revue de FR-007).

    Dans un processus neuf, ANDES est importe AVANT toute compilation, puis
    un reseau a GFL (dont les equations portent un ``Piecewise``) est
    assemble et linearise. Sans :func:`importer_andes`, la compilation
    echouait sur ``NameError: name '__trues' is not defined``. Un processus
    a part, parce que le cache de compilation de ce processus-ci a pu etre
    rempli avant.
    """
    import subprocess
    import sys

    code = (
        "from g2egrid.andes_croise import importer_andes\n"
        "importer_andes()\n"
        "from g2egrid.catalogue import charger\n"
        "from g2egrid.tables import assembler\n"
        "net = assembler(charger('smib_gfl'))\n"
        "net.solve_equilibrium()\n"
        "print(net.build_state_space().A.shape)\n"
    )
    r = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )
    assert r.returncode == 0, r.stderr[-2000:]
