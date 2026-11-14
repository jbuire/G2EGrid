"""Confort de l'atelier : lieu des racines parcouru, equations, lanceur."""

from __future__ import annotations

import pathlib
import socket

import numpy as np
import pytest
import sympy as sp

import g2egrid.catalogue as cat
from g2egrid.equations import classes, documenter, equations
from g2egrid_api import analyse
from g2egrid_api.__main__ import _repond, main

PAGE = pathlib.Path(__file__).resolve().parents[1] / "docs" / "equations_modeles.md"


# ------------------------------------------------------- lieu des racines
def test_le_balayage_dit_qui_participe():
    d = analyse.balayage(
        cat.charger("wscc9"),
        {"genre": "famille", "famille": "SM_PARAMS", "cle": "H"},
        [3.0, 6.0],
        f_max=5.0,
        n_modes=4,
    )
    for pas in d["pas"]:
        for m in pas["modes"]:
            parts = [p for _, p in m["part"]]
            assert parts == sorted(parts, reverse=True)
            assert all(p >= 0.05 for p in parts) and sum(parts) <= 1 + 1e-9
            # sur wscc9, sous 5 Hz, ce sont des modes de machines
            assert m["part"][0][0].startswith("SM")
            assert m["etat"].split(".")[0] in {n for n, _ in m["part"]}


# ------------------------------------------------------------- equations
@pytest.mark.parametrize("nom", sorted(classes()))
def test_chaque_modele_se_documente(nom):
    c = classes()[nom]
    d = equations(nom)
    assert [x["variable"] for x in d["f"]] == list(c.states)
    assert [x["variable"] for x in d["g"]] == list(c.algebs)
    assert [x["variable"] for x in d["h"]] == list(c.output_names())


@pytest.mark.parametrize("nom", ["Line", "SynchronousMachine", "GFM1", "GFL2"])
def test_le_texte_affiche_est_ce_qui_se_calcule(nom):
    """Relu par SymPy et evalue, le texte affiche redonne f et h du modele
    COMPILE, en un point quelconque."""
    c = classes()[nom]
    cm = c.compiled()
    spec = c.spec()
    rng = np.random.default_rng(0)
    x, z = rng.uniform(0.2, 0.9, cm.n_x), rng.uniform(0.2, 0.9, cm.n_z)
    u, p = rng.uniform(0.2, 0.9, cm.n_u), rng.uniform(0.2, 0.9, cm.n_p)
    valeurs = dict(zip(spec.states, x)) | dict(zip(spec.algebs, z))
    valeurs |= dict(zip(spec.inputs, u)) | dict(zip(spec.params, p))
    local = {n: sp.Symbol(n, real=True) for n in valeurs}
    d = equations(nom)
    for cle, attendu in (("f", cm.f(x, z, u, p)), ("h", cm.h(x, z, u, p))):
        lu = [
            float(sp.sympify(e["texte"], locals=local).subs(
                {local[n]: v for n, v in valeurs.items()}
            ))
            for e in d[cle]
        ]  # fmt: skip
        np.testing.assert_allclose(lu, attendu, rtol=1e-9, atol=1e-12)


def test_la_page_de_reference_est_a_jour():
    assert PAGE.read_text(encoding="utf-8") == documenter(), (
        "docs/equations_modeles.md ne decrit plus les modeles : "
        "uv run python -m g2egrid.equations > docs/equations_modeles.md"
    )


def test_api_equations():
    assert "SynchronousMachine" in analyse.equations_modele()["modeles"]
    d = analyse.equations_modele("IB")
    assert d["classe"] == "IB" and d["f"]
    with pytest.raises(ValueError, match="inconnu"):
        analyse.equations_modele("Soucoupe")


# ---------------------------------------------------------------- lanceur
def test_le_lanceur_ne_relance_pas_un_atelier_deja_ouvert():
    """Port deja pris : on n'essaie pas de demarrer un second serveur."""
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        s.listen()
        port = s.getsockname()[1]
        assert _repond("127.0.0.1", port)
        assert main(["--sans-navigateur", "--port", str(port)]) == 0
    assert not _repond("127.0.0.1", port)
