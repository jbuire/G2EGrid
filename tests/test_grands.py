"""La campagne des grands reseaux (``etudes/campagne/grands.py``) : tirages
cibles reproductibles, convertisseurs de distribution gardes ou convertis,
tirages publies inchanges, et statistiques de ``analyse_grands.py``."""

from __future__ import annotations

import gzip
import json

import numpy as np
import pytest
from _etudes import dossier_etude

from g2egrid.analyse.classement import CATEGORIES
from g2egrid.reseau.catalogue import charger
from g2egrid.reseau.tables import assembler

ICI = dossier_etude("campagne")

import analyse_grands
import campagne
import grands


def test_les_grands_reseaux_sont_au_catalogue():
    from g2egrid.multi_niveau.multi_echelle import wscc9_cigre9
    from g2egrid.reseau.zones import aplatir

    for nom in grands.GRANDS:
        t = charger(nom)
        assert len(t.Y_DER) > 100, nom
    a = assembler(charger("wscc9_cigre9")).run_loadflow()
    b = assembler(aplatir(wscc9_cigre9())).run_loadflow()
    assert a.converged and np.abs(a.V - b.V).max() < 1e-12
    assert len(charger("ieee118").T_SMs) == 54


def test_tirages_reproductibles_et_cibles_a_tour_de_role():
    for s in grands.GRANDS:
        n = len(grands.CYCLES.get(s, CATEGORIES))
        tirages = [grands.tirer(i, s) for i in range(2 * n)]
        cycle = [p["cible"] for p in tirages[:n]]
        assert [p["cible"] for p in tirages[n:]] == cycle
        # chaque classe visee au moins une fois par cycle
        assert set(cycle) == set(CATEGORIES)
        assert grands.tirer(3, s) == tirages[3]
        assert tirages[0] != tirages[n]  # meme cible, graines differentes
    assert [grands.cible_de(i, "wscc9_cigre9") for i in range(6)] == list(CATEGORIES)
    # IEEE 118 vise deux fois la frequence et trois fois la resonance
    cycle = [grands.cible_de(i, "ieee118") for i in range(9)]
    assert cycle.count("frequence") == 2 and cycle.count("resonance") == 3
    p = grands.tirer(0, "wscc9_cigre9")
    assert len(p["convertisseurs"]) == 27 and len(p["unites"]) == 3


def test_recette_propre_a_un_reseau():
    """``par_reseau`` remplace, cle par cle, ce qui differe d'un reseau a
    l'autre : un dictionnaire se fond dans celui de la recette commune."""
    commun = grands.poussee("tension", "wscc9_cigre9")
    propre = grands.poussee("tension", "ieee118")
    assert "par_reseau" not in commun and "par_reseau" not in propre
    assert commun["technologies"]["GFL"] == 1.0 and propre["technologies"]["SM"] == 1.0
    # le H de la recette commune reste, le Ka et le K_PSS sont remplaces
    assert propre["SM"]["H"] == commun["SM"]["H"]
    assert propre["SM"]["K_PSS"] == 0.0 and commun["SM"]["K_PSS"] != 0.0
    # la cible tension d'ieee118 : toutes les unites sont des machines
    p = grands.tirer(2, "ieee118")
    assert p["cible"] == "tension"
    assert {u["tech"] for u in p["unites"].values()} == {"SM"}
    assert all(40.0 <= u["Ka"] <= 150.0 for u in p["unites"].values())


def test_construire_garde_et_convertit_les_convertisseurs_de_distribution():
    p = grands.tirer(1, "wscc9_cigre9")  # vise la frequence : des GFM
    t = campagne.construire(p)
    gfm = {b for b, u in p["convertisseurs"].items() if u["tech"] == "GFM"}
    gfm |= {b for b, u in p["unites"].items() if u["tech"] == "GFM"}
    gfl = {b for b, u in p["convertisseurs"].items() if u["tech"] == "GFL"}
    gfl |= {b for b, u in p["unites"].items() if u["tech"] == "GFL"}
    assert gfm and {str(int(r[1])) for r in t.T_GFM1s} == gfm
    assert {str(int(r[1])) for r in t.T_GFL1s} == gfl
    b = int(next(iter(p["convertisseurs"])))
    assert t.PARAMS_PAR_BUS[b]["tau_CL"] == p["convertisseurs"][str(b)]["tau_CL"]
    assert assembler(t).run_loadflow().converged


def test_un_tirage_publie_se_recalcule_a_l_identique():
    """``construire`` garde desormais les convertisseurs d'origine et
    ``evaluer`` accepte des parametres : aucun reseau publie n'en a, et le
    tirage 0 de WSCC9 redonne sa partie reelle publiee."""
    publie = None
    with gzip.open(ICI / "tirages_nuit.jsonl.gz", "rt") as f:
        for ligne in f:
            r = json.loads(ligne)
            if r["systeme"] == "wscc9" and r["i"] == 0:
                publie = r
                break
    assert publie is not None
    r = campagne.evaluer(0, "aucun", "wscc9")
    assert r["max_re"] == pytest.approx(publie["max_re"], rel=1e-7, abs=1e-9)
    assert [m["categorie"] for m in r["modes"]] == [
        m["categorie"] for m in publie["modes"]
    ]


def _tirage(systeme, i, cible, modes=(), statut="ok"):
    return {
        "systeme": systeme,
        "i": i,
        "statut": statut,
        "params": {"cible": cible},
        "max_re": max((m[0] for m in modes), default=-1.0),
        "modes": [
            {"re": re, "f": 1.0, "instable": re > 0, "categorie": c, "groupes": {}}
            for re, c in modes
        ],
    }


def test_statistiques_par_tirage_et_equilibre():
    tirages = [
        # deux modes instables : le tirage compte pour le plus instable
        _tirage("x", 0, "angle", [(0.5, "angle"), (2.0, "tension")]),
        _tirage("x", 1, "frequence", [(0.1, "frequence")]),
        _tirage("x", 2, "tension", [(-0.1, "tension")]),
        _tirage("x", 3, "resonance", statut="point_de_fonctionnement"),
    ]
    r = analyse_grands.resumer(tirages)["x"]
    assert r["instables"] == 2 and r["taux_instables"] == pytest.approx(0.5)
    assert r["classe_dominante"]["tension"] == 1 and r["classe_dominante"]["angle"] == 0
    assert r["classes_presentes"]["angle"] == 1 and r["modes"]["tension"] == 1
    assert r["cible_contre_classe"]["tension"] == {"stable": 1}
    assert r["cible_contre_classe"]["resonance"] == {"point_de_fonctionnement": 1}
    assert analyse_grands.equitabilite(dict.fromkeys(CATEGORIES, 5)) == pytest.approx(1)
    assert analyse_grands.equitabilite({"angle": 9}) == 0.0
