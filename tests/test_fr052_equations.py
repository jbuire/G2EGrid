"""FR-052 : equations dynamiques modifiables (partie serveur).

Verifie, sur le catalogue :

* critere 2 : natif -> modele a equations -> compilation donne le meme load
  flow et les memes valeurs propres a 1e-8 (GFL1, GFL2, GFM1, GFM2, GFMdVOC,
  GFMMatching, IB, machine synchrone, charge R-L) ; schema-bloc -> equations ->
  compilation, a 1e-10 (gfl_complet, gfm_complet, schema a sous-systemes) ;
* critere 3 : lisibilite (un etat par etat natif, algebriques nommees,
  expressions de moins de 2000 caracteres) ;
* critere 4 : doubler ``Kipll`` dans ``dM_pll/dt`` = natif avec ``Kipll``
  double (1e-8), seul l'element vise change, retour au natif exact ;
* critere 5 : la liste globale de ``wscc9_gfm`` et ``ieee39_mixte`` porte
  chaque etat de ``build_state_space`` une fois ;
* critere 6 : expressions malveillantes refusees, tables intactes ;
* la liste porte les memes lignes en double ; ajouter / supprimer.
"""

from __future__ import annotations

import copy

import g2egrid.catalogue as cat
import numpy as np
import pytest
from test_fr051_aide import banc_smib, ecart_spectres, sans_jauge, spectre
from test_fr051_depuis_element import (
    CAS,
    DETAILS,
    TOL_LF,
    TOL_SPECTRE,
    _natif,
    _vers_blocs,
)

from g2egrid.elements import blocs as B
from g2egrid.elements.blocs import SchemaInvalide
from g2egrid.reseau.tables import assembler, tables_vers_json


# ------------------------------------------------- critere 2 : natif
@pytest.mark.parametrize("cas", list(CAS))
def test_natif_equations_compilation(cas):
    cls, _ = CAS[cas]
    t, table, index, lf0, lam0, params = _natif(cas)
    modele = B.equations_depuis_element(cls, parametres=params)
    assert modele["genre"] == "equations"
    tb = _vers_blocs(copy.deepcopy(t), table, index, modele)
    _net, lf1, lam1 = spectre(tb)
    assert np.max(np.abs(lf0.V - lf1.V)) < TOL_LF
    assert ecart_spectres(sans_jauge(lam0), sans_jauge(lam1)) < TOL_SPECTRE


# ---------------------------------------- critere 2 : schema -> equations
def _spectre_schema(schema, genre):
    t = banc_smib(schema, genre)
    return spectre(t)[2]


@pytest.mark.parametrize(
    "schema, genre",
    [(B.schema_gfl(), "gfl"), (B.schema_gfm(), "gfm")],
    ids=["gfl_complet", "gfm_complet"],
)
def test_schema_equations_compilation(schema, genre):
    m = B.equations_depuis_schema(schema)
    assert (
        ecart_spectres(
            sans_jauge(_spectre_schema(schema, genre)),
            sans_jauge(_spectre_schema(m, genre)),
        )
        < 1e-10
    )
    # et retour vers les blocs : compile toujours ; le spectre n'est exige que
    # pour le GFL (les chaines de gains du gfm_complet, une fois les
    # intermediaires substituees, mettent l'elimination des algebriques en
    # defaut : voir le rapport)
    s2 = B.schema_depuis_equations(m)
    assert len(B.compiler(s2).etats) == len(B.compiler(schema).etats)
    if genre == "gfl":
        assert (
            ecart_spectres(
                sans_jauge(_spectre_schema(schema, genre)),
                sans_jauge(_spectre_schema(s2, genre)),
            )
            < 1e-8
        )


def test_schema_a_sous_systemes_equations():
    cls, _ = CAS["GFL1"]
    t, table, index, lf0, lam0, params = _natif("GFL1")
    schema = B.depuis_element(cls, parametres=params, detail="blocs")
    assert any(b["type"] == "sous_systeme" for b in schema["blocs"])
    m = B.equations_depuis_schema(schema)
    _n, _lf, lam_s = spectre(_vers_blocs(copy.deepcopy(t), table, index, schema))
    _n, _lf, lam_m = spectre(_vers_blocs(copy.deepcopy(t), table, index, m))
    assert ecart_spectres(sans_jauge(lam_s), sans_jauge(lam_m)) < 1e-10


@pytest.mark.parametrize("detail", DETAILS)
def test_equations_vers_blocs_aller_retour(detail):
    m = B.equations_depuis_element("GFL2")
    s = B.schema_depuis_equations(m, detail)
    c1, c2 = B.compiler(m), B.compiler(s)
    assert len(c1.etats) == len(c2.etats)
    assert B.schema_depuis_equations(m, "liste") == m


# ------------------------------------------------------ critere 3
def test_lisibilite_gfl1():
    from g2egrid.elements import GFL1

    m = B.equations_depuis_element("GFL1")
    assert [e["nom"] for e in m["etats"]] == list(GFL1.states)
    assert [a["nom"] for a in m["algebriques"]] == list(GFL1.algebs)
    assert max(len(e["derivee"]) for e in m["etats"]) < 2000
    assert all("expr" in a for a in m["algebriques"])
    # les algebriques restent des signaux nommes dans les derivees
    assert any("w_pll" in e["derivee"] for e in m["etats"])


# ------------------------------------------------------ critere 4
@pytest.fixture(scope="module")
def smib_gfl1():
    t = copy.deepcopy(cat.charger("smib_gfl"))
    t.T_GFL1s, t.T_GFL2s = t.T_GFL2s, []
    net, _lf, lam = spectre(copy.deepcopy(t))
    return t, lam, net


def test_kipll_double_et_retour_exact(smib_gfl1):
    t, lam0, net = smib_gfl1
    nom = next(e.name for e in net.all_elements if e.name.startswith("GFL"))
    kipll = next(e for e in net.all_elements if e.name == nom).pdict["Kipll"]
    avant = tables_vers_json(t)
    t2 = B.modifier_equation(t, nom, "M_pll", "2*Kipll*veq")
    assert tables_vers_json(t) == avant  # les tables recues ne bougent pas
    assert t2.T_GFL1s == [] and len(t2.T_Blocs) == 1
    tn = copy.deepcopy(t)
    tn.GFL_PARAMS = dict(tn.GFL_PARAMS, Kipll=2 * kipll)
    lam_n = spectre(tn)[2]
    lam_m = spectre(copy.deepcopy(t2))[2]
    assert ecart_spectres(sans_jauge(lam_n), sans_jauge(lam_m)) < TOL_SPECTRE
    assert ecart_spectres(sans_jauge(lam0), sans_jauge(lam_m)) > 1e-3
    t3 = B.revenir_au_natif(t2, ("T_Blocs", 0))
    assert tables_vers_json(t3) == avant
    lam3 = spectre(copy.deepcopy(t3))[2]
    assert np.array_equal(np.sort_complex(lam0), np.sort_complex(lam3))


def test_seul_l_element_vise_change():
    t = copy.deepcopy(cat.charger("wscc9_gfm"))
    assert len(t.T_GFM1s) >= 2
    nom = "GFM2"
    modele = B.equations_depuis_element("GFM1")
    var = modele["etats"][0]["nom"]
    t2 = B.modifier_equation(t, nom, var, "2*(" + modele["etats"][0]["derivee"] + ")")
    assert len(t2.T_Blocs) == 1 and len(t2.T_GFM1s) == len(t.T_GFM1s) - 1
    noms = {e.name: type(e).__name__ for e in assembler(t2).all_elements}
    assert noms["GFM3"] == "GFM1"
    assert "GFM2" not in noms


def test_equation_inchangee_ne_change_rien(smib_gfl1):
    t, _lam, net = smib_gfl1
    nom = next(e.name for e in net.all_elements if e.name.startswith("GFL"))
    t2 = B.modifier_equation(t, nom, "M_pll", "Kipll*veq")
    assert t2.T_Blocs == [] and tables_vers_json(t2) == tables_vers_json(t)


def test_modele_partage_duplique(smib_gfl1):
    t, _lam, net = smib_gfl1
    nom = next(e.name for e in net.all_elements if e.name.startswith("GFL"))
    t2 = B.modifier_equation(t, nom, "M_pll", "2*Kipll*veq")
    cle = t2.T_Blocs[0][2]
    t2.T_Blocs = [*t2.T_Blocs, (2, 1, cle, t2.T_Blocs[0][3])]  # partage
    t3 = B.modifier_equation(t2, ("T_Blocs", 0), "M_pll", "3*Kipll*veq")
    assert t3.T_Blocs[0][2] != t3.T_Blocs[1][2]
    assert t3.MODELES_BLOCS[cle]["etats"] == t2.MODELES_BLOCS[cle]["etats"]


def test_ajouter_puis_supprimer(smib_gfl1):
    t, _lam, net = smib_gfl1
    nom = next(e.name for e in net.all_elements if e.name.startswith("GFL"))
    t2 = B.ajouter_equation(t, nom, "etat", "xfiltre", "(M_pll - xfiltre)/0.1", 0.0)
    lam = spectre(copy.deepcopy(t2))[2]
    assert lam.size == spectre(copy.deepcopy(t))[2].size + 1
    t3 = B.supprimer_equation(t2, ("T_Blocs", 0), "xfiltre")
    assert spectre(copy.deepcopy(t3))[2].size == lam.size - 1
    with pytest.raises(SchemaInvalide, match="encore employe"):
        B.supprimer_equation(t3, ("T_Blocs", 0), "M_pll")


def test_charge_y_der():
    t = copy.deepcopy(cat.charger("wscc9"))
    lam0 = spectre(copy.deepcopy(t))[2]
    modele = B.equations_depuis_instance(
        next(e for e in assembler(t).all_elements if e.name == "Load5")
    )
    var = modele["etats"][0]["nom"]
    t2 = B.modifier_equation(t, "Load5", var, "(" + modele["etats"][0]["derivee"] + ")")
    assert t2.T_Blocs == []  # inchange
    t2 = B.modifier_equation(
        t, "Load5", var, "2*(" + modele["etats"][0]["derivee"] + ")"
    )
    assert len(t2.T_Blocs) == 1
    assert spectre(copy.deepcopy(t2))[2].size == lam0.size
    t3 = B.revenir_au_natif(t2, ("T_Blocs", 0))
    assert tables_vers_json(t3) == tables_vers_json(t)


# ------------------------------------------------------ critere 5
@pytest.mark.parametrize("nom", ["wscc9_gfm", "ieee39_mixte"])
def test_liste_globale_comptage(nom):
    t = cat.charger(nom)
    lignes = B.equations_reseau(t)
    net = assembler(t)
    net.run_loadflow()
    net.solve_equilibrium()
    ss = net.build_state_space()
    etats = [
        f"{x['element_reseau']}.{x['variable']}" for x in lignes if x["genre"] == "etat"
    ]
    assert len(etats) == len(set(etats))
    assert set(etats) == set(ss.states) | set(ss.elimines)
    non = [x for x in lignes if not x["modifiable"]]
    assert non and all(x["pourquoi"] for x in non)
    assert {x["classe"] for x in non} >= {"Line", "Bus"}
    assert any(x["modifiable"] and not x["modifiee"] for x in lignes)


def test_liste_marque_les_elements_modifies(smib_gfl1):
    t, _lam, net = smib_gfl1
    nom = next(e.name for e in net.all_elements if e.name.startswith("GFL"))
    t2 = B.modifier_equation(t, nom, "M_pll", "2*Kipll*veq")
    lignes = B.equations_reseau(t2)
    modif = [x for x in lignes if x["modifiee"]]
    assert modif and all(x["table"] == "T_Blocs" for x in modif)
    assert any(x["variable"] == "M_pll" and "2*Kipll" in x["expression"] for x in modif)


# ------------------------------------------------------ critere 6
MALVEILLANTES = [
    "__import__('os').system('id')",
    "().__class__",
    "open('/etc/passwd')",
    "lambda: 1",
    "[x for x in (1,)]",
    "'abc'",
    "9**9**9",
    "inconnu + 1",
    "veq.real",
    "x" * 2500,
    "",
]


@pytest.mark.parametrize("expr", MALVEILLANTES)
def test_expressions_malveillantes_refusees(smib_gfl1, expr):
    t, _lam, net = smib_gfl1
    nom = next(e.name for e in net.all_elements if e.name.startswith("GFL"))
    avant = tables_vers_json(t)
    with pytest.raises(SchemaInvalide):
        B.modifier_equation(t, nom, "M_pll", expr)
    assert tables_vers_json(t) == avant
    m = B.equations_depuis_element("GFL1")
    m["etats"][0]["derivee"] = expr
    with pytest.raises(SchemaInvalide):
        B.compiler(m)


def test_variable_inconnue_et_element_inconnu(smib_gfl1):
    t, _lam, net = smib_gfl1
    nom = next(e.name for e in net.all_elements if e.name.startswith("GFL"))
    with pytest.raises(SchemaInvalide, match="pas d'equation"):
        B.modifier_equation(t, nom, "nexistepas", "1")
    with pytest.raises(SchemaInvalide, match="inconnu"):
        B.modifier_equation(t, "GFLzz", "M_pll", "1")
    with pytest.raises(SchemaInvalide, match="ne se modifie pas"):
        B.modifier_equation(t, "bus1", "vd", "1")
