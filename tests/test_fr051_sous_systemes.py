"""FR-051, critere 3 : sous-systemes, aplatissement, regrouper / degrouper.

Tests ecrits A L'AVEUGLE depuis le seul cahier
``.claude/equipe/cahiers/FR-051_editeur_blocs.md`` (contrat : type
``sous_systeme`` de champs ``nom`` et ``schema``, blocs ``port_entree`` /
``port_sortie`` d'``index`` et ``nom``, ``aplatir``, ``regrouper``,
``degrouper``).

Ce qui est verifie, avec des references independantes :

- la valeur d'un sous-systeme calculee a la main (ports lus dans l'ordre des
  ``index``, pas dans l'ordre de la liste) ;
- ``compiler(s)`` et ``compiler(aplatir(s))`` donnent les memes equations
  (empreinte numerique, tolerance 1e-11, voir ``test_fr051_aide``) ; les blocs
  interieurs sont renommes ``<id_ss>__<id>`` recursivement, les ports retires ;
- ``regrouper`` puis ``degrouper`` rend le schema de depart (aux positions
  pres) ; un port par signal entrant, meme s'il alimente plusieurs blocs
  interieurs ; un port par signal interieur utilise dehors ;
- regrouper des blocs quelconques de ``gfl_complet`` / ``gfm_complet`` laisse
  les equations identiques AU RENOMMAGE DOCUMENTE PRES (``x_<id>`` devient
  ``x_<ss>__<id>``, idem ``y_`` et ``mode_``), et le spectre identique sur un
  banc SMIB ;
- imbrication sur trois niveaux, et refus au-dela de la profondeur 8.

TOLERANCE DU SPECTRE. Le cahier dit 1e-12 relatif. Le regroupement permute
l'ordre des etats : la comparaison est donc entre ``eig(A)`` et
``eig(P A P^T)``. Mesure sur le banc SMIB (permutations aleatoires de la
matrice A du GFM, 15 etats, |A| max 5.7e5) : ecart relatif maximal 5.1e-13,
5 a 20 fois plus que sur le GFL (1.3e-14). 1e-12 serait donc dans le bruit de
``numpy.linalg.eigvals`` ; le test emploie 1e-11 (vingt fois l'ecart mesure).
Les EQUATIONS, elles, sont comparees a 1e-11 sur des evaluations
(arrondi ~1e-15) : c'est la voie stricte.

Lectures du cahier retenues :

- le sous-systeme cree par ``regrouper`` est le seul bloc de type
  ``sous_systeme`` du resultat ; ``id`` s'il est donne ;
- sa position ``x, y`` est la MOYENNE des positions des blocs deplaces
  (confondue avec le centre de leur boite englobante quand les blocs ont la
  meme taille, cas des deux blocs testes) ;
- ``degrouper`` remet les blocs interieurs au niveau superieur sous leur
  identifiant interieur ;
- profondeur maximale 8 : 7 niveaux imbriques passent, 10 sont refuses ;
- ``consigne`` et ``sortie`` ne sont pas deplacees dans les regroupements
  aleatoires (leur ordre fixe celui des consignes et des sorties du modele).
"""

from __future__ import annotations

import re

import numpy as np
import pytest
import sympy as sp
from g2egrid.tables import Tables
from test_fr051_aide import (
    banc_smib,
    copie,
    ecart_spectres,
    empreinte,
    memes_empreintes,
    sans_jauge,
    schema_sans_position,
    spectre,
    valeurs_statiques,
)

from g2egrid.elements import blocs as B
from g2egrid.elements.blocs import SchemaInvalide, compiler

TOL_EQ = 1e-11
TOL_SPECTRE = 1e-11  # voir la docstring : le bruit de eig atteint 5e-13


# ------------------------------------------------------------ constructions
def ss_bloc(bid: str, nom: str, blocs: list, liens: list, x=300, y=100) -> dict:
    return {
        "id": bid,
        "type": "sous_systeme",
        "nom": nom,
        "x": x,
        "y": y,
        "schema": {"blocs": blocs, "liens": liens},
    }


def _consignes() -> list[dict]:
    return [
        {"id": "pr", "type": "consigne", "nom": "P_ref", "x": 20, "y": 20},
        {"id": "qr", "type": "consigne", "nom": "Q_ref", "x": 20, "y": 120},
    ]


def schema_diff() -> dict:
    """Un sous-systeme a deux entrees (index 0 : a, 1 : b) et deux sorties
    (index 0 : ``3 (a - b)``, index 1 : ``a - b``). Les ports sont declares
    DANS LE DESORDRE dans la liste : seul ``index`` compte."""
    interieur = [
        {"id": "pb", "type": "port_entree", "index": 1, "nom": "b", "x": 20, "y": 120},
        {"id": "pa", "type": "port_entree", "index": 0, "nom": "a", "x": 20, "y": 20},
        {"id": "s", "type": "somme", "signes": "+-", "x": 170, "y": 70},
        {"id": "g", "type": "gain", "K": 3, "x": 320, "y": 20},
        {"id": "p1", "type": "port_sortie", "index": 1, "nom": "d", "x": 470, "y": 120},
        {"id": "p0", "type": "port_sortie", "index": 0, "nom": "r", "x": 470, "y": 20},
    ]
    liens_int = [
        {"de": "pa", "vers": "s", "entree": 0},
        {"de": "pb", "vers": "s", "entree": 1},
        {"de": "s", "vers": "g"},
        {"de": "g", "vers": "p0"},
        {"de": "s", "vers": "p1"},
    ]
    return {
        "nom": "avec_ss",
        "parametres": {},
        "blocs": _consignes()
        + [
            {"id": "a", "type": "constante", "valeur": 5.0, "x": 20, "y": 220},
            {"id": "b", "type": "constante", "valeur": 2.0, "x": 20, "y": 320},
            ss_bloc("ss", "diff3", interieur, liens_int, x=250, y=250),
            {"id": "oi", "type": "sortie", "signal": "id", "x": 500, "y": 220},
            {"id": "oq", "type": "sortie", "signal": "iq", "x": 500, "y": 320},
        ],
        "liens": [
            {"de": "a", "vers": "ss", "entree": 0},
            {"de": "b", "vers": "ss", "entree": 1},
            {"de": "ss", "vers": "oi", "sortie": 0},
            {"de": "ss", "vers": "oq", "sortie": 1},
        ],
    }


def schema_dynamique() -> dict:
    """``id`` = premier ordre de P_ref, ``iq`` = Q_ref ; P_ref, Q_ref et f
    sont poses a plat (point de depart des regroupements)."""
    return {
        "nom": "dyn",
        "parametres": {"tau": 0.05, "k": 2.0},
        "blocs": _consignes()
        + [
            {"id": "g1", "type": "gain", "K": "k", "x": 150, "y": 20},
            {
                "id": "f",
                "type": "premier_ordre",
                "K": 1,
                "tau": "tau",
                "x": 300,
                "y": 20,
            },
            {"id": "oi", "type": "sortie", "signal": "id", "x": 450, "y": 20},
            {"id": "oq", "type": "sortie", "signal": "iq", "x": 450, "y": 120},
        ],
        "liens": [
            {"de": "pr", "vers": "g1"},
            {"de": "g1", "vers": "f"},
            {"de": "f", "vers": "oi"},
            {"de": "qr", "vers": "oq"},
        ],
    }


def schema_multi() -> dict:
    """Une entree ``vd`` alimente DEUX gains ; leur somme part sur ``id``."""
    return {
        "nom": "multi",
        "parametres": {},
        "blocs": _consignes()
        + [
            {"id": "e", "type": "entree", "signal": "vd", "x": 20, "y": 220},
            {"id": "g1", "type": "gain", "K": 2, "x": 200, "y": 180},
            {"id": "g2", "type": "gain", "K": 3, "x": 200, "y": 280},
            {"id": "s", "type": "somme", "signes": "++", "x": 400, "y": 230},
            {"id": "oi", "type": "sortie", "signal": "id", "x": 600, "y": 230},
            {"id": "oq", "type": "sortie", "signal": "iq", "x": 600, "y": 120},
        ],
        "liens": [
            {"de": "e", "vers": "g1"},
            {"de": "e", "vers": "g2"},
            {"de": "g1", "vers": "s", "entree": 0},
            {"de": "g2", "vers": "s", "entree": 1},
            {"de": "s", "vers": "oi"},
            {"de": "qr", "vers": "oq"},
        ],
    }


def chaine(profondeur: int) -> dict:
    """``profondeur`` sous-systemes emboites ``n1 > n2 > ...`` ; au fond, un
    gain 2 entre un port d'entree et un port de sortie. ``id`` = 2 x 1,5."""
    interieur = {
        "blocs": [
            {
                "id": "pe",
                "type": "port_entree",
                "index": 0,
                "nom": "u",
                "x": 20,
                "y": 20,
            },
            {"id": "g", "type": "gain", "K": 2, "x": 170, "y": 20},
            {
                "id": "ps",
                "type": "port_sortie",
                "index": 0,
                "nom": "y",
                "x": 320,
                "y": 20,
            },
        ],
        "liens": [{"de": "pe", "vers": "g"}, {"de": "g", "vers": "ps"}],
    }
    for k in range(profondeur, 0, -1):
        bloc = {
            "id": f"n{k}",
            "type": "sous_systeme",
            "nom": f"niveau{k}",
            "x": 150,
            "y": 20,
            "schema": interieur,
        }
        interieur = {
            "blocs": [
                {
                    "id": "pe",
                    "type": "port_entree",
                    "index": 0,
                    "nom": "u",
                    "x": 20,
                    "y": 20,
                },
                bloc,
                {
                    "id": "ps",
                    "type": "port_sortie",
                    "index": 0,
                    "nom": "y",
                    "x": 320,
                    "y": 20,
                },
            ],
            "liens": [{"de": "pe", "vers": f"n{k}"}, {"de": f"n{k}", "vers": "ps"}],
        }
    # `interieur` contient maintenant le niveau 0 (ports compris) : on le deplie
    # d'un cran, le niveau superieur du modele n'ayant pas de ports
    n1 = next(b for b in interieur["blocs"] if b["type"] == "sous_systeme")
    return {
        "nom": "chaine",
        "parametres": {},
        "blocs": _consignes()
        + [
            {"id": "c", "type": "constante", "valeur": 1.5, "x": 20, "y": 220},
            n1,
            {"id": "z", "type": "constante", "valeur": 0.0, "x": 20, "y": 320},
            {"id": "oi", "type": "sortie", "signal": "id", "x": 500, "y": 220},
            {"id": "oq", "type": "sortie", "signal": "iq", "x": 500, "y": 320},
        ],
        "liens": [
            {"de": "c", "vers": "n1"},
            {"de": "n1", "vers": "oi"},
            {"de": "z", "vers": "oq"},
        ],
    }


def _bloc(schema: dict, bid: str) -> dict:
    return next(b for b in schema["blocs"] if b["id"] == bid)


# ================================================== valeur d'un sous-systeme
def test_ports_lus_dans_l_ordre_des_index():
    """a = 5, b = 2 : sortie 0 = 3 (5 - 2) = 9, sortie 1 = 5 - 2 = 3. Les
    ports etant declares dans le desordre, une lecture par ordre de liste
    donnerait -9 et -3."""
    assert valeurs_statiques(schema_diff()) == [
        pytest.approx(9.0, abs=1e-12),
        pytest.approx(3.0, abs=1e-12),
    ]


def test_parametres_visibles_dans_les_sous_systemes():
    """Pas de masque en v1 : un gain interieur lit un parametre du modele."""
    s = schema_diff()
    s["parametres"] = {"kk": 4.0}
    _bloc(_bloc(s, "ss")["schema"], "g")["K"] = "kk"
    assert valeurs_statiques(s)[0] == pytest.approx(12.0, abs=1e-12)


def test_sortie_admise_dans_un_sous_systeme():
    """Les blocs ``sortie`` (comme ``entree``, ``constante``, ``consigne``)
    sont admis a tous les niveaux : ``id`` ecrit depuis l'interieur."""
    s = {
        "nom": "sortie_interne",
        "parametres": {},
        "blocs": _consignes()
        + [
            ss_bloc(
                "ss",
                "emet",
                [
                    {"id": "c", "type": "constante", "valeur": 0.7, "x": 20, "y": 20},
                    {"id": "oi", "type": "sortie", "signal": "id", "x": 200, "y": 20},
                ],
                [{"de": "c", "vers": "oi"}],
            ),
            {"id": "z", "type": "constante", "valeur": 0.0, "x": 20, "y": 220},
            {"id": "oq", "type": "sortie", "signal": "iq", "x": 200, "y": 220},
        ],
        "liens": [{"de": "z", "vers": "oq"}],
    }
    assert valeurs_statiques(s) == [pytest.approx(0.7), pytest.approx(0.0)]


def test_unicite_sur_le_modele_aplati():
    """Deux ``sortie`` ``id`` a deux niveaux differents : une seule ecriture
    permise, sur le modele APLATI."""
    compiler(schema_diff())  # temoin : le schema sain compile
    s = schema_diff()
    ss = _bloc(s, "ss")["schema"]
    ss["blocs"].append(
        {"id": "o2", "type": "sortie", "signal": "id", "x": 400, "y": 300}
    )
    ss["liens"].append({"de": "s", "vers": "o2"})
    with pytest.raises(SchemaInvalide):
        compiler(s)


# ================================================================ erreurs
def test_les_ports_sont_refuses_hors_sous_systeme():
    compiler(schema_diff())  # temoin : les ports sont admis DANS un sous-systeme
    s = schema_dynamique()
    s["blocs"].append({"id": "p", "type": "port_entree", "index": 0, "nom": "u"})
    with pytest.raises(SchemaInvalide):
        compiler(s)
    s = schema_dynamique()
    s["blocs"].append({"id": "p", "type": "port_sortie", "index": 0, "nom": "u"})
    with pytest.raises(SchemaInvalide):
        compiler(s)


@pytest.mark.parametrize("niveau", ["haut", "interieur"])
def test_double_souligne_interdit_dans_un_identifiant(niveau):
    """Temoin : le schema sain compile ; un id avec ``__`` est refuse."""
    compiler(schema_diff())
    s = schema_diff()
    if niveau == "haut":
        _bloc(s, "a")["id"] = "a__x"
        s["liens"][0]["de"] = "a__x"
    else:
        ss = _bloc(s, "ss")["schema"]
        _bloc(ss, "g")["id"] = "g__x"
        for lien in ss["liens"]:
            for k in ("de", "vers"):
                if lien[k] == "g":
                    lien[k] = "g__x"
    with pytest.raises(SchemaInvalide):
        compiler(s)


def test_entree_de_sous_systeme_non_reliee_est_refusee():
    compiler(schema_diff())
    s = schema_diff()
    s["liens"] = [
        lien for lien in s["liens"] if lien["vers"] != "ss" or lien["entree"] != 1
    ]
    with pytest.raises(SchemaInvalide):
        compiler(s)


def test_profondeur_maximale():
    """Profondeur <= 8 (lecture : 7 niveaux emboites passent, 10 non)."""
    assert valeurs_statiques(chaine(7))[0] == pytest.approx(3.0, abs=1e-12)
    with pytest.raises(SchemaInvalide):
        compiler(chaine(10))


# ================================================================ aplatir
def test_aplatir_renomme_retire_les_ports_et_raccorde():
    s = schema_diff()
    plat = B.aplatir(copie(s))
    types = {b["type"] for b in plat["blocs"]}
    assert not types & {"sous_systeme", "port_entree", "port_sortie"}
    ids = {b["id"] for b in plat["blocs"]}
    assert {"ss__s", "ss__g"} <= ids
    assert not {"ss", "ss__pa", "ss__pb", "ss__p0", "ss__p1"} & ids
    assert {"pr", "qr", "a", "b", "oi", "oq"} <= ids
    # liens raccordes : a -> ss__s (entree 0), b -> ss__s (entree 1),
    # ss__g -> oi, ss__s -> oq
    liens = {
        (
            lien["de"],
            int(lien.get("sortie", 0)),
            lien["vers"],
            int(lien.get("entree", 0)),
        )
        for lien in plat["liens"]
    }
    assert ("a", 0, "ss__s", 0) in liens
    assert ("b", 0, "ss__s", 1) in liens
    assert ("ss__g", 0, "oi", 0) in liens
    assert ("ss__s", 0, "oq", 0) in liens
    assert ("ss__s", 0, "ss__g", 0) in liens


def test_aplatir_ne_modifie_pas_son_entree_et_est_idempotent():
    s = schema_diff()
    avant = copie(s)
    plat = B.aplatir(s)
    assert s == avant
    assert schema_sans_position(B.aplatir(copie(plat))) == schema_sans_position(plat)


@pytest.mark.parametrize(
    "fabrique",
    [schema_diff, lambda: chaine(3), lambda: chaine(1)],
    ids=["diff", "chaine3", "chaine1"],
)
def test_compiler_et_compiler_aplati_donnent_les_memes_equations(fabrique):
    s = fabrique()
    assert memes_empreintes(
        empreinte(compiler(copie(s))), empreinte(compiler(B.aplatir(copie(s)))), TOL_EQ
    )


def test_imbrication_sur_trois_niveaux_noms_recursifs():
    plat = B.aplatir(chaine(3))
    ids = {b["id"] for b in plat["blocs"]}
    assert "n1__n2__n3__g" in ids
    assert not any(i in ids for i in ("n1", "n1__n2", "n1__n2__n3"))
    assert valeurs_statiques(chaine(3))[0] == pytest.approx(3.0, abs=1e-12)
    c = compiler(chaine(3))
    assert "y_n1__n2__n3__g" in c.algebriques


def test_etat_a_trois_niveaux_dans_un_reseau():
    """Un premier ordre enfoui a trois niveaux garde son pole ``-1/tau`` dans
    le spectre du reseau, et son etat s'appelle ``x_n1__n2__n3__f`` (le nom
    aplati : ``x_`` puis l'identifiant recursif)."""
    tau = 0.05
    interieur = {
        "blocs": [
            {
                "id": "pe",
                "type": "port_entree",
                "index": 0,
                "nom": "u",
                "x": 20,
                "y": 20,
            },
            {
                "id": "f",
                "type": "premier_ordre",
                "K": 1,
                "tau": "tau",
                "x": 170,
                "y": 20,
            },
            {
                "id": "ps",
                "type": "port_sortie",
                "index": 0,
                "nom": "y",
                "x": 320,
                "y": 20,
            },
        ],
        "liens": [{"de": "pe", "vers": "f"}, {"de": "f", "vers": "ps"}],
    }
    for k in (3, 2, 1):
        bloc = ss_bloc(f"n{k}", f"niveau{k}", interieur["blocs"], interieur["liens"])
        interieur = {
            "blocs": [
                {
                    "id": "pe",
                    "type": "port_entree",
                    "index": 0,
                    "nom": "u",
                    "x": 20,
                    "y": 20,
                },
                bloc,
                {
                    "id": "ps",
                    "type": "port_sortie",
                    "index": 0,
                    "nom": "y",
                    "x": 320,
                    "y": 20,
                },
            ],
            "liens": [{"de": "pe", "vers": f"n{k}"}, {"de": f"n{k}", "vers": "ps"}],
        }
    n1 = next(b for b in interieur["blocs"] if b["type"] == "sous_systeme")
    s = {
        "nom": "enfoui",
        "parametres": {"tau": tau},
        "blocs": _consignes()
        + [
            n1,
            {"id": "oi", "type": "sortie", "signal": "id", "x": 400, "y": 20},
            {"id": "oq", "type": "sortie", "signal": "iq", "x": 400, "y": 120},
        ],
        "liens": [
            {"de": "pr", "vers": "n1"},
            {"de": "n1", "vers": "oi"},
            {"de": "qr", "vers": "oq"},
        ],
    }
    assert compiler(s).etats == ["x_n1__n2__n3__f"]
    t = Tables(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3, 20e3, 1.0, 0.0, 20e6, 5e6, 0.0, 0.0),
        ],
        T_Lines_pu=[(1, 1, 2, 1.0, 0.01, 0.05, 0.0, 20e3, 100e6)],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
        T_Shunts=[(1, 2, 0.01)],
        T_Blocs=[(1, 2, "src")],
        MODELES_BLOCS={"src": s},
    )
    _net, _lf, lam = spectre(t)
    assert np.min(np.abs(lam + 1 / tau)) < 1e-6 / tau


# ======================================================= regrouper / degrouper
def _ports(ss: dict, genre: str) -> list[dict]:
    return sorted(
        (b for b in ss["schema"]["blocs"] if b["type"] == genre),
        key=lambda b: b["index"],
    )


def _le_sous_systeme(schema: dict) -> dict:
    ss = [b for b in schema["blocs"] if b["type"] == "sous_systeme"]
    assert len(ss) == 1
    return ss[0]


def test_regrouper_cree_les_ports_et_raccorde():
    """Regrouper e (somme) et c (PI) de l'asservissement : entrent pr et m
    (deux signaux -> deux ports d'entree) ; sort c, utilise par m ET par op
    (un signal -> UN port de sortie)."""
    s = {
        "nom": "asservi",
        "parametres": {"Kp": 0.5, "Ki": 20.0, "tm": 0.02},
        "blocs": _consignes()
        + [
            {"id": "e", "type": "somme", "signes": "+-", "x": 100, "y": 100},
            {"id": "c", "type": "pi", "Kp": "Kp", "Ki": "Ki", "x": 250, "y": 100},
            {
                "id": "m",
                "type": "premier_ordre",
                "K": 1,
                "tau": "tm",
                "x": 400,
                "y": 100,
            },
            {"id": "op", "type": "sortie", "signal": "P", "x": 550, "y": 100},
            {"id": "oq", "type": "sortie", "signal": "Q", "x": 550, "y": 200},
        ],
        "liens": [
            {"de": "pr", "vers": "e", "entree": 0},
            {"de": "m", "vers": "e", "entree": 1},
            {"de": "e", "vers": "c"},
            {"de": "c", "vers": "m"},
            {"de": "c", "vers": "op"},
            {"de": "qr", "vers": "oq"},
        ],
    }
    avant = copie(s)
    r = B.regrouper(s, ["e", "c"], "regulation", id="reg")
    assert s == avant, "l'entree ne doit pas etre modifiee"
    ss = _le_sous_systeme(r)
    assert ss["id"] == "reg" and ss["nom"] == "regulation"
    assert {b["id"] for b in r["blocs"]} == {"pr", "qr", "reg", "m", "op", "oq"}
    assert len(_ports(ss, "port_entree")) == 2
    assert len(_ports(ss, "port_sortie")) == 1
    assert sorted(b["index"] for b in _ports(ss, "port_entree")) == [0, 1]
    # position : moyenne des positions deplacees ((100 + 250)/2, 100)
    assert ss["x"] == pytest.approx(175.0) and ss["y"] == pytest.approx(100.0)
    # les equations sont celles de depart, aux noms pres
    comparer_equations(compiler(avant), compiler(r), {"e", "c"}, "reg", avant)
    # et le retour
    assert schema_sans_position(B.degrouper(r, "reg")) == schema_sans_position(avant)


@pytest.mark.parametrize(
    "ids, entrants, sortants",
    [
        (["g1", "g2", "s"], 1, 1),  # vd alimente g1 ET g2 : un seul port d'entree
        (["g1", "g2"], 1, 2),  # un entrant, deux signaux de sortie
        (["g1"], 1, 1),
    ],
)
def test_un_port_par_signal(ids, entrants, sortants):
    s = schema_multi()
    r = B.regrouper(s, ids, "groupe", id="grp")
    ss = _le_sous_systeme(r)
    assert len(_ports(ss, "port_entree")) == entrants
    assert len(_ports(ss, "port_sortie")) == sortants
    comparer_equations(compiler(s), compiler(r), set(ids), "grp", s)
    assert schema_sans_position(B.degrouper(r, "grp")) == schema_sans_position(s)


def test_regrouper_sans_id_donne_un_id_valide_et_unique():
    s = schema_multi()
    r = B.regrouper(s, ["g1", "g2"], "sans_id")
    ss = _le_sous_systeme(r)
    assert ss["id"] not in {b["id"] for b in s["blocs"]} and "__" not in ss["id"]
    assert schema_sans_position(B.degrouper(r, ss["id"])) == schema_sans_position(s)


@pytest.mark.parametrize(
    "ids",
    [[], ["inconnu"], ["g1", "inconnu"]],
    ids=["vide", "inconnu", "un_inconnu"],
)
def test_regrouper_refuse_les_ids_invalides(ids):
    with pytest.raises(SchemaInvalide):
        B.regrouper(schema_multi(), ids, "x", id="grp")


def test_regrouper_refuse_un_id_deja_pris():
    with pytest.raises(SchemaInvalide):
        B.regrouper(schema_multi(), ["g1"], "x", id="g2")


def test_regrouper_refuse_des_blocs_de_niveaux_differents():
    """Les blocs regroupes sont ``d'un meme niveau`` : ``g`` est a l'interieur
    de ``ss``, ``a`` au niveau superieur."""
    with pytest.raises(SchemaInvalide):
        B.regrouper(schema_diff(), ["a", "g"], "x", id="grp")


def test_degrouper_refuse_un_bloc_qui_n_est_pas_un_sous_systeme():
    with pytest.raises(SchemaInvalide):
        B.degrouper(schema_multi(), "g1")
    with pytest.raises(SchemaInvalide):
        B.degrouper(schema_multi(), "inconnu")


def test_degrouper_un_sous_systeme_ecrit_a_la_main():
    """Le degroupement du sous-systeme ``ss`` de ``schema_diff`` ne change pas
    les equations : meme valeurs (9 et 3), blocs interieurs au niveau haut."""
    r = B.degrouper(schema_diff(), "ss")
    assert not any(
        b["type"] in ("sous_systeme", "port_entree", "port_sortie") for b in r["blocs"]
    )
    assert {"s", "g"} <= {b["id"] for b in r["blocs"]}
    assert valeurs_statiques(r) == [pytest.approx(9.0), pytest.approx(3.0)]


def test_trois_niveaux_par_regroupements_successifs():
    """Regrouper f, puis le resultat, puis le resultat : trois niveaux ;
    equations identiques au renommage pres, et trois degroupements
    redonnent le schema de depart."""
    s = schema_dynamique()
    r1 = B.regrouper(s, ["f"], "A", id="a")
    r2 = B.regrouper(r1, ["a"], "B", id="b")
    r3 = B.regrouper(r2, ["b"], "C", id="c")
    plat = B.aplatir(copie(r3))
    assert "c__b__a__f" in {b["id"] for b in plat["blocs"]}
    assert compiler(r3).etats == ["x_c__b__a__f"]
    comparer_equations(compiler(s), compiler(r3), {"f"}, "c__b__a", s)
    retour = B.degrouper(B.degrouper(B.degrouper(r3, "c"), "b"), "a")
    assert schema_sans_position(retour) == schema_sans_position(s)


def test_regrouper_un_sous_systeme_existant_avec_d_autres_blocs():
    s = schema_diff()
    r = B.regrouper(s, ["ss", "a"], "enveloppe", id="env")
    ss = _le_sous_systeme(r)
    assert ss["id"] == "env"
    inner = {b["id"]: b for b in ss["schema"]["blocs"]}
    assert inner["ss"]["type"] == "sous_systeme" and "a" in inner
    assert valeurs_statiques(r) == [pytest.approx(9.0), pytest.approx(3.0)]


# =========================== regrouper des blocs de gfl_complet / gfm_complet
def _renommeur(ids_tous: set[str], ids: set[str], ss: str):
    motif = re.compile(r"^(x|y|mode)_(.+)$")

    def renom(nom: str) -> str:
        m = motif.match(nom)
        if not m:
            return nom
        p, reste = m.groups()
        # l'identifiant de bloc le plus long dont ``reste`` est le nom ou le
        # prefixe suivi de ``_`` (``y_pk_d`` -> bloc ``pk``, sortie ``d``)
        cand = [i for i in ids_tous if reste == i or reste.startswith(i + "_")]
        if not cand:
            return nom
        bid = max(cand, key=len)
        if bid not in ids:
            return nom
        return f"{p}_{ss}__{bid}{reste[len(bid) :]}"

    return renom


def comparer_equations(c0, c1, ids: set[str], ss: str, schema0: dict) -> None:
    """Les equations de ``c1`` sont celles de ``c0`` au renommage documenté
    pres (``x_<id>`` -> ``x_<ss>__<id>``, idem ``y_`` et ``mode_``), evaluees
    en trois points pseudo-aleatoires : meme valeur a 1e-11 relatif."""
    ids_tous = {b["id"] for b in schema0["blocs"]}
    renom = _renommeur(ids_tous, ids, ss)
    assert {renom(n) for n in c0.etats} == set(c1.etats)
    assert {renom(n) for n in c0.algebriques} == set(c1.algebriques)
    assert {renom(n) for n in c0.parametres} == set(c1.parametres)
    assert c0.consignes == c1.consignes
    assert len(c0.h) == len(c1.h)
    for n, v in c0.parametres.items():
        assert c1.parametres[renom(n)] == v
    sub = {sp.Symbol(n): sp.Symbol(renom(n)) for n in c0.symboles()}
    paires = (
        [
            (e, c1.f[c1.etats.index(renom(n))])
            for n, e in zip(c0.etats, c0.f, strict=True)
        ]
        + [
            (e, c1.g[c1.algebriques.index(renom(n))])
            for n, e in zip(c0.algebriques, c0.g, strict=True)
        ]
        + list(zip(c0.h, c1.h, strict=True))
    )
    libres = sorted(
        {s for a, b in paires for s in (a.xreplace(sub).free_symbols | b.free_symbols)},
        key=lambda s: s.name,
    )
    for graine in (1, 2, 3):
        rng = np.random.default_rng(graine)
        pt = {s: float(rng.uniform(0.6, 1.4)) for s in libres}
        for a, b in paires:
            va = float(sp.N(a.xreplace(sub).xreplace(pt)))
            vb = float(sp.N(b.xreplace(pt)))
            assert abs(va - vb) <= TOL_EQ * max(1.0, abs(va)), (a, b)


def _eligibles(schema: dict) -> list[str]:
    return [b["id"] for b in schema["blocs"] if b["type"] not in ("consigne", "sortie")]


def _sous_ensembles(schema: dict, n: int, graine: int) -> list[list[str]]:
    rng = np.random.default_rng(graine)
    ids = _eligibles(schema)
    out = []
    for _ in range(n):
        k = int(rng.integers(2, 9))
        out.append(sorted(rng.choice(ids, size=k, replace=False).tolist()))
    return out


@pytest.mark.parametrize("genre", ["gfl", "gfm"])
def test_regrouper_des_blocs_quelconques_laisse_les_equations(genre):
    """Six sous-ensembles tires au hasard (graine fixee) : regrouper les
    laisse identiques au renommage documente pres, et ``degrouper`` rend le
    schema."""
    s = B.schema_gfl() if genre == "gfl" else B.schema_gfm()
    c0 = compiler(copie(s))
    for ids in _sous_ensembles(s, 6, graine=20261107):
        r = B.regrouper(copie(s), ids, "groupe", id="grp")
        comparer_equations(c0, compiler(copie(r)), set(ids), "grp", s)
        assert schema_sans_position(
            B.degrouper(copie(r), "grp")
        ) == schema_sans_position(s), ids


GROUPES = {
    # sous-systemes que le cahier cite : PLL et boucle de courant (GFL),
    # filtre LCL (GFM et GFL)
    "gfl": [
        ["pll", "wf", "wp", "dw", "th"],
        ["ed", "cd", "md", "eqc", "cqc", "mq"],
        ["rl1", "cf", "rl2"],
    ],
    "gfm": [["rl1", "cf", "rl2"], ["ep", "inertie", "dw", "amort", "wv", "th"]],
}


@pytest.mark.parametrize(
    "genre, ids",
    [(g, ids) for g, groupes in GROUPES.items() for ids in groupes],
    ids=lambda v: v if isinstance(v, str) else "+".join(v[:2]),
)
def test_regrouper_des_sous_systemes_nommes_spectre_identique_sur_smib(genre, ids):
    """Le spectre du banc SMIB (bus infini, ligne X = 0,15, 10 MVA) ne change
    pas quand on regroupe un sous-systeme du GFL / GFM, ni les equations."""
    s = B.schema_gfl() if genre == "gfl" else B.schema_gfm()
    r = B.regrouper(copie(s), ids, "sous_sys", id="ss")
    comparer_equations(compiler(copie(s)), compiler(copie(r)), set(ids), "ss", s)
    _n0, lf0, lam0 = spectre(banc_smib(s, genre))
    _n1, lf1, lam1 = spectre(banc_smib(r, genre))
    assert np.max(np.abs(lf0.V - lf1.V)) < 1e-12
    assert ecart_spectres(sans_jauge(lam0), sans_jauge(lam1)) < TOL_SPECTRE


@pytest.mark.parametrize("genre", ["gfl", "gfm"])
def test_regrouper_au_hasard_spectre_identique_sur_smib(genre):
    s = B.schema_gfl() if genre == "gfl" else B.schema_gfm()
    _n0, _lf0, lam0 = spectre(banc_smib(s, genre))
    for ids in _sous_ensembles(s, 2, graine=7):
        r = B.regrouper(copie(s), ids, "groupe", id="grp")
        _n1, _lf1, lam1 = spectre(banc_smib(r, genre))
        assert ecart_spectres(sans_jauge(lam0), sans_jauge(lam1)) < TOL_SPECTRE, ids
