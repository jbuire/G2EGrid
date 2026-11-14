"""Elements decrits par SCHEMA-BLOC : un petit Simulink pour definir de
nouveaux equipements du reseau.

Un schema est un dict JSON ::

    {
      "nom": "source_filtree",
      "parametres": {"tau": 0.05},          # parametres de l'element
      "blocs": [
        {"id": "p", "type": "consigne", "nom": "P_ref"},
        {"id": "f", "type": "premier_ordre", "K": 1, "tau": "tau"},
        {"id": "o", "type": "sortie", "signal": "P"},
        ...
      ],
      "liens": [{"de": "p", "vers": "f", "entree": 0}, {"de": "f", "vers": "o"}]
    }

Chaque bloc a un nombre d'entrees et de sorties fixe par son type
(:data:`TYPES`). La plupart n'ont qu'UNE sortie ; les blocs vectoriels du
repere dq (``park``, ``park_inverse``, ``rl_dq``, ``c_dq``) en ont deux, ``d``
et ``q``, et un lien qui part de la seconde le dit par ``"sortie": 1``. Un champ
numerique d'un bloc (``K``, ``tau``...) est un nombre ou une expression des
parametres du schema.

COMPILATION -- chaque sortie d'un bloc qui n'est pas une source (``entree``,
``constante``, ``consigne``) recoit une variable ALGEBRIQUE (``y_<id>``, ou
``y_<id>_d`` et ``y_<id>_q``), avec le residu ``y - (ce que le bloc calcule de
ses entrees)`` ; les blocs dynamiques ajoutent un ETAT par sortie (``x_<id>``,
``x_<id>_d``...). C'est un peu plus de variables qu'il n'en faudrait, mais les
boucles algebriques -- une somme rebouclee sur elle-meme par un gain -- se
traitent alors comme le reste : le modele est un systeme algebro-
differentiel, ce que le moteur sait faire.

LE REPERE dq. Les tensions du bus (``vd``, ``vq``) et les courants injectes
(``id``, ``iq``) sont dans le repere GLOBAL du reseau. Un convertisseur se
regle dans SON repere, d'angle ``theta`` (celui de sa PLL) : ``park`` y passe,
``park_inverse`` en revient. Comme pour les equipements du moteur, l'angle est
compte depuis le repere de reference du reseau : ces deux blocs retranchent
d'eux-memes ``theta_g`` (l'angle de l'element de reference), si bien qu'une PLL
ecrite ``d(theta)/dt = wb.(w_pll - w)`` est exactement celle des GFL du moteur.
``rl_dq`` (branche R-L) et ``c_dq`` (condensateur) ecrivent les equations d'un
filtre dans un repere tournant a la vitesse ``w`` de leur cinquieme entree ;
avec eux, un GFL complet -- filtre LCL, PLL, boucles de courant et de
puissance -- se dessine bloc a bloc (:func:`schema_gfl`).

PARAMETRES RESERVES -- ``wb`` (pulsation de base du reseau, rad/s) et ``kS``
(``Sn/Sb``, rapport de la puissance de l'element a la base du reseau) sont
fournis par l'assembleur : ``wb`` depuis la frequence du reseau, ``kS`` depuis
la colonne ``Sn_VA`` de ``T_Blocs``. Un schema s'en sert sans les declarer.

INTERFACE AVEC LE RESEAU -- l'element a un port alternatif ``t``. Les blocs
``entree`` lisent ``vd``, ``vq`` (tension du bus, repere global), ``V`` (son
module), ``w`` (frequence du repere du reseau) ou ``theta_g`` ; les blocs
``sortie`` ecrivent soit le courant injecte (``id``, ``iq``, en pu base du
reseau), soit la puissance injectee (``P``, ``Q``), d'ou le courant
``i = conj(S/V_mes)``, sur la tension mesuree a travers un premier ordre de
constante ``tau_mes`` (2 ms par defaut, parametre ajoute au schema s'il n'y est
pas). Deux sorties de SERVICE sont facultatives : ``theta`` (l'angle interne)
et ``w`` (la vitesse). Un element qui les ecrit est un FORMEUR pour le moteur :
pose au bus slack, il porte le repere du reseau (``theta_g`` est son angle), et
le centre d'inertie le compte -- c'est ce qui permet a un GFM a schema-bloc de
tenir seul un ilot.

INITIALISATION -- au load flow, l'element prend la puissance programmee de
son bus (``Pg``/``Qg`` de ``Y_DER``, comme un generateur). Ses ``consigne``
-- exactement deux -- sont les reglages qui la lui font produire : ce sont des
entrees de service (comme le ``p_ref`` d'un GFL), calculees avec les etats et
les algebriques par un Newton amorti sur ``f = 0``, ``g = 0`` et ``i = i_lf``.
Les angles qui alimentent un ``park`` partent de celui de la tension du bus :
c'est ce qui place une PLL sur la bonne branche (la PLL calee a 180 degres est
aussi un equilibre, mais instable).
"""

from __future__ import annotations

import ast
import copy
import hashlib
import json
import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import (
    convert_xor,
    parse_expr,
    standard_transformations,
)

from g2egrid.elements.blocs_structure import (
    BLOC_H,
    BLOC_L,
    PROFONDEUR_MAX,
    TYPES_STRUCTURE,
    SchemaInvalide,
    aplatir,
    degrouper,
    noms_entrees,
    noms_sorties,
    regrouper,
    verifier_forme,
)
from g2egrid.elements.limiteurs import _EPS2, HYSTERESIS, MMAX
from g2egrid.noyau.element import Element, InitContext

__all__ = [
    "BLOC_H",
    "BLOC_L",
    "CATEGORIES",
    "PARAMS_RESERVES",
    "PROFONDEUR_MAX",
    "SIGNAUX_ENTREE",
    "SIGNAUX_SORTIE",
    "TYPES",
    "Compilation",
    "SchemaInvalide",
    "aplatir",
    "classe_element",
    "compiler",
    "degrouper",
    "depuis_element",
    "disposer",
    "exemples",
    "regrouper",
    "schema_gfl",
    "schema_gfm",
]

#: les noms de FR-052 : le module qui les definit importe celui-ci (et, pour
#: la liste du reseau, la couche ``reseau``) : on les charge a la demande
_EQUATIONS = {
    "equations_depuis_element": "g2egrid.elements.blocs_equations",
    "equations_depuis_instance": "g2egrid.elements.blocs_equations",
    "equations_depuis_schema": "g2egrid.elements.blocs_equations",
    "schema_depuis_equations": "g2egrid.elements.blocs_equations",
    "equations_reseau": "g2egrid.reseau.equations_reseau",
    "modifier_equation": "g2egrid.reseau.equations_reseau",
    "ajouter_equation": "g2egrid.reseau.equations_reseau",
    "supprimer_equation": "g2egrid.reseau.equations_reseau",
    "revenir_au_natif": "g2egrid.reseau.equations_reseau",
}


def __getattr__(nom: str):
    if nom in _EQUATIONS:
        import importlib

        return getattr(importlib.import_module(_EQUATIONS[nom]), nom)
    raise AttributeError(f"module {__name__!r} has no attribute {nom!r}")


#: type -> nombre d'entrees (``None`` : il se lit dans le bloc, ``signes``
#: d'une somme ou ``n`` d'une fonction), champs du bloc, dynamique ou non,
#: noms des sorties (une seule, ``y``, si rien n'est dit) et des entrees
#: (pour l'affichage)
TYPES: dict[str, dict[str, Any]] = {
    "entree": dict(entrees=0, champs=("signal",), etat=False, titre="entree du reseau"),
    "constante": dict(entrees=0, champs=("valeur",), etat=False, titre="constante"),
    "consigne": dict(
        entrees=0, champs=("nom",), etat=False, titre="consigne (calee au load flow)"
    ),
    "gain": dict(entrees=1, champs=("K",), etat=False, titre="gain K"),
    "somme": dict(entrees=None, champs=("signes",), etat=False, titre="somme"),
    "produit": dict(entrees=2, champs=(), etat=False, titre="produit"),
    "division": dict(entrees=2, champs=(), etat=False, titre="division u1/u2"),
    "fonction": dict(
        entrees=None, champs=("expr", "n"), etat=False, titre="fonction f(u1..un)"
    ),
    "integrateur": dict(entrees=1, champs=("K",), etat=True, titre="integrateur K/s"),
    "premier_ordre": dict(
        entrees=1, champs=("K", "tau"), etat=True, titre="1er ordre K/(1+tau.s)"
    ),
    "pi": dict(entrees=1, champs=("Kp", "Ki"), etat=True, titre="PI Kp + Ki/s"),
    "derivee": dict(
        entrees=1, champs=("K", "tau"), etat=True, titre="derivee filtree K.s/(1+tau.s)"
    ),
    "park": dict(
        entrees=3,
        champs=(),
        etat=False,
        sorties=("d", "q"),
        noms_entrees=("d", "q", "θ"),
        titre="Park : repere global -> repere d'angle theta",
    ),
    "park_inverse": dict(
        entrees=3,
        champs=(),
        etat=False,
        sorties=("d", "q"),
        noms_entrees=("d", "q", "θ"),
        titre="Park inverse : repere d'angle theta -> repere global",
    ),
    "rl_dq": dict(
        entrees=5,
        champs=("R", "L"),
        etat=True,
        sorties=("d", "q"),
        noms_entrees=("v1d", "v1q", "v2d", "v2q", "ω"),
        titre="branche R-L en dq (courant de v1 vers v2)",
    ),
    "c_dq": dict(
        entrees=5,
        champs=("C",),
        etat=True,
        sorties=("d", "q"),
        noms_entrees=("i1d", "i1q", "i2d", "i2q", "ω"),
        titre="condensateur en dq (charge par i1 - i2)",
    ),
    "limiteur_courant": dict(
        entrees=2,
        champs=("Imax", "Taw"),
        etat=False,
        sorties=("d", "q", "aw_d", "aw_q"),
        noms_entrees=("d*", "q*"),
        titre="limitation de courant |i| <= Imax (et retro-calcul aw)",
    ),
    "saturation_tension": dict(
        entrees=2,
        champs=("Mmax", "Taw"),
        etat=False,
        sorties=("d", "q", "aw_d", "aw_q"),
        noms_entrees=("md*", "mq*"),
        titre="saturation de tension |m| <= Mmax (et retro-calcul aw)",
    ),
    "pi_aw": dict(
        entrees=2,
        champs=("Kp", "Ki"),
        etat=True,
        noms_entrees=("e", "aw"),
        titre="PI anti-emballement : Kp + Ki/s, aw ajoute a l'integrale",
    ),
    "sortie": dict(
        entrees=1,
        champs=("signal",),
        etat=False,
        sorties=(),
        titre="sortie vers le reseau",
    ),
    # --- FR-051 : blocs de la bibliotheque de commande
    "pid": dict(
        entrees=1,
        champs=("Kp", "Ki", "Kd", "Tf"),
        etat=True,
        titre="PID Kp + Ki/s + Kd.s/(1+Tf.s)",
    ),
    "avance_retard": dict(
        entrees=1,
        champs=("K", "T1", "T2"),
        etat=True,
        titre="avance-retard K(1+T1.s)/(1+T2.s)",
    ),
    "second_ordre": dict(
        entrees=1,
        champs=("K", "wn", "zeta"),
        etat=True,
        titre="2e ordre K.wn^2/(s^2 + 2.zeta.wn.s + wn^2)",
    ),
    "fonction_transfert": dict(
        entrees=1,
        champs=("num", "den"),
        etat=True,
        titre="fonction de transfert num(s)/den(s)",
    ),
    "saturation": dict(
        entrees=1,
        champs=("min", "max"),
        etat=False,
        titre="saturation min <= y <= max (libre au point de fonctionnement)",
    ),
    "implicite": dict(
        entrees=None,
        champs=("expr", "n"),
        etat=False,
        titre="contrainte implicite expr(u1..un, y) = 0",
    ),
    # --- structure : sous-systemes (voir blocs_structure)
    "sous_systeme": dict(
        entrees=None,
        champs=("nom", "schema"),
        etat=False,
        sorties=(),
        titre="sous-systeme (entrees et sorties = ses ports)",
    ),
    "port_entree": dict(
        entrees=0,
        champs=("index", "nom"),
        etat=False,
        titre="port d'entree du sous-systeme",
    ),
    "port_sortie": dict(
        entrees=1,
        champs=("index", "nom"),
        etat=False,
        sorties=(),
        titre="port de sortie du sous-systeme",
    ),
}
#: les familles de la palette, dans l'ordre d'affichage : (cle, titre)
CATEGORIES: tuple[tuple[str, str], ...] = (
    ("sources", "Sources"),
    ("operations", "Operations"),
    ("dynamique", "Dynamique"),
    ("repere_dq", "Repere dq"),
    ("limiteurs", "Limiteurs"),
    ("structure", "Structure"),
    ("sorties", "Sorties"),
)
_CATEGORIE_DE = {
    **dict.fromkeys(("entree", "constante", "consigne"), "sources"),
    **dict.fromkeys(
        ("gain", "somme", "produit", "division", "fonction", "implicite"), "operations"
    ),
    **dict.fromkeys(
        (
            "integrateur",
            "premier_ordre",
            "pi",
            "derivee",
            "pi_aw",
            "pid",
            "avance_retard",
            "second_ordre",
            "fonction_transfert",
        ),
        "dynamique",
    ),
    **dict.fromkeys(("park", "park_inverse", "rl_dq", "c_dq"), "repere_dq"),
    **dict.fromkeys(
        ("limiteur_courant", "saturation_tension", "saturation"), "limiteurs"
    ),
    **dict.fromkeys(TYPES_STRUCTURE, "structure"),
    "sortie": "sorties",
}
for _t, _d in TYPES.items():
    _d["categorie"] = _CATEGORIE_DE[_t]
SIGNAUX_ENTREE = ("vd", "vq", "V", "w", "theta_g")
SIGNAUX_SORTIE = ("id", "iq", "P", "Q", "theta", "w")
#: les sorties de SERVICE d'un element a schema-bloc, et le nom que le moteur
#: leur connait : un element qui ecrit son angle et sa vitesse peut porter le
#: repere du reseau (bus slack) et compte parmi les formeurs (centre d'inertie)
SORTIES_SERVICE = {"theta": "theta_out", "w": "w"}
#: ordre maximal d'une ``fonction_transfert``
ORDRE_MAX = 20
#: parametres fournis par l'assembleur (valeur par defaut hors reseau)
PARAMS_RESERVES: dict[str, float] = {"wb": 2 * math.pi * 50.0, "kS": 1.0}
#: constante de temps par defaut du filtre de mesure de tension d'une sortie
#: en puissance (parametre ``tau_mes`` du schema), en s
TAU_MESURE = 2e-3
#: fonctions admises dans les expressions
_FONCTIONS = {
    "sin": sp.sin,
    "cos": sp.cos,
    "tan": sp.tan,
    "exp": sp.exp,
    "log": sp.log,
    "sqrt": sp.sqrt,
    "atan2": sp.atan2,
    "tanh": sp.tanh,
    "Abs": sp.Abs,
    "pi": sp.pi,
}


_GLOBALS = {
    "Symbol": sp.Symbol,
    "Integer": sp.Integer,
    "Float": sp.Float,
    "Rational": sp.Rational,
} | _FONCTIONS
_TRANSFORMATIONS = standard_transformations + (convert_xor,)


#: fonctions appelables dans une expression, par leur nom d'ecriture ; les
#: noms ``abs``, ``min`` et ``max`` (minuscules) s'ajoutent a :data:`_FONCTIONS`
_APPELS = {**_FONCTIONS, "abs": sp.Abs, "min": sp.Min, "max": sp.Max}
_CONSTANTES = {"pi"}
#: longueur maximale d'une expression, et profondeur de son arbre
_EXPR_MAX = 1000
_PROFONDEUR_MAX = 40
#: plus grand exposant constant admis (``x**1000000`` ou ``9**9**9`` figeraient
#: SymPy a calculer un entier gigantesque)
_EXPOSANT_MAX = 64.0
_OPERATEURS = (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub, ast.UAdd)


def _valeur_constante(n: ast.AST) -> float | None:
    """La valeur d'un sous-arbre sans symbole (en flottant, jamais en entier
    exact), ou ``None`` s'il en contient un."""
    if isinstance(n, ast.Constant):
        return float(n.value)
    if isinstance(n, ast.UnaryOp):
        v = _valeur_constante(n.operand)
        return None if v is None else (-v if isinstance(n.op, ast.USub) else v)
    if isinstance(n, ast.BinOp):
        a, b = _valeur_constante(n.left), _valeur_constante(n.right)
        if a is None or b is None:
            return None
        try:
            if isinstance(n.op, ast.Add):
                return a + b
            if isinstance(n.op, ast.Sub):
                return a - b
            if isinstance(n.op, ast.Mult):
                return a * b
            if isinstance(n.op, ast.Div):
                return a / b
            return float(a**b) if abs(b) <= _EXPOSANT_MAX else float("inf")
        except (ArithmeticError, ValueError):
            return float("inf")
    return None


def _verifier_expression(
    texte: str, noms: set[str], longueur_max: int = _EXPR_MAX
) -> None:
    """Refuse toute expression qui ne soit pas de l'ARITHMETIQUE sur des noms.

    ``parse_expr`` de SymPy ``eval`` le texte : ``__import__('os').system(..)``
    s'executait. On lit donc d'abord le texte avec :mod:`ast` (qui n'execute
    rien), et on n'admet que des nombres, les noms de ``noms`` et les
    constantes de :data:`_CONSTANTES`, ``+ - * / **`` et l'appel direct des
    fonctions de :data:`_APPELS`. Attributs, indices, lambdas, comprehensions,
    chaines, mots-cles d'appel : tout le reste est refuse avec la raison.
    """
    if len(texte) > longueur_max:
        raise SchemaInvalide(f"expression trop longue ({len(texte)} > {longueur_max})")
    try:
        arbre = ast.parse(texte.strip(), mode="eval")
    except (SyntaxError, ValueError, RecursionError, MemoryError) as err:
        raise SchemaInvalide(f"{texte!r} illisible") from err

    def visiter(n: ast.AST, profondeur: int) -> None:
        if profondeur > _PROFONDEUR_MAX:
            raise SchemaInvalide(f"{texte!r} trop imbriquee")
        profondeur += 1
        if isinstance(n, ast.Expression):
            visiter(n.body, profondeur)
        elif isinstance(n, ast.Constant):
            if isinstance(n.value, bool) or not isinstance(n.value, (int, float)):
                raise SchemaInvalide(f"{texte!r} : seuls les nombres sont admis")
            if isinstance(n.value, int) and abs(n.value) > 10**18:
                raise SchemaInvalide(f"{texte!r} : entier trop grand")
            if isinstance(n.value, float) and not math.isfinite(n.value):
                # 1e999 se lit comme l'infini
                raise SchemaInvalide(f"{texte!r} : nombre hors de portee")
        elif isinstance(n, ast.Name):
            if n.id not in noms and n.id not in _CONSTANTES:
                raise SchemaInvalide(f"{texte!r} emploie {[n.id]}, inconnus")
        elif isinstance(n, ast.UnaryOp) and isinstance(n.op, _OPERATEURS):
            visiter(n.operand, profondeur)
        elif isinstance(n, ast.BinOp) and isinstance(n.op, _OPERATEURS):
            if isinstance(n.op, ast.Pow):
                e = _valeur_constante(n.right)
                if e is not None and abs(e) > _EXPOSANT_MAX:
                    raise SchemaInvalide(
                        f"{texte!r} : exposant constant {e:g} trop grand "
                        f"(au plus {_EXPOSANT_MAX:g})"
                    )
            visiter(n.left, profondeur)
            visiter(n.right, profondeur)
        elif isinstance(n, ast.Call):
            if not isinstance(n.func, ast.Name) or n.func.id not in _APPELS:
                raise SchemaInvalide(
                    f"{texte!r} : appel refuse, les fonctions admises sont "
                    f"{sorted(_APPELS)}"
                )
            if n.keywords or any(isinstance(a, ast.Starred) for a in n.args):
                raise SchemaInvalide(f"{texte!r} : arguments nommes ou etoiles refuses")
            for a in n.args:
                visiter(a, profondeur)
        else:
            raise SchemaInvalide(
                f"{texte!r} : {type(n).__name__} refuse (arithmetique sur les "
                "noms du schema seulement : + - * / ** et fonctions usuelles)"
            )

    visiter(arbre, 0)


def _expr(texte: str, noms: set[str], longueur_max: int = _EXPR_MAX) -> sp.Expr:
    """``texte`` lu en expression sympy : les identifiants qui ne sont pas des
    fonctions de :data:`_FONCTIONS` deviennent des SYMBOLES (et jamais une
    fonction homonyme de sympy, comme ``beta`` ou ``gamma``) ; ceux qui ne
    sont pas dans ``noms`` font lever :class:`SchemaInvalide`.

    Le texte est d'abord verifie par :func:`_verifier_expression` : SymPy ne
    voit jamais autre chose que de l'arithmetique sur des noms connus."""
    _verifier_expression(
        texte.replace("^", "**"), noms, longueur_max
    )  # ``convert_xor``
    try:
        e = parse_expr(
            texte,
            local_dict={n: sp.Symbol(n) for n in noms},
            global_dict=_GLOBALS | _APPELS,
            transformations=_TRANSFORMATIONS,
        )
    except (sp.SympifyError, SyntaxError, TypeError, NameError) as err:
        raise SchemaInvalide(f"{texte!r} illisible") from err
    libres = {str(x) for x in e.free_symbols} - noms
    if libres:
        raise SchemaInvalide(f"{texte!r} emploie {sorted(libres)}, inconnus")
    return e


def _nom_ok(nom: str) -> bool:
    return nom.isidentifier() and not nom.startswith("_")


def _nb_entrees(b: dict) -> int:
    t = TYPES[b["type"]]
    if t["entrees"] is not None:
        return t["entrees"]
    if b["type"] == "somme":
        return len(str(b.get("signes", "++")))
    if b["type"] == "sous_systeme":
        return len(noms_entrees(b))
    return int(b.get("n", 1))


def _noms_sorties(b: dict) -> tuple[str, ...]:
    if b["type"] == "sous_systeme":
        return tuple(noms_sorties(b))
    return TYPES[b["type"]].get("sorties", ("y",))


def _var(prefixe: str, bid: str, sorties: tuple[str, ...], k: int) -> str:
    """Le nom de la variable de la sortie ``k`` d'un bloc : ``y_<id>`` pour un
    bloc a une sortie, ``y_<id>_<sortie>`` sinon."""
    return f"{prefixe}_{bid}" if len(sorties) == 1 else f"{prefixe}_{bid}_{sorties[k]}"


@dataclass
class Compilation:
    """Le schema compile : noms des variables et equations symboliques.

    Les equations sont ecrites sur des SYMBOLES de meme nom que les variables
    de l'element ; :meth:`equations` les reecrit sur l'espace de noms du
    moteur (``s``), et :meth:`latex` les rend lisibles."""

    nom: str
    parametres: dict[str, float]
    etats: list[str]
    algebriques: list[str]
    consignes: list[str]
    #: d(etat)/dt, residus algebriques, et (id, iq) injectes
    f: list[sp.Expr]
    g: list[sp.Expr]
    h: list[sp.Expr]
    #: pour l'affichage : (genre, variable, membre de droite) de chaque bloc
    lignes: list[tuple[str, str, sp.Expr]] = field(default_factory=list)
    #: les variables qui portent un angle de repere (entree theta d'un park) :
    #: l'initialisation les fait partir de l'angle de la tension du bus
    angles: list[str] = field(default_factory=list)
    #: les parametres reserves que le schema emploie (fournis par l'assembleur)
    reserves: list[str] = field(default_factory=list)
    #: sorties de service (``theta_out``, ``w``), apres les deux courants
    services: list[str] = field(default_factory=list)
    #: les limiteurs : (drapeau de mode, module de la consigne libre, Imax),
    #: de quoi ecrire les gardes de la simulation temporelle
    limiteurs: list[tuple[str, sp.Expr, sp.Expr]] = field(default_factory=list)
    #: les saturations libres : (drapeau, entree, min, max), de quoi ecrire
    #: les gardes de la simulation temporelle
    saturations: list[tuple[str, sp.Expr, sp.Expr, sp.Expr]] = field(
        default_factory=list
    )
    #: les consignes CALEES par l'initialisation (exactement deux), et la
    #: valeur des autres (``calee: false``), entrees de service fixees
    calees: list[str] = field(default_factory=list)
    valeurs: dict[str, sp.Expr] = field(default_factory=dict)
    empreinte: str = ""
    #: valeurs de depart du Newton d'initialisation, par variable (``depart``
    #: d'un etat ou d'une variable algebrique d'un modele a equations)
    departs: dict[str, float] = field(default_factory=dict)

    def symboles(self) -> dict[str, sp.Symbol]:
        noms = (
            list(self.parametres)
            + self.etats
            + self.algebriques
            + self.consignes
            + ["vd_t", "vq_t", "w_ref", "theta_g"]
        )
        return {n: sp.Symbol(n) for n in noms}

    def equations(self, s):
        """``(f, g, h)`` sur l'espace de noms ``s`` du moteur."""
        sub = {sym: getattr(s, n) for n, sym in self.symboles().items()}
        return (
            [e.xreplace(sub) for e in self.f],
            [e.xreplace(sub) for e in self.g],
            [e.xreplace(sub) for e in self.h],
        )

    def latex(self) -> list[dict[str, str]]:
        """Les equations du schema, une par bloc, en LaTeX et en texte."""
        out = []
        for genre, var, expr in self.lignes:
            if genre == "implicite":  # 0 = residu(u, y), y = var
                out.append(
                    dict(
                        genre="algebrique",
                        latex=f"0 = {sp.latex(expr)}",
                        texte=f"0 = {sp.sstr(expr)}  (resolue pour {var})",
                    )
                )
                continue
            gauche = (
                rf"\dot{{{sp.latex(sp.Symbol(var))}}}"
                if genre == "etat"
                else sp.latex(sp.Symbol(var))
            )
            texte = f"d{var}/dt" if genre == "etat" else var
            out.append(
                dict(
                    genre=genre,
                    latex=f"{gauche} = {sp.latex(expr)}",
                    texte=f"{texte} = {sp.sstr(expr)}",
                )
            )
        return out


def lire_parametres(parametres: Any) -> dict[str, float]:
    """``{nom: valeur}`` des parametres d'un modele : chaque valeur est un
    nombre FINI (ni booleen, ni NaN, ni infini -- que ``json.loads`` accepte
    pourtant). Leve :class:`SchemaInvalide` en nommant le parametre."""
    if parametres is None:
        return {}
    if not isinstance(parametres, dict):
        raise SchemaInvalide("parametres est un objet {nom: valeur}")
    out: dict[str, float] = {}
    for k, v in parametres.items():
        try:
            if isinstance(v, bool):
                raise ValueError  # noqa: TRY004 -- refuse comme un texte
            x = float(v)  # un nombre ecrit en texte ("0.05") reste admis
        except (TypeError, ValueError):
            raise SchemaInvalide(
                f"parametre {k!r} : {v!r} n'est pas un nombre"
            ) from None
        if not math.isfinite(x):
            raise SchemaInvalide(f"parametre {k!r} : la valeur doit etre finie ({v!r})")
        out[str(k)] = x
    return out


def compiler(schema: dict) -> Compilation:
    """Compile un schema-bloc. Leve :class:`SchemaInvalide` avec la raison.

    Les sous-systemes sont d'abord dissous par :func:`aplatir` : la
    compilation d'un schema et celle de son aplati sont la meme."""
    if isinstance(schema, dict) and schema.get("genre") == "equations":
        from g2egrid.elements.blocs_equations import compiler_equations

        return compiler_equations(schema)
    verifier_forme(schema)
    origine = schema
    schema = aplatir(schema)
    nom = str(schema.get("nom") or "bloc")
    if not _nom_ok(nom):
        raise SchemaInvalide(
            f"nom de modele {nom!r} : lettres, chiffres et _ seulement"
        )
    params = lire_parametres(schema.get("parametres"))
    for k in params:
        if not _nom_ok(k):
            raise SchemaInvalide(f"parametre {k!r} : nom invalide")
        if k in PARAMS_RESERVES:
            raise SchemaInvalide(
                f"parametre {k!r} reserve : il est fourni par le reseau "
                "(wb : pulsation de base ; kS : Sn/Sb, colonne Sn_VA de T_Blocs)"
            )
    admis = set(params) | set(PARAMS_RESERVES)
    blocs = list(schema.get("blocs") or [])
    if not blocs:
        raise SchemaInvalide("schema vide")
    par_id: dict[str, dict] = {}
    for b in blocs:
        if b.get("type") not in TYPES:
            raise SchemaInvalide(
                f"bloc {b.get('id')!r} : type {b.get('type')!r} inconnu"
            )
        bid = str(b.get("id"))
        if not _nom_ok(bid) or bid in par_id:
            raise SchemaInvalide(f"identifiant de bloc {bid!r} invalide ou en double")
        par_id[bid] = b

    # entrees de chaque bloc, depuis les liens : (bloc source, sortie)
    entrees: dict[str, list[tuple[str, int] | None]] = {
        bid: [None] * _nb_entrees(b) for bid, b in par_id.items()
    }
    for lien in schema.get("liens") or []:
        de, vers = str(lien.get("de")), str(lien.get("vers"))
        k, ks = int(lien.get("entree", 0)), int(lien.get("sortie", 0))
        if de not in par_id or vers not in par_id:
            raise SchemaInvalide(f"lien {de} -> {vers} : bloc inconnu")
        n_s = len(_noms_sorties(par_id[de]))
        if n_s == 0:
            raise SchemaInvalide(
                f"lien depuis la sortie {de} : une sortie n'a pas de sortie"
            )
        if not 0 <= ks < n_s:
            raise SchemaInvalide(f"lien {de} -> {vers} : sortie {ks} inexistante")
        if not 0 <= k < len(entrees[vers]):
            raise SchemaInvalide(f"lien {de} -> {vers} : entree {k} inexistante")
        if entrees[vers][k] is not None:
            raise SchemaInvalide(f"bloc {vers} : entree {k + 1} deja reliee")
        entrees[vers][k] = (de, ks)
    for bid, ent in entrees.items():
        for k, e in enumerate(ent):
            if e is None:
                raise SchemaInvalide(f"bloc {bid} : entree {k + 1} non reliee")

    def champ(b: dict, cle: str, defaut: float | None = None) -> sp.Expr:
        return valeur(b, cle, b.get(cle, defaut))

    def valeur(b: dict, cle: str, v: Any) -> sp.Expr:
        if v is None or v == "":
            raise SchemaInvalide(f"bloc {b['id']} : champ {cle!r} manquant")
        if not isinstance(v, str):
            try:
                if isinstance(v, bool):
                    raise ValueError("booleen")
                x = float(v)
                if not math.isfinite(x):
                    raise ValueError("valeur non finie")
                return sp.Integer(int(x)) if x.is_integer() else sp.Float(x)
            except (TypeError, ValueError) as err:
                raise SchemaInvalide(
                    f"bloc {b['id']} : {cle} = {v!r} illisible"
                ) from err
        try:
            return _expr(v, admis)
        except SchemaInvalide as err:
            raise SchemaInvalide(
                f"bloc {b['id']} : {cle} = {err} (seuls les parametres du schema, "
                "wb et kS sont admis)"
            ) from None

    vd, vq, w = sp.Symbol("vd_t"), sp.Symbol("vq_t"), sp.Symbol("w_ref")
    th_g, wb = sp.Symbol("theta_g"), sp.Symbol("wb")
    signaux_entree = {
        "vd": vd,
        "vq": vq,
        "V": sp.sqrt(vd**2 + vq**2),
        "w": w,
        "theta_g": th_g,
    }

    etats, algebs, consignes = [], [], []
    noncalees: dict[str, sp.Expr] = {}
    f, g, lignes = [], [], []
    sortie: dict[str, sp.Expr] = {}
    # la valeur de chaque sortie de chaque bloc : un symbole (algebrique) ou
    # une expression (source)
    val: dict[str, list[sp.Expr]] = {}
    for bid, b in par_id.items():
        t = b["type"]
        if t == "entree":
            s = str(b.get("signal"))
            if s not in SIGNAUX_ENTREE:
                raise SchemaInvalide(
                    f"entree {bid} : signal {s!r} (attendu : {', '.join(SIGNAUX_ENTREE)})"
                )
            val[bid] = [signaux_entree[s]]
        elif t == "constante":
            val[bid] = [champ(b, "valeur")]
        elif t == "consigne":
            n = str(b.get("nom") or bid)
            if not _nom_ok(n) or n in consignes or n in admis:
                raise SchemaInvalide(
                    f"consigne {bid} : nom {n!r} invalide ou en double"
                )
            consignes.append(n)
            calee = b.get("calee", True)
            if not isinstance(calee, bool):
                raise SchemaInvalide(f"consigne {bid} : calee doit etre vrai ou faux")
            if not calee:
                noncalees[n] = champ(b, "valeur", 0)
            val[bid] = [sp.Symbol(n)]
        elif t != "sortie":
            noms_s = _noms_sorties(b)
            val[bid] = []
            for k in range(len(noms_s)):
                y = _var("y", bid, noms_s, k)
                algebs.append(y)
                val[bid].append(sp.Symbol(y))
    tous = algebs + consignes + list(params)
    if len(set(tous)) != len(tous):
        raise SchemaInvalide("deux variables du schema portent le meme nom")

    angles: list[str] = []
    limiteurs: list[tuple[str, sp.Expr, sp.Expr]] = []
    saturations: list[tuple[str, sp.Expr, sp.Expr, sp.Expr]] = []
    valeurs_hors_reseau = {
        sp.Symbol(k): sp.Float(v) for k, v in (params | PARAMS_RESERVES).items()
    }

    def positif(b: dict, cle: str, e: sp.Expr) -> None:
        """Refuse un champ qui doit etre strictement positif et ne l'est pas
        (quand sa valeur se lit : nombre, ou expression des parametres)."""
        try:
            x = float(e.xreplace(valeurs_hors_reseau))
        except (TypeError, ValueError):
            return
        if not x > 0:
            raise SchemaInvalide(
                f"bloc {b['id']} : {cle} = {x:g}, il doit etre strictement positif"
            )

    def non_nul(b: dict, cle: str, e: sp.Expr) -> None:
        """Refuse un champ DIVISEUR qui vaut zero (quand sa valeur se lit) :
        ``tau = 0`` d'un premier ordre compilait, puis divisait par zero a
        l'execution. Une valeur negative reste admise (modele instable voulu)."""
        try:
            x = float(e.xreplace(valeurs_hors_reseau))
        except (TypeError, ValueError):
            return
        if x == 0:
            raise SchemaInvalide(f"bloc {b['id']} : {cle} = 0, on divise par ce champ")

    for bid, b in par_id.items():
        t = b["type"]
        u = [val[de][ks] for de, ks in entrees[bid]]
        if t in ("entree", "constante", "consigne"):
            continue
        if t == "sortie":
            s = str(b.get("signal"))
            if s not in SIGNAUX_SORTIE:
                raise SchemaInvalide(
                    f"sortie {bid} : signal {s!r} (attendu : {', '.join(SIGNAUX_SORTIE)})"
                )
            if s in sortie:
                raise SchemaInvalide(f"sortie {s} ecrite deux fois")
            sortie[s] = u[0]
            continue
        noms_s = _noms_sorties(b)
        ys = [sp.Symbol(_var("y", bid, noms_s, k)) for k in range(len(noms_s))]
        xs = [sp.Symbol(_var("x", bid, noms_s, k)) for k in range(len(noms_s))]

        def etat(x: sp.Symbol, derivee: sp.Expr) -> None:
            etats.append(str(x))
            f.append(derivee)
            lignes.append(("etat", str(x), derivee))

        rhs: list[sp.Expr] | None
        if t == "gain":
            rhs = [champ(b, "K") * u[0]]
        elif t == "somme":
            signes = str(b.get("signes", "++"))
            if set(signes) - {"+", "-"}:
                raise SchemaInvalide(
                    f"somme {bid} : signes {signes!r} (des + et des -)"
                )
            rhs = [sum((1 if c == "+" else -1) * ui for c, ui in zip(signes, u))]
        elif t == "produit":
            rhs = [u[0] * u[1]]
        elif t == "division":
            rhs = [u[0] / u[1]]
        elif t == "fonction":
            n = len(u)
            us = {f"u{k + 1}": ui for k, ui in enumerate(u)}
            if n == 1:
                us["u"] = u[0]
            expr = str(b.get("expr", "u"))
            try:
                e = _expr(expr, admis | set(us))
            except SchemaInvalide as err:
                raise SchemaInvalide(f"fonction {bid} : {err}") from None
            rhs = [e.xreplace({sp.Symbol(k): v for k, v in us.items()})]
        elif t == "integrateur":
            etat(xs[0], champ(b, "K", 1) * u[0])
            rhs = [xs[0]]
        elif t == "premier_ordre":
            K, tau = champ(b, "K", 1), champ(b, "tau")
            non_nul(b, "tau", tau)
            etat(xs[0], (K * u[0] - xs[0]) / tau)
            rhs = [xs[0]]
        elif t == "pi":
            Kp, Ki = champ(b, "Kp"), champ(b, "Ki")
            etat(xs[0], Ki * u[0])
            rhs = [Kp * u[0] + xs[0]]
        elif t == "derivee":
            K, tau = champ(b, "K", 1), champ(b, "tau")
            non_nul(b, "tau", tau)
            etat(xs[0], (K * u[0] - xs[0]) / tau)
            rhs = [(K * u[0] - xs[0]) / tau]
        elif t in ("park", "park_inverse"):
            # l'angle est compte depuis le repere de reference du reseau
            a = u[2] - th_g
            c, s_ = sp.cos(a), sp.sin(a)
            if t == "park":
                rhs = [u[0] * c + u[1] * s_, -u[0] * s_ + u[1] * c]
            else:
                rhs = [u[0] * c - u[1] * s_, u[0] * s_ + u[1] * c]
            de, ks = entrees[bid][2]
            if par_id[de]["type"] not in ("entree", "constante", "consigne"):
                src = _noms_sorties(par_id[de])
                angles.append(_var("y", de, src, ks))
                if (
                    TYPES[par_id[de]["type"]]["etat"]
                    and par_id[de]["type"] not in _ETATS_NUMEROTES
                ):
                    angles.append(_var("x", de, src, ks))
        elif t == "rl_dq":
            R, L = champ(b, "R"), champ(b, "L")
            non_nul(b, "L", L)
            v1d, v1q, v2d, v2q, wr = u
            etat(xs[0], (wb / L) * (v1d - v2d - R * xs[0]) + wb * wr * xs[1])
            etat(xs[1], (wb / L) * (v1q - v2q - R * xs[1]) - wb * wr * xs[0])
            rhs = list(xs)
        elif t in ("limiteur_courant", "saturation_tension"):
            # comme limiter_courant / limiter_tension du moteur : un MODE
            # (drapeau mode_<id>, parametre a 0 en fonctionnement libre)
            # plutot qu'un min dans les equations ; la simulation le bascule
            # sur ses gardes. A drapeau nul, le vecteur passe tel quel et le
            # retro-calcul est nul : le modele linearise est celui du
            # fonctionnement libre. Une borne <= 0 desactive la garde (la
            # saturation de tension du moteur, Mmax = 0, est ainsi « non
            # reglee »).
            drapeau = f"mode_{bid}"
            if drapeau in admis or drapeau in consignes:
                raise SchemaInvalide(f"limiteur {bid} : {drapeau} deja pris")
            m = sp.Symbol(drapeau)
            if t == "limiteur_courant":
                Imax, Taw = champ(b, "Imax", 1.2), champ(b, "Taw", 2e-3)
            else:
                Imax, Taw = champ(b, "Mmax", 1.15), champ(b, "Taw", 20e-3)
            module = sp.sqrt(u[0] ** 2 + u[1] ** 2 + _EPS2)
            k = 1 - m + m * sp.Min(1, Imax / module)
            rhs = [
                u[0] * k,
                u[1] * k,
                m * (u[0] * k - u[0]) / Taw,
                m * (u[1] * k - u[1]) / Taw,
            ]
            limiteurs.append((drapeau, module, Imax))
        elif t == "pi_aw":
            Kp, Ki = champ(b, "Kp"), champ(b, "Ki")
            etat(xs[0], Ki * u[0] + u[1])
            rhs = [Kp * u[0] + xs[0]]
        elif t == "c_dq":
            C = champ(b, "C")
            non_nul(b, "C", C)
            i1d, i1q, i2d, i2q, wr = u
            etat(xs[0], (wb / C) * (i1d - i2d) + wb * wr * xs[1])
            etat(xs[1], (wb / C) * (i1q - i2q) - wb * wr * xs[0])
            rhs = list(xs)
        elif t == "pid":
            # Kp + Ki/s + Kd.s/(1 + Tf.s) : un integrateur (x_1) et le filtre
            # de la derivee (x_2, meme forme que ``derivee``)
            Kp, Ki = champ(b, "Kp", 1), champ(b, "Ki", 0)
            Kd, Tf = champ(b, "Kd", 0), champ(b, "Tf")
            positif(b, "Tf", Tf)
            x1, x2 = sp.Symbol(f"x_{bid}_1"), sp.Symbol(f"x_{bid}_2")
            etat(x1, Ki * u[0])
            etat(x2, (Kd * u[0] - x2) / Tf)
            rhs = [Kp * u[0] + x1 + (Kd * u[0] - x2) / Tf]
        elif t == "avance_retard":
            # K(1 + T1.s)/(1 + T2.s) = K.[T1/T2 + (1 - T1/T2)/(1 + T2.s)]
            K, T1, T2 = champ(b, "K", 1), champ(b, "T1"), champ(b, "T2")
            positif(b, "T2", T2)
            etat(xs[0], (u[0] - xs[0]) / T2)
            rhs = [K * (T1 / T2 * u[0] + (1 - T1 / T2) * xs[0])]
        elif t == "second_ordre":
            # x_1 = y, x_2 = dy/dt : y'' = K.wn^2.u - 2.zeta.wn.y' - wn^2.y
            K, wn, zeta = champ(b, "K", 1), champ(b, "wn"), champ(b, "zeta")
            positif(b, "wn", wn)
            x1, x2 = sp.Symbol(f"x_{bid}_1"), sp.Symbol(f"x_{bid}_2")
            etat(x1, x2)
            etat(x2, K * wn**2 * u[0] - 2 * zeta * wn * x2 - wn**2 * x1)
            rhs = [x1]
        elif t == "fonction_transfert":
            rhs = _fonction_transfert(b, bid, u[0], valeur, etat)
        elif t == "saturation":
            # libre au point de fonctionnement (y = u), comme un limiteur : un
            # MODE (drapeau mode_<id>) que la simulation bascule sur ses gardes
            drapeau = f"mode_{bid}"
            if drapeau in admis or drapeau in consignes:
                raise SchemaInvalide(f"saturation {bid} : {drapeau} deja pris")
            lo, hi = champ(b, "min"), champ(b, "max")
            try:
                bornes = [float(x.xreplace(valeurs_hors_reseau)) for x in (lo, hi)]
            except (TypeError, ValueError):
                bornes = None
            if bornes and bornes[0] > bornes[1]:
                # min > max : la sortie vaudrait max quel que soit l'entree
                raise SchemaInvalide(
                    f"saturation {bid} : min = {bornes[0]:g} depasse max = {bornes[1]:g}"
                )
            m = sp.Symbol(drapeau)
            rhs = [u[0] + m * (sp.Min(sp.Max(u[0], lo), hi) - u[0])]
            saturations.append((drapeau, u[0], lo, hi))
        elif t == "implicite":
            # y tel que expr(u1..un, y) = 0 : une contrainte, pas une affectation
            us = {f"u{k + 1}": ui for k, ui in enumerate(u)}
            if len(u) == 1:
                us["u"] = u[0]
            try:
                e = _expr(str(b.get("expr", "")), admis | set(us) | {"y"})
            except SchemaInvalide as err:
                raise SchemaInvalide(f"implicite {bid} : {err}") from None
            if sp.Symbol("y") not in e.free_symbols:
                raise SchemaInvalide(
                    f"implicite {bid} : l'expression ne contient pas y, "
                    "la contrainte ne definit aucune sortie"
                )
            residu = e.xreplace(
                {sp.Symbol(k): v for k, v in us.items()} | {sp.Symbol("y"): ys[0]}
            )
            g.append(residu)
            lignes.append(("implicite", str(ys[0]), residu))
            rhs = None
        else:  # pragma: no cover - TYPES et ce if vont ensemble
            raise SchemaInvalide(f"bloc {bid} : type {t} sans regle")
        for y, r in zip(ys, rhs or ()):
            g.append(y - r)
            lignes.append(("algebrique", str(y), r))

    return _conclure(
        origine,
        nom,
        params,
        etats,
        algebs,
        consignes,
        noncalees,
        f,
        g,
        lignes,
        sortie,
        (vd, vq),
        angles,
        limiteurs,
        saturations,
    )


def _conclure(
    origine: dict,
    nom: str,
    params: dict[str, float],
    etats: list[str],
    algebs: list[str],
    consignes: list[str],
    noncalees: dict[str, sp.Expr],
    f: list[sp.Expr],
    g: list[sp.Expr],
    lignes: list[tuple[str, str, sp.Expr]],
    sortie: dict[str, sp.Expr],
    tension: tuple[sp.Symbol, sp.Symbol],
    angles: list[str],
    limiteurs: list[tuple[str, sp.Expr, sp.Expr]],
    saturations: list[tuple[str, sp.Expr, sp.Expr, sp.Expr]],
    departs: dict[str, float] | None = None,
) -> Compilation:
    """La fin commune des deux compilateurs (schema-bloc, modele a equations) :
    les sorties vers le reseau (courant injecte, ou puissance -- et alors les
    deux etats de la tension mesuree), les sorties de service, les consignes
    calees, les parametres reserves. ``sortie`` : signal -> expression, sur
    les symboles du moteur (``vd_t``...) ; ``tension`` : ``(vd_t, vq_t)``."""
    vd, vq = tension
    # le courant injecte
    if {"id", "iq"} <= set(sortie) and not ({"P", "Q"} & set(sortie)):
        h = [sortie["id"], sortie["iq"]]
    elif {"P", "Q"} <= set(sortie) and not ({"id", "iq"} & set(sortie)):
        # Le courant se tire de la tension MESUREE, filtree au premier ordre
        # (constante ``tau_mes``) comme dans un vrai convertisseur : calculee
        # sur la tension instantanee du bus, une injection a puissance
        # constante a une conductance incrementale negative qui deferait la
        # moindre resonance du reseau.
        P, Q = sortie["P"], sortie["Q"]
        params.setdefault("tau_mes", TAU_MESURE)
        tm = sp.Symbol("tau_mes")
        vdm, vqm = sp.Symbol("vd_mes"), sp.Symbol("vq_mes")
        for xm, v in ((vdm, vd), (vqm, vq)):
            etats.append(str(xm))
            f.append((v - xm) / tm)
            lignes.append(("etat", str(xm), f[-1]))
        V2 = vdm**2 + vqm**2
        h = [(P * vdm + Q * vqm) / V2, (P * vqm - Q * vdm) / V2]
    else:
        raise SchemaInvalide(
            "il faut deux sorties vers le reseau : id et iq (courant injecte), "
            "ou P et Q (puissance injectee)"
        )
    lignes.append(("algebrique", "i_d", h[0]))
    lignes.append(("algebrique", "i_q", h[1]))
    services = []
    for s_, nom_moteur in SORTIES_SERVICE.items():
        if s_ in sortie:
            h.append(sortie[s_])
            services.append(nom_moteur)
            lignes.append(("algebrique", nom_moteur, sortie[s_]))
    calees = [n for n in consignes if n not in noncalees]
    if len(calees) != 2:
        raise SchemaInvalide(
            f"{len(calees)} consigne(s) calee(s) : il en faut exactement deux, que "
            "l'initialisation cale sur la puissance active et reactive du load flow "
            "(une consigne qui ne doit pas l'etre porte calee = false et une valeur)"
        )
    # les parametres reserves employes rejoignent les parametres, avec leur
    # valeur hors reseau ; l'assembleur les remplace
    libres = set().union(
        *(e.free_symbols for e in f + g + h),
        *(e.free_symbols for e in noncalees.values()),
    )
    reserves = [k for k in PARAMS_RESERVES if sp.Symbol(k) in libres]
    for k in reserves:
        params[k] = PARAMS_RESERVES[k]
    # les drapeaux des limiteurs : des parametres, nuls (fonctionnement libre)
    for drapeau, _m, _i in limiteurs:
        params[drapeau] = 0.0
    for drapeau, *_r in saturations:
        params[drapeau] = 0.0
    _refuser_division_par_zero(params, lignes, noncalees)
    empreinte = hashlib.sha1(
        json.dumps(_sans_positions(origine), sort_keys=True, default=str).encode()
    ).hexdigest()[:10]
    return Compilation(
        nom=nom,
        parametres=params,
        etats=etats,
        algebriques=algebs,
        consignes=consignes,
        f=f,
        g=g,
        h=h,
        lignes=lignes,
        angles=list(dict.fromkeys(angles)),
        reserves=reserves,
        services=services,
        limiteurs=limiteurs,
        saturations=saturations,
        calees=calees,
        valeurs=noncalees,
        empreinte=empreinte,
        departs=dict(departs or {}),
    )


_NON_FINIS = (sp.zoo, sp.oo, -sp.oo, sp.nan)


def _refuser_division_par_zero(
    params: dict[str, float],
    lignes: list[tuple[str, str, sp.Expr]],
    valeurs: dict[str, sp.Expr],
) -> None:
    """Refuse un parametre NUL au denominateur d'une equation (``L = 0`` dans
    ``wb*(...)/L``). Le modele compilait, puis le load flow echouait sur un
    residu infini, avec un diagnostic qui envoyait chercher du cote des
    consignes. Seuls les parametres nuls sont substitues : le controle ne
    coute rien a un modele qui n'en cite aucun."""
    nuls = {sp.Symbol(k): sp.Integer(0) for k, v in params.items() if v == 0}
    if not nuls:
        return
    equations = [(var, e) for _genre, var, e in lignes]
    equations += [(f"la consigne {k}", e) for k, e in valeurs.items()]
    for var, e in equations:
        communs = [s for s in nuls if s in e.free_symbols]
        if not communs:
            continue
        if not e.xreplace({s: nuls[s] for s in communs}).has(*_NON_FINIS):
            continue
        seuls = [s for s in communs if e.xreplace({s: nuls[s]}).has(*_NON_FINIS)]
        noms = ", ".join(sorted(str(s) for s in (seuls or communs)))
        raise SchemaInvalide(
            f"equation de {var} : division par zero, le parametre {noms} vaut 0 "
            "(lui donner une valeur non nulle, ou retirer le terme)"
        )


#: les blocs a plusieurs etats nommes ``x_<id>_<k>`` (et non ``x_<id>``)
_ETATS_NUMEROTES = ("pid", "second_ordre", "fonction_transfert")


def _sans_positions(o: Any) -> Any:
    """Le schema sans les positions ``x``, ``y`` des blocs : deplacer un bloc
    ne change pas l'element, donc ni sa classe ni son empreinte."""
    if isinstance(o, dict):
        bloc = "id" in o and "type" in o
        return {
            k: _sans_positions(v)
            for k, v in o.items()
            if not (bloc and k in ("x", "y"))
        }
    if isinstance(o, list):
        return [_sans_positions(v) for v in o]
    return o


def _fonction_transfert(b: dict, bid: str, u: sp.Expr, valeur, etat) -> list[sp.Expr]:
    """``num(s)/den(s)`` en forme compagne (commandable) : n etats ``x_<id>_k``
    (n = degre du denominateur), ``x_1 = Z, x_2 = s.Z..``, avec ``Z = u/den(s)``.
    ``num`` et ``den`` : listes de coefficients, puissances decroissantes ;
    la fonction doit etre propre (degre du numerateur <= degre du denominateur)."""

    def coeffs(cle: str) -> list[sp.Expr]:
        v = b.get(cle)
        if not isinstance(v, (list, tuple)) or not v:
            raise SchemaInvalide(
                f"fonction_transfert {bid} : {cle} doit etre une liste non vide "
                "de coefficients (puissances decroissantes)"
            )
        return [valeur(b, f"{cle}[{k}]", c) for k, c in enumerate(v)]

    num, den = coeffs("num"), coeffs("den")
    while len(den) > 1 and den[0] == 0:
        den.pop(0)
    while len(num) > 1 and num[0] == 0:
        num.pop(0)
    if den[0] == 0:
        raise SchemaInvalide(f"fonction_transfert {bid} : denominateur nul")
    if len(num) > len(den):
        raise SchemaInvalide(
            f"fonction_transfert {bid} : fonction impropre (degre {len(num) - 1} "
            f"au numerateur, {len(den) - 1} au denominateur)"
        )
    n = len(den) - 1
    if n > ORDRE_MAX:
        raise SchemaInvalide(f"fonction_transfert {bid} : ordre {n} > {ORDRE_MAX}")
    if n == 0:
        return [num[0] / den[0] * u]
    a = [d / den[0] for d in den[1:]]  # a_k, coefficient de s^(n-k)
    bn = [sp.Integer(0)] * (n + 1 - len(num)) + num
    b0 = bn[0] / den[0]
    bk = [c / den[0] for c in bn[1:]]
    z = [sp.Symbol(f"x_{bid}_{k + 1}") for k in range(n)]
    for k in range(n - 1):
        etat(z[k], z[k + 1])
    etat(z[n - 1], u - sum(a[k] * z[n - 1 - k] for k in range(n)))
    return [b0 * u + sum((bk[k] - b0 * a[k]) * z[n - 1 - k] for k in range(n))]


_CLASSES: dict[str, type] = {}


def _newton(F, J, z, p, u, cible) -> tuple[float, np.ndarray]:
    """Newton amorti (moindres carres, recherche lineaire sur la norme du
    residu) sur ``F(z) = cible`` : (residu max, z)."""

    def residu(z):
        r = np.asarray(F(z, p, u), dtype=float).ravel() - cible
        return r, (float(np.linalg.norm(r)) if r.size else 0.0)

    r, n = residu(z)
    for _ in range(100):
        if not np.isfinite(n) or n < 1e-13:
            break
        try:
            Jz = np.asarray(J(z, p, u), dtype=float)
            dz = np.linalg.lstsq(Jz, -r, rcond=None)[0]
        except (np.linalg.LinAlgError, ValueError):
            break
        pas = 1.0
        while pas > 1e-6:
            r2, n2 = residu(z + pas * dz)
            if np.isfinite(n2) and n2 < n:
                break
            pas /= 2
        else:
            break
        z, r, n = z + pas * dz, r2, n2
    m = float(np.abs(r).max()) if r.size else 0.0
    return (m if np.isfinite(m) else np.inf), z


def classe_element(schema: dict) -> type:
    """La classe :class:`~g2egrid.noyau.element.Element` decrite par ``schema``
    (mise en cache par empreinte : deux elements du meme modele partagent
    la meme compilation).

    INITIALISATION -- un Newton amorti sur ``f = 0``, ``g = 0`` et ``i = i_lf``
    dont les inconnues sont les etats, les algebriques et les deux consignes
    CALEES ; les autres consignes sont fixees a leur ``valeur``. Quand le
    schema porte une ``amorce`` (un schema issu de :func:`depuis_element`),
    l'``initialize`` de la classe native sert d'abord de point de depart :
    ses valeurs, reportees par la ``correspondance`` sur les variables du
    schema, puis completees bloc par bloc (chaque sortie intermediaire est
    recalculee de ses entrees), donnent un depart a quelques ulp de la
    solution. Si l'amorce echoue (schema trop modifie, classe native
    absente), le depart d'avant reprend la main."""
    c = compiler(schema)
    cle = f"{c.nom}_{c.empreinte}"
    if cle in _CLASSES:
        return _CLASSES[cle]
    amorce = schema.get("amorce") or None

    noms = c.etats + c.algebriques + c.consignes
    non_calees = list(c.valeurs)
    inconnues_noms = c.etats + c.algebriques + c.calees
    sym = c.symboles()
    residus = c.f + c.g
    entrees = [sym["vd_t"], sym["vq_t"], sym["w_ref"], sym["theta_g"]]
    inconnues = [sym[n] for n in inconnues_noms]
    pars = [sym[n] for n in c.parametres]
    # l'initialisation ne cale que les deux courants injectes
    F = sp.lambdify(
        [inconnues, pars, entrees + [sym[n] for n in non_calees]],
        residus + c.h[:2],
        "numpy",
    )
    J = sp.lambdify(
        [inconnues, pars, entrees + [sym[n] for n in non_calees]],
        sp.Matrix(residus + c.h[:2]).jacobian(inconnues),
        "numpy",
    )
    valeurs_nc = [sp.lambdify(pars, c.valeurs[n], "numpy") for n in non_calees]
    # pour completer le point de depart de l'amorce : chaque sortie de bloc
    # explicite (variable = membre de droite) est recalculee de ses entrees
    explicites = [
        (var, rhs)
        for genre, var, rhs in c.lignes
        if genre == "algebrique" and var in set(c.algebriques)
    ]
    Frhs = (
        sp.lambdify(
            [[sym[n] for n in noms], pars, entrees],
            [r for _v, r in explicites],
            "numpy",
        )
        if amorce and explicites
        else None
    )

    def amorcer(self, ctx: InitContext) -> tuple[dict[str, float], dict[str, float]]:
        """``(connus, consignes recalees)`` depuis l'initialisation native, ou
        ``({}, {})`` si elle n'est pas possible. Les parametres que
        l'initialisation native recalcule (``amorce["recales"]``) sont reportes
        sur l'element, tant qu'ils valent encore leur valeur de depart."""
        from g2egrid.elements.blocs_inverse import initialisation_native

        # Un parametre que l'amorce a deja recale (Kpdc, R et X d'une charge...)
        # redevient LIBRE a l'initialisation suivante, comme dans le natif qui
        # le recalcule a chaque point de fonctionnement : sans cela, un second
        # initialize (continuation, nouveau load flow) gardait la valeur du
        # premier point. Un parametre change par l'utilisateur depuis reste a
        # sa valeur.
        miens = self.__dict__.setdefault("_recales_amorce", {})
        pd = dict(self.pdict)
        for nom, depart in (amorce.get("recales") or {}).items():
            if nom in miens and nom in pd and pd[nom] == miens[nom]:
                pd[nom] = float(depart)
        try:
            init, pnat = initialisation_native(amorce, pd, ctx)
        except Exception:  # noqa: BLE001 - l'amorce n'est qu'un point de depart
            return {}, {}
        recalees: dict[str, float] = {}
        for nom, depart in (amorce.get("recales") or {}).items():
            if nom in c.valeurs:
                fixe = float(
                    valeurs_nc[non_calees.index(nom)](
                        *[self.pdict[k] for k in c.parametres]
                    )
                )
                if fixe == float(depart) and nom in init:
                    recalees[nom] = float(init[nom])
            elif nom in pd and pd[nom] == float(depart):
                v = pnat.get(nom, init.get(nom))
                if v is not None and np.isfinite(v):
                    self.set_param(nom, float(v))
                    miens[nom] = float(v)
        connus = {
            n: float(init[nat])
            for n, nat in (amorce.get("correspondance") or {}).items()
            if nat in init and n in noms
        }
        return connus, recalees

    def completer(connus, p, u, fixes) -> np.ndarray | None:
        """Le point de depart de l'amorce : ``connus`` plus les sorties de
        blocs recalculees, dans l'ordre de ``inconnues_noms``."""
        z = {n: 0.0 for n in noms} | dict(connus) | dict(zip(non_calees, fixes))
        a_calculer = [k for k, (var, _r) in enumerate(explicites) if var not in connus]
        if Frhs is not None and a_calculer:
            tout = np.array([z[n] for n in noms], dtype=float)
            idx = [noms.index(explicites[k][0]) for k in a_calculer]
            with np.errstate(all="ignore"):
                for _ in range(500):
                    r = np.asarray(Frhs(tout, p, u), dtype=float).ravel()
                    nouveau = r[a_calculer]
                    nouveau = np.where(np.isfinite(nouveau), nouveau, 0.0)
                    ecart = np.abs(nouveau - tout[idx]).max()
                    tout[idx] = nouveau
                    if ecart < 1e-15:
                        break
            z = dict(zip(noms, tout))
        depart = np.array([z[n] for n in inconnues_noms], dtype=float)
        return depart if np.all(np.isfinite(depart)) else None

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        V = complex(ctx.v["t"])
        S = complex(ctx.s.get("t", 0j))
        i0 = np.conj(S / V) if abs(V) > 1e-12 else 0j
        cible = np.zeros(len(residus) + 2)
        cible[-2:] = (i0.real, i0.imag)
        connus_nat, recalees = amorcer(self, ctx) if amorce else ({}, {})
        p = [self.pdict[k] for k in c.parametres]
        fixes = [
            recalees.get(n, float(fn(*p))) for n, fn in zip(non_calees, valeurs_nc)
        ]
        u = [V.real, V.imag, ctx.w, ctx.theta_g]
        connus = {"vd_mes": V.real, "vq_mes": V.imag}
        connus |= c.departs
        connus |= {a: float(np.angle(V)) + ctx.theta_g for a in c.angles}
        departs = []
        if connus_nat:
            z0 = completer(connus_nat, p, u, fixes)
            if z0 is not None:
                departs.append(z0)
        for depart in (0.0, 1.0, 0.1, -1.0):
            departs.append(
                np.array([connus.get(n, depart) for n in inconnues_noms], dtype=float)
            )
        meilleur = None
        with np.errstate(all="ignore"):
            for z in departs:
                n = _newton(F, J, z, p, u + fixes, cible)
                if meilleur is None or n[0] < meilleur[0]:
                    meilleur = n
                if n[0] < 1e-9:
                    break
        if meilleur[0] > 1e-6:
            raise RuntimeError(
                f"{self.name}: le schema {c.nom} n'atteint pas le point du load flow "
                f"(residu {meilleur[0]:.2e}) : verifier que ses deux consignes "
                "commandent bien P et Q"
            )
        out = {n: float(v) for n, v in zip(inconnues_noms, meilleur[1])}
        return out | dict(zip(non_calees, fixes))

    # les gardes des limiteurs, comme LimiteurCourant : le module de la
    # consigne LIBRE (les entrees du limiteur) face a Imax, avec hysteresis
    args = [
        [sym[n] for n in c.etats],
        [sym[n] for n in c.algebriques],
        [sp.Symbol(n) for n in (*c.consignes, "w_ref", "theta_g", "vd_t", "vq_t")],
        pars,
    ]
    surveilles = [
        (d, sp.lambdify(args, mod, "numpy"), sp.lambdify(args, imax, "numpy"))
        for d, mod, imax in c.limiteurs
    ]
    # les saturations : l'entree face a ses deux bornes
    saturees = [
        (
            d,
            sp.lambdify(args, ent, "numpy"),
            sp.lambdify(args, lo, "numpy"),
            sp.lambdify(args, hi, "numpy"),
        )
        for d, ent, lo, hi in c.saturations
    ]

    def gardes(self, cm, x, z, u, p) -> list[tuple[str, float, float, int]]:
        if self.bloque:
            return []
        out = []
        for d, module, imax in surveilles:
            courant, Imax = float(module(x, z, u, p)), float(imax(x, z, u, p))
            if Imax <= 0.0:
                continue  # limite non reglee
            if p[self.params.index(d)] == 0.0:
                out.append((d, courant - Imax, 1.0, +1))
            else:
                out.append((d, courant - Imax * (1 - HYSTERESIS), 0.0, -1))
        for d, entree, bas, haut in saturees:
            e, lo, hi = (
                float(entree(x, z, u, p)),
                float(bas(x, z, u, p)),
                float(haut(x, z, u, p)),
            )
            marge = HYSTERESIS * abs(hi - lo)
            if p[self.params.index(d)] == 0.0:
                # libre : on sature quand l'entree franchit une borne
                out.append((d, e - hi, 1.0, +1))
                out.append((d, lo - e, 1.0, +1))
            else:
                # sature : on relache quand elle est revenue dans la bande
                out.append((d, e - (hi - marge), 0.0, -1))
                out.append((d, (lo + marge) - e, 0.0, -1))
        return out

    attrs = dict(
        model_name=f"Blocs_{cle}",
        modes=tuple(d for d, _m, _i in c.limiteurs)
        + tuple(d for d, *_r in c.saturations),
        gardes=gardes,
        bloque=False,
        ports=("t",),
        params=tuple(c.parametres),
        states=tuple(c.etats),
        algebs=tuple(c.algebriques),
        svc_inputs=tuple(c.consignes),
        svc_outputs=tuple(c.services),
        compilation=c,
        equations=staticmethod(c.equations),
        initialize=initialize,
        __doc__=f"Element decrit par le schema-bloc {c.nom!r}.",
    )
    if amorce:
        from g2egrid.elements.blocs_inverse import classe_native

        natif = classe_native(amorce.get("classe"))
        # un element passif suit la pulsation du repere (voir Network.repere)
        attrs["passif"] = bool(getattr(natif, "passif", False))
    classe = type(f"Blocs_{c.nom}", (Element,), attrs)
    _CLASSES[cle] = classe
    return classe


# ------------------------------------------------------------------ exemples
def _lien(de: str, vers: str, entree: int = 0, sortie: int = 0) -> dict:
    lien = {"de": de, "vers": vers, "entree": entree}
    if sortie:
        lien["sortie"] = sortie
    return lien


#: les parametres d'un GFL tels que le moteur les regle (``GFL1``, 50 Hz)
_CLES_GFL = {
    "R1": "R1_pu",
    "L1": "L1_pu",
    "C": "C_pu",
    "R2": "R2_pu",
    "L2": "L2_pu",
    "KpCL": "KpCL",
    "KiCL": "KiCL",
    "Kffv": "Kffv",
    "wff": "wff",
    "Kppll": "Kppll",
    "Kipll": "Kipll",
    "Kpp": "Kpp",
    "Kip": "Kip",
    "Kpq": "Kpq",
    "Kiq": "Kiq",
}


def _inserer_limiteur(
    blocs: list[dict],
    liens: list[dict],
    amont: tuple[str, str],
    pi: tuple[str, str],
    aval: tuple[str, str],
    x: float,
    y: float,
    bid: str = "lim",
    genre: str = "limiteur_courant",
    champs: dict | None = None,
) -> None:
    """Pose le limiteur de courant ``lim`` entre les blocs ``amont`` qui
    font la consigne de courant (d, q) et les sommateurs ``aval`` de la boucle
    de courant ; les PI de la boucle externe (``pi``) deviennent des
    ``pi_aw`` qui recoivent le retro-calcul. Les blocs a droite de ``x`` se decalent pour lui faire
    place. C'est le limiteur des GFL et des GFM du moteur
    (:func:`g2egrid.elements.limiteurs.limiter_courant`), au parametre pres :
    son drapeau s'appelle ``mode_lim``, comme le leur, si bien que les
    manoeuvres de la simulation (``blocage``) le trouvent. La meme pose sert
    a la saturation de tension (``genre="saturation_tension"``, ``bid="vsat"``),
    entre la boucle de courant et le filtre."""
    for b in blocs:
        if b["x"] >= x:
            b["x"] += 160
        if b["id"] in pi:
            b["type"] = "pi_aw"
    champs = champs or {"Imax": "Imax", "Taw": "Taw"}
    blocs.append({"id": bid, "type": genre, **champs, "x": x, "y": y})
    for lien in liens:
        for k, (a, v) in enumerate(zip(amont, aval)):
            if lien["de"] == a and lien["vers"] == v:
                lien["de"] = bid
                if k:
                    lien["sortie"] = k
    for k, (a, c) in enumerate(zip(amont, pi)):
        liens.append(_lien(a, bid, k))
        liens.append(_lien(bid, c, 1, 2 + k))


def schema_gfl(nom: str = "gfl_complet", **reglages) -> dict:
    """Un GFL COMPLET en schema-bloc : filtre LCL (``rl_dq``, ``c_dq``,
    ``rl_dq``), PLL (``park`` -> PI sur ``v_q`` -> integrateur), boucles de
    puissance (PI sur P et Q) et de courant (PI, anticipation de la tension et
    decouplage), retour au repere global (``park_inverse``) et mise a
    l'echelle ``kS``.

    La LIMITATION DE COURANT (bloc ``lim``) ramene la consigne de courant
    dans ``|i| <= Imax`` (pu du convertisseur) ; les PI de puissance suivent
    par retro-calcul (``Taw``). Comme dans le moteur, c'est un mode : libre
    (``mode_lim = 0``) au point de fonctionnement, bascule par la simulation
    temporelle quand la consigne depasse ``Imax``.

    La SATURATION DE TENSION (bloc ``vsat``) ramene la modulation de la boucle
    de courant dans ``|m| <= Mmax`` (1,15 ; ``Mmax = 0`` la desactive, comme
    dans le moteur ou c'est le defaut), et les PI de courant suivent par
    retro-calcul (``Taw_v = tau_CL/2``). Sur un defaut, elle s'enclenche au
    retour de la tension, quand l'onduleur demande plus qu'il ne sait faire.

    Ce sont exactement les equations de :class:`~g2egrid.elements.GFL1` (bus
    continu ideal, ``vdc = 1``), limitation et saturation comprises : les
    valeurs propres d'un reseau sont les memes avec l'un ou l'autre, et une
    simulation de defaut aussi (le moteur regle a ``Mmax = 1,15``). Les
    parametres sont ceux du moteur, reglage automatique compris (``GFL_DEFAULT`` a
    50 Hz), et ``reglages`` en remplace (memes noms que le schema : ``R1``,
    ``L1``, ``KpCL``, ``Kppll``...)."""
    from g2egrid.elements.gfl import GFL1

    p = GFL1("_", 1, Sn=1.0, Sb=1.0, wb=PARAMS_RESERVES["wb"]).pdict
    parametres = {k: float(p[v]) for k, v in _CLES_GFL.items()}
    parametres.update(Imax=float(p["Imax"]), Taw=float(p["Taw"]))
    # saturation de tension reglee (1,15, modulation vectorielle ; le moteur
    # la laisse a 0), retro-calcul du moteur (tau_CL/2)
    parametres.update(Mmax=MMAX, Taw_v=float(p["Taw_v"]))
    parametres.update({k: float(v) for k, v in reglages.items()})
    blocs = [
        # le reseau et les consignes
        {"id": "vd", "type": "entree", "signal": "vd", "x": 20, "y": 30},
        {"id": "vq", "type": "entree", "signal": "vq", "x": 20, "y": 100},
        {"id": "w", "type": "entree", "signal": "w", "x": 600, "y": 800},
        {"id": "pr", "type": "consigne", "nom": "P_ref", "x": 20, "y": 300},
        {"id": "qr", "type": "consigne", "nom": "Q_ref", "x": 20, "y": 480},
        # repere de la PLL
        {"id": "pk", "type": "park", "x": 170, "y": 40},
        # mesure des puissances (repere PLL) : p = ve.ig, q = veq.igd - ved.igq
        {
            "id": "pe",
            "type": "fonction",
            "expr": "u1*u3 + u2*u4",
            "n": 4,
            "x": 170,
            "y": 230,
        },
        {
            "id": "qe",
            "type": "fonction",
            "expr": "-u1*u4 + u2*u3",
            "n": 4,
            "x": 170,
            "y": 430,
        },
        # boucles de puissance
        {"id": "ep", "type": "somme", "signes": "+-", "x": 330, "y": 300},
        {"id": "cp", "type": "pi", "Kp": "Kpp", "Ki": "Kip", "x": 430, "y": 300},
        {"id": "eqp", "type": "somme", "signes": "+-", "x": 330, "y": 470},
        {"id": "cq", "type": "pi", "Kp": "Kpq", "Ki": "Kiq", "x": 430, "y": 470},
        # boucles de courant : PI, anticipation de ve, decouplage w.L1.i
        {"id": "ed", "type": "somme", "signes": "+-", "x": 590, "y": 300},
        {"id": "cd", "type": "pi", "Kp": "KpCL", "Ki": "KiCL", "x": 690, "y": 300},
        {
            "id": "md",
            "type": "fonction",
            "expr": "u1 + Kffv*u2 - wff*L1*u3",
            "n": 3,
            "x": 840,
            "y": 290,
        },
        {"id": "eqc", "type": "somme", "signes": "+-", "x": 590, "y": 470},
        {"id": "cqc", "type": "pi", "Kp": "KpCL", "Ki": "KiCL", "x": 690, "y": 470},
        {
            "id": "mq",
            "type": "fonction",
            "expr": "u1 + Kffv*u2 + wff*L1*u3",
            "n": 3,
            "x": 840,
            "y": 460,
        },
        # filtre LCL, dans le repere de la PLL
        {"id": "rl1", "type": "rl_dq", "R": "R1", "L": "L1", "x": 1000, "y": 350},
        {"id": "cf", "type": "c_dq", "C": "C", "x": 1150, "y": 350},
        {"id": "rl2", "type": "rl_dq", "R": "R2", "L": "L2", "x": 1300, "y": 350},
        # PLL : PI sur veq, + w_ff, integrateur de wb.(w_pll - w)
        {"id": "pll", "type": "pi", "Kp": "Kppll", "Ki": "Kipll", "x": 440, "y": 650},
        {"id": "wf", "type": "constante", "valeur": "wff", "x": 440, "y": 740},
        {"id": "wp", "type": "somme", "signes": "++", "x": 600, "y": 680},
        {"id": "dw", "type": "somme", "signes": "+-", "x": 750, "y": 740},
        {"id": "th", "type": "integrateur", "K": "wb", "x": 860, "y": 740},
        # retour au repere global, a la base du reseau
        {"id": "pinv", "type": "park_inverse", "x": 1450, "y": 120},
        {"id": "kd", "type": "gain", "K": "kS", "x": 1600, "y": 90},
        {"id": "kq", "type": "gain", "K": "kS", "x": 1600, "y": 180},
        {"id": "od", "type": "sortie", "signal": "id", "x": 1750, "y": 90},
        {"id": "oq", "type": "sortie", "signal": "iq", "x": 1750, "y": 180},
    ]
    liens = [
        _lien("vd", "pk", 0),
        _lien("vq", "pk", 1),
        _lien("th", "pk", 2),
        # mesures
        _lien("cf", "pe", 0, 0),
        _lien("cf", "pe", 1, 1),
        _lien("rl2", "pe", 2, 0),
        _lien("rl2", "pe", 3, 1),
        _lien("cf", "qe", 0, 0),
        _lien("cf", "qe", 1, 1),
        _lien("rl2", "qe", 2, 0),
        _lien("rl2", "qe", 3, 1),
        # puissance
        _lien("pr", "ep", 0),
        _lien("pe", "ep", 1),
        _lien("ep", "cp"),
        _lien("qe", "eqp", 0),
        _lien("qr", "eqp", 1),
        _lien("eqp", "cq"),
        # courant
        _lien("cp", "ed", 0),
        _lien("rl1", "ed", 1, 0),
        _lien("ed", "cd"),
        _lien("cd", "md", 0),
        _lien("cf", "md", 1, 0),
        _lien("rl1", "md", 2, 1),
        _lien("cq", "eqc", 0),
        _lien("rl1", "eqc", 1, 1),
        _lien("eqc", "cqc"),
        _lien("cqc", "mq", 0),
        _lien("cf", "mq", 1, 1),
        _lien("rl1", "mq", 2, 0),
        # filtre LCL
        _lien("md", "rl1", 0),
        _lien("mq", "rl1", 1),
        _lien("cf", "rl1", 2, 0),
        _lien("cf", "rl1", 3, 1),
        _lien("wp", "rl1", 4),
        _lien("rl1", "cf", 0, 0),
        _lien("rl1", "cf", 1, 1),
        _lien("rl2", "cf", 2, 0),
        _lien("rl2", "cf", 3, 1),
        _lien("wp", "cf", 4),
        _lien("cf", "rl2", 0, 0),
        _lien("cf", "rl2", 1, 1),
        _lien("pk", "rl2", 2, 0),
        _lien("pk", "rl2", 3, 1),
        _lien("wp", "rl2", 4),
        # PLL
        _lien("cf", "pll", 0, 1),
        _lien("pll", "wp", 0),
        _lien("wf", "wp", 1),
        _lien("wp", "dw", 0),
        _lien("w", "dw", 1),
        _lien("dw", "th"),
        # sortie
        _lien("rl2", "pinv", 0, 0),
        _lien("rl2", "pinv", 1, 1),
        _lien("th", "pinv", 2),
        _lien("pinv", "kd", 0, 0),
        _lien("pinv", "kq", 0, 1),
        _lien("kd", "od"),
        _lien("kq", "oq"),
    ]
    _inserer_limiteur(blocs, liens, ("cp", "cq"), ("cp", "cq"), ("ed", "eqc"), 580, 350)
    # saturation de tension : entre la boucle de courant (md, mq) et le filtre,
    # retro-calcul sur les PI de COURANT -- comme le GFM
    x_rl1 = next(b["x"] for b in blocs if b["id"] == "rl1")
    _inserer_limiteur(
        blocs,
        liens,
        ("md", "mq"),
        ("cd", "cqc"),
        ("rl1", "rl1"),
        x_rl1,
        520,
        bid="vsat",
        genre="saturation_tension",
        champs={"Mmax": "Mmax", "Taw": "Taw_v"},
    )
    return {"nom": nom, "parametres": parametres, "blocs": blocs, "liens": liens}


def schema_gfm(nom: str = "gfm_complet", **reglages) -> dict:
    """Un GFM COMPLET en schema-bloc, a machine synchrone virtuelle (VSM) :

    * la VSM : ``2H.d(dw)/dt = P_ref - p_m - D.dw`` sur la puissance mesuree
      (premier ordre ``wf_p``), puis ``d(theta)/dt = wb.dw`` ; la frequence
      du convertisseur est ``w + dw`` ;
    * le statisme Q-V : ``v_ref = E_ref - nq.q_m`` (``q_m`` filtree a ``wf_q``) ;
    * l'IMPEDANCE VIRTUELLE (``Rv``, ``Lv``, nulle par defaut) : la chute
      ``Rv.i_g + Lv.(d/dt + j.w).i_g`` est retranchee des deux consignes de
      tension ; dynamique (``zv_dyn = 1``, la chute sur ``L2`` mise a
      l'echelle ``Lv/L2``) ou quasi statique (``zv_dyn = 0``) ;
    * la boucle de TENSION : PI sur ``v_ref - v_ed`` et ``0 - v_eq``, plus
      l'anticipation du courant de sortie (``Kffi.i_g``) et le decouplage du
      condensateur (``-+ wff.C.v_e``) : ce sont les consignes de courant ;
    * la boucle de COURANT, le filtre LCL dans le repere de la VSM, le retour
      au repere global et la mise a l'echelle ``kS``, comme pour le GFL.

    La LIMITATION DE COURANT (bloc ``lim``) se place entre la boucle de
    tension et la boucle de courant, comme dans le moteur : la consigne de
    courant est ramenee dans ``|i| <= Imax``, et les PI de TENSION suivent par
    retro-calcul. Libre au point de fonctionnement, elle est basculee par la
    simulation temporelle.

    La SATURATION DE TENSION (bloc ``vsat``) ramene la modulation de la
    boucle de courant dans ``|m| <= Mmax`` (1,15 : modulation vectorielle ;
    ``Mmax = 0`` la desactive, comme dans le moteur ou c'est le defaut), et
    les PI de COURANT suivent par retro-calcul (``Taw_v = tau_CL/2``, rapide,
    le defaut du moteur). Sur un defaut franc, elle s'enclenche avant la
    limitation de courant, quelques millisecondes ; elle demande alors un pas
    d'integration fin (quelques dizaines de us en BDF2).

    Ce sont exactement les equations de :class:`~g2egrid.elements.GFM2`
    (statisme du 2e ordre, bus continu ideal), limitation de courant et
    saturation de tension comprises : ``H`` et
    ``D`` y sont ``1/(2.mp.wc)`` et ``1/mp``, et les deux consignes du GFM
    (``q_ref`` et ``ve_ref``) se fondent en ``E_ref = ve_ref + nq.q_ref`` --
    le load flow ne peut caler que deux consignes, et la dynamique n'en voit
    que la somme. Les parametres par defaut sont ceux du moteur, reglage
    automatique des PI compris (``GFM_DEFAULT`` a 50 Hz), et ``reglages`` en
    remplace (noms du schema : ``H``, ``D``, ``nq``, ``KpVL``...)."""
    from g2egrid.elements.gfm import GFM2

    p = GFM2("_", 1, Sn=1.0, Sb=1.0, wb=PARAMS_RESERVES["wb"]).pdict
    parametres = {
        "R1": p["R1_pu"],
        "L1": p["L1_pu"],
        "C": p["C_pu"],
        "R2": p["R2_pu"],
        "L2": p["L2_pu"],
        "H": 1.0 / (2.0 * p["mp"] * p["wc"]),
        "D": 1.0 / p["mp"],
        "nq": p["nq"],
        "wf_p": p["wf_p"],
        "wf_q": p["wf_q"],
        "KpVL": p["KpVL"],
        "KiVL": p["KiVL"],
        "Kffi": p["Kffi"],
        "KpCL": p["KpCL"],
        "KiCL": p["KiCL"],
        "Kffv": p["Kffv"],
        "wff": p["wff"],
        "Imax": p["Imax"],
        "Taw": p["Taw"],
        # la saturation de tension est REGLEE (le moteur la laisse a Mmax = 0,
        # inactive) : 1,15, la modulation vectorielle ; Mmax = 0 la desactive.
        # Son retro-calcul est celui du moteur, RAPIDE (tau_CL/2) : voir
        # GFM_DEFAULT pour ce que fait un retro-calcul lent.
        "Mmax": MMAX,
        "Taw_v": p["Taw_v"],
        # impedance virtuelle, nulle par defaut (cf. GFM_DEFAULT)
        "Rv": p["Rv_pu"],
        "Lv": p["Lv_pu"],
        "zv_dyn": p["zv_dyn"],
    }
    parametres = {k: float(v) for k, v in parametres.items()}
    parametres.update({k: float(v) for k, v in reglages.items()})

    def f(bid, expr, n, x, y):
        return {"id": bid, "type": "fonction", "expr": expr, "n": n, "x": x, "y": y}

    blocs = [
        # le reseau et les consignes
        {"id": "vd", "type": "entree", "signal": "vd", "x": 20, "y": 30},
        {"id": "vq", "type": "entree", "signal": "vq", "x": 20, "y": 100},
        {"id": "w", "type": "entree", "signal": "w", "x": 560, "y": 380},
        {"id": "pr", "type": "consigne", "nom": "P_ref", "x": 20, "y": 250},
        {"id": "er", "type": "consigne", "nom": "E_ref", "x": 20, "y": 520},
        # repere de la VSM
        {"id": "pk", "type": "park", "x": 170, "y": 40},
        # puissances mesurees, filtrees
        f("pe", "u1*u3 + u2*u4", 4, 170, 190),
        f("qe", "-u1*u4 + u2*u3", 4, 170, 400),
        {
            "id": "pm",
            "type": "premier_ordre",
            "K": 1,
            "tau": "1/wf_p",
            "x": 310,
            "y": 200,
        },
        {
            "id": "qm",
            "type": "premier_ordre",
            "K": 1,
            "tau": "1/wf_q",
            "x": 310,
            "y": 410,
        },
        # VSM : 2H.d(dw)/dt = P_ref - p_m - D.dw ; d(theta)/dt = wb.dw
        {"id": "ep", "type": "somme", "signes": "+--", "x": 460, "y": 230},
        {"id": "inertie", "type": "gain", "K": "1/(2*H)", "x": 550, "y": 230},
        {"id": "dw", "type": "integrateur", "K": 1, "x": 690, "y": 230},
        {"id": "amort", "type": "gain", "K": "D", "x": 550, "y": 310},
        {"id": "wv", "type": "somme", "signes": "++", "x": 840, "y": 330},
        {"id": "th", "type": "integrateur", "K": "wb", "x": 840, "y": 150},
        # statisme Q-V
        f("vr", "u1 - nq*u2", 2, 460, 500),
        # impedance virtuelle : chute Rv.i_g + Lv.(d/dt + j.w).i_g, dynamique
        # (la chute sur L2 a l'echelle Lv/L2) ou quasi statique (j.w.Lv.i_g)
        f(
            "zvd",
            "Rv*u1 - (1-zv_dyn)*Lv*u5*u2 + zv_dyn*Lv/L2*(u3 - u4 - R2*u1)",
            5,
            460,
            690,
        ),
        f(
            "zvq",
            "Rv*u1 + (1-zv_dyn)*Lv*u5*u2 + zv_dyn*Lv/L2*(u3 - u4 - R2*u1)",
            5,
            460,
            800,
        ),
        # boucle de tension : PI, anticipation de i_g, decouplage de C
        {"id": "evd", "type": "somme", "signes": "+--", "x": 600, "y": 500},
        {"id": "cvd", "type": "pi", "Kp": "KpVL", "Ki": "KiVL", "x": 690, "y": 500},
        f("isdr", "u1 + Kffi*u2 - wff*C*u3", 3, 840, 490),
        {"id": "evq", "type": "somme", "signes": "--", "x": 600, "y": 610},
        {"id": "cvq", "type": "pi", "Kp": "KpVL", "Ki": "KiVL", "x": 690, "y": 610},
        f("isqr", "u1 + Kffi*u2 + wff*C*u3", 3, 840, 600),
        # boucle de courant : PI, anticipation de v_e, decouplage de L1
        {"id": "ed", "type": "somme", "signes": "+-", "x": 1000, "y": 500},
        {"id": "cd", "type": "pi", "Kp": "KpCL", "Ki": "KiCL", "x": 1090, "y": 500},
        f("md", "u1 + Kffv*u2 - wff*L1*u3", 3, 1240, 490),
        {"id": "eqc", "type": "somme", "signes": "+-", "x": 1000, "y": 610},
        {"id": "cqc", "type": "pi", "Kp": "KpCL", "Ki": "KiCL", "x": 1090, "y": 610},
        f("mq", "u1 + Kffv*u2 + wff*L1*u3", 3, 1240, 600),
        # filtre LCL, dans le repere de la VSM, a sa frequence
        {"id": "rl1", "type": "rl_dq", "R": "R1", "L": "L1", "x": 1400, "y": 320},
        {"id": "cf", "type": "c_dq", "C": "C", "x": 1560, "y": 320},
        {"id": "rl2", "type": "rl_dq", "R": "R2", "L": "L2", "x": 1720, "y": 320},
        # retour au repere global, a la base du reseau
        {"id": "pinv", "type": "park_inverse", "x": 1880, "y": 120},
        {"id": "kd", "type": "gain", "K": "kS", "x": 2050, "y": 90},
        {"id": "kq", "type": "gain", "K": "kS", "x": 2050, "y": 180},
        {"id": "od", "type": "sortie", "signal": "id", "x": 2200, "y": 90},
        {"id": "oq", "type": "sortie", "signal": "iq", "x": 2200, "y": 180},
        # l'angle et la vitesse de la VSM : de quoi porter le repere du reseau
        {"id": "oth", "type": "sortie", "signal": "theta", "x": 1000, "y": 150},
        {"id": "ow", "type": "sortie", "signal": "w", "x": 1000, "y": 330},
    ]
    liens = [
        _lien("vd", "pk", 0),
        _lien("vq", "pk", 1),
        _lien("th", "pk", 2),
        _lien("th", "oth"),
        _lien("wv", "ow"),
        # mesures
        *(_lien("cf", m, k, k) for m in ("pe", "qe") for k in (0, 1)),
        *(_lien("rl2", m, k + 2, k) for m in ("pe", "qe") for k in (0, 1)),
        _lien("pe", "pm"),
        _lien("qe", "qm"),
        # VSM
        _lien("pr", "ep", 0),
        _lien("pm", "ep", 1),
        _lien("amort", "ep", 2),
        _lien("ep", "inertie"),
        _lien("inertie", "dw"),
        _lien("dw", "amort"),
        _lien("dw", "th"),
        _lien("dw", "wv", 0),
        _lien("w", "wv", 1),
        # statisme Q-V et boucle de tension
        _lien("er", "vr", 0),
        _lien("qm", "vr", 1),
        _lien("vr", "evd", 0),
        _lien("cf", "evd", 1, 0),
        _lien("evd", "cvd"),
        _lien("cvd", "isdr", 0),
        _lien("rl2", "isdr", 1, 0),
        _lien("cf", "isdr", 2, 1),
        _lien("cf", "evq", 0, 1),
        # impedance virtuelle : i_g, v_e, v_g (repere de la VSM) et la vitesse
        _lien("rl2", "zvd", 0, 0),
        _lien("rl2", "zvd", 1, 1),
        _lien("cf", "zvd", 2, 0),
        _lien("pk", "zvd", 3, 0),
        _lien("wv", "zvd", 4),
        _lien("rl2", "zvq", 0, 1),
        _lien("rl2", "zvq", 1, 0),
        _lien("cf", "zvq", 2, 1),
        _lien("pk", "zvq", 3, 1),
        _lien("wv", "zvq", 4),
        _lien("zvd", "evd", 2),
        _lien("zvq", "evq", 1),
        _lien("evq", "cvq"),
        _lien("cvq", "isqr", 0),
        _lien("rl2", "isqr", 1, 1),
        _lien("cf", "isqr", 2, 0),
        # boucle de courant
        _lien("isdr", "ed", 0),
        _lien("rl1", "ed", 1, 0),
        _lien("ed", "cd"),
        _lien("cd", "md", 0),
        _lien("cf", "md", 1, 0),
        _lien("rl1", "md", 2, 1),
        _lien("isqr", "eqc", 0),
        _lien("rl1", "eqc", 1, 1),
        _lien("eqc", "cqc"),
        _lien("cqc", "mq", 0),
        _lien("cf", "mq", 1, 1),
        _lien("rl1", "mq", 2, 0),
        # filtre LCL
        _lien("md", "rl1", 0),
        _lien("mq", "rl1", 1),
        _lien("cf", "rl1", 2, 0),
        _lien("cf", "rl1", 3, 1),
        _lien("wv", "rl1", 4),
        _lien("rl1", "cf", 0, 0),
        _lien("rl1", "cf", 1, 1),
        _lien("rl2", "cf", 2, 0),
        _lien("rl2", "cf", 3, 1),
        _lien("wv", "cf", 4),
        _lien("cf", "rl2", 0, 0),
        _lien("cf", "rl2", 1, 1),
        _lien("pk", "rl2", 2, 0),
        _lien("pk", "rl2", 3, 1),
        _lien("wv", "rl2", 4),
        # sortie
        _lien("rl2", "pinv", 0, 0),
        _lien("rl2", "pinv", 1, 1),
        _lien("th", "pinv", 2),
        _lien("pinv", "kd", 0, 0),
        _lien("pinv", "kq", 0, 1),
        _lien("kd", "od"),
        _lien("kq", "oq"),
    ]
    _inserer_limiteur(
        blocs, liens, ("isdr", "isqr"), ("cvd", "cvq"), ("ed", "eqc"), 990, 520
    )
    # saturation de tension : entre la boucle de courant (md, mq) et le filtre,
    # retro-calcul sur les PI de COURANT
    x_rl1 = next(b["x"] for b in blocs if b["id"] == "rl1")
    _inserer_limiteur(
        blocs,
        liens,
        ("md", "mq"),
        ("cd", "cqc"),
        ("rl1", "rl1"),
        x_rl1,
        520,
        bid="vsat",
        genre="saturation_tension",
        champs={"Mmax": "Mmax", "Taw": "Taw_v"},
    )
    return {"nom": nom, "parametres": parametres, "blocs": blocs, "liens": liens}


_EXEMPLES: dict[str, dict[str, Any]] = {
    "source_filtree": {
        "titre": "Source P/Q filtree (1er ordre)",
        "schema": {
            "nom": "source_filtree",
            "parametres": {"tau": 0.05},
            "blocs": [
                {"id": "pr", "type": "consigne", "nom": "P_ref", "x": 40, "y": 60},
                {"id": "qr", "type": "consigne", "nom": "Q_ref", "x": 40, "y": 200},
                {
                    "id": "fp",
                    "type": "premier_ordre",
                    "K": 1,
                    "tau": "tau",
                    "x": 230,
                    "y": 60,
                },
                {
                    "id": "fq",
                    "type": "premier_ordre",
                    "K": 1,
                    "tau": "tau",
                    "x": 230,
                    "y": 200,
                },
                {"id": "op", "type": "sortie", "signal": "P", "x": 430, "y": 60},
                {"id": "oq", "type": "sortie", "signal": "Q", "x": 430, "y": 200},
            ],
            "liens": [
                _lien("pr", "fp"),
                _lien("qr", "fq"),
                _lien("fp", "op"),
                _lien("fq", "oq"),
            ],
        },
    },
    "asservissement": {
        "titre": "Asservissement PI de P (boucle fermee)",
        "schema": {
            "nom": "source_pi",
            "parametres": {"Kp": 0.5, "Ki": 20, "tm": 0.02},
            "blocs": [
                {"id": "pr", "type": "consigne", "nom": "P_ref", "x": 30, "y": 60},
                {"id": "e", "type": "somme", "signes": "+-", "x": 180, "y": 60},
                {"id": "c", "type": "pi", "Kp": "Kp", "Ki": "Ki", "x": 330, "y": 60},
                {
                    "id": "m",
                    "type": "premier_ordre",
                    "K": 1,
                    "tau": "tm",
                    "x": 330,
                    "y": 190,
                },
                {"id": "op", "type": "sortie", "signal": "P", "x": 530, "y": 60},
                {"id": "qr", "type": "consigne", "nom": "Q_ref", "x": 30, "y": 320},
                {"id": "oq", "type": "sortie", "signal": "Q", "x": 530, "y": 320},
            ],
            "liens": [
                _lien("pr", "e", 0),
                _lien("m", "e", 1),
                _lien("e", "c"),
                _lien("c", "m"),
                _lien("c", "op"),
                _lien("qr", "oq"),
            ],
        },
    },
    "statisme_qv": {
        "titre": "Source a statisme Q(V) mesure",
        "schema": {
            "nom": "statisme_qv",
            "parametres": {"kq": 5, "tv": 0.02},
            "blocs": [
                {"id": "pr", "type": "consigne", "nom": "P_ref", "x": 40, "y": 50},
                {"id": "op", "type": "sortie", "signal": "P", "x": 560, "y": 50},
                {"id": "vr", "type": "consigne", "nom": "V_ref", "x": 40, "y": 170},
                {"id": "v", "type": "entree", "signal": "V", "x": 40, "y": 290},
                {
                    "id": "fv",
                    "type": "premier_ordre",
                    "K": 1,
                    "tau": "tv",
                    "x": 180,
                    "y": 290,
                },
                {"id": "e", "type": "somme", "signes": "+-", "x": 320, "y": 170},
                {"id": "k", "type": "gain", "K": "kq", "x": 440, "y": 170},
                {"id": "oq", "type": "sortie", "signal": "Q", "x": 560, "y": 170},
            ],
            "liens": [
                _lien("pr", "op"),
                _lien("vr", "e", 0),
                _lien("v", "fv"),
                _lien("fv", "e", 1),
                _lien("e", "k"),
                _lien("k", "oq"),
            ],
        },
    },
}


def exemples() -> dict[str, dict[str, Any]]:
    """Les schemas proposes par l'editeur : ``{cle: {titre, schema}}``
    (copies : les modifier ne touche pas aux exemples)."""
    out = copy.deepcopy(_EXEMPLES)
    out["gfl_complet"] = {
        "titre": "GFL complet : LCL, PLL, boucles de courant et de puissance",
        "schema": schema_gfl(),
    }
    out["gfm_complet"] = {
        "titre": "GFM complet : VSM, statisme Q-V, boucles de tension et de courant",
        "schema": schema_gfm(),
    }
    return out


# les fonctions de l'inverse (equations -> blocs) : import en queue de module,
# ``blocs_inverse`` ne charge ``blocs`` que dans ses fonctions
from g2egrid.elements.blocs_inverse import depuis_element, disposer
