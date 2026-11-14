"""Le MODELE A EQUATIONS : un element decrit par ses equations ecrites, en
texte, nom par nom (FR-052).

Un modele a equations est le pendant TEXTE d'un schema-bloc
(:mod:`g2egrid.elements.blocs`) : meme compilation (:class:`Compilation`),
meme classe d'element (:func:`~g2egrid.elements.blocs.classe_element`), meme
pose dans ``T_Blocs`` ; on le reconnait a ``"genre": "equations"``::

    {
      "genre": "equations",
      "nom": "GFL1_GFL3",
      "parametres": {"Kipll": 1200.0, ...},
      "consignes": [{"nom": "p_ref", "calee": true},
                    {"nom": "q_ref", "calee": true}],
      "etats": [{"nom": "theta_pll", "derivee": "wb*(w_pll - w)"}, ...],
      "algebriques": [{"nom": "w_pll", "expr": "1 + Kppll*veq + M_pll"},
                      {"nom": "z", "residu": "...", "implicite": true}, ...],
      "sorties": {"id": "...", "iq": "...", "theta": "...", "w": "..."},
      "amorce": {"classe": "GFL1", "correspondance": {...}, "recales": {...}}
    }

LES NOMS -- une expression s'ecrit avec les noms des etats, des variables
algebriques, des consignes, des parametres, les reserves ``wb`` et ``kS``, et
les entrees du reseau ``vd``, ``vq``, ``V`` (le module de la tension), ``w``
(la vitesse du repere) et ``theta_g``. Le texte passe par l'analyseur SUR de
:mod:`g2egrid.elements.blocs` (liste blanche : arithmetique, fonctions usuelles,
aucun appel, attribut ni indice) : rien d'autre n'est jamais evalue.

``etats`` : ``d<nom>/dt = derivee``. ``algebriques`` : ``nom = expr`` quand la
variable s'explicite ; ``0 = residu`` (``"implicite": true``) sinon, la
variable etant alors determinee par le residu. ``sorties`` : le courant
injecte ``id`` et ``iq`` (pu base reseau), ou la puissance ``P`` et ``Q``
comme pour un schema-bloc ; ``theta`` et ``w`` facultatives font un FORMEUR.
Un etat ou une algebrique peut porter un ``depart`` : la valeur d'ou part le
Newton de l'initialisation.

LES CONVERSIONS -- :func:`equations_depuis_element` et
:func:`equations_depuis_instance` ecrivent un element natif LISIBLEMENT : les
variables algebriques explicites restent des signaux NOMMES (jamais
developpees dans les derivees), les parametres gardent leur nom, une
expression a chaque etat et chaque algebrique de la classe, et le modele porte
l'``amorce`` qui fait partir l'initialisation de celle du natif.
:func:`equations_depuis_schema` ecrit de meme un schema-bloc (les sorties
``y_<bloc>`` purement intermediaires y sont substituees tant que
l'expression reste lisible) ; :func:`schema_depuis_equations` fait le chemin
inverse, vers l'editeur de blocs.

UNE LIMITE -- un modele a equations decrit le fonctionnement LIBRE de
l'element, celui sur lequel on linearise : les drapeaux de limiteur ou de
saturation d'un schema-bloc y deviennent de simples parametres (a 0), et les
gardes de la simulation temporelle ne sont pas reprises.
"""

from __future__ import annotations

import copy
import math
from typing import Any

import sympy as sp
from sympy.printing.str import StrPrinter

from g2egrid.elements.blocs import (
    _APPELS,
    _CONSTANTES,
    PARAMS_RESERVES,
    SORTIES_SERVICE,
    Compilation,
    SchemaInvalide,
    _conclure,
    _expr,
    _nom_ok,
    aplatir,
    compiler,
    lire_parametres,
)

__all__ = [
    "ENTREES_RESEAU",
    "LONGUEUR_MAX",
    "ajouter_equation_modele",
    "compiler_equations",
    "equations_depuis_element",
    "equations_depuis_instance",
    "equations_depuis_schema",
    "est_modele_equations",
    "lignes_modele",
    "modifier_equation_modele",
    "modifier_parametre_modele",
    "schema_depuis_equations",
    "supprimer_equation_modele",
]

#: les entrees du reseau, telles qu'une expression les nomme
ENTREES_RESEAU = ("vd", "vq", "V", "w", "theta_g")
#: longueur maximale d'une expression d'un modele a equations (une derivee
#: native de GFM en fait 600 : la marge d'un schema-bloc aplati)
LONGUEUR_MAX = 2000
#: au-dela, une variable intermediaire d'un schema-bloc n'est plus substituee
#: dans les equations qui l'emploient (elle reste un signal nomme)
_LONGUEUR_SUBSTITUTION = 600
#: ce que les expressions appellent ``vd``, ``vq``, ``w`` : le moteur les
#: nomme ``vd_t``, ``vq_t``, ``w_ref`` (et ``V`` est sqrt(vd_t^2 + vq_t^2))
_VERS_MOTEUR = {
    sp.Symbol("vd"): sp.Symbol("vd_t"),
    sp.Symbol("vq"): sp.Symbol("vq_t"),
    sp.Symbol("w"): sp.Symbol("w_ref"),
    sp.Symbol("V"): sp.sqrt(sp.Symbol("vd_t") ** 2 + sp.Symbol("vq_t") ** 2),
}
_VERS_UTILISATEUR = {"vd_t": "vd", "vq_t": "vq", "w_ref": "w"}
#: noms que le compilateur emploie pour ses propres variables
_NOMS_MOTEUR = {"vd_t", "vq_t", "w_ref"}
_SIGNAUX_SORTIE = ("id", "iq", "P", "Q", "theta", "w")
_INTERDITS = set(ENTREES_RESEAU) | set(_APPELS) | _CONSTANTES | _NOMS_MOTEUR


def est_modele_equations(modele: Any) -> bool:
    """Vrai pour un modele a equations (``"genre": "equations"``)."""
    return isinstance(modele, dict) and modele.get("genre") == "equations"


# ------------------------------------------------------------- ecriture
class _Imprimeur(StrPrinter):
    """L'ecriture d'une expression SymPy que l'analyseur sur relit : les
    flottants au plus court qui redonne le meme double (``repr``), ``min`` et
    ``max`` en minuscules."""

    def _print_Float(self, expr):
        return repr(float(expr))

    def _print_Min(self, expr):
        return "min(" + ", ".join(self._print(a) for a in expr.args) + ")"

    def _print_Max(self, expr):
        return "max(" + ", ".join(self._print(a) for a in expr.args) + ")"


def _ecrire(e: sp.Expr) -> str:
    return _Imprimeur().doprint(e)


def _lisible(e: sp.Expr) -> sp.Expr:
    """``e`` sur les noms de l'utilisateur : symboles simples (sans
    hypothese), ``vd_t`` -> ``vd``..., et ``sqrt(vd^2 + vq^2)`` -> ``V``."""
    sub = {s: sp.Symbol(_VERS_UTILISATEUR.get(s.name, s.name)) for s in e.free_symbols}
    e = sp.sympify(e).xreplace(sub)
    # ``1.0*x`` des natifs : un flottant entier est ecrit en entier (meme valeur)
    entiers = {
        x: sp.Integer(int(x))
        for x in e.atoms(sp.Float)
        if math.isfinite(float(x)) and abs(float(x)) < 1e15 and float(x).is_integer()
    }
    if entiers:
        e = e.xreplace(entiers)
    vd, vq, V = sp.Symbol("vd"), sp.Symbol("vq"), sp.Symbol("V")
    carre = vd**2 + vq**2
    if e.has(vd) and e.has(vq):
        e = e.replace(
            lambda x: (
                x.is_Pow and x.base == carre and x.exp.is_Rational and x.exp.q == 2
            ),
            lambda x: V ** (2 * x.exp),
        )
    return e


def _texte(e: sp.Expr) -> str:
    return _ecrire(_lisible(e))


# ----------------------------------------------------------- compilation
def verifier_forme_equations(modele: Any) -> None:
    """La FORME d'un modele a equations, avant toute lecture (il vient d'un
    fichier ou d'une requete) : ``etats`` et ``algebriques`` sont des listes
    d'objets, ``consignes`` une liste de noms ou d'objets, ``sorties`` et
    ``parametres`` des objets. Leve :class:`SchemaInvalide` avec la raison."""
    if not isinstance(modele, dict):
        raise SchemaInvalide("un modele a equations est un objet {genre, etats, ...}")
    for cle in ("etats", "algebriques", "consignes"):
        v = modele.get(cle)
        if v is None:
            continue
        if not isinstance(v, list):
            raise SchemaInvalide(f"{cle} est une liste")
        for x in v:
            if isinstance(x, dict) or (cle == "consignes" and isinstance(x, str)):
                continue
            attendu = "{nom, ...}" if cle != "consignes" else "un nom ou {nom, calee}"
            raise SchemaInvalide(f"{cle} : chaque entree est {attendu} (recu : {x!r})")
    for cle in ("sorties", "parametres"):
        v = modele.get(cle)
        if v is not None and not isinstance(v, dict):
            raise SchemaInvalide(f"{cle} est un objet {{nom: valeur}}")
    angles = modele.get("angles")
    if angles is not None and not isinstance(angles, list):
        raise SchemaInvalide("angles est une liste de noms")


def _lire(valeur: Any, noms: set[str], contexte: str) -> sp.Expr:
    """L'expression ``valeur`` (texte ou nombre) sur les symboles du moteur."""
    if isinstance(valeur, bool) or valeur is None:
        raise SchemaInvalide(f"{contexte} : expression manquante")
    if isinstance(valeur, (int, float)):
        valeur = repr(valeur)
    if not isinstance(valeur, str):
        raise SchemaInvalide(f"{contexte} : une expression est un texte ou un nombre")
    if not valeur.strip():
        raise SchemaInvalide(f"{contexte} : expression vide")
    try:
        e = _expr(valeur, noms, LONGUEUR_MAX)
    except SchemaInvalide as err:
        raise SchemaInvalide(f"{contexte} : {err}") from None
    if e.has(sp.zoo, sp.nan, sp.oo, -sp.oo):
        raise SchemaInvalide(f"{contexte} : {valeur!r} n'est pas un nombre fini")
    return e.xreplace(_VERS_MOTEUR)


def _angles(
    expressions: list[sp.Expr], variables: set[str], th_g: sp.Symbol
) -> list[str]:
    """Les variables qui portent un angle de repere : celles dont le sinus ou
    le cosinus apparait (``sin(theta - theta_g)``) -- l'initialisation les
    fait partir de l'angle de la tension du bus."""
    out: list[str] = []
    for e in expressions:
        for a in e.atoms(sp.sin, sp.cos):
            arg = a.args[0].xreplace({th_g: 0})
            if arg.could_extract_minus_sign():
                arg = -arg
            if arg.is_Symbol and arg.name in variables and arg.name not in out:
                out.append(arg.name)
    return out


def compiler_equations(modele: dict) -> Compilation:
    """Compile un modele a equations ; leve :class:`SchemaInvalide` avec la
    raison (et le nom de l'equation en cause). Rend la MEME
    :class:`Compilation` qu'un schema-bloc."""
    if not est_modele_equations(modele):
        raise SchemaInvalide('un modele a equations porte "genre": "equations"')
    verifier_forme_equations(modele)
    nom = str(modele.get("nom") or "equations")
    if not _nom_ok(nom):
        raise SchemaInvalide(
            f"nom de modele {nom!r} : lettres, chiffres et _ seulement"
        )

    def nom_variable(n: Any, role: str) -> str:
        n = str(n)
        if not _nom_ok(n):
            raise SchemaInvalide(f"{role} {n!r} : nom invalide (lettres, chiffres, _)")
        if n in _INTERDITS or n in PARAMS_RESERVES:
            raise SchemaInvalide(
                f"{role} {n!r} : nom reserve (entrees du reseau "
                f"{', '.join(ENTREES_RESEAU)}, fonctions, wb, kS)"
            )
        return n

    params: dict[str, float] = {}
    for k, v in lire_parametres(modele.get("parametres")).items():
        if k in PARAMS_RESERVES:
            raise SchemaInvalide(
                f"parametre {k!r} reserve : il est fourni par le reseau "
                "(wb : pulsation de base ; kS : Sn/Sb, colonne Sn_VA de T_Blocs)"
            )
        params[nom_variable(k, "parametre")] = v
    admis_params = set(params) | set(PARAMS_RESERVES)

    consignes: list[str] = []
    noncalees_txt: dict[str, Any] = {}
    for c in modele.get("consignes") or []:
        c = {"nom": c} if isinstance(c, str) else c
        n = nom_variable(c.get("nom"), "consigne")
        calee = c.get("calee", True)
        if not isinstance(calee, bool):
            raise SchemaInvalide(f"consigne {n} : calee doit etre vrai ou faux")
        consignes.append(n)
        if not calee:
            noncalees_txt[n] = c.get("valeur", 0)

    etats: list[str] = []
    algebs: list[str] = []
    for e in modele.get("etats") or []:
        etats.append(nom_variable(e.get("nom"), "etat"))
    for a in modele.get("algebriques") or []:
        algebs.append(nom_variable(a.get("nom"), "variable algebrique"))
    tous = etats + algebs + consignes + list(params)
    doubles = sorted({n for n in tous if tous.count(n) > 1})
    if doubles:
        raise SchemaInvalide(f"nom(s) en double : {', '.join(doubles)}")

    sorties_txt = dict(modele.get("sorties") or {})
    inconnus = sorted(set(sorties_txt) - set(_SIGNAUX_SORTIE))
    if inconnus:
        raise SchemaInvalide(
            f"sortie(s) {inconnus} inconnue(s) (attendu : {', '.join(_SIGNAUX_SORTIE)})"
        )
    if {"P", "Q"} <= set(sorties_txt) and not ({"id", "iq"} & set(sorties_txt)):
        for n in ("vd_mes", "vq_mes", "tau_mes"):
            if n in tous:
                raise SchemaInvalide(
                    f"{n!r} est reserve a la mesure de tension d'une sortie en puissance"
                )

    noms = set(tous) | admis_params | set(ENTREES_RESEAU)
    noncalees = {
        n: _lire(v, admis_params, f"consigne {n}, valeur")
        for n, v in noncalees_txt.items()
    }
    f: list[sp.Expr] = []
    g: list[sp.Expr] = []
    lignes: list[tuple[str, str, sp.Expr]] = []
    departs: dict[str, float] = {}

    def depart(entree: dict, n: str) -> None:
        if entree.get("depart") is None:
            return
        try:
            departs[n] = float(entree["depart"])
        except (TypeError, ValueError):
            raise SchemaInvalide(
                f"{n} : depart {entree['depart']!r} illisible"
            ) from None
        if not math.isfinite(departs[n]):
            raise SchemaInvalide(f"{n} : le depart doit etre fini")

    for n, e in zip(etats, modele.get("etats") or []):
        d = _lire(e.get("derivee"), noms, f"d{n}/dt")
        f.append(d)
        lignes.append(("etat", n, d))
        depart(e, n)
    for n, a in zip(algebs, modele.get("algebriques") or []):
        s = sp.Symbol(n)
        if a.get("implicite"):
            r = _lire(a.get("residu"), noms, f"0 = residu de {n}")
            g.append(r)
            lignes.append(("implicite", n, r))
        else:
            if a.get("expr") is None:
                raise SchemaInvalide(
                    f"{n} : une variable algebrique porte une expression (expr) "
                    'ou un residu (residu, "implicite": true)'
                )
            r = _lire(a.get("expr"), noms, n)
            g.append(s - r)
            lignes.append(("algebrique", n, r))
        depart(a, n)
    sortie = {k: _lire(v, noms, f"sortie {k}") for k, v in sorties_txt.items()}
    vd, vq, th_g = sp.Symbol("vd_t"), sp.Symbol("vq_t"), sp.Symbol("theta_g")
    angles = _angles(f + g + list(sortie.values()), set(etats) | set(algebs), th_g)
    for a in modele.get("angles") or []:
        if a in etats + algebs and a not in angles:
            angles.append(a)
    return _conclure(
        modele,
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
        [],
        [],
        departs,
    )


# --------------------------------------------------- lecture d'un modele
def lignes_modele(modele: dict) -> list[dict]:
    """Les equations d'un modele a equations, une par ligne :
    ``{genre, variable, expression}`` (+ ``implicite`` pour ``0 = residu``).
    ``genre`` : ``etat`` (``d<variable>/dt = expression``), ``algebrique`` ou
    ``sortie``. Les deux etats de la mesure de tension d'une sortie en
    puissance, que la compilation ajoute d'elle-meme, y figurent aussi."""
    out = []
    for e in modele.get("etats") or []:
        out.append(
            {"genre": "etat", "variable": e["nom"], "expression": str(e["derivee"])}
        )
    sorties = modele.get("sorties") or {}
    if {"P", "Q"} <= set(sorties) and not ({"id", "iq"} & set(sorties)):
        for axe in ("d", "q"):
            out.append(
                {
                    "genre": "etat",
                    "variable": f"v{axe}_mes",
                    "expression": f"(v{axe} - v{axe}_mes)/tau_mes",
                    "auto": True,
                }
            )
    for a in modele.get("algebriques") or []:
        if a.get("implicite"):
            out.append(
                {
                    "genre": "algebrique",
                    "variable": a["nom"],
                    "expression": str(a.get("residu")),
                    "implicite": True,
                }
            )
        else:
            out.append(
                {
                    "genre": "algebrique",
                    "variable": a["nom"],
                    "expression": str(a.get("expr")),
                }
            )
    for k, v in sorties.items():
        out.append({"genre": "sortie", "variable": k, "expression": str(v)})
    return out


# ------------------------------------------------------ depuis un natif
def _modele_de_dae(
    nom: str,
    parametres: dict[str, float],
    consignes: list[dict],
    etats: list[tuple[str, sp.Expr]],
    algebs: list[tuple[str, sp.Expr, bool]],
    sorties: dict[str, sp.Expr],
    amorce: dict | None,
) -> dict:
    modele: dict[str, Any] = {
        "genre": "equations",
        "nom": nom,
        "parametres": dict(parametres),
        "consignes": consignes,
        "etats": [{"nom": n, "derivee": _texte(e)} for n, e in etats],
        "algebriques": [
            {"nom": n, "residu": _texte(e), "implicite": True}
            if implicite
            else {"nom": n, "expr": _texte(e)}
            for n, e, implicite in algebs
        ],
        "sorties": {k: _texte(e) for k, e in sorties.items()},
    }
    if amorce:
        modele["amorce"] = amorce
    return modele


def equations_depuis_instance(el, nom: str | None = None) -> dict:
    """Le modele a equations d'un element natif DEJA construit (avec SES
    parametres : ceux d'une ligne de table par exemple).

    Une variable algebrique que son equation resout (lineaire en elle : toutes
    celles de la bibliotheque) s'ecrit ``nom = expr`` ; elle reste un signal
    nomme dans les derivees. Meme perimetre que
    :func:`~g2egrid.elements.blocs.depuis_element` : un seul port alternatif
    ``t`` ; un element a port continu ou a plusieurs ports est refuse."""
    from g2egrid.elements import blocs_inverse as bi

    cls = type(el)
    bi._controler_classe(cls)
    nom_classe = cls.__name__
    reduit = hasattr(cls, "_classe_complete")  # une variante ``reduire(...)``
    pdict = {k: float(v) for k, v in el.pdict.items()}
    wb = float(pdict.get("wb", bi._WB))
    f, g, h, sy = bi._equations(cls, pdict)
    etats, algebs = list(cls.states), list(cls.algebs)
    entrees_svc = list(cls.svc_inputs)
    parametres = {
        k: v for k, v in pdict.items() if k not in PARAMS_RESERVES and math.isfinite(v)
    }
    calees, fixes = bi.consignes_natives(nom_classe, entrees_svc, reduit)
    recales = bi._sonder_recales(el, wb)
    consignes = []
    for n in entrees_svc:
        if n in calees:
            consignes.append({"nom": n, "calee": True})
        else:
            valeur, recalee = fixes.get(n, (0.0, False))
            consignes.append({"nom": n, "calee": False, "valeur": valeur})
            if recalee:
                recales[n] = valeur

    appariement = bi._apparier(g, [sy[a] for a in algebs])
    lignes_alg = []
    for k, a in enumerate(algebs):
        eq = g[appariement[k]]
        sol = bi._resoudre(eq, sy[a])
        lignes_alg.append((a, eq if sol is None else sol, sol is None))

    n_port = cls.n_port_signaux()
    sorties = {"id": h[0], "iq": h[1]}
    svc = list(cls.svc_outputs)
    angle = next((k for k, s in enumerate(svc) if s in bi._SERVICE_ANGLE), None)
    vitesse = next((k for k, s in enumerate(svc) if s in bi._SERVICE_VITESSE), None)
    if angle is not None and vitesse is not None:
        sorties["theta"] = h[n_port + angle]
        sorties["w"] = h[n_port + vitesse]

    amorce = None
    if not reduit:
        amorce = {
            "classe": nom_classe,
            "correspondance": {n: n for n in (*etats, *algebs, *entrees_svc)},
            "recales": {
                k: v for k, v in recales.items() if k in parametres or k in entrees_svc
            },
        }
    modele = _modele_de_dae(
        bi._identifiant(nom) if nom else bi._identifiant(nom_classe),
        parametres,
        consignes,
        list(zip(etats, f)),
        lignes_alg,
        sorties,
        amorce,
    )
    try:
        compiler_equations(modele)
    except SchemaInvalide as err:
        raise SchemaInvalide(f"{nom_classe} : {err}") from None
    return modele


def equations_depuis_element(
    classe_ou_nom: Any,
    parametres: dict[str, float] | None = None,
    f0: float = 50.0,
    nom: str | None = None,
) -> dict:
    """Le modele a equations d'une classe NATIVE (``"GFL1"`` ou la classe),
    ``parametres`` donnes (les defauts sinon) ; ``f0`` : la frequence du
    reseau (Hz) pour les reglages automatiques du constructeur. Voir
    :func:`equations_depuis_instance`."""
    from g2egrid.elements.blocs_inverse import instance_native

    return equations_depuis_instance(
        instance_native(classe_ou_nom, parametres, f0), nom=nom
    )


# ---------------------------------------------------- depuis un schema
def equations_depuis_schema(schema: dict) -> dict:
    """Le modele a equations d'un schema-bloc (sous-systemes compris).

    Chaque equation du schema compile devient une equation du modele :
    ``d<x_bloc>/dt``, ``y_<bloc> = ...``. Les ``y_<bloc>`` purement
    intermediaires -- sortie d'un bloc SANS ETAT, qui s'explicite -- sont
    substituees dans les equations qui les emploient, tant que celles-ci
    restent lisibles (au plus 600 caracteres) ; les autres restent des signaux
    nommes, comme une boucle algebrique."""
    if est_modele_equations(schema):
        return copy.deepcopy(schema)
    from g2egrid.elements.blocs import TYPES, _noms_sorties, _var

    c = compiler(schema)
    plat = aplatir(schema)
    sans_etat: set[str] = set()
    for b in plat["blocs"]:
        t = b["type"]
        if t in ("entree", "constante", "consigne", "sortie", "implicite"):
            continue
        if TYPES[t]["etat"]:
            continue
        noms_s = _noms_sorties(b)
        sans_etat |= {_var("y", str(b["id"]), noms_s, k) for k in range(len(noms_s))}

    explicites: dict[str, sp.Expr] = {}
    residus: dict[str, sp.Expr] = {}
    for genre, var, e in c.lignes:
        if var in c.algebriques:
            (explicites if genre == "algebrique" else residus)[var] = e
    derivees = dict(zip(c.etats, c.f))
    sorties = dict(zip(["id", "iq"], c.h[:2]))
    inverse = {v: k for k, v in SORTIES_SERVICE.items()}
    for nom_service, e in zip(c.services, c.h[2:]):
        sorties[inverse[nom_service]] = e
    conteneurs = [derivees, explicites, residus, sorties]

    def texte_max(d: dict, cles) -> int:
        return max((len(_texte(d[k])) for k in cles), default=0)

    for y in c.algebriques:
        if y not in explicites:
            continue
        rhs = explicites[y]
        # un bloc a etat n'est substitue que si sa sortie est l'etat lui-meme
        if y not in sans_etat and not rhs.is_Symbol:
            continue
        s = sp.Symbol(y)
        if s in rhs.free_symbols:
            continue
        touches = [
            (d, [k for k, e in d.items() if k != y and s in e.free_symbols])
            for d in conteneurs
        ]
        nouveaux = {
            id(d): {k: d[k].xreplace({s: rhs}) for k in cles} for d, cles in touches
        }
        if any(
            len(_texte(e)) > _LONGUEUR_SUBSTITUTION
            for d, cles in touches
            for e in nouveaux[id(d)].values()
        ):
            continue
        for d, _cles in touches:
            d.update(nouveaux[id(d)])
        del explicites[y]
    # la substitution a pu rendre une variable auto-referente : residu
    for y in [y for y, e in explicites.items() if sp.Symbol(y) in e.free_symbols]:
        residus[y] = sp.Symbol(y) - explicites.pop(y)
    algebs = []
    for y in c.algebriques:
        if y in explicites:
            algebs.append((y, explicites[y], False))
        elif y in residus:
            algebs.append((y, residus[y], True))
    consignes = []
    for n in c.consignes:
        if n in c.valeurs:
            v = c.valeurs[n]
            consignes.append(
                {
                    "nom": n,
                    "calee": False,
                    "valeur": float(v) if v.is_Number else _texte(v),
                }
            )
        else:
            consignes.append({"nom": n, "calee": True})
    parametres = {k: v for k, v in c.parametres.items() if k not in PARAMS_RESERVES}
    modele = _modele_de_dae(
        c.nom,
        parametres,
        consignes,
        [(n, derivees[n]) for n in c.etats],
        algebs,
        sorties,
        copy.deepcopy(schema.get("amorce")) if schema.get("amorce") else None,
    )
    compiler_equations(modele)
    return modele


def schema_depuis_equations(modele: dict, detail: str = "blocs") -> dict:
    """Le schema-bloc d'un modele a equations, pour l'editeur de blocs : le
    compilateur inverse de FR-051 lit le modele comme un DAE.

    ``detail`` : ``"blocs"`` (un sous-systeme par equation), ``"equations"``
    (un bloc ``fonction`` par equation) ; ``"liste"`` rend le modele tel quel.
    Les ``depart`` ne sont pas repris (le schema n'en a pas)."""
    from g2egrid.elements import blocs_inverse as bi

    if detail == "liste":
        return copy.deepcopy(modele)
    if detail not in ("blocs", "equations"):
        raise SchemaInvalide(f"detail {detail!r} : 'blocs', 'equations' ou 'liste'")
    c = compiler_equations(modele)
    noms = [*c.etats, *c.algebriques, *c.consignes, *c.parametres]
    reel = {
        n: sp.Symbol(n, real=True) for n in {*noms, "vd_t", "vq_t", "w_ref", "theta_g"}
    }
    sub = {sp.Symbol(n): s for n, s in reel.items()}

    def reel_(e: sp.Expr) -> sp.Expr:
        return e.xreplace(sub)

    pdict = dict(c.parametres)
    amorce = modele.get("amorce") or {}
    fixes = {}
    for n, v in c.valeurs.items():
        val = reel_(v).xreplace({reel[k]: x for k, x in pdict.items()})
        fixes[n] = (float(val), False)
    dae = bi._Dae(
        nom_classe=str(amorce.get("classe") or c.nom),
        reduit=not amorce.get("classe"),
        pdict=pdict,
        params={*pdict, *PARAMS_RESERVES},
        etats=list(c.etats),
        algebs=list(c.algebriques),
        entrees_svc=list(c.consignes),
        sorties_svc=list(c.services),
        n_port=2,
        f=[reel_(e) for e in c.f],
        g=[reel_(e) for e in c.g],
        h=[reel_(e) for e in c.h],
        symboles=reel,
        calees=tuple(c.calees),
        fixes=fixes,
        recales={
            k: float(v)
            for k, v in (amorce.get("recales") or {}).items()
            if k in pdict or k in c.consignes
        },
    )
    return bi.schema_depuis_dae(dae, detail, c.nom)


# ------------------------------------------------------ modifications
def _copie(modele: dict) -> dict:
    if not est_modele_equations(modele):
        raise SchemaInvalide('un modele a equations porte "genre": "equations"')
    return copy.deepcopy(modele)


def _trouver(modele: dict, variable: str, genre: str | None = None):
    """``(liste, entree, cle_expression, genre)`` de l'equation ``variable`` :
    etats, algebriques, puis sorties (``genre`` leve l'ambiguite d'un nom
    comme ``theta``, a la fois etat et sortie)."""
    cherches = (genre,) if genre else ("etat", "algebrique", "sortie")
    for g in cherches:
        if g == "etat":
            for e in modele.get("etats") or []:
                if e["nom"] == variable:
                    return e, "derivee", "etat"
        elif g == "algebrique":
            for a in modele.get("algebriques") or []:
                if a["nom"] == variable:
                    return a, ("residu" if a.get("implicite") else "expr"), g
        elif g == "sortie":
            if variable in (modele.get("sorties") or {}):
                return modele["sorties"], variable, g
        else:
            raise SchemaInvalide(
                f"genre {g!r} : 'etat', 'algebrique' ou 'sortie' attendu"
            )
    raise SchemaInvalide(
        f"le modele {modele.get('nom')!r} n'a pas d'equation pour {variable!r}"
    )


def modifier_equation_modele(
    modele: dict, variable: str, expression: str, genre: str | None = None
) -> dict:
    """Une COPIE du modele ou l'equation de ``variable`` est ``expression``
    (la derivee d'un etat, l'expression d'une algebrique, une sortie). Le
    modele modifie est compile avant d'etre rendu : un nom inconnu, une
    syntaxe refusee, une boucle... leve :class:`SchemaInvalide`."""
    m = _copie(modele)
    entree, cle, _g = _trouver(m, variable, genre)
    entree[cle] = expression
    compiler_equations(m)
    return m


def ajouter_equation_modele(
    modele: dict,
    genre: str,
    variable: str,
    expression: str,
    depart: float | None = None,
) -> dict:
    """Une COPIE du modele avec une equation de plus : ``genre`` ``"etat"``
    (``d<variable>/dt = expression``), ``"algebrique"`` (``variable =
    expression``) ou ``"implicite"`` (``0 = expression``, ``variable`` en
    etant l'inconnue). ``depart`` : la valeur d'ou part le Newton de
    l'initialisation pour cette variable (0 sinon)."""
    m = _copie(modele)
    if genre not in ("etat", "algebrique", "implicite"):
        raise SchemaInvalide(
            f"genre {genre!r} : 'etat', 'algebrique' ou 'implicite' attendu"
        )
    if not isinstance(variable, str) or not _nom_ok(variable):
        raise SchemaInvalide(f"nom de variable {variable!r} invalide")
    pris = {e["nom"] for e in m.get("etats") or []}
    pris |= {a["nom"] for a in m.get("algebriques") or []}
    pris |= {c["nom"] if isinstance(c, dict) else c for c in m.get("consignes") or []}
    pris |= set(m.get("parametres") or {})
    if variable in pris:
        raise SchemaInvalide(f"le nom {variable!r} est deja pris dans le modele")
    if genre == "etat":
        entree: dict[str, Any] = {"nom": variable, "derivee": expression}
        m.setdefault("etats", []).append(entree)
    elif genre == "algebrique":
        entree = {"nom": variable, "expr": expression}
        m.setdefault("algebriques", []).append(entree)
    else:
        entree = {"nom": variable, "residu": expression, "implicite": True}
        m.setdefault("algebriques", []).append(entree)
    if depart is not None:
        entree["depart"] = float(depart)
    compiler_equations(m)
    return m


def _references(modele: dict, variable: str) -> list[str]:
    """Les equations du modele qui emploient ``variable``."""
    import re

    motif = re.compile(rf"(?<![A-Za-z0-9_]){re.escape(variable)}(?![A-Za-z0-9_])")
    out = []
    for ligne in lignes_modele(modele):
        if ligne["variable"] == variable and ligne["genre"] != "sortie":
            continue
        if motif.search(ligne["expression"]):
            out.append(f"{ligne['genre']} {ligne['variable']}")
    return out


def supprimer_equation_modele(modele: dict, variable: str) -> dict:
    """Une COPIE du modele sans l'equation de l'etat ou de la variable
    algebrique ``variable``, qui ne doit plus etre employee ailleurs (sinon
    :class:`SchemaInvalide` dit par qui). Un etat natif d'une sortie ne se
    retire pas : modifier l'equation qui l'emploie d'abord."""
    m = _copie(modele)
    for cle in ("etats", "algebriques"):
        liste = m.get(cle) or []
        for k, e in enumerate(liste):
            if e["nom"] == variable:
                m[cle] = liste[:k] + liste[k + 1 :]
                emplois = _references(m, variable)
                if emplois:
                    raise SchemaInvalide(
                        f"{variable!r} est encore employe par : "
                        f"{', '.join(emplois)} (modifier d'abord ces equations)"
                    )
                compiler_equations(m)
                return m
    raise SchemaInvalide(
        f"le modele {modele.get('nom')!r} n'a ni etat ni variable algebrique "
        f"{variable!r}"
    )


def modifier_parametre_modele(modele: dict, nom: str, valeur: float) -> dict:
    """Une COPIE du modele ou le parametre ``nom`` vaut ``valeur``."""
    m = _copie(modele)
    if nom not in (m.get("parametres") or {}):
        raise SchemaInvalide(f"le modele {m.get('nom')!r} n'a pas de parametre {nom!r}")
    try:
        m["parametres"][nom] = float(valeur)
    except (TypeError, ValueError):
        raise SchemaInvalide(
            f"parametre {nom!r} : {valeur!r} n'est pas un nombre"
        ) from None
    compiler_equations(m)
    return m
