"""Le compilateur INVERSE : equations algebro-differentielles -> blocs.

:func:`depuis_element` ecrit le schema-bloc d'une classe NATIVE de la
bibliotheque (``GFL1``, ``GFM2``, ``SynchronousMachine``, ``IB``, ``RLLoad``...)
a partir de ses equations symboliques ``f, g, h`` : de quoi OUVRIR un
equipement de la bibliotheque dans l'editeur, le modifier bloc a bloc, et le
poser dans un reseau a la place du natif (``T_Blocs``).

CE QUI EST GENERE, pour une classe a un seul port alternatif ``t`` --

* un ``integrateur`` (K = 1) par ETAT, nomme d'apres l'etat ;
* un SOUS-SYSTEME par equation d'etat (``d<etat>/dt``) qui calcule la derivee
  en blocs elementaires : ``somme`` (avec ses signes), ``gain`` (facteur
  numerique ou parametrique), ``produit``, ``division``, ``fonction`` (sin,
  cos, sqrt, min, max...), ``constante`` (expression des parametres).
  ``detail="equations"`` met a la place UN SEUL bloc ``fonction`` par equation ;
* chaque variable ALGEBRIQUE devient un sous-systeme qui la calcule quand son
  equation se resout en elle (lineaire dans cette variable : c'est le cas de
  toutes les classes de la bibliotheque), sinon un bloc ``implicite`` ;
* les entrees du reseau en blocs ``entree`` (``vd``, ``vq``, ``w``,
  ``theta_g``), les courants injectes en blocs ``sortie`` ``id`` et ``iq``,
  l'angle et la vitesse des FORMEURS en sorties de service ``theta`` et ``w``
  (la classe native qui expose ``theta_out`` et ``Omega``/``omega``/``w``).
  Les AUTRES sorties de service (P, Q, Vt, vdc_out...) ne sont PAS reprises :
  ce sont des mesures, qu'aucune equation du modele n'emploie ;
* les entrees de service natives en blocs ``consigne`` : deux CALEES par
  l'initialisation, les autres non calees a une valeur fixe -- le choix, par
  classe, est :data:`CONSIGNES` ;
* les parametres natifs (valeurs par defaut, ou ``parametres`` donnes) en
  ``parametres`` du schema. ``kS`` et ``wb`` restent reserves : l'assembleur
  les fournit (``Sn_VA`` de ``T_Blocs``, frequence du reseau).

LE FONCTIONNEMENT EST LIBRE -- les drapeaux de mode (``mode_lim``,
``mode_vsat``, ``mode_e2``) deviennent des parametres a 0 et les gardes de
limitation du natif ne sont PAS reprises : le modele a blocs genere decrit le
fonctionnement libre du convertisseur, celui sur lequel la linearisation est
faite (a drapeau nul les equations natives sont exactement celles d'un
fonctionnement sans limite). Il en va de meme de l'impedance virtuelle de
limitation (``kRv``) quand ``kRv = 0`` a la generation (le defaut) : le terme
est nul, il n'est pas dessine. Avec ``kRv != 0`` il l'est (le depassement lisse
``ov(|i|)`` s'ecrit alors avec ``min`` et ``max``, a la place du ``Piecewise``
natif ; memes valeurs et meme derivee).

AMORCE DE L'INITIALISATION -- le schema porte ::

    "amorce": {"classe": "GFL2", "correspondance": {"x_isd": "isd", ...},
               "recales": {"Kpdc": 1.0}}

``classe_element`` appelle d'abord l'``initialize`` de la classe native (memes
parametres, meme contexte) et s'en sert comme point de depart du Newton ; si
l'amorce echoue il revient au depart d'avant. ``recales`` (extension du
contrat) liste les parametres -- ou consignes non calees -- que
l'initialisation native RECALCULE au point de fonctionnement (le ``Kpdc`` d'un
GFL2, le ``theta_g0`` d'une machine, le ``R`` et le ``X`` d'une charge, la
frequence ``w_set`` d'un bus infini) avec leur valeur de depart : tant que
l'element les porte encore, ils prennent la valeur du natif, si bien que le
modele a blocs a le spectre du natif ; un parametre que l'utilisateur a change
n'est pas touche (comme un gain donne au constructeur natif).

DISPOSITION -- :func:`disposer` place chaque bloc, a chaque niveau : flux de
gauche a droite (niveau = plus long chemin depuis les sources, les boucles
sont coupees a leur arc de retour), ordre dans la colonne par barycentre, pas
de ``152 x 70`` pour des boites de ``BLOC_L x BLOC_H`` : aucun recouvrement.
"""

from __future__ import annotations

import copy
import inspect
import math
import re
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any

import numpy as np
import sympy as sp

from g2egrid.elements.blocs_structure import (
    BLOC_H,
    BLOC_L,
    SchemaInvalide,
)

__all__ = [
    "CONSIGNES",
    "classe_native",
    "consignes_natives",
    "depuis_element",
    "depuis_instance",
    "disposer",
    "initialisation_native",
    "instance_native",
]

#: pulsation de base par defaut (50 Hz) des natifs construits hors reseau
_WB = 2 * math.pi * 50.0
#: les sorties de service qui font un formeur : (nom natif -> signal du bloc)
_SERVICE_ANGLE = ("theta_out",)
_SERVICE_VITESSE = ("Omega", "omega", "w")
#: entrees du reseau : symbole natif -> signal du bloc ``entree``
_ENTREES = {"vd_t": "vd", "vq_t": "vq", "w_ref": "w", "theta_g": "theta_g"}
#: au-dela, une expression n'est plus ecrite en un seul bloc ``fonction``
_LONGUEUR_FONCTION = 700
_PAS_X = BLOC_L + 60
_PAS_Y = BLOC_H + 20

#: LES CONSIGNES de chaque classe native : ``(calees, non calees)`` ; les non
#: calees sont ``{nom: (valeur, recalee)}`` -- ``recalee`` : la valeur de
#: depart est remplacee par celle de l'initialisation native (voir
#: ``amorce["recales"]``). Deux consignes calees exactement, choisies pour que
#: l'initialisation reussisse ET que le spectre soit celui du natif :
#:
#: * GFL1 (p_ref, q_ref) : les deux commandent P et Q.
#: * GFL2 : ``q_ref`` et ``idc_ref`` (le courant continu que la source doit
#:   fournir, qui boucle le bilan de puissance) ; ``vdc_ref`` est un reglage
#:   (il ne fait que decaler l'integrateur de la boucle de tension).
#: * GFM1/GFM2 : ``p_ref`` et ``ve_ref`` ; la dynamique ne voit que la somme
#:   ``ve_ref + nq.q_ref`` (la f.e.m. de consigne), donc ``q_ref`` est fixee a 0.
#: * GFMdVOC : ``p_ref`` et ``q_ref`` ; ``ve_ref`` est la tension de
#:   reference de l'oscillateur, DONT DEPEND la dynamique (``p_ref/ve_ref^2``) :
#:   elle prend la valeur du natif (recalee).
#: * GFMMatching : ``p_ref`` (qui fixe l'idc de la source) et ``ve_ref`` ;
#:   ``q_ref`` a 0 comme pour le statisme ; ``vdc_ref`` (la frequence est
#:   ``w.vdc/vdc_ref``) prend la valeur du natif.
#: * IB : ``E_ref`` et ``delta_ref`` (la f.e.m. complexe) ; ``w_set``, la
#:   frequence imposee (celle du reseau a l'equilibre), prend celle du natif.
#: * machine synchrone (p_ref, v_ref), charge R-L (dp, dq).
#:
#: Une classe absente : les deux premieres entrees de service sont calees, les
#: autres fixees a 0.
CONSIGNES: dict[str, tuple[tuple[str, ...], dict[str, tuple[float, bool]]]] = {
    "GFL1": (("p_ref", "q_ref"), {}),
    "GFL2": (("q_ref", "idc_ref"), {"vdc_ref": (1.0, True)}),
    "GFM1": (("p_ref", "ve_ref"), {"q_ref": (0.0, False)}),
    "GFM2": (("p_ref", "ve_ref"), {"q_ref": (0.0, False)}),
    "GFMdVOC": (("p_ref", "q_ref"), {"ve_ref": (1.0, True)}),
    "GFMMatching": (
        ("p_ref", "ve_ref"),
        {"q_ref": (0.0, False), "vdc_ref": (1.0, True)},
    ),
    "GFM3": (
        ("p_ref", "ve_ref"),
        {"q_ref": (0.0, False), "vdc_ref": (1.0, True)},
    ),
    "GFM4": (
        ("p_ref", "ve_ref"),
        {"q_ref": (0.0, False), "vdc_ref": (1.0, True)},
    ),
    "IB": (("E_ref", "delta_ref"), {"w_set": (1.0, True)}),
    "SynchronousMachine": (("p_ref", "v_ref"), {}),
    "SynchronousMachineTR": (("p_ref", "v_ref"), {}),
    "RLLoad": (("dp", "dq"), {}),
}


# ------------------------------------------------------------ classes natives
def classe_native(nom: Any) -> type | None:
    """La classe native de ce nom (ou la classe elle-meme), ``None`` si elle
    n'existe pas dans :mod:`g2egrid.elements`."""
    from g2egrid import elements
    from g2egrid.noyau.element import Element

    if isinstance(nom, type):
        return nom if issubclass(nom, Element) else None
    cls = getattr(elements, str(nom), None)
    return cls if isinstance(cls, type) and issubclass(cls, Element) else None


def _pq(ctx) -> dict[str, float]:
    """``P``, ``Q`` consommes que le constructeur d'une charge demande : lus
    sur le contexte d'initialisation (injection = -consommation), ou une
    valeur d'essai."""
    if ctx is None:
        return {"P": 0.5, "Q": 0.2}
    S = complex(ctx.s.get("t", 0j))
    return {"P": -S.real, "Q": -S.imag}


#: arguments de constructeur que ni les parametres ni les bases ne fournissent
_EXTRAS_NATIFS = {
    "RLLoad": _pq,
    "ZLoad": _pq,
    "ILoad": _pq,
    "PLoad": _pq,
}


def _construire_natif(cls: type, params: dict, kS: float, wb: float, ctx=None):
    """Une instance de ``cls`` (hors reseau) : ``Sn = kS``, ``Sb = 1`` (le
    rapport ``kS`` est tout ce qui compte), la pulsation ``wb``, ces
    ``params``, et ce que le constructeur demande en plus pour les charges."""
    sig = inspect.signature(cls.__init__)
    kw: dict[str, Any] = {}
    extras = _EXTRAS_NATIFS.get(cls.__name__)
    donnes = dict(extras(ctx)) if extras else {}
    variadique = False
    for n, prm in sig.parameters.items():
        if n in ("self", "name", "bus", "buses"):
            continue
        if prm.kind is inspect.Parameter.VAR_KEYWORD:
            variadique = True
        elif n == "Sn":
            kw["Sn"] = kS
        elif n == "Sb":
            kw["Sb"] = 1.0
        elif n == "wb":
            kw["wb"] = wb
        elif n in donnes:
            kw[n] = donnes[n]
        elif n in params:
            kw[n] = params[n]
    if variadique:
        for k, v in params.items():
            if k not in kw and k not in ("kS", "wb"):
                kw[k] = v
    try:
        return cls("_", 1, **kw)
    except TypeError as err:
        raise SchemaInvalide(
            f"{cls.__name__} : parametre obligatoire manquant ou invalide ({err})"
        ) from None


def initialisation_native(amorce: dict, pdict: dict, ctx) -> tuple[dict, dict]:
    """``(initialisation, parametres)`` de la classe native de ``amorce``,
    construite avec les parametres ``pdict`` du modele a blocs et initialisee
    dans le MEME contexte. Les parametres ``recales`` encore a leur valeur de
    depart ne sont pas donnes au constructeur : il les regle lui-meme.
    Leve une exception si la classe native n'existe pas ou refuse le point."""
    cls = classe_native(amorce.get("classe"))
    if cls is None:
        raise ValueError(f"classe native {amorce.get('classe')!r} inconnue")
    recales = amorce.get("recales") or {}
    params = {
        k: v
        for k, v in pdict.items()
        if k in cls.params and not (k in recales and v == recales[k])
    }
    el = _construire_natif(
        cls, params, float(pdict.get("kS", 1.0)), float(pdict.get("wb", _WB)), ctx
    )
    init = el.initialize(ctx)
    return init, el.pdict


# -------------------------------------------------------- equations natives
@contextmanager
def _depassement_sans_piecewise():
    """Le depassement lisse ``ov(i)`` de l'impedance virtuelle de limitation,
    ecrit avec ``min`` et ``max`` plutot qu'un ``Piecewise`` (que les blocs ne
    savent pas dessiner) : meme valeur, meme derivee (C1) ::

        ov = min(max(i - I, 0), dI)^2 / (2 dI) + max(i - I - dI, 0)
    """
    from g2egrid.elements import limiteurs

    ancien = limiteurs.depassement_vi

    def ov(i, I_vi, dI_vi):
        a = i - I_vi
        return sp.Min(sp.Max(a, 0), dI_vi) ** 2 / (2 * dI_vi) + sp.Max(a - dI_vi, 0)

    limiteurs.depassement_vi = ov
    try:
        yield
    finally:
        limiteurs.depassement_vi = ancien


def _equations(cls: type, pdict: dict[str, float]):
    """``(f, g, h, symboles)`` de la classe, drapeaux de mode a 0.

    Les symboles portent les noms natifs ; ``symboles`` les donne par nom."""
    from types import SimpleNamespace

    ns = SimpleNamespace()
    symboles: dict[str, sp.Symbol] = {}
    for grp in (cls.params, cls.states, cls.algebs, cls.input_names()):
        for n in grp:
            symboles[n] = sp.Symbol(n, real=True)
            setattr(ns, n, symboles[n])
    kRv = pdict.get("kRv")
    ctx = _depassement_sans_piecewise() if kRv else _nul()
    with ctx:
        f, g, h = cls.equations(ns)
    libre = {
        symboles[n]: sp.Integer(0)
        for n in cls.params
        if n.startswith("mode_") or n in getattr(cls, "modes", ())
    }
    if "kRv" in symboles and not kRv:
        libre[symboles["kRv"]] = sp.Integer(0)
    return (
        [sp.sympify(e).xreplace(libre) for e in f],
        [sp.sympify(e).xreplace(libre) for e in g],
        [sp.sympify(e).xreplace(libre) for e in h],
        symboles,
    )


@contextmanager
def _nul():
    yield


# ---------------------------------------------------- expressions -> blocs
def _identifiant(nom: str) -> str:
    """Un identifiant de bloc valide (lettres, chiffres, ``_``, sans ``__``)."""
    s = re.sub(r"[^0-9a-zA-Z_]", "_", str(nom))
    while "__" in s:
        s = s.replace("__", "_")
    s = s.strip("_")
    return s if s and not s[0].isdigit() else f"v_{s}"


def _nombre(e: sp.Expr) -> str:
    """Un nombre sympy en texte, sans perdre un bit."""
    if e.is_Integer:
        return str(int(e))
    return repr(float(e))


def _texte(e: sp.Expr) -> str:
    """Une expression de parametres en texte relisible par le compilateur."""
    if e.is_Number:
        return _nombre(e)
    t = sp.sstr(e, full_prec=True)
    return t.replace("Min(", "min(").replace("Max(", "max(")


_FONCTIONS = {
    sp.sin: "sin",
    sp.cos: "cos",
    sp.tan: "tan",
    sp.exp: "exp",
    sp.log: "log",
    sp.tanh: "tanh",
    sp.atan2: "atan2",
    sp.Abs: "Abs",
    sp.Min: "min",
    sp.Max: "max",
}


class _Constructeur:
    """Construit les blocs d'UN niveau (un sous-systeme, ou le premier niveau)
    qui calculent des expressions SymPy.

    ``variable(nom)`` donne la source ``(bloc, sortie)`` d'une variable du
    modele (un port d'entree du sous-systeme, cree a la premiere demande, ou
    la source du premier niveau) ; les symboles dont le nom est dans
    ``params`` sont des PARAMETRES, qui n'ont pas de source (ils s'ecrivent
    dans les champs des blocs)."""

    def __init__(self, params: set[str], variable, pris: set[str] | None = None):
        self.params = params
        self.variable = variable
        self.blocs: list[dict] = []
        self.liens: list[dict] = []
        self.pris = pris if pris is not None else set()
        self.memo: dict[sp.Expr, tuple[str, int]] = {}
        self.compteurs: dict[str, int] = {}

    # --- outils
    def id(self, prefixe: str) -> str:
        while True:
            k = self.compteurs.get(prefixe, 0) + 1
            self.compteurs[prefixe] = k
            bid = f"{prefixe}{k}"
            if bid not in self.pris:
                self.pris.add(bid)
                return bid

    def bloc(self, prefixe: str, champs: dict, entrees: list, bid: str | None = None):
        if bid is None:
            bid = self.id(prefixe)
        else:
            self.pris.add(bid)
        self.blocs.append({"id": bid, **champs})
        for k, (src, ks) in enumerate(entrees):
            lien = {"de": src, "vers": bid, "entree": k}
            if ks:
                lien["sortie"] = ks
            self.liens.append(lien)
        return (bid, 0)

    def parametrique(self, e: sp.Expr) -> bool:
        """``e`` ne depend que de nombres et de parametres (et sait s'ecrire
        dans un champ : ni ``min`` ni ``max`` ni ``Piecewise``)."""
        return all(str(s) in self.params for s in e.free_symbols) and not e.has(
            sp.Min, sp.Max, sp.Piecewise
        )

    def constante(self, e: sp.Expr) -> tuple[str, int]:
        return self.bloc("c", {"type": "constante", "valeur": _texte(e)}, [])

    # --- le parcours
    def noeud(self, e: sp.Expr, bid: str | None = None) -> tuple[str, int]:
        """La source ``(bloc, sortie)`` de la valeur de ``e``."""
        if bid is None and e in self.memo:
            return self.memo[e]
        if e.is_Symbol and str(e) not in self.params:
            src = self.variable(str(e))
        elif self.parametrique(e):
            src = self.bloc("c", {"type": "constante", "valeur": _texte(e)}, [], bid)
        elif e.is_Add:
            src = self._somme(e, bid)
        elif e.is_Mul:
            src = self._produit(e, bid)
        elif e.is_Pow:
            src = self._puissance(e, bid)
        elif e.func in _FONCTIONS:
            src = self._fonction_native(e, bid)
        elif isinstance(e, sp.Piecewise):
            raise SchemaInvalide(
                "une expression par morceaux (Piecewise) n'a pas d'equivalent en "
                f"blocs : {sp.sstr(e)[:120]}"
            )
        else:
            raise SchemaInvalide(f"expression non prise en charge : {sp.sstr(e)[:120]}")
        if bid is not None and src[0] != bid:
            # la racine demandee sous un identifiant est une variable ou un
            # sous-resultat deja construit : un gain unitaire la porte
            src = self.bloc("k", {"type": "gain", "K": 1}, [src], bid)
        self.memo.setdefault(e, src)
        return src

    def _scinder(self, t: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
        """``t = coef . reste``, ``coef`` ne dependant que des parametres."""
        if not t.is_Mul:
            return (
                (sp.Integer(1), t) if not self.parametrique(t) else (t, sp.Integer(1))
            )
        coef = [a for a in t.args if self.parametrique(a)]
        reste = [a for a in t.args if not self.parametrique(a)]
        return sp.Mul(*coef), sp.Mul(*reste)

    def _somme(self, e: sp.Expr, bid) -> tuple[str, int]:
        entrees, signes = [], ""
        constantes = []
        for t in e.args:
            if self.parametrique(t):
                constantes.append(t)
                continue
            coef, reste = self._scinder(t)
            neg = bool(coef.could_extract_minus_sign())
            if neg:
                coef = -coef
            src = self.noeud(reste)
            if coef != 1:
                src = self.bloc("k", {"type": "gain", "K": _texte(coef)}, [src])
            entrees.append(src)
            signes += "-" if neg else "+"
        if constantes:
            c = sp.Add(*constantes)
            neg = bool(c.could_extract_minus_sign())
            entrees.append(self.noeud(-c if neg else c))
            signes += "-" if neg else "+"
        return self.bloc("s", {"type": "somme", "signes": signes}, entrees, bid)

    def _produit(self, e: sp.Expr, bid) -> tuple[str, int]:
        coef, reste = self._scinder(e)
        num, den = [], []
        for a in sp.Mul.make_args(reste):
            base, ex = a.as_base_exp()
            if ex.is_Integer and ex < 0:
                den.append(base ** (-ex))
            else:
                num.append(a)
        src = self._chaine(num)
        if den:
            denominateur = self.noeud(sp.Mul(*den))
            src = self.bloc(
                "d",
                {"type": "division"},
                [src, denominateur],
                None if coef != 1 else bid,
            )
        elif coef == 1 and bid is None and len(num) == 1:
            return src
        if coef != 1:
            src = self.bloc("k", {"type": "gain", "K": _texte(coef)}, [src], bid)
        return src

    def _chaine(self, facteurs: list[sp.Expr]) -> tuple[str, int]:
        """Le produit des ``facteurs`` en blocs ``produit`` deux a deux."""
        if not facteurs:
            return self.noeud(sp.Integer(1))
        srcs = [self.noeud(f) for f in facteurs]
        acc = srcs[0]
        for s in srcs[1:]:
            acc = self.bloc("p", {"type": "produit"}, [acc, s])
        return acc

    def _puissance(self, e: sp.Expr, bid) -> tuple[str, int]:
        base, ex = e.as_base_exp()
        if ex == 2:
            s = self.noeud(base)
            return self.bloc("p", {"type": "produit"}, [s, s], bid)
        if ex.is_Integer and ex < 0:
            un = self.noeud(sp.Integer(1))
            return self.bloc(
                "d", {"type": "division"}, [un, self.noeud(base ** (-ex))], bid
            )
        if ex == sp.Rational(1, 2):
            return self.bloc(
                "f",
                {"type": "fonction", "expr": "sqrt(u1)", "n": 1},
                [self.noeud(base)],
                bid,
            )
        s = self.noeud(base)
        return self.bloc(
            "f",
            {"type": "fonction", "expr": f"u1**({_texte(ex)})", "n": 1},
            [s],
            bid,
        )

    def _fonction_native(self, e: sp.Expr, bid) -> tuple[str, int]:
        nom = _FONCTIONS[e.func]
        n = len(e.args)
        args = ", ".join(f"u{k + 1}" for k in range(n))
        return self.bloc(
            "f",
            {"type": "fonction", "expr": f"{nom}({args})", "n": n},
            [self.noeud(a) for a in e.args],
            bid,
        )

    # --- un seul bloc fonction (detail="equations")
    def fonction(self, e: sp.Expr, bid: str | None = None) -> tuple[str, int]:
        """``e`` en UN bloc ``fonction`` dont les entrees sont les variables
        qu'il emploie ; au-dela de :data:`_LONGUEUR_FONCTION` caracteres,
        l'expression est coupee a sa racine (chaque argument devient un bloc
        a son tour, la racine les combine)."""
        if e.is_Symbol and str(e) not in self.params:
            return self.variable(str(e))
        if self.parametrique(e):
            return self.bloc("c", {"type": "constante", "valeur": _texte(e)}, [], bid)
        variables = [
            s
            for s in sorted(e.free_symbols, key=lambda s: str(s))
            if str(s) not in self.params
        ]
        u = {s: sp.Symbol(f"u{k + 1}") for k, s in enumerate(variables)}
        texte = _texte(e.xreplace(u))
        if len(texte) <= _LONGUEUR_FONCTION:
            return self.bloc(
                "f",
                {"type": "fonction", "expr": texte, "n": len(variables)},
                [self.variable(str(s)) for s in variables],
                bid,
            )
        # trop long : on coupe a la racine
        enfants = [self.fonction(a) for a in e.args]
        sub = {a: sp.Symbol(f"u{k + 1}") for k, a in enumerate(e.args)}
        racine = (
            e.func(*[sub[a] for a in e.args])
            if len(set(e.args)) == len(e.args)
            else None
        )
        if racine is None:
            racine = e.func(*[sp.Symbol(f"u{k + 1}") for k in range(len(e.args))])
        return self.bloc(
            "f",
            {"type": "fonction", "expr": _texte(racine), "n": len(e.args)},
            enfants,
            bid,
        )


# ------------------------------------------------------------- disposition
def disposer(schema: dict, chemin: list[str] | tuple[str, ...] = ()) -> dict:
    """Une copie de ``schema`` ou chaque bloc, a chaque niveau (sous-systemes
    compris), a une position ``x``, ``y`` : flux de gauche a droite, aucun
    recouvrement des boites ``BLOC_L x BLOC_H`` (pas de 152 x 70).

    ``chemin`` (identifiants de sous-systemes) limite la disposition a ce
    niveau, et a ce qu'il contient."""
    sch = copy.deepcopy(schema)
    niveau = sch
    for sid in chemin or ():
        b = next((x for x in niveau.get("blocs") or [] if x.get("id") == sid), None)
        if b is None or b.get("type") != "sous_systeme":
            raise SchemaInvalide(f"chemin : {sid!r} n'est pas un sous-systeme")
        niveau = b["schema"]
    _disposer_niveau(niveau)
    return sch


def _disposer_niveau(niveau: dict) -> None:
    blocs = niveau.get("blocs") or []
    for b in blocs:
        if b.get("type") == "sous_systeme" and isinstance(b.get("schema"), dict):
            _disposer_niveau(b["schema"])
    if not blocs:
        return
    ids = [b["id"] for b in blocs]
    pos = {bid: k for k, bid in enumerate(ids)}
    succ: dict[str, list[str]] = {bid: [] for bid in ids}
    pred: dict[str, list[str]] = {bid: [] for bid in ids}
    for lien in niveau.get("liens") or []:
        de, vers = lien.get("de"), lien.get("vers")
        if de in pos and vers in pos and de != vers and vers not in succ[de]:
            succ[de].append(vers)
            pred[vers].append(de)

    # arcs de retour (cycles) par un parcours en profondeur, depuis les sources
    etat: dict[str, int] = {}
    retour: set[tuple[str, str]] = set()
    racines = [i for i in ids if not pred[i]] + ids
    for r in racines:
        if r in etat:
            continue
        pile = [(r, iter(succ[r]))]
        etat[r] = 1
        while pile:
            u, it = pile[-1]
            for v in it:
                if etat.get(v) == 1:
                    retour.add((u, v))
                elif v not in etat:
                    etat[v] = 1
                    pile.append((v, iter(succ[v])))
                    break
            else:
                etat[u] = 2
                pile.pop()
    sans_retour = {u: [v for v in succ[u] if (u, v) not in retour] for u in ids}
    entrants = {i: 0 for i in ids}
    for u in ids:
        for v in sans_retour[u]:
            entrants[v] += 1
    niv = {i: 0 for i in ids}
    file = [i for i in ids if entrants[i] == 0]
    while file:
        u = file.pop(0)
        for v in sans_retour[u]:
            niv[v] = max(niv[v], niv[u] + 1)
            entrants[v] -= 1
            if entrants[v] == 0:
                file.append(v)
    types = {b["id"]: b.get("type") for b in blocs}
    dernier = max(niv.values())
    for i in ids:
        if types[i] in ("entree", "consigne", "port_entree"):
            niv[i] = 0
        elif types[i] in ("sortie", "port_sortie") and dernier > 0:
            niv[i] = dernier
    for i in ids:  # une constante se pose juste avant son premier client
        if types[i] == "constante" and succ[i]:
            niv[i] = max(0, min(niv[v] for v in succ[i]) - 1)

    colonnes: dict[int, list[str]] = {}
    for i in ids:
        colonnes.setdefault(niv[i], []).append(i)
    rang = {i: float(pos[i]) for i in ids}
    for k, i in enumerate(ids):
        rang[i] = float(k)
    for _balayage in range(4):
        for sens in (pred, succ):
            for n in sorted(colonnes):
                col = colonnes[n]

                def bary(i, sens=sens):
                    voisins = [rang[v] for v in sens[i]]
                    return sum(voisins) / len(voisins) if voisins else rang[i]

                col.sort(key=lambda i: (bary(i), pos[i]))
                for r, i in enumerate(col):
                    rang[i] = float(r)
    place = {}
    for n, col in colonnes.items():
        for r, i in enumerate(col):
            place[i] = (20 + n * _PAS_X, 30 + r * _PAS_Y)
    for b in blocs:
        b["x"], b["y"] = place[b["id"]]


# ----------------------------------------------------------- depuis_element
def _sonder_recales(el, wb: float) -> dict[str, float]:
    """Les parametres que l'``initialize`` natif RECALCULE : ``{nom: valeur
    de depart}``. Observes sur une copie de l'element initialisee en un point
    d'essai (tension 1 pu, angle de repere non nul)."""
    from g2egrid.noyau.element import InitContext

    essai = copy.deepcopy(el)
    avant = dict(essai.pdict)
    ctx = InitContext(
        v={"t": complex(1.0, 0.0)},
        s={"t": complex(0.6, 0.1)},
        w=1.0,
        f0=wb / (2 * math.pi),
        theta_g=0.3,
    )
    try:
        essai.initialize(ctx)
    except Exception:  # noqa: BLE001 - une sonde qui echoue ne recale rien
        return {}
    return {
        k: float(avant[k])
        for k, v in essai.pdict.items()
        if np.isfinite(v) and np.isfinite(avant[k]) and v != avant[k]
    }


def _resoudre(g: sp.Expr, z: sp.Symbol) -> sp.Expr | None:
    """``z`` tire de ``g = 0`` quand ``g`` est lineaire en ``z`` (coefficient
    non nul), ``None`` sinon."""
    a = sp.diff(g, z)
    if a == 0 or a.has(z):
        return None
    return -g.xreplace({z: sp.Integer(0)}) / a


def depuis_element(
    classe_ou_nom: Any,
    parametres: dict[str, float] | None = None,
    detail: str = "blocs",
    nom: str | None = None,
    f0: float = 50.0,
) -> dict:
    """Le schema-bloc d'une classe NATIVE de la bibliotheque.

    ``classe_ou_nom`` : la classe (``GFL1``) ou son nom (``"GFL1"``) dans
    :mod:`g2egrid.elements`. ``parametres`` : valeurs qui remplacent les
    defauts (memes noms que ceux de la classe native, y compris les reglages
    de son constructeur : ``tau_CL``...). ``detail`` : ``"blocs"`` (un
    sous-systeme par equation, en blocs elementaires) ou ``"equations"`` (un
    bloc ``fonction`` par equation). ``nom`` : le nom du modele (celui de la
    classe sinon). ``f0`` : la frequence du reseau (Hz) pour laquelle les
    natifs regles automatiquement (PLL, boucles de courant) calculent leurs
    gains : 50 Hz par defaut, 60 Hz pour un reseau a 60 Hz.

    Seule une classe a UN port alternatif ``t`` est prise en charge ; un
    element a port continu (GFM3DC, MMC...) ou a plusieurs ports (branche,
    transformateur) est refuse, avec la raison. Voir le module pour ce qui
    est genere, ce qui ne l'est pas, et l'amorce de l'initialisation."""
    el = instance_native(classe_ou_nom, parametres, f0)
    return depuis_instance(el, detail, nom)


def instance_native(
    classe_ou_nom: Any, parametres: dict[str, float] | None = None, f0: float = 50.0
):
    """Une instance de la classe native (hors reseau), ``parametres`` donnes :
    ``kS = 1``, pulsation de base ``2 pi f0``."""
    cls = classe_native(classe_ou_nom)
    if cls is None:
        raise SchemaInvalide(f"classe native {classe_ou_nom!r} inconnue")
    _controler_classe(cls)
    parametres = dict(parametres or {})
    connus = set(cls.params) | set(inspect.signature(cls.__init__).parameters)
    for k, v in parametres.items():
        if k not in connus:
            raise SchemaInvalide(f"{cls.__name__} : parametre {k!r} inconnu")
        if (
            isinstance(v, bool)
            or not isinstance(v, (int, float))
            or not math.isfinite(v)
        ):
            raise SchemaInvalide(
                f"{cls.__name__} : parametre {k} = {v!r} (un nombre fini est attendu)"
            )
    el = _construire_natif(cls, parametres, 1.0, 2 * math.pi * float(f0))
    # Un parametre DONNE prime sur le reglage automatique du constructeur :
    # le GFM natif recalcule KpVL, KiVL, KpCL, KiCL depuis wb, et des gains
    # lus sur un element d'un reseau a 60 Hz seraient sinon remplaces par
    # ceux de f0 (14 % d'ecart de spectre sur wscc9_gfm).
    for k, v in (parametres or {}).items():
        if k in el.pdict and k not in ("kS", "wb"):
            el.set_param(k, float(v))
    return el


def _controler_classe(cls: type) -> None:
    if cls.ports_dc:
        raise SchemaInvalide(
            f"{cls.__name__} a un port continu ({', '.join(cls.ports_dc)}) : "
            "l'ecriture en blocs ne prend en charge que les elements a un seul "
            "port alternatif t"
        )
    if tuple(cls.ports) != ("t",):
        raise SchemaInvalide(
            f"{cls.__name__} a {len(cls.ports)} port(s) alternatif(s) "
            f"({', '.join(cls.ports) or 'aucun'}) : l'ecriture en blocs ne "
            "prend en charge que les elements a un seul port alternatif t"
        )
    if cls.__name__ not in CONSIGNES and len(cls.svc_inputs) < 2:
        n = len(cls.svc_inputs)
        raise SchemaInvalide(
            f"{cls.__name__} a {n} entree(s) de service "
            f"({', '.join(cls.svc_inputs) or 'aucune'}) : un element a blocs se "
            "cale au load flow par deux consignes (P et Q) ; ce modele ne "
            "s'ecrit pas en blocs en v1"
        )


def consignes_natives(
    nom_classe: str, entrees_svc: list[str], reduit: bool = False
) -> tuple[tuple[str, ...], dict[str, tuple[float, bool]]]:
    """``(calees, fixes)`` des entrees de service d'une classe native : voir
    :data:`CONSIGNES` (les deux premieres sont calees pour une classe
    absente ou une variante ``reduire``)."""
    calees, fixes = CONSIGNES.get(
        nom_classe,
        (tuple(entrees_svc[:2]), {n: (0.0, False) for n in entrees_svc[2:]}),
    )
    if (
        reduit
        or any(n not in entrees_svc for n in calees)
        or any(n not in entrees_svc for n in fixes)
    ):
        calees, fixes = (
            tuple(entrees_svc[:2]),
            {n: (0.0, False) for n in entrees_svc[2:]},
        )
    return tuple(calees), dict(fixes)


@dataclass
class _Dae:
    """Ce dont l'ecriture en blocs a besoin : un element natif, ou un modele a
    equations (:mod:`g2egrid.elements.blocs_equations`) lu comme un DAE."""

    nom_classe: str
    reduit: bool
    pdict: dict[str, float]
    params: set[str]
    etats: list[str]
    algebs: list[str]
    entrees_svc: list[str]
    sorties_svc: list[str]
    #: nombre de signaux de port en tete de ``h`` (les deux courants)
    n_port: int
    f: list[sp.Expr]
    g: list[sp.Expr]
    h: list[sp.Expr]
    symboles: dict[str, sp.Symbol]
    #: consignes calees, et ``{nom: (valeur, recalee)}`` des autres
    calees: tuple[str, ...]
    fixes: dict[str, tuple[float, bool]]
    #: parametres (ou consignes) que l'initialisation native recalcule
    recales: dict[str, float]


def depuis_instance(el, detail: str = "blocs", nom: str | None = None) -> dict:
    """Le schema-bloc d'un element DEJA construit (avec SES parametres, ceux
    d'une ligne de table par exemple). Meme resultat que
    :func:`depuis_element`, qui construit l'element puis l'appelle.

    ``detail="liste"`` rend le MODELE A EQUATIONS de l'element
    (:func:`g2egrid.elements.blocs_equations.equations_depuis_instance`) a la
    place d'un schema-bloc."""
    if detail == "liste":
        from g2egrid.elements.blocs_equations import equations_depuis_instance

        return equations_depuis_instance(el, nom=nom)
    if detail not in ("blocs", "equations"):
        raise SchemaInvalide(f"detail {detail!r} : 'blocs', 'equations' ou 'liste'")
    cls = type(el)
    _controler_classe(cls)
    nom_classe = cls.__name__
    reduit = hasattr(cls, "_classe_complete")  # une variante ``reduire(...)``
    pdict = {k: float(v) for k, v in el.pdict.items()}
    wb = float(pdict.get("wb", _WB))
    f, g, h, symboles = _equations(cls, pdict)
    entrees_svc = list(cls.svc_inputs)
    calees, fixes = consignes_natives(nom_classe, entrees_svc, reduit)
    dae = _Dae(
        nom_classe=nom_classe,
        reduit=reduit,
        pdict=pdict,
        params=set(cls.params),
        etats=list(cls.states),
        algebs=list(cls.algebs),
        entrees_svc=entrees_svc,
        sorties_svc=list(cls.svc_outputs),
        n_port=cls.n_port_signaux(),
        f=f,
        g=g,
        h=h,
        symboles=symboles,
        calees=tuple(calees),
        fixes=dict(fixes),
        recales=_sonder_recales(el, wb),
    )
    return schema_depuis_dae(dae, detail, nom)


def schema_depuis_dae(dae: _Dae, detail: str, nom: str | None) -> dict:
    """Le schema-bloc d'un DAE (voir :class:`_Dae`) : le coeur de
    :func:`depuis_instance`, que le modele a equations emploie aussi."""
    nom_classe, reduit, pdict = dae.nom_classe, dae.reduit, dae.pdict
    f, g, h, symboles = dae.f, dae.g, dae.h, dae.symboles
    params = set(dae.params)
    etats, algebs = list(dae.etats), list(dae.algebs)
    entrees_svc = list(dae.entrees_svc)

    pris: set[str] = set(etats) | set(algebs) | set(entrees_svc) | set(dae.params)

    def libre(base: str) -> str:
        base = _identifiant(base)
        bid, k = base, 1
        while bid in pris:
            k += 1
            bid = f"{base}_{k}"
        pris.add(bid)
        return bid

    blocs: list[dict] = []
    liens: list[dict] = []
    source: dict[str, tuple[str, int]] = {}  # variable native -> source du 1er niveau
    flat: dict[str, str] = {}  # variable native -> variable aplatie du schema

    # --- identifiants fixes d'avance (les equations se referencent entre elles)
    id_etat = {n: _identifiant(n) for n in etats}
    if len(set(id_etat.values())) != len(etats):
        raise SchemaInvalide(f"{nom_classe} : deux etats portent le meme identifiant")
    pris |= set(id_etat.values())
    for n in etats:
        source[n] = (id_etat[n], 0)
        flat[n] = f"y_{id_etat[n]}"

    # --- consignes
    calees, fixes = dae.calees, dae.fixes
    recales = dict(dae.recales)
    for n in entrees_svc:
        bid = libre(f"c_{n}")
        b = {"id": bid, "type": "consigne", "nom": n}
        if n not in calees:
            valeur, recalee = fixes.get(n, (0.0, False))
            b["calee"] = False
            b["valeur"] = valeur
            if recalee:
                recales[n] = valeur
        blocs.append(b)
        source[n] = (bid, 0)
        flat[n] = n

    # --- entrees du reseau (creees a la demande)
    def entree_reseau(sym: str) -> tuple[str, int]:
        if sym not in source:
            bid = libre(f"e_{_ENTREES[sym]}")
            blocs.append({"id": bid, "type": "entree", "signal": _ENTREES[sym]})
            source[sym] = (bid, 0)
        return source[sym]

    # --- les algebriques : explicites (resolues) ou implicites
    solutions: dict[str, sp.Expr | None] = {}
    appariement = _apparier(g, [symboles[a] for a in algebs])
    id_alg: dict[str, str] = {}
    for k, a in enumerate(algebs):
        eq = g[appariement[k]]
        sol = _resoudre(eq, symboles[a])
        solutions[a] = sol
        if sol is None:
            id_alg[a] = libre(f"z_{a}")
            source[a] = (id_alg[a], 0)
            flat[a] = f"y_{id_alg[a]}"
        else:
            id_alg[a] = libre(f"g_{a}")
            source[a] = (id_alg[a], 0)  # la sortie du sous-systeme

    def variable_premier_niveau(nom_var: str) -> tuple[str, int]:
        if nom_var in _ENTREES:
            return entree_reseau(nom_var)
        return source[nom_var]

    # une algebrique qui n'est qu'un autre signal (``a = b``) n'a pas de bloc :
    # elle est ce signal, pour tous ceux qui la lisent
    alias = {
        a: str(sol)
        for a, sol in solutions.items()
        if sol is not None and sol.is_Symbol and str(sol) not in params
    }
    for _ in range(len(alias)):
        for a, cible in list(alias.items()):
            if cible not in alias:
                source[a] = variable_premier_niveau(cible)
                del alias[a]

    def sous_systeme(sid: str, titre: str, e: sp.Expr) -> tuple[str, int] | None:
        """Le sous-systeme ``sid`` qui calcule ``e``, et ses liens d'entree ;
        rend ``(sid, 0)``. Une expression reduite a une variable est rendue
        telle quelle (aucun bloc)."""
        if e.is_Symbol and str(e) not in params:
            return variable_premier_niveau(str(e))
        entrees_ports: list[str] = []

        def variable(nom_var: str) -> tuple[str, int]:
            if nom_var not in entrees_ports:
                entrees_ports.append(nom_var)
            return (f"in_{nom_var}", 0)

        c = _Constructeur(params, variable, pris={"out"})
        sortie = c.noeud(e)
        ports = [
            {
                "id": f"in_{v}",
                "type": "port_entree",
                "index": k,
                "nom": v,
            }
            for k, v in enumerate(entrees_ports)
        ]
        c.liens.append(
            {
                "de": sortie[0],
                "vers": "out",
                "entree": 0,
                **({"sortie": sortie[1]} if sortie[1] else {}),
            }
        )
        interieur = {
            "blocs": ports
            + c.blocs
            + [{"id": "out", "type": "port_sortie", "index": 0, "nom": titre}],
            "liens": c.liens,
        }
        # les ports lus par un bloc (pas un lien) : on relie port -> consommateur
        blocs.append(
            {"id": sid, "type": "sous_systeme", "nom": titre, "schema": interieur}
        )
        for k, v in enumerate(entrees_ports):
            de, ks = variable_premier_niveau(v)
            lien = {"de": de, "vers": sid, "entree": k}
            if ks:
                lien["sortie"] = ks
            liens.append(lien)
        return (sid, 0), sortie

    def bloc_fonction(bid: str, e: sp.Expr) -> tuple[str, int]:
        c = _Constructeur(params, variable_premier_niveau, pris=set(pris))
        src = c.fonction(e, bid)
        blocs.extend(c.blocs)
        liens.extend(c.liens)
        return src

    # une variable resolue est rendue au premier niveau par la source de son
    # calcul ; les sous-systemes ne nous disent pas le bloc final : on le lit

    def calculer(sid: str, titre: str, e: sp.Expr, var: str | None) -> tuple[str, int]:
        """Calcule ``e`` : un bloc ``fonction`` (detail ``equations``) ou un
        sous-systeme (``blocs``) ; ``var`` : la variable native que l'on
        calcule (pour la correspondance de l'amorce)."""
        if detail == "equations":
            src = bloc_fonction(sid, e)
            if var is not None and src[0] == sid:
                flat[var] = f"y_{sid}"
            return src
        res = sous_systeme(sid, titre, e)
        if isinstance(res[0], tuple):
            (s_id, _k), interne = res
            if var is not None and interne[0] != "out":
                flat[var] = f"y_{s_id}__{interne[0]}"
            return (s_id, 0)
        return res

    # sous-systeme(...) rend soit une source (variable), soit (source, sortie)
    # --- construction des variables algebriques
    for a in algebs:
        if solutions[a] is None:
            k = algebs.index(a)
            eq = g[appariement[k]]
            _ecrire_implicite(
                blocs,
                liens,
                id_alg[a],
                eq,
                symboles[a],
                params,
                variable_premier_niveau,
            )
        else:
            calculer(id_alg[a], a, solutions[a], a)

    # --- integrateurs et equations d'etat
    id_eq = {n: libre(f"f_{n}") for n in etats}
    for k, n in enumerate(etats):
        # le facteur parametrique de l'equation (wb/L1, KiCL...) devient le
        # gain K de l'integrateur : la derivee n'est alors calculee que de
        # l'ecart, ce qui est lisible, et garde des coefficients moderes dans
        # les equations algebriques du modele (un KiCL de 6e5 y rendrait la
        # matrice mal conditionnee pour l'elimination des algebriques)
        coef, reste = _facteur_parametrique(f[k], params)
        src = calculer(id_eq[n], f"d{n}/dt", reste, None)
        blocs.append(
            {
                "id": id_etat[n],
                "type": "integrateur",
                "K": 1 if coef == 1 else _texte(coef),
            }
        )
        lien = {"de": src[0], "vers": id_etat[n], "entree": 0}
        if src[1]:
            lien["sortie"] = src[1]
        liens.append(lien)

    # --- sorties vers le reseau
    n_port = dae.n_port
    sorties = [("id", h[0]), ("iq", h[1])]
    svc = list(dae.sorties_svc)
    angle = next((k for k, s in enumerate(svc) if s in _SERVICE_ANGLE), None)
    vitesse = next((k for k, s in enumerate(svc) if s in _SERVICE_VITESSE), None)
    if angle is not None and vitesse is not None:
        sorties += [("theta", h[n_port + angle]), ("w", h[n_port + vitesse])]
    for signal, e in sorties:
        titre = {"id": "id injecte", "iq": "iq injecte"}.get(signal, signal)
        sid = libre(f"h_{signal}")
        src = calculer(sid, titre, e, None)
        oid = libre(f"o_{signal}")
        blocs.append({"id": oid, "type": "sortie", "signal": signal})
        lien = {"de": src[0], "vers": oid, "entree": 0}
        if src[1]:
            lien["sortie"] = src[1]
        liens.append(lien)

    parametres = {
        k: v for k, v in pdict.items() if k not in ("kS", "wb") and np.isfinite(v)
    }
    schema: dict[str, Any] = {
        "nom": _identifiant(nom) if nom else _identifiant(nom_classe),
        "parametres": parametres,
        "blocs": blocs,
        "liens": liens,
    }
    if not reduit:
        correspondance = {}
        for n in etats:
            correspondance[f"x_{id_etat[n]}"] = n
            correspondance[f"y_{id_etat[n]}"] = n
        for n in algebs:
            if n in flat:
                correspondance[flat[n]] = n
        for n in entrees_svc:
            correspondance[n] = n
        schema["amorce"] = {
            "classe": nom_classe,
            "correspondance": correspondance,
            "recales": {
                k: v for k, v in recales.items() if k in parametres or k in entrees_svc
            },
        }
    return disposer(schema)


def _facteur_parametrique(e: sp.Expr, params: set[str]) -> tuple[sp.Expr, sp.Expr]:
    """``e = coef . reste`` ou ``coef`` est un produit de facteurs qui ne
    dependent que des parametres (au moins un symbole) ; ``(1, e)`` sinon."""
    if not e.is_Mul:
        return sp.Integer(1), e
    coef = [
        a
        for a in e.args
        if a.free_symbols
        and all(str(x) in params for x in a.free_symbols)
        and not a.has(sp.Min, sp.Max, sp.Piecewise)
    ]
    if not coef:
        return sp.Integer(1), e
    reste = sp.Mul(*[a for a in e.args if a not in coef])
    return sp.Mul(*coef), reste


def _apparier(g: list[sp.Expr], z: list[sp.Symbol]) -> list[int]:
    """Pour chaque variable algebrique ``k``, l'indice de l'equation qui la
    definit : ``g[k]`` quand elle la contient (le cas des classes de la
    bibliotheque), un appariement structurel sinon."""
    n = len(z)
    if all(g[k].has(z[k]) for k in range(n)):
        return list(range(n))
    from scipy.optimize import linear_sum_assignment

    cout = np.array(
        [[0.0 if g[i].has(z[j]) else 1e6 for j in range(n)] for i in range(n)]
    )
    _lignes, colonnes = linear_sum_assignment(cout.T)
    return [int(i) for i in colonnes]


def _ecrire_implicite(blocs, liens, bid, eq, z, params, variable) -> None:
    """Le bloc ``implicite`` qui definit ``z`` par ``eq = 0`` : ses entrees
    sont les autres variables de ``eq``, ``y`` est ``z``."""
    variables = [
        s for s in sorted(eq.free_symbols, key=str) if str(s) not in params and s != z
    ]
    u = {s: sp.Symbol(f"u{k + 1}") for k, s in enumerate(variables)}
    texte = _texte(eq.xreplace(u | {z: sp.Symbol("y")}))
    blocs.append({"id": bid, "type": "implicite", "expr": texte, "n": len(variables)})
    for k, s in enumerate(variables):
        de, ks = variable(str(s))
        lien = {"de": de, "vers": bid, "entree": k}
        if ks:
            lien["sortie"] = ks
        liens.append(lien)
