"""Sous-systemes d'un schema-bloc : aplatir, regrouper, degrouper.

Un SOUS-SYSTEME est un bloc de type ``sous_systeme`` qui contient un schema
(``{"blocs": [...], "liens": [...]}``, de la meme forme que le schema de
premier niveau, positions ``x``, ``y`` comprises). A l'interieur, les blocs
``port_entree`` et ``port_sortie`` (champs ``index`` a partir de 0 et ``nom``)
font l'interface : le sous-systeme a autant d'entrees que de ``port_entree``
et autant de sorties que de ``port_sortie``, dans l'ordre des ``index``,
nommees par leurs ``nom``. L'imbrication est admise jusqu'a
:data:`PROFONDEUR_MAX` niveaux.

Les sous-systemes sont de l'ORGANISATION, pas de la physique : il n'y a ni
masque de parametres (les parametres du modele sont visibles partout) ni
portee des noms. Le compilateur ne les voit jamais : :func:`aplatir` les
dissout d'abord, chaque bloc interieur etant renomme
``<id_sous_systeme>__<id_interieur>`` (recursivement), les ports retires et
les liens raccordes. ``compiler(s)`` et ``compiler(aplatir(s))`` sont donc la
meme compilation, par construction. C'est pourquoi les identifiants saisis ne
contiennent pas ``__`` : le double souligne est le separateur du renommage
(une fois aplati, un schema n'a plus de structure et ses ``__`` sont admis).

:func:`regrouper` et :func:`degrouper` sont les deux operations d'edition :
l'interface les appelle (``/api/blocs/regrouper``, ``/degrouper``) au lieu
d'en recopier les regles en JavaScript. Elles travaillent sur une COPIE,
n'altere jamais le schema recu, et acceptent ``chemin`` (la liste des
identifiants de sous-systemes qui mene au niveau a editer) pour regrouper ou
degrouper a l'interieur d'un sous-systeme.

Convention des listes : les liens sont ecrits comme le fait
:func:`g2egrid.elements.blocs.schema_gfl` (``entree`` toujours, ``sortie`` si
elle est non nulle). Regrouper puis degrouper rend le schema de depart, aux
positions pres et -- quand les blocs regroupes se suivent dans la liste -- a
l'ordre des listes pres ; sinon seul l'ordre des blocs et des liens peut
changer, jamais leur contenu.
"""

from __future__ import annotations

import copy
import re
from typing import Any

__all__ = [
    "BLOC_H",
    "BLOC_L",
    "PROFONDEUR_MAX",
    "TYPES_STRUCTURE",
    "SchemaInvalide",
    "aplatir",
    "degrouper",
    "regrouper",
    "verifier_forme",
]


class SchemaInvalide(ValueError):
    """Un schema qu'on ne sait pas compiler, avec la raison."""


#: dimensions (px) de la boite d'un bloc, communes a l'interface et aux
#: dispositions automatiques : deux blocs ne se recouvrent pas si leurs boites
#: ``BLOC_L x BLOC_H`` sont disjointes
BLOC_L = 92
BLOC_H = 50
#: profondeur maximale d'imbrication des sous-systemes
PROFONDEUR_MAX = 8
#: les types de blocs de structure
TYPES_STRUCTURE = ("sous_systeme", "port_entree", "port_sortie")
#: pas de la grille de placement des ports et ecart au bloc le plus proche
_PAS_Y = 70
_ECART_X = 48


def _types() -> dict[str, dict[str, Any]]:
    """``blocs.TYPES`` (import tardif : ``blocs`` importe ce module)."""
    from g2egrid.elements import blocs

    return blocs.TYPES


def _lien(de: str, vers: str, entree: int = 0, sortie: int = 0) -> dict:
    lien = {"de": de, "vers": vers, "entree": entree}
    if sortie:
        lien["sortie"] = sortie
    return lien


def _id_valide(bid: Any) -> bool:
    return isinstance(bid, str) and bid.isidentifier() and not bid.startswith("_")


def verifier_forme(schema: Any, ou: str = "le schema", profondeur: int = 0) -> None:
    """La FORME d'un schema-bloc, avant toute lecture : un objet, dont
    ``parametres`` est un objet, ``blocs`` et ``liens`` des listes d'objets,
    les numeros d'entree et de sortie des liens des entiers -- a tous les
    niveaux. Un schema vient d'un fichier ou d'une requete : mal forme, il
    doit etre REFUSE avec la raison (:class:`SchemaInvalide`), pas faire
    echouer le code qui le lit (``'str' object has no attribute 'get'``)."""
    if not isinstance(schema, dict):
        raise SchemaInvalide(
            f"{ou} : un schema est un objet {{nom, parametres, blocs, liens}}"
        )
    if profondeur > PROFONDEUR_MAX:
        raise SchemaInvalide(
            f"sous-systemes imbriques sur plus de {PROFONDEUR_MAX} niveaux ({ou})"
        )
    p = schema.get("parametres")
    if p is not None and not isinstance(p, dict):
        raise SchemaInvalide(f"{ou} : parametres est un objet {{nom: valeur}}")
    for cle in ("blocs", "liens"):
        v = schema.get(cle)
        if v is not None and not isinstance(v, list):
            raise SchemaInvalide(f"{ou} : {cle} est une liste")
    for b in schema.get("blocs") or []:
        if not isinstance(b, dict):
            raise SchemaInvalide(
                f"{ou} : chaque bloc est un objet {{id, type, ...}} (recu : {b!r})"
            )
        if b.get("type") == "sous_systeme" and b.get("schema") is not None:
            if not isinstance(b["schema"], dict):
                raise SchemaInvalide(f"sous-systeme {b.get('id')!r} : schema manquant")
            verifier_forme(b["schema"], f"sous-systeme {b.get('id')!r}", profondeur + 1)
    for lien in schema.get("liens") or []:
        if not isinstance(lien, dict):
            raise SchemaInvalide(
                f"{ou} : chaque lien est un objet {{de, vers, entree, sortie}} "
                f"(recu : {lien!r})"
            )
        for cle in ("entree", "sortie"):
            v = lien.get(cle, 0)
            try:
                if isinstance(v, bool):
                    raise ValueError  # noqa: TRY004 -- refuse comme un texte
                int(v)
            except (TypeError, ValueError):
                raise SchemaInvalide(
                    f"lien {lien.get('de')} -> {lien.get('vers')} : {cle} = {v!r}, "
                    "un numero (entier a partir de 0) est attendu"
                ) from None


def contient_structure(schema: dict) -> bool:
    """Le schema (a n'importe quel niveau) a-t-il un bloc de structure ?"""
    return any(b.get("type") in TYPES_STRUCTURE for b in schema.get("blocs") or [])


def _ports(schema: dict, type_port: str, ou: str) -> list[dict]:
    """Les ports ``type_port`` d'un niveau, tries par ``index``."""
    out = [b for b in schema.get("blocs") or [] if b.get("type") == type_port]
    vus: set[int] = set()
    for b in out:
        k = b.get("index")
        if isinstance(k, bool) or not isinstance(k, int) or k < 0:
            raise SchemaInvalide(
                f"{ou} : le port {b.get('id')!r} a un index {k!r} "
                "(un entier a partir de 0)"
            )
        if k in vus:
            raise SchemaInvalide(f"{ou} : deux ports {type_port} d'index {k}")
        vus.add(k)
    return sorted(out, key=lambda b: b["index"])


def noms_entrees(bloc: dict) -> list[str]:
    """Les noms des entrees d'un bloc ``sous_systeme`` (ordre des ``index``)."""
    sch = bloc.get("schema") or {}
    return [
        str(p.get("nom") or f"e{k + 1}")
        for k, p in enumerate(
            _ports(sch, "port_entree", f"sous-systeme {bloc.get('id')}")
        )
    ]


def noms_sorties(bloc: dict) -> list[str]:
    """Les noms des sorties d'un bloc ``sous_systeme`` (ordre des ``index``)."""
    sch = bloc.get("schema") or {}
    return [
        str(p.get("nom") or f"s{k + 1}")
        for k, p in enumerate(
            _ports(sch, "port_sortie", f"sous-systeme {bloc.get('id')}")
        )
    ]


# ------------------------------------------------------------------ aplatir
class _Niveau:
    """Un niveau du schema (le premier, ou l'interieur d'un sous-systeme)."""

    def __init__(self, schema, prefixe, parent, sid, chemin, profondeur):
        if profondeur > PROFONDEUR_MAX:
            raise SchemaInvalide(
                f"sous-systemes imbriques sur plus de {PROFONDEUR_MAX} niveaux "
                f"({chemin})"
            )
        ou = chemin or "le schema"
        self.prefixe, self.parent, self.sid, self.chemin = prefixe, parent, sid, chemin
        self.blocs: dict[str, dict] = {}
        for b in schema.get("blocs") or []:
            bid = b.get("id")
            if not _id_valide(bid):
                raise SchemaInvalide(f"{ou} : identifiant de bloc {bid!r} invalide")
            if "__" in bid:
                raise SchemaInvalide(
                    f"{ou} : identifiant {bid!r} : le double souligne (__) est "
                    "reserve au renommage des sous-systemes"
                )
            if bid in self.blocs:
                raise SchemaInvalide(f"{ou} : identifiant de bloc {bid!r} en double")
            if b.get("type") not in _types():
                raise SchemaInvalide(f"bloc {bid!r} : type {b.get('type')!r} inconnu")
            self.blocs[bid] = b
        if parent is None and any(
            b["type"] in ("port_entree", "port_sortie") for b in self.blocs.values()
        ):
            bad = next(
                b["id"]
                for b in self.blocs.values()
                if b["type"] in ("port_entree", "port_sortie")
            )
            raise SchemaInvalide(
                f"bloc {bad!r} : port_entree et port_sortie ne sont admis que "
                "dans un sous-systeme"
            )
        self.liens: list[dict] = []
        self.vers: dict[tuple[str, int], tuple[str, int]] = {}
        for lien in schema.get("liens") or []:
            de, vers = str(lien.get("de")), str(lien.get("vers"))
            if de not in self.blocs or vers not in self.blocs:
                raise SchemaInvalide(f"lien {de} -> {vers} : bloc inconnu")
            k, ks = int(lien.get("entree", 0)), int(lien.get("sortie", 0))
            cle = (vers, k)
            if cle in self.vers:
                raise SchemaInvalide(f"bloc {vers} : entree {k + 1} deja reliee")
            self.vers[cle] = (de, ks)
            self.liens.append(dict(de=de, vers=vers, entree=k, sortie=ks))
        self.entrees = _ports(schema, "port_entree", ou)
        self.sorties = _ports(schema, "port_sortie", ou)
        self.enfants: dict[str, _Niveau] = {}
        for bid, b in self.blocs.items():
            if b["type"] != "sous_systeme":
                continue
            interieur = b.get("schema")
            if not isinstance(interieur, dict):
                raise SchemaInvalide(f"sous-systeme {bid!r} : schema manquant")
            ch = f"{chemin} > {bid}" if chemin else bid
            self.enfants[bid] = _Niveau(
                interieur, f"{prefixe}{bid}__", self, bid, ch, profondeur + 1
            )
        self._verifier_liens()

    def _verifier_liens(self) -> None:
        for lien in self.liens:
            de, vers, k, ks = lien["de"], lien["vers"], lien["entree"], lien["sortie"]
            tv, td = self.blocs[vers]["type"], self.blocs[de]["type"]
            if tv == "sous_systeme" and not 0 <= k < len(self.enfants[vers].entrees):
                raise SchemaInvalide(f"lien {de} -> {vers} : entree {k} inexistante")
            if tv == "port_sortie" and k != 0:
                raise SchemaInvalide(f"lien {de} -> {vers} : entree {k} inexistante")
            if tv == "port_entree":
                raise SchemaInvalide(
                    f"lien {de} -> {vers} : un port_entree n'a pas d'entree"
                )
            if td == "sous_systeme" and not 0 <= ks < len(self.enfants[de].sorties):
                raise SchemaInvalide(f"lien {de} -> {vers} : sortie {ks} inexistante")
            if td == "port_entree" and ks != 0:
                raise SchemaInvalide(f"lien {de} -> {vers} : sortie {ks} inexistante")
            if td == "port_sortie":
                raise SchemaInvalide(
                    f"lien depuis {de} : un port_sortie n'a pas de sortie"
                )

    def resoudre(self, bid: str, sortie: int, vus: frozenset = frozenset()):
        """La source REELLE ``(id aplati, sortie)`` du signal ``(bid, sortie)``
        de ce niveau, en traversant les sous-systemes et les ports ; ``None``
        si le signal vient d'un port d'entree que rien n'alimente."""
        b = self.blocs[bid]
        t = b["type"]
        if t == "sous_systeme":
            enfant = self.enfants[bid]
            port = enfant.sorties[sortie]
            cle = (id(enfant), port["id"])
            if cle in vus:
                raise SchemaInvalide(
                    f"boucle de ports sans bloc dans le sous-systeme {enfant.chemin}"
                )
            src = enfant.vers.get((port["id"], 0))
            if src is None:
                nom = port.get("nom") or port["id"]
                raise SchemaInvalide(
                    f"sous-systeme {enfant.chemin} : la sortie {nom!r} n'est reliee a rien"
                )
            return enfant.resoudre(src[0], src[1], vus | {cle})
        if t == "port_entree":
            k = next(j for j, p in enumerate(self.entrees) if p["id"] == bid)
            cle = (id(self), bid)
            if cle in vus:
                raise SchemaInvalide(
                    f"boucle de ports sans bloc dans le sous-systeme {self.chemin}"
                )
            src = self.parent.vers.get((self.sid, k))
            if src is None:
                return None
            return self.parent.resoudre(src[0], src[1], vus | {cle})
        return (self.prefixe + bid, sortie)

    def emettre(self, blocs: list[dict], liens: list[dict]) -> None:
        """Verse ce niveau, aplati, dans ``blocs`` et ``liens``."""
        for bid, b in self.blocs.items():
            if b["type"] == "sous_systeme":
                self.enfants[bid].emettre(blocs, liens)
            elif b["type"] not in ("port_entree", "port_sortie"):
                copie = copy.deepcopy(b)
                copie["id"] = self.prefixe + bid
                blocs.append(copie)
        for lien in self.liens:
            if self.blocs[lien["vers"]]["type"] in TYPES_STRUCTURE:
                continue  # raccorde depuis l'autre cote, par ``resoudre``
            src = self.resoudre(lien["de"], lien["sortie"])
            if src is None:
                continue
            liens.append(
                _lien(src[0], self.prefixe + lien["vers"], lien["entree"], src[1])
            )


def aplatir(schema: dict) -> dict:
    """Le schema sans sous-systeme : les blocs interieurs renommes
    ``<id_sous_systeme>__<id_interieur>`` (recursivement), les ports retires,
    les liens raccordes.

    Un schema sans structure est rendu tel quel (copie). Le schema recu n'est
    jamais modifie. Leve :class:`SchemaInvalide` pour un ``port_entree`` ou
    un ``port_sortie`` hors d'un sous-systeme, un identifiant qui contient
    ``__``, un sous-systeme dont une sortie ne mene a rien, une imbrication de
    plus de :data:`PROFONDEUR_MAX` niveaux."""
    verifier_forme(schema)
    if not contient_structure(schema):
        return copy.deepcopy(schema)
    racine = _Niveau(schema, "", None, None, "", 0)
    blocs: list[dict] = []
    liens: list[dict] = []
    racine.emettre(blocs, liens)
    out = {
        k: copy.deepcopy(v) for k, v in schema.items() if k not in ("blocs", "liens")
    }
    out["blocs"], out["liens"] = blocs, liens
    return out


# ----------------------------------------------------------------- edition
def _acceder(schema: dict, chemin) -> dict:
    """Le niveau (dict ``blocs``/``liens``) que designe ``chemin`` dans
    ``schema`` (modifie en place : l'appelant a deja copie)."""
    niveau = schema
    for sid in chemin or ():
        b = next(
            (x for x in niveau.get("blocs") or [] if x.get("id") == sid),
            None,
        )
        if b is None or b.get("type") != "sous_systeme":
            raise SchemaInvalide(f"chemin : {sid!r} n'est pas un sous-systeme")
        niveau = b["schema"]
    niveau.setdefault("blocs", [])
    niveau.setdefault("liens", [])
    return niveau


def _slug(nom: str) -> str:
    s = re.sub(r"[^0-9a-zA-Z]+", "_", str(nom).strip().lower()).strip("_")
    s = re.sub(r"_+", "_", s)
    if not s or s[0].isdigit():
        s = "ss_" + s if s else "ss"
    return s


def _id_libre(base: str, pris: set[str]) -> str:
    base = base.replace("__", "_")
    if base not in pris:
        return base
    k = 2
    while f"{base}_{k}" in pris:
        k += 1
    return f"{base}_{k}"


def _nom_sortie(bloc: dict, k: int) -> str:
    """Le nom d'un signal : celui du bloc, suffixe du nom de la sortie pour un
    bloc a plusieurs sorties."""
    from g2egrid.elements import blocs

    noms = blocs._noms_sorties(bloc)
    bid = str(bloc["id"])
    return bid if len(noms) <= 1 else f"{bid}_{noms[k]}"


def _boite(blocs: list[dict]) -> tuple[float, float, float, float] | None:
    """``(xmin, xmax, ymin, ymax)`` des positions des blocs, si tous en ont."""
    if not blocs or any("x" not in b or "y" not in b for b in blocs):
        return None
    xs = [b["x"] for b in blocs]
    ys = [b["y"] for b in blocs]
    return min(xs), max(xs), min(ys), max(ys)


def _demi(v: float) -> float | int:
    """``v`` divise par 2, entier quand c'est possible (JSON propre)."""
    r = v / 2
    return int(r) if float(r).is_integer() else r


def regrouper(
    schema: dict,
    ids: list[str],
    nom: str,
    id: str | None = None,
    chemin: list[str] | tuple[str, ...] = (),
) -> dict:
    """Range les blocs ``ids`` (d'un meme niveau) dans un nouveau sous-systeme.

    * un ``port_entree`` par SIGNAL qui entre (un seul port par sortie
      source, meme si elle alimente plusieurs blocs interieurs), un
      ``port_sortie`` par signal interieur utilise dehors ; les liens sont
      raccordes de part et d'autre ;
    * le sous-systeme se place au centre des blocs deplaces, qui gardent leurs
      positions relatives (le schema interieur est recale a partir de
      ``(160, 30)``, les ports de part et d'autre) ;
    * ``id`` : l'identifiant du sous-systeme (derive de ``nom`` si absent) ;
    * ``chemin`` : les identifiants de sous-systemes qui menent au niveau a
      editer (le premier niveau si vide).

    Les equations du modele ne changent pas : le schema aplati est le meme,
    aux noms des variables pres (``<id_sous_systeme>__...``)."""
    verifier_forme(schema)
    sch = copy.deepcopy(schema)
    niveau = _acceder(sch, chemin)
    blocs, liens = niveau["blocs"], niveau["liens"]
    aplatis = sorted(str(b.get("id")) for b in blocs if "__" in str(b.get("id")))
    if aplatis:
        # un schema APLATI porte des identifiants <sous_systeme>__<bloc> : un
        # sous-systeme de plus donnerait un schema que plus rien ne relit
        raise SchemaInvalide(
            f"regrouper : ce niveau porte des identifiants en __ ({aplatis[0]}...) : "
            "c'est un schema aplati ; regrouper le schema d'origine, ou renommer "
            "ces blocs"
        )
    ids = [str(i) for i in (ids or [])]
    if not ids:
        raise SchemaInvalide("regrouper : aucun bloc a regrouper")
    if len(set(ids)) != len(ids):
        raise SchemaInvalide("regrouper : un bloc est cite deux fois")
    par_id = {str(b.get("id")): b for b in blocs}
    for i in ids:
        if i not in par_id:
            raise SchemaInvalide(f"regrouper : bloc {i!r} inconnu")
        if par_id[i].get("type") in ("port_entree", "port_sortie"):
            raise SchemaInvalide(
                f"regrouper : {i!r} est un port du sous-systeme, il ne se regroupe pas"
            )
    nom = str(nom).strip() if nom is not None else ""
    if not nom:
        raise SchemaInvalide("regrouper : le sous-systeme veut un nom")
    pris = set(par_id)
    if id is None:
        sid = _id_libre(_slug(nom), pris)
    else:
        sid = str(id)
        if not _id_valide(sid) or "__" in sid:
            raise SchemaInvalide(f"regrouper : identifiant {sid!r} invalide")
        if sid in pris:
            raise SchemaInvalide(f"regrouper : l'identifiant {sid!r} existe deja")
    dedans = set(ids)
    deplaces = [b for b in blocs if str(b.get("id")) in dedans]

    entrees_ext: dict[tuple[str, int], int] = {}  # source externe -> index du port
    sorties_int: dict[tuple[str, int], int] = {}  # source interne -> index du port
    liens_int: list[dict] = []
    ports_e: list[tuple[str, int, str]] = []  # (id du port, index, nom)
    ports_s: list[tuple[str, int, str]] = []
    pris_int = set(dedans)

    def port(prefixe, k):
        pid = _id_libre(f"{prefixe}{k + 1}", pris_int)
        pris_int.add(pid)
        return pid

    nouveaux_ext: dict[tuple[str, int], dict] = {}
    for lien in liens:
        de, vers = str(lien["de"]), str(lien["vers"])
        k, ks = int(lien.get("entree", 0)), int(lien.get("sortie", 0))
        d_in, v_in = de in dedans, vers in dedans
        if d_in and v_in:
            liens_int.append(_lien(de, vers, k, ks))
        elif v_in:  # un signal entre
            cle = (de, ks)
            if cle not in entrees_ext:
                idx = len(entrees_ext)
                entrees_ext[cle] = idx
                ports_e.append((port("e", idx), idx, _nom_sortie(par_id[de], ks)))
                nouveaux_ext[cle] = _lien(de, sid, idx, ks)
            pid = ports_e[entrees_ext[cle]][0]
            liens_int.append(_lien(pid, vers, k))
        elif d_in:  # un signal sort
            cle = (de, ks)
            if cle not in sorties_int:
                idx = len(sorties_int)
                sorties_int[cle] = idx
                ports_s.append((port("s", idx), idx, _nom_sortie(par_id[de], ks)))
                liens_int.append(_lien(de, ports_s[idx][0], 0, ks))

    boite = _boite(deplaces)
    interieur_blocs = copy.deepcopy(deplaces)
    bloc_ss: dict[str, Any] = {"id": sid, "type": "sous_systeme", "nom": nom}
    if boite is not None:
        xmin, xmax, ymin, ymax = boite
        dx, dy = 160 - xmin, 30 - ymin
        for b in interieur_blocs:
            b["x"] += dx
            b["y"] += dy
        bloc_ss["x"] = _demi(xmin + xmax)
        bloc_ss["y"] = _demi(ymin + ymax)
    ports_blocs: list[dict] = []
    for pid, idx, pnom in ports_e:
        p = {"id": pid, "type": "port_entree", "index": idx, "nom": pnom}
        if boite is not None:
            p["x"], p["y"] = 20, 30 + idx * _PAS_Y
        ports_blocs.append(p)
    sorties_blocs: list[dict] = []
    for pid, idx, pnom in ports_s:
        p = {"id": pid, "type": "port_sortie", "index": idx, "nom": pnom}
        if boite is not None:
            p["x"] = xmax + dx + BLOC_L + _ECART_X
            p["y"] = 30 + idx * _PAS_Y
        sorties_blocs.append(p)
    bloc_ss["schema"] = {
        "blocs": ports_blocs + interieur_blocs + sorties_blocs,
        "liens": liens_int,
    }

    # la liste des blocs : le sous-systeme a la place du premier bloc deplace
    premier = next(i for i, b in enumerate(blocs) if str(b.get("id")) in dedans)
    nouveaux = []
    for i, b in enumerate(blocs):
        if i == premier:
            nouveaux.append(bloc_ss)
        if str(b.get("id")) not in dedans:
            nouveaux.append(b)
    niveau["blocs"] = nouveaux
    # les liens : ceux du dehors gardent leur ordre relatif
    ordre = []
    vus_ext = set()
    for lien in liens:
        de, vers = str(lien["de"]), str(lien["vers"])
        k, ks = int(lien.get("entree", 0)), int(lien.get("sortie", 0))
        if de in dedans and vers in dedans:
            continue
        if vers in dedans:
            cle = (de, ks)
            if cle in vus_ext:
                continue
            vus_ext.add(cle)
            ordre.append(nouveaux_ext[cle])
        elif de in dedans:
            ordre.append(_lien(sid, vers, k, sorties_int[(de, ks)]))
        else:
            ordre.append(_copie_lien(lien))
    niveau["liens"] = ordre
    return sch


def _copie_lien(lien: dict) -> dict:
    """Un lien dans la forme canonique (``entree`` toujours, ``sortie`` si
    elle est non nulle), les autres champs conserves."""
    out = copy.deepcopy(lien)
    out.setdefault("entree", 0)
    if not out.get("sortie"):
        out.pop("sortie", None)
    return out


def degrouper(schema: dict, id: str, chemin: list[str] | tuple[str, ...] = ()) -> dict:
    """L'inverse de :func:`regrouper` : dissout le sous-systeme ``id`` dans son
    niveau. Les blocs interieurs reprennent leur identifiant (prefixe
    ``<id>_`` repete jusqu'a etre libre si un autre bloc le porte deja), les
    ports disparaissent et les liens se raccordent ; le schema interieur est
    recentre sur la position du sous-systeme. Les sous-systemes imbriques
    restent des sous-systemes."""
    verifier_forme(schema)
    sch = copy.deepcopy(schema)
    niveau = _acceder(sch, chemin)
    blocs, liens = niveau["blocs"], niveau["liens"]
    ss = next((b for b in blocs if b.get("id") == id), None)
    if ss is None or ss.get("type") != "sous_systeme":
        raise SchemaInvalide(f"degrouper : {id!r} n'est pas un sous-systeme")
    inter = ss.get("schema") or {}
    i_blocs = inter.get("blocs") or []
    i_liens = [_copie_lien(lien) for lien in inter.get("liens") or []]
    ports_e = _ports(inter, "port_entree", f"sous-systeme {id}")
    ports_s = _ports(inter, "port_sortie", f"sous-systeme {id}")
    id_port_e = {p["id"]: k for k, p in enumerate(ports_e)}
    id_port_s = {p["id"]: k for k, p in enumerate(ports_s)}
    reels = [b for b in i_blocs if b.get("type") not in ("port_entree", "port_sortie")]

    # les identifiants : libres dans le niveau (hors le sous-systeme lui-meme)
    pris = {str(b.get("id")) for b in blocs} - {id}
    renom: dict[str, str] = {}
    for b in reels:
        bid = str(b["id"])
        nouveau = bid
        while nouveau in pris:
            nouveau = f"{id}_{nouveau}"
        pris.add(nouveau)
        renom[bid] = nouveau

    # position : le centre du schema interieur sur le centre du sous-systeme
    boite = _boite(reels)
    dx = dy = 0
    if boite is not None and "x" in ss and "y" in ss:
        xmin, xmax, ymin, ymax = boite
        dx = ss["x"] - _demi(xmin + xmax)
        dy = ss["y"] - _demi(ymin + ymax)
    sortis = []
    for b in reels:
        c = copy.deepcopy(b)
        c["id"] = renom[str(b["id"])]
        if boite is not None and "x" in ss and "y" in ss:
            c["x"] = b["x"] + dx
            c["y"] = b["y"] + dy
        sortis.append(c)

    # sources exterieures de chaque entree du sous-systeme
    ext_vers_ss: dict[int, tuple[str, int]] = {}
    for lien in liens:
        if str(lien["vers"]) == id:
            ext_vers_ss[int(lien.get("entree", 0))] = (
                str(lien["de"]),
                int(lien.get("sortie", 0)),
            )
    interne_vers: dict[tuple[str, int], tuple[str, int]] = {
        (str(lien["vers"]), lien["entree"]): (str(lien["de"]), lien.get("sortie", 0))
        for lien in i_liens
    }

    def source_interne(bid: str, ks: int):
        """``(id, sortie)`` du signal interieur, ``None`` s'il vient d'un port
        d'entree qu'aucun lien n'alimente ; traverse les ports d'entree."""
        if bid in id_port_e:
            return ext_vers_ss.get(id_port_e[bid])
        return (renom[bid], ks)

    nouveaux_liens: list[dict] = []
    ou_inserer: int | None = None  # les liens interieurs, a la place du premier
    for lien in liens:
        de, vers = str(lien["de"]), str(lien["vers"])
        if ou_inserer is None and id in (de, vers):
            ou_inserer = len(nouveaux_liens)
        if vers == id:
            continue  # remplace par les liens interieurs
        k, ks = int(lien.get("entree", 0)), int(lien.get("sortie", 0))
        if de == id:
            if not 0 <= ks < len(ports_s):
                raise SchemaInvalide(f"degrouper : sortie {ks} de {id!r} inexistante")
            src = interne_vers.get((ports_s[ks]["id"], 0))
            if src is None:
                continue
            res = source_interne(src[0], src[1])
            if res is not None:
                nouveaux_liens.append(_lien(res[0], vers, k, res[1]))
        else:
            nouveaux_liens.append(_copie_lien(lien))
    internes = []
    for lien in i_liens:
        de, vers = lien["de"], lien["vers"]
        if vers in id_port_s or vers in id_port_e:
            continue
        res = source_interne(de, lien.get("sortie", 0))
        if res is not None:
            internes.append(_lien(res[0], renom[vers], lien["entree"], res[1]))
    k = len(nouveaux_liens) if ou_inserer is None else ou_inserer
    nouveaux_liens[k:k] = internes

    pos = next(i for i, b in enumerate(blocs) if b.get("id") == id)
    niveau["blocs"] = blocs[:pos] + sortis + blocs[pos + 1 :]
    niveau["liens"] = nouveaux_liens
    return sch
