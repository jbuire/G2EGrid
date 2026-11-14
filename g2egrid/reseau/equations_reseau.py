"""Les EQUATIONS DYNAMIQUES du reseau, par element et en liste globale, et leur
edition (FR-052).

``dx/dt = f(x, u)`` de chaque equipement est lisible et MODIFIABLE :

* :func:`equations_reseau` : la liste globale, une ligne par equation de tous
  les elements assembles (equipements, branches, noeuds), avec ce qui se peut
  modifier et pourquoi le reste ne se peut pas ;
* :func:`modele_element` : le modele a equations
  (:mod:`g2egrid.elements.blocs_equations`) d'UN element ;
* :func:`modifier_equation`, :func:`ajouter_equation`,
  :func:`supprimer_equation`, :func:`modifier_parametre` : l'edition, qui
  rend de NOUVELLES tables (les tables recues ne bougent jamais, un refus non
  plus) ;
* :func:`revenir_au_natif` : l'element retrouve sa ligne d'origine.

COMMENT UN ELEMENT NATIF EST MODIFIE -- il est d'abord REMPLACE par son modele
a equations, ecrit avec SES parametres (ceux de la ligne de table et de
``PARAMS_PAR_BUS``), range dans ``MODELES_BLOCS`` sous le nom libre
``<classe>_<element>`` et pose dans ``T_Blocs`` a sa place (meme bus, meme
``Sn_VA``) ; l'equation est ensuite remplacee dans ce modele. La ligne native
est retiree de sa table ; le modele garde en ``"origine"`` ce qu'il faut pour
la lui rendre (:func:`revenir_au_natif`). Un modele partage par plusieurs
elements est d'abord duplique : seul l'element vise change.

Une CHARGE R-L de ``Y_DER`` (ou de ``T_Loads``) n'est pas une ligne a part :
elle est la consommation du bus. Elle est remplacee comme dans FR-051 : sa
puissance est declaree en generation negative au meme bus (meme injection
nette pour le load flow) avant d'y poser l'element a equations.

CE QUI NE SE MODIFIE PAS (et la raison s'affiche) : les noeuds (loi des
noeuds, capacite du bus), les branches (deux ports alternatifs), les
elements a port continu, les condensateurs et reactances shunt (admittances du
Y-bus). Un element ouvert par un disjoncteur ou cite par ``SIMPLIFICATIONS``
n'est pas remplace non plus : le message dit lequel.
"""

from __future__ import annotations

import copy
import math
from typing import Any

import sympy as sp

from g2egrid.elements.blocs import SchemaInvalide, compiler
from g2egrid.elements.blocs_equations import (
    _Imprimeur,
    ajouter_equation_modele,
    compiler_equations,
    equations_depuis_instance,
    equations_depuis_schema,
    est_modele_equations,
    lignes_modele,
    modifier_equation_modele,
    modifier_parametre_modele,
    supprimer_equation_modele,
)
from g2egrid.elements.bus import Bus
from g2egrid.noyau.element import Element
from g2egrid.reseau.tables import TABLES, Tables, assembler

__all__ = [
    "MODIFIABLES",
    "ajouter_equation",
    "equations_reseau",
    "modele_element",
    "modifier_equation",
    "modifier_parametre",
    "noms_element",
    "remplacer_par_modele",
    "retour_au_natif",
    "revenir_au_natif",
    "supprimer_equation",
]

#: les tables dont une ligne se remplace par un element a equations
MODIFIABLES = (
    "T_GFM1s",
    "T_GFM2s",
    "T_GFM3s",
    "T_GFM4s",
    "T_GFL1s",
    "T_GFL2s",
    "T_SMs",
    "T_IBs",
    "T_Loads",
    "Y_DER",
    "T_Blocs",
)
_RAISONS_TABLE = {
    "T_Shunts": (
        "reactance ou condensateur shunt : une admittance du Y-bus, sans "
        "equation dynamique de l'equipement a modifier"
    ),
    "T_Lines_SI": "branche a deux ports alternatifs : equations imposees par le reseau",
    "T_Lines_pu": "branche a deux ports alternatifs : equations imposees par le reseau",
    "T_Transfos": "branche a deux ports alternatifs : equations imposees par le reseau",
}


# ------------------------------------------------------------ les noms
def noms_element(t: Tables, table: str, index: int) -> list[str]:
    """Les noms que l'assembleur donne a l'element de ``table[index]`` (le
    premier present dans le reseau est le bon). Une ligne de ``Y_DER`` donne
    la charge (``Load<bus>``, si le bus en porte une) puis le noeud
    (``bus<bus>``)."""
    r = list(getattr(t, table)[index])
    if table.startswith(("T_GFM", "T_GFL")):
        return [f"{table[2:5]}{int(r[1])}"]
    if table in ("T_SMs", "T_IBs"):
        return [f"{table[2:4]}{int(r[1])}"]
    if table in ("T_Lines_SI", "T_Lines_pu"):
        return [f"L{int(r[0])}_{int(r[1])}_{int(r[2])}"]
    if table == "T_Transfos":
        return [f"T{int(r[0])}_{int(r[1])}_{int(r[2])}"]
    if table == "T_LinesDC":
        return [f"CDC{int(r[0])}_{int(r[1])}_{int(r[2])}"]
    if table == "T_Loads":
        return [f"Load{r[4]}{int(r[0])}_{int(r[1])}"]
    if table == "T_Shunts":
        return [f"Xsh{int(r[0])}_{int(r[1])}", f"Csh{int(r[0])}_{int(r[1])}"]
    if table == "T_Blocs":
        return [f"Blocs{int(r[0])}_{int(r[1])}"]
    if table == "T_DCDC":
        return [f"DCDC{int(r[0])}_{int(r[1])}"]
    if table == "T_LoadsDC":
        return [f"LoadDC{int(r[0])}_{int(r[1])}"]
    if table == "Y_DER":
        return [f"Load{int(r[0])}", f"bus{int(r[0])}"]
    if table == "T_BusDC":
        return [f"bus{int(r[0])}"]
    if table == "T_VSCs":
        return [f"GFL{int(r[1])}", f"GFM{int(r[1])}"]
    return []


def _correspondance(t: Tables) -> dict[str, tuple[str, int]]:
    """Nom d'element assemble -> ``(table, indice de ligne)``. Le premier
    gagne : les tables de reseau alternatif d'abord, le continu ensuite."""
    out: dict[str, tuple[str, int]] = {}
    for table in TABLES:
        for i in range(len(getattr(t, table))):
            for nom in noms_element(t, table, i):
                out.setdefault(nom, (table, i))
    return out


def _localiser(t: Tables, element: Any) -> tuple[str, str, int]:
    """``(nom, table, indice)`` de l'element designe par son nom
    (``"GFL3"``), par ``(table, indice)`` ou par ``{"table", "index"}``."""
    if isinstance(element, dict):
        element = (element.get("table"), element.get("index"))
    if isinstance(element, (tuple, list)) and len(element) == 2:
        table, index = element
        if table not in TABLES:
            raise SchemaInvalide(
                f"table inconnue : {table!r} (tables connues : {', '.join(TABLES)})"
            )
        n = len(getattr(t, table))
        if not isinstance(index, int) or not 0 <= index < n:
            raise SchemaInvalide(
                f"{table} n'a pas de ligne d'indice {index} : "
                + (f"l'indice va de 0 a {n - 1}" if n else "la table est vide")
            )
        noms = noms_element(t, table, index)
        if not noms:
            raise SchemaInvalide(f"{table}[{index}] : aucun element dynamique")
        nom = noms[0]
        if table == "Y_DER":
            r = t.Y_DER[index]
            if not (float(r[7]) or float(r[8])):
                nom = noms[1]  # un bus sans charge : le noeud
        return nom, table, index
    corr = _correspondance(t)
    nom = str(element)
    if nom not in corr:
        raise SchemaInvalide(f"element {nom!r} inconnu des tables")
    return (nom, *corr[nom])


# ----------------------------------------------- modeles des natifs (cache)
_CACHE: dict[Any, dict] = {}
_CACHE_MAX = 256


def _modele_natif(el: Element) -> dict:
    """Le modele a equations d'un element natif assemble (copie : le cache
    est partage, par classe et par valeur de parametres)."""
    cle = (type(el), tuple(sorted((k, float(v)) for k, v in el.pdict.items())))
    if cle not in _CACHE:
        if len(_CACHE) >= _CACHE_MAX:
            _CACHE.pop(next(iter(_CACHE)))
        _CACHE[cle] = equations_depuis_instance(el)
    return copy.deepcopy(_CACHE[cle])


def _modele_en_equations(modele: dict) -> dict:
    """Le modele de ``MODELES_BLOCS`` en equations (un schema-bloc est converti)."""
    return modele if est_modele_equations(modele) else equations_depuis_schema(modele)


# --------------------------------------------------------- liste globale
def _raison(el: Element, table: str | None) -> str | None:
    """Pourquoi l'element ne se modifie pas, ``None`` s'il se modifie."""
    if isinstance(el, Bus):
        return (
            "noeud du reseau : la loi des noeuds et la capacite du bus sont "
            "imposees par le reseau"
        )
    if el.ports_dc:
        return (
            f"{type(el).__name__} a un port continu ({', '.join(el.ports_dc)}) : "
            "le modele a equations ne prend en charge que les elements a un "
            "seul port alternatif"
        )
    if tuple(el.ports) != ("t",):
        return (
            f"{type(el).__name__} a {len(el.ports)} port(s) alternatif(s) : "
            "le modele a equations ne prend en charge que les elements a un "
            "seul port alternatif"
        )
    if table is None:
        return "element sans ligne de table : on ne sait pas le remplacer"
    if table in _RAISONS_TABLE:
        return _RAISONS_TABLE[table]
    if table not in MODIFIABLES:
        return f"element de {table} : pas de remplacement par un modele a equations"
    return None


def _lignes_spec(cls: type) -> list[dict]:
    """Les equations (etats, algebriques, courants de port) d'une classe,
    lues sur son modele compile."""
    from types import SimpleNamespace

    spec = cls.compiled().spec
    ns = SimpleNamespace()
    for noms in (spec.states, spec.algebs, spec.inputs, spec.params):
        for n in noms:
            setattr(ns, n, sp.Symbol(n))
    f, g, h = spec.equations(ns)
    out = []
    imp = _Imprimeur()
    for n, e in zip(spec.states, f):
        out.append(
            {"genre": "etat", "variable": n, "expression": imp.doprint(sp.sympify(e))}
        )
    for n, e in zip(spec.algebs, g):
        out.append(
            {
                "genre": "algebrique",
                "variable": n,
                "expression": imp.doprint(sp.sympify(e)),
                "implicite": True,
            }
        )
    for n, e in list(zip(spec.outputs, h))[: cls.n_port_signaux()]:
        out.append(
            {"genre": "sortie", "variable": n, "expression": imp.doprint(sp.sympify(e))}
        )
    return out


def equations_reseau(tables: Tables) -> list[dict]:
    """La liste globale des equations du reseau : une ligne par equation de
    chaque element assemble (equipements, branches et noeuds compris) ::

        {"element": "GFL3", "table": "T_GFL1s", "index": 0, "bus": 3,
         "classe": "GFL1", "genre": "etat" | "algebrique" | "sortie",
         "variable": "theta_pll", "expression": "wb*(w_pll - w)",
         "modifiable": True, "modifiee": False}

    ``modifiable`` dit si la ligne se modifie (:func:`modifier_equation`) ;
    sinon ``pourquoi`` donne la raison. ``modifiee`` : l'element est deja
    remplace par un modele a equations ou a blocs (on montre SON modele).
    Une algebrique qui n'est pas explicitee (``0 = residu``) porte
    ``implicite``. Les noeuds, branches et elements a port continu sont
    ecrits comme ils sont compiles (noms de signaux du modele), sans les
    entrees renommees des equipements. Chaque etat de ``build_state_space``
    apparait une fois, ``<element_reseau>.<variable>``. ``element`` est le nom
    a montrer : pour un element deja converti (``T_Blocs``), son nom natif
    d'origine ; ``element_reseau``, son nom dans le reseau assemble."""
    net = assembler(tables)
    corr = _correspondance(tables)
    t_blocs = {f"Blocs{int(r[0])}_{int(r[1])}": r for r in tables.T_Blocs}
    modeles_blocs: dict[str, dict] = {}
    cache_spec: dict[type, list[dict]] = {}
    lignes: list[dict] = []
    for el in net.all_elements:
        table, index = corr.get(el.name, (None, None))
        raison = _raison(el, table)
        modifiee = el.name in t_blocs
        modele = None
        if raison is None:
            try:
                if modifiee:
                    nom_m = t_blocs[el.name][2]
                    if nom_m not in modeles_blocs:
                        modeles_blocs[nom_m] = _modele_en_equations(
                            tables.MODELES_BLOCS[nom_m]
                        )
                    modele = modeles_blocs[nom_m]
                else:
                    modele = _modele_natif(el)
            except SchemaInvalide as err:
                raison = str(err)
        # Un element converti s'appelle Blocs<id>_<bus> dans le reseau
        # assemble ; la liste garde le nom NATIF d'origine (celui que
        # l'utilisateur connait) et donne le nom assemble a part.
        nom_affiche = el.name
        if modifiee:
            ligne_b = t_blocs[el.name]
            org = (tables.MODELES_BLOCS.get(ligne_b[2]) or {}).get("origine") or {}
            if org.get("element") and org.get("id_bloc") == int(ligne_b[0]):
                nom_affiche = str(org["element"])
        base = {
            "element": nom_affiche,
            "element_reseau": el.name,
            "table": table,
            "index": index,
            "bus": int(el.buses[0])
            if getattr(el, "buses", None)
            else getattr(el, "bus", None),
            "classe": type(el).__name__,
        }
        if modele is not None:
            for ligne in lignes_modele(modele):
                lignes.append(base | ligne | {"modifiable": True, "modifiee": modifiee})
            continue
        cls = type(el)
        if cls not in cache_spec:
            cache_spec[cls] = _lignes_spec(cls)
        for ligne in cache_spec[cls]:
            lignes.append(
                base
                | ligne
                | {"modifiable": False, "modifiee": False, "pourquoi": raison}
            )
    return lignes


# --------------------------------------------------------- edition
def _element_assemble(net, nom: str) -> Element:
    for el in net.all_elements:
        if el.name == nom:
            return el
    raise SchemaInvalide(
        f"l'element {nom!r} n'existe pas dans le reseau assemble (ligne ouverte "
        "par un disjoncteur, ou sans puissance ?)"
    )


def _modele(t: Tables, nom: str, table: str, index: int, net=None):
    """``(modele, element, reseau)`` : le modele a equations de l'element
    (le natif converti, ou le modele de ``T_Blocs`` converti si c'est un
    schema-bloc) ; ``element`` et ``reseau`` (assemble) sont ``None`` pour un
    element deja a modele."""
    if table == "T_Blocs":
        return (
            copy.deepcopy(_modele_en_equations(t.MODELES_BLOCS[t.T_Blocs[index][2]])),
            None,
            None,
        )
    net = net or assembler(t)
    el = _element_assemble(net, nom)
    raison = _raison(el, table)
    if raison:
        raise SchemaInvalide(f"{nom} ne se modifie pas : {raison}")
    return _modele_natif(el), el, net


def modele_element(tables: Tables, element: Any) -> dict:
    """Le modele a equations de l'element : le natif converti (avec ses
    parametres), ou le modele deja pose dans ``T_Blocs`` (un schema-bloc est
    converti). Leve :class:`SchemaInvalide` si l'element ne se modifie pas
    (la raison est dans le message)."""
    nom, table, index = _localiser(tables, element)
    return _modele(tables, nom, table, index)[0]


def _equivalents(a: dict, b: dict) -> bool:
    """Les deux modeles compilent en les memes equations (meme structure
    SymPy) : une modification qui ne change rien ne change rien."""
    ca, cb = compiler_equations(a), compiler_equations(b)
    return (
        ca.f == cb.f
        and ca.g == cb.g
        and ca.h == cb.h
        and ca.parametres == cb.parametres
        and ca.etats == cb.etats
        and ca.algebriques == cb.algebriques
        and ca.valeurs == cb.valeurs
        and ca.departs == cb.departs
    )


def _nom_libre(t: Tables, base: str) -> str:
    nom, k = base, 1
    while nom in t.MODELES_BLOCS:
        k += 1
        nom = f"{base}_{k}"
    return nom


def _poser_natif(
    t: Tables, nom: str, table: str, index: int, el, net, modele: dict
) -> str:
    """Remplace, DANS ``t`` (une copie privee), l'element natif ``nom`` par
    ``modele`` dans ``T_Blocs`` ; rend le nom du modele."""
    raison = _raison(el, table)
    if raison:
        raise SchemaInvalide(f"{nom} ne se modifie pas : {raison}")
    if t.SIMPLIFICATIONS and nom in t.SIMPLIFICATIONS:
        raise SchemaInvalide(
            f"{nom} est cite par SIMPLIFICATIONS : retirer cette hypothese avant "
            "de remplacer l'element par un modele a equations"
        )
    row = list(getattr(t, table)[index])
    if (table, int(row[0])) in {(a, int(b)) for a, b in t.OUVERTS}:
        raise SchemaInvalide(
            f"{nom} est ouvert par un disjoncteur : le refermer d'abord"
        )
    bus = int(row[1]) if table != "Y_DER" else int(row[0])
    retraits: list[dict] = []
    remplacements: list[dict] = []
    if table == "Y_DER":
        Sn = float(t.Sb)
        nouveau = list(row)
        nouveau[5] = float(row[5]) - float(row[7])
        nouveau[6] = float(row[6]) - float(row[8])
        nouveau[7] = nouveau[8] = 0.0
        remplacements.append(
            {"table": "Y_DER", "index": index, "ligne": list(row), "posee": nouveau}
        )
        t.Y_DER[index] = tuple(nouveau)
    elif table == "T_Loads":
        Sn = float(t.Sb)
        i_bus = next((k for k, r in enumerate(t.Y_DER) if int(r[0]) == bus), None)
        if i_bus is None:
            raise SchemaInvalide(f"{nom} : le bus {bus} n'est pas dans Y_DER")
        old = list(t.Y_DER[i_bus])
        nouveau = list(old)
        nouveau[5] = float(old[5]) - float(row[2])
        nouveau[6] = float(old[6]) - float(row[3])
        remplacements.append(
            {"table": "Y_DER", "index": i_bus, "ligne": old, "posee": nouveau}
        )
        t.Y_DER[i_bus] = tuple(nouveau)
        retraits.append({"table": table, "index": index, "ligne": list(row)})
        getattr(t, table).pop(index)
    else:
        Sn = float(row[2])
        retraits.append({"table": table, "index": index, "ligne": list(row)})
        getattr(t, table).pop(index)
    if table in ("Y_DER", "T_Loads"):
        autres = [
            e.name
            for e in net.all_elements
            if e.name != nom
            and len(e.ports) == 1
            and not e.ports_dc
            and e.buses[0] == bus
            and not isinstance(e, Bus)
            and e.admittance_ybus() is None
            and getattr(e, "scheduled_injection", lambda: None)() is None
        ]
        if autres:
            raise SchemaInvalide(
                f"{nom} : le bus {bus} porte deja {', '.join(autres)}, dont "
                "l'injection n'est pas programmee ; une charge remplacee par un "
                "modele a equations prend la puissance du load flow, il lui "
                "faut un bus a elle"
            )
    rid = max((int(r[0]) for r in t.T_Blocs), default=0) + 1
    cle = _nom_libre(t, f"{type(el).__name__}_{nom}")
    modele = copy.deepcopy(modele)
    modele["nom"] = cle
    modele["origine"] = {
        "element": nom,
        "id_bloc": rid,
        "retraits": retraits,
        "remplacements": remplacements,
    }
    if modele.get("amorce"):
        # la ligne native d'origine, a portee de la lecture de l'amorce
        premiere = (retraits or remplacements)[0]
        modele["amorce"]["natif"] = {
            "table": premiere["table"],
            "ligne": list(premiere["ligne"]),
        }
    t.MODELES_BLOCS = {**t.MODELES_BLOCS, cle: modele}
    t.T_Blocs = [*t.T_Blocs, (rid, bus, cle, Sn)]
    return cle


def _remplacer_modele(t: Tables, index: int, modele: dict) -> str:
    """Dans ``t`` (copie privee) : le modele de la ligne ``T_Blocs[index]`` est
    ``modele`` (un modele partage est d'abord duplique) ; rend sa cle."""
    row = t.T_Blocs[index]
    rid, bus, cle = int(row[0]), int(row[1]), row[2]
    partage = sum(1 for r in t.T_Blocs if r[2] == cle) > 1
    modele = copy.deepcopy(modele)
    ancien = t.MODELES_BLOCS[cle]
    origine = (ancien.get("origine") if isinstance(ancien, dict) else None) or None
    # Un SCHEMA-BLOC dont on modifie une equation devient un modele a
    # equations : le dessin (blocs, positions, sous-systemes) ne doit pas
    # disparaitre sous le meme nom. Il reste dans MODELES_BLOCS, et l'element
    # recoit un modele a part, ``<nom>_eq``.
    dessin = not est_modele_equations(ancien) and est_modele_equations(modele)
    if partage or dessin:
        nouvelle = _nom_libre(t, f"{cle}_Blocs{rid}_{bus}" if partage else f"{cle}_eq")
        if origine and origine.get("id_bloc") != rid:
            origine = None
        elif origine:
            # le proprietaire de l'origine part avec sa copie
            ancien = dict(ancien)
            ancien.pop("origine", None)
            t.MODELES_BLOCS = {**t.MODELES_BLOCS, cle: ancien}
        cle = nouvelle
        modele["nom"] = cle
        t.T_Blocs = [
            (*r[:2], cle, *r[3:]) if k == index else r for k, r in enumerate(t.T_Blocs)
        ]
    if origine:
        modele["origine"] = origine
    else:
        modele.pop("origine", None)
    t.MODELES_BLOCS = {**t.MODELES_BLOCS, cle: modele}
    return cle


def _editer(tables: Tables, element: Any, transformer) -> tuple[Tables, dict, int]:
    """Applique ``transformer(modele) -> modele`` au modele de l'element, dans
    de NOUVELLES tables ; rend ``(tables, modele, bus)``. Les tables recues ne
    changent pas, un refus non plus. Si le modele reste le meme, les tables
    sont rendues telles quelles (une copie)."""
    t = tables.copie()
    nom, table, index = _localiser(t, element)
    ligne = t.Y_DER[index] if table == "Y_DER" else getattr(t, table)[index]
    bus = int(ligne[0]) if table == "Y_DER" else int(ligne[1])
    courant, el, net = _modele(t, nom, table, index)
    nouveau = transformer(copy.deepcopy(courant))
    if _equivalents(courant, nouveau):
        return t, courant, bus
    if table == "T_Blocs":
        cle = _remplacer_modele(t, index, nouveau)
    else:
        cle = _poser_natif(t, nom, table, index, el, net, nouveau)
    return t, t.MODELES_BLOCS[cle], bus


def modifier_equation(
    tables: Tables,
    element: Any,
    variable: str,
    expression: str,
    genre: str | None = None,
) -> Tables:
    """De nouvelles tables ou l'equation de ``variable`` de ``element`` est
    ``expression`` : la derivee d'un etat (``d<variable>/dt``), l'expression
    d'une variable algebrique, une sortie (``id``, ``iq``...). ``genre``
    leve l'ambiguite d'un nom a la fois etat et sortie (``theta``).

    L'expression est verifiee AVANT toute ecriture (analyseur sur, noms
    connus, compilation) : une erreur leve :class:`SchemaInvalide` et les
    tables recues ne bougent pas. Une equation inchangee ne change rien (pas de
    conversion inutile). Voir le module pour le remplacement d'un natif."""
    t, _m, _bus = _editer(
        tables,
        element,
        lambda m: modifier_equation_modele(m, variable, expression, genre),
    )
    return t


def ajouter_equation(
    tables: Tables,
    element: Any,
    genre: str,
    variable: str,
    expression: str,
    depart: float | None = None,
) -> Tables:
    """De nouvelles tables ou ``element`` a une equation de plus : un ETAT
    (``genre="etat"``, ``d<variable>/dt = expression``), une variable
    ALGEBRIQUE (``"algebrique"``, ``variable = expression``) ou implicite
    (``"implicite"``, ``0 = expression``). ``depart`` : la valeur d'ou part le
    Newton de l'initialisation pour cette variable."""
    t, _m, _bus = _editer(
        tables,
        element,
        lambda m: ajouter_equation_modele(m, genre, variable, expression, depart),
    )
    return t


def supprimer_equation(tables: Tables, element: Any, variable: str) -> Tables:
    """De nouvelles tables sans l'etat ou la variable algebrique ``variable``
    de ``element`` ; elle ne doit plus etre employee par une autre equation
    (le message dit par qui)."""
    t, _m, _bus = _editer(
        tables, element, lambda m: supprimer_equation_modele(m, variable)
    )
    return t


def _schema_bloc_de(tables: Tables, element: Any):
    """``(copie des tables, indice, schema)`` si l'element est une ligne de
    ``T_Blocs`` decrite par un schema-bloc ; ``None`` sinon."""
    t = tables.copie()
    _nom, table, index = _localiser(t, element)
    if table != "T_Blocs":
        return None
    modele = t.MODELES_BLOCS[t.T_Blocs[index][2]]
    return None if est_modele_equations(modele) else (t, index, modele)


def modifier_parametre(tables: Tables, element: Any, nom: str, valeur: float) -> Tables:
    """De nouvelles tables ou le parametre ``nom`` du modele de ``element``
    vaut ``valeur``. ``PARAMS_PAR_BUS`` surcharge les parametres d'un element
    a equations ; un reglage de ce bus pour ``nom`` est donc retire (sinon il
    l'emporterait sur la valeur donnee)."""
    if not math.isfinite(float(valeur)):
        raise SchemaInvalide(f"parametre {nom!r} : la valeur doit etre finie")
    dessin = _schema_bloc_de(tables, element)
    if dessin is not None:
        # un SCHEMA-BLOC reste un schema-bloc : seul son parametre change (le
        # passer en equations pour une valeur ferait perdre le dessin)
        t, index, ancien = dessin
        if nom not in (ancien.get("parametres") or {}):
            raise SchemaInvalide(
                f"le modele {ancien.get('nom')!r} n'a pas de parametre {nom!r}"
            )
        nouveau = copy.deepcopy(ancien)
        nouveau["parametres"][nom] = float(valeur)
        compiler(nouveau)  # un refus (diviseur nul...) avant d'ecrire
        _remplacer_modele(t, index, nouveau)
        bus = int(t.T_Blocs[index][1])
    else:
        t, _m, bus = _editer(
            tables, element, lambda m: modifier_parametre_modele(m, nom, valeur)
        )
    par_bus = t.PARAMS_PAR_BUS.get(bus)
    if par_bus and nom in par_bus:
        t.PARAMS_PAR_BUS = {
            **t.PARAMS_PAR_BUS,
            bus: {k: v for k, v in par_bus.items() if k != nom},
        }
    return t


def _rendre_ligne(t: Tables, r: dict, nom: str) -> int:
    """Defait un remplacement de ligne (``Y_DER`` d'une charge). La ligne est
    retrouvee par son BUS, pas par son indice (des bus ont pu etre ajoutes ou
    retires depuis). Si elle est restee telle que le remplacement l'a posee,
    la ligne d'origine revient EXACTEMENT ; si l'utilisateur l'a retouchee
    entre-temps (consigne de tension, type du bus...), seules les puissances
    deplacees par le remplacement reviennent, et ses retouches restent."""
    liste = getattr(t, r["table"])
    gardee = list(r["ligne"])
    k = next((i for i, x in enumerate(liste) if int(x[0]) == int(gardee[0])), None)
    if k is None:
        raise SchemaInvalide(
            f"{nom} : le bus {int(gardee[0])} n'est plus dans {r['table']} ; "
            "le retour au modele natif n'est plus possible (Ctrl+Z annule les "
            "dernieres editions)"
        )
    posee = r.get("posee")
    actuelle = list(liste[k])
    if posee is None or actuelle == list(posee):
        liste[k] = tuple(gardee)
        return k
    for c in (5, 6, 7, 8):  # Pg, Qg, Pl, Ql
        actuelle[c] = float(actuelle[c]) + float(gardee[c]) - float(posee[c])
    liste[k] = tuple(actuelle)
    return k


def _reinserer(t: Tables, r: dict) -> int:
    """Rend une ligne retiree a sa table, a sa place si elle existe encore ;
    si un autre element a pris son identifiant depuis, elle en recoit un neuf."""
    liste = getattr(t, r["table"])
    ligne = list(r["ligne"])
    pris = {int(x[0]) for x in liste}
    if int(ligne[0]) in pris:
        ligne[0] = max(pris) + 1
    k = min(int(r["index"]), len(liste))
    liste.insert(k, tuple(ligne))
    return k


def remplacer_par_modele(
    tables: Tables, element: Any, genre: str = "blocs", detail: str = "blocs"
) -> tuple[Tables, str]:
    """De nouvelles tables ou l'element NATIF est remplace par son modele :
    un schema-bloc (``genre="blocs"``, ``detail`` : ``"blocs"`` ou
    ``"equations"``) ou un modele a equations (``genre="equations"``), ecrit
    avec SES parametres. Rend ``(tables, nom du modele)``.

    C'est le meme remplacement que celui d'une equation modifiee (voir le
    module) : la puissance d'une charge passe en generation negative a son
    bus, un element ouvert ou cite par ``SIMPLIFICATIONS`` est refuse, et le
    modele garde son ``origine`` -- :func:`revenir_au_natif` rend l'element."""
    from g2egrid.elements.blocs_inverse import depuis_instance

    if genre not in ("blocs", "equations"):
        raise SchemaInvalide(f"genre {genre!r} : 'blocs' ou 'equations'")
    t = tables.copie()
    nom, table, index = _localiser(t, element)
    if table == "T_Blocs":
        raise SchemaInvalide(
            f"{nom} est deja decrit par un modele (T_Blocs) : l'ouvrir dans "
            "l'onglet Blocs"
        )
    net = assembler(t)
    el = _element_assemble(net, nom)
    raison = _raison(el, table)
    if raison:
        raise SchemaInvalide(f"{nom} ne se remplace pas : {raison}")
    modele = _modele_natif(el) if genre == "equations" else depuis_instance(el, detail)
    cle = _poser_natif(t, nom, table, index, el, net, modele)
    return t, cle


def revenir_au_natif(tables: Tables, element: Any) -> Tables:
    """De nouvelles tables ou l'element, remplace par un modele a equations
    (ou un schema-bloc issu d'un natif), retrouve sa ligne d'origine : les
    tables sont celles d'avant le remplacement (le spectre aussi, exactement).
    Leve :class:`SchemaInvalide` si l'element n'a pas d'origine native
    (modele ecrit a la main, ou partage)."""
    return retour_au_natif(tables, element)[0]


def retour_au_natif(tables: Tables, element: Any) -> tuple[Tables, tuple[str, int]]:
    """:func:`revenir_au_natif`, avec le LIEU de l'element rendu :
    ``(tables, (table, indice))`` -- la ligne reinseree (``T_GFL2s``,
    ``T_Loads``...), ou celle de ``Y_DER`` pour la charge d'un bus."""
    t = tables.copie()
    nom, table, index = _localiser(t, element)
    if table != "T_Blocs":
        raise SchemaInvalide(f"{nom} est deja natif")
    row = t.T_Blocs[index]
    cle = row[2]
    modele = t.MODELES_BLOCS[cle]
    origine = modele.get("origine") if isinstance(modele, dict) else None
    if not origine or origine.get("id_bloc") != int(row[0]):
        raise SchemaInvalide(
            f"{nom} n'a pas d'origine native connue : il ne vient pas du "
            "remplacement d'un element de la bibliotheque"
        )
    lieu = None
    for r in origine.get("remplacements") or []:
        k = _rendre_ligne(t, r, nom)
        lieu = lieu or (r["table"], k)
    retraits = sorted(origine.get("retraits") or [], key=lambda r: r["index"])
    for r in retraits:
        k = _reinserer(t, r)
        if r is retraits[0]:
            lieu = (r["table"], k)  # la ligne de l'element prime sur celle du bus
    if lieu is None:
        raise SchemaInvalide(
            f"{nom} : l'origine native du modele est vide, rien a rendre"
        )
    t.T_Blocs = [r for k, r in enumerate(t.T_Blocs) if k != index]
    if not any(r[2] == cle for r in t.T_Blocs):
        t.MODELES_BLOCS = {k: v for k, v in t.MODELES_BLOCS.items() if k != cle}
    return t, lieu
