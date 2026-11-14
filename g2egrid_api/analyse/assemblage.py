"""L'assemblage d'un reseau edite et son cache.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import numpy as np
from g2egrid.reseau.tables import TABLES, Tables, assembler, tables_depuis_json
from ._commun import Assemble, _CACHE, _CACHE_MAX, _empreinte


# ------------------------------------------------------------- assemblage
#: les colonnes qui portent un numero de bus, et le reseau (alternatif ou
#: continu) auquel il doit appartenir, par table. Une table absente d'ici
#: n'en porte pas (``T_VSCs`` a deux : voir :func:`references_inconnues`).
_BUS_ALTERNATIF = {
    "T_Lines_SI": ("From", "To"),
    "T_Lines_pu": ("From", "To"),
    "T_Transfos": ("From", "To"),
    "T_Shunts": ("BusID",),
    "T_IBs": ("BusID",),
    "T_SMs": ("BusID",),
    "T_GFM1s": ("BusID",),
    "T_GFM2s": ("BusID",),
    "T_GFM3s": ("BusID",),
    "T_GFM4s": ("BusID",),
    "T_GFL1s": ("BusID",),
    "T_GFL2s": ("BusID",),
    "T_Loads": ("BusID",),
    "T_Blocs": ("BusID",),
    "T_VSCs": ("BusAC",),
}
_BUS_CONTINU = {
    "T_LinesDC": ("From", "To"),
    "T_ShuntsDC": ("BusID",),
    "T_DCDC": ("BusDC",),
    "T_LoadsDC": ("BusDC",),
    "T_VSCs": ("BusDC",),
}


def references_inconnues(t: Tables) -> list[str]:
    """Les lignes de tables qui nomment un bus qui n'existe pas.

    Un numero de bus faux (une ligne vers le bus 99 quand le reseau s'arrete
    a 9) levait un ``KeyError`` que l'API rendait sous la forme ``"99"``. On
    le dit avant d'assembler : quelle table, quelle ligne, quel bus, et ou le
    declarer. Rend une liste vide si tout est coherent.
    """
    ac = {int(r[0]) for r in t.Y_DER}
    dc = {int(r[0]) for r in t.T_BusDC}
    out: list[str] = []
    for tables, connus, ou in (
        (_BUS_ALTERNATIF, ac, "Y_DER"),
        (_BUS_CONTINU, dc, "T_BusDC"),
    ):
        for nom, colonnes in tables.items():
            noms = TABLES[nom]
            for i, r in enumerate(getattr(t, nom)):
                for col in colonnes:
                    try:
                        b = int(r[noms.index(col)])
                    except (ValueError, TypeError, IndexError):
                        continue  # la forme de la ligne est controlee ailleurs
                    if b not in connus:
                        out.append(
                            f"{nom}, ligne n. {i + 1} (colonne {col}) : le bus {b} "
                            f"n'existe pas dans {ou} ; corrigez le numero ou "
                            f"ajoutez ce bus a {ou}"
                        )
    return out


def controler_entree(t: Tables) -> None:
    """Refuse, en francais, un reseau que l'assemblage ne saurait pas lire :
    vide, ou qui nomme des bus inexistants."""
    if not t.Y_DER:
        raise ValueError(
            "reseau vide : la table Y_DER n'a aucun bus. Ajoutez au moins un bus "
            "et une source (un bus slack) avant de calculer."
        )
    pbs = references_inconnues(t)
    if pbs:
        plus = f" (et {len(pbs) - 3} autre(s))" if len(pbs) > 3 else ""
        raise ValueError(" ; ".join(pbs[:3]) + plus)


def construire(t: Tables) -> Assemble:
    """Assemble et resout, en gardant le resultat en cache."""
    cle = _empreinte(t)
    if cle in _CACHE:
        return _CACHE[cle]
    controler_entree(t)
    net = assembler(t)
    a = Assemble(tables=t, net=net)
    a.lf = net.run_loadflow()
    if len(_CACHE) >= _CACHE_MAX:
        _CACHE.pop(next(iter(_CACHE)))
    _CACHE[cle] = a
    return a


#: colonnes qui portent un numero de bus ou d'element : des entiers
_COLONNES_ENTIERES = {"ID", "BusID", "From", "To", "BusAC", "BusDC"}
#: suffixes des colonnes qui portent une grandeur physique : des nombres finis
_SUFFIXES_NOMBRES = (
    "_pu", "_V", "_VA", "_W", "_var", "_km", "_rad", "_deg", "_S", "_mH_km", "_uF_km",
)  # fmt: skip


def lire_tables(d: dict) -> Tables:
    """Les tables du corps d'une requete, avec leur FORME verifiee.

    ``tables_depuis_json`` est permissif (un reseau en cours de dessin est
    incomplet). Mais une ligne trop courte, un numero de bus qui n'est pas un
    entier ou une reactance qui est du texte ne sont pas un reseau incomplet,
    c'est une erreur de saisie : on dit laquelle (table, ligne, colonne), au
    lieu de ``invalid literal for int()`` ou d'un ``IndexError``.
    """
    if not isinstance(d, dict):
        raise ValueError("les tables doivent etre un objet {nom_de_table: lignes}")
    for nom in TABLES:
        v = d.get(nom)
        if v is not None and not isinstance(v, (list, tuple)):
            raise ValueError(f"{nom} doit etre une liste de lignes")
        for i, r in enumerate(v or []):
            if not isinstance(r, (list, tuple)):
                raise ValueError(
                    f"{nom}, ligne n. {i + 1} doit etre une liste de valeurs"
                )
    t = tables_depuis_json(d)
    for nom, colonnes in TABLES.items():
        for i, r in enumerate(getattr(t, nom)):
            if len(r) < len(colonnes):
                raise ValueError(
                    f"{nom}, ligne n. {i + 1} : {len(r)} valeur(s) au lieu de "
                    f"{len(colonnes)} ({', '.join(colonnes)})"
                )
            for col, v in zip(colonnes, r):
                if col in _COLONNES_ENTIERES:
                    ok = _est_entier(v)
                    attendu = "un entier (numero de bus ou d'element)"
                elif col.endswith(_SUFFIXES_NOMBRES):
                    ok = _est_nombre(v)
                    attendu = "un nombre fini"
                else:
                    continue
                if not ok:
                    raise ValueError(
                        f"{nom}, ligne n. {i + 1}, colonne {col} : {v!r} n'est pas "
                        f"{attendu}"
                    )
    return t


def _est_nombre(v) -> bool:
    if isinstance(v, bool):
        return False
    try:
        return bool(np.isfinite(float(v)))
    except (TypeError, ValueError):
        return False


def _est_entier(v) -> bool:
    return _est_nombre(v) and float(v).is_integer()


def vider_cache() -> int:
    n = len(_CACHE)
    _CACHE.clear()
    return n
