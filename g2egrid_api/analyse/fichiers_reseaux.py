"""Les fichiers de reseaux : liste, lecture, enregistrement, import d'autres outils.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import json
import pathlib
import re
import g2egrid.reseau.catalogue
from g2egrid.reseau.tables import TABLES
from g2egrid.reseau.zones import (
    FORMAT,
    VERSION,
    Reseau,
    aplatir,
    reseau_depuis_json,
    reseau_vers_json,
)
from .catalogue_api import _reponse_reseau


# ----------------------------------------------------------------- fichiers
def fichiers() -> list[dict]:
    """Les reseaux au format natif -- exemples et sauvegardes --, SANS les assembler.

    On lit l'en-tete de chaque fichier et rien de plus : une liste doit
    s'afficher tout de suite, et aplatir wscc3sm_zones pour compter ses bus
    couterait une seconde par ligne. ``origine`` dit d'ou vient chacun :
    un exemple livre avec l'outil, ou une sauvegarde de l'utilisateur.

    UN FICHIER ILLISIBLE EST LISTE QUAND MEME, avec son erreur. C'est
    exactement celui qu'on cherche quand on ouvre cette liste, et le taire le
    rendrait invisible sans le rendre valide.
    """
    suffixe = g2egrid.reseau.catalogue.SUFFIXE
    sources = (
        ("exemple", g2egrid.reseau.catalogue.DOSSIER_RESEAUX),
        ("sauvegarde", g2egrid.reseau.catalogue.dossier_sauvegardes()),
    )
    out = []
    for origine, dossier in sources:
        if not dossier.is_dir():
            continue
        for p in sorted(dossier.glob(f"*{suffixe}")):
            out.append(_fiche_fichier(p, suffixe, origine))
    return out


def _fiche_fichier(p: pathlib.Path, suffixe: str, origine: str) -> dict:
    etat = p.stat()
    fic = {
        "nom": p.name[: -len(suffixe)],
        "origine": origine,
        "octets": etat.st_size,
        "modifie": etat.st_mtime,
    }
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
        fic.update(
            titre=str(d.get("nom", fic["nom"])),
            f_zoom=float(d.get("f_zoom", 5.0)),
            zones=list(d.get("zones") or {}),
            instances=len(d.get("instances") or []),
            bus_transport=len((d.get("transport") or {}).get("Y_DER") or []),
        )
    except Exception as err:  # noqa: BLE001 -- fichier tronque, JSON casse
        fic["erreur"] = f"{type(err).__name__}: {err}"
    return fic


#: Ce qu'on accepte comme nom de reseau enregistre. Volontairement etroit :
#: voir :func:`_nom_de_fichier`.
_NOM_FICHIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}")


def _nom_de_fichier(nom: str) -> str:
    """Un nom de reseau dont on puisse faire un nom de fichier.

    LE NOM DEVIENT UN CHEMIN (``g2egrid/reseaux/<nom>.g2egrid.json``). Sans ce
    filtre, un nom portant ``..`` ou une barre ecrirait ou il veut sur le
    disque de qui sert l'interface, et un nom vide ecraserait ``.g2egrid.json``.
    On refuse plutot que d'assainir en silence : un reseau enregistre sous un
    autre nom que celui qu'on a tape est un reseau qu'on ne retrouve pas.
    """
    nom = str(nom or "").strip()
    if not _NOM_FICHIER.fullmatch(nom):
        raise ValueError(
            f"nom de reseau invalide : {nom!r}. Lettres, chiffres, '_' et '-' "
            "seulement, 64 caracteres au plus, en commencant par une lettre ou "
            "un chiffre -- ce nom devient un nom de fichier."
        )
    return nom


_REGLAGES_JSON = (
    "SM_PARAMS",
    "GFM_PARAMS",
    "GFL_PARAMS",
    "PARAMS_PAR_BUS",
    "OUVERTS",
    # le traitement des noeuds sans capacite (FR-005), quand ce n'est pas le defaut
    "NOEUDS",
    "R_NOEUDS",
)
_CHAMPS_INSTANCE = (
    "zone",
    "prefixe",
    "bus_amont",
    "K",
    "transfo",
    "charge_residuelle",
    "shunt_pcc",
)


def _vers_fichier(d: dict) -> dict:
    """La forme que rend :func:`depuis_catalogue`, remise a celle du fichier.

    L'interface ouvre un reseau, l'edite, et le renvoie -- pour l'aplatir, ou
    pour l'enregistrer. Encore faut-il que ce qu'elle renvoie soit ce que ces
    routes attendent. Or la reponse de lecture et le fichier n'ont pas la meme
    forme : le fichier porte une enveloppe (format, version), met f0, Sb et
    les reglages en tete, range le transport en tables nues, et reserve son
    ``positions`` de tete au TRANSPORT ; la lecture, elle, donne ``tables`` et
    ``positions`` du reseau APLATI -- c'est ce que lisent le schema et toutes
    les routes de calcul -- et le transport sous sa propre cle. Renvoyer la
    lecture telle quelle a ``aplatir`` echouait donc sur "format 'None'
    inconnu" : on pouvait ouvrir un reseau, pas le sauver.

    On convertit ici plutot que d'ecrire un second lecteur : il n'y a qu'un
    analyseur du format, ``g2egrid.reseau.zones.reseau_depuis_json``, et le dedoubler
    serait exactement la derive que ce projet a deja payee ailleurs.

    Trois formes sont acceptees : le fichier lui-meme (rendu tel quel), la
    reponse de lecture (transport, zones, instances), et un reseau dessine a
    la main, qui n'a que ``tables`` et ``positions`` -- son transport est
    alors le reseau entier.
    """
    # Un fichier natif se reconnait a son ``transport``, pas a son seul champ
    # ``format`` : l'ancien bouton Enregistrer du constructeur ecrivait
    # {"format": "g2egrid", "version": 1, "tables": ..., "positions": ...} --
    # il SE DISAIT natif en rangeant ses tables a plat. Pris pour natif, ce
    # fichier donnait un transport vide : 0 bus, et aucune erreur.
    if d.get("format") and "transport" in d:
        return d
    lecture = d.get("transport")
    tr = (lecture or {}).get("tables") or d.get("tables") or {}
    positions = (lecture or {}).get("positions") if lecture else d.get("positions")

    def reglages(tj: dict) -> dict:
        # les tables de lecture portent leurs reglages A COTE des tables, la
        # ou le fichier les range sous une cle a part
        return {k: tj[k] for k in _REGLAGES_JSON if tj.get(k)}

    return {
        "format": FORMAT,
        "version": VERSION,
        "nom": d.get("nom") or tr.get("nom") or "reseau",
        "f_zoom": d.get("f_zoom", 5.0),
        "f0": tr.get("f0", 50.0),
        "Sb": tr.get("Sb", 100e6),
        "reglages": reglages(tr),
        "transport": {k: tr.get(k) or [] for k in TABLES},
        "zones": {
            nom: {
                "tables": {k: (z.get("tables") or {}).get(k) or [] for k in TABLES},
                "bus_pcc": z.get("bus_pcc", 0),
                "reglages": reglages(z.get("tables") or {}),
                "positions": z.get("positions") or {},
            }
            for nom, z in (d.get("zones") or {}).items()
        },
        # ``bus`` est un champ DERIVE de la lecture -- les numeros qu'occupe
        # l'instance a plat -- et n'appartient pas au fichier
        "instances": [
            {k: i[k] for k in _CHAMPS_INSTANCE if k in i}
            for i in (d.get("instances") or [])
        ],
        "positions": positions or {},
    }


def _reseau_depuis(d: dict) -> Reseau:
    """Un :class:`~g2egrid.reseau.zones.Reseau` depuis le JSON du fichier.

    ``reseau_depuis_json`` suppose son entree bien formee : une instance sans
    ``zone``, un transfo de poste a quatre valeurs au lieu de cinq lui font
    lever un ``KeyError`` ou un ``TypeError`` nus -- le premier donnerait une
    422 dont le message serait le seul mot ``'zone'``, le second une 500. On
    les rattrape ici pour que le corps de requete mal forme soit une reponse a
    rendre, comme le reste.
    """
    try:
        return reseau_depuis_json(_vers_fichier(d))
    except ValueError:
        raise  # deja redige pour un lecteur (format, version)
    except (KeyError, TypeError, IndexError, AttributeError) as err:
        raise ValueError(
            f"description de reseau illisible ({type(err).__name__} : {err}). "
            "Attendu le format d'un fichier .g2egrid.json : 'format', 'version', "
            "'transport', 'zones', 'instances'."
        ) from err


def enregistrer_fichier(reseau: dict, nom: str | None = None) -> dict:
    """Enregistre un reseau STRUCTURE au format natif, et rend sa fiche.

    LE RESEAU EST APLATI AVANT D'ETRE ECRIT, et on jette le resultat. Ce n'est
    pas une precaution de style : ``aplatir`` porte tous les garde-fous de
    structure -- prefixes qui se chevauchent, bus amont absent, slack en zone,
    ``Un`` inegaux a la fusion, reglage de famille qui deborde -- et ils sont
    rediges pour un utilisateur. Les faire jouer AVANT l'ecriture evite de
    laisser sur le disque un fichier que rien ne pourra plus relire, et rend
    l'erreur au moment ou l'on peut encore la corriger.
    """
    r = _reseau_depuis(reseau)
    nom = _nom_de_fichier(nom or r.nom)
    aplatir(r)
    chemin = g2egrid.reseau.catalogue.enregistrer(r, nom)
    return {
        "nom": nom,
        "chemin": str(chemin),
        "octets": chemin.stat().st_size,
        "zones": list(r.zones),
        "instances": len(r.instances),
    }


def lire_fichier(reseau: dict) -> dict:
    """N'importe quelle forme de reseau, rendue sous la forme de LECTURE.

    Le navigateur recoit un fichier -- ouvert depuis le disque de
    l'utilisateur -- dont il ne doit pas avoir a connaitre le format : un
    ``.g2egrid.json`` natif, la reponse d'une lecture enregistree telle quelle,
    ou l'ancien JSON a plat du constructeur. Tous passent par l'unique
    analyseur du format, et reviennent sous la forme que rend
    :func:`depuis_catalogue` -- celle que le constructeur sait editer.
    """
    return _reponse_reseau(_reseau_depuis(reseau))


def convertir_fichier(reseau: dict) -> dict:
    """Le contenu d'un vrai ``.g2egrid.json``, a telecharger.

    Le constructeur tient le reseau sous la forme de LECTURE ; telecharger
    celle-la donnerait un fichier que le catalogue ne sait pas relire, et
    qu'on ne pourrait donc pas deposer dans ``g2egrid/reseaux/``. On passe par
    le meme chemin qu'un enregistrement -- aplatissement compris, pour que
    ses garde-fous jouent AVANT qu'un fichier illisible parte sur un disque.
    """
    r = _reseau_depuis(reseau)
    aplatir(r)
    return reseau_vers_json(r)
