"""Les reseaux de test, charges par programme.

DEUX SOURCES, UN SEUL CONTRAT. Un reseau se decrit soit par un notebook de
``reseaux_test/``, dont la cellule de donnees ne contient que des tables, soit
par un fichier ``g2egrid/reseaux/<nom>.g2egrid.json``, qui sait en plus decrire
des ZONES instanciees plusieurs fois. :func:`charger` rend dans les deux cas
des tables APLATIES : c'est ce dont dependent l'assembleur, les etudes et
l'interface, et ce contrat ne bouge pas. Qui veut la structure en zones
demande :func:`charger_reseau`.

Les fichiers ``.g2egrid.json`` passent AVANT les notebooks. Un reseau qu'on a
converti est donc lu par sa forme native, sans executer de code ; le notebook
du meme nom reste la comme document, et le jour ou les deux divergent, c'est le
fichier qui fait foi.
"""

from __future__ import annotations

import json
import math
import os
import pathlib
import re
from dataclasses import dataclass

import numpy as np

from g2egrid.reseau.network import Network
from g2egrid.reseau.tables import TABLES, Tables, assembler
from g2egrid.reseau.zones import Reseau, aplatir, reseau_depuis_json, reseau_vers_json

__all__ = [
    "RESEAUX",
    "Fiche",
    "charger",
    "charger_reseau",
    "construire",
    "enregistrer",
    "fiche",
    "lister",
]

#: La racine des donnees : le dossier ``g2egrid/`` du depot. Une installation
#: non editable (``pipx install``) n'emporte pas les exemples, qui ne font pas
#: partie du paquet : ``G2EGRID_DONNEES`` designe alors le ``g2egrid/`` d'un
#: depot clone.
RACINE = pathlib.Path(
    os.environ.get("G2EGRID_DONNEES") or pathlib.Path(__file__).resolve().parents[2]
)

DOSSIER = RACINE / "reseaux_test"

#: Les reseaux au format natif : exemples livres et sauvegardes de l'interface.
DOSSIER_RESEAUX = RACINE / "reseaux"

SUFFIXE = ".g2egrid.json"

#: Les tables attendues dans la cellule de donnees d'un notebook. Elles sont
#: LUES dans le contrat de colonnes, et non recopiees : une table ajoutee a
#: ``g2egrid.reseau.tables.TABLES`` (les quatre tables continues du lot 4, par exemple)
#: devient lisible par les notebooks sans qu'on ait a y penser ici. Une
#: cellule qui ne la declare pas prend une table vide, comme avant.
_TABLES = tuple(TABLES)
_REGLAGES = ("SM_PARAMS", "GFM_PARAMS", "GFL_PARAMS")

#: L'ordre de la documentation : du plus simple au plus gros. Le balayage de
#: ``g2egrid/reseaux/`` s'y ajoute a la fin, sans doublon : les noms deja cites
#: gardent leur rang, et un fichier natif neuf arrive apres.
_RESEAUX_NOTEBOOKS = (
    "smib",
    "smib_sm",
    "smib_gfm",
    "smib_gfl",
    "three_machine_ib",
    "wscc9",
    "wscc9_gfm",
    "kundur2area",
    "ieee14",
    "ieee33",
    "ieee39",
    "cigre_mv",
    "wscc9_cigre",
    "wscc3sm_cigre3",
    # liaisons a courant continu : deux zones alternatives reliees par un
    # VSC-HVDC, puis trois terminaux. Ce ne sont pas des jeux de parametres
    # generiques mais des POINTS DE FONCTIONNEMENT choisis -- voir les
    # notebooks : lignes alternatives de 150 km et cables continus a
    # 0,05 ohm/km, en deca desquels le reseau devient instable pour des
    # raisons physiques (resonance LC desamortie par des terminaux a
    # puissance constante).
    "hvdc2t",
    "mtdc3t",
)


def _natifs() -> tuple[str, ...]:
    """Les noms des fichiers ``.g2egrid.json`` presents, tries."""
    if not DOSSIER_RESEAUX.is_dir():
        return ()
    return tuple(
        sorted(p.name[: -len(SUFFIXE)] for p in DOSSIER_RESEAUX.glob(f"*{SUFFIXE}"))
    )


RESEAUX: tuple[str, ...] = tuple(dict.fromkeys((*_RESEAUX_NOTEBOOKS, *_natifs())))


# LES SAUVEGARDES NE VONT PAS AVEC LES EXEMPLES. Le bouton "Au catalogue" de
# l'interface ecrivait dans g2egrid/reseaux/, a cote des exemples natifs, et
# cela cassait deux choses : un fichier non suivi apparaissait dans un dossier
# que git suit, et -- plus grave -- la suite de tests tombait, parce que la
# non-regression bit a bit enumere les reseaux alternatifs du catalogue et
# exige une reference pour chacun. Enregistrer un reseau depuis l'interface
# faisait donc echouer les tests. Les sauvegardes ont leur sous-dossier,
# ignore par git, lu a chaque appel et jamais fige dans RESEAUX.
def dossier_sauvegardes() -> pathlib.Path:
    """Calcule a l'appel, pour suivre un DOSSIER_RESEAUX redirige (tests)."""
    return DOSSIER_RESEAUX / "sauvegardes"


def sauvegardes() -> tuple[str, ...]:
    """Les noms des reseaux enregistres par l'utilisateur, lus maintenant."""
    d = dossier_sauvegardes()
    if not d.is_dir():
        return ()
    return tuple(sorted(p.name[: -len(SUFFIXE)] for p in d.glob(f"*{SUFFIXE}")))


@dataclass
class Fiche:
    """Ce qu'on peut dire d'un reseau sans l'assembler."""

    nom: str
    titre: str
    bus: int
    lignes: int
    transfos: int
    sm: int
    gfm: int
    gfl: int
    charges: int
    f0: float
    Sb: float
    f_zoom: float


# ------------------------------------------------------------- les deux formes
#: Un nom de reseau est un nom de fichier nu : ni separateur ni ``..``. Le nom
#: arrive de l'atelier (``/api/reseau/importer``) et finit dans un chemin, puis,
#: pour un notebook, dans un ``exec`` : sans ce filtre, ``../../x`` faisait
#: executer n'importe quel ``.ipynb`` du disque.
_NOM_RESEAU = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}")


def _verifier_nom(nom: str) -> str:
    if not isinstance(nom, str) or not _NOM_RESEAU.fullmatch(nom):
        raise ValueError(f"nom de reseau invalide : {nom!r}")
    return nom


def chemin_natif(nom: str) -> pathlib.Path | None:
    """Le fichier ``.g2egrid.json`` du reseau ``nom``, s'il existe.

    Les exemples d'abord, les sauvegardes ensuite : une sauvegarde qui
    porterait le nom d'un exemple ne doit pas le masquer.
    """
    _verifier_nom(nom)
    for dossier in (DOSSIER_RESEAUX, dossier_sauvegardes()):
        p = dossier / f"{nom}{SUFFIXE}"
        if p.exists():
            return p
    return None


def _cellule_donnees(nom: str) -> dict:
    """Execute la seule cellule de donnees du notebook et rend son espace de noms."""
    chemin = DOSSIER / f"{_verifier_nom(nom)}.ipynb"
    if not chemin.exists():
        if not DOSSIER_RESEAUX.is_dir() and not DOSSIER.is_dir():
            # installation sans les exemples : le dire, plutot que de laisser
            # croire a une faute de frappe dans le nom
            raise FileNotFoundError(
                f"reseau inconnu : {nom} -- le catalogue est introuvable sous "
                f"{RACINE} (les exemples ne font pas partie du paquet installe) : "
                "donner dans la variable G2EGRID_DONNEES le dossier d'un depot "
                "clone qui contient reseaux/ et reseaux_test/"
            )
        raise FileNotFoundError(f"reseau inconnu : {nom}")
    nb = json.loads(chemin.read_text(encoding="utf-8"))
    code = [c for c in nb["cells"] if c["cell_type"] == "code"]
    for cellule in code[:3]:  # imports, donnees, assemblage
        src = "".join(cellule["source"])
        if "Y_DER = [" not in src and "Y_DER=[" not in src:
            continue
        # quelques cellules de donnees calculent leurs valeurs (np.pi, sqrt...)
        g: dict = {"np": np, "math": math}
        exec(compile(src, f"{nom}:donnees", "exec"), g)  # noqa: S102
        if "Y_DER" in g:
            return g
    raise ValueError(f"{nom} : aucune cellule de donnees exploitable")


def charger_reseau(nom: str) -> Reseau:
    """Le reseau ``nom`` AVEC sa structure : transport, zones, instances.

    Un notebook ne decrit pas de zone -- il est deja a plat. On lui rend donc
    un :class:`~g2egrid.reseau.zones.Reseau` a une seule piece, dont le transport porte
    toutes les tables : le code qui manipule des reseaux structures marche
    ainsi sur les deux sources, sans avoir a demander laquelle il tient.
    """
    p = chemin_natif(nom)
    if p is not None:
        return reseau_depuis_json(json.loads(p.read_text(encoding="utf-8")))
    t = charger(nom)
    return Reseau(transport=t, nom=t.nom, f_zoom=_zoom_notebook(nom))


def charger(nom: str) -> Tables:
    """Les tables APLATIES du reseau ``nom``.

    Le contrat de cette fonction ne bouge pas : elle rend des tables pretes
    pour ``tables.assembler``, que le reseau soit ecrit a plat dans un notebook
    ou decoupe en zones dans un fichier natif.
    """
    p = chemin_natif(nom)
    if p is not None:
        return aplatir(charger_reseau(nom))
    g = _cellule_donnees(nom)
    return Tables(
        **{k: list(g.get(k, [])) for k in _TABLES},
        **{k: dict(g.get(k, {})) for k in _REGLAGES},
        PARAMS_PAR_BUS=dict(g.get("PARAMS_PAR_BUS", {})),
        f0=float(g.get("f0", 50.0)),
        Sb=float(g.get("Sb", 100e6)),
        nom=str(g.get("NOM_RESEAU", nom)),
        # la cellule de donnees peut fixer le traitement des noeuds sans
        # capacite : les reseaux de parite MATLAB gardent la resistance de Node.m
        **({"NOEUDS": str(g["NOEUDS"])} if "NOEUDS" in g else {}),
        **({"R_NOEUDS": float(g["R_NOEUDS"])} if "R_NOEUDS" in g else {}),
    )


def enregistrer(r: Reseau, nom: str | None = None) -> pathlib.Path:
    """Ecrit ``r`` dans ``g2egrid/reseaux/sauvegardes/`` et rend le chemin.

    Pas dans ``g2egrid/reseaux/`` : ce dossier-la porte les exemples natifs,
    suivis par git et couverts par les tests. Un exemple s'y depose a la main,
    deliberement ; une sauvegarde n'a rien a y faire.
    """
    dossier = dossier_sauvegardes()
    dossier.mkdir(parents=True, exist_ok=True)
    chemin = dossier / f"{nom or r.nom}{SUFFIXE}"
    chemin.write_text(
        json.dumps(reseau_vers_json(r), indent=1, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return chemin


# ----------------------------------------------------------------- les fiches
def _zoom_notebook(nom: str) -> float:
    try:
        return float(_cellule_donnees(nom).get("F_ZOOM", 5.0))
    except (FileNotFoundError, ValueError):
        return 5.0


def fiche(nom: str) -> Fiche:
    """Le resume d'un reseau, sans l'assembler ni le resoudre."""
    p = chemin_natif(nom)
    if p is not None:
        d = json.loads(p.read_text(encoding="utf-8"))
        titre, zoom = str(d.get("nom", nom)), float(d.get("f_zoom", 5.0))
        t = charger(nom)
    else:
        g = _cellule_donnees(nom)
        titre, zoom = str(g.get("NOM_RESEAU", nom)), float(g.get("F_ZOOM", 5.0))
        t = charger(nom)
    n_gfm = sum(len(getattr(t, f"T_GFM{i}s")) for i in (1, 2, 3, 4))
    n_charges = sum(1 for r in t.Y_DER if r[7] or r[8]) + len(t.T_Loads)
    return Fiche(
        nom=nom,
        titre=titre,
        bus=len(t.Y_DER),
        lignes=len(t.T_Lines_SI) + len(t.T_Lines_pu),
        transfos=len(t.T_Transfos),
        sm=len(t.T_SMs),
        gfm=n_gfm,
        gfl=len(t.T_GFL1s) + len(t.T_GFL2s),
        charges=n_charges,
        f0=t.f0,
        Sb=t.Sb,
        f_zoom=zoom,
    )


def lister() -> list[Fiche]:
    """Les fiches de tous les reseaux disponibles."""
    out = []
    for nom in RESEAUX:
        try:
            out.append(fiche(nom))
        except (FileNotFoundError, ValueError):  # notebook absent ou illisible
            continue
    return out


def construire(nom: str, resoudre: bool = True) -> Network:
    """Assemble le reseau ``nom``, load flow et equilibre compris."""
    net = assembler(charger(nom))
    if resoudre:
        net.run_loadflow()
        net.solve_equilibrium()
    return net
