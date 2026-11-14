"""Les exports (MATLAB, notebook, tables).

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import pathlib
import re
import tempfile
from g2egrid.formats import notebook
from g2egrid.formats.export import exporter, vers_matlab
from g2egrid.reseau.tables import Tables
from g2egrid.reseau.zones import Reseau
from .catalogue_api import _reponse_reseau


# ------------------------------------------------------------------ exports
#: Les trois sorties, et ce qu'elles valent. ``donnees`` est la cellule de
#: tables a coller dans un notebook, ``notebook`` le notebook complet,
#: ``matlab`` le fichier reseau de buildFromMatrices.
EXPORTS = ("donnees", "notebook", "matlab", "matpower")


def exporter_vers(t: Tables, quoi: str) -> dict:
    """Le reseau dans l'un des trois formats, ET CE QU'IL Y PERD.

    Le champ ``avertissements`` n'est pas un ornement. L'export ``.m`` laisse
    en commentaire ce que ``buildFromMatrices`` ne sait pas prendre -- les
    shunts, les charges a courant ou a puissance constante, les reglages fins
    des machines : le fichier produit reste executable, mais il decrit un
    reseau plus pauvre que celui qu'on a dessine. Avaler ces messages
    donnerait un export qui a l'air complet et qui ne l'est pas, ce qui est
    exactement la panne qu'on ne diagnostique pas.

    Le refus du continu, lui, n'est pas un avertissement mais une
    ``ValueError`` de :func:`g2egrid.formats.export.vers_matlab` : un reseau ampute de
    sa liaison continue n'est pas une approximation du reseau, c'est un autre
    reseau qui converge et qui ment. Elle remonte en 422 avec son message.
    """
    if quoi not in EXPORTS:
        raise ValueError(f"export '{quoi}' inconnu ; attendu l'un de {EXPORTS}")
    if quoi == "donnees":
        return {
            "quoi": quoi,
            "fichier": f"{_nom_de_sortie(t)}.py",
            "type": "text/x-python",
            "contenu": exporter(t),
            "avertissements": [],
        }
    if quoi == "matlab":
        source, avertissements = vers_matlab(t)
        return {
            "quoi": quoi,
            "fichier": f"{_nom_matlab(source, _nom_de_sortie(t))}.m",
            "type": "text/x-matlab",
            "contenu": source,
            "avertissements": avertissements,
        }
    if quoi == "matpower":
        from g2egrid.formats.interop import vers_matpower

        avert = []
        if t.T_BusDC:
            avert.append("le reseau continu n'a pas d'equivalent MATPOWER : omis")
        avert.append(
            "MATPOWER ne porte que le load flow : la dynamique (reglages des "
            "machines et convertisseurs) n'est pas exportee"
        )
        return {
            "quoi": quoi,
            "fichier": f"{_nom_de_sortie(t)}_mpc.m",
            "type": "text/x-matlab",
            "contenu": vers_matpower(t),
            "avertissements": avert,
        }
    return _export_notebook(t)


def importer_externe(nom: str, contenu: str, charges: str = "Z") -> dict:
    """Un fichier MATPOWER (``.m``), PSS/E (``.raw``) ou pandapower (JSON),
    lu par :mod:`g2egrid.formats.interop`, rendu sous la forme de lecture, avec des
    positions pour le schema et ce que l'import n'a pas pu reprendre."""
    from g2egrid.formats import interop

    ext = pathlib.Path(nom).suffix.lower()
    base = pathlib.Path(nom).name.split(".")[0] or "importe"
    if ext == ".m":
        t = interop.depuis_matpower(contenu, nom=base, charges=charges)
    elif ext == ".raw":
        t = interop.depuis_psse(contenu, nom=base, charges=charges)
    elif ext == ".json" and '"pandapower' in contenu[:2000]:
        t = interop.depuis_pandapower(contenu, nom=base, charges=charges)
    else:
        raise ValueError(
            f"{nom} : format non reconnu (MATPOWER .m, PSS/E .raw, pandapower .json)"
        )
    aretes = [(int(r[1]), int(r[2])) for k in TABLES_BRANCHES for r in getattr(t, k)]
    pos = interop.disposition([int(r[0]) for r in t.Y_DER], aretes)
    r = Reseau(transport=t, nom=t.nom, positions=pos)
    return {**_reponse_reseau(r), "avertissements": list(t.avertissements)}


TABLES_BRANCHES = ("T_Lines_SI", "T_Lines_pu", "T_Transfos")


def _nom_de_sortie(t: Tables) -> str:
    """Le nom du reseau, ramene a ce qu'un nom de fichier accepte."""
    ident = re.sub(r"[^0-9A-Za-z_-]", "_", t.nom or "reseau").strip("_")
    return ident or "reseau"


def _nom_matlab(source: str, defaut: str) -> str:
    """Le nom de la fonction MATLAB, LU dans le source plutot que refait.

    ``export._nom_fonction`` est prive, et en recopier la regle donnerait deux
    assainissements a tenir d'accord. Or MATLAB exige que le fichier porte le
    nom de la fonction qu'il contient : les voir diverger donnerait un
    ``.m`` que MATLAB refuse d'appeler.
    """
    m = re.match(r"function\s+net\s*=\s*([0-9A-Za-z_]+)", source)
    return m.group(1) if m else defaut


def _export_notebook(t: Tables) -> dict:
    """Le notebook complet d'un reseau, et ses cellules une a une.

    Les deux, et pas l'un ou l'autre. ``contenu`` est le ``.ipynb`` tel que
    :func:`g2egrid.formats.notebook.ecrire` le pose sur le disque -- donc exactement ce
    qu'on obtiendrait en l'ecrivant soi-meme, sans que l'interface ait a
    connaitre le format nbformat. Mais ``nbformat`` n'est installe que pour le
    developpement : sur une installation qui ne l'a pas, rendre une 500 pour un
    export serait absurde alors que les CELLULES, elles, se fabriquent sans
    lui. On rend donc toujours les cellules, le fichier quand c'est possible,
    et un avertissement quand ca ne l'est pas.
    """
    nom = _nom_de_sortie(t)
    donnees = notebook.cellule_donnees(t)
    intro = _intro_notebook(t)
    cellules = [
        {"genre": genre, "source": source}
        for genre, source in notebook.cellules(intro, donnees)
    ]
    contenu, avertissements = "", []
    try:
        with tempfile.TemporaryDirectory() as tmp:
            chemin = notebook.ecrire(pathlib.Path(tmp) / f"{nom}.ipynb", intro, donnees)
            contenu = chemin.read_text(encoding="utf-8")
    except ImportError as err:
        avertissements.append(
            f"notebook non assemble : {err}. Les cellules sont rendues telles "
            "quelles ; installer nbformat (uv sync) rend aussi le .ipynb."
        )
    return {
        "quoi": "notebook",
        "fichier": f"{nom}.ipynb",
        "type": "application/x-ipynb+json",
        "contenu": contenu,
        "cellules": cellules,
        "avertissements": avertissements,
    }


def _intro_notebook(t: Tables) -> str:
    """Le markdown de tete : ce que le reseau contient, en un coup d'oeil.

    Sans date ni horodatage, volontairement : un notebook reexporte sans
    changer le reseau doit redonner le MEME fichier, sinon tout reexport
    apparait comme une modification.
    """
    gfm = sum(len(getattr(t, f"T_GFM{i}s")) for i in (1, 2, 3, 4))
    gfl = len(t.T_GFL1s) + len(t.T_GFL2s)
    lignes = [
        f"# {t.nom}",
        "",
        "Reseau exporte de l'interface g2egrid.",
        "",
        f"- {len(t.Y_DER)} bus, "
        f"{len(t.T_Lines_SI) + len(t.T_Lines_pu)} ligne(s), "
        f"{len(t.T_Transfos)} transfo(s)",
        f"- {len(t.T_SMs)} machine(s) synchrone(s), {gfm} grid-forming, "
        f"{gfl} grid-following, {len(t.T_IBs)} bus infini(s)",
        f"- f0 = {t.f0:g} Hz, Sb = {t.Sb / 1e6:g} MVA",
    ]
    if t.T_BusDC:
        lignes.append(
            f"- partie continue : {len(t.T_BusDC)} bus, "
            f"{len(t.T_LinesDC)} cable(s), {len(t.T_VSCs)} convertisseur(s)"
        )
    return "\n".join(lignes)


# ------------------------------------------------------------------ export
# ``exporter`` vit dans ``g2egrid.formats.export`` et n'est plus que reexporte ici.
# Le sens de la dependance comptait : ``g2egrid.formats.notebook`` avait besoin de la
# cellule de donnees et allait la chercher dans ``g2egrid_api``, si bien que
# le PAQUET dependait de l'INTERFACE -- et devenait inutilisable sur une
# installation sans l'extra ``api``. Une couche de transport peut dependre
# du modele ; l'inverse est une inversion qui se paie tot ou tard.
