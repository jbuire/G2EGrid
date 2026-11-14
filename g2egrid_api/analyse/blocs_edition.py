"""L'edition des schemas-blocs, sans HTTP : sous-systemes, element de la
bibliotheque en blocs, bibliotheque de modeles.

Ce que l'interface fait d'un schema (regrouper, degrouper, aplatir, ouvrir un
GFL en blocs, ranger un modele) se fait aussi en code : ces fonctions ne
connaissent ni FastAPI ni les codes d'erreur, elles levent
:class:`~g2egrid.elements.blocs.SchemaInvalide` (une ``ValueError``) avec un
message en francais, que le transport rend en 422.

LA BIBLIOTHEQUE de modeles est un dossier de fichiers JSON, un par modele,
dans ``modeles/`` sous le dossier des sauvegardes du catalogue
(:func:`g2egrid.reseau.catalogue.dossier_sauvegardes`) : partagee entre les
reseaux, elle survit a la session. Le nom du fichier est le ``nom`` du schema,
valide comme celui d'un reseau (nom de fichier nu : lettres, chiffres, ``_``,
``-``, 64 caracteres au plus, pas de chemin).
"""

from __future__ import annotations

import datetime
import json
import os
import pathlib
import tempfile

from g2egrid.elements import blocs, blocs_equations
from g2egrid.elements.blocs import SchemaInvalide
from g2egrid.reseau import catalogue
from g2egrid.reseau.tables import Tables, assembler

from .resultats import _nom_element


def blocs_regrouper(schema: dict, ids: list[str], nom: str, chemin=()) -> dict:
    """Le schema ou les blocs ``ids`` sont ranges dans un sous-systeme ``nom``
    (voir :func:`g2egrid.elements.blocs.regrouper`) : ``{"schema": ...}``."""
    return {"schema": blocs.regrouper(schema, ids, nom, chemin=list(chemin or ()))}


def blocs_degrouper(schema: dict, id: str, chemin=()) -> dict:
    """Le schema ou le sous-systeme ``id`` est dissous : ``{"schema": ...}``."""
    return {"schema": blocs.degrouper(schema, id, chemin=list(chemin or ()))}


def blocs_aplatir(schema: dict) -> dict:
    """Le schema sans sous-systeme : ``{"schema": ...}``."""
    return {"schema": blocs.aplatir(schema)}


def blocs_depuis_element(
    classe: str | None = None,
    parametres: dict | None = None,
    detail: str = "blocs",
    t: Tables | None = None,
    table: str | None = None,
    index: int | None = None,
    f0: float | None = None,
) -> dict:
    """Le schema-bloc d'un element de la bibliotheque : ``{"schema": ...}``.

    Soit une CLASSE native (``classe``, ``parametres``, ``f0``), soit un
    element d'un reseau (``t``, ``table``, ``index``) avec SES parametres."""
    if table is not None:
        if t is None or index is None:
            raise SchemaInvalide(
                "un element de reseau demande ses tables et son indice"
            )
        m = _schema_element_du_reseau(t, table, index, detail)
        return {"schema": m, "modele": m} if detail == "liste" else {"schema": m}
    if not classe:
        raise SchemaInvalide(
            "il faut une classe native (classe) ou un element de reseau "
            "(tables, table, index)"
        )
    kw = {} if f0 is None else {"f0": float(f0)}
    m = blocs.depuis_element(classe, parametres, detail, **kw)
    return {"schema": m, "modele": m} if detail == "liste" else {"schema": m}


def _schema_element_du_reseau(t: Tables, table: str, index: int, detail: str) -> dict:
    from g2egrid.elements.blocs_inverse import depuis_instance
    from g2egrid.reseau.tables import TABLES

    if table not in TABLES:
        raise SchemaInvalide(
            f"table inconnue : {table!r} (tables connues : {', '.join(TABLES)})"
        )
    if table == "T_Blocs":
        raise SchemaInvalide(
            "cet element est deja decrit par un schema-bloc (T_Blocs) : ouvrir "
            "son modele dans l'onglet Blocs"
        )
    n = len(getattr(t, table))
    if not 0 <= index < n:
        raise SchemaInvalide(
            f"{table} n'a pas de ligne d'indice {index} : "
            + (f"l'indice va de 0 a {n - 1}" if n else "la table est vide")
        )
    net = assembler(t)
    noms = _nom_element(t, table, index)
    par_nom = {e.name: e for e in net.all_elements}
    el = next((par_nom[n] for n in noms if n in par_nom), None)
    if el is None:
        raise SchemaInvalide(f"{table}[{index}] : aucun modele dynamique a decrire")
    return depuis_instance(el, detail)


# ------------------------------------------------------------ bibliotheque
def dossier_bibliotheque() -> pathlib.Path:
    """Le dossier des modeles, calcule a l'appel (il suit un catalogue redirige)."""
    return catalogue.dossier_sauvegardes() / "modeles"


def _nom_modele(nom) -> str:
    try:
        return catalogue._verifier_nom(nom)
    except ValueError:
        raise SchemaInvalide(
            f"nom de modele invalide : {nom!r} (lettres, chiffres, _ et - ; 64 "
            "caracteres au plus, sans chemin)"
        ) from None


def _chemin(nom: str) -> pathlib.Path:
    return dossier_bibliotheque() / f"{_nom_modele(nom)}.json"


#: taille maximale d'un modele range (octets de son JSON) : un GFM complet en
#: blocs en fait 60 000 ; au-dela, ce n'est plus un modele dessine a la main
TAILLE_MODELE_MAX = 2_000_000


class ModeleAbsent(FileNotFoundError):
    """Un modele qui n'est pas dans la bibliotheque (le transport en fait une 404)."""


def bibliotheque_liste() -> list[dict]:
    """``[{nom, titre, modifie}]``, par nom : ``titre`` est le titre du schema
    (son nom a defaut), ``modifie`` la date de derniere modification (ISO UTC,
    a la seconde)."""
    d = dossier_bibliotheque()
    if not d.is_dir():
        return []
    out = []
    for p in sorted(d.glob("*.json")):
        nom = p.name[: -len(".json")]
        try:
            _nom_modele(nom)
            schema = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue  # un fichier etranger ou abime ne casse pas la liste
        titre = schema.get("titre") if isinstance(schema, dict) else None
        out.append(
            {
                "nom": nom,
                "titre": str(titre or nom),
                "modifie": datetime.datetime.fromtimestamp(
                    p.stat().st_mtime, tz=datetime.timezone.utc
                )
                .replace(microsecond=0)
                .isoformat(),
            }
        )
    return out


def bibliotheque_lire(nom: str) -> dict:
    """``{"schema": ...}`` du modele ``nom``."""
    p = _chemin(nom)
    if not p.is_file():
        raise ModeleAbsent(f"modele {nom!r} absent de la bibliotheque")
    return {"schema": json.loads(p.read_text(encoding="utf-8"))}


def bibliotheque_enregistrer(schema: dict) -> dict:
    """Range ``schema`` sous son ``nom`` (remplace le modele du meme nom) :
    ``{"nom": ...}``. Le schema n'a pas a compiler (un brouillon se range),
    mais il doit avoir un nom valide et la forme d'un schema."""
    if blocs_equations.est_modele_equations(schema):
        # un modele a equations se range aussi (il s'ouvre dans l'onglet Blocs)
        blocs_equations.verifier_forme_equations(schema)
    else:
        if not isinstance(schema, dict) or not isinstance(schema.get("blocs"), list):
            raise SchemaInvalide(
                "ce n'est pas un modele : un schema-bloc (liste de blocs) ou un "
                'modele a equations ("genre": "equations") est attendu'
            )
        # un fichier mal forme ne doit pas entrer : il casserait a l'ouverture
        blocs.verifier_forme(schema)
    nom = _nom_modele(schema.get("nom"))
    d = dossier_bibliotheque()
    d.mkdir(parents=True, exist_ok=True)
    try:
        # NaN et inf ne sont pas du JSON : le fichier ne se relirait pas ailleurs
        texte = json.dumps(
            schema, ensure_ascii=False, indent=1, sort_keys=False, allow_nan=False
        )
    except ValueError:
        raise SchemaInvalide(
            "le schema contient une valeur non finie (NaN ou infini)"
        ) from None
    if len(texte.encode("utf-8")) > TAILLE_MODELE_MAX:
        raise SchemaInvalide(
            f"modele trop gros pour la bibliotheque ({len(texte) // 1000} ko, "
            f"au plus {TAILLE_MODELE_MAX // 1000} ko)"
        )
    # ecriture atomique : un lecteur ne voit jamais un fichier a moitie ecrit
    fd, tmp = tempfile.mkstemp(dir=d, suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(texte)
        os.replace(tmp, d / f"{nom}.json")
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return {"nom": nom}


def bibliotheque_supprimer(nom: str) -> dict:
    """Retire le modele ``nom`` de la bibliotheque : ``{"nom": ...}``."""
    p = _chemin(nom)
    if not p.is_file():
        raise ModeleAbsent(f"modele {nom!r} absent de la bibliotheque")
    p.unlink()
    return {"nom": nom}
