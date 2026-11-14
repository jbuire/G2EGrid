"""Les balayages d'un parametre et le suivi des modes.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

from collections.abc import Sequence
import numpy as np
from scipy.optimize import linear_sum_assignment
from g2egrid.analyse import modal
from g2egrid.reseau.tables import TABLES, Tables
from ._commun import JAUGE
from .assemblage import construire


# ------------------------------------------------------------------ balayage
#: Les trois endroits ou vit un parametre, et ce qu'il faut pour le designer.
#: ``famille`` vaut pour tous les equipements d'un type, ``bus`` pour un seul
#: (``PARAMS_PAR_BUS`` l'emporte sur la famille), ``table`` pour une cellule
#: d'une table de reseau -- une reactance de ligne, une puissance de charge.
CIBLES = ("famille", "bus", "table")


def _appliquer(t: Tables, cible: dict, valeur: float) -> Tables:
    """Une COPIE des tables ou ``cible`` vaut ``valeur``.

    Une copie, et non une modification en place : le cache des reseaux
    assembles est indexe par une empreinte des tables, et muter l'original
    ferait repondre le cache pour un reseau qu'il n'a pas assemble.
    """
    genre = cible.get("genre")
    if genre not in CIBLES:
        raise ValueError(f"cible de genre {genre!r} ; attendu l'un de {CIBLES}")
    t2 = t.copie()
    cle = cible.get("cle")
    if genre == "famille":
        fam = cible.get("famille")
        if fam not in ("SM_PARAMS", "GFM_PARAMS", "GFL_PARAMS"):
            raise ValueError(f"famille inconnue : {fam!r}")
        getattr(t2, fam)[cle] = valeur
    elif genre == "bus":
        b = int(cible["bus"])
        t2.PARAMS_PAR_BUS.setdefault(b, {})[cle] = valeur
    else:
        table, ligne, col = cible["table"], int(cible["ligne"]), cible["colonne"]
        rows = getattr(t2, table, None)
        if rows is None:
            raise ValueError(f"table inconnue : {table!r}")
        if not 0 <= ligne < len(rows):
            raise ValueError(f"{table} n'a pas de ligne {ligne}")
        cols = TABLES.get(table, ())
        if col not in cols:
            raise ValueError(f"{table} n'a pas de colonne {col!r} (parmi {list(cols)})")
        # on rend une ligne du MEME type que l'originale : les tables sont
        # tantot des listes de tuples, tantot des listes de listes, et
        # l'assembleur comme l'empreinte du cache s'appuient dessus
        origine = rows[ligne]
        r = list(origine)
        r[cols.index(col)] = valeur
        rows[ligne] = tuple(r) if isinstance(origine, tuple) else r
    return t2


def _suivre(precedents: np.ndarray, courants: np.ndarray) -> np.ndarray:
    """Apparie deux spectres voisins ; rend l'ordre a donner a ``courants``.

    Sans appariement, un lieu des racines n'est qu'un nuage : ``numpy.eig``
    ne promet aucun ordre stable, et deux pas voisins d'un balayage rendent
    les memes modes dans un ordre different des que deux d'entre eux se
    croisent. On relie donc chaque mode a son plus proche voisin du pas
    precedent, par un appariement GLOBAL de cout minimal -- le glouton se
    trompe precisement la ou c'est interessant, quand deux modes se
    rencontrent puis repartent.
    """
    cout = np.abs(precedents[:, None] - courants[None, :])
    lignes, cols = linear_sum_assignment(cout)
    ordre = np.empty(len(courants), dtype=int)
    ordre[lignes] = cols
    return ordre


def _participations_par_element(P: np.ndarray, etats: list[str], seuil=0.05):
    """Les participations de chaque mode (lignes de ``P``), sommees par
    element -- le prefixe ``SM1`` de ``SM1.delta`` -- et reduites a celles
    qui depassent ``seuil``, par ordre decroissant."""
    noms = list(dict.fromkeys(e.split(".", 1)[0] for e in etats))
    colonne = {n: j for j, n in enumerate(noms)}
    M = np.zeros((len(etats), len(noms)))
    for i, e in enumerate(etats):
        M[i, colonne[e.split(".", 1)[0]]] = 1.0
    out = []
    for ligne in P @ M:
        garde = [j for j in np.argsort(-ligne) if ligne[j] >= seuil][:4]
        out.append([[noms[j], round(float(ligne[j]), 4)] for j in garde])
    return out


def balayage(
    t: Tables,
    cible: dict,
    valeurs: Sequence[float],
    f_max: float | None = None,
    n_modes: int = 12,
) -> dict:
    """Le lieu des modes quand un parametre varie.

    Rend, pour chaque valeur, le spectre reduit aux ``n_modes`` les moins
    amortis de la bande, **suivis d'un pas a l'autre** pour que chaque lieu
    soit une courbe et non un nuage (cf. :func:`_suivre`).

    Le premier pas fixe l'identite des lieux : les modes y sont tries du
    moins amorti au plus amorti, et les pas suivants s'y raccrochent. Un
    mode qui sort de la bande en cours de route laisse un trou dans son
    lieu plutot que d'aspirer le voisin le plus proche.

    Chaque mode porte aussi ses PARTICIPATIONS, regroupees par element
    (``part`` : les elements a plus de 5 %, du plus au moins participant) et
    son etat dominant : le long d'un lieu, c'est ce qui dit QUI oscille, et
    quand un mode change de nature -- un mode de machine qui devient un mode
    de convertisseur a l'approche d'une instabilite.

    Le nombre de valeurs est borne : chaque pas assemble le reseau et
    resout son equilibre, ce qui coute de l'ordre de la seconde sur un gros
    reseau. Au-dela d'une soixantaine de pas, c'est un script qu'il faut,
    pas une interface qui attend.
    """
    valeurs = [float(v) for v in valeurs]
    if not 1 <= len(valeurs) <= 60:
        raise ValueError(f"{len(valeurs)} valeurs ; attendu entre 1 et 60")

    pas: list[dict] = []
    reference: np.ndarray | None = None
    for v in valeurs:
        entree = {"valeur": v}
        try:
            a = construire(_appliquer(t, cible, v))
            ss = a.etat_stationnaire()
            tout, _, P = modal.spectre(ss)
            garde = np.abs(tout) > JAUGE
            entree["jauge"] = int((~garde).sum())
            if f_max is not None:
                garde &= modal.frequence(tout) <= f_max
            k = np.flatnonzero(garde)
            # les moins amortis d'abord : ce sont eux qu'on surveille
            k = k[np.argsort(modal.amortissement(tout[k]))][:n_modes]
            lam = tout[k]
            if reference is None or len(lam) != len(reference):
                reference = lam
            else:
                ordre = _suivre(reference, lam)
                k, lam = k[ordre], lam[ordre]
                reference = lam
            parts = _participations_par_element(P[k], list(ss.states))
            entree["modes"] = [
                {
                    "re": float(z.real),
                    "im": float(z.imag),
                    "f": float(f),
                    "zeta": float(zt),
                    "part": part,
                    "etat": ss.states[int(np.argmax(P[i]))],
                }
                for z, f, zt, part, i in zip(
                    lam,
                    modal.frequence(lam),
                    modal.amortissement(lam),
                    parts,
                    k,
                )
            ]
            entree["max_re"] = float(max(z.real for z in lam)) if len(lam) else None
        except (ValueError, KeyError, RuntimeError, np.linalg.LinAlgError) as err:
            # Un point de fonctionnement qui n'existe pas EST un resultat :
            # c'est souvent la fin du domaine de validite, et l'effacer
            # ferait croire que le balayage s'arrete sans raison.
            entree["erreur"] = str(err) or type(err).__name__
        pas.append(entree)

    return {"cible": cible, "valeurs": valeurs, "pas": pas, "n_modes": n_modes}
