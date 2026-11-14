"""Transport et distribution : impact, comparaisons, penetration.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

from collections.abc import Sequence
import numpy as np
from g2egrid.analyse import modal
from g2egrid.multi_niveau.impact import impact_distribution
from g2egrid.scenarios.robustesse import etudier as etudier_robustesse
from g2egrid.scenarios import simplification as simp
from g2egrid.reseau.tables import Tables, tables_vers_json
from g2egrid.reseau.zones import aplatir
from ._commun import JAUGE
from .catalogue_api import _positions, _reponse_reseau
from .fichiers_reseaux import _reseau_depuis


# --------------------------------------------------- transport / distribution
def robustesse(
    t: Tables,
    f_min: float = 0.05,
    f_max: float = 5.0,
    zeta_max: float | None = 30.0,
    frac: float = 0.1,
    puissance: bool = True,
    deconnexion: bool = True,
) -> dict:
    """Les modes les plus sensibles aux variations de puissance et au N-1.

    Voir :mod:`g2egrid.scenarios.robustesse`.
    """
    e = etudier_robustesse(t, f_min, f_max, zeta_max, frac, puissance, deconnexion)
    pert = lambda p: {"genre": p.genre, "libelle": p.libelle}  # noqa: E731
    effet = lambda x: {  # noqa: E731
        **pert(x.perturbation),
        "disparu": x.disparu,
        "re": None if x.disparu else float(x.lam.real),
        "im": None if x.disparu else float(x.lam.imag),
        "f": None if x.disparu else float(modal.frequence(x.lam)[0]),
        "zeta": None if x.disparu else float(modal.amortissement(x.lam)[0]),
        "df": x.df,
        "dzeta": x.dzeta,
        "dsigma": x.dsigma,
        "score": x.score,
        "hors_limite": x.hors_limite,
    }
    return {
        "nom": e.nom,
        "frac": e.frac,
        "n_variantes": e.n_variantes,
        # FR-009 : le cas de base et les variantes hors des limites des
        # convertisseurs (leurs effets restent dans les scores, marques)
        "depassements": [d.en_dict() for d in e.depassements],
        "hors_limite": [
            {**pert(p), "depassements": [d.en_dict() for d in ds]}
            for p, ds in e.hors_limite
        ],
        "modes": [
            {
                "re": m.lam.real,
                "im": m.lam.imag,
                "f": m.f,
                "zeta": m.zeta,
                "dominant": m.dominant,
                "elements": [{"nom": n, "part": v} for n, v in m.elements],
                "score_puissance": m.score("puissance"),
                "score_deconnexion": m.score("deconnexion"),
                "zeta_min": m.zeta_min,
                "pire": effet(m.pire) if m.pire else None,
                "disparitions": m.disparitions,
                "effets": [effet(x) for x in m.effets],
            }
            for m in e.modes
        ],
        "instables": [
            {
                **pert(p),
                "max_re": re,
                "elements": [{"nom": n, "part": v} for n, v in qui],
            }
            for p, re, qui in e.instables
        ],
        "echecs": [{**pert(p), "raison": r} for p, r in e.echecs],
    }


def simplification_inventaire(t: Tables) -> dict:
    """Ce que l'onglet des hypotheses de simplification propose : chaque type
    d'element du reseau, ses etats et ses elements, et les preregles (ceux
    qui s'appliquent a ce reseau). Voir :mod:`g2egrid.scenarios.simplification`."""
    inv = simp.inventaire(t)
    preregles = []
    for nom, (description, _) in simp.PREREGLES.items():
        s = simp.preregle(t, nom)
        if s:
            preregles.append(
                {"nom": nom, "description": description, "simplifications": s}
            )
    return {
        "types": [
            {
                "type": typ,
                "etats": d["etats"],
                "elements": [
                    {"nom": n, "dynamiques": dyn} for n, dyn in d["elements"].items()
                ],
            }
            for typ, d in inv.items()
        ],
        "preregles": preregles,
        "simplifications": t.SIMPLIFICATIONS,
    }


def simplification_comparer(
    t: Tables, simplifications: dict, preregles: Sequence[str] = ()
) -> dict:
    """Les modes du modele COMPLET (sans aucune simplification) contre ceux
    du modele simplifie par ``simplifications`` et les ``preregles``."""
    s = simp.fusionner(simp.preregle(t, *preregles), simplifications)
    c = simp.comparer(t.copie(SIMPLIFICATIONS={}), s)
    cplx = lambda lam: (  # noqa: E731
        (None, None) if lam is None else (float(lam.real), float(lam.imag))
    )
    modes = []
    for m in c.modes:
        (re, im), (re_r, im_r) = cplx(m.lam), cplx(m.lam_reduit)
        modes.append(
            {
                "statut": m.statut,
                "re": re,
                "im": im,
                "f": m.f,
                "zeta": m.zeta,
                "re_r": re_r,
                "im_r": im_r,
                "f_r": m.f_reduit,
                "zeta_r": m.zeta_reduit,
                "ecart": m.ecart,
                "dzeta": m.dzeta,
                "part_retiree": m.part_retiree,
                "critique": m.critique(),
                "elements": [{"nom": n, "part": v} for n, v in m.elements],
            }
        )
    return {
        "simplifications": c.simplifications,
        "n_etats": c.n_etats,
        "n_etats_reduit": c.n_etats_reduit,
        "max_re": c.max_re,
        "max_re_reduit": c.max_re_reduit,
        "stable": c.stable,
        "stable_reduit": c.stable_reduit,
        "ecart_max": c.ecart_max(),
        "ecart_max_em": c.ecart_max(f_max=5.0),
        "dzeta_max": c.dzeta_max(),
        "comptes": {
            k: len(getattr(c, k))
            for k in ("retrouves", "retires", "perdus", "nouveaux")
        },
        "alertes": c.alertes(),
        "sure": c.sure(),
        "modes": modes,
    }


def impact(reseau: dict, f_max: float = 5.0, par_poste: bool = True) -> dict:
    """L'impact de la dynamique de la distribution sur le transport.

    Voir :mod:`g2egrid.multi_niveau.impact`. Le reseau est pris sous n'importe laquelle des
    formes qu'accepte :func:`lire_fichier` ; il doit avoir des zones.
    """
    r = _reseau_depuis(reseau)
    if not r.instances:
        raise ValueError(
            "ce reseau n'a pas de zone : l'impact de la distribution se mesure "
            "sur un reseau decrit par zones (transport + instances)"
        )
    im = impact_distribution(r, f_max=f_max, par_poste=par_poste)
    mode = lambda m: {  # noqa: E731
        "f_eq": m.f_equivalent,
        "zeta_eq": m.zeta_equivalent,
        "f": m.f_complet,
        "zeta": m.zeta_complet,
        "ecart_f": m.ecart_f,
        "ecart_zeta": m.ecart_zeta,
        "part_distribution": m.part_distribution,
        "dominant": m.dominant,
        "re_eq": m.lam_equivalent.real,
        "im_eq": m.lam_equivalent.imag,
        "re": m.lam_complet.real,
        "im": m.lam_complet.imag,
    }
    return {
        "ecart_tension": im.ecart_tension,
        "max_re_equivalent": im.max_re_equivalent,
        "max_re_complet": im.max_re_complet,
        "soutirage": {k: {"P": v.real, "Q": v.imag} for k, v in im.soutirage.items()},
        "paires": [mode(m) for m in im.paires],
        "degrades": [mode(m) for m in im.degrades],
        "crees": [
            {
                "f": c.f,
                "zeta": c.zeta,
                "part_distribution": c.part_distribution,
                "dominant": c.dominant,
                "re": c.lam.real,
                "im": c.lam.imag,
            }
            for c in im.crees
        ],
        "par_poste": {
            k: [{"f_eq": a, "ecart_f": b, "ecart_zeta": c} for a, b, c in v]
            for k, v in im.par_poste.items()
        },
    }


# ------------------------------------------------ transport et distribution
def _comparaison_json(c, f_max: float = 100.0) -> dict:
    modes = [m for m in c.modes if m.dans_bande(f_max)]
    pire = c.pire(f_max)
    return {
        "etude": c.etude,
        "ecart_max": c.ecart_max(f_max),
        "dzeta_max": c.dzeta_max(f_max),
        "manquants": len(c.absents(f_max)),
        "pire": None if pire is None else _mode_compare_json(pire),
        "modes": [_mode_compare_json(m) for m in modes],
    }


def _mode_compare_json(m) -> dict:
    return {
        "f": m.f,
        "zeta": m.zeta,
        "re": m.lam_complet.real,
        "im": m.lam_complet.imag,
        "f_reduit": m.f_reduit,
        "zeta_reduit": m.zeta_reduit,
        "re_reduit": None if m.lam_reduit is None else m.lam_reduit.real,
        "im_reduit": None if m.lam_reduit is None else m.lam_reduit.imag,
        "ecart": m.ecart,
        "dzeta": m.dzeta,
        "part": m.part,
        "elements": [{"nom": n, "part": v} for n, v in m.elements],
    }


def multi_ajouter_gfl(
    reseau: dict, n: int = 5, graine: int | None = None, part: float = 0.05
) -> dict:
    """Ajoute ``n`` GFL a des bus tires au hasard dans les distributions
    (voir :func:`g2egrid.multi_niveau.multi_echelle.ajouter_gfl`)."""
    from g2egrid.multi_niveau.multi_echelle import ajouter_gfl

    r = _reseau_depuis(reseau)
    if not r.instances:
        raise ValueError("il faut un reseau a zones (des distributions instanciees)")
    r2, ajouts = ajouter_gfl(r, n, graine, part)
    return _reponse_reseau(r2) | {
        "ajouts": [{"zone": z, "bus": b, "P_W": p} for z, b, p in ajouts]
    }


def multi_comparer(reseau: dict, instances: list[int] | None = None) -> dict:
    """Les deux modeles reduits confrontes au modele complet : le transport
    avec ses distributions en charges (R-L, puis hybride), chaque
    distribution seule derriere un bus infini et son impedance de Thevenin,
    et les modes d'interaction. Voir :mod:`g2egrid.multi_niveau.multi_echelle`."""
    from g2egrid.multi_niveau.multi_echelle import (
        _production_gfl,
        comparer_distribution,
        comparer_partition,
        comparer_transport,
        modes_interaction,
        partition,
    )

    r = _reseau_depuis(reseau)
    if not r.instances:
        raise ValueError("il faut un reseau a zones (des distributions instanciees)")
    instances = list(range(len(r.instances))) if instances is None else instances
    plat = aplatir(r)
    P, complet = partition(r, plat)
    lam = complet.lam[np.abs(complet.lam) > JAUGE]
    return {
        # dx/dt = Ax + Bu partage en x_T, x_D : l'autre cote dans les entrees,
        # fige ou quasi statique (complement de Schur)
        "partition": {
            f"{cote}_{'quasi_statique' if qs else 'fige'}": _comparaison_json(
                comparer_partition(r, cote, qs, part=(P, complet))
            )
            for cote in ("transport", "distribution")
            for qs in (False, True)
        },
        "p_gfl": _production_gfl(r),
        "max_re": float(lam.real.max()) if lam.size else 0.0,
        "transport_rl": _comparaison_json(
            comparer_transport(r, plat=plat, complet=complet)
        ),
        "transport_hybride": _comparaison_json(
            comparer_transport(r, plat=plat, complet=complet, hybride=True)
        ),
        "distributions": [
            {
                "instance": f"{r.instances[k].zone}@{r.instances[k].bus_amont}",
                **_comparaison_json(
                    comparer_distribution(r, k, plat=plat, complet=complet)
                ),
            }
            for k in instances
        ],
        "interaction": [
            {
                "f": m.f,
                "zeta": m.zeta,
                "re": m.lam.real,
                "im": m.lam.imag,
                "part_transport": m.part_transport,
                "elements": [{"nom": n, "part": v} for n, v in m.elements],
            }
            for m in modes_interaction(r, complet, plat)
        ],
    }


def multi_penetration(
    reseau: dict, n_max: int = 30, pas: int = 10, graine: int = 1, part: float = 0.05
) -> dict:
    """Le balayage de penetration : des GFL ajoutes par paquets, et a chaque
    palier l'ecart des modeles reduits au complet."""
    from g2egrid.multi_niveau.multi_echelle import penetration

    r = _reseau_depuis(reseau)
    if not r.instances:
        raise ValueError("il faut un reseau a zones (des distributions instanciees)")
    if n_max > 200 or pas < 1:
        raise ValueError("n_max <= 200 et pas >= 1")
    out = []
    for pal in penetration(r, n_max, pas, graine, part):
        e = {"n_gfl": pal.n_gfl, "p_gfl": pal.p_gfl, "erreur": pal.erreur}
        if not pal.erreur:
            e |= {
                "max_re": pal.max_re,
                "transport_rl": {
                    "ecart": pal.transport.ecart_max(),
                    "dzeta": pal.transport.dzeta_max(),
                },
                "transport_hybride": {
                    "ecart": pal.transport_hybride.ecart_max(),
                    "dzeta": pal.transport_hybride.dzeta_max(),
                },
                "distribution": {
                    "ecart": max(
                        (c.ecart_max() for c in pal.distribution),
                        default=0.0,
                    ),
                    "dzeta": max(
                        (c.dzeta_max() for c in pal.distribution),
                        default=0.0,
                    ),
                },
                "interaction": len([m for m in pal.interaction if m.zeta < 30]),
            }
        out.append(e)
    return {"paliers": out}


def definir_zone(reseau: dict, bus: list[int], nom: str) -> dict:
    """Les bus ``bus`` du transport deviennent une zone (voir
    :func:`g2egrid.reseau.zones.extraire_zone`), rendue sous la forme de LECTURE pour
    que le constructeur entre aussitot en mode zones."""
    from g2egrid.reseau.zones import extraire_zone

    r = extraire_zone(_reseau_depuis(reseau), bus, nom)
    aplatir(r)  # ses garde-fous jouent avant de rendre la main
    return _reponse_reseau(r)


def aplatir_reseau(reseau: dict) -> dict:
    """Un reseau structure, rendu A PLAT : les tables que l'assembleur lit.

    C'est le meme aplatissement que celui du catalogue, offert sur un reseau
    qui n'est pas enregistre : le constructeur peut donc montrer ce que ses
    zones donnent -- et le faire calculer -- sans avoir a ecrire un fichier
    d'abord.
    """
    r = _reseau_depuis(reseau)
    t = aplatir(r)
    return {"nom": t.nom, "tables": tables_vers_json(t), "positions": _positions(t)}
