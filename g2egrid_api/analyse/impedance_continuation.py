"""Impedance, boucle mineure, continuation et statismes.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import numpy as np
from g2egrid.reseau.tables import Tables, assembler
from .assemblage import construire


# ------------------------------------------------- impedance et continuation
def elements_reseau(t: Tables) -> list[dict]:
    """Les elements du reseau assemble et leurs parametres, pour les listes
    de choix (continuation, boucle mineure)."""
    net = construire(t).net
    out = []
    for el_ in net.all_elements:
        if not getattr(el_, "params", None):
            continue
        out.append(
            {
                "nom": el_.name,
                "classe": type(el_).__name__,
                "bus": [int(b) for b in (getattr(el_, "buses", None) or [])],
                "parametres": {k: float(v) for k, v in zip(el_.params, el_.p)},
                "convertisseur": el_.name.startswith(("GFM", "GFL")),
            }
        )
    return out


def _bifurcations(c) -> list[dict]:
    return [
        {
            "genre": b.genre,
            "mu": b.mu,
            "re": b.valeur.real,
            "im": b.valeur.imag,
            "f": b.f,
            "sens": b.sens,
            "elements": [{"nom": n, "part": p} for n, p in b.elements[:4]],
            "hors_limite": b.hors_limite,
            "texte": str(b),
        }
        for b in c.bifurcations
    ]


def _limites_suite(c) -> dict:
    """Les points d'une suite hors des limites des convertisseurs (FR-009) :
    masque par ``mu``, premier ``mu`` hors limite (raffine) et detail."""
    return {
        "hors_limite": c.hors_limite.tolist(),
        "mu_hors_limite": c.mu_hors_limite,
        "depassements": [[d.en_dict() for d in ds] for ds in c.depassements],
    }


def _trajectoires(c, n_modes: int = 8) -> list[dict]:
    """Les modes que le parametre DEPLACE (plus de 1 % de leur module), les
    moins amortis d'abord, un par paire ; a defaut, les moins amortis."""
    lam = c.lam
    if not lam.size:
        return []
    echelle = max(float(np.nanmax(np.abs(lam))), 1.0)
    garde = [
        j
        for j in range(lam.shape[1])
        if np.nanmin(np.abs(lam[:, j])) > 1e-9 * echelle
        and np.nanmax(lam[:, j].imag) >= -1e-9 * echelle
    ]
    bouge = [
        j
        for j in garde
        if np.nanmax(np.abs(lam[:, j] - lam[0, j])) > 0.01 * abs(lam[0, j])
    ]
    garde = bouge or garde
    garde.sort(
        key=lambda j: -np.nanmax(lam[:, j].real / np.maximum(np.abs(lam[:, j]), 1e-9))
    )
    return [
        {"re": lam[:, j].real.tolist(), "im": lam[:, j].imag.tolist()}
        for j in garde[:n_modes]
    ]


def balayage_statismes(
    t: Tables,
    debut: float = 0.25,
    fin: float = 4.0,
    pas: int = 9,
    stations: list[int] | None = None,
) -> dict:
    """Le statisme Vdc-P des stations en ``droop`` (toutes, ou ``stations``)
    multiplie par ``mu``, de ``debut`` a ``fin`` (echelle logarithmique) :
    partage de puissance, tensions continues et modes a chaque pas."""
    from g2egrid.scenarios import continuation as cont

    if debut <= 0 or fin <= 0:
        raise ValueError("facteurs de statisme strictement positifs")
    fab = cont.par_statisme(t, stations)
    c = cont.continuer(fab, np.geomspace(debut, fin, pas))
    mus = [float(m) for m in c.mu]
    mes = [fab.mesures.get(m) for m in mus]
    convertisseurs = sorted({b for m in mes if m for b in m["p_mw"]})
    noeuds = sorted({b for m in mes if m for b in m["vdc"]})
    k0 = cont.stations_statisme(t)
    return {
        "mu": mus,
        "stations": fab.stations,
        "kdroop": {str(b): k0[b] for b in fab.stations},
        "p_mw": {
            str(b): [m["p_mw"].get(b) if m else None for m in mes]
            for b in convertisseurs
        },
        "vdc": {str(b): [m["vdc"].get(b) if m else None for m in mes] for b in noeuds},
        "max_re": [float(np.nanmax(x.real[np.abs(x) > 1e-6])) for x in c.lam],
        "trajectoires": _trajectoires(c),
        "bifurcations": _bifurcations(c),
        **_limites_suite(c),
        "resume": c.resume(),
    }


def continuation_parametre(
    t: Tables,
    element: str,
    parametre: str,
    debut: float,
    fin: float,
    pas: int = 15,
    log: bool = False,
    equilibre: bool = False,
) -> dict:
    """Les modes le long d'un parametre d'element, de ``debut`` a ``fin``."""
    from g2egrid.scenarios import continuation as cont

    net = assembler(t)  # neuf : la continuation change ses parametres
    net.solve_equilibrium()
    if log:
        if debut <= 0 or fin <= 0:
            raise ValueError("echelle logarithmique : bornes strictement positives")
        valeurs = np.geomspace(debut, fin, pas)
    else:
        valeurs = np.linspace(debut, fin, pas)
    c = cont.continuer(cont.par_parametre(net, element, parametre, equilibre), valeurs)
    return {
        "mu": c.mu.tolist(),
        "max_re": [float(np.nanmax(x.real[np.abs(x) > 1e-6])) for x in c.lam],
        "trajectoires": _trajectoires(c),
        "bifurcations": _bifurcations(c),
        **_limites_suite(c),
        "resume": c.resume(),
    }


def boucle_mineure(
    t: Tables, element: str, f_min: float = 0.1, f_max: float = 5000.0, n: int = 800
) -> dict:
    """Le critere de la boucle mineure d'un convertisseur contre le reste du
    reseau (voir :mod:`g2egrid.analyse.impedance`)."""
    from g2egrid.analyse import impedance

    net = assembler(t)
    net.solve_equilibrium()
    f = np.geomspace(f_min, f_max, n)
    b = impedance.boucle_mineure(net, element, f)
    # eigvals ne garde aucun ordre d'une frequence a l'autre : on suit les
    # deux lieux par continuite, sinon le trace saute d'une branche a l'autre
    lieux = b.nyquist["lieux"].copy()
    for k in range(1, len(f)):
        a, c = lieux[k - 1], lieux[k]
        if abs(a[0] - c[0]) + abs(a[1] - c[1]) > abs(a[0] - c[1]) + abs(a[1] - c[0]):
            lieux[k] = c[::-1]
    sv = lambda M: np.linalg.svd(M, compute_uv=False)  # noqa: E731
    return {
        "element": element,
        "bus": b.bus,
        "forme": b.forme,
        "f": f.tolist(),
        "lieux": [
            {"re": lieux[:, i].real.tolist(), "im": lieux[:, i].imag.tolist()}
            for i in range(lieux.shape[1])
        ],
        "z_reste": sv(b.Z_r).max(axis=1).tolist(),
        "z_equipement": (1 / sv(b.Y_e).min(axis=1)).tolist(),
        "passivite": impedance.passivite(b.Y_e).tolist(),
        "non_passif": b.non_passif,
        "tours": b.nyquist["tours"],
        "instables": b.instables,
        "applicable": b.applicable,
        "stable": b.stable,
        "marge": b.nyquist["marge"],
        "f_marge": b.nyquist["f_marge"],
        "resume": b.resume(),
    }
