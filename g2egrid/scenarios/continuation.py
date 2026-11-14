"""Continuation : suivre les modes le long d'un parametre, trouver ou ils
traversent l'axe imaginaire.

Une etude modale repond pour UN point de fonctionnement ; la question qu'on se
pose en pratique est « jusqu'ou ? » -- jusqu'a quel gain de PLL, quelle
charge, quelle puissance d'un parc le reseau reste stable, et quel mode le
lache en premier. On fait donc varier un parametre ``mu`` le long d'une suite
de valeurs, et a chaque pas :

* on reconstruit le modele (:func:`par_parametre`, :func:`par_statisme`, ou
  toute fonction ``mu -> Network`` resolu) et son spectre ;
* on APPARIE chaque mode a celui du pas precedent (algorithme hongrois sur la
  distance des valeurs propres et la ressemblance des participations, comme
  :mod:`g2egrid.scenarios.robustesse`) : un mode garde son indice d'un bout a l'autre,
  meme quand deux modes se croisent ;
* un mode dont la partie reelle change de signe entre deux pas est une
  BIFURCATION, raffinee par dichotomie : de HOPF pour une paire complexe (une
  oscillation nait, a la frequence ``Im lam / 2 pi``), NOEUD-COL pour une
  valeur propre reelle (le point d'equilibre disparait ou change de nature) ;
* un pas qui n'a plus de point d'equilibre (le load flow ne converge plus)
  ARRETE la suite au dernier point converge : ``Continuation.arret_a`` garde
  le ``mu`` refuse, et :meth:`Continuation.resume` le dit. La suite ne
  cherche pas la frontiere : elle s'arrete la ou le modele cesse d'exister.

Les tensions des bus sont gardees a chaque pas.

Le mode de jauge (valeur propre nulle d'un reseau sans bus infini) n'est pas
une bifurcation : on l'ecarte, comme partout ailleurs.

LES LIMITES DES CONVERTISSEURS sont verifiees a chaque pas
(:mod:`g2egrid.noyau.limites`). Un point ou un convertisseur depasse son courant
maximal (ou sa modulation maximale) est MARQUE, pas retire : la suite
continue, ``Continuation.depassements`` dit lesquels a chaque ``mu``,
``Continuation.mu_hors_limite`` donne le premier ``mu`` hors limite (raffine
par dichotomie), et chaque bifurcation porte ``hors_limite``. Pourquoi ne pas
s'arreter ? Parce que sans point d'equilibre la suite n'a plus rien a suivre,
alors qu'au-dela de ``Imax`` il en reste un -- celui du modele libre, qui dit ce que
ferait un convertisseur plus gros, et ou la commande perdrait sa stabilite
s'il l'etait : c'est une information de dimensionnement, que l'arret
cacherait. Ce point ne decrit plus le materiel installe, et le resultat le
dit. Un seul avertissement :class:`~g2egrid.noyau.limites.DepassementLimite` est
emis, a la fin de la suite, au lieu d'un par pas.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field, replace
from typing import Callable

import numpy as np
from scipy.optimize import linear_sum_assignment

from g2egrid.noyau import limites
from g2egrid.analyse import modal
from g2egrid.noyau.progres import signaler

__all__ = [
    "Bifurcation",
    "Continuation",
    "continuer",
    "par_parametre",
]

#: En dessous (relatif a max|lam|), une valeur propre est le mode de jauge.
JAUGE = 1e-9


@dataclass
class Bifurcation:
    genre: str  # "hopf" ou "noeud-col"
    mu: float
    valeur: complex  # valeur propre a la traversee
    mode: int  # indice du mode suivi (colonne de Continuation.lam)
    elements: list[tuple[str, float]]  # participations cumulees par element [%]
    sens: int  # +1 : devient instable quand mu croit ; -1 : se stabilise
    #: un des deux points qui encadrent la traversee (au pas de dichotomie
    #: pres) a un convertisseur au-dela de sa limite : la bifurcation est
    #: celle du modele libre, que le materiel installe n'atteint pas tel quel
    hors_limite: bool = False

    @property
    def f(self) -> float:
        return abs(self.valeur.imag) / (2 * np.pi)

    def __str__(self) -> str:
        quoi = "destabilise" if self.sens > 0 else "stabilise"
        el = ", ".join(f"{n} {p:.0f} %" for n, p in self.elements[:3])
        hors = " [hors limite des convertisseurs]" if self.hors_limite else ""
        if self.genre == "hopf":
            return f"Hopf a mu = {self.mu:.6g} ({self.f:.3g} Hz, {quoi}) : {el}{hors}"
        return f"noeud-col a mu = {self.mu:.6g} ({quoi}) : {el}{hors}"


@dataclass
class Continuation:
    mu: np.ndarray  # valeurs du parametre effectivement calculees
    lam: np.ndarray  # (n_mu, n_modes), modes suivis ; NaN si absent
    etats: list[str]
    tensions: dict[str, np.ndarray] = field(default_factory=dict)  # |V| par bus
    bifurcations: list[Bifurcation] = field(default_factory=list)
    #: le ``mu`` qui n'a pas de point d'equilibre et a arrete la suite (le
    #: dernier point calcule est ``mu[-1]``) ; None si la suite est allee au bout
    arret_a: float | None = None
    #: les convertisseurs au-dela d'une limite (:class:`g2egrid.noyau.limites.Depassement`),
    #: une liste par valeur de ``mu`` (vide : tout est dans les limites)
    depassements: list[list[limites.Depassement]] = field(default_factory=list)
    #: le premier ``mu`` de la suite ou un convertisseur sort de ses limites,
    #: raffine par dichotomie (``mu[0]`` si le point de depart y est deja) ;
    #: None si aucun point n'est hors limite
    mu_hors_limite: float | None = None

    @property
    def hors_limite(self) -> np.ndarray:
        """Masque des ``mu`` dont le point a un convertisseur hors limite."""
        return np.array([bool(d) for d in self.depassements], dtype=bool)

    @property
    def marge(self) -> float | None:
        """Le premier ``mu`` qui destabilise, s'il y en a un (sans egard aux
        limites des convertisseurs : voir :attr:`mu_hors_limite`)."""
        d = [b.mu for b in self.bifurcations if b.sens > 0]
        return min(d) if d else None

    def resume(self) -> str:
        lignes = [str(b) for b in self.bifurcations] or ["aucune bifurcation"]
        if self.arret_a is not None:
            lignes.append(
                f"suite arretee : pas de point d'equilibre a mu = {self.arret_a:.6g} "
                f"(dernier point calcule : mu = {self.mu[-1]:.6g})"
            )
        if self.mu_hors_limite is not None:
            premier = next(d for d in self.depassements if d)
            qui = ", ".join(f"{d.element} ({d.grandeur})" for d in premier)
            lignes.append(
                f"hors limite des convertisseurs a partir de mu = "
                f"{self.mu_hors_limite:.6g} ({int(self.hors_limite.sum())} point(s) "
                f"sur {len(self.depassements)}) : {qui}"
            )
        return "\n".join(lignes)


# ------------------------------------------------------------- fabriques
def par_parametre(net, element: str, nom: str, equilibre: bool = False):
    """``mu -> net`` ou ``element.nom = mu``. ``equilibre`` : refaire load flow
    et equilibre (parametre qui deplace le point de fonctionnement) ; sinon
    seule la dynamique change (gain de regulation)."""
    el = next((e for e in net.all_elements if e.name == element), None)
    if el is None:
        raise ValueError(f"element {element!r} introuvable")
    if nom not in el.params:
        raise KeyError(f"{element} : parametre inconnu {nom!r}")
    if any(e.x0 is None for e in net.all_elements):
        net.solve_equilibrium()

    def fabrique(mu: float):
        el.set_param(nom, mu)
        if equilibre:
            net.run_loadflow()
            net.solve_equilibrium()
        return net

    return fabrique


def stations_statisme(t) -> dict[int, float]:
    """Les stations en statisme Vdc-P de ``t`` : bus alternatif -> ``kdroop``."""
    return {
        int(r[1]): float((t.PARAMS_PAR_BUS.get(int(r[1])) or {}).get("kdroop") or 0.0)
        for r in t.T_VSCs
        if str(r[3]) == "droop"
    }


def par_statisme(t, stations=None):
    """``mu -> reseau`` ou le statisme ``kdroop`` des stations en ``droop``
    (toutes, ou les bus alternatifs ``stations``) vaut ``mu`` fois celui de
    ``t`` : ``mu < 1`` raidit les stations, ``mu > 1`` les assouplit.

    La fabrique garde, pour chaque ``mu`` calcule, ce que le load flow dit des
    stations dans ``fabrique.mesures[mu]`` : la puissance alternative de
    chaque convertisseur a port continu (MW) et la tension de chaque noeud
    continu (pu) -- le partage de puissance que le statisme decide."""
    from g2egrid.reseau.tables import assembler

    k0 = stations_statisme(t)
    if not k0:
        raise ValueError("aucune station en statisme (mode droop) dans T_VSCs")
    choisies = set(k0) if stations is None else {int(b) for b in stations}
    inconnues = choisies - set(k0)
    if inconnues:
        raise ValueError(f"pas de station en statisme aux bus {sorted(inconnues)}")

    def fabrique(mu: float):
        par_bus = {int(b): dict(v) for b, v in t.PARAMS_PAR_BUS.items()}
        for b in choisies:
            par_bus.setdefault(b, {})["kdroop"] = k0[b] * float(mu)
        net = assembler(t.copie(PARAMS_PAR_BUS=par_bus))
        net.solve_equilibrium()
        idx = net._bus_index()
        prog = net._injections_programmees()
        fabrique.mesures[float(mu)] = {
            "p_mw": {
                int(el.buses[0]): float(
                    (net.lf.S[idx[el.buses[0]]].real - prog[el.buses[0]][0])
                    * net.Sb
                    / 1e6
                )
                for el in net._elements_vsc()
            },
            "vdc": {
                int(b): float(v) for b, v in zip(sorted(net.bus_specs_dc), net.lf.Vdc)
            },
        }
        return net

    fabrique.mesures = {}
    fabrique.stations = sorted(choisies)
    return fabrique


# ------------------------------------------------------------- spectre
@dataclass
class _Point:
    mu: float
    lam: np.ndarray
    p: np.ndarray
    etats: list[str]
    tensions: dict[str, float]
    depassements: list = field(default_factory=list)


def _construire(fabrique, mu):
    """``fabrique(mu)``, l'avertissement de depassement tu : la suite en
    emet un seul, a la fin (voir :func:`continuer`)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", limites.DepassementLimite)
        return fabrique(mu)


def _depassements(net) -> list:
    """RECALCULES, pas lus sur ``net.depassements`` : une fabrique qui ne
    refait pas l'equilibre (:func:`par_parametre` sur un gain, ou sur
    ``Imax`` lui-meme) laisse l'attribut a sa valeur d'avant."""
    if not hasattr(net, "all_elements"):
        return []
    return limites.depassements(net)


def _point(fabrique, mu) -> _Point | None:
    try:
        net = _construire(fabrique, mu)
        ss = net.build_state_space()
    except (RuntimeError, np.linalg.LinAlgError, ValueError):
        return None
    lam, _, p = modal.spectre(ss)
    tensions = {}
    if getattr(net, "lf", None) is not None:
        for b, i in net._bus_index().items():
            tensions[f"bus{b}"] = float(abs(net.lf.V[i]))
    return _Point(float(mu), lam, p, list(ss.states), tensions, _depassements(net))


def _hors_limite(fabrique, mu) -> bool | None:
    """Le point ``mu`` a-t-il un convertisseur hors limite ? ``None`` s'il n'a
    pas de point d'equilibre. Sans spectre : load flow et equilibre seuls."""
    try:
        net = _construire(fabrique, mu)
    except (RuntimeError, np.linalg.LinAlgError, ValueError):
        return None
    return bool(_depassements(net))


def _raffiner_hors_limite(fabrique, dedans: float, dehors: float, tol, iter_max):
    """Dichotomie sur ``[dedans, dehors]`` : le premier ``mu`` hors limite."""
    for _ in range(iter_max):
        if abs(dehors - dedans) <= tol * max(abs(dedans), 1.0):
            break
        milieu = 0.5 * (dedans + dehors)
        h = _hors_limite(fabrique, milieu)
        if h is None:
            break
        if h:
            dehors = milieu
        else:
            dedans = milieu
    return dehors


def _apparier(a: _Point, b: _Point) -> np.ndarray:
    """``j = perm[i]`` : le mode ``j`` de ``b`` continue le mode ``i`` de
    ``a`` (memes etats : la suite ne change pas la structure du modele)."""
    if a.etats != b.etats:
        raise ValueError("la structure du modele change le long de la continuation")
    ressemblance = np.sqrt(a.p) @ np.sqrt(b.p).T
    echelle = np.maximum(np.abs(a.lam[:, None]), 1.0)
    distance = np.abs(a.lam[:, None] - b.lam[None, :]) / echelle
    cout = (1 - ressemblance) + np.minimum(distance, 10.0)
    i, j = linear_sum_assignment(cout)
    perm = np.empty(len(a.lam), dtype=int)
    perm[i] = j
    return perm


def _signe(lam: complex, echelle: float) -> int:
    return 1 if lam.real > 1e-9 * echelle else -1


def _elements(p: np.ndarray, etats: list[str]) -> list[tuple[str, float]]:
    s: dict[str, float] = {}
    for v, e in zip(p, etats):
        k = e.split(".", 1)[0]
        s[k] = s.get(k, 0.0) + 100 * float(v)
    return sorted(s.items(), key=lambda t: -t[1])


def _raffiner(fabrique, a: _Point, b: _Point, i: int, tol: float, iter_max: int):
    """Dichotomie sur ``[a.mu, b.mu]`` pour le mode ``i`` de ``a`` dont la
    partie reelle change de signe. Rend ``(mu, lam, p, hors_limite)``, le
    dernier vrai si l'une des deux bornes finales a un convertisseur hors
    limite."""
    echelle = max(float(np.abs(a.lam).max()), 1.0)
    s0 = _signe(a.lam[i], echelle)
    cour = a
    fin = b
    for _ in range(iter_max):
        if abs(fin.mu - cour.mu) <= tol * max(abs(cour.mu), 1.0):
            break
        m = _point(fabrique, 0.5 * (cour.mu + fin.mu))
        if m is None:
            break
        # reindexer m dans l'ordre de cour, pour garder l'indice i
        perm = _apparier(cour, m)
        m = replace(m, lam=m.lam[perm], p=m.p[perm])
        if _signe(m.lam[i], echelle) == s0:
            cour = m
        else:
            fin = m
    # la traversee, par interpolation lineaire de Re lam entre les bornes
    perm = _apparier(cour, fin)
    l1, l2 = cour.lam[i], fin.lam[perm[i]]
    t = 0.0 if l2.real == l1.real else -l1.real / (l2.real - l1.real)
    t = float(np.clip(t, 0.0, 1.0))
    mu = cour.mu + t * (fin.mu - cour.mu)
    lam = complex(l1 + t * (l2 - l1))
    return mu, lam, cour.p[i], bool(cour.depassements or fin.depassements)


def continuer(
    fabrique: Callable[[float], object],
    valeurs,
    tol: float = 1e-4,
    iter_max: int = 40,
    raffiner: bool = True,
) -> Continuation:
    """Suit le spectre de ``fabrique(mu)`` pour ``mu`` dans ``valeurs``
    (monotones), repere et raffine les bifurcations, et marque les points ou
    un convertisseur sort de ses limites (voir l'en-tete du module). La suite
    s'arrete au dernier point converge si un pas n'a plus de point
    d'equilibre. ``raffiner`` regle le raffinement du premier ``mu`` hors
    limite des convertisseurs."""
    valeurs = [float(v) for v in valeurs]
    if len(valeurs) < 2:
        raise ValueError("il faut au moins deux valeurs de mu")
    pts: list[_Point] = []
    arret_a = None
    mu_hors = None
    bifs: list[Bifurcation] = []
    for n, mu in enumerate(valeurs):
        signaler(n / len(valeurs), f"mu = {mu:.6g}")
        pt = _point(fabrique, mu)
        if pt is None:
            if not pts:
                raise RuntimeError(f"aucun point d'equilibre a mu = {mu}")
            arret_a = mu
            break
        if pt.depassements and mu_hors is None:
            # le premier point hors limite ; entre lui et le precedent, la
            # frontiere, raffinee par dichotomie
            mu_hors = pt.mu
            if pts and raffiner:
                mu_hors = _raffiner_hors_limite(
                    fabrique, pts[-1].mu, pt.mu, tol, iter_max
                )
        if pts:
            perm = _apparier(pts[-1], pt)
            pt = replace(pt, lam=pt.lam[perm], p=pt.p[perm])
            a = pts[-1]
            echelle = max(float(np.abs(a.lam).max()), 1.0)
            for i in range(len(a.lam)):
                la, lb = a.lam[i], pt.lam[i]
                if min(abs(la), abs(lb)) < JAUGE * echelle:
                    continue
                if la.imag < -1e-9 * echelle:
                    continue  # une paire n'est comptee qu'une fois
                sa, sb = _signe(la, echelle), _signe(lb, echelle)
                if sa == sb:
                    continue
                mu_c, lam_c, p_c, hors = _raffiner(fabrique, a, pt, i, tol, iter_max)
                oscillant = abs(lam_c.imag) > 1e-6 * echelle
                sens = sb if pt.mu > a.mu else -sb
                bifs.append(
                    Bifurcation(
                        "hopf" if oscillant else "noeud-col",
                        mu_c,
                        lam_c,
                        i,
                        _elements(p_c, a.etats),
                        sens,
                        hors,
                    )
                )
        pts.append(pt)
    # les modes de chaque pas sont dans l'ordre du premier : un tableau suffit
    lam = np.array([p.lam for p in pts])
    noms = sorted({k for p in pts for k in p.tensions}, key=lambda s: int(s[3:]))
    tensions = {k: np.array([p.tensions.get(k, np.nan) for p in pts]) for k in noms}
    bifs.sort(key=lambda b: b.mu if valeurs[-1] > valeurs[0] else -b.mu)
    hors = [p for p in pts if p.depassements]
    if hors:
        limites.signaler(
            hors[0].depassements,
            quoi=(
                f"continuation, {len(hors)} point(s) sur {len(pts)} hors limite, "
                f"a partir de mu = {mu_hors:.6g}"
            ),
        )
    return Continuation(
        mu=np.array([p.mu for p in pts]),
        lam=lam,
        etats=pts[0].etats,
        tensions=tensions,
        bifurcations=bifs,
        arret_a=arret_a,
        depassements=[p.depassements for p in pts],
        mu_hors_limite=mu_hors,
    )
