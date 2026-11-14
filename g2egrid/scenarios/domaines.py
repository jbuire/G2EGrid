"""Le domaine d'un mode : alternatif, continu, ou couple AC/DC.

Dans un reseau a liaisons continues, la premiere question qu'on pose a un mode
n'est pas sa classe (angle, tension, convertisseur...) mais son COTE : vit-il
dans un reseau alternatif, dans le reseau continu, ou traverse-t-il les
convertisseurs ? Un mode couple AC/DC est celui ou une perturbation d'un cote
se lit de l'autre -- la raison d'etre d'une etude AC/DC unifiee.

LE DOMAINE D'UN ETAT se lit sur son element :

* noeuds continus, cables continus, DC/DC, charges continues : ``DC`` ;
* convertisseurs a port continu : ``AC`` pour le filtre, la PLL, les boucles
  de courant et de puissance ; ``DC`` pour la boucle de tension continue
  (``M_d`` du GFL a port continu) ;
* convertisseurs a bus continu interne (GFM3, GFM4, GFL2) : ``DC`` pour leur
  condensateur et leur courant continus (``vdc``, ``idc``, et ``M_d`` du
  GFL2) ;
* tout le reste -- machines, reseau alternatif, charges : ``AC``.

LE DOMAINE D'UN MODE vient de la participation de ses etats (somme 1) :
``AC`` si la part continue est sous ``seuil``, ``DC`` si la part alternative
l'est, ``AC/DC`` sinon. Chaque etat appartient aussi a une ILE -- composante
connexe du reseau alternatif ou du reseau continu -- et un mode dont deux iles
portent chacune au moins ``seuil`` est INTER-ILES : il traverse une liaison
continue (un parc offshore et la cote, deux reseaux asynchrones).

::

    from g2egrid import catalogue, domaines

    net = catalogue.construire("cigre_b4")
    for m in domaines.classer_modes(net)[:5]:
        print(m.texte())
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from g2egrid.analyse import modal

#: les deux cotes
DOMAINES = ("AC", "DC")
#: part en dessous de laquelle un cote (ou une ile) ne compte pas dans un mode
SEUIL = 0.1
#: un mode de jauge (zero exact) n'a pas de domaine
JAUGE = 1e-6

#: etats DC des convertisseurs, par classe d'element
_ETATS_DC_CONVERTISSEUR = {
    "GFL2DC": ("M_d",),
    # le MMC : son cote continu, l'energie des bras et la boucle de tension
    "MMC": ("idc", "W", "M_d", "vdc_f"),
    "MMCEnergie": ("idc", "W", "M_d", "vdc_f", "M_z", "M_W"),
    "GFL2": ("vdc", "idc", "M_d"),
    "GFM3": ("vdc", "idc"),
    "GFM4": ("vdc", "idc"),
    "GFMMatching": ("vdc", "idc"),
}


def _classe(el) -> str:
    return getattr(type(el), "_classe_complete", type(el)).__name__


def domaine_etat(el, etat: str) -> str:
    """``"AC"`` ou ``"DC"`` pour l'etat ``etat`` de l'element ``el``."""
    from g2egrid.elements import BusDC

    if isinstance(el, BusDC) or (el.ports_dc and not el.ports):
        return "DC"
    if etat in _ETATS_DC_CONVERTISSEUR.get(_classe(el), ()):
        return "DC"
    return "AC"


def _composantes(noeuds, aretes) -> dict[int, int]:
    """Noeud -> plus petit noeud de sa composante connexe."""
    parent = {n: n for n in noeuds}

    def racine(n):
        while parent[n] != n:
            parent[n] = parent[parent[n]]
            n = parent[n]
        return n

    for a, b in aretes:
        if a in parent and b in parent:
            ra, rb = racine(a), racine(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    return {n: racine(n) for n in noeuds}


def iles(net) -> dict[str, dict[int, int]]:
    """Les iles du reseau : ``{"AC": {bus: representant}, "DC": {...}}``.

    Une ile alternative est une composante connexe des bus alternatifs par les
    branches alternatives (lignes, transformateurs) ; une ile continue, celle
    des noeuds continus par les cables. Le representant est le plus petit
    numero de bus de l'ile."""
    ac = set(net.bus_specs)
    dc = set(net.bus_specs_dc)
    aretes_ac, aretes_dc = [], []
    for el in net.elements:
        if len(el.ports) == 2 and not el.ports_dc:
            aretes_ac.append((int(el.buses[0]), int(el.buses[1])))
        elif len(el.ports_dc) == 2 and not el.ports:
            aretes_dc.append((int(el.buses[0]), int(el.buses[1])))
    return {"AC": _composantes(ac, aretes_ac), "DC": _composantes(dc, aretes_dc)}


def etiquettes_etats(net, etats: list[str]) -> tuple[list[str], list[str]]:
    """Pour chaque etat ``"ELEMENT.x"`` : son domaine et le nom de son ile
    (``"AC 1"``, ``"DC 201"``, du plus petit bus de l'ile)."""
    par_nom = {el.name: el for el in net.all_elements}
    ile = iles(net)
    doms, noms = [], []
    for e in etats:
        nom, _, x = e.partition(".")
        el = par_nom.get(nom)
        if el is None:
            doms.append("AC")
            noms.append("?")
            continue
        d = domaine_etat(el, x)
        doms.append(d)
        if d == "DC":
            b = next((int(b) for b in el.buses if int(b) in ile["DC"]), None)
        else:
            b = next((int(b) for b in el.buses if int(b) in ile["AC"]), None)
        noms.append(f"{d} {ile[d][b]}" if b is not None else d)
    return doms, noms


@dataclass
class ModeDomaine:
    """Un mode, vu par cotes et par iles."""

    indice: int
    lam: complex
    part_ac: float
    part_dc: float
    domaine: str  # "AC", "DC" ou "AC/DC"
    iles: list[tuple[str, float]] = field(default_factory=list)

    @property
    def f(self) -> float:
        return abs(self.lam.imag) / (2 * np.pi)

    @property
    def zeta(self) -> float:
        return 100.0 * -self.lam.real / max(abs(self.lam), 1e-30)

    @property
    def inter_iles(self) -> bool:
        """Deux iles au moins portent chacune ``SEUIL`` du mode."""
        return sum(p >= SEUIL for _, p in self.iles) >= 2

    @property
    def couple(self) -> bool:
        """Couple AC/DC, ou alternatif d'une ile a l'autre par le continu."""
        return self.domaine == "AC/DC" or self.inter_iles

    def texte(self) -> str:
        iles = ", ".join(f"{n} {100 * p:.0f} %" for n, p in self.iles[:3])
        return (
            f"{self.lam.real:9.3f} {self.f:8.2f} Hz  zeta {self.zeta:5.1f} %  "
            f"{self.domaine:5s} (DC {100 * self.part_dc:.0f} %)  {iles}"
        )


def classer_modes(net, ss=None, seuil: float = SEUIL, oscillants: bool = False):
    """Les modes du reseau ``net`` (equilibre resolu), du moins amorti au plus
    amorti, avec leur domaine et leurs iles. Les modes de jauge sont ecartes ;
    ``oscillants`` ne garde qu'un mode par paire conjuguee (Im > 0)."""
    if ss is None:
        ss = net.build_state_space()
    lam, _, p = modal.spectre(ss)
    doms, noms_iles = etiquettes_etats(net, list(ss.states))
    est_dc = np.array([d == "DC" for d in doms])
    ile_index: dict[str, list[int]] = {}
    for j, n in enumerate(noms_iles):
        ile_index.setdefault(n, []).append(j)
    out = []
    for k, v in enumerate(lam):
        if abs(v) <= JAUGE or (oscillants and v.imag < 0):
            continue
        part_dc = float(p[k][est_dc].sum())
        part_ac = 1.0 - part_dc
        if part_dc < seuil:
            dom = "AC"
        elif part_ac < seuil:
            dom = "DC"
        else:
            dom = "AC/DC"
        parts = sorted(
            ((n, float(p[k][idx].sum())) for n, idx in ile_index.items()),
            key=lambda x: -x[1],
        )
        out.append(ModeDomaine(k, complex(v), part_ac, part_dc, dom, parts))
    out.sort(key=lambda m: -m.lam.real)
    return out


def bilan(modes: list[ModeDomaine]) -> dict[str, int]:
    """Combien de modes de chaque domaine, et combien sont inter-iles."""
    out = {d: 0 for d in ("AC", "DC", "AC/DC")}
    for m in modes:
        out[m.domaine] += 1
    out["inter-iles"] = sum(m.inter_iles for m in modes)
    return out


#: le continu rendu RIGIDE : tensions des noeuds figees, cables instantanes.
#: Comparer a ce modele dit ce que la dynamique du reseau continu fait aux
#: modes alternatifs (:func:`impact_continu`).
CONTINU_RIGIDE = {"BusDC": {"v": "fige"}, "LineDC": {"i": "statique"}}


def impact_continu(t, f_max: float = 100.0):
    """Ce que la dynamique du reseau continu fait aux modes : comparaison du
    reseau ``t`` (tables) a son double ou le continu est rigide
    (:data:`CONTINU_RIGIDE`), par :func:`g2egrid.scenarios.simplification.comparer`.
    Un mode alternatif RETROUVE mais deplace est un mode que le continu
    influence ; un mode RETIRE est un mode du continu lui-meme."""
    from g2egrid.scenarios.simplification import comparer

    if not t.T_BusDC:
        raise ValueError("pas de reseau continu (T_BusDC vide)")
    simpl = {k: dict(v) for k, v in CONTINU_RIGIDE.items()}
    if not t.T_LinesDC:
        simpl.pop("LineDC")
    return comparer(t, simpl)
