"""Transport et distribution : faut-il modeliser l'un pour etudier l'autre ?

LA QUESTION -- les etudes de stabilite se font d'ordinaire par morceaux.
L'exploitant de transport voit la distribution comme une CHARGE ; celui de
distribution voit le transport comme un RESEAU INFINI derriere une
IMPEDANCE. Les deux simplifications etaient sures tant que la distribution
etait passive. Avec des convertisseurs grid-following (GFL) dans la
distribution, le sont-elles encore ? Ce module construit de quoi y repondre :

* **des reseaux** : WSCC9 dont chaque charge est remplacee par deux reseaux
  CIGRE MV a GFL (:func:`wscc9_cigre6`), et IEEE14 dont quatre charges
  sont remplacees par des departs IEEE33 a GFL (:func:`ieee14_ieee33`).
  Chaque distribution est une ZONE a elle, pour qu'on puisse y ajouter des
  GFL independamment ;
* **des GFL ajoutes au hasard** dans la distribution
  (:func:`ajouter_gfl`), un par un, pour faire croitre la penetration ;
* **les deux modeles reduits** :

  - la distribution etudiee SEULE (:func:`distribution_seule`) : tout le
    reste du reseau devient un bus infini derriere son impedance de
    Thevenin, vue du poste ;
  - le transport etudie seul : chaque distribution devient la charge
    statique qu'elle soutire (:func:`g2egrid.multi_niveau.impact.equivalent_statique`) ;

* **les comparaisons** au modele complet (:func:`comparer_distribution`,
  :func:`comparer_transport`) et le **balayage de penetration**
  (:func:`penetration`).

Les modes sont apparies par la ressemblance de leurs participations sur
les etats communs, comme dans :mod:`g2egrid.scenarios.robustesse` : les deux modeles
compares ont la MEME numerotation (on decoupe le reseau aplati), donc les
memes noms d'etats.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from g2egrid.analyse import modal
from g2egrid.reseau.catalogue import charger
from g2egrid.elements.bus import NOEUDS_DEFAUT, R_CAUSALITE
from g2egrid.multi_niveau.impact import bus_d_instance, equivalent_statique
from g2egrid.scenarios.robustesse import (
    SANS_POINT_DE_FONCTIONNEMENT,
    _apparier,
    _elements,
    _Spectre,
    _spectre,
)
from g2egrid.reseau.tables import TABLES, Tables, assembler
from g2egrid.reseau.zones import FAMILLES, Instance, Reseau, TransfoPoste, Zone, aplatir

__all__ = [
    "ComparaisonModes",
    "ModeCompare",
    "Partition",
    "ajouter_gfl",
    "ajouter_unites",
    "comparer_partition",
    "partition",
    "regler_gfm",
    "remplacer_machines",
    "comparer_distribution",
    "comparer_transport",
    "poste_faible",
    "regler_gfl",
    "equivalent_hybride",
    "evaluer",
    "modes_interaction",
    "distribution_seule",
    "ieee14_ieee33",
    "ieee39_mixte",
    "impedance_thevenin",
    "penetration",
    "wscc9_cigre6",
    "wscc9_cigre9",
]

SLACK, PV, PQ = 1, 2, 3
JAUGE = 1e-6
#: au-dela, les modes sont les resonances electromagnetiques des cables et
#: des filtres : hors du propos (et identiques dans tous les modeles)
F_MAX = 100.0
#: au-dela, un mode ne boucle pas une periode avant de s'eteindre : il est
#: aperiodique, et son ecart n'a pas de sens oscillatoire
ZETA_MAX = 95.0
#: la part de la charge d'un bus du transport que ses distributions reprennent
PART_REMPLACEE = 0.98


# ---------------------------------------------------------------- les zones
def zone_depuis(
    t: Tables,
    bus_pcc: int,
    f0: float,
    Sb: float,
    echelle_charge: float = 1.0,
    nom: str = "zone",
    echelle_q: float | None = None,
    graine: int | None = None,
    dispersion: float = 0.0,
) -> Zone:
    """Une ZONE faite d'un reseau de distribution du catalogue.

    Le bus infini est retire (le transport le remplace), le slack devient un
    bus PQ, les charges sont mises a l'echelle, et la base Sb comme la
    frequence sont celles du transport (``echelle_q`` : une autre echelle
    pour le reactif, sinon celle de l'actif). ``dispersion`` tire, bus par
    bus, un facteur dans ``[1 - d, 1 + d]`` (graine ``graine``), puis
    renormalise : le total reste celui demande, mais deux zones faites du
    meme reseau ne sont plus identiques -- deux copies exactes au meme poste
    donneraient des paires de modes symetriques, partagees a 50/50. Les GFL restent en PQ (ce sont des
    PV : ils suivent, a leur puissance MPPT).
    """
    eq = echelle_charge if echelle_q is None else echelle_q
    z = t.copie(f0=f0, Sb=Sb, nom=nom)
    k = t.Sb / Sb
    z.T_Shunts = [(a, b, float(c) * k) for a, b, c in t.T_Shunts]
    z.T_ShuntsDC = [(a, b, float(c) * k) for a, b, c in t.T_ShuntsDC]
    z.T_IBs = []
    rng = np.random.default_rng(graine)
    fac = (
        1 + dispersion * rng.uniform(-1, 1, len(t.Y_DER))
        if dispersion
        else np.ones(len(t.Y_DER))
    )
    pl = np.array([float(r[7]) for r in t.Y_DER])
    ql = np.array([float(r[8]) for r in t.Y_DER])
    fp = fac * (pl.sum() / max((fac * pl).sum(), 1e-30))
    fq = fac * (ql.sum() / max((fac * ql).sum(), 1e-30))
    y = []
    for i, r in enumerate(t.Y_DER):
        r = list(r)
        if int(r[1]) == SLACK:
            r[1] = PQ
            r[5] = r[6] = 0.0
        r[7] = float(r[7]) * echelle_charge * fp[i]
        r[8] = float(r[8]) * eq * fq[i]
        y.append(tuple(r))
    z.Y_DER = y
    z.T_Loads = [
        (a, b, float(p) * echelle_charge, float(q) * eq, ty)
        for a, b, p, q, ty in t.T_Loads
    ]
    return Zone(tables=z, bus_pcc=int(bus_pcc), nom=nom)


def _raccorder(
    transport: Tables,
    dist: Tables,
    bus_pcc: int,
    departs: list[tuple[int, int]],
    nom: str,
    transfo: tuple[float, float, float] = (0.003, 0.12, 50e6),
    prise: float = 1.0,
    dispersion: float = 0.2,
    premier: int = 1,
) -> tuple[Tables, dict[str, Zone], list[Instance]]:
    """Des distributions a la place des charges du transport.

    ``departs`` : ``(bus_amont, K)``, deux zones par bus si le bus y figure
    deux fois. ``prise`` : rapport de la prise du transfo de poste (1,05
    releve la tension de la distribution de 5 %, comme le regleur en charge
    d'un poste source). Chaque zone est une copie de ``dist`` dont les charges sont
    mises a l'echelle pour que ``K`` departs remplacent exactement leur part
    de la charge du bus : ``charge_residuelle`` retire ensuite du bus ce que
    la distribution consomme, et il n'y reste (presque) rien.
    """
    p_dist, q_dist = dist.charge_totale()
    un = {int(r[0]): float(r[2]) for r in transport.Y_DER}
    charge = {int(r[0]): float(r[7]) for r in transport.Y_DER}
    reactif = {int(r[0]): float(r[8]) for r in transport.Y_DER}
    parts: dict[int, int] = {}
    for b, _k in departs:
        parts[b] = parts.get(b, 0) + 1
    zones, instances = {}, []
    un_pcc = {int(r[0]): float(r[2]) for r in dist.Y_DER}[bus_pcc]
    R, X, Sn = transfo
    for i, (b, K) in enumerate(departs):
        # ce que ce groupe de K departs remplace, en actif et en reactif : 98 %
        # de la charge du bus, le reste etant la consommation propre du poste
        # (une charge residuelle exactement nulle tomberait a -1e-8 var par
        # arrondi, que l'assembleur refuse)
        echelle = PART_REMPLACEE * charge[b] / parts[b] / (K * p_dist)
        # Le reactif est plafonne a l'echelle de l'actif : une distribution
        # garde son propre facteur de puissance. Charger un CIGRE (Q/P = 0,06)
        # du Q/P de 0,4 d'un bus de transport lui imposait sept fois son
        # reactif nominal, et des tensions de 0,79 pu ; le reactif en trop
        # reste au bus de transport (compensation, charges HT).
        echelle_q = min(PART_REMPLACEE * reactif[b] / parts[b] / (K * q_dist), echelle)
        zn = f"{nom}{i + premier}"
        zones[zn] = zone_depuis(
            dist,
            bus_pcc,
            transport.f0,
            transport.Sb,
            echelle,
            zn,
            echelle_q,
            graine=i + premier,
            dispersion=dispersion,
        )
        instances.append(
            Instance(
                zone=zn,
                prefixe=100 * (i + premier),
                bus_amont=b,
                K=K,
                transfo=TransfoPoste(Sn, un[b], un_pcc * prise, R, X),
                charge_residuelle=True,
            )
        )
    return transport, zones, instances


def _wscc9_cigre(par_bus: int, K: int, prise: float, nom: str) -> Reseau:
    """WSCC9 dont chaque charge (bus 5, 6, 8) est remplacee par ``par_bus``
    reseaux CIGRE MV a GFL, chacun une zone a lui, de ``K`` departs."""
    # wscc9 garde la resistance de Node.m, pour la parite MATLAB ; le reseau
    # multi-echelle est un reseau d'etude, aux noeuds par defaut (FR-005)
    tr = charger("wscc9").copie(NOEUDS=NOEUDS_DEFAUT, R_NOEUDS=R_CAUSALITE)
    cigre = charger("cigre_mv")
    departs = [(b, K) for b in (5, 6, 8) for _ in range(par_bus)]
    tr, zones, inst = _raccorder(
        tr, cigre, 0, departs, "cigre", (0.003, 0.12, 50e6), prise
    )
    return Reseau(transport=tr, zones=zones, instances=inst, nom=nom)


def wscc9_cigre6(K: int = 3, prise: float = 1.05) -> Reseau:
    """WSCC9, dont chaque charge (bus 5, 6, 8) est remplacee par DEUX
    reseaux CIGRE MV a GFL -- six distributions, chacune une zone a elle
    (99 bus).

    Chaque instance agrege ``K`` departs identiques ; les charges CIGRE
    sont mises a l'echelle pour que les deux instances d'un bus y remplacent
    exactement sa charge.
    """
    return _wscc9_cigre(2, K, prise, "wscc9_cigre6")


def wscc9_cigre9(K: int = 2, prise: float = 1.05) -> Reseau:
    """WSCC9, dont chaque charge est remplacee par TROIS reseaux CIGRE MV :
    neuf distributions, 144 bus. ``K = 2`` garde six departs par bus, comme
    :func:`wscc9_cigre6`."""
    return _wscc9_cigre(3, K, prise, "wscc9_cigre9")


def ieee14_ieee33(prise: float = 1.08) -> Reseau:
    """IEEE14, dont les charges des bus 9, 10, 13 et 14 sont remplacees par
    des departs IEEE33 -- quatre distributions radiales, chacune avec trois
    GFL (PV) en bout de ligne, aux bus 18, 25 et 33."""
    tr = charger("ieee14")
    d33 = charger("ieee33")
    p = d33.charge_totale()[0]
    # trois PV en bout d'antenne, 10 % de la charge du depart chacun
    gfl = []
    y = []
    for r in d33.Y_DER:
        r = list(r)
        if int(r[0]) in (18, 25, 33):
            gfl.append((len(gfl) + 1, int(r[0]), 0.15 * p, 40e3))
            r[5] = 0.1 * p
        y.append(tuple(r))
    d33 = d33.copie(Y_DER=y, T_GFL1s=gfl)
    departs = [(9, 8), (10, 3), (13, 4), (14, 4)]
    tr, zones, inst = _raccorder(
        tr, d33, 1, departs, "ieee33_", (0.003, 0.08, 10e6), prise
    )
    return Reseau(transport=tr, zones=zones, instances=inst, nom="ieee14_ieee33")


def ieee39_mixte(seuil: float = 100e6) -> Reseau:
    """Le grand exemple : IEEE39 (New England, 10 machines), dont chaque
    charge d'au moins ``seuil`` (17 bus) est remplacee par une distribution,
    alternativement un CIGRE MV (trois GFL) et un depart IEEE33 (trois GFL en
    bout d'antenne). ``K`` est l'entier qui evite de charger une
    distribution au-dela de sa charge nominale."""
    import math

    tr = charger("ieee39")
    cigre = charger("cigre_mv")
    d33 = (
        ieee14_ieee33().zones["ieee33_1"].tables
    )  # l'IEEE33 a GFL, non mis a l'echelle
    d33 = charger("ieee33").copie(
        T_GFL1s=d33.T_GFL1s,
        Y_DER=[
            (*r[:5], float(g[5]), *r[6:])
            for r, g in zip(charger("ieee33").Y_DER, d33.Y_DER)
        ],
    )
    charges = [
        (int(r[0]), float(r[7]))
        for r in tr.Y_DER
        if int(r[1]) == PQ and float(r[7]) >= seuil
    ]
    p_c, p_i = cigre.charge_totale()[0], d33.charge_totale()[0]
    dep_c = [
        (b, math.ceil(PART_REMPLACEE * p / p_c))
        for k, (b, p) in enumerate(charges)
        if k % 2 == 0
    ]
    dep_i = [
        (b, math.ceil(PART_REMPLACEE * p / p_i))
        for k, (b, p) in enumerate(charges)
        if k % 2 == 1
    ]
    tr, zc, ic = _raccorder(tr, cigre, 0, dep_c, "cigre", (0.003, 0.12, 50e6), 1.05)
    tr, zi, ii = _raccorder(
        tr, d33, 1, dep_i, "ieee33_", (0.003, 0.08, 10e6), 1.08, premier=len(dep_c) + 1
    )
    return Reseau(
        transport=tr, zones={**zc, **zi}, instances=[*ic, *ii], nom="ieee39_mixte"
    )


# ------------------------------------------------------ des GFL au hasard
def _unites(t: Tables) -> set[int]:
    return {int(r[1]) for nom in (*FAMILLES, "T_IBs") for r in getattr(t, nom)}


def ajouter_unites(
    r: Reseau,
    n: int = 1,
    graine: int | None = None,
    part: float = 0.05,
    zones: list[str] | None = None,
    genre: str = "GFL",
) -> tuple[Reseau, list[tuple[str, int, float]]]:
    """Ajoute ``n`` convertisseurs a des bus tires au hasard dans la
    distribution : des GFL (PV) ou des GFM (batteries), ``genre``.

    Le tirage est uniforme sur les bus libres de toutes les zones (ou de
    ``zones``) : ni le PCC, ni un bus qui porte deja une unite. Chaque unite
    produit ``part`` de la charge de sa zone (5 %), et son bus garde son
    type PQ : le point de fonctionnement est impose, le GFM regle ensuite
    autour de lui. Le reseau rendu est une COPIE ; rend aussi ``(zone, bus
    local, P [W])`` par unite ajoutee.
    """
    if genre not in ("GFL", "GFM"):
        raise ValueError(f"genre {genre!r} : attendu 'GFL' ou 'GFM'")
    table = "T_GFL1s" if genre == "GFL" else "T_GFM1s"
    rng = np.random.default_rng(graine)
    out = _copie_zones(r)
    ajouts = []
    for _ in range(n):
        libres = [
            (zn, int(row[0]))
            for zn, z in out.zones.items()
            if zones is None or zn in zones
            for row in z.tables.Y_DER
            if int(row[0]) != z.bus_pcc and int(row[0]) not in _unites(z.tables)
        ]
        if not libres:
            break
        zn, b = libres[int(rng.integers(len(libres)))]
        t = out.zones[zn].tables
        p = part * t.charge_totale()[0]
        lignes = getattr(t, table)
        ident = max((int(x[0]) for x in lignes), default=0) + 1
        setattr(t, table, [*lignes, (ident, b, 1.2 * p, 40e3)])
        t.Y_DER = [
            (*x[:5], float(x[5]) + p, *x[6:]) if int(x[0]) == b else x for x in t.Y_DER
        ]
        ajouts.append((zn, b, p))
    return out, ajouts


def ajouter_gfl(
    r: Reseau,
    n: int = 1,
    graine: int | None = None,
    part: float = 0.05,
    zones: list[str] | None = None,
) -> tuple[Reseau, list[tuple[str, int, float]]]:
    """:func:`ajouter_unites` pour des GFL."""
    return ajouter_unites(r, n, graine, part, zones, "GFL")


def _copie_zones(r: Reseau, **tables_modifs) -> Reseau:
    zones = {}
    for k, z in r.zones.items():
        t = z.tables.copie()
        for champ, reglages in tables_modifs.items():
            setattr(t, champ, {**getattr(t, champ), **reglages})
        zones[k] = Zone(
            tables=t, bus_pcc=z.bus_pcc, positions=dict(z.positions), nom=z.nom
        )
    return Reseau(
        transport=r.transport,
        zones=zones,
        instances=list(r.instances),
        nom=r.nom,
        positions=dict(r.positions),
        f_zoom=r.f_zoom,
    )


def regler_gfl(r: Reseau, **reglages) -> Reseau:
    """Une copie ou les GFL de TOUTES les zones prennent ``reglages``
    (reglage de famille de chaque zone : ``Kppll=..., Kip=...``)."""
    return _copie_zones(r, GFL_PARAMS=reglages)


def regler_gfm(r: Reseau, **reglages) -> Reseau:
    """Une copie ou les GFM de TOUTES les zones prennent ``reglages``."""
    return _copie_zones(r, GFM_PARAMS=reglages)


def remplacer_machines(r: Reseau, bus: dict[int, str]) -> Reseau:
    """Une copie ou la machine synchrone de chaque bus du TRANSPORT de
    ``bus`` devient un convertisseur de meme puissance, ``"GFM"`` (batterie)
    ou ``"GFL"`` (PV) ; le bus reste PV."""
    t = r.transport.copie()
    for b, genre in bus.items():
        sm = next((x for x in t.T_SMs if int(x[1]) == b), None)
        if sm is None:
            raise ValueError(f"pas de machine synchrone au bus {b} du transport")
        t.T_SMs = [x for x in t.T_SMs if int(x[1]) != b]
        table = {"GFL": "T_GFL1s", "GFM": "T_GFM1s"}[genre]
        lignes = getattr(t, table)
        ident = max((int(x[0]) for x in lignes), default=0) + 1
        setattr(t, table, [*lignes, (ident, b, sm[2], 40e3)])
    return Reseau(
        transport=t,
        zones=r.zones,
        instances=list(r.instances),
        nom=r.nom,
        positions=dict(r.positions),
        f_zoom=r.f_zoom,
    )


def poste_faible(r: Reseau, X: float) -> Reseau:
    """Une copie ou chaque transfo de poste a la reactance ``X`` [pu, base
    du transfo] : un raccordement plus faible."""
    inst = [
        Instance(
            zone=i.zone,
            prefixe=i.prefixe,
            bus_amont=i.bus_amont,
            K=i.K,
            transfo=None
            if i.transfo is None
            else TransfoPoste(*tuple(i.transfo)[:4], X),
            charge_residuelle=i.charge_residuelle,
            shunt_pcc=i.shunt_pcc,
        )
        for i in r.instances
    ]
    return Reseau(
        transport=r.transport,
        zones=r.zones,
        instances=inst,
        nom=r.nom,
        positions=dict(r.positions),
        f_zoom=r.f_zoom,
    )


# -------------------------------------------- le transport vu du poste
def _admittances_unites(t: Tables, net, lf) -> dict[int, complex]:
    """Ce que chaque bus voit a la terre, pour un court-circuit : les
    charges en impedance constante (au point du load flow), les machines
    derriere leur reactance transitoire, les GFM derriere leur filtre, les
    bus infinis derriere leur impedance. Les GFL sont des sources de
    courant : ils n'ajoutent rien. En pu base systeme."""
    Sb = t.Sb
    idx = net._bus_index()
    y: dict[int, complex] = {}

    def ajoute(b, v):
        y[b] = y.get(b, 0j) + v

    for r in t.Y_DER:
        b = int(r[0])
        s = complex(float(r[7]), float(r[8])) / Sb
        v = abs(lf.V[idx[b]])
        if abs(s) > 0:
            ajoute(b, s.conjugate() / v**2)
    for r in t.T_Loads:
        b = int(r[1])
        s = complex(float(r[2]), float(r[3])) / Sb
        ajoute(b, s.conjugate() / abs(lf.V[idx[b]]) ** 2)
    from g2egrid.elements.gfm import GFM_DEFAULT
    from g2egrid.elements.sync_machine import SM_DEFAULT

    for r in t.T_SMs:
        b, Sn = int(r[1]), float(r[2])
        p = SM_DEFAULT | t.SM_PARAMS | t.PARAMS_PAR_BUS.get(b, {})
        xd1 = p["Ll"] + p["Lad"] * p["Lfd"] / (p["Lad"] + p["Lfd"])
        ajoute(b, (Sn / Sb) / complex(p["Ra"], xd1))
    for table in ("T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s"):
        for r in getattr(t, table):
            b, Sn = int(r[1]), float(r[2])
            p = GFM_DEFAULT | t.GFM_PARAMS | t.PARAMS_PAR_BUS.get(b, {})
            z = complex(p["R1_pu"] + p["R2_pu"], p["L1_pu"] + p["L2_pu"])
            ajoute(b, (Sn / Sb) / z)
    for r in t.T_IBs:
        b, Sn = int(r[1]), float(r[2])
        ajoute(b, (Sn / Sb) / complex(float(r[5]), float(r[6])))
    return y


def impedance_thevenin(t: Tables, bus: int) -> complex:
    """L'impedance de Thevenin du reseau ``t`` vue du bus ``bus`` [pu base
    systeme] : Y-bus des branches et des capacites, plus les admittances a
    la terre de :func:`_admittances_unites`."""
    net = assembler(t)
    lf = net.run_loadflow()
    Y = net.build_ybus().copy()
    idx = net._bus_index()
    for b, y in _admittances_unites(t, net, lf).items():
        Y[idx[b], idx[b]] += y
    Z = np.linalg.inv(Y)
    return complex(Z[idx[bus], idx[bus]])


def _garder(t: Tables, bus: set[int]) -> Tables:
    """Une copie des tables reduite aux lignes dont TOUS les bus sont dans
    ``bus``."""
    out = t.copie()
    for nom, cols in TABLES.items():
        ic = [
            j
            for j, c in enumerate(cols)
            if c in ("BusID", "From", "To", "BusAC", "BusDC")
        ]
        if nom in ("Y_DER", "T_BusDC"):
            ic = [0]
        if not ic:
            continue
        setattr(
            out, nom, [x for x in getattr(t, nom) if all(int(x[j]) in bus for j in ic)]
        )
    out.PARAMS_PAR_BUS = {b: v for b, v in t.PARAMS_PAR_BUS.items() if int(b) in bus}
    out.OUVERTS = [o for o in t.OUVERTS if o[0] != "Y_DER" or int(o[1]) in bus]
    return out


def distribution_seule(
    r: Reseau, k: int, zth: bool = True, plat: Tables | None = None
) -> tuple[Tables, complex]:
    """L'instance ``k`` etudiee SEULE : son poste et sa distribution, et au
    bus amont un bus infini derriere l'impedance de Thevenin de TOUT le
    reste du reseau (transport et autres distributions), vue de ce bus.

    Le bus amont devient le slack, a la tension et a l'angle du modele
    complet : la distribution voit le meme point de fonctionnement. Ses
    charges et unites de transport sont dans l'equivalent, pas dans le
    modele reduit. ``zth=False`` : bus infini raide (impedance nulle).
    Rend les tables et l'impedance employee [pu base systeme].
    """
    plat = plat if plat is not None else aplatir(r)
    inst = r.instances[k]
    amont = int(inst.bus_amont)
    bus = bus_d_instance(r, inst)
    # le reste du reseau, vu du bus amont
    reste = _garder(plat, {int(x[0]) for x in plat.Y_DER} - bus)
    z = impedance_thevenin(reste, amont) if zth else complex(0, 1e-4)
    net = assembler(plat)
    lf = net.run_loadflow()
    v = lf.V[net._bus_index()[amont]]

    red = _garder(plat, bus | {amont})
    # au bus amont : rien du transport, sinon la reference et l'equivalent
    for nom in (*FAMILLES, "T_IBs", "T_Loads", "T_Shunts"):
        setattr(red, nom, [x for x in getattr(red, nom) if int(x[1]) != amont])
    red.PARAMS_PAR_BUS.pop(amont, None)
    red.Y_DER = [
        (b, SLACK, un, abs(v), float(np.angle(v)), 0.0, 0.0, 0.0, 0.0)
        if int(b) == amont
        else (b, ty, un, V, th, pg, qg, pl, ql)
        for b, ty, un, V, th, pg, qg, pl, ql in red.Y_DER
    ]
    red.T_IBs = [
        (1, amont, red.Sb, abs(v), float(np.degrees(np.angle(v))), z.real, z.imag)
    ]
    return red, z


# ------------------------------------------------------------ comparaisons
@dataclass
class ModeCompare:
    lam_complet: complex
    lam_reduit: complex | None  # None : absent du modele reduit
    part: float  # participation de la partie etudiee au mode [0-1]
    elements: list[tuple[str, float]]

    @property
    def f(self) -> float:
        return float(modal.frequence(self.lam_complet)[0])

    @property
    def zeta(self) -> float:
        return float(modal.amortissement(self.lam_complet)[0])

    @property
    def f_reduit(self) -> float | None:
        return (
            None
            if self.lam_reduit is None
            else float(modal.frequence(self.lam_reduit)[0])
        )

    @property
    def zeta_reduit(self) -> float | None:
        return (
            None
            if self.lam_reduit is None
            else float(modal.amortissement(self.lam_reduit)[0])
        )

    def dans_bande(self, f_max: float = F_MAX) -> bool:
        """Oscillant (f > 0,05 Hz) et de frequence PROPRE abs(lambda)/2.pi au plus
        ``f_max`` : la frequence amortie seule laisserait passer les modes
        electromagnetiques tres rapides (abs(Re) ~ 1e4) qui portent, en repere dq,
        la rotation a f0 -- 60 Hz amortis a 100 %, hors du propos."""
        return self.f > 0.05 and abs(self.lam_complet) / (2 * np.pi) <= f_max

    @property
    def ecart(self) -> float:
        """abs(d lambda) / abs(lambda) ; 1 si le mode manque au modele reduit."""
        if self.lam_reduit is None:
            return 1.0
        return abs(self.lam_reduit - self.lam_complet) / abs(self.lam_complet)

    @property
    def dzeta(self) -> float | None:
        return None if self.lam_reduit is None else self.zeta_reduit - self.zeta


@dataclass
class ComparaisonModes:
    """Les modes de la partie etudiee, modele complet contre modele reduit."""

    etude: str  # "distribution" ou "transport"
    modes: list[ModeCompare] = field(default_factory=list)
    #: modes du modele COMPLET que la partie etudiee porte a moins de
    #: ``part_min`` mais que le modele reduit rate : ils n'existent que par
    #: l'interaction
    manques: int = 0

    def _candidats(self, f_max: float, zeta_max: float) -> list[ModeCompare]:
        return [m for m in self.modes if m.dans_bande(f_max) and m.zeta <= zeta_max]

    def pire(
        self, f_max: float = F_MAX, zeta_max: float = ZETA_MAX
    ) -> ModeCompare | None:
        """Le mode oscillant le plus deplace, parmi ceux que le modele reduit
        RETROUVE (les absents se comptent a part : :meth:`absents`)."""
        cand = [m for m in self._candidats(f_max, zeta_max) if m.lam_reduit is not None]
        return max(cand, key=lambda m: m.ecart, default=None)

    def absents(
        self, f_max: float = F_MAX, zeta_max: float = ZETA_MAX
    ) -> list[ModeCompare]:
        """Les modes que le modele reduit n'a pas : cote distribution, ce sont
        typiquement des modes PARTAGES entre deux departs du meme poste, dont
        un depart seul ne porte que la moitie."""
        return [m for m in self._candidats(f_max, zeta_max) if m.lam_reduit is None]

    def ecart_max(self, f_max: float = F_MAX, zeta_max: float = ZETA_MAX) -> float:
        m = self.pire(f_max, zeta_max)
        return m.ecart if m else 0.0

    def dzeta_max(self, f_max: float = F_MAX, zeta_max: float = ZETA_MAX) -> float:
        cand = [
            abs(m.dzeta)
            for m in self._candidats(f_max, zeta_max)
            if m.dzeta is not None
        ]
        return max(cand, default=0.0)


def _comparer(complet, reduit, etats_etudies: set[str], part_min: float, etude: str):
    """Apparie les modes de ``complet`` portes a plus de ``part_min`` par
    ``etats_etudies`` a ceux de ``reduit``."""
    part = _part(complet, etats_etudies)
    garde = np.where((part >= part_min) & (np.abs(complet.lam) > JAUGE))[0]
    paires = _apparier(complet, reduit, garde) if garde.size else {}
    c = ComparaisonModes(etude=etude)
    for i in garde:
        j = paires.get(int(i))
        c.modes.append(
            ModeCompare(
                lam_complet=complex(complet.lam[i]),
                lam_reduit=None if j is None else complex(reduit.lam[j]),
                part=float(part[i]),
                elements=_elements(complet.p[i], complet.etats),
            )
        )
    c.modes.sort(key=lambda m: -m.ecart)
    return c


def _part(complet, etats: set[str]) -> np.ndarray:
    idx = [i for i, e in enumerate(complet.etats) if e in etats]
    return complet.p[:, idx].sum(axis=1)


def _etats_de(etats: list[str], noms_elements: set[str]) -> set[str]:
    return {e for e in etats if e.split(".", 1)[0] in noms_elements}


def _noms_elements(t: Tables) -> set[str]:
    return {el.name for el in assembler(t).all_elements}


def comparer_distribution(
    r: Reseau,
    k: int,
    zth: bool = True,
    part_min: float = 0.5,
    plat: Tables | None = None,
    complet=None,
) -> ComparaisonModes:
    """Les modes de la distribution ``k`` : modele complet contre
    distribution seule derriere un bus infini (et l'impedance de Thevenin
    si ``zth``). Les modes compares sont ceux que la distribution porte a
    plus de ``part_min`` dans le modele complet."""
    plat = plat if plat is not None else aplatir(r)
    complet = complet if complet is not None else _spectre(plat)
    red, _ = distribution_seule(r, k, zth, plat)
    reduit = _spectre(red)
    noms = _noms_elements(red) - {"IB" + str(r.instances[k].bus_amont)}
    return _comparer(
        complet, reduit, _etats_de(complet.etats, noms), part_min, "distribution"
    )


def equivalent_hybride(r: Reseau, plat: Tables | None = None) -> Tables:
    """Comme :func:`g2egrid.multi_niveau.impact.equivalent_statique`, mais la production
    des GFL de chaque distribution y est une injection a PUISSANCE CONSTANTE
    (``T_Loads`` de type ``P``, puissance negative), le reste du soutirage
    restant en R-L. Un GFL suit sa consigne de puissance quelle que soit la
    tension : en R-L, l'equivalent statique lui prete un comportement
    d'impedance qu'il n'a pas. Le point de fonctionnement est le meme."""
    plat = plat if plat is not None else aplatir(r)
    eq, _ = equivalent_statique(r, plat=plat)
    p_gfl: dict[int, float] = {}
    for inst in r.instances:
        t = r.zones[inst.zone].tables
        gfl = {int(x[1]) for x in (*t.T_GFL1s, *t.T_GFL2s)}
        p = inst.K * sum(float(x[5]) for x in t.Y_DER if int(x[0]) in gfl)
        p_gfl[int(inst.bus_amont)] = p_gfl.get(int(inst.bus_amont), 0.0) + p
    eq.Y_DER = [
        (*x[:7], float(x[7]) + p_gfl.get(int(x[0]), 0.0), x[8]) for x in eq.Y_DER
    ]
    ident = max((int(x[0]) for x in eq.T_Loads), default=0)
    eq.T_Loads = [
        *eq.T_Loads,
        *(
            (ident + i + 1, b, -p, 0.0, "P")
            for i, (b, p) in enumerate(sorted(p_gfl.items()))
            if p
        ),
    ]
    return eq


def _sans_injections(t: Tables) -> Tables:
    """Les injections a puissance constante de l'equivalent hybride ne sont
    pas des elements du transport : leurs etats ne comptent pas comme tels."""
    return t.copie(T_Loads=[x for x in t.T_Loads if x[4] != "P" or float(x[2]) >= 0])


def comparer_transport(
    r: Reseau,
    part_min: float = 0.5,
    plat: Tables | None = None,
    complet=None,
    hybride: bool = False,
) -> ComparaisonModes:
    """Les modes du transport : modele complet contre transport ou chaque
    distribution est la charge statique qu'elle soutire -- en R-L, ou, avec
    ``hybride``, sa production GFL en puissance constante
    (:func:`equivalent_hybride`)."""
    plat = plat if plat is not None else aplatir(r)
    complet = complet if complet is not None else _spectre(plat)
    if hybride:
        eq = equivalent_hybride(r, plat=plat)
    else:
        eq, _ = equivalent_statique(r, plat=plat)
    reduit = _spectre(eq)
    noms = _noms_elements(_sans_injections(eq))
    return _comparer(
        complet, reduit, _etats_de(complet.etats, noms), part_min, "transport"
    )


# ------------------------------------------------ partition de l'etat x_T, x_D
@dataclass
class Partition:
    """Le modele d'etat du reseau COMPLET, ``dx/dt = A x + B u``, et le
    partage de ses etats entre transport et distribution ::

        d/dt [x_T]   [A_TT  A_TD] [x_T]
             [x_D] = [A_DT  A_DD] [x_D] + B u

    Etudier le transport seul, c'est ne garder que ``dx_T/dt = A_TT x_T +
    A_TD x_D + B_T u`` en faisant passer ``x_D`` du cote des ENTREES : les
    modes du transport sont alors ceux de ``A_TT`` (:meth:`transport`). Le
    couplage ``A_TD``, ``A_DT`` est exactement ce que cette separation
    neglige. Entre les deux, une distribution QUASI STATIQUE
    (``dx_D/dt = 0``) donne ``A_TT - A_TD A_DD^-1 A_DT`` : le complement de
    Schur, dont les equivalents de charge (R-L, hybride) sont des
    approximations. Symetriquement pour la distribution."""

    A: np.ndarray
    etats: list[str]
    transport: np.ndarray  # masque des etats du transport

    def _blocs(self, cote: str):
        t = self.transport if cote == "transport" else ~self.transport
        a, b = np.where(t)[0], np.where(~t)[0]
        A = self.A
        return a, A[np.ix_(a, a)], A[np.ix_(a, b)], A[np.ix_(b, a)], A[np.ix_(b, b)]

    def reduit(self, cote: str = "transport", quasi_statique: bool = False) -> _Spectre:
        """Les modes d'un cote, l'autre FIGE dans les entrees (``A_TT``) ou
        QUASI STATIQUE (complement de Schur)."""
        idx, Aaa, Aab, Aba, Abb = self._blocs(cote)
        M = Aaa
        if quasi_statique and Abb.size:
            try:
                M = Aaa - Aab @ np.linalg.solve(Abb, Aba)
            except np.linalg.LinAlgError:
                M = Aaa - Aab @ np.linalg.lstsq(Abb, Aba, rcond=None)[0]
        return _spectre_matrice(M, [self.etats[i] for i in idx])

    def couplage(self) -> float:
        """La force relative du couplage : ``||A_TD|| ||A_DT|| / (||A_TT||
        ||A_DD||)`` (normes de Frobenius)."""
        _, Att, Atd, Adt, Add = self._blocs("transport")
        n = np.linalg.norm
        return float(n(Atd) * n(Adt) / max(n(Att) * n(Add), 1e-300))


def _spectre_matrice(A: np.ndarray, etats: list[str]) -> _Spectre:
    lam, Vr = np.linalg.eig(A)
    p = np.abs(Vr.T) * np.abs(np.linalg.inv(Vr))
    p = p / np.maximum(p.sum(axis=1, keepdims=True), 1e-30)
    k = np.where((np.abs(lam) > JAUGE) & (lam.imag >= -1e-9))[0]
    return _Spectre(lam[k], p[k], etats)


def partition(r: Reseau, plat: Tables | None = None) -> tuple[Partition, _Spectre]:
    """La partition ``x_T`` / ``x_D`` du modele complet, et son spectre."""
    plat = plat if plat is not None else aplatir(r)
    net = assembler(plat)
    net.solve_equilibrium()
    ss = net.build_state_space()
    eq, _ = equivalent_statique(r, plat=plat)
    noms = _noms_elements(_sans_injections(eq))
    etats = list(ss.states)
    masque = np.array([e.split(".", 1)[0] in noms for e in etats])
    return Partition(np.asarray(ss.A), etats, masque), _spectre_matrice(ss.A, etats)


def comparer_partition(
    r: Reseau,
    cote: str = "transport",
    quasi_statique: bool = False,
    part_min: float = 0.5,
    plat: Tables | None = None,
    part: tuple[Partition, _Spectre] | None = None,
) -> ComparaisonModes:
    """Les modes d'un cote, modele complet contre bloc diagonal (l'autre
    cote dans les entrees) ou complement de Schur (l'autre cote quasi
    statique). Voir :class:`Partition`."""
    P, complet = part if part is not None else partition(r, plat)
    reduit = P.reduit(cote, quasi_statique)
    garde = set(reduit.etats)
    return _comparer(complet, reduit, garde, part_min, cote)


@dataclass
class ModeInteraction:
    """Un mode que transport et distribution portent TOUS DEUX notablement :
    aucun des deux modeles reduits ne peut le representer."""

    lam: complex
    part_transport: float
    elements: list[tuple[str, float]]

    @property
    def f(self) -> float:
        return float(modal.frequence(self.lam)[0])

    @property
    def zeta(self) -> float:
        return float(modal.amortissement(self.lam)[0])


def modes_interaction(
    r: Reseau, complet, plat: Tables, seuil: float = 0.1
) -> list[ModeInteraction]:
    """Les modes oscillants (0,05-``F_MAX`` Hz, une fois par paire) dont la
    part du transport est entre ``seuil`` et ``1 - seuil``."""
    eq, _ = equivalent_statique(r, plat=plat)
    part = _part(complet, _etats_de(complet.etats, _noms_elements(eq)))
    f = modal.frequence(complet.lam)
    f_propre = np.abs(complet.lam) / (2 * np.pi)
    out = []
    for i in np.where(
        (part > seuil) & (part < 1 - seuil) & (f > 0.05) & (f_propre <= F_MAX)
    )[0]:
        out.append(
            ModeInteraction(
                complex(complet.lam[i]),
                float(part[i]),
                _elements(complet.p[i], complet.etats),
            )
        )
    return sorted(out, key=lambda m: m.zeta)


# --------------------------------------------------- balayage de penetration
@dataclass
class Palier:
    n_gfl: int
    p_gfl: float  # [W], production GFL totale dans la distribution (K compris)
    transport: ComparaisonModes | None = None
    #: le transport avec l'equivalent hybride (GFL en puissance constante)
    transport_hybride: ComparaisonModes | None = None
    distribution: list[ComparaisonModes] = field(default_factory=list)
    interaction: list[ModeInteraction] = field(default_factory=list)
    #: la partition x_T / x_D du modele complet : ``(cote, quasi_statique)``
    #: -> comparaison (voir :class:`Partition`)
    partition: dict[tuple[str, bool], ComparaisonModes] = field(default_factory=dict)
    max_re: float = 0.0  # marge du modele complet (modes physiques)
    ajouts: list[tuple[str, int, float]] = field(default_factory=list)
    erreur: str = ""


def _production_gfl(r: Reseau) -> float:
    """La production des GFL de la distribution [W], departs K compris."""
    tot = 0.0
    for inst in r.instances:
        t = r.zones[inst.zone].tables
        gfl = {int(x[1]) for x in (*t.T_GFL1s, *t.T_GFL2s)}
        tot += inst.K * sum(float(x[5]) for x in t.Y_DER if int(x[0]) in gfl)
    return tot


def evaluer(r: Reseau, instances: list[int] | None = None) -> Palier:
    """Un palier : les deux modeles reduits confrontes au complet, et les
    modes d'interaction."""
    instances = list(range(len(r.instances))) if instances is None else instances
    pal = Palier(n_gfl=0, p_gfl=_production_gfl(r))
    plat = aplatir(r)
    P, complet = partition(r, plat)
    lam = complet.lam[np.abs(complet.lam) > JAUGE]
    pal.max_re = float(lam.real.max()) if lam.size else 0.0
    for cote in ("transport", "distribution"):
        for qs in (False, True):
            pal.partition[cote, qs] = comparer_partition(r, cote, qs, part=(P, complet))
    pal.transport = comparer_transport(r, plat=plat, complet=complet)
    pal.transport_hybride = comparer_transport(
        r, plat=plat, complet=complet, hybride=True
    )
    for k in instances:
        # une distribution seule sans point de fonctionnement (load flow du
        # modele reduit) est ecartee et signalee, sans arreter le palier
        try:
            pal.distribution.append(
                comparer_distribution(r, k, plat=plat, complet=complet)
            )
        except RuntimeError as err:
            pal.erreur += f"distribution {k} : {err} ; "
    pal.interaction = modes_interaction(r, complet, plat)
    return pal


def penetration(
    r: Reseau,
    n_max: int = 30,
    pas: int = 5,
    graine: int = 1,
    part: float = 0.05,
    instances: list[int] | None = None,
) -> list[Palier]:
    """Ajoute des GFL par paquets de ``pas`` (graine fixe : les memes a
    chaque appel) et, a chaque palier, confronte les deux modeles reduits au
    modele complet. ``instances`` : les distributions a etudier seules
    (toutes par defaut). Un palier sans point de fonctionnement porte son
    ``erreur`` et arrete le balayage."""
    paliers = []
    courant, n, ajouts = r, 0, []
    while True:
        try:
            pal = evaluer(courant, instances)
        except SANS_POINT_DE_FONCTIONNEMENT as err:
            pal = Palier(
                n_gfl=n,
                p_gfl=_production_gfl(courant),
                erreur=str(err).splitlines()[0][:200],
            )
        pal.n_gfl, pal.ajouts = n, ajouts
        paliers.append(pal)
        if pal.erreur or n >= n_max:
            break
        courant, ajouts = ajouter_gfl(courant, pas, graine=graine + n, part=part)
        n += pas
    return paliers
