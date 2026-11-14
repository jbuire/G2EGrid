"""L'impact des reseaux de distribution sur le reseau de transport.

LA QUESTION -- un exploitant de transport modelise d'ordinaire la
distribution comme une CHARGE : une puissance soutiree au poste, sans
dynamique. Tant que la distribution etait passive, l'erreur etait faible.
Avec des convertisseurs dans la distribution, elle ne l'est plus forcement :
la distribution peut deplacer les modes electromecaniques du transport, les
amortir ou les desamortir, et en creer de nouveaux -- dont certains
dangereux. Ce module mesure tout cela, sur un reseau decrit par zones.

LA METHODE, et pourquoi celle-la -- on compare le modele COMPLET a un
EQUIVALENT STATIQUE ou chaque instance de zone est remplacee par la puissance
qu'elle soutire reellement au poste, lue sur le load flow du modele complet.
Le point de fonctionnement du transport est alors IDENTIQUE dans les deux
modeles (verifie : memes tensions aux bus du transport) : la seule difference
est la dynamique de la distribution. Un equivalent construit sur les charges
nominales melangerait deux effets -- le deplacement du point de
fonctionnement, et la dynamique -- et on ne saurait pas lequel on mesure.

Ce qui est rendu, pour chaque mode du transport :
  - son deplacement en frequence et en amortissement quand la distribution
    devient dynamique ;
  - la PART de la distribution dans ce mode (participations de ses etats) :
    un mode ou la distribution pese 40 % n'est plus un mode du transport
    seul, et un correcteur du transport ne suffira pas a l'amortir ;

et, a part, les modes que la distribution CREE -- ceux que l'equivalent n'a
pas, domines par des etats de distribution -- avec une alerte pour tout mode
peu amorti ou instable.

``par_poste`` refait la comparaison en ne remplacant qu'une instance a la
fois : c'est ce qui dit QUEL poste deplace quel mode.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import linear_sum_assignment

from g2egrid.analyse import modal
from g2egrid.reseau.flux import flux_branches
from g2egrid.reseau.tables import TABLES, Tables, assembler
from g2egrid.reseau.zones import Instance, Reseau, aplatir

__all__ = [
    "ModeCompare",
    "ModeCree",
    "Impact",
    "bus_d_instance",
    "equivalent_statique",
    "impact_distribution",
]

#: En dessous, une valeur propre est un mode de jauge (reference d'angle
#: flottante), pas un mode physique.
JAUGE = 1e-6

#: Colonnes qui designent un bus, table par table.
_COLONNES_BUS = ("BusID", "From", "To", "BusAC", "BusDC")


@dataclass
class ModeCompare:
    """Un mode du transport, sans puis avec la dynamique de la distribution."""

    f_equivalent: float  # Hz
    zeta_equivalent: float  # %
    f_complet: float
    zeta_complet: float
    part_distribution: float  # % de participation des etats de distribution
    dominant: str  # l'etat le plus participant, dans le modele complet
    lam_equivalent: complex = 0j
    lam_complet: complex = 0j

    @property
    def ecart_f(self) -> float:
        """[%] de la frequence de l'equivalent."""
        return (
            100 * (self.f_complet - self.f_equivalent) / max(self.f_equivalent, 1e-12)
        )

    @property
    def ecart_zeta(self) -> float:
        """[points] d'amortissement : ce qu'un cahier des charges regarde."""
        return self.zeta_complet - self.zeta_equivalent


@dataclass
class ModeCree:
    """Un mode du modele complet que l'equivalent statique n'a pas."""

    f: float
    zeta: float
    part_distribution: float
    dominant: str
    lam: complex = 0j

    @property
    def alerte(self) -> bool:
        return self.zeta < 5.0


@dataclass
class Impact:
    paires: list[ModeCompare]
    crees: list[ModeCree]
    #: des modes qui existaient sans la distribution, mais qu'elle rend
    #: instables ou peu amortis -- l'autre moitie des alertes
    degrades: list[ModeCompare]
    #: la puissance soutiree a chaque poste, en pu systeme, par instance
    soutirage: dict[str, complex]
    #: ecart maximal de tension aux bus du transport, equivalent / complet :
    #: la preuve que les deux modeles ont le MEME point de fonctionnement
    ecart_tension: float
    max_re_equivalent: float
    max_re_complet: float
    #: instance -> [(f de l'equivalent, ecart f %, ecart zeta pts)] quand SEULE
    #: cette instance est rendue dynamique
    par_poste: dict[str, list[tuple[float, float, float]]] = field(default_factory=dict)
    equivalent: Tables | None = None


def _nom_instance(inst: Instance) -> str:
    return f"{inst.zone}@{inst.bus_amont}"


def bus_d_instance(r: Reseau, inst: Instance) -> set[int]:
    """Les numeros qu'occupe une instance dans le reseau APLATI.

    Sans transfo de poste, le PCC de la zone fusionne avec le bus amont : il
    appartient alors au transport, et n'est pas compte ici.
    """
    z = r.zones[inst.zone]
    locaux = [int(x[0]) for x in z.tables.Y_DER] + [int(x[0]) for x in z.tables.T_BusDC]
    out = set()
    for loc in locaux:
        if inst.transfo is None and loc == int(z.bus_pcc):
            continue
        out.add(int(inst.prefixe) + loc)
    return out


def _retirer(t: Tables, bus: set[int]) -> Tables:
    """Une copie des tables sans les lignes qui touchent ``bus``."""
    t2 = t.copie()
    for nom, cols in TABLES.items():
        ic = [j for j, c in enumerate(cols) if c in _COLONNES_BUS]
        if nom in ("Y_DER", "T_BusDC"):
            ic = [0]
        rows = getattr(t2, nom)
        setattr(t2, nom, [r for r in rows if not any(int(r[j]) in bus for j in ic)])
    t2.PARAMS_PAR_BUS = {
        b: v for b, v in t2.PARAMS_PAR_BUS.items() if int(b) not in bus
    }
    return t2


def equivalent_statique(
    r: Reseau, remplacer: list[int] | None = None, plat: Tables | None = None
) -> tuple[Tables, dict[str, complex]]:
    """Le reseau ou les instances ``remplacer`` deviennent des charges.

    Chaque instance remplacee disparait -- ses bus, ses branches, ses
    equipements, son transfo de poste -- et le bus amont recoit en charge
    la puissance que l'instance y soutirait dans le modele COMPLET. C'est ce
    qui garde le point de fonctionnement du transport inchange.

    Rend les tables et, par instance remplacee, la puissance soutiree (pu).
    """
    plat = plat if plat is not None else aplatir(r)
    remplacer = list(range(len(r.instances))) if remplacer is None else list(remplacer)
    net = assembler(plat)
    lf = net.run_loadflow()
    transits = flux_branches(net, lf)

    a_retirer: set[int] = set()
    soutirage: dict[str, complex] = {}
    ajout: dict[int, complex] = {}
    for k in remplacer:
        inst = r.instances[k]
        bus = bus_d_instance(r, inst)
        amont = int(inst.bus_amont)
        s = 0j
        for tr in transits:
            if tr.de == amont and tr.vers in bus:
                s += complex(tr.P_de, tr.Q_de)
            elif tr.vers == amont and tr.de in bus:
                s += complex(tr.P_vers, tr.Q_vers)
        soutirage[_nom_instance(inst)] = s
        ajout[amont] = ajout.get(amont, 0j) + s
        a_retirer |= bus

    eq = _retirer(plat, a_retirer)
    lignes = []
    for row in eq.Y_DER:
        row = list(row)
        s = ajout.get(int(row[0]))
        if s is not None:
            row[7] = row[7] + s.real * plat.Sb
            row[8] = row[8] + s.imag * plat.Sb
        lignes.append(tuple(row))
    eq.Y_DER = lignes
    return eq, soutirage


def _spectre(t: Tables):
    net = assembler(t)
    net.solve_equilibrium()
    ss = net.build_state_space()
    lam, _, p = modal.spectre(ss)
    p = p / np.maximum(p.sum(axis=1, keepdims=True), 1e-30)
    return net, ss, lam, p


def _etats_distribution(
    net, ss, bus_distribution: set[int], amonts: set[int]
) -> np.ndarray:
    """Masque des etats portes par la distribution.

    Un element est de la distribution si tous ses bus sont dans la
    distribution ou au poste amont, et au moins un dans la distribution --
    ce qui range le transfo de poste cote distribution : il disparait avec
    l'instance dans l'equivalent.
    """
    a_elle = {}
    for el in net.all_elements:
        bs = set(int(b) for b in el.buses)
        a_elle[el.name] = bool(bs & bus_distribution) and bs <= (
            bus_distribution | amonts
        )
    return np.array([a_elle.get(n.split(".")[0], False) for n in ss.states])


def _garde(lam, f_max, zeta_max):
    f = modal.frequence(lam)
    z = modal.amortissement(lam)
    return (np.abs(lam) > JAUGE) & (lam.imag > 1e-6) & (f <= f_max) & (z <= zeta_max)


def _apparier(lam_eq, lam_c, seuil=0.3):
    cout = np.abs(lam_eq[:, None] - lam_c[None, :]) / np.maximum(
        np.abs(lam_eq[:, None]), 1e-9
    )
    i, j = linear_sum_assignment(cout)
    ok = cout[i, j] < seuil
    return i[ok], j[ok]


def impact_distribution(
    r: Reseau,
    f_max: float = 5.0,
    zeta_max: float = 95.0,
    par_poste: bool = True,
    f_alerte: float = 1000.0,
) -> Impact:
    """Mesure ce que la dynamique de la distribution change au transport.

    ``f_alerte`` borne les alertes : au-dela, les modes peu amortis sont les
    resonances electromagnetiques des cables (dizaines de kHz, ~0 %
    d'amortissement), presentes dans tout modele de reseau et hors du propos
    ici. Un mode INSTABLE, lui, est signale a toute frequence.
    """
    plat = aplatir(r)
    net_c, ss_c, lam_c, p_c = _spectre(plat)
    eq, soutirage = equivalent_statique(r, plat=plat)
    net_e, ss_e, lam_e, _ = _spectre(eq)

    # preuve : le meme point de fonctionnement aux bus du transport
    idx_c, idx_e = net_c._bus_index(), net_e._bus_index()
    commun = [b for b in idx_e if b in idx_c]
    ecart_v = max(abs(net_c.lf.V[idx_c[b]] - net_e.lf.V[idx_e[b]]) for b in commun)

    bus_dist = (
        set().union(*(bus_d_instance(r, i) for i in r.instances))
        if r.instances
        else set()
    )
    amonts = {int(i.bus_amont) for i in r.instances}
    dist = _etats_distribution(net_c, ss_c, bus_dist, amonts)
    part = 100 * (p_c[:, dist].sum(axis=1) if dist.any() else np.zeros(len(lam_c)))

    # On apparie les spectres ENTIERS, pas seulement la bande etudiee : un
    # mode du complet n'est "cree" par la distribution que s'il n'a de
    # partenaire nulle part dans l'equivalent. N'apparier que 0-5 Hz faisait
    # passer pour crees les resonances electromagnetiques du TRANSPORT
    # (bus4.vq a 993 Hz...), qui existent dans les deux modeles.
    osc_e = np.where((np.abs(lam_e) > JAUGE) & (lam_e.imag > 1e-6))[0]
    osc_c = np.where((np.abs(lam_c) > JAUGE) & (lam_c.imag > 1e-6))[0]
    ie, ic = _apparier(lam_e[osc_e], lam_c[osc_c])
    fe, ze = modal.frequence(lam_e), modal.amortissement(lam_e)
    fc, zc = modal.frequence(lam_c), modal.amortissement(lam_c)
    noms = list(ss_c.states)
    dans_bande_e = _garde(lam_e, f_max, zeta_max)
    dangereux = lambda lam, f, z: lam.real > 0 or (z < 5.0 and f <= f_alerte)  # noqa: E731

    paires, degrades = [], []
    for a, b in zip(ie, ic):
        k_e, k_c = osc_e[a], osc_c[b]
        m = ModeCompare(
            float(fe[k_e]),
            float(ze[k_e]),
            float(fc[k_c]),
            float(zc[k_c]),
            float(part[k_c]),
            noms[int(np.argmax(p_c[k_c]))],
            complex(lam_e[k_e]),
            complex(lam_c[k_c]),
        )
        if dans_bande_e[k_e]:
            paires.append(m)
        # un mode qui existait deja mais que la distribution rend dangereux
        if dangereux(lam_c[k_c], fc[k_c], zc[k_c]) and not dangereux(
            lam_e[k_e], fe[k_e], ze[k_e]
        ):
            degrades.append(m)
    paires.sort(key=lambda m: m.zeta_equivalent)

    apparies = set(osc_c[ic])
    crees = []
    for k in osc_c:
        if k in apparies:
            continue
        dans_bande = fc[k] <= f_max and zc[k] <= zeta_max and part[k] > 50.0
        if dangereux(lam_c[k], fc[k], zc[k]) or dans_bande:
            crees.append(
                ModeCree(
                    float(fc[k]),
                    float(zc[k]),
                    float(part[k]),
                    noms[int(np.argmax(p_c[k]))],
                    complex(lam_c[k]),
                )
            )
    crees.sort(key=lambda m: m.zeta)
    ge = np.where(dans_bande_e)[0]

    phys = lambda lam: float(lam.real[np.abs(lam) > JAUGE].max())  # noqa: E731
    out = Impact(
        paires=paires,
        crees=crees,
        degrades=degrades,
        soutirage=soutirage,
        ecart_tension=float(ecart_v),
        max_re_equivalent=phys(lam_e),
        max_re_complet=phys(lam_c),
        equivalent=eq,
    )

    if par_poste and len(r.instances) > 1:
        tous = list(range(len(r.instances)))
        for k, inst in enumerate(r.instances):
            # SEULE l'instance k reste dynamique : on remplace toutes les autres
            eq_k, _ = equivalent_statique(r, [j for j in tous if j != k], plat=plat)
            _, _, lam_k, _ = _spectre(eq_k)
            gk = np.where((np.abs(lam_k) > JAUGE) & (lam_k.imag > 1e-6))[0]
            a, b = _apparier(lam_e[ge], lam_k[gk])
            fk, zk = modal.frequence(lam_k), modal.amortissement(lam_k)
            out.par_poste[_nom_instance(inst)] = [
                (
                    float(fe[ge[x]]),
                    float(100 * (fk[gk[y]] - fe[ge[x]]) / max(fe[ge[x]], 1e-12)),
                    float(zk[gk[y]] - ze[ge[x]]),
                )
                for x, y in zip(a, b)
            ]
    return out
