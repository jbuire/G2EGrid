"""Fatigue mecanique d'un arbre : comptage rainflow, courbe S-N, cumul de Miner.

POURQUOI -- un groupe electrogene alimente par un redresseur et un onduleur voit
la charge informatique TELLE QU'ELLE EST : des variations de puissance de
plusieurs dizaines de pour cent en quelques secondes, repetees des milliers de
fois. Chaque variation tord l'accouplement entre le moteur et l'alternateur ;
la tenue de l'arbre se juge sur le NOMBRE et l'AMPLITUDE de ces cycles, pas sur
leur maximum. Le couple d'arbre est la sortie ``C_arbre`` du modele
(:class:`~g2egrid.elements.GroupeElectrogene`) ; ce module en tire un
endommagement.

LA CHAINE DE CALCUL ::

    couple d'arbre (pu)  ->  couple [N.m]  ->  contrainte de cisaillement [Pa]
        ->  comptage rainflow : (amplitude, moyenne, compte) de chaque cycle
        ->  correction de contrainte moyenne (Goodman)
        ->  courbe S-N (Basquin, coude de Haibach) : cycles a rupture N_i
        ->  regle de Miner : D = somme(n_i / N_i)

D = 1 est la rupture. La duree de vie en repetant la meme sollicitation vaut
``duree / D``.

CE QUE C'EST ET CE QUE CE N'EST PAS -- la methode est celle qu'emploient les
normes (rainflow ASTM E1049-85, regle lineaire de Miner), avec ses limites
connues : la regle de Miner ignore l'ordre des cycles, et le coude de la
courbe S-N est un modele. Les DONNEES MATERIAU par defaut
(:data:`COURBE_ACIER_ARBRE`, :data:`ARBRE_DEFAUT`) sont des ORDRES DE GRANDEUR
d'un acier d'arbre allie, pas ceux d'un arbre precis : le dommage absolu en
depend a la puissance ``m`` (``m`` ~ 7 : un facteur 2 sur la contrainte change
le dommage d'un facteur 130). Ce que le calcul donne de robuste, c'est la
COMPARAISON entre deux sollicitations sur le meme arbre (avec et sans lissage
par batterie, avec deux reprises de charge). Pour un chiffre absolu, remplacer
les donnees par la fiche de l'accouplement.

Le modele dynamique est un modele MOYEN : le couple d'arbre qu'il calcule ne
contient ni la pulsation d'allumage des cylindres ni les harmoniques 6 du
redresseur. Ces composantes, a haute frequence, ajoutent des cycles en grand
nombre mais de faible amplitude ; elles se superposent par l'entree
``dT_exc`` du groupe, ou se traitent a part.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import pairwise

import numpy as np

__all__ = [
    "ARBRE_DEFAUT",
    "COURBE_ACIER_ARBRE",
    "Arbre",
    "CourbeSN",
    "Cycles",
    "ResultatFatigue",
    "couple_base",
    "evaluer_fatigue",
    "limite_non_resonante_epri",
    "points_de_retournement",
    "rainflow",
    "securite_goodman",
    "tension_cisaillement",
]


# ----------------------------------------------------------------------
# unites
# ----------------------------------------------------------------------
def couple_base(Sn_VA: float, f0: float = 50.0, paires_poles: int = 2) -> float:
    """Le couple de base [N.m] d'une machine : ``Sn / omega_mecanique``.

    ``omega_mecanique = 2.pi.f0 / paires_poles`` (157 rad/s, 1500 tr/min pour
    deux paires de poles a 50 Hz). Un couple de 1 pu est celui de la puissance
    apparente nominale a la vitesse nominale.
    """
    if paires_poles < 1:
        raise ValueError("paires_poles >= 1")
    return float(Sn_VA) / (2.0 * np.pi * f0 / paires_poles)


def tension_cisaillement(
    couple_Nm, diametre_m: float, diametre_interieur_m: float = 0.0
):
    """Contrainte de cisaillement maximale [Pa] d'un arbre en torsion, en peau :
    ``tau = C.(D/2)/J`` avec ``J = pi.(D^4 - d^4)/32`` ; arbre plein (``d = 0``) :
    ``tau = 16.C/(pi.D^3)``. Les gros arbres d'alternateurs hydrauliques sont
    creux (Grand Coulee GC-19 : 2540 / 2134 mm)."""
    D, d = float(diametre_m), float(diametre_interieur_m)
    if not 0.0 <= d < D:
        raise ValueError("diametres : 0 <= interieur < exterieur")
    J = np.pi * (D**4 - d**4) / 32.0
    return np.asarray(couple_Nm, dtype=float) * (D / 2.0) / J


def securite_goodman(tau_a, tau_m, S_se: float, S_su: float):
    """Coefficient de securite de Goodman en cisaillement,
    ``1/n = tau_a/S_se + tau_m/S_su``.

    ``S_se`` : limite d'endurance en cisaillement de la piece (facteurs de Marin
    compris), ``S_su`` : resistance ultime en cisaillement. ``n > 1`` : vie
    infinie au sens du critere ; c'est la metrique de criblage du PNNL
    (arXiv 2607.14412, eq. 12)."""
    return 1.0 / (
        np.asarray(tau_a, dtype=float) / S_se + np.asarray(tau_m, dtype=float) / S_su
    )


def limite_non_resonante_epri(
    UTS_Pa: float,
    tau_m_Pa: float = 50e6,
    rapport_cisaillement: float = 0.29,
    marin: float = 0.5,
    Kt: float = 4.0,
    SF: float = 2.0,
) -> float:
    """Oscillation de puissance mecanique admissible hors resonance, en fraction
    de la puissance nominale, par la methode simplifiee de l'EPRI (rapport
    3002036374, annexe A) ::

        tau'_e = rapport.UTS        limite d'endurance en cisaillement, eprouvette polie
        tau_a  = tau'_e.marin/(Kt.SF)   admissible nominale, piece reelle
        dP/P   = dT/T = tau_a / tau_m   tau_m : cisaillement moyen au couple nominal

    Avec les hypotheses de l'EPRI (0,29 UTS, Marin 0,5, Kt = 4 pour les clavetages
    frettes, SF = 2, tau_m = 50 MPa), un acier de 80 ksi (552 MPa) donne 20 %."""
    return rapport_cisaillement * UTS_Pa * marin / (Kt * SF) / tau_m_Pa


# ----------------------------------------------------------------------
# comptage rainflow
# ----------------------------------------------------------------------
def points_de_retournement(serie, seuil: float = 0.0) -> np.ndarray:
    """Les extrema successifs d'une serie, en ignorant les retournements
    de moins de ``seuil``.

    Le filtre est une HYSTERESIS : un extremum n'est retenu que si le signal
    s'en eloigne d'au moins ``seuil`` ; les oscillations plus petites sont du
    bruit numerique ou de mesure, qui gonflerait le comptage de milliers de
    cycles sans dommage. Le premier et le dernier point sont toujours gardes.
    """
    x = np.asarray(serie, dtype=float)
    if x.ndim != 1:
        raise ValueError("serie : un vecteur")
    if len(x) < 2:
        return x.copy()
    if seuil < 0.0:
        raise ValueError("seuil >= 0")
    pts = [x[0]]
    ext = x[0]  # extremum courant, dans le sens courant
    sens = 0  # 0 : indetermine, +1 : on monte, -1 : on descend
    for v in x[1:]:
        if sens == 0:
            if v - ext > seuil:
                sens, ext = +1, v
            elif ext - v > seuil:
                sens, ext = -1, v
        elif sens > 0:
            if v > ext:
                ext = v
            elif ext - v > seuil:  # le sommet est confirme
                pts.append(ext)
                sens, ext = -1, v
        else:
            if v < ext:
                ext = v
            elif v - ext > seuil:  # le creux est confirme
                pts.append(ext)
                sens, ext = +1, v
    pts.append(ext)
    pts = np.array(pts)
    # supprime les doublons consecutifs
    garde = np.concatenate(([True], np.diff(pts) != 0.0))
    return pts[garde]


@dataclass
class Cycles:
    """Le resultat du comptage : un cycle par ligne.

    ``plage`` : amplitude crete a crete ; ``moyenne`` : valeur moyenne du cycle ;
    ``compte`` : 1 pour un cycle complet, 0,5 pour un demi-cycle.
    """

    plage: np.ndarray
    moyenne: np.ndarray
    compte: np.ndarray

    @property
    def n(self) -> float:
        """Nombre de cycles equivalents (un demi-cycle compte 0,5)."""
        return float(self.compte.sum())

    @property
    def amplitude(self) -> np.ndarray:
        return self.plage / 2.0

    def __len__(self) -> int:
        return len(self.plage)


def rainflow(serie, seuil: float = 0.0) -> Cycles:
    """Comptage rainflow a trois points (ASTM E1049-85, article 5.4.4).

    Les cycles fermes comptent 1, ceux du residu (les montees et descentes qui
    ne se referment pas) 0,5. A seuil nul, la somme des comptes est
    ``(n_retournements - 1) / 2``.
    """
    pts = points_de_retournement(serie, seuil)
    pile: list[float] = []
    plage: list[float] = []
    moy: list[float] = []
    cpt: list[float] = []
    for p in pts:
        pile.append(float(p))
        while len(pile) >= 3:
            X = abs(pile[-1] - pile[-2])
            Y = abs(pile[-2] - pile[-3])
            if X < Y:
                break
            m = 0.5 * (pile[-2] + pile[-3])
            if len(pile) == 3:
                # Y contient le point de depart : demi-cycle, on l'ecarte
                plage.append(Y), moy.append(m), cpt.append(0.5)
                del pile[0]
            else:
                plage.append(Y), moy.append(m), cpt.append(1.0)
                del pile[-3:-1]
    for a, b in pairwise(pile):
        plage.append(abs(b - a)), moy.append(0.5 * (a + b)), cpt.append(0.5)
    return Cycles(np.array(plage), np.array(moy), np.array(cpt))


# ----------------------------------------------------------------------
# courbe S-N
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class CourbeSN:
    """Courbe de Wohler en cisaillement, loi de Basquin avec coude.

    Au-dessus du coude, ``N = N_D.(S_D/S)^m`` ; sous le coude (``S < S_D``),
    la pente passe a ``m2``. ``m2 = 2m - 1`` est la modification de Haibach :
    les cycles sous la limite d'endurance ne sont pas gratuits sous une charge
    d'amplitude variable, car les gros cycles fissurent et abaissent la limite.
    ``m2 = inf`` est la regle elementaire de Miner : aucun dommage sous ``S_D``.

    ``S_D`` est l'AMPLITUDE de contrainte [Pa] tenue ``N_D`` cycles, pour une
    contrainte moyenne nulle.
    """

    S_D: float
    N_D: float
    m: float
    m2: float = np.inf
    nom: str = ""

    def __post_init__(self) -> None:
        if not (self.S_D > 0.0 and self.N_D > 0.0 and self.m > 1.0 and self.m2 > 1.0):
            raise ValueError("courbe S-N : S_D, N_D > 0 ; pentes > 1")

    def cycles_a_rupture(self, amplitude) -> np.ndarray:
        """``N`` pour chaque amplitude de contrainte [Pa] ; ``inf`` si elle est
        nulle, ou si elle est sous le coude avec ``m2 = inf``."""
        S = np.abs(np.asarray(amplitude, dtype=float))
        N = np.full(S.shape, np.inf)
        pos = S > 0.0
        if np.isinf(self.m2):
            haut = pos & (S >= self.S_D)
            N[haut] = self.N_D * (self.S_D / S[haut]) ** self.m
            return N
        pente = np.where(S >= self.S_D, self.m, self.m2)
        N[pos] = self.N_D * (self.S_D / S[pos]) ** pente[pos]
        return N


#: Acier d'arbre allie traite (type 42CrMo4), cisaillement, avec effets d'entaille,
#: de taille et d'etat de surface deja deduits : amplitude de 90 MPa tenue 2.10^6
#: cycles, pente 7 (pente de Haibach 13 sous le coude). ORDRE DE GRANDEUR, a
#: remplacer par les donnees de l'arbre etudie.
COURBE_ACIER_ARBRE = CourbeSN(
    S_D=90e6, N_D=2e6, m=7.0, m2=13.0, nom="acier allie d'arbre (ordre de grandeur)"
)


@dataclass(frozen=True)
class Arbre:
    """L'arbre (ou l'accouplement) dont on juge la tenue."""

    diametre_m: float = 0.11
    #: diametre interieur [m] d'un arbre creux (0 : plein)
    diametre_interieur_m: float = 0.0
    #: resistance ultime au cisaillement [Pa], pour la correction de Goodman
    tau_u: float = 600e6
    #: paires de poles de la machine (vitesse = 2.pi.f0 / paires)
    paires_poles: int = 2


ARBRE_DEFAUT = Arbre()


# ----------------------------------------------------------------------
# endommagement
# ----------------------------------------------------------------------
@dataclass
class ResultatFatigue:
    """Le bilan de fatigue d'une serie de couple."""

    duree_s: float
    #: nombre de cycles equivalents comptes (demi-cycles pour 0,5)
    n_cycles: float
    #: dommage de Miner de la serie (1 = rupture)
    dommage: float
    #: contrainte de cisaillement maximale et moyenne [Pa]
    tau_max: float
    tau_moyen: float
    #: amplitude [Pa] du plus gros cycle, et amplitude equivalente de Miner
    amplitude_max: float
    amplitude_equivalente: float
    #: les cycles : amplitude de contrainte EQUIVALENTE apres Goodman [Pa], compte,
    #: dommage de chaque cycle
    cycles: Cycles
    amplitude_goodman: np.ndarray = field(
        repr=False, default_factory=lambda: np.zeros(0)
    )
    dommage_par_cycle: np.ndarray = field(
        repr=False, default_factory=lambda: np.zeros(0)
    )
    courbe: CourbeSN = COURBE_ACIER_ARBRE
    #: contrainte moyenne [Pa] de chaque cycle et resistance ultime, pour refaire
    #: le calcul a une autre amplitude (:meth:`marge_amplitude`)
    moyenne: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(0))
    tau_u: float = 600e6
    goodman: bool = True

    @property
    def dommage_par_heure(self) -> float:
        """Dommage accumule par heure si la serie se repete sans fin."""
        return self.dommage * 3600.0 / self.duree_s

    @property
    def vie_heures(self) -> float:
        """Heures de fonctionnement jusqu'a D = 1 en repetant la serie."""
        d = self.dommage_par_heure
        return float(np.inf) if d <= 0.0 else 1.0 / d

    def dommage_a_l_echelle(self, k: float) -> float:
        """Le dommage si TOUTES les amplitudes de contrainte etaient multipliees
        par ``k`` (les moyennes restant ce qu'elles sont)."""
        S_a = self.cycles.amplitude * k
        if self.goodman:
            S_m = np.maximum(self.moyenne, 0.0)
            with np.errstate(divide="ignore"):
                S_eq = np.where(
                    S_m < self.tau_u, S_a / (1.0 - S_m / self.tau_u), np.inf
                )
        else:
            S_eq = S_a
        with np.errstate(divide="ignore", invalid="ignore"):
            d = self.cycles.compte / self.courbe.cycles_a_rupture(S_eq)
        return float(np.sum(np.where(np.isnan(d), 0.0, d)))

    def marge_amplitude(self, vie_cible_h: float = 25 * 8760.0) -> float:
        """Le facteur ``k`` par lequel on peut multiplier les amplitudes de
        contrainte de la serie avant que sa vie, repetee sans fin, tombe a
        ``vie_cible_h`` heures (25 ans par defaut).

        ``k > 1`` : de la marge ; ``k < 1`` : la serie use l'arbre plus vite que la
        cible. C'est l'indicateur robuste du calcul, bien plus que le dommage
        absolu, qui varie comme ``k^m`` (``m ~ 7``) avec les donnees materiau.
        Borne a 1e3. Obtenu par bissection sur ``log k`` : le dommage croit avec ``k``.
        """
        cible = self.duree_s / (3600.0 * vie_cible_h)  # dommage de la serie a la cible
        if self.dommage_a_l_echelle(1e3) < cible:
            return 1e3
        lo, hi = -7.0, float(np.log(1e3))
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if self.dommage_a_l_echelle(float(np.exp(mid))) < cible:
                lo = mid
            else:
                hi = mid
        return float(np.exp(0.5 * (lo + hi)))

    def histogramme(
        self, n_classes: int = 12
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Le spectre de charge : ``(bornes, comptes, dommages)`` par classe
        d'amplitude Goodman. Les bornes sont celles de ``n_classes`` classes de
        meme largeur entre 0 et l'amplitude maximale."""
        a = self.amplitude_goodman
        if not len(a) or a.max() <= 0.0:
            b = np.linspace(0.0, 1.0, n_classes + 1)
            return b, np.zeros(n_classes), np.zeros(n_classes)
        b = np.linspace(0.0, a.max() * (1 + 1e-12), n_classes + 1)
        k = np.clip(np.digitize(a, b) - 1, 0, n_classes - 1)
        comptes = np.bincount(k, weights=self.cycles.compte, minlength=n_classes)
        dom = np.bincount(k, weights=self.dommage_par_cycle, minlength=n_classes)
        return b, comptes, dom

    def __repr__(self) -> str:
        return (
            f"<ResultatFatigue {self.duree_s:g} s, {self.n_cycles:.1f} cycles, "
            f"D = {self.dommage:.3e}, vie {self.vie_heures:.3g} h>"
        )


def evaluer_fatigue(
    t,
    couple_pu,
    Sn_VA: float,
    f0: float = 50.0,
    arbre: Arbre = ARBRE_DEFAUT,
    courbe: CourbeSN = COURBE_ACIER_ARBRE,
    seuil_pu: float = 1e-3,
    goodman: bool = True,
) -> ResultatFatigue:
    """Evalue le dommage de fatigue d'un arbre soumis a un couple ``couple_pu(t)``.

    ``couple_pu`` : couple d'arbre en pu de la machine (``C_arbre``). ``seuil_pu`` :
    les retournements de moins de ce couple sont ignores (hysteresis du
    comptage ; 1e-3 pu = 0,1 % du couple nominal, sous le bruit de calcul).
    ``goodman`` : ramene chaque cycle a une contrainte moyenne nulle,
    ``S_eq = S_a / (1 - S_m/tau_u)`` ; un cycle dont la moyenne atteint
    ``tau_u`` est declare destructif (dommage infini).
    """
    t = np.asarray(t, dtype=float)
    c = np.asarray(couple_pu, dtype=float)
    if t.shape != c.shape or len(t) < 2:
        raise ValueError("t et couple_pu : meme longueur, >= 2 points")
    if not np.all(np.isfinite(c)):
        raise ValueError("couple_pu contient des valeurs non finies")
    cb = couple_base(Sn_VA, f0, arbre.paires_poles)
    tau = tension_cisaillement(c * cb, arbre.diametre_m, arbre.diametre_interieur_m)
    cycles = rainflow(
        tau,
        float(
            tension_cisaillement(
                seuil_pu * cb, arbre.diametre_m, arbre.diametre_interieur_m
            )
        ),
    )
    S_a = cycles.amplitude
    if goodman:
        S_m = np.maximum(cycles.moyenne, 0.0)  # une moyenne negative ne protege pas
        with np.errstate(divide="ignore"):
            S_eq = np.where(S_m < arbre.tau_u, S_a / (1.0 - S_m / arbre.tau_u), np.inf)
    else:
        S_eq = S_a
    with np.errstate(divide="ignore", invalid="ignore"):
        d_i = cycles.compte / courbe.cycles_a_rupture(S_eq)
    d_i = np.where(np.isnan(d_i), 0.0, d_i)
    dommage = float(d_i.sum())
    m = courbe.m
    fin = np.isfinite(S_eq)
    amp_eq = (
        float(
            (np.sum(cycles.compte[fin] * S_eq[fin] ** m) / cycles.compte[fin].sum())
            ** (1.0 / m)
        )
        if fin.any() and cycles.compte[fin].sum() > 0.0
        else 0.0
    )
    return ResultatFatigue(
        duree_s=float(t[-1] - t[0]),
        n_cycles=cycles.n,
        dommage=dommage,
        tau_max=float(np.abs(tau).max()),
        tau_moyen=float(tau.mean()),
        amplitude_max=float(S_a.max()) if len(S_a) else 0.0,
        amplitude_equivalente=amp_eq,
        cycles=cycles,
        amplitude_goodman=S_eq,
        dommage_par_cycle=d_i,
        courbe=courbe,
        moyenne=cycles.moyenne,
        tau_u=arbre.tau_u,
        goodman=goodman,
    )
