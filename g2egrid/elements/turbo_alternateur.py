"""Alternateur synchrone a ligne d'arbre MULTI-MASSES : turbo-alternateurs et
groupes hydrauliques.

L'ELECTRIQUE EST CELUI DE ``SynchronousMachine`` (17 etats, inchanges, regulateur
de vitesse et AVR compris). La masse unique du rotor est remplacee par une
chaine de masses reliees par des troncons d'arbre elastiques -- le modele a masses
concentrees des etudes de resonance sous-synchrone (IEEE First Benchmark Model,
1977) et des papiers recents sur la fatigue des lignes d'arbres face aux centres
de donnees (Hossain et al., arXiv 2605.01173, eq. 1 ; Chatterjee et al., arXiv
2607.14412) ::

    dd_i/dt          = wb.(Dw_i - Dw_(i+1))                       torsion du troncon i
    C_i              = K_i.d_i + Dm_i.(Dw_i - Dw_(i+1))           couple transmis de i a i+1
    2.H_r.dDw_r/dt   = Cm_r + C_(r-1) - C_r - D_r.Dw_r - [Ce si r = alternateur]

``d_i`` en radians ELECTRIQUES, ``K`` en pu de couple par radian electrique,
``H`` en s, tout en base de la machine (``Sn``). La masse de l'alternateur porte
l'inertie ``H`` et l'amortissement ``KD`` de ``SynchronousMachine``, que la ligne
fixe par defaut a ``H`` et ``D`` de sa masse d'alternateur ; le couple
moteur ``Pm`` (sortie du regulateur de vitesse) se repartit entre les etages de
turbine selon les fractions ``F_r`` (somme 1) ::

    Cm_r = F_r.Pm / (w_ref + kCm.Dw_r)

``kCm = 0`` (defaut) : couple moteur INDEPENDANT de la vitesse, comme dans les
etudes de resonance sous-synchrone (couple mecanique constant du FBM) et chez le
PNNL, dont l'amortissement ``D`` de la roue contient deja l'autoreglage de
l'eau. ``kCm = 1`` : ``Cm = Pm/w``, la convention de ``SynchronousMachine`` ; sa
pente ``-Pm/w^2`` est un amortissement de ``F_r.Pm`` pu sur chaque etage (0,24 pu
sur l'etage HP du FBM a 0,8 pu de charge), du meme ordre que l'amortissement
mecanique : a ne pas compter deux fois avec un ``D`` mesure.

Les couples de troncon ``C_<troncon>`` sont des sorties de service : ce sont les
grandeurs de la fatigue (:mod:`g2egrid.fatigue`).

La ligne d'arbre se decrit par :class:`LigneArbre` ; :data:`ARBRE_FBM` est celle du
First Benchmark Model (892,4 MVA, 3600 tr/min, 6 masses : HP, IP, LPA, LPB, GEN,
EXC), dont les frequences de torsion a rotor libre, 15,71 / 20,21 / 25,55 / 32,28 /
47,46 Hz, sont retrouvees exactement (``tests/test_turbo_alternateur.py``).
:func:`ligne_deux_masses` construit une ligne turbine-alternateur a partir de
grandeurs physiques (inerties en kg.m2, raideur en N.m/rad), comme les groupes
hydrauliques du PNNL ; :func:`ligne_depuis_physique` fait de meme pour un nombre
quelconque de masses.

GROUPES NUCLEAIRES. Un groupe nucleaire a eau legere tourne a DEMI-VITESSE
(quatre poles : 1500 tr/min a 50 Hz, 1800 a 60 Hz), sans corps moyenne pression,
avec un corps haute pression et plusieurs corps basse pression a double flux :

- :data:`ARBRE_NUCLEAIRE_EPRI` : le groupe nucleaire reel (anonyme) du rapport
  EPRI 1013460 (2006, annexe A), 1500 MVA, 1800 tr/min, 60 Hz, neuf masses (un
  corps HP et trois corps BP, chacun en deux demi-corps, et l'alternateur ;
  excitation statique, donc pas de masse d'excitatrice) ;
- :data:`ARBRE_AM600` : la ligne a trois masses (HP, BP, alternateur) du projet
  AM600, 600 MWe, 1500 tr/min, 50 Hz (Mudau et Field, KNS 2018, tableau 3).

Les machines dont on ne connait que l'inertie totale ``H`` (c'est le cas de
toutes celles des theses de Boudot et de Mazauric) recoivent une de ces lignes
MISE A L'ECHELLE par :func:`ligne_a_l_echelle` : memes deformees modales, memes
frequences de torsion, inertie totale imposee. C'est une hypothese, et elle est
dite comme telle la ou elle est employee. :data:`LIGNES_ARBRE` nomme les lignes
pour les tables (``SM_PARAMS["arbre"]``, voir :mod:`g2egrid.reseau.tables`).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from g2egrid.elements.sync_machine import SynchronousMachine
from g2egrid.noyau.element import InitContext

__all__ = [
    "ARBRE_AM600",
    "ARBRE_FBM",
    "ARBRE_NUCLEAIRE_EPRI",
    "LIGNES_ARBRE",
    "TAU_TORSION",
    "LigneArbre",
    "classe_multi_masses",
    "ligne_a_l_echelle",
    "ligne_depuis_physique",
    "ligne_deux_masses",
    "turbo_alternateur",
    "turbo_alternateur_nomme",
]


@dataclass(frozen=True)
class LigneArbre:
    """Une ligne d'arbre a masses concentrees, en pu de la machine.

    ``masses`` : noms dans l'ordre de la ligne ; ``generateur`` : nom de la masse
    de l'alternateur (son inertie est le ``H`` de la machine, son amortissement
    son ``KD``). ``H``, ``D`` (amortissement propre), ``F`` (fraction du couple
    moteur) : un par masse, celles de l'alternateur ignorees. ``K``, ``Dm``
    (raideur et amortissement mutuel) : un par troncon.
    """

    masses: tuple[str, ...]
    generateur: str
    H: tuple[float, ...]
    K: tuple[float, ...]
    D: tuple[float, ...]
    Dm: tuple[float, ...]
    F: tuple[float, ...]
    nom: str = ""
    #: les troncons qui sont des ARBRES PHYSIQUES (accouplements entre rotors),
    #: quand le modele decoupe un rotor en plusieurs masses : c'est sur eux que
    #: se juge la fatigue. Vide : tous les troncons.
    accouplements: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        n = len(self.masses)
        if self.generateur not in self.masses:
            raise ValueError("la masse de l'alternateur doit etre dans la ligne")
        if not (len(self.H) == len(self.D) == len(self.F) == n):
            raise ValueError("H, D, F : une valeur par masse")
        if not (len(self.K) == len(self.Dm) == n - 1):
            raise ValueError("K, Dm : une valeur par troncon")
        if min(self.K) <= 0.0:
            raise ValueError("raideurs > 0")
        g = self.masses.index(self.generateur)
        if abs(sum(f for i, f in enumerate(self.F) if i != g) - 1.0) > 1e-9:
            raise ValueError("les fractions du couple moteur doivent faire 1")
        inconnus = set(self.accouplements) - set(self.troncons)
        if inconnus:
            raise ValueError(f"accouplements inconnus : {sorted(inconnus)}")

    @property
    def arbres(self) -> tuple[str, ...]:
        """Les troncons a juger en fatigue : :attr:`accouplements`, ou tous."""
        return self.accouplements or self.troncons

    @property
    def H_total(self) -> float:
        """L'inertie de toute la ligne [s], alternateur compris."""
        return float(sum(self.H))

    @property
    def troncons(self) -> tuple[str, ...]:
        return tuple(f"{a}_{b}" for a, b in zip(self.masses[:-1], self.masses[1:]))

    def frequences_libres(self, wb: float, H_gen: float | None = None) -> np.ndarray:
        """Frequences de torsion [Hz] de la ligne seule, alternateur libre (sans
        raideur electrique), amortissement ignore. Le mode rigide (0 Hz) est rendu."""
        n = len(self.masses)
        H = np.array(self.H, dtype=float)
        if H_gen is not None:
            H[self.masses.index(self.generateur)] = H_gen
        Km = np.zeros((n, n))
        for i, k in enumerate(self.K):
            Km[i, i] += k
            Km[i + 1, i + 1] += k
            Km[i, i + 1] -= k
            Km[i + 1, i] -= k
        w2 = np.linalg.eigvals(np.diag(wb / (2.0 * H)) @ Km).real
        return np.sort(np.sqrt(np.abs(w2))) / (2.0 * np.pi)

    def modele_mecanique(
        self, wb: float, H_gen: float | None = None, D_gen: float = 0.0
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Le modele d'etat LINEAIRE de la ligne seule, ``(A, B, C)`` :
        entree le couple electrique ``Ce`` [pu] applique a l'alternateur (couple
        moteur fige), sorties les couples des troncons [pu], etats
        ``(Dw_1..Dw_n, d_1..d_(n-1))``. C'est la fonction ``H_plant`` du PNNL
        (arXiv 2607.14412, eq. 4 : couple d'arbre sur couple electrique) ; avec le
        reseau, la machine ajoute la raideur et l'amortissement electriques
        (``Network.build_state_space`` les contient)."""
        n = len(self.masses)
        g = self.masses.index(self.generateur)
        H = np.array(self.H, dtype=float)
        D = np.array(self.D, dtype=float)
        if H_gen is not None:
            H[g] = H_gen
        D[g] = D_gen
        K = np.array(self.K, dtype=float)
        Dm = np.array(self.Dm, dtype=float)
        m = n - 1
        # couples de troncon : C = K.d + Dm.(Dw_i - Dw_(i+1))
        Cw = np.zeros((m, n))
        for i in range(m):
            Cw[i, i], Cw[i, i + 1] = Dm[i], -Dm[i]
        Cd = np.diag(K)
        # accelerations : 2H.dDw_r/dt = C_(r-1) - C_r - D_r.Dw_r - Ce.[r = g]
        T = np.zeros((n, m))  # couple de troncon recu par chaque masse
        for i in range(m):
            T[i, i], T[i + 1, i] = -1.0, 1.0
        Minv = np.diag(1.0 / (2.0 * H))
        A = np.zeros((n + m, n + m))
        A[:n, :n] = Minv @ (T @ Cw - np.diag(D))
        A[:n, n:] = Minv @ T @ Cd
        for i in range(m):
            A[n + i, i], A[n + i, i + 1] = wb, -wb
        B = np.zeros((n + m, 1))
        B[g, 0] = -1.0 / (2.0 * H[g])
        C = np.hstack([Cw, Cd])
        return A, B, C


#: Amortissement mecanique par defaut : PROPORTIONNEL A L'INERTIE (amortissement
#: de Rayleigh en masse), ``D_r = 4.H_r/tau``. Avec ``2.H.dDw/dt = -D.Dw + ...``,
#: chaque mode de torsion a rotor libre decroit alors en ``exp(-t/tau)``
#: (``sigma = -D/(4.H)``) ; ``tau = 10 s`` est au milieu de la plage de 4 a 30 s
#: que Kundur donne pour les modes de torsion et que retiennent Hossain et al.
#: (arXiv 2605.01173, p. 3). Le modele de reference n'a pas d'amortissement :
#: c'est l'hypothese a reprendre avec les donnees du constructeur. Le reseau
#: ajoute son amortissement (amortisseurs de l'alternateur), different par mode.
TAU_TORSION = 10.0

_H_FBM = (0.092897, 0.155589, 0.858670, 0.884215, 0.868495, 0.0342165)

#: IEEE First Benchmark Model (IEEE SSR Working Group, Trans. PAS-96, 1977) :
#: 892,4 MVA, 60 Hz, 3600 tr/min, six masses, fractions du couple moteur
#: 0,30 / 0,26 / 0,22 / 0,22.
ARBRE_FBM = LigneArbre(
    masses=("HP", "IP", "LPA", "LPB", "GEN", "EXC"),
    generateur="GEN",
    H=_H_FBM,
    K=(19.303, 34.929, 52.038, 70.858, 2.822),
    D=tuple(4.0 * h / TAU_TORSION for h in _H_FBM),
    Dm=(0.0, 0.0, 0.0, 0.0, 0.0),
    F=(0.30, 0.26, 0.22, 0.22, 0.0, 0.0),
    nom="IEEE First Benchmark Model",
)


def ligne_deux_masses(
    J_turbine: float,
    J_generateur: float,
    K: float,
    D_turbine: float,
    Sn_VA: float,
    vitesse_tr_min: float,
    f0: float = 50.0,
    nom: str = "",
) -> tuple[LigneArbre, float]:
    """Une ligne turbine - alternateur a partir de grandeurs PHYSIQUES.

    ``J`` en kg.m2, ``K`` en N.m/rad (mecanique), ``D_turbine`` en N.m.s/rad
    (amortissement visqueux de l'eau sur la roue). Conversion en pu de la machine
    (``Sn``, vitesse ``w_m``, ``p = 60 f0 / vitesse`` paires de poles) ::

        H = J.w_m^2/(2.Sn)     K_pu = K/(p.C_base)     D_pu = D.w_m/C_base

    avec ``C_base = Sn/w_m`` ; le ``1/p`` de ``K_pu`` vient de l'angle en radians
    ELECTRIQUES du modele (``theta_m = theta_e/p``). Rend ``(ligne, H_generateur)``,
    a passer comme ``H`` de la machine."""
    wm = 2.0 * np.pi * vitesse_tr_min / 60.0
    p = 60.0 * f0 / vitesse_tr_min
    cb = Sn_VA / wm
    Ht = J_turbine * wm**2 / (2.0 * Sn_VA)
    Hg = J_generateur * wm**2 / (2.0 * Sn_VA)
    ligne = LigneArbre(
        masses=("TURB", "GEN"),
        generateur="GEN",
        H=(Ht, Hg),
        K=(K / (p * cb),),
        D=(D_turbine * wm / cb, 0.0),
        Dm=(0.0,),
        F=(1.0, 0.0),
        nom=nom,
    )
    return ligne, Hg


def ligne_depuis_physique(
    masses: tuple[str, ...],
    generateur: str,
    J: tuple[float, ...],
    K: tuple[float, ...],
    F: tuple[float, ...],
    Sn_VA: float,
    vitesse_tr_min: float,
    f0: float = 50.0,
    tau: float = TAU_TORSION,
    nom: str = "",
    accouplements: tuple[str, ...] = (),
) -> LigneArbre:
    """Une ligne a ``n`` masses a partir de grandeurs PHYSIQUES : ``J`` [kg.m2]
    par masse, ``K`` [N.m/rad mecanique] par troncon, ``F`` la fraction du
    couple moteur de chaque masse (0 pour l'alternateur).

    Memes conversions que :func:`ligne_deux_masses` ::

        H = J.w_m^2/(2.Sn)      K_pu = K/(p.C_base),  C_base = Sn/w_m

    L'amortissement, que les donnees publiees ne donnent pas, est celui de
    :data:`TAU_TORSION` (proportionnel a l'inertie, ``D = 4.H/tau``)."""
    wm = 2.0 * np.pi * vitesse_tr_min / 60.0
    p = 60.0 * f0 / vitesse_tr_min
    cb = Sn_VA / wm
    H = tuple(float(j) * wm**2 / (2.0 * Sn_VA) for j in J)
    return LigneArbre(
        masses=tuple(masses),
        generateur=generateur,
        H=H,
        K=tuple(float(k) / (p * cb) for k in K),
        D=tuple(4.0 * h / tau for h in H),
        Dm=(0.0,) * (len(masses) - 1),
        F=tuple(F),
        nom=nom,
        accouplements=tuple(accouplements),
    )


def ligne_a_l_echelle(
    ligne: LigneArbre,
    H_total: float,
    rapport_f0: float = 1.0,
    nom: str | None = None,
) -> LigneArbre:
    """La ligne ``ligne`` ramenee a l'inertie totale ``H_total`` [s] d'une
    autre machine, SANS CHANGER ses frequences de torsion ni ses deformees.

    Toutes les inerties sont multipliees par ``k = H_total / ligne.H_total`` ;
    les raideurs et les amortissements le sont aussi, de sorte que chaque
    pulsation ``w^2 ~ wb.K/(2.H)`` est conservee. ``rapport_f0 = f0 d'origine /
    f0 de la machine cible`` (1,2 pour porter une ligne a 60 Hz sur un reseau a
    50 Hz) corrige les raideurs du changement de ``wb``, pour garder les
    frequences en hertz.

    C'EST UNE HYPOTHESE DE TRAVAIL, pas une donnee : elle suppose que la machine
    cible a la meme repartition d'inertie et de raideur que la machine d'origine.
    Elle sert quand on ne connait d'une machine que son ``H`` total."""
    if H_total <= 0.0 or rapport_f0 <= 0.0:
        raise ValueError("H_total et rapport_f0 > 0")
    k = H_total / ligne.H_total
    return LigneArbre(
        masses=ligne.masses,
        generateur=ligne.generateur,
        H=tuple(h * k for h in ligne.H),
        K=tuple(x * k * rapport_f0 for x in ligne.K),
        D=tuple(d * k for d in ligne.D),
        Dm=tuple(d * k for d in ligne.Dm),
        F=ligne.F,
        nom=ligne.nom if nom is None else nom,
        accouplements=ligne.accouplements,
    )


#: Groupe NUCLEAIRE de 1500 MVA, 1800 tr/min, 60 Hz : "the model of an actual
#: nuclear turbine-generator" (anonyme) du rapport EPRI 1013460, *Torsional
#: Interaction Between Electrical Network Phenomena and Turbine-Generator
#: Shafts : Plant Vulnerability* (2006), annexe A. Un corps haute pression et
#: trois corps basse pression a double flux, chacun en DEUX masses (une par
#: flux), puis l'alternateur ; excitation statique. Inerties en kg.m2 et
#: raideurs en N.m/rad telles que publiees ; puissance repartie egalement entre
#: les huit demi-corps (12,5 %), comme dans le rapport, qui previent que le
#: corps HP d'un groupe reel en prend un peu plus. Les ARBRES PHYSIQUES sont
#: les quatre accouplements entre rotors (les "segments 2, 4, 6 et 8" du
#: rapport) ; les quatre autres troncons sont internes a un rotor.
#: Inertie totale 3,75 s ; modes a rotor libre : 8,15 / 15,83 / 22,05 / 23,71 Hz
#: (ligne entiere), 50,6 / 52,5 / 54,5 Hz (les deux flux d'un corps BP l'un
#: contre l'autre) et 133,9 Hz (les deux flux du corps HP).
ARBRE_NUCLEAIRE_EPRI = ligne_depuis_physique(
    masses=("HPA", "HPB", "LP1A", "LP1B", "LP2A", "LP2B", "LP3A", "LP3B", "GEN"),
    generateur="GEN",
    J=(
        7335.88,
        7057.94,
        38515.45,
        38544.93,
        38544.93,
        38574.41,
        38574.41,
        38603.89,
        70832.05,
    ),
    K=(
        2484192000.0,
        235266000.0,
        1883484000.0,
        363679200.0,
        1883484000.0,
        448293600.0,
        1883484000.0,
        603826800.0,
    ),
    F=(0.125,) * 8 + (0.0,),
    Sn_VA=1500e6,
    vitesse_tr_min=1800.0,
    f0=60.0,
    nom="groupe nucleaire 1500 MVA (EPRI 1013460)",
    accouplements=("HPB_LP1A", "LP1B_LP2A", "LP2B_LP3A", "LP3B_GEN"),
)

#: Ligne d'arbre du projet AM600 (cycle vapeur de 600 MWe pour un reacteur a eau
#: pressurisee de taille moyenne, 1800 a 2000 MWt), demi-vitesse, 50 Hz : trois
#: masses, corps HP a simple flux, UN SEUL corps BP a double flux, alternateur.
#: Mudau et Field, "Rotordynamic Analysis of the AM600 Turbine-Generator
#: Shaftline", Trans. Korean Nuclear Society, automne 2018, tableau 3 (inerties
#: 59 042 / 140 625 / 74 864 kg.m2, raideurs 7,29e8 et 1,43e9 N.m/rad).
#:
#: DEUX RESERVES. (1) L'article ne donne ni la puissance apparente ni le partage
#: de puissance : on prend 648,6 MW bruts (Field, Nucl. Eng. Technol. 49, 2017,
#: tableau du bilan thermique) a cos phi 0,9, soit 720 MVA et H = 4,70 s, et
#: 30 % / 70 % du couple moteur pour HP / BP (ordre de grandeur d'un cycle a
#: vapeur saturee). (2) Les donnees du tableau 3 donnent 18,87 et 28,80 Hz, alors
#: que l'article annonce 18,2 et 23,9 Hz pour ce meme modele (18,6 et 24,2 Hz en
#: elements finis) : le premier mode est retrouve a 4 %, le second ne l'est pas.
#: A prendre comme un ORDRE DE GRANDEUR d'une ligne courte et raide.
ARBRE_AM600 = ligne_depuis_physique(
    masses=("HP", "LP", "GEN"),
    generateur="GEN",
    J=(59042.0, 140625.0, 74864.0),
    K=(7.29e8, 1.43e9),
    F=(0.30, 0.70, 0.0),
    Sn_VA=648.6e6 / 0.9,
    vitesse_tr_min=1500.0,
    f0=50.0,
    nom="AM600, 600 MWe (Mudau et Field, KNS 2018)",
)

#: Les lignes d'arbre NOMMEES, pour les tables et le catalogue
#: (``SM_PARAMS["arbre"]``) : nom -> (ligne, frequence du reseau d'origine).
LIGNES_ARBRE: dict[str, tuple[LigneArbre, float]] = {
    "fbm": (ARBRE_FBM, 60.0),
    "nucleaire_epri": (ARBRE_NUCLEAIRE_EPRI, 60.0),
    "am600": (ARBRE_AM600, 50.0),
}


_CLASSES: dict[tuple, type] = {}


def classe_multi_masses(ligne: LigneArbre) -> type:
    """La classe d'element (mise en cache) d'une machine a cette TOPOLOGIE de
    ligne d'arbre : memes masses, meme alternateur. Les valeurs numeriques
    (H, K, D...) sont des parametres de l'instance."""
    cle = (ligne.masses, ligne.generateur)
    if cle in _CLASSES:
        return _CLASSES[cle]
    masses, gen = ligne.masses, ligne.generateur
    g = masses.index(gen)
    autres = tuple(m for m in masses if m != gen)
    troncons = ligne.troncons
    p_masses = tuple(f"{x}_{m}" for m in autres for x in ("H", "D", "F"))
    p_troncons = tuple(f"{x}_{t}" for t in troncons for x in ("K", "Dm"))
    etats = tuple(f"Dw_{m}" for m in autres) + tuple(f"d_{t}" for t in troncons)

    def dw(s, i):
        """Ecart de vitesse de la masse i (l'alternateur : ``Dwr``)."""
        return s.Dwr if i == g else getattr(s, f"Dw_{masses[i]}")

    def equations(s):
        f, gg, h = SynchronousMachine.equations(s)
        Ce = s.phi_d * s.i_q - s.phi_q * s.i_d
        C = []
        for i, t in enumerate(troncons):
            C.append(
                getattr(s, f"K_{t}") * getattr(s, f"d_{t}")
                + getattr(s, f"Dm_{t}") * (dw(s, i) - dw(s, i + 1))
            )
        acc = {}
        for i, m in enumerate(masses):
            couple = (C[i - 1] if i > 0 else 0) - (C[i] if i < len(troncons) else 0)
            if i == g:
                acc[m] = (couple - Ce - s.KD * s.Dwr) / (2 * s.H)
            else:
                Cm = getattr(s, f"F_{m}") * s.Pm / (s.kCm * dw(s, i) + s.w_ref)
                acc[m] = (Cm + couple - getattr(s, f"D_{m}") * dw(s, i)) / (
                    2 * getattr(s, f"H_{m}")
                )
        f[6] = acc[gen]
        f = (
            f
            + [acc[m] for m in autres]
            + [s.wb * (dw(s, i) - dw(s, i + 1)) for i in range(len(troncons))]
        )
        return f, gg, list(h) + C

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = SynchronousMachine.initialize(self, ctx)
        p = self.pdict
        Pm = init["Pm"]
        Ce = init["phi_d"] * init["i_q"] - init["phi_q"] * init["i_d"]
        cumul = 0.0
        for i, t in enumerate(troncons):
            m = masses[i]
            cumul += Ce * -1.0 if i == g else p[f"F_{m}"] * Pm / ctx.w
            init[f"d_{t}"] = cumul / p[f"K_{t}"]
        for m in autres:
            init[f"Dw_{m}"] = 0.0
        return init

    def __init__(self, name, bus, Sn, Sb, wb, ligne_=None, **params):
        lg = ligne if ligne_ is None else ligne_
        if (lg.masses, lg.generateur) != cle:
            raise ValueError(f"{name}: ligne d'arbre d'une autre topologie")
        for i, m in enumerate(lg.masses):
            if m == gen:
                continue
            params.setdefault(f"H_{m}", lg.H[i])
            params.setdefault(f"D_{m}", lg.D[i])
            params.setdefault(f"F_{m}", lg.F[i])
        for i, t in enumerate(lg.troncons):
            params.setdefault(f"K_{t}", lg.K[i])
            params.setdefault(f"Dm_{t}", lg.Dm[i])
        params.setdefault("kCm", 0.0)
        params.setdefault("H", lg.H[g])
        # l'amortissement mecanique est porte par la ligne (D, Dm) : pas de KD
        # global sur l'alternateur, qui amortirait aussi les modes de torsion
        params.setdefault("KD", lg.D[g])
        SynchronousMachine.__init__(self, name, bus, Sn=Sn, Sb=Sb, wb=wb, **params)
        self.ligne = lg

    nom = "MultiMasses_" + "_".join(masses)
    cls = type(
        nom,
        (SynchronousMachine,),
        {
            "__doc__": f"Machine synchrone a ligne d'arbre {', '.join(masses)}.",
            "__module__": __name__,
            "model_name": nom,
            "params": SynchronousMachine.params + ("kCm",) + p_masses + p_troncons,
            "states": SynchronousMachine.states + etats,
            "svc_outputs": SynchronousMachine.svc_outputs
            + tuple(f"C_{t}" for t in troncons),
            "equations": staticmethod(equations),
            "initialize": initialize,
            "__init__": __init__,
            "troncons": troncons,
        },
    )
    _CLASSES[cle] = cls
    return cls


def turbo_alternateur(
    name: str,
    bus: int,
    Sn: float,
    Sb: float,
    wb: float,
    ligne: LigneArbre = ARBRE_FBM,
    **params,
):
    """Une machine synchrone a ligne d'arbre ``ligne`` (FBM par defaut).

    Les parametres electriques et de regulation sont ceux de
    ``SynchronousMachine`` (``SM_DEFAULT``), surchargeables ; ``H`` de
    l'alternateur et les valeurs de la ligne viennent de ``ligne``."""
    return classe_multi_masses(ligne)(name, bus, Sn, Sb, wb, ligne_=ligne, **params)


def turbo_alternateur_nomme(
    name: str,
    bus: int,
    Sn: float,
    Sb: float,
    wb: float,
    arbre: str,
    **params,
):
    """Une machine a la ligne d'arbre NOMMEE ``arbre`` (cle de
    :data:`LIGNES_ARBRE`), mise a l'echelle de SON inertie.

    C'est la porte des tables et du catalogue (``SM_PARAMS["arbre"]`` ou
    ``PARAMS_PAR_BUS[bus]["arbre"]``). Le ``H`` de la machine y garde son sens
    de fiche constructeur -- l'inertie de TOUTE la ligne, turbines comprises -- et
    la ligne nommee le repartit entre ses masses par :func:`ligne_a_l_echelle`,
    frequences de torsion conservees en hertz."""
    from g2egrid.elements.sync_machine import SM_DEFAULT

    cle = str(arbre).lower()
    if cle not in LIGNES_ARBRE:
        raise ValueError(
            f"{name}: ligne d'arbre {arbre!r} inconnue ; {', '.join(LIGNES_ARBRE)}"
        )
    ligne, f0_origine = LIGNES_ARBRE[cle]
    H_total = float(params.pop("H", SM_DEFAULT["H"]))
    f0 = wb / (2.0 * np.pi)
    lg = ligne_a_l_echelle(ligne, H_total, rapport_f0=f0_origine / f0)
    return turbo_alternateur(name, bus, Sn, Sb, wb, ligne=lg, **params)
