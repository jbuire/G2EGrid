"""Simulation temporelle du modele NON LINEAIRE.

Le modele d'etat de :meth:`~g2egrid.reseau.network.Network.build_state_space` est une
*linearisation* : les jacobiens sont evalues une fois au point de
fonctionnement, et les variables algebriques eliminees par complement de
Schur. Ce module prend le meme reseau, le meme cablage
(:class:`~g2egrid.reseau.network.Cablage`) et les memes equations, mais les garde
**non lineaires** et les integre dans le temps.

La difficulte est que le systeme ne se met pas sous forme d'ODE
explicitement. Chaque element porte ses propres equations algebriques
``g_i(x_i, z_i, u_i) = 0``, mais ses entrees sont cablees aux SORTIES des
autres, et une sortie ``h_j`` peut elle-meme dependre de ``u_j`` :

.. math::

    0 = g(x, z, u), \\qquad u = R_{uext}\\,u_{ext} + R_y\\,h(x, z, u)

On ne peut donc pas resoudre un ``z_i`` isolement puis le reporter chez le
voisin : tous les ``z`` et tout le vecteur ``u`` se resolvent **ensemble**,
par Newton, a chaque etat ``x`` que l'integrateur demande. Le jacobien de ce
systeme est analytique, assemble a partir des jacobiens par element et de la
matrice de cablage ``Ry`` ::

        [   dg/dz         dg/du      ]
    J = [                             ]
        [ -Ry.dh/dz   I - Ry.dh/du   ]

C'est la meme combinaison que l'elimination A/B/C/D du modele lineaire,
gardee sous forme de matrice explicite : Newton la redemande a chaque itere,
au lieu de la resoudre une fois pour toutes.

La methode est celle du portage de F. Kelada (``timedomain/emt.py`` de
``FKELADA/G2ELin``, branche ``python-migration``), adaptee au cablage de
``g2egrid``.

Ce que ce module simule reste un modele **dq moyenne** : pas de commutation,
pas de repere abc ni de desequilibre, lignes a constantes localisees. Ce
n'est donc pas un EMT au sens d'un EMTP ; c'est le meme jeu d'equations que
l'analyse modale, resolu sans linearisation.

CE QUI S'ARTICULE AUTOUR
    * :class:`Evenement` et :meth:`ModeleNonLineaire.scenario` decrivent la
      perturbation ; le scenario porte ses instants de discontinuite, dont
      :meth:`ModeleNonLineaire.simuler` fait des bornes d'intervalles ;
    * :func:`mesures` reconstruit le long d'une trajectoire les grandeurs
      qu'un exploitant lit -- module et angle de tension, frequence
      instantanee, P et Q par element -- et non les etats bruts ;
    * :func:`linearisation` rend la matrice ``A`` INSTANTANEE en un point de
      la trajectoire, ce que le modele linearise a l'equilibre ne sait pas
      donner : un mode peut n'etre dangereux que pendant le transitoire.
"""

from __future__ import annotations

import copy
import math
import warnings
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.integrate import solve_ivp
from scipy.optimize import root

from g2egrid.elements.bus import Bus, ResistiveBus
from g2egrid.elements.limiteurs import DUREE_MIN_MODE, IMAX_BLOQUE
from g2egrid.elements.load import Defaut, DefautDC
from g2egrid.elements.repere import Repere
from g2egrid.noyau.integrateur import integrer_bdf2, reinitialiser
from g2egrid.noyau.limites import TOLERANCE as TOLERANCE_GARDE
from g2egrid.reseau.network import Cablage, Network
from g2egrid.noyau.progres import signaler

__all__ = [
    "Courbe",
    "Evenement",
    "Manoeuvre",
    "ModeleNonLineaire",
    "Protection",
    "Trajectoire",
    "linearisation",
    "mesures",
    "protections_ac",
    "protections_dc",
    "reponse_indicielle",
    "simuler_lineaire",
]

#: Nom lisible des deux bouts d'une branche, pour les mesures de transit.
_BOUT = {"f": "from", "t": "to"}


@dataclass
class Trajectoire:
    """Le resultat d'une simulation : trajectoires des etats, et des sorties
    si on les a redemandees."""

    t: np.ndarray
    x: np.ndarray  # (n_t, n_x)
    etats: list[str]
    sorties: dict[str, np.ndarray] = field(default_factory=dict)
    #: nombre de resolutions algebriques (une par evaluation du second membre)
    n_solves: int = 0
    message: str = ""
    #: les morceaux integres chacun avec SON modele, quand une manoeuvre
    #: (ouverture, echelon de charge...) a change le reseau en cours de route :
    #: ``(i0, i1, modele, colonnes)`` -- les instants ``t[i0:i1]``, le modele
    #: qui les a produits, et les colonnes de ``x`` qui portent ses etats. Un
    #: etat d'un element disparu vaut ``NaN`` apres sa disparition.
    segments: list = field(default_factory=list)
    #: les basculements de mode des limiteurs et butees, ``(t, element,
    #: drapeau, valeur)`` -- "GFM5 entre en limitation de courant a 0,102 s"
    basculements: list = field(default_factory=list)
    #: les declenchements de protection, ``(t, protection, [(genre, cible)])``
    #: -- ``t`` est l'instant ou la mesure a franchi le seuil ; les manoeuvres
    #: sont posees ``retard`` plus tard
    declenchements: list = field(default_factory=list)
    #: relais a retombee (:attr:`Protection.retombee`) : ``(t, protection)``
    #: des DEMARRAGES (seuil franchi, temporisation lancee) et des RETOMBEES
    #: (mesure revenue sous le seuil de retombee avant la fin de la
    #: temporisation : pas de declenchement)
    demarrages: list = field(default_factory=list)
    retombees: list = field(default_factory=list)
    #: les ``(z, u)`` resolus par l'integrateur, un par instant (``None`` la
    #: ou ils ne sont pas gardes) -- seulement pour un modele a noeuds
    #: ideaux, qu'on ne peut pas relire en re-resolvant l'algebrique a ``x``
    #: fige (indice 2). Vide sinon.
    algebriques: list = field(default_factory=list)

    def etat(self, nom: str) -> np.ndarray:
        return self.x[:, self.etats.index(nom)]

    def __repr__(self) -> str:
        return (
            f"<Trajectoire {len(self.t)} points sur [{self.t[0]:g}, {self.t[-1]:g}] s, "
            f"{self.x.shape[1]} etats, {self.n_solves} resolutions algebriques>"
        )


# ----------------------------------------------------------------------
# perturbations
# ----------------------------------------------------------------------
@dataclass
class Evenement:
    """Une perturbation datee sur une entree externe.

    ``genre`` vaut ``"echelon"`` (saut de ``amplitude`` a ``instant``) ou
    ``"rampe"`` (montee lineaire de 0 a ``amplitude`` entre ``instant`` et
    ``instant + duree``, puis palier). ``amplitude`` s'AJOUTE a la consigne
    d'equilibre : un evenement n'ecrase jamais le point de fonctionnement, il
    l'ecarte.

    La rampe n'est pas cosmetique. Un echelon excite tous les modes du
    modele, jusqu'aux 1e7 rad/s des troncons de cable, avec la meme energie ;
    une rampe d'un dixieme de seconde ne sollicite que ce qui est lent devant
    elle. Quand on regarde un mode electromecanique, c'est la rampe qui donne
    une trajectoire lisible -- et c'est aussi ce que fait un automate de
    consigne reel, qui ne saute pas.
    """

    instant: float
    entree: str
    amplitude: float
    genre: str = "echelon"
    duree: float = 0.0

    def __post_init__(self) -> None:
        if self.genre not in ("echelon", "rampe"):
            raise ValueError(f"genre {self.genre!r} inconnu ('echelon' ou 'rampe')")
        if self.genre == "rampe" and self.duree <= 0.0:
            raise ValueError(f"{self.entree}: une rampe demande duree > 0")

    def discontinuites(self) -> tuple[float, ...]:
        """Les instants ou la consigne n'est pas derivable.

        Les deux bornes d'une rampe en font partie : la consigne y est
        continue, mais sa PENTE saute, et un integrateur a pas et a ordre
        variables se fait reprendre de la meme facon sur une cassure de pente
        que sur un saut.
        """
        if self.genre == "rampe":
            return (float(self.instant), float(self.instant + self.duree))
        return (float(self.instant),)

    def valeur(self, t: float) -> float:
        """L'ecart de consigne apporte par cet evenement a l'instant ``t``."""
        if t < self.instant:
            return 0.0
        if self.genre == "echelon" or t >= self.instant + self.duree:
            return self.amplitude
        return self.amplitude * (t - self.instant) / self.duree


@dataclass
class Courbe:
    """Une entree externe pilotee par une COURBE tabulee ``(t, valeurs)``.

    C'est l'entree d'une simulation de plusieurs minutes : la courbe de charge
    d'un centre de donnees, un profil de production, une consigne de puissance
    mesuree. Contrairement a un :class:`Evenement`, qui ECARTE la consigne
    d'equilibre, la courbe en donne la valeur ABSOLUE (``absolue=True``, le
    defaut) : ``valeurs`` est directement la puissance demandee. Avec
    ``absolue=False``, elle s'ajoute a l'equilibre comme un evenement.

    ``interpolation`` :

    * ``"pchip"`` (defaut) -- interpolation cubique monotone par morceaux,
      derivable une fois, sans depassement entre deux points. Elle ne declare
      AUCUNE discontinuite : l'integrateur n'est pas redemarre aux points de la
      courbe, ce qui compte pour une courbe de plusieurs centaines de points ;
    * ``"lineaire"`` -- segments de droite ; chaque point est une cassure de
      pente, donc un instant de decoupage de l'integration (cher si la courbe
      est dense) ;
    * ``"palier"`` -- la valeur reste constante jusqu'au point suivant ; chaque
      point est un saut.

    Avant le premier point la courbe garde sa premiere valeur, apres le dernier
    sa derniere. Pour que la simulation parte de l'equilibre, la valeur en
    ``t[0]`` doit etre celle que l'element a retenue a l'equilibre ;
    :meth:`ModeleNonLineaire.scenario` le verifie.
    """

    entree: str
    t: np.ndarray
    valeurs: np.ndarray
    interpolation: str = "pchip"
    absolue: bool = True

    def __post_init__(self) -> None:
        from scipy.interpolate import PchipInterpolator

        self.t = np.asarray(self.t, dtype=float)
        self.valeurs = np.asarray(self.valeurs, dtype=float)
        if self.t.ndim != 1 or self.t.shape != self.valeurs.shape or len(self.t) < 2:
            raise ValueError(
                f"{self.entree}: t et valeurs, deux vecteurs de >= 2 points"
            )
        if np.any(np.diff(self.t) <= 0.0):
            raise ValueError(f"{self.entree}: les instants doivent croitre strictement")
        if not np.all(np.isfinite(self.valeurs)):
            raise ValueError(f"{self.entree}: valeurs non finies")
        if self.interpolation not in ("pchip", "lineaire", "palier"):
            raise ValueError(
                f"interpolation {self.interpolation!r} inconnue "
                "('pchip', 'lineaire' ou 'palier')"
            )
        self._pchip = (
            PchipInterpolator(self.t, self.valeurs)
            if self.interpolation == "pchip"
            else None
        )

    def discontinuites(self) -> tuple[float, ...]:
        if self.interpolation == "pchip":
            return ()
        return tuple(float(x) for x in self.t)

    def valeur(self, t: float) -> float:
        """La valeur de la courbe a l'instant ``t`` (absolue ou ecart)."""
        if t <= self.t[0]:
            return float(self.valeurs[0])
        if t >= self.t[-1]:
            return float(self.valeurs[-1])
        if self._pchip is not None:
            return float(self._pchip(t))
        if self.interpolation == "palier":
            return float(self.valeurs[np.searchsorted(self.t, t, side="right") - 1])
        return float(np.interp(t, self.t, self.valeurs))


#: Les manoeuvres qui changent le RESEAU lui-meme, et non une consigne.
GENRES_MANOEUVRE = (
    "ouverture",
    "charge",
    "saut_de_phase",
    "defaut",
    "elimination",
    "defaut_dc",
    "elimination_dc",
    "blocage",
)
#: impedance de defaut par defaut (pu systeme) : alternatif, continu
IMPEDANCE_DEFAUT = {"defaut": (1e-3, 1e-2), "defaut_dc": (1e-4, 1e-3)}


@dataclass
class Manoeuvre:
    """Une manoeuvre datee sur le reseau : ce qu'un :class:`Evenement` ne peut
    pas decrire, parce qu'elle change le modele et pas une consigne.

    - ``"ouverture"`` : ``cible`` est le nom d'un element assemble (``L3_5_7``,
      ``SM2``, ``Load5``, ``Xsh4_7``...) ou d'un condensateur (``Csh2_5``, cf.
      :func:`g2egrid.reseau.disjoncteurs.nom_element`). L'element disparait du
      modele ; une ligne reprend a ses deux noeuds la demi-charge capacitive
      qu'elle leur avait versee. Declencher la machine de reference passe le
      repere a un :class:`~g2egrid.elements.Repere` qui suit la plus grosse
      machine restante : il n'y a pas de slack en dynamique.
    - ``"charge"`` : la charge ``cible`` varie de ``dp`` (P) et ``dq`` (Q),
      en valeur relative (``dp = 0.1`` : +10 %). Une charge a impedance
      constante change d'impedance a la tension du point de fonctionnement --
      c'est l'enclenchement d'un depart, pas une consigne de puissance.
    - ``"saut_de_phase"`` : la tension du bus ``cible`` (``"bus5"`` ou ``5``)
      tourne instantanement de ``angle`` degres. Pour un saut de phase du
      reseau amont, un bus infini a deja une entree ``delta_ref`` : c'est un
      :class:`Evenement`.
    - ``"defaut"`` : un court-circuit a la terre au bus ``cible``,
      d'impedance ``R + jX`` (pu systeme) -- un element
      :class:`~g2egrid.elements.Defaut` pose le temps du defaut.
      ``"elimination"`` le retire ; l'elimination par ouverture d'une ligne
      s'ecrit avec une ``"ouverture"`` au meme instant.
    - ``"defaut_dc"`` : un court-circuit POLE A POLE au noeud continu
      ``cible`` (``"bus213"`` ou ``213``), branche R-L
      (:class:`~g2egrid.elements.DefautDC`), ``R = 1e-4``, ``X = 1e-3`` pu par
      defaut ; ``"elimination_dc"`` le retire. Ouvrir un cable continu --
      l'action d'un disjoncteur a courant continu -- est une ``"ouverture"``
      du cable (``CDC7_213_214``).
    - ``"blocage"`` : le convertisseur ``cible`` (``GFL113``, ``GFM121``) se
      BLOQUE -- ses interrupteurs s'ouvrent, il ne debite plus ni cote
      alternatif ni cote continu. Represente par sa propre limitation de
      courant, ramenee a ``IMAX_BLOQUE`` (presque zero) et figee : la boucle
      de courant eteint le courant en quelques millisecondes, le filtre reste
      en place. C'est un convertisseur a capacite de blocage (pont complet),
      pas les diodes de roue libre d'un demi-pont, qui continueraient
      d'alimenter un defaut continu depuis l'alternatif. Le retirer d'un coup
      (``"ouverture"``) couperait en zero seconde le courant inductif de son
      transformateur, et la surtension effondrerait le poste.

    RIEN N'EST RECALCULE PAR UN LOAD FLOW : les elements gardent leurs
    consignes, leurs etats sont repris par leur nom, les variables algebriques
    se resolvent de nouveau a l'instant de la manoeuvre. C'est ce qui separe
    une manoeuvre en simulation d'un disjoncteur ouvert en regime permanent
    (:mod:`g2egrid.reseau.disjoncteurs`), ou le reseau est re-equilibre.
    """

    instant: float
    genre: str
    cible: str
    dp: float = 0.0
    dq: float = 0.0
    angle: float = 0.0  # degres
    R: float | None = None  # impedance de defaut, pu systeme (IMPEDANCE_DEFAUT)
    X: float | None = None

    def __post_init__(self) -> None:
        if self.genre not in GENRES_MANOEUVRE:
            raise ValueError(
                f"manoeuvre {self.genre!r} inconnue ; {', '.join(GENRES_MANOEUVRE)}"
            )
        self.cible = str(self.cible)
        r, x = IMPEDANCE_DEFAUT.get(self.genre, (1e-3, 1e-2))
        self.R = r if self.R is None else float(self.R)
        self.X = x if self.X is None else float(self.X)

    def discontinuites(self) -> tuple[float, ...]:
        return (float(self.instant),)


#: Les protections que :meth:`ModeleNonLineaire.simuler` sait surveiller.
GENRES_PROTECTION = (
    "surintensite_dc",
    "sous_tension_dc",
    "surintensite_convertisseur",
    "surtension_ac",
)


@dataclass
class Protection:
    """Une protection surveillee PENDANT la simulation : une mesure, un seuil,
    un retard, et les manoeuvres qu'elle declenche.

    - ``"surintensite_dc"`` : ``|i|`` du cable continu ``cible``
      (``CDC7_213_214``, pu systeme) depasse ``seuil`` -> le disjoncteur a
      courant continu l'OUVRE (``"ouverture"`` du cable) ;
    - ``"sous_tension_dc"`` : la tension du noeud continu ``cible``
      (``bus213``) passe sous ``seuil`` (pu) -> les convertisseurs du noeud se
      BLOQUENT (``"blocage"``) ;
    - ``"surintensite_convertisseur"`` : ``|is|`` du convertisseur ``cible``
      (pu de l'equipement) depasse ``seuil`` -> il se bloque ;
    - ``"surtension_ac"`` : la tension ``|ve|`` aux bornes du convertisseur
      ``cible`` depasse ``seuil`` (pu) -> il se bloque. C'est la protection
      qui fait decrocher un parc eolien dont l'ile n'a plus de formeur.

    Le seuil est une GARDE de l'integrateur, comme celles des limiteurs :
    l'integration s'arrete a l'instant exact ou la mesure le franchit. Les
    manoeuvres sont alors posees ``retard`` secondes plus tard -- le temps de
    detection et d'ouverture (quelques millisecondes pour un disjoncteur
    continu hybride, une fraction pour le blocage d'un convertisseur). Une
    protection ne declenche qu'une fois ; ``actions`` (``[(genre, cible)]``)
    remplace les manoeuvres par defaut.

    RETOMBEE -- par defaut (``retombee = None``) la protection est VERROUILLEE
    des que le seuil est franchi : elle declenche ``retard`` plus tard, meme
    si la mesure est redescendue entre-temps. Avec ``retombee`` (le rapport
    de retombee d'un relais a temps constant, 0,95 a 0,98 pour un relais a
    maximum, plus de 1 pour un relais a minimum), le franchissement ne fait
    que DEMARRER la temporisation : si la mesure repasse le seuil de retombee
    ``retombee.seuil`` avant son terme, le relais retombe et se reamorce ; il
    ne declenche que si le defaut a TENU toute la temporisation. C'est ce
    qu'il faut pour un relais de surtension temporise : un pic transitoire
    apres une manoeuvre ne doit pas le faire declencher.
    """

    genre: str
    cible: str
    seuil: float
    retard: float = 0.0
    actions: list | None = None
    retombee: float | None = None

    def __post_init__(self) -> None:
        if self.genre not in GENRES_PROTECTION:
            raise ValueError(
                f"protection {self.genre!r} inconnue ; {', '.join(GENRES_PROTECTION)}"
            )
        if self.retard < 0:
            raise ValueError(f"protection {self.cible} : retard negatif")
        if self.retombee is not None:
            minimum = self.genre == "sous_tension_dc"
            if self.retombee <= 0 or (
                self.retombee < 1 if minimum else self.retombee > 1
            ):
                raise ValueError(
                    f"protection {self.cible} : rapport de retombee {self.retombee} "
                    f"({'>= 1 pour un relais a minimum' if minimum else 'dans ]0, 1] pour un relais a maximum'})"
                )
        self.cible = str(self.cible)

    @property
    def nom(self) -> str:
        return f"{self.genre}:{self.cible}"

    def _etats(self) -> tuple[str, ...]:
        if self.genre == "surintensite_dc":
            return (f"{self.cible}.i",)
        if self.genre == "sous_tension_dc":
            c = self.cible if self.cible.startswith("bus") else f"bus{self.cible}"
            return (f"{c}.v",)
        if self.genre == "surtension_ac":
            return (f"{self.cible}.ved", f"{self.cible}.veq")
        return (f"{self.cible}.isd", f"{self.cible}.isq")

    def indices(self, m: "ModeleNonLineaire") -> list[int] | None:
        """Les indices des etats mesures dans le modele ``m`` ; ``None`` si
        l'equipement n'y est plus (ouvert, bloque)."""
        noms = {n: i for i, n in enumerate(m.c.state_names)}
        idx = [noms.get(e) for e in self._etats()]
        return None if any(i is None for i in idx) else idx

    def garde(self, x: np.ndarray, idx: list[int], seuil: float | None = None) -> float:
        """Positive quand la protection doit declencher (``seuil`` remplace
        celui de la protection : c'est la garde de retombee)."""
        s = self.seuil if seuil is None else seuil
        if self.genre == "sous_tension_dc":
            return s - float(x[idx[0]])
        return float(np.sqrt(sum(float(x[i]) ** 2 for i in idx))) - s

    def garde_retombee(self, x: np.ndarray, idx: list[int]) -> float:
        """Pendant la temporisation : NEGATIVE quand le relais retombe."""
        return self.garde(x, idx, self.retombee * self.seuil)

    def manoeuvres(self, t: float, m: "ModeleNonLineaire") -> list[Manoeuvre]:
        """Les manoeuvres declenchees, datees de ``t`` (deja retardees)."""
        if self.actions:
            return [Manoeuvre(t, g, c) for g, c in self.actions]
        if self.genre == "surintensite_dc":
            return [Manoeuvre(t, "ouverture", self.cible)]
        if self.genre in ("surintensite_convertisseur", "surtension_ac"):
            return [Manoeuvre(t, "blocage", self.cible)]
        b = int(self.cible[3:] if self.cible.startswith("bus") else self.cible)
        return [
            Manoeuvre(t, "blocage", el.name)
            for el in m.net.elements
            if el.ports and el.ports_dc and int(el.buses[1]) == b
        ]


def protections_dc(
    net,
    k_courant: float = 3.0,
    v_min: float = 0.8,
    i_convertisseur: float = 2.0,
    v_max_ac: float = 1.3,
    retard_disjoncteur: float = 3e-3,
    retard_blocage: float = 2e-4,
    retard_surtension: float = 1e-3,
) -> list[Protection]:
    """Un plan de protection type du reseau continu de ``net`` (equilibre
    resolu) :

    - chaque cable continu a son disjoncteur, qui s'ouvre ``retard_disjoncteur``
      apres que ``|i|`` a depasse ``k_courant`` fois son courant d'equilibre --
      ou le courant MOYEN des cables, pour un cable presque a vide, qu'un
      simple report de charge ne doit pas faire declencher ;
    - chaque convertisseur a port continu se bloque ``retard_blocage`` apres
      que la tension de son noeud continu est passee sous ``v_min``, ou que
      son courant a depasse ``i_convertisseur`` pu de sa puissance ;
    - chaque source SANS port continu (parc eolien ou photovoltaique, GFM a
      batterie) se bloque ``retard_surtension`` apres que la tension a ses
      bornes a depasse ``v_max_ac`` : un parc dont l'ile a perdu son formeur
      decroche. Les stations HVDC, elles, tiennent le creux et le rebond qui
      suivent un blocage voisin.
    """
    from g2egrid.elements.branch import LineDC

    cables = [el for el in net.elements if isinstance(el, LineDC)]
    i0 = {
        el.name: abs(float(el.x0[0])) if el.x0 is not None and len(el.x0) else 0.0
        for el in cables
    }
    plancher = float(np.mean(list(i0.values()))) if i0 else 0.0
    out = [
        Protection(
            "surintensite_dc",
            el.name,
            k_courant * max(i0[el.name], plancher, 1e-3),
            retard_disjoncteur,
        )
        for el in cables
    ]
    convs = [el for el in net.elements if el.ports and el.ports_dc]
    noeuds = sorted({int(el.buses[1]) for el in convs})
    out += [
        Protection("sous_tension_dc", f"bus{b}", v_min, retard_blocage) for b in noeuds
    ]
    out += [
        Protection(
            "surintensite_convertisseur", el.name, i_convertisseur, retard_blocage
        )
        for el in convs
    ]
    out += [
        Protection("surtension_ac", el.name, v_max_ac, retard_surtension)
        for el in net.elements
        if "mode_lim" in el.params and "ved" in el.states and not el.ports_dc
    ]
    return out


#: familles de convertisseurs visees par :func:`protections_ac`
FAMILLES_AC = ("GFL", "GFM")


def protections_ac(
    net,
    v_max: float = 1.10,
    retard: float = 0.2,
    retombee: float | None = 0.98,
    action: str = "ouverture",
    familles: Sequence[str] = ("GFL",),
    cibles: Sequence[str] | None = None,
) -> list[Protection]:
    """Le relais de SURTENSION de chaque convertisseur sans port continu :
    il declenche si la tension a ses bornes (``|v_e|``, pu) depasse ``v_max``
    pendant ``retard`` secondes d'affilee (relais a retombee, rapport
    ``retombee`` ; ``None`` : relais verrouille au premier franchissement).

    C'est la protection de decouplage d'un producteur raccorde en
    distribution (seuil de 1,10 pu en regime permanent, quelques centaines de
    millisecondes) : elle sert a etudier les CASCADES de surtension, ou le
    declenchement d'un producteur qui absorbait du reactif fait monter la
    tension de ses voisins jusqu'a leur propre seuil.

    ``action`` : ``"ouverture"`` (le disjoncteur du producteur s'ouvre, filtre
    compris -- le cas d'un onduleur de distribution) ou ``"blocage"`` (les
    interrupteurs se bloquent, le filtre reste en place). ``familles`` :
    ``"GFL"``, ``"GFM"`` ; ``cibles`` restreint a des equipements nommes.
    """
    from g2egrid.elements.gfl import _GFLBase
    from g2egrid.elements.gfm import _GFMBase

    classes = {"GFL": _GFLBase, "GFM": _GFMBase}
    inconnues = [f for f in familles if f not in classes]
    if inconnues:
        raise ValueError(f"famille(s) inconnue(s) : {inconnues} ; {FAMILLES_AC}")
    if action not in ("ouverture", "blocage"):
        raise ValueError(f"action {action!r} : 'ouverture' ou 'blocage'")
    voulues = tuple(classes[f] for f in familles)
    out = []
    for el in net.elements:
        if not isinstance(el, voulues) or el.ports_dc or "ved" not in el.states:
            continue
        if cibles is not None and el.name not in cibles:
            continue
        out.append(
            Protection(
                "surtension_ac",
                el.name,
                v_max,
                retard,
                actions=[(action, el.name)],
                retombee=retombee,
            )
        )
    return out


def _copie_element(el):
    """Copie d'un element dont on peut changer les parametres sans toucher
    l'original : le modele d'avant la manoeuvre doit rester intact, il sert a
    relire la premiere partie de la trajectoire."""
    e = copy.copy(el)
    e.p = np.array(el.p, dtype=float)
    e.pdict = dict(el.pdict)
    return e


class ModeleNonLineaire:
    """Le reseau, garde non lineaire, integrable dans le temps.

    Construit a partir d'un ``Network`` deja equilibre (``run_loadflow`` puis
    ``solve_equilibrium``) ::

        m = ModeleNonLineaire(net)
        tr = m.simuler((0, 1.0), uext_fn=m.echelon("SM1.p_ref", 0.05, t0=0.1))
    """

    def __init__(self, net: Network):
        from g2egrid.elements.bus import IdealBus

        #: les noeuds ideaux (``noeuds_sans_capacite="ideal"``) : leur loi des
        #: noeuds peut ne porter que sur des etats (courants en serie), une
        #: contrainte d'indice 2. Le systeme algebrique n'a alors pas de
        #: solution a ``x`` fige : seul BDF2, qui resout etats et contraintes
        #: ENSEMBLE, l'integre (FR-005 (c)).
        self.ideaux = [el.name for el in net.all_elements if isinstance(el, IdealBus)]
        if any(el.x0 is None for el in net.all_elements):
            net.solve_equilibrium()
        self.net = net
        self.c: Cablage = net.cablage()
        self.p = [el.p for el in self.c.els]
        self._n_solves = 0
        # le cablage est une matrice d'incidence : creux, et fixe
        self._Ry = sp.csr_matrix(self.c.Ry)
        # la derniere factorisation du jacobien algebrique (Newton modifie)
        self._lu = None

    # ------------------------------------------------------------------ noms
    @property
    def etats(self) -> list[str]:
        return self.c.state_names

    @property
    def entrees(self) -> list[str]:
        return self.c.uext_names

    @property
    def nom_sorties(self) -> list[str]:
        return self.c.output_names

    @property
    def noms_z(self) -> list[str]:
        """Les variables algebriques, nommees ``element.variable``."""
        return [
            f"{el.name}.{a}"
            for el, cm in zip(self.c.els, self.c.cms)
            for a in cm.spec.algebs
        ]

    @property
    def noms_u(self) -> list[str]:
        """Les entrees de tous les elements, nommees ``element.entree``."""
        return [
            f"{el.name}.{i}"
            for el, cm in zip(self.c.els, self.c.cms)
            for i in cm.spec.inputs
        ]

    # ------------------------------------------------------- decoupage
    def _par_element(self, v: np.ndarray, offsets: list[int], tailles) -> list:
        return [v[o : o + n] for o, n in zip(offsets, tailles)]

    def _x(self, x):
        return self._par_element(x, self.c.ox, [cm.n_x for cm in self.c.cms])

    def _z(self, z):
        return self._par_element(z, self.c.oz, [cm.n_z for cm in self.c.cms])

    def _u(self, u):
        return self._par_element(u, self.c.ou, [cm.n_u for cm in self.c.cms])

    # ------------------------------------------------------- systeme couple
    def residu(self, zu: np.ndarray, xs: list, uext: np.ndarray) -> np.ndarray:
        """``[g(x, z, u) ; u - Ruext.uext - Ry.h(x, z, u)]``, nul a la solution."""
        c = self.c
        z, u = zu[: c.nz], zu[c.nz :]
        zs, us = self._z(z), self._u(u)
        g, y = [], []
        for cm, pi, xi, zi, ui in zip(c.cms, self.p, xs, zs, us):
            if cm.n_z:
                g.append(cm.g(xi, zi, ui, pi))
            y.append(cm.h(xi, zi, ui, pi))
        g = np.concatenate(g) if g else np.zeros(0)
        y = np.concatenate(y) if y else np.zeros(0)
        return np.concatenate([g, u - (c.Ruext @ uext + c.Ry @ y)])

    def jacobien_creux(self, zu: np.ndarray, xs: list, uext: np.ndarray):
        """Le jacobien analytique du residu, en ``(z, u)``, assemble CREUX.

        ``x`` ne varie pas pendant cette resolution : ``gx`` et ``hx`` n'y
        entrent donc pas, seuls ``gz``, ``gu``, ``hz`` et ``hu``.

        POURQUOI CREUX -- sur wscc9 le systeme couple a 144 inconnues (18
        variables algebriques, 126 signaux de cablage) et son jacobien n'est
        plein qu'a 1,1 % : chaque element ne voit que ses propres variables, et
        le cablage ``Ry`` n'est qu'une matrice d'incidence. L'assembler en
        matrice pleine et le confier a une factorisation dense coutait, a lui
        seul, plus de la moitie d'une simulation EMT -- bien plus que
        l'evaluation des equations elles-memes.

        STRUCTURE FIXE -- construire a chaque appel des matrices
        ``scipy.sparse`` (verification des index, ``hstack``, ``vstack``)
        coutait plus que les equations. Or ce qui change d'un appel a l'autre
        n'est que la VALEUR des blocs ``gz, gu, hz, hu`` des elements ; leur
        place, et la facon dont ``Ry`` les recombine, sont fixees par le
        cablage. :meth:`_structure` les precalcule une fois : un appel ne fait
        plus qu'evaluer les blocs, un produit creux fixe, et poser les valeurs
        dans une matrice CSC de structure connue.
        """
        st = self._structure()
        c = self.c
        z, u = zu[: c.nz], zu[c.nz :]
        zs, us = self._z(z), self._u(u)
        vals = []
        for cm, pi, xi, zi, ui in zip(c.cms, self.p, xs, zs, us):
            if cm.n_z:
                vals.append(np.ravel(cm.gz(xi, zi, ui, pi)))
                vals.append(np.ravel(cm.gu(xi, zi, ui, pi)))
                vals.append(np.ravel(cm.hz(xi, zi, ui, pi)))
            vals.append(np.ravel(cm.hu(xi, zi, ui, pi)))
        v = np.concatenate(vals) if vals else np.zeros(0)
        data = (st["M"] @ v + st["const"])[st["perm"]]
        J = sp.csc_matrix(st["forme"])
        J.data, J.indices, J.indptr = data, st["indices"], st["indptr"]
        return J

    def _structure(self) -> dict:
        """La structure du jacobien, calculee une fois (voir :meth:`jacobien_creux`).

        On numerote chaque coefficient des blocs d'elements dans l'ordre ou
        :meth:`jacobien_creux` les concatene (``gz, gu, hz`` puis ``hu``,
        element par element, blocs pleins aplatis ligne a ligne). Chaque
        coefficient va a une ou plusieurs places du jacobien :

        - ``gz``, ``gu`` : une seule, en haut ;
        - ``hz``, ``hu`` : autant que d'entrees ``p`` cablees sur la sortie
          ``y``, avec le poids ``-Ry[p, y]`` -- c'est ``-Ry.hz`` et
          ``-Ry.hu`` ;

        plus l'identite du bloc ``u``, constante. ``M`` porte ces poids :
        ``valeurs = M.v + const``, dans l'ordre des places uniques, que
        ``perm`` range enfin dans l'ordre CSC.
        """
        if getattr(self, "_st", None) is not None:
            return self._st
        c = self.c
        Ry = sp.csc_matrix(c.Ry)
        lignes, colonnes, sources, poids = [], [], [], []
        n = 0

        def haut(r0, c0, nr, nc):
            nonlocal n
            i, j = np.divmod(np.arange(nr * nc), nc)
            lignes.append(r0 + i)
            colonnes.append(c0 + j)
            sources.append(n + np.arange(nr * nc))
            poids.append(np.ones(nr * nc))
            n += nr * nc

        def bas(y0, c0, nr, nc):
            nonlocal n
            for i in range(nr):
                col = Ry[:, y0 + i]  # les entrees cablees sur la sortie y0 + i
                for pp, w in zip(col.indices, col.data):
                    lignes.append(np.full(nc, c.nz + pp))
                    colonnes.append(c0 + np.arange(nc))
                    sources.append(n + i * nc + np.arange(nc))
                    poids.append(np.full(nc, -w))
            n += nr * nc

        for k, cm in enumerate(c.cms):
            sz, su, sy = c.oz[k], c.ou[k], c.oy[k]
            if cm.n_z:
                haut(sz, sz, cm.n_z, cm.n_z)
                haut(sz, c.nz + su, cm.n_z, cm.n_u)
                bas(sy, sz, cm.n_y, cm.n_z)
            bas(sy, c.nz + su, cm.n_y, cm.n_u)

        N = c.nz + c.nu
        cat = lambda L, t: np.concatenate(L).astype(t) if L else np.zeros(0, t)  # noqa: E731
        li, co = cat(lignes, np.int64), cat(colonnes, np.int64)
        so, po = cat(sources, np.int64), cat(poids, float)
        # l'identite du bloc u : des places qui existent toujours
        diag = c.nz + np.arange(c.nu)
        cle_places = np.concatenate([li * N + co, diag * N + diag])
        places, inv = np.unique(cle_places, return_inverse=True)
        M = sp.csr_matrix((po, (inv[: len(li)], so)), shape=(len(places), n))
        const = np.zeros(len(places))
        np.add.at(const, inv[len(li) :], 1.0)
        # l'ordre CSC des places : par colonne, puis par ligne
        r, cc = np.divmod(places, N)
        perm = np.lexsort((r, cc))
        indptr = np.concatenate([[0], np.cumsum(np.bincount(cc, minlength=N))])
        self._st = {
            "M": M,
            "const": const,
            "perm": perm,
            "indices": r[perm].astype(np.int32),
            "indptr": indptr.astype(np.int32),
            "forme": (N, N),
        }
        return self._st

    def jacobien(self, zu: np.ndarray, xs: list, uext: np.ndarray) -> np.ndarray:
        """Le meme, en matrice pleine -- pour les diagnostics et les tests."""
        return self.jacobien_creux(zu, xs, uext).toarray()

    def resoudre_algebrique(
        self,
        x: np.ndarray,
        uext: np.ndarray,
        z_essai: np.ndarray | None = None,
        u_essai: np.ndarray | None = None,
        t: float | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Newton sur le systeme couple, a ``x`` et ``uext`` figes.

        ``hybr`` (dogleg de Powell, ``hybrj`` de MINPACK des qu'on lui donne
        un jacobien) d'abord ; repli sur Levenberg-Marquardt, plus lent mais
        plus robuste loin de la solution, s'il n'y arrive pas.

        ``t`` ne sert qu'au diagnostic d'echec -- c'est la premiere chose
        qu'on veut savoir d'une simulation qui casse.
        """
        c = self.c
        if z_essai is None:
            z_essai = c.z0()
        if u_essai is None:
            u_essai = c.u0()
        xs = self._x(x)
        zu0 = np.concatenate([z_essai, u_essai])

        # D'ABORD UN NEWTON CREUX, MODIFIE. Partant de la solution du pas
        # precedent, il converge en quelques iterations, chacune ne coutant
        # qu'une resolution par une factorisation LU creuse -- et cette
        # factorisation est GARDEE d'une resolution a l'autre (methode des
        # cordes) : d'un pas d'integration au suivant, le jacobien du systeme
        # algebrique bouge a peine, et l'evaluer coutait autant que tout le
        # reste. On ne le refait que si la convergence ralentit (contraction
        # du residu moins bonne que 1/2). MINPACK (``hybrj``) reste le repli :
        # plus robuste loin de la solution, mais en matrice PLEINE.
        zu = zu0.copy()
        r = self.residu(zu, xs, uext)
        lu, frais = self._lu, False
        for _ in range(25):
            nr = np.abs(r).max()
            # Seul un residu au niveau de l'arrondi dispense d'un pas : un
            # residu de courant de 1e-13 devient 3e-8 de dv/dt sur un noeud
            # de C = 1e-3 (wb/C = 3e5). On converge donc sur le PAS.
            if nr < 1e-14:
                self._lu = lu
                self._n_solves += 1
                return zu[: c.nz], zu[c.nz :]
            if lu is None:
                try:
                    lu = spla.splu(self.jacobien_creux(zu, xs, uext))
                except RuntimeError:  # singulier : on laisse MINPACK faire
                    break
                frais = True
            dzu = lu.solve(r)
            if not np.all(np.isfinite(dzu)):
                if frais:
                    break
                lu = None
                continue
            zu_n = zu - dzu
            r_n = self.residu(zu_n, xs, uext)
            nr_n = np.abs(r_n).max()
            if not nr_n <= 0.5 * nr and not frais:
                # la corde ne tient plus : jacobien au point courant, et on
                # ne garde le pas que s'il a fait baisser le residu
                lu = None
                if nr_n < nr:
                    zu, r = zu_n, r_n
                continue
            zu, r = zu_n, r_n
            echelle = 1.0 + np.abs(zu).max()
            if np.abs(dzu).max() < 1e-12 * echelle and nr_n <= 1e-8 * echelle:
                self._lu = lu
                self._n_solves += 1
                return zu[: c.nz], zu[c.nz :]
            frais = False
        self._lu = None

        for methode in ("hybr", "lm"):
            sol = root(
                self.residu, zu0, args=(xs, uext), jac=self.jacobien, method=methode
            )
            self._n_solves += 1
            if sol.success and self._resolu(sol.x, xs, uext):
                return sol.x[: c.nz], sol.x[c.nz :]
        raise RuntimeError(self._diagnostic(sol, xs, uext, t))

    def _resolu(self, zu: np.ndarray, xs: list, uext: np.ndarray) -> bool:
        """Le residu est-il VRAIMENT nul ?

        Le ``success`` de MINPACK -- de Levenberg-Marquardt surtout -- ne dit
        que "les deux derniers iteres sont proches" : sur un systeme sans
        solution, ``lm`` s'arrete sur place et se declare reussi avec un
        residu de 1. On juge donc sur le residu lui-meme, relativement a la
        taille des inconnues (des courants de 1e8 ne se ferment pas a 1e-12).
        """
        r = self.residu(zu, xs, uext)
        return bool(np.abs(r).max() <= 1e-8 * (1.0 + np.abs(zu).max()))

    def _diagnostic(self, sol, xs: list, uext: np.ndarray, t: float | None) -> str:
        """Le message d'echec de Newton, rendu exploitable.

        Un ``sol.message`` de MINPACK ne dit rien du reseau ; il faut l'y
        raccrocher. Le residu se lit en deux blocs -- les equations ``g`` des
        elements, puis la loi de cablage -- et ``Cablage`` sait de quel
        element vient chaque ligne, par ses offsets ``oz`` et ``ou``. Savoir
        que le pire terme est une equation de flux de machine, ou au
        contraire une tension de port de convertisseur, oriente d'emblee vers
        le bon modele : dans le premier cas l'etat est parti trop loin, dans
        le second c'est la fermeture du reseau qui ne se boucle plus.
        """
        c = self.c
        r = self.residu(sol.x, xs, uext)
        quand = "" if t is None else f" a t = {t:g} s"
        i = int(np.argmax(np.abs(r)))
        # searchsorted a droite, moins un : sur des offsets repetes (un
        # element sans variable algebrique ne fait pas avancer oz), c'est bien
        # le dernier element de meme offset qui possede la ligne.
        if i < c.nz:
            k = int(np.searchsorted(c.oz, i, side="right")) - 1
            nom = c.cms[k].spec.algebs[i - c.oz[k]]
            pire = f"equation algebrique '{nom}' de {c.els[k].name}"
        else:
            j = i - c.nz
            k = int(np.searchsorted(c.ou, j, side="right")) - 1
            nom = c.els[k].input_names()[j - c.ou[k]]
            pire = f"loi de cablage sur l'entree '{nom}' de {c.els[k].name}"
        ng = float(np.abs(r[: c.nz]).max()) if c.nz else 0.0
        nc = float(np.abs(r[c.nz :]).max()) if c.nu else 0.0
        return (
            f"resolution algebrique en echec{quand} : |residu| = "
            f"{float(np.abs(r).max()):.3e} (equations g : {ng:.3e}, "
            f"loi de cablage : {nc:.3e})\n"
            f"  pire terme : {pire}, residu {r[i]:+.3e}\n"
            f"  {sol.message}"
        )

    def amorcer(self, x, uext, z, u, t: float | None = None):
        """L'amorce ``(x, z, u)`` d'un intervalle d'integration.

        D'ordinaire ``x`` inchange et la solution algebrique a ``x`` fige. Un
        modele a noeuds ideaux n'en a pas toujours (indice 2) : on le
        reinitialise (:func:`~g2egrid.noyau.integrateur.reinitialiser`), ce qui
        remet aussi ``x`` sur la loi des noeuds apres une manoeuvre qui a
        retire une branche parcourue par un courant."""
        x = np.asarray(x, dtype=float)
        if self.ideaux:
            try:
                return reinitialiser(self, x, z, u, uext)
            except RuntimeError as err:
                quand = "" if t is None else f" a t = {t:g} s"
                raise RuntimeError(f"{err}{quand}") from None
        return (x.copy(), *self.resoudre_algebrique(x, uext, z, u, t))

    def derivee(
        self,
        x: np.ndarray,
        uext: np.ndarray,
        z_essai: np.ndarray | None = None,
        u_essai: np.ndarray | None = None,
        t: float | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(dx/dt, z, u)`` : resout l'algebrique, puis evalue les ``f``."""
        z, u = self.resoudre_algebrique(x, uext, z_essai, u_essai, t)
        xs, zs, us = self._x(x), self._z(z), self._u(u)
        parts = [
            cm.f(xi, zi, ui, pi)
            for cm, pi, xi, zi, ui in zip(self.c.cms, self.p, xs, zs, us)
            if cm.n_x
        ]
        return (np.concatenate(parts) if parts else np.zeros(0)), z, u

    def jacobien_etats(
        self, x: np.ndarray, uext: np.ndarray, z: np.ndarray, u: np.ndarray
    ) -> np.ndarray:
        """``d(dx/dt)/dx`` a ``(x, z, u)`` deja resolus : la matrice A instantanee.

        Par le theoreme des fonctions implicites, ``z`` et ``u`` sont des
        fonctions de ``x`` le long de la solution du systeme couple ::

            J_zu . [dz/dx ; du/dx] = -[gx ; -Ry.hx]

        et ``d(dx/dt)/dx = fx + fz.dz/dx + fu.du/dx``. C'est exactement
        l'elimination de :meth:`~g2egrid.reseau.network.Network.build_state_space`, evaluee au point
        courant au lieu du point de fonctionnement. Sans elle, ``Radau``
        estimerait ce jacobien par differences finies : ``n_x`` resolutions
        de Newton par jacobien, ce qui rend l'integration cent fois plus lente.
        """
        c = self.c
        xs, zs, us = self._x(x), self._z(z), self._u(u)
        Fx = np.zeros((c.nx, c.nx))
        Fz = np.zeros((c.nx, c.nz))
        Fu = np.zeros((c.nx, c.nu))
        Gx = np.zeros((c.nz, c.nx))
        Hx = np.zeros((c.ny, c.nx))
        for k, (cm, pi, xi, zi, ui) in enumerate(zip(c.cms, self.p, xs, zs, us)):
            sx, sz, su, sy = c.ox[k], c.oz[k], c.ou[k], c.oy[k]
            fx_ = slice(sx, sx + cm.n_x)
            if cm.n_x:
                Fx[fx_, fx_] = cm.fx(xi, zi, ui, pi)
                Fu[fx_, su : su + cm.n_u] = cm.fu(xi, zi, ui, pi)
                if cm.n_z:
                    Fz[fx_, sz : sz + cm.n_z] = cm.fz(xi, zi, ui, pi)
            if cm.n_z:
                Gx[sz : sz + cm.n_z, fx_] = cm.gx(xi, zi, ui, pi)
            Hx[sy : sy + cm.n_y, fx_] = cm.hx(xi, zi, ui, pi)
        J = self.jacobien_creux(np.concatenate([z, u]), xs, uext)
        rhs = -np.vstack([Gx, -(self._Ry @ Hx)])
        S = spla.splu(J).solve(rhs)  # [dz/dx ; du/dx]
        return Fx + Fz @ S[: c.nz] + Fu @ S[c.nz :]

    # ------------------------------------------------------------ manoeuvres
    def _angle_repere(self, u: np.ndarray) -> float:
        """La valeur courante de ``theta_g``, lue sur l'entree d'un element."""
        c = self.c
        for k, el in enumerate(c.els):
            return float(u[c.ou[k] + el.thetaref_idx()])
        return 0.0

    def appliquer(self, man: Manoeuvre, x, z, u):
        """Le modele APRES la manoeuvre, et l'etat repris : ``(m, x, z, u)``.

        Le modele d'avant n'est pas touche : les elements modifies sont des
        copies. Les etats sont repris par leur nom ; ceux d'un element retire
        disparaissent, celui d'un repere cree prend la valeur de ``theta_g``.
        ``z`` et ``u`` ne sont que des amorces, repris par nom eux aussi.
        """
        net = self.net
        x = np.asarray(x, dtype=float).copy()
        if man.genre == "saut_de_phase":
            return self, self._saut_de_phase(man, x), z, u

        net2 = copy.copy(net)
        net2.elements = [_copie_element(e) for e in net.elements]
        net2.buses = {b: _copie_element(e) for b, e in net.buses.items()}
        net2.condensateurs = dict(getattr(net, "condensateurs", {}))
        par_nom = {e.name: e for e in net2.elements}

        if man.genre == "charge":
            self._varier_charge(par_nom, man)
        elif man.genre == "defaut":
            b = self._bus_de(man.cible)
            if b not in net2.buses:
                raise ValueError(f"defaut : bus {b} inconnu")
            nom = f"Defaut_bus{b}"
            if nom in par_nom:
                raise ValueError(f"defaut : il y a deja un defaut au bus {b}")
            d = Defaut(nom, b, R=man.R, X=man.X, wb=net2.wb)
            d.x0, d.z0 = np.zeros(2), np.zeros(0)
            d.u0 = np.zeros(len(d.input_names()))
            net2.elements.append(d)
        elif man.genre == "elimination":
            nom = f"Defaut_bus{self._bus_de(man.cible)}"
            if nom not in par_nom:
                raise ValueError(f"elimination : pas de defaut en cours ({nom})")
            net2.elements = [e for e in net2.elements if e.name != nom]
        elif man.genre == "defaut_dc":
            b = self._bus_de(man.cible)
            if b not in getattr(net2, "bus_specs_dc", {}):
                raise ValueError(f"defaut_dc : noeud continu {b} inconnu")
            nom = f"DefautDC_bus{b}"
            if nom in par_nom:
                raise ValueError(f"defaut_dc : il y a deja un defaut au noeud {b}")
            d = DefautDC(nom, b, R=man.R, X=man.X, wb=net2.wb)
            d.x0, d.z0 = np.zeros(1), np.zeros(0)
            d.u0 = np.zeros(len(d.input_names()))
            net2.elements.append(d)
        elif man.genre == "elimination_dc":
            nom = f"DefautDC_bus{self._bus_de(man.cible)}"
            if nom not in par_nom:
                raise ValueError(f"elimination_dc : pas de defaut en cours ({nom})")
            net2.elements = [e for e in net2.elements if e.name != nom]
        elif man.genre == "blocage":
            el = par_nom.get(man.cible)
            if el is None or "mode_lim" not in el.params:
                convs = sorted(n for n, e in par_nom.items() if "mode_lim" in e.params)
                raise ValueError(
                    f"blocage : {man.cible!r} n'est pas un convertisseur ; "
                    f"convertisseurs : {convs}"
                )
            if el.bloque:
                return self, x, z, u  # deja bloque
            el.set_param("Imax", IMAX_BLOQUE)
            if "bloque" in el.params:
                el.set_param("bloque", 1.0)  # decouple du noeud continu
            el.bloque = True
        else:  # ouverture
            self._ouvrir(net2, par_nom, man, self._angle_repere(u))

        m2 = self._reporter_modes(ModeleNonLineaire(net2))
        if man.genre == "blocage":
            # le convertisseur bloque passe en mode limite, et y reste
            kb = next(i for i, e in enumerate(m2.c.els) if e.name == man.cible)
            m2 = m2.avec_mode(kb, "mode_lim", 1.0)
        return m2, *self._reprendre(m2, x, z, u)

    @staticmethod
    def _bus_de(cible) -> int:
        c = str(cible)
        return int(c[3:] if c.startswith("bus") else c)

    # ------------------------------------------------------------ modes
    @property
    def _k_modes(self) -> list[int]:
        """Les elements qui ont des modes discrets (limiteurs, butees)."""
        return [k for k, el in enumerate(self.c.els) if el.modes]

    def modes(self) -> dict[str, dict[str, float]]:
        """L'etat courant des drapeaux, par element."""
        return {
            self.c.els[k].name: {
                d: float(self.p[k][self.c.els[k].params.index(d)])
                for d in self.c.els[k].modes
            }
            for k in self._k_modes
        }

    def avec_mode(self, k: int, drapeau: str, valeur: float) -> "ModeleNonLineaire":
        """Le meme modele, un drapeau change : memes equations, meme cablage,
        un seul parametre different. Rien a reconstruire ; seule la
        factorisation gardee par le Newton modifie est a refaire."""
        m2 = copy.copy(self)
        m2.p = list(self.p)
        pk = np.array(self.p[k], dtype=float)
        pk[self.c.els[k].params.index(drapeau)] = valeur
        m2.p[k] = pk
        m2._lu = None
        m2._n_solves = 0
        return m2

    def _reporter_modes(self, m2: "ModeleNonLineaire") -> "ModeleNonLineaire":
        """Les drapeaux de ``self`` reportes sur ``m2``, element par element
        (une manoeuvre ne remet pas un limiteur a zero)."""
        courants = self.modes()
        if not courants:
            return m2
        m2.p = list(m2.p)
        for k in m2._k_modes:
            el = m2.c.els[k]
            if el.name not in courants:
                continue
            pk = np.array(m2.p[k], dtype=float)
            for d, v in courants[el.name].items():
                pk[el.params.index(d)] = v
            m2.p[k] = pk
        return m2

    def _gardes(self, x, z, u) -> list[tuple[int, str, float, float, int]]:
        """``(k, drapeau, g, valeur_apres, sens)`` de toutes les gardes."""
        out = []
        ks = self._k_modes
        if not ks:
            return out
        xs, zs, us = self._x(x), self._z(z), self._u(u)
        for k in ks:
            el, cm = self.c.els[k], self.c.cms[k]
            for d, g, v, sens in el.gardes(cm, xs[k], zs[k], us[k], self.p[k]):
                out.append((k, d, float(g), float(v), int(sens)))
        return out

    def _saut_de_phase(self, man: Manoeuvre, x: np.ndarray) -> np.ndarray:
        nom = man.cible if man.cible.startswith("bus") else f"bus{int(man.cible)}"
        noms = self.c.state_names
        if f"{nom}.vd" not in noms:
            raise ValueError(
                f"saut de phase : {nom} n'a pas de tension d'etat (noeud sans "
                "capacite, ou bus inconnu)"
            )
        i, j = noms.index(f"{nom}.vd"), noms.index(f"{nom}.vq")
        a = math.radians(man.angle)
        vd, vq = x[i], x[j]
        x[i] = vd * math.cos(a) - vq * math.sin(a)
        x[j] = vd * math.sin(a) + vq * math.cos(a)
        return x

    @staticmethod
    def _varier_charge(par_nom: dict, man: Manoeuvre) -> None:
        el = par_nom.get(man.cible)
        if el is None:
            raise ValueError(f"charge {man.cible!r} inconnue")
        if "P0" in el.pdict:  # courant ou puissance constante
            el.set_param("P0", el.pdict["P0"] * (1.0 + man.dp))
            el.set_param("Q0", el.pdict["Q0"] * (1.0 + man.dq))
            return
        if not ({"R", "X"} <= set(el.pdict) and hasattr(el, "P")):
            raise ValueError(f"{man.cible} n'est pas une charge")
        S0 = complex(el.P, el.Q)
        S1 = complex(el.P * (1.0 + man.dp), el.Q * (1.0 + man.dq))
        if abs(S1) < 1e-9 * max(abs(S0), 1e-12):
            raise ValueError(
                f"{man.cible} : charge ramenee a zero ; c'est une ouverture"
            )
        # Z = |V0|^2 / conj(S) a la tension du point de fonctionnement
        Z1 = complex(el.pdict["R"], el.pdict["X"]) * np.conj(S0) / np.conj(S1)
        el.set_param("R", Z1.real)
        el.set_param("X", Z1.imag)
        el.P, el.Q = S1.real, S1.imag

    def _ouvrir(self, net2, par_nom: dict, man: Manoeuvre, theta_g: float) -> None:
        cible = man.cible
        if cible in net2.condensateurs:
            b, B = net2.condensateurs.pop(cible)
            self._reprendre_capacite(net2, b, B)
            return
        el = par_nom.get(cible)
        if el is None:
            connus = sorted([*par_nom, *net2.condensateurs])
            raise ValueError(f"ouverture : {cible!r} inconnu ; elements : {connus}")
        net2.elements = [e for e in net2.elements if e.name != cible]
        for b in el.buses if getattr(el, "b_demi", 0.0) else ():
            self._reprendre_capacite(net2, b, el.b_demi)

        k_ref, _ty, _tw = self.net._reference(self.c.els)
        ref = self.c.els[k_ref].name if k_ref is not None else None
        if cible == ref:
            rep_ = Repere("Repere", wb=net2.wb)
            rep_.theta0 = theta_g
            rep_.x0 = np.array([theta_g])
            rep_.z0 = np.zeros(0)
            rep_.u0 = np.array([1.0, 1.0, theta_g])  # w_in, w_ref, theta_g
            net2.elements.append(rep_)
            net2.reference = rep_.name
            net2.repere_suit = self._plus_grosse_machine(net2)
        elif isinstance(net2.repere_suit, dict) and cible in net2.repere_suit:
            # centre d'inertie : la machine perdue en sort, les autres restent
            reste = {n: w for n, w in net2.repere_suit.items() if n != cible}
            net2.repere_suit = reste or None
        elif cible == net2.repere_suit:
            net2.repere_suit = self._plus_grosse_machine(net2)

    @staticmethod
    def _plus_grosse_machine(net2) -> str | None:
        """La machine que le repere suivra : le plus gros formeur restant."""
        candidats = [
            (e.pdict.get("kS", 0.0), e.name)
            for e in net2.elements
            if "theta_out" in e.svc_outputs
            and Network.sortie_vitesse(e) is not None
            and not isinstance(e, Repere)
        ]
        return max(candidats)[1] if candidats else None

    @staticmethod
    def _reprendre_capacite(net2, b: int, dB: float) -> None:
        """Le noeud ``b`` perd la capacite ``dB`` (pu) ; sans capacite, il
        devient un noeud resistif, comme a l'assemblage."""
        bus = net2.buses[b]
        C = bus.pdict.get("C", 0.0) - dB
        if C > 1e-9:
            bus.set_param("C", C)
            return
        if getattr(net2, "noeuds_sans_capacite", "resistance") == "ideal":
            from g2egrid.elements.bus import IdealBus

            rb = IdealBus(bus.bus, bus.Un, name=bus.name)
            rb.spec_si = dict(bus.spec_si)
            rb.finaliser(net2.wb)
            # la tension : reprise de l'etat du noeud capacitif (_reprendre)
            rb.x0, rb.z0 = np.zeros(0), np.zeros(2)
        else:
            rb = ResistiveBus(bus.bus, bus.Un, name=bus.name)
            rb.spec_si = dict(bus.spec_si)
            if rb.spec_si.get("R") is None:
                rb.spec_si["R"] = net2.r_causalite
            rb.finaliser(net2.wb)
            rb.x0, rb.z0 = np.zeros(0), np.zeros(0)
        rb.u0 = np.zeros(len(rb.input_names()))
        net2.buses[b] = rb

    def _reprendre(self, m2: "ModeleNonLineaire", x, z, u):
        """L'etat du modele ``m2`` repris par nom sur celui de ``self``."""

        def par_nom(noms_src, v_src, noms_dst, defaut):
            src = dict(zip(noms_src, v_src))
            return np.array([src.get(n, d) for n, d in zip(noms_dst, defaut)])

        x2 = par_nom(self.c.state_names, x, m2.c.state_names, m2.c.x0())
        # une variable algebrique nouvelle reprend l'etat de meme nom : la
        # tension d'un noeud capacitif qui perd sa capacite et devient ideal
        z2 = par_nom(
            [*self.c.state_names, *self.noms_z],
            np.concatenate([np.asarray(x, dtype=float), np.asarray(z, dtype=float)]),
            m2.noms_z,
            m2.c.z0(),
        )
        u2 = par_nom(self.noms_u, u, m2.noms_u, m2.c.u0())
        return x2, z2, u2

    def traducteur_uext(self, m2: "ModeleNonLineaire") -> Callable:
        """``uext`` de ``self`` -> ``uext`` de ``m2``, par nom. Une entree d'un
        element retire est oubliee ; une entree nouvelle prend sa valeur
        d'equilibre."""
        src = {n: i for i, n in enumerate(self.c.uext_names)}
        base = m2.uext0()
        paires = [(j, src[n]) for j, n in enumerate(m2.c.uext_names) if n in src]
        dst = np.array([j for j, _ in paires], dtype=int)
        ori = np.array([i for _, i in paires], dtype=int)

        def traduire(v: np.ndarray) -> np.ndarray:
            out = base.copy()
            out[dst] = np.asarray(v)[ori]
            return out

        return traduire

    # ------------------------------------------------------------- entrees
    def uext0(self) -> np.ndarray:
        return self.c.uext0()

    def scenario(self, *evenements) -> Callable[[float], np.ndarray]:
        """Une fonction ``uext(t)`` composant plusieurs :class:`Evenement`.

        Les evenements s'AJOUTENT : plusieurs sur la meme entree se cumulent
        -- un creux de consigne se decrit par deux echelons de signes opposes
        -- et rien n'oblige a ce qu'ils portent sur la meme entree.

        La fonction rendue porte un attribut ``instants``, la liste triee de
        ses discontinuites. :meth:`simuler` la lit toute seule : l'appelant
        n'a jamais a redonner des dates qu'il a deja ecrites dans ses
        evenements, ni a se souvenir que l'oublier coute du temps de calcul.
        """
        u0 = self.uext0()
        manoeuvres = tuple(
            sorted(
                (ev for ev in evenements if isinstance(ev, Manoeuvre)),
                key=lambda m: m.instant,
            )
        )
        evenements = [ev for ev in evenements if not isinstance(ev, Manoeuvre)]
        plan = [(self.c.iu(ev.entree), ev) for ev in evenements]
        # une courbe ABSOLUE remplace la consigne d'equilibre ; les ecarts
        # (evenements, courbes relatives) s'y ajoutent ensuite
        absolues = [(i, ev) for i, ev in plan if getattr(ev, "absolue", False)]
        ecarts = [(i, ev) for i, ev in plan if not getattr(ev, "absolue", False)]
        if len({i for i, _ in absolues}) < len(absolues):
            raise ValueError("deux courbes absolues sur la meme entree")
        for i, ev in absolues:
            depart = ev.valeur(float(ev.t[0]))
            if abs(depart - u0[i]) > 1e-6 * max(1.0, abs(u0[i])):
                warnings.warn(
                    f"{ev.entree}: la courbe part de {depart:.6g} alors que "
                    f"l'equilibre vaut {u0[i]:.6g} : la simulation demarre par un "
                    "saut de consigne",
                    stacklevel=2,
                )

        def uext(t: float) -> np.ndarray:
            u = u0.copy()
            for i, ev in absolues:
                u[i] = ev.valeur(t)
            for i, ev in ecarts:
                u[i] += ev.valeur(t)
            return u

        uext.instants = tuple(
            sorted(
                {s for ev in (*evenements, *manoeuvres) for s in ev.discontinuites()}
            )
        )
        uext.manoeuvres = manoeuvres
        return uext

    def echelon(self, nom: str, amplitude: float, t0: float = 0.0):
        """Une fonction ``uext(t)`` : l'entree ``nom`` saute de ``amplitude``
        a l'instant ``t0``, les autres restent a leur valeur d'equilibre.

        Cas particulier de :meth:`scenario` a un seul evenement ; la fonction
        rendue porte donc, elle aussi, son ``instants``.
        """
        return self.scenario(Evenement(t0, nom, amplitude))

    # ------------------------------------------------------------ simulation
    def simuler(
        self,
        t_span: tuple[float, float],
        x0: np.ndarray | None = None,
        uext_fn=None,
        t_eval: np.ndarray | None = None,
        methode: str = "auto",
        rtol: float = 1e-6,
        atol: float = 1e-8,
        premier_pas: float | None = None,
        sorties: list[str] | None = None,
        instants: Sequence[float] | None = None,
        manoeuvres: Sequence[Manoeuvre] | None = None,
        limiteurs: bool = True,
        pas: float | None = None,
        protections: Sequence[Protection] | None = None,
    ) -> Trajectoire:
        """Integre le modele non lineaire.

        ``x0`` vaut le point de fonctionnement par defaut ; ``uext_fn(t)``
        rend le vecteur d'entrees externes a l'instant ``t`` (cf.
        :meth:`echelon` et :meth:`scenario`), et vaut la consigne d'equilibre
        par defaut.

        Les modes les plus raides du modele montent a 10^6-10^7 rad/s
        (resistances de causalite, troncons de cable courts) : l'integrateur
        par defaut est donc ``Radau``, implicite et L-stable. Chaque Newton
        repart de la solution du precedent.

        DECOUPAGE AUX DISCONTINUITES -- ``instants`` donne les dates ou la
        consigne saute ou casse de pente ; a defaut, elles sont lues sur
        l'attribut ``instants`` du scenario. L'integration se fait alors par
        intervalles CONTINUS enchaines, chacun repartant de l'etat final du
        precedent.

        Sans ce decoupage, un integrateur a pas adaptatif ne decouvre
        l'echelon qu'APRES l'avoir franchi, par un rejet d'erreur, et ramene
        son pas par bissections successives jusqu'a le border a la tolerance
        pres. Chaque bissection coute ici une resolution de Newton couple et
        une factorisation du jacobien, et rien de ce travail ne sert a la
        trajectoire. Lui donner la date le supprime entierement.
        ``instants=[]`` retablit l'integration d'un bloc, ce qui n'a d'autre
        interet que de mesurer l'ecart.

        MANOEUVRES -- une :class:`Manoeuvre` (ouverture, echelon de charge,
        saut de phase), lue sur l'attribut ``manoeuvres`` du scenario ou passee
        ici, change le MODELE a son instant : l'integration continue avec le
        modele d'apres (:meth:`appliquer`), et la trajectoire rendue garde ses
        segments (:attr:`Trajectoire.segments`) pour la relecture. Ses colonnes
        sont celles du modele de depart, suivies des etats apparus ; un etat
        d'element retire vaut ``NaN`` apres son retrait.

        INTEGRATEUR -- ``methode="auto"`` (le defaut) prend BDF2 pour un
        reseau a noeuds ideaux (le defaut des reseaux, FR-005), que seul BDF2
        integre, et Radau sinon. ``methode="BDF2"`` remplace ``solve_ivp`` par
        l'integrateur DAE monolithique de :mod:`g2egrid.noyau.integrateur` : un pas
        FIXE ``pas`` (``1/(400.f0)`` par defaut), etats et contraintes resolus
        ensemble, les resonances plus rapides que le pas AMORTIES au lieu
        d'etre suivies. C'est le choix pour un defaut ou une manoeuvre ; Radau,
        a pas adaptatif et d'ordre 5, reste celui de la validation fine du
        petit signal (ecart au modele lineaire ∝ a^2).

        LIMITEURS ET BUTEES -- avec ``limiteurs=True`` (le defaut), les gardes
        des modes discrets (:mod:`g2egrid.elements.limiteurs`) sont surveillees
        comme des evenements de ``solve_ivp`` : a l'instant precis ou l'une
        franchit zero, l'integration s'arrete, le drapeau bascule
        (:meth:`avec_mode`) et elle reprend. Les basculements sont rendus dans
        :attr:`Trajectoire.basculements`. Une garde deja franchie au depart
        (un GFM surcharge a l'equilibre) fait partir en mode limite.

        PROTECTIONS -- chaque :class:`Protection` armee est une garde de plus :
        a l'instant ou sa mesure franchit le seuil, l'integration s'arrete, ses
        manoeuvres sont posees ``retard`` plus tard, et l'integration reprend.
        Les declenchements sont rendus dans :attr:`Trajectoire.declenchements`
        (``(instant, protection, manoeuvres)``).

        AVANCEMENT -- la simulation signale ou elle en est
        (:func:`g2egrid.noyau.progres.signaler`) au debut de chaque intervalle et,
        en BDF2, tous les ``PAS_RAPPEL`` pas : une tache de l'atelier la suit
        et peut l'interrompre.

        Sur l'intervalle qui PRECEDE une discontinuite, la consigne est
        evaluee juste avant la date de celle-ci : sans ce recul d'un ulp, le
        dernier pas de l'intervalle verrait deja la valeur d'apres, et
        l'echelon serait avance d'un pas.
        """
        x0 = self.c.x0() if x0 is None else np.asarray(x0, dtype=float)
        self._n_solves = 0
        if methode == "auto":
            methode = "BDF2" if self.ideaux else "Radau"
        if self.ideaux and methode != "BDF2":
            raise ValueError(
                f"noeud(s) ideal(aux) {', '.join(self.ideaux[:5])}"
                f"{'...' if len(self.ideaux) > 5 else ''} : leur loi des noeuds "
                "peut etre une contrainte d'indice 2, que seul BDF2 integre "
                "(methode='BDF2') ; Radau resout l'algebrique a x fige, ce qui "
                "n'a pas de solution."
            )

        if instants is None:
            instants = getattr(uext_fn, "instants", ())
        manoeuvres = [*getattr(uext_fn, "manoeuvres", ()), *(manoeuvres or ())]
        t_ini, t_fin = float(t_span[0]), float(t_span[1])
        dates = {float(s) for s in instants} | {float(m.instant) for m in manoeuvres}
        coupures = sorted(s for s in dates if t_ini < s < t_fin)
        bornes = [t_ini, *coupures, t_fin]
        if t_eval is not None:
            t_eval = np.asarray(t_eval, dtype=float)

        implicite = methode in ("Radau", "BDF", "LSODA")
        modele, traduire = self, None
        utilises = [self]
        cache = {"z": self.c.z0(), "u": self.c.u0()}
        noms_colonnes = list(self.c.state_names)
        morceaux: list[tuple[np.ndarray, np.ndarray, ModeleNonLineaire]] = []
        basculements: list[tuple[float, str, str, float]] = []
        #: element -> instant avant lequel il ne peut pas ressortir de son mode
        #: (element, drapeau) -> instant ou ce drapeau peut de nouveau changer
        verrou: dict[tuple[int, str], float] = {}
        x_depart, message = x0, ""
        armees = list(protections or [])
        declenchements: list[tuple[float, str, list]] = []
        appliquees: set[int] = set()
        #: relais a retombee dont la temporisation court : (protection, terme)
        en_attente: list[tuple[Protection, float]] = []
        demarrages: list[tuple[float, str]] = []
        retombees: list[tuple[float, str]] = []

        def source(t):
            return self.uext0() if uext_fn is None else uext_fn(t)

        def declencher(p, t, m):
            """``p`` declenche a ``t`` : ses manoeuvres entrent au plan, et les
            bornes de l'integration avec elles. Rend l'indice de ``t``. Un
            relais a retombee, lui, ne fait que DEMARRER sa temporisation : son
            terme devient une borne, ou il declenchera s'il n'est pas retombe."""
            nonlocal bornes
            armees.remove(p)
            if p.retombee is not None and p.retard > 0:
                terme = t + p.retard
                en_attente.append((p, terme))
                demarrages.append((t, p.nom))
                bornes = sorted(
                    set(bornes) | {t} | ({terme} if terme < t_fin else set())
                )
                return bornes.index(t)
            mans = p.manoeuvres(t + p.retard, m)
            manoeuvres.extend(mans)
            declenchements.append((t, p.nom, [(mm.genre, mm.cible) for mm in mans]))
            nouvelles = {t} | {mm.instant for mm in mans if mm.instant < t_fin}
            bornes = sorted(set(bornes) | nouvelles)
            return bornes.index(t)

        def avancement(tt):
            if t_fin > t_ini:
                signaler((tt - t_ini) / (t_fin - t_ini), f"t = {tt:.2f} s")

        k = -1
        while k < len(bornes) - 2:
            k += 1
            a, b = bornes[k], bornes[k + 1]
            avancement(a)

            # --- les relais a retombee dont la temporisation s'acheve ici :
            #     le defaut a tenu, ils declenchent ---
            for pp, terme in [e for e in en_attente if e[1] == a]:
                en_attente.remove((pp, terme))
                if pp.indices(modele) is None:
                    continue  # l'equipement est deja parti
                mans = pp.manoeuvres(a, modele)
                manoeuvres.extend(mans)
                declenchements.append(
                    (a - pp.retard, pp.nom, [(mm.genre, mm.cible) for mm in mans])
                )

            # --- les manoeuvres datees du debut de cet intervalle ---
            for man in [
                m_
                for m_ in manoeuvres
                if float(m_.instant) == a and id(m_) not in appliquees
            ]:
                appliquees.add(id(man))
                uext_a = source(a) if traduire is None else traduire(source(a))
                x_depart, z_a, u_a = modele.amorcer(
                    x_depart, uext_a, cache["z"], cache["u"], a
                )
                apres, x_depart, z2, u2 = modele.appliquer(man, x_depart, z_a, u_a)
                cache = {"z": z2, "u": u2}
                if apres is not modele:
                    modele = apres
                    utilises.append(modele)
                    traduire = self.traducteur_uext(modele)
                    noms_colonnes += [
                        n for n in modele.c.state_names if n not in noms_colonnes
                    ]

            avant = np.nextafter(b, a)  # dernier flottant avant la coupure

            if t_eval is None:
                demandes = None
            else:
                dernier = k == len(bornes) - 2
                dedans = t_eval >= a
                dedans &= (t_eval <= b) if dernier else (t_eval < b)
                demandes = t_eval[dedans]

            # --- l'intervalle, coupe aux basculements de mode ---
            t_a = a
            while True:
                u_defaut = modele.uext0()

                def consigne(t, _b=b, _avant=avant, _tr=traduire, _ud=u_defaut):
                    if uext_fn is None:
                        return _ud
                    v = uext_fn(_avant if t >= _b else t)
                    return v if _tr is None else _tr(v)

                def verrouille(kk, d, t, _m):
                    # un drapeau ACTIF depuis moins de DUREE_MIN_MODE ne peut
                    # pas retomber ; les autres drapeaux du meme element (la
                    # limitation de courant quand la tension sature) restent
                    # libres de basculer
                    actif = _m.p[kk][_m.c.els[kk].params.index(d)] != 0.0
                    return actif and verrou.get((kk, d), -np.inf) > t

                # une garde deja franchie au depart : on bascule avant d'integrer
                gardes = []
                if limiteurs and modele._k_modes:
                    for _essai in range(20):
                        x_depart, z_a, u_a = modele.amorcer(
                            x_depart, consigne(t_a), cache["z"], cache["u"], t_a
                        )
                        cache = {"z": z_a, "u": u_a}
                        gardes = [
                            g
                            for g in modele._gardes(x_depart, z_a, u_a)
                            if not verrouille(g[0], g[1], t_a, modele)
                        ]
                        franchie = next(
                            (g for g in gardes if g[4] * g[2] > TOLERANCE_GARDE), None
                        )
                        if franchie is None:
                            break
                        kk, d, _g, v, _s = franchie
                        modele = modele.avec_mode(kk, d, v)
                        utilises.append(modele)
                        basculements.append((t_a, modele.c.els[kk].name, d, v))
                        verrou[(kk, d)] = t_a + DUREE_MIN_MODE
                    else:
                        raise RuntimeError(
                            f"modes discrets instables a t = {t_a:g} s : les "
                            "limiteurs basculent sans fin"
                        )

                # la fin de ce sous-intervalle : b, ou la levee d'un verrou,
                # pour reprendre avec la garde de sortie reamorcee
                leves = [
                    verrou[(kk, d)]
                    for kk in modele._k_modes
                    for d in modele.c.els[kk].modes
                    if verrouille(kk, d, t_a, modele) and verrou[(kk, d)] < b
                ]
                fin = min([b, *leves])
                # les gardes surveillees, reperees par leur rang dans la liste
                # COMPLETE que rend _gardes (celle de ``evaluer``)
                rangs = []
                if gardes:
                    complete = modele._gardes(x_depart, cache["z"], cache["u"])
                    rangs = [
                        i
                        for i, g in enumerate(complete)
                        if not verrouille(g[0], g[1], t_a, modele)
                    ]
                    gardes = [complete[i] for i in rangs]

                # --- les relais dont la temporisation court : retombes des le
                #     depart, ou surveilles par une garde de retombee ---
                att = []
                for pp, terme in list(en_attente):
                    ix = pp.indices(modele)
                    if ix is None:  # l'equipement est parti : plus rien a temporiser
                        en_attente.remove((pp, terme))
                    elif pp.garde_retombee(x_depart, ix) < 0:
                        en_attente.remove((pp, terme))
                        armees.append(pp)
                        retombees.append((t_a, pp.nom))
                    else:
                        att.append((pp, ix))

                # --- les protections armees dont l'equipement est la ---
                prot = [
                    (pp, ix) for pp in armees if (ix := pp.indices(modele)) is not None
                ]
                deja = next((pp for pp, ix in prot if pp.garde(x_depart, ix) > 0), None)
                if deja is not None:
                    # franchie des le depart : elle declenche maintenant
                    k = declencher(deja, t_a, modele) - 1
                    break

                evenements_modes = []
                if gardes:
                    memo: dict = {}

                    def evaluer(t, x, _m=modele, _u=consigne, _memo=memo):
                        cle = (t, x.tobytes())
                        if _memo.get("cle") != cle:
                            z_e, u_e = _m.resoudre_algebrique(
                                x, _u(t), cache["z"], cache["u"], t
                            )
                            _memo["cle"], _memo["g"] = cle, _m._gardes(x, z_e, u_e)
                        return _memo["g"]

                    for i_g, (_k, _d, _g, _v, sens) in zip(rangs, gardes):

                        def ev(t, x, _i=i_g, _ev=evaluer):
                            return _ev(t, x)[_i][2]

                        ev.terminal = True
                        ev.direction = sens
                        evenements_modes.append(ev)

                def second_membre(t, x, _u=consigne, _m=modele):
                    xdot, z, u = _m.derivee(x, _u(t), cache["z"], cache["u"], t)
                    cache["z"], cache["u"] = z, u
                    return xdot

                def jacobien(t, x, _u=consigne, _m=modele):
                    uext = _u(t)
                    z, u = _m.resoudre_algebrique(x, uext, cache["z"], cache["u"], t)
                    cache["z"], cache["u"] = z, u
                    return _m.jacobien_etats(x, uext, z, u)

                if demandes is None:
                    rendus, ajout = None, False
                else:
                    reste = demandes[demandes >= t_a] if t_a > a else demandes
                    if fin < b:
                        reste = reste[reste < fin]
                    # Avec ``t_eval``, ``sol.y`` ne contient QUE ces instants :
                    # l'etat final, dont repart la suite, n'y serait pas. On
                    # ajoute donc ``fin``, quitte a le retirer de ce qu'on garde.
                    ajout = not len(reste) or reste[-1] < fin
                    rendus = np.append(reste, fin) if ajout else reste

                evenements_protections = []
                for pp, ix in prot:

                    def ev_p(t, x, _p=pp, _ix=ix):
                        return _p.garde(x, _ix)

                    ev_p.terminal = True
                    ev_p.direction = 1
                    evenements_protections.append(ev_p)
                for pp, ix in att:

                    def ev_r(t, x, _p=pp, _ix=ix):
                        return _p.garde_retombee(x, _ix)

                    ev_r.terminal = True
                    ev_r.direction = -1
                    evenements_protections.append(ev_r)

                if methode == "BDF2":
                    x_depart, z_a, u_a = modele.amorcer(
                        x_depart, consigne(t_a), cache["z"], cache["u"], t_a
                    )

                    def valeurs(x, z, u, _m=modele, _r=rangs, _p=prot, _a=att):
                        tout = _m._gardes(x, z, u) if _r else []
                        return (
                            [tout[i][2] for i in _r]
                            + [pp.garde(x, ix) for pp, ix in _p]
                            + [pp.garde_retombee(x, ix) for pp, ix in _a]
                        )

                    sens_tous = (
                        [g[4] for g in gardes] + [1] * len(prot) + [-1] * len(att)
                    )
                    sol = integrer_bdf2(
                        modele,
                        (t_a, fin),
                        x_depart,
                        z_a,
                        u_a,
                        consigne,
                        pas or 1.0 / (400.0 * self.net.f0),
                        t_rendus=rendus,
                        gardes=(valeurs, sens_tous) if sens_tous else None,
                        rappel=avancement,
                        # un modele d'indice 2 ne se relit pas en re-resolvant
                        # l'algebrique : on garde celui de l'integration
                        garder_algebriques=bool(modele.ideaux),
                    )
                    cache = {"z": sol.z, "u": sol.u}
                else:
                    sol = solve_ivp(
                        second_membre,
                        (t_a, fin),
                        x_depart,
                        method=methode,
                        t_eval=rendus,
                        rtol=rtol,
                        atol=atol,
                        first_step=premier_pas,
                        jac=jacobien if implicite else None,
                        events=(evenements_modes + evenements_protections) or None,
                    )
                if not sol.success:
                    raise RuntimeError(
                        f"integration en echec sur [{t_a:g}, {fin:g}] s : {sol.message}"
                    )
                modele._verifier_fini(sol)
                message = sol.message
                # une garde franchie avant le premier instant demande : aucun
                # point rendu, et ``sol.y`` est alors une liste vide
                tk = np.asarray(sol.t, dtype=float)
                xk = np.asarray(sol.y, dtype=float).reshape(len(x_depart), -1).T
                zuk = getattr(sol, "zu", None)
                arrete = sol.status == 1  # une garde a franchi zero
                if ajout and not arrete:
                    tk, xk = tk[:-1], xk[:-1]
                    zuk = None if zuk is None else zuk[:-1]
                if t_eval is None and morceaux:
                    tk, xk = tk[1:], xk[1:]  # la jonction appartient au precedent
                    zuk = None if zuk is None else zuk[1:]
                morceaux.append((tk, xk, modele, zuk))

                if not arrete:
                    x_depart = sol.y[:, -1]
                    if fin >= b:
                        break
                    t_a = fin  # un verrou se leve : on reamorce sa garde
                    continue
                i_ev = min(
                    (i for i, te in enumerate(sol.t_events) if len(te)),
                    key=lambda i: sol.t_events[i][0],
                )
                t_a = float(sol.t_events[i_ev][0])
                x_depart = sol.y_events[i_ev][0]
                if i_ev >= len(gardes) + len(prot):
                    # un relais retombe avant le terme : il se reamorce, et
                    # l'integration repart de cet instant
                    pp = att[i_ev - len(gardes) - len(prot)][0]
                    en_attente[:] = [e for e in en_attente if e[0] is not pp]
                    armees.append(pp)
                    retombees.append((t_a, pp.nom))
                    bornes = sorted(set(bornes) | {t_a})
                    k = bornes.index(t_a) - 1
                    break
                if i_ev >= len(gardes):
                    # une protection : ses manoeuvres entrent au plan, et
                    # l'integration repart de son instant
                    k = declencher(prot[i_ev - len(gardes)][0], t_a, modele) - 1
                    break
                kk, d, _g, v, _s = gardes[i_ev]
                modele = modele.avec_mode(kk, d, v)
                utilises.append(modele)
                basculements.append((t_a, modele.c.els[kk].name, d, v))
                verrou[(kk, d)] = t_a + DUREE_MIN_MODE
                # sans t_eval, ``sol`` s'arrete SUR l'instant du basculement :
                # il appartient au morceau d'avant, et la suite ne le redonne pas

        # --- la trajectoire : colonnes du depart, puis celles apparues ---
        col = {n: i for i, n in enumerate(noms_colonnes)}
        blocs_t, blocs_x, segments, i0 = [], [], [], 0
        algebriques: list = []
        if any(zk is not None for *_r, zk in morceaux):
            for tk, _xk, mm, zk in morceaux:
                if zk is None:
                    algebriques += [None] * len(tk)
                else:
                    nz = mm.c.nz
                    algebriques += [(r[:nz].copy(), r[nz:].copy()) for r in zk]
        for tk, xk, mm, _zk in morceaux:
            cols = np.array([col[n] for n in mm.c.state_names], dtype=int)
            if len(cols) == len(noms_colonnes) and (cols == np.arange(len(cols))).all():
                X = xk
            else:
                X = np.full((len(tk), len(noms_colonnes)), np.nan)
                X[:, cols] = xk
            blocs_t.append(tk)
            blocs_x.append(X)
            if segments and segments[-1][2] is mm:
                segments[-1] = (segments[-1][0], i0 + len(tk), mm, cols)
            else:
                segments.append((i0, i0 + len(tk), mm, cols))
            i0 += len(tk)

        tr = Trajectoire(
            t=np.concatenate(blocs_t),
            x=np.vstack(blocs_x),
            etats=noms_colonnes,
            n_solves=sum(m_._n_solves for m_ in utilises),
            message=message,
            segments=segments,
            basculements=basculements,
            declenchements=declenchements,
            demarrages=demarrages,
            retombees=retombees,
            algebriques=algebriques,
        )
        if sorties:
            tr.sorties = self.sorties(tr, uext_fn, sorties)
        return tr

    def _verifier_fini(self, sol) -> None:
        """Refuse une trajectoire qui contient un ``NaN`` ou un infini.

        Un reseau instable finit par faire deborder ses grandeurs, et
        ``solve_ivp`` peut rendre ``success = True`` avec des ``NaN`` dedans :
        Newton, lui, converge parfaitement sur des ``NaN``, puisque aucune
        comparaison a une tolerance n'echoue. Mieux vaut lever ici, en nommant
        l'etat et l'instant, que laisser la pollution atteindre les mesures.
        """
        fini = np.isfinite(sol.y)
        if fini.all():
            return
        i = int(np.flatnonzero(~fini.all(axis=1))[0])
        j = int(np.flatnonzero(~fini.all(axis=0))[0])
        raise RuntimeError(
            f"integration divergente : {self.c.state_names[i]} n'est plus fini "
            f"a t = {sol.t[j]:g} s"
        )

    # --------------------------------------------- relecture d'une trajectoire
    def _segments(self, tr: Trajectoire) -> list:
        """Les segments d'une trajectoire, ou un seul -- ce modele -- pour une
        trajectoire construite a la main."""
        return tr.segments or [(0, len(tr.t), self, np.arange(len(tr.etats)))]

    def _le_long(self, tr: Trajectoire, uext_fn=None) -> Iterator[tuple]:
        """Iterateur ``(j, t, x, z, u, y, modele)`` sur les instants rendus.

        ``solve_ivp`` ne rend que les etats : les ``z`` et ``u`` qu'il a
        resolus vivent a SES points de pas, pas aux instants rendus. On
        re-resout donc le systeme couple une fois par instant rendu, ce qui
        coute a peu pres autant que l'integration elle-meme -- d'ou le fait
        que ce soit optionnel, et que :meth:`sorties` et :func:`mesures`
        partagent cette boucle-ci plutot que d'en avoir chacune une.

        Chaque instant est relu avec le MODELE de son segment (cf.
        :attr:`Trajectoire.segments`) : ``x``, ``z``, ``u`` et ``y`` sont les
        siens, et ``uext_fn`` -- ecrite pour le modele de depart -- lui est
        traduite par nom. ``y`` est le vecteur COMPLET des sorties d'element,
        courants de port compris ; ``Cablage.rows_sorties`` dit lesquelles
        sont exposees.
        """
        for i0, i1, mm, cols in self._segments(tr):
            traduire = None if mm is self else self.traducteur_uext(mm)
            u_defaut = mm.uext0()
            z_e, u_e = mm.c.z0(), mm.c.u0()
            c = mm.c
            for j in range(i0, i1):
                t = tr.t[j]
                if uext_fn is None:
                    uext = u_defaut
                else:
                    uext = uext_fn(t) if traduire is None else traduire(uext_fn(t))
                x = tr.x[j, cols]
                gardes = tr.algebriques[j] if tr.algebriques else None
                if gardes is not None:
                    z_e, u_e = gardes
                elif mm.ideaux:
                    raise ValueError(
                        "trajectoire d'un modele a noeuds ideaux sans ses variables "
                        "algebriques : elle ne se relit pas en re-resolvant a x fige "
                        "(indice 2) ; la simuler en BDF2 les garde"
                    )
                else:
                    z_e, u_e = mm.resoudre_algebrique(x, uext, z_e, u_e, t)
                xs, zs, us = mm._x(x), mm._z(z_e), mm._u(u_e)
                y = np.concatenate(
                    [
                        cm.h(xi, zi, ui, pi)
                        for cm, pi, xi, zi, ui in zip(c.cms, mm.p, xs, zs, us)
                    ]
                )
                yield j, t, x, z_e, u_e, y, mm

    def sorties(
        self, tr: Trajectoire, uext_fn=None, noms: list[str] | None = None
    ) -> dict[str, np.ndarray]:
        """Les trajectoires des sorties nommees, le long d'une trajectoire.

        Une sortie d'un element retire par une manoeuvre vaut ``NaN`` apres
        son retrait.
        """
        segs = self._segments(tr)
        connues: list[str] = []
        for _i0, _i1, mm, _cols in segs:
            connues += [n for n in mm.c.output_names if n not in connues]
        noms = noms or connues
        inconnues = [n for n in noms if n not in connues]
        if inconnues:
            raise ValueError(f"sortie(s) inconnue(s) : {inconnues}")
        lignes: dict[int, list[tuple[int, int]]] = {}
        for _i0, _i1, mm, _cols in segs:
            idx = {n: i for i, n in enumerate(mm.c.output_names)}
            lignes[id(mm)] = [
                (j, int(mm.c.rows_sorties[idx[n]]))
                for j, n in enumerate(noms)
                if n in idx
            ]
        out = np.full((len(tr.t), len(noms)), np.nan)
        for j, _t, _x, _z, _u, y, mm in self._le_long(tr, uext_fn):
            for jj, r in lignes[id(mm)]:
                out[j, jj] = y[r]
        return {n: out[:, j] for j, n in enumerate(noms)}

    # ------------------------------------------------------------ controle
    def residu_equilibre(self) -> tuple[float, float]:
        """``(|dx/dt|, |z - z0| + |u - u0|)`` au point de fonctionnement.

        Les deux doivent etre nuls : c'est la verification croisee entre le
        modele non lineaire et le point autour duquel on linearise. Si le
        second ne l'est pas, les points de fonctionnement par element ne
        forment pas un point fixe du reseau entier -- invisible en modal,
        fatal en temporel.

        Avec des noeuds ideaux, l'algebrique a ``x`` fige n'a pas toujours de
        solution unique (indice 2) : on ne le resout pas, on evalue au point
        de fonctionnement les ``f`` et le residu du systeme couple, nuls
        ensemble en un point fixe.
        """
        c = self.c
        if self.ideaux:
            x, z, u = c.x0(), c.z0(), c.u0()
            xs, zs, us = self._x(x), self._z(z), self._u(u)
            parts = [
                cm.f(xi, zi, ui, pi)
                for cm, pi, xi, zi, ui in zip(c.cms, self.p, xs, zs, us)
                if cm.n_x
            ]
            xdot = np.concatenate(parts) if parts else np.zeros(0)
            r = self.residu(np.concatenate([z, u]), xs, self.uext0())
            return float(np.abs(xdot).max()), float(np.abs(r).max())
        xdot, z, u = self.derivee(c.x0(), self.uext0(), c.z0(), c.u0())
        ecart = np.abs(z - c.z0()).max() if c.nz else 0.0
        return float(np.abs(xdot).max()), float(max(ecart, np.abs(u - c.u0()).max()))


# ----------------------------------------------------------------------
# mesures physiques le long d'une trajectoire
# ----------------------------------------------------------------------
def _plan_mesures(m: ModeleNonLineaire) -> dict[str, tuple[str, int, int]]:
    """Le catalogue des mesures du reseau : nom -> ``(genre, element, port)``.

    Il sert deux fois : a valider le ``quoi`` de :func:`mesures`, et a
    n'evaluer dans la boucle par instant que ce qui a ete demande. Sur un
    reseau de quelques centaines d'etats, le catalogue complet compte
    plusieurs centaines de grandeurs ; les calculer toutes quand on ne veut
    qu'une frequence multiplierait le cout de la relecture.

    Un noeud sans capacite (:class:`~g2egrid.elements.bus.ResistiveBus`) n'a
    pas d'etat de tension, donc pas de ``dv/dt`` analytique : il n'offre ni
    frequence brute ni frequence filtree. C'est assume -- sa tension est une
    relation algebrique ``R.isum``, et sa derivee ne s'obtiendrait qu'en
    derivant tout le reste du reseau.

    UN NOEUD CONTINU N'EST PAS UN NOEUD ALTERNATIF, et la distinction ne se
    voit pas a l'heritage : ``BusDC`` derive de ``Bus``, pour que le cablage
    le retrouve comme n'importe quel noeud. Mais il n'a qu'UNE sortie, sa
    tension, la ou un noeud alternatif en a deux (``vd_o``, ``vq_o``). Le
    lire comme un noeud alternatif prendrait la premiere sortie de l'element
    SUIVANT pour son ``vq`` : pas une erreur bruyante, un module de tension
    faux. Un noeud continu offre donc ``Vdc_busN``, et ni angle ni frequence
    -- un continu n'a pas de phase, et sa "frequence" n'existe pas.
    """
    plan: dict[str, tuple[str, int, int]] = {}
    for k, el in enumerate(m.c.els):
        if isinstance(el, Bus):
            n = el.bus
            if getattr(el, "dc", False):
                plan[f"Vdc_bus{n}"] = ("Vdc", k, 0)
                continue
            plan[f"V_bus{n}"] = ("V", k, 0)
            plan[f"angle_bus{n}"] = ("angle", k, 0)
            if m.c.cms[k].n_x:
                plan[f"f_bus{n}"] = ("f", k, 0)
                plan[f"fm_bus{n}"] = ("fm", k, 0)
            continue
        for pk, pt in enumerate(el.ports):
            bout = "" if len(el.ports) == 1 else "_" + _BOUT.get(pt, pt)
            plan[f"P_{el.name}{bout}"] = ("P", k, pk)
            plan[f"Q_{el.name}{bout}"] = ("Q", k, pk)
        # un port continu ne porte qu'une puissance, et pas de reactif
        for pk, pt in enumerate(getattr(el, "ports_dc", ())):
            bout = "_dc" if len(el.ports_dc) == 1 else "_" + _BOUT.get(pt, pt) + "dc"
            plan[f"P_{el.name}{bout}"] = ("Pdc", k, pk)
    return plan


def _premier_ordre(v: np.ndarray, t: np.ndarray, tau: float) -> np.ndarray:
    """Filtre passe-bas du premier ordre, EXACT a pas irregulier.

    Une trajectoire n'est a pas constant que si on lui a impose ``t_eval`` ;
    autrement l'integrateur rend ses propres pas, resserres precisement dans
    les transitoires. Le coefficient usuel ``a = dt/tau`` d'un filtre
    numerique y donnerait une constante de temps differente a chaque pas,
    donc un filtre plus mordant la ou on le regarde. La reponse exacte du
    premier ordre a une entree tenue sur ``dt`` vaut ``a = 1 - exp(-dt/tau)``,
    recalcule pas par pas : le resultat ne depend alors plus du maillage.
    """
    y = np.empty_like(v)
    if not len(v):
        return y
    y[0] = v[0]
    for k in range(1, len(v)):
        a = 1.0 - math.exp(-(t[k] - t[k - 1]) / tau)
        y[k] = y[k - 1] + a * (v[k] - y[k - 1])
    return y


def mesures(
    m: ModeleNonLineaire,
    tr: Trajectoire,
    uext_fn=None,
    quoi: Sequence[str] | None = None,
) -> dict[str, np.ndarray]:
    """Les grandeurs PHYSIQUES le long d'une trajectoire, pas les etats bruts.

    Ce que rend :meth:`ModeleNonLineaire.simuler` est un vecteur d'etats :
    des flux de machine, des courants de branche, des tensions dq de noeud.
    Ce qu'on compare a un enregistrement, c'est un module de tension, une
    frequence, un transit. D'ou cette couche, qui les reconstruit a partir
    des ``(x, z, u, y)`` re-resolus par :meth:`ModeleNonLineaire._le_long`.

    Par bus ``N`` : ``V_busN``, ``angle_busN`` (degres), ``f_busN`` et
    ``fm_busN`` (Hz) -- ou, pour un bus CONTINU, le seul ``Vdc_busN`` (pu),
    un continu n'ayant ni phase ni frequence. Par port continu,
    ``P_<nom>_dc``, sans reactif. Par element a un port : ``P_<nom>``,
    ``Q_<nom>``, en
    convention INJECTEE au bus, celle de tout le paquet -- un generateur
    injecte du positif, une charge du negatif. Par branche, les quatre
    ``P_<nom>_from``, ``Q_<nom>_from``, ``P_<nom>_to``, ``Q_<nom>_to``, avec
    la meme convention : la puissance injectee AU bus par la branche, donc
    negative du cote d'ou part le transit.

    RESERVE SUR LES TRANSITS DE LIGNE -- un load flow compte dans le transit
    d'une ligne la moitie de sa charge capacitive, ``|v|^2.b/2`` a chaque
    bout. Ici la susceptance des lignes n'est PAS dans la branche :
    :mod:`g2egrid.reseau.tables` la repartit par moities dans la capacite ``C`` des
    deux bus, ou elle se confond avec celle des autres lignes aboutissant au
    meme noeud et avec les condensateurs de ``T_Shunts``. Elle n'est donc pas
    reconstituable ligne par ligne depuis le reseau assemble, et les quatre
    transits rendus ici sont ceux de la BRANCHE R-L seule, sans correction :
    ils different d'un transit de load flow de ``|v|^2.b/2`` a chaque bout.
    La verification interne, elle, tient : ``P_from + P_to = -R.|i|^2``, les
    pertes joule de la branche.

    La frequence ``f_busN`` vient de la derivee ANALYTIQUE de l'angle ::

        dphi/dt = (vd.vq' - vq.vd') / |v|^2

    ou ``vd'`` et ``vq'`` sont les equations du noeud, et non une difference
    finie sur la trajectoire : a pas adaptatif, deriver numeriquement un
    angle qui bouge de 1e-6 rad rendrait un bruit de plusieurs hertz. Elle
    est rendue en Hz comme ``f0.(w_repere + dphi/dt / wb)``, ``w_repere``
    etant la pulsation du repere dq vue par le noeud (cf.
    ``Network.repere``) : a l'equilibre ``dphi/dt`` est nul et la mesure vaut
    exactement ``f0``.

    ``fm_busN`` est la meme, passee dans un premier ordre de constante de
    temps une periode -- le reglage d'un relais de frequence, et ce qui rend
    la grandeur lisible, la frequence instantanee d'un modele dq portant tous
    les modes rapides du reseau.

    ``quoi`` restreint aux mesures nommees ; a defaut, tout le catalogue.
    """
    segs = m._segments(tr)
    plans = {id(mm): _plan_mesures(mm) for _i0, _i1, mm, _c in segs}
    plan: dict[str, tuple] = {}
    for _i0, _i1, mm, _c in segs:
        for nom, v in plans[id(mm)].items():
            plan.setdefault(nom, v)
    if quoi is None:
        garde = list(plan)
    else:
        inconnues = [n for n in quoi if n not in plan]
        if inconnues:
            raise KeyError(
                f"mesure(s) inconnue(s) : {inconnues} ; disponibles : {sorted(plan)}"
            )
        garde = list(quoi)

    # une frequence filtree se calcule a partir de la brute : on l'ajoute au
    # travail meme si l'appelant ne l'a pas demandee.
    travail = set(garde) | {"f_" + n[3:] for n in garde if n.startswith("fm_")}
    # le meme capteur peut designer un autre element d'un segment a l'autre
    # (les indices changent quand un element disparait) : un plan par modele
    par_modele: dict[int, dict[int, list[tuple[str, str, int]]]] = {}
    for _i0, _i1, mm, _c in segs:
        par_element: dict[int, list[tuple[str, str, int]]] = {}
        for nom in travail:
            if nom not in plans[id(mm)]:
                continue
            genre, k, pk = plans[id(mm)][nom]
            if genre != "fm":  # calculee apres coup, sur toute la trajectoire
                par_element.setdefault(k, []).append((nom, genre, pk))
        par_modele[id(mm)] = par_element

    out = {nom: np.full(len(tr.t), np.nan) for nom in travail}
    for j, _t, x, z, u, y, mm in m._le_long(tr, uext_fn):
        c = mm.c
        for k, capteurs in par_modele[id(mm)].items():
            el, cm = c.els[k], c.cms[k]
            sx, sz, su, sy = c.ox[k], c.oz[k], c.ou[k], c.oy[k]
            vpoint = None  # dv/dt du noeud, evalue au plus une fois par instant
            for nom, genre, pk in capteurs:
                if genre == "Pdc":
                    iv = el.port_dc_input_idx(pk)
                    ii = el.port_dc_output_idx(pk)
                    out[nom][j] = u[su + iv] * y[sy + ii]
                    continue
                if genre == "Vdc":
                    out[nom][j] = y[sy]  # v_o du noeud continu
                    continue
                if genre in ("P", "Q"):
                    ivd, ivq = el.port_input_idx(pk)
                    vd, vq = u[su + ivd], u[su + ivq]
                    # les courants de port sont les PREMIERES sorties (cf.
                    # Element : h rend d'abord les ports, puis le service)
                    iid, iiq = el.port_output_idx(pk)
                    cd, cq = y[sy + iid], y[sy + iiq]
                    out[nom][j] = (
                        vd * cd + vq * cq if genre == "P" else vq * cd - vd * cq
                    )
                    continue
                vd, vq = y[sy], y[sy + 1]  # vd_o, vq_o du noeud
                if genre == "V":
                    out[nom][j] = math.hypot(vd, vq)
                elif genre == "angle":
                    out[nom][j] = math.degrees(math.atan2(vq, vd))
                else:  # "f"
                    if vpoint is None:
                        vpoint = cm.f(
                            x[sx : sx + cm.n_x],
                            z[sz : sz + cm.n_z],
                            u[su : su + cm.n_u],
                            mm.p[k],
                        )
                    dphi = (vd * vpoint[1] - vq * vpoint[0]) / (vd * vd + vq * vq)
                    out[nom][j] = m.net.f0 * (
                        u[su + el.wref_idx()] + dphi / el.pdict["wb"]
                    )

    tau = 1.0 / m.net.f0  # une periode
    for nom in garde:
        if plan[nom][0] == "fm":
            out[nom] = _premier_ordre(out["f_" + nom[3:]], tr.t, tau)
    return {nom: out[nom] for nom in garde}


# ----------------------------------------------------------------------
# linearisation le long d'une trajectoire
# ----------------------------------------------------------------------
def linearisation(
    m: ModeleNonLineaire, tr: Trajectoire, t: float, uext_fn=None
) -> np.ndarray:
    """La matrice ``A`` INSTANTANEE a l'instant ``t`` d'une trajectoire simulee.

    Le modele d'etat linearise a l'equilibre repond a "ce mode est-il amorti
    autour du point de fonctionnement". Il ne repond pas a "ce mode est-il
    deja dangereux PENDANT le transitoire", parce qu'un grand ecart deplace
    le point de linearisation : angle rotor loin de sa valeur nominale,
    tension effondree, excitation en butee. Cette fonction rend la matrice
    ``A`` la ou la trajectoire se trouve vraiment ; son spectre se lit comme
    celui de :meth:`~g2egrid.reseau.network.Network.build_state_space`, et sur une trajectoire non
    perturbee elle redonne exactement ``ss.A``.

    INSTANT RETENU -- celui des instants rendus le plus proche de ``t``, sans
    interpolation. Interpoler entre deux pas donnerait un etat par lequel le
    systeme n'est jamais passe ; la matrice, elle, varie lentement devant le
    pas d'integration, et imposer ``t_eval`` a
    :meth:`ModeleNonLineaire.simuler` suffit a choisir ou l'on veut pouvoir
    lineariser.
    """
    if not len(tr.t):
        raise ValueError("trajectoire vide : rien a lineariser")
    k = int(np.argmin(np.abs(np.asarray(tr.t, dtype=float) - float(t))))
    t_k = float(tr.t[k])
    # le modele du segment qui contient cet instant (apres une manoeuvre, la
    # matrice A n'a plus la meme taille : c'est celle du reseau d'apres)
    _i0, _i1, mm, cols = next(sg for sg in m._segments(tr) if sg[0] <= k < sg[1])
    uext = mm.uext0() if uext_fn is None else uext_fn(t_k)
    if mm is not m and uext_fn is not None:
        uext = m.traducteur_uext(mm)(uext)
    x = tr.x[k, cols]
    if mm.ideaux:
        raise ValueError(
            "linearisation le long d'une trajectoire : le modele a des noeuds "
            "ideaux, dont la contrainte d'indice 2 n'a pas de matrice d'etat "
            "locale par elimination de z ; lineariser le reseau "
            "(Network.build_state_space), qui reduit les courants en serie"
        )
    z, u = mm.resoudre_algebrique(x, uext, mm.c.z0(), mm.c.u0(), t_k)
    return mm.jacobien_etats(x, uext, z, u)


def reponse_indicielle(
    net: Network,
    entree: str,
    amplitude: float = 0.01,
    duree: float = 1.0,
    n: int = 401,
    t0: float = 0.0,
    **kw,
) -> Trajectoire:
    """Reponse du modele NON LINEAIRE a un echelon de consigne.

    Le pendant non lineaire de
    :func:`g2egrid.multi_niveau.separation.reponse_indicielle`, qui rend la reponse du
    modele d'etat lineaire au meme echelon. Pour une petite amplitude les
    deux doivent coincider ; l'ecart entre elles mesure la non-linearite.
    """
    m = ModeleNonLineaire(net)
    return m.simuler(
        (0.0, duree),
        uext_fn=m.echelon(entree, amplitude, t0),
        t_eval=np.linspace(0.0, duree, n),
        **kw,
    )


def simuler_lineaire(
    m: ModeleNonLineaire,
    t: np.ndarray,
    uext_fn,
    sorties: list[str] | None = None,
    ss=None,
) -> Trajectoire:
    """Simule le modele d'etat LINEARISE sous la meme entree que :meth:`simuler`.

    Le pendant lineaire de :meth:`ModeleNonLineaire.simuler` : memes entrees
    ``uext_fn(t)`` (un scenario, une :class:`Courbe`), meme point de
    fonctionnement, mais la dynamique est celle de ``dx = A.dx + B.du`` --
    ``(A, B)`` de :meth:`Network.build_state_space` -- et non celle des equations
    non lineaires. Ce qu'il rend est donc ``x0 + dx`` et ``y0 + dy``.

    DISCRETISATION EXACTE -- pour un pas ``dt`` constant, l'entree est tenue
    LINEAIRE entre deux instants (bloqueur d'ordre un) et l'etat avance par
    l'exponentielle de la matrice augmentee ``[[A, B, 0], [0, 0, I], [0, 0, 0]]``.
    C'est exact pour une entree affine par morceaux, quelle que soit la raideur
    de ``A`` (ses modes montent a 10^4 s^-1 ici) : aucun critere de stabilite
    n'impose le pas, qui ne sert qu'a echantillonner la courbe. Une simulation
    de dix minutes a 10 ms coute une seconde, contre plusieurs minutes pour le
    modele non lineaire -- c'est ce qui permet de l'utiliser en balayage.

    Les manoeuvres du scenario (ouvertures, defauts) changent le reseau : le
    modele lineaire, lui, reste celui du point de depart. Elles sont ignorees
    ici, et c'est ce que la comparaison avec le non lineaire mesure.
    ``sorties`` : noms de sorties exposees (toutes par defaut).
    """
    from scipy.linalg import expm

    t = np.asarray(t, dtype=float)
    if t.ndim != 1 or len(t) < 2:
        raise ValueError("t : au moins deux instants")
    dts = np.diff(t)
    if np.any(dts <= 0.0) or np.ptp(dts) > 1e-9 * dts.mean():
        raise ValueError("simuler_lineaire demande un pas constant (t = np.linspace)")
    dt = float(dts.mean())
    if ss is None:
        ss = m.net.build_state_space()
    if list(ss.inputs) != list(m.entrees) or list(ss.states) != list(m.etats):
        raise ValueError(
            "le modele d'etat et le modele non lineaire n'ont pas les memes "
            "entrees ou etats (noeuds ideaux, injections ?)"
        )
    n = ss.B.shape[0]
    u0 = m.uext0()
    du = np.array([uext_fn(float(ti)) for ti in t]) - u0  # (nt, nu)
    # on ne garde que les entrees qui bougent : la matrice augmentee reste petite
    actives = np.flatnonzero(np.abs(du).max(axis=0) > 0.0)
    nt = len(t)
    dx = np.zeros((nt, n))
    if len(actives):
        B = ss.B[:, actives]
        k = len(actives)
        M = np.zeros((n + 2 * k, n + 2 * k))
        M[:n, :n] = ss.A
        M[:n, n : n + k] = B
        M[n : n + k, n + k :] = np.eye(k)
        E = expm(M * dt)
        Ad, G0, G1 = E[:n, :n], E[:n, n : n + k], E[:n, n + k :]
        ua = du[:, actives]
        pente = np.diff(ua, axis=0) / dt
        for j in range(nt - 1):
            dx[j + 1] = Ad @ dx[j] + G0 @ ua[j] + G1 @ pente[j]

    x0 = m.c.x0()
    tr = Trajectoire(t=t, x=x0[None, :] + dx, etats=list(m.etats))
    noms = sorties if sorties is not None else list(ss.outputs)
    inconnues = [nm for nm in noms if nm not in ss.outputs]
    if inconnues:
        raise ValueError(
            f"sortie(s) inconnue(s) du modele d'etat : {inconnues} "
            "(un modele d'etat ne porte que les sorties de service)"
        )
    y0 = m.sorties(
        Trajectoire(t=t[:1], x=x0[None, :], etats=list(m.etats)), None, list(noms)
    )
    for nm in noms:
        j = ss.outputs.index(nm)
        tr.sorties[nm] = y0[nm][0] + dx @ ss.C[j] + du @ ss.D[j]
    return tr
