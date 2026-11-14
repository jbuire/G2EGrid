"""Noeud (bus) avec capacite shunt.

C'est ici que se joue la difference d'architecture avec le framework MATLAB.

MATLAB : chaque noeud recevait une RESISTANCE DE CAUSALITE fictive
(``setResistor(100)``) uniquement pour donner au bloc une forme causale
entree->sortie exploitable par l'assemblage LFT. Consequences : shunt non
physique de 0,01 pu a chaque bus, boucle algebrique raide (solveur ode23t
obligatoire) et impedance de source saturant a 100 pu en haute frequence.

Ici : le bus porte sa CAPACITE SHUNT REELLE (charge capacitive des lignes,
capacite de jeu de barres). La tension devient alors un ETAT :

    dvd/dt = (wb/C) . ( isum_d + w.C.vq )
    dvq/dt = (wb/C) . ( isum_q - w.C.vd )

``isum`` est la somme des courants injectes par tous les elements du bus :
c'est la loi des noeuds, exacte, sans element fictif. Le modele global
devient une EDO propre (aucune boucle algebrique tension/courant).

``C`` est la susceptance shunt en pu a la frequence nominale (B = w0.C).
Une valeur faible rend le bus "raide" et introduit des modes rapides : c'est
le comportement EMT physique, pas un artefact.

COHERENCE AVEC LE LOAD FLOW : cette capacite etant un vrai element shunt,
elle DOIT figurer dans le Y-bus (Y[k,k] += j.B). Sinon la loi des noeuds du
load flow (somme des courants = 0) contredirait celle du modele dynamique
(somme des courants = j.w.C.V) et l'equilibre ne se refermerait pas.
``Network`` s'en charge automatiquement.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from g2egrid.noyau.element import Element, InitContext

__all__ = [
    "NOEUDS_DEFAUT",
    "R_CAUSALITE",
    "R_NODE_M",
    "Bus",
    "BusDC",
    "IdealBus",
    "ResistiveBus",
]

#: ce que devient un noeud sans capacite par defaut : ``"ideal"``, la loi des
#: noeuds exacte (FR-005). La resistance de causalite de ``Node.m`` faussait
#: le point de fonctionnement (jusqu'a 77 % de la charge en pertes fictives
#: sur ``ieee33``) ; une resistance assez grande pour ne plus le fausser
#: (1e5 pu) rendait le modele mille fois plus raide et bruitait les valeurs
#: propres proches de zero.
NOEUDS_DEFAUT = "ideal"
#: la resistance de causalite de ``Node.m`` (``setResistor(100)``), pu
#: SYSTEME : celle du mode ``"resistance"``, que gardent les reseaux compares
#: a MATLAB
R_NODE_M = 100.0
#: la resistance du mode ``"resistance"`` quand rien d'autre n'est dit
R_CAUSALITE = R_NODE_M


class Bus(Element):
    """Jeu de barres avec capacite shunt : la tension est un etat."""

    passif = True

    ports = ()  # le bus n'a pas de port : il EST le noeud
    params = ("C", "wb")
    states = ("vd", "vq")
    svc_inputs = ("isum_d", "isum_q")
    svc_outputs = ("vd_o", "vq_o", "vmag")

    def __init__(
        self,
        bus_id: int,
        Un: float = 20e3,
        *,
        btype: int = 3,
        V: float = 1.0,
        theta: float = 0.0,
        Pg: float = 0.0,
        Qg: float = 0.0,
        Pl: float = 0.0,
        Ql: float = 0.0,
        C: float = 0.0,
        R: float | None = None,
        name: str | None = None,
        wb: float | None = None,
    ):
        """Declare un bus, sur le modele d'une ligne de ``Y_DER`` cote MATLAB.

        ``[BusID, Type, V_base, V_pu, Theta_rad, Pg_W, Qg_var, Pl_W, Ql_var]``

        UNITES -- ``Un`` en volts et ``Pg``/``Qg``/``Pl``/``Ql`` en **SI**
        (W et var), comme dans ``Y_DER``. La conversion en pu est faite par
        :meth:`Network.add`, qui seul connait la base ``Sb``. Seule exception :
        ``C``, susceptance shunt deja **en pu** (charge capacitive des lignes
        aboutissant au bus, plus une eventuelle compensation).

        ``btype`` : 1 = slack, 2 = PV, 3 = PQ (convention de ``Y_DER``).

        ``R`` : la resistance de causalite (pu systeme) si ``C`` est nul ;
        ``None`` prend celle du reseau (:attr:`Network.r_causalite`).

        Le bus n'a pas de port au sens des autres elements : il EST le noeud.
        Tant que ``wb`` n'est pas fourni, l'objet n'est qu'une DECLARATION ;
        ``Network.add`` l'acheve en lui donnant la pulsation de base et, si
        ``C`` est nul, en lui substituant un :class:`ResistiveBus`.
        """
        self.name = name if name is not None else f"bus{int(bus_id)}"
        self.buses = [int(bus_id)]
        self.Un = float(Un)
        self.spec_si = {
            "btype": int(btype),
            "V": float(V),
            "theta": float(theta),
            "Pg": float(Pg),
            "Qg": float(Qg),
            "Pl": float(Pl),
            "Ql": float(Ql),
            "C": float(C),
            "R": None if R is None else float(R),
        }
        self.x0 = self.z0 = self.u0 = None
        if wb is not None:
            self.finaliser(wb)

    def finaliser(self, wb: float) -> None:
        """Fige les parametres du modele une fois ``wb`` connu."""
        C = self.spec_si["C"]
        self.pdict = {"C": float(C), "wb": float(wb)}
        self.p = np.array([float(C), float(wb)])

    @property
    def bus(self) -> int:
        return self.buses[0]

    @staticmethod
    def equations(s):
        f = [
            (s.wb / s.C) * (s.isum_d + s.w_ref * s.C * s.vq),
            (s.wb / s.C) * (s.isum_q - s.w_ref * s.C * s.vd),
        ]
        h = [s.vd, s.vq, sp.sqrt(s.vd**2 + s.vq**2)]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        V = ctx.v["_self"]
        # a l'equilibre : isum = +j.w.C.V (dv/dt = 0 dans les equations ci-dessus :
        # isum_d = -w.C.vq, isum_q = +w.C.vd) -- les courants injectes par les
        # elements nourrissent la capacite du bus
        C, w = self.pdict["C"], ctx.w
        isum = 1j * w * C * V
        return {"vd": V.real, "vq": V.imag, "isum_d": isum.real, "isum_q": isum.imag}


class BusDC(Bus):
    """Jeu de barres CONTINU : une capacite, donc un etat de tension.

    Le pendant continu de :class:`Bus`, et le meme raisonnement : la tension
    du noeud est l'etat de sa capacite, et la loi des noeuds ``somme(i) = 0``
    l'alimente. Un seul signal de chaque cote au lieu de deux -- pas d'axe q
    sur un reseau continu.

    HERITER DE ``Bus`` N'EST PAS COSMETIQUE : le cablage repere les noeuds par
    ``isinstance(el, Bus)`` pour construire ``bus_pos`` et pour brancher la loi
    des noeuds. Un noeud continu qui descendrait d'``Element`` obligerait a
    doubler ces deux passages ; ici il les traverse tels quels, et seul le
    GENRE de chaque port decide ensuite combien de signaux y transitent.

    C'EST AUSSI LE CONDENSATEUR DU CONVERTISSEUR. Un VSC a port DC ne porte
    plus son propre etat ``vdc`` : il lit la tension du bus et y renvoie son
    courant. La capacite du bus vaut donc la somme des ``Cdc_pu`` des VSC qui
    s'y raccordent, des demi-capacites des cables et des shunts -- le decalque
    exact de ce que ``Cbus`` fait deja du cote alternatif.

    ``C = 0`` EST INTERDIT : sans capacite, ``dv/dt = (wb/C).isum`` n'est plus
    une equation d'etat mais une equation algebrique ``isum = 0`` SANS
    variable pour la satisfaire -- le cote alternatif s'en tire avec une
    resistance de causalite (:class:`ResistiveBus`), qui n'a pas de sens ici
    puisque rien n'oblige un noeud continu a exister sans capacite.
    """

    model_name = "BusDC"
    dc = True
    passif = True

    #: la tension du noeud ne se deduit pas de ses PROPRES equations : celles-ci
    #: disent seulement que ``isum`` est nul a l'equilibre. C'est la loi des
    #: noeuds du reseau continu, resolue par ``Network._equilibre_dc``, qui la
    #: fixe. Un Newton par element porterait ici sur un jacobien identiquement
    #: nul ; il n'a rien a raffiner.
    raffiner_equilibre = False

    ports = ()
    ports_dc = ()
    params = ("C", "wb")
    states = ("v",)
    svc_inputs = ("isum",)
    svc_outputs = ("v_o",)

    def __init__(
        self,
        bus_id: int,
        Un: float = 640e3,
        *,
        C: float = 0.0,
        V: float = 1.0,
        name: str | None = None,
        wb: float | None = None,
    ):
        """Declare un noeud continu. ``C`` est en pu, ``V`` est le point de depart.

        ``V`` n'est pas une consigne : le reseau continu n'a pas de load flow
        au lot 3, et cette valeur sert d'amorce au Newton de
        ``Network._equilibre_dc``.
        """
        self.name = name if name is not None else f"bus{int(bus_id)}"
        self.buses = [int(bus_id)]
        self.Un = float(Un)
        self.spec_si = {"C": float(C), "V": float(V)}
        self.x0 = self.z0 = self.u0 = None
        if wb is not None:
            self.finaliser(wb)

    @staticmethod
    def equations(s):
        return [(s.wb / s.C) * s.isum], [], [s.v]

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        # isum est nul a l'equilibre : le noeud continu n'a pas de terme de
        # couplage dq a compenser, contrairement a son homologue alternatif.
        return {"v": float(np.real(ctx.v["_self"])), "isum": 0.0}


class ResistiveBus(Bus):
    """Noeud sans capacite : RESISTANCE DE CAUSALITE, comme ``Node.m``.

    Lorsqu'aucune ligne n'amene de charge capacitive au noeud, ``Network`` y
    pose une resistance vers la terre et la tension devient une simple
    relation algebrique

        vd = R . isum_d      vq = R . isum_q

    Le noeud n'a alors AUCUN etat. Ce shunt de 1/R n'est pas physique : il
    n'existe que pour rendre le bloc causal. Contrairement a MATLAB, il est
    aussi injecte dans le Y-bus du load flow (``shunt_admittance``), sans
    quoi l'equilibre ne se refermerait pas et le modele serait linearise hors
    point de fonctionnement.

    C'est le mode ``"resistance"``, celui de MATLAB, que gardent les
    reseaux compares a MATLAB ; le defaut est le noeud ideal
    (:class:`IdealBus`, :data:`NOEUDS_DEFAUT`). R vaut les 100 pu de
    ``setResistor(100)`` (:data:`R_NODE_M`) : 0,01 Sb consommes a chaque
    noeud, jusqu'a 77 % de la charge sur ``ieee33`` (constat FR-005). Une
    autre valeur se regle par ``Network(r_causalite=...)``,
    ``Tables.R_NOEUDS`` ou ``Bus(R=...)`` ; R grand rend le noeud raide
    (valeurs propres en ``-R``).
    """

    model_name = "ResistiveBus"
    params = ("R", "wb")
    states = ()

    def finaliser(self, wb: float) -> None:
        R = self.spec_si["R"]
        self.pdict = {"R": float(R), "wb": float(wb)}
        self.p = np.array([float(R), float(wb)])

    @staticmethod
    def equations(s):
        vd = s.R * s.isum_d
        vq = s.R * s.isum_q
        return [], [], [vd, vq, sp.sqrt(vd**2 + vq**2)]

    def shunt_admittance(self) -> complex:
        return complex(1.0 / self.pdict["R"], 0.0)

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        V = ctx.v["_self"]
        R = self.pdict["R"]
        return {"isum_d": V.real / R, "isum_q": V.imag / R}


class IdealBus(Bus):
    """Noeud SANS capacite ni resistance de causalite : la loi des noeuds exacte.

    La tension du noeud est une variable ALGEBRIQUE, et son equation est la
    loi des noeuds elle-meme ::

        0 = isum_d        0 = isum_q

    Aucun element fictif. Le prix : quand tous les elements du noeud ont leur
    courant pour ETAT -- deux branches R-L en serie, une ligne et un
    transformateur --, la contrainte ne porte que sur des etats (``i1 = i2``)
    et ne dit rien de la tension. Le systeme est d'indice 2 : ``gz`` est
    singulier. C'est :meth:`g2egrid.reseau.network.Network.build_state_space` qui le
    resout, globalement, en derivant la contrainte (la tension du noeud est
    alors celle qui garde les deux courants egaux) puis en eliminant l'etat
    devenu redondant -- la REDUCTION ALGEBRIQUE des courants en serie. Deux
    branches ``R1 + jX1`` et ``R2 + jX2`` en serie donnent exactement le
    spectre d'une branche ``(R1 + R2) + j(X1 + X2)``.

    En temporel, seul BDF2 l'integre (FR-005 (c)) : il resout etats et
    contraintes ENSEMBLE, et sa matrice de Newton reste inversible pour un pas
    fini, indice 2 compris. Radau, qui resout ``g = 0`` a ``x`` fige, ne le
    peut pas : a ``x`` fige, la contrainte ``i1 = i2`` ne dit rien de la
    tension.
    """

    model_name = "IdealBus"
    params = ("wb",)
    states = ()
    algebs = ("vd", "vq")
    #: la tension d'equilibre vient du load flow ; ``0 = isum`` n'a rien a raffiner
    raffiner_equilibre = False

    def finaliser(self, wb: float) -> None:
        self.pdict = {"wb": float(wb)}
        self.p = np.array([float(wb)])

    @staticmethod
    def equations(s):
        g = [s.isum_d, s.isum_q]
        return [], g, [s.vd, s.vq, sp.sqrt(s.vd**2 + s.vq**2)]

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        V = ctx.v["_self"]
        return {"vd": V.real, "vq": V.imag, "isum_d": 0.0, "isum_q": 0.0}
