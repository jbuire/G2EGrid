"""Charge R-L dynamique avec echelon de puissance.

Circuit R-L serie (impedance equivalente calculee a la tension du load flow)
PLUS une injection de courant commandee en puissance, filtree au premier
ordre. Cette seconde voie fournit les entrees ``dp`` / ``dq`` : un echelon de
puissance active/reactive propre, sans changer l'impedance.

Pour une demande dS = dP + j.dQ sous la tension V = vd + j.vq, le courant
appele vaut I = conj(dS/V), soit en dq :

    i_pd = ( dP.vd + dQ.vq ) / (vd^2 + vq^2)
    i_pq = ( dP.vq - dQ.vd ) / (vd^2 + vq^2)

(on verifie P = vd.i_pd + vq.i_pq = dP et Q = vq.i_pd - vd.i_pq = dQ.)
"""

from __future__ import annotations

import numpy as np

from g2egrid.noyau.element import Element, InitContext

__all__ = [
    "ChargeDC_P",
    "ChargeDC_R",
    "Defaut",
    "ILoad",
    "PLoad",
    "RAPPORT_XR_REACTANCE",
    "RLLoad",
    "ReactanceShunt",
    "ZLoad",
]


def _exiger_inductif(classe: str, name: str, Q: float) -> None:
    """Refuse une charge R-L serie dont le reactif n'est pas strictement positif.

    L'impedance equivalente ``Z = V^2/conj(S)`` a pour reactance ``X = V^2.Q /
    |S|^2``, et les equations de la branche s'ecrivent avec ``wb/X``. A ``Q = 0``
    la reactance est nulle (division par zero : NaN, matrice singuliere) ; a
    ``Q < 0`` c'est une inductance NEGATIVE, dont le modele est un systeme
    instable a la place d'une charge capacitive -- une fausse instabilite. Le
    pendant MATLAB (``AdvancedLoad`` mode 'Z', meme ``X = V^2.Q/|S|^2``) n'a pas
    de garde ; ``tables.py`` refuse deja ``Ql <= 0``. Pour un reactif nul ou
    negatif, utiliser :class:`PLoad` ou :class:`ILoad` (pas d'inductance).
    """
    if not Q > 0.0:
        raise ValueError(
            f"{name} ({classe}) : Q = {Q:g} pu. Une charge R-L serie demande Q > 0 "
            "(sa reactance vaut V^2.Q/|S|^2 : nulle a Q = 0, negative donc "
            "instable a Q < 0). Utiliser PLoad ou ILoad, qui acceptent un "
            "reactif nul ou negatif."
        )


class RLLoad(Element):
    """Charge R-L serie + echelon de puissance filtre."""

    passif = True

    ports = ("t",)
    params = ("R", "X", "Tf", "wb")
    states = ("i_d", "i_q", "dp_f", "dq_f")
    svc_inputs = ("dp", "dq")
    svc_outputs = ("p", "q")

    def __init__(
        self, name: str, bus: int, P: float, Q: float, wb: float, Tf: float = 1e-3
    ):
        # R, X sont recalcules dans initialize() a la tension du load flow.
        _exiger_inductif("RLLoad", name, Q)
        super().__init__(name, bus, R=1.0, X=1.0, Tf=Tf, wb=wb)
        self.P, self.Q = float(P), float(Q)  # consommation nominale (pu)

    def scheduled_injection(self) -> complex:
        return -complex(self.P, self.Q)  # une charge injecte du negatif

    @staticmethod
    def equations(s):
        f = [
            (s.wb / s.X) * (s.vd_t - s.R * s.i_d + s.w_ref * s.X * s.i_q),
            (s.wb / s.X) * (s.vq_t - s.R * s.i_q - s.w_ref * s.X * s.i_d),
            (s.dp - s.dp_f) / s.Tf,
            (s.dq - s.dq_f) / s.Tf,
        ]
        v2 = s.vd_t**2 + s.vq_t**2
        i_pd = (s.dp_f * s.vd_t + s.dq_f * s.vq_t) / v2
        i_pq = (s.dp_f * s.vq_t - s.dq_f * s.vd_t) / v2
        # courant total APPELE par la charge
        i_d, i_q = s.i_d + i_pd, s.i_q + i_pq
        h = [
            -i_d,
            -i_q,  # injecte = -appele
            s.vd_t * i_d + s.vq_t * i_q,  # p consomme
            s.vq_t * i_d - s.vd_t * i_q,  # q consomme
        ]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        V = ctx.v["t"]
        S = complex(self.P, self.Q)  # consommee
        if abs(S) < 1e-12:
            raise ValueError(f"{self.name}: charge nulle, impedance indefinie")
        Z = (abs(V) ** 2) / np.conj(S)  # impedance equivalente
        self.set_param("R", Z.real)
        self.set_param("X", Z.imag)
        courant = V / Z  # courant appele
        return {
            "i_d": courant.real,
            "i_q": courant.imag,
            "dp_f": 0.0,
            "dq_f": 0.0,
            "dp": 0.0,
            "dq": 0.0,
        }


class ZLoad(Element):
    """Charge a impedance constante, fidele a ``symLoad`` de G2Elin V1.

    Branche R-L serie pure : deux etats ``i_d``/``i_q``, et rien d'autre.
    L'impedance equivalente est figee a la tension du load flow, comme dans
    la reference V1 ::

        Z = V^2 / conj(S)   ->   R = Re(Z),  X = Im(Z)

    ce qui donne, a V = 1 pu et en base systeme,
    ``R = P/(P^2+Q^2)`` et ``X = Q/(P^2+Q^2)``.

    Difference avec :class:`RLLoad` : celui-ci ajoute une seconde voie
    d'injection commandee en puissance (entrees ``dp``/``dq``), filtree au
    premier ordre, soit **deux etats de plus** qui n'existent pas dans V1.
    Utiliser ``ZLoad`` quand on veut comparer au modele de reference ;
    ``RLLoad`` quand on veut appliquer un echelon de puissance.

    Le pendant MATLAB est ``AdvancedLoad`` en mode 'Z'.
    """

    passif = True

    ports = ("t",)
    params = ("R", "X", "wb")
    states = ("i_d", "i_q")
    svc_outputs = ("p", "q")

    def __init__(self, name: str, bus: int, P: float, Q: float, wb: float):
        # R, X sont recalcules dans initialize() a la tension du load flow.
        _exiger_inductif("ZLoad", name, Q)
        super().__init__(name, bus, R=1.0, X=1.0, wb=wb)
        self.P, self.Q = float(P), float(Q)  # consommation nominale (pu)

    def scheduled_injection(self) -> complex:
        return -complex(self.P, self.Q)  # une charge injecte du negatif

    @staticmethod
    def equations(s):
        f = [
            (s.wb / s.X) * (s.vd_t - s.R * s.i_d + s.w_ref * s.X * s.i_q),
            (s.wb / s.X) * (s.vq_t - s.R * s.i_q - s.w_ref * s.X * s.i_d),
        ]
        h = [
            -s.i_d,
            -s.i_q,  # injecte = -appele
            s.vd_t * s.i_d + s.vq_t * s.i_q,  # p consomme
            s.vq_t * s.i_d - s.vd_t * s.i_q,  # q consomme
        ]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        V = ctx.v["t"]
        S = complex(self.P, self.Q)
        if abs(S) < 1e-12:
            raise ValueError(f"{self.name}: charge nulle, impedance indefinie")
        Z = (abs(V) ** 2) / np.conj(S)
        self.set_param("R", Z.real)
        self.set_param("X", Z.imag)
        courant = V / Z
        return {"i_d": courant.real, "i_q": courant.imag}


#: Rapport X/R d'une reactance shunt dont la table ne donne que B. 50 est
#: l'ordre de grandeur d'une inductance de compensation reelle (facteur de
#: qualite de quelques dizaines a quelques centaines) : assez pour que ses
#: pertes restent negligeables au load flow, assez peu pour que son courant ait
#: un amortissement propre en simulation.
RAPPORT_XR_REACTANCE = 50.0


class ReactanceShunt(Element):
    """Reactance shunt : une branche R-L vers la terre, a admittance constante.

    POURQUOI PAS UNE CAPACITE NEGATIVE -- un condensateur de compensation
    s'ajoute simplement a la capacite du noeud (``T_Shunts`` a ``B > 0``) : son
    courant est deja dans l'equation d'etat du bus, ``C.dv/dt = isum - j.w.C.v``.
    Une reactance ne peut pas s'y fondre. En capacite NEGATIVE, son courant
    croitrait avec la frequence au lieu de decroitre (ce n'est juste qu'a
    50/60 Hz), et des qu'elle depasse la charge des lignes du noeud, la
    capacite totale devient negative : l'equation du bus change de signe, et
    l'assemblage d'avant basculait sans le dire sur un noeud resistif ou la
    reactance disparaissait du modele dynamique, laissant un point de
    fonctionnement qui n'etait plus un equilibre (``|dx/dt| ~ 2.10^4`` sur wscc9).
    F. Kelada fait le meme constat dans son portage : la reactance y est une
    branche R-L dont l'autre bout est la terre, et c'est ce qu'on fait ici.

    ``B`` est la susceptance nominale en pu systeme (negative), ``Q = |B|.V^2``
    a la tension nominale ; ``X = 1/|B|`` et ``R = X/RAPPORT_XR_REACTANCE``.
    Au load flow, c'est une admittance du Y-bus (:meth:`admittance_ybus`) : sa
    consommation suit ``|V|^2``, comme celle d'un condensateur.
    """

    passif = True

    model_name = "ReactanceShunt"
    ports = ("t",)
    params = ("R", "X", "wb")
    states = ("i_d", "i_q")
    svc_outputs = ("p", "q")

    def __init__(
        self, name: str, bus: int, B: float, wb: float, XR: float = RAPPORT_XR_REACTANCE
    ):
        if float(B) >= 0.0:
            raise ValueError(
                f"{name}: B = {B:g} pu ; une reactance a une susceptance NEGATIVE "
                "(un condensateur s'ajoute a la capacite du noeud)"
            )
        X = 1.0 / abs(float(B))
        super().__init__(name, bus, R=X / float(XR), X=X, wb=wb)
        self.B = float(B)

    def admittance_ybus(self) -> complex:
        return 1.0 / complex(self.pdict["R"], self.pdict["X"])

    @staticmethod
    def equations(s):
        f = [
            (s.wb / s.X) * (s.vd_t - s.R * s.i_d + s.w_ref * s.X * s.i_q),
            (s.wb / s.X) * (s.vq_t - s.R * s.i_q - s.w_ref * s.X * s.i_d),
        ]
        h = [
            -s.i_d,
            -s.i_q,  # injecte = -appele
            s.vd_t * s.i_d + s.vq_t * s.i_q,  # p consomme
            s.vq_t * s.i_d - s.vd_t * s.i_q,  # q consomme
        ]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        courant = ctx.v["t"] / complex(self.pdict["R"], self.pdict["X"])  # appele
        return {"i_d": courant.real, "i_q": courant.imag}


class Defaut(Element):
    """Court-circuit a la terre : une branche R-L d'impedance ``R + jX`` (pu).

    Il n'existe qu'en simulation, pose et retire par une manoeuvre
    (``Manoeuvre(genre="defaut")`` puis ``"elimination"``) : il n'a pas de
    point de fonctionnement, et le load flow ne le voit jamais. Son courant
    part de ZERO a l'instant du defaut, comme celui de toute inductance. Un
    defaut franc se modelise par une petite impedance, pas nulle : l'inductance
    donne au courant de defaut sa dynamique, la resistance son amortissement.
    """

    passif = True

    model_name = "Defaut"
    ports = ("t",)
    params = ("R", "X", "wb")
    states = ("i_d", "i_q")
    svc_outputs = ("p", "q")

    def __init__(self, name: str, bus: int, R: float, X: float, wb: float):
        if float(X) <= 0.0:
            raise ValueError(f"{name}: X = {X} ; il faut une inductance de defaut > 0")
        super().__init__(name, bus, R=R, X=X, wb=wb)

    equations = ReactanceShunt.equations

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        return {"i_d": 0.0, "i_q": 0.0}


class DefautDC(Element):
    """Court-circuit POLE A POLE d'un noeud continu : une branche R-L de
    ``R`` et ``X`` (pu systeme, base du noeud), courant appele ``i`` ::

        (X/wb) . di/dt = v - R.i

    Comme :class:`Defaut`, il n'existe qu'en simulation, pose et retire par
    une manoeuvre (``defaut_dc`` puis ``elimination_dc``), et son courant part
    de zero : l'inductance des cables et la sienne donnent au courant de
    defaut sa montee, la decharge des condensateurs du reseau continu son
    premier pic."""

    passif = True

    model_name = "DefautDC"
    ports = ()
    ports_dc = ("t",)
    params = ("R", "X", "wb")
    states = ("i",)
    svc_outputs = ("p",)

    def __init__(self, name: str, bus: int, R: float, X: float, wb: float):
        if float(X) <= 0.0:
            raise ValueError(f"{name}: X = {X} ; il faut une inductance de defaut > 0")
        if float(R) < 0.0:
            raise ValueError(f"{name}: R = {R} ; resistance de defaut negative")
        super().__init__(name, bus, R=R, X=X, wb=wb)

    @staticmethod
    def equations(s):
        f = [(s.wb / s.X) * (s.v_t - s.R * s.i)]
        return f, [], [-s.i, s.v_t * s.i]

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        return {"i": 0.0}


class _CurrentSourceLoad(Element):
    """Base des charges a courant impose (courant constant, puissance constante).

    La charge appelle un courant de consigne ``i_ref(v)`` qui ne depend que de
    la tension a ses bornes ; le courant reel le suit au premier ordre, de
    constante de temps ``Tf``, comme ``AdvancedLoad`` en modes 'I' et 'P' ::

        Tf . di/dt = i_ref(v) - i

    Les entrees ``dp`` / ``dq`` ajoutent un echelon a la consommation nominale.
    Contrairement a une charge R-L, un reactif negatif ne pose aucun probleme :
    le modele n'a pas d'inductance.

    SOUS ``Vbt`` (0,7 pu), LA CHARGE DEVIENT UNE IMPEDANCE -- la convention
    des outils de reseau. Une puissance ou un courant constants n'ont pas de
    sens a tension nulle : le courant de consigne y perd sa direction, et le
    Newton sa solution, des qu'un defaut ou une ile sans formeur fait
    s'effondrer la tension. C'est un MODE discret (``mode_bt``, cf.
    :mod:`g2egrid.elements.limiteurs`) : sous ``Vbt``, la consigne est celle de
    l'impedance qui appelle la meme puissance a ``Vbt`` -- continue a la
    frontiere --, et la charge revient a son modele au-dessus de
    ``Vbt.(1 + HYSTERESIS)``. A drapeau nul, les equations sont exactement
    celles d'avant.
    """

    passif = True

    ports = ("t",)
    params = ("P0", "Q0", "V0", "Tf", "mode_bt", "Vbt")
    modes = ("mode_bt",)
    states = ("i_d", "i_q")
    svc_inputs = ("dp", "dq")
    svc_outputs = ("p", "q")

    def __init__(self, name: str, bus: int, P: float, Q: float, Tf: float = 1e-2):
        # V0 est recale dans initialize() sur la tension du load flow.
        super().__init__(
            name,
            bus,
            P0=float(P),
            Q0=float(Q),
            V0=1.0,
            Tf=Tf,
            mode_bt=0.0,
            Vbt=0.7,
            **self._params_supplementaires(),
        )
        self.P, self.Q = float(P), float(Q)  # consommation nominale (pu)

    def scheduled_injection(self) -> complex:
        return -complex(self.P, self.Q)  # une charge injecte du negatif

    def _params_supplementaires(self) -> dict[str, float]:
        """Parametres propres a une sous-classe (aucun par defaut)."""
        return {}

    @classmethod
    def _consigne_libre(cls, s, P, Q, v2):
        """Le courant de consigne ``(i_ref_d, i_ref_q)`` en fonctionnement
        libre : ``k(v).conj(S/V)``."""
        k = cls._echelle(s)
        # I = conj(S / V) = (P.vd + Q.vq + j.(P.vq - Q.vd)) / |V|^2
        return k * (P * s.vd_t + Q * s.vq_t) / v2, k * (P * s.vq_t - Q * s.vd_t) / v2

    @staticmethod
    def _echelle(s):
        """Facteur k(v) tel que i_ref = conj(S/V) . k(v) ; 1 a la tension V0."""
        raise NotImplementedError

    @staticmethod
    def _echelle_bt(s):
        """Le meme facteur, a ``|V| = Vbt`` : celui de l'impedance de
        repli."""
        raise NotImplementedError

    @staticmethod
    def _puissances(s):
        """La puissance ``(P, Q)`` demandee : nominale plus echelon. Une charge
        pilotee par une consigne ABSOLUE (:class:`ChargeACPref`) la surcharge."""
        return s.P0 + s.dp, s.Q0 + s.dq

    @classmethod
    def equations(cls, s):
        P, Q = cls._puissances(s)
        v2 = s.vd_t**2 + s.vq_t**2
        i_ref_d, i_ref_q = cls._consigne_libre(s, P, Q, v2)
        # sous Vbt : l'impedance qui appelle la meme puissance a Vbt
        kz = cls._echelle_bt(s) / s.Vbt**2
        i_ref_d = (1 - s.mode_bt) * i_ref_d + s.mode_bt * kz * (P * s.vd_t + Q * s.vq_t)
        i_ref_q = (1 - s.mode_bt) * i_ref_q + s.mode_bt * kz * (P * s.vq_t - Q * s.vd_t)
        f = [(i_ref_d - s.i_d) / s.Tf, (i_ref_q - s.i_q) / s.Tf]
        h = [
            -s.i_d,
            -s.i_q,  # injecte = -appele
            s.vd_t * s.i_d + s.vq_t * s.i_q,  # p consomme
            s.vq_t * s.i_d - s.vd_t * s.i_q,  # q consomme
        ]
        return f, [], h

    def _recaler_direction(self, V: complex) -> None:
        """Hook : recaler sur ``V`` les parametres propres a la sous-classe."""

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        V = ctx.v["t"]
        if abs(V) < 1e-9:
            raise ValueError(f"{self.name}: tension nulle au bus {self.buses[0]}")
        self.set_param("V0", abs(V))
        self._recaler_direction(V)
        courant = np.conj(complex(self.P, self.Q) / V)
        return {"i_d": courant.real, "i_q": courant.imag, "dp": 0.0, "dq": 0.0}

    def gardes(self, cm, x, z, u, p) -> list[tuple[str, float, float, int]]:
        """Entree en impedance sous ``Vbt``, sortie au-dessus de
        ``Vbt.(1 + HYSTERESIS)``."""
        from g2egrid.elements.limiteurs import HYSTERESIS

        noms = self.input_names()
        v = float(np.hypot(u[noms.index("vd_t")], u[noms.index("vq_t")]))
        vbt = p[self.params.index("Vbt")]
        if p[self.params.index("mode_bt")] == 0.0:
            return [("mode_bt", vbt - v, 1.0, +1)]
        return [("mode_bt", v - vbt * (1 + HYSTERESIS), 0.0, +1)]


class ILoad(_CurrentSourceLoad):
    """Charge a courant constant : ``|I|`` et facteur de puissance fixes.

    La puissance appelee varie comme ``|V|`` : ``S(V) = S0 . |V| / V0``, ou
    ``V0`` est la tension du load flow. Pendant MATLAB : ``AdvancedLoad`` en
    mode 'I'.

    COURANT REACTIF DE DIRECTION FIXE. Un courant de module ET de facteur de
    puissance fixes, ``I = conj(S0/V)|V|/V0``, tourne avec l'angle de V : son
    admittance incrementale ``dI/dV = [[0, Q], [0, P]]/V0^2`` n'est pas
    symetrique -- sa partie symetrique ``[[0, Q/2], [Q/2, P]]`` a un
    determinant ``-Q^2/4 < 0``, donc une direction de conductance NEGATIVE
    (active), et sa partie antisymetrique est une susceptance. Filtree par
    ``Tf`` avec le reste du courant, cette admittance retarde la susceptance :
    ``j.a/(1 + j.w.Tf)`` a pour partie reelle ``a.w.Tf/(1 + (w.Tf)^2)``, de
    signe de ``w``. Le seul amortissement d'un reseau etant ``R/(2L)`` (de
    l'ordre de la seconde^-1), tout mode LC passait instable (Re +2 a
    +6000 s^-1 selon ``C`` et ``Q``) des que ``Q != 0``. Retarder moins ne
    corrige rien -- une direction qui suit V sans retard est PIRE (Re jusqu'a
    +4e6 s^-1, essaye) : c'est l'admittance non passive elle-meme qui est en
    cause, pas le filtre. ``PLoad`` a une admittance symetrique et n'a pas ce
    defaut.

    Correction : seule la partie ACTIVE du courant suit l'angle de V ; la partie
    REACTIVE garde la direction de V au load flow (``cphi0``, ``sphi0``) ::

        I = (P.v/|v| - j.Q.e^{j.phi0}) / V0 ,    dI/dV = [[0, 0], [0, P]]/V0^2

    admittance symetrique semi-definie positive : passive. Inchange : le
    point d'equilibre, la puissance ``P(|V|) = P0.|V|/V0``, et tout a
    ``Q = 0``. Change : pour un ecart d'ANGLE de V, le reactif appele n'est
    plus a facteur de puissance constant (``Q ~ Q0.cos(phi - phi0)``) --
    c'est le prix d'un modele passif.
    """

    params = _CurrentSourceLoad.params + ("cphi0", "sphi0")

    def _params_supplementaires(self) -> dict[str, float]:
        return dict(cphi0=1.0, sphi0=0.0)  # recales par _recaler_direction

    def _recaler_direction(self, V: complex) -> None:
        self.set_param("cphi0", V.real / abs(V))
        self.set_param("sphi0", V.imag / abs(V))

    @classmethod
    def _consigne_libre(cls, s, P, Q, v2):
        """Module constant, partie ACTIVE qui suit l'angle de V, partie
        REACTIVE de direction fixe (celle de V au load flow) ::

            I = (P.v/|v| - j.Q.e^{j.phi0}) / V0

        A l'equilibre, c'est ``conj(S/V)`` ; la puissance appelee reste
        ``S0.|V|/V0`` en module de V."""
        mod = v2**0.5
        return (P * s.vd_t / mod + Q * s.sphi0) / s.V0, (
            P * s.vq_t / mod - Q * s.cphi0
        ) / s.V0

    @staticmethod
    def _echelle(s):
        return (s.vd_t**2 + s.vq_t**2) ** 0.5 / s.V0

    @staticmethod
    def _echelle_bt(s):
        return s.Vbt / s.V0


class PLoad(_CurrentSourceLoad):
    """Charge a puissance constante : ``S = P + jQ`` quelle que soit la tension.

    Le courant appele varie comme ``1/|V|``. Pendant MATLAB : ``AdvancedLoad``
    en mode 'P'.
    """

    @staticmethod
    def _echelle(s):
        return 1.0

    @staticmethod
    def _echelle_bt(s):
        return 1.0


class ChargeDC_R(Element):
    """Charge continue resistive : ``i = v/R``, sans aucun etat.

    Le seul element continu PUREMENT algebrique du lot. Il n'a besoin ni
    d'etat ni de variable algebrique : sa sortie est une fonction explicite de
    son entree, et ``Network.solve_equilibrium`` la saute sans Newton.

    Attention -- une charge resistive ne declare PAS ``shunt_admittance`` :
    ce shunt-la part dans le Y-bus du load flow ALTERNATIF, ou son bus n'a
    aucun indice.
    """

    passif = True

    model_name = "ChargeDC_R"
    ports = ()
    ports_dc = ("t",)
    params = ("R",)
    svc_outputs = ("p",)

    def __init__(self, name: str, bus: int, R: float):
        if float(R) <= 0.0:
            raise ValueError(f"{name}: R = {R} ; resistance de charge non physique")
        super().__init__(name, bus, R=R)

    @staticmethod
    def equations(s):
        i = s.v_t / s.R  # courant APPELE
        return [], [], [-i, s.v_t * i]

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        return {}


class ChargeDC_P(Element):
    """Charge continue a puissance constante, courant filtre au premier ordre.

    Le pendant continu de :class:`PLoad` ::

        Tf . di/dt = (P0 + dp)/v_t - i

    ``dp`` est l'entree d'echelon : elle permet de perturber la liaison du
    cote continu sans toucher au reseau alternatif, ce qui est precisement ce
    qu'on veut pour exciter les modes d'une MTDC.

    Le filtre n'est pas cosmetique : sans lui, ``i = P/v`` serait une relation
    algebrique a transmission directe de la tension vers le courant, et la
    fermeture LFT du reseau (``I - Ry.Hu``) verrait une boucle algebrique
    tension/courant la ou le reste du paquet n'en a aucune.
    """

    passif = True

    model_name = "ChargeDC_P"
    ports = ()
    ports_dc = ("t",)
    params = ("P0", "Tf")
    states = ("i",)
    svc_inputs = ("dp",)
    svc_outputs = ("p",)

    def __init__(self, name: str, bus: int, P: float, Tf: float = 1e-2):
        super().__init__(name, bus, P0=float(P), Tf=Tf)
        self.P = float(P)  # consommation nominale (pu base systeme)

    def scheduled_injection_dc(self) -> float:
        return -self.P  # une charge injecte du negatif

    @staticmethod
    def equations(s):
        i_ref = (s.P0 + s.dp) / s.v_t
        f = [(i_ref - s.i) / s.Tf]
        h = [-s.i, s.v_t * s.i]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        v = float(np.real(ctx.v["t"]))
        if abs(v) < 1e-9:
            raise ValueError(f"{self.name}: tension nulle au bus {self.buses[0]}")
        return {"i": self.P / v, "dp": 0.0}
