"""Blocs elementaires d'un centre de donnees (data center).

La chaine electrique d'un centre de donnees alimente par le reseau MT est ::

    reseau -- transfo -- REDRESSEUR A DIODES -- bus continu -- ONDULEUR V/f -- CHARGE
                                                   |
                                          DC/DC BATTERIE (en parallele)

et, en secours, un groupe electrogene (:mod:`g2egrid.elements.groupe_electrogene`)
se raccorde au meme jeu de barres MT. Chaque bloc est un :class:`~g2egrid.Element`
ordinaire : il declare ses equations symboliques, et c'est le meme modele qui
sert a la simulation NON LINEAIRE (``ModeleNonLineaire``) et au modele d'etat
LINEARISE (``Network.build_state_space``). Le transformateur est le
``Transformer`` du noyau, le condensateur alternatif est la capacite du bus
(``add_bus(..., C=...)``) et le condensateur continu celle du ``BusDC``.

====================  =====================================================
Classe                Role
====================  =====================================================
RedresseurDiode       pont de diodes triphase, modele MOYEN dq (1 etat)
OnduleurVF            onduleur de tension a V et f CONSTANTS, bus continu
BatterieDCDC          DC/DC + batterie (etat de charge), statisme de tension
ChargeACPref          charge alternative a puissance constante, entree p_ref
====================  =====================================================

Les blocs se placent en base systeme ``Sb`` ; un bloc de puissance ``Sn``
ramene son courant par ``kS = Sn/Sb``.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from g2egrid.elements.gfm import GFM3DC, _GFMBase, _GFMDCBase
from g2egrid.elements.load import PLoad
from g2egrid.noyau.element import Element, InitContext

__all__ = [
    "BatterieDCDC",
    "ChargeACPref",
    "OnduleurVF",
    "RedresseurDiode",
    "reactif_commutation",
]


# ----------------------------------------------------------------------
# redresseur a diodes
# ----------------------------------------------------------------------
class RedresseurDiode(Element):
    """Pont de diodes triphase : modele MOYEN dans le repere dq.

    POURQUOI UN MODELE MOYEN -- un pont de diodes commute six fois par
    periode ; le repere dq, lui, ne voit que le FONDAMENTAL. Deux relations
    exactes du pont ideal a conduction continue le ferment :

    * cote continu, la tension moyenne de sortie est proportionnelle a
      l'amplitude de la tension alternative, diminuee de la chute de
      commutation, qui se comporte comme une resistance (non dissipative :
      c'est de l'energie echangee avec l'inductance de source) ::

          Vr = kv.|V| ,    v_pont = Vr - Rcom.idc

      En grandeurs physiques ``Rcom = 3.Xc/pi`` [ohm]. En pu, avec pour base
      continue la tension redressee a vide ``3.sqrt(2)/pi.U`` et la meme
      puissance de base que cote alternatif, l'impedance de base continue vaut
      ``18/pi^2`` fois l'alternative, d'ou ``Rcom = pi.Xcom/6`` avec ``Xcom``
      en pu alternatif (une premiere version prenait ``3.Xcom/pi``, la formule
      physique appliquee a des pu : 1,8 fois trop). ``Xcom`` est exprime dans
      la base alternative DU PONT, celle qui fait de ``3.sqrt(2)/pi`` fois sa
      tension la base continue ; ``kv`` joue le role d'un transformateur ideal
      entre le bus et le pont. Exprime dans la base du bus, la reactance
      vaudrait ``Xcom/kv^2`` (``kv`` est cale a 1,02 - 1,03 dans le reseau
      type : 5 % d'ecart) ;

    * cote alternatif, le fondamental du courant est EN PHASE avec la tension
      (au deplacement ``tphi`` pres, du a la commutation) et la puissance est
      celle que la diode laisse passer sans perte : ``P = v_pont.idc``.

    Le courant continu est un ETAT, celui de l'inductance de lissage ``Ld`` en
    serie avec les deux inductances de commutation qui conduisent (``2.Lc``,
    soit ``pi^2.Xcom/9`` en pu continu) ::

        (Lt/wb) . didc/dt = Vr - (Rdc + Rcom).idc - v_dc ,   Lt = Ld + pi^2.Xcom/9

    et le courant alternatif en est une fonction ALGEBRIQUE de ``idc`` et de la
    tension du bus, sans filtre : les diodes commutent sur la tension
    instantanee, leur courant en suit la phase sans retard. Un filtre du
    premier ordre sur ce courant (essaye d'abord) retarde ce suivi et rend
    instable, vers 500 Hz, la resonance entre la capacite du bus alternatif et
    l'inductance du transformateur : de l'energie est renvoyee au bon moment
    pour l'entretenir. Cette transmission directe tension vers courant
    suppose donc un bus alternatif AVEC capacite (la capacite alternative du
    redresseur, ``add_bus(..., C=...)``) ; sur un noeud ideal sans capacite,
    elle ferait une boucle algebrique (cf. ``ChargeDC_P``).

    Le deplacement du fondamental du courant vient de la commutation : pour
    un pont a diodes (angle de retard nul), l'empietement ``mu`` verifie
    ``v_pont/Vr = Vd/Vd0 = (1 + cos(mu))/2`` et ``tan(phi) = (2.mu -
    sin(2.mu))/(1 - cos(2.mu))`` (forme exacte, Kimbark ; l'approximation
    ``cos(phi) ~ Vd/Vd0`` surestime ``tan(phi)`` de 6 %). Le modele garde
    ``Q = tphi.P`` avec ``tphi`` cale au point de fonctionnement ;
    :func:`reactif_commutation` donne la valeur physique, et
    ``g2egrid.datacenter.reseau_datacenter`` l'impose au load flow. Le
    linearise a donc ``dQ/dP = tphi`` et non la pente physique (plus forte
    d'environ 40 % : l'empietement croit avec le courant).

    DIODES BLOQUEES -- un pont ne conduit que dans un sens : si la tension
    continue depasse ``Vr`` le courant s'annule. Plutot qu'un mode discret
    (blocage/conduction, comme les limiteurs), dont les basculements repetes
    quand ``Vr`` et ``v_dc`` sont proches ont ete mesures couteux -- des
    dizaines de milliers de resolutions en 30 ms --, la conduction est une
    caracteristique LISSE : le courant continu cherche ``max(0, Vr - v_dc)/Rd``,
    le maximum etant lisse sur ``eps_v`` ::

        (Lt/wb) . didc/dt = g(Vr - v_dc) - Rd.idc ,   g(x) = (x + sqrt(x^2 + eps_v^2))/2

    Pont conducteur (``x >> eps_v``), ``g(x) = x + eps_v^2/(4x)`` : le modele
    est celui d'une inductance en serie avec ``Rd`` a ``eps_v^2/(4x)`` pres
    (0,1 % de la chute ``Rd.idc`` a pleine charge pour ``eps_v = 1e-3`` : le
    rayon est choisi pour cela), et son linearise celui du pont conducteur. Pont bloque, ``g``
    tombe en ``eps_v^2/(4|x|)`` et le courant s'eteint en ``Lt/(wb.Rd)`` (une
    dizaine de ms), sans evenement. Pas de conduction discontinue (courant
    intermittent) : c'est l'hypothese du modele moyen, valable tant que le courant
    continu reste au-dessus de quelques % du nominal.

    CALAGE -- ``kv`` et ``tphi`` ne sont pas des donnees de depart mais le
    resultat du dimensionnement : ``initialize`` les cale pour que le pont
    delivre au point du load flow exactement la tension continue du bus et le
    reactif qu'il consomme (comme ``RLLoad`` recale son impedance). Le pont est
    donc un VSC pour le load flow mixte (mode ``VdcQ`` : il tient la tension
    continue et son reactif ; sa puissance active sort du bilan du noeud).

    Pas d'angle de commande, pas de PLL : un redresseur a diodes est un
    element PASSIF du point de vue de la commande, ce qui en fait un
    consommateur de puissance purement reseau.
    """

    ports = ("t",)
    ports_dc = ("dc",)
    modes_dc = ("PQ", "VdcQ")
    params = (
        "kS",
        "kv",
        "Rdc",
        "Xcom",
        "Ld",
        "tphi",
        "wb",
        "eps_v",
    )
    states = ("idc",)
    svc_outputs = ("p", "q", "vr", "idc_out")

    def __init__(
        self,
        name: str,
        bus: int,
        bus_dc: int,
        Sn: float,
        Sb: float,
        wb: float,
        Rdc: float = 0.005,
        Xcom: float = 0.03,
        Ld: float = 0.1,
        eps_v: float = 1e-3,
    ):
        """``Rdc`` : resistance continue (pertes cuivre + diodes), ``Xcom`` :
        reactance de commutation (inductance de source vue du pont), ``Ld`` :
        inductance de lissage, tous en pu de ``Sn`` ; ``eps_v`` : rayon de
        lissage de la conduction [pu de tension continue]."""
        if float(Rdc) <= 0.0:
            raise ValueError(f"{name}: Rdc = {Rdc} ; resistance continue non physique")
        if float(Ld) <= 0.0 or float(eps_v) <= 0.0:
            raise ValueError(f"{name}: Ld et eps_v doivent etre > 0")
        super().__init__(
            name,
            [bus, int(bus_dc)],
            kS=Sn / Sb,
            kv=1.0,
            Rdc=Rdc,
            Xcom=Xcom,
            Ld=Ld,
            tphi=0.0,
            wb=wb,
            eps_v=eps_v,
        )
        self.Sn, self.Sb = float(Sn), float(Sb)

    @staticmethod
    def equations(s):
        V2 = s.vd_t**2 + s.vq_t**2
        V = sp.sqrt(V2)
        Rcom = sp.pi * s.Xcom / 6
        Lt = s.Ld + sp.pi**2 * s.Xcom / 9
        Vr = s.kv * V
        x = Vr - s.v_dc
        didc = (s.wb / Lt) * (
            (x + sp.sqrt(x**2 + s.eps_v**2)) / 2 - (s.Rdc + Rcom) * s.idc
        )
        P = (Vr - Rcom * s.idc) * s.idc  # puissance alternative consommee
        Q = s.tphi * P
        ref_d = (P * s.vd_t + Q * s.vq_t) / V2
        ref_q = (P * s.vq_t - Q * s.vd_t) / V2
        h = [
            -s.kS * ref_d,  # alternatif : injecte = -appele
            -s.kS * ref_q,
            s.kS * s.idc,  # continu : le pont ALIMENTE le noeud
            s.kS * P,  # p, q consommes (base Sb)
            s.kS * Q,
            Vr,
            s.idc,
        ]
        f = [didc]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        p = self.pdict
        V = ctx.v["t"]
        vdc = float(np.real(ctx.v["dc"]))
        S = ctx.s.get("t", 0j)  # puissance alternative INJECTEE (< 0 : consommee)
        P = -complex(S).real / p["kS"]
        Q = -complex(S).imag / p["kS"]
        if abs(V) < 1e-9 or vdc < 1e-9:
            raise ValueError(f"{self.name}: tension nulle a un port du redresseur")
        # P = (Vr - Rcom.idc).idc, avec Vr = vdc + x et x = g^-1(Rd.idc) :
        # g(x) = s => x = s - eps^2/(4 s), d'ou (vdc + Rdc.idc).idc - eps^2/(4 Rd) = P.
        # Un P negatif n'a pas de sens pour un pont de diodes ; il n'apparait que
        # dans les essais du Newton du load flow (differences finies autour de
        # zero), ou la continuation de la formule est ce qu'il lui faut.
        Rdc, eps = p["Rdc"], p["eps_v"]
        Rcom = np.pi * p["Xcom"] / 6
        Rd = Rdc + Rcom
        discr = vdc**2 + 4.0 * Rdc * (P + eps**2 / (4.0 * Rd))
        if discr < 0.0:
            raise ValueError(
                f"{self.name}: puissance {P:.4g} pu impossible a la tension "
                f"continue {vdc:.4g} pu (resistance Rdc trop forte)"
            )
        idc = (-vdc + np.sqrt(discr)) / (2.0 * Rdc)
        g = max(Rd * idc, 1e-12)
        Vr = vdc + g - eps**2 / (4.0 * g)
        self.set_param("kv", Vr / abs(V))
        self.set_param("tphi", Q / P if P > 1e-9 else 0.0)
        return {"idc": idc}


def reactif_commutation(red: RedresseurDiode) -> float:
    """Le reactif [pu systeme] que le pont consomme a son point d'equilibre,
    d'apres le deplacement du fondamental du courant du a la commutation.

    Pont a diodes (angle de retard nul), empietement ``mu`` ::

        Vd = Vd0.(1 + cos(mu))/2     tan(phi1) = (2.mu - sin(2.mu))/(1 - cos(2.mu))

    (Kimbark, *Direct current transmission*, ch. 3 ; la forme exacte du
    fondamental, que ``tests/test_datacenter.py`` verifie par serie de Fourier
    du courant de phase). ``red`` doit avoir son equilibre."""
    p = red.pdict
    idc = float(red.x0[red.states.index("idc")])
    noms = red.input_names()
    V = float(np.hypot(red.u0[noms.index("vd_t")], red.u0[noms.index("vq_t")]))
    Vr = p["kv"] * V
    v_pont = Vr - np.pi * p["Xcom"] / 6 * idc
    mu = float(np.arccos(np.clip(2.0 * v_pont / Vr - 1.0, -1.0, 1.0)))
    if mu < 1e-9:
        return 0.0
    tphi = (2.0 * mu - np.sin(2.0 * mu)) / (1.0 - np.cos(2.0 * mu))
    return float(p["kS"] * v_pont * idc * tphi)


# ----------------------------------------------------------------------
# DC/DC + batterie
# ----------------------------------------------------------------------
class BatterieDCDC(Element):
    """Convertisseur DC/DC et batterie, en parallele sur un bus continu.

    Vu du bus continu, le convertisseur est une source de COURANT dont la
    boucle interne est du premier ordre (``Tc``). Sa consigne combine ::

        i_libre = p_ref/v + Kv.g(v_ref - v) + Kh.(vf - v) ,   i_ref = sat(i_libre)

    * ``p_ref`` : puissance programmee (pu systeme, > 0 = decharge) ;
    * ``Kv.g(v_ref - v)`` : STATISME UNILATERAL de tension. ``g`` est le
      maximum lisse de zero, ``g(x) = (x + sqrt(x^2 + eps^2))/2`` : la batterie
      ne soutient le bus que quand ``v`` passe SOUS ``v_ref``, et reste inerte
      au-dessus. Un statisme bilateral partagerait durablement la charge avec le
      redresseur et videlait la batterie. ``v_ref`` se place sous la tension
      de service du bus : batterie en veille tant que le reseau est la, prise de
      charge a la disparition du redresseur. RELEVER ``v_ref`` au-dessus de la
      tension du redresseur fait prendre la charge a la batterie, l'ABAISSER
      progressivement la rend au redresseur : c'est la reprise progressive du
      groupe electrogene ;
    * ``Kh.(vf - v)`` : statisme sur la composante RAPIDE de la tension, ``vf``
      etant la tension continue filtree (constante ``Thp``). A l'equilibre
      ``vf = v`` : le terme s'annule, la batterie ne fournit aucune puissance
      moyenne et n'absorbe que les fluctuations rapides de la charge. C'est le
      LISSAGE.

    A l'equilibre le courant vaut ``P / v`` exactement : ``initialize`` cale
    ``p_ref`` pour annuler le courant residu de ``g`` (de l'ordre de
    ``Kv.eps^2/(4|x|)``) et pour compenser la saturation douce (la consigne
    libre est ``sat^-1(P/v)``), de sorte que la batterie ne se decharge pas en
    veille et que son etat de charge est un vrai point d'equilibre.

    La batterie porte son ETAT DE CHARGE ``soc`` (0 a 1) ::

        dsoc/dt = - P_bat / Eb ,   P_bat = p_dc + Rint.p_dc^2/Voc^2

    ``Eb`` est l'energie utile en pu systeme x seconde (``E[Wh].3600/Sb``),
    ``Rint`` la resistance interne en pu, ``Voc`` la tension a vide
    ``Voc0.(1 + dVoc.(soc - soc0))``. La tension de la batterie n'entre pas
    dans le bus continu : le DC/DC la decouple. ``soc`` est un integrateur pur :
    le modele linearise a un mode nul (``lambda ~ 0``), l'etat de charge n'a pas
    de point d'equilibre propre.

    Le courant est borne par une saturation DOUCE, ``Imax.r/(1 + r^4)^(1/4)`` avec
    ``r = i_libre/Imax`` : quasi lineaire jusqu'a 50 % de ``Imax`` (1,5 %
    d'ecart), derivable partout, de pente 1 a l'origine. Quand le courant du point de
    fonctionnement s'approche de ``Imax``, le gain de la boucle de tension du modele
    linearise est donc inferieur a ``Kv``.
    """

    ports = ()
    ports_dc = ("dc",)
    params = (
        "Kv",
        "Kh",
        "Thp",
        "Tc",
        "Imax",
        "Eb",
        "Rint",
        "Voc0",
        "dVoc",
        "soc0",
        "eps_v",
    )
    states = ("i", "vf", "soc")
    svc_inputs = ("p_ref", "v_ref")
    svc_outputs = ("p", "soc_out", "p_bat")

    def __init__(
        self,
        name: str,
        bus_dc: int,
        Kv: float = 20.0,
        Kh: float = 0.0,
        Thp: float = 1.0,
        Tc: float = 5e-4,
        Imax: float = 0.5,
        Eb: float = 600.0,
        Rint: float = 0.01,
        Voc0: float = 1.0,
        dVoc: float = 0.1,
        soc0: float = 0.9,
        eps_v: float = 3e-3,
        P: float = 0.0,
        vref: float = 0.92,
    ):
        """``Eb`` : energie utile [pu systeme . s] ; ``P`` : puissance
        programmee a l'etat de depart ; ``vref`` : consigne de tension ;
        ``eps_v`` : rayon de lissage de ``g``."""
        if float(Tc) <= 0.0 or float(Thp) <= 0.0 or float(eps_v) <= 0.0:
            raise ValueError(f"{name}: Tc, Thp et eps_v doivent etre > 0")
        if float(Eb) <= 0.0 or float(Imax) <= 0.0:
            raise ValueError(f"{name}: Eb et Imax doivent etre > 0")
        if not 0.0 < float(soc0) < 1.0:
            raise ValueError(f"{name}: soc0 = {soc0} hors de ]0, 1[")
        super().__init__(
            name,
            [int(bus_dc)],
            Kv=Kv,
            Kh=Kh,
            Thp=Thp,
            Tc=Tc,
            Imax=Imax,
            Eb=Eb,
            Rint=Rint,
            Voc0=Voc0,
            dVoc=dVoc,
            soc0=soc0,
            eps_v=eps_v,
        )
        self.P, self.vref = float(P), float(vref)

    #: la batterie ne tient pas la tension de son ile continue (statisme)
    reference_tension_dc = False

    @staticmethod
    def equations(s):
        v = s.v_dc
        x = s.v_ref - v
        g = (x + sp.sqrt(x**2 + s.eps_v**2)) / 2  # max(0, x) lisse
        i_libre = s.p_ref / v + s.Kv * g + s.Kh * (s.vf - v)
        r = i_libre / s.Imax
        i_ref = s.Imax * r / (1 + r**4) ** sp.Rational(1, 4)  # saturation douce
        p_dc = v * s.i
        voc = s.Voc0 * (1 + s.dVoc * (s.soc - s.soc0))
        p_bat = p_dc + s.Rint * p_dc**2 / voc**2
        f = [
            (i_ref - s.i) / s.Tc,
            (v - s.vf) / s.Thp,
            -p_bat / s.Eb,
        ]
        return f, [], [s.i, p_dc, s.soc, p_bat]

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        v = float(np.real(ctx.v["dc"]))
        p = self.pdict
        x = self.vref - v
        g = 0.5 * (x + np.hypot(x, p["eps_v"]))
        # p_ref annule le courant residu du statisme et compense la saturation
        # douce : la consigne libre est sat^-1(P/v), le courant vaut P/v exactement
        # (sat^-1 : r -> r/(1 - r^4)^(1/4), pour |r| < 1)
        r = (self.P / v) / p["Imax"]
        if abs(r) >= 1.0:
            raise ValueError(
                f"{self.name}: puissance {self.P:.4g} pu au-dela de Imax.v "
                f"({p['Imax'] * v:.4g} pu)"
            )
        i_libre = p["Imax"] * r / (1.0 - r**4) ** 0.25
        p_ref = v * i_libre - v * p["Kv"] * g
        return {
            "i": self.P / v,
            "vf": v,
            "soc": p["soc0"],
            "p_ref": p_ref,
            "v_ref": self.vref,
        }


# ----------------------------------------------------------------------
# onduleur V/f constants
# ----------------------------------------------------------------------
class OnduleurVF(GFM3DC):
    """Onduleur de tension a TENSION ET FREQUENCE CONSTANTES, bus continu.

    C'est un grid-forming sans statisme : la frequence est une consigne
    ``w_set`` (1 pu = la frequence nominale) et le module de la tension aux
    bornes du condensateur une consigne ``ve_ref``. Aucun terme P-f ni Q-V
    (``mp = nq = 0``) ; la puissance qu'il fournit est celle que sa charge
    appelle, et elle est tiree du bus continu -- c'est pour le bus continu une
    charge a puissance constante, le comportement qui rend ce bus fragile.

    Meme filtre LCL, memes boucles de tension et de courant que ``GFM3DC``
    (13 etats), dont il herite aussi les limitations de courant et de
    modulation (qui s'activent des que la tension continue s'effondre) ; seule
    differe l'equation de l'angle ::

        dtheta/dt = wb.(w_set - w_ref)

    L'angle ne depend de rien d'autre que du repere : l'ilot alimente par
    l'onduleur est INDEPENDANT de la frequence du reseau amont, et le modele a
    un mode de jauge (``lambda = 0``) tant que rien d'autre ne l'ancre, comme
    tout reseau sans bus infini.

    Au load flow mixte, l'onduleur forme son ilot alternatif : il en est le
    slack (mode ``Vtheta``).
    """

    svc_inputs = GFM3DC.svc_inputs + ("w_set",)

    def __init__(self, name, bus, bus_dc, Sn, Sb, wb, **params):
        for k in ("mp", "nq"):
            if params.get(k, 0.0) != 0.0:
                raise ValueError(f"{name}: un onduleur V/f n'a pas de statisme ({k})")
            params[k] = 0.0
        super().__init__(name, bus, bus_dc, Sn, Sb, wb, **params)

    @staticmethod
    def equations(s):
        w = s.w_set
        lcl, inner, (p, q, md, mq, iD, iQ, Vt) = _GFMBase._common(
            s, w, _GFMDCBase._vdc_modulation(s)
        )
        f = (
            lcl
            + [s.wf_p * (p - s.pm), s.wb * (w - s.w_ref), s.wf_q * (q - s.qm)]
            + inner
        )
        h = [iD, iQ, _GFMDCBase._courant_dc(s, md, mq), p, q, w, Vt, s.theta, s.v_dc]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = super().initialize(ctx)
        init["w_set"] = ctx.w
        return init


# ----------------------------------------------------------------------
# charge alternative pilotee en puissance
# ----------------------------------------------------------------------
class ChargeACPref(PLoad):
    """Charge alternative a puissance constante, pilotee par ``p_ref``.

    La charge informatique d'un centre de donnees : elle appelle la puissance
    ``p_ref`` (pu systeme, ENTREE ABSOLUE, pas un echelon) sous un facteur de
    puissance donne par ``q_ref``. Hors ``p_ref``/``q_ref`` c'est exactement
    un :class:`PLoad` : courant filtre au premier ordre (``Tf``), puissance
    constante jusqu'a ``Vbt`` (0,7 pu), impedance constante en dessous.

    ``p_ref`` se pilote par une courbe de charge de plusieurs minutes
    (:class:`g2egrid.temporel.Courbe`).
    """

    svc_inputs = ("p_ref", "q_ref")

    @staticmethod
    def _puissances(s):
        return s.p_ref, s.q_ref

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = super().initialize(ctx)
        init.pop("dp"), init.pop("dq")
        init["p_ref"], init["q_ref"] = self.P, self.Q
        return init
