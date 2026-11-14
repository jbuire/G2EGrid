"""Convertisseurs grid-following : portage FIDELE de ``matlab/G2Elin_V2/core/GFL1..2.m``.

====== ================================================= ======
Classe Description                                       Etats
====== ================================================= ======
GFL1   boucles externes P/Q, bus DC ideal                  12
GFL2   boucle externe vdc/Q, bus DC dynamique              14
====== ================================================= ======

Structure commune (identique a MATLAB) :

* **filtre LCL** (6 etats) : ``isd/isq``, ``igd/igq``, ``ved/veq`` ;
* **PLL** : ``dM_pll = Kipll.veq``, ``w_pll = M_pll + Kppll.veq + wff``,
  ``dtheta_pll = wb.w_pll`` -- la PLL se verrouille sur la tension du
  CONDENSATEUR (``veq = 0`` a l'equilibre) ;
* **boucles externes** : PI sur P (GFL1) ou sur ``vdc`` (GFL2), et PI sur Q ;
* **boucle de courant interne** : ``M_CLd/q``, decouplage ``wff``,
  feedforward ``Kffv`` ;
* **bus DC** (GFL2) : ``vdc``, ``idc``, gain ``Kpdc`` recale au point de
  fonctionnement (correction CPL, comme MATLAB).

Comme cote MATLAB, les gains sont auto-regles depuis les temps de reponse
(``tau_CL``, ``tau_PL``, ``tau_QL``, ``tau_PLL``, ``tau_VL``).

Repere : l'axe d de la PLL est en avance de ``theta_pll`` sur l'axe D global.
``dtheta_pll = wb.w_pll`` (MATLAB, angle absolu) devient
``wb.(w_pll - w_ref)`` en repere global tournant : meme matrice A, seule la
colonne ``w_ref`` de B change.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from g2egrid.noyau.element import Element, InitContext
from g2egrid.elements.limiteurs import (
    LimiteurCourantTension,
    limiter_courant,
    limiter_tension,
)

__all__ = ["GFL1", "GFL2", "GFL2DC", "GFL_DEFAULT", "GFL2_DEFAULT"]

# Valeurs par defaut : reprise LITTERALE des proprietes de GFL1.m.
GFL_DEFAULT = dict(
    R1_pu=0.015,
    L1_pu=0.1,
    C_pu=0.11,
    R2_pu=0.003,
    L2_pu=0.01,
    Kffv=1.0,
    wff=1.0,
    tau_CL=0.01,
    zeta_CL=0.7,
    tau_PL=0.1,
    zeta_PL=0.7,
    tau_QL=0.1,
    zeta_QL=0.7,
    tau_PLL=0.05,
    zeta_PLL=0.7,
    # limitation de courant (cf. g2egrid.elements.limiteurs)
    Imax=1.2,
    Taw=2e-3,
    mode_lim=0.0,
    # saturation de tension : |m| <= Mmax sur la modulation de la boucle de
    # courant, retro-calcul sur ses integrateurs -- comme le GFM. Mmax = 0 :
    # pas de saturation (defaut, y compris sur bus continu dynamique, voir
    # limiteurs.REGLAGES_BUS_CONTINU) ; MMAX = 1,15 pour une modulation
    # vectorielle. Taw_v = None : tau_CL/2, fixe au constructeur (voir
    # GFM_DEFAULT pour ce que fait un retro-calcul lent).
    Mmax=0.0,
    Taw_v=None,
    mode_vsat=0.0,
)

# GFL2 : cote continu + boucle de tension DC (GFL2.m).
GFL2_DEFAULT = dict(Cdc_pu=0.2, Gdc_pu=17e-6, Tdc=0.001, tau_VL=0.1, zeta_VL=0.7)


class _GFLBase(LimiteurCourantTension, Element):
    """Tronc commun des GFL : filtre LCL + PLL + boucle de courant,
    limitation de courant et saturation de tension."""

    ports = ("t",)
    #: variables algebriques, dans l'ordre de ``algVec`` de symGFL
    algebs = ("md", "mq", "w_pll")
    _base_params = (
        "R1_pu",
        "L1_pu",
        "C_pu",
        "R2_pu",
        "L2_pu",
        "KpCL",
        "KiCL",
        "Kffv",
        "wff",
        "Kppll",
        "Kipll",
        "kS",
        "wb",
        "Imax",
        "Taw",
        "mode_lim",
        "Mmax",
        "Taw_v",
        "mode_vsat",
    )

    def _tune(self, p: dict, wb: float, given: dict | None = None) -> None:
        """Auto-reglage des PI — ``computeControllerGains()`` de MATLAB.

        Un gain passe explicitement au constructeur n'est jamais ecrase.
        """
        given = given or {}
        K0_CL = 1.0 / p["R1_pu"]  # vdc0 = 1
        tau0_CL = p["L1_pu"] / (wb * p["R1_pu"])
        auto = dict(
            KpCL=(6 * tau0_CL - p["tau_CL"]) / (K0_CL * p["tau_CL"]),
            KiCL=(3 / (p["zeta_CL"] * p["tau_CL"])) ** 2 * tau0_CL / K0_CL,
        )
        wn_pll = 3 / (p["zeta_PLL"] * p["tau_PLL"])  # V0 = 1
        auto["Kipll"] = wn_pll**2 / wb
        auto["Kppll"] = 2 * p["zeta_PLL"] * wn_pll / wb
        for k, v in auto.items():
            if k not in given:
                p[k] = v
        if p.get("Taw_v") is None:
            p["Taw_v"] = p["tau_CL"] / 2  # cf. GFL_DEFAULT

    @staticmethod
    def _common(s, vdc, isd_ref, isq_ref):
        """Rotation, LCL, PLL et modulation — communs aux deux GFL.

        Rend aussi ``(awv_d, awv_q)``, le retro-calcul de la saturation de
        tension, a ajouter (fois ``vdc``) aux integrateurs de la boucle de
        courant. Les ``md``, ``mq`` rendus sont SATURES : c'est ce que
        l'onduleur applique, au filtre comme au bus continu ; les algebriques
        ``md``, ``mq`` restent la modulation LIBRE, que lit la garde."""
        # Angle repere par rapport au repere GLOBAL, comme symGFL (theta_g
        # est l'angle de l'element de reference d'angle du reseau).
        ct, st = sp.cos(s.theta_pll - s.theta_g), sp.sin(s.theta_pll - s.theta_g)
        vgd = s.vd_t * ct + s.vq_t * st  # global -> repere PLL
        vgq = -s.vd_t * st + s.vq_t * ct

        # md, mq et w_pll sont ALGEBRIQUES, comme dans symGFL (algVec = [md mq
        # w_pll]). Ce n'est pas cosmetique : la reference substitue dans son
        # jacobien les valeurs d'equilibre md_0, mq_0 et w_pll_0 issues de ses
        # propres formules, et non ce que rendraient alg1..alg3 au meme point.
        # Son initialisation des integrateurs etant incoherente (cf. le mode
        # init='symGFL'), les deux jeux different, et cela se lit sur les
        # coefficients en vdc et sur les termes de couplage en w_pll.
        w_pll, md, mq = s.w_pll, s.md, s.mq
        alg = [
            md
            - (1 / vdc)
            * (
                s.KpCL * (isd_ref - s.isd)
                + s.Kffv * s.ved
                - s.wff * s.L1_pu * s.isq
                + s.M_CLd
            ),
            mq
            - (1 / vdc)
            * (
                s.KpCL * (isq_ref - s.isq)
                + s.Kffv * s.veq
                + s.wff * s.L1_pu * s.isd
                + s.M_CLq
            ),
            w_pll - s.M_pll - s.Kppll * s.veq - s.wff,
        ]
        # saturation de tension (cf. limiter_tension) : a mode_vsat = 0, la
        # modulation passe telle quelle et le retro-calcul est nul
        md, mq, awv_d, awv_q = limiter_tension(s, md, mq)

        lcl = [
            (s.wb / s.L1_pu)
            * (md * vdc - s.ved - s.R1_pu * s.isd + w_pll * s.L1_pu * s.isq),
            (s.wb / s.L1_pu)
            * (mq * vdc - s.veq - s.R1_pu * s.isq - w_pll * s.L1_pu * s.isd),
            (s.wb / s.L2_pu)
            * (s.ved - vgd - s.R2_pu * s.igd + w_pll * s.L2_pu * s.igq),
            (s.wb / s.L2_pu)
            * (s.veq - vgq - s.R2_pu * s.igq - w_pll * s.L2_pu * s.igd),
            (s.wb / s.C_pu) * (s.isd - s.igd + w_pll * s.C_pu * s.veq),
            (s.wb / s.C_pu) * (s.isq - s.igq - w_pll * s.C_pu * s.ved),
        ]
        iD = s.kS * (s.igd * ct - s.igq * st)
        iQ = s.kS * (s.igd * st + s.igq * ct)
        Vt = sp.sqrt(s.ved**2 + s.veq**2)
        return lcl, alg, (w_pll, md, mq, iD, iQ, Vt), (awv_d, awv_q)

    @staticmethod
    def _retro_calcul(Ki, aw):
        """Le retro-calcul ``aw`` de la limitation de courant sur un
        integrateur de gain ``Ki`` -- nul si ``Ki = 0``.

        La boucle de tension continue des GFL2 est PROPORTIONNELLE par defaut
        (``Kidc = 0``, comme GFL2.m) : son etat ``M_d`` n'integre rien, c'est
        le courant du point de fonctionnement, un biais. Il n'y a pas
        d'emballement a combattre, et le retro-calcul ne ferait que deplacer
        ce biais pendant la limitation : la tension continue garderait l'ecart
        ensuite, aucune integrale ne la ramenant (``cigre_mv_dc``, defaut au
        bus 15 : 1,041 pu au lieu de 1 apres 5 ms de limitation)."""
        return sp.Piecewise((0, sp.Eq(Ki, 0)), (aw, True))

    def _init_common(self, ctx: InitContext, vdc: float = 1.0) -> dict[str, float]:
        """Equilibre analytique EXACT (la PLL s'aligne sur ``ve``).

        ``vdc`` est la tension continue AU POINT DE DEPART : 1 pour les
        variantes qui portent leur propre condensateur, la tension du bus
        continu pour celles a port DC (cf. :class:`GFL2DC`).
        """
        p = self.pdict
        w0 = ctx.w
        V = ctx.v["t"]
        S = ctx.s.get("t", 0j)
        Ig = (np.conj(S / V) if abs(V) > 1e-12 else 0j) / p["kS"]

        Ve = V + complex(p["R2_pu"], w0 * p["L2_pu"]) * Ig
        # veq = 0 a l'equilibre : la PLL s'aligne sur ve. Son angle est
        # compte depuis le REPERE GLOBAL, qui tourne de theta_g (l'angle rotor
        # de la machine de reference) : oublier theta_g laissait au Newton le
        # soin de le rattraper, ce qu'il fait tant que |theta_g| < pi/2 -- au-
        # dela, il pouvait tomber sur la branche ou la PLL est calee a 180
        # degres (ved < 0), un equilibre instable.
        angle_ve = float(np.angle(Ve))
        theta_pll = angle_ve + ctx.theta_g
        ved, veq = abs(Ve), 0.0
        ig = Ig * np.exp(-1j * angle_ve)
        igd, igq = ig.real, ig.imag
        isd = igd - w0 * p["C_pu"] * veq
        isq = igq + w0 * p["C_pu"] * ved

        md = (ved + p["R1_pu"] * isd - w0 * p["L1_pu"] * isq) / vdc
        mq = (veq + p["R1_pu"] * isq + w0 * p["L1_pu"] * isd) / vdc

        p0 = ved * igd + veq * igq
        q0 = -ved * igq + veq * igd
        return {
            "isd": isd,
            "isq": isq,
            "igd": igd,
            "igq": igq,
            "ved": ved,
            "veq": veq,
            "M_d": isd,
            "M_q": isq,
            "M_CLd": md * vdc - p["Kffv"] * ved + p["wff"] * p["L1_pu"] * isq,
            "M_CLq": mq * vdc - p["Kffv"] * veq - p["wff"] * p["L1_pu"] * isd,
            "M_pll": w0 - p["wff"],
            "theta_pll": theta_pll,
            "md": md,
            "mq": mq,
            "w_pll": w0,
            "_p0": p0,
            "_q0": q0,
            "_md": md,
            "_mq": mq,
            "_vdc": vdc,
        }


class GFL1(_GFLBase):
    """GFL a boucles P/Q, bus DC ideal (12 etats) — cf ``GFL1.m``."""

    params = _GFLBase._base_params + ("Kpp", "Kip", "Kpq", "Kiq")
    states = (
        "isd",
        "isq",
        "igd",
        "igq",
        "ved",
        "veq",
        "M_d",
        "M_q",
        "M_CLd",
        "M_CLq",
        "M_pll",
        "theta_pll",
    )
    svc_inputs = ("p_ref", "q_ref")
    # La SORTIE de service se nomme w_pll_out : symGFL declare aussi un
    # ALGEBRIQUE w_pll, et le compilateur refuse le meme nom des deux cotes.
    svc_outputs = ("P", "Q", "w_pll_out", "Vt", "theta_pll_out")

    def __init__(self, name: str, bus: int, Sn: float, Sb: float, wb: float, **params):
        p = dict(GFL_DEFAULT)
        p.update(params)
        self._tune(p, wb, params)
        # GFL1.m : reglage du 1er ordre sur les boucles externes, sauf si
        # l'utilisateur impose directement les gains.
        p.setdefault("Kip", 3 / p["tau_PL"])
        p.setdefault("Kpp", 0.0)
        p.setdefault("Kiq", 3 / p["tau_QL"])
        p.setdefault("Kpq", 0.0)
        p["kS"] = Sn / Sb
        p["wb"] = wb
        super().__init__(name, bus, **{k: p[k] for k in self.params})
        self.Sn, self.Sb = float(Sn), float(Sb)

    @staticmethod
    def equations(s):
        pe = s.ved * s.igd + s.veq * s.igq
        qe = -s.ved * s.igq + s.veq * s.igd
        isd_ref = s.Kpp * (s.p_ref - pe) + s.M_d
        isq_ref = s.Kpq * (qe - s.q_ref) + s.M_q
        isd_ref, isq_ref, aw_d, aw_q = limiter_courant(s, isd_ref, isq_ref)
        lcl, alg, (w_pll, md, mq, iD, iQ, Vt), (awv_d, awv_q) = _GFLBase._common(
            s, 1, isd_ref, isq_ref
        )
        f = lcl + [
            s.Kip * (s.p_ref - pe) + aw_d,
            s.Kiq * (qe - s.q_ref) + aw_q,
            s.KiCL * (isd_ref - s.isd) + awv_d,
            s.KiCL * (isq_ref - s.isq) + awv_q,
            s.Kipll * s.veq,
            s.wb * (w_pll - s.w_ref),
        ]
        h = [iD, iQ, pe, qe, w_pll, Vt, s.theta_pll]
        return f, alg, h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = self._init_common(ctx)
        init["p_ref"] = init.pop("_p0")
        init["q_ref"] = init.pop("_q0")
        return init


class GFL2(_GFLBase):
    """GFL a boucle vdc/Q, bus DC dynamique (14 etats) — cf ``GFL2.m``."""

    params = _GFLBase._base_params + (
        "Kpdc",
        "Kidc",
        "Kpq",
        "Kiq",
        "Cdc_pu",
        "Gdc_pu",
        "Tdc",
    )
    states = (
        "isd",
        "isq",
        "igd",
        "igq",
        "ved",
        "veq",
        "vdc",
        "idc",
        "M_d",
        "M_q",
        "M_CLd",
        "M_CLq",
        "M_pll",
        "theta_pll",
    )
    svc_inputs = ("vdc_ref", "q_ref", "idc_ref")
    svc_outputs = ("P", "Q", "w_pll_out", "Vt", "vdc_out", "theta_pll_out")

    #: conventions d'initialisation acceptees par ``init``
    INITS = ("physique", "symGFL")

    def __init__(
        self,
        name: str,
        bus: int,
        Sn: float,
        Sb: float,
        wb: float,
        init: str = "physique",
        **params,
    ):
        """``init`` choisit la convention de calcul du point d'equilibre.

        ``'physique'`` (defaut) pose chaque integrateur a la valeur qui annule
        sa propre derivee : ``M_d = isd``, ``M_q = isq``, ``M_pll`` tel que
        ``w_pll = w``.

        ``'symGFL'`` rejoue celle de G2Elin V1, qui DIVISE ces trois valeurs par
        le gain de l'integrateur ::

            M_d_0   = Is0_d / Kid        (Kid = -77783 ici, donc M_d_0 ~ 0)
            M_q_0   = Is0_q / Kiq
            M_pll_0 = (-1 - veq0.Kppll) / Kipll

        C'est une initialisation de schema-bloc Simulink, ou le gain est
        applique APRES l'integrateur, transposee telle quelle au modele d'etat
        ou il est applique avant. L'equilibre obtenu ne verifie donc pas les
        equations -- ``isd_ref`` y vaut 0 au lieu de ``isd`` -- mais c'est ce
        que la reference linearise, et cela se lit sur deux coefficients de A
        (``isd <- vdc`` et ``isq <- vdc``). Dans ce mode,
        ``raffiner_equilibre`` passe a False.
        """
        if init not in self.INITS:
            raise ValueError(
                f"init={init!r} inconnu (attendu : {', '.join(self.INITS)})"
            )
        self.init_mode = init
        self.raffiner_equilibre = init != "symGFL"
        p = dict(GFL_DEFAULT)
        p.update(GFL2_DEFAULT)
        p.update(params)
        self._tune(p, wb, params)
        p.setdefault("Kiq", 100.0)  # GFL2.m : pas d'auto-reglage
        p.setdefault("Kpq", 0.0)
        # Kpdc/Kidc : recales au point de fonctionnement, sauf si imposes.
        self._auto_kpdc = "Kpdc" not in params
        p.setdefault("Kpdc", 1.0)
        p.setdefault("Kidc", 0.0)
        p["kS"] = Sn / Sb
        p["wb"] = wb
        self.tuning = {k: p[k] for k in ("tau_VL", "zeta_VL")}
        super().__init__(name, bus, **{k: p[k] for k in self.params})
        self.Sn, self.Sb = float(Sn), float(Sb)

    @staticmethod
    def equations(s):
        pe = s.ved * s.igd + s.veq * s.igq
        qe = -s.ved * s.igq + s.veq * s.igd
        isd_ref = s.Kpdc * (s.vdc - s.vdc_ref) + s.M_d
        isq_ref = s.Kpq * (qe - s.q_ref) + s.M_q
        isd_ref, isq_ref, aw_d, aw_q = limiter_courant(s, isd_ref, isq_ref)
        lcl, alg, (w_pll, md, mq, iD, iQ, Vt), (awv_d, awv_q) = _GFLBase._common(
            s, s.vdc, isd_ref, isq_ref
        )
        f = lcl + [
            (s.wb / s.Cdc_pu) * (s.idc - s.Gdc_pu * s.vdc - md * s.isd - mq * s.isq),
            (s.idc_ref - s.idc) / s.Tdc,
            s.Kidc * (s.vdc - s.vdc_ref) + _GFLBase._retro_calcul(s.Kidc, aw_d),
            s.Kiq * (qe - s.q_ref) + aw_q,
            s.KiCL * (isd_ref - s.isd) + awv_d * s.vdc,
            s.KiCL * (isq_ref - s.isq) + awv_q * s.vdc,
            s.Kipll * s.veq,
            s.wb * (w_pll - s.w_ref),
        ]
        h = [iD, iQ, pe, qe, w_pll, Vt, s.vdc, s.theta_pll]
        return f, alg, h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = self._init_common(ctx)
        p = self.pdict
        md, mq, vdc = init.pop("_md"), init.pop("_mq"), init.pop("_vdc")
        init.pop("_p0")
        q0 = init.pop("_q0")
        idc = p["Gdc_pu"] * vdc + md * init["isd"] + mq * init["isq"]
        init.update(vdc=vdc, idc=idc, vdc_ref=vdc, q_ref=q0, idc_ref=idc)
        if self.init_mode == "symGFL":
            # Les gains de V1 sont ceux de la convention (ref - mesure) ; les
            # notres sont ceux de (mesure - ref), d'ou les signes opposes.
            init["M_d"] = init["isd"] / -p["Kidc"] if p["Kidc"] else 0.0
            init["M_q"] = init["isq"] / -p["Kiq"] if p["Kiq"] else 0.0
            init["M_pll"] = (-1.0 - init["veq"] * p["Kppll"]) / p["Kipll"]
        # Reglage CPL-aware de Kpdc au point de fonctionnement (GFL2.m)
        if self._auto_kpdc:
            P0_dc = md * init["isd"] + mq * init["isq"]
            Geff = p["Gdc_pu"] - P0_dc / vdc**2
            self.set_param(
                "Kpdc",
                (3 * p["Cdc_pu"] / (self.tuning["tau_VL"] * p["wb"]) - Geff) / md,
            )
        return init


class GFL2DC(_GFLBase):
    """GFL a boucle vdc/Q, PORT continu (12 etats) — cf :class:`GFL2`.

    Meme decoupage que :class:`~g2egrid.elements.gfm._GFMDCBase` : les deux
    etats continus ``vdc`` et ``idc`` sortent du convertisseur. ``vdc`` devient
    l'entree de port ``v_dc``, et le condensateur devient la capacite du
    :class:`~g2egrid.elements.bus.BusDC`.

    DIFFERENCE AVEC LE GFM : ``vdc_ref`` RESTE une entree de service. Un
    grid-following a port continu tient lui-meme sa tension DC -- c'est sa
    boucle externe, ``isd_ref = Kpdc.(v_dc - vdc_ref) + M_d`` -- la ou un
    grid-forming laisse cette tension a la charge du reseau continu. Ce qui
    disparait, c'est seulement ``idc_ref`` : le courant continu n'est plus un
    etat du convertisseur mais ce que le noeud lui apporte.

    ``Cdc_pu`` RESTE un parametre, bien que le condensateur soit maintenant
    dans le bus : il ne sert plus a une equation d'etat mais au REGLAGE de
    ``Kpdc``, cale au point de fonctionnement comme dans ``GFL2.m``.
    L'assembleur le renseigne depuis la capacite du noeud auquel le
    convertisseur est raccorde.
    """

    ports_dc = ("dc",)
    params = _GFLBase._base_params + (
        "Kpdc",
        "Kidc",
        "Kpq",
        "Kiq",
        "Cdc_pu",
        "Gdc_pu",
        "kVdc",
        "bloque",
    )
    states = (
        "isd",
        "isq",
        "igd",
        "igq",
        "ved",
        "veq",
        "M_d",
        "M_q",
        "M_CLd",
        "M_CLq",
        "M_pll",
        "theta_pll",
    )
    svc_inputs = ("vdc_ref", "q_ref")
    svc_outputs = ("P", "Q", "w_pll_out", "Vt", "vdc_out", "theta_pll_out")

    #: Un grid-following ne tient pas la tension de son bus alternatif : son
    #: bus est PQ, et le reactif est programme (``q_ref``). Du cote continu, sa
    #: boucle externe porte sur ``vdc`` -- a l'equilibre ``v_dc = vdc_ref`` --
    #: donc il SAIT tenir la tension continue : c'est ``VdcQ``, et sa puissance
    #: alternative devient l'inconnue. En ``PQ``, c'est l'inverse : les deux
    #: puissances alternatives sont programmees, la tension continue sort du
    #: bilan de l'ile, et ``vdc_ref`` se cale dessus a l'initialisation.
    #:
    #: ``droop`` : statisme tension continue / puissance, ``kdroop`` pu de
    #: tension par pu de puissance (base Sn). Sa boucle externe est deja
    #: proportionnelle (``Kidc = 0``) : c'est un statisme, dont la pente est
    #: ``Kpdc``. En ``droop``, ``Kpdc`` n'est donc plus recale sur le
    #: condensateur mais sur le statisme, ``Kpdc = 1/(kdroop.ved0)``, de sorte
    #: que le modele dynamique et le load flow suivent la meme pente.
    modes_dc = ("VdcQ", "PQ", "droop", "droopI")
    #: statisme (``droop``), pose par l'assembleur ; ``None`` sinon
    kdroop: float | None = None
    #: statisme COURANT (``droopI``) : consignes ``vdc_ref`` et ``p_ref`` (pu
    #: du convertisseur) de la loi du load flow, posees par l'assembleur
    droop_courant: dict | None = None

    #: conventions d'initialisation acceptees par ``init``
    INITS = ("physique", "symGFL")

    def __init__(
        self,
        name: str,
        bus: int,
        bus_dc: int,
        Sn: float,
        Sb: float,
        wb: float,
        init: str = "physique",
        **params,
    ):
        """``bus`` est le noeud alternatif, ``bus_dc`` le noeud continu.

        Les deux arrivent separes plutot que dans une liste : l'ordre des bus
        d'un element a port DC (alternatifs d'abord, continus ensuite) est une
        convention du noyau, pas quelque chose que l'appelant doit retenir.
        """
        if init not in self.INITS:
            raise ValueError(
                f"init={init!r} inconnu (attendu : {', '.join(self.INITS)})"
            )
        self.init_mode = init
        self.raffiner_equilibre = init != "symGFL"
        p = dict(GFL_DEFAULT)
        p.update(GFL2_DEFAULT)
        p.update(params)
        self._tune(p, wb, params)
        p.setdefault("Kiq", 100.0)  # GFL2.m : pas d'auto-reglage
        p.setdefault("Kpq", 0.0)
        # Kpdc/Kidc : recales au point de fonctionnement, sauf si imposes.
        self._auto_kpdc = "Kpdc" not in params
        p.setdefault("Kpdc", 1.0)
        p.setdefault("Kidc", 0.0)
        p.setdefault("kVdc", 1.0)
        p.setdefault("bloque", 0.0)
        p["kS"] = Sn / Sb
        p["wb"] = wb
        self.tuning = {k: p[k] for k in ("tau_VL", "zeta_VL")}
        Element.__init__(self, name, [bus, bus_dc], **{k: p[k] for k in self.params})
        self.Sn, self.Sb = float(Sn), float(Sb)

    @staticmethod
    def equations(s):
        pe = s.ved * s.igd + s.veq * s.igq
        qe = -s.ved * s.igq + s.veq * s.igd
        isd_ref = s.Kpdc * (s.v_dc - s.vdc_ref) + s.M_d
        isq_ref = s.Kpq * (qe - s.q_ref) + s.M_q
        isd_ref, isq_ref, aw_d, aw_q = limiter_courant(s, isd_ref, isq_ref)
        # bloque : la modulation ne lit plus le noeud continu (cf. GFM3DC)
        vdc_m = s.v_dc * (1 - s.bloque) + s.bloque
        lcl, alg, (w_pll, md, mq, iD, iQ, Vt), (awv_d, awv_q) = _GFLBase._common(
            s, vdc_m, isd_ref, isq_ref
        )
        f = lcl + [
            s.Kidc * (s.v_dc - s.vdc_ref) + _GFLBase._retro_calcul(s.Kidc, aw_d),
            s.Kiq * (qe - s.q_ref) + aw_q,
            s.KiCL * (isd_ref - s.isd) + awv_d * vdc_m,
            s.KiCL * (isq_ref - s.isq) + awv_q * vdc_m,
            s.Kipll * s.veq,
            s.wb * (w_pll - s.w_ref),
        ]
        i_dc = (
            -s.kS
            * s.kVdc
            * (s.Gdc_pu * s.v_dc + md * s.isd + mq * s.isq)
            * (1 - s.bloque)
        )
        h = [iD, iQ, i_dc, pe, qe, w_pll, Vt, s.v_dc, s.theta_pll]
        return f, alg, h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        vdc = float(np.real(ctx.v["dc"]))
        init = self._init_common(ctx, vdc=vdc)
        p = self.pdict
        md, mq = init.pop("_md"), init.pop("_mq")
        init.pop("_vdc")
        init.pop("_p0")
        q0 = init.pop("_q0")
        init.update(vdc_ref=vdc, q_ref=q0)
        if self.init_mode == "symGFL":
            # Les gains de V1 sont ceux de la convention (ref - mesure) ; les
            # notres sont ceux de (mesure - ref), d'ou les signes opposes.
            init["M_d"] = init["isd"] / -p["Kidc"] if p["Kidc"] else 0.0
            init["M_q"] = init["isq"] / -p["Kiq"] if p["Kiq"] else 0.0
            init["M_pll"] = (-1.0 - init["veq"] * p["Kppll"]) / p["Kipll"]
        if self.kdroop and self.droop_courant:
            # statisme courant (droopI) : p = v_dc.(i_ref + (v_dc - v_ref)/k) ;
            # Kpdc donne la meme pente dp/dv_dc au point de fonctionnement
            d = self.droop_courant
            i_ref = d["p_ref"] / d["vdc_ref"]
            pente = i_ref + (2 * vdc - d["vdc_ref"]) / self.kdroop
            self.set_param("Kpdc", pente / init["ved"])
        elif self.kdroop:
            # statisme : dp = ved.disd = dv_dc/kdroop
            self.set_param("Kpdc", 1.0 / (self.kdroop * init["ved"]))
        # Reglage CPL-aware de Kpdc au point de fonctionnement (GFL2.m)
        elif self._auto_kpdc:
            P0_dc = md * init["isd"] + mq * init["isq"]
            Geff = p["Gdc_pu"] - P0_dc / vdc**2
            self.set_param(
                "Kpdc",
                (3 * p["Cdc_pu"] / (self.tuning["tau_VL"] * p["wb"]) - Geff) / md,
            )
        return init
