"""Convertisseur modulaire multiniveaux : deux MODELES MOYENS a energie.

Le modele des stations des reseaux continus a haute tension. Ce n'est ni un
modele de sous-modules ni un modele de bras : c'est le modele MOYEN dont se
servent les etudes de stabilite petit signal des reseaux MTDC -- le
convertisseur y est une source de tension commandee de chaque cote, et ce
qui le distingue d'un VSC a deux niveaux, c'est que SON BUS CONTINU N'EXISTE
PAS. Il n'y a pas de condensateur entre les deux cotes : il y a l'ENERGIE
``W`` stockee dans les condensateurs des sous-modules, que les deux cotes
alimentent.

====================  ==================================================
cote alternatif        ``isd/isq`` a travers ``L1 = Lc + Larm/2``, ``R1``
cote continu           ``idc`` a travers ``Ldc = 2.Larm/3``, ``Rdc``
energie                ``dW/dt = v_z.idc - (e_d.isd + e_q.isq)``
commande               PI de courant alternatif (``M_CLd/q``), boucles
                       externes vdc/Q (``M_d/M_q``) sur la tension continue
                       mesuree et filtree (``vdc_f``, ``T_vdc``), et PLL
                       (``M_pll``, ``theta_pll``)
====================  ==================================================

La tension alternative que les bras fabriquent est IMPOSEE par la commande
(modulation compensee, ``e = e_ref``, tant que l'energie des bras le permet :
la saturation de la modulation n'est pas representee). Ce qui distingue les
deux modeles, c'est ce qui fixe ``v_z``, la tension continue des bras :

* :class:`MMC` (11 etats) -- ENERGIE LIBRE, le modele moyen de reference du
  CIGRE (brochures B4.57/B4.58) : ``v_z`` est la somme des tensions des
  condensateurs, ``v_z = v_c = sqrt(W/W0)``. Vu du continu, la station est
  une inductance ``Ldc`` en serie avec la capacite EQUIVALENTE de ses
  sous-modules, ``C_eq = 2.W0`` (en secondes) : c'est elle qui tient la
  tension du noeud entre deux actions de la commande.
* :class:`MMCEnergie` (13 etats) -- ENERGIE REGULEE : ``v_z`` est fabrique
  par une boucle de courant continu (``M_z``), dont la consigne vient d'une
  boucle d'energie (``M_W``),
  ``idc_ref = p_t + KpW.(W0 - W) + M_W`` (anticipation de la puissance
  MESUREE aux bornes ; l'integrateur reprend les pertes), et
  ``v_z = v_dc - [KpZ.(idc_ref - idc) + M_z]``. L'energie est decouplee de
  la tension continue : la station ne montre plus sa capacite au reseau
  continu, dont les noeuds ne gardent que celle des cables. Un statisme
  raide sur un noeud de faible capacite devient alors une boucle tres rapide,
  que la moindre dynamique interne destabilise -- le resultat connu des
  reseaux MTDC a MMC regules en energie, que le cas CIGRE B4 montre
  (``tests/test_mmc.py``).

LE PORT CONTINU est celui de :class:`~g2egrid.elements.gfl.GFL2DC` : entree
``v_dc``, sortie ``i_dc = -kS.kVdc.idc`` (le courant que la station tire de
son noeud). Les modes sont les siens (``VdcQ``, ``PQ``, ``droop``) et sa
boucle externe aussi, ``isd_ref = Kpdc.(v_dc - vdc_ref) + M_d``, reglee de
meme sur la capacite vue du noeud (celle du noeud, plus ``C_eq`` en energie
libre) ou sur le statisme. Un MMC n'a pas de condensateur continu : la
capacite qu'il apporte a son noeud (``Cdc_pu``, faible) est celle de ses
equipements de pole.

Le BLOCAGE (``bloque = 1``, manoeuvre ``"blocage"``) coupe la tension qui
entraine le courant continu : les bras bloques l'eteignent en
:data:`T_BLOCAGE`, et le courant, celui d'une inductance, reste continu au
noeud. Les diodes de roue libre (qui empecheraient la tension continue de
s'inverser) ne sont pas representees.

LA MESURE DE TENSION CONTINUE EST FILTREE (``T_vdc``, 5 ms). Sans filtre,
le proportionnel de la boucle de tension agit encore a la frequence de
resonance du noeud continu (sa capacite entre le cable et ``Ldc``, 150 a
250 Hz) ; a travers l'integrateur d'energie, il y donne a la station une
conductance NEGATIVE, et le cas CIGRE DCS1 part instable vers 200 Hz. Le
filtre tourne cette action en reactance. De meme, l'anticipation de tension
alternative est nulle par defaut (``Kffv = 0``) : a 1, le systeme CIGRE
complet est instable vers 195 Hz.

Les gains sont auto-regles depuis les temps de reponse, comme ceux des GFL :
``tau_CL`` (courant alternatif), ``tau_Z`` (courant continu), ``tau_W``
(energie), ``tau_VL`` (tension continue), ``tau_PLL``.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from g2egrid.noyau.element import Element, InitContext
from g2egrid.elements.limiteurs import LimiteurCourant, limiter_courant

__all__ = ["MMC", "MMCEnergie", "MMC_DEFAULT"]

#: Reglages par defaut, en pu de la station (base Sn, tension continue pole a
#: pole). ``W0`` est l'energie stockee nominale en secondes de puissance
#: assignee (30 kJ/MVA, l'ordre de grandeur des stations a demi-pont).
MMC_DEFAULT = dict(
    L1_pu=0.15,
    R1_pu=0.003,
    Ldc_pu=0.1,
    Rdc_pu=0.002,
    W0=0.03,
    Kffv=0.0,
    wff=1.0,
    tau_CL=0.002,
    zeta_CL=0.7,
    tau_Z=0.005,
    zeta_Z=0.7,
    tau_W=0.1,
    zeta_W=0.7,
    tau_VL=0.02,
    # filtre de la mesure de tension continue [s]
    T_vdc=0.005,
    tau_PLL=0.05,
    zeta_PLL=0.7,
    # capacite apportee au noeud continu (equipements de pole), pu station
    Cdc_pu=0.02,
    Imax=1.2,
    Taw=2e-3,
    mode_lim=0.0,
)


#: Constante de temps d'extinction du courant continu d'une station BLOQUEE
#: [s] : les bras bloques inserent leurs condensateurs dans le chemin du
#: courant, qui s'eteint en une milliseconde environ.
T_BLOCAGE = 1e-3


def _pi_second_ordre(L: float, R: float, wb: float, tau: float, zeta: float):
    """PI d'un courant dans ``(L/wb).di/dt = u - R.i`` : poles a
    ``wn = 3/(zeta.tau)``, amortissement ``zeta`` (la regle des GFL)."""
    wn = 3.0 / (zeta * tau)
    return 2 * zeta * wn * L / wb - R, wn**2 * L / wb


class MMC(LimiteurCourant, Element):
    """MMC, modele moyen a ENERGIE LIBRE, port continu (11 etats).

    Se declare comme une station de ``T_GFL2s`` raccordee par ``T_VSCs``, avec
    le reglage ``modele = "mmc"`` (famille ``GFL_PARAMS`` ou
    ``PARAMS_PAR_BUS``) ; ``"mmc_energie"`` donne :class:`MMCEnergie`.
    """

    #: l'energie est-elle regulee (boucles M_z, M_W) ?
    regule = False
    #: l'initialisation analytique est EXACTE (voir :meth:`initialize`) : un
    #: Newton par element n'a rien a raffiner, et le Powell de scipy, qui
    #: explore loin de son depart, y essayait des energies negatives
    raffiner_equilibre = False
    ports = ("t",)
    ports_dc = ("dc",)
    params = (
        "L1_pu",
        "R1_pu",
        "Ldc_pu",
        "Rdc_pu",
        "W0",
        "KpCL",
        "KiCL",
        "Kffv",
        "wff",
        "Kppll",
        "Kipll",
        "Kpdc",
        "Kidc",
        "Kpq",
        "Kiq",
        "Cdc_pu",
        "T_vdc",
        "kS",
        "kVdc",
        "wb",
        "Imax",
        "Taw",
        "mode_lim",
        "bloque",
    )
    states = (
        "isd",
        "isq",
        "idc",
        "W",
        "M_d",
        "M_q",
        "M_CLd",
        "M_CLq",
        "vdc_f",
        "M_pll",
        "theta_pll",
    )
    svc_inputs = ("vdc_ref", "q_ref")
    svc_outputs = ("P", "Q", "w_pll_out", "Vt", "vdc_out", "theta_pll_out", "W_out")

    #: les modes de :class:`~g2egrid.elements.gfl.GFL2DC`, pour les memes raisons
    modes_dc = ("VdcQ", "PQ", "droop", "droopI")
    #: statisme (``droop``), pose par l'assembleur ; ``None`` sinon
    kdroop: float | None = None
    #: statisme COURANT (``droopI``), cf. :class:`~g2egrid.elements.gfl.GFL2DC`
    droop_courant: dict | None = None

    def __init__(
        self,
        name: str,
        bus: int,
        bus_dc: int,
        Sn: float,
        Sb: float,
        wb: float,
        **params,
    ):
        p = dict(MMC_DEFAULT)
        p.update(params)
        auto = {}
        auto["KpCL"], auto["KiCL"] = _pi_second_ordre(
            p["L1_pu"], p["R1_pu"], wb, p["tau_CL"], p["zeta_CL"]
        )
        auto["KpZ"], auto["KiZ"] = _pi_second_ordre(
            p["Ldc_pu"], p["Rdc_pu"], wb, p["tau_Z"], p["zeta_Z"]
        )
        # energie : dW/dt ~ KpW.(W0 - W) + M_W, un second ordre en W
        wn_w = 3.0 / (p["zeta_W"] * p["tau_W"])
        auto["KpW"], auto["KiW"] = 2 * p["zeta_W"] * wn_w, wn_w**2
        wn_pll = 3.0 / (p["zeta_PLL"] * p["tau_PLL"])
        auto["Kipll"], auto["Kppll"] = wn_pll**2 / wb, 2 * p["zeta_PLL"] * wn_pll / wb
        for k, v in auto.items():
            p.setdefault(k, v)
        self._auto_kpdc = "Kpdc" not in params
        p.setdefault("Kpdc", 1.0)
        p.setdefault("Kidc", 0.0)
        p.setdefault("Kiq", 100.0)
        p.setdefault("Kpq", 0.0)
        p.setdefault("kVdc", 1.0)
        if abs(float(p["kVdc"]) - 1.0) > 1e-12:
            # kVdc ne pondere que le COURANT rendu au noeud, pas la tension
            # que lisent les boucles internes (``vdc_m``) : les deux ne sont
            # homogenes que pour kVdc = 1. La garde tables._verifier_bases_dc
            # l'impose deja a l'assemblage ; on la tient aussi ici pour qu'un
            # MMC construit a la main ne fabrique pas un modele faux. Le
            # parametre reste expose, comme celui de GFM3DC/GFL2DC (meme
            # famille de ports, memes fiches).
            raise ValueError(
                f"{name}: kVdc = {p['kVdc']:g}, mais le modele moyen ne pese "
                "que le courant continu, pas la tension de ses boucles : "
                "seul kVdc = 1 (base de tension continue du convertisseur "
                "egale a celle du noeud) est represente."
            )
        p.setdefault("bloque", 0.0)
        p["kS"] = Sn / Sb
        p["wb"] = wb
        self.tuning = {"tau_VL": p["tau_VL"]}
        Element.__init__(self, name, [bus, bus_dc], **{k: p[k] for k in self.params})
        self.Sn, self.Sb = float(Sn), float(Sb)

    @staticmethod
    def _equations(s, regule: bool):
        ct, st = sp.cos(s.theta_pll - s.theta_g), sp.sin(s.theta_pll - s.theta_g)
        vgd = s.vd_t * ct + s.vq_t * st  # global -> repere PLL
        vgq = -s.vd_t * st + s.vq_t * ct
        w_pll = s.M_pll + s.Kppll * vgq + s.wff

        # boucles externes : tension continue (ou statisme) et reactif
        qe = vgq * s.isd - vgd * s.isq
        pe = vgd * s.isd + vgq * s.isq
        isd_ref = s.Kpdc * (s.vdc_f - s.vdc_ref) + s.M_d
        isq_ref = s.Kpq * (qe - s.q_ref) + s.M_q
        isd_ref, isq_ref, aw_d, aw_q = limiter_courant(s, isd_ref, isq_ref)

        # cote alternatif : la tension que fabriquent les bras
        e_d = (
            s.KpCL * (isd_ref - s.isd)
            + s.M_CLd
            + s.Kffv * vgd
            - s.wff * s.L1_pu * s.isq
        )
        e_q = (
            s.KpCL * (isq_ref - s.isq)
            + s.M_CLq
            + s.Kffv * vgq
            + s.wff * s.L1_pu * s.isd
        )
        p_ac = e_d * s.isd + e_q * s.isq

        # cote continu : la tension des bras, et ce qui la fixe ; bloque, la
        # station est decouplee de son noeud
        vdc_m = s.v_dc * (1 - s.bloque) + s.bloque
        integrateurs = []
        if regule:
            idc_ref = pe + s.KpW * (s.W0 - s.W) + s.M_W
            u_z = s.KpZ * (idc_ref - s.idc) + s.M_z
            v_z = vdc_m - u_z
            integrateurs = [s.KiZ * (idc_ref - s.idc), s.KiW * (s.W0 - s.W)]
        else:
            v_z = sp.sqrt(s.W / s.W0)
        ouvert = 1 - s.bloque

        f = [
            (s.wb / s.L1_pu) * (e_d - vgd - s.R1_pu * s.isd + w_pll * s.L1_pu * s.isq),
            (s.wb / s.L1_pu) * (e_q - vgq - s.R1_pu * s.isq - w_pll * s.L1_pu * s.isd),
            (s.wb / s.Ldc_pu) * ((vdc_m - v_z) * ouvert - s.Rdc_pu * s.idc)
            - s.bloque * s.idc / T_BLOCAGE,
            v_z * s.idc * ouvert - p_ac,
            s.Kidc * (s.vdc_f - s.vdc_ref) + aw_d,
            s.Kiq * (qe - s.q_ref) + aw_q,
            s.KiCL * (isd_ref - s.isd),
            s.KiCL * (isq_ref - s.isq),
            *integrateurs,
            (s.v_dc - s.vdc_f) / s.T_vdc,
            s.Kipll * vgq,
            s.wb * (w_pll - s.w_ref),
        ]
        iD = s.kS * (s.isd * ct - s.isq * st)
        iQ = s.kS * (s.isd * st + s.isq * ct)
        # le courant d'une inductance : continu, meme au blocage
        i_dc = -s.kS * s.kVdc * s.idc
        Vt = sp.sqrt(vgd**2 + vgq**2)
        h = [iD, iQ, i_dc, pe, qe, w_pll, Vt, s.v_dc, s.theta_pll, s.W]
        return f, [], h

    @staticmethod
    def equations(s):
        return MMC._equations(s, regule=False)

    def _capacite_vue(self) -> float:
        """La capacite que la boucle de tension continue voit au noeud, pu
        station : celle du noeud, plus celle des sous-modules si l'energie
        est libre (``C_eq/wb = 2.W0``)."""
        p = self.pdict
        return p["Cdc_pu"] + (0.0 if self.regule else 2 * p["W0"] * p["wb"])

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        """Equilibre analytique exact : PLL calee sur la tension du bus, et le
        courant continu qui equilibre l'energie, ``v_z.idc = p_ac`` avec
        ``v_z = v_dc - Rdc.idc`` ; ``W = W0`` si l'energie est regulee,
        ``W = W0.v_z^2`` sinon."""
        p = self.pdict
        w0 = ctx.w
        V = ctx.v["t"]
        vdc = float(np.real(ctx.v["dc"]))
        S = ctx.s.get("t", 0j)
        Ig = (np.conj(S / V) if abs(V) > 1e-12 else 0j) / p["kS"]
        angle_v = float(np.angle(V))
        vgd = abs(V)
        i = Ig * np.exp(-1j * angle_v)
        isd, isq = i.real, i.imag
        L, R = p["L1_pu"], p["R1_pu"]
        e_d = vgd + R * isd - w0 * L * isq
        e_q = R * isq + w0 * L * isd
        p_ac = e_d * isd + e_q * isq
        # Rdc.idc^2 - vdc.idc + p_ac = 0, la racine proche de p_ac/vdc
        Rdc = p["Rdc_pu"]
        idc = 2 * p_ac / (vdc + np.sqrt(max(vdc**2 - 4 * Rdc * p_ac, 0.0)))
        v_z = vdc - Rdc * idc
        q0 = -vgd * isq
        if self.kdroop and self.droop_courant:
            # statisme courant (droopI) : meme pente dp/dv_dc qu'au load flow
            d = self.droop_courant
            pente = d["p_ref"] / d["vdc_ref"] + (2 * vdc - d["vdc_ref"]) / self.kdroop
            self.set_param("Kpdc", pente / vgd)
        elif self.kdroop:
            # statisme : dp = vgd.disd = dv_dc/kdroop
            self.set_param("Kpdc", 1.0 / (self.kdroop * vgd))
        elif self._auto_kpdc:
            # la capacite vue, et la conductance negative d'une puissance
            # constante (idc = p/v_dc)
            Geff = -idc / vdc
            self.set_param(
                "Kpdc",
                (3 * self._capacite_vue() / (self.tuning["tau_VL"] * p["wb"]) - Geff)
                / (vgd / vdc),
            )
        init = {
            "isd": isd,
            "isq": isq,
            "idc": idc,
            "W": p["W0"] if self.regule else p["W0"] * v_z**2,
            "M_d": isd,
            "M_q": isq,
            "M_CLd": e_d - p["Kffv"] * vgd + p["wff"] * L * isq,
            "M_CLq": e_q - p["wff"] * L * isd,
            "vdc_f": vdc,
            "M_pll": w0 - p["wff"],
            "theta_pll": angle_v + ctx.theta_g,
            "vdc_ref": vdc,
            "q_ref": q0,
        }
        if self.regule:
            init.update(M_z=Rdc * idc, M_W=idc - vgd * isd)
        return init


class MMCEnergie(MMC):
    """MMC, modele moyen a ENERGIE REGULEE, port continu (13 etats).

    Les boucles de courant continu (``M_z``) et d'energie (``M_W``) en plus :
    voir la docstring du module. ``modele = "mmc_energie"``.
    """

    regule = True
    params = MMC.params[:7] + ("KpZ", "KiZ", "KpW", "KiW") + MMC.params[7:]
    states = MMC.states[:8] + ("M_z", "M_W") + MMC.states[8:]

    @staticmethod
    def equations(s):
        return MMC._equations(s, regule=True)
