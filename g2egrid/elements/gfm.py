"""Convertisseurs grid-forming : portage FIDELE de ``matlab/G2Elin_V2/core/GFM1..4.m``.

Quatre variantes, memes equations et memes valeurs par defaut que MATLAB,
leurs deux decoupages a port continu, et deux autres lois de synchronisation :

=========== ============================================ ======
Classe      Description                                  Etats
=========== ============================================ ======
GFM1        statisme 1er ordre, bus DC ideal               13
GFM2        statisme 2e ordre (etat ``Dw``), DC ideal      14
GFM3        statisme 1er ordre, bus DC dynamique           15
GFM4        statisme 2e ordre, bus DC dynamique            16
GFM3DC      statisme 1er ordre, PORT continu               13
GFM4DC      statisme 2e ordre, PORT continu                14
GFMdVOC     oscillateur virtuel dispatchable, DC ideal     14
GFMMatching appariement ``w = vdc``, DC dynamique          14
=========== ============================================ ======

Les deux dernieres se choisissent par le reglage ``loi`` (``"dvoc"`` sur une
ligne de T_GFM1s/T_GFM2s, ``"matching"`` sur T_GFM3s/T_GFM4s), voir
:data:`LOIS_GFM`. ``GFM3DC`` et ``GFM4DC`` sont les memes que ``GFM3``/``GFM4``, mais leurs deux
etats continus sont sortis dans le reseau : le condensateur devient un
``BusDC``, la boucle de courant DC une ``SourceDC``. Voir
:class:`_GFMDCBase`.

Structure commune (identique a MATLAB) :

* **filtre LCL** (6 etats) : ``isd/isq`` cote convertisseur, ``igd/igq`` cote
  reseau, ``ved/veq`` aux bornes du condensateur ;
* **statisme** : ``w = w_ref + mp.(p_ref - pm)`` (1er ordre) ou etat ``Dw``
  filtre a ``wc`` (2e ordre) ; ``ved_ref = ve_ref + nq.(q_ref - qm)`` ;
* **boucles internes** : PI de tension (``M_VLd/q``) puis PI de courant
  (``M_CLd/q``), avec decouplage ``wff`` et feedforward ``Kffi``/``Kffv`` ;
* **impedance virtuelle** ``Rv_pu``, ``Lv_pu`` (nulle par defaut) : la
  consigne de tension perd la chute virtuelle, dynamique (``zv_dyn = 1``) ou
  quasi statique (``zv_dyn = 0``) -- de quoi stabiliser un GFM face a un
  reseau fort ;
* **bus DC** (GFM3/4) : ``vdc``, ``idc`` avec regulation ``Kpdc`` calee au
  point de fonctionnement (correction CPL, comme MATLAB).

Les gains PI sont auto-regles a partir des temps de reponse
(``tau_VL``, ``tau_CL``, ``tau_DC``), exactement comme
``computeControllerGains()`` cote MATLAB.

Repere : l'axe d local est en avance de ``theta`` sur l'axe D global. MATLAB
ecrit la meme rotation avec ``theta_g - theta`` (interface (V,0) + angle).
Le modele MATLAB integre l'angle ABSOLU (``dtheta = wb.w``) et exclut cette
ligne de l'equilibre ; en repere global tournant a ``w_ref`` la meme equation
s'ecrit ``dtheta = wb.(w - w_ref)``, qui s'annule bien a l'equilibre. Les deux
donnent la MEME matrice A (seule la colonne ``w_ref`` de B differe).
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from g2egrid.noyau.element import Element, InitContext
from g2egrid.elements.limiteurs import (
    _EPS2,
    REGLAGES_BUS_CONTINU,
    LimiteurCourantTension,
    depassement_vi,
    impedance_virtuelle,
    limiter_courant,
    limiter_tension,
)

__all__ = [
    "GFM1",
    "GFM2",
    "GFM3",
    "GFM4",
    "GFM3DC",
    "GFM4DC",
    "GFMdVOC",
    "GFMMatching",
    "LOIS_GFM",
    "GFM_DEFAULT",
    "GFM_DC_DEFAULT",
]

# Valeurs par defaut : reprise LITTERALE des proprietes de GFM1.m / GFM2.m.
GFM_DEFAULT = dict(
    R1_pu=0.01,
    L1_pu=0.1,
    C_pu=0.1,
    R2_pu=0.003,
    L2_pu=0.01,
    mp=0.001,
    nq=0.001,
    wf_p=0.1 * 2 * np.pi * 50,
    wf_q=0.1 * 2 * np.pi * 50,
    wff=1.0,
    Kffi=1.0,
    Kffv=1.0,
    tau_VL=15e-3,
    zeta_VL=0.707,
    tau_CL=0.1e-3,
    zeta_CL=0.707,
    wc=None,  # GFM2/4 : defaut 0.1*wb (fixe au constructeur)
    # limitation de courant (pu de l'equipement) et son drapeau -- voir
    # g2egrid.elements.limiteurs ; a mode_lim = 0 les equations sont celles de
    # GFM1.m, au bit pres
    Imax=1.2,
    Taw=2e-3,  # constante de temps du retro-calcul (anti-emballement)
    mode_lim=0.0,
    # saturation de tension : |m| <= Mmax (modulation), retro-calcul sur les
    # integrateurs de la boucle de courant. Mmax = 0 : pas de saturation de
    # tension (defaut sur bus continu ideal) ; MMAX = 1,15 pour une
    # modulation vectorielle, le defaut des GFM a bus continu dynamique (voir
    # limiteurs.REGLAGES_BUS_CONTINU). Elle s'active par reglage (GFM_PARAMS,
    # PARAMS_PAR_BUS) : sur un defaut franc, elle s'enclenche avant la
    # limitation de courant et demande un pas d'integration fin (quelques
    # dizaines de us en BDF2).
    Mmax=0.0,
    # Taw_v = None : tau_CL/2, de l'ordre de la boucle de courant (fixe au
    # constructeur). Un retro-calcul lent (20 ms, l'ancien defaut) laisse
    # l'integrateur de courant (KiCL ~ 4e5 /s) s'emballer pendant un defaut :
    # la modulation reste collee a Mmax dans une direction qui ne commande
    # plus le courant, la saturation ne relache plus et le courant monte a
    # plusieurs pu -- quand le Newton ne lache pas avant. Les valeurs
    # intermediaires (2 ms, 0,2 ms) le font lacher : la modulation libre colle
    # au coin du min.
    Taw_v=None,
    mode_vsat=0.0,
    # IMPEDANCE VIRTUELLE Rv + j.Lv (pu de l'equipement) : la consigne de la
    # boucle de tension perd la chute virtuelle Rv.i_g + Lv.(d/dt + j.w).i_g.
    # Vu du reseau, le GFM est une source derriere une impedance de plus : face
    # a un reseau FORT (peu d'impedance), elle adoucit la boucle de puissance
    # (gain P-delta en 1/X) et amortit la synchronisation. Rv ~ Lv/5 amortit
    # en plus le mode du flux de l'inductance totale, vers 50 Hz dans le
    # repere dq. Rv = Lv = 0 : les equations de GFM1.m au bit pres.
    Rv_pu=0.0,
    Lv_pu=0.0,
    # zv_dyn = 1 (defaut) : inductance virtuelle DYNAMIQUE, realisee sur la
    # chute mesuree aux bornes de L2 -- Lv.(d/dt + j.w).i_g vaut exactement
    # (Lv/L2).(v_e - v_g - R2.i_g) --, sans derivee ni etat de plus ; une
    # inductance reelle de Lv en serie ne ferait pas autrement. zv_dyn = 0 :
    # forme QUASI STATIQUE, (Rv + j.w.Lv).i_g, la plus repandue ; elle ignore
    # le d/dt et, avec la boucle de tension par defaut, la DESTABILISE (mode
    # a 60-130 Hz des Lv = 0,05 pu sur un reseau fort) -- a reserver aux
    # boucles de tension lentes. Meme equilibre dans les deux cas.
    zv_dyn=1.0,
    # limitation par IMPEDANCE VIRTUELLE (voir limiteurs.impedance_virtuelle) :
    # au-dela de I_vi (pu de l'equipement), une impedance Zv = Rv + j.Xv,
    # Rv = kRv.ov(|is|), Xv = XsurR_vi.Rv, abaisse la consigne de tension.
    # kRv = 0 : inactive (defaut), equations de GFM1.m au bit pres. dI_vi :
    # largeur du raccord parabolique du seuil (jacobien continu). XsurR_vi
    # est BAS (0,1) et non "eleve" : avec les boucles de tension du MATLAB
    # (tau_VL = 15 ms) une reactance virtuelle algebrique dans le repere dq
    # instabilise le mode LC du filtre (voir limiteurs.impedance_virtuelle).
    kRv=0.0,
    XsurR_vi=0.1,
    I_vi=1.0,
    dI_vi=0.05,
)

# GFM3/GFM4 : cote continu.
GFM_DC_DEFAULT = dict(
    Cdc_pu=0.1, Gdc_pu=17.28e-6, Tdc=0.0001, tau_DC=5e-3, zeta_DC=0.707
)


class _GFMBase(LimiteurCourantTension, Element):
    """Tronc commun des GFM : filtre LCL + boucles internes + statisme."""

    ports = ("t",)
    _base_params = (
        "R1_pu",
        "L1_pu",
        "C_pu",
        "R2_pu",
        "L2_pu",
        "mp",
        "nq",
        "wf_p",
        "wf_q",
        "wff",
        "Kffi",
        "Kffv",
        "KpVL",
        "KiVL",
        "KpCL",
        "KiCL",
        "kS",
        "wb",
        "Imax",
        "Taw",
        "mode_lim",
        "Mmax",
        "Taw_v",
        "mode_vsat",
        "Rv_pu",
        "Lv_pu",
        "zv_dyn",
        "kRv",
        "XsurR_vi",
        "I_vi",
        "dI_vi",
    )
    svc_inputs = ("p_ref", "q_ref", "ve_ref")
    svc_outputs = ("P", "Q", "w", "Vt", "theta_out")

    has_dc = False
    second_order_droop = False

    def __init__(self, name: str, bus: int, Sn: float, Sb: float, wb: float, **params):
        p = dict(GFM_DEFAULT)
        if self.has_dc:
            p.update(GFM_DC_DEFAULT)
            p.update(REGLAGES_BUS_CONTINU)  # saturation de tension active
        p.update(params)
        if p.get("wc") is None:
            p["wc"] = 0.1 * wb  # GFM2.m : wc = 0.1*wb
        if p.get("Taw_v") is None:
            p["Taw_v"] = p["tau_CL"] / 2  # cf. GFM_DEFAULT
        # --- auto-reglage des PI (computeControllerGains de MATLAB) ---
        K0_VL = wb / p["C_pu"]
        p["KpVL"] = 6 / (K0_VL * p["tau_VL"])
        p["KiVL"] = (3 / (p["zeta_VL"] * p["tau_VL"])) ** 2 / K0_VL
        K0_CL = 1.0 / p["R1_pu"]  # vdc0 = 1
        tau0_CL = p["L1_pu"] / (wb * p["R1_pu"])
        p["KpCL"] = (6 * tau0_CL - p["tau_CL"]) / (K0_CL * p["tau_CL"])
        p["KiCL"] = (3 / (p["zeta_CL"] * p["tau_CL"])) ** 2 * tau0_CL / K0_CL
        # Kpdc : recale au point de fonctionnement, SAUF si l'utilisateur en
        # impose un (cas d'un reglage repris d'un autre jeu de donnees).
        self._auto_kpdc = "Kpdc" not in params
        p.setdefault("Kpdc", 0.0)
        p["kS"] = Sn / Sb
        p["wb"] = wb
        if self.ports_dc:
            # Rapport des bases de tension continue, bus / convertisseur. Il ne
            # sert qu'a la variante a port DC : les GFM monolithiques portent
            # leur propre condensateur, donc leur propre base.
            p.setdefault("kVdc", 1.0)
            # convertisseur BLOQUE (manoeuvre "blocage") : 0 en service
            p.setdefault("bloque", 0.0)
        self.tuning = {k: p[k] for k in ("tau_VL", "zeta_VL", "tau_CL", "zeta_CL")}
        if self.has_dc:
            self.tuning.update(tau_DC=p["tau_DC"], zeta_DC=p["zeta_DC"])
        super().__init__(name, bus, **{k: p[k] for k in self.params})
        self.Sn, self.Sb = float(Sn), float(Sb)

    # ------------------------------------------------------------------
    # Bloc commun aux quatre variantes : rotation, LCL, boucles internes.
    # Renvoie (equations LCL, grandeurs utiles) pour la sous-classe.
    # ------------------------------------------------------------------
    @staticmethod
    def _common(s, w, vdc, ved_ref=None):
        # L'angle est repere par rapport au repere GLOBAL, comme symGFM qui
        # tourne de theta_g - theta : theta_g est l'angle de l'element de
        # reference, nul a l'equilibre mais pas sa derivee.
        ct, st = sp.cos(s.theta - s.theta_g), sp.sin(s.theta - s.theta_g)
        vgd = s.vd_t * ct + s.vq_t * st  # repere global -> local
        vgq = -s.vd_t * st + s.vq_t * ct

        p = s.ved * s.igd + s.veq * s.igq
        q = -s.ved * s.igq + s.veq * s.igd
        if ved_ref is None:  # statisme Q-V ; dVOC fournit sa propre consigne
            ved_ref = s.ve_ref + (s.q_ref - s.qm) * s.nq
        # impedance virtuelle (cf. GFM_DEFAULT) : la consigne est la f.e.m.
        # virtuelle, sur l'axe d, moins la chute Rv.i_g + Lv.(d/dt + j.w).i_g.
        # Dynamique (zv_dyn = 1) : Lv.(d/dt + j.w).i_g = (Lv/L2).(v_e - v_g -
        # R2.i_g), la chute sur L2 mise a l'echelle. Quasi statique (0) :
        # j.w.Lv.i_g, a la frequence du GFM.
        kqs = (1 - s.zv_dyn) * w * s.Lv_pu
        kdyn = s.zv_dyn * s.Lv_pu / s.L2_pu
        ved_ref = ved_ref - (
            s.Rv_pu * s.igd - kqs * s.igq + kdyn * (s.ved - vgd - s.R2_pu * s.igd)
        )
        veq_ref = -(
            s.Rv_pu * s.igq + kqs * s.igd + kdyn * (s.veq - vgq - s.R2_pu * s.igq)
        )
        # limitation par impedance virtuelle : la consigne perd en plus Zv.is
        # (nulle sous le seuil I_vi, et partout si kRv = 0)
        dv_d, dv_q = impedance_virtuelle(s, s.isd, s.isq)
        ved_ref = ved_ref - dv_d
        veq_ref = veq_ref - dv_q

        isd_ref = (
            s.KpVL * (ved_ref - s.ved)
            + s.M_VLd
            + s.Kffi * s.igd
            - s.wff * s.C_pu * s.veq
        )
        isq_ref = (
            s.KpVL * (veq_ref - s.veq)
            + s.M_VLq
            + s.Kffi * s.igq
            + s.wff * s.C_pu * s.ved
        )
        # limitation de courant : consigne ramenee dans |i| <= Imax, et les
        # integrateurs de la boucle de TENSION la suivent (retro-calcul)
        isd_ref, isq_ref, aw_d, aw_q = limiter_courant(s, isd_ref, isq_ref)
        md = (1 / vdc) * (
            s.KpCL * (isd_ref - s.isd)
            + s.Kffv * s.ved
            - s.wff * s.L1_pu * s.isq
            + s.M_CLd
        )
        mq = (1 / vdc) * (
            s.KpCL * (isq_ref - s.isq)
            + s.Kffv * s.veq
            + s.wff * s.L1_pu * s.isd
            + s.M_CLq
        )
        # saturation de tension : la modulation ramenee dans |m| <= Mmax, et
        # les integrateurs de la boucle de COURANT la suivent (retro-calcul)
        md, mq, awv_d, awv_q = limiter_tension(s, md, mq)

        lcl = [
            (s.wb / s.L1_pu)
            * (md * vdc - s.ved - s.R1_pu * s.isd + w * s.L1_pu * s.isq),
            (s.wb / s.L1_pu)
            * (mq * vdc - s.veq - s.R1_pu * s.isq - w * s.L1_pu * s.isd),
            (s.wb / s.L2_pu) * (s.ved - vgd - s.R2_pu * s.igd + w * s.L2_pu * s.igq),
            (s.wb / s.L2_pu) * (s.veq - vgq - s.R2_pu * s.igq - w * s.L2_pu * s.igd),
            (s.wb / s.C_pu) * (s.isd - s.igd + w * s.C_pu * s.veq),
            (s.wb / s.C_pu) * (s.isq - s.igq - w * s.C_pu * s.ved),
        ]
        inner = [
            s.KiVL * (ved_ref - s.ved) + aw_d,
            s.KiVL * (veq_ref - s.veq) + aw_q,
            s.KiCL * (isd_ref - s.isd) + awv_d * vdc,
            s.KiCL * (isq_ref - s.isq) + awv_q * vdc,
        ]
        # courant injecte au bus : repere local -> global, base systeme
        iD = s.kS * (s.igd * ct - s.igq * st)
        iQ = s.kS * (s.igd * st + s.igq * ct)
        Vt = sp.sqrt(s.ved**2 + s.veq**2)
        return lcl, inner, (p, q, md, mq, iD, iQ, Vt)

    # ------------------------------------------------------------------
    def _init_common(self, ctx: InitContext, vdc: float = 1.0) -> dict[str, float]:
        """Equilibre analytique EXACT (le repere local s'aligne sur ``ve``).

        ``vdc`` est la tension continue AU POINT DE DEPART. Elle vaut 1 pour
        les variantes qui portent leur propre condensateur -- c'est ce que fait
        MATLAB -- et la tension du bus continu pour celles a port DC, chez qui
        la valeur ne vient plus du modele mais du reseau.
        """
        p = self.pdict
        w0 = ctx.w
        V = ctx.v["t"]  # tension bus (pu systeme)
        S = ctx.s.get("t", 0j)
        Ig = (np.conj(S / V) if abs(V) > 1e-12 else 0j) / p["kS"]  # pu conv.

        # tension aux bornes du condensateur, repere global
        Ve = V + complex(p["R2_pu"], w0 * p["L2_pu"]) * Ig
        # Le repere local s'aligne sur la f.e.m. virtuelle E, derriere les
        # deux impedances virtuelles : la constante (Rv_pu, Lv_pu), dont la
        # chute a l'equilibre vaut (Rv + j.w.Lv).i_g sous ses deux formes, et
        # celle de LIMITATION, Zv.is, nulle sous le seuil I_vi. C'est E que la
        # boucle de tension met sur l'axe d. Sans impedance virtuelle, E = Ve
        # et veq = 0 exactement (delta = 0). |is| ne depend pas du repere : la
        # chute se calcule dans le repere global, sans iteration. L'angle est
        # compte depuis le REPERE GLOBAL, qui tourne de theta_g (cf.
        # GFL._init_common) : sans theta_g, le Newton devait le rattraper, et
        # pouvait tomber sur une autre branche.
        Is = Ig + 1j * w0 * p["C_pu"] * Ve
        Zv = 0j
        if p["kRv"] > 0.0:
            i_s = float(np.sqrt(abs(Is) ** 2 + _EPS2))  # meme regularisation
            ov = float(depassement_vi(i_s, p["I_vi"], p["dI_vi"]))
            Rv = p["kRv"] * ov
            Zv = complex(Rv, p["XsurR_vi"] * Rv)
        E = Ve + complex(p["Rv_pu"], w0 * p["Lv_pu"]) * Ig
        if Zv != 0:
            E = E + Zv * Is
        angle_e = float(np.angle(E))
        theta = angle_e + ctx.theta_g
        rot = np.exp(-1j * angle_e)
        delta = float(np.angle(Ve)) - angle_e
        ved, veq = abs(Ve) * np.cos(delta), abs(Ve) * np.sin(delta)
        ig = Ig * rot
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
            "pm": p0,
            "qm": q0,
            "theta": theta,
            "M_VLd": isd - p["Kffi"] * igd + p["wff"] * p["C_pu"] * veq,
            "M_VLq": isq - p["Kffi"] * igq - p["wff"] * p["C_pu"] * ved,
            "M_CLd": md * vdc - p["Kffv"] * ved + p["wff"] * p["L1_pu"] * isq,
            "M_CLq": mq * vdc - p["Kffv"] * veq - p["wff"] * p["L1_pu"] * isd,
            "p_ref": p0,
            "q_ref": q0,
            # f.e.m. virtuelle : la consigne de la boucle de tension
            "ve_ref": abs(E),
            "_md": md,
            "_mq": mq,
            "_vdc": vdc,
        }


class GFM1(_GFMBase):
    """GFM a statisme 1er ordre, bus DC ideal (13 etats) — cf ``GFM1.m``."""

    params = _GFMBase._base_params
    states = (
        "isd",
        "isq",
        "igd",
        "igq",
        "ved",
        "veq",
        "pm",
        "qm",
        "theta",
        "M_VLd",
        "M_VLq",
        "M_CLd",
        "M_CLq",
    )

    @staticmethod
    def equations(s):
        w = s.mp * (s.p_ref - s.pm) + s.w_ref
        lcl, inner, (p, q, md, mq, iD, iQ, Vt) = _GFMBase._common(s, w, 1)
        f = (
            lcl
            + [s.wf_p * (p - s.pm), s.wf_q * (q - s.qm), s.wb * (w - s.w_ref)]
            + inner
        )
        h = [iD, iQ, p, q, w, Vt, s.theta]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        return self._init_common(ctx)


class GFM2(_GFMBase):
    """GFM a statisme 2e ordre, bus DC ideal (14 etats) — cf ``GFM2.m``."""

    second_order_droop = True
    params = _GFMBase._base_params + ("wc",)
    states = (
        "isd",
        "isq",
        "igd",
        "igq",
        "ved",
        "veq",
        "pm",
        "Dw",
        "theta",
        "qm",
        "M_VLd",
        "M_VLq",
        "M_CLd",
        "M_CLq",
    )

    @staticmethod
    def equations(s):
        w = s.Dw + s.w_ref
        lcl, inner, (p, q, md, mq, iD, iQ, Vt) = _GFMBase._common(s, w, 1)
        f = (
            lcl
            + [
                s.wf_p * (p - s.pm),
                s.mp * s.wc * (s.p_ref - s.pm) - s.wc * s.Dw,
                s.wb * (w - s.w_ref),
                s.wf_q * (q - s.qm),
            ]
            + inner
        )
        h = [iD, iQ, p, q, w, Vt, s.theta]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = self._init_common(ctx)
        init["Dw"] = 0.0
        return init


class GFM3(_GFMBase):
    """GFM a statisme 1er ordre, bus DC dynamique (15 etats) — cf ``GFM3.m``."""

    has_dc = True
    params = _GFMBase._base_params + ("Cdc_pu", "Gdc_pu", "Tdc", "Kpdc")
    states = (
        "isd",
        "isq",
        "igd",
        "igq",
        "ved",
        "veq",
        "vdc",
        "idc",
        "pm",
        "theta",
        "qm",
        "M_VLd",
        "M_VLq",
        "M_CLd",
        "M_CLq",
    )
    svc_inputs = ("p_ref", "q_ref", "ve_ref", "vdc_ref")
    svc_outputs = ("P", "Q", "w", "Vt", "theta_out", "vdc_out")

    @staticmethod
    def equations(s):
        w = s.mp * (s.p_ref - s.pm) + s.w_ref
        lcl, inner, (p, q, md, mq, iD, iQ, Vt) = _GFMBase._common(s, w, s.vdc)
        idc_ref = s.p_ref / s.vdc_ref + s.Kpdc * (s.vdc_ref - s.vdc)
        dc = [
            (s.wb / s.Cdc_pu) * (s.idc - s.Gdc_pu * s.vdc - md * s.isd - mq * s.isq),
            (idc_ref - s.idc) / s.Tdc,
        ]
        f = (
            lcl
            + dc
            + [s.wf_p * (p - s.pm), s.wb * (w - s.w_ref), s.wf_q * (q - s.qm)]
            + inner
        )
        h = [iD, iQ, p, q, w, Vt, s.theta, s.vdc]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = self._init_common(ctx)
        p = self.pdict
        md, mq, vdc = init.pop("_md"), init.pop("_mq"), init.pop("_vdc")
        init["vdc"] = vdc
        init["idc"] = p["Gdc_pu"] * vdc + md * init["isd"] + mq * init["isq"]
        init["vdc_ref"] = vdc
        # Reglage CPL-aware de Kpdc au point de fonctionnement (GFM3.m)
        if self._auto_kpdc:
            Geff = p["Gdc_pu"] - init["pm"] / vdc**2
            self.set_param(
                "Kpdc", 3 * p["Cdc_pu"] / (self.tuning["tau_DC"] * p["wb"]) - Geff
            )
        # La source continue doit fournir idc, pertes du filtre comprises, et
        # non p_ref/vdc (la puissance au condensateur alternatif) : vdc_ref se
        # cale pour que idc_ref = p_ref/vdc_ref + Kpdc.(vdc_ref - vdc) = idc.
        # Sans ce calage, l'ecart (les pertes R1, R2) etait laisse au Newton
        # local, qui ne le resorbe pas toujours sous 1e-6.
        Kpdc, idc, pref = self.pdict["Kpdc"], init["idc"], init["p_ref"]
        # Kpdc.vr^2 - b.vr + pref = 0, b = Kpdc.vdc + idc. Deux racines : on
        # garde celle qui est au voisinage de vdc (la branche physique). Avec
        # Kpdc > 0 c'est la racine "+" ; avec Kpdc < 0 -- batterie en charge
        # (P < -0,19 pu), ou la correction CPL de Kpdc change de signe -- c'est
        # la racine "-" : l'ancien code ne calait pas vdc_ref dans ce cas, et
        # laissait au Newton local un residu de 7e-10.
        b = Kpdc * vdc + idc
        if Kpdc > 1e-12:  # inchange : racine "+"
            init["vdc_ref"] = (b + np.sqrt(max(b * b - 4 * Kpdc * pref, 0.0))) / (
                2 * Kpdc
            )
        elif Kpdc < -1e-12:
            disc = np.sqrt(max(b * b - 4 * Kpdc * pref, 0.0))
            racines = ((b + disc) / (2 * Kpdc), (b - disc) / (2 * Kpdc))
            init["vdc_ref"] = min(racines, key=lambda r: abs(r - vdc))
        elif abs(idc) > 1e-12:
            init["vdc_ref"] = pref / idc
        return init


class GFM4(_GFMBase):
    """GFM a statisme 2e ordre, bus DC dynamique (16 etats) — cf ``GFM4.m``."""

    has_dc = True
    second_order_droop = True
    params = _GFMBase._base_params + ("wc", "Cdc_pu", "Gdc_pu", "Tdc", "Kpdc")
    states = (
        "isd",
        "isq",
        "igd",
        "igq",
        "ved",
        "veq",
        "vdc",
        "idc",
        "pm",
        "Dw",
        "theta",
        "qm",
        "M_VLd",
        "M_VLq",
        "M_CLd",
        "M_CLq",
    )
    svc_inputs = ("p_ref", "q_ref", "ve_ref", "vdc_ref")
    svc_outputs = ("P", "Q", "w", "Vt", "theta_out", "vdc_out")

    @staticmethod
    def equations(s):
        w = s.Dw + s.w_ref
        lcl, inner, (p, q, md, mq, iD, iQ, Vt) = _GFMBase._common(s, w, s.vdc)
        idc_ref = s.p_ref / s.vdc_ref + s.Kpdc * (s.vdc_ref - s.vdc)
        dc = [
            (s.wb / s.Cdc_pu) * (s.idc - s.Gdc_pu * s.vdc - md * s.isd - mq * s.isq),
            (idc_ref - s.idc) / s.Tdc,
        ]
        f = (
            lcl
            + dc
            + [
                s.wf_p * (p - s.pm),
                s.mp * s.wc * (s.p_ref - s.pm) - s.wc * s.Dw,
                s.wb * (w - s.w_ref),
                s.wf_q * (q - s.qm),
            ]
            + inner
        )
        h = [iD, iQ, p, q, w, Vt, s.theta, s.vdc]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = GFM3.initialize(self, ctx)
        init["Dw"] = 0.0
        return init


class GFMdVOC(_GFMBase):
    """GFM a OSCILLATEUR VIRTUEL dispatchable (dVOC), bus DC ideal (14 etats).

    La boucle de synchronisation n'est plus un statisme mais l'oscillateur
    de Colombino, Gross et Dorfler (dVOC), ecrit en coordonnees polaires dans
    sa forme pour reseau inductif (rotation ``kappa = pi/2``) ::

        w      = w_ref + eta.(p_ref/v_ref^2 - p/Vr^2)
        dVr/dt = wb.eta.[(q_ref/v_ref^2 - q/Vr^2) + alpha.(v_ref^2 - Vr^2)/v_ref^2].Vr

    ``Vr`` est l'amplitude de l'oscillateur ; elle devient la consigne de la
    boucle de tension (``ved_ref = Vr``), le reste -- filtre LCL, boucles de
    tension et de courant, limitation -- est celui des autres GFM. Pres du
    nominal, ``eta`` est la pente du statisme P-f (pu/pu) ; ``alpha`` fixe la
    raideur de l'amplitude : le statisme Q-V equivalent vaut ``1/(2.alpha)``.
    A la difference du statisme, la dynamique de tension est INTRINSEQUE a
    l'oscillateur, et les puissances sont normalisees par ``Vr^2`` : c'est ce
    qui donne a dVOC ses garanties de stabilite presque globale.

    ``p`` et ``q`` sont les puissances filtrees (``wf_p``, ``wf_q``), comme
    pour les autres variantes.
    """

    params = _GFMBase._base_params + ("eta", "alpha")
    states = (
        "isd",
        "isq",
        "igd",
        "igq",
        "ved",
        "veq",
        "pm",
        "qm",
        "theta",
        "Vr",
        "M_VLd",
        "M_VLq",
        "M_CLd",
        "M_CLq",
    )

    def __init__(self, name, bus, Sn, Sb, wb, **params):
        params.setdefault("eta", 0.02)
        params.setdefault("alpha", 5.0)
        super().__init__(name, bus, Sn, Sb, wb, **params)

    @staticmethod
    def equations(s):
        v2 = s.ve_ref**2
        w = s.w_ref + s.eta * (s.p_ref / v2 - s.pm / s.Vr**2)
        lcl, inner, (p, q, md, mq, iD, iQ, Vt) = _GFMBase._common(s, w, 1, s.Vr)
        dVr = (
            s.wb
            * s.eta
            * ((s.q_ref / v2 - s.qm / s.Vr**2) + s.alpha * (v2 - s.Vr**2) / v2)
            * s.Vr
        )
        f = (
            lcl
            + [
                s.wf_p * (p - s.pm),
                s.wf_q * (q - s.qm),
                s.wb * (w - s.w_ref),
                dVr,
            ]
            + inner
        )
        h = [iD, iQ, p, q, w, Vt, s.theta]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = self._init_common(ctx)
        init["Vr"] = init["ve_ref"]  # la f.e.m. virtuelle (= ved sans Zv)
        return init


class GFMMatching(_GFMBase):
    """GFM a COMMANDE PAR APPARIEMENT (matching control), bus DC dynamique
    (14 etats).

    Arghir, Jouini et Dorfler : la frequence du convertisseur est la tension
    continue, a un gain pres ::

        w = w_ref . vdc / vdc_ref

    Le condensateur continu joue le role de l'inertie d'une machine
    synchrone -- son energie ``Cdc.vdc^2/2`` celui de l'energie cinetique --,
    et le regulateur de la source continue celui du regulateur de vitesse :
    ``idc_ref = p_ref/vdc_ref + Kpdc.(vdc_ref - vdc)``, dont le gain fixe le
    statisme P-f (``Delta w / w = Delta P / (Kpdc . vdc_ref)`` en regime
    etabli). Une baisse de la tension continue -- la source primaire n'arrive
    plus a suivre -- se lit donc immediatement en frequence : c'est la
    propriete que le statisme n'a pas.

    La tension alternative reste reglee par statisme Q-V (``nq``) et les
    boucles internes communes. ``Kpdc`` vaut par defaut ``1/mdc``, ``mdc``
    etant le statisme P-f voulu (2 %).
    """

    has_dc = True
    params = _GFMBase._base_params + ("Cdc_pu", "Gdc_pu", "Tdc", "Kpdc")
    states = (
        "isd",
        "isq",
        "igd",
        "igq",
        "ved",
        "veq",
        "vdc",
        "idc",
        "theta",
        "qm",
        "M_VLd",
        "M_VLq",
        "M_CLd",
        "M_CLq",
    )
    svc_inputs = ("p_ref", "q_ref", "ve_ref", "vdc_ref")
    svc_outputs = ("P", "Q", "w", "Vt", "theta_out", "vdc_out")

    def __init__(self, name, bus, Sn, Sb, wb, mdc: float = 0.02, **params):
        params.setdefault("Kpdc", 1.0 / mdc)
        super().__init__(name, bus, Sn, Sb, wb, **params)

    @staticmethod
    def equations(s):
        w = s.w_ref * s.vdc / s.vdc_ref
        lcl, inner, (p, q, md, mq, iD, iQ, Vt) = _GFMBase._common(s, w, s.vdc)
        idc_ref = s.p_ref / s.vdc_ref + s.Kpdc * (s.vdc_ref - s.vdc)
        dc = [
            (s.wb / s.Cdc_pu) * (s.idc - s.Gdc_pu * s.vdc - md * s.isd - mq * s.isq),
            (idc_ref - s.idc) / s.Tdc,
        ]
        f = lcl + dc + [s.wb * (w - s.w_ref), s.wf_q * (q - s.qm)] + inner
        h = [iD, iQ, p, q, w, Vt, s.theta, s.vdc]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = self._init_common(ctx)
        p = self.pdict
        md, mq, vdc = init.pop("_md"), init.pop("_mq"), init.pop("_vdc")
        init.pop("pm")
        init["vdc"] = vdc
        init["idc"] = p["Gdc_pu"] * vdc + md * init["isd"] + mq * init["isq"]
        init["vdc_ref"] = vdc
        # la source continue fournit exactement la puissance appelee
        init["p_ref"] = init["idc"] * vdc
        return init


#: loi de synchronisation -> classe, selon que le bus continu est ideal
#: (tables T_GFM1s/T_GFM2s) ou dynamique (T_GFM3s/T_GFM4s). Choisie par le
#: reglage ``loi`` (GFM_PARAMS ou PARAMS_PAR_BUS) ; absente : statisme.
LOIS_GFM = {
    ("dvoc", False): GFMdVOC,
    ("matching", True): GFMMatching,
}


class _GFMDCBase(_GFMBase):
    """Tronc commun des GFM A PORT CONTINU : le condensateur est dans le bus.

    C'EST LE MEME MODELE QUE ``GFM3``/``GFM4``, AUTREMENT DECOUPE. La version
    monolithique porte deux etats continus de plus, ``vdc`` et ``idc``, et
    referme sur eux-memes le condensateur et la boucle de courant DC ::

        dvdc/dt = (wb/Cdc).(idc - Gdc.vdc - md.isd - mq.isq)
        didc/dt = (idc_ref - idc)/Tdc

    Ici ces deux etats sortent du convertisseur :

    * ``vdc`` DEVIENT l'entree de port ``v_dc``. Le condensateur n'a pas
      disparu, il est devenu le :class:`~g2egrid.elements.bus.BusDC` auquel le
      convertisseur se raccorde -- ``C = Cdc_pu`` -- et la loi des noeuds de ce
      bus redonne la premiere equation a l'identique, a condition que le
      courant renvoye au bus soit bien ``-(Gdc.v_dc + md.isd + mq.isq)``, le
      courant APPELE change de signe ;
    * la boucle de courant DC devient un element de reseau, une
      :class:`~g2egrid.elements.source.SourceDC` dont la docstring donne la
      correspondance exacte des trois parametres.

    Ce qu'on y gagne : plusieurs convertisseurs peuvent partager un meme noeud
    continu, et un cable continu peut les separer. Ce qu'il faut verifier :
    que le decoupage n'a rien change au spectre -- c'est le test d'equivalence
    de ``tests/test_dc.py``.

    ``kVdc = Un_dc / Vdc_n`` rapporte la base de tension du noeud continu a
    celle du convertisseur. RESERVE ASSUMEE : il ne pondere ici que le
    COURANT rendu au bus, pas la tension lue par les boucles internes, qui
    prend ``v_dc`` telle quelle. Les deux ne sont homogenes que pour
    ``kVdc = 1``, la valeur par defaut et la seule que les tests exercent ; le
    cas des bases distinctes reste a trancher avec le lot 4, qui apportera les
    tensions nominales continues.
    """

    has_dc = True
    ports_dc = ("dc",)
    svc_outputs = ("P", "Q", "w", "Vt", "theta_out", "vdc_out")

    #: Un grid-forming tient la TENSION de son bus alternatif -- ``ve_ref`` est
    #: une entree de service, et ``initialize`` la cale sur le module lu au bus.
    #: Son bus est donc PV, jamais PQ, et ``PQ``/``VdcQ`` n'ont pas de sens ici.
    #: Reste ce qu'il fait du cote continu : soit il suit son statisme et sa
    #: puissance est programmee (``PV``), soit il devient la reference de
    #: tension de son ile continue et c'est sa puissance qui devient
    #: l'inconnue (``VdcV``). Il n'a pas de regulateur de ``vdc`` : en
    #: ``VdcV``, la tension continue est une DONNEE du point de fonctionnement,
    #: et l'ile doit etre tenue dynamiquement par ailleurs. Seul formeur de son
    #: ile alternative -- la station d'un parc offshore --, il en est le slack
    #: (``Vtheta``) : sa puissance est celle que l'ile appelle, et le reseau
    #: continu, tenu ailleurs, la transporte.
    modes_dc = ("PV", "VdcV", "Vtheta")

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
        """``bus`` est le noeud alternatif, ``bus_dc`` le noeud continu.

        Les deux arrivent separes plutot que dans une liste : l'ordre des bus
        d'un element a port DC (alternatifs d'abord, continus ensuite) est une
        convention du noyau, pas quelque chose que l'appelant doit retenir.
        """
        super().__init__(name, [bus, int(bus_dc)], Sn, Sb, wb, **params)

    @staticmethod
    def _courant_dc(s, md, mq):
        """Courant INJECTE au noeud continu : l'oppose du courant appele. Nul
        quand le convertisseur est bloque."""
        return (
            -s.kS
            * s.kVdc
            * (s.Gdc_pu * s.v_dc + md * s.isd + mq * s.isq)
            * (1 - s.bloque)
        )

    @staticmethod
    def _vdc_modulation(s):
        """La tension continue que lit la modulation : ``v_dc`` en service,
        1 pu une fois BLOQUE -- le convertisseur ne voit plus son noeud, et sa
        modulation ne divise plus par une tension qui peut passer par zero
        pendant un defaut continu. A ``bloque = 0``, c'est ``v_dc`` au bit
        pres."""
        return s.v_dc * (1 - s.bloque) + s.bloque

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        return self._init_common(ctx, vdc=float(np.real(ctx.v["dc"])))


class GFM3DC(_GFMDCBase):
    """GFM a statisme 1er ordre, port continu (13 etats) — cf :class:`GFM3`."""

    params = _GFMBase._base_params + ("Gdc_pu", "kVdc", "bloque")
    states = (
        "isd",
        "isq",
        "igd",
        "igq",
        "ved",
        "veq",
        "pm",
        "theta",
        "qm",
        "M_VLd",
        "M_VLq",
        "M_CLd",
        "M_CLq",
    )

    @staticmethod
    def equations(s):
        w = s.mp * (s.p_ref - s.pm) + s.w_ref
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


class GFM4DC(_GFMDCBase):
    """GFM a statisme 2e ordre, port continu (14 etats) — cf :class:`GFM4`."""

    second_order_droop = True
    params = _GFMBase._base_params + ("wc", "Gdc_pu", "kVdc", "bloque")
    states = (
        "isd",
        "isq",
        "igd",
        "igq",
        "ved",
        "veq",
        "pm",
        "Dw",
        "theta",
        "qm",
        "M_VLd",
        "M_VLq",
        "M_CLd",
        "M_CLq",
    )

    @staticmethod
    def equations(s):
        w = s.Dw + s.w_ref
        lcl, inner, (p, q, md, mq, iD, iQ, Vt) = _GFMBase._common(
            s, w, _GFMDCBase._vdc_modulation(s)
        )
        f = (
            lcl
            + [
                s.wf_p * (p - s.pm),
                s.mp * s.wc * (s.p_ref - s.pm) - s.wc * s.Dw,
                s.wb * (w - s.w_ref),
                s.wf_q * (q - s.qm),
            ]
            + inner
        )
        h = [iD, iQ, _GFMDCBase._courant_dc(s, md, mq), p, q, w, Vt, s.theta, s.v_dc]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = super().initialize(ctx)
        init["Dw"] = 0.0
        return init
