"""Bus infini : portage FIDELE de ``matlab/G2Elin_V2/core/IB.m``.

F.e.m. constante ``E/_delta`` derriere une impedance serie R-L, plus un etat
d'ANGLE qui integre l'ecart de frequence de consigne :

    di_d/dt = (wb/X).(e_d - v_gd - R.i_d) + wb.w.i_q
    di_q/dt = (wb/X).(e_q - v_gq - R.i_q) - wb.w.i_d
    dtheta/dt = wb.(w_set - w_ref)
    e_d = E.cos(delta + theta) ,  e_q = E.sin(delta + theta)

ETATS (3) : i_d, i_q, theta        -- comme IB.m
ENTREES de service : E_ref, delta_ref, w_set (frequence imposee par la source)

Cote MATLAB, ``theta`` est sorti comme ``theta_slack`` et sert de reference
d'angle a tout le reseau ; ici le repere global est fixe, et c'est l'angle de
la f.e.m. qui porte l'integrale de la frequence. C'est la meme physique : dans
les deux cas l'angle de la source relatif au reste du reseau est
``integrale(w_set - w_reseau)``, et le mode marginal correspondant (lambda = 0)
apparait des deux cotes.

ECART ASSUME AVEC IB.m : le fichier MATLAB pose ``wb = obj.W_ref`` (= 1 pu)
la ou la physique demande ``wb = 2.pi.f``. Les commentaires du code MATLAB
("couplage ... en rad/s", "Delta_omega_ref en p.u. * wb = rad/s") montrent que
c'est bien ``2.pi.f`` qui etait vise. Le portage utilise la valeur physique :
avec ``wb = 1`` la dynamique electrique de la source serait 314 fois trop
lente (pole a -R/X = -0,1 s^-1 au lieu de -31,4 s^-1).
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from g2egrid.noyau.element import Element, InitContext

__all__ = ["IB", "SourceDC", "ConvertisseurDCDC", "VoltageSource", "IB_DEFAULT"]

# Valeurs par defaut de la table T_IBs des reseaux MATLAB.
IB_DEFAULT = dict(R=0.001, X=0.01)


class IB(Element):
    """Bus infini : f.e.m. constante derriere R-L, avec etat d'angle."""

    ports = ("t",)
    params = ("R", "X", "kS", "wb")
    states = ("i_d", "i_q", "theta")
    svc_inputs = ("E_ref", "delta_ref", "w_set")
    svc_outputs = ("P", "Q", "Vt", "theta_out", "omega")

    def __init__(
        self,
        name: str,
        bus: int,
        R: float = IB_DEFAULT["R"],
        X: float = IB_DEFAULT["X"],
        wb: float = 2 * np.pi * 50,
        Sn: float | None = None,
        Sb: float | None = None,
    ):
        """``R``/``X`` en pu base ``Sn`` (base systeme si ``Sn``/``Sb`` omis)."""
        kS = 1.0 if (Sn is None or Sb is None) else Sn / Sb
        super().__init__(name, bus, R=R, X=X, kS=kS, wb=wb)
        self.Sn = float(Sn) if Sn is not None else None
        self.Sb = float(Sb) if Sb is not None else None

    @staticmethod
    def equations(s):
        # L'angle de la f.e.m. est repere par rapport au repere GLOBAL : quand
        # la source est elle-meme la reference d'angle, theta_g vaut theta et
        # la f.e.m. ne tourne pas, comme dans symIB (V1). Quand elle ne l'est
        # pas, elle tourne de theta - theta_g par rapport au reste du reseau.
        ang = s.delta_ref + s.theta - s.theta_g
        ed = s.E_ref * sp.cos(ang)
        eq = s.E_ref * sp.sin(ang)
        f = [
            (s.wb / s.X) * (ed - s.vd_t - s.R * s.i_d) + s.wb * s.w_ref * s.i_q,
            (s.wb / s.X) * (eq - s.vq_t - s.R * s.i_q) - s.wb * s.w_ref * s.i_d,
            s.wb * (s.w_set - s.w_ref),
        ]
        iD, iQ = s.kS * s.i_d, s.kS * s.i_q
        h = [
            iD,
            iQ,
            s.vd_t * iD + s.vq_t * iQ,  # P injecte (base systeme)
            s.vq_t * iD - s.vd_t * iQ,  # Q injecte
            sp.sqrt(s.vd_t**2 + s.vq_t**2),
            s.theta,
            s.w_set,
        ]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        V = ctx.v["t"]
        S = ctx.s.get("t", 0j)  # puissance INJECTEE
        courant = (np.conj(S / V) if abs(V) > 1e-12 else 0j) / self.pdict["kS"]
        E = V + complex(self.pdict["R"], self.pdict["X"]) * courant
        return {
            "i_d": courant.real,
            "i_q": courant.imag,
            "theta": 0.0,
            "E_ref": abs(E),
            "delta_ref": float(np.angle(E)),
            "w_set": ctx.w,
        }


class SourceDC(Element):
    """Source continue : f.e.m. constante derriere R-L. L'analogue DC de :class:`IB`.

        di/dt = (wb/X).(E_ref - v_t - R.i)

    Un seul etat, un seul port : ni axe q, ni etat d'angle -- il n'y a pas de
    frequence a integrer sur un reseau continu.

    LA F.E.M. EST IMPOSEE, pas deduite d'un load flow. Le lot 3 s'arrete a un
    Newton sur les seules tensions continues ; c'est le Newton unifie AC/DC du
    lot 4 qui saura attribuer une puissance a cette source et en tirer ``E``.
    D'ici la, ``E`` est une donnee, et c'est ``v_t`` que le reseau resout.

    Ce modele est aussi l'EQUIVALENT EXACT de la boucle de courant DC que
    ``GFM3``/``GFM4`` portaient en interne. Leur seconde equation continue ::

        didc/dt = (idc_ref - idc)/Tdc ,  idc_ref = p_ref/vdc_ref + Kpdc.(vdc_ref - vdc)

    se factorise par ``Kpdc`` en ::

        didc/dt = (Kpdc/Tdc).[ p_ref/(Kpdc.vdc_ref) + vdc_ref - vdc - idc/Kpdc ]

    terme a terme la forme ci-dessus, avec ``R = 1/Kpdc``,
    ``X = wb.Tdc/Kpdc`` et ``E = vdc_ref + p_ref/(Kpdc.vdc_ref)``. C'est ce
    qui rend le decoupage "VSC + bus DC + source" spectralement identique au
    ``GFM3`` monolithique, et ``tests/test_dc.py`` le verifie.
    """

    ports = ()
    ports_dc = ("t",)
    params = ("R", "X", "kS", "wb")
    states = ("i",)
    svc_inputs = ("E_ref",)
    svc_outputs = ("P", "Vt")

    def __init__(
        self,
        name: str,
        bus: int,
        R: float = IB_DEFAULT["R"],
        X: float = IB_DEFAULT["X"],
        wb: float = 2 * np.pi * 50,
        E: float = 1.0,
        Sn: float | None = None,
        Sb: float | None = None,
    ):
        """``R``/``X`` en pu base ``Sn`` (base systeme si ``Sn``/``Sb`` omis)."""
        if float(R) <= 0.0:
            raise ValueError(
                f"{name}: R = {R} ; une f.e.m. sans resistance serie laisse son "
                "courant indetermine a l'equilibre"
            )
        kS = 1.0 if (Sn is None or Sb is None) else Sn / Sb
        super().__init__(name, bus, R=R, X=X, kS=kS, wb=wb)
        self.E = float(E)
        self.Sn = float(Sn) if Sn is not None else None
        self.Sb = float(Sb) if Sb is not None else None

    @staticmethod
    def equations(s):
        f = [(s.wb / s.X) * (s.E_ref - s.v_t - s.R * s.i)]
        iT = s.kS * s.i
        h = [iT, s.v_t * iT, s.v_t]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        v = float(np.real(ctx.v["t"]))
        return {"i": (self.E - v) / self.pdict["R"], "E_ref": self.E}


class ConvertisseurDCDC(Element):
    """Convertisseur DC/DC raccordant une batterie ou un champ PV a un bus continu.

    Le convertisseur est vu du bus continu par sa boucle de courant, du
    premier ordre (``Tc``), et par la consigne que sa commande lui donne ::

        Tc . di/dt = i_ref - i ,   i_ref = (P0 + dp)/v_t + Kv.(v_ref - v_t)

    Deux cas, par le seul gain ``Kv`` :

    * **PV** (``Kv = 0``) -- le suivi du point de puissance maximale tient la
      puissance ``P0`` quelle que soit la tension : une source a puissance
      constante, le pendant producteur de :class:`~g2egrid.elements.ChargeDC_P` ;
    * **batterie** (``Kv > 0``) -- le convertisseur SOUTIENT la tension du bus
      par un statisme ``Kv`` [pu de courant par pu de tension] autour de sa
      consigne ``P0`` : il participe au reglage de la tension continue, et
      peut servir de reference a une ile continue.

    ``P0`` et le courant sont en pu base systeme. La batterie elle-meme (etat
    de charge) n'est pas modelisee : a l'echelle des modes etudies, elle est
    une source de tension ideale derriere le convertisseur.
    """

    ports = ()
    ports_dc = ("t",)
    params = ("P0", "Kv", "vref", "Tc")
    states = ("i",)
    svc_inputs = ("dp",)
    svc_outputs = ("p",)

    def __init__(
        self,
        name: str,
        bus: int,
        P: float,
        Kv: float = 0.0,
        Tc: float = 5e-3,
        vref: float = 1.0,
    ):
        if float(Tc) <= 0.0:
            raise ValueError(f"{name}: Tc = {Tc} ; constante de temps non physique")
        super().__init__(
            name, bus, P0=float(P), Kv=float(Kv), vref=float(vref), Tc=float(Tc)
        )

    @property
    def batterie(self) -> bool:
        return self.pdict["Kv"] > 0.0

    #: une batterie a statisme tient la tension de son ile continue
    @property
    def reference_tension_dc(self) -> bool:
        return self.batterie

    @staticmethod
    def equations(s):
        i_ref = (s.P0 + s.dp) / s.v_t + s.Kv * (s.vref - s.v_t)
        return [(i_ref - s.i) / s.Tc], [], [s.i, s.v_t * s.i]

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        v = float(np.real(ctx.v["t"]))
        p = self.pdict
        return {"i": p["P0"] / v + p["Kv"] * (p["vref"] - v), "dp": 0.0}

    def scheduled_injection_dc(self) -> float | None:
        # le PV connait sa puissance ; celle de la batterie depend de la tension
        return None if self.batterie else self.pdict["P0"]

    def caler_dc(self, v: float, p: float) -> None:
        """Au slack continu, la batterie fournit le reliquat ``p`` a la tension
        ``v`` imposee : sa consigne ``P0`` s'en deduit (``i = p/v``)."""
        self.set_param("P0", p - self.pdict["Kv"] * (self.pdict["vref"] - v) * v)


# Nom historique du portage.
VoltageSource = IB
