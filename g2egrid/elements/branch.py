"""Branches serie R-L dynamiques : lignes et transformateurs.

Le courant de branche est un ETAT (di/dt), pas une relation algebrique :
c'est ce qui rend le modele EMT/dq et donne les modes electriques rapides.

Equations dans le repere dq global (pu, X = w0.L en pu) :

    di_d/dt = (wb/X) . ( vd_f/n - vd_t - R.i_d + w.X.i_q )
    di_q/dt = (wb/X) . ( vq_f/n - vq_t - R.i_q - w.X.i_d )

``i`` circule de ``f`` vers ``t`` (cote t). Courants INJECTES :
port f : -i/n   (le transfo ideal conserve la puissance)   port t : +i
"""

from __future__ import annotations

import numpy as np

from g2egrid.noyau.element import Element, InitContext

__all__ = ["Branch", "Line", "LineDC", "Transformer"]


class Branch(Element):
    """Branche R-L dynamique entre deux bus, avec rapport de transformation."""

    passif = True

    model_name = "Branch"
    ports = ("f", "t")
    params = ("R", "X", "n", "wb")
    states = ("i_d", "i_q")
    svc_outputs = ("p_f", "q_f")

    @staticmethod
    def equations(s):
        vf_d, vf_q = s.vd_f / s.n, s.vq_f / s.n
        f = [
            (s.wb / s.X) * (vf_d - s.vd_t - s.R * s.i_d + s.w_ref * s.X * s.i_q),
            (s.wb / s.X) * (vf_q - s.vq_t - s.R * s.i_q - s.w_ref * s.X * s.i_d),
        ]
        # courants injectes : -i/n cote f, +i cote t
        ifd, ifq = -s.i_d / s.n, -s.i_q / s.n
        h = [
            ifd,
            ifq,  # port f
            s.i_d,
            s.i_q,  # port t
            s.vd_f * ifd
            + s.vq_f * ifq,  # p_f (injecte au bus f, donc < 0 si transit f->t)
            s.vq_f * ifd - s.vd_f * ifq,  # q_f
        ]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        R, X, n = self.pdict["R"], self.pdict["X"], self.pdict["n"]
        i = (ctx.v["f"] / n - ctx.v["t"]) / complex(R, X)
        return {"i_d": i.real, "i_q": i.imag}


class Line(Branch):
    """Ligne de transmission R-L (rapport unitaire)."""

    model_name = "Branch"

    def __init__(
        self, name: str, bus_from: int, bus_to: int, R: float, X: float, wb: float
    ):
        super().__init__(name, [bus_from, bus_to], R=R, X=X, n=1.0, wb=wb)


class LineDC(Element):
    """Cable continu R-L : le courant est un etat, comme cote alternatif.

    Un seul axe, donc une seule equation et un seul etat -- pas de terme de
    couplage ``w.X.i``, qui n'existe que parce qu'un repere dq tourne.

    La capacite du cable N'EST PAS ici : elle se repartit par moities dans la
    capacite des deux noeuds, exactement comme ``tables`` le fait deja des
    charges capacitives des lignes alternatives. Un troncon de cable se lit
    donc comme un pi dont les deux demi-capacites ont ete versees aux bus.
    """

    passif = True

    model_name = "LineDC"
    ports = ()
    ports_dc = ("f", "t")
    params = ("R", "X", "wb")
    states = ("i",)
    svc_outputs = ("p_f",)

    def __init__(
        self, name: str, bus_from: int, bus_to: int, R: float, X: float, wb: float
    ):
        if float(R) <= 0.0:
            raise ValueError(
                f"{name}: R = {R} ; un cable continu sans resistance laisse son "
                "courant indetermine a l'equilibre"
            )
        super().__init__(name, [bus_from, bus_to], R=R, X=X, wb=wb)

    @staticmethod
    def equations(s):
        f = [(s.wb / s.X) * (s.v_f - s.v_t - s.R * s.i)]
        # courants injectes : -i cote f, +i cote t (i circule de f vers t)
        h = [-s.i, s.i, -s.v_f * s.i]
        return f, [], h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        vf, vt = float(np.real(ctx.v["f"])), float(np.real(ctx.v["t"]))
        return {"i": (vf - vt) / self.pdict["R"]}


class Transformer(Branch):
    """Transformateur R-L avec rapport de prise ``n`` (pu/pu, cote ``f``)."""

    model_name = "Branch"

    def __init__(
        self,
        name: str,
        bus_from: int,
        bus_to: int,
        R: float,
        X: float,
        wb: float,
        n: float = 1.0,
    ):
        super().__init__(name, [bus_from, bus_to], R=R, X=X, n=n, wb=wb)
