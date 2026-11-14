"""Les transits du load flow, branche par branche.

Le load flow rend des tensions et des injections nettes aux bus ; ce qu'on
lit sur un schema de reseau, ce sont des TRANSITS -- combien passe dans cette
ligne, dans quel sens, et ce qu'elle en perd. On les recalcule ici a partir
des tensions, avec la meme admittance serie et le meme rapport de
transformation que ``Network.build_ybus`` : deux formules differentes pour la
meme branche donneraient un transit qui ne boucle pas avec les injections.

RESERVE -- la susceptance des lignes n'est PAS dans la branche : g2egrid la
range par moities dans la capacite des bus (``Cbus``), ou elle se confond avec
celle des autres lignes du meme noeud. Les transits rendus sont donc ceux de
la branche R-L seule, comme ceux de :func:`g2egrid.temporel.mesures` ; ils
different d'un transit de load flow classique de ``|v|^2.b/2`` a chaque bout.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from g2egrid.reseau.network import Network

__all__ = ["Transit", "flux_branches"]


@dataclass
class Transit:
    """Le transit d'une branche, en pu systeme, compte a partir de chaque bout.

    ``P_de`` est la puissance qui ENTRE dans la branche par le bus ``de`` ;
    ``P_de + P_vers`` vaut donc les pertes joule de la branche.
    """

    nom: str
    de: int
    vers: int
    P_de: float
    Q_de: float
    P_vers: float
    Q_vers: float

    @property
    def pertes(self) -> float:
        return self.P_de + self.P_vers


def flux_branches(net: Network, lf=None) -> list[Transit]:
    """Les transits de toutes les branches alternatives du reseau."""
    lf = lf if lf is not None else net.lf
    if lf is None:
        raise ValueError("load flow non resolu : appeler net.run_loadflow() d'abord")
    idx = net._bus_index()
    out = []
    for el in net.elements:
        if len(el.ports) != 2 or not all(b in idx for b in el.buses):
            continue
        f, t = el.buses
        Vf, Vt = lf.V[idx[f]], lf.V[idx[t]]
        y = 1.0 / complex(el.pdict["R"], el.pdict["X"])
        n = el.pdict.get("n", 1.0)
        i_f = y / (n * n) * Vf - y / n * Vt
        i_t = y * Vt - y / n * Vf
        Sf, St = Vf * np.conj(i_f), Vt * np.conj(i_t)
        out.append(Transit(el.name, f, t, Sf.real, Sf.imag, St.real, St.imag))
    return out
