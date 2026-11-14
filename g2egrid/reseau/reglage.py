"""Regler par INTENTION : temps de reponse et amortissement, donnees constructeur.

DEUX CONVERSIONS, dans les deux sens.

**Boucles des convertisseurs.** Les gains des PI des GFM et GFL ne se donnent
pas a la main : ils se deduisent, au constructeur de chaque element, d'un
temps de reponse ``tau`` et d'un amortissement ``zeta`` par placement de poles
(``wn = 3/(zeta.tau)``, les formules de ``computeControllerGains`` du MATLAB).
:func:`boucles` rend, pour un jeu de reglages, chaque boucle avec son
intention ET les gains qui en resultent -- LUS sur l'element construit, pas
recalcules ici : c'est le modele qui fait foi. :func:`intention_depuis_gains`
fait le chemin inverse, pour un jeu de donnees qui arrive avec des gains.

**Machines synchrones.** Le modele est ecrit en inductances et resistances
d'enroulements (``Lad``, ``Lfd``, ``Rfd``...), la forme de Park ; une fiche
constructeur, ou un livre, donne des reactances et des constantes de temps
(``Xd``, ``X'd``, ``T'd0``...). :func:`machine_depuis_standard` convertit l'une
en l'autre avec les definitions CLASSIQUES de Kundur (chap. 4) ::

    X'd  = Ll + Lad.Lfd/(Lad + Lfd)          T'd0  = (Lad + Lfd)/(wb.Rfd)
    X''d = Ll + 1/(1/Lad + 1/Lfd + 1/L1d)    T''d0 = (L1d + Lad.Lfd/(Lad+Lfd))/(wb.R1d)

(et de meme sur l'axe q avec ``L1q``, ``L2q``), et :func:`standard_depuis_machine`
fait l'inverse. Verification : le jeu par defaut du paquet est celui de
l'exemple 3.1 de Kundur, et redonne les valeurs de son exemple 4.1 (X'd =
0,30, X''d = 0,23, T'd0 = 8,07 s, T''d0 = 0,03 s, X'q = 0,65, X''q = 0,25).

Idee reprise du portage de F. Kelada (``tuning.py``, et la conversion de son
preset Kundur) ; les formules sont celles du livre.
"""

from __future__ import annotations

import math

__all__ = [
    "BOUCLES",
    "STANDARD",
    "boucles",
    "intention_depuis_gains",
    "machine_depuis_standard",
    "standard_depuis_machine",
]

#: Les boucles reglables, famille par famille : (nom, cle de tau, cle de zeta
#: ou None, gains qu'elle produit).
BOUCLES: dict[str, list[tuple[str, str, str | None, tuple[str, ...]]]] = {
    "GFM": [
        ("tension", "tau_VL", "zeta_VL", ("KpVL", "KiVL")),
        ("courant", "tau_CL", "zeta_CL", ("KpCL", "KiCL")),
    ],
    "GFL": [
        ("courant", "tau_CL", "zeta_CL", ("KpCL", "KiCL")),
        ("PLL", "tau_PLL", "zeta_PLL", ("Kppll", "Kipll")),
        ("puissance active", "tau_PL", None, ("Kip",)),
        ("puissance reactive", "tau_QL", None, ("Kiq",)),
    ],
}

#: Les donnees standard d'une machine, dans l'ordre d'une fiche.
STANDARD = (
    "Xd",
    "Xq",
    "Xl",
    "Xd_p",
    "Xd_pp",
    "Xq_p",
    "Xq_pp",
    "Td0_p",
    "Td0_pp",
    "Tq0_p",
    "Tq0_pp",
    "Ra",
)


def _element(famille: str, reglages: dict, f0: float):
    """Un element de la famille, construit avec ces reglages (Sn = Sb)."""
    from g2egrid.elements import GFL1, GFM1

    classe = {"GFM": GFM1, "GFL": GFL1}[famille]
    wb = 2 * math.pi * f0
    return classe("reglage", 1, Sn=1.0, Sb=1.0, wb=wb, **reglages)


def boucles(famille: str, reglages: dict | None = None, f0: float = 50.0) -> list[dict]:
    """Chaque boucle : son intention (tau, zeta) et les gains qui en sortent.

    ``reglages`` : ce qu'on donnerait a ``GFM_PARAMS`` ou ``PARAMS_PAR_BUS``
    (les cles absentes prennent les valeurs par defaut du modele).
    """
    if famille not in BOUCLES:
        raise ValueError(f"famille {famille!r} : {sorted(BOUCLES)}")
    reglages = dict(reglages or {})
    el = _element(famille, reglages, f0)
    from g2egrid.elements.gfl import GFL_DEFAULT
    from g2egrid.elements.gfm import GFM_DEFAULT

    defauts = GFM_DEFAULT if famille == "GFM" else GFL_DEFAULT
    out = []
    for nom, cle_tau, cle_zeta, gains in BOUCLES[famille]:
        tau = float(reglages.get(cle_tau, defauts.get(cle_tau)))
        zeta = (
            float(reglages.get(cle_zeta, defauts.get(cle_zeta))) if cle_zeta else None
        )
        # la PLL et le courant ont un zeta ; les boucles de puissance du GFL
        # sont du premier ordre (Ki = 3/tau), sans amortissement a regler
        out.append(
            {
                "nom": nom,
                "cle_tau": cle_tau,
                "cle_zeta": cle_zeta,
                "tau": tau,
                "zeta": zeta,
                "wn": 3.0 / (zeta * tau) if zeta else 3.0 / tau,
                "gains": {g: float(el.pdict[g]) for g in gains},
            }
        )
    return out


def intention_depuis_gains(
    famille: str, gains: dict, reglages: dict | None = None, f0: float = 50.0
) -> dict:
    """``tau`` et ``zeta`` de chaque boucle dont on donne les deux gains.

    L'inverse exact des formules du constructeur, pour les plantes que le
    modele suppose (filtre et bus continu de ``reglages``). Rend un dict
    ``{cle_tau: ..., cle_zeta: ...}`` pret a remettre dans les reglages.
    """
    reglages = dict(reglages or {})
    from g2egrid.elements.gfl import GFL_DEFAULT
    from g2egrid.elements.gfm import GFM_DEFAULT

    p = {**(GFM_DEFAULT if famille == "GFM" else GFL_DEFAULT), **reglages}
    wb = 2 * math.pi * f0
    out: dict[str, float] = {}

    def pi_ordre1(Kp, Ki, K0, tau0):
        # Kp = (6 tau0 - tau)/(K0 tau) ; Ki = (3/(zeta tau))^2 . tau0/K0
        tau = 6 * tau0 / (K0 * Kp + 1)
        zeta = 3 / (tau * math.sqrt(Ki * K0 / tau0))
        return tau, zeta

    if {"KpCL", "KiCL"} <= set(gains):
        K0, tau0 = 1.0 / p["R1_pu"], p["L1_pu"] / (wb * p["R1_pu"])
        out["tau_CL"], out["zeta_CL"] = pi_ordre1(
            gains["KpCL"], gains["KiCL"], K0, tau0
        )
    if famille == "GFM" and {"KpVL", "KiVL"} <= set(gains):
        # Kp = 6/(K0 tau) ; Ki = (3/(zeta tau))^2 / K0, K0 = wb/C
        K0 = wb / p["C_pu"]
        tau = 6 / (K0 * gains["KpVL"])
        out["tau_VL"] = tau
        out["zeta_VL"] = 3 / (tau * math.sqrt(gains["KiVL"] * K0))
    if famille == "GFL" and {"Kppll", "Kipll"} <= set(gains):
        # Kipll = wn^2/wb ; Kppll = 2 zeta wn / wb
        wn = math.sqrt(gains["Kipll"] * wb)
        zeta = gains["Kppll"] * wb / (2 * wn)
        out["tau_PLL"], out["zeta_PLL"] = 3 / (zeta * wn), zeta
    if famille == "GFL":
        for cle_g, cle_t in (("Kip", "tau_PL"), ("Kiq", "tau_QL")):
            if cle_g in gains:
                out[cle_t] = 3.0 / gains[cle_g]
    return out


# ------------------------------------------------------------------ machines
def standard_depuis_machine(p: dict, f0: float = 50.0) -> dict:
    """Les donnees standard (reactances, constantes de temps en s) d'un jeu
    de parametres de :class:`~g2egrid.elements.SynchronousMachine`."""
    wb = 2 * math.pi * f0
    Ll, Lad, Laq = p["Ll"], p["Lad"], p["Laq"]
    Lfd, L1d, L1q, L2q = p["Lfd"], p["L1d"], p["L1q"], p["L2q"]
    par_d = Lad * Lfd / (Lad + Lfd)
    par_q = Laq * L1q / (Laq + L1q)
    return {
        "Xd": Ll + Lad,
        "Xq": Ll + Laq,
        "Xl": Ll,
        "Xd_p": Ll + par_d,
        "Xd_pp": Ll + 1 / (1 / Lad + 1 / Lfd + 1 / L1d),
        "Xq_p": Ll + par_q,
        "Xq_pp": Ll + 1 / (1 / Laq + 1 / L1q + 1 / L2q),
        "Td0_p": (Lad + Lfd) / (wb * p["Rfd"]),
        "Td0_pp": (L1d + par_d) / (wb * p["R1d"]),
        "Tq0_p": (Laq + L1q) / (wb * p["R1q"]),
        "Tq0_pp": (L2q + par_q) / (wb * p["R2q"]),
        "Ra": p["Ra"],
    }


def machine_depuis_standard(d: dict, f0: float = 50.0) -> dict:
    """Les parametres du modele de Park depuis une fiche standard.

    ``d`` porte les cles de :data:`STANDARD` (reactances en pu base machine,
    constantes de temps a circuit ouvert en secondes). Leve si la fiche n'est
    pas physique : il faut ``Xd > X'd > X''d > Xl`` et de meme sur l'axe q.
    """
    manquantes = [k for k in STANDARD if k not in d]
    if manquantes:
        raise ValueError(f"donnees manquantes : {manquantes}")
    wb = 2 * math.pi * f0
    Ll = d["Xl"]
    for axe, X, Xp, Xpp in (
        ("d", d["Xd"], d["Xd_p"], d["Xd_pp"]),
        ("q", d["Xq"], d["Xq_p"], d["Xq_pp"]),
    ):
        if not X > Xp > Xpp > Ll > 0:
            raise ValueError(
                f"axe {axe} : il faut X > X' > X'' > Xl > 0 (recu {X:g}, {Xp:g}, "
                f"{Xpp:g}, {Ll:g})"
            )
    for k in ("Td0_p", "Td0_pp", "Tq0_p", "Tq0_pp"):
        if d[k] <= 0:
            raise ValueError(f"{k} = {d[k]:g} s : une constante de temps est > 0")
    Lad, Laq = d["Xd"] - Ll, d["Xq"] - Ll
    # X' = Ll + Lad.Lf/(Lad + Lf)  =>  Lf = Lad.(X' - Ll)/(Lad - (X' - Ll))
    a = d["Xd_p"] - Ll
    Lfd = Lad * a / (Lad - a)
    b = d["Xq_p"] - Ll
    L1q = Laq * b / (Laq - b)
    # X'' = Ll + 1/(1/Lad + 1/Lf + 1/L1)  =>  1/L1 = 1/(X'' - Ll) - 1/Lad - 1/Lf
    L1d = 1 / (1 / (d["Xd_pp"] - Ll) - 1 / Lad - 1 / Lfd)
    L2q = 1 / (1 / (d["Xq_pp"] - Ll) - 1 / Laq - 1 / L1q)
    par_d = Lad * Lfd / (Lad + Lfd)
    par_q = Laq * L1q / (Laq + L1q)
    return {
        "Ra": d["Ra"],
        "Ll": Ll,
        "Lad": Lad,
        "Laq": Laq,
        "Lfd": Lfd,
        "L1d": L1d,
        "L1q": L1q,
        "L2q": L2q,
        "Rfd": (Lad + Lfd) / (wb * d["Td0_p"]),
        "R1d": (L1d + par_d) / (wb * d["Td0_pp"]),
        "R1q": (Laq + L1q) / (wb * d["Tq0_p"]),
        "R2q": (L2q + par_q) / (wb * d["Tq0_pp"]),
    }
