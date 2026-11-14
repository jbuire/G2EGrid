"""L'inertie d'un reseau : celle que les reglages promettent, et celle que la
simulation montre.

1. INERTIE EQUIVALENTE D'UNE BOUCLE DE PUISSANCE -- Kelada, Buire et
   Hadjsaid [1] ramenent les lois de synchronisation des grid-forming (statisme
   filtre, synchronverter/VSM, dVOC, appariement) a la forme canonique de
   l'equation d'oscillation d'une machine ::

        2H . d(dw)/dt = P* - p - D . dw          (dw en pu de w_b)

   et en deduisent, par identification, l'inertie ``H`` et l'amortissement
   ``D`` que chaque reglage represente (leur table 1) :

   ===================  ==========================  ===================
   loi                  2H                          D
   ===================  ==========================  ===================
   statisme filtre      1/(m_p.w_f)                 1/m_p
   VSM                  J.w*                        D_p.w*
   dVOC                 e*^2/(eta.w_f)              e*^2/eta
   appariement          v_dc*.K_p,dc/(w_f.K_theta)  v_dc*.K_p,dc/K_theta
   ===================  ==========================  ===================

   :func:`inertie_equivalente` applique ces formules aux modeles de g2egrid
   (en pu, ``w`` en pu de ``w_b``) ; :func:`reglage_equivalent` fait
   l'inverse -- le reglage de chaque loi qui donne un couple ``(H, D)`` --,
   ce qui permet de comparer les lois A INERTIE EGALE.

   Deux adaptations aux modeles de g2egrid. Le statisme du SECOND ordre (GFM2,
   GFM4) filtre deux fois -- la mesure de puissance (``w_fp``) puis l'ecart
   de pulsation (``w_c``) : au premier ordre, ses deux constantes de temps
   s'ajoutent, ``2H = (1/w_c + 1/w_fp)/m_p``. L'appariement de g2egrid n'a
   pas de filtre de mesure sur ``v_dc`` : c'est le condensateur continu qui en
   tient lieu, ``2H = C_dc/w_b`` (``C_dc`` en pu « w_b.C »), ``D = K_pdc``.

2. INERTIE EFFECTIVE, FENETRE PAR FENETRE -- :func:`inertie_glissante`
   estime sur une trajectoire l'inertie que chaque groupe d'unites fournit
   reellement, sur des fenetres successives (0,1 s par defaut) ::

        H_k = f0 . integrale(P_inertielle dt) / (2 . (f(t_a) - f(t_b)))   [MW.s]

   ``f`` est la frequence du centre d'inertie des machines (celle du reseau),
   et la puissance INERTIELLE d'une unite est ce qu'elle fournit au-dela de
   sa reponse primaire : ``P_e - P_m`` pour une machine (``P_m`` : sa
   turbine), ``P_e - (p* - D.dw)`` pour un GFM (son statisme). Pour un GFM
   libre, ``H_k`` vaut son inertie equivalente ; LIMITE en courant, il ne
   suit plus sa propre loi -- sa puissance est plafonnee -- et ``H_k`` dit
   ce qui reste de son inertie virtuelle.

[1] F. Kelada, J. Buire, N. Hadjsaid, « Comparative study of decentralized
    grid-forming converter controls for inverter-based microgrids », Archivos
    de Ciencia e Investigacion, 2025 (hal-05262196).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = [
    "Inertie",
    "inertie_equivalente",
    "inertie_reseau",
    "reglage_equivalent",
    "inertie_glissante",
    "sorties_inertie",
]


@dataclass
class Inertie:
    """L'inertie ``H`` [s] et l'amortissement ``D`` [pu/pu] equivalents d'une
    unite, sur sa puissance ``S`` [MVA], et la formule qui les donne."""

    nom: str
    loi: str
    H: float
    D: float
    S: float
    formule: str

    @property
    def HS(self) -> float:
        """Energie inertielle [MW.s]."""
        return self.H * self.S


def _classe(el) -> str:
    return type(el).__name__


def inertie_equivalente(el, v_ref: float = 1.0) -> Inertie | None:
    """L'inertie equivalente d'un element (``None`` si ce n'est pas une
    unite de production). ``v_ref`` : amplitude de consigne du dVOC [pu]."""
    c, p = _classe(el), el.pdict
    S = getattr(el, "Sn", 0.0) / 1e6
    if c.startswith("SynchronousMachine"):
        return Inertie(el.name, "machine", p["H"], p.get("KD", 0.0), S, "H")
    if c in ("GFM1", "GFM3", "GFM3DC"):
        H = 1.0 / (2 * p["mp"] * p["wf_p"])
        return Inertie(el.name, "statisme", H, 1 / p["mp"], S, "2H = 1/(m_p.w_fp)")
    if c in ("GFM2", "GFM4", "GFM4DC"):
        H = (1.0 / p["wc"] + 1.0 / p["wf_p"]) / (2 * p["mp"])
        return Inertie(
            el.name, "statisme 2e ordre", H, 1 / p["mp"], S, "2H = (1/w_c + 1/w_fp)/m_p"
        )
    if c == "GFMdVOC":
        v2 = v_ref**2
        H = v2 / (2 * p["eta"] * p["wf_p"])
        return Inertie(el.name, "dVOC", H, v2 / p["eta"], S, "2H = v*^2/(eta.w_fp)")
    if c == "GFMMatching":
        H = p["Cdc_pu"] / (2 * p["wb"])
        return Inertie(el.name, "appariement", H, p["Kpdc"], S, "2H = C_dc/w_b")
    if c.startswith(("GFL", "MMC")):
        return Inertie(el.name, "suiveur", 0.0, 0.0, S, "aucune (suit le reseau)")
    return None


def inertie_reseau(net) -> dict:
    """L'inertie equivalente de chaque unite du reseau, et les totaux [MW.s]
    par genre (machines, grid-forming)."""
    unites = [i for el in net.elements if (i := inertie_equivalente(el)) is not None]
    sm = sum(i.HS for i in unites if i.loi == "machine")
    gfm = sum(i.HS for i in unites if i.loi not in ("machine", "suiveur"))
    return {"unites": unites, "H_sm": sm, "H_gfm": gfm, "H_total": sm + gfm}


def reglage_equivalent(
    loi: str, H: float, D: float, wb: float = 2 * math.pi * 50, **fixes
) -> dict:
    """Le reglage de la boucle de puissance ``loi`` qui donne l'inertie ``H``
    [s] et l'amortissement ``D`` [pu/pu] -- la table 2 de [1], en pu.

    - ``statisme`` (GFM1) : ``m_p = 1/D``, ``w_fp = D/(2H)`` ;
    - ``statisme2`` (GFM2) : ``m_p = 1/D``, ``w_fp`` fixe (``fixes``, 0,5 w_b
      par defaut), ``w_c`` complete l'inertie ;
    - ``dvoc`` : ``eta = v*^2/D``, ``w_fp = D/(2H)`` ;
    - ``appariement`` : ``K_pdc = D``, ``C_dc = 2H.w_b``.
    """
    if H <= 0 or D <= 0:
        raise ValueError("H et D doivent etre positifs")
    if loi == "statisme":
        return {"mp": 1 / D, "wf_p": D / (2 * H)}
    if loi == "statisme2":
        wf_p = fixes.get("wf_p", 0.5 * wb)
        reste = 2 * H / D - 1 / wf_p
        if reste <= 0:
            raise ValueError(
                f"H = {H:g} s inaccessible avec w_fp = {wf_p:g} rad/s (minimum {1 / (2 * wf_p) * D:g} s)"
            )
        return {"mp": 1 / D, "wf_p": wf_p, "wc": 1 / reste}
    if loi == "dvoc":
        v2 = fixes.get("v_ref", 1.0) ** 2
        return {"loi": "dvoc", "eta": v2 / D, "wf_p": D / (2 * H)}
    if loi == "appariement":
        return {"loi": "matching", "Kpdc": D, "Cdc_pu": 2 * H * wb}
    raise ValueError(f"loi {loi!r} inconnue (statisme, statisme2, dvoc, appariement)")


# ---------------------------------------------------------------- mesure
def sorties_inertie(net) -> list[str]:
    """Les sorties a demander a la simulation pour :func:`inertie_glissante`."""
    out = []
    for el in net.elements:
        i = inertie_equivalente(el)
        if i is None or i.loi == "suiveur":
            continue
        # une machine : la puissance d'ENTREFER (couple electrique), pas celle
        # aux bornes, que les pertes et la resistance de causalite faussent
        out.append(f"{el.name}.Ce" if i.loi == "machine" else f"{el.name}.P")
    return out


def _consigne(el, nom: str) -> float:
    return float(el.u0[list(type(el).svc_inputs).index(nom)])


def inertie_glissante(
    tr,
    net,
    t0: float,
    fenetre: float = 0.1,
    duree: float = 3.0,
    f0: float | None = None,
    seuil_hz: float = 1e-3,
) -> dict:
    """L'inertie effective par fenetres de ``fenetre`` secondes apres
    l'evenement ``t0``, par groupe (machines, grid-forming, total) et par
    unite [MW.s]. Une fenetre ou la frequence varie de moins de
    ``seuil_hz`` (autour du nadir) n'a pas d'inertie mesurable : NaN.

    ``tr`` doit porter les sorties de :func:`sorties_inertie` ; ``net`` est
    le reseau simule, a son point d'equilibre (``el.u0``)."""
    f0 = f0 or net.f0
    t = tr.t
    unites = {
        i.nom: (i, el)
        for el in net.elements
        if (i := inertie_equivalente(el)) is not None
    }
    # frequence du centre d'inertie des machines [Hz] ; un groupe declenche
    # (NaN) en sort
    num = np.zeros_like(t)
    den = np.zeros_like(t)
    for nom, (i, el) in unites.items():
        if i.loi == "machine" and f"{nom}.Dwr" in tr.etats:
            w = tr.etat(f"{nom}.Dwr")
            ok = np.isfinite(w)
            num += np.where(ok, i.H * i.S * np.nan_to_num(w), 0.0)
            den += np.where(ok, i.H * i.S, 0.0)
    if not np.any(den > 0):
        raise ValueError("pas de machine synchrone : frequence du reseau indefinie")
    f = num / den * f0

    def puissance_inertielle(nom, i, el):
        if i.loi == "suiveur":
            return None
        if i.loi == "machine":
            # 2H.dw/dt = Pm/w - Ce - KD.dw : la puissance inertielle est
            # Ce.w - Pm + KD.dw (w = 1 + dw)
            Ce = tr.sorties.get(f"{nom}.Ce")
            if Ce is None:
                return None
            dw = np.nan_to_num(tr.etat(f"{nom}.Dwr"))
            Pm = np.nan_to_num(tr.etat(f"{nom}.Pm"))
            return (np.nan_to_num(Ce) * (1 + dw) - Pm + i.D * dw) * i.S
        P = tr.sorties.get(f"{nom}.P")
        if P is None:
            return None
        P = np.nan_to_num(P)
        p_ref = _consigne(el, "p_ref")
        return (P - (p_ref - i.D * f / f0)) * i.S

    pin = {
        n: v
        for n, (i, el) in unites.items()
        if (v := puissance_inertielle(n, i, el)) is not None
    }
    groupes = {
        "machines": [n for n in pin if unites[n][0].loi == "machine"],
        "gfm": [n for n in pin if unites[n][0].loi != "machine"],
    }
    bornes = np.arange(t0, t0 + duree + 1e-9, fenetre)
    out = {
        "t": [],
        "f": [],
        "machines": [],
        "gfm": [],
        "total": [],
        "unites": {n: [] for n in pin},
    }
    # l'inertie CUMULEE depuis l'evenement : l'energie inertielle fournie
    # depuis t0, rapportee a la chute de frequence depuis t0 -- elle lisse
    # les oscillations electromecaniques qu'une fenetre de 0,1 s attrape
    cumul = {"machines": [], "gfm": [], "total": []}
    i0 = int(np.searchsorted(t, t0))
    for a, b in zip(bornes[:-1], bornes[1:]):
        sel = (t >= a) & (t <= b)
        if sel.sum() < 2:
            continue
        ia, ib = np.flatnonzero(sel)[[0, -1]]
        df = f[ia] - f[ib]
        out["t"].append(float(0.5 * (a + b) - t0))
        out["f"].append(float(f[ib]))

        def H(noms):
            if abs(df) < seuil_hz:
                return float("nan")
            e = sum(np.trapezoid(pin[n][sel], t[sel]) for n in noms)
            return float(f0 * e / (2 * df))

        for g, noms in groupes.items():
            out[g].append(H(noms))
        out["total"].append(H(list(pin)))
        sel0 = (t >= t0) & (t <= b)
        df0 = f[i0] - f[ib]
        for g, noms in (*groupes.items(), ("total", list(pin))):
            e = sum(np.trapezoid(pin[n][sel0], t[sel0]) for n in noms)
            cumul[g].append(
                float(f0 * e / (2 * df0)) if abs(df0) >= seuil_hz else float("nan")
            )
        for n in pin:
            out["unites"][n].append(H([n]))
    out["cumul"] = cumul
    out["nominal"] = {
        "machines": sum(unites[n][0].HS for n in groupes["machines"]),
        "gfm": sum(unites[n][0].HS for n in groupes["gfm"]),
    }
    out["nominal"]["total"] = out["nominal"]["machines"] + out["nominal"]["gfm"]
    return out
