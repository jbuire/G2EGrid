"""Les fiches d'equipement : les equations d'un modele, LISIBLES et ordonnees.

POURQUOI -- :mod:`g2egrid.elements.equations` rend ce que SymPy compile : c'est exact,
mais tout y est developpe. Une consigne de courant limitee apparait
recopiee, avec son ``Min(1, Imax/|i*|)``, dans chacune des equations qui
la lisent ; ``d isd/dt`` d'un GFM occupe dix lignes. Une fiche ecrit le meme
modele comme on le dessine : par SOUS-SYSTEME (repere, synchronisation,
boucles externes, boucle de courant, onduleur et filtre, bus continu), avec
ses signaux intermediaires nommes (``i_sd*``, ``m_d``, ``omega``...).

LA FICHE NE PEUT PAS DERIVER DU MODELE. Ses equations sont ecrites ici, a la
main, mais :func:`verifier` substitue chaque signal intermediaire par sa
definition et compare le resultat au modele compile, etat par etat, variable
algebrique par variable algebrique, sortie par sortie : la difference doit
etre nulle, exactement (``tests/test_fiches.py``). La fiche decrit le
fonctionnement LIBRE -- limitation de courant inactive (``mode_lim = 0``),
convertisseur en service (``bloque = 0``) -- et le dit ; les modes discrets
sont decrits en texte.

Un modele sans fiche ecrite a la main en recoit une generique : ses
equations compilees, rangees en etats, variables algebriques et sorties.

- :func:`fiche` : la fiche d'une classe, en MathML (affichage navigateur),
  LaTeX et texte ;
- :func:`verifier` : l'ecart entre la fiche et le modele compile.
"""

from __future__ import annotations

import contextlib
import inspect
import re
from dataclasses import dataclass, field

import sympy as sp

from g2egrid.elements.equations import _classe, _symboles

__all__ = ["FICHES", "fiche", "mathml_symbole", "symbole_affiche", "verifier"]


# ------------------------------------------------------------------ symboles
#: noms de modele -> nom d'affichage (sympy le met en indices/exposants)
_NOMS = {
    "isd": "i_sd",
    "isq": "i_sq",
    "igd": "i_gd",
    "igq": "i_gq",
    "ved": "v_ed",
    "veq": "v_eq",
    "vd_t": "v_D",
    "vq_t": "v_Q",
    "vdc": "v_dc",
    "v_dc": "v_dc",
    "vdc_f": "v_dc,f",
    "idc": "i_dc",
    "pm": "p_m",
    "qm": "q_m",
    "Dw": "Δω",
    "w": "ω",
    "w_ref": "ω_ref",
    "w_pll": "ω_pll",
    "wb": "ω_b",
    "wff": "ω_ff",
    "wf_p": "ω_fp",
    "wf_q": "ω_fq",
    "wc": "ω_c",
    "theta": "θ",
    "theta_g": "θ_g",
    "theta_pll": "θ_pll",
    "C_pu": "C_f",
    "kS": "k_S",
    "kVdc": "k_Vdc",
    "mp": "m_p",
    "nq": "n_q",
    "md": "m_d",
    "mq": "m_q",
    "Vr": "V_r",
    "W0": "W_0",
    "Imax": "I_max",
    "Taw": "T_aw",
    "Tdc": "T_dc",
    "T_vdc": "T_vdc",
    "id_t": "i_D",
    "iq_t": "i_Q",
    "i_dc": "i_dc,port",
    "eta": "η",
    "alpha": "α",
    "mode_lim": "mode_lim",
}

#: ce que designe chaque symbole (etats, parametres, signaux des fiches)
DESCRIPTIONS = {
    "isd": "courant cote onduleur (axe d)",
    "isq": "courant cote onduleur (axe q)",
    "igd": "courant injecte au reseau (axe d)",
    "igq": "courant injecte au reseau (axe q)",
    "ved": "tension du condensateur de filtre (axe d)",
    "veq": "tension du condensateur de filtre (axe q)",
    "vd_t": "tension du bus, repere global (axe D)",
    "vq_t": "tension du bus, repere global (axe Q)",
    "v_gd": "tension du bus dans le repere du convertisseur (d)",
    "v_gq": "tension du bus dans le repere du convertisseur (q)",
    "vdc": "tension continue",
    "v_dc": "tension du noeud continu",
    "vdc_f": "tension continue mesuree, filtree",
    "idc": "courant continu",
    "pm": "puissance active mesuree (filtree)",
    "qm": "puissance reactive mesuree (filtree)",
    "p_e": "puissance active instantanee",
    "q_e": "puissance reactive instantanee",
    "Dw": "ecart de pulsation (statisme du 2e ordre)",
    "w": "pulsation du formeur",
    "w_ref": "pulsation du repere (systeme)",
    "w_pll": "pulsation estimee par la PLL",
    "theta": "angle du formeur",
    "theta_g": "angle du repere de reference",
    "theta_pll": "angle de la PLL",
    "M_d": "integrateur de la boucle externe (d)",
    "M_q": "integrateur de la boucle externe (q)",
    "M_VLd": "integrateur de la boucle de tension (d)",
    "M_VLq": "integrateur de la boucle de tension (q)",
    "M_CLd": "integrateur de la boucle de courant (d)",
    "M_CLq": "integrateur de la boucle de courant (q)",
    "M_pll": "integrateur de la PLL",
    "M_z": "integrateur de la boucle de courant continu",
    "M_W": "integrateur de la boucle d'energie",
    "W": "energie stockee dans les sous-modules",
    "Vr": "amplitude de l'oscillateur virtuel",
    "md": "indice de modulation (d)",
    "mq": "indice de modulation (q)",
    "e_d": "tension fabriquee par l'onduleur (d)",
    "e_q": "tension fabriquee par l'onduleur (q)",
    "v_z": "tension continue des bras",
    "R1_pu": "resistance cote onduleur",
    "L1_pu": "inductance cote onduleur",
    "C_pu": "capacite du filtre",
    "R2_pu": "resistance cote reseau",
    "L2_pu": "inductance cote reseau",
    "Ldc_pu": "inductance continue des bras",
    "Rdc_pu": "resistance continue des bras",
    "Cdc_pu": "capacite continue",
    "Gdc_pu": "conductance de fuite continue",
    "W0": "energie nominale des sous-modules",
    "mp": "statisme P-f",
    "nq": "statisme Q-V",
    "wf_p": "pulsation de coupure du filtre de P",
    "wf_q": "pulsation de coupure du filtre de Q",
    "wc": "pulsation de coupure du statisme du 2e ordre",
    "wff": "pulsation de decouplage",
    "Kffv": "anticipation de tension",
    "Kffi": "anticipation de courant",
    "kS": "rapport des puissances Sn/Sb",
    "wb": "pulsation de base",
    "p_ref": "consigne de puissance active",
    "q_ref": "consigne de puissance reactive",
    "ve_ref": "consigne de tension",
    "vdc_ref": "consigne de tension continue",
    "idc_ref": "consigne de courant continu",
}


def symbole_affiche(nom: str) -> sp.Symbol:
    """Le symbole d'affichage d'un nom de modele : ``KpCL`` -> ``K_pCL``,
    ``R1_pu`` -> ``R_1``, ``p_ref`` -> ``p^*``, ``isd`` -> ``i_sd``."""
    if nom in _NOMS:
        return sp.Symbol(_NOMS[nom])
    m = re.fullmatch(r"(.+)_ref", nom)
    if m:
        base = _NOMS.get(m.group(1), m.group(1))
        if base == m.group(1) and len(base) > 1 and "_" not in base:
            base = base[0] + "_" + base[1:]
        return sp.Symbol(base + "^*")
    m = re.fullmatch(r"([A-Za-z]+?)(\d*)_pu", nom)
    if m:
        lettres, chiffres = m.groups()
        if chiffres:
            return sp.Symbol(f"{lettres}_{chiffres}")
        return sp.Symbol(f"{lettres[0]}_{lettres[1:]}" if len(lettres) > 1 else lettres)
    m = re.fullmatch(r"K([a-z]+[A-Za-z]*)", nom)
    if m:
        return sp.Symbol("K_" + m.group(1))
    m = re.fullmatch(r"M_(\w+)", nom)
    if m:
        return sp.Symbol("M_" + m.group(1))
    return sp.Symbol(nom)


def mathml_symbole(nom: str) -> str:
    """Le symbole d'affichage de ``nom``, en MathML."""
    return (
        "<math>" + sp.mathml(symbole_affiche(nom), printer="presentation") + "</math>"
    )


# --------------------------------------------------------------------- fiche
@dataclass
class Ligne:
    """Une equation de la fiche.

    ``genre`` : ``"def"`` (signal intermediaire, ``nom = expr``), ``"alg"``
    (variable algebrique du modele, ``nom = expr``), ``"etat"``
    (``d nom/dt = expr``) ou ``"sortie"`` (``nom = expr``).
    """

    genre: str
    nom: str
    expr: sp.Expr


@dataclass
class Bloc:
    titre: str
    texte: str = ""
    lignes: list[Ligne] = field(default_factory=list)


class _Redaction:
    """Ce qu'une fiche ecrit : ses blocs, et la definition de chaque signal
    intermediaire, pour la verification."""

    def __init__(self, s):
        self.s = s
        self.blocs: list[Bloc] = []
        self.defs: dict[sp.Symbol, sp.Expr] = {}

    def bloc(self, titre: str, texte: str = "") -> Bloc:
        b = Bloc(titre, texte)
        self.blocs.append(b)
        return b

    def signal(self, b: Bloc, nom: str, expr) -> sp.Symbol:
        x = sp.Symbol(f"§{nom}", real=True)
        self.defs[x] = sp.sympify(expr)
        b.lignes.append(Ligne("def", nom, sp.sympify(expr)))
        return x

    @staticmethod
    def etat(b: Bloc, nom: str, expr) -> None:
        b.lignes.append(Ligne("etat", nom, sp.sympify(expr)))

    @staticmethod
    def alg(b: Bloc, nom: str, expr) -> None:
        b.lignes.append(Ligne("alg", nom, sp.sympify(expr)))

    @staticmethod
    def sortie(b: Bloc, nom: str, expr) -> None:
        b.lignes.append(Ligne("sortie", nom, sp.sympify(expr)))


#: le texte des hypotheses communes aux convertisseurs
LIBRE = (
    "Fonctionnement libre : limitation de courant inactive (mode_lim = 0), "
    "saturation de tension inactive (mode_vsat = 0, GFM), impedance virtuelle "
    "inactive (kRv = 0, GFM) et convertisseur en service (bloque = 0). En "
    "limitation, la consigne de courant est ramenee "
    "dans le disque |i*| <= I_max, direction conservee, et les integrateurs de "
    "la boucle externe la suivent par retro-calcul (constante T_aw) ; en "
    "saturation de tension, la modulation est ramenee dans |m| <= M_max et ce "
    "sont les integrateurs de la boucle de courant qui la suivent."
)


def _sur(e, v):
    """``e / v``, sans le ``/ 1`` d'un bus continu ideal."""
    return e if v == 1 else e / v


def _park(r: _Redaction, angle, titre: str, texte: str):
    s = r.s
    b = r.bloc(titre, texte)
    ct, st = sp.cos(angle - s.theta_g), sp.sin(angle - s.theta_g)
    vgd = r.signal(b, "v_gd", s.vd_t * ct + s.vq_t * st)
    vgq = r.signal(b, "v_gq", -s.vd_t * st + s.vq_t * ct)
    return vgd, vgq, ct, st


def _lcl(r: _Redaction, w, md, mq, vdc, vgd, vgq) -> Bloc:
    s = r.s
    b = r.bloc(
        "Onduleur et filtre LCL",
        "L'onduleur moyen fabrique m.v_dc ; le filtre R1-L1 / C_f / R2-L2 "
        "relie cette tension au bus, dans le repere tournant a la pulsation "
        "du convertisseur (termes de couplage en ω).",
    )
    e_d, e_q = (md, mq) if vdc == 1 else (md * vdc, mq * vdc)
    r.etat(
        b,
        "isd",
        (s.wb / s.L1_pu) * (e_d - s.ved - s.R1_pu * s.isd + w * s.L1_pu * s.isq),
    )
    r.etat(
        b,
        "isq",
        (s.wb / s.L1_pu) * (e_q - s.veq - s.R1_pu * s.isq - w * s.L1_pu * s.isd),
    )
    r.etat(b, "ved", (s.wb / s.C_pu) * (s.isd - s.igd + w * s.C_pu * s.veq))
    r.etat(b, "veq", (s.wb / s.C_pu) * (s.isq - s.igq - w * s.C_pu * s.ved))
    r.etat(
        b,
        "igd",
        (s.wb / s.L2_pu) * (s.ved - vgd - s.R2_pu * s.igd + w * s.L2_pu * s.igq),
    )
    r.etat(
        b,
        "igq",
        (s.wb / s.L2_pu) * (s.veq - vgq - s.R2_pu * s.igq - w * s.L2_pu * s.igd),
    )
    return b


def _raccordement(r: _Redaction, ct, st, port_dc=None, i=None) -> None:
    s = r.s
    i_d, i_q = i if i is not None else (s.igd, s.igq)
    b = r.bloc(
        "Raccordement au reseau",
        "Le courant injecte au bus, ramene du repere du convertisseur au "
        "repere global et de la base convertisseur a la base systeme (k_S).",
    )
    r.sortie(b, "id_t", s.kS * (i_d * ct - i_q * st))
    r.sortie(b, "iq_t", s.kS * (i_d * st + i_q * ct))
    if port_dc is not None:
        r.sortie(b, "i_dc", port_dc)


# ----------------------------------------------------------------------- GFL
def _gfl(r: _Redaction, variante: str) -> None:
    s = r.s
    vgd, vgq, ct, st = _park(
        r,
        s.theta_pll,
        "Changement de repere (angle de la PLL)",
        "La tension du bus est lue dans le repere dq de la PLL.",
    )
    b = r.bloc("Puissances mesurees", "Aux bornes du condensateur de filtre.")
    pe = r.signal(b, "p_e", s.ved * s.igd + s.veq * s.igq)
    qe = r.signal(b, "q_e", -s.ved * s.igq + s.veq * s.igd)

    b = r.bloc(
        "PLL",
        "Synchronisation : la PLL annule v_eq, la tension du condensateur sur "
        "l'axe q, par un correcteur PI dont la sortie est la pulsation.",
    )
    r.alg(b, "w_pll", s.M_pll + s.Kppll * s.veq + s.wff)
    r.etat(b, "M_pll", s.Kipll * s.veq)
    r.etat(b, "theta_pll", s.wb * (s.w_pll - s.w_ref))

    if variante == "GFL1":
        b = r.bloc(
            "Boucles externes P et Q",
            "Deux correcteurs PI fixent les consignes de courant a partir des "
            "consignes de puissance.",
        )
        isd_r = r.signal(b, "isd_ref", s.Kpp * (s.p_ref - pe) + s.M_d)
        r.etat(b, "M_d", s.Kip * (s.p_ref - pe))
        vdc = 1
    else:
        v = s.vdc if variante == "GFL2" else s.v_dc
        b = r.bloc(
            "Boucles externes v_dc et Q",
            "La boucle de tension continue fixe le courant actif : c'est "
            "elle qui renvoie au reseau ce que le bus continu recoit. Sans "
            "integrateur (K_idc = 0), c'est un statisme de pente K_pdc.",
        )
        isd_r = r.signal(b, "isd_ref", s.Kpdc * (v - s.vdc_ref) + s.M_d)
        r.etat(b, "M_d", s.Kidc * (v - s.vdc_ref))
        vdc = v
    isq_r = r.signal(b, "isq_ref", s.Kpq * (qe - s.q_ref) + s.M_q)
    r.etat(b, "M_q", s.Kiq * (qe - s.q_ref))

    b = r.bloc(
        "Boucle de courant",
        "Correcteur PI, decouplage des axes (ω_ff.L1) et anticipation de la "
        "tension du condensateur (K_ffv) ; la modulation divise par v_dc.",
    )
    r.alg(
        b,
        "md",
        _sur(
            s.KpCL * (isd_r - s.isd)
            + s.Kffv * s.ved
            - s.wff * s.L1_pu * s.isq
            + s.M_CLd,
            vdc,
        ),
    )
    r.alg(
        b,
        "mq",
        _sur(
            s.KpCL * (isq_r - s.isq)
            + s.Kffv * s.veq
            + s.wff * s.L1_pu * s.isd
            + s.M_CLq,
            vdc,
        ),
    )
    r.etat(b, "M_CLd", s.KiCL * (isd_r - s.isd))
    r.etat(b, "M_CLq", s.KiCL * (isq_r - s.isq))

    _lcl(r, s.w_pll, s.md, s.mq, vdc, vgd, vgq)

    port = None
    if variante == "GFL2":
        b = r.bloc(
            "Bus continu",
            "Le condensateur continu recoit le courant de la source, moins "
            "celui qu'appelle l'onduleur (m.i_s) et la fuite G_dc ; la source "
            "suit sa consigne avec la constante T_dc.",
        )
        r.etat(
            b,
            "vdc",
            (s.wb / s.Cdc_pu)
            * (s.idc - s.Gdc_pu * s.vdc - s.md * s.isd - s.mq * s.isq),
        )
        r.etat(b, "idc", (s.idc_ref - s.idc) / s.Tdc)
    elif variante == "GFL2DC":
        port = -s.kS * s.kVdc * (s.Gdc_pu * s.v_dc + s.md * s.isd + s.mq * s.isq)
    _raccordement(r, ct, st, port)


# ----------------------------------------------------------------------- GFM
def _gfm(r: _Redaction, c) -> None:
    s = r.s
    nom = c.__name__
    dc_interne = nom in ("GFM3", "GFM4", "GFMMatching")
    port = nom in ("GFM3DC", "GFM4DC")
    vgd, vgq, ct, st = _park(
        r,
        s.theta,
        "Changement de repere (angle du formeur)",
        "Le formeur impose son propre angle θ : la tension du bus est lue "
        "dans SON repere.",
    )
    b = r.bloc(
        "Mesure des puissances",
        "Puissances aux bornes du condensateur, filtrees au premier ordre.",
    )
    p = r.signal(b, "p_e", s.ved * s.igd + s.veq * s.igq)
    q = r.signal(b, "q_e", -s.ved * s.igq + s.veq * s.igd)
    if nom != "GFMMatching":
        r.etat(b, "pm", s.wf_p * (p - s.pm))
    r.etat(b, "qm", s.wf_q * (q - s.qm))

    if nom == "GFMdVOC":
        b = r.bloc(
            "Synchronisation : oscillateur virtuel (dVOC)",
            "La pulsation et l'amplitude V_r de l'oscillateur suivent les "
            "puissances normalisees par V_r^2.",
        )
        v2 = s.ve_ref**2
        w = r.signal(b, "w", s.w_ref + s.eta * (s.p_ref / v2 - s.pm / s.Vr**2))
        r.etat(
            b,
            "Vr",
            s.wb
            * s.eta
            * ((s.q_ref / v2 - s.qm / s.Vr**2) + s.alpha * (v2 - s.Vr**2) / v2)
            * s.Vr,
        )
    elif nom == "GFMMatching":
        b = r.bloc(
            "Synchronisation : appariement (matching)",
            "La pulsation est l'image de la tension continue : le "
            "condensateur continu tient le role de l'inertie.",
        )
        w = r.signal(b, "w", s.w_ref * s.vdc / s.vdc_ref)
    elif c.second_order_droop:
        b = r.bloc(
            "Synchronisation : statisme P-f du 2e ordre",
            "L'ecart de pulsation suit le statisme m_p avec la coupure ω_c "
            "(inertie virtuelle).",
        )
        w = r.signal(b, "w", s.Dw + s.w_ref)
        r.etat(b, "Dw", s.mp * s.wc * (s.p_ref - s.pm) - s.wc * s.Dw)
    else:
        b = r.bloc(
            "Synchronisation : statisme P-f",
            "La pulsation du formeur baisse quand sa puissance depasse sa consigne.",
        )
        w = r.signal(b, "w", s.mp * (s.p_ref - s.pm) + s.w_ref)
    r.etat(b, "theta", s.wb * (w - s.w_ref))

    b = r.bloc(
        "Boucle de tension",
        "Consigne de tension du condensateur : la f.e.m. virtuelle (statisme "
        "Q-V ; dVOC : l'amplitude V_r) moins la chute dans l'impedance "
        "virtuelle R_v + j.L_v -- dynamique (zv_dyn = 1 : la chute sur L2 a "
        "l'echelle L_v/L2) ou quasi statique (zv_dyn = 0) ; nulle par defaut. "
        "Correcteur PI, anticipation du courant reseau (K_ffi) et decouplage "
        "(ω_ff.C_f).",
    )
    if nom == "GFMdVOC":
        e_r = s.Vr
    else:
        e_r = r.signal(b, "e_ref", s.ve_ref + (s.q_ref - s.qm) * s.nq)
    kqs = (1 - s.zv_dyn) * w * s.Lv_pu
    kdyn = s.zv_dyn * s.Lv_pu / s.L2_pu
    dzd = r.signal(
        b,
        "dv_vd",
        s.Rv_pu * s.igd - kqs * s.igq + kdyn * (s.ved - vgd - s.R2_pu * s.igd),
    )
    dzq = r.signal(
        b,
        "dv_vq",
        s.Rv_pu * s.igq + kqs * s.igd + kdyn * (s.veq - vgq - s.R2_pu * s.igq),
    )
    ved_r = r.signal(b, "ved_ref", e_r - dzd)
    veq_r = r.signal(b, "veq_ref", -dzq)
    isd_r = r.signal(
        b,
        "isd_ref",
        s.KpVL * (ved_r - s.ved) + s.M_VLd + s.Kffi * s.igd - s.wff * s.C_pu * s.veq,
    )
    isq_r = r.signal(
        b,
        "isq_ref",
        s.KpVL * (veq_r - s.veq) + s.M_VLq + s.Kffi * s.igq + s.wff * s.C_pu * s.ved,
    )
    r.etat(b, "M_VLd", s.KiVL * (ved_r - s.ved))
    r.etat(b, "M_VLq", s.KiVL * (veq_r - s.veq))

    vdc = s.vdc if dc_interne else (s.v_dc if port else 1)
    b = r.bloc(
        "Boucle de courant",
        "Correcteur PI, decouplage (ω_ff.L1) et anticipation de la tension "
        "du condensateur (K_ffv) ; la modulation divise par v_dc.",
    )
    md = r.signal(
        b,
        "md",
        _sur(
            s.KpCL * (isd_r - s.isd)
            + s.Kffv * s.ved
            - s.wff * s.L1_pu * s.isq
            + s.M_CLd,
            vdc,
        ),
    )
    mq = r.signal(
        b,
        "mq",
        _sur(
            s.KpCL * (isq_r - s.isq)
            + s.Kffv * s.veq
            + s.wff * s.L1_pu * s.isd
            + s.M_CLq,
            vdc,
        ),
    )
    r.etat(b, "M_CLd", s.KiCL * (isd_r - s.isd))
    r.etat(b, "M_CLq", s.KiCL * (isq_r - s.isq))

    _lcl(r, w, md, mq, vdc, vgd, vgq)

    courant_port = None
    if dc_interne:
        b = r.bloc(
            "Bus continu",
            "La source continue suit la consigne i_dc* (puissance programmee "
            "et regulation proportionnelle de v_dc) avec la constante T_dc ; "
            "le condensateur recoit sa difference avec le courant appele.",
        )
        idc_r = r.signal(
            b, "idc_ref", s.p_ref / s.vdc_ref + s.Kpdc * (s.vdc_ref - s.vdc)
        )
        r.etat(
            b,
            "vdc",
            (s.wb / s.Cdc_pu) * (s.idc - s.Gdc_pu * s.vdc - md * s.isd - mq * s.isq),
        )
        r.etat(b, "idc", (idc_r - s.idc) / s.Tdc)
    elif port:
        courant_port = -s.kS * s.kVdc * (s.Gdc_pu * s.v_dc + md * s.isd + mq * s.isq)
    _raccordement(r, ct, st, courant_port)


# ----------------------------------------------------------------------- MMC
def _mmc(r: _Redaction, c) -> None:
    s = r.s
    vgd, vgq, ct, st = _park(
        r,
        s.theta_pll,
        "Changement de repere (angle de la PLL)",
        "La tension du bus est lue dans le repere dq de la PLL.",
    )
    b = r.bloc("PLL", "La PLL annule v_gq, la tension du bus sur l'axe q.")
    w_pll = r.signal(b, "w_pll", s.M_pll + s.Kppll * vgq + s.wff)
    r.etat(b, "M_pll", s.Kipll * vgq)
    r.etat(b, "theta_pll", s.wb * (w_pll - s.w_ref))

    b = r.bloc(
        "Boucles externes v_dc et Q",
        "Sur la tension continue MESUREE ET FILTREE (T_vdc) ; sans "
        "integrateur, la boucle de tension continue est un statisme.",
    )
    pe = r.signal(b, "p_e", vgd * s.isd + vgq * s.isq)
    qe = r.signal(b, "q_e", vgq * s.isd - vgd * s.isq)
    r.etat(b, "vdc_f", (s.v_dc - s.vdc_f) / s.T_vdc)
    isd_r = r.signal(b, "isd_ref", s.Kpdc * (s.vdc_f - s.vdc_ref) + s.M_d)
    isq_r = r.signal(b, "isq_ref", s.Kpq * (qe - s.q_ref) + s.M_q)
    r.etat(b, "M_d", s.Kidc * (s.vdc_f - s.vdc_ref))
    r.etat(b, "M_q", s.Kiq * (qe - s.q_ref))

    b = r.bloc(
        "Boucle de courant alternatif",
        "Les bras fabriquent la tension e commandee (modulation compensee).",
    )
    e_d = r.signal(
        b,
        "e_d",
        s.KpCL * (isd_r - s.isd) + s.M_CLd + s.Kffv * vgd - s.wff * s.L1_pu * s.isq,
    )
    e_q = r.signal(
        b,
        "e_q",
        s.KpCL * (isq_r - s.isq) + s.M_CLq + s.Kffv * vgq + s.wff * s.L1_pu * s.isd,
    )
    r.etat(b, "M_CLd", s.KiCL * (isd_r - s.isd))
    r.etat(b, "M_CLq", s.KiCL * (isq_r - s.isq))

    b = r.bloc("Cote alternatif", "L1 = Lc + Larm/2, R1, jusqu'au bus.")
    r.etat(
        b,
        "isd",
        (s.wb / s.L1_pu) * (e_d - vgd - s.R1_pu * s.isd + w_pll * s.L1_pu * s.isq),
    )
    r.etat(
        b,
        "isq",
        (s.wb / s.L1_pu) * (e_q - vgq - s.R1_pu * s.isq - w_pll * s.L1_pu * s.isd),
    )

    if c.regule:
        b = r.bloc(
            "Cote continu et energie (regulee)",
            "La boucle d'energie fixe la consigne de courant continu, "
            "anticipee par la puissance mesuree ; la boucle de courant "
            "continu fabrique la tension des bras v_z.",
        )
        idc_r = r.signal(b, "idc_ref", pe + s.KpW * (s.W0 - s.W) + s.M_W)
        v_z = r.signal(b, "v_z", s.v_dc - (s.KpZ * (idc_r - s.idc) + s.M_z))
        r.etat(b, "M_z", s.KiZ * (idc_r - s.idc))
        r.etat(b, "M_W", s.KiW * (s.W0 - s.W))
    else:
        b = r.bloc(
            "Cote continu et energie (libre)",
            "La tension continue des bras est la somme des tensions des "
            "condensateurs des sous-modules : vue du continu, la station est "
            "L_dc en serie avec leur capacite equivalente (2.W_0).",
        )
        v_z = r.signal(b, "v_z", sp.sqrt(s.W / s.W0))
    p_ac = e_d * s.isd + e_q * s.isq
    r.etat(b, "idc", (s.wb / s.Ldc_pu) * (s.v_dc - v_z - s.Rdc_pu * s.idc))
    r.etat(b, "W", v_z * s.idc - p_ac)
    _raccordement(r, ct, st, -s.kS * s.kVdc * s.idc, (s.isd, s.isq))


#: les classes qui ont une fiche ecrite, et sa redaction
FICHES = {
    "GFL1": lambda r, c: _gfl(r, "GFL1"),
    "GFL2": lambda r, c: _gfl(r, "GFL2"),
    "GFL2DC": lambda r, c: _gfl(r, "GFL2DC"),
    **{
        n: _gfm
        for n in (
            "GFM1",
            "GFM2",
            "GFM3",
            "GFM4",
            "GFM3DC",
            "GFM4DC",
            "GFMdVOC",
            "GFMMatching",
        )
    },
    "MMC": _mmc,
    "MMCEnergie": _mmc,
}

#: le schema (electrique et de commande) que l'interface dessine
SCHEMAS = {
    "GFL1": "gfl_pq",
    "GFL2": "gfl_vdc",
    "GFL2DC": "gfl_vdc_port",
    "GFM1": "gfm",
    "GFM2": "gfm2",
    "GFM3": "gfm_dc",
    "GFM4": "gfm2_dc",
    "GFM3DC": "gfm_port",
    "GFM4DC": "gfm2_port",
    "GFMdVOC": "gfm_dvoc",
    "GFMMatching": "gfm_matching",
    "MMC": "mmc",
    "MMCEnergie": "mmc_energie",
}


def _base(c) -> type:
    """La classe COMPLETE d'une variante reduite (meme nom)."""
    while c.__bases__ and c.__bases__[0].__name__ == c.__name__:
        c = c.__bases__[0]
    return c


def _rediger(c) -> _Redaction:
    spec = c.spec()
    s = _symboles(spec)
    r = _Redaction(s)
    redac = FICHES.get(c.__name__)
    if redac is not None:
        redac(r, c)
        return r
    f, g, h = spec.equations(s)
    b = r.bloc("Equations d'etat", "dx/dt = f(x, z, u)")
    for n, e in zip(spec.states, f):
        r.etat(b, n, e)
    if spec.algebs:
        b = r.bloc("Equations algebriques", "0 = g(x, z, u)")
        for n, e in zip(spec.algebs, g):
            b.lignes.append(Ligne("residu", n, sp.sympify(e)))
    b = r.bloc("Sorties", "y = h(x, z, u)")
    for n, e in zip(spec.outputs, h):
        r.sortie(b, n, e)
    return r


def _deplier(e: sp.Expr, defs: dict) -> sp.Expr:
    for _ in range(len(defs) + 1):
        if not (e.free_symbols & defs.keys()):
            return e
        e = e.xreplace(defs)
    return e


def verifier(nom_ou_classe) -> dict[str, sp.Expr]:
    """``{equation: ecart}`` entre la fiche et le modele compile, en
    fonctionnement libre (``mode_lim = 0``, ``bloque = 0``). Une fiche
    exacte rend des ecarts tous nuls."""
    c = _base(_classe(nom_ou_classe))
    spec = c.spec()
    s = _symboles(spec)
    f, g, h = spec.equations(s)
    libre = {
        getattr(s, n): 0
        for n in ("mode_lim", "mode_vsat", "bloque", "kRv")
        if hasattr(s, n)
    }
    r = _rediger(c)
    if c.__name__ not in FICHES:
        return {}
    modele = {("etat", n): e for n, e in zip(spec.states, f)}
    modele |= {("sortie", n): e for n, e in zip(spec.outputs, h)}
    algs = dict(zip(spec.algebs, g))
    ecarts: dict[str, sp.Expr] = {}
    vus = set()
    for b in r.blocs:
        for li in b.lignes:
            if li.genre == "def":
                continue
            fiche_e = _deplier(li.expr, r.defs)
            if li.genre == "alg":
                # le residu du modele s'annule quand z vaut sa definition
                ref = sp.sympify(algs[li.nom]).xreplace({getattr(s, li.nom): fiche_e})
                d = ref.xreplace(libre)
            else:
                vus.add((li.genre, li.nom))
                d = sp.sympify(modele[(li.genre, li.nom)]).xreplace(libre) - fiche_e
            d = sp.simplify(sp.expand(d))
            ecarts[f"{li.genre}:{li.nom}"] = d
    for n in spec.states:
        if ("etat", n) not in vus:
            ecarts[f"etat:{n}"] = sp.nan  # un etat que la fiche oublie
    return ecarts


# ------------------------------------------------------------------- rendus
def _aplatir(e):
    """Les sommes imbriquees ``((a + b) - c)`` d'une expression non evaluee,
    remises a plat ``a + b - c``, sans changer l'ordre des termes."""
    if not getattr(e, "args", ()) or isinstance(e, sp.Symbol):
        return e
    args = [_aplatir(a) for a in e.args]
    if isinstance(e, sp.Add):
        plats = []
        for a in args:
            plats.extend(a.args if isinstance(a, sp.Add) else (a,))
        return sp.Add(*plats, evaluate=False)
    if isinstance(e, sp.Derivative) or isinstance(e, sp.Function) and not args:
        return e
    try:
        return e.func(*args, evaluate=False)
    except TypeError:
        return e.func(*args)


def _affiche(e: sp.Expr, ecrite: bool = True) -> sp.Expr:
    rempl = {}
    for x in e.free_symbols:
        nom = x.name
        rempl[x] = symbole_affiche(nom[1:] if nom.startswith("§") else nom)
    if not ecrite:
        return e.xreplace(rempl)
    with sp.evaluate(False):
        return _aplatir(e.xreplace(rempl))


def _gauche(li: Ligne) -> sp.Expr:
    x = symbole_affiche(li.nom)
    if li.genre == "etat":
        return sp.Derivative(sp.Function(x.name)(sp.Symbol("t")), sp.Symbol("t"))
    return x


def _mathml(gauche, droite, residu: bool = False) -> str:
    g = "<mn>0</mn>" if residu else _mathml_gauche(gauche)
    d = sp.mathml(droite, printer="presentation", order="none")
    return f'<math display="block"><mrow>{g}<mo>=</mo>{d}</mrow></math>'


def _mathml_gauche(x) -> str:
    if isinstance(x, sp.Derivative):
        nom = sp.mathml(sp.Symbol(x.expr.func.__name__), printer="presentation")
        return f"<mfrac><mrow><mi>d</mi>{nom}</mrow><mrow><mi>d</mi><mi>t</mi></mrow></mfrac>"
    return sp.mathml(x, printer="presentation")


def _latex_gauche(x) -> str:
    if isinstance(x, sp.Derivative):
        return rf"\frac{{d\,{sp.latex(sp.Symbol(x.expr.func.__name__))}}}{{dt}}"
    return sp.latex(x)


def fiche(nom_ou_classe) -> dict:
    """La fiche d'un modele : blocs ordonnes d'equations (MathML, LaTeX,
    texte), le schema a dessiner et la legende de ses symboles."""
    c = _base(_classe(nom_ou_classe))
    spec = c.spec()
    # L'affichage garde l'ordre REDIGE des termes (i* - i, et non -i + i*) :
    # les expressions sont construites sans evaluation. La verification, elle,
    # travaille sur la forme evaluee (:func:`verifier`).
    ecrite = c.__name__ in FICHES
    with sp.evaluate(False) if ecrite else contextlib.nullcontext():
        r = _rediger(c)
    blocs = []
    noms: set[str] = set()
    for b in r.blocs:
        lignes = []
        for li in b.lignes:
            droite = _affiche(li.expr, ecrite)
            gauche = _gauche(li)
            residu = li.genre == "residu"
            lignes.append(
                {
                    "genre": li.genre,
                    "nom": li.nom,
                    "mathml": _mathml(gauche, droite, residu),
                    "latex": ("0" if residu else _latex_gauche(gauche))
                    + " = "
                    + sp.latex(droite, order="none"),
                    "texte": sp.sstr(li.expr).replace("§", ""),
                }
            )
            noms.add(li.nom)
            noms |= {x.name.lstrip("§") for x in li.expr.free_symbols}
        blocs.append({"titre": b.titre, "texte": b.texte, "lignes": lignes})
    legende = [
        {
            "nom": n,
            "mathml": mathml_symbole(n),
            "description": DESCRIPTIONS.get(n, ""),
            "genre": (
                "etat"
                if n in spec.states
                else "parametre"
                if n in spec.params
                else "entree"
                if n in spec.inputs
                else "algebrique"
                if n in spec.algebs
                else "sortie"
                if n in spec.outputs
                else "signal"
            ),
        }
        for n in sorted(noms)
    ]
    doc = inspect.getdoc(c) or ""
    resume = doc.split("\n\n", 1)[0].replace("\n", " ")
    resume = re.sub(r":\w+:`~?([^`]+)`", lambda m: m.group(1).split(".")[-1], resume)
    return {
        "classe": c.__name__,
        "resume": resume.replace("``", ""),
        "ecrite": ecrite,
        "hypotheses": LIBRE if ecrite else "",
        "schema": SCHEMAS.get(c.__name__),
        "etats": list(spec.states),
        "blocs": blocs,
        "legende": legende,
    }
