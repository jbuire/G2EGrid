"""Residus, grammiens et rapport par mode critique.

Les participations (:mod:`g2egrid.analyse.modal`) disent quels ETATS font un mode. Elles
ne disent ni par ou l'exciter, ni ou le voir, ni ce qui en ferait bouger la
valeur propre -- c'est pourtant ce qui sert a placer un stabilisateur, un
amortisseur de convertisseur ou une mesure. Ce module repond a ces questions
depuis le meme modele d'etat ``(A, B, C, D)``.

RESIDUS. La fonction de transfert se decompose sur les modes ::

    H(s) = D + sum_k  R_k / (s - lam_k),     R_k = (C v_k) (w_k B)

``v_k`` vecteur propre a droite, ``w_k`` a gauche, ``w_k v_k = 1``. Le residu
``R_k[i, j]`` est le poids du mode ``k`` dans le transfert de l'entree ``j``
vers la sortie ``i`` : grand, le mode est a la fois commandable par ``j`` et
observable sur ``i``, et c'est la que se boucle un amortisseur (la valeur
propre se deplace, au premier ordre, de ``R_k[i, j] . K`` pour un gain ``K``
de ``i`` vers ``j``). ``w_k B`` (commandabilite modale) et ``C v_k``
(observabilite modale) dependent chacun de la normalisation des vecteurs
propres ; leur produit n'en depend pas. Les deux premiers ne s'interpretent
donc que RELATIVEMENT, entree contre entree ou sortie contre sortie dans un
meme mode : on les rend en % du plus grand.

Les residus portent les unites des entrees et des sorties : comparer un
residu en ``pu/pu`` a un residu en ``rad/pu`` n'a pas de sens. Les rangs se
lisent donc par paire d'unites homogenes.

INDICE DE DOMINANCE. Au voisinage de sa resonance, le mode ``k`` donne un
gain de crete d'environ ``sigma_max(R_k) / |Re lam_k|`` : c'est l'indice
qu'on retient pour classer les modes critiques, il combine le faible
amortissement et la presence aux entrees-sorties.

GRAMMIENS. Le grammien de commandabilite ``W_c`` et celui d'observabilite
``W_o`` resolvent les equations de Lyapunov ::

    A W_c + W_c A^H + B B^H = 0          A^H W_o + W_o A + C^H C = 0

Ils n'existent que pour un ``A`` stable ; or un reseau sans bus infini porte un
mode de jauge nul, et un tirage de campagne peut etre instable. On les calcule
donc dans la base modale, sur la seule partie STABLE du spectre, ou ils ont
une forme close ::

    W_c[i, j] = -(B_m B_m^H)[i, j] / (lam_i + conj(lam_j)),    B_m = W B
    W_o[i, j] = -(C_m^H C_m)[i, j] / (conj(lam_i) + lam_j),    C_m = C V

Les modes ecartes sont rendus a part. Les valeurs singulieres de Hankel,
racines des valeurs propres de ``W_c W_o``, ne dependent pas de la base : elles
mesurent l'energie que chaque direction d'etat transmet des entrees aux
sorties, et leur decroissance dit combien d'etats un modele reduit equilibre
devrait garder.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from g2egrid.analyse.modal import _decomposition, amortissement, forme_modale, frequence
from g2egrid.analyse.modal import sensibilite_top

__all__ = [
    "Grammiens",
    "RapportMode",
    "commandabilite",
    "dominance",
    "grammiens",
    "observabilite",
    "rapport",
    "rapport_mode",
    "residus",
]


def _io(ss, entrees, sorties):
    """Indices et noms des entrees et sorties retenues (toutes par defaut)."""

    def choix(noms, demande):
        if demande is None:
            return np.arange(len(noms)), list(noms)
        idx = [noms.index(d) if isinstance(d, str) else int(d) for d in demande]
        return np.array(idx, dtype=int), [noms[i] for i in idx]

    ie, ne = choix(list(ss.inputs), entrees)
    io, no = choix(list(ss.outputs), sorties)
    return ie, ne, io, no


def _instable(lam, tol=1e-6):
    return lam.real >= -tol * np.maximum(1.0, np.abs(lam))


def residus(ss, entrees=None, sorties=None) -> np.ndarray:
    """Les residus ``R_k = (C v_k)(w_k B)`` de tous les modes, forme
    ``(n_modes, n_sorties, n_entrees)``, dans l'ordre de ``eig(A)``."""
    ie, _, io, _ = _io(ss, entrees, sorties)
    _, Vr, Wl = _decomposition(ss)
    Cv = ss.C[io] @ Vr  # (ny, n)
    wB = Wl @ ss.B[:, ie]  # (n, nu)
    return np.einsum("ik,kj->kij", Cv, wB)


def commandabilite(ss, k: int, entrees=None) -> list[tuple[str, float]]:
    """``(entree, %)`` : ``|w_k B_j|`` en % de la plus grande, decroissant."""
    ie, ne, _, _ = _io(ss, entrees, None)
    _, _, Wl = _decomposition(ss)
    c = np.abs(Wl[k] @ ss.B[:, ie])
    c = 100 * c / max(c.max(initial=0.0), 1e-300)
    return [(ne[j], float(c[j])) for j in np.argsort(-c)]


def observabilite(ss, k: int, sorties=None) -> list[tuple[str, float]]:
    """``(sortie, %)`` : ``|C_i v_k|`` en % de la plus grande, decroissant."""
    _, _, io, no = _io(ss, None, sorties)
    _, Vr, _ = _decomposition(ss)
    o = np.abs(ss.C[io] @ Vr[:, k])
    o = 100 * o / max(o.max(initial=0.0), 1e-300)
    return [(no[i], float(o[i])) for i in np.argsort(-o)]


def dominance(ss, entrees=None, sorties=None) -> np.ndarray:
    """``sigma_max(R_k) / |Re lam_k|`` pour chaque mode ; infini pour un mode
    non amorti (``Re lam_k = 0``)."""
    lam = _decomposition(ss)[0]
    R = residus(ss, entrees, sorties)
    s = np.linalg.svd(R, compute_uv=False).max(axis=1) if R.size else np.zeros(len(lam))
    with np.errstate(divide="ignore"):
        return np.where(np.abs(lam.real) > 0, s / np.abs(lam.real), np.inf)


@dataclass
class Grammiens:
    """Grammiens de la partie stable, en base modale."""

    lam: np.ndarray  # valeurs propres stables retenues
    Wc: np.ndarray  # commandabilite (base modale)
    Wo: np.ndarray  # observabilite (base modale)
    V: np.ndarray  # vecteurs propres a droite correspondants (x = V z)
    exclus: np.ndarray  # valeurs propres ecartees (instables ou nulles)

    @property
    def hankel(self) -> np.ndarray:
        """Valeurs singulieres de Hankel, decroissantes."""
        h = np.linalg.eigvals(self.Wc @ self.Wo).real
        return np.sort(np.sqrt(np.clip(h, 0.0, None)))[::-1]

    def physique(self) -> tuple[np.ndarray, np.ndarray]:
        """``(W_c, W_o)`` dans la base des etats : ``V W_c V^H`` et
        ``V^-H W_o V^-1`` (pseudo-inverse : ``V`` n'a que les colonnes stables).
        Reels a l'arrondi pres."""
        Vp = np.linalg.pinv(self.V)
        Wc = self.V @ self.Wc @ self.V.conj().T
        Wo = Vp.conj().T @ self.Wo @ Vp
        return Wc.real, Wo.real

    def ordre(self, tol: float = 1e-3) -> int:
        """Nombre d'etats a garder pour que la plus grande valeur de Hankel
        ecartee soit sous ``tol`` fois la plus grande."""
        h = self.hankel
        if not h.size or h[0] == 0:
            return 0
        return int(np.sum(h > tol * h[0]))


def grammiens(ss, entrees=None, sorties=None, tol: float = 1e-6) -> Grammiens:
    """Grammiens de commandabilite et d'observabilite de la partie stable du
    spectre (voir le module)."""
    ie, _, io, _ = _io(ss, entrees, sorties)
    lam, Vr, Wl = _decomposition(ss)
    garde = ~_instable(lam, tol)
    ls = lam[garde]
    Bm = Wl[garde] @ ss.B[:, ie]
    Cm = ss.C[io] @ Vr[:, garde]
    Wc = -(Bm @ Bm.conj().T) / (ls[:, None] + ls.conj()[None, :])
    Wo = -(Cm.conj().T @ Cm) / (ls.conj()[:, None] + ls[None, :])
    return Grammiens(lam=ls, Wc=Wc, Wo=Wo, V=Vr[:, garde], exclus=lam[~garde])


# ----------------------------------------------------------------------
# rapport par mode critique
# ----------------------------------------------------------------------
@dataclass
class RapportMode:
    """Tout ce qu'on sait d'un mode, pour decider quoi en faire."""

    indice: int
    valeur: complex
    f: float  # frequence amortie [Hz]
    zeta: float  # amortissement [%]
    dominance: float
    etats: list[tuple[str, float]]  # participations [%]
    elements: list[tuple[str, float]]  # participations cumulees par element [%]
    forme: list[tuple[str, float, float]]  # (etat, module, phase [deg])
    entrees: list[tuple[str, float]]  # commandabilite [% du max]
    sorties: list[tuple[str, float]]  # observabilite [% du max]
    residus: list[tuple[str, str, float, float]]  # (sortie, entree, |R|, phase)
    sensibilites: list[tuple[str, str, float]] = field(default_factory=list)

    @property
    def stable(self) -> bool:
        return self.valeur.real < 0

    def texte(self) -> str:
        """Le rapport en Markdown."""
        v = self.valeur
        lignes = [
            f"### Mode {self.indice} : {v.real:.4g} {'+' if v.imag >= 0 else '-'} "
            f"{abs(v.imag):.4g}j  ({self.f:.3g} Hz, zeta {self.zeta:.2f} %"
            f"{'' if self.stable else ', INSTABLE'})",
            "",
            "Elements : "
            + ", ".join(f"{n} {p:.0f} %" for n, p in self.elements if p >= 1),
            "",
            "Etats : " + ", ".join(f"{n} {p:.0f} %" for n, p in self.etats),
            "",
            "Forme modale : "
            + ", ".join(f"{n} {m:.2f}/{a:+.0f} deg" for n, m, a in self.forme),
            "",
            "Commandable par : "
            + ", ".join(f"{n} {p:.0f} %" for n, p in self.entrees[:4]),
            "",
            "Observable sur : "
            + ", ".join(f"{n} {p:.0f} %" for n, p in self.sorties[:4]),
            "",
            f"Indice de dominance : {self.dominance:.3g}",
        ]
        if self.residus:
            lignes += ["", "| sortie | entree | abs(R) | phase |", "|---|---|---|---|"]
            lignes += [
                f"| {o} | {i} | {r:.3g} | {a:+.0f} |" for o, i, r, a in self.residus
            ]
        if self.sensibilites:
            lignes += [
                "",
                "Coefficients de A les plus influents : "
                + ", ".join(f"d{a}/d{b} ({s:.2g})" for a, b, s in self.sensibilites),
            ]
        return "\n".join(lignes)


def _element(etat: str) -> str:
    return etat.split(".", 1)[0]


def rapport_mode(ss, k: int, entrees=None, sorties=None, n: int = 5) -> RapportMode:
    """Le rapport du mode ``k`` (indice dans ``eig(A)``)."""
    lam, Vr, Wl = _decomposition(ss)
    if not -len(lam) <= k < len(lam):
        raise IndexError(f"mode {k} hors du spectre ({len(lam)} modes)")
    k = k % len(lam)
    p = np.abs(Vr[:, k]) * np.abs(Wl[k])
    p = 100 * p / max(p.sum(), 1e-300)
    noms = list(ss.states)
    par_el: dict[str, float] = {}
    for nom, pk in zip(noms, p):
        par_el[_element(nom)] = par_el.get(_element(nom), 0.0) + float(pk)
    ie, ne, io, no = _io(ss, entrees, sorties)
    R = np.outer(ss.C[io] @ Vr[:, k], Wl[k] @ ss.B[:, ie])
    plats = np.argsort(-np.abs(R), axis=None)[:n]
    ii, jj = np.unravel_index(plats, R.shape)
    sigma = np.linalg.svd(R, compute_uv=False).max() if R.size else 0.0
    re = abs(lam[k].real)
    return RapportMode(
        indice=int(k),
        valeur=complex(lam[k]),
        f=float(frequence(lam[k])[0]),
        zeta=float(amortissement(lam[k])[0]),
        dominance=float(sigma / re) if re > 0 else float("inf"),
        etats=[(noms[i], float(p[i])) for i in np.argsort(-p)[:n]],
        elements=sorted(par_el.items(), key=lambda t: -t[1]),
        forme=forme_modale(ss, k, n),
        entrees=commandabilite(ss, k, entrees),
        sorties=observabilite(ss, k, sorties),
        residus=[
            (no[i], ne[j], float(abs(R[i, j])), float(np.degrees(np.angle(R[i, j]))))
            for i, j in zip(ii, jj)
        ],
        sensibilites=sensibilite_top(ss, k, n=min(n, 4)),
    )


def rapport(
    ss,
    n_modes: int = 3,
    entrees=None,
    sorties=None,
    f_min: float = 0.0,
    f_max: float | None = None,
    critere: str = "amortissement",
) -> list[RapportMode]:
    """Les rapports des ``n_modes`` modes critiques, un par paire conjuguee
    (``Im >= 0``), mode de jauge exclu. ``critere`` : ``"amortissement"``
    (les moins amortis d'abord) ou ``"dominance"`` (les plus forts aux
    entrees-sorties retenues d'abord)."""
    lam = _decomposition(ss)[0]
    f = frequence(lam)
    ok = (lam.imag >= 0) & (np.abs(lam) > 1e-6 * max(1.0, np.abs(lam).max()))
    ok &= f >= f_min
    if f_max is not None:
        ok &= f <= f_max
    idx = np.flatnonzero(ok)
    if critere == "amortissement":
        idx = idx[np.argsort(amortissement(lam[idx]))]
    elif critere == "dominance":
        idx = idx[np.argsort(-dominance(ss, entrees, sorties)[idx])]
    else:
        raise ValueError(f"critere {critere!r} : 'amortissement' ou 'dominance'")
    return [rapport_mode(ss, int(k), entrees, sorties) for k in idx[:n_modes]]
