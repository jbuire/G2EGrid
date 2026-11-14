"""Comparer deux versions d'un reseau : ce qui a change, et ce que ca change.

LA QUESTION -- on retouche un reseau (un gain, une ligne, un statisme, un
convertisseur de plus) et on veut savoir, d'un coup d'oeil, ce que la
retouche a fait aux modes. Relire deux T_eig cote a cote ne suffit pas : ils
n'ont pas le meme nombre de lignes, et l'ordre (par partie reelle) change des
qu'un mode en double un autre.

LA METHODE -- les deux spectres sont apparies EN ENTIER, par affectation
optimale (algorithme hongrois) sur la distance relative des valeurs propres,
comme :mod:`g2egrid.multi_niveau.impact` apparie l'equivalent et le modele complet. Une
paire dont les deux valeurs propres sont a plus de ``seuil`` (30 %) l'une de
l'autre n'en est pas une : les deux modes sont alors declares l'un disparu,
l'autre apparu. Apparier tout le spectre, et pas seulement la bande regardee,
evite de faire passer pour "apparu" un mode qui n'a fait que glisser a
travers la borne de la bande.

A cote des modes, on rend ce qui a change dans les tables -- lignes ajoutees,
retirees, cellules modifiees, reglages -- et l'ecart de point de
fonctionnement : la retouche a-t-elle deplace le load flow, ou seulement la
dynamique ?
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment

from g2egrid.analyse import modal
from g2egrid.reseau.tables import TABLES, Tables, assembler

__all__ = [
    "Comparaison",
    "DiffTable",
    "ModeApparie",
    "ModeSeul",
    "comparer",
    "differences_tables",
]

#: En dessous, une valeur propre est un mode de jauge (reference d'angle
#: flottante), pas un mode physique.
JAUGE = 1e-6

#: Reglages par famille : des dictionnaires, compares cle par cle.
_REGLAGES = ("SM_PARAMS", "GFM_PARAMS", "GFL_PARAMS", "PARAMS_PAR_BUS")


@dataclass
class ModeApparie:
    lam_a: complex
    lam_b: complex
    dominant_a: str
    dominant_b: str

    @property
    def f_a(self) -> float:
        return float(modal.frequence(np.array([self.lam_a]))[0])

    @property
    def f_b(self) -> float:
        return float(modal.frequence(np.array([self.lam_b]))[0])

    @property
    def zeta_a(self) -> float:
        return float(modal.amortissement(np.array([self.lam_a]))[0])

    @property
    def zeta_b(self) -> float:
        return float(modal.amortissement(np.array([self.lam_b]))[0])

    @property
    def ecart_f(self) -> float:
        """[%] de la frequence de A (0 pour un mode aperiodique)."""
        return 100 * (self.f_b - self.f_a) / self.f_a if self.f_a > 1e-9 else 0.0

    @property
    def ecart_zeta(self) -> float:
        """[points] d'amortissement."""
        return self.zeta_b - self.zeta_a

    @property
    def degrade(self) -> bool:
        """Devenu instable, ou passe sous 5 % d'amortissement."""
        mauvais = lambda lam, z: lam.real > 0 or z < 5.0  # noqa: E731
        return mauvais(self.lam_b, self.zeta_b) and not mauvais(self.lam_a, self.zeta_a)


@dataclass
class ModeSeul:
    """Un mode present dans une seule des deux versions."""

    lam: complex
    dominant: str

    @property
    def f(self) -> float:
        return float(modal.frequence(np.array([self.lam]))[0])

    @property
    def zeta(self) -> float:
        return float(modal.amortissement(np.array([self.lam]))[0])


@dataclass
class DiffTable:
    table: str
    ajoutees: list[tuple] = field(default_factory=list)
    retirees: list[tuple] = field(default_factory=list)
    #: (cle de ligne, colonne, avant, apres)
    modifiees: list[tuple[Any, str, Any, Any]] = field(default_factory=list)

    @property
    def vide(self) -> bool:
        return not (self.ajoutees or self.retirees or self.modifiees)


@dataclass
class Comparaison:
    paires: list[ModeApparie]
    disparus: list[ModeSeul]  # dans A seulement
    apparus: list[ModeSeul]  # dans B seulement
    tables: list[DiffTable]
    #: bus -> (V_a, V_b, theta_a, theta_b), bus communs aux deux versions
    point: dict[int, tuple[float, float, float, float]]
    n_etats: tuple[int, int]
    max_re: tuple[float, float]

    @property
    def ecart_tension(self) -> float:
        return max((abs(v[1] - v[0]) for v in self.point.values()), default=0.0)

    @property
    def degrades(self) -> list[ModeApparie]:
        return [p for p in self.paires if p.degrade]


# ---------------------------------------------------------------- tables
def _egal(a, b) -> bool:
    try:
        return bool(np.isclose(float(a), float(b), rtol=1e-12, atol=0.0))
    except (TypeError, ValueError):
        return a == b


def _par_cle(rows: list[tuple]) -> dict:
    """Les lignes indexees par leur premiere colonne (ID ou BusID) ; par leur
    rang si cette colonne ne les distingue pas."""
    cles = [r[0] for r in rows]
    if len(set(cles)) == len(cles):
        return {r[0]: tuple(r) for r in rows}
    return {("rang", i): tuple(r) for i, r in enumerate(rows)}


def differences_tables(a: Tables, b: Tables) -> list[DiffTable]:
    """Ce qui a change d'une version a l'autre, table par table.

    Les lignes sont reconnues par leur premiere colonne -- ``BusID`` pour
    ``Y_DER``, ``ID`` ailleurs : un equipement dont on change un parametre
    est une ligne MODIFIEE, pas une ligne retiree puis ajoutee. Les
    reglages (``SM_PARAMS``...) sont compares cle par cle ; ``f0`` et ``Sb``
    vont dans une table ``reseau``. Seules les tables qui different sont
    rendues.
    """
    out = []
    for nom, cols in TABLES.items():
        da, db = _par_cle(getattr(a, nom)), _par_cle(getattr(b, nom))
        d = DiffTable(nom)
        d.retirees = [da[k] for k in da if k not in db]
        d.ajoutees = [db[k] for k in db if k not in da]
        for k in da:
            if k in db:
                for j, c in enumerate(cols):
                    va, vb = da[k][j], db[k][j]
                    if not _egal(va, vb):
                        d.modifiees.append((k, c, va, vb))
        if not d.vide:
            out.append(d)
    for nom in _REGLAGES:
        ra, rb = getattr(a, nom) or {}, getattr(b, nom) or {}
        d = DiffTable(nom)
        for k in sorted(set(ra) | set(rb), key=str):
            va, vb = ra.get(k), rb.get(k)
            if isinstance(va, dict) or isinstance(vb, dict):
                va, vb = va or {}, vb or {}
                for c in sorted(set(va) | set(vb), key=str):
                    if not _egal(va.get(c), vb.get(c)):
                        d.modifiees.append((k, str(c), va.get(c), vb.get(c)))
            elif not _egal(va, vb):
                d.modifiees.append((k, str(k), va, vb))
        if not d.vide:
            out.append(d)
    d = DiffTable("reseau")
    for c in ("f0", "Sb"):
        if not _egal(getattr(a, c), getattr(b, c)):
            d.modifiees.append(("reseau", c, getattr(a, c), getattr(b, c)))
    if not d.vide:
        out.append(d)
    return out


# ---------------------------------------------------------------- modes
def _spectre(t: Tables):
    net = assembler(t)
    net.solve_equilibrium()
    ss = net.build_state_space()
    lam, _, p = modal.spectre(ss)
    return net, ss, lam, p


def _physiques(lam: np.ndarray) -> np.ndarray:
    """Les modes physiques, une fois chacun : les paires oscillantes par leur
    membre de frequence positive, les modes reels tels quels."""
    return np.where((np.abs(lam) > JAUGE) & (lam.imag >= -1e-9))[0]


def comparer(a: Tables, b: Tables, seuil: float = 0.3) -> Comparaison:
    """Compare deux versions d'un reseau (voir l'en-tete du module)."""
    net_a, ss_a, lam_a, p_a = _spectre(a)
    net_b, ss_b, lam_b, p_b = _spectre(b)
    ka, kb = _physiques(lam_a), _physiques(lam_b)
    la, lb = lam_a[ka], lam_b[kb]
    # Un mode reel ne s'apparie pas a une paire oscillante : on le pousse loin.
    cout = np.abs(la[:, None] - lb[None, :]) / np.maximum(np.abs(la[:, None]), 1e-9)
    genre = (np.abs(la.imag) > 1e-6)[:, None] != (np.abs(lb.imag) > 1e-6)[None, :]
    cout = np.where(genre, 1e6, cout)
    i, j = linear_sum_assignment(cout)
    ok = cout[i, j] < seuil
    i, j = i[ok], j[ok]

    noms_a, noms_b = list(ss_a.states), list(ss_b.states)
    dom = lambda p, k, noms: noms[int(np.argmax(p[k]))]  # noqa: E731
    paires = [
        ModeApparie(
            complex(la[x]),
            complex(lb[y]),
            dom(p_a, ka[x], noms_a),
            dom(p_b, kb[y], noms_b),
        )
        for x, y in zip(i, j)
    ]
    paires.sort(key=lambda m: -m.lam_b.real)
    seuls_a = sorted(set(range(len(ka))) - set(i.tolist()))
    seuls_b = sorted(set(range(len(kb))) - set(j.tolist()))
    disparus = [ModeSeul(complex(la[x]), dom(p_a, ka[x], noms_a)) for x in seuls_a]
    apparus = [ModeSeul(complex(lb[y]), dom(p_b, kb[y], noms_b)) for y in seuls_b]

    ia, ib = net_a._bus_index(), net_b._bus_index()
    point = {
        int(bus): (
            float(abs(net_a.lf.V[ia[bus]])),
            float(abs(net_b.lf.V[ib[bus]])),
            float(np.angle(net_a.lf.V[ia[bus]])),
            float(np.angle(net_b.lf.V[ib[bus]])),
        )
        for bus in ia
        if bus in ib
    }
    return Comparaison(
        paires=paires,
        disparus=disparus,
        apparus=apparus,
        tables=differences_tables(a, b),
        point=point,
        n_etats=(len(noms_a), len(noms_b)),
        max_re=(
            float(la.real.max()) if len(la) else 0.0,
            float(lb.real.max()) if len(lb) else 0.0,
        ),
    )
