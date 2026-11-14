"""Hypotheses de simplification du calcul des modes.

Chaque etat ``x`` d'un element peut etre rendu ALGEBRIQUE (``dx/dt = 0`` :
``"statique"``) ou FIGE (``x = x0`` : ``"fige"``), par type d'element
(``"Line"``, ``"SynchronousMachine"``, ``"GFL1"``...) ou par element
(``"SM2"``) -- c'est le champ :attr:`g2egrid.reseau.tables.Tables.SIMPLIFICATIONS`.

Ce module dit ce qu'on PERD a le faire : il apparie les modes du modele
complet a ceux du modele simplifie et separe

- les modes RETROUVES, avec leur deplacement ;
- les modes RETIRES, portes surtout par les etats supprimes : c'est le prix
  attendu de l'hypothese, et il ne coute rien tant que ces modes sont bien
  amortis ;
- les modes PERDUS, portes par les etats gardes mais que le modele simplifie
  n'a plus ;
- les modes NOUVEAUX, que le modele simplifie invente.

Une simplification est SURE pour la stabilite quand elle ne retire ni ne perd
aucun mode instable ou peu amorti, n'en invente aucun, et deplace peu les
modes retrouves.

::

    from g2egrid.reseau.catalogue import charger
    from g2egrid.scenarios.simplification import comparer, preregle

    t = charger("kundur")
    c = comparer(t, preregle(t, "reseau_quasi_statique"))
    print(c.resume())
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from g2egrid.analyse import modal
from g2egrid.scenarios.robustesse import JAUGE, _apparier, _elements, _spectre
from g2egrid.reseau.tables import Tables, assembler

#: au-dela, le mode est tenu pour mal amorti [%]
ZETA_ALERTE = 5.0
#: un mode mal amorti n'alerte qu'en dessous de cette frequence [Hz] ; un
#: mode INSTABLE alerte toujours
F_ALERTE = 100.0

#: les hypotheses classiques, par etat de chaque type d'element. Un
#: preregle ne vaut que pour les types presents dans le reseau.
_RESEAU = {"i_d": "statique", "i_q": "statique"}
_BC = ("M_CLd", "M_CLq", "isd", "isq")
_CONV = (
    "GFL",
    "GFL1",
    "GFL2",
    "GFL2DC",
    "GFM",
    "GFM1",
    "GFM2",
    "GFM3",
    "GFM3DC",
    "GFM4",
    "GFM4DC",
    "MMC",
    "MMCEnergie",
)
PREREGLES: dict[str, tuple[str, dict[str, dict[str, str]]]] = {
    "reseau_quasi_statique": (
        "Réseau en régime sinusoïdal établi (phaseurs) : courants des lignes, "
        "transformateurs, charges et tensions des bus algébriques",
        {
            **{
                c: dict(_RESEAU)
                for c in (
                    "Line",
                    "Transformer",
                    "Branch",
                    "ReactanceShunt",
                    "RLLoad",
                    "GenericLoad",
                    "ZLoad",
                    "PLoad",
                    "ILoad",
                    "IB",
                    "VoltageSource",
                    "Defaut",
                )
            },
            "Bus": {"vd": "statique", "vq": "statique"},
            "SynchronousMachineTR": {"i_gd": "statique", "i_gq": "statique"},
        },
    ),
    "flux_stator": (
        "Flux statoriques des machines synchrones algébriques "
        "(modèle de Park réduit, dit « RMS »)",
        {
            c: {"phi_d": "statique", "phi_q": "statique"}
            for c in ("SynchronousMachine", "SM", "SynchronousMachineTR")
        },
    ),
    "amortisseurs": (
        "Enroulements amortisseurs algébriques (plus de subtransitoire)",
        {
            c: {"phi_1d": "statique", "phi_1q": "statique", "phi_2q": "statique"}
            for c in ("SynchronousMachine", "SM", "SynchronousMachineTR")
        },
    ),
    "excitation": (
        "Mesure de tension et excitatrice des machines instantanées, retour "
        "de stabilisation de l'AVR supprimé (l'amplificateur reste dynamique)",
        {
            c: {"e1": "statique", "efd": "statique", "e3": "fige"}
            for c in ("SynchronousMachine", "SM", "SynchronousMachineTR")
        },
    ),
    "pss": (
        "Stabilisateurs (PSS) instantanés",
        {
            c: {s: "statique" for s in ("v1", "v2", "vpss")}
            for c in ("SynchronousMachine", "SM", "SynchronousMachineTR")
        },
    ),
    "turbine": (
        "Turbine et régulateur de vitesse instantanés",
        {
            c: {"Pm": "statique", "Dw1": "statique"}
            for c in ("SynchronousMachine", "SM", "SynchronousMachineTR")
        },
    ),
    "filtres_convertisseurs": (
        "Filtres LC des convertisseurs algébriques",
        {
            c: {s: "statique" for s in ("isd", "isq", "igd", "igq", "ved", "veq")}
            for c in _CONV
            if not c.startswith("MMC")  # un MMC n'a pas de filtre
        },
    ),
    "boucles_courant": (
        "Boucles de courant des convertisseurs instantanées : le courant de "
        "l'onduleur suit sa consigne (intégrateur et inductance côté "
        "convertisseur algébriques ensemble, sans quoi le système est d'indice 2)",
        {c: dict.fromkeys(_BC, "statique") for c in _CONV},
    ),
    "boucles_tension_gfm": (
        "Boucles de tension ET de courant des GFM instantanées : la tension "
        "du condensateur suit sa consigne",
        {
            c: dict.fromkeys((*_BC, "M_VLd", "M_VLq", "ved", "veq"), "statique")
            for c in _CONV
            if "GFM" in c
        },
    ),
    "pll": (
        "Correcteur PI des PLL des GFL instantané (l'angle reste dynamique)",
        {c: {"M_pll": "statique"} for c in _CONV if "GFL" in c or "MMC" in c},
    ),
    "mesures_puissance": (
        "Filtres de mesure de puissance instantanés (GFM, charges)",
        {
            **{c: {"pm": "statique", "qm": "statique"} for c in _CONV if "GFM" in c},
            **{
                c: {"dp_f": "statique", "dq_f": "statique"}
                for c in ("RLLoad", "GenericLoad")
            },
        },
    ),
    "bus_dc": (
        "Bus continu des convertisseurs parfaitement régulé",
        {
            c: {"idc": "statique"}
            if c.startswith("MMC")
            else {"vdc": "statique", "idc": "statique"}
            for c in _CONV
        },
    ),
}


def classe(el) -> str:
    """Le nom du type COMPLET d'un element, meme deja reduit."""
    return getattr(type(el), "_classe_complete", type(el)).__name__


def inventaire(t: Tables) -> dict[str, dict]:
    """Les types d'elements du reseau, leurs etats et leurs elements :
    ``{type: {"etats": [...], "elements": {nom: {etat: sort}}}}``. Le sort
    est celui des dynamiques de famille, SANS les simplifications."""
    net = assembler(t.copie(SIMPLIFICATIONS={}))
    out: dict[str, dict] = {}
    for el in net.all_elements:
        dyn = type(el).dynamiques()
        if not dyn:
            continue
        d = out.setdefault(classe(el), {"etats": list(dyn), "elements": {}})
        d["elements"][el.name] = dyn
    return dict(sorted(out.items()))


def preregle(
    t: Tables, *noms: str, inv: dict | None = None
) -> dict[str, dict[str, str]]:
    """Les simplifications des preregles ``noms``, restreintes aux types et
    aux etats presents dans ``t`` (``inv`` : son :func:`inventaire`, s'il est
    deja fait -- il assemble le reseau)."""
    inv = inv if inv is not None else inventaire(t)
    out: dict[str, dict[str, str]] = {}
    for nom in noms:
        if nom not in PREREGLES:
            raise ValueError(
                f"preregle {nom!r} inconnu ; attendu : {', '.join(PREREGLES)}"
            )
        for typ, etats in PREREGLES[nom][1].items():
            if typ not in inv:
                continue
            garde = {e: s for e, s in etats.items() if e in inv[typ]["etats"]}
            if garde:
                out.setdefault(typ, {}).update(garde)
    return out


def fusionner(*simplifications: dict) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for s in simplifications:
        for k, v in s.items():
            out.setdefault(k, {}).update(v)
    return out


# --------------------------------------------------------------- comparaison
@dataclass
class ModeSimplifie:
    lam: complex  # modele complet (None seulement pour un mode NOUVEAU)
    lam_reduit: complex | None
    #: part du mode portee par les etats supprimes [0-1]
    part_retiree: float
    elements: list[tuple[str, float]]
    statut: str  # "retrouve", "retire", "perdu", "nouveau"
    #: indice du mode dans le spectre reduit (:attr:`ComparaisonSimplification.reduit`)
    j_reduit: int | None = None

    @staticmethod
    def _fz(lam):
        return (
            (None, None)
            if lam is None
            else (
                float(modal.frequence(lam)[0]),
                float(modal.amortissement(lam)[0]),
            )
        )

    @property
    def f(self):
        return self._fz(self.lam)[0]

    @property
    def zeta(self):
        return self._fz(self.lam)[1]

    @property
    def f_reduit(self):
        return self._fz(self.lam_reduit)[0]

    @property
    def zeta_reduit(self):
        return self._fz(self.lam_reduit)[1]

    @property
    def ecart(self) -> float | None:
        if self.lam is None or self.lam_reduit is None:
            return None
        return float(abs(self.lam_reduit - self.lam) / max(abs(self.lam), 1e-9))

    @property
    def dzeta(self) -> float | None:
        if self.lam is None or self.lam_reduit is None:
            return None
        return self.zeta_reduit - self.zeta

    def critique(
        self, zeta_alerte: float = ZETA_ALERTE, f_alerte: float = F_ALERTE
    ) -> bool:
        """Instable, ou mal amorti sous ``f_alerte``, dans l'un OU l'autre
        modele."""
        for lam in (self.lam, self.lam_reduit):
            if lam is None:
                continue
            f, z = modal.frequence(lam)[0], modal.amortissement(lam)[0]
            if lam.real > 0 or (z < zeta_alerte and f <= f_alerte):
                return True
        return False


@dataclass
class ComparaisonSimplification:
    simplifications: dict[str, dict]
    n_etats: int
    n_etats_reduit: int
    max_re: float
    max_re_reduit: float
    modes: list[ModeSimplifie] = field(default_factory=list)
    #: le spectre du modele reduit (valeurs propres, participations, etats,
    #: et vecteurs propres si :func:`comparer` les a gardes)
    reduit: object | None = None

    def _statut(self, s: str) -> list[ModeSimplifie]:
        return [m for m in self.modes if m.statut == s]

    @property
    def retrouves(self):
        return self._statut("retrouve")

    @property
    def retires(self):
        return self._statut("retire")

    @property
    def perdus(self):
        return self._statut("perdu")

    @property
    def nouveaux(self):
        return self._statut("nouveau")

    @property
    def stable(self) -> bool:
        return self.max_re < 0

    @property
    def stable_reduit(self) -> bool:
        return self.max_re_reduit < 0

    def ecart_max(self, f_max: float = 100.0, zeta_max: float = 95.0) -> float:
        return max(
            (
                m.ecart
                for m in self.retrouves
                if m.lam.imag > 0 and m.f <= f_max and m.zeta <= zeta_max
            ),
            default=0.0,
        )

    def dzeta_max(self, f_max: float = 100.0, zeta_max: float = 95.0) -> float:
        return max(
            (
                abs(m.dzeta)
                for m in self.retrouves
                if m.lam.imag > 0 and m.f <= f_max and m.zeta <= zeta_max
            ),
            default=0.0,
        )

    def alertes(
        self, zeta_alerte: float = ZETA_ALERTE, f_alerte: float = F_ALERTE
    ) -> list[str]:
        """Ce qui fait perdre (ou inventer) une information de stabilite."""
        out = []
        if self.stable != self.stable_reduit:
            out.append(
                "verdict de stabilité changé : "
                f"{'stable' if self.stable else 'instable'} → "
                f"{'stable' if self.stable_reduit else 'instable'}"
            )
        for m in self.retires + self.perdus:
            if m.critique(zeta_alerte, f_alerte):
                out.append(
                    f"mode {m.statut} critique : {m.f:.2f} Hz, ζ = {m.zeta:.1f} % "
                    f"({', '.join(n for n, _ in m.elements[:2])})"
                )
        for m in self.nouveaux:
            if m.critique(zeta_alerte, f_alerte):
                out.append(
                    f"mode nouveau critique : {m.f_reduit:.2f} Hz, "
                    f"ζ = {m.zeta_reduit:.1f} %"
                )
        for m in self.retrouves:
            if m.lam.real * m.lam_reduit.real < 0:
                out.append(
                    f"mode {m.f:.2f} Hz : Re change de signe "
                    f"({m.lam.real:.3g} → {m.lam_reduit.real:.3g})"
                )
        return out

    def sure(
        self, zeta_alerte: float = ZETA_ALERTE, f_alerte: float = F_ALERTE
    ) -> bool:
        return not self.alertes(zeta_alerte, f_alerte)

    def resume(self, n: int = 8) -> str:
        lignes = [
            f"{self.n_etats} → {self.n_etats_reduit} etats ; max Re "
            f"{self.max_re:.4g} → {self.max_re_reduit:.4g}",
            f"retrouves {len(self.retrouves)}, retires {len(self.retires)}, "
            f"perdus {len(self.perdus)}, nouveaux {len(self.nouveaux)} ; "
            f"ecart max {100 * self.ecart_max():.2f} %, "
            f"dzeta max {self.dzeta_max():.2f} pts",
        ]
        lignes += [f"  ! {a}" for a in self.alertes()]
        osc = sorted(
            (m for m in self.retrouves if m.lam.imag > 0),
            key=lambda m: m.zeta,
        )[:n]
        for m in osc:
            lignes.append(
                f"  {m.f:8.3f} Hz {m.zeta:6.2f} % -> {m.f_reduit:8.3f} Hz "
                f"{m.zeta_reduit:6.2f} %   {', '.join(x for x, _ in m.elements)}"
            )
        return "\n".join(lignes)


def comparer(
    t: Tables,
    simplifications: dict[str, dict],
    complet=None,
    part_retiree: float = 0.5,
    vecteurs: bool = False,
) -> ComparaisonSimplification:
    """Les modes de ``t`` sans puis avec ``simplifications`` (qui s'ajoutent
    a celles que ``t`` porte deja). ``vecteurs`` garde les vecteurs propres
    du modele reduit (pour classer ses modes, :mod:`g2egrid.analyse.classement`)."""
    complet = complet if complet is not None else _spectre(t)
    simp = fusionner(t.SIMPLIFICATIONS, simplifications)
    reduit = _spectre(t.copie(SIMPLIFICATIONS=simp), vecteurs=vecteurs)
    garde = np.arange(complet.lam.size)
    paires = _apparier(complet, reduit, garde) if garde.size else {}
    retires = [i for i, e in enumerate(complet.etats) if e not in set(reduit.etats)]
    part_ret = complet.p[:, retires].sum(axis=1)
    c = ComparaisonSimplification(
        simplifications=simp,
        n_etats=len(complet.etats),
        n_etats_reduit=len(reduit.etats),
        max_re=float(complet.lam.real.max(initial=-np.inf)),
        max_re_reduit=float(reduit.lam.real.max(initial=-np.inf)),
        reduit=reduit,
    )
    for i in range(complet.lam.size):
        j = paires.get(i)
        c.modes.append(
            ModeSimplifie(
                lam=complex(complet.lam[i]),
                lam_reduit=None if j is None else complex(reduit.lam[j]),
                part_retiree=float(part_ret[i]),
                elements=_elements(complet.p[i], complet.etats),
                statut=(
                    "retrouve"
                    if j is not None
                    else "retire"
                    if part_ret[i] >= part_retiree
                    else "perdu"
                ),
                j_reduit=j,
            )
        )
    pris = set(paires.values())
    for j in range(reduit.lam.size):
        if j not in pris and abs(reduit.lam[j]) > JAUGE:
            c.modes.append(
                ModeSimplifie(
                    lam=None,
                    lam_reduit=complex(reduit.lam[j]),
                    part_retiree=0.0,
                    elements=_elements(reduit.p[j], reduit.etats),
                    statut="nouveau",
                    j_reduit=j,
                )
            )
    return c
