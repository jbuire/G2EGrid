"""Quels modes bougent le plus quand le reseau change d'etat ?

LA QUESTION -- un mode lu sur un seul point de fonctionnement ne dit pas
s'il est ROBUSTE. Un mode bien amorti au point nominal peut se desamortir
des qu'une unite pousse sa production, que la charge monte, ou qu'une ligne
declenche. On cherche donc les modes les plus SENSIBLES :

  - aux variations de puissance : chaque unite non slack pousse sa
    production de ``frac`` (10 % par defaut), chaque charge monte de
    ``frac``, et toute la charge monte ensemble de ``frac`` ;
  - aux deconnexions : chaque element manoeuvrable (ligne, transformateur,
    unite, charge, shunt) est ouvert a son tour, comme un N-1.

LA METHODE -- chaque variante est reassemblee, son load flow et son
equilibre recalcules, son spectre recalcule. On apparie ensuite chaque mode
du cas de base a un mode de la variante par affectation optimale
(algorithme hongrois) sur un cout qui melange :

  - la RESSEMBLANCE des participations, calculee sur les etats communs aux
    deux modeles (``sum sqrt(p_a . p_b)``, 1 pour deux modes identiques) :
    c'est elle qui suit un mode quand deux modes voisins se croisent, la
    ou la seule distance des valeurs propres se trompe ;
  - la distance relative des valeurs propres.

Les etats communs sont compares par leur NOM : une deconnexion retire des
etats (ceux de l'element ouvert), mais les autres gardent le leur.

Pour chaque mode, on rend l'effet de chaque perturbation (ecarts de
frequence, d'amortissement et de partie reelle) et un classement ; pour le
reseau, le classement des modes par sensibilite. Un mode que la variante ne
retrouve pas (cout trop grand) est dit DISPARU : c'est typiquement le mode
porte par l'element qu'on vient d'ouvrir.

LES LIMITES DES CONVERTISSEURS (:mod:`g2egrid.noyau.limites`) -- une variante qui
pousse un convertisseur au-dela de son courant maximal (ou de sa modulation
maximale) est MARQUEE : ``Etude.hors_limite`` la nomme avec ses
depassements, ses effets portent ``hors_limite``, et le compte rendu la
signale en tete. Elle n'est ni retiree ni sortie des scores. Qu'une
perturbation de 10 % ou un N-1 sorte un convertisseur de ses limites est en
soi un resultat de robustesse, que l'arret ou le retrait cacheraient ; le
spectre de la variante reste celui du modele libre, juste tant que le
materiel suit, et chacun peut ecarter ces effets en le sachant. Les retirer
des scores changerait en silence les classements des etudes deja faites.
Le cas de base est verifie de meme (``Etude.depassements``). Un seul
avertissement :class:`~g2egrid.noyau.limites.DepassementLimite` resume l'etude.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import linear_sum_assignment

from g2egrid.noyau import limites
from g2egrid.analyse import modal
from g2egrid.noyau.progres import signaler
from g2egrid.reseau.disjoncteurs import (
    OUVRABLES,
    UNITES,
    basculer,
    en_service,
    nom_element,
)
from g2egrid.noyau.loadflow import SLACK
from g2egrid.reseau.tables import Tables, assembler

__all__ = [
    "Effet",
    "Etude",
    "ModeSensible",
    "Perturbation",
    "etudier",
    "perturbations",
    "resume",
]

#: En dessous, une valeur propre est un mode de jauge (reference flottante).
JAUGE = 1e-6

#: Ce que leve la chaine load flow -> equilibre -> modele quand une variante n'a
#: PAS de point de fonctionnement : Newton ou equilibre non converge
#: (``RuntimeError``), jacobien singulier (``LinAlgError``, sous ``ValueError``),
#: reseau incoherent apres une deconnexion (``ValueError``), division par zero
#: (``ArithmeticError``). Tout le reste (``TypeError``, ``KeyError``,
#: ``AttributeError``...) est une erreur de programmation et doit remonter.
SANS_POINT_DE_FONCTIONNEMENT = (RuntimeError, ValueError, ArithmeticError)
#: Au-dela, l'appariement est refuse : le mode a disparu de la variante.
COUT_MAX = 1.2
#: Part minimale d'un mode sur les etats que la variante garde : en dessous,
#: le mode etait porte par l'element ouvert, et il disparait avec lui.
PART_MIN = 0.5


@dataclass
class Perturbation:
    #: "puissance", "deconnexion", ou "ilotage" : une deconnexion qui coupe
    #: le reseau en plusieurs iles vivantes (une unite restee seule sur son
    #: transformateur ouvert, par exemple). Ses effets sont rendus, mais
    #: n'entrent pas dans les scores : le reseau n'est plus le meme.
    genre: str
    libelle: str
    table: str = ""
    ident: int = 0


@dataclass
class Effet:
    perturbation: Perturbation
    lam: complex | None  # None : le mode a disparu
    df: float = 0.0  # [Hz]
    dzeta: float = 0.0  # [points de %]
    dsigma: float = 0.0  # [1/s], > 0 : moins stable
    score: float = 0.0  # |d lam| / |lam0|, ou 1 pour un mode disparu
    #: la variante a un convertisseur au-dela de sa limite (voir
    #: ``Etude.hors_limite``) : l'effet est celui du modele libre
    hors_limite: bool = False

    @property
    def disparu(self) -> bool:
        return self.lam is None


@dataclass
class ModeSensible:
    lam: complex
    f: float
    zeta: float
    dominant: str
    elements: list[tuple[str, float]]  # element -> participation cumulee
    effets: list[Effet] = field(default_factory=list)  # du plus fort au plus faible

    def suivis(self, genre: str | None = None) -> list[Effet]:
        """Les effets ou le mode est RETROUVE (ni disparu, ni ilotage)."""
        return [
            e
            for e in self.effets
            if not e.disparu
            and e.perturbation.genre != "ilotage"
            and (genre is None or e.perturbation.genre == genre)
        ]

    def score(self, genre: str | None = None) -> float:
        """Le plus fort deplacement relatif ``|d lam| / |lam|`` du mode.

        Un mode qui DISPARAIT n'y compte pas : c'est presque toujours le mode
        porte par l'element qu'on vient d'ouvrir, ce qui ne dit rien de sa
        sensibilite (voir :attr:`disparitions`).
        """
        return max((e.score for e in self.suivis(genre)), default=0.0)

    @property
    def disparitions(self) -> list[str]:
        return [
            e.perturbation.libelle
            for e in self.effets
            if e.disparu and e.perturbation.genre != "ilotage"
        ]

    @property
    def zeta_min(self) -> float:
        """L'amortissement le plus bas atteint parmi les variantes [%]."""
        return min([self.zeta, *(self.zeta + e.dzeta for e in self.suivis())])

    @property
    def pire(self) -> Effet | None:
        """La perturbation qui desamortit le plus ce mode."""
        reste = self.suivis()
        return min(reste, key=lambda e: e.dzeta) if reste else None


@dataclass
class Etude:
    nom: str
    frac: float
    modes: list[ModeSensible]  # du plus sensible au moins sensible
    #: perturbations dont la variante n'a pas de point de fonctionnement
    echecs: list[tuple[Perturbation, str]] = field(default_factory=list)
    #: perturbations qui rendent le reseau instable : max Re, et les elements
    #: qui portent le mode instable
    instables: list[tuple[Perturbation, float, list[tuple[str, float]]]] = field(
        default_factory=list
    )
    n_variantes: int = 0
    #: les convertisseurs hors limite au point de BASE (:class:`g2egrid.noyau.limites.Depassement`)
    depassements: list[limites.Depassement] = field(default_factory=list)
    #: les variantes dont le point a un convertisseur hors limite, avec ses
    #: depassements
    hors_limite: list[tuple[Perturbation, list[limites.Depassement]]] = field(
        default_factory=list
    )


# ------------------------------------------------------------ perturbations
def _charger(t: Tables, unites: dict, frac: float, y_der, loads) -> Tables:
    """Une copie ou les charges designees montent de ``frac``, et ou la
    production suit : le surplus est reparti sur les unites reglantes (SM et
    GFM a batterie, pas les GFL, qui sont des PV) au prorata de leur
    puissance nominale, comme le ferait le reglage primaire. La part
    du slack reste au slack (le load flow la prend). Sans ce partage, le
    slack porterait tout, et une charge eloignee du slack ferait surtout
    croitre le transit vers elle -- ce qui mesure le reseau, plus la charge.
    """
    t2 = t.copie()
    dp = 0.0
    for i in y_der:
        r = t.Y_DER[i]
        dp += float(r[7]) * frac
        _mettre(t2, "Y_DER", i, 7, float(r[7]) * (1 + frac))
        _mettre(t2, "Y_DER", i, 8, float(r[8]) * (1 + frac))
    for i in loads:
        r = t.T_Loads[i]
        dp += float(r[2]) * frac
        _mettre(t2, "T_Loads", i, 2, float(r[2]) * (1 + frac))
        _mettre(t2, "T_Loads", i, 3, float(r[3]) * (1 + frac))
    _repartir(t2, unites, dp)
    return t2


def _production_slack(t: Tables, bus: int) -> float:
    """La production [W] du slack au point de base, lue sur le load flow."""
    net = assembler(t)
    lf = net.run_loadflow()
    i = net._bus_index()[bus]
    pl = next(float(r[7]) for r in t.Y_DER if int(r[0]) == bus)
    return float(lf.S[i].real) * net.Sb + pl


def _repartir(t: Tables, unites: dict, dp: float, sauf: int | None = None) -> None:
    """Repartit ``dp`` [W] sur les unites REGLANTES (sauf celle du bus
    ``sauf``), au prorata de leur puissance nominale ; la part du slack reste
    au slack.

    Les GFL n'y prennent pas part : par convention un GFL est un PV, qui
    produit a sa puissance MPPT et n'a pas de reserve. Ce sont les machines
    synchrones et les GFM -- a batterie -- qui reprennent un ecart.
    """
    reglantes = {
        b: u for b, u in unites.items() if b != sauf and not u[0].startswith("T_GFL")
    }
    sn_total = sum(u[2] for u in reglantes.values())
    if not sn_total:
        return
    for i, r in enumerate(t.Y_DER):
        b = int(r[0])
        if b in reglantes and int(r[1]) != SLACK:
            _mettre(t, "Y_DER", i, 5, float(r[5]) + dp * reglantes[b][2] / sn_total)


def _mettre(t: Tables, table: str, i: int, col: int, valeur: float) -> None:
    rows = getattr(t, table)
    r = list(rows[i])
    r[col] = valeur
    rows[i] = tuple(r) if isinstance(rows[i], tuple) else r


def perturbations(
    t: Tables, frac: float = 0.1, puissance: bool = True, deconnexion: bool = True
) -> list[tuple[Perturbation, Tables]]:
    """Les variantes du reseau ``t``, chacune avec la perturbation qui la fait.

    Toute puissance qui change -- charge qui monte, unite ou charge qui
    declenche -- est reprise par les unites reglantes (SM et GFM a batterie ;
    les GFL sont des PV, a leur puissance MPPT) au prorata de leur puissance
    nominale (voir :func:`_repartir`) ; seule la hausse de
    production d'une unite est laissee au slack, puisque c'est elle qu'on
    etudie. La perte du slack lui-meme est reprise par la nouvelle
    reference que designe :func:`g2egrid.reseau.disjoncteurs.en_service`.
    """
    out: list[tuple[Perturbation, Tables]] = []
    unites = {
        int(r[1]): (table, int(r[0]), float(r[2]))
        for table in UNITES
        for r in getattr(t, table)
    }
    y_der = {int(r[0]): r for r in t.Y_DER}
    if puissance:
        for i, r in enumerate(t.Y_DER):
            b = int(r[0])
            if int(r[1]) == SLACK or not float(r[5]) or b not in unites:
                continue
            t2 = t.copie()
            _mettre(t2, "Y_DER", i, 5, float(r[5]) * (1 + frac))
            nom = nom_element(t, *unites[b][:2])
            out.append((Perturbation("puissance", f"{nom} P +{frac:.0%}"), t2))
        for i, r in enumerate(t.Y_DER):
            if float(r[7]) or float(r[8]):
                t2 = _charger(t, unites, frac, y_der=[i], loads=[])
                libelle = f"charge bus {int(r[0])} +{frac:.0%}"
                out.append((Perturbation("puissance", libelle), t2))
        for i, r in enumerate(t.T_Loads):
            t2 = _charger(t, unites, frac, y_der=[], loads=[i])
            nom = nom_element(t, "T_Loads", int(r[0]))
            out.append((Perturbation("puissance", f"{nom} +{frac:.0%}"), t2))
        tout = _charger(t, unites, frac, range(len(t.Y_DER)), range(len(t.T_Loads)))
        out.append((Perturbation("puissance", f"charge totale +{frac:.0%}"), tout))
    if deconnexion:
        deja = {(a, int(b)) for a, b in t.OUVERTS}
        for table in OUVRABLES:
            for r in getattr(t, table):
                ident = int(r[0])
                if (table, ident) in deja:
                    continue
                if table == "Y_DER" and not (float(r[7]) or float(r[8])):
                    continue
                nom = nom_element(t, table, ident)
                t2 = basculer(t, table, ident, True)
                # la puissance perdue est reprise par les autres unites
                if table in UNITES:
                    bus = int(r[1])
                    ligne = y_der.get(bus)
                    if ligne is not None:
                        pg = float(ligne[5])
                        if int(ligne[1]) == SLACK:
                            pg = _production_slack(t, bus)
                        _repartir(t2, unites, pg, sauf=bus)
                elif table == "Y_DER":
                    _repartir(t2, unites, -float(r[7]))
                elif table == "T_Loads":
                    _repartir(t2, unites, -float(r[2]))
                try:
                    iles = len(en_service(t2)[1].iles)
                except ValueError:
                    iles = 1  # l'echec se lira a l'assemblage
                genre = "ilotage" if iles > 1 else "deconnexion"
                out.append((Perturbation(genre, f"ouverture {nom}", table, ident), t2))
    return out


# ------------------------------------------------------------------ spectre
@dataclass
class _Spectre:
    lam: np.ndarray  # modes physiques, un par paire (Im >= 0)
    p: np.ndarray  # participations de ces modes, sur ``etats``
    etats: list[str]
    #: vecteurs propres a droite de ces modes (colonnes), s'ils sont gardes
    V: np.ndarray | None = None
    #: convertisseurs hors limite au point de fonctionnement
    depassements: list = field(default_factory=list)


def _spectre(t: Tables, vecteurs: bool = False) -> _Spectre:
    net = assembler(t)
    net.solve_equilibrium()
    ss = net.build_state_space()
    # dans la base complete : un courant elimine aux noeuds ideaux garde sa
    # participation, et deux modeles se comparent etat par etat (FR-044)
    lam, V, p, etats = modal.spectre_complet(ss)
    k = np.where((np.abs(lam) > JAUGE) & (lam.imag >= -1e-9))[0]
    return _Spectre(
        lam[k],
        p[k],
        etats,
        V[:, k] if vecteurs else None,
        list(net.depassements),
    )


def _apparier(base: _Spectre, var: _Spectre, garde: np.ndarray) -> dict[int, int]:
    """Mode ``i`` de la base (parmi ``garde``) -> mode ``j`` de la variante."""
    commun = [e for e in base.etats if e in set(var.etats)]
    ia = [base.etats.index(e) for e in commun]
    jb = {e: j for j, e in enumerate(var.etats)}
    ib = [jb[e] for e in commun]
    pa = base.p[garde][:, ia]
    pb = var.p[:, ib]
    ressemblance = np.sqrt(pa) @ np.sqrt(pb).T
    # un mode porte surtout par des etats retires disparait avec eux
    porte = pa.sum(axis=1) >= PART_MIN
    la, lb = base.lam[garde], var.lam
    distance = np.abs(la[:, None] - lb[None, :]) / np.maximum(np.abs(la[:, None]), 1e-9)
    cout = (1 - ressemblance) + np.minimum(distance, 10.0)
    genre = (np.abs(la.imag) > 1e-6)[:, None] != (np.abs(lb.imag) > 1e-6)[None, :]
    cout = np.where(genre, 1e6, cout)
    i, j = linear_sum_assignment(cout)
    return {
        int(garde[x]): int(y)
        for x, y in zip(i, j)
        if cout[x, y] < COUT_MAX and porte[x]
    }


def _elements(p: np.ndarray, etats: list[str], n: int = 3) -> list[tuple[str, float]]:
    somme: dict[str, float] = {}
    for v, e in zip(p, etats):
        cle = e.split(".", 1)[0]
        somme[cle] = somme.get(cle, 0.0) + float(v)
    return sorted(somme.items(), key=lambda kv: -kv[1])[:n]


# -------------------------------------------------------------------- etude
def etudier(
    t: Tables,
    f_min: float = 0.05,
    f_max: float = 5.0,
    zeta_max: float | None = None,
    frac: float = 0.1,
    puissance: bool = True,
    deconnexion: bool = True,
) -> Etude:
    """La sensibilite des modes oscillants de ``[f_min, f_max]`` Hz, et
    d'amortissement au plus ``zeta_max`` % s'il est donne.

    Voir l'en-tete du module. Les modes sont rendus du plus sensible au moins
    sensible (score = plus fort deplacement relatif de la valeur propre) ;
    pour chacun, les effets des perturbations, du plus fort au plus faible.
    Les variantes hors limite des convertisseurs sont marquees
    (``Etude.hors_limite``) et resumees par un seul avertissement.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", limites.DepassementLimite)
        etude = _etudier(t, f_min, f_max, zeta_max, frac, puissance, deconnexion)
    hors = [(None, etude.depassements)] if etude.depassements else []
    hors += etude.hors_limite
    if hors:
        quoi = (
            f"robustesse de {etude.nom} : cas de base et "
            if etude.depassements
            else f"robustesse de {etude.nom} : "
        ) + f"{len(etude.hors_limite)} variante(s) sur {etude.n_variantes} hors limite"
        limites.signaler(hors[0][1], quoi=quoi)
    return etude


def _etudier(t, f_min, f_max, zeta_max, frac, puissance, deconnexion) -> Etude:
    base = _spectre(t)
    f = modal.frequence(base.lam)
    zeta = modal.amortissement(base.lam)
    dedans = (f >= f_min) & (f <= f_max)
    if zeta_max is not None:
        dedans &= zeta <= zeta_max
    garde = np.where(dedans)[0]
    modes = {
        int(k): ModeSensible(
            lam=complex(base.lam[k]),
            f=float(f[k]),
            zeta=float(zeta[k]),
            dominant=base.etats[int(np.argmax(base.p[k]))],
            elements=_elements(base.p[k], base.etats),
        )
        for k in garde
    }
    etude = Etude(nom=t.nom, frac=frac, modes=[], depassements=base.depassements)
    variantes = perturbations(t, frac, puissance, deconnexion)
    etude.n_variantes = len(variantes)
    if not garde.size:
        return etude
    for n, (pert, t2) in enumerate(variantes):
        signaler(n / len(variantes), pert.libelle)
        try:
            var = _spectre(t2)
        except SANS_POINT_DE_FONCTIONNEMENT as exc:
            etude.echecs.append((pert, str(exc).splitlines()[0][:200]))
            continue
        hors = bool(var.depassements)
        if hors:
            etude.hors_limite.append((pert, var.depassements))
        if var.lam.size and var.lam.real.max() > 0:
            k = int(np.argmax(var.lam.real))
            etude.instables.append(
                (pert, float(var.lam[k].real), _elements(var.p[k], var.etats))
            )
        paires = _apparier(base, var, garde)
        for k, m in modes.items():
            j = paires.get(k)
            if j is None:
                m.effets.append(Effet(pert, None, score=1.0, hors_limite=hors))
                continue
            l1 = complex(var.lam[j])
            m.effets.append(
                Effet(
                    pert,
                    l1,
                    df=float(modal.frequence(l1)[0]) - m.f,
                    dzeta=float(modal.amortissement(l1)[0]) - m.zeta,
                    dsigma=l1.real - m.lam.real,
                    score=abs(l1 - m.lam) / abs(m.lam),
                    hors_limite=hors,
                )
            )
    for m in modes.values():
        m.effets.sort(key=lambda e: -e.score)
    etude.modes = sorted(modes.values(), key=lambda m: -m.score())
    return etude


def resume(etude: Etude, n_modes: int = 8, n_effets: int = 3) -> str:
    """Un compte rendu texte de l'etude."""
    lignes = [
        f"{etude.nom} : {len(etude.modes)} modes, {etude.n_variantes} variantes "
        f"(puissance +{etude.frac:.0%}, N-1)"
    ]
    for d in etude.depassements:
        lignes.append(f"  CAS DE BASE HORS LIMITE : {d}")
    for p, deps in etude.hors_limite:
        lignes.append(
            f"  HORS LIMITE apres {p.libelle} : {' ; '.join(str(d) for d in deps)}"
        )
    hors = {p.libelle for p, _ in etude.hors_limite}
    marque = lambda libelle: " [hors limite]" if libelle in hors else ""
    for p, re, qui in etude.instables:
        acteurs = ", ".join(f"{n} {100 * v:.0f}%" for n, v in qui)
        lignes.append(
            f"  INSTABLE apres {p.libelle} (max Re = {re:+.3g} ; {acteurs})"
            f"{marque(p.libelle)}"
        )
    ilots = sorted(
        {
            e.perturbation.libelle
            for m in etude.modes
            for e in m.effets
            if e.perturbation.genre == "ilotage"
        }
    )
    if ilots:
        lignes.append(f"  ilotages (hors score) : {', '.join(ilots)}")
    for p, raison in etude.echecs:
        quoi = " (ilotage)" if p.genre == "ilotage" else ""
        lignes.append(
            f"  sans point de fonctionnement apres {p.libelle}{quoi} : {raison}"
        )
    for m in etude.modes[:n_modes]:
        acteurs = ", ".join(f"{n} {100 * v:.0f}%" for n, v in m.elements)
        lignes.append(
            f"- {m.f:.3f} Hz, zeta {m.zeta:.1f} % [{acteurs}] : "
            f"score P {m.score('puissance'):.3f}, N-1 {m.score('deconnexion'):.3f}, "
            f"zeta min {m.zeta_min:.1f} %"
        )
        if m.disparitions:
            lignes.append(f"    disparait avec : {', '.join(m.disparitions)}")
        for genre in ("puissance", "deconnexion"):
            for e in m.suivis(genre)[:n_effets]:
                lignes.append(
                    f"    {e.perturbation.libelle} : df {e.df:+.3f} Hz, "
                    f"dzeta {e.dzeta:+.2f} pts, score {e.score:.3f}"
                    f"{' [hors limite]' if e.hors_limite else ''}"
                )
    return "\n".join(lignes)
