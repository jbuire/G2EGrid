"""Le catalogue des MACHINES : fiches d'alternateurs publiees, pretes a l'emploi.

POURQUOI -- le catalogue des reseaux (:mod:`g2egrid.reseau.catalogue`) decrit
des reseaux entiers ; il manquait de quoi remplacer UNE machine par une autre
sans recopier douze reactances. Chaque :class:`FicheMachine` porte les donnees
STANDARD d'un alternateur (``Xd``, ``X'd``, ``T'd0``..., celles d'une fiche
constructeur), son inertie, sa filiere et sa source, et sait se convertir dans
les parametres de Park du modele (:func:`g2egrid.reseau.reglage.machine_depuis_standard`).

D'OU VIENNENT LES FICHES

- ``smr_*`` : les six turbo-alternateurs de petits reacteurs modulaires (SMR) de
  la these de Boudot, tableau IV.1 (25 a 245,5 MVA, d'apres Anderson et Fouad,
  *Power System Control and Stability*, 2002) ;
- ``ne39_g1`` a ``ne39_g10`` : les dix machines du reseau IEEE 39 noeuds de la
  meme these, tableau A.1 (cinq nucleaires, trois thermiques, une hydraulique,
  une interconnexion), avec leurs stabilisateurs (tableau A.5) ;
- ``kundur_900_a``, ``kundur_900_b`` : les machines de 900 MVA du reseau a quatre
  machines de la these de Mazauric, tableau III-4 (Kundur, exemple 12.6) ;
- ``nucleaire_epri_1500`` : le groupe nucleaire de 1500 MVA, 1800 tr/min du
  rapport EPRI 1013460 (2006, annexe A), le seul de la liste publie AVEC sa
  ligne d'arbre (:data:`g2egrid.elements.ARBRE_NUCLEAIRE_EPRI`).

CE QUE LES THESES NE DONNENT PAS. Aucune des deux ne decrit la ligne d'arbre
de ses machines : elles en donnent l'inertie totale ``H``, pas la repartition
entre turbines et alternateur ni les raideurs. Le champ ``arbre`` d'une fiche
nomme donc la ligne de :data:`g2egrid.elements.LIGNES_ARBRE` qu'on lui PRETE
(mise a l'echelle de son ``H``, voir
:func:`g2egrid.elements.turbo_alternateur_nomme`) : la ligne EPRI pour les
groupes nucleaires a 60 Hz, celle de l'AM600 a 50 Hz. C'est une hypothese ;
seule ``nucleaire_epri_1500`` porte sa propre ligne.

LES REGULATIONS. Le champ ``regulation`` d'une fiche traduit celles de la
source dans le modele de machine du paquet :

- regulateur de vitesse : le modele a un premier ordre (statisme ``mp`` en %,
  constante ``Tg``), la ou les theses emploient l'IEEEG1 (vapeur) et l'IEEEG3
  (hydraulique). On garde le statisme, et ``Tg`` est le PREMIER MOMENT de la
  reponse du regulateur d'origine (son retard moyen) -- :func:`retard_moyen_ieeeg1`,
  :func:`retard_moyen_ieeeg3`. Les etages de resurchauffe et la reponse inverse
  de l'eau ne sont donc pas representes un par un ;
- regulateur de tension : l'avance-retard ``(1 + s.Tc)/(1 + s.Tb)`` d'un ST1
  (reduction du gain transitoire) s'ecrit EXACTEMENT avec la boucle de
  stabilisation du modele, ``Tfd = Tc`` et ``Kfd = (Tb - Tc)/Ka`` ; la boucle
  ``s.Kf/(1 + s.Tf)`` d'un EXST1 est cette boucle meme (``Kfd = Kf``,
  ``Tfd = Tf``) ;
- stabilisateur : gain, filtre passe-haut et deux avances de phase, terme a terme.

Usage ::

    from g2egrid.reseau.machines import fiche_machine
    f = fiche_machine("smr_245")
    reglages = f.reglages(f0=50.0)            # Park + H + regulation
    t.PARAMS_PAR_BUS[37] = reglages           # la machine du bus 37 devient un SMR
"""

from __future__ import annotations

from dataclasses import dataclass, field

from g2egrid.reseau.reglage import machine_depuis_standard

__all__ = [
    "MACHINES",
    "FicheMachine",
    "fiche_machine",
    "lister_machines",
    "retard_moyen_ieeeg1",
    "retard_moyen_ieeeg3",
]

BOUDOT = (
    "these de Boudot sur l'integration de parcs de SMR dans les reseaux electriques"
)
MAZAURIC = (
    "these d'A.-L. Mazauric sur la flexibilite du nucleaire et la stabilite du reseau"
)
EPRI = (
    "EPRI 1013460, Torsional Interaction Between Electrical Network Phenomena "
    "and Turbine-Generator Shafts : Plant Vulnerability, 2006, annexe A"
)


def retard_moyen_ieeeg1(
    T1: float,
    T2: float,
    T3: float,
    etages: tuple[tuple[float, float], ...],
) -> float:
    """Le retard moyen [s] d'un regulateur IEEEG1 : avance-retard
    ``(1 + s.T2)/(1 + s.T1)``, servomoteur ``T3``, puis les etages de la turbine
    ``(K_i, T_i)`` en cascade (sorties haute pression seules, somme des ``K_i``
    egale a 1).

    C'est le premier moment de sa reponse impulsionnelle : la somme des retards
    moins l'avance, chaque etage pese par sa part de puissance ::

        tau = T1 - T2 + T3 + somme_i K_i.(T_4 + ... + T_i)

    Un premier ordre de constante ``tau`` a le meme gain statique et le meme
    retard moyen : c'est la reduction retenue pour le ``Tg`` du modele."""
    cumul, moment, somme = 0.0, 0.0, 0.0
    for K, T in etages:
        cumul += T
        moment += K * cumul
        somme += K
    return T1 - T2 + T3 + moment / somme


def retard_moyen_ieeeg3(
    Tg: float,
    s_perm: float,
    s_trans: float,
    Tr: float,
    Tw: float,
    a11: float,
    a13: float,
    a21: float,
    a23: float,
) -> float:
    """Le retard moyen [s] d'un regulateur hydraulique IEEEG3 : vannage a
    statisme transitoire (``s_trans``, ``Tr``) et permanent (``s_perm``), puis
    conduite forcee ``a23.(1 + (a11 - a13.a21/a23).Tw.s)/(1 + a11.Tw.s)``.

    Au premier ordre en ``s`` : ``(Tg + s_trans.Tr)/s_perm`` pour le vannage, et
    ``a13.a21/a23 . Tw`` pour l'eau (la reponse inverse compte comme un
    retard)."""
    return (Tg + s_trans * Tr) / s_perm + (a13 * a21 / a23) * Tw


#: IEEEG1 des machines 1 a 9 du reseau 39 noeuds (Boudot, tableau A.2) :
#: K = 20 (statisme 5 %), T1 = 0,2, T2 = 1, T3 = 0,6, etages (K, T) =
#: (0,3 ; 0,6), (0,25 ; 0,5), (0,3 ; 0,8), (0,15 ; 1). Retard moyen : 1,26 s.
_REG_VAPEUR_39 = {
    "mp": 5.0,
    "Tg": retard_moyen_ieeeg1(
        0.2, 1.0, 0.6, ((0.3, 0.6), (0.25, 0.5), (0.3, 0.8), (0.15, 1.0))
    ),
}
#: IEEEG3 de la machine 10 (tableau A.3) : statisme permanent 4 %, transitoire
#: 20 % sur 10 s, eau 0,75 s. Retard moyen : 52,4 s.
_REG_HYDRAULIQUE_39 = {
    "mp": 4.0,
    "Tg": retard_moyen_ieeeg3(0.05, 0.04, 0.2, 10.0, 0.75, 0.5, 1.0, 1.5, 1.0),
}
#: AVR ST1 de toutes les machines (tableau A.4) : Ka = 200, Ta = 0,015 s,
#: reduction du gain transitoire Tc/Tb = 1/10 -> Tfd = 1 s, Kfd = 9/200.
_AVR_ST1_39 = {"Ka": 200.0, "Ta": 0.015, "Tfd": 1.0, "Kfd": (10.0 - 1.0) / 200.0}


@dataclass(frozen=True)
class FicheMachine:
    """La fiche d'un alternateur : donnees standard, inertie, origine."""

    nom: str
    titre: str
    #: puissance apparente [VA], tension nominale [V], facteur de puissance
    Sn_VA: float
    Un_V: float
    cos_phi: float
    #: inertie de TOUTE la ligne d'arbre [s], en base ``Sn_VA``
    H: float
    #: reactances [pu, base machine] et constantes de temps a circuit ouvert [s],
    #: aux cles de :data:`g2egrid.reseau.reglage.STANDARD`
    standard: dict
    #: ``"nucleaire"``, ``"thermique"``, ``"hydraulique"``, ``"interconnexion"``
    filiere: str
    source: str
    #: frequence du reseau d'origine [Hz] (``None`` : la source ne la fixe pas)
    f0: float | None = None
    #: ligne d'arbre pretee (cle de ``LIGNES_ARBRE``), ou ``None``
    arbre: dict | str | None = None
    #: reglages du modele equivalents aux regulations de la source (``mp``, ``Tg``,
    #: PSS...) ; vide : la source n'en donne pas
    regulation: dict = field(default_factory=dict)
    remarque: str = ""

    def park(self, f0: float | None = None) -> dict:
        """Les inductances et resistances du modele de Park, a la frequence
        ``f0`` (celle de la fiche par defaut ; les resistances de rotor en
        dependent par ``wb``)."""
        f = self.f0 if f0 is None else f0
        if f is None:
            raise ValueError(f"{self.nom} : preciser f0 (la source ne la fixe pas)")
        return machine_depuis_standard(self.standard, f0=float(f))

    def nom_arbre(self, f0: float | None = None) -> str | None:
        """La ligne d'arbre pretee a la frequence ``f0``."""
        if isinstance(self.arbre, dict):
            f = self.f0 if f0 is None else f0
            return self.arbre.get(round(f)) if f is not None else None
        return self.arbre

    def reglages(
        self, f0: float | None = None, arbre: bool = False, regulation: bool = True
    ) -> dict:
        """Les reglages a poser dans ``SM_PARAMS`` ou ``PARAMS_PAR_BUS`` : Park,
        ``H``, et les regulations de la source si ``regulation``. ``arbre`` :
        ajoute la cle ``arbre`` (ligne multi-masses pretee, voir le module)."""
        reg = self.park(f0) | {"H": self.H}
        if regulation:
            reg |= self.regulation
        if arbre:
            ligne = self.nom_arbre(f0)
            if ligne is None:
                raise ValueError(f"{self.nom} : pas de ligne d'arbre a cette frequence")
            reg["arbre"] = ligne
        return reg


def _std(Xd, Xq, Xd_p, Xq_p, Xd_pp, Xq_pp, Td0_p, Tq0_p, Td0_pp, Tq0_pp, Xl, Ra):
    """Une fiche standard, dans l'ordre des colonnes des tableaux des theses."""
    return {
        "Xd": Xd,
        "Xq": Xq,
        "Xl": Xl,
        "Xd_p": Xd_p,
        "Xd_pp": Xd_pp,
        "Xq_p": Xq_p,
        "Xq_pp": Xq_pp,
        "Td0_p": Td0_p,
        "Td0_pp": Td0_pp,
        "Tq0_p": Tq0_p,
        "Tq0_pp": Tq0_pp,
        "Ra": Ra,
    }


#: ligne d'arbre pretee a un groupe nucleaire, selon la frequence du reseau
_ARBRE_NUCLEAIRE = {60: "nucleaire_epri", 50: "am600"}

# ---------------------------------------------------------------- SMR (Boudot)
#   Sn [MVA], cos phi, Vn [kV], H [s], xd, xq, x'd, x'q, x''d, x''q,
#   T'd0, T'q0, T''d0, T''q0, xl, rs                  (tableau IV.1)
_SMR = {
    "smr_25": (25.0, 0.8, 13.8, 5.016, 1.25, 1.22, 0.232, 0.715, 0.12, 0.12,
               4.75, 1.5, 0.059, 0.21, 0.119, 0.0014),
    "smr_35": (35.29, 0.85, 13.8, 4.389, 1.4, 1.372, 0.231, 0.8015, 0.118, 0.118,
               5.5, 1.25, 0.04, 0.035, 0.117, 0.0),
    "smr_51": (51.2, 0.8, 13.8, 5.078, 1.27, 1.24, 0.209, 0.85, 0.105, 0.116,
               6.6, 1.25, 0.04, 0.035, 0.104, 0.0),
    "smr_75": (75.0, 0.8, 13.8, 6.187, 1.05, 0.98, 0.185, 0.36, 0.13, 0.13,
               6.1, 0.3, 0.038, 0.099, 0.07, 0.0031),
    "smr_100": (100.0, 0.8, 13.8, 4.985, 1.18, 1.05, 0.22, 0.38, 0.145, 0.145,
                5.9, 0.3, 0.042, 0.092, 0.075, 0.0035),
    "smr_245": (245.5, 0.85, 14.4, 4.624, 1.71, 1.63, 0.32, 0.51, 0.21, 0.21,
                7.1, 0.38, 0.038, 0.073, 0.125, 0.0032),
}  # fmt: skip

# ------------------------------------------------- IEEE 39 noeuds (Boudot, A.1)
#   bus, filiere, Sn [MVA], H, xd, xq, x'd, x'q, x''d, x''q, T'd0, T'q0, T''d0,
#   T''q0, xl, rs ; PSS (tableau A.5) : Kpss, Tw, T1, T2, T3, T4
_NE39 = {
    1: (39, "interconnexion", 10000.0, 5.0, 2.0, 1.9, 0.6, 0.8, 0.4, 0.4,
        7.0, 0.7, 0.05, 0.035, 0.3, 0.0, (1.0, 10.0, 5.0, 0.6, 3.0, 0.5)),
    2: (31, "nucleaire", 700.0, 4.329, 2.065, 1.974, 0.4879, 1.19, 0.35, 0.35,
        6.56, 1.5, 0.05, 0.035, 0.245, 0.0, (0.5, 10.0, 5.0, 0.4, 1.0, 0.1)),
    3: (32, "nucleaire", 800.0, 4.475, 1.996, 1.896, 0.4248, 0.7008, 0.36, 0.36,
        5.7, 1.5, 0.05, 0.035, 0.2432, 0.0, (0.5, 10.0, 3.0, 0.2, 2.0, 0.2)),
    4: (33, "thermique", 800.0, 3.575, 2.096, 2.064, 0.3488, 1.328, 0.28, 0.28,
        5.69, 1.5, 0.05, 0.035, 0.236, 0.0, (2.0, 10.0, 1.0, 0.1, 1.0, 0.3)),
    5: (34, "thermique", 300.0, 4.333, 2.01, 1.86, 0.396, 0.498, 0.267, 0.267,
        5.4, 0.44, 0.05, 0.035, 0.162, 0.0, (1.0, 10.0, 1.5, 0.2, 1.0, 0.1)),
    6: (35, "nucleaire", 800.0, 4.35, 2.032, 1.928, 0.4, 0.6512, 0.32, 0.32,
        7.3, 0.4, 0.05, 0.035, 0.1792, 0.0, (4.0, 10.0, 0.5, 0.1, 0.5, 0.05)),
    7: (36, "thermique", 700.0, 3.771, 2.065, 2.044, 0.343, 1.302, 0.308, 0.308,
        5.66, 1.5, 0.05, 0.035, 0.2254, 0.0, (7.5, 10.0, 0.2, 0.02, 0.5, 0.1)),
    8: (37, "nucleaire", 700.0, 3.471, 2.03, 1.96, 0.399, 0.6377, 0.315, 0.315,
        6.7, 0.41, 0.05, 0.035, 0.196, 0.0, (2.0, 10.0, 1.0, 0.2, 0.5, 0.1)),
    9: (38, "nucleaire", 1000.0, 3.45, 2.106, 2.05, 0.57, 0.587, 0.45, 0.45,
        4.79, 1.96, 0.05, 0.035, 0.298, 0.0, (2.0, 10.0, 1.0, 0.5, 2.0, 0.1)),
    # machine hydraulique a poles saillants : la these ne donne ni x'q ni T'q0
    # (un seul amortisseur d'axe q). Le modele a DEUX enroulements d'axe q : le
    # premier est rendu inoperant (x'q = 0,98.xq, soit une inductance de fuite
    # de 22 pu ; T'q0 = 1 s, sans effet), le second porte x''q et T''q0.
    10: (30, "hydraulique", 1000.0, 4.2, 1.0, 0.69, 0.31, 0.98 * 0.69, 0.25, 0.25,
         10.2, 1.0, 0.05, 0.035, 0.125, 0.0, (1.0, 10.0, 1.0, 0.05, 3.0, 0.5)),
}  # fmt: skip

#: bus du reseau 39 noeuds -> nom de la fiche de sa machine
BUS_NE39 = {v[0]: f"ne39_g{k}" for k, v in _NE39.items()}


def _fiches() -> dict[str, FicheMachine]:
    out: dict[str, FicheMachine] = {}
    for nom, (Sn, cosphi, Vn, H, *x) in _SMR.items():
        out[nom] = FicheMachine(
            nom=nom,
            titre=f"turbo-alternateur de SMR, {Sn:g} MVA",
            Sn_VA=Sn * 1e6,
            Un_V=Vn * 1e3,
            cos_phi=cosphi,
            H=H,
            standard=_std(*x),
            filiere="nucleaire",
            source=BOUDOT + ", tableau IV.1 (d'apres Anderson et Fouad, 2002)",
            arbre=_ARBRE_NUCLEAIRE,
            remarque=(
                "certaines valeurs sont extrapolees par l'auteur (cases grisees du "
                "tableau) ; pour les trois plus petites, x''d depasse xl de 0,001 pu "
                "seulement : l'amortisseur d'axe d est presque sans fuite"
            ),
        )
    for k, (bus, filiere, Sn, H, *x, pss) in _NE39.items():
        Kp, Tw, T1, T2, T3, T4 = pss
        reg = {
            # PSS conventionnel : gain, filtre passe-haut, deux avances de phase
            "K_PSS": Kp,
            "T_HP": Tw,
            "T1n": T1,
            "T1d": T2,
            "T2n": T3,
            "T2d": T4,
        } | _AVR_ST1_39
        reg |= _REG_HYDRAULIQUE_39 if filiere == "hydraulique" else _REG_VAPEUR_39
        out[f"ne39_g{k}"] = FicheMachine(
            nom=f"ne39_g{k}",
            titre=f"IEEE 39 noeuds, machine {k} (bus {bus}), {filiere}, {Sn:g} MVA",
            Sn_VA=Sn * 1e6,
            Un_V=16.5e3,
            cos_phi=0.85,
            H=H,
            standard=_std(*x),
            filiere=filiere,
            source=BOUDOT + ", tableaux A.1 a A.5 (benchmark IEEE PES)",
            f0=60.0,
            arbre="nucleaire_epri" if filiere == "nucleaire" else None,
            regulation=reg,
            remarque=(
                "x'q et T'q0 absents de la these (poles saillants) : premier "
                "amortisseur d'axe q rendu inoperant"
                if filiere == "hydraulique"
                else ""
            ),
        )
    kundur = _std(1.8, 1.7, 0.3, 0.55, 0.25, 0.25, 8.0, 0.4, 0.03, 0.05, 0.2, 0.0025)
    for lettre, H in (("a", 6.5), ("b", 6.175)):
        out[f"kundur_900_{lettre}"] = FicheMachine(
            nom=f"kundur_900_{lettre}",
            titre=f"reseau 4 machines, zone {lettre.upper()}, 900 MVA",
            Sn_VA=900e6,
            Un_V=20e3,
            cos_phi=0.9,
            H=H,
            standard=kundur,
            filiere="nucleaire",
            source=MAZAURIC + ", tableau III-4 (Kundur, exemple 12.6)",
            f0=50.0,
            arbre=_ARBRE_NUCLEAIRE,
            # regulation de vitesse "nucleaire" de la these : statisme 4 %,
            # T4 = 10 s pour tenir la rampe de 5 %Pn/min (tableau III-6)
            regulation={"mp": 4.0, "Tg": 10.0},
            remarque=(
                "la filiere depend du reseau : dans la these, GEN_A2 et GEN_B2 "
                "sont nucleaires, GEN_A1 thermique, GEN_B1 hydraulique"
            ),
        )
    out["nucleaire_epri_1500"] = FicheMachine(
        nom="nucleaire_epri_1500",
        titre="groupe nucleaire reel (anonyme), 1500 MVA, 1800 tr/min",
        Sn_VA=1500e6,
        Un_V=25e3,
        cos_phi=0.9,
        H=3.7495,
        standard=_std(1.99, 1.86, 0.37, 0.564, 0.28, 0.28, 10.0, 0.48, 0.033, 0.055, 0.22, 0.002),
        filiere="nucleaire",
        source=EPRI,
        f0=60.0,
        arbre="nucleaire_epri",
        # excitation statique EXST1 du rapport : Ka = 250, Ta = 0,01 s, Tr = 0,02 s,
        # boucle de stabilisation Kf = 0,05, Tf = 1 s (Tc = Tb : pas d'avance-retard)
        regulation={"Ka": 250.0, "Ta": 0.01, "Tr": 0.02, "Kfd": 0.05, "Tfd": 1.0},
        remarque=(
            "H est la somme des inerties publiees de la ligne d'arbre "
            "(316 584 kg.m2 a 1800 tr/min)"
        ),
    )  # fmt: skip
    return out


#: Toutes les fiches, par nom.
MACHINES: dict[str, FicheMachine] = _fiches()


def fiche_machine(nom: str) -> FicheMachine:
    """La fiche ``nom`` ; leve avec la liste des noms si elle n'existe pas."""
    try:
        return MACHINES[nom]
    except KeyError:
        raise KeyError(
            f"machine {nom!r} inconnue ; fiches : {', '.join(MACHINES)}"
        ) from None


def lister_machines(filiere: str | None = None) -> list[FicheMachine]:
    """Les fiches, eventuellement d'une seule ``filiere``."""
    return [f for f in MACHINES.values() if filiere in (None, f.filiere)]
