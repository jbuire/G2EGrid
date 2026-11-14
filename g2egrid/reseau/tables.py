"""Assemblage d'un reseau a partir des tables des notebooks.

Les notebooks de ``reseaux_test/`` decrivent un reseau par des tables
(``Y_DER``, ``T_Lines_SI``, ``T_Transfos``, ``T_SMs``...) puis l'assemblent dans
une cellule. Ce module refait exactement ce que fait cette cellule, mais comme
une fonction : c'est ce qu'il faut des qu'on veut construire plusieurs variantes
d'un meme reseau par programme, comme dans
``etudes/transition_ibr.ipynb``.

Les notebooks gardent leur cellule d'assemblage : elle se lit, elle s'edite, et
c'est leur interet pedagogique. ``test_assembleur_identique_au_notebook`` verifie
que les deux chemins donnent le meme reseau.

LE CONTRAT DE COLONNES VIT ICI, PAS DANS L'INTERFACE. :data:`TABLES` et les
deux codecs :func:`tables_vers_json` / :func:`tables_depuis_json` ont ete
ramenes de ``g2egrid_api.analyse``, ou ils vivaient a cote de l'assembleur qui
les lit. La raison est arithmetique : des l'instant ou les tables continues
existent, deux copies d'un meme contrat se mettent a deriver, et c'est
l'assembleur qui perd -- il lit des colonnes que l'interface n'edite plus.
Une seule definition, la ou on l'utilise.

LE CONTINU : quatre tables de plus, et aucune colonne cachee.

``T_BusDC``     un noeud continu, son type (1 = slack, 2 = P impose) et son
                injection NETTE. ``Pinj_W`` s'appelle ainsi et pas ``Pg``/``Pl``
                pour qu'on ne la confonde pas avec le couple generation/charge
                du cote alternatif, qui, lui, se compense.

``T_LinesDC``   un cable R-L a etat, dont la capacite va aux bus par moities,
                comme le ``B/2`` des lignes alternatives. Deux colonnes
                optionnelles, ``Config`` et ``Circuits``, le decrivent par
                pole (:func:`cable_dc`).

``T_ShuntsDC``  le condensateur propre d'un noeud continu.

``T_VSCs``      quelle machine de ``T_GFM3s``/``T_GFM4s``/``T_GFL2s`` devient un
                convertisseur A PORT CONTINU, sur quel noeud, et dans quel mode.
                La cle est ``BusAC`` : les ``ID`` sont locaux a chaque table,
                le numero de bus alternatif, lui, est unique.

LA CAPACITE D'UN NOEUD CONTINU EST CALCULEE ICI, comme ``Cbus`` du cote
alternatif : demi-capacites des cables, plus les ``Cdc_pu`` des convertisseurs
raccordes, plus les shunts. Le lot 3 laissait l'appelant renseigner ``BusDC.C``
a la main et l'a signale comme une dette ; c'est cette dette qu'on solde. Un
convertisseur a port continu n'a plus de condensateur a lui : le noeud EST son
condensateur, et deux convertisseurs sur un meme noeud y mettent les deux.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field, replace
from typing import Any

from g2egrid.elements import (
    GFL1,
    GFL2,
    GFL2DC,
    GFM1,
    GFM2,
    GFM3,
    GFM3DC,
    GFM4,
    GFM4DC,
    IB,
    Bus,
    BusDC,
    ChargeDC_P,
    ChargeDC_R,
    ConvertisseurDCDC,
    ILoad,
    Line,
    LineDC,
    MMC,
    MMCEnergie,
    PLoad,
    ReactanceShunt,
    RLLoad,
    SynchronousMachine,
    SynchronousMachineTR,
    turbo_alternateur_nomme,
    Transformer,
)
from g2egrid.elements.bus import NOEUDS_DEFAUT, R_CAUSALITE
from g2egrid.elements.gfl import GFL2_DEFAULT
from g2egrid.elements.gfm import GFM_DC_DEFAULT
from g2egrid.elements.mmc import MMC_DEFAULT
from g2egrid.elements.blocs import classe_element
from g2egrid.noyau.loadflow import DC_P, DC_SLACK
from g2egrid.reseau.network import Network

CHARGES_SPECIALES = {"I": ILoad, "P": PLoad}

#: cle des reglages (``SM_PARAMS``, ``GFM_PARAMS``, ``GFL_PARAMS`` ou une
#: entree de ``PARAMS_PAR_BUS``) qui selectionne les dynamiques des equipements
#: : ``{"Dw": "statique", "vdc": "fige", "idc": "fige"}``. Voir
#: :func:`appliquer_dynamiques`.
DYNAMIQUES = "dynamiques"
#: cle des reglages d'une machine synchrone qui lui donne une ligne d'arbre
#: multi-masses : un nom de :data:`g2egrid.elements.LIGNES_ARBRE`
ARBRE = "arbre"
#: les deux cas du convertisseur DC/DC (``T_DCDC.Type``)
TYPES_DCDC = ("batterie", "pv")
#: statisme d'une batterie, en pu de courant par pu de tension et par pu de
#: puissance nominale (voir :class:`~g2egrid.elements.ConvertisseurDCDC`)
KV_BATTERIE = 20.0
#: les reperes de reference admis par ``Tables.REFERENCE``
REFERENCES = ("slack", "coi", "nominale")

#: Les colonnes de chaque table, dans l'ordre. C'est le contrat entre le
#: constructeur graphique et l'assembleur : l'interface edite exactement ces
#: colonnes, et :func:`assembler` les lit dans cet ordre.
TABLES: dict[str, list[str]] = {
    "Y_DER": [
        "BusID",
        "Type",
        "Un_V",
        "V_pu",
        "Theta_rad",
        "Pg_W",
        "Qg_var",
        "Pl_W",
        "Ql_var",
    ],
    "T_Lines_SI": ["ID", "From", "To", "Length_km", "R_ohm_km", "X_ohm_km", "B_S_km"],
    "T_Lines_pu": [
        "ID",
        "From",
        "To",
        "Length_km",
        "R_pu",
        "X_pu",
        "B_pu",
        "Un_V",
        "Sb_VA",
    ],
    "T_Shunts": ["ID", "BusID", "B_pu"],
    "T_Transfos": ["ID", "From", "To", "Sn_VA", "U1_V", "U2_V", "R_pu", "X_pu"],
    "T_IBs": ["ID", "BusID", "Sn_VA", "E_pu", "theta_deg", "R_pu", "X_pu"],
    "T_SMs": ["ID", "BusID", "Sn_VA"],
    "T_GFM1s": ["ID", "BusID", "Sn_VA", "Vdc_V"],
    "T_GFM2s": ["ID", "BusID", "Sn_VA", "Vdc_V"],
    "T_GFM3s": ["ID", "BusID", "Sn_VA", "Vdc_V"],
    "T_GFM4s": ["ID", "BusID", "Sn_VA", "Vdc_V"],
    "T_GFL1s": ["ID", "BusID", "Sn_VA", "Vdc_V"],
    "T_GFL2s": ["ID", "BusID", "Sn_VA", "Vdc_V"],
    "T_Loads": ["ID", "BusID", "P_W", "Q_var", "Type"],
    # --- le reseau continu -------------------------------------------------
    "T_BusDC": ["BusID", "Type", "Un_V", "V_pu", "Pinj_W"],
    "T_LinesDC": ["ID", "From", "To", "Length_km", "R_ohm_km", "L_mH_km", "C_uF_km"],
    "T_ShuntsDC": ["ID", "BusID", "C_pu"],
    "T_VSCs": ["ID", "BusAC", "BusDC", "Mode"],
    # convertisseur DC/DC d'une batterie (statisme de tension) ou d'un PV
    # (puissance constante), et charge continue a puissance ("P") ou
    # resistance ("R") constante
    "T_DCDC": ["ID", "BusDC", "Pn_W", "Type", "P_W"],
    "T_LoadsDC": ["ID", "BusDC", "P_W", "Type"],
    # element decrit par un SCHEMA-BLOC (voir :mod:`g2egrid.elements.blocs`) : ``Modele``
    # est une cle de ``Tables.MODELES_BLOCS`` ; ``Sn_VA`` donne le parametre
    # reserve ``kS = Sn/Sb`` (une ligne sans elle : ``kS = 1``)
    "T_Blocs": ["ID", "BusID", "Modele", "Sn_VA"],
}

#: De la machine monolithique a sa variante A PORT CONTINU. Une entree absente
#: n'est pas un oubli : ``GFM1``, ``GFM2`` et ``GFL1`` ont un bus continu IDEAL,
#: ils n'ont pas d'etat continu a sortir dans le reseau, donc pas de port.
VSC_DC: dict[str, tuple[str, Any]] = {
    "T_GFM3s": ("GFM3DC", GFM3DC),
    "T_GFM4s": ("GFM4DC", GFM4DC),
    "T_GFL2s": ("GFL2DC", GFL2DC),
}

#: Colonnes OPTIONNELLES, en fin de ligne : une ligne qui ne les porte pas
#: garde le sens qu'elle avait. ``T_LinesDC.Config`` dit comment lire ses
#: parametres lineiques, ``Circuits`` combien de circuits identiques sont en
#: parallele (voir :func:`cable_dc`).
COLONNES_OPTIONNELLES: dict[str, list[str]] = {
    "T_LinesDC": ["Config", "Circuits"],
}

#: Les configurations d'un cable continu, et ce qu'elles font des parametres
#: lineiques PAR POLE (r, l, c) pour en tirer le conducteur equivalent
#: unique du modele, a la tension POLE A POLE (celle du noeud continu) :
#: facteurs ``(R et L, C)``.
#:
#: * ``equivalent`` -- les colonnes sont DEJA celles du conducteur equivalent
#:   (le sens historique de la table) ;
#: * ``monopole_sym`` -- monopole symetrique +/-U : la boucle traverse les deux
#:   poles (``2r``, ``2l``) et les deux capacites pole-terre sont en serie
#:   entre les poles (``c/2``) ;
#: * ``bipole`` -- bipole equilibre a retour metallique sans courant : meme
#:   equivalent que le monopole symetrique, la station de chaque extremite
#:   agregeant ses deux poles ;
#: * ``monopole_metallique`` -- monopole asymetrique a retour metallique : la
#:   boucle traverse le pole et le retour (``2r``, ``2l``), la seule capacite
#:   sous tension est celle du pole (``c``) ;
#: * ``monopole_terre`` -- monopole a retour par la terre : le pole seul.
CONFIGS_DC: dict[str, tuple[float, float]] = {
    "equivalent": (1.0, 1.0),
    "monopole_sym": (2.0, 0.5),
    "bipole": (2.0, 0.5),
    "monopole_metallique": (2.0, 1.0),
    "monopole_terre": (1.0, 1.0),
}


def cable_dc(ligne: tuple) -> tuple:
    """La ligne de ``T_LinesDC`` ramenee au conducteur EQUIVALENT :
    ``(ID, From, To, Length_km, R_ohm_km, L_mH_km, C_uF_km)``.

    Les colonnes optionnelles ``Config`` (:data:`CONFIGS_DC`, defaut
    ``equivalent``) et ``Circuits`` (defaut 1) sont appliquees : ``n``
    circuits identiques en parallele divisent R et L par ``n`` et multiplient
    C par ``n``. Une ligne a sept colonnes est rendue telle quelle.
    """
    lid, i, j, km, r, ind, c, *opt = ligne
    config = str(opt[0]).strip().lower() if opt and opt[0] not in (None, "") else ""
    config = config or "equivalent"
    if config not in CONFIGS_DC:
        raise ValueError(
            f"T_LinesDC, cable {lid} : configuration {config!r} inconnue "
            f"({', '.join(CONFIGS_DC)})"
        )
    n = opt[1] if len(opt) > 1 and opt[1] not in (None, "") else 1
    n = int(n)
    if n < 1:
        raise ValueError(f"T_LinesDC, cable {lid} : Circuits = {n} (au moins 1)")
    k_rl, k_c = CONFIGS_DC[config]
    return (
        lid,
        i,
        j,
        float(km),
        float(r) * k_rl / n,
        float(ind) * k_rl / n,
        float(c) * k_c * n,
    )


#: Le modele d'une station a port continu de ``T_GFL2s`` (reglage ``modele``
#: de ``GFL_PARAMS`` ou de ``PARAMS_PAR_BUS``) : ``"vsc"`` (defaut), le VSC a
#: deux niveaux et filtre LCL (:class:`~g2egrid.elements.GFL2DC`) ; ``"mmc"``,
#: le modele moyen du convertisseur modulaire multiniveaux, energie libre
#: (:class:`~g2egrid.elements.MMC`) ; ``"mmc_energie"``, le meme a energie
#: regulee (:class:`~g2egrid.elements.MMCEnergie`). Un MMC n'existe qu'avec
#: son port continu.
MODELES_STATION = ("vsc", "mmc", "mmc_energie")
CLASSES_MMC = {"mmc": MMC, "mmc_energie": MMCEnergie}
#: les reglages qu'un MMC comprend : ceux du VSC qui n'en sont pas (filtre
#: LCL, Gdc...) sont ecartes, pour qu'une famille GFL_PARAMS commune aux deux
#: modeles reste valable
_CLES_MMC = frozenset(MMCEnergie.params) | frozenset(MMC_DEFAULT)


def _modele_station(reg: dict, bus) -> str:
    """Le modele de station du reglage ``reg`` (la cle est retiree)."""
    modele = str(reg.pop("modele", "vsc") or "vsc").lower()
    if modele not in MODELES_STATION:
        raise ValueError(
            f"GFL du bus {int(bus)} : modele {modele!r} inconnu "
            f"({', '.join(MODELES_STATION)})"
        )
    return modele


#: Capacite continue nominale par famille, en pu base CONVERTISSEUR. C'est ce
#: que le convertisseur apporte au noeud auquel il se raccorde.
CDC_DEFAUT: dict[str, float] = {
    "T_GFM3s": GFM_DC_DEFAULT["Cdc_pu"],
    "T_GFM4s": GFM_DC_DEFAULT["Cdc_pu"],
    "T_GFL2s": GFL2_DEFAULT["Cdc_pu"],
}


@dataclass
class Tables:
    """Les tables d'un reseau, dans l'ordre et le format des notebooks.

    Chaque champ est la table du meme nom. ``copie()`` en rend un exemplaire
    independant : c'est la base des variantes d'un meme reseau.
    """

    Y_DER: list[tuple] = field(default_factory=list)
    T_Lines_SI: list[tuple] = field(default_factory=list)
    T_Lines_pu: list[tuple] = field(default_factory=list)
    T_Shunts: list[tuple] = field(default_factory=list)
    T_Transfos: list[tuple] = field(default_factory=list)
    T_IBs: list[tuple] = field(default_factory=list)
    T_SMs: list[tuple] = field(default_factory=list)
    T_GFM1s: list[tuple] = field(default_factory=list)
    T_GFM2s: list[tuple] = field(default_factory=list)
    T_GFM3s: list[tuple] = field(default_factory=list)
    T_GFM4s: list[tuple] = field(default_factory=list)
    T_GFL1s: list[tuple] = field(default_factory=list)
    T_GFL2s: list[tuple] = field(default_factory=list)
    T_Loads: list[tuple] = field(default_factory=list)
    #: le reseau CONTINU. Vides, ils laissent le reseau strictement identique a
    #: ce qu'il etait : ``Network.run_loadflow`` ne change de chemin que si
    #: ``bus_specs_dc`` est non vide.
    T_BusDC: list[tuple] = field(default_factory=list)
    T_LinesDC: list[tuple] = field(default_factory=list)
    T_ShuntsDC: list[tuple] = field(default_factory=list)
    T_VSCs: list[tuple] = field(default_factory=list)
    T_DCDC: list[tuple] = field(default_factory=list)
    T_LoadsDC: list[tuple] = field(default_factory=list)
    T_Blocs: list[tuple] = field(default_factory=list)
    f0: float = 50.0
    Sb: float = 100e6
    SM_PARAMS: dict[str, Any] = field(default_factory=dict)
    GFM_PARAMS: dict[str, Any] = field(default_factory=dict)
    GFL_PARAMS: dict[str, Any] = field(default_factory=dict)
    #: reglages propres a un equipement, par numero de bus. Ils l'emportent sur
    #: ``SM_PARAMS``/``GFM_PARAMS``/``GFL_PARAMS``, qui valent pour toute la
    #: famille : c'est ce qui permet de regler autrement les convertisseurs de
    #: la distribution et ceux du transport.
    PARAMS_PAR_BUS: dict[int, dict[str, Any]] = field(default_factory=dict)
    nom: str = "reseau"
    #: disjoncteurs OUVERTS, en couples ``(table, ID)`` -- voir
    #: :mod:`g2egrid.reseau.disjoncteurs`. Les tables elles-memes restent intactes :
    #: c'est l'assembleur qui retire ce qui est ouvert, et ce qui en devient
    #: hors tension.
    OUVERTS: list[tuple[str, int]] = field(default_factory=list)
    #: le repere de reference du modele dynamique -- voir
    #: :meth:`g2egrid.reseau.network.Network.reference_flottante`. ``"slack"`` : l'angle
    #: de l'element du bus slack (regle de V1) ; ``"coi"`` : un repere flottant
    #: qui suit le centre d'inertie des formeurs ; ``"nominale"`` : un repere
    #: flottant a la vitesse nominale.
    REFERENCE: str = "slack"
    #: ce que devient un bus sans capacite : ``"ideal"`` (defaut, FR-005 : loi
    #: des noeuds exacte, courants en serie reduits algebriquement ; en
    #: temporel, BDF2) ou ``"resistance"`` (resistance de causalite, comme
    #: Node.m : le mode des reseaux compares a MATLAB)
    NOEUDS: str = NOEUDS_DEFAUT
    #: la resistance de causalite du mode ``"resistance"``, pu systeme (les
    #: 100 pu de Node.m par defaut)
    R_NOEUDS: float = R_CAUSALITE
    #: hypotheses de simplification du calcul des modes : par TYPE d'element
    #: (nom de classe : ``"Line"``, ``"SynchronousMachine"``, ``"GFL1"``...)
    #: ou par element (``"L3_5_7"``, ``"SM2"``), ``{etat: "statique" |
    #: "fige" | "dynamique"}``. ``"statique"`` : dx/dt = 0, l'etat devient
    #: algebrique ; ``"fige"`` : x = x0. Le reglage d'un element l'emporte sur
    #: celui de son type, et tous deux s'ajoutent aux dynamiques de famille
    #: (``..._PARAMS["dynamiques"]``). Voir :func:`appliquer_dynamiques`.
    SIMPLIFICATIONS: dict[str, dict[str, Any]] = field(default_factory=dict)
    #: les modeles des elements de ``T_Blocs`` : nom -> schema-bloc (le dict
    #: JSON de :func:`g2egrid.elements.blocs.compiler`)
    MODELES_BLOCS: dict[str, dict[str, Any]] = field(default_factory=dict)

    def copie(self, **modifs) -> "Tables":
        """Copie independante, avec les champs de ``modifs`` remplaces."""
        return replace(deepcopy(self), **modifs)

    def bus(self, bid: int) -> tuple:
        """La ligne de ``Y_DER`` du bus ``bid``."""
        for r in self.Y_DER:
            if int(r[0]) == bid:
                return r
        raise KeyError(f"bus {bid} absent de Y_DER")

    def charge_totale(self) -> tuple[float, float]:
        """Somme des charges de ``Y_DER`` et de ``T_Loads``, en W et var."""
        P = sum(r[7] for r in self.Y_DER) + sum(r[2] for r in self.T_Loads)
        Q = sum(r[8] for r in self.Y_DER) + sum(r[3] for r in self.T_Loads)
        return P, Q


# --------------------------------------------------------------- conversions
def _liste(v) -> list:
    return [list(r) for r in v]


def tables_vers_json(t: Tables) -> dict:
    """Les tables, en JSON pur (listes de listes)."""
    d = {nom: _liste(getattr(t, nom)) for nom in TABLES}
    d.update(
        f0=t.f0,
        Sb=t.Sb,
        nom=t.nom,
        SM_PARAMS=dict(t.SM_PARAMS),
        GFM_PARAMS=dict(t.GFM_PARAMS),
        GFL_PARAMS=dict(t.GFL_PARAMS),
        PARAMS_PAR_BUS={str(k): dict(v) for k, v in t.PARAMS_PAR_BUS.items()},
    )
    if t.OUVERTS:
        d["OUVERTS"] = [[str(a), int(b)] for a, b in t.OUVERTS]
    if t.REFERENCE != "slack":
        d["REFERENCE"] = t.REFERENCE
    if t.NOEUDS != NOEUDS_DEFAUT:
        d["NOEUDS"] = t.NOEUDS
    if t.R_NOEUDS != R_CAUSALITE:
        d["R_NOEUDS"] = float(t.R_NOEUDS)
    if t.SIMPLIFICATIONS:
        d["SIMPLIFICATIONS"] = {k: dict(v) for k, v in t.SIMPLIFICATIONS.items()}
    if t.MODELES_BLOCS:
        d["MODELES_BLOCS"] = deepcopy(t.MODELES_BLOCS)
    return d


def tables_depuis_json(d: dict) -> Tables:
    """L'inverse. Les cles absentes prennent une table vide."""
    t = Tables(
        f0=float(d.get("f0", 50.0)),
        Sb=float(d.get("Sb", 100e6)),
        nom=str(d.get("nom", "reseau")),
        SM_PARAMS=dict(d.get("SM_PARAMS") or {}),
        GFM_PARAMS=dict(d.get("GFM_PARAMS") or {}),
        GFL_PARAMS=dict(d.get("GFL_PARAMS") or {}),
        PARAMS_PAR_BUS={
            int(k): dict(v) for k, v in (d.get("PARAMS_PAR_BUS") or {}).items()
        },
        OUVERTS=[(str(a), int(b)) for a, b in (d.get("OUVERTS") or [])],
        REFERENCE=str(d.get("REFERENCE") or "slack"),
        NOEUDS=str(d.get("NOEUDS") or NOEUDS_DEFAUT),
        R_NOEUDS=float(d.get("R_NOEUDS") or R_CAUSALITE),
        SIMPLIFICATIONS={
            str(k): dict(v) for k, v in (d.get("SIMPLIFICATIONS") or {}).items()
        },
        MODELES_BLOCS=deepcopy(dict(d.get("MODELES_BLOCS") or {})),
    )
    for nom in TABLES:
        setattr(t, nom, [tuple(r) for r in (d.get(nom) or [])])
    return t


# ----------------------------------------------------------------- le continu
def _vsc_par_bus(t: Tables) -> dict[int, tuple]:
    """``T_VSCs`` indexee par son BUS ALTERNATIF, et la table d'ou vient le VSC.

    Rend ``{bus_ac: (id, bus_dc, mode, nom_table, classe)}``. C'est ici que se
    fait la seule resolution de nom du lot : ``T_VSCs`` ne dit pas quelle
    machine poser, elle dit que celle qui est DEJA declaree sur ce bus devient
    un convertisseur a port continu.
    """
    out: dict[int, tuple] = {}
    for ligne in t.T_VSCs:
        vid, bus_ac, bus_dc, mode = ligne[0], int(ligne[1]), int(ligne[2]), ligne[3]
        if bus_ac in out:
            raise ValueError(
                f"T_VSCs, VSC {vid} : le bus alternatif {bus_ac} porte deja un "
                "convertisseur. La cle de cette table est BusAC, pas ID."
            )
        trouves = [
            nom
            for nom in (
                "T_GFM1s",
                "T_GFM2s",
                "T_GFM3s",
                "T_GFM4s",
                "T_GFL1s",
                "T_GFL2s",
            )
            if any(int(r[1]) == bus_ac for r in getattr(t, nom))
        ]
        if not trouves:
            raise ValueError(
                f"T_VSCs, VSC {vid} : aucun convertisseur declare au bus "
                f"alternatif {bus_ac} (cherche dans T_GFM1s..T_GFM4s, "
                "T_GFL1s, T_GFL2s)"
            )
        if len(trouves) > 1:
            raise ValueError(
                f"T_VSCs, VSC {vid} : le bus alternatif {bus_ac} porte "
                f"plusieurs convertisseurs ({', '.join(trouves)})"
            )
        table = trouves[0]
        if table not in VSC_DC:
            raise ValueError(
                f"T_VSCs, VSC {vid} : {table} n'a pas de variante a port "
                "continu. Son bus continu est IDEAL : il n'a aucun etat a "
                "sortir dans le reseau. Utiliser T_GFM3s, T_GFM4s ou T_GFL2s."
            )
        nom_classe, classe = VSC_DC[table]
        modes = getattr(classe, "modes_dc", ())
        if mode not in modes:
            raise ValueError(
                f"T_VSCs, VSC {vid} (bus {bus_ac}) : mode '{mode}' impossible "
                f"pour {nom_classe} ; modes admis : {', '.join(modes)}"
            )
        out[bus_ac] = (vid, bus_dc, mode, table, classe)
    return out


def _capacites_dc(t: Tables, vscs: dict[int, tuple], UN_DC: dict[int, float], wb):
    """La capacite de chaque noeud continu, en pu base SYSTEME.

    Le decalque de ``Cbus`` : demi-capacites des cables (schema en pi), shunts,
    et les condensateurs des convertisseurs raccordes. Un ``Cdc_pu`` est donne
    en pu base CONVERTISSEUR ; le ramener a la base systeme demande
    ``kS = Sn/Sb`` et le rapport des bases de tension continue, d'ou le
    ``kS . kVdc^2``.
    """
    Sb = t.Sb
    C = {b: 0.0 for b in UN_DC}
    for _id, i, j, km, _r, _l, c_km in map(cable_dc, t.T_LinesDC):
        i, j = int(i), int(j)
        Zb = UN_DC[i] ** 2 / Sb
        c_pu = wb * (c_km * 1e-6 * km) * Zb
        C[i] += c_pu / 2.0
        C[j] += c_pu / 2.0
    for _id, b, c_pu in t.T_ShuntsDC:
        C[int(b)] += c_pu
    for bus_ac, (_vid, bus_dc, _mode, table, _cls) in vscs.items():
        ligne = next(r for r in getattr(t, table) if int(r[1]) == bus_ac)
        Sn = float(ligne[2])
        famille = _famille_params(t, table)
        reg = _reglages(t, famille, bus_ac)
        defaut = CDC_DEFAUT[table]
        if table == "T_GFL2s" and _modele_station(dict(reg), bus_ac) in CLASSES_MMC:
            defaut = MMC_DEFAULT["Cdc_pu"]  # pas de condensateur continu
        cdc = float(reg.get("Cdc_pu", defaut))
        C[bus_dc] += cdc * (Sn / Sb)  # kVdc = 1, cf. _verifier_bases_dc
    return C


def _famille_params(t: Tables, table: str) -> dict:
    return t.GFL_PARAMS if table.startswith("T_GFL") else t.GFM_PARAMS


def _reglages(t: Tables, famille: dict, bus: int) -> dict:
    return famille | t.PARAMS_PAR_BUS.get(int(bus), {})


def _verifier_bases_dc(t: Tables, vscs: dict[int, tuple], UN_DC: dict[int, float]):
    """La base de tension continue du convertisseur est celle de son noeud.

    RESERVE DU LOT 3, TRANCHEE ICI. ``kVdc = Un_dc / Vdc_n`` ne pondere que le
    COURANT rendu au noeud, pas la tension que lisent les boucles internes du
    convertisseur : les deux ne sont homogenes que pour ``kVdc = 1``. Plutot
    que de laisser l'assembleur fabriquer silencieusement un modele faux des
    que les deux tensions different, on exige qu'elles coincident et on le dit.
    Le jour ou le cote tension sera pondere lui aussi, cette garde tombe et
    ``kVdc`` se calcule ici.
    """
    for bus_ac, (vid, bus_dc, _mode, table, _cls) in vscs.items():
        ligne = next(r for r in getattr(t, table) if int(r[1]) == bus_ac)
        vdc_n = float(ligne[3])
        un_dc = float(UN_DC[bus_dc])
        if abs(vdc_n - un_dc) > 1e-9 * max(un_dc, 1.0):
            raise ValueError(
                f"T_VSCs, VSC {vid} : la colonne Vdc_V de {table} vaut "
                f"{vdc_n:g} V et Un_V du bus continu {bus_dc} vaut {un_dc:g} V. "
                "Les deux bases de tension continue doivent coincider "
                "(kVdc = 1) : la ponderation du cote tension n'est pas portee."
            )


def _loi_gfm(reg: dict, GFM: type, bus) -> type:
    """La classe d'un GFM selon le reglage ``loi`` (retire de ``reg``) :
    ``"statisme"`` (defaut, la classe de la table), ``"dvoc"`` (bus continu
    ideal : T_GFM1s/T_GFM2s) ou ``"matching"`` (bus continu dynamique :
    T_GFM3s/T_GFM4s). Voir :data:`g2egrid.elements.gfm.LOIS_GFM`."""
    from g2egrid.elements.gfm import LOIS_GFM

    loi = str(reg.pop("loi", "statisme") or "statisme").lower()
    if loi == "statisme":
        return GFM
    classe = LOIS_GFM.get((loi, GFM.has_dc))
    if classe is None:
        possibles = sorted({k for k, dc in LOIS_GFM if dc == GFM.has_dc})
        raise ValueError(
            f"GFM du bus {int(bus)} : loi {loi!r} impossible avec {GFM.__name__} "
            f"(bus continu {'dynamique' if GFM.has_dc else 'ideal'}) ; "
            f"possibles : statisme, {', '.join(possibles)}"
        )
    return classe


def assembler(t: Tables) -> Network:
    """Construit le ``Network`` decrit par les tables.

    Reprend, dans le meme ordre, la cellule d'assemblage des notebooks, plus
    les reglages par bus de ``PARAMS_PAR_BUS`` et, depuis le lot 4, le reseau
    continu. L'ordre compte : les noeuds continus sont declares AVANT les
    convertisseurs, qui s'y raccordent, et apres les bus alternatifs, pour que
    l'ordre des elements d'un reseau sans continu ne bouge pas d'un cran.

    DISJONCTEURS -- un reseau qui declare des disjoncteurs ouverts
    (``t.OUVERTS``) est d'abord ramene a sa partie EN SERVICE par
    :func:`g2egrid.reseau.disjoncteurs.en_service` : elements ouverts retires, iles sans
    formeur hors tension, une reference par ile vivante. Sans disjoncteur
    ouvert, rien ne change.
    """
    if t.OUVERTS:
        from g2egrid.reseau.disjoncteurs import (
            en_service,
        )  # import tardif : il importe ce module

        t, _rapport = en_service(t)
    Sb = t.Sb
    if t.REFERENCE not in REFERENCES:
        raise ValueError(
            f"REFERENCE = {t.REFERENCE!r} : {', '.join(map(repr, REFERENCES))}"
        )

    #: element -> selection de dynamiques, appliquee une fois le reseau pose
    selections: dict[str, dict] = {}

    def reglages(famille: dict, bus: int, nom: str | None = None) -> dict:
        """Reglages de famille, surcharges par ceux du bus. La cle
        ``dynamiques`` n'est pas un parametre : elle est mise de cote pour
        l'element ``nom``, fusionnee cle a cle (le bus complete la famille)."""
        par_bus = t.PARAMS_PAR_BUS.get(int(bus), {})
        reg = famille | par_bus
        dyn = dict(famille.get(DYNAMIQUES) or {}) | dict(par_bus.get(DYNAMIQUES) or {})
        reg.pop(DYNAMIQUES, None)
        if dyn and nom is not None:
            selections[nom] = dyn
        return reg

    # --- bases de tension par bus, lues dans Y_DER ---
    UN = {int(r[0]): r[2] for r in t.Y_DER}
    UN_DC = {int(r[0]): float(r[2]) for r in t.T_BusDC}

    # --- le reseau continu : quel convertisseur devient un VSC, et ou ---
    vscs = _vsc_par_bus(t)
    for bus_ac, (vid, bus_dc, _mode, _table, _cls) in vscs.items():
        if bus_dc not in UN_DC:
            raise ValueError(
                f"T_VSCs, VSC {vid} (bus alternatif {bus_ac}) : le noeud "
                f"continu {bus_dc} n'est pas declare dans T_BusDC"
            )
    _verifier_bases_dc(t, vscs, UN_DC)

    # --- capacite shunt de chaque bus : moitie de la charge des lignes qui y
    #     aboutissent (schema en pi), comme C_half = B_pu/2 dans Network.m. Un
    #     bus qui n'en recoit aucune devient un noeud ALGEBRIQUE a resistance de
    #     causalite R = R_NOEUDS pu (1e5 ; Node.m / setResistor(100) : 100). ---
    Cbus = {int(r[0]): 0.0 for r in t.Y_DER}
    for _id, i, j, km, _r, _x, b_km in t.T_Lines_SI:
        b_pu = b_km * km * (UN[i] ** 2 / Sb)
        Cbus[i] += b_pu / 2.0
        Cbus[j] += b_pu / 2.0
    for _id, i, j, _km, _r, _x, b_pu, Un, Sb_l in t.T_Lines_pu:
        b_pu = b_pu * (UN[i] ** 2 / Sb) / (Un**2 / Sb_l)  # base ligne -> systeme
        Cbus[i] += b_pu / 2.0
        Cbus[j] += b_pu / 2.0
    # Un condensateur (B > 0) est de la capacite de noeud de plus ; une
    # reactance (B < 0) devient une branche R-L vers la terre, plus bas
    # (ReactanceShunt) : fondue dans Cbus, elle serait une capacite NEGATIVE.
    for _id, b, B_pu in t.T_Shunts:
        if B_pu > 0:
            Cbus[b] += B_pu
    for b, c in Cbus.items():
        if c < 0:
            raise ValueError(
                f"bus {b} : capacite shunt negative ({c:g} pu). Une susceptance "
                "de ligne negative n'a pas de sens physique ; une reactance se "
                "declare dans T_Shunts avec B_pu < 0."
            )

    # --- charges speciales de T_Loads : leur consommation s'ajoute a celle de
    #     Y_DER dans le load flow ---
    S_spec = {int(r[0]): [0.0, 0.0] for r in t.Y_DER}
    for _id, b, P, Q, typ in t.T_Loads:
        if typ not in CHARGES_SPECIALES:
            raise ValueError(
                f"T_Loads, charge {_id} : type '{typ}' inconnu ('I' ou 'P')"
            )
        S_spec[b][0] += P
        S_spec[b][1] += Q

    net = Network(f0=t.f0, Sb=Sb, noeuds_sans_capacite=t.NOEUDS, r_causalite=t.R_NOEUDS)

    for bid, btype, Un, Vpu, theta, Pg, Qg, Pl, Ql in t.Y_DER:
        net.add(
            Bus(
                int(bid),
                Un,
                btype=int(btype),
                V=Vpu,
                theta=theta,
                Pg=Pg,
                Qg=Qg,
                Pl=Pl + S_spec[int(bid)][0],
                Ql=Ql + S_spec[int(bid)][1],
                C=Cbus[int(bid)],
            )
        )

    # --- les noeuds continus, et la capacite qu'ils portent ---
    Cdc = _capacites_dc(t, vscs, UN_DC, net.wb)
    for bid, dctype, Un, Vpu, Pinj in t.T_BusDC:
        bid, dctype = int(bid), int(dctype)
        if dctype not in (DC_SLACK, DC_P):
            raise ValueError(
                f"T_BusDC, bus {bid} : Type = {dctype} inconnu "
                f"({DC_SLACK} = slack continu, {DC_P} = P impose)"
            )
        net.add(BusDC(bid, Un, C=Cdc[bid], V=float(Vpu)))
        spec = net.bus_specs_dc[bid]
        spec.btype = dctype
        spec.Pinj = float(Pinj) / Sb

    # Chaque ligne garde la demi-charge capacitive qu'elle a versee a chacun de
    # ses deux noeuds (``b_demi``, pu systeme) : l'ouvrir en cours de simulation
    # doit la leur reprendre (g2egrid.temporel.Manoeuvre).
    for lid, i, j, km, r_km, x_km, b_km in t.T_Lines_SI:
        ZB = UN[i] ** 2 / Sb
        ligne = Line(
            f"L{int(lid)}_{i}_{j}",
            i,
            j,
            R=r_km * km / ZB,
            X=x_km * km / ZB,
            wb=net.wb,
        )
        ligne.b_demi = b_km * km * (UN[i] ** 2 / Sb) / 2.0
        net.add(ligne)

    for lid, i, j, _km, r_pu, x_pu, b_pu, Un, Sb_l in t.T_Lines_pu:
        k = (Un**2 / Sb_l) / (UN[i] ** 2 / Sb)  # base de la ligne -> base systeme
        ligne = Line(f"L{int(lid)}_{i}_{j}", i, j, R=r_pu * k, X=x_pu * k, wb=net.wb)
        ligne.b_demi = b_pu * (UN[i] ** 2 / Sb) / (Un**2 / Sb_l) / 2.0
        net.add(ligne)

    # Cables continus : R-L a etat, capacite deja versee aux deux noeuds ; la
    # configuration (monopole, bipole...) est ramenee au conducteur equivalent.
    for lid, i, j, km, r_km, l_km, _c in map(cable_dc, t.T_LinesDC):
        i, j = int(i), int(j)
        ZB = UN_DC[i] ** 2 / Sb
        net.add(
            LineDC(
                f"CDC{int(lid)}_{i}_{j}",
                i,
                j,
                R=r_km * km / ZB,
                X=net.wb * (l_km * 1e-3 * km) / ZB,
                wb=net.wb,
            )
        )

    # Convertisseurs DC/DC et charges continues, sur leur noeud continu.
    for did, b, pn, typ, p in t.T_DCDC:
        b, typ = int(b), str(typ).lower()
        if typ not in TYPES_DCDC:
            raise ValueError(
                f"T_DCDC, convertisseur {did} : type {typ!r} inconnu "
                f"({', '.join(TYPES_DCDC)})"
            )
        kS = float(pn) / Sb
        net.add(
            ConvertisseurDCDC(
                f"DCDC{int(did)}_{b}",
                b,
                P=float(p) / Sb,
                # une batterie soutient la tension : 20 pu de courant par pu de
                # tension, rapportes a sa puissance ; un PV tient sa puissance
                Kv=KV_BATTERIE * kS if typ == "batterie" else 0.0,
            )
        )
    for lid, b, p, typ in t.T_LoadsDC:
        b, typ = int(b), str(typ).upper()
        if typ == "P":
            net.add(ChargeDC_P(f"LoadDC{int(lid)}_{b}", b, P=float(p) / Sb))
        elif typ == "R":
            if p <= 0:
                raise ValueError(
                    f"T_LoadsDC, charge {lid} : P = {p} W pour une resistance"
                )
            net.add(ChargeDC_R(f"LoadDC{int(lid)}_{b}", b, R=1.0 / (float(p) / Sb)))
        else:
            raise ValueError(
                f"T_LoadsDC, charge {lid} : type {typ!r} inconnu ('P' ou 'R')"
            )

    # Une prise hors nominal se declare par U1_V/U2_V, comme dans Network.m :
    # n = (U1/U2)/(Un_i/Un_j), egal a 1 quand U1, U2 sont les tensions des bus.
    for tid, i, j, Sn, U1, U2, R, X in t.T_Transfos:
        k = (Sb / Sn) * (U1 / UN[i]) ** 2  # base du transfo -> base systeme
        n = (U1 / U2) / (UN[i] / UN[j])
        net.add(
            Transformer(f"T{int(tid)}_{i}_{j}", i, j, R=R * k, X=X * k, wb=net.wb, n=n)
        )

    for _id, b, Sn, _E, _th, R, X in t.T_IBs:
        net.add(IB(f"IB{b}", b, R=R, X=X, wb=net.wb, Sn=Sn, Sb=Sb))

    # Colonnes optionnelles Rt_pu, Xt_pu des tables de DER : transfo de
    # raccordement INTEGRE a l'equipement. Sans elles, la machine n'a pas de
    # transfo (a declarer dans T_Transfos) et le convertisseur garde son filtre
    # R2/L2 par defaut.
    #
    # Cle ``arbre`` des reglages d'une machine (famille ou bus) : le nom d'une
    # ligne d'arbre de ``g2egrid.elements.LIGNES_ARBRE``. La masse unique du
    # rotor devient cette ligne MULTI-MASSES, mise a l'echelle du ``H`` de la
    # machine ; les couples des troncons sont ses sorties ``C_<troncon>``.
    for _id, b, Sn, *Tr in t.T_SMs:
        reg = reglages(t.SM_PARAMS, b, f"SM{b}")
        arbre = reg.pop(ARBRE, None)
        if arbre and Tr:
            raise ValueError(
                f"T_SMs, bus {b} : une ligne d'arbre ({arbre!r}) ne se combine pas "
                "au transformateur integre ; declarer le transformateur dans "
                "T_Transfos"
            )
        if arbre:
            net.add(
                turbo_alternateur_nomme(
                    f"SM{b}", b, Sn=Sn, Sb=Sb, wb=net.wb, arbre=arbre, **reg
                )
            )
        elif Tr:
            net.add(
                SynchronousMachineTR(
                    f"SM{b}",
                    b,
                    Sn=Sn,
                    Sb=Sb,
                    wb=net.wb,
                    Rt=Tr[0],
                    Lt=Tr[1],
                    **reg,
                )
            )
        else:
            net.add(
                SynchronousMachine(
                    f"SM{b}",
                    b,
                    Sn=Sn,
                    Sb=Sb,
                    wb=net.wb,
                    **reg,
                )
            )

    for T_GFM, GFM in (
        (t.T_GFM1s, GFM1),
        (t.T_GFM2s, GFM2),
        (t.T_GFM3s, GFM3),
        (t.T_GFM4s, GFM4),
    ):
        for _id, b, Sn, _Vdc, *Tr in T_GFM:
            tr = dict(R2_pu=Tr[0], L2_pu=Tr[1]) if Tr else {}
            reg = reglages(t.GFM_PARAMS, b, f"GFM{b}") | tr
            classe_gfm = _loi_gfm(reg, GFM, b)
            if int(b) in vscs:
                _vid, bus_dc, mode, _table, classe = vscs[int(b)]
                nq_lf = reg.pop("statisme_qv", None)
                net.add(classe(f"GFM{b}", b, bus_dc, Sn=Sn, Sb=Sb, wb=net.wb, **reg))
                net.vsc_modes[int(b)] = mode
                if nq_lf:
                    net.vsc_qv[int(b)] = float(nq_lf)
                continue
            net.add(classe_gfm(f"GFM{b}", b, Sn=Sn, Sb=Sb, wb=net.wb, **reg))

    for T_GFL, GFL in ((t.T_GFL1s, GFL1), (t.T_GFL2s, GFL2)):
        for _id, b, Sn, _Vdc, *Tr in T_GFL:
            tr = dict(R2_pu=Tr[0], L2_pu=Tr[1]) if Tr else {}
            reg = reglages(t.GFL_PARAMS, b, f"GFL{b}") | tr
            modele = _modele_station(reg, b)
            # Pose sur la famille, le modele vaut pour les STATIONS (les GFL a
            # port continu) ; pose sur un bus sans port continu, c'est une
            # erreur de description.
            par_bus = t.PARAMS_PAR_BUS.get(int(b), {})
            mmc = modele in CLASSES_MMC
            if mmc and int(b) not in vscs and "modele" in par_bus:
                raise ValueError(
                    f"GFL du bus {int(b)} : un MMC a forcement un port continu "
                    "(T_GFL2s raccorde par T_VSCs)"
                )
            if int(b) in vscs:
                _vid, bus_dc, mode, _table, classe = vscs[int(b)]
                kdroop = reg.pop("kdroop", None)
                nq_lf = reg.pop("statisme_qv", None)
                if mmc:
                    # pas de filtre LCL : les reglages du VSC n'ont pas cours
                    classe = CLASSES_MMC[modele]
                    reg = {k: v for k, v in reg.items() if k in _CLES_MMC}
                # Le condensateur n'est plus dans le convertisseur, mais le
                # reglage de son Kpdc en a toujours besoin : on lui rend la
                # capacite REELLE de son noeud, ramenee a sa base (GFL2DC).
                reg = reg | dict(Cdc_pu=net.bus_specs_dc[bus_dc].C * (Sb / float(Sn)))
                el = net.add(
                    classe(f"GFL{b}", b, bus_dc, Sn=Sn, Sb=Sb, wb=net.wb, **reg)
                )
                net.vsc_modes[int(b)] = mode
                if nq_lf:
                    net.vsc_qv[int(b)] = float(nq_lf)
                if mode in ("droop", "droopI"):
                    net.vsc_statisme[int(b)] = kdroop
                    el.kdroop = kdroop
                    if mode == "droopI":
                        # le modele dynamique cale son Kpdc sur la pente du
                        # statisme au point de fonctionnement (cf. GFL2DC)
                        el.droop_courant = dict(
                            vdc_ref=float(net.bus_specs_dc[bus_dc].V),
                            p_ref=float(
                                net.bus_specs[int(b)].Pg - net.bus_specs[int(b)].Pl
                            )
                            * Sb
                            / float(Sn),
                        )
                continue
            net.add(GFL(f"GFL{b}", b, Sn=Sn, Sb=Sb, wb=net.wb, **reg))

    # Elements decrits par un schema-bloc : la classe est compilee (et mise en
    # cache) depuis le modele ; les parametres du schema valent par defaut, le
    # reseau fournit les reserves (wb, kS = Sn/Sb), PARAMS_PAR_BUS surcharge.
    for bid_, b, modele, *Sn in t.T_Blocs:
        if modele not in t.MODELES_BLOCS:
            raise ValueError(
                f"T_Blocs, element {int(bid_)} : modele {modele!r} absent de "
                "MODELES_BLOCS"
            )
        classe = classe_element(t.MODELES_BLOCS[modele])
        reg = dict(classe.compilation.parametres)
        reserves = dict(wb=net.wb, kS=float(Sn[0]) / Sb if Sn and Sn[0] else 1.0)
        reg.update({k: v for k, v in reserves.items() if k in reg})
        reg.update(
            {
                k: float(v)
                for k, v in t.PARAMS_PAR_BUS.get(int(b), {}).items()
                if k in reg
            }
        )
        net.add(classe(f"Blocs{int(bid_)}_{int(b)}", int(b), **reg))

    # Charges classiques : colonnes Pl_W, Ql_var de Y_DER -> charge R-L a
    # impedance constante (GenericLoad.m). Une inductance ne peut pas consommer
    # un reactif negatif ou nul : condensateur dans T_Shunts, ou charge speciale.
    for bid, *_, Pl, Ql in t.Y_DER:
        if Pl == 0 and Ql == 0:
            continue
        if Ql <= 0:
            raise ValueError(
                f"bus {int(bid)} : Ql = {Ql:g} var. Une charge R-L demande Ql > 0 ; "
                "declarer le condensateur dans T_Shunts ou la charge dans T_Loads "
                "('P' ou 'I')"
            )
        net.add(RLLoad(f"Load{int(bid)}", int(bid), P=Pl / Sb, Q=Ql / Sb, wb=net.wb))

    for lid, b, P, Q, typ in t.T_Loads:
        net.add(
            CHARGES_SPECIALES[typ](f"Load{typ}{int(lid)}_{b}", b, P=P / Sb, Q=Q / Sb)
        )

    # Reactances shunt, en dernier : un reseau sans reactance garde la liste
    # d'elements -- et donc l'ordre des etats -- qu'il avait avant elles. Les
    # condensateurs ne sont pas des elements (ils sont dans la capacite du
    # noeud) ; on garde leur nom et leur B pour pouvoir les manoeuvrer.
    net.condensateurs = {}
    for sid, b, B_pu in t.T_Shunts:
        if B_pu < 0:
            net.add(
                ReactanceShunt(f"Xsh{int(sid)}_{int(b)}", int(b), B=B_pu, wb=net.wb)
            )
        else:
            net.condensateurs[f"Csh{int(sid)}_{int(b)}"] = (int(b), float(B_pu))

    appliquer_dynamiques(net, _avec_simplifications(net, selections, t.SIMPLIFICATIONS))
    if t.REFERENCE == "coi":
        net.reference_flottante("coi")
    elif t.REFERENCE == "nominale":
        net.reference_flottante(None)
    return net


def _avec_simplifications(
    net: Network, selections: dict[str, dict], simplifications: dict[str, dict]
) -> dict[str, dict]:
    """Les selections de dynamiques de chaque element, hypotheses de
    simplification comprises : type, puis famille ou bus, puis l'element
    lui-meme -- le plus precis l'emporte, etat par etat."""
    if not simplifications:
        return selections
    noms = {el.name for el in net.all_elements}
    classes = {type(el).__name__ for el in net.all_elements}
    inconnus = sorted(set(simplifications) - noms - classes)
    if inconnus:
        raise ValueError(
            f"SIMPLIFICATIONS : {', '.join(inconnus)} n'est ni un element ni un "
            "type d'element de ce reseau"
        )
    out = {k: dict(v) for k, v in selections.items()}
    for el in net.all_elements:
        par_type = simplifications.get(type(el).__name__, {})
        par_nom = simplifications.get(el.name, {})
        if par_type or par_nom:
            out[el.name] = {**par_type, **out.get(el.name, {}), **par_nom}
    return out


def appliquer_dynamiques(net: Network, selections: dict[str, dict]) -> None:
    """Applique a chaque element sa selection de dynamiques.

    ``selections`` : ``{nom d'element: {etat: sort}}``, ou ``sort`` vaut
    ``"statique"`` (``dx/dt = 0``), ``"fige"`` (``x = x0``), un nombre (``x``
    fige a cette valeur) ou ``"dynamique"`` (rien a faire : c'est ce qui permet
    a un bus de rendre une dynamique que sa famille supprime). Voir
    :meth:`g2egrid.noyau.element.Element.reduire`.

    Une selection de FAMILLE (``GFM_PARAMS``) vaut pour des variantes qui n'ont
    pas toutes les memes etats -- ``Dw`` n'existe pas dans ``GFM1`` --, donc
    chaque element ignore les noms qu'il n'a pas. Mais un nom qu'AUCUN element
    ne connait est une faute de frappe, et elle est signalee.
    """
    par_nom = {el.name: el for el in net.all_elements}
    vus: set[str] = set()
    demandes: set[str] = set()
    for nom, dyn in selections.items():
        el = par_nom.get(nom)
        if el is None:
            continue  # element retire (disjoncteur ouvert, ile morte)
        statique, fige = [], {}
        for etat, sort in dyn.items():
            demandes.add(etat)
            if etat in el.states:
                vus.add(etat)
            if sort == "statique":
                statique.append(etat)
            elif sort == "fige":
                fige[etat] = None
            elif sort == "dynamique":
                continue
            elif isinstance(sort, (int, float)) and not isinstance(sort, bool):
                fige[etat] = float(sort)
            else:
                raise ValueError(
                    f"{nom}: dynamique {etat!r} -> {sort!r} ; attendu 'statique', "
                    "'fige', 'dynamique' ou une valeur numerique"
                )
        el.reduire(statique=statique, fige=fige, strict=False)
    inconnus = sorted(demandes - vus)
    if inconnus:
        raise ValueError(
            f"dynamiques : {', '.join(map(repr, inconnus))} n'est l'etat d'aucun "
            "element du reseau"
        )
