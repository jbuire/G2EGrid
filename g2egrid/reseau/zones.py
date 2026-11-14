"""Zones interconnectees : un reseau decrit UNE FOIS, instancie plusieurs fois.

POURQUOI. Un reseau de transport auquel se raccordent plusieurs reseaux de
distribution s'ecrivait jusqu'ici A PLAT : une seule grande table de bus,
numerotee a la main a partir de 100, 200, 300. C'est ce qu'a demande
``wscc3sm_cigre3`` (WSCC 3 machines + 3 distributions CIGRE, 54 bus), et c'est
intenable : la meme distribution y est recopiee trois fois, a trois
numerotations differentes. Corriger une impedance de depart demandait trois
corrections, et rien ne disait que les trois copies etaient encore identiques.

Une :class:`Zone` est ce reseau decrit une fois, en numerotation LOCALE et
SANS SLACK -- elle n'est pas un reseau autonome, elle est un morceau qui pend
au bout d'un bus de transport. Une :class:`Instance` la raccorde a un bus
precis, sous un prefixe de numerotation, et peut representer ``K`` departs
identiques.

CE QUI APPARTIENT A QUI. Le transfo de poste appartient a l'INSTANCE, pas a la
zone : sa puissance depend de ``K`` et son primaire est le bus amont, deux
choses que la zone ne connait pas. Le banc de condensateurs du poste source
(:attr:`Instance.shunt_pcc`) est dans le meme cas -- il est dimensionne site
par site, et les valeurs de ``wscc3sm_cigre3`` (0,35 pu derriere K = 4 ;
0,25 pu derriere K = 3) ne sont justement PAS dans un rapport K. Les transfos
de depart, eux, appartiennent a la zone : ils se repetent avec elle.

L'ECHELLE K. Agreger ``K`` departs identiques revient a multiplier les
puissances par ``K`` et a diviser les impedances par ``K``, si bien que les
tensions, les pertes relatives et les dynamiques EN PU restent celles d'un
seul depart. La regle n'est donc pas la meme d'une table a l'autre, et
:data:`ECHELLE` la donne colonne par colonne. Deux pieges y sont explicites :
``T_Lines_pu`` ne voit multiplier que sa colonne ``Sb_VA``, parce que ses pu
sont deja rapportees a cette base ; et ``T_VSCs`` ne bouge pas, parce qu'elle
ne porte aucune grandeur physique.

L'ORDRE D'EMISSION EST UN CONTRAT. :func:`aplatir` emet le transport, puis les
instances dans l'ordre, TABLE PAR TABLE. L'ordre d'insertion des elements fixe
l'ordre des etats, donc l'ordre des lignes et des colonnes de ``A`` : s'en
ecarter donnerait le meme reseau avec une matrice permutee, et la comparaison
avec le reseau ecrit a plat n'aurait plus de sens. C'est ce qui permet a
``test_zones`` de comparer ``A`` a 1e-12 pres.

LES REGLAGES PAR FAMILLE D'UNE ZONE SONT REPLIES PAR BUS. ``GFM_PARAMS`` vaut
pour toute la famille ; deux zones peuvent avoir des GFM regles differemment,
et un seul dictionnaire global ne saurait pas les distinguer. A
l'aplatissement, les reglages de famille d'une zone deviennent donc des
entrees de ``PARAMS_PAR_BUS``, qui, elle, est indexee par bus et survit au
melange.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, NamedTuple

from g2egrid.elements.bus import NOEUDS_DEFAUT, R_CAUSALITE
from g2egrid.reseau.tables import TABLES, Tables

__all__ = [
    "extraire_zone",
    "importer_reseau",
    "Instance",
    "Reseau",
    "TransfoPoste",
    "Zone",
    "aplatir",
    "renumeroter",
    "reseau_depuis_json",
    "reseau_vers_json",
]

#: Convention conservee des reseaux ecrits a plat : une zone commence a 100.
#: Plus bas, ses bus se confondraient avec ceux d'un transport ordinaire.
PREFIXE_MIN = 100

#: ``Y_DER.Type`` du slack. Une zone n'en a pas : elle pend au transport.
SLACK = 1
PV, PQ = 2, 3

#: Ce que ``K`` multiplie (``"*"``) ou divise (``"/"``), table par table et
#: colonne par colonne. Une colonne absente ne bouge pas.
ECHELLE: dict[str, dict[int, str]] = {
    "Y_DER": {5: "*", 6: "*", 7: "*", 8: "*"},  # Pg, Qg, Pl, Ql
    "T_Lines_SI": {4: "/", 5: "/", 6: "*"},  # R/km, X/km, B/km
    "T_Lines_pu": {8: "*"},  # Sb_VA SEULEMENT : les pu y sont deja rapportees
    "T_Shunts": {2: "*"},  # B_pu
    "T_Transfos": {3: "*"},  # Sn_VA
    "T_IBs": {2: "*"},  # Sn_VA
    "T_SMs": {2: "*"},  # Sn_VA
    "T_GFM1s": {2: "*"},
    "T_GFM2s": {2: "*"},
    "T_GFM3s": {2: "*"},
    "T_GFM4s": {2: "*"},
    "T_GFL1s": {2: "*"},
    "T_GFL2s": {2: "*"},
    "T_Loads": {2: "*", 3: "*"},  # P_W, Q_var
    "T_BusDC": {4: "*"},  # Pinj_W
    "T_LinesDC": {4: "/", 5: "/", 6: "*"},  # R/km, L/km, C/km
    "T_ShuntsDC": {2: "*"},  # C_pu
    "T_VSCs": {},  # aucune grandeur physique
    "T_DCDC": {2: "*", 4: "*"},  # Pn_W, P_W
    "T_LoadsDC": {2: "*"},  # P_W
    "T_Blocs": {3: "*"},  # Sn_VA (la puissance, elle, est celle du bus)
}

#: Les colonnes qui designent un BUS ALTERNATIF. Elles prennent le prefixe, et
#: ce sont les seules que la fusion du PCC avec son bus amont peut detourner.
BUS_AC: dict[str, tuple[int, ...]] = {
    "Y_DER": (0,),
    "T_Lines_SI": (1, 2),
    "T_Lines_pu": (1, 2),
    "T_Shunts": (1,),
    "T_Transfos": (1, 2),
    "T_IBs": (1,),
    "T_SMs": (1,),
    "T_GFM1s": (1,),
    "T_GFM2s": (1,),
    "T_GFM3s": (1,),
    "T_GFM4s": (1,),
    "T_GFL1s": (1,),
    "T_GFL2s": (1,),
    "T_Loads": (1,),
    "T_VSCs": (1,),
    "T_Blocs": (1,),
}

#: Les colonnes qui designent un BUS CONTINU. Numerotation distincte de celle
#: des bus alternatifs, mais meme prefixe : une instance garde ses deux
#: numerotations groupees.
BUS_DC: dict[str, tuple[int, ...]] = {
    "T_BusDC": (0,),
    "T_LinesDC": (1, 2),
    "T_ShuntsDC": (1,),
    "T_VSCs": (2,),
    "T_DCDC": (1,),
    "T_LoadsDC": (1,),
}

#: Les colonnes qui portent l'IDENTIFIANT de l'element. Elles prennent le
#: prefixe elles aussi, parce que le nom de l'element assemble en est bati
#: (``L{ID}_{From}_{To}``) : sans cela, trois instances d'une meme zone
#: donneraient trois elements de meme nom.
IDS: dict[str, tuple[int, ...]] = {
    nom: (0,) for nom in TABLES if nom not in ("Y_DER", "T_BusDC")
}

#: Les tables de DER, et le dictionnaire de reglages de famille dont elles
#: relevent. Sert au repliement des reglages d'une zone dans PARAMS_PAR_BUS.
FAMILLES: dict[str, str] = {
    "T_SMs": "SM_PARAMS",
    "T_GFM1s": "GFM_PARAMS",
    "T_GFM2s": "GFM_PARAMS",
    "T_GFM3s": "GFM_PARAMS",
    "T_GFM4s": "GFM_PARAMS",
    "T_GFL1s": "GFL_PARAMS",
    "T_GFL2s": "GFL_PARAMS",
}


class TransfoPoste(NamedTuple):
    """Le transfo de poste d'une instance, pour UN depart.

    ``Sn_VA`` est la puissance d'un seul depart : :func:`aplatir` la multiplie
    par ``K``. Les deux bus, eux, ne sont pas ici -- le primaire est le bus
    amont de l'instance et le secondaire le PCC de la zone, et c'est l'instance
    qui les connait.
    """

    Sn_VA: float
    U1_V: float
    U2_V: float
    R_pu: float
    X_pu: float


@dataclass
class Zone:
    """Un reseau decrit en numerotation LOCALE, sans slack, pret a instancier.

    ``tables`` est un :class:`~g2egrid.reseau.tables.Tables` ordinaire, a ceci pres
    qu'aucun bus n'y est slack : une zone ne tient pas sa propre tension, elle
    la recoit du transport par son ``bus_pcc``. ``bus_pcc`` est le point de
    raccordement, cote zone -- le jeu de barres du poste source quand il y a un
    transfo de poste, le bus a fusionner avec le bus amont quand il n'y en a
    pas.
    """

    tables: Tables
    bus_pcc: int = 0
    positions: dict[int, tuple[float, float]] = field(default_factory=dict)
    nom: str = "zone"

    def charge_totale(self) -> tuple[float, float]:
        """La charge d'UN depart, en W et var (cf. ``Tables.charge_totale``)."""
        return self.tables.charge_totale()


@dataclass
class Instance:
    """Un raccordement d'une zone a un bus du transport.

    ``zone`` est le NOM d'une zone de :attr:`Reseau.zones`, pas l'objet : c'est
    ce nom qui se lit dans le fichier ``.g2egrid.json``, et le garder evite deux
    verites sur la meme zone. ``K`` est le nombre de departs identiques que
    cette instance agrege.

    ``transfo`` a ``None`` fait FUSIONNER le PCC de la zone avec ``bus_amont``
    (memes charges, memes shunts, un seul bus) ; les deux ``Un`` doivent alors
    coincider. ``charge_residuelle`` retranche ``K x charge_totale()`` de la
    zone a la charge du bus amont : le bus de transport garde exactement la
    consommation qu'il avait avant qu'on lui accroche une distribution, et le
    point de fonctionnement du transport ne bouge pas.
    """

    zone: str
    prefixe: int
    bus_amont: int
    K: int = 1
    transfo: TransfoPoste | tuple | None = None
    charge_residuelle: bool = False
    #: banc de condensateurs du poste source, en pu base SYSTEME. Valeur
    #: FINALE, non multipliee par K : un gradin se dimensionne pour le poste,
    #: pas par depart.
    shunt_pcc: float | None = None


@dataclass
class Reseau:
    """Un transport, des zones, et les instances qui les y raccordent."""

    transport: Tables
    zones: dict[str, Zone] = field(default_factory=dict)
    instances: list[Instance] = field(default_factory=list)
    nom: str = "reseau"
    positions: dict[int, tuple[float, float]] = field(default_factory=dict)
    f_zoom: float = 5.0

    def tables(self) -> Tables:
        """Le reseau a plat (raccourci de :func:`aplatir`)."""
        return aplatir(self)


# ---------------------------------------------------------------- numerotation
def _decale(v, prefixe: int, fusion: dict[int, int] | None):
    """Un numero local, vu du reseau complet."""
    n = int(v)
    if fusion is not None and n in fusion:
        return fusion[n]
    return n + prefixe


def renumeroter(
    t: Tables,
    prefixe: int,
    K: int = 1,
    fusion: dict[int, int] | None = None,
) -> Tables:
    """Les tables d'une zone, vues du reseau complet.

    Deux operations, et rien d'autre : les numeros locaux (bus et identifiants
    d'elements) prennent le ``prefixe``, et les grandeurs physiques prennent
    l'echelle ``K`` de :data:`ECHELLE`. ``fusion`` donne les bus ALTERNATIFS
    locaux qui, au lieu du prefixe, prennent un numero impose -- c'est ainsi
    que le PCC d'une zone sans transfo de poste devient le bus amont.

    La fonction ne fusionne PAS les lignes de ``Y_DER`` qu'elle amene a porter
    le meme numero : renommer et fusionner sont deux gestes distincts, et
    fusionner demande de savoir quelle ligne fait foi, ce que seul
    :func:`aplatir` sait. ``PARAMS_PAR_BUS`` voit ses cles decalees et ses
    VALEURS INCHANGEES : un reglage ne depend pas du nombre de departs.
    """
    if K < 1:
        raise ValueError(f"K = {K} : il faut au moins un depart")
    out = t.copie()
    for nom in TABLES:
        echelle = ECHELLE[nom]
        ac, dc, ids = BUS_AC.get(nom, ()), BUS_DC.get(nom, ()), IDS.get(nom, ())
        lignes = []
        for ligne in getattr(t, nom):
            r = list(ligne)
            for c in ids:
                r[c] = int(r[c]) + prefixe
            for c in ac:
                r[c] = _decale(r[c], prefixe, fusion)
            for c in dc:
                r[c] = int(r[c]) + prefixe
            if K != 1:
                for c, op in echelle.items():
                    if c < len(r):
                        r[c] = r[c] * K if op == "*" else r[c] / K
            lignes.append(tuple(r))
        setattr(out, nom, lignes)
    out.PARAMS_PAR_BUS = {
        _decale(b, prefixe, fusion): dict(v) for b, v in t.PARAMS_PAR_BUS.items()
    }
    # un disjoncteur ouvert dans la zone l'est dans chaque instance : ses ID
    # prennent le prefixe comme ceux des elements qu'ils designent (un
    # delestage de Y_DER designe un BUS, qui peut fusionner avec le bus amont)
    out.OUVERTS = [
        (a, _decale(b, prefixe, fusion) if a == "Y_DER" else int(b) + prefixe)
        for a, b in t.OUVERTS
    ]
    return out


def _replier_reglages(t: Tables) -> Tables:
    """Les reglages de famille d'une zone, portes sur ses bus.

    ``GFM_PARAMS`` vaut pour toute la famille : une fois deux zones melees dans
    une meme table, il ne sait plus dire de quels GFM il parle. On le reporte
    donc sur chaque bus porteur, dans ``PARAMS_PAR_BUS``, que l'assembleur lit
    bus par bus et qui survit au melange. Un reglage deja present par bus
    l'emporte, comme dans ``tables.assembler``.
    """
    out = t.copie()
    par_bus = {int(b): dict(v) for b, v in t.PARAMS_PAR_BUS.items()}
    for table, famille in FAMILLES.items():
        reglages = getattr(t, famille)
        if not reglages:
            continue
        for ligne in getattr(t, table):
            b = int(ligne[1])
            par_bus[b] = dict(reglages) | par_bus.get(b, {})
    out.PARAMS_PAR_BUS = par_bus
    out.SM_PARAMS, out.GFM_PARAMS, out.GFL_PARAMS = {}, {}, {}
    return out


# ----------------------------------------------------------------- validations
def _numeros_ac(t: Tables) -> set[int]:
    return {int(r[0]) for r in t.Y_DER}


def _verifier_zone(z: Zone, nom: str, transport: Tables) -> None:
    slacks = [int(r[0]) for r in z.tables.Y_DER if int(r[1]) == SLACK]
    if slacks:
        raise ValueError(
            f"zone '{nom}' : le bus {slacks[0]} y est slack (Type = {SLACK}). "
            "Une zone pend au transport, elle ne tient pas sa propre tension ; "
            "le slack est un bus du transport."
        )
    if z.bus_pcc not in _numeros_ac(z.tables):
        raise ValueError(
            f"zone '{nom}' : bus_pcc = {z.bus_pcc} absent de son Y_DER local"
        )
    if abs(z.tables.f0 - transport.f0) > 1e-12:
        raise ValueError(
            f"zone '{nom}' : f0 = {z.tables.f0:g} Hz contre {transport.f0:g} Hz "
            "au transport. Un reseau n'a qu'une frequence."
        )
    if abs(z.tables.Sb - transport.Sb) > 1e-12 * max(transport.Sb, 1.0):
        raise ValueError(
            f"zone '{nom}' : Sb = {z.tables.Sb:g} VA contre {transport.Sb:g} VA "
            "au transport. Les colonnes en pu base SYSTEME (T_Shunts.B_pu, "
            "T_ShuntsDC.C_pu) changeraient de sens a l'aplatissement."
        )


def _verifier_familles(r: Reseau) -> None:
    """Un reglage de famille du transport ne doit pas deborder sur une zone.

    ``tables.assembler`` lit ``famille | PARAMS_PAR_BUS[bus]`` : une cle que le
    transport regle et que la zone ne regle pas traverserait donc le
    repliement et s'appliquerait aux machines de la zone, silencieusement. On
    le refuse plutot que de fabriquer un modele dont personne ne peut dire d'ou
    vient le reglage.
    """
    for nom, z in r.zones.items():
        utilisees = {FAMILLES[table] for table in FAMILLES if getattr(z.tables, table)}
        for famille in sorted(utilisees):
            cles_t = set(getattr(r.transport, famille))
            cles_z = set(getattr(z.tables, famille))
            manquantes = sorted(cles_t - cles_z)
            if manquantes:
                raise ValueError(
                    f"zone '{nom}' : le transport regle {famille}"
                    f"{manquantes} et la zone non. Ce reglage s'appliquerait "
                    "aussi aux machines de la zone, sans que rien ne le dise. "
                    f"Donner ces cles dans le {famille} de la zone, ou les "
                    "retirer du transport."
                )


def _verifier_instances(r: Reseau) -> None:
    bus_transport = _numeros_ac(r.transport)
    occupes: dict[int, str] = {}
    vus: set[int] = set()
    for k, inst in enumerate(r.instances):
        etiquette = f"instance {k} (zone '{inst.zone}', prefixe {inst.prefixe})"
        if inst.zone not in r.zones:
            raise ValueError(f"{etiquette} : zone inconnue ; {sorted(r.zones)}")
        if inst.prefixe < PREFIXE_MIN:
            raise ValueError(
                f"{etiquette} : prefixe < {PREFIXE_MIN}. La convention garde "
                "les deux premiers chiffres pour le transport."
            )
        if inst.prefixe in vus:
            raise ValueError(f"{etiquette} : ce prefixe est deja pris")
        vus.add(inst.prefixe)
        if inst.K < 1:
            raise ValueError(f"{etiquette} : K = {inst.K}, il faut K >= 1")
        if inst.bus_amont not in bus_transport:
            raise ValueError(
                f"{etiquette} : bus amont {inst.bus_amont} absent du Y_DER du transport"
            )
        z = r.zones[inst.zone]
        for b in _numeros_ac(z.tables):
            if inst.transfo is None and b == z.bus_pcc:
                continue  # ce bus-la fusionne, il n'occupe pas de numero neuf
            n = b + inst.prefixe
            if n in bus_transport:
                raise ValueError(
                    f"{etiquette} : le bus {n} existe deja dans le transport"
                )
            if n in occupes:
                raise ValueError(
                    f"{etiquette} : le bus {n} est deja pris par "
                    f"{occupes[n]}. Les plages de prefixes se chevauchent."
                )
            occupes[n] = etiquette
        if inst.transfo is None:
            un_zone = float(z.tables.bus(z.bus_pcc)[2])
            un_amont = float(r.transport.bus(inst.bus_amont)[2])
            if abs(un_zone - un_amont) > 1e-9 * max(un_amont, 1.0):
                raise ValueError(
                    f"{etiquette} : sans transfo de poste, le PCC "
                    f"{z.bus_pcc} ({un_zone:g} V) fusionne avec le bus amont "
                    f"{inst.bus_amont} ({un_amont:g} V) ; les deux Un doivent "
                    "etre egaux, ou il faut un transfo."
                )
        if inst.transfo is not None and any(
            int(row[0]) == 0 for row in z.tables.T_Transfos
        ):
            raise ValueError(
                f"{etiquette} : la zone a un transfo d'ID 0, qui est celui que "
                "le transfo de poste de l'instance occupe."
            )
        if inst.shunt_pcc is not None and any(
            int(row[0]) == 0 for row in z.tables.T_Shunts
        ):
            raise ValueError(
                f"{etiquette} : la zone a un shunt d'ID 0, qui est celui que "
                "le banc de condensateurs du poste occupe."
            )


# ---------------------------------------------------------------- aplatissement
def aplatir(r: Reseau) -> Tables:
    """Le reseau complet, en une seule table, prete pour ``tables.assembler``.

    Emet le transport, puis les instances dans l'ordre, TABLE PAR TABLE.
    L'ordre d'insertion fixe l'ordre des etats du modele assemble : le meme
    reseau emis dans un autre ordre donnerait la meme physique avec une matrice
    ``A`` permutee, et plus aucune comparaison ligne a ligne avec un reseau
    ecrit a plat ne tiendrait.

    Chaque instance emet, dans cet ordre : son transfo de poste puis les
    transfos de la zone, son banc de condensateurs puis les shunts de la zone.
    Le poste vient avant ce qu'il alimente, comme on le lit sur un schema.
    """
    for nom, z in r.zones.items():
        _verifier_zone(z, nom, r.transport)
    _verifier_familles(r)
    _verifier_instances(r)

    t = r.transport.copie(nom=r.nom)
    morceaux: list[Tables] = []
    for inst in r.instances:
        z = r.zones[inst.zone]
        fusion = None if inst.transfo is not None else {z.bus_pcc: inst.bus_amont}
        zt = renumeroter(_replier_reglages(z.tables), inst.prefixe, inst.K, fusion)
        pcc = inst.bus_amont if inst.transfo is None else z.bus_pcc + inst.prefixe

        # --- le poste : transfo et gradin, qui appartiennent a l'instance ---
        if inst.transfo is not None:
            Sn, U1, U2, R, X = tuple(inst.transfo)
            zt.T_Transfos = [
                (inst.prefixe, inst.bus_amont, pcc, Sn * inst.K, U1, U2, R, X),
                *zt.T_Transfos,
            ]
        if inst.shunt_pcc is not None:
            zt.T_Shunts = [(inst.prefixe, pcc, inst.shunt_pcc), *zt.T_Shunts]

        # --- fusion du PCC avec le bus amont : un seul bus, charges ajoutees ---
        if inst.transfo is None:
            ligne_pcc = next(r_ for r_ in zt.Y_DER if int(r_[0]) == inst.bus_amont)
            zt.Y_DER = [r_ for r_ in zt.Y_DER if int(r_[0]) != inst.bus_amont]
            t.Y_DER = [
                (
                    b,
                    ty,
                    Un,
                    V,
                    th,
                    Pg + ligne_pcc[5],
                    Qg + ligne_pcc[6],
                    Pl + ligne_pcc[7],
                    Ql + ligne_pcc[8],
                )
                if int(b) == inst.bus_amont
                else (b, ty, Un, V, th, Pg, Qg, Pl, Ql)
                for b, ty, Un, V, th, Pg, Qg, Pl, Ql in t.Y_DER
            ]

        # --- la charge que la distribution reprend au bus amont ---
        if inst.charge_residuelle:
            P, Q = z.tables.charge_totale()
            K = inst.K
            t.Y_DER = [
                (b, ty, Un, V, th, Pg, Qg, Pl - K * P, Ql - K * Q)
                if int(b) == inst.bus_amont
                else (b, ty, Un, V, th, Pg, Qg, Pl, Ql)
                for b, ty, Un, V, th, Pg, Qg, Pl, Ql in t.Y_DER
            ]
        morceaux.append(zt)

    for nom in TABLES:
        lignes = list(getattr(t, nom))
        for zt in morceaux:
            lignes.extend(getattr(zt, nom))
        setattr(t, nom, lignes)
    for zt in morceaux:
        t.PARAMS_PAR_BUS.update(zt.PARAMS_PAR_BUS)
        t.MODELES_BLOCS = {**t.MODELES_BLOCS, **zt.MODELES_BLOCS}
        t.OUVERTS = [*t.OUVERTS, *zt.OUVERTS]
    return t


# ------------------------------------------------------------------- le fichier
FORMAT = "g2egrid"
#: les marques de format LUES : celle de l'outil, et les deux noms qu'il a
#: portes (g2egrid en developpement, g2egrid une fois publie) -- un fichier
#: enregistre par l'un s'ouvre dans l'autre. Ecrites en deux morceaux pour
#: que le renommage de la publication ne les confonde pas.
FORMATS_LUS = (FORMAT, "g2" + "elin", "g2" + "egrid")
VERSION = 1
_REGLAGES = ("SM_PARAMS", "GFM_PARAMS", "GFL_PARAMS")


def _tables_seules(t: Tables) -> dict:
    return {nom: [list(ligne) for ligne in getattr(t, nom)] for nom in TABLES}


def _reglages_vers_json(t: Tables) -> dict:
    d: dict[str, Any] = {k: dict(getattr(t, k)) for k in _REGLAGES}
    d["PARAMS_PAR_BUS"] = {str(k): dict(v) for k, v in t.PARAMS_PAR_BUS.items()}
    # l'etat des disjoncteurs voyage avec les reglages : ce n'est pas une table
    # de construction, et le bloc des reglages traverse deja tous les chemins
    # du fichier
    d["OUVERTS"] = [[str(a), int(b)] for a, b in t.OUVERTS]
    d["SIMPLIFICATIONS"] = {k: dict(v) for k, v in t.SIMPLIFICATIONS.items()}
    d["MODELES_BLOCS"] = deepcopy(t.MODELES_BLOCS)
    # le traitement des noeuds sans capacite, quand ce n'est pas le defaut
    # (les reseaux de parite MATLAB gardent la resistance de Node.m)
    if t.NOEUDS != NOEUDS_DEFAUT:
        d["NOEUDS"] = t.NOEUDS
    if t.R_NOEUDS != R_CAUSALITE:
        d["R_NOEUDS"] = float(t.R_NOEUDS)
    return {k: v for k, v in d.items() if v}


def _tables_depuis_json(d: dict, reglages: dict, f0: float, Sb: float, nom: str):
    t = Tables(
        f0=float(f0),
        Sb=float(Sb),
        nom=nom,
        SM_PARAMS=dict(reglages.get("SM_PARAMS") or {}),
        GFM_PARAMS=dict(reglages.get("GFM_PARAMS") or {}),
        GFL_PARAMS=dict(reglages.get("GFL_PARAMS") or {}),
        PARAMS_PAR_BUS={
            int(k): dict(v) for k, v in (reglages.get("PARAMS_PAR_BUS") or {}).items()
        },
        OUVERTS=[(str(a), int(b)) for a, b in (reglages.get("OUVERTS") or [])],
        SIMPLIFICATIONS={
            str(k): dict(v) for k, v in (reglages.get("SIMPLIFICATIONS") or {}).items()
        },
        MODELES_BLOCS=deepcopy(dict(reglages.get("MODELES_BLOCS") or {})),
        NOEUDS=str(reglages.get("NOEUDS") or NOEUDS_DEFAUT),
        R_NOEUDS=float(reglages.get("R_NOEUDS") or R_CAUSALITE),
    )
    for nom_table in TABLES:
        setattr(t, nom_table, [tuple(r) for r in (d.get(nom_table) or [])])
    return t


def _positions_vers_json(p: dict) -> dict:
    return {str(k): [float(v[0]), float(v[1])] for k, v in p.items()}


def _positions_depuis_json(p: dict | None) -> dict:
    return {int(k): (float(v[0]), float(v[1])) for k, v in (p or {}).items()}


def reseau_vers_json(r: Reseau) -> dict:
    """Le reseau, en JSON pur -- le contenu d'un fichier ``.g2egrid.json``.

    ``f0`` et ``Sb`` sont en tete et non dans chaque bloc de tables : le
    fichier decrit UN reseau, et les validations de :func:`aplatir` exigent de
    toute facon que transport et zones s'accordent dessus. Les rendre une fois
    evite un fichier qui se contredit lui-meme.
    """
    return {
        "format": FORMAT,
        "version": VERSION,
        "nom": r.nom,
        "f_zoom": r.f_zoom,
        "f0": r.transport.f0,
        "Sb": r.transport.Sb,
        "reglages": _reglages_vers_json(r.transport),
        "transport": _tables_seules(r.transport),
        "zones": {
            nom: {
                "tables": _tables_seules(z.tables),
                "bus_pcc": z.bus_pcc,
                "reglages": _reglages_vers_json(z.tables),
                "positions": _positions_vers_json(z.positions),
            }
            for nom, z in r.zones.items()
        },
        "instances": [
            {
                "zone": i.zone,
                "prefixe": i.prefixe,
                "bus_amont": i.bus_amont,
                "K": i.K,
                "transfo": None if i.transfo is None else list(i.transfo),
                "charge_residuelle": bool(i.charge_residuelle),
                "shunt_pcc": i.shunt_pcc,
            }
            for i in r.instances
        ],
        "positions": _positions_vers_json(r.positions),
    }


def reseau_depuis_json(d: dict) -> Reseau:
    """L'inverse de :func:`reseau_vers_json`."""
    if d.get("format") not in FORMATS_LUS:
        raise ValueError(f"format '{d.get('format')}' inconnu ; attendu '{FORMAT}'")
    if int(d.get("version", 0)) > VERSION:
        raise ValueError(
            f"version {d.get('version')} : ce paquet lit jusqu'a la {VERSION}"
        )
    f0, Sb = d.get("f0", 50.0), d.get("Sb", 100e6)
    nom = str(d.get("nom", "reseau"))
    transport = _tables_depuis_json(
        d.get("transport") or {}, d.get("reglages") or {}, f0, Sb, nom
    )
    zones = {}
    for nom_zone, z in (d.get("zones") or {}).items():
        zones[nom_zone] = Zone(
            tables=_tables_depuis_json(
                z.get("tables") or {}, z.get("reglages") or {}, f0, Sb, nom_zone
            ),
            bus_pcc=int(z.get("bus_pcc", 0)),
            positions=_positions_depuis_json(z.get("positions")),
            nom=nom_zone,
        )
    instances = []
    for i in d.get("instances") or []:
        tr = i.get("transfo")
        instances.append(
            Instance(
                zone=str(i["zone"]),
                prefixe=int(i["prefixe"]),
                bus_amont=int(i["bus_amont"]),
                K=int(i.get("K", 1)),
                transfo=None if tr is None else TransfoPoste(*tr),
                charge_residuelle=bool(i.get("charge_residuelle", False)),
                shunt_pcc=None if i.get("shunt_pcc") is None else float(i["shunt_pcc"]),
            )
        )
    return Reseau(
        transport=transport,
        zones=zones,
        instances=instances,
        nom=nom,
        positions=_positions_depuis_json(d.get("positions")),
        f_zoom=float(d.get("f_zoom", 5.0)),
    )


# ------------------------------------------------------------------ extraction
def extraire_zone(r: Reseau, bus: list[int] | set[int], nom: str = "zone") -> Reseau:
    """Fait des bus ``bus`` du transport une ZONE, instanciee a la meme place.

    Le geste de l'atelier "definir une zone" : on selectionne un ensemble de
    bus, et ils deviennent une distribution, decrite une fois et raccordee au
    transport par une :class:`Instance`. Le reseau APLATI est le meme qu'avant
    -- memes elements, memes valeurs --, seuls les numeros des bus de la zone
    changent (``prefixe + numero local``).

    La zone se raccorde au reste par UNE SEULE branche :

    * un **transformateur** devient le transfo de poste de l'instance (le
      secondaire, cote zone, est son PCC) ;
    * une **ligne** reste dans la zone ; son extremite cote transport y entre
      comme PCC, un bus sans charge ni production que l'aplatissement FUSIONNE
      avec le bus amont.

    Refuse une zone qui contient le slack, qui touche au reseau continu ou qui
    se raccorde par plusieurs branches : aucune ne se decrirait comme une
    instance.
    """
    t = r.transport
    zone = {int(b) for b in bus}
    tous = _numeros_ac(t)
    inconnus = zone - tous
    if inconnus:
        raise ValueError(f"bus {sorted(inconnus)} absents du transport")
    if not zone:
        raise ValueError("aucun bus selectionne")
    if nom in r.zones:
        raise ValueError(f"une zone s'appelle deja '{nom}'")
    slack = [int(x[0]) for x in t.Y_DER if int(x[0]) in zone and int(x[1]) == SLACK]
    if slack:
        raise ValueError(
            f"le bus {slack[0]} est le slack : une zone pend au transport, "
            "elle ne peut pas porter la reference"
        )
    if any(int(x[1]) in zone for x in t.T_VSCs):
        raise ValueError("la zone porte un convertisseur raccorde au reseau continu")

    # --- la frontiere : branches a cheval ---
    frontiere = []
    for nom_t in ("T_Lines_SI", "T_Lines_pu", "T_Transfos"):
        for k, x in enumerate(getattr(t, nom_t)):
            dedans = (int(x[1]) in zone, int(x[2]) in zone)
            if dedans[0] != dedans[1]:
                frontiere.append((nom_t, k, x, dedans))
    if len(frontiere) != 1:
        raise ValueError(
            f"la zone se raccorde au transport par {len(frontiere)} branche(s) ; "
            "il en faut exactement une (un transformateur de poste ou une ligne)"
        )
    table_f, k_f, br, dedans = frontiere[0]
    pcc_plat = int(br[2]) if dedans[1] else int(br[1])
    amont = int(br[1]) if dedans[1] else int(br[2])

    # --- numerotation locale et prefixe ---
    locaux = sorted(zone)
    local = {b: i + 1 for i, b in enumerate(locaux)}
    avec_ligne = table_f != "T_Transfos"
    if avec_ligne:
        local[amont] = len(locaux) + 1  # le PCC, copie du bus amont
    plus_grand = max([*tous, *(i.prefixe + 1000 for i in r.instances)], default=0)
    prefixe = max(PREFIXE_MIN, 100 * ((plus_grand // 100) + 1))

    zt = Tables(f0=t.f0, Sb=t.Sb, nom=nom)
    reste = t.copie()
    for nom_t in TABLES:
        cols_ac = BUS_AC.get(nom_t, ())
        a_garder, a_zone = [], []
        for k, x in enumerate(getattr(t, nom_t)):
            if nom_t == table_f and k == k_f:
                if avec_ligne:  # la ligne frontiere entre dans la zone
                    a_zone.append(x)
                continue  # le transfo frontiere devient celui du poste
            buses = [int(x[c]) for c in cols_ac]
            if buses and all(b in zone for b in buses):
                a_zone.append(x)
            elif any(b in zone for b in buses):
                raise ValueError(
                    f"{nom_t} : l'element {x[0]} est a cheval sur la frontiere"
                )
            else:
                a_garder.append(x)
        # renumerotation locale : bus et identifiants
        ids = IDS.get(nom_t, ())
        zlignes = []
        for n_id, x in enumerate(a_zone, start=1):
            y = list(x)
            for c in cols_ac:
                y[c] = local[int(y[c])]
            for c in ids:
                y[c] = n_id
            zlignes.append(tuple(y))
        setattr(zt, nom_t, zlignes)
        setattr(reste, nom_t, a_garder)

    if avec_ligne:
        ligne_amont = next(x for x in t.Y_DER if int(x[0]) == amont)
        zt.Y_DER = [
            *zt.Y_DER,
            (local[amont], 3, ligne_amont[2], ligne_amont[3], 0, 0, 0, 0, 0),
        ]
        bus_pcc, transfo = local[amont], None
    else:
        _id, i, j, Sn, U1, U2, R, X = br
        if int(i) != amont:  # le poste se lit de l'amont vers la zone
            U1, U2 = U2, U1
        bus_pcc, transfo = local[pcc_plat], TransfoPoste(Sn, U1, U2, R, X)

    zt.PARAMS_PAR_BUS = {
        local[b]: dict(v) for b, v in t.PARAMS_PAR_BUS.items() if int(b) in zone
    }
    reste.PARAMS_PAR_BUS = {
        b: v for b, v in t.PARAMS_PAR_BUS.items() if int(b) not in zone
    }
    for fam in _REGLAGES:  # la zone herite des reglages de famille du transport
        setattr(zt, fam, dict(getattr(t, fam)))
    reste.OUVERTS = [
        o for o in t.OUVERTS if not (o[0] == "Y_DER" and int(o[1]) in zone)
    ]

    pos_zone = {local[b]: p for b, p in r.positions.items() if int(b) in zone}
    if avec_ligne and amont in r.positions:
        pos_zone[local[amont]] = r.positions[amont]
    positions = {b: p for b, p in r.positions.items() if int(b) not in zone}

    zones = dict(r.zones)
    zones[nom] = Zone(tables=zt, bus_pcc=bus_pcc, positions=pos_zone, nom=nom)
    instances = [
        *r.instances,
        Instance(zone=nom, prefixe=prefixe, bus_amont=amont, K=1, transfo=transfo),
    ]
    return Reseau(
        transport=reste,
        zones=zones,
        instances=instances,
        nom=r.nom,
        positions=positions,
        f_zoom=r.f_zoom,
    )


# ------------------------------------------------------------ importer
def _numeros(t: Tables) -> list[int]:
    """Tous les numeros qu'occupe un reseau a plat : bus (alternatifs et
    continus) et identifiants d'elements."""
    out = []
    for nom in TABLES:
        cols = (*BUS_AC.get(nom, ()), *BUS_DC.get(nom, ()), *IDS.get(nom, ()))
        for r in getattr(t, nom):
            out.extend(int(r[c]) for c in cols)
    return out


def _production_des_slacks(t: Tables) -> dict[int, float] | None:
    """La production [W] de chaque slack, lue sur le load flow du reseau
    seul ; ``None`` s'il ne converge pas."""
    from g2egrid.reseau.tables import assembler

    try:
        net = assembler(t)
        lf = net.run_loadflow()
    except Exception:
        return None
    idx = net._bus_index()
    return {
        int(r[0]): (float(lf.S[idx[int(r[0])]].real) + net.bus_specs[int(r[0])].Pl)
        * net.Sb
        for r in t.Y_DER
        if int(r[1]) == SLACK
    }


def importer_reseau(
    r: Reseau,
    source: Reseau | Tables,
    prefixe: int | None = None,
    raccord: dict | None = None,
    garder_slack: bool = False,
    garder_ib: bool = False,
) -> tuple[Reseau, dict]:
    """Importe ``source`` dans le TRANSPORT de ``r``, et l'y raccorde.

    Le geste de l'atelier "importer un reseau" : un second reseau vient se
    poser a cote du premier, pour y etre raccorde. Ce qui se passe :

    * **numerotation** -- bus et identifiants de ``source`` prennent un
      ``prefixe`` (par defaut la centaine au-dessus de tout numero deja pris,
      zones instanciees comprises) : rien ne se chevauche ;
    * **reference** -- un reseau raccorde n'a qu'une reference d'angle, celle
      de ``r``. Les **bus infinis** importes (``T_IBs``), qui representaient
      justement le reseau amont que ``r`` devient, sont retires, et leur bus
      devient un bus ordinaire (PQ) ; les autres **slacks** passent en PV,
      avec la production que leur donnait le load flow du reseau importe
      seul -- le point de fonctionnement importe ne bouge donc pas. ``garder_ib``
      et ``garder_slack`` desactivent l'un et l'autre ;
    * **reglages** -- ceux de famille de ``source`` sont reportes bus par bus
      (``PARAMS_PAR_BUS``) : ils ne se melent pas a ceux de ``r`` ;
    * **base** -- les seules grandeurs en pu de la base systeme (``T_Shunts``,
      ``T_ShuntsDC``) sont converties si ``Sb`` differe ; un reseau importe a une autre frequence prend celle de ``r``,
      reactances et pu inchangees ;
    * **raccordement** -- ``raccord = {"mode": ..., "bus_importe": b (numero
      d'ORIGINE), "bus_existant": B}`` : ``"fusion"`` fait de ``b`` et ``B``
      un seul bus (memes tensions nominales ; charges et productions
      s'ajoutent), ``"ligne"`` pose une ligne, ``"transfo"`` un
      transformateur, avec les valeurs par defaut de la palette. Sans
      ``raccord``, le reseau importe reste une ile, a raccorder a la main.

    Un ``source`` a zones est aplati d'abord : il arrive a plat dans le
    transport. Rend le nouveau reseau et ``{"prefixe", "bus", "messages"}``,
    ou ``bus`` donne le nouveau numero de chaque bus importe.
    """
    src = aplatir(source) if isinstance(source, Reseau) else source.copie()
    t = r.transport
    messages: list[str] = []
    if isinstance(source, Reseau) and source.instances:
        messages.append(
            f"{len(source.instances)} instance(s) de zone du reseau importe aplaties"
        )
    if abs(src.f0 - t.f0) > 1e-9:
        messages.append(
            f"f0 = {src.f0:g} Hz dans le reseau importe : il prend les {t.f0:g} Hz "
            "d'ici, reactances et pu conservees (materiel dimensionne pour ce reseau)"
        )
        src.f0 = t.f0

    # --- base de puissance : seuls les shunts sont en pu systeme ---
    if abs(src.Sb - t.Sb) > 1e-12 * max(t.Sb, 1.0):
        k = src.Sb / t.Sb
        src.T_Shunts = [(a, b, float(c) * k) for a, b, c in src.T_Shunts]
        src.T_ShuntsDC = [(a, b, float(c) * k) for a, b, c in src.T_ShuntsDC]
        messages.append(f"Sb {src.Sb:g} VA ramene a {t.Sb:g} VA (shunts convertis)")
        src.Sb = t.Sb

    # --- une seule reference : celle du reseau existant ---
    ib = {int(x[1]) for x in src.T_IBs}
    retires = set() if garder_ib else ib
    production = None if garder_slack else _production_des_slacks(src)
    if not garder_slack and production is None:
        messages.append(
            "le load flow du reseau importe seul ne converge pas : ses slacks "
            "passent en PV avec la production de leur table"
        )
    if retires:
        src.T_IBs = []
        messages.append(f"bus infini retire aux bus {sorted(retires)} (d'origine)")
    y = []
    for ligne in src.Y_DER:
        ligne = list(ligne)
        b = int(ligne[0])
        if b in retires:
            ligne[1] = PQ
        elif int(ligne[1]) == SLACK and not garder_slack:
            ligne[1] = PV
            if production is not None:
                ligne[5] = production[b]
            messages.append(
                f"slack {b} (d'origine) passe en PV a {ligne[5] / 1e6:.1f} MW"
            )
        y.append(tuple(ligne))
    src.Y_DER = y

    # --- numerotation ---
    deja = _numeros(aplatir(r))
    haut = max(deja, default=0)
    if prefixe is None:
        prefixe = 100 * (haut // 100 + 1)
    elif prefixe <= haut:
        raise ValueError(f"prefixe {prefixe} : les numeros jusqu'a {haut} sont pris")

    fusion = None
    mode = (raccord or {}).get("mode")
    if mode not in (None, "", "aucun", "fusion", "ligne", "transfo"):
        raise ValueError(f"raccordement {mode!r} : fusion, ligne ou transfo")
    if mode in ("fusion", "ligne", "transfo"):
        b_imp, b_ex = int(raccord["bus_importe"]), int(raccord["bus_existant"])
        un_imp = {int(x[0]): float(x[2]) for x in src.Y_DER}.get(b_imp)
        un_ex = {int(x[0]): float(x[2]) for x in t.Y_DER}.get(b_ex)
        if un_imp is None:
            raise ValueError(f"le reseau importe n'a pas de bus {b_imp}")
        if un_ex is None:
            raise ValueError(f"le transport n'a pas de bus {b_ex}")
        if mode == "fusion":
            if abs(un_imp - un_ex) > 1e-6 * un_ex:
                raise ValueError(
                    f"fusion impossible : {un_imp / 1e3:g} kV contre {un_ex / 1e3:g} kV "
                    "-- raccorder par un transformateur"
                )
            fusion = {b_imp: b_ex}

    nt = renumeroter(_replier_reglages(src), prefixe, 1, fusion)
    correspondance = {int(x[0]): _decale(x[0], prefixe, fusion) for x in src.Y_DER}
    correspondance |= {int(x[0]): int(x[0]) + prefixe for x in src.T_BusDC}

    out = t.copie()
    if fusion is not None:
        (b_ex,) = fusion.values()
        imp = next(x for x in nt.Y_DER if int(x[0]) == b_ex)
        nt.Y_DER = [x for x in nt.Y_DER if int(x[0]) != b_ex]
        rows = []
        for x in out.Y_DER:
            x = list(x)
            if int(x[0]) == b_ex:
                for c in (5, 6, 7, 8):
                    x[c] = float(x[c]) + float(imp[c])
                if int(x[1]) == PQ and int(imp[1]) == PV:
                    x[1], x[3] = PV, imp[3]  # une unite importee tient la tension
            rows.append(tuple(x))
        out.Y_DER = rows
        messages.append(f"bus {b_imp} (d'origine) fusionne avec le bus {b_ex}")

    for nom in TABLES:
        setattr(out, nom, [*getattr(out, nom), *getattr(nt, nom)])
    out.PARAMS_PAR_BUS = {**out.PARAMS_PAR_BUS, **nt.PARAMS_PAR_BUS}
    out.OUVERTS = [*out.OUVERTS, *nt.OUVERTS]

    if fusion is not None:
        porteurs = [
            int(x[1]) for nom in (*FAMILLES, "T_IBs") for x in getattr(out, nom)
        ]
        if porteurs.count(b_ex) > 1:
            raise ValueError(
                f"fusion impossible : le bus {b_ex} porterait deux unites "
                "-- raccorder par une ligne ou un transformateur"
            )

    if mode in ("ligne", "transfo"):
        b_new = correspondance[b_imp]
        table = "T_Lines_pu" if mode == "ligne" else "T_Transfos"
        ident = max((int(x[0]) for x in getattr(out, table)), default=0) + 1
        if mode == "ligne":
            ligne = (ident, b_ex, b_new, 1.0, 0.01, 0.1, 0.0, un_ex, t.Sb)
        else:
            ligne = (ident, b_ex, b_new, t.Sb, un_ex, un_imp, 0.0, 0.05)
        setattr(out, table, [*getattr(out, table), ligne])
        messages.append(
            f"{'ligne' if mode == 'ligne' else 'transformateur'} {ident} pose entre "
            f"le bus {b_ex} et le bus {b_new}"
        )

    # --- le dessin : l'existant a gauche, l'importe a droite ---
    positions = {b: (0.02 + 0.6 * p[0], p[1]) for b, p in r.positions.items()}
    for b, p in (source.positions if isinstance(source, Reseau) else {}).items():
        nb = correspondance.get(int(b))
        if nb is not None and nb not in positions:
            positions[nb] = (0.66 + 0.32 * p[0], p[1])

    nouveau = Reseau(
        transport=out,
        zones=r.zones,
        instances=r.instances,
        nom=r.nom,
        positions=positions,
        f_zoom=r.f_zoom,
    )
    if mode in (None, "", "aucun"):
        messages.append(
            "le reseau importe reste une ile sans reference tant qu'il n'est pas "
            "raccorde : tracer une ligne ou un transformateur vers le reseau existant"
        )
    aplatir(nouveau)  # ses garde-fous jouent avant de rendre la main
    return nouveau, {"prefixe": prefixe, "bus": correspondance, "messages": messages}
