"""Disjoncteurs et ilotage : le reseau EN SERVICE, pour les calculs de regime permanent.

Un disjoncteur ouvert se declare dans ``Tables.OUVERTS``, une liste de couples
``(table, ID)`` -- ``("T_Lines_pu", 3)`` ouvre la ligne 3, ``("T_SMs", 2)``
declenche la machine 2, ``("Y_DER", 5)`` deleste la charge du bus 5 (une
charge de ``Y_DER`` n'a pas d'ID propre : c'est son bus). Les tables elles-memes
ne bougent pas : on peut refermer, et le reseau d'origine est intact.

:func:`en_service` rend les tables que les calculs doivent voir, en trois temps.

1. **Retirer ce qui est ouvert.** Une unite declenchee laisse son bus sans
   production : ``Pg = Qg = 0``, et un bus PV ou slack redevient PQ.
2. **Mettre hors tension les iles mortes.** Les branches restantes decoupent
   le reseau alternatif en iles. Une ile qui n'a plus aucun FORMEUR de reseau
   -- machine synchrone, grid-forming, bus infini -- n'a rien pour tenir sa
   tension et sa frequence : un grid-following ne fait que suivre une tension
   (et sa protection anti-ilotage le declencherait). Elle est retiree, bus,
   charges et equipements compris.
3. **Une reference par ile vivante.** Chaque ile garde son slack s'il est
   encore en service ; sinon, son plus gros formeur le devient. C'est ce que
   fait tout outil de load flow, et c'est deja ce que le paquet fait pour les
   deux cotes d'une liaison a courant continu. Dans le modele dynamique,
   chaque ile de plus apporte un mode de jauge de plus (son angle absolu est
   libre) : il est filtre comme les autres.

Pour les evenements EN COURS de simulation (ouverture a t = 0,1 s), voir
:class:`g2egrid.temporel.Manoeuvre` : la, rien n'est recalcule par un load
flow -- les elements gardent leurs consignes, et l'ile qui se separe evolue
d'elle-meme.

Idee reprise du portage de F. Kelada (``network/breakers.py``), reecrite sur
les tables de ce paquet.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from g2egrid.reseau.tables import TABLES, Tables

__all__ = [
    "BRANCHES",
    "CABLES_DC",
    "FORMEURS",
    "OUVRABLES",
    "UNITES",
    "Rapport",
    "basculer",
    "en_service",
    "est_ouvert",
    "nom_element",
    "ouvrables",
]

#: Les branches alternatives : ce sont elles qui decoupent le reseau en iles.
BRANCHES = ("T_Lines_SI", "T_Lines_pu", "T_Transfos")
#: Les unites qui TIENNENT une tension et une frequence.
FORMEURS = ("T_IBs", "T_SMs", "T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s")
#: Toutes les unites (les grid-following suivent, ne forment pas).
UNITES = (*FORMEURS, "T_GFL1s", "T_GFL2s")
#: Les cables continus : leur ouverture ne decoupe pas le reseau alternatif
#: mais peut decouper le reseau continu en iles (voir ``Rapport.iles_dc``).
CABLES_DC = ("T_LinesDC",)
#: Ce qu'on peut ouvrir, et comment le nommer.
OUVRABLES = (*BRANCHES, *UNITES, "Y_DER", "T_Loads", "T_Shunts", *CABLES_DC)

#: Les colonnes qui designent un bus, table par table (celles de
#: ``zones.BUS_AC``, redites ici pour ne pas dependre des zones).
_COLONNES_BUS = {
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

SLACK, PV, PQ = 1, 2, 3


@dataclass
class Rapport:
    """Ce que :func:`en_service` a fait, pour le dire a l'utilisateur."""

    ouverts: list[tuple[str, int]] = field(default_factory=list)
    #: bus retires parce que leur ile n'a plus de formeur
    hors_tension: list[int] = field(default_factory=list)
    #: une entree par ile vivante : ses bus et son bus de reference
    iles: list[dict] = field(default_factory=list)
    messages: list[str] = field(default_factory=list)
    #: une entree par ile CONTINUE : ses noeuds, ses convertisseurs et ses
    #: references de tension (vide si le reseau n'a pas de continu)
    iles_dc: list[dict] = field(default_factory=list)


def _ligne(t: Tables, table: str, ident: int):
    for r in getattr(t, table):
        if int(r[0]) == int(ident):
            return r
    return None


def est_ouvert(t: Tables, table: str, ident: int) -> bool:
    return (table, int(ident)) in {(a, int(b)) for a, b in t.OUVERTS}


def basculer(t: Tables, table: str, ident: int, ouvert: bool) -> Tables:
    """Une copie des tables avec ce disjoncteur ouvert (ou referme)."""
    if table not in OUVRABLES:
        raise ValueError(f"{table} : pas un element manoeuvrable ; {OUVRABLES}")
    if _ligne(t, table, ident) is None:
        raise ValueError(f"{table} : pas d'element d'ID {ident}")
    out = t.copie()
    reste = [(a, int(b)) for a, b in t.OUVERTS if (a, int(b)) != (table, int(ident))]
    out.OUVERTS = [*reste, (table, int(ident))] if ouvert else reste
    return out


def nom_element(t: Tables, table: str, ident: int) -> str:
    """Le nom de l'element ASSEMBLE, tel que le modele dynamique le connait.

    C'est celui que :class:`~g2egrid.temporel.Manoeuvre` vise. Les noms sont ceux
    de :func:`g2egrid.reseau.tables.assembler` ; un condensateur de ``T_Shunts`` n'y
    est pas un element (il est dans la capacite du noeud) et se nomme
    ``Csh<ID>_<bus>``, comme la reactance ``Xsh<ID>_<bus>``.
    """
    r = _ligne(t, table, ident)
    if r is None:
        raise ValueError(f"{table} : pas d'element d'ID {ident}")
    if table in ("T_Lines_SI", "T_Lines_pu"):
        return f"L{int(r[0])}_{int(r[1])}_{int(r[2])}"
    if table == "T_Transfos":
        return f"T{int(r[0])}_{int(r[1])}_{int(r[2])}"
    if table == "T_IBs":
        return f"IB{int(r[1])}"
    if table == "T_SMs":
        return f"SM{int(r[1])}"
    if table.startswith("T_GFM"):
        return f"GFM{int(r[1])}"
    if table.startswith("T_GFL"):
        return f"GFL{int(r[1])}"
    if table == "Y_DER":
        return f"Load{int(r[0])}"
    if table == "T_Loads":
        return f"Load{r[4]}{int(r[0])}_{int(r[1])}"
    if table == "T_Shunts":
        genre = "Xsh" if float(r[2]) < 0 else "Csh"
        return f"{genre}{int(r[0])}_{int(r[1])}"
    if table == "T_LinesDC":
        return f"CDC{int(r[0])}_{int(r[1])}_{int(r[2])}"
    raise ValueError(f"{table} : pas un element manoeuvrable")


def ouvrables(t: Tables) -> list[dict]:
    """Tous les elements qu'on peut ouvrir, avec leur nom et leur etat."""
    out = []
    for table in OUVRABLES:
        for r in getattr(t, table):
            if table == "Y_DER" and not (float(r[7]) or float(r[8])):
                continue  # un bus sans charge n'a rien a delester
            ident = int(r[0])
            out.append(
                {
                    "table": table,
                    "id": ident,
                    "nom": nom_element(t, table, ident),
                    "ouvert": est_ouvert(t, table, ident),
                    "genre": _genre(table),
                }
            )
    return out


def _genre(table: str) -> str:
    if table in BRANCHES:
        return "branche"
    if table in UNITES:
        return "unite"
    if table == "T_Shunts":
        return "shunt"
    if table in CABLES_DC:
        return "cable"
    return "charge"


def _iles(bus: list[int], liens: list[tuple[int, int]]) -> list[list[int]]:
    voisins: dict[int, set[int]] = {b: set() for b in bus}
    for a, b in liens:
        if a in voisins and b in voisins:
            voisins[a].add(b)
            voisins[b].add(a)
    vus: set[int] = set()
    iles = []
    for b in bus:
        if b in vus:
            continue
        pile, ile = [b], []
        vus.add(b)
        while pile:
            c = pile.pop()
            ile.append(c)
            for d in voisins[c]:
                if d not in vus:
                    vus.add(d)
                    pile.append(d)
        iles.append(sorted(ile))
    return iles


def _nom_station(t: Tables, bus_ac: int) -> str:
    """Le nom assemble du convertisseur du bus alternatif ``bus_ac``."""
    for table in ("T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s", "T_GFL1s", "T_GFL2s"):
        for r in getattr(t, table):
            if int(r[1]) == bus_ac:
                return nom_element(t, table, int(r[0]))
    return f"VSC@{bus_ac}"


def _iles_continues(out: Tables) -> list[list[int]]:
    bus = [int(r[0]) for r in out.T_BusDC]
    liens = [(int(r[1]), int(r[2])) for r in out.T_LinesDC]
    return _iles(bus, liens)


def _references_dc(out: Tables, ile: list[int]) -> list[str]:
    """Ce qui tient la tension d'une ile continue : un VSC en mode ``Vdc...``
    ou a statisme, une batterie (DC/DC a statisme)."""
    refs = []
    for r in out.T_VSCs:
        if int(r[2]) in ile and (
            str(r[3]).startswith("Vdc") or str(r[3]) in ("droop", "droopI")
        ):
            refs.append(f"{_nom_station(out, int(r[1]))} ({r[3]})")
    for r in out.T_DCDC:
        if int(r[1]) in ile and str(r[3]).lower() == "batterie":
            refs.append(f"DCDC{int(r[0])}_{int(r[1])} (batterie)")
    return refs


def _continu(
    t: Tables,
    out: Tables,
    perdus: list[tuple],
    mortes: set[int],
    rap: Rapport,
) -> None:
    """Le cote continu apres les manoeuvres : orphelins, references, iles.

    * un noeud continu dont le SEUL convertisseur etait dans une ile
      alternative morte est orphelin : il est retire avec ses shunts et ses
      cables (les cables pendants, dont l'autre extremite reste, aussi) ;
    * un noeud de type 1 dont le convertisseur ``Vdc...`` a disparu (declenche
      ou orphelin) n'a plus de tenant : il passe au type 2 (passif), faute de
      quoi il resterait un slack ideal que le modele dynamique ne realise pas ;
    * les iles continues sont rapportees ; une ile qui contient encore des
      convertisseurs mais aucune reference de tension est REFUSEE ;
    * chaque station perdue est dite : ce qu'elle tenait, ce qui reste.
    """
    if not t.T_BusDC:
        return
    avant = len(_iles_continues(t))
    vivants = {int(r[0]) for r in out.T_BusDC}

    def occupe(b: int) -> bool:
        return (
            any(int(r[2]) == b for r in out.T_VSCs)
            or any(int(r[1]) == b for r in out.T_DCDC)
            or any(int(r[1]) == b for r in out.T_LoadsDC)
        )

    # --- 1. les noeuds orphelins d'une ile alternative morte ---
    orphelins = {
        int(r[2]) for r, _c in perdus if int(r[1]) in mortes and not occupe(int(r[2]))
    }
    a_voir = set(orphelins)
    while a_voir:
        b = a_voir.pop()
        pendants = [r for r in out.T_LinesDC if b in (int(r[1]), int(r[2]))]
        for r in pendants:
            autre = int(r[2]) if int(r[1]) == b else int(r[1])
            rap.messages.append(
                f"cable continu {int(r[0])} ({int(r[1])}-{int(r[2])}) retire : "
                f"le noeud continu {b} est orphelin (sa station est dans une "
                "ile alternative hors tension)"
            )
            if autre in vivants and autre not in orphelins:
                reste = [
                    x
                    for x in out.T_LinesDC
                    if x not in pendants and autre in (int(x[1]), int(x[2]))
                ]
                if not reste and not occupe(autre):
                    orphelins.add(autre)
                    a_voir.add(autre)
        out.T_LinesDC = [r for r in out.T_LinesDC if r not in pendants]
        vivants.discard(b)
        rap.messages.append(f"noeud continu {b} retire (orphelin)")
    if orphelins:
        out.T_BusDC = [r for r in out.T_BusDC if int(r[0]) not in orphelins]
        out.T_ShuntsDC = [r for r in out.T_ShuntsDC if int(r[1]) not in orphelins]

    # --- 2. un noeud de type 1 qui perd son tenant devient passif ---
    for r, _cause in perdus:
        b = int(r[2])
        if b in orphelins or b not in vivants or not str(r[3]).startswith("Vdc"):
            continue
        autres = [
            x for x in out.T_VSCs if int(x[2]) == b and str(x[3]).startswith("Vdc")
        ]
        autres += [
            x for x in out.T_DCDC if int(x[1]) == b and str(x[3]).lower() == "batterie"
        ]
        if autres:
            continue
        lignes = []
        for x in out.T_BusDC:
            if int(x[0]) == b and int(x[1]) == 1:
                x = (x[0], 2, *x[2:])
                rap.messages.append(
                    f"noeud continu {b} : plus aucun tenant de tension, il passe "
                    "de type 1 (slack) a type 2 (passif)"
                )
            lignes.append(x)
        out.T_BusDC = lignes

    # --- 3. les iles continues, et leurs references ---
    iles = _iles_continues(out)
    for ile in iles:
        conv = [_nom_station(t, int(r[1])) for r in out.T_VSCs if int(r[2]) in ile]
        refs = _references_dc(out, ile)
        rap.iles_dc.append({"bus": ile, "convertisseurs": conv, "references": refs})
        if conv and not refs:
            raise ValueError(
                f"ile continue {ile} sans reference de tension : ses "
                f"convertisseur(s) {', '.join(conv)} sont tous a puissance "
                "imposee. Un convertisseur en mode Vdc... (ou a statisme), une "
                "batterie ou une source doit tenir la tension de chaque ile "
                "continue."
            )
    if len(iles) > avant:
        rap.messages.append(
            f"reseau continu coupe en {len(iles)} iles : "
            + " ; ".join(str(i["bus"]) for i in rap.iles_dc)
        )

    # --- 4. ce que les stations perdues tenaient ---
    for r, cause in perdus:
        b_ac, b_dc, mode = int(r[1]), int(r[2]), str(r[3])
        nom = _nom_station(t, b_ac)
        if mode.startswith("Vdc"):
            role = f"tenait la tension du noeud continu {b_dc}"
        elif mode in ("droop", "droopI"):
            role = f"participait au reglage de la tension (statisme, noeud {b_dc})"
        else:
            role = f"imposait sa puissance au noeud continu {b_dc}"
        msg = f"station {nom} (bus alternatif {b_ac}, mode {mode}) perdue ({cause}) : elle {role}."
        ile = next((i for i in rap.iles_dc if b_dc in i["bus"]), None)
        if ile is not None and ile["convertisseurs"]:
            msg += (
                f" L'ile continue {ile['bus']} n'est plus tenue que par "
                f"{', '.join(ile['references'])}, qui reprend sa puissance."
                " A verifier : le point de fonctionnement change et le spectre"
                " peut devenir instable (effet physique du convertisseur"
                " restant, pas un artefact d'assemblage)."
            )
        rap.messages.append(msg)


def en_service(t: Tables) -> tuple[Tables, Rapport]:
    """Les tables du reseau en service, et le rapport de ce qui a ete fait."""
    rap = Rapport(ouverts=[(a, int(b)) for a, b in t.OUVERTS])
    out = t.copie()
    out.OUVERTS = []
    if not t.OUVERTS:
        return out, rap

    # --- 1. retirer ce qui est ouvert ---
    a_retirer: dict[str, set[int]] = {}
    for table, ident in rap.ouverts:
        if table not in OUVRABLES:
            raise ValueError(f"OUVERTS : {table} n'est pas un element manoeuvrable")
        if _ligne(t, table, ident) is None:
            raise ValueError(f"OUVERTS : {table} n'a pas d'element d'ID {ident}")
        a_retirer.setdefault(table, set()).add(int(ident))

    sans_prod: set[int] = set()
    for table in UNITES:
        for r in getattr(t, table):
            if int(r[0]) in a_retirer.get(table, ()):
                sans_prod.add(int(r[1]))
    delestes = a_retirer.get("Y_DER", set())
    ydr = []
    for r in out.Y_DER:
        r = list(r)
        b = int(r[0])
        if b in sans_prod:
            r[5] = r[6] = 0.0
            if int(r[1]) in (SLACK, PV):
                r[1] = PQ
        if b in delestes:
            r[7] = r[8] = 0.0
        ydr.append(tuple(r))
    out.Y_DER = ydr
    for table, ids in a_retirer.items():
        if table != "Y_DER":
            setattr(
                out, table, [r for r in getattr(out, table) if int(r[0]) not in ids]
            )
    # un convertisseur declenche emporte sa declaration de port continu
    out.T_VSCs = [r for r in out.T_VSCs if int(r[1]) not in sans_prod]

    # --- 2. les iles, et celles qui n'ont plus de formeur ---
    bus = [int(r[0]) for r in out.Y_DER]
    liens = [(int(r[1]), int(r[2])) for table in BRANCHES for r in getattr(out, table)]
    formeurs: dict[int, list[tuple[float, int]]] = {}
    for table in FORMEURS:
        for r in getattr(out, table):
            formeurs.setdefault(int(r[1]), []).append((float(r[2]), int(r[1])))
    mortes: set[int] = set()
    for ile in _iles(bus, liens):
        candidats = [f for b in ile for f in formeurs.get(b, [])]
        if not candidats:
            mortes |= set(ile)
            rap.messages.append(
                f"ile {ile} hors tension : aucune machine, aucun grid-forming, "
                "aucun bus infini pour tenir sa tension"
            )
            continue
        slacks = [b for b in ile if int(out.bus(b)[1]) == SLACK]
        if slacks:
            ref = slacks[0]
        else:
            ref = max(candidats)[1]  # le plus gros formeur de l'ile
            out.Y_DER = [
                (r[0], SLACK, *r[2:]) if int(r[0]) == ref else r for r in out.Y_DER
            ]
            rap.messages.append(
                f"ile {ile} : le bus {ref} (plus gros formeur) en devient le slack"
            )
        rap.iles.append({"bus": ile, "reference": ref})
    if mortes:
        rap.hors_tension = sorted(mortes)
        for table, cols in _COLONNES_BUS.items():
            setattr(
                out,
                table,
                [
                    r
                    for r in getattr(out, table)
                    if not any(int(r[c]) in mortes for c in cols)
                ],
            )
        out.PARAMS_PAR_BUS = {
            b: v for b, v in out.PARAMS_PAR_BUS.items() if int(b) not in mortes
        }
    perdus = [
        (r, "declenchee" if int(r[1]) in sans_prod else "ile alternative hors tension")
        for r in t.T_VSCs
        if int(r[1]) in sans_prod or int(r[1]) in mortes
    ]
    _continu(t, out, perdus, mortes, rap)
    if not rap.iles:
        raise ValueError(
            "plus aucune ile sous tension : tous les formeurs de reseau sont "
            "declenches ou isoles"
        )
    return out, rap


# ``TABLES`` est importe pour que ``_COLONNES_BUS`` reste verifiable : chaque
# table qui y figure doit exister dans le contrat des colonnes.
assert set(_COLONNES_BUS) <= set(TABLES)
