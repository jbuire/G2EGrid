"""La validation d'un reseau edite et ses reparations proposees.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import numpy as np
from g2egrid.noyau.loadflow import PQ
from g2egrid.reseau.tables import CHARGES_SPECIALES, Tables
from g2egrid.reseau import disjoncteurs as disj
from ._commun import EQUIPEMENTS
from .assemblage import construire
from .catalogue_api import _iles


# ----------------------------------------------------------------- validation
# ------------------------------------------------------------- reparations
#: genre des unites qui tiennent une reference d'angle, par table
_FORMEURS = ("T_IBs", "T_SMs", "T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s")

_LIBELLES = {
    "slack": "Prendre {nom} comme slack",
    "un_slack": "Garder {nom} comme seul slack",
    "pv": "Passer le bus {bus} en PV",
    "pq": "Passer le bus {bus} en PQ",
    "raccorder": "Raccorder le bus {bus} au bus {vers}",
    "charge_p": "Charge du bus {bus} a puissance constante",
    "cdc": "Ajouter un condensateur au bus continu {bus}",
}


def _formeur_principal(t: Tables, ile) -> tuple[str, int] | None:
    """L'unite formatrice de plus grande puissance de l'ile : ``(table, bus)``."""
    ile = set(ile)
    cands = [
        (float(r[2]), table, int(r[1]))
        for table in _FORMEURS
        for r in getattr(t, table)
        if int(r[1]) in ile
    ]
    if not cands:
        return None
    _sn, table, b = max(cands)
    return table, b


def _plus_proche(t: Tables, b: int, positions: dict | None) -> int | None:
    """Le bus de meme tension le plus proche de ``b`` (sur le schema s'il est
    donne, par numero sinon), hors de l'ile de ``b``."""
    un = {int(r[0]): float(r[2]) for r in t.Y_DER}
    autres = [x for x in un if x != b and un[x] == un.get(b)] or [
        x for x in un if x != b
    ]
    if not autres:
        return None
    pos = {int(k): v for k, v in (positions or {}).items()}
    if b in pos and all(x in pos for x in autres):
        pb = np.asarray(pos[b], dtype=float)
        return min(autres, key=lambda x: float(np.linalg.norm(np.asarray(pos[x]) - pb)))
    return min(autres, key=lambda x: abs(x - b))


def _proposer(t: Tables, action: str, positions=None, **args) -> dict:
    """La reparation a joindre a un probleme, si elle est possible :
    ``{"reparation": {"action", "args", "libelle"}}`` (vide sinon)."""
    if action in ("slack", "un_slack"):
        f = _formeur_principal(t, args["ile"])
        if f is None:
            return {}
        args = {**args, "bus": f[1]}
        nom = f"{f[0][2:-1]} du bus {f[1]}"
    elif action == "raccorder":
        vers = _plus_proche(t, args["bus"], positions)
        if vers is None:
            return {}
        args = {**args, "vers": vers}
        nom = ""
    else:
        nom = ""
    return {
        "reparation": {
            "action": action,
            "args": args,
            "libelle": _LIBELLES[action].format(nom=nom, **args),
        }
    }


def _mettre_type(t: Tables, b: int, ty: int) -> None:
    t.Y_DER = [(r[0], ty, *r[2:]) if int(r[0]) == b else r for r in t.Y_DER]


def reparer(
    t: Tables, action: str, args: dict, positions: dict | None = None
) -> Tables:
    """Applique une reparation proposee par :func:`valider` ; rend une copie."""
    t = t.copie()
    b = int(args.get("bus", 0))
    if action == "slack":
        _mettre_type(t, b, 1)
    elif action == "un_slack":
        ile = {int(x) for x in args["ile"]}
        unites = {int(r[1]) for tb in _FORMEURS for r in getattr(t, tb)}
        for r in list(t.Y_DER):
            x = int(r[0])
            if x in ile and int(r[1]) == 1 and x != b:
                _mettre_type(t, x, 2 if x in unites else PQ)
        _mettre_type(t, b, 1)
    elif action in ("pv", "pq"):
        _mettre_type(t, b, 2 if action == "pv" else PQ)
    elif action == "raccorder":
        vers = int(args.get("vers") or _plus_proche(t, b, positions) or 0)
        if not vers:
            raise ValueError(f"bus {b} : aucun bus ou le raccorder")
        un = {int(r[0]): float(r[2]) for r in t.Y_DER}
        # une ligne courte, comme celles du reseau a cette tension s'il y en a
        modele = next(
            (r for r in t.T_Lines_pu if float(r[7]) == un[b]),
            None,
        )
        R, X, B = (modele[4], modele[5], modele[6]) if modele else (0.01, 0.1, 0.0)
        lid = max((int(r[0]) for r in t.T_Lines_pu), default=0) + 1
        t.T_Lines_pu = [*t.T_Lines_pu, (lid, b, vers, 1, R, X, B, un[b], t.Sb)]
    elif action == "charge_p":
        r = next(r for r in t.Y_DER if int(r[0]) == b)
        P, Q = float(r[7]), float(r[8])
        t.Y_DER = [(*x[:7], 0.0, 0.0) if int(x[0]) == b else x for x in t.Y_DER]
        lid = max((int(x[0]) for x in t.T_Loads), default=0) + 1
        t.T_Loads = [*t.T_Loads, (lid, b, P, Q, "P")]
    elif action == "cdc":
        sid = max((int(x[0]) for x in t.T_ShuntsDC), default=0) + 1
        t.T_ShuntsDC = [*t.T_ShuntsDC, (sid, b, float(args.get("C_pu", 0.01)))]
    else:
        raise ValueError(f"reparation inconnue : {action!r}")
    return t


def valider(t: Tables, positions: dict | None = None) -> list[dict]:
    """Ce qui empechera l'assemblage, ou le rendra faux, avant d'essayer.

    Le but n'est pas de dupliquer les controles de ``assembler`` (qui leve
    des ``ValueError`` explicites), mais de les rendre visibles PENDANT
    qu'on dessine, et d'y ajouter ce qu'un assembleur ne peut pas voir :
    un bus isole, l'absence de slack, deux equipements au meme endroit.

    Et un controle qui demande le point de fonctionnement : un convertisseur
    au-dela de sa limite de courant (ou de modulation) a l'equilibre, que le
    modele petit signal ne represente plus (:func:`_avis_limites`).
    """
    origine = t
    p: list[dict] = []
    isoles_voulus: set[int] = set()
    # Des disjoncteurs ouverts : on valide ce qui reste EN SERVICE (une ile
    # hors tension n'est pas un "bus isole", c'est l'effet voulu de
    # l'ouverture), et on dit ce que l'ouverture a fait.
    if t.OUVERTS:
        try:
            t, rap = disj.en_service(t)
        except ValueError as err:
            return [{"gravite": "erreur", "message": str(err)}]
        p += [{"gravite": "avis", "message": m} for m in rap.messages]
        # une ile d'un seul bus qui garde son formeur : une unite a vide,
        # isolee par l'ouverture -- voulu, donc un avis et non une erreur
        isoles_voulus = {i["bus"][0] for i in rap.iles if len(i["bus"]) == 1}
    ids = [int(r[0]) for r in t.Y_DER]
    if not ids:
        return [*p, {"gravite": "erreur", "message": "aucun bus"}]
    if len(set(ids)) != len(ids):
        p.append({"gravite": "erreur", "message": "deux bus portent le meme numero"})

    # --- connexite alternative, ile par ile ---------------------------
    voisins: dict[int, set[int]] = {b: set() for b in ids}
    for nom in ("T_Lines_SI", "T_Lines_pu", "T_Transfos"):
        for r in getattr(t, nom):
            a, b = int(r[1]), int(r[2])
            for x in (a, b):
                if x not in voisins:
                    p.append(
                        {
                            "gravite": "erreur",
                            "message": f"{nom} {r[0]} : bus {x} inconnu",
                        }
                    )
            if a in voisins and b in voisins:
                voisins[a].add(b)
                voisins[b].add(a)

    iles = _iles(ids, voisins)

    # UNE REFERENCE D'ANGLE PAR ILE ALTERNATIVE, ET NON UNE POUR TOUT LE
    # RESEAU. Une liaison a courant continu separe les iles alternatives :
    # elle transmet de la puissance mais aucun angle, et chaque cote garde
    # donc sa propre reference. Exiger un slack unique ferait declarer faux
    # un HVDC a deux terminaux, qui est le cas d'usage meme du continu ; ne
    # rien exiger laisserait passer une ile sans reference, dont le modele
    # d'etat n'aurait pas de repere.
    types = {int(r[0]): int(r[1]) for r in t.Y_DER}
    for ile in iles:
        n_slack = sum(1 for b in ile if types.get(b) == 1)
        ou = f"l'ile {{{', '.join(str(b) for b in sorted(ile))}}}"
        quoi = ou if len(iles) > 1 else "le reseau"
        if n_slack == 0:
            p.append(
                {
                    "gravite": "erreur",
                    "message": f"aucun bus slack (type 1) dans {quoi}",
                    **_proposer(t, "slack", ile=sorted(ile)),
                }
            )
        elif n_slack > 1:
            p.append(
                {
                    "gravite": "erreur",
                    "message": f"{n_slack} bus slack dans {quoi} : il en faut un seul",
                    **_proposer(t, "un_slack", ile=sorted(ile)),
                }
            )

    # UN BUS N'EST PAS ISOLE S'IL PORTE UN CONVERTISSEUR : il est relie par
    # le cote continu, que les tables de branches alternatives ne montrent
    # pas. Le signaler serait accuser un terminal HVDC d'etre en l'air.
    # Un bus isole est exactement une ile d'un seul element : ``_iles``
    # couvre tous les numeros, il n'y a donc rien a rattraper apres.
    sur_vsc = {int(r[1]) for r in t.T_VSCs} | {int(r[2]) for r in t.T_VSCs}
    for ile in iles:
        if len(ile) != 1:
            continue
        (b,) = ile
        if b in sur_vsc:
            continue
        if b in isoles_voulus:
            p.append(
                {
                    "gravite": "avis",
                    "message": f"bus {b} isole par les disjoncteurs : son unite "
                    "tourne a vide, dans son propre ilot",
                }
            )
            continue
        p.append(
            {
                "gravite": "erreur",
                "message": f"bus {b} isole : aucune branche n'y arrive",
                **_proposer(t, "raccorder", positions=positions, bus=b),
            }
        )

    # convention : un GFM est a batterie (il tient la tension, bus PV ou
    # slack), un GFL est un PV (il suit a sa puissance MPPT, bus PQ). Un
    # convertisseur raccorde au continu suit, lui, le mode de T_VSCs.
    genre_bus = {int(r[0]): int(r[1]) for r in t.Y_DER}
    raccordes = {int(r[1]) for r in t.T_VSCs}
    for table in ("T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s", "T_GFL1s", "T_GFL2s"):
        formeur = table.startswith("T_GFM")
        for r in getattr(t, table):
            b = int(r[1])
            if b in raccordes or b not in genre_bus:
                continue
            if formeur and genre_bus[b] == PQ:
                p.append(
                    {
                        "gravite": "avis",
                        "message": f"GFM du bus {b} en PQ : un GFM a batterie tient "
                        "sa tension, on le declare d'ordinaire en PV",
                        **_proposer(t, "pv", bus=b),
                    }
                )
            elif not formeur and genre_bus[b] != PQ:
                p.append(
                    {
                        "gravite": "avis",
                        "message": f"GFL du bus {b} en PV ou slack : un GFL (PV) suit "
                        "a sa puissance sans tenir la tension, on le declare "
                        "d'ordinaire en PQ",
                        **_proposer(t, "pq", bus=b),
                    }
                )

    # equipements
    occupe: dict[int, str] = {}
    for fam, spec in EQUIPEMENTS.items():
        if spec["table"] in ("T_Loads", "T_Shunts"):
            continue
        for r in getattr(t, spec["table"]):
            b = int(r[1])
            if b not in voisins:
                p.append(
                    {"gravite": "erreur", "message": f"{fam} {r[0]} : bus {b} inconnu"}
                )
            elif b in occupe:
                p.append(
                    {
                        "gravite": "erreur",
                        "message": f"bus {b} : {occupe[b]} et {fam} au meme bus",
                    }
                )
            else:
                occupe[b] = fam

    # charges R-L : une inductance ne consomme pas un reactif negatif
    for r in t.Y_DER:
        if (r[7] or r[8]) and r[8] <= 0:
            p.append(
                {
                    "gravite": "erreur",
                    "message": (
                        f"bus {int(r[0])} : Ql = {r[8]:g} var. Une charge R-L demande Ql > 0 ; "
                        "declarer le condensateur en gradin, ou la charge en type I ou P"
                    ),
                    **_proposer(t, "charge_p", bus=int(r[0])),
                }
            )
    for r in t.T_Loads:
        if str(r[4]) not in CHARGES_SPECIALES:
            p.append(
                {
                    "gravite": "erreur",
                    "message": f"charge {r[0]} : type '{r[4]}' inconnu ('I' ou 'P')",
                }
            )

    # avertissements
    for b in ids:
        if types.get(b) == 2 and b not in occupe:
            p.append(
                {
                    "gravite": "avis",
                    "message": f"bus {b} declare PV mais sans equipement",
                    **_proposer(t, "pq", bus=b),
                }
            )
    if not occupe:
        p.append(
            {
                "gravite": "avis",
                "message": "aucun equipement : le reseau n'aura aucune dynamique",
            }
        )

    # --- raccordements continus ---------------------------------------
    # L'assembleur refuse deja tout ceci, avec de bons messages ; mais il ne
    # parle qu'au moment ou l'on demande un calcul. Les repeter ici, c'est
    # les montrer PENDANT qu'on dessine, au lieu de laisser decouvrir par un
    # load flow en echec qu'un trait pose trois gestes plus tot etait faux.
    bus_dc = {int(r[0]) for r in t.T_BusDC}
    porteurs = {
        int(r[1])
        for table in ("T_GFM3s", "T_GFM4s", "T_GFL2s")
        for r in getattr(t, table)
    }
    vus_ac: set[int] = set()
    for r in t.T_VSCs:
        vid, ac, dc = r[0], int(r[1]), int(r[2])
        if dc not in bus_dc:
            p.append(
                {
                    "gravite": "erreur",
                    "message": f"VSC {vid} : le bus continu {dc} n'existe pas",
                }
            )
        if ac not in porteurs:
            p.append(
                {
                    "gravite": "erreur",
                    "message": (
                        f"VSC {vid} : le bus {ac} ne porte ni GFM3, ni GFM4, ni "
                        "GFL2 -- seules ces familles ont une variante a port continu"
                    ),
                }
            )
        if ac in vus_ac:
            p.append(
                {
                    "gravite": "erreur",
                    "message": f"VSC {vid} : le bus {ac} est deja raccorde",
                }
            )
        vus_ac.add(ac)
        # UN VSC NE PEUT PAS ETRE SUR LE SLACK : il tient exactement deux
        # grandeurs, et sur le bus de reference les deux puissances sont
        # libres. C'est la regle qui surprend le plus en dessinant, parce
        # que le premier equipement pose devient justement la reference.
        # ...SAUF s'il FORME son ile (mode Vtheta) : il en est alors le slack,
        # et ses deux puissances sortent du bilan de l'ile -- la station d'un
        # parc offshore (CIGRE B4)
        if types.get(ac) == 1 and str(r[3]) != "Vtheta":
            p.append(
                {
                    "gravite": "erreur",
                    "message": (
                        f"VSC {vid} : le bus {ac} est le slack. Un convertisseur "
                        "tient exactement deux grandeurs, et sur la reference "
                        "d'angle les deux puissances sont libres -- il faut une "
                        "autre source pour tenir cette ile."
                    ),
                }
            )
    for b in bus_dc:
        relie = any(int(r[2]) == b for r in t.T_VSCs)
        relie = relie or any(b in (int(r[1]), int(r[2])) for r in t.T_LinesDC)
        if not relie:
            p.append(
                {
                    "gravite": "erreur",
                    "message": f"bus continu {b} isole : ni cable ni convertisseur",
                }
            )
    # un noeud continu sans capacite : l'assembleur le refuse (equation
    # algebrique sans variable), autant le dire en dessinant
    try:
        from g2egrid.reseau.tables import _capacites_dc, _vsc_par_bus

        un_dc = {int(r[0]): float(r[2]) for r in t.T_BusDC}
        cdc = _capacites_dc(t, _vsc_par_bus(t), un_dc, 2 * np.pi * t.f0)
    except Exception:  # noqa: BLE001 -- d'autres erreurs le disent deja
        cdc = {}
    for b, c in sorted(cdc.items()):
        if c <= 0:
            p.append(
                {
                    "gravite": "erreur",
                    "message": f"bus continu {b} sans capacite : ni cable, ni "
                    "condensateur, ni convertisseur qui en porte un",
                    **_proposer(t, "cdc", bus=b),
                }
            )
    if not any(x["gravite"] == "erreur" for x in p):
        p += _avis_convergence(origine)
    if not any(x["gravite"] == "erreur" for x in p):
        p += _avis_limites(origine)
    return p


def _avis_convergence(t: Tables) -> list[dict]:
    """Une ERREUR si le load flow du reseau ne converge pas.

    Sans point de fonctionnement, ni les modes, ni le temporel, ni aucun
    calcul de l'atelier ne peut partir : ils refusent tous avec ce meme
    message. La validation est l'endroit ou le dire, PENDANT le dessin, plutot
    que de laisser un reseau « valide » (``ok: true``) que tout calcul refuse.
    C'est une erreur et non un avis pour la meme raison : ``ok`` doit dire si
    on peut calculer. Le calcul est celui que les autres routes font ensuite,
    et son cache les sert. Seule la non-convergence est signalee ici : les
    autres refus de l'assemblage ont leurs propres messages.
    """
    try:
        construire(t)
    except RuntimeError as err:
        return [
            {
                "gravite": "erreur",
                "message": (
                    f"{err} : aucun calcul n'est "
                    "possible. Verifiez les charges et productions, les "
                    "impedances de ligne et la tension des bus de controle."
                ),
            }
        ]
    except Exception:  # noqa: BLE001 -- les autres refus ont leurs messages
        return []
    return []


#: les tables dont les equipements portent un limiteur de courant (GFM, GFL,
#: et, par ``T_VSCs``, leurs variantes a port continu et les MMC)
_TABLES_LIMITEES = (
    "T_GFM1s",
    "T_GFM2s",
    "T_GFM3s",
    "T_GFM4s",
    "T_GFL1s",
    "T_GFL2s",
)


def _avis_limites(t: Tables) -> list[dict]:
    """Un avis par convertisseur au-dela de sa limite au point d'equilibre.

    Le seul controle de :func:`valider` qui calcule : load flow et equilibre,
    par :func:`construire`, dont le cache les garde pour les calculs qui
    suivent (modes, temporel...). Il n'est fait que sur un reseau sans erreur
    et qui porte des convertisseurs. Un reseau qui n'a pas encore de point
    de fonctionnement (load flow qui ne converge pas, reseau en cours de
    dessin) n'a pas d'avis ici : c'est le calcul, quand on le lance, qui
    dira pourquoi. Chaque avis porte, sous ``depassement``, le detail chiffre
    (:meth:`g2egrid.noyau.limites.Depassement.en_dict`).
    """
    if not any(getattr(t, n, None) for n in _TABLES_LIMITEES):
        return []
    try:
        net = construire(t).equilibre()
    except Exception:  # noqa: BLE001 -- le calcul, lance, dira pourquoi
        return []
    return [
        {
            "gravite": "avis",
            "message": (
                f"{d} : point de fonctionnement au-dela de la limite ; le "
                "modele petit signal (mode libre) ne represente pas ce "
                "convertisseur, qui serait en limitation"
            ),
            "depassement": d.en_dict(),
        }
        for d in net.depassements
    ]
