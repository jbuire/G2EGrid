"""Les resultats d'un point de fonctionnement : topologie, load flow, modes, detail d'un mode, signaux, equations, fiches.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import numpy as np

from g2egrid.analyse import modal
from g2egrid.analyse.sensibilites import sensibilites as sensibilites_parametriques
from g2egrid.multi_niveau.separation import reponse_indicielle as reponse_lineaire
from g2egrid.reseau import disjoncteurs as disj
from g2egrid.reseau.flux import flux_branches
from g2egrid.reseau.network import Network
from g2egrid.reseau.tables import TABLES, Tables, assembler
from g2egrid.temporel import (
    GENRES_MANOEUVRE,
    Evenement,
    Manoeuvre,
    ModeleNonLineaire,
    Trajectoire,
)
from g2egrid.temporel import mesures as mesures_physiques

from ._commun import EQUIPEMENTS, JAUGE, Assemble
from .assemblage import construire
from .catalogue_api import _positions


# --------------------------------------------------------------- resultats
def topologie(t: Tables, positions: dict | None = None) -> dict:
    """Bus, branches et equipements, pour le schema."""
    pos = positions or _positions(t)
    equip: dict[int, list[dict]] = {}
    for fam, spec in EQUIPEMENTS.items():
        for r in getattr(t, spec["table"]):
            b = int(r[1])
            equip.setdefault(b, []).append(
                {
                    "famille": fam,
                    "id": r[0],
                    "couleur": spec["couleur"],
                    "ligne": list(r),
                }
            )
    bus = []
    for r in t.Y_DER:
        b = int(r[0])
        bus.append(
            {
                "id": b,
                "type": int(r[1]),
                "Un": float(r[2]),
                "V": float(r[3]),
                "Pg": float(r[5]),
                "Pl": float(r[7]),
                "Ql": float(r[8]),
                "equipements": equip.get(b, []),
                "x": pos.get(str(b), pos.get(b, (0.5, 0.5)))[0],
                "y": pos.get(str(b), pos.get(b, (0.5, 0.5)))[1],
            }
        )
    branches = []
    for nom, genre in (
        ("T_Lines_SI", "ligne"),
        ("T_Lines_pu", "ligne"),
        ("T_Transfos", "transfo"),
    ):
        for r in getattr(t, nom):
            branches.append(
                {
                    "table": nom,
                    "genre": genre,
                    "id": r[0],
                    "de": int(r[1]),
                    "vers": int(r[2]),
                    "ligne": list(r),
                }
            )
    return {
        "bus": bus,
        "branches": branches,
        "positions": {str(k): list(v) for k, v in pos.items()},
    }


def loadflow(t: Tables) -> dict:
    a = construire(t)
    lf, net = a.lf, a.net
    idx = net._bus_index()
    ydr = {int(r[0]): r for r in t.Y_DER}
    bus = []
    for b, i in sorted(idx.items()):
        r = ydr.get(int(b))
        bus.append(
            {
                "id": b,
                "type": int(r[1]) if r is not None else 3,
                "V": float(abs(lf.V[i])),
                "angle": float(np.degrees(np.angle(lf.V[i]))),
                "P": float(lf.S[i].real),
                "Q": float(lf.S[i].imag),
                # la charge du bus, pour que l'onglet puisse rendre la
                # PRODUCTION (P net + charge) et non la seule injection nette
                "Pl": float(r[7]) / t.Sb if r is not None else 0.0,
                "Ql": float(r[8]) / t.Sb if r is not None else 0.0,
            }
        )
    # les transits, pour les dessiner sur le schema : ce qu'on lit d'abord sur
    # un plan de reseau, ce n'est pas une tension, c'est ce qui passe ou
    transits = [
        {
            "nom": tr.nom,
            "de": tr.de,
            "vers": tr.vers,
            "P_de": tr.P_de,
            "Q_de": tr.Q_de,
            "P_vers": tr.P_vers,
            "Q_vers": tr.Q_vers,
            "pertes": tr.pertes,
        }
        for tr in flux_branches(net, lf)
    ]
    return {
        "branches": transits,
        "converge": bool(lf.converged),
        "iterations": int(lf.iterations),
        "mismatch": float(lf.mismatch),
        "bus": bus,
        # les noeuds continus : le point de fonctionnement d'un statisme Vdc-P
        # se lit sur sa caracteristique (inspecteur de l'atelier)
        "bus_dc": [
            {"id": int(b), "V": float(v)}
            for b, v in zip(sorted(net.bus_specs_dc), getattr(lf, "Vdc", []))
        ],
        "pertes_P": float(lf.S.real.sum()),
        "pertes_Q": float(lf.S.imag.sum()),
    }


def _domaines_modes(net, ss):
    """Le DOMAINE de chaque mode (alternatif, continu, couple), des qu'il y a
    un reseau continu : c'est la premiere question qu'on pose a un mode d'un
    reseau AC/DC (:mod:`g2egrid.scenarios.domaines`). Rend ``(champs, bilan)`` :
    ``champs(indice)`` donne les champs a ajouter au mode d'indice ``indice``
    (vide sans continu), ``bilan`` le nombre de modes par domaine."""
    if not net.bus_specs_dc:
        return (lambda _k: {}), None
    from g2egrid.scenarios import domaines

    classes = domaines.classer_modes(net, ss)
    dom = {m.indice: m for m in classes}

    def champs(k):
        d = dom.get(int(k))
        if d is None:
            return {}
        return {
            "domaine": d.domaine,
            "part_dc": 100.0 * d.part_dc,
            "inter_iles": d.inter_iles,
            "iles": [{"nom": n, "part": 100.0 * x} for n, x in d.iles[:3] if x >= 0.01],
        }

    return champs, domaines.bilan(classes)


def modes(t: Tables, f_max: float | None = None, n_etats: int = 3) -> dict:
    """Les modes, et la marge de stabilite lue sur les modes PHYSIQUES.

    LE MODE DE JAUGE N'EST PAS UN MODE, et il fausse justement le nombre
    qu'on lit en premier. Un reseau dont la reference d'angle n'est tenue
    par aucun bus infini porte un zero exact, rendu a 1e-12 pres par
    l'algebre ; il y en a un par ile dans ce cas, et un reseau a liaison
    continue en compte plusieurs -- trois sur ``hvdc2t``, cinq sur
    ``mtdc3t``. Prendre le maximum de ``Re`` sur tout le spectre affichait
    donc "max Re = +0,00e+00" sur des reseaux dont la marge reelle est
    -14,2 et -0,880 : impossible d'y distinguer un reseau juste stable d'un
    reseau confortable.

    Les modes de jauge sont donc ecartes du classement et de la marge, mais
    leur nombre est RENDU plutot que taire : il compte, puisqu'il dit
    combien d'iles flottent sans reference absolue.
    """
    a = construire(t)
    ss = a.etat_stationnaire()
    liste = [
        m
        for m in modal.modes(ss, f_max=f_max, n_etats=n_etats)
        if abs(m.valeur) > JAUGE
    ]
    lam = ss.eig()
    physiques = lam[np.abs(lam) > JAUGE]
    domaine, bilan_dom = _domaines_modes(a.net, ss)

    return {
        "domaines": bilan_dom,
        "n_etats": len(ss.states),
        "etats": list(ss.states),
        "max_re": float(physiques.real.max()) if physiques.size else 0.0,
        "jauge": int((np.abs(lam) <= JAUGE).sum()),
        "stable": bool(physiques.size and physiques.real.max() < 0),
        "modes": [
            {
                "indice": m.indice,
                "re": m.valeur.real,
                "im": m.valeur.imag,
                "f": m.f,
                "zeta": m.zeta,
                "etats": [{"nom": n, "part": p} for n, p in m.etats],
                **domaine(m.indice),
            }
            for m in liste
        ],
    }


def impact_continu(t: Tables, n: int = 12) -> dict:
    """Ce que la dynamique du reseau continu fait aux modes : le reseau
    compare a son double ou le continu est rigide (tensions figees, cables
    instantanes ; :func:`g2egrid.scenarios.domaines.impact_continu`). Rend les ``n``
    modes retrouves les plus deplaces, et les modes retires (ceux du continu
    lui-meme), perdus et nouveaux."""
    from g2egrid.scenarios import domaines

    c = domaines.impact_continu(t)

    def ligne(m):
        return {
            "f": m.f,
            "zeta": m.zeta,
            "re": float(m.lam.real) if m.lam is not None else None,
            "f_rigide": m.f_reduit,
            "zeta_rigide": m.zeta_reduit,
            "re_rigide": float(m.lam_reduit.real) if m.lam_reduit is not None else None,
            "ecart": m.ecart,
            "dzeta": m.dzeta,
            "part_dc": 100.0 * m.part_retiree,
            "elements": [{"nom": e, "part": 100.0 * x} for e, x in m.elements[:3]],
        }

    deplaces = sorted(
        (m for m in c.retrouves if m.ecart is not None and m.f < 1000),
        key=lambda m: -m.ecart,
    )
    return {
        "resume": c.resume(),
        "n_etats": c.n_etats,
        "n_etats_rigide": c.n_etats_reduit,
        "max_re": c.max_re,
        "max_re_rigide": c.max_re_reduit,
        "deplaces": [ligne(m) for m in deplaces[:n]],
        "retires": [ligne(m) for m in sorted(c.retires, key=lambda m: -m.lam.real)[:n]],
        "perdus": [ligne(m) for m in c.perdus[:n]],
        "nouveaux": len(c.nouveaux),
    }


def detail_mode(
    t: Tables,
    indice: int,
    n_etats: int = 8,
    n_sensibilites: int = 8,
    n_parametres: int = 10,
) -> dict:
    """La forme modale et les sensibilites d'UN mode, designe par son indice.

    La table des modes dit QUI oscille ; elle ne dit ni comment les etats se
    placent les uns par rapport aux autres, ni quel coefficient du modele
    commande l'amortissement. C'est pourtant la que se joue la lecture d'un
    mode : deux vitesses rotor en phase, c'est un mode local ou tout le
    groupe bouge ensemble ; a 180 degres, c'est un mode inter-machines. Et
    la sensibilite designe l'equation ou agir -- un statisme, un gain de
    boucle -- au lieu de laisser tatonner.

    ``indice`` est celui que rend :func:`modes`, c'est-a-dire le rang dans le
    spectre rendu par ``ss.eig()``, pas le rang dans la table affichee (qui
    est triee par amortissement et filtree en frequence).

    ``parametres`` va un cran plus loin que ``sensibilites`` : non plus le
    coefficient de ``A`` mais le PARAMETRE d'equipement qui deplace le mode
    (:mod:`g2egrid.analyse.sensibilites`), normalise -- l'effet d'une variation de
    100 % du parametre, point de fonctionnement fixe.
    """
    a = construire(t)
    ss = a.etat_stationnaire()
    lam = ss.eig()
    if not 0 <= indice < len(lam):
        raise ValueError(
            f"mode {indice} hors du spectre : l'indice va de 0 a {len(lam) - 1} "
            f"({len(lam)} modes)"
        )
    v = complex(lam[indice])
    return {
        "indice": int(indice),
        "re": v.real,
        "im": v.imag,
        "f": float(modal.frequence(np.array([v]))[0]),
        "zeta": float(modal.amortissement(np.array([v]))[0]),
        "forme": [
            {"nom": nom, "module": mod, "phase": ph}
            for nom, mod, ph in modal.forme_modale(ss, indice, n=n_etats)
        ],
        "sensibilites": [
            {"ligne": i, "colonne": j, "valeur": v_}
            for i, j, v_ in modal.sensibilite_top(ss, indice, n=n_sensibilites)
        ],
        "projection": _projection(a.net, ss, v),
        **_sensibilites_parametres(a, v, n_parametres),
    }


def _projection(net: Network, ss, valeur: complex) -> dict:
    """Le mode vu sur le schema : par bus, la participation cumulee des
    elements qui y sont raccordes (le bus lui-meme compris) et le phaseur --
    module relatif, phase -- de leur etat le plus participant ; par branche,
    sa participation. Les phases sont relatives a l'etat le plus participant
    de tout le mode, comme la forme modale."""
    lam, Vr, p = modal.spectre(ss)
    k = int(np.argmin(np.abs(lam - valeur)))
    v, pk = Vr[:, k], p[k]
    v = v / v[int(np.argmax(pk))]
    rang = {nom: i for i, nom in enumerate(ss.states)}
    par_el: dict[str, list[int]] = {}
    for nom, i in rang.items():
        par_el.setdefault(nom.split(".", 1)[0], []).append(i)
    bus: dict[int, dict] = {}
    branches = []
    vmax = float(np.abs(v).max()) or 1.0
    for el_ in net.all_elements:
        idx = par_el.get(el_.name)
        if not idx:
            continue
        part = float(pk[idx].sum())
        lies = [int(b) for b in (getattr(el_, "buses", None) or [el_.bus])]
        if len(lies) >= 2 and len(set(lies)) >= 2 and not getattr(el_, "ports_dc", ()):
            branches.append(
                {"element": el_.name, "de": lies[0], "vers": lies[1], "part": part}
            )
            continue
        b = lies[0]
        i = max(idx, key=lambda j: pk[j])
        d = bus.setdefault(
            b, {"part": 0.0, "pmax": -1.0, "module": 0.0, "phase": 0.0, "etat": ""}
        )
        d["part"] += part
        if pk[i] > d["pmax"]:
            d.update(
                pmax=float(pk[i]),
                module=float(abs(v[i]) / vmax),
                phase=float(np.degrees(np.angle(v[i]))),
                etat=ss.states[i],
            )
    for d in bus.values():
        d.pop("pmax")
    return {"bus": {str(b): d for b, d in bus.items()}, "branches": branches}


def _sensibilites_parametres(a: Assemble, lam: complex, n: int) -> dict:
    """Les ``n`` parametres qui deplacent le plus le mode ``lam``.

    Le mode est designe par sa VALEUR, pas par son indice : le calcul refait le
    modele d'etat, et rien ne garantit que ``eig`` y range le spectre dans le
    meme ordre au dernier bit pres. Un reseau reduit en indice (noeuds ideaux)
    n'a pas de sensibilite analytique : on le dit plutot que de lancer des
    differences finies sur tout le reseau depuis une requete interactive.
    """
    if n <= 0:
        return {"parametres": []}
    try:
        sens = sensibilites_parametriques(a.net, lam)
    except NotImplementedError as err:
        return {"parametres": [], "parametres_erreur": str(err)}
    return {
        "parametres": [
            {
                "element": s.element,
                "parametre": s.parametre,
                "valeur": s.valeur,
                "re": s.normalisee.real,
                "im": s.normalisee.imag,
                "df": s.df,
                "dzeta": s.dzeta,
            }
            for s in sens[:n]
        ]
    }


def signaux(t: Tables) -> dict:
    a = construire(t)
    ss = a.etat_stationnaire()
    return {
        "entrees": list(ss.inputs),
        "etats": list(ss.states),
        "sorties": list(ss.outputs),
        "mesures": noms_mesures(t),
        "manoeuvres": manoeuvrables(a.net),
    }


def manoeuvrables(net: Network) -> dict:
    """Ce qu'une :class:`~g2egrid.temporel.Manoeuvre` peut viser dans ce reseau.

    Les noms sont ceux du modele ASSEMBLE : ce sont eux que la simulation
    connait. ``bus`` : les noeuds qui ont une tension d'etat, les seuls ou un
    saut de phase a un sens.
    """
    from g2egrid.elements import (
        IB,
        Branch,
        ILoad,
        LineDC,
        PLoad,
        ReactanceShunt,
        RLLoad,
        ZLoad,
    )

    elements = []
    for el in net.elements:
        if isinstance(el, Branch):
            genre = "branche"
        elif isinstance(el, LineDC):
            genre = "cable_dc"
        elif "mode_lim" in el.params:
            genre = "convertisseur"
        elif isinstance(el, (RLLoad, ZLoad, ILoad, PLoad)):
            genre = "charge"
        elif isinstance(el, ReactanceShunt):
            genre = "shunt"
        elif (
            "theta_out" in el.svc_outputs
            or "theta_pll_out" in el.svc_outputs
            or isinstance(el, IB)
        ):
            genre = "unite"
        else:
            continue
        elements.append({"nom": el.name, "genre": genre})
    elements += [
        {"nom": n, "genre": "shunt"} for n in getattr(net, "condensateurs", {})
    ]
    bus = [b for b, e in sorted(net.buses.items()) if "vd" in e.states]
    return {
        "elements": elements,
        "bus": bus,
        "bus_dc": sorted(int(b) for b in getattr(net, "bus_specs_dc", {})),
        "genres": list(GENRES_MANOEUVRE),
    }


def reglage_boucles(
    famille: str, reglages: dict | None = None, f0: float = 50.0
) -> dict:
    """Les boucles d'un GFM ou d'un GFL : intention (tau, zeta) et gains.

    Voir :mod:`g2egrid.reseau.reglage`. ``reglages`` : les reglages EFFECTIFS de
    l'unite (famille puis bus), tels que l'assembleur les passera.
    """
    from g2egrid.reseau import reglage

    propres = {k: v for k, v in (reglages or {}).items() if isinstance(v, (int, float))}
    return {"boucles": reglage.boucles(famille, propres, f0)}


def machine_standard(
    valeurs: dict, f0: float = 50.0, sens: str = "vers_standard"
) -> dict:
    """Donnees constructeur <-> parametres du modele de Park d'une machine.

    ``sens="vers_standard"`` : ``valeurs`` porte les parametres du modele
    (defauts completes) et l'on rend la fiche (Xd, X'd, T'd0...) ;
    ``sens="depuis_standard"`` : l'inverse, pour PARAMS_PAR_BUS.
    """
    from g2egrid.elements.sync_machine import SM_DEFAULT
    from g2egrid.reseau import reglage

    if sens == "vers_standard":
        p = {**SM_DEFAULT, **{k: v for k, v in valeurs.items() if k in SM_DEFAULT}}
        return {"standard": reglage.standard_depuis_machine(p, f0)}
    if sens == "depuis_standard":
        return {"parametres": reglage.machine_depuis_standard(valeurs, f0)}
    raise ValueError(f"sens {sens!r} : 'vers_standard' ou 'depuis_standard'")


def equations_modele(nom: str | None = None) -> dict:
    """Les equations d'un modele, ecrites par SymPy (:mod:`g2egrid.elements.equations`).

    Sans ``nom`` : la liste des modeles documentables.
    """
    from g2egrid.elements import equations as eq

    if nom is None:
        return {"modeles": list(eq.classes())}
    return eq.equations(nom)


def court_circuit(
    t: Tables,
    gfm: str = "tension",
    gfl: str = "aucun",
    charges: bool = False,
    dynamique: bool = False,
    f: float = 1.0,
    machine: str = "soustransitoire",
) -> dict:
    """La puissance de court-circuit et le SCR de chaque noeud alternatif
    (:func:`g2egrid.analyse.court_circuit.puissance_court_circuit`) ; ``machine`` :
    ``"soustransitoire"`` (x''d et facteurs CEI 60909, defaut) ou
    ``"transitoire"`` (x'd, ancien calcul)."""
    from g2egrid.analyse.court_circuit import puissance_court_circuit

    return puissance_court_circuit(
        t,
        gfm=gfm,
        gfl=gfl,
        charges=charges,
        dynamique=dynamique,
        f=f,
        machine=machine,
    ).vers_json()


def blocs_types() -> dict:
    """Les types de blocs de l'editeur, leurs champs et les signaux du reseau."""
    from g2egrid.elements import blocs

    return {
        "types": {
            k: {
                "entrees": v["entrees"],
                "champs": list(v["champs"]),
                "etat": v["etat"],
                "titre": v["titre"],
                "categorie": v["categorie"],
                "sorties": list(v.get("sorties", ("y",))),
                "noms_entrees": list(v.get("noms_entrees", ())),
            }
            for k, v in blocs.TYPES.items()
        },
        "categories": [{"id": k, "titre": t} for k, t in blocs.CATEGORIES],
        "entrees": list(blocs.SIGNAUX_ENTREE),
        "sorties": list(blocs.SIGNAUX_SORTIE),
        "reserves": dict(blocs.PARAMS_RESERVES),
        "exemples": blocs.exemples(),
    }


def blocs_compiler(schema: dict) -> dict:
    """Compile un schema-bloc : ses variables et ses equations en LaTeX.

    Un schema invalide leve :class:`g2egrid.elements.blocs.SchemaInvalide` (une
    ``ValueError``, donc une 422 avec la raison)."""
    import sympy as sp

    from g2egrid.elements import blocs

    c = blocs.compiler(schema)
    libres = set().union(*(e.free_symbols for e in c.f + c.g + c.h))
    return {
        "nom": c.nom,
        "parametres": c.parametres,
        "etats": c.etats,
        "algebriques": c.algebriques,
        "consignes": c.consignes,
        "equations": c.latex(),
        "sorties": [
            rf"i_d = {sp.latex(c.h[0])}",
            rf"i_q = {sp.latex(c.h[1])}",
        ],
        "sorties_texte": [f"id = {sp.sstr(c.h[0])}", f"iq = {sp.sstr(c.h[1])}"],
        "entrees_reseau": sorted(
            str(x) for x in libres if str(x) in ("vd_t", "vq_t", "w_ref", "theta_g")
        ),
        "reserves": c.reserves,
    }


def _nom_element(t: Tables, table: str, index: int) -> list[str]:
    """Les noms que l'assembleur donne a l'element de ``table[index]`` (le
    premier present dans le reseau est le bon)."""
    r = list(getattr(t, table)[index])
    if table.startswith(("T_GFM", "T_GFL")):
        return [f"{table[2:5]}{int(r[1])}"]
    if table in ("T_SMs", "T_IBs"):
        return [f"{table[2:4]}{int(r[1])}"]
    if table in ("T_Lines_SI", "T_Lines_pu"):
        return [f"L{int(r[0])}_{int(r[1])}_{int(r[2])}"]
    if table == "T_Transfos":
        return [f"T{int(r[0])}_{int(r[1])}_{int(r[2])}"]
    if table == "T_LinesDC":
        return [f"CDC{int(r[0])}_{int(r[1])}_{int(r[2])}"]
    if table == "T_Loads":
        return [f"Load{r[4]}{int(r[0])}_{int(r[1])}"]
    if table == "T_Shunts":
        return [f"Xsh{int(r[0])}_{int(r[1])}"]
    if table == "T_Blocs":
        return [f"Blocs{int(r[0])}_{int(r[1])}"]
    if table == "T_DCDC":
        return [f"DCDC{int(r[0])}_{int(r[1])}"]
    if table == "T_LoadsDC":
        return [f"LoadDC{int(r[0])}_{int(r[1])}"]
    if table == "Y_DER":
        return [f"Load{int(r[0])}", f"bus{int(r[0])}"]
    if table == "T_BusDC":
        return [f"bus{int(r[0])}"]
    if table == "T_VSCs":
        return [f"GFL{int(r[1])}", f"GFM{int(r[1])}"]
    return []


def fiche_element(t: Tables, table: str, index: int) -> dict:
    """La fiche d'un element du reseau (:mod:`g2egrid.elements.fiches`) : ses equations
    ordonnees par sous-systeme, la valeur de ses parametres et, si le load
    flow converge, son point d'equilibre."""
    from g2egrid.elements import fiches

    if table not in TABLES:
        raise ValueError(
            f"table inconnue : {table!r} (tables connues : {', '.join(TABLES)})"
        )
    n = len(getattr(t, table))
    if not 0 <= index < n:
        raise ValueError(
            f"{table} n'a pas de ligne d'indice {index} : "
            + (f"l'indice va de 0 a {n - 1}" if n else "la table est vide")
        )
    net = assembler(t)
    noms = _nom_element(t, table, index)
    par_nom = {e.name: e for e in net.all_elements}
    el = next((par_nom[n] for n in noms if n in par_nom), None)
    if el is None:
        raise ValueError(f"{table}[{index}] : aucun modele dynamique a decrire")
    equilibre = {}
    try:
        net.run_loadflow()
        net.solve_equilibrium()
        equilibre = {n: float(v) for n, v in zip(el.states, el.x0)}
    except Exception:  # noqa: BLE001 -- la fiche reste utile sans equilibre
        equilibre = {}
    f = fiches.fiche(type(el))
    p = el.pdict
    f["element"] = {
        "nom": el.name,
        "modele": type(el).model_name or type(el).__name__,
        "bus": [int(b) for b in el.buses],
        "Sn": getattr(el, "Sn", None),
    }
    f["parametres"] = [
        {
            "nom": k,
            "mathml": fiches.mathml_symbole(k),
            "valeur": float(v),
            "description": fiches.DESCRIPTIONS.get(k, ""),
        }
        for k, v in p.items()
        if isinstance(v, (int, float, np.floating))
    ]
    f["equilibre"] = equilibre
    return f


def disjoncteurs(t: Tables) -> dict:
    """L'etat des disjoncteurs : ce qu'on peut ouvrir, et ce que l'ouverture
    de ceux qui le sont a fait (iles, bus hors tension, nouveau slack)."""
    rap = disj.en_service(t)[1] if t.OUVERTS else disj.Rapport()
    return {
        "ouvrables": disj.ouvrables(t),
        "ouverts": [[a, int(b)] for a, b in t.OUVERTS],
        "hors_tension": rap.hors_tension,
        "iles": rap.iles,
        "messages": rap.messages,
    }


def noms_mesures(t: Tables) -> list[str]:
    """Les mesures physiques disponibles sur ce reseau.

    On les obtient en mesurant le point de fonctionnement lui-meme, sur une
    trajectoire d'UN point : le catalogue des mesures depend du reseau
    assemble (un bus de plus, quatre mesures de plus), et le deduire d'une
    liste ecrite a la main dans l'interface serait un doublon qui derive.
    Le cout est d'une resolution algebrique, et le modele est deja en cache.
    """
    a = construire(t)
    m = a.non_lineaire()
    return sorted(_mesures_equilibre(m))


def _mesures_equilibre(m: ModeleNonLineaire) -> dict:
    """Les mesures au point de fonctionnement, sans integrer. Les variables
    algebriques sont celles de l'equilibre : un modele a noeuds ideaux ne
    se relit pas en re-resolvant a ``x`` fige."""
    c = m.c
    tr = Trajectoire(
        t=np.array([0.0]),
        x=c.x0()[None, :],
        etats=list(c.state_names),
        algebriques=[(c.z0(), c.u0())],
    )
    return mesures_physiques(m, tr)


#: Les etats qu'on regarde en premier, par ordre de preference : la vitesse
#: d'un rotor, la frequence d'un grid-forming, la puissance filtree, la PLL
#: d'un grid-following. Prendre les premiers etats venus donnerait les
#: tensions des bus, qui viennent en tete de l'assemblage et ne disent rien
#: d'un transitoire electromecanique.
_PREFERES = (".Dwr", ".Dw", ".w", ".pm", ".theta_pll", ".theta", ".vd")


def _etats_par_defaut(etats: list[str], n: int = 4) -> list[str]:
    out: list[str] = []
    for suffixe in _PREFERES:
        for e in etats:
            if e.endswith(suffixe) and e not in out:
                out.append(e)
                if len(out) >= n:
                    return out
    return out or list(etats[:n])


def _evenements(brut: list[dict] | None, entree: str, amplitude: float, t0: float):
    """Les evenements du scenario, ou l'echelon unique a defaut.

    L'interface a commence par un echelon et garde ce cas simple ; un
    scenario a plusieurs evenements n'est que sa generalisation, et les deux
    doivent aboutir au meme objet pour ne pas faire deux chemins de code.
    """
    if not brut:
        return [Evenement(t0, entree, amplitude)]
    out = []
    for e in brut:
        genre = str(e.get("genre", "echelon"))
        if genre in GENRES_MANOEUVRE:
            out.append(
                Manoeuvre(
                    instant=float(e["instant"]),
                    genre=genre,
                    cible=str(e["cible"]),
                    dp=float(e.get("dp", 0.0)),
                    dq=float(e.get("dq", 0.0)),
                    angle=float(e.get("angle", 0.0)),
                    R=float(e["R"]) if e.get("R") is not None else None,
                    X=float(e["X"]) if e.get("X") is not None else None,
                )
            )
            continue
        out.append(
            Evenement(
                instant=float(e["instant"]),
                entree=str(e["entree"]),
                amplitude=float(e["amplitude"]),
                genre=genre,
                duree=float(e.get("duree", 0.0)),
            )
        )
    return out


def _sans_nan(v) -> list:
    """Une serie en JSON : un ``NaN`` (grandeur d'un element retire par une
    manoeuvre) devient ``null``, que JSON sait dire et qu'un trace sait sauter."""
    return [
        None if not np.isfinite(x) else float(x) for x in np.asarray(v, dtype=float)
    ]


def temporel(
    t: Tables,
    entree: str,
    amplitude: float = 0.01,
    duree: float = 1.0,
    t0: float = 0.0,
    n: int = 401,
    etats: list[str] | None = None,
    voie: str = "les_deux",
    evenements: list[dict] | None = None,
    mesures: list[str] | None = None,
    methode: str = "auto",
    pas: float | None = None,
    limiteurs: bool = True,
    protections_dc: bool = False,
) -> dict:
    """La reponse a un scenario, par l'une des deux voies ou les deux.

    ``voie`` vaut ``lineaire``, ``non_lineaire`` ou ``les_deux``. Dans ce
    dernier cas on rend aussi l'ecart maximal entre les deux, rapporte a
    l'amplitude de la reponse : c'est la mesure de la non-linearite du point
    de fonctionnement, et le controle croise des deux chemins de code.

    ``evenements`` decrit un scenario a plusieurs perturbations ; a defaut
    c'est l'echelon unique ``(entree, amplitude, t0)``. Une MANOEUVRE
    (ouverture, charge, saut de phase, defaut, elimination) y figure avec son
    ``genre`` et sa ``cible`` (cf. :class:`g2egrid.temporel.Manoeuvre`).

    ``methode`` choisit l'integrateur de la voie non lineaire : ``Radau``
    (pas adaptatif, pour la validation fine), ``BDF2`` (pas fixe ``pas``,
    DAE monolithique : defauts et manoeuvres), ou ``auto`` -- BDF2 des que le
    scenario contient une manoeuvre, Radau sinon. ``limiteurs`` active les
    limitations de courant et les butees d'excitation ; leurs basculements
    sont rendus dans ``basculements``. ``protections_dc`` arme le plan de
    protection type du reseau continu (:func:`g2egrid.temporel.protections_dc`)
    ; ses declenchements sont rendus dans ``declenchements``.

    ``mesures`` demande des grandeurs PHYSIQUES (module de tension,
    frequence, transits) en plus des etats bruts. Elles ne sont rendues que
    sur la voie NON LINEAIRE : elles se reconstruisent a partir des
    ``(x, z, u, y)`` re-resolus le long de la trajectoire, ce que la reponse
    lineaire -- une exponentielle de matrice sur les etats -- ne fournit pas.

    LA VOIE LINEAIRE NE SUIT QU'UN ECHELON. Sa reponse est calculee en forme
    close ; un scenario a plusieurs evenements, ou une rampe, sortirait de ce
    cadre. Quand le scenario n'est pas un echelon unique, la voie lineaire
    est donc silencieusement omise plutot que rendue fausse, et le champ
    ``voies`` du resultat dit ce qui a reellement ete calcule.
    """
    if not (np.isfinite(duree) and duree > 0):
        raise ValueError(f"duree = {duree} : la duree simulee doit etre > 0 (secondes)")
    if int(n) < 2:
        raise ValueError(f"n = {n} : il faut au moins 2 instants de calcul")
    if pas is not None and not (np.isfinite(pas) and pas > 0):
        raise ValueError(f"pas = {pas} : le pas de temps doit etre > 0 (secondes)")
    a = construire(t)
    ss = a.etat_stationnaire()
    evs = _evenements(evenements, entree, amplitude, t0)
    for e in evs:
        if isinstance(e, Evenement) and e.entree not in ss.inputs:
            raise ValueError(f"entree inconnue : {e.entree}")
    consignes = [e for e in evs if isinstance(e, Evenement)]
    echelon_simple = (
        len(evs) == 1
        and consignes
        and evs[0].genre == "echelon"
        and evs[0].instant >= 0.0
    )
    if voie in ("lineaire", "les_deux") and not echelon_simple:
        voie = "non_lineaire"
    if consignes:
        entree, amplitude = consignes[0].entree, consignes[0].amplitude
    t0 = evs[0].instant
    etats = etats or _etats_par_defaut(ss.states)
    manquants = [e for e in etats if e not in ss.states]
    if manquants:
        raise ValueError(f"etat(s) inconnu(s) : {', '.join(manquants)}")
    idx = [ss.states.index(e) for e in etats]
    temps = np.linspace(0.0, duree, int(n))

    series: dict[str, dict[str, list[float]]] = {e: {} for e in etats}
    if voie in ("lineaire", "les_deux"):
        u = np.zeros(len(ss.inputs))
        u[ss.iu(entree)] = amplitude
        xl = reponse_lineaire(ss.A, ss.B, np.maximum(temps - t0, 0.0), u)
        for e, i in zip(etats, idx):
            series[e]["lineaire"] = xl[:, i].tolist()
    releves: dict[str, list[float]] = {}
    basculements: list[dict] = []
    declenchements: list[dict] = []
    if voie in ("non_lineaire", "les_deux"):
        m = a.non_lineaire()
        uext_fn = m.scenario(*evs)
        if methode == "auto":
            # BDF2 pour une manoeuvre ou une protection, et toujours pour un
            # reseau a noeuds ideaux (indice 2, que Radau n'integre pas)
            methode = (
                "BDF2"
                if protections_dc
                or m.ideaux
                or any(isinstance(e, Manoeuvre) for e in evs)
                else "Radau"
            )
        plan = None
        if protections_dc:
            from g2egrid.temporel import protections_dc as plan_dc

            if not a.net.bus_specs_dc:
                raise ValueError("protections DC : le reseau n'a pas de reseau continu")
            plan = plan_dc(a.net)
        if methode not in ("Radau", "BDF2"):
            raise ValueError(
                f"integrateur {methode!r} inconnu ('auto', 'Radau', 'BDF2')"
            )
        tr = m.simuler(
            (0.0, duree),
            uext_fn=uext_fn,
            t_eval=temps,
            methode=methode,
            pas=pas,
            limiteurs=limiteurs,
            protections=plan,
        )
        declenchements = [
            {
                "t": float(t_),
                "protection": nom_,
                "manoeuvres": [{"genre": g_, "cible": c_} for g_, c_ in mans_],
            }
            for t_, nom_, mans_ in tr.declenchements
        ]
        basculements = [
            {"t": float(t_), "element": el_, "drapeau": d_, "valeur": float(v_)}
            for t_, el_, d_, v_ in tr.basculements
        ]
        # les premieres colonnes sont celles du modele de depart, dans son
        # ordre ; une manoeuvre peut en AJOUTER (repere) et en vider (NaN)
        x0 = m.c.x0()
        xn = tr.x[:, : len(x0)] - x0
        for e, i in zip(etats, idx):
            series[e]["non_lineaire"] = _sans_nan(xn[:, i])
        if mesures:
            dispo = _mesures_equilibre(m)
            inconnues = [q for q in mesures if q not in dispo]
            if inconnues:
                raise ValueError(f"mesure(s) inconnue(s) : {', '.join(inconnues)}")
            vals = mesures_physiques(m, tr, uext_fn=uext_fn, quoi=mesures)
            releves = {k: _sans_nan(v) for k, v in vals.items()}

    ecart = None
    if voie == "les_deux":
        num = max(
            float(
                np.abs(
                    np.array(series[e]["non_lineaire"])
                    - np.array(series[e]["lineaire"])
                ).max()
            )
            for e in etats
        )
        den = max(float(np.abs(np.array(series[e]["lineaire"])).max()) for e in etats)
        ecart = {"absolu": num, "relatif": num / den if den > 1e-15 else float("nan")}
    return {
        "t": temps.tolist(),
        "series": series,
        "ecart": ecart,
        "entree": entree,
        "voies": (["lineaire", "non_lineaire"] if voie == "les_deux" else [voie]),
        "mesures": releves,
        "integrateur": methode if voie != "lineaire" else None,
        "basculements": basculements,
        "declenchements": declenchements,
        "evenements": [
            {
                "instant": e.instant,
                "entree": e.entree,
                "amplitude": e.amplitude,
                "genre": e.genre,
                "duree": e.duree,
            }
            if isinstance(e, Evenement)
            else {
                "instant": e.instant,
                "genre": e.genre,
                "cible": e.cible,
                "dp": e.dp,
                "dq": e.dq,
                "angle": e.angle,
                "R": e.R,
                "X": e.X,
            }
            for e in evs
        ],
    }
