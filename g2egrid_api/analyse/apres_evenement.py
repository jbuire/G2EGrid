"""Apres un evenement : simulation d'un defaut, d'un declenchement ou d'une manoeuvre, et le reseau qui en resulte.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import re
from collections.abc import Sequence
import numpy as np
from g2egrid.reseau.network import Network
from g2egrid.reseau.tables import Tables
from g2egrid.temporel import Evenement, Manoeuvre
from g2egrid.reseau import disjoncteurs as disj
from g2egrid.temporel import mesures as mesures_physiques
from .assemblage import construire
from .resultats import _evenements, _mesures_equilibre, _sans_nan, manoeuvrables


# --------------------------------------------------------- apres un evenement
def _nombre(x) -> float | None:
    """Un flottant en JSON : ``NaN`` (element retire, bus hors tension) -> null."""
    x = float(x)
    return x if np.isfinite(x) else None


def _relais_ac(net: Network, relais: dict) -> list:
    """Le plan :func:`g2egrid.temporel.protections_ac` decrit par l'interface :
    ``{"v_max", "retard", "retombee" (None : verrouille), "action",
    "familles"}``."""
    from g2egrid.temporel import protections_ac

    retombee = relais.get("retombee", 0.98)
    return protections_ac(
        net,
        v_max=float(relais.get("v_max", 1.10)),
        retard=float(relais.get("retard", 0.2)),
        retombee=None if retombee in (None, "", 0) else float(retombee),
        action=str(relais.get("action", "ouverture")),
        familles=tuple(relais.get("familles") or ("GFL",)),
    )


def _valeur_finale(v: np.ndarray, temps: np.ndarray, fenetre: float) -> float:
    """La moyenne sur la derniere ``fenetre`` secondes : le regime atteint,
    sans l'oscillation residuelle. ``NaN`` si l'element n'est plus la."""
    if not np.isfinite(v[-1]):
        return float("nan")
    w = v[(temps >= temps[-1] - fenetre) & np.isfinite(v)]
    return float(w.mean())


def _regime_lf(t: Tables, evs: list, ouvertures: list[str]) -> dict:
    """Le regime permanent du LOAD FLOW apres les memes ouvertures : celui
    vers lequel la simulation devrait converger, aux regulations pres
    (statisme des machines, charges a impedance constante, frequence).

    Il n'existe que si le scenario se ramene a des ouvertures : une consigne,
    un echelon de charge ou un defaut tenu n'ont pas d'equivalent dans les
    tables du load flow."""
    from g2egrid.scenarios.cascade import tensions_bus

    if any(isinstance(e, Evenement) for e in evs):
        return {"ok": False, "raison": "le scenario change une consigne"}
    tous = {o["nom"]: (o["table"], o["id"]) for o in disj.ouvrables(t)}
    defauts: set[str] = set()
    for e in evs:
        if e.genre == "defaut":
            defauts.add(e.cible)
        elif e.genre == "elimination":
            defauts.discard(e.cible)
        elif e.genre not in ("ouverture", "blocage"):
            return {"ok": False, "raison": f"{e.genre} : pas d'equivalent en load flow"}
    if defauts:
        return {"ok": False, "raison": "un defaut n'est pas elimine"}
    noms = list(dict.fromkeys(ouvertures))
    inconnus = [c for c in noms if c not in tous]
    if inconnus:
        return {"ok": False, "raison": f"{', '.join(inconnus)} : pas de disjoncteur"}
    t2 = t.copie()
    t2.OUVERTS = list(
        dict.fromkeys([(a, int(b)) for a, b in t.OUVERTS] + [tous[c] for c in noms])
    )
    try:
        V, ok = tensions_bus(t2)
    except (RuntimeError, ValueError, np.linalg.LinAlgError) as err:
        return {"ok": False, "raison": f"load flow : {err}"}
    if not ok:
        return {"ok": False, "raison": "le load flow ne converge pas"}
    return {"ok": True, "V": {str(b): v for b, v in V.items()}, "ouverts": noms}


def apres_evenement(
    t: Tables,
    evenements: list[dict],
    duree: float = 5.0,
    n: int = 501,
    pas: float | None = 1e-3,
    limiteurs: bool = True,
    relais: dict | None = None,
    protections_dc: bool = False,
    v_bas: float = 0.9,
    v_haut: float = 1.1,
    fenetre: float = 0.2,
) -> dict:
    """Les grandeurs electriques du reseau apres un evenement transitoire
    (ouverture, defaut, blocage, saut de phase) ou de consigne (echelon,
    rampe) : simulation non lineaire (BDF2) de tout le reseau, puis, pour
    CHAQUE bus et CHAQUE unite, ce qu'il faut pour juger l'evenement.

    - ``bus`` : tension avant (``V0``), extremes du transitoire (``Vmin``,
      ``Vmax`` et leurs instants), regime atteint (``Vfin``, moyenne sur la
      derniere ``fenetre`` s), ecart ``dV``, temps passe hors de
      ``[v_bas, v_haut]`` et, quand le scenario se ramene a des ouvertures,
      la tension du load flow apres les memes ouvertures (``V_lf``) ;
    - ``unites`` : P et Q avant et apres de chaque unite (machine,
      convertisseur, charge), et l'instant ou elle est sortie ;
    - ``frequence`` : enveloppe min / max / moyenne des frequences des bus
      des sources de tension (machines, GFM) ;
    - ``relais`` arme des relais de surtension sur les convertisseurs
      (:func:`g2egrid.temporel.protections_ac`) : leurs demarrages,
      retombees et declenchements sont rendus -- c'est ainsi qu'on suit une
      cascade de surtension.

    Le PAS de BDF2 vaut 1 ms par defaut, pas le ``1/(400.f0)`` du temporel :
    on regarde ici des secondes (regulations, temporisations des relais) ;
    BDF2 amortit les resonances plus rapides que le pas au lieu de les
    suivre. Un defaut proche ou une resonance de filtre demandent un pas
    plus fin (``pas=None`` revient a ``1/(400.f0)``).
    """
    from g2egrid.temporel import protections_dc as plan_dc

    if not evenements:
        raise ValueError("aucun evenement : le scenario est vide")
    if duree <= 0:
        raise ValueError("la duree doit etre positive")
    if not v_bas < v_haut:
        raise ValueError("bornes de tension : v_bas < v_haut")
    a = construire(t)
    ss = a.etat_stationnaire()
    evs = _evenements(evenements, "", 0.0, 0.0)
    for e in evs:
        if isinstance(e, Evenement) and e.entree not in ss.inputs:
            raise ValueError(f"entree inconnue : {e.entree}")
    m = a.non_lineaire()
    plan = _relais_ac(a.net, relais) if relais else []
    if protections_dc:
        if not a.net.bus_specs_dc:
            raise ValueError("protections DC : le reseau n'a pas de reseau continu")
        plan += plan_dc(a.net)
    retards = {p.nom: p.retard for p in plan}
    sc = m.scenario(*evs)
    temps = np.linspace(0.0, float(duree), int(n))
    tr = m.simuler(
        (0.0, float(duree)),
        uext_fn=sc,
        t_eval=temps,
        methode="BDF2",
        pas=pas,
        limiteurs=limiteurs,
        protections=plan or None,
    )

    # les sorties d'unites : evenements du scenario puis declenchements
    sortie: dict[str, float] = {}
    for e in evs:
        if isinstance(e, Manoeuvre) and e.genre in ("ouverture", "blocage"):
            sortie.setdefault(e.cible, float(e.instant))
    declenchements = []
    for t_dem, nom, mans in tr.declenchements:
        t_act = float(t_dem) + retards.get(nom, 0.0)
        declenchements.append(
            {
                "t_demarrage": float(t_dem),
                "t": t_act,
                "protection": nom,
                "manoeuvres": [{"genre": g, "cible": c} for g, c in mans],
            }
        )
        for g, c in mans:
            if g in ("ouverture", "blocage"):
                sortie.setdefault(c, t_act)

    dispo = _mesures_equilibre(m)
    bus = sorted(int(k[5:]) for k in dispo if re.fullmatch(r"V_bus\d+", k))
    unites = [
        e
        for e in manoeuvrables(a.net)["elements"]
        if e["genre"] in ("unite", "convertisseur", "charge")
        and f"P_{e['nom']}" in dispo
        and f"Q_{e['nom']}" in dispo
    ]
    # la frequence du reseau : celle des bus des SOURCES DE TENSION (machines,
    # GFM, bus infini) quand ils en portent une, sinon de tous les bus SAUF
    # ceux d'une unite sortie ou d'un defaut : a l'instant de l'ouverture,
    # l'angle d'un tel bus saute, et sa derivee vaut des centaines de hertz,
    # ce qui n'est pas une frequence.
    els = {e.name: e for e in a.net.elements}

    def bus_de(nom):
        el = els.get(nom)
        return int(el.buses[0]) if el is not None and len(el.buses) else None

    sources = {bus_de(u["nom"]) for u in unites if u["genre"] == "unite"}
    exclus = {bus_de(n) for n in sortie} | {
        int(re.sub(r"\D", "", str(e.cible)) or -1)
        for e in evs
        if isinstance(e, Manoeuvre) and e.genre == "defaut"
    }
    tous_f = sorted(
        (int(k[6:]) for k in dispo if re.fullmatch(r"fm_bus\d+", k)),
    )
    choix = [b for b in tous_f if b in sources] or [
        b for b in tous_f if b not in exclus
    ]
    freqs = [f"fm_bus{b}" for b in choix]
    quoi = (
        [f"V_bus{b}" for b in bus]
        + [f"{g}_{u['nom']}" for u in unites for g in ("P", "Q")]
        + freqs
    )
    vals = {
        k: np.asarray(v, dtype=float)
        for k, v in mesures_physiques(m, tr, uext_fn=sc, quoi=quoi).items()
    }

    lf = _regime_lf(t, evs, list(sortie))
    un = {int(r[0]): float(r[2]) for r in t.Y_DER}
    dt = np.diff(temps)
    lignes_bus = []
    for b in bus:
        v = vals[f"V_bus{b}"]
        fini = np.isfinite(v)
        if not fini.any():
            lignes_bus.append({"bus": b, "Un_kV": un.get(b, 0.0) / 1e3, "V0": None})
            continue
        i_min, i_max = int(np.nanargmin(v)), int(np.nanargmax(v))
        hors = (fini & ((v > v_haut) | (v < v_bas)))[:-1]
        v_fin = _valeur_finale(v, temps, fenetre)
        v_lf = lf["V"].get(str(b)) if lf["ok"] else None
        lignes_bus.append(
            {
                "bus": b,
                "Un_kV": un.get(b, 0.0) / 1e3,
                "V0": _nombre(v[0]),
                "Vmin": _nombre(v[i_min]),
                "t_min": float(temps[i_min]),
                "Vmax": _nombre(v[i_max]),
                "t_max": float(temps[i_max]),
                "Vfin": _nombre(v_fin),
                "dV": _nombre(v_fin - v[0]),
                "t_hors": float(np.sum(dt[hors])),
                "V_lf": v_lf,
                "ecart_lf": _nombre(v_fin - v_lf) if v_lf is not None else None,
            }
        )
    lignes_unites = []
    for u in unites:
        p, q = vals[f"P_{u['nom']}"], vals[f"Q_{u['nom']}"]
        p_fin, q_fin = (
            _valeur_finale(p, temps, fenetre),
            _valeur_finale(q, temps, fenetre),
        )
        if u["nom"] in sortie or not np.isfinite(p_fin):
            p_fin, q_fin = 0.0, 0.0
        lignes_unites.append(
            {
                "nom": u["nom"],
                "genre": u["genre"],
                "P0": _nombre(p[0]),
                "Q0": _nombre(q[0]),
                "Pfin": _nombre(p_fin),
                "Qfin": _nombre(q_fin),
                "dP": _nombre(p_fin - p[0]),
                "dQ": _nombre(q_fin - q[0]),
                "sortie": sortie.get(u["nom"]),
            }
        )
    if freqs:
        F = np.vstack([vals[k] for k in freqs])
        with np.errstate(all="ignore"):
            enveloppe = {
                "min": np.nanmin(F, axis=0),
                "max": np.nanmax(F, axis=0),
                "moy": np.nanmean(F, axis=0),
            }
        frequence = {
            "min": _nombre(np.nanmin(enveloppe["min"])),
            "max": _nombre(np.nanmax(enveloppe["max"])),
            "fin": _nombre(_valeur_finale(enveloppe["moy"], temps, fenetre)),
            "bus": [int(k[6:]) for k in freqs],
            "series": {k: _sans_nan(np.round(v, 5)) for k, v in enveloppe.items()},
        }
    else:
        frequence = {"min": None, "max": None, "fin": None, "bus": [], "series": {}}

    def serie_unite(nom: str, g: str) -> list:
        # une unite sortie ne produit plus rien : 0, plutot qu'un trou
        v = vals[f"{g}_{nom}"].copy()
        if nom in sortie:
            v[(temps >= sortie[nom]) & ~np.isfinite(v)] = 0.0
        return _sans_nan(np.round(v, 6))

    valides = [r for r in lignes_bus if r["V0"] is not None]
    haut = max(valides, key=lambda r: r["Vmax"], default=None)
    bas = min(valides, key=lambda r: r["Vmin"], default=None)
    return {
        "t": temps.tolist(),
        "bus": lignes_bus,
        "unites": lignes_unites,
        "frequence": frequence,
        "series": {
            "V": {str(b): _sans_nan(np.round(vals[f"V_bus{b}"], 6)) for b in bus},
            "P": {u["nom"]: serie_unite(u["nom"], "P") for u in unites},
            "Q": {u["nom"]: serie_unite(u["nom"], "Q") for u in unites},
        },
        "declenchements": declenchements,
        "demarrages": [{"t": float(x), "protection": p} for x, p in tr.demarrages],
        "retombees": [{"t": float(x), "protection": p} for x, p in tr.retombees],
        "basculements": len(tr.basculements),
        "regime_lf": {k: v for k, v in lf.items() if k != "V"},
        "resume": {
            "v_max": haut
            and {"bus": haut["bus"], "V": haut["Vmax"], "t": haut["t_max"]},
            "v_min": bas and {"bus": bas["bus"], "V": bas["Vmin"], "t": bas["t_min"]},
            "n_hors_bornes": sum(1 for r in valides if r["t_hors"] > 0),
            "n_sorties": len(sortie),
            "n_declenchements": len(declenchements),
        },
        "bornes": [v_bas, v_haut],
        "relais": relais or None,
        "evenements": evenements,
    }


def cascade(
    t: Tables,
    initiales: list[str],
    v_max: float = 1.10,
    familles: Sequence[str] = ("GFL",),
    ordre: str = "vague",
) -> dict:
    """La cascade de surtension en regimes permanents successifs
    (:func:`g2egrid.scenarios.cascade.cascade_surtension`) : les ``initiales`` s'ouvrent,
    puis les unites protegees dont le bus depasse ``v_max`` declenchent,
    toutes ensemble (``ordre="vague"``) ou la plus haute d'abord
    (``"une_a_une"``), jusqu'a ce que plus rien ne depasse."""
    from g2egrid.scenarios.cascade import cascade_surtension, unites_protegees

    if not initiales:
        raise ValueError("aucune unite initiale")
    familles = tuple(familles)
    c = cascade_surtension(
        t, list(initiales), float(v_max), familles=familles, ordre=ordre
    )
    return {
        "seuil": c.seuil,
        "ordre": ordre,
        "issue": c.issue,
        "initiales": c.initiales,
        "declenchees": c.declenchees,
        "bus_unite": {n: b for _tab, _i, b, n in unites_protegees(t, familles)},
        "V0": {str(b): v for b, v in c.V0.items()},
        "etapes": [
            {
                "rang": e.rang,
                "ouvertes": e.ouvertes,
                "hors_service": e.hors_service,
                "au_dessus": e.au_dessus,
                "converge": e.converge,
                "v_max": _nombre(e.v_max),
                "bus_v_max": e.bus_v_max,
                "V": {str(b): v for b, v in e.V.items()},
            }
            for e in c.etapes
        ],
    }
