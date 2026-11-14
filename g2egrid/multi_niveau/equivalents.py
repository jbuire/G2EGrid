"""Equivalents de centrales et de zones : agreger des unites en une seule.

Une UNITE est une machine synchrone (SM), un grid-forming (GFM, batterie) ou
un grid-following (GFL, PV) sur son bus basse tension, raccordee par son
transformateur a un bus haute tension. Un GROUPE d'unites est remplace, au
bus haute tension ``cible``, par :

``modele="technologie"``
    une unite equivalente PAR TECHNOLOGIE (une SM, un GFM, un GFL), chacune
    derriere le transformateur equivalent (les transformateurs du groupe en
    parallele) ;

``modele="machine"``
    une SEULE machine synchrone equivalente, qui porte aussi les GFM (leur
    inertie virtuelle et leur statisme) ; les GFL deviennent une injection a
    puissance constante (``gfl="charge"``) ou sont absorbes dans la machine
    (``gfl="absorbe"`` : sa puissance, pas d'inertie).

Les regles d'agregation sont celles de l'equivalent d'une flotte de
machines en parallele (Min, Paganini, Mallada) : les inerties et les
amortissements s'ajoutent (en MW.s et en MW/pu), les statismes aussi (en
inverse : les gains s'ajoutent). La constante de temps de la turbine
equivalente n'a pas d'expression exacte : elle est ici la moyenne des
constantes, ponderee par les gains de regulation. Un GFM a statisme
``mp`` et filtre de puissance ``wf_p`` est une machine virtuelle
d'inertie ``H = 1 / (2 mp wf_p)`` et d'amortissement ``D = 1 / mp`` : dans
une machine unique, il apporte ``H`` et, au choix, ``D`` en amortissement
(``gfm="kd"``) ou son statisme au regulateur de vitesse
(``gfm="gouverneur"``).

Le point de fonctionnement du reseau HAUTE tension est garde : le bus BT de
l'equivalent recoit la tension qui fait passer, par le transformateur
equivalent, exactement le courant que le groupe injectait au bus cible dans
le modele complet.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from g2egrid.reseau.tables import Tables, assembler

GENRES = {"T_SMs": "SM", "T_GFM1s": "GFM", "T_GFL1s": "GFL"}
TABLE = {v: k for k, v in GENRES.items()}


@dataclass
class Unite:
    genre: str  # "SM", "GFM", "GFL"
    bus: int  # bus BT
    bus_ht: int
    Sn: float  # [VA]
    transfo: tuple  # la ligne de T_Transfos
    params: dict  # parametres numeriques de l'element (pdict)
    reglages: dict  # famille | bus : ce que les tables disent
    S: complex  # puissance injectee au bus BT [pu base systeme]
    I: complex  # noqa: E741 -- courant injecte au bus HT [pu base systeme]
    V: complex  # tension BT [pu]
    vdc: float | None = None  # colonne Vdc_V des convertisseurs


@dataclass
class Groupe:
    bus: list[int]  # les bus BT des unites du groupe
    cible: int  # le bus HT ou raccorder l'equivalent
    modele: str = "technologie"  # ou "machine"
    gfl: str = "charge"  # modele "machine" : "charge" ou "absorbe"
    gfm: str = "kd"  # modele "machine" : "kd" ou "gouverneur"
    nom: str = ""
    #: ajustements de l'equivalent, par technologie : un FACTEUR sur chaque
    #: parametre de l'unite equivalente (``"H"``, ``"Ka"``, ``"mp"``...) et
    #: sur l'impedance de son transformateur (``"X"``) -- ce qu'une
    #: optimisation fait varier autour de l'agregation directe
    ajustements: dict[str, dict[str, float]] = field(default_factory=dict)


def unites(t: Tables) -> dict[int, Unite]:
    """Les unites de ``t`` raccordees par un transformateur, par bus BT."""
    net = assembler(t)
    lf = net.run_loadflow()
    if not lf.converged:
        raise ValueError("load flow du modele complet non converge")
    idx = net._bus_index()
    par_nom = {el.name: el for el in net.all_elements}
    out = {}
    for table, genre in GENRES.items():
        famille = {"SM": t.SM_PARAMS, "GFM": t.GFM_PARAMS, "GFL": t.GFL_PARAMS}[genre]
        for row in getattr(t, table):
            b = int(row[1])
            tr = next((x for x in t.T_Transfos if b in (int(x[1]), int(x[2]))), None)
            if tr is None:
                continue
            ht = int(tr[2]) if int(tr[1]) == b else int(tr[1])
            el = par_nom[f"{genre}{b}"]
            V = lf.V[idx[b]]
            S = lf.S[idx[b]]
            out[b] = Unite(
                genre=genre,
                bus=b,
                bus_ht=ht,
                Sn=float(row[2]),
                transfo=tuple(tr),
                params=dict(el.pdict),
                reglages={**famille, **t.PARAMS_PAR_BUS.get(b, {})},
                S=complex(S),
                I=complex(np.conj(S / V)),
                V=complex(V),
                vdc=float(row[3]) if len(row) > 3 else None,
            )
    return out


def _moyenne(us: list[Unite], cle: str, poids=None) -> float:
    w = np.array([u.Sn if poids is None else poids(u) for u in us], float)
    v = np.array([float(u.params[cle]) for u in us])
    return float((w * v).sum() / w.sum())


def _params_sm(sms: list[Unite], gfms: list[Unite], gfm: str) -> dict:
    """Les parametres de la machine equivalente (base : la somme des Sn,
    GFM compris)."""
    S = sum(u.Sn for u in sms) + sum(u.Sn for u in gfms)
    p = {}
    # les parametres electriques : moyenne ponderee (exacte a pu egaux)
    for cle in sms[0].params:
        if cle in ("H", "KD", "mp", "Tg", "theta_g0", "mode_e2"):
            continue
        vals = {float(u.params[cle]) for u in sms}
        if len(vals) > 1:
            p[cle] = _moyenne(sms, cle)
    # inertie [MW.s] et amortissement [MW/pu] : ils s'ajoutent
    h = sum(u.params["H"] * u.Sn for u in sms)
    kd = sum(u.params["KD"] * u.Sn for u in sms)
    gain = sum(u.Sn / (u.params["mp"] / 100) for u in sms)  # [MW/pu]
    tg = sum(u.params["Tg"] * u.Sn / u.params["mp"] for u in sms) / sum(
        u.Sn / u.params["mp"] for u in sms
    )
    for u in gfms:
        mp, wf = u.params["mp"], u.params["wf_p"]
        h += u.Sn / (2 * mp * wf)
        if gfm == "kd":
            kd += u.Sn / mp
        else:
            gain += u.Sn / mp
    p.update(H=h / S, KD=kd / S, mp=100 * S / gain, Tg=tg)
    return p


def _params_gfm(us: list[Unite]) -> dict:
    S = sum(u.Sn for u in us)
    p = {}
    for cle in us[0].params:
        if cle in ("mp", "wf_p", "nq", "wf_q", "mode_lim"):
            continue
        if len({float(u.params[cle]) for u in us}) > 1:
            p[cle] = _moyenne(us, cle)
    # statismes : les gains s'ajoutent ; inertie virtuelle : elle s'ajoute
    mp = S / sum(u.Sn / u.params["mp"] for u in us)
    nq = S / sum(u.Sn / u.params["nq"] for u in us)
    h = sum(u.Sn / (2 * u.params["mp"] * u.params["wf_p"]) for u in us) / S
    hq = sum(u.Sn / (2 * u.params["nq"] * u.params["wf_q"]) for u in us) / S
    p.update(mp=mp, nq=nq, wf_p=1 / (2 * h * mp), wf_q=1 / (2 * hq * nq))
    return p


def _params_gfl(us: list[Unite]) -> dict:
    p = {}
    for cle in us[0].params:
        if cle == "mode_lim":
            continue
        if len({float(u.params[cle]) for u in us}) > 1:
            p[cle] = _moyenne(us, cle)
    return p


def _intentions(us: list[Unite]) -> dict:
    """Les reglages par intention (``tau_CL``...) du bus, s'ils sont
    communs a tout le groupe : ils ne se moyennent pas, ils se gardent."""
    cles = set.intersection(*(set(u.reglages) for u in us)) if us else set()
    return {
        k: us[0].reglages[k]
        for k in cles
        if k.startswith(("tau_", "zeta_"))
        and len({repr(u.reglages[k]) for u in us}) == 1
    }


def _transfo_equivalent(us: list[Unite], Sb: float, ident: int, bus: int, cible: int):
    """Les transformateurs en parallele, sur la base de la somme des Sn."""
    Sn = sum(float(u.transfo[3]) for u in us)
    y = sum(float(u.transfo[3]) / Sb / complex(u.transfo[6], u.transfo[7]) for u in us)
    z = Sn / Sb / y
    tr = us[0].transfo
    u_bt, u_ht = (tr[4], tr[5]) if int(tr[1]) == us[0].bus else (tr[5], tr[4])
    return (
        ident,
        bus,
        cible,
        Sn,
        u_bt,
        u_ht,
        float(z.real),
        float(z.imag),
    ), z * Sb / Sn


def agreger(t: Tables, groupes: list[Groupe], cache: dict | None = None) -> Tables:
    """Une copie de ``t`` ou chaque groupe est remplace par son equivalent.
    ``cache`` (un dict, vide au premier appel) garde le load flow du modele
    complet d'un appel a l'autre : une optimisation agrege cent fois le
    meme reseau."""
    cache = {} if cache is None else cache
    if "unites" not in cache:
        cache["unites"], cache["tensions"] = unites(t), _tensions(t)
    tout, V_ht = cache["unites"], cache["tensions"]
    eq = t.copie()
    bus_libre = max(int(r[0]) for r in t.Y_DER) + 1
    ids = {
        k: max((int(r[0]) for r in getattr(t, k)), default=0) for k in TABLE.values()
    }
    id_tr = max((int(r[0]) for r in t.T_Transfos), default=0)
    id_ch = max((int(r[0]) for r in t.T_Loads), default=0)
    ydict = {int(r[0]): tuple(r) for r in t.Y_DER}
    ppb = dict(eq.PARAMS_PAR_BUS)
    retires: set[int] = set()
    nouveaux_y, nouveaux_tr, charges = [], [], []
    ajouts = {k: [] for k in TABLE.values()}
    for g in groupes:
        us = [tout[b] for b in g.bus]
        retires |= set(g.bus)
        slack = any(int(ydict[b][1]) == 1 for b in g.bus)
        par = {k: [u for u in us if u.genre == k] for k in ("SM", "GFM", "GFL")}
        if g.modele == "technologie":
            lots = [(k, par[k]) for k in ("SM", "GFM", "GFL") if par[k]]
        elif g.modele == "machine":
            if not par["SM"]:
                raise ValueError(f"groupe {g.nom} : une machine unique demande une SM")
            dans = par["SM"] + par["GFM"] + (par["GFL"] if g.gfl == "absorbe" else [])
            lots = [("SM", dans)]
            if g.gfl == "charge" and par["GFL"]:
                s = sum(u.S for u in par["GFL"]) * t.Sb
                id_ch += 1
                charges.append((id_ch, g.cible, -s.real, -s.imag, "P"))
        else:
            raise ValueError(f"modele {g.modele!r} : 'technologie' ou 'machine'")
        for k_lot, (genre, membres) in enumerate(lots):
            b = bus_libre
            bus_libre += 1
            id_tr += 1
            tr, z_sys = _transfo_equivalent(membres, t.Sb, id_tr, b, g.cible)
            adj = dict(g.ajustements.get(genre, {}))
            kx = adj.pop("X", 1.0)
            if kx != 1.0:
                tr = (*tr[:6], tr[6] * kx, tr[7] * kx)
                z_sys = z_sys * kx
            nouveaux_tr.append(tr)
            # la tension BT qui fait passer le meme courant au bus HT
            i_ht = sum(u.I for u in membres)
            V = V_ht[g.cible] + z_sys * i_ht
            S = V * np.conj(i_ht)
            ref = ydict[membres[0].bus]
            typ = (
                1
                if (slack and k_lot == 0)
                else (3 if genre == "GFL" and int(ref[1]) == 3 else 2)
            )
            nouveaux_y.append(
                (
                    b,
                    typ,
                    ref[2],
                    float(abs(V)),
                    float(np.angle(V)),
                    float(S.real * t.Sb),
                    float(S.imag * t.Sb) if typ == 3 else 0.0,
                    0.0,
                    0.0,
                )
            )
            Sn = sum(u.Sn for u in membres)
            if genre == "SM":
                sms = [u for u in membres if u.genre == "SM"]
                gfms = [u for u in membres if u.genre == "GFM"]
                reg = _params_sm(sms, gfms, g.gfm)
                reg = {**reg, **_intentions(sms)}
                # la puissance des GFL absorbes : sans inertie ni statisme
                gfls = [u for u in membres if u.genre == "GFL"]
                if gfls:
                    s_sm = sum(u.Sn for u in sms) + sum(u.Sn for u in gfms)
                    reg["H"] *= s_sm / Sn
                    reg["KD"] *= s_sm / Sn
                    reg["mp"] *= Sn / s_sm
            elif genre == "GFM":
                reg = {**_params_gfm(membres), **_intentions(membres)}
            else:
                reg = {**_params_gfl(membres), **_intentions(membres)}
            for k, f in adj.items():
                reg[k] = float(reg.get(k, membres[0].params[k])) * f
            ppb[b] = reg
            ids[TABLE[genre]] += 1
            if genre == "SM":
                ajouts["T_SMs"].append((ids["T_SMs"], b, Sn))
            else:
                vdc = next((float(u.vdc) for u in membres if u.vdc), 40e3)
                ajouts[TABLE[genre]].append((ids[TABLE[genre]], b, Sn, vdc))
    for b in retires:
        ppb.pop(b, None)
    eq.Y_DER = [r for r in t.Y_DER if int(r[0]) not in retires] + nouveaux_y
    eq.T_Transfos = [
        r for r in t.T_Transfos if int(r[1]) not in retires and int(r[2]) not in retires
    ] + nouveaux_tr
    for table in TABLE.values():
        setattr(
            eq,
            table,
            [r for r in getattr(t, table) if int(r[1]) not in retires] + ajouts[table],
        )
    eq.T_Loads = [*t.T_Loads, *charges]
    eq.PARAMS_PAR_BUS = ppb
    return eq


def _tensions(t: Tables) -> dict[int, complex]:
    net = assembler(t)
    lf = net.run_loadflow()
    idx = net._bus_index()
    return {b: complex(lf.V[i]) for b, i in idx.items()}


@dataclass
class Plan:
    """Un plan d'agregation nomme, pour les etudes."""

    nom: str
    groupes: list[Groupe] = field(default_factory=list)
