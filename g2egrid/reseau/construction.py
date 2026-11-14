"""Construire un reseau EN CODE, avec les conventions du constructeur graphique.

Le constructeur de l'interface et ce module font la meme chose -- ajouter des
lignes aux tables -- et suivent les memes conventions, pour qu'un reseau
dessine a la souris et un reseau ecrit en code soient interchangeables :

- un bus nu est PQ ;
- poser une machine ou un convertisseur sur un bus lui cree SON bus terminal
  et SON transformateur de raccordement, comme dans tous les reseaux de test ;
- une charge "classique" (impedance constante) va dans les colonnes Pl/Ql de
  ``Y_DER`` ; une charge a courant ou puissance constante va dans ``T_Loads`` ;
- les numeros de bus alternatifs et continus partagent le meme espace.

Tout agit EN PLACE sur un :class:`~g2egrid.reseau.tables.Tables` et rend le numero ou
l'identifiant cree, pour enchainer ::

    from g2egrid.reseau.construction import *
    t = reseau_vide("demo")
    a = poser_bus(t, 230e3)
    b = poser_bus(t, 230e3, Pl_W=100e6, Ql_var=30e6)
    poser_ligne(t, a, b, R_pu=0.01, X_pu=0.1)
    poser_equipement(t, "IB", a, slack=True)
    poser_equipement(t, "GFM1", b, Sn_VA=50e6, P_W=40e6)

Le point de fonctionnement (qui est slack, PV, PQ, avec quelles consignes) se
regle ensuite par :func:`fixer_point` -- c'est l'onglet "Point de
fonctionnement" de l'interface.
"""

from __future__ import annotations

from g2egrid.reseau.tables import TABLES, Tables

__all__ = [
    "FAMILLES",
    "SLACK",
    "PV",
    "PQ",
    "reseau_vide",
    "nouveau_numero",
    "poser_bus",
    "poser_ligne",
    "poser_transfo",
    "poser_equipement",
    "poser_charge",
    "poser_shunt",
    "equipements",
    "fixer_point",
    "un_slack_tout_pq",
]

SLACK, PV, PQ = 1, 2, 3

#: famille -> table ; ce sont les memes noms que la palette de l'interface
FAMILLES = {
    "SM": "T_SMs",
    "GFM1": "T_GFM1s",
    "GFM2": "T_GFM2s",
    "GFM3": "T_GFM3s",
    "GFM4": "T_GFM4s",
    "GFL1": "T_GFL1s",
    "GFL2": "T_GFL2s",
    "IB": "T_IBs",
}


def reseau_vide(nom: str = "reseau", f0: float = 50.0, Sb: float = 100e6) -> Tables:
    t = Tables(f0=f0, Sb=Sb, nom=nom)
    for k in TABLES:
        setattr(t, k, [])
    return t


def nouveau_numero(t: Tables) -> int:
    ids = [int(r[0]) for r in t.Y_DER] + [int(r[0]) for r in t.T_BusDC]
    return max(ids) + 1 if ids else 1


def _nouvel_id(rows) -> int:
    return max((int(r[0]) for r in rows), default=0) + 1


def _ligne_bus(t: Tables, b: int) -> int:
    for i, r in enumerate(t.Y_DER):
        if int(r[0]) == int(b):
            return i
    raise KeyError(f"bus {b} absent de Y_DER")


def poser_bus(
    t: Tables,
    Un_V: float,
    type: int = PQ,
    V_pu: float = 1.0,
    Pl_W: float = 0.0,
    Ql_var: float = 0.0,
    numero: int | None = None,
) -> int:
    b = nouveau_numero(t) if numero is None else int(numero)
    t.Y_DER.append(
        (b, type, float(Un_V), float(V_pu), 0, 0, 0, float(Pl_W), float(Ql_var))
    )
    return b


def poser_ligne(
    t: Tables,
    de: int,
    vers: int,
    R_pu: float,
    X_pu: float,
    B_pu: float = 0.0,
    longueur_km: float = 1.0,
) -> int:
    """Une ligne en pu, sur la base du reseau (``Un`` du bus de depart, ``Sb``)."""
    un = t.Y_DER[_ligne_bus(t, de)][2]
    i = _nouvel_id(t.T_Lines_pu)
    t.T_Lines_pu.append(
        (i, int(de), int(vers), float(longueur_km), R_pu, X_pu, B_pu, un, t.Sb)
    )
    return i


def poser_transfo(
    t: Tables, de: int, vers: int, Sn_VA: float, R_pu: float = 0.0, X_pu: float = 0.05
) -> int:
    u1 = t.Y_DER[_ligne_bus(t, de)][2]
    u2 = t.Y_DER[_ligne_bus(t, vers)][2]
    i = _nouvel_id(t.T_Transfos)
    t.T_Transfos.append((i, int(de), int(vers), float(Sn_VA), u1, u2, R_pu, X_pu))
    return i


def poser_equipement(
    t: Tables,
    famille: str,
    bus: int,
    Sn_VA: float | None = None,
    P_W: float = 0.0,
    Q_var: float = 0.0,
    slack: bool = False,
    type: int | None = None,
    U_V: float = 20e3,
    X_transfo: float = 0.05,
    Vdc_V: float = 40e3,
) -> int:
    """Pose une unite derriere SON transfo, sur SON bus terminal ; rend ce bus.

    Un bus infini (``IB``) garde la tension de son hote. Le type du bus
    terminal vaut ``type`` s'il est donne, ``SLACK`` si ``slack``, ``PV`` sinon.
    """
    if famille not in FAMILLES:
        raise ValueError(
            f"famille {famille!r} inconnue ; attendu l'une de {list(FAMILLES)}"
        )
    Sn = float(Sn_VA if Sn_VA is not None else t.Sb)
    hote = t.Y_DER[_ligne_bus(t, bus)]
    un = hote[2] if famille == "IB" else float(U_V)
    ty = type if type is not None else (SLACK if slack else PV)
    b = poser_bus(t, un, type=ty)
    i = _ligne_bus(t, b)
    r = list(t.Y_DER[i])
    r[5], r[6] = float(P_W), float(Q_var)
    t.Y_DER[i] = tuple(r)
    table = FAMILLES[famille]
    rows = getattr(t, table)
    k = _nouvel_id(rows)
    if famille == "IB":
        rows.append((k, b, Sn, 1.0, 0.0, 0.001, 0.05))
    elif famille == "SM":
        rows.append((k, b, Sn))
    else:
        rows.append((k, b, Sn, float(Vdc_V)))
    poser_transfo(t, bus, b, Sn, 0.0, X_transfo)
    return b


def poser_charge(
    t: Tables, bus: int, P_W: float, Q_var: float = 0.0, genre: str = "Z"
) -> None:
    """``genre`` : ``Z`` impedance constante (dans Y_DER), ``P`` ou ``I`` (T_Loads)."""
    if genre == "Z":
        i = _ligne_bus(t, bus)
        r = list(t.Y_DER[i])
        r[7] += float(P_W)
        r[8] += float(Q_var)
        t.Y_DER[i] = tuple(r)
    elif genre in ("P", "I"):
        t.T_Loads.append(
            (_nouvel_id(t.T_Loads), int(bus), float(P_W), float(Q_var), genre)
        )
    else:
        raise ValueError(f"genre de charge {genre!r} : attendu 'Z', 'P' ou 'I'")


def poser_shunt(t: Tables, bus: int, B_pu: float) -> int:
    i = _nouvel_id(t.T_Shunts)
    t.T_Shunts.append((i, int(bus), float(B_pu)))
    return i


def equipements(t: Tables) -> list[dict]:
    """Les unites du reseau : famille, identifiant, bus, et consignes de leur bus."""
    out = []
    for fam, table in FAMILLES.items():
        for r in getattr(t, table):
            b = int(r[1])
            y = t.Y_DER[_ligne_bus(t, b)]
            out.append(
                {
                    "famille": fam,
                    "table": table,
                    "id": int(r[0]),
                    "bus": b,
                    "Sn_VA": float(r[2]),
                    "type": int(y[1]),
                    "V_pu": float(y[3]),
                    "P_W": float(y[5]),
                    "Q_var": float(y[6]),
                }
            )
    return out


def fixer_point(
    t: Tables,
    bus: int,
    type: int | None = None,
    V_pu: float | None = None,
    P_W: float | None = None,
    Q_var: float | None = None,
) -> None:
    """Regle le type et les consignes d'un bus : le point de fonctionnement."""
    i = _ligne_bus(t, bus)
    r = list(t.Y_DER[i])
    if type is not None:
        r[1] = int(type)
    if V_pu is not None:
        r[3] = float(V_pu)
    if P_W is not None:
        r[5] = float(P_W)
    if Q_var is not None:
        r[6] = float(Q_var)
    t.Y_DER[i] = tuple(r)


def un_slack_tout_pq(
    t: Tables, slack_bus: int, lf_bus: dict[int, complex] | None = None
) -> None:
    """Le cas le plus courant : UN equipement en slack, tous les autres en PQ.

    Passer un PV en PQ demande de lui donner un reactif. Avec ``lf_bus`` (bus
    -> injection nette en pu, typiquement le resultat d'un load flow), le
    reactif est celui que le bus fournissait : le point de fonctionnement ne
    bouge pas, seule change la facon de le decrire. Sans, il vaut zero.

    L'ancien slack, lui, n'avait pas de consigne de puissance active : il
    recoit aussi celle qu'il fournissait au load flow. Sans ``lf_bus``, sa
    consigne ``P`` reste celle des tables.
    """
    for r in list(t.Y_DER):
        b = int(r[0])
        if b == int(slack_bus):
            fixer_point(t, b, type=SLACK)
        elif int(r[1]) in (SLACK, PV):
            q = p = None
            if lf_bus is not None and b in lf_bus:
                q = lf_bus[b].imag * t.Sb + r[8]  # Qg = Q net + Ql
                if int(r[1]) == SLACK:
                    p = lf_bus[b].real * t.Sb + r[7]  # Pg = P net + Pl
            fixer_point(t, b, type=PQ, P_W=p, Q_var=q if q is not None else 0.0)
