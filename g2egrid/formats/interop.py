"""Interoperabilite : lire (et ecrire) les formats des autres outils.

========================== =========================== =====================
Format                     Lecture                     Ecriture
========================== =========================== =====================
MATPOWER (``.m``, dict)    :func:`depuis_matpower`     :func:`vers_matpower`
PSS/E RAW (v33, v34)       :func:`depuis_psse`
pandapower (JSON)          :func:`depuis_pandapower`
========================== =========================== =====================

Aucun de ces lecteurs n'a besoin de l'outil d'origine : le ``.m`` de
MATPOWER, le RAW de PSS/E et le JSON de pandapower sont lus tels quels.

CE QUI EST REPRIS -- le reseau du load flow : bus, lignes, transformateurs a
deux enroulements (rapport hors nominal compris), charges, shunts, unites et
leurs consignes. La DYNAMIQUE ne l'est pas : aucun de ces formats ne la
porte (PSS/E la met dans un ``.dyr`` a part, dont les modeles ne sont pas
ceux de G2Egrid). Les unites prennent les reglages par defaut du catalogue,
comme ieee39 ; ``unites`` choisit leur nature. Ce qui n'a pas d'equivalent
-- dephaseur, transformateur a trois enroulements -- est ignore et SIGNALE
dans ``t.avertissements`` : un import muet qui change le reseau serait pire
qu'un import qui refuse. La CONDUCTANCE SHUNT (``Gs`` de MATPOWER) est reprise
APPROCHEE : les tables n'ont pas de conductance a impedance constante (``T_Shunts``
ne porte qu'une susceptance), elle devient une charge a PUISSANCE constante
``Gs`` MW, exacte a 1 pu et fausse de ``Gs.(V^2 - 1)`` ailleurs ; l'import le
signale, avec le plus gros ``Gs`` en cause.

Conventions de G2Egrid rappelees : ``Y_DER.Type`` vaut 1 (slack), 2 (PV),
3 (PQ) ; les impedances sont en pu base systeme (100 MVA).

LES CHARGES -- ``charges="Z"`` (defaut, la convention du catalogue) : une
charge a ``P >= 0`` et ``Q > 0`` devient une R-L a impedance constante
(``Y_DER``), les autres une charge a puissance constante (``T_Loads``, type
``"P"``). ``charges="P"`` : toutes a puissance constante, ce que supposent
MATPOWER, PSS/E (part PL, QL) et pandapower ; le load flow est alors
exactement le leur. Avec ``"Z"``, il s'en ecarte d'autant que les tensions
s'ecartent de 1 pu.

Le PLANCHER ``R >= X / x_sur_r_max`` (100 par defaut) amortit les branches
sans resistance -- deux transformateurs sans pertes en parallele font une
boucle de courant non amortie a f0 exactement --, au prix d'un ecart de
load flow ; ``x_sur_r_max=None`` le retire.

CGMES (CIM) et l'export FMU ne sont pas couverts.
"""

from __future__ import annotations

import json
import math
import pathlib
import re
from collections import defaultdict

import numpy as np

from g2egrid.reseau.tables import Tables

__all__ = [
    "depuis_matpower",
    "depuis_pandapower",
    "depuis_psse",
    "disposition",
    "importer",
    "vers_matpower",
]

SB = 100e6
#: le rapport X/R le plus eleve admis : une branche sans resistance
#: n'amortit rien, et deux transformateurs sans pertes en parallele font une
#: boucle de courant non amortie a f0 exactement
X_SUR_R_MAX = 100.0
SLACK, PV, PQ = 1, 2, 3


# ----------------------------------------------------------------------
# outils communs
# ----------------------------------------------------------------------
def sn_unite(p_w: float, qmax_var: float = 0.0, mbase_va: float = 0.0) -> float:
    """La puissance nominale d'une unite : sa base si elle est donnee et
    suffit, sinon de quoi porter sa production a facteur de puissance 0,85
    (100 MVA au moins pour une machine sans production, un compensateur)."""
    if mbase_va > 0 and mbase_va >= abs(p_w):
        return float(mbase_va)
    return float(max(abs(p_w) / 0.85, abs(qmax_var), 1e6))


def disposition(buses, aretes, graine: int = 0) -> dict[int, tuple[float, float]]:
    """Des positions pour l'interface : placement par ressorts
    (Fruchterman-Reingold), ramene dans [0,05 ; 0,95]."""
    buses = list(buses)
    rng = np.random.default_rng(graine)
    n = len(buses)
    if n == 0:
        return {}
    ix = {b: k for k, b in enumerate(buses)}
    x = rng.uniform(0, 1, (n, 2))
    k = 1.0 / np.sqrt(n)
    E = [(ix[a], ix[b]) for a, b in aretes if a in ix and b in ix and a != b]
    t = 0.1
    for _ in range(500):
        d = x[:, None, :] - x[None, :, :]
        r = np.linalg.norm(d, axis=2) + 1e-9
        f = (d / r[..., None] * (k * k / r)[..., None]).sum(axis=1)
        for a, b in E:
            v = x[a] - x[b]
            ra = np.linalg.norm(v) + 1e-9
            g = v / ra * ra * ra / k
            f[a] -= g
            f[b] += g
        lf = np.linalg.norm(f, axis=1) + 1e-9
        x += f / lf[:, None] * np.minimum(lf, t)[:, None]
        t *= 0.99
    x -= x.min(axis=0)
    x /= x.max(axis=0) + 1e-12
    x = 0.05 + 0.9 * x
    return {b: (float(x[ix[b], 0]), float(x[ix[b], 1])) for b in buses}


class _Reseau:
    """Le reseau intermediaire, commun a tous les lecteurs (unites SI,
    impedances en pu base ``SB``)."""

    def __init__(self, nom: str, f0: float):
        self.nom, self.f0 = nom, f0
        self.bus: dict[int, list] = {}  # id -> [type, Un_V, V_pu, theta_rad]
        self.unites: dict[int, list] = {}  # id -> [P_W, Q_var, Sn_VA, table]
        self.charges: dict[int, list] = defaultdict(lambda: [0.0, 0.0])
        self.shunts: dict[int, float] = defaultdict(float)  # B_pu
        self.branches: list[tuple] = []  # (de, vers, R, X, B, tap)
        self.lignes_si: list[tuple] = []  # (de, vers, km, R, X, B) par km
        self.transfos: list[tuple] = []  # (de, vers, Sn, U1, U2, R, X) base propre
        self.ibs: dict[int, tuple] = {}  # bus -> (E, theta_deg)
        self.avertissements: list[str] = []

    def unite(self, b, p, q, sn, table):
        u = self.unites.get(b)
        if u is None:
            self.unites[b] = [p, q, sn, table]
        else:  # plusieurs groupes au meme bus : une unite agregee
            u[0] += p
            u[1] += q
            u[2] += sn
            if u[3] != table:
                self.avertissements.append(
                    f"bus {b} : {u[3]} et {table} agreges en {u[3]}"
                )

    def avertir(self, texte: str):
        self.avertissements.append(texte)

    def tables(
        self,
        x_min: float = 0.0,
        x_sur_r_max: float | None = X_SUR_R_MAX,
        charges: str = "Z",
    ) -> Tables:
        t = Tables(f0=self.f0, Sb=SB, nom=self.nom)
        # la reference d'angle a zero : le load flow part d'angles nuls, et
        # un slack a 50 deg (pandapower, apres un couplage Yd) le fait diverger
        slacks = [b for b, v in sorted(self.bus.items()) if v[0] == SLACK]
        th0 = self.bus[slacks[0]][3] if slacks else 0.0
        if th0:
            self.avertir(
                f"angles decales de {-math.degrees(th0):g} deg (slack {slacks[0]} a 0)"
            )
            for v in self.bus.values():
                if v[3]:  # un angle connu ; 0 = pas de donnee, depart a plat
                    v[3] -= th0
            self.ibs = {
                b: (e, th - math.degrees(th0)) for b, (e, th) in self.ibs.items()
            }
        for k in (
            "T_Lines_SI",
            "T_Lines_pu",
            "T_Shunts",
            "T_Transfos",
            "T_IBs",
            "T_SMs",
            "T_GFL1s",
            "T_Loads",
        ):
            setattr(t, k, [])
        if charges not in ("Z", "P"):
            raise ValueError(f"charges {charges!r} : 'Z' (R-L) ou 'P' (puissance)")
        rl = {
            b: c
            for b, c in self.charges.items()
            if charges == "Z" and c[0] >= 0 and c[1] > 0
        }  # R-L : P >= 0 et Q > 0
        y = []
        for b, (typ, un, v, th) in sorted(self.bus.items()):
            u = self.unites.get(b, [0.0, 0.0, 0.0, None])
            pg = u[0]  # au slack : valeur de depart ; au PV : consigne
            qg = u[1] if typ == PQ else 0.0  # un PQ tient P et Q
            pl, ql = rl.get(b, (0.0, 0.0))
            y.append((b, typ, un, v, th, pg, qg, pl, ql))
        t.Y_DER = y
        t.T_Loads = [
            (k + 1, b, c[0], c[1], "P")
            for k, (b, c) in enumerate(
                sorted(
                    (b, c)
                    for b, c in self.charges.items()
                    if b not in rl and (c[0] or c[1])
                )
            )
        ]
        for b, (p, q, sn, table) in sorted(self.unites.items()):
            rows = getattr(t, table)
            if table == "T_GFL1s":
                rows.append((len(rows) + 1, b, sn, 0.0))
            else:
                rows.append((len(rows) + 1, b, sn))
        t.T_IBs = [
            (k + 1, b, 1e10, e, th, 0.001, 0.01)
            for k, (b, (e, th)) in enumerate(sorted(self.ibs.items()))
        ]
        t.T_Shunts = [
            (k + 1, b, s) for k, (b, s) in enumerate(sorted(self.shunts.items())) if s
        ]
        lignes, transfos = [], []
        for de, vers, R, X, B, tap in self.branches:
            X = max(X, x_min)
            if x_sur_r_max:
                R = max(R, abs(X) / x_sur_r_max)
            ui, uj = self.bus[de][1], self.bus[vers][1]
            if tap == 0 and ui == uj:
                lignes.append((len(lignes) + 1, de, vers, 1, R, X, B, ui, SB))
            else:
                n = tap or 1.0
                # prise ideale cote origine, impedance cote extremite (MATPOWER
                # et PSS/E) ; l'assembleur rapporte R, X a U1 = n.Un_i : 1/n^2
                if B:
                    self.avertir(
                        f"branche {de}-{vers} : susceptance d'un transformateur ignoree"
                    )
                transfos.append(
                    (len(transfos) + 1, de, vers, SB, n * ui, uj, R / n**2, X / n**2)
                )
        for de, vers, sn, u1, u2, R, X in self.transfos:
            if x_sur_r_max:
                R = max(R, abs(X) / x_sur_r_max)
            transfos.append((len(transfos) + 1, de, vers, sn, u1, u2, R, X))
        t.T_Lines_pu, t.T_Transfos = lignes, transfos
        t.T_Lines_SI = [
            (k + 1, de, vers, km, r, x, b)
            for k, (de, vers, km, r, x, b) in enumerate(self.lignes_si)
        ]
        # un bus PV sans unite n'a rien pour tenir sa tension
        for i, r in enumerate(t.Y_DER):
            if r[1] == PV and r[0] not in self.unites:
                t.Y_DER[i] = (r[0], PQ, *r[2:])
                self.avertir(f"bus {r[0]} : PV sans unite, devient PQ")
        t.avertissements = list(self.avertissements)
        return t


# ----------------------------------------------------------------------
# MATPOWER
# ----------------------------------------------------------------------
_MATRICE = re.compile(r"mpc\.(\w+)\s*=\s*\[(.*?)\]\s*;", re.S)
_SCALAIRE = re.compile(r"mpc\.(\w+)\s*=\s*([-+0-9.eEInfa]+)\s*;")


def _lire_m(texte: str) -> dict:
    texte = "\n".join(ln.split("%", 1)[0] for ln in texte.splitlines())
    mpc: dict = {}
    for nom, v in _SCALAIRE.findall(texte):
        mpc[nom] = float(v)
    for nom, corps in _MATRICE.findall(texte):
        lignes = [ln for ln in re.split(r"[;\n]", corps) if ln.strip()]
        rows = []
        for ln in lignes:
            vals = [
                float(x.replace("Inf", "inf"))
                for x in re.split(r"[\s,]+", ln.strip())
                if x
            ]
            if vals:
                rows.append(vals)
        if rows:
            mpc[nom] = np.array(rows, dtype=float)
    return mpc


def depuis_matpower(
    source,
    nom: str | None = None,
    f0: float = 50.0,
    un_defaut: float = 138e3,
    unites: str = "T_SMs",
    x_min: float = 0.0,
    x_sur_r_max: float | None = X_SUR_R_MAX,
    charges: str = "Z",
) -> Tables:
    """Un cas MATPOWER : chemin ou texte d'un ``.m``, ou dict ``mpc`` (cles
    ``baseMVA``, ``bus``, ``gen``, ``branch`` -- celui de PYPOWER convient).
    ``un_defaut`` : tension de base des bus dont ``BASE_KV`` vaut 0 (IEEE57).
    ``unites`` : table des generateurs (``"T_SMs"`` par defaut)."""
    if isinstance(source, dict):
        mpc = {k: np.asarray(v, dtype=float) for k, v in source.items() if k in (
            "bus", "gen", "branch")}  # fmt: skip
        mpc["baseMVA"] = float(source.get("baseMVA", 100.0))
    else:
        texte = str(source)
        p = pathlib.Path(texte) if len(texte) < 4096 and "\n" not in texte else None
        if p is not None and p.exists():
            nom = nom or p.stem
            texte = p.read_text(encoding="utf-8", errors="replace")
        mpc = _lire_m(texte)
    for cle in ("bus", "branch"):
        if cle not in mpc:
            raise ValueError(f"MATPOWER : mpc.{cle} absent")
    base = float(mpc.get("baseMVA", 100.0))
    kz = SB / (base * 1e6)  # pu base du cas -> pu base 100 MVA
    r = _Reseau(nom or "matpower", f0)
    typ = {3: SLACK, 2: PV, 1: PQ}
    sans_base = 0
    avec_gs: list[tuple[float, int]] = []
    for row in mpc["bus"]:
        b, t_, pd, qd, gs, bs = int(row[0]), int(row[1]), *row[2:6]
        if t_ == 4:
            continue  # bus isole
        un = row[9] * 1e3 if len(row) > 9 and row[9] > 0 else un_defaut
        sans_base += not (len(row) > 9 and row[9] > 0)
        r.bus[b] = [typ.get(t_, PQ), un, float(row[7]), math.radians(row[8])]
        if pd or qd:
            r.charges[b][0] += pd * 1e6
            r.charges[b][1] += qd * 1e6
        if gs:  # conductance shunt : P constante (approchee, voir l'en-tete)
            r.charges[b][0] += gs * 1e6
            avec_gs.append((abs(gs), b))
        if bs:
            r.shunts[b] += bs / base
    if avec_gs:
        gmax, bmax = max(avec_gs)
        r.avertir(
            f"{len(avec_gs)} bus avec conductance shunt Gs (le plus gros : "
            f"{gmax:g} MW au bus {bmax}) : reprise comme charge a puissance "
            "constante, MATPOWER la traite en impedance constante (ecart "
            "Gs.(V^2 - 1))"
        )
    if sans_base:
        r.avertir(f"{sans_base} bus sans BASE_KV : {un_defaut / 1e3:g} kV pris")
    for row in mpc.get("gen", np.zeros((0, 10))):
        b = int(row[0])
        if len(row) > 7 and row[7] <= 0 or b not in r.bus:
            continue
        p, q = row[1] * 1e6, row[2] * 1e6
        mb = row[6] * 1e6 if len(row) > 6 else 0.0
        r.unite(b, p, q, sn_unite(p, row[3] * 1e6, mb), unites)
        r.bus[b][2] = float(row[5])
    dephases = 0
    for row in mpc["branch"]:
        if len(row) > 10 and row[10] <= 0:
            continue
        de, vers = int(row[0]), int(row[1])
        if de not in r.bus or vers not in r.bus:
            continue
        if len(row) > 9 and row[9]:
            dephases += 1
        r.branches.append(
            (
                de,
                vers,
                row[2] * kz,
                row[3] * kz,
                row[4] / kz,
                row[8] if len(row) > 8 else 0,
            )
        )
    if dephases:
        r.avertir(f"{dephases} dephaseur(s) : angle de dephasage ignore")
    return r.tables(x_min, x_sur_r_max, charges)


def vers_matpower(t: Tables, chemin=None) -> str:
    """Le texte d'un cas MATPOWER (``mpc`` version 2, base ``t.Sb``) ; ecrit
    dans ``chemin`` s'il est donne. Toutes les unites deviennent des
    generateurs, les bus infinis aussi ; les lignes en unites SI sont
    ramenees en pu."""
    Sb = float(t.Sb)
    base = Sb / 1e6
    UN = {int(r[0]): float(r[2]) for r in t.Y_DER}
    typ = {SLACK: 3, PV: 2, PQ: 1}
    charges = {int(r[0]): [r[7], r[8]] for r in t.Y_DER}
    for _id, b, p, q, _ty in getattr(t, "T_Loads", []):
        charges[int(b)][0] += p
        charges[int(b)][1] += q
    shunts = defaultdict(float)
    for _id, b, B in getattr(t, "T_Shunts", []):
        shunts[int(b)] += B
    bus = []
    for b, ty, un, v, th, *_ in t.Y_DER:
        b = int(b)
        pd, qd = charges[b]
        bus.append(
            [b, typ[int(ty)], pd / 1e6, qd / 1e6, 0, shunts[b] * base, 1, v,
             math.degrees(th), un / 1e3, 1, 1.1, 0.9]
        )  # fmt: skip
    pg = {int(r[0]): (r[5], r[6], r[3]) for r in t.Y_DER}
    gen = []
    for table in ("T_IBs", "T_SMs", "T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s",
                  "T_GFL1s", "T_GFL2s"):  # fmt: skip
        for row in getattr(t, table, []):
            b, sn = int(row[1]), float(row[2])
            p, q, v = pg.get(b, (0, 0, 1))
            gen.append(
                [b, p / 1e6, q / 1e6, 9999, -9999, v, min(sn, 1e5 * Sb) / 1e6, 1,
                 max(sn, abs(p)) / 1e6, 0]
            )  # fmt: skip
    branch = []
    for _id, i, j, km, r_, x_, b_, un, sbl in getattr(t, "T_Lines_pu", []):
        k = (Sb / sbl) * (un / UN[int(i)]) ** 2
        branch.append(
            [int(i), int(j), r_ * k, x_ * k, b_ / k, 0, 0, 0, 0, 0, 1, -360, 360]
        )
    for _id, i, j, km, r_, x_, b_ in getattr(t, "T_Lines_SI", []):
        zb = UN[int(i)] ** 2 / Sb
        branch.append(
            [int(i), int(j), r_ * km / zb, x_ * km / zb, b_ * km * zb, 0, 0, 0, 0,
             0, 1, -360, 360]
        )  # fmt: skip
    for _id, i, j, sn, u1, u2, r_, x_ in getattr(t, "T_Transfos", []):
        i, j = int(i), int(j)
        k = (Sb / sn) * (u1 / UN[i]) ** 2
        n = (u1 / u2) / (UN[i] / UN[j])
        # l'assembleur place R.k cote extremite, apres la prise ideale n cote
        # origine : c'est exactement la branche MATPOWER (ratio = n)
        branch.append([i, j, r_ * k, x_ * k, 0, 0, 0, 0, n, 0, 1, -360, 360])

    def matrice(nomm, rows, entete):
        corps = "\n".join("\t" + "\t".join(f"{v:.10g}" for v in r) + ";" for r in rows)
        return f"%% {entete}\nmpc.{nomm} = [\n{corps}\n];\n"

    fonction = re.sub(r"\W", "_", t.nom or "cas")
    texte = (
        f"function mpc = {fonction}\n"
        f"%% exporte par G2Egrid (f0 = {t.f0:g} Hz)\n"
        "mpc.version = '2';\n"
        f"mpc.baseMVA = {base:g};\n\n"
        + matrice("bus", bus, "bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin")
        + "\n"
        + matrice("gen", gen, "bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin")
        + "\n"
        + matrice("branch", branch, "fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax")
    )  # fmt: skip
    if chemin is not None:
        pathlib.Path(chemin).write_text(texte, encoding="utf-8")
    return texte


# ----------------------------------------------------------------------
# PSS/E RAW
# ----------------------------------------------------------------------
_SECTIONS = {
    "BUS": "bus",
    "LOAD": "load",
    "FIXED SHUNT": "shunt",
    "GENERATOR": "gen",
    "BRANCH": "branch",
    "NON-TRANSFORMER BRANCH": "branch",
    "SYSTEM SWITCHING DEVICE": "switch",
    "TRANSFORMER": "transfo",
}


def _champs(ligne: str) -> list[str]:
    """Une ligne RAW decoupee : virgules hors apostrophes, commentaire ``/``
    hors apostrophes retire."""
    out, cur, dedans = [], [], False
    for c in ligne:
        if c == "'":
            dedans = not dedans
            continue
        if c == "/" and not dedans:
            break
        if c == "," and not dedans:
            out.append("".join(cur).strip())
            cur = []
            continue
        cur.append(c)
    if "".join(cur).strip():
        out.append("".join(cur).strip())
    return out


def _sections_raw(lignes: list[str], rev: int) -> dict[str, list[list[str]]]:
    ordre = (
        ["bus", "load", "shunt", "gen", "branch", "transfo"]
        if rev < 34
        else ["bus", "load", "shunt", "gen", "branch", "switch", "transfo"]
    )
    out: dict[str, list[list[str]]] = defaultdict(list)
    courant, k = "bus", 0
    for ln in lignes:
        corps = ln.split("/")[0].strip()
        if corps == "Q":
            break
        if corps == "0" or re.fullmatch(r"0\s*", corps) or corps.startswith("0 "):
            m = re.search(r"BEGIN\s+(.+?)\s+DATA", ln.upper())
            k += 1
            if m:
                courant = _SECTIONS.get(m.group(1).strip(), m.group(1).strip().lower())
            else:
                courant = ordre[k] if k < len(ordre) else f"section{k}"
            continue
        if corps:
            out[courant].append(_champs(ln))
    return out


def depuis_psse(
    source,
    nom: str | None = None,
    unites: str = "T_SMs",
    x_min: float = 0.0,
    x_sur_r_max: float | None = X_SUR_R_MAX,
    charges: str = "Z",
) -> Tables:
    """Un fichier PSS/E RAW (versions 33 et 34) : chemin ou texte. Bus,
    charges, shunts fixes, generateurs, branches et transformateurs a deux
    enroulements (codes CW 1-3, CZ 1-3). ``x_min`` [pu base 100 MVA] : la
    reactance minimale des branches (bretelles)."""
    texte = str(source)
    p = pathlib.Path(texte) if len(texte) < 4096 and "\n" not in texte else None
    if p is not None and p.exists():
        nom = nom or p.stem
        texte = p.read_text(encoding="utf-8", errors="replace")
    lignes = texte.splitlines()
    tete = _champs(lignes[0])
    sbase = float(tete[1]) * 1e6
    rev = int(float(tete[2])) if len(tete) > 2 and tete[2] else 33
    f0 = float(tete[5]) if len(tete) > 5 and tete[5] else 50.0
    if rev not in (33, 34):
        raise ValueError(f"PSS/E RAW version {rev} : seules 33 et 34 sont lues")
    kz = SB / sbase
    S = _sections_raw(lignes[3:], rev)
    r = _Reseau(nom or "psse", f0)
    typ = {1: PQ, 2: PV, 3: SLACK}
    for row in S["bus"]:
        if int(row[3]) == 4:
            continue
        r.bus[int(row[0])] = [
            typ.get(int(row[3]), PQ),
            float(row[2]) * 1e3,
            float(row[7]),
            math.radians(float(row[8])),
        ]
    for row in S["load"]:
        b = int(row[0])
        if int(row[2]) and b in r.bus:
            r.charges[b][0] += float(row[5]) * 1e6
            r.charges[b][1] += float(row[6]) * 1e6
            # charges a courant (IP, IQ) et admittance (YP, YQ) : prises a V = 1
            if len(row) > 10:
                r.charges[b][0] += (float(row[7] or 0) + float(row[9] or 0)) * 1e6
                r.charges[b][1] += (float(row[8] or 0) - float(row[10] or 0)) * 1e6
    for row in S["shunt"]:
        b = int(row[0])
        if int(row[2]) and b in r.bus:
            if float(row[3]):
                r.charges[b][0] += float(row[3]) * 1e6
            r.shunts[b] += float(row[4]) * 1e6 / SB
    for row in S["gen"]:
        b = int(row[0])
        if not int(row[14]) or b not in r.bus:
            continue
        pg, qg, qmax, mb = (float(row[i]) * 1e6 for i in (2, 3, 4, 8))
        r.unite(b, pg, qg, sn_unite(pg, qmax, mb), unites)
        r.bus[b][2] = float(row[6])
    ist = 13 if rev < 34 else 23
    for row in S["branch"]:
        i, j = int(row[0]), abs(int(row[1]))
        if len(row) > ist and not int(row[ist]):
            continue
        if i in r.bus and j in r.bus:
            R, X, B = float(row[3]) * kz, float(row[4]) * kz, float(row[5]) / kz
            r.branches.append((i, j, R, X, B, 0.0))
    tr = S["transfo"]
    k, trois = 0, 0
    while k < len(tr):
        l1 = tr[k]
        if int(l1[2] or 0):  # trois enroulements : 5 lignes
            trois += 1
            k += 5
            continue
        l2, l3, l4 = tr[k + 1 : k + 4]
        k += 4
        if not int(l1[11]):
            continue
        i, j = int(l1[0]), int(l1[1])
        if i not in r.bus or j not in r.bus:
            continue
        cw, cz = int(l1[4]), int(l1[5])
        R, X, sb12 = float(l2[0]), float(l2[1]), float(l2[2] or 0) * 1e6
        if cz == 3:  # pertes cuivre [W] et |Z| [pu base propre]
            R, X = R / sb12, math.sqrt(max(X**2 - (R / sb12) ** 2, 0.0))
        if cz in (2, 3):
            R, X = R * SB / sb12, X * SB / sb12
        else:
            R, X = R * kz, X * kz
        w1, w2 = float(l3[0]), float(l4[0])
        if cw == 2:  # en kV
            w1, w2 = w1 * 1e3 / r.bus[i][1], w2 * 1e3 / r.bus[j][1]
        elif cw == 3:  # en pu de la tension nominale de l'enroulement
            n1 = float(l3[1] or 0) * 1e3 or r.bus[i][1]
            n2 = float(l4[1] or 0) * 1e3 or r.bus[j][1]
            w1, w2 = w1 * n1 / r.bus[i][1], w2 * n2 / r.bus[j][1]
        if float(l3[2] or 0):
            r.avertir(f"transformateur {i}-{j} : dephasage ignore")
        r.branches.append((i, j, R, X, 0.0, w1 / w2))
    if trois:
        r.avertir(f"{trois} transformateur(s) a trois enroulements ignore(s)")
    for b in r.bus:
        if r.bus[b][0] != PQ and b not in r.unites:
            r.bus[b][0] = PQ
    return r.tables(x_min, x_sur_r_max, charges)


# ----------------------------------------------------------------------
# pandapower
# ----------------------------------------------------------------------
def _tableau_pp(obj) -> list[dict]:
    """Un DataFrame serialise par pandapower (``orient = split``) en lignes."""
    if obj is None:
        return []
    if isinstance(obj, dict) and "_object" in obj:
        obj = obj["_object"]
    if isinstance(obj, str):
        obj = json.loads(obj)
    cols = obj.get("columns", [])
    return [
        dict(zip(cols, row), index=idx)
        for idx, row in zip(obj.get("index", []), obj.get("data", []))
    ]


def _vrai(v) -> bool:
    return bool(v) if v is not None else True


def _num(v, defaut=0.0) -> float:
    try:
        x = float(v)
    except (TypeError, ValueError):
        return defaut
    return defaut if math.isnan(x) else x


def depuis_pandapower(
    source,
    nom: str | None = None,
    sgen: str = "T_GFL1s",
    x_sur_r_max: float | None = X_SUR_R_MAX,
    charges: str = "Z",
) -> Tables:
    """Un reseau pandapower : chemin ou texte de son JSON (``pp.to_json``),
    ou l'objet ``pandapowerNet`` lui-meme. Bus, lignes, transformateurs a
    deux enroulements, charges, shunts, generateurs (``gen`` : machines
    synchrones PV), generateurs statiques (``sgen`` : ``T_GFL1s`` par
    defaut, a puissance imposee), reseau amont (``ext_grid`` : bus infini,
    slack) et interrupteurs (bus-bus fermes : bus fusionnes ; ouverts sur
    ligne ou transformateur : branche hors service)."""
    if hasattr(source, "bus") and hasattr(source, "line"):  # pandapowerNet
        import pandapower as pp  # noqa: PLC0415

        source = pp.to_json(source)
    texte = str(source)
    p = pathlib.Path(texte) if len(texte) < 4096 and "\n" not in texte else None
    if p is not None and p.exists():
        nom = nom or p.stem.replace(".json", "")
        texte = p.read_text(encoding="utf-8")
    d = json.loads(texte)
    o = d.get("_object", d)
    T = {k: _tableau_pp(o.get(k)) for k in (
        "bus", "line", "trafo", "trafo3w", "load", "sgen", "gen", "ext_grid",
        "shunt", "switch", "impedance", "xward", "ward")}  # fmt: skip
    f0 = _num(o.get("f_hz"), 50.0)
    r = _Reseau(nom or str(o.get("name") or "pandapower"), f0)
    sn_mva = _num(o.get("sn_mva"), 1.0)
    del sn_mva  # les puissances sont lues en MW, la base est celle de G2Egrid

    # interrupteurs : bus-bus fermes -> bus fusionnes ; ouverts -> hors service
    parent: dict[int, int] = {}

    def racine(b):
        while parent.get(b, b) != b:
            b = parent[b]
        return b

    # element -> bus ou il est ouvert (un interrupteur de ligne ne coupe que
    # SON extremite : la ligne reste sous tension depuis l'autre)
    ouverts: dict[str, dict[int, set]] = {"l": defaultdict(set), "t": defaultdict(set)}
    for s in T["switch"]:
        ferme = _vrai(s.get("closed"))
        if s.get("et") == "b" and ferme:
            a, c = racine(int(s["bus"])), racine(int(s["element"]))
            if a != c:
                parent[max(a, c)] = min(a, c)
        elif s.get("et") in ("l", "t") and not ferme:
            ouverts[s["et"]][int(s["element"])].add(int(s["bus"]))

    # identifiants G2Egrid : index pandapower + 1
    def ident(b):
        return racine(int(b)) + 1

    for b in T["bus"]:
        if not _vrai(b.get("in_service")):
            continue
        i = ident(b["index"])
        if i not in r.bus:
            r.bus[i] = [PQ, _num(b.get("vn_kv")) * 1e3, 1.0, 0.0]
    w = 2 * math.pi * f0
    for ln in T["line"]:
        if not _vrai(ln.get("in_service")):
            continue
        i, j = ident(ln["from_bus"]), ident(ln["to_bus"])
        par = max(_num(ln.get("parallel"), 1.0), 1.0)
        km = _num(ln.get("length_km"), 1.0)
        r_, x_ = _num(ln.get("r_ohm_per_km")) / par, _num(ln.get("x_ohm_per_km")) / par
        b_ = w * _num(ln.get("c_nf_per_km")) * 1e-9 * par
        coupe = ouverts["l"].get(ln["index"], set())
        if coupe:
            bouts = {int(ln["from_bus"]), int(ln["to_bus"])} - coupe
            if len(bouts) == 1:
                # ligne a vide : son admittance d'entree (schema en pi exact)
                # devient un shunt au bus qui l'alimente
                k = ident(bouts.pop())
                zb = r.bus[k][1] ** 2 / SB
                z = complex(r_, x_) * km / zb
                yb = 0.5j * b_ * km * zb
                if yb:
                    r.shunts[k] += (yb + 1 / (z + 1 / yb)).imag
            continue
        if _num(ln.get("g_us_per_km")):
            r.avertir(f"ligne {ln['index']} : conductance ignoree")
        if i == j:
            continue
        if r.bus[i][1] != r.bus[j][1]:
            r.avertir(f"ligne {ln['index']} entre deux tensions de base : base de {i}")
        r.lignes_si.append((i, j, km, r_, x_, b_))
    for tr in T["trafo"]:
        if not _vrai(tr.get("in_service")) or tr["index"] in ouverts["t"]:
            # un transformateur ouvert d'un cote ne transmet rien (sa
            # magnetisation, ignoree ici, est tout ce qui resterait)
            continue
        hv, lv = ident(tr["hv_bus"]), ident(tr["lv_bus"])
        sn = _num(tr.get("sn_mva")) * 1e6
        par = max(_num(tr.get("parallel"), 1.0), 1.0)
        z, rr = _num(tr.get("vk_percent")) / 100, _num(tr.get("vkr_percent")) / 100
        x = math.sqrt(max(z * z - rr * rr, 0.0))
        u1, u2 = _num(tr.get("vn_hv_kv")) * 1e3, _num(tr.get("vn_lv_kv")) * 1e3
        pas = _num(tr.get("tap_step_percent")) / 100
        prise = _num(tr.get("tap_pos")) - _num(tr.get("tap_neutral"))
        if pas and prise:
            if tr.get("tap_side") == "lv":
                u2 *= 1 + prise * pas
            else:
                u1 *= 1 + prise * pas
        if _num(tr.get("shift_degree")) % 360:
            r.avertir(
                f"transformateur {tr['index']} : couplage "
                f"{_num(tr.get('shift_degree')):g} deg ignore (angles decales)"
            )
        # pandapower rapporte l'impedance au cote BT (prise BT comprise) ;
        # l'assembleur, a U1 (prise HT comprise) : on la ramene
        k = ((u2 / r.bus[lv][1]) / (u1 / r.bus[hv][1])) ** 2
        r.transfos.append((hv, lv, sn * par, u1, u2, rr * k, x * k))
    if T["trafo3w"]:
        r.avertir(
            f"{len(T['trafo3w'])} transformateur(s) a trois enroulements ignore(s)"
        )
    for k in ("impedance", "ward", "xward"):
        if T[k]:
            r.avertir(f"{len(T[k])} element(s) '{k}' ignore(s)")
    for c in T["load"]:
        if _vrai(c.get("in_service")):
            i, s = ident(c["bus"]), _num(c.get("scaling"), 1.0)
            r.charges[i][0] += _num(c.get("p_mw")) * s * 1e6
            r.charges[i][1] += _num(c.get("q_mvar")) * s * 1e6
    for sh in T["shunt"]:
        if not _vrai(sh.get("in_service")):
            continue
        i = ident(sh["bus"])
        k = (
            _num(sh.get("step"), 1.0)
            * (r.bus[i][1] / (_num(sh.get("vn_kv")) * 1e3 or r.bus[i][1])) ** 2
        )
        r.shunts[i] += -_num(sh.get("q_mvar")) * k * 1e6 / SB
        if _num(sh.get("p_mw")):
            r.charges[i][0] += _num(sh.get("p_mw")) * k * 1e6
    for e in T["ext_grid"]:
        if _vrai(e.get("in_service")):
            i = ident(e["bus"])
            v, th = _num(e.get("vm_pu"), 1.0), _num(e.get("va_degree"))
            r.bus[i] = [SLACK, r.bus[i][1], v, math.radians(th)]
            r.ibs[i] = (v, th)
    for g in T["gen"]:
        if not _vrai(g.get("in_service")):
            continue
        i = ident(g["bus"])
        p = _num(g.get("p_mw")) * _num(g.get("scaling"), 1.0) * 1e6
        sn = _num(g.get("sn_mva")) * 1e6
        r.unite(i, p, 0.0, sn_unite(p, _num(g.get("max_q_mvar")) * 1e6, sn), "T_SMs")
        if r.bus[i][0] != SLACK:
            r.bus[i][0] = SLACK if _vrai(g.get("slack")) and g.get("slack") else PV
        r.bus[i][2] = _num(g.get("vm_pu"), 1.0)
    for g in T["sgen"]:
        if not _vrai(g.get("in_service")):
            continue
        i = ident(g["bus"])
        s = _num(g.get("scaling"), 1.0)
        p, q = _num(g.get("p_mw")) * s * 1e6, _num(g.get("q_mvar")) * s * 1e6
        sn = _num(g.get("sn_mva")) * 1e6 or abs(complex(p, q)) / 0.9
        if r.bus[i][0] != PQ or i in r.unites:
            # sur un bus deja tenu, une injection fixe : une charge negative
            r.charges[i][0] -= p
            r.charges[i][1] -= q
            continue
        r.unite(i, p, q, max(sn, 1e3), sgen)
    if not r.ibs and not any(b[0] == SLACK for b in r.bus.values()):
        raise ValueError("pandapower : ni ext_grid ni gen slack")
    return r.tables(x_sur_r_max=x_sur_r_max, charges=charges)


# ----------------------------------------------------------------------
def importer(chemin, **kw) -> Tables:
    """Choisit le lecteur d'apres l'extension : ``.m`` (MATPOWER), ``.raw``
    (PSS/E), ``.json`` (pandapower ; un ``.g2egrid.json`` se lit avec
    :func:`g2egrid.reseau.catalogue.charger`)."""
    p = pathlib.Path(chemin)
    ext = p.suffix.lower()
    if ext == ".m":
        return depuis_matpower(p, **kw)
    if ext == ".raw":
        return depuis_psse(p, **kw)
    if ext == ".json":
        return depuis_pandapower(p, **kw)
    raise ValueError(f"{p.name} : format inconnu (.m, .raw, .json)")
