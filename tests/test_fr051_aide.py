"""Aides communes des tests FR-051 (editeur de blocs), sans test propre.

Ecrites a partir du SEUL cahier des charges
(``.claude/equipe/cahiers/FR-051_editeur_blocs.md``) : les tests visent le
contrat, pas l'implementation. Ce module ne contient que des outils de mesure
INDEPENDANTS du code teste :

- la fonction de transfert d'un schema, calculee sur son modele LINEARISE
  (jacobiens symboliques des equations compilees, resolution en precision
  etendue : la tolerance de 1e-12 du cahier se juge sur la formule, pas sur le
  conditionnement de la resolution) ;
- la valeur statique d'un schema sans etat ;
- l'empreinte numerique d'un schema compile (equivalence symbolique testee en
  evaluant f, g, h en points pseudo-aleatoires fixes) ;
- la pose d'un element a blocs sur un reseau, le calcul d'un spectre, la
  comparaison de deux spectres.
"""

from __future__ import annotations

import copy
import hashlib

import mpmath as mp
import numpy as np
import sympy as sp
from g2egrid.tables import Tables, assembler
from scipy.optimize import linear_sum_assignment

from g2egrid.elements import blocs as B

mp.mp.dps = 40


# ------------------------------------------------------------------ schemas
def banc_bloc(bloc: dict, parametres: dict | None = None, entrees: int = 1) -> dict:
    """Un banc minimal : une entree ``vd`` (symbole ``vd_t``), le bloc ``b`` a
    tester, sa sortie sur ``id`` ; ``iq`` vaut 0 (constante)."""
    blocs = [
        # le compilateur exige deux consignes calees (v1 comme v2), inutilisees ici
        {"id": "pr", "type": "consigne", "nom": "P_ref"},
        {"id": "qr", "type": "consigne", "nom": "Q_ref"},
        {"id": "e", "type": "entree", "signal": "vd"},
        dict(bloc, id="b"),
        {"id": "oi", "type": "sortie", "signal": "id"},
        {"id": "z", "type": "constante", "valeur": 0},
        {"id": "oq", "type": "sortie", "signal": "iq"},
    ]
    liens = [
        {"de": "e", "vers": "b", "entree": 0},
        {"de": "b", "vers": "oi"},
        {"de": "z", "vers": "oq"},
    ]
    for k in range(1, entrees):  # entrees supplementaires : constantes nulles
        blocs.append({"id": f"c{k}", "type": "constante", "valeur": 0})
        liens.append({"de": f"c{k}", "vers": "b", "entree": k})
    return {
        "nom": "banc",
        "parametres": dict(parametres or {}),
        "blocs": blocs,
        "liens": liens,
    }


def _mat(lignes: int, colonnes: int, entrees) -> mp.matrix:
    m = mp.matrix(lignes, colonnes)
    for (i, j), v in entrees.items():
        m[i, j] = v
    return m


def _reel(e: sp.Expr, valeurs: dict) -> mp.mpf:
    return mp.mpf(str(sp.N(e.xreplace(valeurs), 40)))


def fdt(
    schema: dict, w: float, sortie: int = 0, entree: str = "vd_t", point=None
) -> complex:
    """G(j w) de ``entree`` vers la sortie ``sortie`` (rang dans h) du schema,
    sur son modele linearise autour de ``point`` (dict de symboles, nul par
    defaut -- un bloc lineaire n'a pas de point privilegie).

    ``s x = Jfx x + Jfz z + Jfu u`` ; ``0 = Jgx x + Jgz z + Jgu u`` ;
    ``y = Jhx x + Jhz z + Jhu u``. Resolution en 40 chiffres (mpmath)."""
    c = B.compiler(schema)
    sym = c.symboles()
    X = [sym[n] for n in c.etats]
    Z = [sym[n] for n in c.algebriques]
    u = sym[entree]
    val = {sym[n]: sp.Float(v) for n, v in c.parametres.items()}
    for n in c.etats + c.algebriques + c.consignes + [entree]:
        val[sym[n]] = sp.Float(0.0)
    val.update(point or {})
    nx, nz = len(X), len(Z)

    def jac(exprs, vars_):
        return {
            (i, j): _reel(sp.diff(e, v), val)
            for i, e in enumerate(exprs)
            for j, v in enumerate(vars_)
        }

    Jfx, Jfz = jac(c.f, X), jac(c.f, Z)
    Jgx, Jgz = jac(c.g, X), jac(c.g, Z)
    Jfu = {(i, 0): _reel(sp.diff(e, u), val) for i, e in enumerate(c.f)}
    Jgu = {(i, 0): _reel(sp.diff(e, u), val) for i, e in enumerate(c.g)}
    h = c.h[sortie]
    jhx = [_reel(sp.diff(h, v), val) for v in X]
    jhz = [_reel(sp.diff(h, v), val) for v in Z]
    jhu = _reel(sp.diff(h, u), val)
    s = mp.mpc(0, w)
    n = nx + nz
    M = mp.matrix(n, n)
    rhs = mp.matrix(n, 1)
    for i in range(nx):
        for j in range(nx):
            M[i, j] = (s if i == j else 0) - Jfx[(i, j)]
        for j in range(nz):
            M[i, nx + j] = -Jfz[(i, j)]
        rhs[i] = Jfu[(i, 0)]
    for i in range(nz):
        for j in range(nx):
            M[nx + i, j] = Jgx[(i, j)]
        for j in range(nz):
            M[nx + i, nx + j] = Jgz[(i, j)]
        rhs[nx + i] = -Jgu[(i, 0)]
    sol = mp.lu_solve(M, rhs) if n else mp.matrix(0, 1)
    y = jhu + sum(jhx[j] * sol[j] for j in range(nx))
    y += sum(jhz[j] * sol[nx + j] for j in range(nz))
    return complex(y)


def valeurs_statiques(schema: dict) -> list[float]:
    """Les sorties ``h`` d'un schema SANS etat ni entree variable : on resout
    ``g = 0`` pour les algebriques puis on evalue ``h``."""
    c = B.compiler(schema)
    assert not c.etats, "schema statique attendu"
    inc = [sp.Symbol(n) for n in c.algebriques]
    if inc:
        sol = sp.solve(c.g, inc, dict=True)
        assert len(sol) == 1, "solution algebrique unique attendue"
        sub = sol[0]
    else:
        sub = {}
    # les parametres du modele prennent leur valeur (un gain K = "kk")
    valeurs = {sp.Symbol(k): v for k, v in c.parametres.items()}
    return [float(sp.N(e.xreplace(sub).xreplace(valeurs))) for e in c.h]


def erreur_relative(a: complex, b: complex) -> float:
    return abs(a - b) / abs(b)


# ---------------------------------------------------- empreinte d'un schema
def empreinte(c) -> dict:
    """Empreinte d'une compilation : noms (hash), tailles, et pour deux points
    pseudo-aleatoires fixes la somme ponderee et la somme des carres de chaque
    vecteur f, g, h. Deux compilations symboliquement egales ont la meme
    empreinte (a l'arrondi pres) ; l'inverse est vrai avec une probabilite
    pratiquement egale a 1."""
    noms = "|".join(
        list(c.etats)
        + ["#"]
        + list(c.algebriques)
        + ["#"]
        + list(c.consignes)
        + ["#"]
        + sorted(c.parametres)
    )
    libres = sorted(
        {s for e in (*c.f, *c.g, *c.h) for s in e.free_symbols}, key=lambda s: s.name
    )
    out = {
        "noms": hashlib.sha1(noms.encode()).hexdigest()[:12],
        "taille": (len(c.etats), len(c.algebriques), len(c.h)),
        "pts": [],
    }
    for graine in (12345, 6789):
        rng = np.random.default_rng(graine)
        pt = {s: float(rng.uniform(0.6, 1.4)) for s in libres}
        ligne = []
        for vec in (c.f, c.g, c.h):
            v = np.array([float(sp.N(e.xreplace(pt))) for e in vec])
            w = 1.0 + 0.37 * np.arange(v.size)
            ligne += [float(w @ v), float(v @ v)]
        out["pts"].append(ligne)
    return out


def memes_empreintes(a: dict, b: dict, tol: float = 1e-11) -> bool:
    if a["noms"] != b["noms"] or a["taille"] != b["taille"]:
        return False
    for la, lb in zip(a["pts"], b["pts"], strict=True):
        for x, y in zip(la, lb, strict=True):
            if abs(x - y) > tol * max(1.0, abs(y)):
                return False
    return True


# ------------------------------------------------------------ reseaux, spectres
def banc_smib(schema: dict, genre: str, X: float = 0.15) -> Tables:
    """Bus infini, ligne de reactance ``X``, element a blocs de 10 MVA au
    bus 2 (``genre`` : ``"gfl"`` PQ a 8 MW / 1 Mvar, ``"gfm"`` PV a 8 MW),
    comme ``tests/test_blocs.py``."""
    t = Tables(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3 if genre == "gfl" else 2, 20e3, 1.01, 0.0, 8e6, 1e6, 0.0, 0.0),
        ],
        T_Lines_pu=[(1, 1, 2, 1.0, 0.01, X, 0.0, 20e3, 100e6)],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
        T_Blocs=[(1, 2, "m", 10e6)],
        MODELES_BLOCS={"m": schema},
    )
    return t


def spectre(t: Tables):
    """(net, lf, valeurs propres de A) ; le load flow doit converger."""
    net = assembler(t)
    lf = net.run_loadflow()
    assert lf.converged, "le load flow ne converge pas"
    net.solve_equilibrium()
    ss = net.build_state_space()
    return net, lf, np.linalg.eigvals(ss.A)


def sans_jauge(lam: np.ndarray, seuil: float = 1e-6) -> np.ndarray:
    """Retire le(s) mode(s) de jauge (valeur propre nulle)."""
    return lam[np.abs(lam) > seuil]


def ecart_spectres(a: np.ndarray, b: np.ndarray) -> float:
    """Ecart relatif maximal entre deux spectres apres appariement optimal
    (probleme d'affectation) : trier par partie reelle serait fragile quand
    deux modes ont des parties reelles voisines. Relatif a ``max(1, |a|)`` :
    les modes de module inferieur a 1 rad/s sont juges en absolu (le bruit
    numerique de ``eig`` est de l'ordre de 1e-11 sur ces matrices)."""
    assert a.size == b.size, f"{a.size} modes contre {b.size}"
    cout = np.abs(a[:, None] - b[None, :])
    i, j = linear_sum_assignment(cout)
    return float(np.max(cout[i, j] / np.maximum(1.0, np.abs(a[i]))))


def copie(schema: dict) -> dict:
    return copy.deepcopy(schema)


# --------------------------------------------------- structure d'un schema
def liens_canoniques(schema: dict) -> set[tuple]:
    return {
        (
            lien["de"],
            int(lien.get("sortie", 0)),
            lien["vers"],
            int(lien.get("entree", 0)),
        )
        for lien in schema["liens"]
    }


def sans_position(bloc: dict) -> dict:
    return {k: v for k, v in bloc.items() if k not in ("x", "y")}


def schema_sans_position(schema: dict) -> dict:
    """Forme comparable d'un schema (niveau courant, recursif) : blocs par id
    sans positions, liens canoniques, sous-schemas traites de meme."""
    blocs = {}
    for b in schema["blocs"]:
        b = sans_position(b)
        if b.get("type") == "sous_systeme":
            b["schema"] = schema_sans_position(b["schema"])
        blocs[b["id"]] = b
    return {
        "blocs": blocs,
        "liens": liens_canoniques(schema),
        "parametres": schema.get("parametres"),
    }


def tous_les_niveaux(schema: dict):
    """Parcourt le schema : rend (chemin, liste de blocs) pour chaque niveau."""
    yield "", schema["blocs"]
    for b in schema["blocs"]:
        if b.get("type") == "sous_systeme":
            for chemin, blocs in tous_les_niveaux(b["schema"]):
                yield f"{b['id']}/{chemin}", blocs


BLOC_L, BLOC_H = 92, 50


def recouvrements(blocs: list[dict]) -> list[tuple[str, str]]:
    """Les paires de blocs dont les boites ``BLOC_L x BLOC_H`` (coin haut
    gauche en ``x, y``) ont une intersection d'aire non nulle."""
    out = []
    for i, a in enumerate(blocs):
        for b in blocs[i + 1 :]:
            # chevauchement = dimension de la boite - ecart des coins
            ox = BLOC_L - abs(a["x"] - b["x"])
            oy = BLOC_H - abs(a["y"] - b["y"])
            if ox > 1e-9 and oy > 1e-9:
                out.append((a["id"], b["id"]))
    return out
