"""Integrateur DAE MONOLITHIQUE : etats, variables algebriques et cablage, ensemble.

POURQUOI -- :meth:`ModeleNonLineaire.simuler` confie les etats a ``Radau`` et
resout l'algebrique DANS chaque appel du second membre : un Newton complet par
evaluation, et Radau en demande des milliers. Surtout, Radau est un
integrateur a PAS ADAPTATIF d'ordre 5 : il suit a la tolerance pres TOUS les
modes du modele, y compris les resonances electromagnetiques de quelques
dizaines de kHz qu'un defaut franc excite (inductance de defaut contre
capacite de bus) -- mesure : 17 800 pas pour 50 ms de defaut sur ``smib_gfm``,
plusieurs minutes de calcul.

Les outils EMT font autrement, et VeraGrid aussi : un pas FIXE, un schema
implicite qui AMORTIT ce qu'il ne resout pas (trapezes pour VeraGrid-EMT,
Euler implicite pour VeraGrid-RMS), et un seul systeme non lineaire par pas,
etats et contraintes ensemble. C'est ce que fait ce module, en BDF2 :

.. math::

    \\frac{\\alpha_0 x_{n+1} + \\alpha_1 x_n + \\alpha_2 x_{n-1}}{h} = f(x_{n+1}, z_{n+1}, u_{n+1})
    \\qquad 0 = g(\\ldots) \\qquad u = R_{uext} u_{ext} + R_y h(\\ldots)

resolu par un Newton sur ``y = [x, z, u]`` dont le jacobien ::

    [ a.I - fx     -fz      -fu     ]
    [    gx         gz       gu     ]      a = alpha_0 / h
    [ -Ry.hx     -Ry.hz   I - Ry.hu ]

est creux, a STRUCTURE FIXE (calculee une fois, cf. :func:`structure_dae`), et
factorise en LU gardee d'un pas a l'autre (Newton modifie) : elle ne se refait
que si le pas change ou si la convergence ralentit.

BDF2 et non les trapezes : les deux sont d'ordre 2 et A-stables, mais BDF2 est
L-stable -- une resonance bien plus rapide que le pas est AMORTIE, au lieu de
rebondir d'un pas a l'autre (le "ringing" des trapezes). Ce qu'il en coute :
une erreur d'amortissement sur les modes proches de la frequence du pas, qui
se regle en choisissant le pas (``1/(400.f0)`` par defaut, ~40 us a 60 Hz :
les modes electromecaniques et de controle jusqu'au kHz sont justes a mieux
que le pourcent).

:func:`equilibre_pseudo_transitoire` reprend le meme pas implicite pour
trouver un point d'equilibre quand l'initialisation explicite n'en donne pas
un (VeraGrid, "pseudo-transient continuation") : un Euler implicite a pas
croissant sur ``dx/dtau = f``, puis un Newton de finition, et le residu
verifie avant d'accepter.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

__all__ = [
    "equilibre_pseudo_transitoire",
    "integrer_bdf2",
    "jacobien_dae",
    "reinitialiser",
    "residu_dae",
    "structure_dae",
]


# ------------------------------------------------------------------ structure
def structure_dae(m) -> dict:
    """La structure du jacobien DAE complet de ``m``, calculee une fois.

    Meme principe que ``ModeleNonLineaire._structure`` (le jacobien de la seule
    partie algebrique) : on numerote les coefficients des blocs d'elements dans
    l'ordre ou :func:`jacobien_dae` les concatene, chacun va a une ou plusieurs
    places (via ``-Ry`` pour les blocs ``h``), et ``M`` porte ces poids. S'y
    ajoutent l'identite du bloc ``u`` (constante) et la diagonale du bloc ``x``,
    ou se pose ``a = alpha_0/h`` -- qui change avec le pas, sans que la
    structure change.
    """
    st = getattr(m, "_st_dae", None)
    if st is not None:
        return st
    c = m.c
    nx, nz, nu = c.nx, c.nz, c.nu
    N = nx + nz + nu
    oz, ou = nx, nx + nz
    Ry = sp.csc_matrix(c.Ry)
    lignes, colonnes, sources, poids = [], [], [], []
    ordre: list[tuple[int, str]] = []
    n = 0

    def plein(r0, c0, nr, nc, w):
        nonlocal n
        if nr * nc == 0:
            return
        i, j = np.divmod(np.arange(nr * nc), nc)
        lignes.append(r0 + i)
        colonnes.append(c0 + j)
        sources.append(n + np.arange(nr * nc))
        poids.append(np.full(nr * nc, float(w)))
        n += nr * nc

    def cable(y0, c0, nr, nc):
        nonlocal n
        if nr * nc == 0:
            return
        for i in range(nr):
            col = Ry[:, y0 + i]
            for pp, wgt in zip(col.indices, col.data):
                lignes.append(np.full(nc, ou + pp))
                colonnes.append(c0 + np.arange(nc))
                sources.append(n + i * nc + np.arange(nc))
                poids.append(np.full(nc, -wgt))
        n += nr * nc

    position: dict[tuple[int, str], tuple[int, int]] = {}
    for k, cm in enumerate(c.cms):
        sx, sz, su, sy = c.ox[k], c.oz[k], c.ou[k], c.oy[k]
        blocs = [
            ("fx", sx, sx, cm.n_x, cm.n_x, -1.0, plein),
            ("fz", sx, oz + sz, cm.n_x, cm.n_z, -1.0, plein),
            ("fu", sx, ou + su, cm.n_x, cm.n_u, -1.0, plein),
            ("gx", oz + sz, sx, cm.n_z, cm.n_x, 1.0, plein),
            ("gz", oz + sz, oz + sz, cm.n_z, cm.n_z, 1.0, plein),
            ("gu", oz + sz, ou + su, cm.n_z, cm.n_u, 1.0, plein),
            ("hx", sy, sx, cm.n_y, cm.n_x, None, cable),
            ("hz", sy, oz + sz, cm.n_y, cm.n_z, None, cable),
            ("hu", sy, ou + su, cm.n_y, cm.n_u, None, cable),
        ]
        for nom, r0, c0, nr, nc, w, pose in blocs:
            if nr * nc == 0:
                continue
            position[(k, nom)] = (n, nr * nc)
            if pose is plein:
                plein(r0, c0, nr, nc, w)
            else:
                cable(r0, c0, nr, nc)
            ordre.append((k, nom))

    cat = lambda L, t: np.concatenate(L).astype(t) if L else np.zeros(0, t)  # noqa: E731
    li, co = cat(lignes, np.int64), cat(colonnes, np.int64)
    so, po = cat(sources, np.int64), cat(poids, float)
    diag_u = ou + np.arange(nu)
    diag_x = np.arange(nx)
    cles = np.concatenate([li * N + co, diag_u * N + diag_u, diag_x * N + diag_x])
    places, inv = np.unique(cles, return_inverse=True)
    n_li = len(li)
    M = sp.csr_matrix((po, (inv[:n_li], so)), shape=(len(places), n))
    const = np.zeros(len(places))
    np.add.at(const, inv[n_li : n_li + nu], 1.0)
    diag = np.zeros(len(places))
    diag[inv[n_li + nu :]] = 1.0
    r, cc = np.divmod(places, N)
    perm = np.lexsort((r, cc))
    indptr = np.concatenate([[0], np.cumsum(np.bincount(cc, minlength=N))])
    st = {
        "M": M,
        "const": const[perm],
        "diag": diag[perm],
        "perm": perm,
        "indices": r[perm].astype(np.int32),
        "indptr": indptr.astype(np.int32),
        "N": N,
        "ordre": ordre,
        "n_valeurs": n,
        # par lot et par bloc : ou vont ses valeurs dans le vecteur ``v``
        "lots": [
            (
                i_lot,
                nom,
                np.array([position[(k, nom)][0] for k in lot.ks])[None, :]
                + np.arange(position[(lot.ks[0], nom)][1])[:, None],
            )
            for i_lot, lot in enumerate(lots(m))
            for nom in ("fx", "fz", "fu", "gx", "gz", "gu", "hx", "hz", "hu")
            if (lot.ks[0], nom) in position
        ],
    }
    m._st_dae = st
    return st


def lots(m) -> list[SimpleNamespace]:
    """Les elements de ``m`` groupes par MODELE COMPILE, avec les indices qui
    rangent leurs variables en colonnes : ``X = x[lot.ix]`` est ``(n_x,
    n_instances)``, et ``f[lot.ix] = ...`` y reverse le resultat."""
    L = getattr(m, "_lots_dae", None)
    if L is not None:
        return L
    c = m.c
    groupes: dict[int, list[int]] = {}
    for k, cm in enumerate(c.cms):
        groupes.setdefault(id(cm), []).append(k)
    L = []
    for ks in groupes.values():
        cm = c.cms[ks[0]]
        ks = np.array(ks)

        def rangs(offsets, n, _ks=ks):
            return np.array(offsets)[_ks][None, :] + np.arange(n)[:, None]

        L.append(
            SimpleNamespace(
                cm=cm,
                ks=ks,
                ix=rangs(c.ox, cm.n_x),
                iz=rangs(c.oz, cm.n_z),
                iu=rangs(c.ou, cm.n_u),
                iy=rangs(c.oy, cm.n_y),
            )
        )
    m._lots_dae = L
    return L


def _parametres(m, lot) -> np.ndarray:
    """Les parametres du lot en colonnes, empiles une fois par MODELE : un
    clone de mode (``avec_mode``) a sa propre liste ``p``, donc son propre
    cache."""
    cache = m.__dict__.get("_P_lots")
    if cache is None or cache[0] is not m.p:
        cache = m._P_lots = (m.p, {})
    P = cache[1].get(id(lot))
    if P is None:
        P = cache[1][id(lot)] = np.column_stack([m.p[k] for k in lot.ks])
    return P


def jacobien_dae(m, x, z, u, a: float):
    """Le jacobien de :func:`residu_dae` en ``y = [x, z, u]``, creux (CSC)."""
    st = structure_dae(m)
    v = np.empty(st["n_valeurs"])
    L = lots(m)
    cache = {}
    for i_lot, nom, pos in st["lots"]:
        lot = L[i_lot]
        if i_lot not in cache:
            cache[i_lot] = (x[lot.ix], z[lot.iz], u[lot.iu], _parametres(m, lot))
        v[pos] = lot.cm.lot(nom)(*cache[i_lot])
    data = (st["M"] @ v)[st["perm"]] + st["const"] + a * st["diag"]
    N = st["N"]
    J = sp.csc_matrix((N, N))
    J.data, J.indices, J.indptr = data, st["indices"], st["indptr"]
    return J


def residu_dae(m, x, z, u, uext, xpoint) -> np.ndarray:
    """``[xpoint - f ; g ; u - Ruext.uext - Ry.h]`` au point ``(x, z, u)``."""
    c = m.c
    f, g, y = np.empty(c.nx), np.empty(c.nz), np.empty(c.ny)
    for lot in lots(m):
        cm = lot.cm
        X, Z, U, P = x[lot.ix], z[lot.iz], u[lot.iu], _parametres(m, lot)
        if cm.n_x:
            f[lot.ix] = cm.lot("f")(X, Z, U, P)
        if cm.n_z:
            g[lot.iz] = cm.lot("g")(X, Z, U, P)
        y[lot.iy] = cm.lot("h")(X, Z, U, P)
    Ruext = getattr(m, "_Ruext", None)
    if Ruext is None:
        Ruext = m._Ruext = sp.csr_matrix(c.Ruext)
    return np.concatenate([xpoint - f, g, u - (Ruext @ uext + m._Ry @ y)])


# ------------------------------------------------------------------ un pas
class _Newton:
    """Le Newton modifie d'un pas implicite, et sa factorisation gardee."""

    def __init__(self, m):
        self.m = m
        self.lu = None
        self.a = None
        self.n_factorisations = 0

    def resoudre(self, x, z, u, uext, a, hist, tol=1e-8):
        """Resout ``a.x + hist = f(x, z, u)`` et les contraintes ; rend
        ``(x, z, u)`` ou ``None`` si Newton n'y arrive pas."""
        m = self.m
        nx, nz = len(x), len(z)
        y0 = np.concatenate([x, z, u])
        y = y0.copy()
        for essai in range(2):
            if self.lu is None or self.a is None or abs(a - self.a) > 1e-12 * abs(a):
                J = jacobien_dae(m, y[:nx], y[nx : nx + nz], y[nx + nz :], a)
                try:
                    self.lu = spla.splu(J)
                except RuntimeError:
                    self.lu = None
                    return None
                self.a = a
                self.n_factorisations += 1
            pas_prec = np.inf
            for it in range(10):
                xk, zk, uk = y[:nx], y[nx : nx + nz], y[nx + nz :]
                F = residu_dae(m, xk, zk, uk, uext, a * xk + hist)
                dy = self.lu.solve(F)
                if not np.all(np.isfinite(dy)):
                    break
                y = y - dy
                n_dy = np.abs(dy).max()
                if n_dy <= tol * (1.0 + np.abs(y).max()):
                    # La corde a converge, mais lentement : la LU gardee s'est
                    # eloignee du jacobien (l'etat a derive depuis). On la
                    # refera au pas suivant, avant que chaque pas ne coute dix
                    # residus au lieu de deux.
                    if it >= 3:
                        self.lu = None
                    return y[:nx], y[nx : nx + nz], y[nx + nz :]
                if n_dy > 0.8 * pas_prec:  # la corde ne tient plus
                    break
                pas_prec = n_dy
            # jacobien au point courant, puis on reessaie une fois
            self.lu = None
            if essai == 0:
                continue
        # Dernier recours, un Newton PLEIN depuis l'amorce : la corde part du
        # jacobien de l'amorce, et s'egare quand la solution du pas en est
        # loin -- un retro-calcul tres rapide qui s'enclenche avec la
        # limitation (des milliers sur ``y_vsat_aw`` en un pas), le premier
        # pas apres une reinitialisation. Un jacobien a chaque iteration y
        # converge en trois ou quatre.
        y = y0.copy()
        for _it in range(25):
            xk, zk, uk = y[:nx], y[nx : nx + nz], y[nx + nz :]
            F = residu_dae(m, xk, zk, uk, uext, a * xk + hist)
            try:
                dy = spla.splu(jacobien_dae(m, xk, zk, uk, a)).solve(F)
            except RuntimeError:
                return None
            self.n_factorisations += 1
            if not np.all(np.isfinite(dy)):
                return None
            y = y - dy
            if np.abs(dy).max() <= tol * (1.0 + np.abs(y).max()):
                return y[:nx], y[nx : nx + nz], y[nx + nz :]
        return None


# ------------------------------------------------------------ reinitialiser
#: le pas des deux Euler implicites de :func:`reinitialiser`, en secondes
H_REINIT = 1e-8


def reinitialiser(m, x, z, u, uext, h: float = H_REINIT, iter_max: int = 20):
    """Un point ``(x, z, u)`` coherent pour repartir, a ``x`` presque fige.

    A noeuds ideaux, l'algebrique a ``x`` fige n'a pas de solution unique
    (indice 2) : la tension d'un noeud ou ne se rencontrent que des courants
    d'etat est celle qui annule ``d(isum)/dt``, pas une solution de
    ``isum = 0``. Deux Euler implicites, de pas ``+h`` puis ``-h``, resolvent
    ``x``, ``z`` et ``u`` ensemble :

    - le premier remet ``x`` sur la contrainte. Si une manoeuvre a retire une
      branche qui portait du courant (elimination d'un defaut, ouverture),
      les courants restants violent la loi des noeuds ; la tension du pas,
      ``~ Phi/h``, est l'impulsion de la limite ideale, et le saut de ``x``
      qu'elle produit conserve le flux, quel que soit ``h`` ;
    - le second, de ``x`` desormais coherent (sans impulsion), defait le
      ``h`` de dynamique reguliere du premier a ``O(h^2)`` pres, et rend la
      tension qui annule ``d(isum)/dt`` a l'instant de depart.

    Sans lui, un intervalle repart de la tension d'AVANT l'evenement : les
    gardes des limiteurs la voient au depart, et le premier pas stocke une
    tension ``Phi/pas`` qu'aucun reseau ne porte. Le pas retrograde compte :
    deux pas en avant decalent l'etat de ``2h`` a chaque reprise, et les
    basculements a repetition d'un limiteur verrouille amplifient ce
    decalage jusqu'a changer la sequence -- un schema-bloc, dont le
    retro-calcul a un mode vers 1e6 1/s, y est sensible des 1e-9 s. La
    tolerance porte sur ``x`` : la tension impulsive du premier pas ferait
    une tolerance relative a tout ``y`` bien trop lache sur les etats.
    """
    nx, nz = len(x), len(z)
    y = np.concatenate([x, z, u]).astype(float)
    for a in (1.0 / h, -1.0 / h):
        x_prec = y[:nx].copy()
        for _it in range(iter_max):
            xk, zk, uk = y[:nx], y[nx : nx + nz], y[nx + nz :]
            F = residu_dae(m, xk, zk, uk, uext, a * (xk - x_prec))
            dy = spla.splu(jacobien_dae(m, xk, zk, uk, a)).solve(F)
            if not np.all(np.isfinite(dy)):
                raise RuntimeError("reinitialisation : pas de Newton non fini")
            y = y - dy
            if np.abs(dy[:nx]).max() <= 1e-12 * (1.0 + np.abs(y[:nx]).max()) and (
                np.abs(dy[nx:]).max() <= 1e-9 * (1.0 + np.abs(y[nx:]).max())
            ):
                break
        else:
            raise RuntimeError(
                f"reinitialisation : Newton non converge en {iter_max} iterations"
            )
    return y[:nx], y[nx : nx + nz], y[nx + nz :]


# ------------------------------------------------------------------ integrer
#: le nombre de pas entre deux appels du ``rappel`` d'avancement
PAS_RAPPEL = 200


def integrer_bdf2(
    m,
    t_span: tuple[float, float],
    x0: np.ndarray,
    z0: np.ndarray,
    u0: np.ndarray,
    consigne,
    pas: float,
    t_rendus: np.ndarray | None = None,
    gardes=None,
    rappel=None,
    garder_algebriques: bool = False,
) -> SimpleNamespace:
    """Integre ``m`` en BDF2 a pas ``pas`` sur ``t_span``.

    ``consigne(t)`` rend ``uext`` ; ``gardes`` est ``(valeurs, sens)`` : une
    fonction ``(x, z, u) -> liste de g`` et la liste des sens de franchissement
    a surveiller. La premiere garde franchie ARRETE l'integration, a l'instant
    interpole lineairement entre les deux pas qui l'encadrent (a ``pas`` pres,
    le pas etant petit devant la dynamique qu'un limiteur regarde).

    Rend un objet a la maniere de ``solve_ivp`` : ``t``, ``y`` (etats, une
    colonne par instant rendu), ``status`` (0 fin, 1 garde franchie),
    ``t_events``, ``y_events``, ``success``, ``message`` -- plus ``z`` et ``u``,
    la solution algebrique au dernier pas, dont la suite peut repartir.
    ``rappel(t)``, s'il est donne, est appele tous les ``PAS_RAPPEL`` pas :
    l'avancement d'une longue simulation, et son point d'interruption.

    ``garder_algebriques`` : rendre aussi ``zu``, les ``(z, u)`` resolus aux
    instants rendus (``(n_t, nz + nu)``, interpoles comme les etats). Un
    modele d'indice 2 (noeuds ideaux) n'a pas de solution algebrique a ``x``
    fige : c'est la seule facon de relire ses tensions de noeud.
    """
    t0, t1 = float(t_span[0]), float(t_span[1])
    newton = _Newton(m)
    x_n = np.asarray(x0, dtype=float).copy()
    z, u = np.asarray(z0, dtype=float).copy(), np.asarray(u0, dtype=float).copy()
    x_prec = None
    h_prec = None
    ts, xs_ = [t0], [x_n.copy()]
    zus = [np.concatenate([z, u])] if garder_algebriques else None
    g_prec = None
    n_gardes = 0
    if gardes is not None:
        valeurs, sens = gardes
        g_prec = np.asarray(valeurs(x_n, z, u), dtype=float)
        n_gardes = len(sens)
    t = t0
    status, t_ev, y_ev = 0, None, None
    i_ev = None
    fin = t1 - 1e-12 * max(1.0, abs(t1))
    while t < fin:
        h = min(pas, t1 - t)
        if x_prec is None:
            a0, a1, a2 = 1.0, -1.0, 0.0  # Euler implicite au premier pas
            predit = x_n.copy()
        else:
            w = h / h_prec
            a0, a1, a2 = (1 + 2 * w) / (1 + w), -(1 + w), w * w / (1 + w)
            predit = x_n + (x_n - x_prec) * w
        hist = (a1 * x_n + a2 * (x_prec if x_prec is not None else 0.0)) / h
        sol = newton.resoudre(predit, z, u, consigne(t + h), a0 / h, hist)
        if sol is None:
            raise RuntimeError(
                f"BDF2 : Newton en echec a t = {t + h:g} s (pas {h:g} s) ; "
                "reduire le pas, ou le probleme n'a plus de solution"
            )
        x, z, u = sol
        if not np.all(np.isfinite(x)):
            raise RuntimeError(f"BDF2 : etat non fini a t = {t + h:g} s")
        if gardes is not None:
            g = np.asarray(valeurs(x, z, u), dtype=float)
            franchies = [
                i
                for i in range(n_gardes)
                if (sens[i] > 0 and g_prec[i] < 0 <= g[i])
                or (sens[i] < 0 and g_prec[i] > 0 >= g[i])
            ]
            if franchies:
                # l'instant de franchissement, par interpolation lineaire
                frac = {
                    i: g_prec[i] / (g_prec[i] - g[i]) if g_prec[i] != g[i] else 1.0
                    for i in franchies
                }
                i_ev = min(franchies, key=lambda i: frac[i])
                th = frac[i_ev]
                t_ev = t + th * h
                y_ev = x_n + th * (x - x_n)
                ts.append(t_ev)
                xs_.append(y_ev.copy())
                if zus is not None:
                    zus.append(zus[-1] + th * (np.concatenate([z, u]) - zus[-1]))
                status = 1
                break
            g_prec = g
        x_prec, x_n, h_prec = x_n, x, h
        t = t + h
        ts.append(t)
        xs_.append(x.copy())
        if zus is not None:
            zus.append(np.concatenate([z, u]))
        if rappel is not None and len(ts) % PAS_RAPPEL == 0:
            rappel(t)

    T = np.array(ts)
    X = np.array(xs_).T  # (n_x, n_t)
    ZU = np.array(zus).T if zus is not None else None
    if t_rendus is not None:
        tr = np.asarray(t_rendus, dtype=float)
        tr = tr[tr <= T[-1] + 1e-12 * max(1.0, abs(T[-1]))]
        if status == 1:
            tr = tr[tr < t_ev]
        X = np.array([np.interp(tr, T, X[i]) for i in range(X.shape[0])]).reshape(
            X.shape[0], len(tr)
        )
        if ZU is not None:
            ZU = np.array(
                [np.interp(tr, T, ZU[i]) for i in range(ZU.shape[0])]
            ).reshape(ZU.shape[0], len(tr))
        T = tr
    t_events = [[] for _ in range(n_gardes)]
    y_events = [[] for _ in range(n_gardes)]
    if status == 1:
        t_events[i_ev] = [t_ev]
        y_events[i_ev] = [y_ev]
    return SimpleNamespace(
        t=T,
        y=X,
        status=status,
        t_events=[np.array(e) for e in t_events],
        y_events=[np.array(e) for e in y_events],
        success=True,
        message="BDF2 : fin de l'intervalle"
        if status == 0
        else "BDF2 : garde franchie",
        z=z,
        u=u,
        zu=None if ZU is None else ZU.T,
        n_factorisations=newton.n_factorisations,
    )


# ------------------------------------------------------------------ equilibre
def equilibre_pseudo_transitoire(
    m,
    uext=None,
    x0=None,
    tol: float = 1e-9,
    dtau0: float = 1e-3,
    iter_max: int = 400,
):
    """Un point d'equilibre de ``m`` par continuation pseudo-transitoire.

    Euler implicite sur ``dx/dtau = f(x, z, u)``, contraintes tenues a chaque
    pas, avec un pas ``dtau`` qui croit comme le residu decroit (regle SER :
    ``dtau <- dtau . |f|_prec / |f|``, borne a x10 par pas) ; le pas fini
    regularise le jacobien la ou un Newton direct buterait sur la direction de
    jauge (angle absolu libre). Rend ``(x, z, u, residu)`` ; leve si le residu
    ``|f|`` ne descend pas sous ``tol``.
    """
    uext = m.uext0() if uext is None else uext
    x = m.c.x0() if x0 is None else np.asarray(x0, dtype=float)
    # pas la solution algebrique a x fige, qu'un noeud ideal n'a pas toujours
    x, z, u = m.amorcer(x, uext, m.c.z0(), m.c.u0())
    newton = _Newton(m)

    def residu_f(x, z, u):
        return np.abs(residu_dae(m, x, z, u, uext, np.zeros(len(x)))[: len(x)]).max()

    r = residu_f(x, z, u)
    dtau = dtau0
    for _k in range(iter_max):
        if r < tol:
            return x, z, u, r
        sol = newton.resoudre(x, z, u, uext, 1.0 / dtau, -x / dtau)
        if sol is None:
            dtau /= 4.0
            newton.lu = None
            continue
        x2, z2, u2 = sol
        r2 = residu_f(x2, z2, u2)
        x, z, u = x2, z2, u2
        dtau = min(dtau * min(10.0, max(1.0, r / max(r2, 1e-300))), 1e8)
        r = r2
    if r < tol:
        return x, z, u, r
    raise RuntimeError(
        f"equilibre pseudo-transitoire : residu {r:.2e} apres {iter_max} pas "
        "(le point cherche n'existe peut-etre pas)"
    )
