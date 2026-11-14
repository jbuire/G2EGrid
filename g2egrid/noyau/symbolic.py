"""Couche symbolique-numerique hybride (approche ANDES).

Un modele d'element est declare une seule fois, symboliquement, sous forme DAE :

    dx/dt = f(x, z, u, p)      equations differentielles
    0     = g(x, z, u, p)      equations algebriques
    y     = h(x, z, u, p)      equations de sortie

SymPy derive les jacobiens analytiquement, puis `lambdify` produit des
fonctions numpy rapides. La compilation est faite **une fois par CLASSE**
(les parametres restent symboliques), pas par instance : tous les SM du
reseau partagent le meme code compile.

C'est la correction du defaut du framework MATLAB, ou `buildSymbolicModel`
+ `jacobian` etaient re-executes a chaque appel de `getLinearModel`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from types import SimpleNamespace

import numpy as np
import sympy as sp

__all__ = ["CompiledModel", "ModelSpec", "clear_cache", "compile_spec"]


# --------------------------------------------------------------------------
# Declaration
# --------------------------------------------------------------------------
@dataclass
class ModelSpec:
    """Declaration symbolique d'un modele DAE d'element.

    Les listes contiennent des NOMS ; les symboles sont crees par
    :func:`compile_spec`, qui les passe a ``equations`` dans un namespace.
    """

    name: str
    params: tuple[str, ...] = ()
    states: tuple[str, ...] = ()
    algebs: tuple[str, ...] = ()
    inputs: tuple[str, ...] = ()
    outputs: tuple[str, ...] = ()
    # equations(s) -> (f, g, h) : listes d'expressions sympy, memes longueurs
    # que states / algebs / outputs respectivement.
    equations: callable = None

    def __post_init__(self) -> None:
        for grp in ("params", "states", "algebs", "inputs", "outputs"):
            setattr(self, grp, tuple(getattr(self, grp)))
        seen: dict[str, str] = {}
        for grp in ("params", "states", "algebs", "inputs", "outputs"):
            for nm in getattr(self, grp):
                if nm in seen:
                    raise ValueError(
                        f"{self.name}: nom '{nm}' declare deux fois "
                        f"({seen[nm]} et {grp})"
                    )
                seen[nm] = grp


# --------------------------------------------------------------------------
# Compilation
# --------------------------------------------------------------------------
@dataclass
class CompiledModel:
    """Modele DAE compile : residus + jacobiens, en fonctions numpy."""

    spec: ModelSpec
    n_x: int
    n_z: int
    n_u: int
    n_y: int
    n_p: int
    # residus : (x, z, u, p) -> ndarray
    f: callable
    g: callable
    h: callable
    # jacobiens : (x, z, u, p) -> ndarray 2D
    fx: callable
    fz: callable
    fu: callable
    gx: callable
    gz: callable
    gu: callable
    hx: callable
    hz: callable
    hu: callable
    #: les expressions symboliques (arguments, f, g, h et leurs jacobiens),
    #: gardees pour compiler a la demande les versions EN LOT (:meth:`lot`)
    expressions: dict | None = field(default=None, repr=False)
    _lots: dict = field(default_factory=dict, repr=False)
    _dparams: dict = field(default_factory=dict, repr=False)

    # -- derivees par rapport aux parametres -------------------------------
    def derivees_parametres(self, figes: dict[str, float] | None = None):
        """``d J / d p`` du jacobien local complet, sous forme creuse.

        ``J`` est le jacobien EMPILE de l'element ::

                    x    z    u
            J = [ fx   fz   fu ]   f
                [ gx   gz   gu ]   g
                [ hx   hz   hu ]   h

        Rend ``(ip, lignes, colonnes, valeurs)`` : pour chaque coefficient non
        nul de ``dJ/dp``, l'indice du parametre, sa position dans ``J``, et une
        fonction ``valeurs(x, z, u, p)`` qui les evalue TOUS d'un appel. Les
        derivees sont SYMBOLIQUES -- les equations sont deja la, il suffit de
        les deriver une fois de plus -- et compilees a la premiere demande,
        pour la classe entiere : c'est ce qui rend la sensibilite d'un mode a
        chaque parametre (:mod:`g2egrid.analyse.sensibilites`) exacte et bon marche.

        ``figes`` : ``{nom: valeur}`` des parametres DRAPEAUX des modes
        discrets (:attr:`g2egrid.noyau.element.Element.modes`), a leur valeur au
        point de fonctionnement. Les derivees sont prises a ces valeurs, et
        compilees une fois par classe ET par jeu de valeurs.

        POURQUOI FIGER LES DRAPEAUX -- un limiteur ecrit ``(1 - m) + m.Min(1,
        a)``, et un GFM en empile deux (courant, puis tension sur la
        modulation qui en depend) autour de l'impedance virtuelle, un
        ``Piecewise`` : un seul coefficient de ``fx`` y compte 15 000
        operations SymPy, sa derivee 120 000, et ``dJ/dp`` d'un GFM ne rendait
        plus (plus de dix minutes). Les drapeaux, eux, ne varient pas autour
        d'un point de fonctionnement : remplaces par leur valeur (0, en mode
        libre), les ``Min`` disparaissent et il reste les equations du
        convertisseur libre. La derivee par rapport a un drapeau se prend les
        AUTRES figes, puis a sa valeur : elle reste exacte. Un drapeau a 1
        (un point de fonctionnement pose a la main en mode limite : le modele
        linearise est d'ordinaire en mode libre, et la simulation garde ses
        drapeaux pour elle) ramene son ``Min`` : exact encore, mais de l'ordre
        de la minute pour un GFM.

        POURQUOI ``d/dv (dR/dp)`` -- les derivees croisees sont egales, et le
        residu ``R = (f, g, h)`` est bien plus court que son jacobien : on le
        derive par rapport a ``p``, puis par rapport aux variables
        ``v = (x, z, u)``. Le resultat est le meme, au dernier chiffre pres.
        """
        figes = figes or {}
        cle = tuple(sorted((nom, float(v)) for nom, v in figes.items()))
        if cle in self._dparams:
            return self._dparams[cle]
        ex = self.expressions
        if ex is None:
            raise RuntimeError(f"{self.spec.name}: expressions non conservees")
        xs, zs, us, ps = ex["args"]
        vs = [*xs, *zs, *us]  # colonnes de J, dans l'ordre x, z, u
        residus = [*ex["f"], *ex["g"], *ex["h"]]  # lignes de J, f, g, h
        valeur = {ps[self.spec.params.index(nom)]: _nombre(v) for nom, v in cle}
        fige = [r.xreplace(valeur) for r in residus]
        libres = [r.free_symbols for r in fige]
        ip, li, co, exprs = [], [], [], []
        for k, sym in enumerate(ps):
            if sym in valeur:  # un drapeau : les autres figes, lui a sa valeur
                autres = {s: v for s, v in valeur.items() if s != sym}
                drs = [
                    sp.diff(r.xreplace(autres), sym).xreplace({sym: valeur[sym]})
                    if sym in r.free_symbols
                    else sp.S.Zero
                    for r in residus
                ]
            else:
                drs = [
                    sp.diff(r, sym) if sym in lib else sp.S.Zero
                    for r, lib in zip(fige, libres)
                ]
            for i, dr in enumerate(drs):
                if dr == 0:
                    continue
                lib = dr.free_symbols
                for j, v in enumerate(vs):
                    if v not in lib:
                        continue
                    # la derivee d'un seuil (Heaviside des limiteurs) est un
                    # Dirac : nulle partout sauf SUR le seuil, ou la
                    # linearisation n'a de toute facon pas de sens
                    d = sp.diff(dr, v).replace(
                        lambda a: isinstance(a, sp.DiracDelta),
                        lambda a: sp.S.Zero,
                    )
                    if d == 0:
                        continue
                    ip.append(k)
                    li.append(i)
                    co.append(j)
                    exprs.append(d)
        n = len(exprs)
        if n:
            brut = sp.lambdify(ex["args"], exprs, "numpy", cse=True)

            def valeurs(x, z, u, p, _b=brut, _n=n):
                return np.asarray(_b(x, z, u, p), dtype=float).reshape(_n)

        else:

            def valeurs(x, z, u, p):
                return np.zeros(0)

        self._dparams[cle] = out = (
            np.array(ip, dtype=int),
            np.array(li, dtype=int),
            np.array(co, dtype=int),
            valeurs,
        )
        return out

    # -- evaluation en lot --------------------------------------------------
    def lot(self, nom: str):
        """La fonction ``nom`` (``f``, ``g``, ``h``, ``fx``...) EN LOT.

        Elle prend ``X, Z, U, P`` de forme ``(n, n_instances)`` -- une colonne
        par element de la meme classe -- et rend ``(n_valeurs, n_instances)``,
        les blocs de jacobien aplatis ligne a ligne comme ``np.ravel``.

        POURQUOI -- un pas d'integration evalue ``f``, ``g`` et ``h`` de chaque
        element ; appeles un par un, c'est une cinquantaine d'appels Python
        par residu sur wscc9, et le cout d'appel l'emporte sur le calcul. Les
        memes expressions, lambdifiees avec des tableaux pour arguments,
        s'evaluent pour toutes les lignes (ou toutes les machines) d'un seul
        appel NumPy. Compilees a la premiere demande seulement : la voie
        lineaire n'en a pas besoin.
        """
        fn = self._lots.get(nom)
        if fn is not None:
            return fn
        exprs = self.expressions
        if exprs is None:
            raise RuntimeError(f"{self.spec.name}: expressions non conservees")
        elements = list(exprs[nom])  # ligne a ligne, comme np.ravel
        # Une expression CONSTANTE (un zero de jacobien, le plus souvent)
        # sortirait en scalaire au milieu de tableaux : on les met de cote une
        # fois pour toutes, et les autres sortent d'un seul tableau NumPy.
        variables = [i for i, e in enumerate(elements) if e.free_symbols]
        constantes = np.array(
            [float(e) if not e.free_symbols else 0.0 for e in elements]
        )
        if not variables:

            def fn(X, Z, U, P, _c=constantes):
                return np.repeat(_c[:, None], P.shape[1], axis=1)

        else:
            brut = sp.lambdify(
                exprs["args"],
                sp.Matrix([elements[i] for i in variables]),
                "numpy",
                cse=True,
            )
            iv = np.array(variables)

            def fn(X, Z, U, P, _b=brut, _iv=iv, _c=constantes):
                n = P.shape[1]
                out = np.repeat(_c[:, None], n, axis=1)
                out[_iv] = np.asarray(_b(X, Z, U, P), dtype=float).reshape(len(_iv), n)
                return out

        self._lots[nom] = fn
        return fn

    # -- acces par nom ----------------------------------------------------
    def ix(self, name: str) -> int:
        return self.spec.states.index(name)

    def iz(self, name: str) -> int:
        return self.spec.algebs.index(name)

    def iu(self, name: str) -> int:
        return self.spec.inputs.index(name)

    def iy(self, name: str) -> int:
        return self.spec.outputs.index(name)

    # -- evaluation groupee ------------------------------------------------
    def residual(self, x, z, u, p) -> np.ndarray:
        """Residu DAE complet [f; g] (nul a l'equilibre)."""
        return np.concatenate((self.f(x, z, u, p), self.g(x, z, u, p)))

    def linearize(self, x, z, u, p):
        """Linearisation autour de (x, z, u) : elimine z par complement de
        Schur et renvoie (A, B, C, D) explicites.

            dx/dt = fx.x + fz.z + fu.u ,  0 = gx.x + gz.z + gu.u
            =>  z = -gz^-1 (gx.x + gu.u)
            A = fx - fz gz^-1 gx        B = fu - fz gz^-1 gu
            C = hx - hz gz^-1 gx        D = hu - hz gz^-1 gu
        """
        fx, fz, fu = self.fx(x, z, u, p), self.fz(x, z, u, p), self.fu(x, z, u, p)
        gx, gz, gu = self.gx(x, z, u, p), self.gz(x, z, u, p), self.gu(x, z, u, p)
        hx, hz, hu = self.hx(x, z, u, p), self.hz(x, z, u, p), self.hu(x, z, u, p)
        if self.n_z:
            sol_x = np.linalg.solve(gz, gx)
            sol_u = np.linalg.solve(gz, gu)
            A = fx - fz @ sol_x
            B = fu - fz @ sol_u
            C = hx - hz @ sol_x
            D = hu - hz @ sol_u
        else:
            A, B, C, D = fx, fu, hx, hu
        return A, B, C, D


def _nombre(v: float) -> sp.Expr:
    """``v`` en nombre SymPy, ENTIER s'il l'est : un drapeau nul doit annuler
    son terme (``0.m = 0``), et ``0.0`` flottant ne le ferait pas toujours."""
    return sp.Integer(int(v)) if float(v).is_integer() else sp.Float(v)


_CACHE: dict[str, CompiledModel] = {}


def clear_cache() -> None:
    """Vide le cache de compilation (utile en test / rechargement)."""
    _CACHE.clear()


def compile_spec(spec: ModelSpec, use_cache: bool = True) -> CompiledModel:
    """Derive les jacobiens symboliquement et lambdifie le tout.

    Le resultat est mis en cache sous ``spec.name`` : la deuxieme demande
    pour la meme classe de modele est instantanee.
    """
    if use_cache and spec.name in _CACHE:
        return _CACHE[spec.name]

    xs = sp.symbols(spec.states, real=True) if spec.states else ()
    zs = sp.symbols(spec.algebs, real=True) if spec.algebs else ()
    us = sp.symbols(spec.inputs, real=True) if spec.inputs else ()
    ps = sp.symbols(spec.params, real=True) if spec.params else ()
    # sp.symbols renvoie un scalaire si un seul nom -> normaliser en tuple
    xs, zs, us, ps = (
        (v,) if isinstance(v, sp.Symbol) else tuple(v) for v in (xs, zs, us, ps)
    )

    ns = SimpleNamespace()
    for grp, syms in (
        (spec.states, xs),
        (spec.algebs, zs),
        (spec.inputs, us),
        (spec.params, ps),
    ):
        for nm, s in zip(grp, syms):
            setattr(ns, nm, s)

    f_e, g_e, h_e = spec.equations(ns)
    for label, got, want in (
        ("f", len(f_e), len(spec.states)),
        ("g", len(g_e), len(spec.algebs)),
        ("h", len(h_e), len(spec.outputs)),
    ):
        if got != want:
            raise ValueError(
                f"{spec.name}: {label} a {got} equation(s) pour {want} variable(s)"
            )

    F, G, H = (sp.Matrix(e) if e else sp.zeros(0, 1) for e in (f_e, g_e, h_e))
    X, Z, U = (sp.Matrix(v) if v else sp.zeros(0, 1) for v in (xs, zs, us))

    args = (list(xs), list(zs), list(us), list(ps))

    def _vec(expr, n):
        fn = sp.lambdify(args, expr, "numpy", cse=True)

        def wrapped(x, z, u, p, _fn=fn, _n=n):
            return np.asarray(_fn(x, z, u, p), dtype=float).reshape(_n)

        return wrapped

    def _mat(expr, nr, nc):
        fn = sp.lambdify(args, expr, "numpy", cse=True)

        def wrapped(x, z, u, p, _fn=fn, _r=nr, _c=nc):
            return np.asarray(_fn(x, z, u, p), dtype=float).reshape(_r, _c)

        return wrapped

    nx, nz, nu, ny, npar = (
        len(spec.states),
        len(spec.algebs),
        len(spec.inputs),
        len(spec.outputs),
        len(spec.params),
    )

    jac = dict(
        fx=F.jacobian(X) if nx and nx else sp.zeros(nx, nx),
        fz=F.jacobian(Z) if nx and nz else sp.zeros(nx, nz),
        fu=F.jacobian(U) if nx and nu else sp.zeros(nx, nu),
        gx=G.jacobian(X) if nz and nx else sp.zeros(nz, nx),
        gz=G.jacobian(Z) if nz and nz else sp.zeros(nz, nz),
        gu=G.jacobian(U) if nz and nu else sp.zeros(nz, nu),
        hx=H.jacobian(X) if ny and nx else sp.zeros(ny, nx),
        hz=H.jacobian(Z) if ny and nz else sp.zeros(ny, nz),
        hu=H.jacobian(U) if ny and nu else sp.zeros(ny, nu),
    )
    cm = CompiledModel(
        spec=spec,
        n_x=nx,
        n_z=nz,
        n_u=nu,
        n_y=ny,
        n_p=npar,
        f=_vec(F, nx),
        g=_vec(G, nz),
        h=_vec(H, ny),
        fx=_mat(jac["fx"], nx, nx),
        fz=_mat(jac["fz"], nx, nz),
        fu=_mat(jac["fu"], nx, nu),
        gx=_mat(jac["gx"], nz, nx),
        gz=_mat(jac["gz"], nz, nz),
        gu=_mat(jac["gu"], nz, nu),
        hx=_mat(jac["hx"], ny, nx),
        hz=_mat(jac["hz"], ny, nz),
        hu=_mat(jac["hu"], ny, nu),
        expressions=dict(args=args, f=F, g=G, h=H, **jac),
    )
    if use_cache:
        _CACHE[spec.name] = cm
    return cm
