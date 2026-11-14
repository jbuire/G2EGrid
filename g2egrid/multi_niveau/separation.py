"""Separation d'un modele d'etat global en transport et distribution.

On part du modele lineaire global ``dx/dt = A x + B u`` et d'une partition de
ses etats en ``x_T`` (transport) et ``x_D`` (distribution) ::

    [dx_T]   [A_TT  A_TD] [x_T]   [B_T]
    [dx_D] = [A_DT  A_DD] [x_D] + [B_D] u

``A_TD`` est l'action de la distribution sur le transport, ``A_DT`` celle du
transport sur la distribution. Trois modeles du transport en decoulent :

* ``A_TT`` seul : on neglige toute dynamique de la distribution ;
* l'elimination exacte, qui depend de la frequence ::

      A_eff(s) = A_TT + A_TD (sI - A_DD)^-1 A_DT

* son approximation quasi statique, pour une distribution rapide ::

      A_qs = A_TT - A_TD A_DD^-1 A_DT

Point de vue electrique : la distribution ne voit le transport qu'a travers
quelques etats (la tension du poste de raccordement, et la pulsation du repere
global que lui impose la machine de reference), et le transport ne voit la
distribution qu'a travers le courant qu'elle injecte au poste. ``A_DT`` et
``A_TD`` n'ont donc que quelques colonnes et lignes non nulles, ce qui donne
un modele d'interface ::

    dx_D/dt = A_D x_D + B_D y_T          y_T : etats de transport vus par la distribution
    z_D     = C_D x_D                    z_D : contribution de la distribution a dx_T/dt

C'est ce sous-modele que l'on reduit (troncature equilibree), en conservant
ses entrees et ses sorties d'interface, avant de le reconnecter au transport.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np
from scipy import linalg

__all__ = [
    "SousModele",
    "Separation",
    "separer",
    "troncature_equilibree",
    "reponse_indicielle",
]


@dataclass
class SousModele:
    """Modele d'etat ``dx = A x + B u, y = C x + D u`` avec des noms."""

    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    entrees: list[str]
    sorties: list[str]
    etats: list[str] = field(default_factory=list)

    @property
    def ordre(self) -> int:
        return self.A.shape[0]

    def transfert(self, f: np.ndarray) -> np.ndarray:
        """Matrice de transfert a ``2.pi.f``, forme ``(len(f), ny, nu)``."""
        f = np.atleast_1d(np.asarray(f, dtype=float))
        n = self.ordre
        G = np.empty((f.size, self.C.shape[0], self.B.shape[1]), dtype=complex)
        for k, fk in enumerate(f):
            G[k] = self.C @ np.linalg.solve(
                2j * np.pi * fk * np.eye(n) - self.A, self.B
            )
            G[k] += self.D
        return G


@dataclass
class Separation:
    """Partition d'un modele global et modeles qui en decoulent."""

    etats: list[str]
    iT: np.ndarray  # indices des etats de transport
    iD: np.ndarray  # indices des etats de distribution
    A: np.ndarray
    B: np.ndarray
    entrees: list[str]
    i_vus: np.ndarray  # etats de transport qui excitent la distribution
    i_excites: np.ndarray  # etats de transport excites par la distribution
    distribution: SousModele  # modele d'interface de la distribution

    # --- blocs ------------------------------------------------------------
    @property
    def ATT(self) -> np.ndarray:
        return self.A[np.ix_(self.iT, self.iT)]

    @property
    def ATD(self) -> np.ndarray:
        return self.A[np.ix_(self.iT, self.iD)]

    @property
    def ADT(self) -> np.ndarray:
        return self.A[np.ix_(self.iD, self.iT)]

    @property
    def ADD(self) -> np.ndarray:
        return self.A[np.ix_(self.iD, self.iD)]

    @property
    def etats_T(self) -> list[str]:
        return [self.etats[i] for i in self.iT]

    @property
    def etats_D(self) -> list[str]:
        return [self.etats[i] for i in self.iD]

    # --- modeles du transport --------------------------------------------
    def a_effectif(self, s: complex) -> np.ndarray:
        """Elimination exacte : ``A_TT + A_TD (sI - A_DD)^-1 A_DT``."""
        n = len(self.iD)
        return self.ATT + self.ATD @ np.linalg.solve(s * np.eye(n) - self.ADD, self.ADT)

    def a_quasi_statique(self) -> np.ndarray:
        """Distribution supposee infiniment rapide : ``A_TT - A_TD A_DD^-1 A_DT``."""
        return self.ATT - self.ATD @ np.linalg.solve(self.ADD, self.ADT)

    def modele_global(
        self, dist: SousModele | None = None, entrees: Sequence[str] = ()
    ) -> tuple[np.ndarray, np.ndarray, list[str]]:
        """Transport reconnecte a un modele de distribution (complet par defaut).

        Rend ``(A, B, etats)`` : les etats de transport d'abord, puis ceux du
        modele de distribution. ``B`` porte les entrees globales ``entrees`` :
        leur action sur le transport vient du modele global, leur action sur
        la distribution du modele ``dist``, qui doit alors les avoir parmi ses
        entrees (voir ``entrees`` de :func:`separer`).
        """
        dist = dist or self.distribution
        nT, nD = len(self.iT), dist.ordre
        n_vus = len(self.i_vus)
        # entrees de dist : d'abord les etats de transport vus, puis les
        # entrees globales gardees ; sorties : les derivees excitees
        B_vus, B_ext = dist.B[:, :n_vus], dist.B[:, n_vus:]
        D_vus = dist.D[:, :n_vus]
        E_vus = np.zeros((n_vus, nT))
        E_vus[np.arange(n_vus), self.i_vus] = 1.0
        E_exc = np.zeros((nT, len(self.i_excites)))
        E_exc[self.i_excites, np.arange(len(self.i_excites))] = 1.0

        A = np.block(
            [
                [self.ATT + E_exc @ D_vus @ E_vus, E_exc @ dist.C],
                [B_vus @ E_vus, dist.A],
            ]
        )
        cols = [self.entrees.index(n) for n in entrees]
        ext = dist.entrees[n_vus:]
        B = np.zeros((nT + nD, len(entrees)))
        B[:nT] = self.B[np.ix_(self.iT, cols)]
        for j, nom in enumerate(entrees):
            if nom in ext:
                k = ext.index(nom)
                B[nT:, j] = B_ext[:, k]
                B[:nT, j] += E_exc @ dist.D[:, n_vus + k]
        etats = self.etats_T + (dist.etats or [f"xr{k}" for k in range(nD)])
        return A, B, etats


def separer(
    ss,
    distribution: Callable[[str], bool],
    entrees: Sequence[str] = (),
    tol: float = 1e-12,
) -> Separation:
    """Partitionne le modele d'etat ``ss`` (``StateSpace`` de ``g2egrid``).

    ``distribution(nom_etat)`` dit si un etat appartient a la distribution.
    ``entrees`` : entrees globales (``ss.inputs``) que le modele de
    distribution doit garder, par exemple les consignes de ses GFL.
    """
    etats = list(ss.states)
    masque = np.array([distribution(e) for e in etats])
    iD, iT = np.flatnonzero(masque), np.flatnonzero(~masque)
    A, B = np.asarray(ss.A), np.asarray(ss.B)
    ATD, ADT = A[np.ix_(iT, iD)], A[np.ix_(iD, iT)]
    seuil = tol * np.abs(A).max()
    i_vus = np.flatnonzero(np.abs(ADT).max(axis=0) > seuil)
    i_excites = np.flatnonzero(np.abs(ATD).max(axis=1) > seuil)

    noms_T = [etats[i] for i in iT]
    cols = [list(ss.inputs).index(n) for n in entrees]
    dist = SousModele(
        A=A[np.ix_(iD, iD)],
        B=np.hstack([ADT[:, i_vus], B[np.ix_(iD, cols)]]),
        C=ATD[i_excites, :],
        D=np.zeros((len(i_excites), len(i_vus) + len(cols))),
        entrees=[noms_T[i] for i in i_vus] + list(entrees),
        sorties=["d/dt " + noms_T[i] for i in i_excites],
        etats=[etats[i] for i in iD],
    )
    return Separation(
        etats=etats,
        iT=iT,
        iD=iD,
        A=A,
        B=B,
        entrees=list(ss.inputs),
        i_vus=i_vus,
        i_excites=i_excites,
        distribution=dist,
    )


# ----------------------------------------------------------------------
# reduction
# ----------------------------------------------------------------------
def _separer_stable(A, B, C, marge):
    """Decompose ``(A, B, C)`` en partie stable et partie lente ou instable.

    Forme de Schur reelle ordonnee, puis decouplage par une equation de
    Sylvester : ``A = blocdiag(A_s, A_i)`` dans une base ou les deux parties
    ne se voient plus. Les modes de partie reelle superieure a ``-marge``
    (modes nuls, integrateurs, instables) sont gardes tels quels.
    """
    T, Z, ns = linalg.schur(A, output="real", sort=lambda re, im: re < -marge)
    if ns == A.shape[0]:
        return (A, B, C), None
    T11, T12, T22 = T[:ns, :ns], T[:ns, ns:], T[ns:, ns:]
    X = linalg.solve_sylvester(T11, -T22, -T12)  # T11 X - X T22 = -T12
    n = A.shape[0]
    W = np.eye(n)
    W[:ns, ns:] = X
    Winv = np.eye(n)
    Winv[:ns, ns:] = -X
    Bt, Ct = Winv @ Z.T @ B, C @ Z @ W
    stable = (T11, Bt[:ns], Ct[:, :ns])
    lente = (T22, Bt[ns:], Ct[:, ns:])
    return stable, lente


def troncature_equilibree(
    m: SousModele, ordre: int, methode: str = "residualisation", marge: float = 1e-6
) -> tuple[SousModele, np.ndarray]:
    """Reduction equilibree, entrees et sorties conservees.

    Les modes lents ou instables (partie reelle > ``-marge``) sont gardes
    exactement ; la partie stable est mise sous forme equilibree (methode de la
    racine carree), puis ses etats les moins commandables et observables sont
    elimines, de l'une des deux facons :

    * ``"troncature"`` : on les supprime. Exact en haute frequence, mais le
      gain statique est faux ;
    * ``"residualisation"`` : on les suppose a l'equilibre (``dx = 0``). Le
      gain statique est exact, ce qui compte ici : les modes du transport qui
      nous interessent sont lents devant la distribution. A l'ordre 0, c'est
      l'approximation quasi statique.

    ``ordre`` est l'ordre TOTAL du modele reduit. Rend le modele reduit et les
    valeurs singulieres de Hankel de la partie stable.
    """
    if methode not in ("troncature", "residualisation"):
        raise ValueError(f"methode {methode!r} inconnue")
    stable, lente = _separer_stable(m.A, m.B, m.C, marge)
    As, Bs, Cs = stable
    n_lent = 0 if lente is None else lente[0].shape[0]
    r = ordre - n_lent
    if r < 0:
        raise ValueError(f"ordre {ordre} < {n_lent} modes lents a conserver")

    P = linalg.solve_continuous_lyapunov(As, -Bs @ Bs.T)
    Q = linalg.solve_continuous_lyapunov(As.T, -Cs.T @ Cs)
    Lp, Lq = _facteur(P), _facteur(Q)
    U, hsv, Vt = np.linalg.svd(Lq.T @ Lp)
    k = int(np.sum(hsv > hsv[0] * 1e-14))  # partie minimale, numeriquement
    S = np.diag(hsv[:k] ** -0.5)
    Tg = Lp @ Vt[:k].T @ S  # base equilibree (droite)
    Tl = S @ U[:, :k].T @ Lq.T  # projection (gauche)
    Ab, Bb, Cb = Tl @ As @ Tg, Tl @ Bs, Cs @ Tg
    Dr = m.D.copy()
    if methode == "troncature" or r >= k:
        Ar, Br, Cr = Ab[:r, :r], Bb[:r], Cb[:, :r]
    else:
        A11, A12, A21, A22 = Ab[:r, :r], Ab[:r, r:], Ab[r:, :r], Ab[r:, r:]
        X = np.linalg.solve(A22, np.hstack([A21, Bb[r:]]))
        Ar = A11 - A12 @ X[:, :r]
        Br = Bb[:r] - A12 @ X[:, r:]
        Cr = Cb[:, :r] - Cb[:, r:] @ X[:, :r]
        Dr = Dr - Cb[:, r:] @ X[:, r:]

    if lente is not None:
        Al, Bl, Cl = lente
        Ar = linalg.block_diag(Ar, Al)
        Br = np.vstack([Br, Bl])
        Cr = np.hstack([Cr, Cl])
    reduit = SousModele(
        A=Ar,
        B=Br,
        C=Cr,
        D=Dr,
        entrees=list(m.entrees),
        sorties=list(m.sorties),
        etats=[f"xr{k}" for k in range(Ar.shape[0])],
    )
    return reduit, hsv


def _facteur(M: np.ndarray) -> np.ndarray:
    """Facteur ``L`` tel que ``M = L L'``, robuste a une matrice semi-definie."""
    M = 0.5 * (M + M.T)
    w, V = np.linalg.eigh(M)
    return V * np.sqrt(np.clip(w, 0.0, None))


# ----------------------------------------------------------------------
# simulation
# ----------------------------------------------------------------------
def reponse_indicielle(
    A: np.ndarray, B: np.ndarray, t: np.ndarray, u: np.ndarray | None = None
) -> np.ndarray:
    """Reponse de ``dx = A x + B u`` a un echelon ``u`` depuis ``x = 0``.

    Discretisation exacte (bloqueur d'ordre zero) par l'exponentielle de la
    matrice augmentee ``[[A, B], [0, 0]]`` : valable pour des modes aussi
    rapides que l'on veut et pour une matrice ``A`` singuliere. ``t`` doit
    etre a pas constant. Rend ``x`` de forme ``(len(t), n)``.
    """
    n, m = B.shape
    u = np.ones(m) if u is None else np.asarray(u, dtype=float)
    dt = t[1] - t[0]
    M = np.zeros((n + m, n + m))
    M[:n, :n], M[:n, n:] = A, B
    E = linalg.expm(M * dt)
    Ad, Bd = E[:n, :n], E[:n, n:] @ u
    x = np.zeros((len(t), n))
    for k in range(1, len(t)):
        x[k] = Ad @ x[k - 1] + Bd
    return x
