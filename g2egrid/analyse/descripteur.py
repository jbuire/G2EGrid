"""Forme descripteur creuse ``(E, J)`` et valeurs propres par shift-invert.

Le modele d'etat ``A`` (:meth:`Network.build_state_space`) s'obtient en
eliminant les variables algebriques -- complement de Schur, ``BZ^-1``. Cette
elimination REMPLIT : une matrice de cablage a quelques termes par ligne
devient une ``A`` largement pleine, et ``eig(A)`` coute ``O(n^3)`` sur tout le
spectre quand on ne veut que la poignee de modes proches de l'axe
imaginaire. On garde donc le DAE lineaire TEL QUEL ::

    E . dx/dt = J . x,      E = [I 0; 0 0],      J = [AX AZ; BX BZ]

dont les valeurs propres FINIES sont exactement celles de ``A`` (les autres
sont infinies : ce sont les contraintes). ``J`` reste creuse : chaque
equipement ne touche que ses propres lignes et celles de ses bus.

Les valeurs propres proches d'un point ``sigma`` du plan complexe sortent
d'Arnoldi sur l'operateur ``(J - sigma.E)^-1 . E``, dont les valeurs propres
sont ``1/(lam - sigma)`` : les ``lam`` proches de ``sigma`` deviennent les plus
grandes, les infinies tombent a zero, ou Arnoldi ne les cherche pas. Une seule
factorisation LU creuse de ``J - sigma.E`` par point.

:func:`balayage` repete l'operation le long d'une droite ``Re = sigma_r`` a
plusieurs frequences, pour couvrir une bande sans rien manquer : c'est la
recherche des modes critiques d'un grand reseau.

Les vecteurs propres a GAUCHE, pour les participations, s'obtiennent par
iteration inverse sur la transposee conjuguee, avec la meme LU. Pour la
partie differentielle, ce sont exactement les vecteurs a gauche de ``A`` (la
partie algebrique ``w_z = -BZ^-H.AZ^H.w_x`` les complete) : les participations
sont celles de :func:`g2egrid.analyse.modal.spectre`, sans avoir forme ``A``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

__all__ = [
    "Descripteur",
    "ModesDescripteur",
    "balayage",
    "descripteur",
    "modes_proches",
]


@dataclass
class Descripteur:
    """Le faisceau ``J - lam.E`` d'un reseau, creux."""

    J: sp.csc_matrix
    E: sp.csc_matrix
    etats: list[str]  # les nx premieres inconnues
    algebriques: list[str]  # les nz suivantes

    @property
    def nx(self) -> int:
        return len(self.etats)

    @property
    def n(self) -> int:
        return self.J.shape[0]

    @property
    def densite(self) -> float:
        return self.J.nnz / max(self.n**2, 1)


def descripteur(net, tol: float = 0.0) -> Descripteur:
    """Le faisceau du reseau ``net`` au point d'equilibre (resolu si besoin).
    ``tol`` : coefficients de ``J`` en dessous de ``tol.max|J|`` ecartes."""
    if any(el.x0 is None for el in net.all_elements):
        net.solve_equilibrium()
    d = net.dae_lineaire()
    nx, nz = d.AX.shape[0], d.BZ.shape[0]

    def creux(M):
        M = sp.csr_matrix(M)
        if tol > 0 and M.nnz:
            M.data[np.abs(M.data) < tol * np.abs(M.data).max()] = 0.0
            M.eliminate_zeros()
        return M

    J = (
        sp.bmat([[creux(d.AX), creux(d.AZ)], [creux(d.BX), creux(d.BZ)]], format="csc")
        if nz
        else creux(d.AX).tocsc()
    )
    E = sp.diags(np.r_[np.ones(nx), np.zeros(nz)], format="csc")
    return Descripteur(J=J, E=E, etats=list(d.etats), algebriques=list(d.algebriques))


@dataclass
class ModesDescripteur:
    lam: np.ndarray  # valeurs propres finies trouvees
    V: np.ndarray  # vecteurs a droite, partie differentielle (nx, k)
    p: np.ndarray  # participations (k, nx), somme 1 par mode
    etats: list[str]
    #: toutes les valeurs propres finies a moins de ``rayon`` du sigma demande
    #: ont ete trouvees (inf : le faisceau n'en a pas d'autres)
    rayon: float = np.inf

    def __len__(self) -> int:
        return len(self.lam)


def _lu(desc: Descripteur, sigma: complex):
    M = (desc.J - sigma * desc.E).astype(complex).tocsc()
    echelle = max(float(np.abs(desc.J.data).max(initial=1.0)), 1.0)
    for decalage in (0.0, 1e-9 * echelle, 1e-6 * echelle):
        try:
            return spla.splu(M - decalage * desc.E if decalage else M), sigma + decalage
        except RuntimeError as err:  # "Factor is exactly singular"
            if "singular" not in str(err).lower():
                raise
            derniere = err
    raise RuntimeError(f"J - sigma.E singuliere autour de sigma = {sigma} : {derniere}")


def _gauche(lu, desc: Descripteur, v: np.ndarray, iterations: int = 3) -> np.ndarray:
    """Vecteur a gauche par iteration inverse : ``(J - mu.E)^H w = E^H w``."""
    w = desc.E @ v  # depart non orthogonal au vecteur cherche, en general
    if not np.any(w):
        w = np.ones(desc.n, dtype=complex)
    for _ in range(iterations):
        w = lu.solve(desc.E.T @ w, trans="H")
        w /= np.linalg.norm(w)
    return w


def modes_proches(
    desc: Descripteur,
    sigma: complex = 0.0,
    k: int = 6,
    gauche: bool = True,
    residu_max: float = 1e-8,
) -> ModesDescripteur:
    """Les ``k`` valeurs propres finies les plus proches de ``sigma`` (moins
    si certaines n'atteignent pas un residu relatif ``residu_max``)."""
    nx = desc.nx
    if not 0 < k < nx - 1:
        raise ValueError(f"k = {k} hors de ]0, nx-1[ avec nx = {nx}")
    # sigma = 0 tombe EXACTEMENT sur le mode de jauge d'un reseau sans bus
    # infini : (J - sigma.E)^-1 y prend une valeur propre de 1e12, et Arnoldi,
    # dont la precision est relative a la plus grande, rend alors du bruit pour
    # tous les autres modes. On regarde un peu a cote ; le tri final se fait
    # toujours par distance au sigma demande.
    s0 = sigma + (0.02 + 0.03j) if abs(sigma) < 0.05 else sigma
    lu, s = _lu(desc, s0)
    E = desc.E.astype(complex)
    op = spla.LinearOperator(
        (desc.n, desc.n), matvec=lambda x: lu.solve(E @ x), dtype=complex
    )
    nu, X = spla.eigs(op, k=k, which="LM")
    garde = np.abs(nu) > 1e-12  # infini : nu = 0
    lam = s + 1.0 / nu[garde]
    X = X[:, garde]
    # Arnoldi rend les k valeurs les plus proches de s : tout ce qui est plus
    # pres de s que la plus lointaine est la. Une valeur infinie parmi elles
    # veut dire qu'il n'y a plus de valeur finie au-dela.
    rayon = np.inf if not garde.all() else float(np.abs(lam - s).max())
    # garde-fou : le residu relatif de chaque paire (J - lam.E) x
    normeJ = spla.norm(desc.J, 1)
    res = np.array(
        [
            np.linalg.norm(desc.J @ X[:, j] - lam[j] * (desc.E @ X[:, j]))
            / ((normeJ + abs(lam[j])) * np.linalg.norm(X[:, j]))
            for j in range(len(lam))
        ]
    )
    ok = res < residu_max
    if not ok.all():  # une valeur rejetee : la couverture s'arrete avant elle
        rayon = min(rayon, float(np.abs(lam[~ok] - s).min()))
    rayon = max(rayon - abs(s - sigma), 0.0)
    lam, X = lam[ok], X[:, ok]
    ordre = np.argsort(np.abs(lam - sigma))
    lam, X = lam[ordre], X[:, ordre]
    V = X[:nx]
    p = np.zeros((len(lam), nx))
    if gauche:
        for j in range(len(lam)):
            luj, _ = _lu(desc, lam[j] + 1e-8 * max(abs(lam[j]), 1.0))
            w = _gauche(luj, desc, X[:, j])
            pj = np.abs(np.conj(w[:nx]) * V[:, j])
            p[j] = pj / max(pj.sum(), 1e-300)
    return ModesDescripteur(lam=lam, V=V, p=p, etats=list(desc.etats), rayon=rayon)


def balayage(
    desc: Descripteur,
    f_max: float,
    f_min: float = 0.0,
    sigma_r: float = 0.0,
    profondeur: float = 0.5,
    k: int = 8,
    pas: int | None = None,
    gauche: bool = True,
) -> ModesDescripteur:
    """TOUS les modes du rectangle ``|Re lam - sigma_r| <= profondeur``,
    ``f_min <= Im lam / 2.pi <= f_max``, un par paire conjuguee (``Im >= 0``),
    tries par partie reelle decroissante (les plus critiques d'abord).

    Les points ``sigma_r + 2j.pi.f`` sont espaces de ``df`` (``pas`` points ;
    par defaut un tous les 0,5 Hz). Chacun doit couvrir le disque de rayon
    ``sqrt((pi.df)^2 + profondeur^2)`` : tant que les valeurs trouvees ne
    l'atteignent pas, ``k`` double. Le rectangle est donc couvert quel que soit
    le nombre de modes qu'il contient -- une grappe de modes de regulateurs
    vers -0,2 s^-1 en compte une dizaine sur un grand reseau."""
    pas = pas or max(2, int(np.ceil((f_max - f_min) / 0.5)) + 1)
    freqs = np.linspace(f_min, f_max, pas)
    df = (f_max - f_min) / (pas - 1)
    besoin = float(np.hypot(np.pi * df, profondeur))
    trouves: list[ModesDescripteur] = []
    for f in freqs:
        kk = k
        while True:
            m = modes_proches(desc, sigma_r + 2j * np.pi * f, kk, gauche)
            if m.rayon >= besoin or kk >= desc.nx - 2:
                break
            kk = min(2 * kk, desc.nx - 2)
        trouves.append(m)
    # dedoublonnage : une valeur trouvee par plusieurs points ne compte qu'une
    # fois, mais une valeur MULTIPLE garde sa multiplicite (le plus grand
    # nombre de copies vues par un meme point)
    groupes: list[tuple[complex, int, int]] = []  # (lam, point, indice)
    for a, m in enumerate(trouves):
        for i, lam in enumerate(m.lam):
            if lam.imag < -1e-9 * max(abs(lam), 1.0):
                continue
            if abs(lam.real - sigma_r) > profondeur:
                continue
            if not f_min - 1e-9 <= lam.imag / (2 * np.pi) <= f_max + 1e-9:
                continue
            groupes.append((lam, a, i))
    retenus: list[tuple[int, int]] = []
    vus: list[complex] = []
    for lam, _a, _i in groupes:
        tol = 1e-6 * max(abs(lam), 1.0)
        if any(abs(lam - v) <= tol for v in vus):
            continue
        vus.append(lam)
        meme = [(a, i) for l2, a, i in groupes if abs(l2 - lam) <= tol]
        par_point: dict[int, list[int]] = {}
        for a, i in meme:
            par_point.setdefault(a, []).append(i)
        a, idx = max(par_point.items(), key=lambda t: len(t[1]))
        retenus += [(a, i) for i in idx]
    lam = np.array([trouves[a].lam[i] for a, i in retenus], dtype=complex)
    V = np.array([trouves[a].V[:, i] for a, i in retenus]).T.reshape(desc.nx, -1)
    p = np.array([trouves[a].p[i] for a, i in retenus]).reshape(-1, desc.nx)
    ordre = np.argsort(-lam.real)
    return ModesDescripteur(
        lam=lam[ordre], V=V[:, ordre], p=p[ordre], etats=list(desc.etats)
    )
