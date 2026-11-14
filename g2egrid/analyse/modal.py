"""Analyse modale : valeurs propres, participations, et ce qui s'en deduit.

Les notebooks calculent tous la meme chose a la main. Ce module en fait des
fonctions, pour que les notebooks, les etudes, les tests et l'interface web
partent du meme code.

Convention des facteurs de participation, celle des notebooks ::

    lam, Vr = eig(A)        # Vr : vecteurs propres a DROITE, en colonnes
    Wl = inv(Vr)            # ses lignes sont les vecteurs propres a GAUCHE
    p[k, i] = |Vr[i, k]| . |Wl[k, i]|      puis normalisation par mode

Prendre ``Wl = inv(Vr)`` apparie de lui-meme chaque vecteur a gauche avec
celui de droite du meme mode. Les calculer par deux ``eig`` separes, sur ``A``
puis sur ``A^H``, ne garantit aucun ordre commun : c'est un piege classique,
qui produit des participations fausses -- silencieusement -- des que deux
modes sont proches.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

__all__ = [
    "Mode",
    "amortissement",
    "forme_modale",
    "frequence",
    "modes",
    "participations",
    "reponse_indicielle_modale",
    "reponse_libre",
    "sensibilite",
    "sensibilite_top",
    "spectre",
    "valeurs_propres_partielles",
]


@dataclass
class Mode:
    """Un mode, avec les etats qui y participent le plus."""

    indice: int
    valeur: complex
    f: float  # frequence amortie [Hz]
    zeta: float  # amortissement [%]
    etats: list[
        tuple[str, float]
    ]  # (nom, participation en %), du plus fort au plus faible

    @property
    def stable(self) -> bool:
        return self.valeur.real < 0

    @property
    def oscillant(self) -> bool:
        return abs(self.valeur.imag) > 1e-6

    @property
    def dominant(self) -> str:
        return self.etats[0][0] if self.etats else ""


def _A(ss) -> np.ndarray:
    """Accepte une matrice ou un ``StateSpace``."""
    return ss if isinstance(ss, np.ndarray) else ss.A


def _noms(ss, etats: list[str] | None) -> list[str]:
    if etats is not None:
        return etats
    if not isinstance(ss, np.ndarray):
        return list(ss.states)
    return [f"x{i}" for i in range(len(_A(ss)))]


def _decomposition(ss) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(lam, Vr, Wl)`` avec ``Wl = inv(Vr)``, donc ``Wl @ Vr = I``.

    Tout ce module part de la : la bi-orthonormalisation est acquise, ce qui
    dispense de diviser par ``w_k . v_k`` partout (le ``coef`` du MATLAB de
    G2Elin V1) et evite surtout d'apparier a la main deux ``eig`` separes.
    """
    lam, Vr = np.linalg.eig(_A(ss))
    return lam, Vr, np.linalg.inv(Vr)


def amortissement(lam) -> np.ndarray:
    """Amortissement ``zeta = -100 * Re(lam) / |lam|`` [%], vectoriel.

    Une valeur propre nulle (le mode de jauge d'un reseau sans bus infini)
    n'a pas d'amortissement defini : on la compte a 100 %, ce qui la range
    avec les modes non oscillants plutot qu'avec les modes critiques.
    """
    lam = np.atleast_1d(np.asarray(lam, dtype=complex))
    mag = np.abs(lam)
    ok = mag > 1e-12
    return np.where(ok, -100.0 * lam.real / np.where(ok, mag, 1.0), 100.0)


def frequence(lam) -> np.ndarray:
    """Frequence AMORTIE ``|Im(lam)| / 2.pi`` [Hz], vectoriel.

    C'est la frequence que l'on lit sur un enregistrement, pas la frequence
    propre ``|lam| / 2.pi`` : les deux ne se confondent qu'a faible
    amortissement.
    """
    lam = np.atleast_1d(np.asarray(lam, dtype=complex))
    return np.abs(lam.imag) / (2 * np.pi)


def spectre(ss) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(valeurs propres, vecteurs a droite, participations normalisees)``.

    ``p[k, i]`` est la participation de l'etat ``i`` au mode ``k``, en
    fraction (somme 1 par mode).
    """
    lam, Vr, Wl = _decomposition(ss)
    p = np.abs(Vr.T) * np.abs(Wl)
    p = p / np.maximum(p.sum(axis=1, keepdims=True), 1e-30)
    return lam, Vr, p


def spectre_complet(ss) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """:func:`spectre` lu dans la base COMPLETE des etats, et ses noms.

    Aux noeuds ideaux, la reduction d'indice elimine un courant par
    contrainte de serie (``x = T.x_r``, :attr:`StateSpace.elimines`). Lues sur
    les seuls etats gardes, les participations donnent au courant garde la
    part de celui qu'il represente : le reseau perd la moitie de la sienne
    quand deux branches sont en serie, et un classement par groupes d'etats
    (:mod:`g2egrid.classement`) change de categorie a valeur propre identique
    (FR-044 : 17 modes instables sur 288 dans la campagne).

    Dans la base complete, ``v = T.v_r`` et ``w = w_r.T_g``, ou ``T_g`` est
    l'inverse a gauche de ``T`` dont le noyau porte les directions RAPIDES
    des noeuds ideaux (:attr:`StateSpace.T_gauche`) : ``w.v = w_r.v_r = 1``
    reste vrai, et ``p`` est la limite, quand ``R`` tend vers l'infini, des
    participations du modele a resistance de causalite (un mode lent y est
    orthogonal au mode parasite en ``-R``). La pseudo-inverse de ``T``,
    elle, repartit la part d'un courant elimine sur TOUS ceux de sa
    contrainte : au noeud d'un GFM, une part de ligne passait a sa boucle de
    courant, et une resonance devenait une interaction de convertisseur.
    Sans reduction, c'est :func:`spectre`.
    """
    lam, Vr, Wl = _decomposition(ss)
    T = getattr(ss, "T", None)
    if T is None:
        p = np.abs(Vr.T) * np.abs(Wl)
        noms = _noms(ss, None)
    else:
        Tg = getattr(ss, "T_gauche", None)
        if Tg is None:  # a defaut : la pseudo-inverse
            Tg = np.linalg.solve(T.T @ T, T.T)
        Vr = T @ Vr
        Wl = Wl @ Tg
        p = np.abs(Vr.T) * np.abs(Wl)
        noms = list(ss.etats_complets)
    p = p / np.maximum(p.sum(axis=1, keepdims=True), 1e-30)
    return lam, Vr, p, noms


def participations(ss) -> np.ndarray:
    """Seulement la matrice des participations (cf. :func:`spectre`)."""
    return spectre(ss)[2]


def modes(
    ss,
    etats: list[str] | None = None,
    f_max: float | None = None,
    zeta_max: float | None = None,
    n_etats: int = 3,
    oscillants: bool = False,
) -> list[Mode]:
    """Les modes, du moins amorti au plus amorti.

    ``f_max`` borne la frequence amortie retenue [Hz], ``zeta_max``
    l'amortissement [%] (au-dela d'environ 95 %, une paire complexe ne fait
    plus une oscillation), ``oscillants`` ne garde que les paires complexes
    de partie imaginaire positive.
    """
    lam, _, p = spectre(ss)
    noms = _noms(ss, etats)
    zeta = amortissement(lam)
    f = frequence(lam)

    garde = np.ones(len(lam), dtype=bool)
    if f_max is not None:
        garde &= f <= f_max
    if zeta_max is not None:
        garde &= zeta <= zeta_max
    if oscillants:
        garde &= lam.imag > 1e-6

    out = []
    for k in np.flatnonzero(garde):
        top = np.argsort(-p[k])[:n_etats]
        out.append(
            Mode(
                indice=int(k),
                valeur=complex(lam[k]),
                f=float(f[k]),
                zeta=float(zeta[k]),
                etats=[(noms[i], float(100 * p[k, i])) for i in top],
            )
        )
    out.sort(key=lambda m: m.zeta)
    return out


# ----------------------------------------------------------------------
# sensibilite aux coefficients de A
# ----------------------------------------------------------------------
def sensibilite(ss, k: int) -> np.ndarray:
    """``dlambda_k / dA[i, j] = Wl[k, i]*Vr[j, k]``, masquee aux zeros de ``A``.

    Le masque n'est pas de la cosmetique : un coefficient nul de ``A`` l'est
    parce que les deux etats ne sont relies par aucune equation du modele, et
    aucun parametre physique ne le fera bouger. Sa sensibilite est
    parfaitement definie mathematiquement et parfaitement inutilisable : sans
    le masque, les plus grandes valeurs du tableau designent regulierement des
    couples d'etats qui n'ont aucun chemin entre eux. C'est ce que fait le
    ``A./A`` de ``modal_analysis_sens.m`` de G2Elin V1.

    Rend une matrice complexe ``(n, n)`` : la partie reelle deplace
    l'amortissement du mode, la partie imaginaire sa frequence.
    """
    lam, Vr, Wl = _decomposition(ss)
    if not -len(lam) <= k < len(lam):
        raise IndexError(f"mode {k} hors du spectre ({len(lam)} modes)")
    S = np.outer(Wl[k], Vr[:, k])
    return S * (_A(ss) != 0)


def sensibilite_top(
    ss, k: int, n: int = 8, etats: list[str] | None = None
) -> list[tuple[str, str, float]]:
    """Les ``n`` couples ``(etat ligne, etat colonne, |sensibilite|)`` dominants.

    La lecture est celle du modele : ``A[i, j]`` est l'effet de l'etat ``j``
    sur la derivee de l'etat ``i``. Le couple rendu designe donc un
    coefficient d'equation, pas une paire d'etats symetrique -- ``(i, j)`` et
    ``(j, i)`` n'ont aucune raison de sortir ensemble.
    """
    S = np.abs(sensibilite(ss, k))
    noms = _noms(ss, etats)
    plats = np.argsort(-S, axis=None)[: max(n, 0)]
    i, j = np.unravel_index(plats, S.shape)
    return [(noms[a], noms[b], float(S[a, b])) for a, b in zip(i, j)]


# ----------------------------------------------------------------------
# forme modale
# ----------------------------------------------------------------------
def forme_modale(
    ss, k: int, n: int = 5, etats: list[str] | None = None
) -> list[tuple[str, float, float]]:
    """``(nom, module, phase [deg])`` des ``n`` etats les plus participants.

    Le vecteur propre a droite n'est defini qu'a un facteur complexe pres :
    ``eig`` le rend de norme 1, mais de phase arbitraire, qui change d'une
    machine a l'autre et d'une version de LAPACK a l'autre. On divise donc
    tout par la composante de l'etat le plus participant : celle-ci vaut
    ``(1, 0 deg)`` et les autres se lisent RELATIVEMENT a elle. C'est ce qui
    rend le diagramme polaire interpretable -- deux vitesses rotor a 180 deg
    l'une de l'autre, c'est un mode inter-machines.

    Les phases sont dans ``]-180, 180]``.
    """
    lam, Vr, p = spectre(ss)
    if not -len(lam) <= k < len(lam):
        raise IndexError(f"mode {k} hors du spectre ({len(lam)} modes)")
    noms = _noms(ss, etats)
    top = np.argsort(-p[k])[: max(n, 0)]
    v = Vr[top, k]
    v = v / v[0]  # l'etat dominant sert de reference de phase
    return [
        (noms[i], float(np.abs(c)), float(np.degrees(np.angle(c))))
        for i, c in zip(top, v)
    ]


# ----------------------------------------------------------------------
# reponses en forme close
# ----------------------------------------------------------------------
def reponse_libre(ss, x0: np.ndarray, t: np.ndarray) -> np.ndarray:
    """``x(t) = Vr @ diag(exp(lam*t)) @ Wl @ x0``, sans integration numerique.

    Le developpement modal donne la reponse a un instant quelconque au prix
    d'une exponentielle scalaire : pas de pas de temps a choisir, donc pas de
    raideur. Un reseau dq porte des modes a 1e7 rad/s a cote de modes
    electromecaniques a 10 rad/s ; un integrateur explicite s'y etrangle, la
    forme close l'ignore.

    ``t`` n'a pas a etre a pas constant. Rend un tableau REEL
    ``(len(t), n_x)`` : les parties imaginaires se compensent deux a deux
    entre paires conjuguees, et ce qui reste est du bruit d'arrondi.
    """
    lam, Vr, Wl = _decomposition(ss)
    x0 = np.asarray(x0, dtype=float).ravel()
    t = np.atleast_1d(np.asarray(t, dtype=float))
    if x0.size != len(lam):
        raise ValueError(f"x0 de taille {x0.size}, attendu {len(lam)}")
    c = Wl @ x0
    return (Vr @ (c[:, None] * np.exp(np.outer(lam, t)))).T.real


def reponse_indicielle_modale(
    ss,
    u: np.ndarray,
    t: np.ndarray,
    par_mode: bool = False,
    etats: list[str] | list[int] | None = None,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Reponse a un echelon d'entree depuis ``x = 0``, mode par mode ::

        x(t) = somme_k  v_k . (w_k . B . u) . (exp(lam_k.t) - 1) / lam_k

    L'interet n'est pas la reponse elle-meme -- :func:`g2egrid.multi_niveau.separation.
    reponse_indicielle` la donne par discretisation exacte, et plus vite --
    mais sa decomposition : chaque terme de la somme est la part d'un mode, ce
    qui dit lequel porte le depassement ou l'oscillation residuelle.

    Un mode de valeur propre nulle, le mode de jauge d'un reseau sans bus
    infini, ferait diverger ``1 / lam_k`` ; sa contribution est le passage a
    la limite, ``t``, c'est-a-dire une derive lineaire. Elle est bien reelle :
    l'echelon deplace l'angle de reference sans retour.

    ``etats`` (noms ou indices) restreint les colonnes rendues. Avec
    ``par_mode``, s'en passer coute ``n_modes x len(t) x n_x`` flottants, soit
    plusieurs centaines de Mo des la centaine d'etats : sur un vrai reseau, on
    selectionne.

    Rend ``x`` de forme ``(len(t), n_sel)``, et avec ``par_mode`` le couple
    ``(x, contributions)``, ``contributions`` de forme
    ``(n_modes, len(t), n_sel)`` et ``x = contributions.sum(axis=0)``.
    """
    if isinstance(ss, np.ndarray):
        raise TypeError("reponse_indicielle_modale veut un StateSpace : il lui faut B")
    lam, Vr, Wl = _decomposition(ss)
    u = np.atleast_1d(np.asarray(u, dtype=float)).ravel()
    t = np.atleast_1d(np.asarray(t, dtype=float))
    B = np.asarray(ss.B)
    if u.size != B.shape[1]:
        raise ValueError(f"u de taille {u.size}, attendu {B.shape[1]}")

    cols = _colonnes(ss, etats)
    c = Wl @ (B @ u)
    petit = np.abs(lam) <= 1e-12 * max(float(np.abs(lam).max()), 1.0)
    denom = np.where(petit, 1.0, lam)[:, None]
    psi = np.where(petit[:, None], t[None, :], (np.exp(np.outer(lam, t)) - 1.0) / denom)
    modal = c[:, None] * psi  # (n_modes, n_t)

    if not par_mode:
        return (Vr[cols] @ modal).T.real
    # (n_modes, n_t, n_sel) : chaque mode garde sa propre trajectoire
    contributions = np.einsum("ik,kt->kti", Vr[cols], modal).real
    return contributions.sum(axis=0), contributions


def _colonnes(ss, etats) -> np.ndarray:
    """Indices des etats demandes (noms ou indices), tous par defaut."""
    if etats is None:
        return np.arange(len(_A(ss)))
    noms = _noms(ss, None)
    return np.array([noms.index(e) if isinstance(e, str) else int(e) for e in etats])


# ----------------------------------------------------------------------
# spectre partiel
# ----------------------------------------------------------------------
def valeurs_propres_partielles(
    ss, k: int, sigma: complex = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """Les ``k`` valeurs propres les plus proches de ``sigma``, shift-and-invert.

    ``numpy.linalg.eig`` calcule TOUT le spectre, en ``O(n^3)`` et avec une
    matrice de vecteurs propres pleine : passe quelques milliers d'etats -- un
    reseau de distribution complet raccorde a un reseau de transport -- ce
    n'est plus tenable, alors qu'on ne veut jamais que la poignee de modes
    electromecaniques autour de ``sigma``. Arnoldi sur ``(A - sigma.I)^-1``
    les sort directement : cette transformation envoie a l'infini les valeurs
    propres proches de ``sigma``, c'est-a-dire la ou Arnoldi converge vite.

    ``sigma`` se choisit dans le plan complexe : ``0`` pour les modes lents,
    ``-1 + 2j.pi.f`` pour ceux autour de ``f`` Hz.

    ARPACK impose ``k < n - 1`` : il lui faut un sous-espace de Krylov plus
    grand que le nombre de modes demandes. Pour un spectre complet, c'est
    ``numpy.linalg.eig`` qu'il faut, pas ce shift-and-invert.

    Rend ``(lam, Vr)`` trie par distance croissante a ``sigma``.
    """
    A = _A(ss)
    n = A.shape[0]
    if not 0 < k < n - 1:
        raise ValueError(f"k = {k} hors de ]0, n-1[ avec n = {n} (limite ARPACK)")

    # A est REELLE, et c'est un piege : devant une matrice reelle, ARPACK
    # emploie son shift-invert reel, qui factorise Re[(A - sigma.I)^-1]. Cet
    # operateur a les memes vecteurs propres, mais pas le meme classement :
    # il rend les paires conjuguees les plus resonnantes, pas les k valeurs
    # propres les plus proches de sigma. Sur SMIB, a sigma = -1 + 2j.pi.1,5,
    # il rendait des modes a 15 et 27 de distance au lieu de 9,4 -- et c'est
    # pourtant le cas d'usage que cette fonction recommande. Passer A en
    # complexe force le shift-invert complexe, qui classe correctement.
    M = sp.csc_matrix(A).astype(complex)

    # sigma exactement sur une valeur propre rend (A - sigma.I) singuliere.
    # Ce n'est pas un cas tordu ici : un reseau sans bus infini porte un mode
    # de jauge EXACTEMENT nul, et sigma = 0 est justement ce qu'on prend pour
    # chercher les modes lents. Un decalage infime suffit a conditionner la
    # factorisation, et il ne change pas les valeurs propres rendues : le
    # shift-invert converge vers celles de A, sigma ne fixe que l'endroit ou
    # il regarde.
    echelle = max(float(np.abs(A).max()), 1.0)
    for decalage in (0.0, 1e-9 * echelle, 1e-6 * echelle):
        try:
            lam, Vr = spla.eigs(M, k=k, sigma=sigma + decalage)
            break
        except RuntimeError as err:  # "Factor is exactly singular"
            if "singular" not in str(err).lower():
                raise
            derniere = err
    else:
        raise RuntimeError(
            f"shift-and-invert impossible autour de sigma = {sigma} : "
            f"{derniere}. Choisir un sigma qui ne soit pas exactement une "
            "valeur propre."
        )

    ordre = np.argsort(np.abs(lam - sigma))
    return lam[ordre], Vr[:, ordre]
