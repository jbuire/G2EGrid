"""Reponse frequentielle d'un equipement vue de ses bornes.

Chaque equipement est un systeme dq a deux entrees de port (``vd_t``, ``vq_t``)
et deux sorties de port (``id_t``, ``iq_t``, courant INJECTE au noeud). En
figeant ses autres entrees (consignes, ``w_ref``, ``theta_g``), sa linearisation
donne une matrice d'admittance 2x2 ::

    [di_d]          [dv_d]
    [di_q] = -Y(s) . [dv_q]        Y(s) = -(C (sI - A)^-1 B + D)

Le signe fait de ``Y`` l'admittance d'ENTREE de l'equipement (courant compte
positif vers l'equipement), pour retrouver les conventions usuelles : une source
de tension parfaite a une impedance nulle, une source de courant parfaite une
impedance infinie. C'est ce que montre le trace de ``Z = Y^-1`` :

* une **source de tension** (machine synchrone, grid-forming) tient sa tension :
  un creux de tension lui fait debiter beaucoup de courant, donc ``|Z|`` est
  faible en basse frequence ;
* une **source de courant** (grid-following) tient son courant : une variation
  de tension ne le change presque pas, donc ``|Z|`` est eleve en basse
  frequence.

Une matrice dq n'est pas un scalaire : ses deux valeurs singulieres donnent la
raideur de l'equipement dans sa direction la plus raide et dans la plus molle.
Toutes les sources ont une direction molle en basse frequence, celle de l'angle,
que rien ne fixe quand on regarde un equipement seul ; c'est dans la direction du
module de tension que se lit la nature de la source.
"""

from __future__ import annotations

import numpy as np

from g2egrid.noyau.element import Element

__all__ = [
    "reponse",
    "admittance_dq",
    "impedance_dq",
    "axes_tension",
    "gain",
    "valeurs_singulieres",
    "vecteur_complexe",
    "module_db",
]


def reponse(
    el: Element,
    f: np.ndarray,
    entrees: tuple[str, ...],
    sorties: tuple[str, ...],
) -> np.ndarray:
    """Matrice de transfert ``sorties``/``entrees`` de l'element, a ``2.pi.f``.

    ``el`` doit avoir ete initialise (``Network.solve_equilibrium``). Rend un
    tableau complexe de forme ``(len(f), len(sorties), len(entrees))``.
    """
    if el.x0 is None:
        raise ValueError(f"{el.name}: equilibre non calcule (solve_equilibrium)")
    cm = el.compiled()
    A, B, C, D = cm.linearize(el.x0, el.z0, el.u0, el.p)
    iu = [cm.iu(n) for n in entrees]
    iy = [cm.iy(n) for n in sorties]
    B, C, D = B[:, iu], C[iy, :], D[np.ix_(iy, iu)]
    f = np.atleast_1d(np.asarray(f, dtype=float))
    n = A.shape[0]
    G = np.empty((f.size, len(sorties), len(entrees)), dtype=complex)
    for k, fk in enumerate(f):
        G[k] = C @ np.linalg.solve(2j * np.pi * fk * np.eye(n) - A, B) + D
    return G


def admittance_dq(el: Element, f: np.ndarray, port: int = 0) -> np.ndarray:
    """Admittance d'entree dq de l'element.

    ``Y[k]`` est la matrice 2x2 a la frequence ``f[k]``, en pu base systeme.

    Le shunt de ``shunt_admittance`` n'est pas ajoute : il ne sert qu'au Y-bus
    du load flow. La sortie de courant du modele le porte deja -- la machine
    synchrone injecte ``ig = i - v/Rg``, resistance de causalite comprise.
    """
    vd, vq = (el.input_names()[i] for i in el.port_input_idx(port))
    id_, iq = (el.output_names()[i] for i in el.port_output_idx(port))
    return -reponse(el, f, (vd, vq), (id_, iq))


def impedance_dq(el: Element, f: np.ndarray, port: int = 0) -> np.ndarray:
    """Impedance d'entree dq de l'element : inverse de :func:`admittance_dq`."""
    return np.linalg.inv(admittance_dq(el, f, port))


def axes_tension(el: Element, port: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Les deux directions dq du point de fonctionnement : module et angle.

    ``e_mod`` est le vecteur unitaire porte par la tension de bornes : une
    perturbation dans cette direction est une variation de MODULE. ``e_ang`` lui
    est orthogonal : une perturbation dans cette direction est, au premier ordre,
    une variation d'ANGLE de la tension.
    """
    if el.u0 is None:
        raise ValueError(f"{el.name}: equilibre non calcule (solve_equilibrium)")
    iv, iq = el.port_input_idx(port)
    v = np.array([el.u0[iv], el.u0[iq]])
    e_mod = v / np.linalg.norm(v)
    return e_mod, np.array([-e_mod[1], e_mod[0]])


def gain(M: np.ndarray, e: np.ndarray) -> np.ndarray:
    """Gain ``||M e||`` : module de la reponse a une excitation unitaire selon ``e``.

    C'est la norme du vecteur ``M e``, et non la projection ``e' M e`` : sous une
    impedance inductive, une variation de MODULE de tension appelle un courant en
    QUADRATURE, que la projection annulerait.
    """
    return np.linalg.norm(np.einsum("kij,j->ki", M, e), axis=1)


def valeurs_singulieres(M: np.ndarray) -> np.ndarray:
    """Valeurs singulieres de chaque matrice de ``M``, forme ``(len(f), 2)``."""
    return np.linalg.svd(M, compute_uv=False)


def vecteur_complexe(M: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Transferts en vecteur complexe d'une matrice dq 2x2.

    Pour une grandeur complexe ``x = x_d + j.x_q``, une matrice dq se decompose
    en une composante directe et une composante inverse (couplage) ::

        M_direct = ((M_dd + M_qq) + j (M_qd - M_dq)) / 2
        M_inverse = ((M_dd - M_qq) + j (M_qd + M_dq)) / 2

    Un equipement symetrique (``M_dd = M_qq``, ``M_qd = -M_dq``) a une
    composante inverse nulle : sa reponse se lit alors sur le seul ``M_direct``,
    qui est le Bode habituel d'une impedance ou d'une admittance.
    """
    dd, dq, qd, qq = M[..., 0, 0], M[..., 0, 1], M[..., 1, 0], M[..., 1, 1]
    return 0.5 * ((dd + qq) + 1j * (qd - dq)), 0.5 * ((dd - qq) + 1j * (qd + dq))


def module_db(x: np.ndarray) -> np.ndarray:
    """Module en dB, avec un plancher qui evite ``log10(0)``."""
    return 20 * np.log10(np.maximum(np.abs(x), 1e-300))
