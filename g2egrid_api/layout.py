"""Placement automatique des bus, pour le schema du constructeur.

Un reseau n'a pas de coordonnees : il n'a qu'une topologie. Il faut donc en
inventer, et qu'elles soient LISIBLES -- les bus relies proches, les autres
ecartes -- et surtout **stables** : deux appels sur le meme reseau doivent
rendre le meme dessin, sinon le schema saute a chaque rechargement.

D'ou un Fruchterman-Reingold maison, en numpy, avec un tirage initial a
graine fixe : une trentaine de lignes, contre une dependance de plus
(networkx) pour le meme resultat.
"""

from __future__ import annotations

import numpy as np

__all__ = ["placer"]


def placer(
    bus: list[int],
    aretes: list[tuple[int, int]],
    iterations: int = 300,
    graine: int = 0,
) -> dict[int, tuple[float, float]]:
    """Place les bus dans le carre [0, 1] x [0, 1].

    Modele a ressorts : les aretes attirent en ``d^2 / k``, tous les couples
    se repoussent en ``k^2 / d``, et le deplacement admissible decroit a
    chaque iteration (recuit) pour que le dessin se fige au lieu d'osciller.
    """
    n = len(bus)
    if n == 0:
        return {}
    if n == 1:
        return {bus[0]: (0.5, 0.5)}

    idx = {b: i for i, b in enumerate(bus)}
    liens = np.array(
        [(idx[a], idx[b]) for a, b in aretes if a in idx and b in idx and a != b],
        dtype=int,
    ).reshape(-1, 2)

    rng = np.random.default_rng(graine)
    pos = rng.random((n, 2)) - 0.5
    k = 1.0 / np.sqrt(n)  # distance d'equilibre
    pas = 0.1

    for _ in range(iterations):
        delta = pos[:, None, :] - pos[None, :, :]
        d = np.maximum(np.linalg.norm(delta, axis=2), 1e-3)
        # repulsion de tous contre tous
        force = (delta / d[:, :, None]) * (k * k / d)[:, :, None]
        np.einsum("iij->ij", force)[...] = 0.0
        f = force.sum(axis=1)
        # attraction le long des aretes
        if len(liens):
            dl = pos[liens[:, 0]] - pos[liens[:, 1]]
            dd = np.maximum(np.linalg.norm(dl, axis=1), 1e-3)
            fa = (dl / dd[:, None]) * (dd * dd / k)[:, None]
            np.add.at(f, liens[:, 0], -fa)
            np.add.at(f, liens[:, 1], +fa)
        norme = np.maximum(np.linalg.norm(f, axis=1), 1e-9)
        pos += (f / norme[:, None]) * np.minimum(norme, pas)[:, None]
        pas *= 0.99

    # recadrage dans [0, 1], en gardant les proportions
    mini, maxi = pos.min(axis=0), pos.max(axis=0)
    etendue = np.maximum(maxi - mini, 1e-6).max()
    pos = (pos - mini) / etendue
    pos += (1.0 - pos.max(axis=0)) / 2.0  # centre le dessin
    return {b: (float(pos[i, 0]), float(pos[i, 1])) for b, i in idx.items()}
