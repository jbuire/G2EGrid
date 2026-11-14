"""L'avancement d'un calcul long, et son interruption.

Une etude de robustesse refait load flow et spectre pour chaque variante ;
une continuation, pour chaque pas. Depuis l'interface, on veut voir ou elle
en est et pouvoir l'arreter. Mais le moteur ne doit rien savoir de HTTP ni
des fils d'execution : il se contente d'appeler :func:`signaler` dans ses
boucles. Qui veut suivre installe un rappel par :func:`suivi` ; sans rappel,
:func:`signaler` ne coute qu'une lecture de variable.

Le rappel est porte par une ``contextvars.ContextVar`` : chaque fil (chaque
tache de l'atelier) a le sien, sans etat global partage. Pour interrompre, le
rappel leve :class:`Interrompu` -- l'exception remonte de la boucle comme
n'importe quelle erreur, et le calcul s'arrete au prochain point de controle.
"""

from __future__ import annotations

import contextlib
import contextvars
from collections.abc import Callable, Iterator

__all__ = ["Interrompu", "signaler", "suivi"]

_RAPPEL: contextvars.ContextVar[Callable[[float, str], None] | None] = (
    contextvars.ContextVar("g2egrid_progres", default=None)
)


class Interrompu(Exception):
    """Le calcul a ete interrompu a la demande."""


def signaler(fraction: float, message: str = "") -> None:
    """Le calcul en est a ``fraction`` (0 a 1) ; ``message`` dit ou."""
    rappel = _RAPPEL.get()
    if rappel is not None:
        rappel(max(0.0, min(1.0, float(fraction))), message)


@contextlib.contextmanager
def suivi(rappel: Callable[[float, str], None]) -> Iterator[None]:
    """Installe ``rappel(fraction, message)`` pour la duree du bloc."""
    jeton = _RAPPEL.set(rappel)
    try:
        yield
    finally:
        _RAPPEL.reset(jeton)
