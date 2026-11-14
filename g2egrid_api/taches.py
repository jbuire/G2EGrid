"""Les calculs longs de l'atelier, en tache de fond.

Une requete HTTP qui dure deux minutes fige l'interface, ne dit rien de son
avancement et ne s'arrete pas quand on change d'avis. Les calculs longs
(robustesse, continuation, boucle mineure, rapport) passent donc par une
TACHE : la requete de lancement rend tout de suite un identifiant ; le
navigateur suit l'avancement par un flux SSE (``text/event-stream``), et peut
demander l'arret.

L'avancement vient du moteur lui-meme (:mod:`g2egrid.noyau.progres`) : chaque tache
installe un rappel qui met a jour son etat et leve
:class:`~g2egrid.noyau.progres.Interrompu` si l'arret a ete demande. Le calcul
s'arrete ainsi au prochain point de controle de sa boucle, proprement.

Les taches tournent dans des fils (le calcul numerique libere le GIL dans
LAPACK) ; deux au plus a la fois, les suivantes attendent -- mais pas sans
limite : au-dela de :data:`MAX_ACTIVES` taches non terminees (en cours ou en
attente), :func:`lancer` refuse avec :class:`TropDeTaches` (l'API repond
429), au lieu d'empiler des calculs de plusieurs minutes dont chacun garde
les tables d'un reseau. Les vingt dernieres terminees sont gardees pour
qu'on puisse relire leur resultat.
"""

from __future__ import annotations

import contextvars
import logging
import os
import threading
import time
import uuid
from collections import OrderedDict
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any

from g2egrid.noyau.progres import Interrompu, suivi

__all__ = [
    "MAX_ACTIVES",
    "MAX_GARDEES",
    "Tache",
    "TropDeTaches",
    "annuler",
    "lancer",
    "lire",
]

MAX_GARDEES = 20
#: taches non terminees admises (deux tournent, les autres attendent) ;
#: ``G2EGRID_TACHES_MAX`` le change
MAX_ACTIVES = 8
journal = logging.getLogger("g2egrid_api")


class TropDeTaches(RuntimeError):
    """Trop de calculs longs deja lances : on refuse le suivant."""


def _max_actives() -> int:
    try:
        v = int(os.environ.get("G2EGRID_TACHES_MAX", ""))
    except ValueError:
        return MAX_ACTIVES
    return v if v > 0 else MAX_ACTIVES


_POOL = ThreadPoolExecutor(max_workers=2, thread_name_prefix="g2egrid-tache")
_TACHES: OrderedDict[str, "Tache"] = OrderedDict()
_VERROU = threading.Lock()


@dataclass
class Tache:
    id: str
    genre: str
    etat: str = "en_attente"  # en_attente, en_cours, fini, echec, annule
    fraction: float = 0.0
    message: str = ""
    resultat: Any = None
    erreur: str = ""
    debut: float = field(default_factory=time.time)
    fin: float | None = None
    arret: threading.Event = field(default_factory=threading.Event)
    #: incremente a chaque changement : le flux SSE n'envoie que du neuf
    revision: int = 0

    def resume(self, avec_resultat: bool = False) -> dict:
        d = {
            "id": self.id,
            "genre": self.genre,
            "etat": self.etat,
            "fraction": self.fraction,
            "message": self.message,
            "erreur": self.erreur,
            "duree": (self.fin or time.time()) - self.debut,
            "revision": self.revision,
        }
        if avec_resultat and self.etat == "fini":
            d["resultat"] = self.resultat
        return d

    @property
    def terminee(self) -> bool:
        return self.etat in ("fini", "echec", "annule")


def _toucher(t: Tache, **champs) -> None:
    with _VERROU:
        for k, v in champs.items():
            setattr(t, k, v)
        t.revision += 1


def lancer(genre: str, calcul: Callable[[], Any]) -> Tache:
    """Lance ``calcul()`` en tache de fond ; rend la tache tout de suite."""
    t = Tache(id=uuid.uuid4().hex[:12], genre=genre)
    with _VERROU:
        actives = sum(1 for x in _TACHES.values() if not x.terminee)
        if actives >= _max_actives():
            raise TropDeTaches(
                f"trop de calculs en cours ({actives} taches non terminees, "
                f"au plus {_max_actives()}) : attendez la fin d'un calcul ou "
                "arretez-en un (bouton Arreter), puis relancez."
            )
        _TACHES[t.id] = t
        while len(_TACHES) > MAX_GARDEES:
            vieille = next(iter(_TACHES.values()))
            if not vieille.terminee:
                break
            _TACHES.popitem(last=False)

    def rappel(fraction: float, message: str) -> None:
        if t.arret.is_set():
            raise Interrompu("arret demande")
        _toucher(t, fraction=fraction, message=message)

    def executer() -> None:
        if t.arret.is_set():
            _toucher(t, etat="annule", fin=time.time())
            return
        _toucher(t, etat="en_cours")
        try:
            with suivi(rappel):
                r = calcul()
        except Interrompu:
            _toucher(t, etat="annule", fin=time.time(), message="interrompu")
        except Exception as err:  # noqa: BLE001 -- rendu a l'utilisateur
            # une HTTPException des fonctions de route porte deja son message
            # (et, pour une erreur interne, sa reference : la trace est au
            # journal, jamais dans la reponse)
            detail = getattr(err, "detail", None) or str(err) or type(err).__name__
            if getattr(err, "status_code", None) is None and not isinstance(
                err, (ValueError, KeyError, RuntimeError)
            ):
                ident = uuid.uuid4().hex[:8]
                journal.error("erreur interne de tache %s", ident, exc_info=err)
                detail = (
                    f"erreur interne du serveur (reference {ident}) : la trace "
                    "est dans le journal du serveur."
                )
            _toucher(t, etat="echec", fin=time.time(), erreur=detail)
        else:
            _toucher(t, etat="fini", fin=time.time(), fraction=1.0, resultat=r)

    _POOL.submit(contextvars.copy_context().run, executer)
    return t


def lire(id_: str) -> Tache | None:
    with _VERROU:
        return _TACHES.get(id_)


def annuler(id_: str) -> Tache | None:
    t = lire(id_)
    if t is not None and not t.terminee:
        t.arret.set()
        _toucher(t, message="arret demande...")
    return t


def liste() -> list[dict]:
    with _VERROU:
        return [t.resume() for t in reversed(_TACHES.values())]
