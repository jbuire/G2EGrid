"""Les zones : aplatir, definir.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import g2egrid.reseau.catalogue
from g2egrid.reseau.zones import Reseau, aplatir
from .catalogue_api import _positions, _positions_connues, _reponse_reseau
from .fichiers_reseaux import _reseau_depuis


# -------------------------------------------------------------------- zones
def _source_import(source: dict) -> Reseau:
    """Le reseau a importer, avec des positions pour TOUS ses noeuds."""
    if source.get("catalogue"):
        r = g2egrid.reseau.catalogue.charger_reseau(str(source["catalogue"]))
    else:
        r = _reseau_depuis(source)
    if r.instances:
        pos = _positions(aplatir(r))
    else:
        pos = _positions_connues(r.positions, r.transport)
    r.positions = {int(b): tuple(p) for b, p in pos.items()}
    return r


def importer(
    reseau: dict,
    source: dict,
    raccord: dict | None = None,
    garder_slack: bool = False,
    garder_ib: bool = False,
) -> dict:
    """Importe un reseau dans ``reseau`` et l'y raccorde (voir
    :func:`g2egrid.reseau.zones.importer_reseau`) ; rend la forme de LECTURE, plus
    ``import`` : prefixe, correspondance des bus, messages."""
    from g2egrid.reseau.zones import importer_reseau

    r = _reseau_depuis(reseau)
    pos = _positions_connues(r.positions, r.transport)
    r.positions = {int(b): tuple(p) for b, p in pos.items()}
    nouveau, info = importer_reseau(
        r, _source_import(source), None, raccord, garder_slack, garder_ib
    )
    return _reponse_reseau(nouveau) | {
        "import": {
            "prefixe": info["prefixe"],
            "bus": {str(a): b for a, b in info["bus"].items()},
            "messages": info["messages"],
        }
    }
