"""Les valeurs par defaut proposees par l'interface.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import numpy as np
from g2egrid.reseau.tables import CDC_DEFAUT, CONFIGS_DC, MODELES_STATION, TABLES
from g2egrid.reseau.zones import (
    ECHELLE,
    FAMILLES,
    FORMAT,
    PREFIXE_MIN,
    SLACK,
    VERSION,
    Instance,
    TransfoPoste,
)
from ._commun import EQUIPEMENTS, GROUPES_DYNAMIQUES


# ------------------------------------------------------------------ defauts
def defauts() -> dict:
    """Les colonnes des tables et les reglages par defaut de chaque famille.

    C'est ce qui permet au constructeur de proposer un equipement deja
    parametre : les valeurs viennent des classes elles-memes, pas d'une
    copie dans l'interface qui deriverait.
    """
    from g2egrid.elements.gfl import GFL2_DEFAULT, GFL_DEFAULT
    from g2egrid.elements.mmc import MMC_DEFAULT
    from g2egrid.elements.gfm import GFM_DC_DEFAULT, GFM_DEFAULT
    from g2egrid.elements.limiteurs import REGLAGES_BUS_CONTINU
    from g2egrid.elements.sync_machine import SM_DEFAULT

    def propre(d: dict) -> dict:
        # les drapeaux ``mode_*`` des limiteurs sont un etat de simulation,
        # pas un reglage : on ne les propose pas
        return {
            k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
            for k, v in d.items()
            if not k.startswith("mode_")
        }

    return {
        "tables": TABLES,
        "configs_dc": list(CONFIGS_DC),
        "equipements": {
            k: {"table": v["table"], "titre": v["titre"], "couleur": v["couleur"]}
            for k, v in EQUIPEMENTS.items()
        },
        "reglages": {
            "SM": propre(SM_DEFAULT),
            "GFM": propre({**GFM_DEFAULT, **GFM_DC_DEFAULT}),
            "GFL": propre({**GFL_DEFAULT, **GFL2_DEFAULT}),
            "MMC": propre(MMC_DEFAULT),
        },
        # ce qui change pour un GFM dont la tension continue est une VARIABLE
        # (GFM3/GFM4, ports continus) : la saturation de tension y est active
        "reglages_bus_continu": {
            "tables": sorted(t for t in CDC_DEFAUT if t.startswith("T_GFM")),
            "valeurs": propre(REGLAGES_BUS_CONTINU),
        },
        "modeles_station": list(MODELES_STATION),
        "types_bus": {"1": "slack", "2": "PV", "3": "PQ"},
        "dynamiques": GROUPES_DYNAMIQUES,
        # Les modes admissibles d'un VSC viennent des CLASSES, par leur
        # attribut modes_dc, et non d'une liste recopiee ici. Un GFM tient
        # la tension alternative, un GFL le reactif : ils n'ont pas les
        # memes couples possibles, et l'assembleur refuse les autres. Les
        # proposer depuis la source evite de laisser choisir un mode que le
        # modele rejettera.
        "modes_vsc": _modes_vsc(),
        # Ce qu'il faut pour EDITER des zones, et non seulement les afficher.
        "zones": _defauts_zones(),
    }


def _modes_vsc() -> dict[str, list[str]]:
    from g2egrid.elements.gfl import GFL2DC
    from g2egrid.elements.gfm import GFM3DC, GFM4DC

    return {
        "GFM3": list(getattr(GFM3DC, "modes_dc", ())),
        "GFM4": list(getattr(GFM4DC, "modes_dc", ())),
        "GFL2": list(getattr(GFL2DC, "modes_dc", ())),
    }


def _defauts_zones() -> dict:
    """Les regles que le constructeur doit connaitre pour editer des zones.

    Tout vient de :mod:`g2egrid.reseau.zones`, y compris l'echelle ``K``. La recopier
    dans l'interface serait la garantie qu'elle derive, et une echelle fausse
    NE SE VOIT PAS : elle ne fait pas echouer l'assemblage, elle rend un reseau
    plausible dont les puissances ne sont pas celles qu'on a saisies. Elle est
    rendue par NOM de colonne et non par indice, parce que c'est le nom que
    l'interface affiche et qu'un indice se decale quand une table gagne une
    colonne.
    """
    neuve = Instance(zone="", prefixe=PREFIXE_MIN, bus_amont=0)
    return {
        "format": {"format": FORMAT, "version": VERSION},
        "prefixe_min": PREFIXE_MIN,
        #: le ``Y_DER.Type`` qu'une zone ne peut pas porter : elle pend au
        #: transport et ne tient pas sa propre tension.
        "type_slack": SLACK,
        "instance": {
            "zone": neuve.zone,
            "prefixe": neuve.prefixe,
            "bus_amont": neuve.bus_amont,
            "K": neuve.K,
            "transfo": None,
            "charge_residuelle": neuve.charge_residuelle,
            "shunt_pcc": neuve.shunt_pcc,
        },
        "transfo_poste": list(TransfoPoste._fields),
        "echelle": {
            table: {
                TABLES[table][c]: op
                for c, op in colonnes.items()
                if c < len(TABLES[table])
            }
            for table, colonnes in ECHELLE.items()
        },
        # Quelle table de DER releve de quel dictionnaire de reglages : c'est
        # ce qui dit a l'interface que le GFM_PARAMS d'une zone ne vaut que
        # pour les GFM de cette zone, et sera replie par bus a l'aplatissement.
        "familles": dict(FAMILLES),
    }
