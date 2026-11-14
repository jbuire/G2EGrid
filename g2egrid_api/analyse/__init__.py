"""Le coeur de l'interface : des tables a des resultats.

Une seule regle d'architecture : **rien ici ne sait ce qu'est HTTP**. Ce
module prend des tables de reseau (le meme :class:`g2egrid.reseau.tables.Tables` que
les notebooks) et rend des structures python ; ``app.py`` ne fait que les
convertir en JSON. C'est ce qui permet de tester ces fonctions comme le
reste du paquet, et de les appeler depuis un script sans lancer de serveur.

Le reseau construit est **mis en cache**, indexe par une empreinte de ses
tables : assembler wscc3sm_cigre3 et resoudre son equilibre demande une
dizaine de secondes, et l'interface redemande le meme reseau a chaque
onglet.


Le module est un PAQUET depuis le 2026-10-28 : une section de l'ancien ``analyse.py`` par sous-module. Ce fichier reexporte tous les noms d'avant, prives compris.
"""

# ruff: noqa: F401  (reexportation de tout l ancien module)
from __future__ import annotations

import hashlib
import json
import pathlib
import re
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment

import g2egrid.reseau.catalogue
from g2egrid.analyse import identification, modal
from g2egrid.analyse.comparaison import comparer as comparer_versions
from g2egrid.analyse.sensibilites import (
    sensibilites as sensibilites_parametriques,
)
from g2egrid.formats import notebook
from g2egrid.formats.export import exporter, vers_matlab
from g2egrid.multi_niveau.impact import impact_distribution
from g2egrid.multi_niveau.separation import (
    reponse_indicielle as reponse_lineaire,
)
from g2egrid.noyau.loadflow import PQ
from g2egrid.reseau import construction
from g2egrid.reseau import disjoncteurs as disj
from g2egrid.reseau.catalogue import RESEAUX, fiche
from g2egrid.reseau.flux import flux_branches
from g2egrid.reseau.network import Network
from g2egrid.reseau.tables import (
    CDC_DEFAUT,
    CHARGES_SPECIALES,
    CONFIGS_DC,
    MODELES_STATION,
    TABLES,
    Tables,
    assembler,
    tables_depuis_json,
    tables_vers_json,
)
from g2egrid.reseau.zones import (
    ECHELLE,
    FAMILLES,
    FORMAT,
    PREFIXE_MIN,
    SLACK,
    VERSION,
    Instance,
    Reseau,
    TransfoPoste,
    Zone,
    aplatir,
    reseau_depuis_json,
    reseau_vers_json,
)
from g2egrid.scenarios import simplification as simp
from g2egrid.scenarios.robustesse import etudier as etudier_robustesse
from g2egrid.temporel import (
    GENRES_MANOEUVRE,
    Evenement,
    Manoeuvre,
    ModeleNonLineaire,
    Trajectoire,
)
from g2egrid.temporel import mesures as mesures_physiques

from ..layout import placer
from ._commun import (
    _CACHE,
    _CACHE_MAX,
    EQUIPEMENTS,
    GROUPES_DYNAMIQUES,
    JAUGE,
    ZETA_OSCILLANT,
    Assemble,
    _empreinte,
)
from .apres_evenement import (
    _nombre,
    _regime_lf,
    _relais_ac,
    _valeur_finale,
    apres_evenement,
    cascade,
)
from .assemblage import (
    _BUS_ALTERNATIF,
    _BUS_CONTINU,
    _COLONNES_ENTIERES,
    _SUFFIXES_NOMBRES,
    _est_entier,
    _est_nombre,
    construire,
    controler_entree,
    lire_tables,
    references_inconnues,
    vider_cache,
)
from .balayages import (
    CIBLES,
    _appliquer,
    _participations_par_element,
    _suivre,
    balayage,
)
from .blocs_edition import (
    ModeleAbsent,
    bibliotheque_enregistrer,
    bibliotheque_lire,
    bibliotheque_liste,
    bibliotheque_supprimer,
    blocs_aplatir,
    blocs_degrouper,
    blocs_depuis_element,
    blocs_regrouper,
)
from .equations_edition import (
    blocs_remplacer,
    equations_ajouter,
    equations_depuis_schema,
    equations_element,
    equations_liste_reseau,
    equations_modifier,
    equations_parametre,
    equations_revenir,
    equations_supprimer,
    equations_verifier,
    equations_vers_blocs,
)
from .catalogue_api import (
    _bus_d_instance,
    _iles,
    _noeuds,
    _positions,
    _positions_connues,
    _reponse_reseau,
    catalogue,
    depuis_catalogue,
)
from .controles import (
    t_eig,
    un_slack_tout_pq,
    validation,
)
from .exports import (
    EXPORTS,
    TABLES_BRANCHES,
    _export_notebook,
    _intro_notebook,
    _nom_de_sortie,
    _nom_matlab,
    exporter_vers,
    importer_externe,
)
from .fichiers_reseaux import (
    _CHAMPS_INSTANCE,
    _NOM_FICHIER,
    _REGLAGES_JSON,
    _fiche_fichier,
    _nom_de_fichier,
    _reseau_depuis,
    _vers_fichier,
    convertir_fichier,
    enregistrer_fichier,
    fichiers,
    lire_fichier,
)
from .generation import (
    _tables_de,
    code_python,
    comparer,
)
from .impedance_continuation import (
    _bifurcations,
    _limites_suite,
    _trajectoires,
    balayage_statismes,
    boucle_mineure,
    continuation_parametre,
    elements_reseau,
)
from .reglages_defaut import (
    _defauts_zones,
    _modes_vsc,
    defauts,
)
from .reparations import (
    _FORMEURS,
    _LIBELLES,
    _TABLES_LIMITEES,
    _avis_convergence,
    _avis_limites,
    _formeur_principal,
    _mettre_type,
    _plus_proche,
    _proposer,
    reparer,
    valider,
)
from .resultats import (
    _PREFERES,
    _domaines_modes,
    _etats_par_defaut,
    _evenements,
    _mesures_equilibre,
    _nom_element,
    _projection,
    _sans_nan,
    _sensibilites_parametres,
    blocs_compiler,
    blocs_types,
    court_circuit,
    detail_mode,
    disjoncteurs,
    equations_modele,
    fiche_element,
    impact_continu,
    loadflow,
    machine_standard,
    manoeuvrables,
    modes,
    noms_mesures,
    reglage_boucles,
    signaux,
    temporel,
    topologie,
)
from .transport_distribution import (
    _comparaison_json,
    _mode_compare_json,
    aplatir_reseau,
    definir_zone,
    impact,
    multi_ajouter_gfl,
    multi_comparer,
    multi_penetration,
    robustesse,
    simplification_comparer,
    simplification_inventaire,
)
from .zones_api import (
    _source_import,
    importer,
)

__all__ = [
    "EXPORTS",
    "TABLES",
    "Assemble",
    "aplatir_reseau",
    "construire",
    "defauts",
    "depuis_catalogue",
    "enregistrer_fichier",
    "exporter",
    "exporter_vers",
    "fichiers",
    "loadflow",
    "modes",
    "signaux",
    "tables_depuis_json",
    "temporel",
    "topologie",
    "valider",
    "vider_cache",
]
