"""Bibliotheque d'equipements (portage fidele des classes MATLAB)."""

from g2egrid.elements.branch import Branch, Line, LineDC, Transformer
from g2egrid.elements.bus import Bus, BusDC, IdealBus, ResistiveBus
from g2egrid.elements.datacenter import (
    BatterieDCDC,
    ChargeACPref,
    OnduleurVF,
    RedresseurDiode,
)
from g2egrid.elements.gfl import GFL1, GFL2, GFL2_DEFAULT, GFL2DC, GFL_DEFAULT
from g2egrid.elements.gfm import (
    GFM1,
    GFM2,
    GFM3,
    GFM3DC,
    GFM4,
    GFM4DC,
    GFM_DC_DEFAULT,
    GFM_DEFAULT,
    LOIS_GFM,
    GFMdVOC,
    GFMMatching,
)
from g2egrid.elements.groupe_electrogene import GROUPE_DEFAULT, GroupeElectrogene
from g2egrid.elements.load import (
    ChargeDC_P,
    ChargeDC_R,
    Defaut,
    DefautDC,
    ILoad,
    PLoad,
    ReactanceShunt,
    RLLoad,
    ZLoad,
)
from g2egrid.elements.mmc import MMC, MMC_DEFAULT, MMCEnergie
from g2egrid.elements.repere import Repere
from g2egrid.elements.source import (
    IB,
    IB_DEFAULT,
    ConvertisseurDCDC,
    SourceDC,
    VoltageSource,
)
from g2egrid.elements.sync_machine import (
    KUNDUR,
    SM_DEFAULT,
    SynchronousMachine,
    SynchronousMachineTR,
)
from g2egrid.elements.turbo_alternateur import (
    ARBRE_AM600,
    ARBRE_FBM,
    ARBRE_NUCLEAIRE_EPRI,
    LIGNES_ARBRE,
    TAU_TORSION,
    LigneArbre,
    classe_multi_masses,
    ligne_a_l_echelle,
    ligne_depuis_physique,
    ligne_deux_masses,
    turbo_alternateur,
    turbo_alternateur_nomme,
)

# Alias historiques du portage Python.
GFM = GFM1
GFL = GFL1
SM = SynchronousMachine
GenericLoad = RLLoad

__all__ = [
    "ARBRE_AM600",
    "ARBRE_FBM",
    "ARBRE_NUCLEAIRE_EPRI",
    "LIGNES_ARBRE",
    "TAU_TORSION",
    "LigneArbre",
    "classe_multi_masses",
    "ligne_a_l_echelle",
    "ligne_depuis_physique",
    "ligne_deux_masses",
    "turbo_alternateur",
    "turbo_alternateur_nomme",
    "BatterieDCDC",
    "ChargeACPref",
    "GROUPE_DEFAULT",
    "GroupeElectrogene",
    "OnduleurVF",
    "RedresseurDiode",
    "Branch",
    "Line",
    "Transformer",
    "Bus",
    "BusDC",
    "ResistiveBus",
    "IdealBus",
    "LineDC",
    "SourceDC",
    "ConvertisseurDCDC",
    "ChargeDC_R",
    "ChargeDC_P",
    "GFM3DC",
    "GFM4DC",
    "GFMdVOC",
    "GFMMatching",
    "LOIS_GFM",
    "GFL2DC",
    "MMC",
    "MMC_DEFAULT",
    "MMCEnergie",
    "RLLoad",
    "ZLoad",
    "ReactanceShunt",
    "Defaut",
    "DefautDC",
    "Repere",
    "ILoad",
    "PLoad",
    "GenericLoad",
    "IB",
    "VoltageSource",
    "IB_DEFAULT",
    "SynchronousMachine",
    "SynchronousMachineTR",
    "SM",
    "SM_DEFAULT",
    "KUNDUR",
    "GFM1",
    "GFM2",
    "GFM3",
    "GFM4",
    "GFM",
    "GFM_DEFAULT",
    "GFM_DC_DEFAULT",
    "GFL1",
    "GFL2",
    "GFL",
    "GFL_DEFAULT",
    "GFL2_DEFAULT",
]
