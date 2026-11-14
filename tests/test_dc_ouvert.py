"""Les regles du GFM/GFL avec sa source de l'editeur (outils ``gfm_batterie``,
``gfm_pv``, ``gfm_dc``, ``gfl_*`` dans ``web/app.js``, fonctions
``poserConvertisseurSource`` et ``accorderIleDC``).

L'editeur pose un convertisseur avec son bus continu et un condensateur de
0,2 pu, puis accorde sa commande a ce qu'on branche cote continu :

* batterie -> le convertisseur regle sa puissance (GFM : PV ; GFL : PQ), la
  batterie tient la tension continue ;
* PV ou rien -> le convertisseur tient la tension continue (VdcV / VdcQ).

Ces tests verifient que ces combinaisons ont un point de fonctionnement et
sont stables, et que les autres -- un PV sans convertisseur a Vdc, ou le bus
continu sans condensateur -- sont bien celles que l'editeur evite.
"""

import numpy as np
import pytest

from g2egrid.robustesse import _spectre
from g2egrid.tables import Tables, assembler


def reseau(conv, mode, dcdc=None, pinj=10e6, c_pu=0.2):
    t = Tables(f0=50, Sb=100e6, nom="ouvert")
    t.Y_DER = [
        (1, 1, 20e3, 1.0, 0, 0, 0, 0, 0),
        (2, 2 if conv == "GFM" else 3, 20e3, 1.0, 0, pinj, 0, 0, 0),
        (3, 3, 20e3, 1.0, 0, 0, 0, 20e6, 5e6),
    ]
    t.T_IBs = [(1, 1, 100e6, 1.0, 0, 0.001, 0.05)]
    t.T_Lines_pu = [
        (1, 1, 3, 1, 0.01, 0.05, 0.0, 20e3, 100e6),
        (2, 3, 2, 1, 0.01, 0.05, 0.0, 20e3, 100e6),
    ]
    if conv == "GFM":
        t.T_GFM3s = [(1, 2, 100e6, 40e3)]
    else:
        t.T_GFL2s = [(1, 2, 100e6, 40e3)]
    t.T_BusDC = [(10, 1, 40e3, 1.0, 0.0)]
    if c_pu:
        t.T_ShuntsDC = [(1, 10, c_pu)]
    t.T_VSCs = [(1, 2, 10, mode)]
    if dcdc:
        t.T_DCDC = [(1, 10, 50e6, dcdc, 10e6)]
    return t


def re_max(t):
    return float(_spectre(t).lam.real.max())


# la commande que l'editeur choisit pour chaque branchement
ACCORD = [
    ("GFM", "VdcV", None),
    ("GFM", "VdcV", "PV"),
    ("GFM", "PV", "batterie"),
    ("GFL", "VdcQ", None),
    ("GFL", "VdcQ", "PV"),
    ("GFL", "PQ", "batterie"),
]


@pytest.mark.parametrize("conv,mode,dcdc", ACCORD)
def test_commande_accordee_est_stable(conv, mode, dcdc):
    for pinj in (10e6, -10e6):
        assert re_max(reseau(conv, mode, dcdc, pinj)) < 1e-3


@pytest.mark.parametrize("conv,mode", [("GFM", "PV"), ("GFL", "PQ")])
def test_pv_sans_convertisseur_a_vdc_n_a_pas_de_reference(conv, mode):
    with pytest.raises(Exception, match="reference de tension"):
        assembler(reseau(conv, mode, "PV")).run_loadflow()


def test_sans_condensateur_le_gfm_a_batterie_oscille():
    """Le mode bus continu / courant du DC/DC : instable sans condensateur,
    d'ou celui que l'editeur pose d'office."""
    sp = _spectre(reseau("GFM", "PV", "batterie", c_pu=0.0))
    k = np.argmax(sp.lam.real)
    assert sp.lam[k].real > 0
    porteurs = {sp.etats[j] for j in np.argsort(-sp.p[k])[:2]}
    assert porteurs == {"bus10.v", "DCDC1_10.i"}
