"""Les reseaux faibles du catalogue (``etudes/reseaux_faibles``) : SCR < 3 au
raccordement de chaque centrale GFL, load flow, plan de tension et stabilite
au point livre ; et le fichier du catalogue est bien celui que reconstruit le
script.

Le SCR de ces reseaux est celui de leur construction : machines en reactance
TRANSITOIRE (``machine="transitoire"``), comme dans l'etude publiee. Le calcul
par defaut de :mod:`g2egrid.court_circuit` est passe en x''d avec les facteurs
de la CEI 60909 (revue du 2026-10-25) ; il releve le SCR de 2 a 13 % (bus 32
d'``ieee39_faible`` : 2,89 -> 3,06). La convention est dite ici plutot que de
changer le reseau ou le critere.
"""

from __future__ import annotations

import importlib.util

import numpy as np
import pytest
from g2egrid.catalogue import RESEAUX, charger
from g2egrid.court_circuit import puissance_court_circuit
from g2egrid.tables import assembler, tables_vers_json

FAIBLES = ("wscc9_faible", "ieee39_faible", "cigre_mv_faible")
from _etudes import dossier_etude

SCRIPT = dossier_etude("reseaux_faibles", sur_le_chemin=False) / "construire.py"


def _script():
    spec = importlib.util.spec_from_file_location("construire_faibles", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.mark.parametrize("nom", FAIBLES)
def test_reseau_faible(nom):
    assert nom in RESEAUX
    t = charger(nom)
    r = puissance_court_circuit(t, machine="transitoire")
    centrales = [n for n in r.noeuds if n.S_gfl > 2e6]
    assert centrales, "pas de centrale GFL"
    assert all(n.SCR < 3 for n in centrales), [(n.bus, n.SCR) for n in centrales]
    net = assembler(t)
    lf = net.run_loadflow()
    assert lf.converged
    assert 0.93 < np.abs(lf.V).min() and np.abs(lf.V).max() < 1.08
    net.solve_equilibrium()
    lam = np.linalg.eigvals(net.build_state_space().A)
    assert lam[np.abs(lam) > 1e-6].real.max() < 0


@pytest.mark.parametrize("nom", FAIBLES)
def test_le_fichier_est_celui_du_script(nom):
    m = _script()
    construire = {f.__name__: f for f in m.CONSTRUCTEURS}[nom]
    _nom, t = construire()
    a, b = tables_vers_json(t), tables_vers_json(charger(nom))
    for cle in ("Y_DER", "T_SMs", "T_GFL1s", "T_IBs", "T_Lines_pu", "T_Lines_SI"):
        assert np.allclose(
            np.array([[float(x) for x in r] for r in a[cle]] or [[0]]),
            np.array([[float(x) for x in r] for r in b[cle]] or [[0]]),
        ), cle
