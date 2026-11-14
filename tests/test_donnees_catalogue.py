"""Les donnees des reseaux publies du catalogue contre leur source.

* ``ieee33`` contre ``case33bw`` (MATPOWER, Baran et Wu 1989) : la solution de
  reference est celle de ``pandapower.networks.case33bw()`` gardee dans
  ``tests/donnees/interop/references.json`` (FR-021 : les impedances de ligne du
  catalogue etaient 1,157 fois trop petites). Load flow en noeuds ideaux, sans
  la resistance de causalite, qui est une autre question (FR-005) ;
* ``ieee14`` contre ``case14`` (MATPOWER), solution de
  ``pandapower.networks.case14()`` : les prises des transformateurs 4-7, 4-9 et
  5-6 (0,978, 0,969, 0,932) etaient ignorees (1,6e-2 pu d'ecart).
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
import pytest

from g2egrid.catalogue import charger
from g2egrid.tables import assembler

REF = json.loads(
    (
        pathlib.Path(__file__).resolve().parent
        / "donnees"
        / "interop"
        / "references.json"
    ).read_text()
)


def _lf_ideal(nom: str):
    t = charger(nom)
    t.NOEUDS = "ideal"
    net = assembler(t)
    return net.run_loadflow(), net._bus_index()


def _ecarts(nom: str, ref: str, decalage: int = 0):
    lf, ix = _lf_ideal(nom)
    r = REF[ref]
    dv, da = [], []
    for b, vm, va in zip(r["bus"], r["vm"], r["va"]):
        V = lf.V[ix[b + decalage]]
        dv.append(abs(abs(V) - vm))
        da.append(np.degrees(np.angle(V)) - va)
    return max(dv), max(np.abs(da))


def test_ieee33_ligne_1_2_de_case33bw():
    """Ligne 1-2 : 0,0922 + j0,0470 ohm (Baran et Wu), soit sur la base du
    reseau (12,66 kV, 10 MVA, 16,03 ohm) R = 0,005753 pu. Le catalogue portait
    0,004967 pu (0,0796 ohm)."""
    t = charger("ieee33")
    _id, i, j, _km, r_pu, x_pu, _b, Un, Sb_l = t.T_Lines_pu[0]
    zb = Un**2 / Sb_l
    assert (i, j) == (1, 2)
    assert r_pu * zb == pytest.approx(0.0922, rel=1e-12)
    assert x_pu * zb == pytest.approx(0.0470, rel=1e-12)


def test_ieee33_retrouve_case33bw():
    """Load flow du catalogue (noeuds ideaux) = case33bw de MATPOWER, a 1e-12
    pres en module et en angle (mesure : 2e-14) ; tension minimale 0,91309 au
    bus 18, la valeur publiee."""
    dv, da = _ecarts("ieee33", "case33bw")
    assert dv < 1e-12 and da < 1e-10
    assert min(REF["case33bw"]["vm"]) == pytest.approx(0.91309, abs=1e-5)


def test_ieee33_pertes_de_case33bw():
    """Pertes de la base de Baran et Wu : 202,68 kW (valeur courante de la
    litterature pour ce cas, retrouvee par l'ecart au slack de pandapower)."""
    lf, _ = _lf_ideal("ieee33")
    charge = sum(r[7] for r in charger("ieee33").Y_DER)
    pertes = lf.S.real.sum() * charger("ieee33").Sb
    assert pertes == pytest.approx(202.68e3, abs=0.05e3)
    assert charge == pytest.approx(3.715e6, rel=1e-6)


def test_ieee14_prises_des_transformateurs():
    """Les trois transformateurs a prise de case14 sont des ``T_Transfos`` dont
    le rapport (U1/U2 aux bases de tension egales) est la prise de MATPOWER."""
    t = charger("ieee14")
    prises = {(int(r[1]), int(r[2])): r[4] / r[5] for r in t.T_Transfos}
    assert prises == pytest.approx({(4, 7): 0.978, (4, 9): 0.969, (5, 6): 0.932})
    assert not {(4, 7), (4, 9), (5, 6)} & {(int(r[1]), int(r[2])) for r in t.T_Lines_pu}


def test_ieee14_retrouve_case14():
    """Load flow du catalogue (noeuds ideaux) = case14 a 1e-4 pu pres (mesure :
    7,3e-5 en module, au bus 4). Ce qui reste : le bus 4 de case14 a une demande
    reactive de -3,9 Mvar, qu'une charge R-L ne peut pas porter ; le catalogue la
    remplace par 0,478 Mvar et un shunt de 0,04378 pu (admittance constante :
    -4,06 Mvar a V4 = 1,018 pu), plus la resistance de 1e-4 pu des branches que
    MATPOWER a sans resistance (6,5e-5 sans elle, 7,3e-5 avec)."""
    dv, da = _ecarts("ieee14", "case14")
    assert dv < 1e-4 and da < 0.01  # angle en degres
