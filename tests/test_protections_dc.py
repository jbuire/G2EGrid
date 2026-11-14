"""Defauts continus, blocage des convertisseurs et protections DC
(:mod:`g2egrid.temporel`) : manoeuvres ``defaut_dc``, ``elimination_dc``,
``blocage`` ; protections surveillees pendant la simulation ; repli en
impedance des charges sous 0,7 pu."""

import numpy as np
import pytest

from g2egrid import catalogue
from g2egrid.temporel import (
    Manoeuvre,
    ModeleNonLineaire,
    Protection,
    protections_dc,
)


@pytest.fixture(scope="module")
def dcs1():
    net = catalogue.construire("cigre_b4_dcs1")
    return net, ModeleNonLineaire(net)


def _v(tr, nom, t):
    j = tr.etats.index(nom)
    return float(tr.x[min(np.searchsorted(tr.t, t), len(tr.t) - 1), j])


def test_defaut_franc_bloque_les_stations(dcs1):
    """Defaut franc pole a pole au noeud C1 : Cm-C1 se bloque d'abord (la
    tension de son noeud s'effondre en 0,2 ms), puis Cm-A1, quand l'onde a
    traverse les 200 km de cable."""
    net, m = dcs1
    sc = m.scenario(Manoeuvre(0.02, "defaut_dc", "bus221"))
    tr = m.simuler((0, 0.06), uext_fn=sc, methode="BDF2", pas=2.5e-5,
                   protections=protections_dc(net))  # fmt: skip
    noms = [d[1] for d in tr.declenchements]
    assert noms[:2] == ["sous_tension_dc:bus221", "sous_tension_dc:bus201"]
    t_c1, t_a1 = tr.declenchements[0][0], tr.declenchements[1][0]
    assert 0.02 < t_c1 < 0.0205 < t_a1 < 0.025
    assert tr.declenchements[0][2] == [("blocage", "GFL121")]


def test_blocage_eteint_le_courant(dcs1):
    """Un convertisseur bloque ramene son courant a zero en quelques
    millisecondes, par sa propre boucle de courant, et ne prend plus rien a son
    noeud continu."""
    net, m = dcs1
    sc = m.scenario(Manoeuvre(0.02, "blocage", "GFL121"))
    tr = m.simuler((0, 0.06), uext_fn=sc, methode="BDF2", pas=5e-5)
    i0 = np.hypot(_v(tr, "GFL121.isd", 0.019), _v(tr, "GFL121.isq", 0.019))
    i1 = np.hypot(_v(tr, "GFL121.isd", 0.05), _v(tr, "GFL121.isq", 0.05))
    assert i0 > 0.4 and i1 < 0.01
    # plus de transit : le cable se decharge
    assert abs(_v(tr, "CDC1_201_221.i", 0.059)) < 0.05 * abs(
        _v(tr, "CDC1_201_221.i", 0.019)
    )


def test_defaut_puis_elimination(dcs1):
    """Un defaut resistif, elimine 10 ms plus tard : la tension continue
    revient, sans qu'aucune station ne se bloque."""
    net, m = dcs1
    sc = m.scenario(
        Manoeuvre(0.02, "defaut_dc", "bus221", R=0.5),
        Manoeuvre(0.03, "elimination_dc", "bus221"),
    )
    tr = m.simuler((0, 0.2), uext_fn=sc, methode="BDF2", pas=5e-5,
                   protections=protections_dc(net))  # fmt: skip
    assert not tr.declenchements
    assert _v(tr, "bus221.v", 0.025) < 0.97
    assert abs(_v(tr, "bus221.v", 0.199) - _v(tr, "bus221.v", 0.019)) < 0.02


def test_disjoncteur_dc_ouvre_le_cable(dcs1):
    """Une protection de surintensite a seuil bas sur le cable : le
    disjoncteur l'ouvre ``retard`` apres le franchissement."""
    net, m = dcs1
    p = Protection("surintensite_dc", "CDC1_201_221", seuil=4.5, retard=2e-3)
    sc = m.scenario(Manoeuvre(0.02, "defaut_dc", "bus221", R=1e-3))
    tr = m.simuler((0, 0.05), uext_fn=sc, methode="BDF2", pas=2.5e-5, protections=[p])
    (t, nom, mans), *_ = tr.declenchements
    assert nom == "surintensite_dc:CDC1_201_221" and mans == [
        ("ouverture", "CDC1_201_221")
    ]
    j = tr.etats.index("CDC1_201_221.i")
    assert np.isnan(tr.x[-1, j])  # le cable n'est plus dans le modele
    assert np.isfinite(
        tr.x[np.searchsorted(tr.t, t + 1e-3), j]
    )  # ouvert apres le retard


def test_manoeuvres_et_protections_refusees(dcs1):
    net, m = dcs1
    with pytest.raises(ValueError, match="protection"):
        Protection("surchauffe", "x", 1.0)
    for man, motif in (
        (Manoeuvre(0.01, "defaut_dc", "bus999"), "noeud continu"),
        (Manoeuvre(0.01, "elimination_dc", "bus221"), "pas de defaut"),
        (Manoeuvre(0.01, "blocage", "CDC1_201_221"), "convertisseur"),
    ):
        with pytest.raises(ValueError, match=motif):
            m.simuler((0, 0.02), uext_fn=m.scenario(man), methode="BDF2", pas=1e-4)


def test_plan_type(dcs1):
    """Un disjoncteur par cable, sous-tension et surintensite par station."""
    net, _ = dcs1
    genres = sorted(p.nom for p in protections_dc(net))
    assert genres == [
        "sous_tension_dc:bus201",
        "sous_tension_dc:bus221",
        "surintensite_convertisseur:GFL101",
        "surintensite_convertisseur:GFL121",
        "surintensite_dc:CDC1_201_221",
    ]


def test_blocage_d_une_station_du_cigre_b4():
    """Bloquer Cb-B2 (1,5 GW) : le poste B2 plonge puis se rétablit, les
    charges passent en impedance sous 0,7 pu le temps du creux, et le reseau
    continu reprend un regime sans cascade de declenchements."""
    net = catalogue.construire("cigre_b4")
    m = ModeleNonLineaire(net)
    sc = m.scenario(Manoeuvre(0.02, "blocage", "GFL113"))
    tr = m.simuler((0, 0.08), uext_fn=sc, methode="BDF2", pas=5e-5,
                   protections=protections_dc(net))  # fmt: skip
    assert {d[1] for d in tr.declenchements} <= {"sous_tension_dc:bus213"}
    assert any(b[2] == "mode_bt" for b in tr.basculements)
    j = tr.etats.index("bus12.vd")
    creux = tr.x[(tr.t > 0.02) & (tr.t < 0.03), j]
    assert creux.min() < 0.5 < _v(tr, "bus12.vd", 0.079)
    assert 0.9 < _v(tr, "bus202.v", 0.079) < 1.3


def test_api_temporel_avec_protections():
    analyse = pytest.importorskip("g2egrid_api.analyse")
    d = analyse.temporel(
        catalogue.charger("cigre_b4_dcs1"),
        "",
        duree=0.04,
        n=81,
        voie="non_lineaire",
        evenements=[{"instant": 0.02, "genre": "defaut_dc", "cible": "bus221"}],
        etats=["bus221.v"],
        pas=2.5e-5,
        protections_dc=True,
    )
    assert d["integrateur"] == "BDF2"
    assert d["declenchements"][0]["protection"] == "sous_tension_dc:bus221"
    assert d["declenchements"][0]["manoeuvres"] == [
        {"genre": "blocage", "cible": "GFL121"}
    ]
    s = analyse.manoeuvrables(catalogue.construire("cigre_b4_dcs1"))
    assert s["bus_dc"] == [201, 221]
    assert {"nom": "GFL121", "genre": "convertisseur"} in s["elements"]
