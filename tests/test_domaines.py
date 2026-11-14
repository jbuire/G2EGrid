"""Le domaine des modes (:mod:`g2egrid.domaines`) : alternatif, continu,
couple AC/DC, et les iles qu'ils traversent."""

import numpy as np
import pytest

from g2egrid import catalogue, domaines


@pytest.fixture(scope="module")
def b4():
    return catalogue.construire("cigre_b4")


def test_reseau_alternatif_tout_ac():
    net = catalogue.construire("wscc9")
    ms = domaines.classer_modes(net)
    assert ms and all(m.domaine == "AC" and m.part_dc == 0 for m in ms)
    assert domaines.bilan(ms)["DC"] == 0


def test_domaine_des_etats(b4):
    """Noeuds et cables continus, boucle Vdc des GFL a port continu : DC ;
    filtre, PLL, machines, reseau alternatif : AC."""
    ss = b4.build_state_space()
    doms, _ = domaines.etiquettes_etats(b4, list(ss.states))
    d = dict(zip(ss.states, doms))
    assert d["bus213.v"] == "DC"
    assert d["CDC7_213_214.i"] == "DC"
    assert d["GFL113.M_d"] == "DC"
    assert d["GFL113.isd"] == "AC"
    assert d["GFL113.M_pll"] == "AC"
    assert d["SM11.Dwr"] == "AC"
    assert d["bus11.vd"] == "AC"


def test_iles_du_cigre_b4(b4):
    """Six iles alternatives (A, B, C1+C2, D, E, F, postes de conversion compris) et
    deux reseaux continus (DCS1, et DCS2+DCS3 relies par Cd-E1)."""
    ile = domaines.iles(b4)
    assert len(set(ile["AC"].values())) == 6
    assert set(ile["DC"].values()) == {201, 202}
    # la station offshore est dans l'ile de son parc
    assert ile["AC"][121] == ile["AC"][21] == ile["AC"][22]


def test_classement_complet(b4):
    """Chaque mode physique a un domaine ; les trois domaines sont presents
    et les parts somment a 1."""
    ss = b4.build_state_space()
    ms = domaines.classer_modes(b4, ss)
    lam = ss.eig()
    assert len(ms) == int((np.abs(lam) > domaines.JAUGE).sum())
    b = domaines.bilan(ms)
    assert b["AC"] > 0 and b["DC"] > 0 and b["AC/DC"] > 0
    for m in ms:
        assert m.part_ac + m.part_dc == pytest.approx(1.0)
        assert sum(p for _, p in m.iles) == pytest.approx(1.0)


def test_mode_electromecanique_inter_iles(b4):
    """Le mode electromecanique de A vers 1,9 Hz est en partie porte par le
    reseau B, que seul le continu relie a A : un mode inter-iles."""
    ms = domaines.classer_modes(b4, oscillants=True)
    m = next(m for m in ms if 1.8 < m.f < 2.0 and m.iles[0][0] == "AC 1")
    assert m.domaine == "AC" and m.inter_iles and m.couple
    assert dict(m.iles)["AC 10"] > domaines.SEUIL


def test_resonances_des_cables_sont_dc(b4):
    """Les resonances du reseau continu (dizaines a centaines de hertz) sont
    portees a plus de 90 % par les noeuds et les cables continus."""
    ms = domaines.classer_modes(b4, oscillants=True)
    dc = [m for m in ms if m.domaine == "DC" and m.f > 10]
    assert len(dc) >= 3
    assert all(m.part_dc > 0.9 for m in dc)


def test_impact_du_continu():
    """Continu rigide : les modes du continu sont retires, et des modes
    alternatifs sont deplaces -- la marque du couplage."""
    c = domaines.impact_continu(catalogue.charger("cigre_b4"))
    assert c.n_etats_reduit < c.n_etats
    assert c.retires
    assert max(m.ecart for m in c.retrouves if m.ecart is not None) > 0.05


def test_impact_sans_continu_refuse():
    with pytest.raises(ValueError, match="pas de reseau continu"):
        domaines.impact_continu(catalogue.charger("wscc9"))


def test_api_t_eig_porte_les_domaines():
    analyse = pytest.importorskip("g2egrid_api.analyse")
    d = analyse.t_eig(catalogue.charger("hvdc2t"))
    assert d["domaines"]["DC"] >= 1
    assert all(m["domaine"] in ("AC", "DC", "AC/DC") for m in d["modes"])
    assert analyse.t_eig(catalogue.charger("wscc9"))["domaines"] is None
