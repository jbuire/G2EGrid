"""La topologie d'un cable continu : colonnes optionnelles ``Config`` et
``Circuits`` de ``T_LinesDC`` (:func:`g2egrid.tables.cable_dc`)."""

import pytest

from g2egrid.catalogue import charger
from g2egrid.tables import CONFIGS_DC, assembler, cable_dc

R, L, C = 0.011, 2.615, 0.2185  # cable 200 kV du CIGRE B4, par pole


def test_sept_colonnes_inchangees():
    ligne = (1, 201, 221, 200.0, R, L, C)
    assert cable_dc(ligne) == ligne


@pytest.mark.parametrize(
    "config, k_rl, k_c",
    [
        ("equivalent", 1.0, 1.0),
        ("monopole_sym", 2.0, 0.5),
        ("bipole", 2.0, 0.5),
        ("monopole_metallique", 2.0, 1.0),
        ("monopole_terre", 1.0, 1.0),
    ],
)
def test_configurations(config, k_rl, k_c):
    _, _, _, km, r, ind, c = cable_dc((1, 201, 221, 200, R, L, C, config, 1))
    assert km == 200.0
    assert (r, ind, c) == pytest.approx((k_rl * R, k_rl * L, k_c * C))


def test_circuits_en_parallele():
    _, _, _, _, r, ind, c = cable_dc((1, 1, 2, 100, R, L, C, "bipole", 2))
    assert (r, ind, c) == pytest.approx((R, L, C / 2 * 2))


def test_refus():
    with pytest.raises(ValueError, match="inconnue"):
        cable_dc((1, 1, 2, 100, R, L, C, "tripole"))
    with pytest.raises(ValueError, match="au moins 1"):
        cable_dc((1, 1, 2, 100, R, L, C, "bipole", 0))
    assert set(CONFIGS_DC) >= {"equivalent", "monopole_sym", "bipole"}


def test_cigre_decrit_par_pole():
    """Le CIGRE B4 du catalogue decrit ses lignes PAR POLE avec leur
    configuration ; l'assembleur en tire le conducteur equivalent du
    tableau 6 : boucle 2r, capacite c/2, circuits en parallele."""
    t = charger("cigre_b4")
    configs = {r[7] for r in t.T_LinesDC}
    assert configs == {"monopole_sym", "bipole"}
    net = assembler(t)
    zb = 400e3**2 / t.Sb
    cable = next(e for e in net.elements if e.name == "CDC1_201_221")
    assert cable.pdict["R"] == pytest.approx(2 * R * 200 / zb)
    # DCS3, Bb-A1 - Bb-B1 : 400 km, deux circuits
    zb = 800e3**2 / t.Sb
    ligne = next(e for e in net.elements if e.name == "CDC2_202_211")
    assert ligne.pdict["R"] == pytest.approx(2 * 0.0114 * 400 / 2 / zb)
