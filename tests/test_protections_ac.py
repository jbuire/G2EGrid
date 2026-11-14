"""Relais de surtension des convertisseurs (:func:`g2egrid.temporel.protections_ac`)
et rapport de retombee (:attr:`g2egrid.temporel.Protection.retombee`) : un pic
plus court que la temporisation ne fait pas declencher un relais a retombee,
une surtension tenue le fait ; un relais verrouille declenche sur le pic."""

import numpy as np
import pytest

from g2egrid.catalogue import charger
from g2egrid.tables import assembler
from g2egrid.temporel import (
    Evenement,
    ModeleNonLineaire,
    Protection,
    protections_ac,
)

SEUIL, RETARD = 1.10, 0.3


@pytest.fixture(scope="module")
def modele():
    net = assembler(charger("smib_gfl"))
    net.solve_equilibrium()
    return ModeleNonLineaire(net)


def _simuler(m, evenements, retombee, duree=0.8):
    prot = protections_ac(m.net, v_max=SEUIL, retard=RETARD, retombee=retombee)
    sc = m.scenario(*evenements)
    return m.simuler(
        (0.0, duree),
        uext_fn=sc,
        methode="BDF2",
        pas=1e-4,
        t_eval=np.linspace(0.0, duree, 161),
        protections=prot,
    )


# la tension amont monte de 0,12 pu : le GFL voit ~1,12 pu
PIC = [Evenement(0.1, "IB1.E_ref", 0.12), Evenement(0.2, "IB1.E_ref", -0.12)]
TENUE = [Evenement(0.1, "IB1.E_ref", 0.12)]


def test_le_plan_vise_les_gfl(modele):
    prot = protections_ac(modele.net)
    assert [p.cible for p in prot] == ["GFL2"]
    p = prot[0]
    assert p.genre == "surtension_ac" and p.actions == [("ouverture", "GFL2")]
    assert p.seuil == 1.10 and p.retombee == 0.98
    assert protections_ac(modele.net, familles=("GFM",)) == []
    with pytest.raises(ValueError, match="famille"):
        protections_ac(modele.net, familles=("SM",))
    with pytest.raises(ValueError, match="action"):
        protections_ac(modele.net, action="fusion")


def test_rapport_de_retombee_valide():
    with pytest.raises(ValueError, match="retombee"):
        Protection("surtension_ac", "GFL2", 1.1, 0.2, retombee=1.05)
    with pytest.raises(ValueError, match="retombee"):
        Protection("sous_tension_dc", "bus3", 0.8, 0.2, retombee=0.9)
    Protection("sous_tension_dc", "bus3", 0.8, 0.2, retombee=1.05)


def test_un_pic_plus_court_que_la_temporisation_ne_declenche_pas(modele):
    tr = _simuler(modele, PIC, retombee=0.98)
    assert tr.declenchements == []
    # le transitoire de l'echelon peut faire demarrer et retomber le relais
    # plusieurs fois ; a la fin, chaque demarrage a sa retombee
    assert tr.demarrages and len(tr.retombees) == len(tr.demarrages)
    assert {n for _t, n in tr.demarrages} == {"surtension_ac:GFL2"}
    t_dem, t_ret = tr.demarrages[-1][0], tr.retombees[-1][0]
    assert 0.1 <= t_dem < 0.12 and 0.2 <= t_ret < t_dem + RETARD
    # le GFL est toujours la a la fin
    assert np.isfinite(tr.etat("GFL2.ved")[-1])


def test_un_relais_verrouille_declenche_sur_le_pic(modele):
    tr = _simuler(modele, PIC, retombee=None)
    assert [d[1] for d in tr.declenchements] == ["surtension_ac:GFL2"]
    t_dem = tr.declenchements[0][0]
    ved = tr.etat("GFL2.ved")
    assert np.isnan(ved[tr.t > t_dem + RETARD + 1e-6]).all()


def test_une_surtension_tenue_declenche_au_terme(modele):
    tr = _simuler(modele, TENUE, retombee=0.98)
    assert len(tr.declenchements) == 1
    t_dem, nom, mans = tr.declenchements[0]
    assert nom == "surtension_ac:GFL2" and mans == [("ouverture", "GFL2")]
    # le declenchement part du DERNIER demarrage ; les retombees eventuelles
    # du transitoire le precedent
    assert tr.demarrages[-1] == (t_dem, nom)
    assert all(t < t_dem for t, _n in tr.retombees)
    assert 0.1 <= t_dem < 0.12
    ved = tr.etat("GFL2.ved")
    assert np.isfinite(ved[tr.t < t_dem + RETARD - 1e-6]).all()
    assert np.isnan(ved[tr.t > t_dem + RETARD + 1e-6]).all()
