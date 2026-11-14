"""L'onglet Evenements cote serveur : :func:`g2egrid_api.analyse.apres_evenement`
(les grandeurs de chaque bus et de chaque unite apres un evenement) et
:func:`g2egrid_api.analyse.cascade` (la cascade quasi statique)."""

from __future__ import annotations

import pytest

from g2egrid.catalogue import charger
from g2egrid_api import analyse

#: la tension amont monte de 0,12 pu a 0,1 s : le GFL voit ~1,12 pu
ECHELON = [
    {"genre": "echelon", "instant": 0.1, "entree": "IB1.E_ref", "amplitude": 0.12}
]
RELAIS = {"v_max": 1.10, "retard": 0.3, "retombee": 0.98, "familles": ["GFL"]}


@pytest.fixture(scope="module")
def smib():
    return charger("smib_gfl")


@pytest.fixture(scope="module")
def surtension(smib):
    return analyse.apres_evenement(
        smib, ECHELON, duree=0.8, n=161, pas=1e-4, relais=RELAIS
    )


def test_le_relais_declenche_et_le_tableau_le_dit(surtension):
    r = surtension
    assert [d["protection"] for d in r["declenchements"]] == ["surtension_ac:GFL2"]
    d = r["declenchements"][0]
    assert d["t"] == pytest.approx(d["t_demarrage"] + 0.3)
    assert 0.1 <= d["t_demarrage"] < 0.12
    u = {x["nom"]: x for x in r["unites"]}
    assert u["GFL2"]["sortie"] == pytest.approx(d["t"])
    assert u["GFL2"]["Pfin"] == 0.0 and u["GFL2"]["P0"] > 0.4
    assert u["IB1"]["sortie"] is None
    # la serie d'une unite sortie vaut 0 apres sa sortie, pas null
    assert r["series"]["P"]["GFL2"][-1] == 0.0
    assert r["resume"]["n_declenchements"] == 1


def test_les_bus_extremes_et_hors_bornes(surtension):
    r = surtension
    b = {x["bus"]: x for x in r["bus"]}
    assert set(b) == {1, 2}
    assert b[2]["Vmax"] > 1.10 and b[2]["t_max"] >= 0.1
    assert b[2]["t_hors"] > 0.3  # au-dessus de 1,10 pu toute la temporisation
    assert b[2]["V0"] == pytest.approx(1.0, abs=1e-3)
    assert b[2]["dV"] == pytest.approx(b[2]["Vfin"] - b[2]["V0"])
    assert r["resume"]["v_max"]["bus"] == 2
    # une consigne n'a pas d'equivalent en load flow
    assert r["regime_lf"]["ok"] is False and "consigne" in r["regime_lf"]["raison"]
    assert b[2]["V_lf"] is None
    assert len(r["series"]["V"]["2"]) == len(r["t"]) == 161


def test_une_ouverture_se_compare_au_load_flow(smib):
    r = analyse.apres_evenement(
        smib,
        [{"genre": "ouverture", "instant": 0.1, "cible": "GFL2"}],
        duree=1.0,
        n=101,
    )
    assert r["regime_lf"] == {"ok": True, "ouverts": ["GFL2"]}
    for b in r["bus"]:
        assert abs(b["ecart_lf"]) < 1e-3
    # la frequence est celle du bus infini
    assert r["frequence"]["bus"] == [1]
    assert r["frequence"]["min"] == pytest.approx(50.0, abs=0.01)
    assert r["relais"] is None and r["demarrages"] == []


def test_les_erreurs_sont_dites(smib):
    with pytest.raises(ValueError, match="vide"):
        analyse.apres_evenement(smib, [])
    with pytest.raises(ValueError, match="entree inconnue"):
        analyse.apres_evenement(
            smib,
            [{"genre": "echelon", "instant": 0.1, "entree": "X.y", "amplitude": 1}],
        )
    with pytest.raises(ValueError, match="bornes"):
        analyse.apres_evenement(smib, ECHELON, v_bas=1.1, v_haut=0.9)


def test_la_cascade_quasi_statique():
    t = charger("wscc3sm_cigre3_surtension")
    c = analyse.cascade(t, ["GFL303", "GFL304", "GFL305"], ordre="une_a_une")
    assert c["issue"] == "arretee"
    assert sorted(c["declenchees"]) == ["GFL307", "GFL311", "GFL313", "GFL314"]
    assert c["bus_unite"]["GFL314"] == 314
    e0 = c["etapes"][0]
    assert e0["ouvertes"] == ["GFL303", "GFL304", "GFL305"]
    assert e0["v_max"] > c["seuil"] and "GFL314" in e0["au_dessus"]
    assert set(c["V0"]) == set(e0["V"])
    with pytest.raises(ValueError, match="aucune"):
        analyse.cascade(t, [])


def test_les_routes(smib):
    fastapi = pytest.importorskip("fastapi")
    from g2egrid_api import app as routes
    from g2egrid_api.schemas import CorpsApresEvenement, CorpsCascade

    tables = analyse.tables_vers_json(smib)
    with pytest.raises(fastapi.HTTPException) as err:
        routes.apres_evenement(CorpsApresEvenement(tables=tables, evenements=[]))
    assert err.value.status_code == 422
    with pytest.raises(fastapi.HTTPException) as err:
        routes.cascade(CorpsCascade(tables=tables, initiales=["GFL2"], ordre="hasard"))
    assert err.value.status_code == 422
    assert {"apres_evenement", "cascade"} <= set(routes.GENRES_TACHE)
