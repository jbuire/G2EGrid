"""Cascade de surtension en regimes permanents (:mod:`g2egrid.cascade`), sur
le cas d'etude ``wscc3sm_cigre3_surtension`` (``etudes/surtension_cascade``) :
33 GFL a tan(phi) = -0,3, tension haute sans depassement.

Ce que le cas affirme, et que ces tests tiennent :

- la perte de GFL qui ABSORBENT du reactif fait monter la tension du poste
  110 kV : c'est leur Q qui manque ; leur P, seul, la ferait baisser ;
- la cascade depend de l'ordre des declenchements : ensemble (vague) elle
  gagne les trois zones, une a une elle s'arrete dans la zone 3.
"""

from __future__ import annotations

import pytest

from g2egrid.cascade import cascade_surtension, tensions_bus, unites_protegees
from g2egrid.catalogue import charger

NOM = "wscc3sm_cigre3_surtension"
INITIALES = ["GFL303", "GFL304", "GFL305"]


@pytest.fixture(scope="module")
def reseau():
    return charger(NOM)


def test_le_depart_est_haut_sans_depasser(reseau):
    V, ok = tensions_bus(reseau)
    assert ok
    hta = [v for b, v in V.items() if b >= 100 and b % 100 not in (0, 50)]
    assert 1.095 < max(hta) < 1.10
    assert len(unites_protegees(reseau)) == 33


def test_perdre_le_reactif_absorbe_fait_monter_la_tension(reseau):
    V0, _ = tensions_bus(reseau)
    bus = {n: b for _t, _i, b, n in unites_protegees(reseau)}

    def sans(p, q):
        t = reseau.copie()
        cibles = {bus[n] for n in INITIALES}
        t.Y_DER = [
            tuple(
                list(r[:5])
                + [0.0 if (p and int(r[0]) in cibles) else r[5]]
                + [0.0 if (q and int(r[0]) in cibles) else r[6]]
                + list(r[7:])
            )
            for r in reseau.Y_DER
        ]
        return tensions_bus(t)[0]

    t = reseau.copie()
    t.OUVERTS = [
        (tab, i) for tab, i, _b, n in unites_protegees(reseau) if n in INITIALES
    ]
    total, perte_q, perte_p = tensions_bus(t)[0], sans(False, True), sans(True, False)
    # au poste 110 kV de la zone 3 : le Q fait tout (X/R eleve du poste)
    assert total[300] - V0[300] > 5e-3
    assert perte_q[300] - V0[300] > 5e-3
    assert abs(perte_p[300] - V0[300]) < 0.1 * (perte_q[300] - V0[300])
    # le long du depart qui les portait (cables a R/X = 0,7) : les deux
    # effets s'opposent, de 25 mpu chacun, et le Q l'emporte de peu
    assert perte_p[305] - V0[305] < -0.02
    assert perte_q[305] - V0[305] > 0.02
    assert 0 < total[305] - V0[305] < 0.005
    # l'autre zone monte aussi, par le transport, mais bien moins
    assert 0 < total[200] - V0[200] < 0.3 * (total[300] - V0[300])


def test_la_cascade_depend_de_l_ordre(reseau):
    vague = cascade_surtension(reseau, INITIALES)
    une = cascade_surtension(reseau, INITIALES, ordre="une_a_une")
    assert vague.issue == une.issue == "arretee"
    assert len(vague.declenchees) == 21
    zones = {n[3] for n in vague.declenchees}
    assert zones == {"1", "2", "3"}
    assert sorted(une.declenchees) == ["GFL307", "GFL311", "GFL313", "GFL314"]
    # une a une : un declenchement par regime, le plus haut d'abord
    assert [e.ouvertes for e in une.etapes[1:]] == [
        ["GFL314"],
        ["GFL311"],
        ["GFL313"],
        ["GFL307"],
    ]
    assert une.etapes[-1].au_dessus == []
    # en vague, chaque regime ouvre tout ce qui depassait au precedent
    for avant, apres in zip(vague.etapes, vague.etapes[1:]):
        assert apres.ouvertes == avant.au_dessus


def test_un_seuil_haut_arrete_tout(reseau):
    c = cascade_surtension(reseau, INITIALES, v_max=1.2)
    assert c.issue == "arretee" and c.declenchees == []
    assert len(c.etapes) == 1 and c.etapes[0].ouvertes == INITIALES


def test_les_erreurs_sont_dites(reseau):
    with pytest.raises(ValueError, match="inconnue"):
        cascade_surtension(reseau, ["GFL999"])
    with pytest.raises(ValueError, match="ordre"):
        cascade_surtension(reseau, INITIALES, ordre="hasard")
    with pytest.raises(ValueError, match="famille"):
        cascade_surtension(reseau, INITIALES, familles=("SM",))
