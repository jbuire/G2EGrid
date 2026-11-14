"""Construction en code, transits, et impact de la distribution sur le transport."""

from __future__ import annotations

import numpy as np
import pytest

from g2egrid.catalogue import charger, charger_reseau
from g2egrid.construction import (
    PQ,
    SLACK,
    equipements,
    poser_bus,
    poser_charge,
    poser_equipement,
    poser_ligne,
    reseau_vide,
    un_slack_tout_pq,
)
from g2egrid.flux import flux_branches
from g2egrid.impact import equivalent_statique, impact_distribution
from g2egrid.tables import assembler


# ----------------------------------------------------------------------
# transits
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nom", ["wscc9", "wscc3sm_cigre3"])
def test_les_transits_bouclent_avec_les_injections(nom):
    """Injection nette d'un bus = transits sortants + ses shunts.

    Si ce bilan ne tombe pas a l'arrondi, les transits affiches sur le
    schema ne sont pas ceux du load flow -- ils racontent un autre reseau.
    """
    net = assembler(charger(nom))
    lf = net.run_loadflow()
    tr = flux_branches(net, lf)
    idx = net._bus_index()
    ysh = {b: 1j * net.bus_specs[b].C for b in idx}
    for el in net.all_elements:
        y = el.shunt_admittance()
        if y is not None:
            ysh[el.buses[0]] += y
    for b, i in idx.items():
        s = sum(complex(t.P_de, t.Q_de) for t in tr if t.de == b)
        s += sum(complex(t.P_vers, t.Q_vers) for t in tr if t.vers == b)
        s += abs(lf.V[i]) ** 2 * np.conj(ysh[b])
        assert abs(lf.S[i] - s) < 1e-10, f"bus {b}"


# ----------------------------------------------------------------------
# construction en code
# ----------------------------------------------------------------------
def test_construire_un_reseau_en_code():
    """Les memes conventions que le constructeur graphique, en quatre lignes."""
    t = reseau_vide("demo")
    a = poser_bus(t, 230e3)
    b = poser_bus(t, 230e3, Pl_W=80e6, Ql_var=20e6)
    poser_ligne(t, a, b, R_pu=0.01, X_pu=0.1)
    poser_equipement(t, "IB", a, slack=True)
    g = poser_equipement(t, "GFM1", b, Sn_VA=100e6, P_W=50e6)
    # l'unite a SON bus terminal et SON transfo
    assert g not in (a, b)
    assert any(r[1] == b and r[2] == g for r in t.T_Transfos)
    net = assembler(t)
    net.solve_equilibrium()
    ss = net.build_state_space()
    assert ss.A.shape[0] > 10
    assert net.lf.mismatch < 1e-8


def test_un_slack_tout_pq_ne_bouge_pas_le_point_de_fonctionnement():
    """Decrire le meme point autrement : les tensions ne changent pas.

    Les PV passent en PQ avec le reactif qu'ils fournissaient ; le load flow
    doit retomber exactement sur les memes tensions.
    """
    t = charger("wscc9")
    net = assembler(t)
    lf = net.run_loadflow()
    idx = net._bus_index()
    injections = {b: lf.S[i] for b, i in idx.items()}
    t2 = t.copie()
    un_slack_tout_pq(t2, slack_bus=1, lf_bus=injections)
    assert [int(r[1]) for r in t2.Y_DER].count(SLACK) == 1
    assert all(int(r[1]) == PQ for r in t2.Y_DER if int(r[0]) != 1)
    lf2 = assembler(t2).run_loadflow()
    assert np.abs(lf2.V - lf.V).max() < 1e-9


def test_changer_de_slack_ne_bouge_pas_le_point_de_fonctionnement():
    """Le slack passe au bus 2 : l'ancien slack garde la puissance qu'il donnait.

    Non-regression : l'ancien slack passait en PQ avec P = 0, et le nouveau
    slack reprenait toute sa production (|V| deplace de 9,2e-2 pu sur wscc9).
    La reference d'angle change de bus : on compare les modules et les ecarts
    d'angle au bus 1, a la tolerance du load flow.
    """
    t = charger("wscc9")
    net = assembler(t)
    lf = net.run_loadflow()
    idx = net._bus_index()
    injections = {b: lf.S[i] for b, i in idx.items()}
    t2 = t.copie()
    un_slack_tout_pq(t2, slack_bus=2, lf_bus=injections)
    assert [int(r[0]) for r in t2.Y_DER if int(r[1]) == SLACK] == [2]
    lf2 = assembler(t2).run_loadflow()
    assert np.abs(np.abs(lf2.V) - np.abs(lf.V)).max() < 1e-9
    d1 = np.angle(lf.V / lf.V[idx[1]])
    d2 = np.angle(lf2.V / lf2.V[idx[1]])
    assert np.abs(d2 - d1).max() < 1e-9


def test_equipements_liste_les_unites():
    t = charger("wscc9")
    eq = equipements(t)
    assert [e["famille"] for e in eq] == ["SM", "SM", "SM"]
    assert {e["bus"] for e in eq} == {1, 2, 3}


def test_charge_speciale_dans_t_loads():
    t = reseau_vide()
    b = poser_bus(t, 20e3)
    poser_charge(t, b, 1e6, 0.2e6, genre="P")
    poser_charge(t, b, 2e6, 0.0, genre="Z")
    assert t.T_Loads[0][4] == "P"
    assert t.Y_DER[0][7] == 2e6
    with pytest.raises(ValueError):
        poser_charge(t, b, 1e6, genre="X")


# ----------------------------------------------------------------------
# impact de la distribution
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def zones():
    return charger_reseau("wscc3sm_zones")


def test_l_equivalent_garde_le_point_de_fonctionnement(zones):
    """C'est ce qui rend la comparaison propre : seule la dynamique change.

    Un equivalent construit sur les charges nominales deplacerait aussi le
    point de fonctionnement, et l'on ne saurait plus ce qu'on mesure.
    """
    im = impact_distribution(zones, par_poste=False)
    assert im.ecart_tension < 1e-10
    assert len(im.soutirage) == 3
    # la distribution SOUTIRE de la puissance active au poste
    assert all(s.real > 0 for s in im.soutirage.values())


def test_une_distribution_passive_ne_change_rien_d_important(zones):
    """Distribution passive : la modeliser en charge est legitime.

    Les modes du transport bougent de moins de 0,5 % en frequence et de
    moins d'un point d'amortissement, la distribution n'y participe pas, et
    rien n'est cree ni degrade dans la bande electromecanique.

    MAIS pas "rien du tout" : la capacite des cables de distribution decale
    une resonance electromagnetique du transport vers 500 Hz, dont
    l'amortissement passe de 5,6 % a ~4 %. L'outil a raison de la signaler ;
    c'est une propriete du reseau, pas un artefact, et elle est hors de la
    bande des modes electromecaniques.
    """
    im = impact_distribution(zones, par_poste=False)
    assert im.paires, "aucun mode du transport apparie"
    for m in im.paires:
        assert abs(m.ecart_f) < 0.5
        assert abs(m.ecart_zeta) < 1.0
        assert m.part_distribution < 5.0
    assert all(m.f_complet > 100.0 for m in im.degrades), "degrade en basse frequence"
    assert all(c.f > 5.0 for c in im.crees), (
        "un mode basse frequence cree par du passif"
    )


def test_equivalent_retire_les_instances(zones):
    eq, soutirage = equivalent_statique(zones)
    assert len(eq.Y_DER) == 9, "il ne doit rester que les bus du transport"
    assert not eq.T_Lines_SI or all(
        int(r[1]) < 100 and int(r[2]) < 100 for r in eq.T_Lines_SI
    )


def test_des_gfm_de_distribution_deplacent_les_modes_du_transport():
    """Avec des convertisseurs, la distribution n'est plus une charge.

    Trois GFM par depart suffisent a deplacer un mode inter-machines de
    plusieurs pourcents et de plusieurs points d'amortissement -- l'erreur
    que commet un modele de transport qui la reduit a une charge.
    """
    r = charger_reseau("wscc3sm_zones")
    z = list(r.zones.values())[0]
    for b in (5, 9, 12):
        poser_equipement(z.tables, "GFM1", b, Sn_VA=4e6, P_W=3e6)
    z.tables.GFM_PARAMS = {"mp": 0.005}
    im = impact_distribution(r, par_poste=False)
    assert im.ecart_tension < 1e-10
    assert max(abs(m.ecart_zeta) for m in im.paires) > 2.0
    assert max(abs(m.ecart_f) for m in im.paires) > 2.0
