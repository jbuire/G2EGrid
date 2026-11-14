"""Transport et distribution : les reseaux, les GFL au hasard, les modeles
reduits (voir g2egrid.multi_echelle)."""

from __future__ import annotations

import numpy as np
import pytest

from g2egrid.multi_echelle import (
    PART_REMPLACEE,
    ajouter_gfl,
    comparer_distribution,
    comparer_transport,
    distribution_seule,
    equivalent_hybride,
    ieee14_ieee33,
    impedance_thevenin,
    wscc9_cigre6,
)
from g2egrid.robustesse import _spectre
from g2egrid.tables import assembler
from g2egrid.zones import aplatir


@pytest.fixture(scope="module")
def wscc():
    return wscc9_cigre6()


@pytest.fixture(scope="module")
def plat(wscc):
    return aplatir(wscc)


@pytest.fixture(scope="module")
def complet(plat):
    return _spectre(plat)


def test_six_cigre_remplacent_les_charges(wscc, plat):
    assert len(wscc.zones) == 6
    assert sorted(i.bus_amont for i in wscc.instances) == [5, 5, 6, 6, 8, 8]
    # il ne reste au poste que la part non reprise par les distributions
    charge = {5: 125e6, 6: 90e6, 8: 100e6}
    for r in plat.Y_DER:
        if int(r[0]) in charge:
            assert r[7] == pytest.approx(
                (1 - PART_REMPLACEE) * charge[int(r[0])], rel=1e-9
            )
    # les six zones ne sont pas identiques (pas de paires symetriques)
    charges = [tuple(x[7] for x in z.tables.Y_DER) for z in wscc.zones.values()]
    assert len(set(charges)) == 6


def test_les_reseaux_d_etude_suivent_les_noeuds_par_defaut(plat):
    """wscc9 garde la resistance de Node.m (parite MATLAB) ; le reseau
    multi-echelle bati dessus est un reseau d'etude : noeuds par defaut,
    comme son fichier du catalogue (FR-005, FR-008)."""
    from g2egrid.catalogue import charger
    from g2egrid.elements.bus import NOEUDS_DEFAUT

    assert charger("wscc9").NOEUDS != NOEUDS_DEFAUT
    assert plat.NOEUDS == charger("wscc9_cigre6").NOEUDS == NOEUDS_DEFAUT
    assert aplatir(ieee14_ieee33()).NOEUDS == NOEUDS_DEFAUT


def test_les_deux_reseaux_sont_stables_et_convergent():
    for r in (wscc9_cigre6(), ieee14_ieee33()):
        net = assembler(aplatir(r))
        lf = net.run_loadflow()
        assert lf.converged
        assert abs(lf.V).min() > 0.8
        net.solve_equilibrium()
        lam = net.build_state_space().eig()
        assert lam[np.abs(lam) > 1e-6].real.max() < 0


def test_gfl_au_hasard_reproductibles_et_sur_bus_libres(wscc):
    r1, a1 = ajouter_gfl(wscc, 8, graine=4)
    r2, a2 = ajouter_gfl(wscc, 8, graine=4)
    assert a1 == a2 and len(a1) == 8
    assert len({(z, b) for z, b, _ in a1}) == 8  # jamais deux fois le meme bus
    for z, b, p in a1:
        zone = r1.zones[z]
        assert b != zone.bus_pcc
        assert any(int(x[1]) == b for x in zone.tables.T_GFL1s)
        pg = next(float(x[5]) for x in zone.tables.Y_DER if int(x[0]) == b)
        assert pg == pytest.approx(p)
    # le reseau d'origine n'a pas bouge
    assert sum(len(z.tables.T_GFL1s) for z in wscc.zones.values()) == 18


def test_distribution_seule_bus_infini_au_poste(wscc, plat):
    red, z = distribution_seule(wscc, 0, plat=plat)
    amont = wscc.instances[0].bus_amont
    assert [int(r[0]) for r in red.Y_DER if int(r[1]) == 1] == [amont]
    assert [int(x[1]) for x in red.T_IBs] == [amont]
    assert not red.T_SMs  # rien du transport
    assert 0 < z.imag < 1 and z.real >= 0
    assert assembler(red).run_loadflow().converged


def test_thevenin_d_un_bus_infini_est_son_impedance():
    from g2egrid.catalogue import charger

    t = charger("smib")
    ib = t.T_IBs[0]
    zib = complex(ib[5], ib[6]) * t.Sb / ib[2]
    z = impedance_thevenin(t, int(ib[1]))
    # le bus infini en parallele avec le reste : |Z| ne peut que baisser
    assert abs(z) <= abs(zib) * (1 + 1e-9)


def test_les_modeles_reduits_retrouvent_les_modes(wscc, plat, complet):
    ct = comparer_transport(wscc, plat=plat, complet=complet)
    ch = comparer_transport(wscc, plat=plat, complet=complet, hybride=True)
    cd = comparer_distribution(wscc, 0, plat=plat, complet=complet)
    assert ct.modes and cd.modes
    # l'equivalent hybride fait mieux que les charges R-L
    assert ch.ecart_max() < ct.ecart_max() < 0.1
    assert cd.ecart_max() < 0.02


def test_equivalent_hybride_garde_le_point_de_fonctionnement(wscc, plat):
    eq = equivalent_hybride(wscc, plat=plat)
    net_c, net_e = assembler(plat), assembler(eq)
    lf_c, lf_e = net_c.run_loadflow(), net_e.run_loadflow()
    ic, ie = net_c._bus_index(), net_e._bus_index()
    ecart = max(abs(lf_c.V[ic[b]] - lf_e.V[ie[b]]) for b in ie)
    assert ecart < 1e-6
    assert any(x[4] == "P" and x[2] < 0 for x in eq.T_Loads)


def test_partition_perturbation_singuliere(wscc, plat):
    """Le transport voit la distribution RAPIDE : quasi statique, il la
    rend mieux que figee. La distribution voit le transport LENT : fige, il
    le rend mieux que quasi statique."""
    from g2egrid.multi_echelle import comparer_partition, partition

    P = partition(wscc, plat)
    n = len(P[0].etats)
    assert 0 < P[0].transport.sum() < n
    ecart = {
        (c, q): comparer_partition(wscc, c, q, part=P).ecart_max()
        for c in ("transport", "distribution")
        for q in (False, True)
    }
    assert ecart["transport", True] < ecart["transport", False]
    assert ecart["transport", True] < 0.02
    assert ecart["distribution", False] < ecart["distribution", True]
    assert ecart["distribution", False] < 0.05


def test_penetration_ne_masque_que_l_absence_de_point_de_fonctionnement(
    monkeypatch, wscc
):
    """Un palier sans point de fonctionnement porte son ``erreur`` ; une
    erreur de programmation remonte."""
    from g2egrid import multi_echelle

    def sans_pf(r, instances=None):
        raise RuntimeError("load flow non converge\nsuite")

    monkeypatch.setattr(multi_echelle, "evaluer", sans_pf)
    paliers = multi_echelle.penetration(wscc, n_max=5, pas=5)
    assert len(paliers) == 1 and paliers[0].erreur == "load flow non converge"

    def bug(r, instances=None):
        raise AttributeError("faute de frappe")

    monkeypatch.setattr(multi_echelle, "evaluer", bug)
    with pytest.raises(AttributeError):
        multi_echelle.penetration(wscc, n_max=5, pas=5)
