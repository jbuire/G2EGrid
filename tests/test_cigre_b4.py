"""Le CIGRE B4 DC Grid Test System du catalogue, contre le load flow de
reference de DIgSILENT PowerFactory (« CIGRE WG B4.57 and B4.58 DC grid test
system, Part I », tableaux 13 et 17), et les deux modes de convertisseur
qu'il a fait entrer : ``Vtheta`` (une station forme son ile offshore) et
``droop`` (statisme tension continue / puissance).

Le reseau est construit par ``etudes/cigre_b4/construire.py`` ; ses ecarts
au modele de reference (bipoles agreges, DC/DC ideaux, stations GFM/GFL) y
sont listes."""

import numpy as np
import pytest

from g2egrid.catalogue import charger, construire
from g2egrid.robustesse import _spectre
from g2egrid.tables import assembler
from g2egrid.temporel import ModeleNonLineaire

#: tensions continues du tableau 17, noeuds positifs (pu), par noeud g2egrid
REFERENCE_COMPLET = {
    201: 1.0000,  # Bm-A1
    221: 1.0108,  # Bm-C1
    202: 1.0100,  # Bb-A1
    211: 1.0002,  # Bb-B1 (= Bb-B1x)
    213: 0.9861,  # Bb-B2
    214: 0.9947,  # Bb-B4
    222: 1.0138,  # Bb-C2
    231: 1.0135,  # Bb-D1
    241: 1.0065,  # Bb-E1 = Bm-E1
    212: 1.0000,  # Bm-B2
    215: 0.9922,  # Bm-B3
    216: 1.0002,  # Bm-B5
    251: 1.0068,  # Bm-F1
}


@pytest.fixture(scope="module")
def complet():
    net = assembler(charger("cigre_b4"))
    net.run_loadflow()
    return net


def _vdc(net):
    return dict(zip(sorted(net.bus_specs_dc), net.lf.Vdc))


def test_dcs1_comme_powerfactory():
    """DCS1 (tableau 13) : 400 MW de C vers A, Bm-C1 a 1,0108 pu."""
    net = assembler(charger("cigre_b4_dcs1"))
    lf = net.run_loadflow()
    assert lf.converged
    assert _vdc(net)[221] == pytest.approx(1.0108, abs=2e-4)


def test_systeme_complet_comme_powerfactory(complet):
    """Les treize tensions continues du tableau 17, a 1e-3 pu pres (7e-4 au
    pire, en Bb-B2 : les pertes des stations ne sont pas celles du modele de
    reference) : les
    deux references (Cb-A1 1,01 ; Cm-B2 1,00), les trois statismes, les
    quatre iles offshore formees par leur station."""
    assert complet.lf.converged
    v = _vdc(complet)
    ecarts = {b: abs(v[b] - ref) for b, ref in REFERENCE_COMPLET.items()}
    assert max(ecarts.values()) < 1e-3, ecarts


def test_statisme_suivi(complet):
    """Chaque station en statisme injecte p_ref + s.(v_dc - 1)."""
    idx = complet._bus_index()
    prog = complet._injections_programmees()
    v = _vdc(complet)
    for b_ac, b_dc, p_ref, sn in (
        (111, 211, 1500e6, 2400e6),
        (113, 213, 1700e6, 1600e6),
        (114, 215, 800e6, 1200e6),
    ):
        k = complet.vsc_statisme[b_ac]
        p = complet.lf.S[idx[b_ac]].real - prog[b_ac][0]
        attendu = p_ref / complet.Sb + (sn / complet.Sb) / k * (v[b_dc] - 1.0)
        assert p == pytest.approx(attendu, abs=1e-9)


def test_formeurs_offshore_portent_leur_ile(complet):
    """Une station Vtheta fournit ce que son ile appelle : en E1, la charge
    de 100 MW et les pertes, rien de plus."""
    idx = complet._bus_index()
    p_e1 = complet.lf.S[idx[141]].real * complet.Sb
    assert -110e6 < -p_e1 < -100e6 or 100e6 < p_e1 < 110e6


def test_jacobien_des_modes_vtheta_et_droop():
    """Le jacobien analytique, regles de chaine du formeur (dP du slack) et
    du statisme (dp/dv_dc) comprises, contre les differences finies."""
    net = assembler(charger("cigre_b4"))
    pb = net._probleme_acdc()
    x0 = pb.depart()
    x = x0 + 0.01 * np.random.default_rng(7).standard_normal(pb.n_x)
    J = pb.jacobien(x)
    Jn = np.zeros_like(J)
    for k in range(pb.n_x):
        h = 1e-7 * max(1.0, abs(x[k]))
        xp, xm = x.copy(), x.copy()
        xp[k] += h
        xm[k] -= h
        Jn[:, k] = -(pb.residu(xp) - pb.residu(xm)) / (2 * h)
    assert np.abs(J - Jn).max() / np.abs(Jn).max() < 1e-6


def test_point_fixe_et_stable():
    """L'equilibre est un point fixe du modele non lineaire, et le reseau,
    avec ses stations reglees en MMC, est stable en petit signal."""
    net = construire("cigre_b4")
    dx, dz = ModeleNonLineaire(net).residu_equilibre()
    assert dx < 1e-5 and dz < 1e-8
    sp = _spectre(charger("cigre_b4"))
    assert float(sp.lam.real.max()) < 0


def test_dcs1_lcl_par_defaut_instable():
    """Ce que le reglage MMC evite : avec le filtre LCL et l'anticipation de
    tension par defaut, la resonance du filtre de Cm-C1 avec le
    transformateur et le reseau C part instable vers 230 Hz."""
    t = charger("cigre_b4_dcs1")
    sp = _spectre(t.copie(GFL_PARAMS={}))
    k = int(np.argmax(sp.lam.real))
    assert sp.lam[k].real > 0
    assert 200 < abs(sp.lam[k].imag) / (2 * np.pi) < 260
    # le mode se partage a parts presque egales entre les quatre etats du
    # filtre (ve, ig) ; lequel domine depend des etats que la reduction
    # d'indice garde aux noeuds ideaux
    assert sp.etats[int(np.argmax(sp.p[k]))].startswith(("GFL121.ve", "GFL121.ig"))


def test_valide_dans_l_atelier():
    """Le validateur de l'atelier admet un formeur offshore sur le slack de
    son ile (mode Vtheta) : aucune erreur sur le reseau du catalogue."""
    analyse = pytest.importorskip("g2egrid_api.analyse")
    erreurs = [
        p for p in analyse.valider(charger("cigre_b4")) if p["gravite"] == "erreur"
    ]
    assert not erreurs, erreurs


def test_balayage_des_statismes():
    """Plus raides (facteur < 1), les statismes reprennent les ecarts de
    tension continue : chaque station s'ecarte davantage de son P_ref, et le
    reseau reste stable sur tout l'intervalle."""
    from g2egrid import continuation as cont

    fab = cont.par_statisme(charger("cigre_b4"))
    c = cont.continuer(fab, [0.5, 1.0, 2.0])
    p = {mu: fab.mesures[mu]["p_mw"] for mu in (0.5, 1.0, 2.0)}
    for b, pref in ((111, 1500.0), (113, 1700.0), (114, 800.0)):
        ecarts = [abs(p[mu][b] - pref) for mu in (0.5, 1.0, 2.0)]
        assert ecarts[0] > ecarts[1] > ecarts[2], (b, ecarts)
    assert not c.bifurcations
    assert max(float(x.real[abs(x) > 1e-6].max()) for x in c.lam) < 0


def test_la_reduction_d_indice_garde_ses_etats_le_long_du_balayage():
    """Aux noeuds ideaux, la reduction d'indice elimine un etat par
    contrainte. Le choix basculait entre GFL113.igq et GFL51.igq de mu = 0,25
    a 1 (pivot decide a l'arrondi pres) : memes valeurs propres, mais la
    continuation ne retrouvait plus ses etats et refusait de suivre."""
    from g2egrid import continuation as cont

    fab = cont.par_statisme(charger("cigre_b4"))
    mus = (0.25, 1.0, 4.0)
    elimines = [cont._construire(fab, mu).build_state_space().elimines for mu in mus]
    assert elimines[0] and elimines[0] == elimines[1] == elimines[2]
    assert len(cont.continuer(fab, mus).mu) == 3


def test_balayage_sans_statisme_refuse():
    from g2egrid import continuation as cont

    with pytest.raises(ValueError, match="aucune station en statisme"):
        cont.par_statisme(charger("mtdc3t"))
