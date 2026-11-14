"""Tests du module de fatigue : rainflow, courbe S-N, regle de Miner.

Chaque test compare a une REFERENCE CALCULABLE A LA MAIN : l'exemple du comptage
de la norme ASTM E1049-85, un dommage de Miner qu'on sait ecrire en une ligne
pour une sinusoide pure, la pente de Basquin lue entre deux points de la courbe.
Aucune tolerance n'est reglee sur le resultat du code.
"""

import numpy as np
import pytest

from g2egrid.fatigue import (
    Arbre,
    CourbeSN,
    couple_base,
    evaluer_fatigue,
    points_de_retournement,
    rainflow,
    tension_cisaillement,
)


def test_rainflow_exemple_astm_e1049():
    """La serie -2, 1, -3, 5, -1, 3, -4, 4, -2 de la norme donne : un cycle
    complet de plage 4 (moyenne 1) et six demi-cycles de plages 3, 4, 8, 9, 8, 6."""
    c = rainflow([-2, 1, -3, 5, -1, 3, -4, 4, -2])
    complets = [(r, m) for r, m, n in zip(c.plage, c.moyenne, c.compte) if n == 1.0]
    demis = sorted(r for r, n in zip(c.plage, c.compte) if n == 0.5)
    assert complets == [(4.0, 1.0)]
    assert demis == [3.0, 4.0, 6.0, 8.0, 8.0, 9.0]
    # (retournements - 1) / 2 cycles equivalents
    assert c.n == pytest.approx((9 - 1) / 2)


def test_rainflow_sinusoide_pure():
    """Une sinusoide d'amplitude a et de N periodes : N cycles (a un demi pres)
    de plage 2a, tous de meme moyenne."""
    N, a, m0 = 50, 3.0, 10.0
    t = np.linspace(0.0, N, 200 * N + 1)
    c = rainflow(m0 + a * np.sin(2 * np.pi * t + 0.3))
    assert c.n == pytest.approx(N, abs=1.0)
    # les cycles de pleine plage (ex aequo a l'arrondi pres : fermes ou en
    # residu selon le sens du ``<`` flottant) portent toute l'energie
    gros = c.plage > 0.99 * 2 * a
    assert c.compte[gros].sum() == pytest.approx(N, abs=1.0)
    assert np.allclose(c.plage[gros], 2 * a, atol=2e-3)
    assert np.allclose(c.moyenne[gros], m0, atol=2e-3)


def test_hysteresis_ignore_le_bruit():
    """Un bruit de 0,01 (ecart-type) sur un signal de plage 2 n'ajoute aucun cycle
    au-dessus du seuil 0,1 (10 sigma), et en ajoute des centaines sans seuil."""
    rng = np.random.default_rng(1)
    t = np.linspace(0, 10, 4001)
    x = np.sin(2 * np.pi * 0.5 * t) + 0.01 * rng.standard_normal(len(t))
    assert rainflow(x, seuil=0.1).n == pytest.approx(5.5, abs=0.5)
    assert rainflow(x, seuil=0.0).n > 100


def test_serie_constante_sans_cycle():
    c = rainflow(np.full(20, 0.7))
    assert c.n == 0.0
    assert len(points_de_retournement(np.full(20, 0.7))) == 1


def test_courbe_sn_pentes_et_coude():
    """``N`` vaut ``N_D`` a ``S_D``, suit la pente ``m`` au-dessus, ``m2`` en dessous."""
    sn = CourbeSN(S_D=100e6, N_D=1e6, m=5.0, m2=9.0)
    assert sn.cycles_a_rupture(100e6) == pytest.approx(1e6)
    # doubler la contrainte au-dessus du coude divise N par 2^m
    assert sn.cycles_a_rupture(200e6) == pytest.approx(1e6 / 2**5)
    # la diviser par deux sous le coude multiplie N par 2^m2
    assert sn.cycles_a_rupture(50e6) == pytest.approx(1e6 * 2**9)
    # Miner elementaire : aucun dommage sous la limite d'endurance
    assert np.isinf(CourbeSN(100e6, 1e6, 5.0).cycles_a_rupture(50e6))
    assert np.isinf(sn.cycles_a_rupture(0.0))


def test_unites_couple_et_contrainte():
    """Couple de base d'un groupe de 2,5 MVA a 1500 tr/min, et ``16.C/(pi.d^3)``."""
    cb = couple_base(2.5e6, 50.0, 2)
    assert cb == pytest.approx(2.5e6 / (2 * np.pi * 25.0))  # 15 915 N.m
    assert tension_cisaillement(1000.0, 0.1) == pytest.approx(
        16 * 1000 / (np.pi * 1e-3)
    )


def test_miner_sinusoide_dommage_analytique():
    """Un couple sinusoidal de plage donnee, N periodes : ``D = N / N_f(S)`` avec
    ``S`` l'amplitude de contrainte apres Goodman. Meme calcul ecrit a la main."""
    Sn, f0 = 2.5e6, 50.0
    arbre = Arbre(diametre_m=0.11, tau_u=600e6, paires_poles=2)
    sn = CourbeSN(S_D=5e6, N_D=2e6, m=7.0)
    moyen, ampl, N = 0.5, 0.2, 100  # pu
    t = np.linspace(0.0, 100.0, 100 * 400 + 1)
    r = evaluer_fatigue(
        t, moyen + ampl * np.sin(2 * np.pi * t / 1.0), Sn, f0, arbre, sn, seuil_pu=0.0
    )
    cb = couple_base(Sn, f0, 2)
    tau_a = tension_cisaillement(ampl * cb, 0.11)
    tau_m = tension_cisaillement(moyen * cb, 0.11)
    s_eq = tau_a / (1 - tau_m / arbre.tau_u)
    attendu = N / (2e6 * (5e6 / s_eq) ** 7.0)
    assert tau_a > sn.S_D  # on est au-dessus du coude : formule de Basquin simple
    assert r.dommage == pytest.approx(
        attendu, rel=2e-2
    )  # le comptage est a un cycle pres
    # 100 s de sollicitation donnent `attendu` : la vie en heures est 100 / (3600 D)
    assert r.vie_heures == pytest.approx(100.0 / (3600.0 * attendu), rel=2e-2)


def test_dommage_croit_avec_la_puissance_de_l_amplitude():
    """Doubler l'amplitude du couple multiplie le dommage par 2^m (loi de Basquin)."""
    t = np.linspace(0.0, 60.0, 6001)
    sn = CourbeSN(S_D=3e6, N_D=2e6, m=7.0)
    d = [
        evaluer_fatigue(
            t, 0.3 + a * np.sin(2 * np.pi * t), 2.5e6, courbe=sn, seuil_pu=0.0
        ).dommage
        for a in (0.1, 0.2)
    ]
    # meme moyenne pour tous les cycles : le facteur de Goodman se simplifie et le
    # rapport est 2^7, aux extrema echantillonnes pres
    assert d[1] / d[0] == pytest.approx(2.0**7, rel=1e-4)


def test_goodman_aggrave_et_moyenne_negative_inoffensive():
    """Une moyenne positive augmente le dommage ; une moyenne negative ne le
    diminue pas (hypothese prudente)."""
    t = np.linspace(0.0, 30.0, 3001)
    sn = CourbeSN(S_D=3e6, N_D=2e6, m=7.0)

    def dom(moyen):
        return evaluer_fatigue(
            t, moyen + 0.2 * np.sin(2 * np.pi * t), 2.5e6, courbe=sn, seuil_pu=0.0
        ).dommage

    assert dom(0.8) > dom(0.0)
    assert dom(-0.8) == pytest.approx(dom(0.0), rel=1e-3)


def test_entrees_invalides():
    t = np.linspace(0, 1, 11)
    with pytest.raises(ValueError):
        evaluer_fatigue(t, np.ones(10), 1e6)
    with pytest.raises(ValueError):
        evaluer_fatigue(t, np.where(t > 0.5, np.nan, 1.0), 1e6)
    with pytest.raises(ValueError):
        CourbeSN(S_D=-1.0, N_D=1e6, m=5.0)


def test_histogramme_somme_des_dommages():
    t = np.linspace(0.0, 60.0, 6001)
    c = 0.4 + 0.3 * np.sin(2 * np.pi * t / 3) + 0.1 * np.sin(2 * np.pi * t / 0.7)
    r = evaluer_fatigue(t, c, 2.5e6, seuil_pu=1e-3)
    _b, comptes, dom = r.histogramme(10)
    assert comptes.sum() == pytest.approx(r.n_cycles)
    assert dom.sum() == pytest.approx(r.dommage, rel=1e-12)


def test_marge_amplitude_inverse_de_la_loi_de_basquin():
    """La marge ``k`` qui ramene la vie a la cible est l'inverse de la loi de
    Basquin, branche haute ou branche sous le coude. Les deux sont verifiees
    contre la formule fermee."""
    t = np.linspace(0.0, 60.0, 6001)
    sn = CourbeSN(S_D=1e6, N_D=2e6, m=3.0, m2=5.0)
    r = evaluer_fatigue(
        t,
        0.3 + 0.2 * np.sin(2 * np.pi * t),
        2.5e6,
        courbe=sn,
        seuil_pu=0.0,
        goodman=False,
    )
    S = r.amplitude_max
    assert S > sn.S_D  # k = 1 est au-dessus du coude

    # cible courte (3,6 s) : la solution est au-dessus du coude, pente m
    cible_h = 1e-3
    d_cible = r.duree_s / (3600.0 * cible_h)
    k = r.marge_amplitude(cible_h)
    assert k == pytest.approx((d_cible / r.dommage) ** (1.0 / 3.0), rel=1e-6)
    assert k * S > sn.S_D and k > 1.0

    # cible longue (25 ans) : la solution passe sous le coude, pente m2
    d_long = r.duree_s / (3600.0 * 25 * 8760.0)
    k2 = r.marge_amplitude()
    somme = np.sum(r.cycles.compte * r.cycles.amplitude**5.0)
    attendu = sn.S_D * (d_long * sn.N_D / somme) ** (1.0 / 5.0)
    assert k2 * S < sn.S_D
    assert k2 == pytest.approx(attendu, rel=1e-6)
    # a la marge trouvee, la vie vaut la cible
    assert r.dommage_a_l_echelle(k2) * 3600.0 / r.duree_s == pytest.approx(
        1.0 / (25 * 8760.0), rel=1e-6
    )
    # une serie qui use trop vite a une marge < 1
    assert r.marge_amplitude(1e9) < 1.0
