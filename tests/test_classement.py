"""Classement des instabilites par type de stabilite (g2egrid.classement)."""

from __future__ import annotations

import numpy as np

from g2egrid.classement import CATEGORIE_DE, classer, decider, groupe, part_commune


def test_groupes():
    assert groupe("SM31.Dwr") == "SM_rotor"
    assert groupe("SM31.Pm") == "SM_gouverneur"
    assert groupe("SM31.efd") == "SM_excitation"
    assert groupe("SM31.phi_d") == "SM_stator"
    assert groupe("GFM32.M_CLd") == "GFM_courant"
    assert groupe("GFM32.theta") == "GFM_synchronisation"
    assert groupe("GFL33.M_pll") == "GFL_pll"
    assert groupe("GFL33.M_q") == "GFL_puissance_reactive"
    assert groupe("L4_6_9.i_d") == "reseau"
    assert groupe("bus5.vd") == "reseau"
    assert groupe("Load5.dq_f") == "charge"
    assert set(CATEGORIE_DE.values()) <= {
        "angle",
        "frequence",
        "tension",
        "convertisseur_lent",
        "convertisseur_rapide",
        "resonance",
    }


def test_categorie_dominante_et_seconde():
    etats = ["SM1.Dwr", "SM1.theta", "SM1.efd", "L1.i_d"]
    c = classer(np.array([0.3, 0.2, 0.4, 0.1]), etats)
    assert c["categorie"] == "angle" and c["seconde"] == "tension"


def test_mode_commun_est_de_la_frequence():
    etats = ["SM1.theta", "SM2.theta", "SM1.Dwr", "SM2.Dwr"]
    M = {"SM1": 5.0, "SM2": 3.0}
    p = np.array([0.25, 0.25, 0.25, 0.25])
    ensemble = np.array([1.0, 1.0, 0.1, 0.1], complex)
    oppose = np.array([1.0, -5 / 3, 0.1, -0.1], complex)
    assert part_commune(ensemble, etats, M) > 0.99
    assert part_commune(oppose, etats, M) < 0.01
    assert (
        classer(p, etats, commun=part_commune(ensemble, etats, M))["categorie"]
        == "frequence"
    )
    assert (
        classer(p, etats, commun=part_commune(oppose, etats, M))["categorie"] == "angle"
    )


def test_convertisseur_lent_ou_rapide_selon_la_frequence():
    """§B5 de la classification 2021 : la FREQUENCE du phenomene separe
    l'interaction lente (< 10 Hz) de la rapide, pas la boucle qui y
    participe -- une PLL peut porter une interaction rapide."""
    etats = ["GFL3.M_pll", "GFL3.theta_pll", "L1_1_2.i_d"]
    p = np.array([0.4, 0.3, 0.3])
    assert classer(p, etats, f=3.0)["categorie"] == "convertisseur_lent"
    assert classer(p, etats, f=45.0)["categorie"] == "convertisseur_rapide"


def test_boucle_de_puissance_gfm_compte_comme_un_regulateur_de_vitesse():
    """La boucle de puissance du GFM gere la puissance, comme le regulateur
    de vitesse d'une machine : elle releve de la frequence, pas des
    commandes de convertisseur."""
    etats = ["GFM2.theta", "GFM2.pm", "SM1.Dwr"]
    c = classer(np.array([0.4, 0.3, 0.3]), etats, f=0.8)
    assert c["categorie"] == "frequence"
    assert c["seconde"] == "angle"
    assert c["part_commande"] == 0


def test_resonance_electrique_sans_commande_ni_machine():
    etats = ["L1_1_2.i_d", "L1_1_2.i_q", "SM1.phi_d", "GFM2.M_CLd"]
    c = classer(np.array([0.5, 0.3, 0.15, 0.05]), etats, f=120.0)
    assert c["categorie"] == "resonance" and c["sous_categorie"] == "electrique"


def test_angle_local_ou_inter_zones():
    etats = ["SM1.Dwr", "SM1.theta", "SM2.Dwr", "SM2.theta"]
    p = np.full(4, 0.25)
    assert classer(p, etats, f=1.6)["sous_categorie"] == "local"
    assert classer(p, etats, f=0.5)["sous_categorie"] == "inter_zones"


def test_reclasser_depuis_les_groupes_enregistres():
    """La campagne n'enregistre que les parts par groupe : l'arbre doit
    pouvoir s'appliquer sur elles seules."""
    d = decider({"reseau": 0.47, "GFM_courant": 0.33, "GFM_tension": 0.19}, 47.7)
    assert d["categorie"] == "convertisseur_rapide"
    assert d["seconde"] == "resonance"


# ------------------------------------------------- procedure complete
from g2egrid.classement import procedure  # noqa: E402


def test_procedure_commandes_par_la_sensibilite():
    """Les commandes reglent la marge : pilotee par convertisseur, meme si le
    reseau porte le mode."""
    c = procedure(
        {"reseau": 0.7, "GFM_courant": 0.3},
        f=60.0,
        sensibilites={"commandes": 5.0, "angle": 1.0, "reseau": 9.0},
    )
    assert c["classe"] == "convertisseur_rapide"
    assert c["critere"] == "sensibilité"
    assert c["mecanisme"] == "boucle de courant-réseau"


def test_procedure_resonance_et_paire():
    c = procedure(
        {"reseau": 0.8, "GFM_courant": 0.2},
        f=60.0,
        sensibilites={"commandes": 1.0, "angle": 2.0, "reseau": 5.0},
    )
    assert c["classe"] == "resonance"
    assert c["mecanisme"] == "réseau-courant GFM"


def test_procedure_pas_de_resonance_aperiodique():
    """Une derive aperiodique portee par le reseau n'est pas une resonance."""
    c = procedure({"reseau": 0.8, "SM_rotor": 0.2}, f=0.0)
    assert c["classe"] != "resonance"
    assert c["nature"] == "apériodique"


def test_procedure_frequence_et_angle():
    g = {"SM_rotor": 0.7, "SM_excitation": 0.3}
    assert procedure(g, f=0.5, commun=0.9)["classe"] == "frequence"
    a = procedure(g, f=1.5, commun=0.1)
    assert (a["classe"], a["sous_classe"]) == ("angle", "local")


def test_procedure_synchronisme_gfm():
    g = {"GFM_synchronisation": 0.6, "SM_rotor": 0.2, "reseau": 0.2}
    loi = procedure(
        g, f=1.0, commun=0.1, sensibilites={"frequence": 3.0, "reseau": 1.0}
    )
    assert (loi["classe"], loi["sous_classe"], loi["seconde"]) == (
        "angle",
        "synchronisme GFM",
        "convertisseur_lent",
    )
    reseau = procedure(
        g, f=1.0, commun=0.1, sensibilites={"frequence": 1.0, "reseau": 3.0}
    )
    assert reseau["mecanisme"] == "couple synchronisant"
