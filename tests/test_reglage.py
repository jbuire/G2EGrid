"""Regler par intention : boucles des convertisseurs, fiches de machines."""

from __future__ import annotations

import pytest

import g2egrid.catalogue as cat
from g2egrid.elements import SM_DEFAULT
from g2egrid.reglage import (
    boucles,
    intention_depuis_gains,
    machine_depuis_standard,
    standard_depuis_machine,
)
from g2egrid.tables import assembler
from g2egrid_api import analyse

#: Kundur, exemple 12.6 (reseau a deux zones) : la fiche des quatre machines.
KUNDUR_12_6 = dict(
    Xd=1.8, Xq=1.7, Xl=0.2, Xd_p=0.3, Xq_p=0.55, Xd_pp=0.25, Xq_pp=0.25,
    Ra=0.0025, Td0_p=8.0, Tq0_p=0.4, Td0_pp=0.03, Tq0_pp=0.05,
)  # fmt: skip


def test_le_jeu_par_defaut_redonne_l_exemple_4_1_de_kundur():
    """Le jeu par defaut est l'exemple 3.1 de Kundur ; ses donnees standard
    sont celles que le livre en tire a l'exemple 4.1 (60 Hz)."""
    s = standard_depuis_machine(SM_DEFAULT, f0=60.0)
    attendu = dict(
        Xd=1.81, Xq=1.76, Xd_p=0.30, Xd_pp=0.23, Xq_p=0.65, Xq_pp=0.25,
        Td0_p=8.07, Td0_pp=0.03, Tq0_p=1.0, Tq0_pp=0.07,
    )  # fmt: skip
    for k, v in attendu.items():
        assert s[k] == pytest.approx(v, rel=1e-2), k


def test_aller_retour_exact():
    s = standard_depuis_machine(SM_DEFAULT, f0=50.0)
    p = machine_depuis_standard(s, f0=50.0)
    for k, v in p.items():
        assert v == pytest.approx(SM_DEFAULT[k], rel=1e-12), k
    assert standard_depuis_machine(
        {**SM_DEFAULT, **machine_depuis_standard(KUNDUR_12_6, 60.0)}, 60.0
    ) == pytest.approx(KUNDUR_12_6, rel=1e-12)


def test_une_fiche_non_physique_est_refusee():
    with pytest.raises(ValueError, match="axe d"):
        machine_depuis_standard({**KUNDUR_12_6, "Xd_pp": 0.35})  # X''d > X'd
    with pytest.raises(ValueError, match="Td0_p"):
        machine_depuis_standard({**KUNDUR_12_6, "Td0_p": 0.0})


def test_une_machine_saisie_par_sa_fiche():
    """Les quatre machines de Kundur, posees dans wscc9 par leur fiche : le
    reseau s'assemble, reste stable, et la fiche se relit."""
    t = cat.charger("wscc9")
    p = machine_depuis_standard(KUNDUR_12_6, f0=t.f0)
    t.PARAMS_PAR_BUS = {int(r[1]): dict(p) for r in t.T_SMs}
    net = assembler(t)
    assert net.build_state_space().is_stable()
    sm = next(e for e in net.elements if e.name == "SM1")
    assert standard_depuis_machine(sm.pdict, t.f0)["Xd_p"] == pytest.approx(0.3)


@pytest.mark.parametrize("famille", ["GFM", "GFL"])
def test_gains_lus_sur_le_modele_et_inverse(famille):
    reglages = (
        {"tau_VL": 0.03, "zeta_VL": 1.0, "tau_CL": 2e-4}
        if famille == "GFM"
        else {"tau_CL": 0.02, "tau_PLL": 0.08, "zeta_PLL": 0.9, "tau_PL": 0.2}
    )
    bs = boucles(famille, reglages, f0=50.0)
    for b in bs:
        if b["cle_tau"] in reglages:
            assert b["tau"] == reglages[b["cle_tau"]]
    gains = {k: v for b in bs for k, v in b["gains"].items()}
    retour = intention_depuis_gains(famille, gains, reglages, f0=50.0)
    for k, v in retour.items():
        attendu = next(
            b["tau"] if b["cle_tau"] == k else b["zeta"]
            for b in bs
            if k in (b["cle_tau"], b["cle_zeta"])
        )
        assert v == pytest.approx(attendu, rel=1e-10), k
    if famille == "GFM":
        (tension,) = [b for b in bs if b["nom"] == "tension"]
        assert tension["wn"] == pytest.approx(100.0)  # 3/(1 x 0,03)


def test_balayage_d_un_temps_de_reponse():
    """tau_VL se balaie comme n'importe quel reglage de famille."""
    d = analyse.balayage(
        cat.charger("wscc9_gfm"),
        {"genre": "famille", "famille": "GFM_PARAMS", "cle": "tau_VL"},
        [0.01, 0.015, 0.03],
        f_max=50.0,
    )
    assert len(d["pas"]) == 3 and all(p.get("modes") for p in d["pas"])


def test_api():
    b = analyse.reglage_boucles("GFM", {"tau_VL": 0.03}, 50.0)["boucles"]
    assert b[0]["tau"] == 0.03
    s = analyse.machine_standard(dict(SM_DEFAULT), 60.0)["standard"]
    assert s["Xd_p"] == pytest.approx(0.30, rel=1e-2)
    p = analyse.machine_standard(KUNDUR_12_6, 60.0, "depuis_standard")["parametres"]
    assert set(p) >= {"Lad", "Lfd", "Rfd", "L1d", "R1d"}
