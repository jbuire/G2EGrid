"""La recherche des modes sensibles aux variations de puissance et au N-1."""

from __future__ import annotations

import numpy as np
import pytest

from g2egrid import modal
from g2egrid.catalogue import charger
from g2egrid.robustesse import _apparier, _spectre, etudier, perturbations, resume
from g2egrid.tables import assembler


@pytest.fixture(scope="module")
def wscc9():
    return charger("wscc9")


@pytest.fixture(scope="module")
def etude(wscc9):
    return etudier(wscc9, zeta_max=30)


def test_les_perturbations_couvrent_unites_charges_et_ouvrables(wscc9):
    libelles = [p.libelle for p, _ in perturbations(wscc9)]
    assert "SM2 P +10%" in libelles and "SM3 P +10%" in libelles
    assert "SM1 P +10%" not in libelles  # le slack n'a pas de consigne
    assert "charge totale +10%" in libelles
    assert "ouverture L6_8_9" in libelles and "ouverture SM2" in libelles


def test_une_hausse_de_charge_est_reprise_par_les_unites(wscc9):
    (_, t2), *_ = [
        (p, t) for p, t in perturbations(wscc9) if p.libelle == "charge bus 5 +10%"
    ]
    avant = {int(r[0]): r for r in wscc9.Y_DER}
    apres = {int(r[0]): r for r in t2.Y_DER}
    assert apres[5][7] == pytest.approx(1.1 * avant[5][7])
    hausse = sum(apres[b][5] - avant[b][5] for b in (2, 3))
    sn = {int(r[1]): float(r[2]) for r in wscc9.T_SMs}
    # part des unites non slack, au prorata des puissances nominales
    attendu = 0.1 * avant[5][7] * (sn[2] + sn[3]) / sum(sn.values())
    assert hausse == pytest.approx(attendu)


def test_ouvrir_un_transfo_de_groupe_est_un_ilotage(wscc9):
    genres = {p.libelle: p.genre for p, _ in perturbations(wscc9)}
    assert genres["ouverture T3_3_9"] == "ilotage"
    assert genres["ouverture L6_8_9"] == "deconnexion"


def test_une_variante_identique_se_retrouve_exactement(wscc9):
    base = _spectre(wscc9)
    garde = np.arange(len(base.lam))
    paires = _apparier(base, _spectre(wscc9.copie()), garde)
    assert all(i == j for i, j in paires.items())
    assert len(paires) == len(base.lam)


def test_le_mode_d_une_unite_ouverte_disparait(etude):
    # le mode a 2,2 Hz est porte a 72 % par SM3
    m = min(etude.modes, key=lambda m: abs(m.f - 2.2))
    assert m.elements[0][0] == "SM3"
    assert "ouverture SM3" in m.disparitions


def test_les_ecarts_sont_ceux_des_spectres(wscc9, etude):
    m = etude.modes[0]
    e = next(e for e in m.suivis("deconnexion"))
    t2 = next(t for p, t in perturbations(wscc9) if p.libelle == e.perturbation.libelle)
    net = assembler(t2)
    net.solve_equilibrium()
    lam = net.build_state_space().eig()
    assert np.min(np.abs(lam - e.lam)) < 1e-6 * abs(e.lam)
    assert e.df == pytest.approx(float(modal.frequence(e.lam)[0]) - m.f)
    assert e.score == pytest.approx(abs(e.lam - m.lam) / abs(m.lam))


def test_le_classement_et_le_resume(etude):
    scores = [m.score() for m in etude.modes]
    assert scores == sorted(scores, reverse=True)
    assert all(m.zeta_min <= m.zeta for m in etude.modes)
    texte = resume(etude)
    assert "WSCC9" in texte and "ilotages" in texte


def test_la_puissance_bouge_peu_les_modes_electromecaniques(etude):
    """Le constat de l'etude : +10 % de puissance deplace les modes de
    moins de 3 %, une deconnexion jusqu'a plusieurs dizaines de %."""
    assert max(m.score("puissance") for m in etude.modes) < 0.03
    assert max(m.score("deconnexion") for m in etude.modes) > 0.2


def test_les_gfl_pv_ne_reprennent_pas_la_charge():
    """Un GFL est un PV : a sa puissance MPPT, il ne reprend pas d'ecart ;
    le GFM (a batterie) et les machines, si."""
    t = charger("wscc9_cigre")
    (t2,) = [v for p, v in perturbations(t) if p.libelle == "charge totale +10%"]
    avant = {int(r[0]): float(r[5]) for r in t.Y_DER}
    apres = {int(r[0]): float(r[5]) for r in t2.Y_DER}
    gfl = [int(r[1]) for r in t.T_GFL1s]
    gfm = [int(r[1]) for r in t.T_GFM1s]
    assert gfl and all(apres[b] == avant[b] for b in gfl)
    assert gfm and all(apres[b] > avant[b] for b in gfm)


@pytest.mark.parametrize(
    ("erreur", "attendu"),
    [
        (RuntimeError("load flow non converge"), "echec"),
        (np.linalg.LinAlgError("jacobien singulier"), "echec"),
        (ValueError("reseau sans slack"), "echec"),
        (TypeError("erreur de programmation"), "remonte"),
        (KeyError("bus"), "remonte"),
    ],
)
def test_seule_l_absence_de_point_de_fonctionnement_est_un_echec(
    monkeypatch, wscc9, erreur, attendu
):
    """``etudier`` range une variante sans point de fonctionnement dans
    ``echecs`` ; une erreur de programmation (TypeError, KeyError...) remonte
    au lieu d'etre avalee par un ``except Exception``."""
    from g2egrid import robustesse

    vrai = robustesse._spectre
    appels = []

    def spectre(t, vecteurs=False):
        appels.append(1)
        if len(appels) == 1:  # le cas de base doit passer
            return vrai(t, vecteurs)
        raise erreur

    monkeypatch.setattr(robustesse, "_spectre", spectre)
    if attendu == "remonte":
        with pytest.raises(type(erreur)):
            etudier(wscc9, zeta_max=30)
    else:
        e = etudier(wscc9, zeta_max=30)
        assert len(e.echecs) == e.n_variantes > 0
