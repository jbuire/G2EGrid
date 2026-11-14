"""FR-005 : un reseau a noeuds ideaux (``NOEUDS = "ideal"``, le defaut) se
simule en BDF2.

La loi des noeuds d'un noeud ideal peut ne porter que sur des etats
(courants en serie) : contrainte d'indice 2, sans solution algebrique a ``x``
fige. BDF2 resout etats et contraintes ensemble et l'integre ; Radau est
refuse. Le noeud ideal est la limite de la resistance de causalite quand
``R`` tend vers l'infini : les reponses doivent s'accorder avec ``R = 1e5``
a l'ordre de ``1/R`` pres.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2egrid.catalogue import charger
from g2egrid.tables import assembler
from g2egrid.temporel import (
    Evenement,
    Manoeuvre,
    ModeleNonLineaire,
    linearisation,
    mesures,
)

TE = np.linspace(0.0, 0.1, 51)


def _modele(nom, **reglages):
    net = assembler(charger(nom).copie(**reglages))
    net.run_loadflow()
    net.solve_equilibrium()
    return ModeleNonLineaire(net)


@pytest.fixture(scope="module")
def ideal():
    return _modele("wscc9", NOEUDS="ideal")


@pytest.fixture(scope="module")
def resistif():
    return _modele("wscc9", NOEUDS="resistance", R_NOEUDS=1e5)


def _simuler(m, *evenements):
    sc = m.scenario(*evenements)
    tr = m.simuler((0.0, 0.1), uext_fn=sc, methode="BDF2", pas=1e-4, t_eval=TE)
    return tr, mesures(m, tr, sc)


def test_le_modele_connait_ses_noeuds_ideaux(ideal, resistif):
    assert sorted(ideal.ideaux) == ["bus1", "bus2", "bus3"]
    assert resistif.ideaux == []


def test_radau_est_refuse_et_bdf2_choisi(ideal):
    with pytest.raises(ValueError, match="BDF2"):
        ideal.simuler((0.0, 0.01), methode="Radau")
    # sans methode : BDF2, choisi de lui-meme
    tr = ideal.simuler((0.0, 0.002), pas=1e-4)
    assert tr.algebriques and np.isfinite(tr.x).all()


def test_echelon_comme_a_resistance_infinie(ideal, resistif):
    ev = Evenement(0.01, "SM1.p_ref", 0.05)
    tr_i, v_i = _simuler(ideal, ev)
    tr_r, v_r = _simuler(resistif, ev)
    assert tr_i.algebriques and len(tr_i.algebriques) == len(tr_i.t)
    assert not tr_r.algebriques  # rien a garder : l'algebrique se re-resout
    assert tr_i.etats == tr_r.etats
    rep = tr_r.x - tr_r.x[0]
    ecart = np.abs((tr_i.x - tr_i.x[0]) - rep).max()
    assert ecart < 1e-4 * np.abs(rep).max()
    # la tension des noeuds ideaux (bus des machines) se relit
    for b in (1, 2, 3):
        assert np.isfinite(v_i[f"V_bus{b}"]).all()
        assert np.abs(v_i[f"V_bus{b}"] - v_r[f"V_bus{b}"]).max() < 1e-5


@pytest.mark.parametrize(
    "evenements",
    [
        [Manoeuvre(0.02, "ouverture", "L3_5_7")],
        [
            Manoeuvre(0.02, "defaut", "bus5", R=1e-2, X=1e-2),
            Manoeuvre(0.05, "elimination", "bus5"),
        ],
        [Manoeuvre(0.02, "charge", "Load5", dp=0.2, dq=0.1)],
    ],
    ids=["ouverture", "defaut", "charge"],
)
def test_manoeuvres(ideal, resistif, evenements):
    _tr_i, v_i = _simuler(ideal, *evenements)
    _tr_r, v_r = _simuler(resistif, *evenements)
    for b in range(1, 10):
        vi, vr = v_i[f"V_bus{b}"], v_r[f"V_bus{b}"]
        assert np.isfinite(vi).all()
        assert np.abs(vi - vr).max() < 2e-3


def _inductif(**reglages):
    """Bus infini, ligne sans capacite, charge au bus 2 : deux noeuds ideaux
    ou ne se rencontrent que des courants d'etat (indice 2)."""
    from g2egrid.tables import Tables

    t = Tables(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3, 20e3, 1.0, 0.0, 0.0, 0.0, 5e6, 1e6),
        ],
        T_Lines_pu=[(1, 1, 2, 1.0, 0.01, 0.1, 0.0, 20e3, 100e6)],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
        **reglages,
    )
    net = assembler(t)
    net.run_loadflow()
    net.solve_equilibrium()
    return net


def test_l_elimination_a_un_noeud_inductif_est_reinitialisee(monkeypatch):
    """Retirer le defaut, qui porte ~10 pu, viole la loi du noeud 2. La
    reprise (``reinitialiser``) applique l'impulsion de la limite ideale
    avant le premier pas : aucune tension ``Phi/pas`` (~30 pu ici) n'est
    stockee, et la reponse ne depend pas du pas ``h`` de la reprise."""
    import functools

    import g2egrid.integrateur as integrateur
    import g2egrid.temporel as temporel

    m = ModeleNonLineaire(_inductif())
    assert sorted(m.ideaux) == ["bus1", "bus2"]
    sc = m.scenario(
        Manoeuvre(0.02, "defaut", "bus2", R=1e-3, X=1e-2),
        Manoeuvre(0.05, "elimination", "bus2"),
    )
    v = {}
    for h in (1e-8, 1e-9):
        monkeypatch.setattr(
            temporel, "reinitialiser", functools.partial(integrateur.reinitialiser, h=h)
        )
        tr = m.simuler((0.0, 0.08), uext_fn=sc, methode="BDF2", pas=1e-4)
        v[h] = mesures(m, tr, sc)["V_bus2"]
    pendant = v[1e-8][(tr.t > 0.021) & (tr.t < 0.05)]
    assert pendant.max() < 0.2  # le defaut tient le noeud bas
    assert v[1e-8].max() < 1.01 and v[1e-8][-1] > 0.95  # pas d'impulsion stockee
    assert np.abs(v[1e-8] - v[1e-9]).max() < 1e-8


def test_l_impedance_d_un_noeud_inductif_se_lit_sur_le_dae():
    """Un courant injecte au noeud 2 y voit ``Z ~ s.L`` : impropre, sans
    modele d'etat. Lue sur le DAE, elle est la limite de la resistance de
    causalite (ecart en ``|Z|/R``)."""
    from g2egrid.impedance import impedance_bus

    f = np.array([1.0, 50.0, 1000.0])
    z = impedance_bus(_inductif(), 2, f)
    assert abs(z[-1, 0, 0]) > 1.5  # elle croit avec la frequence
    for r in (1e2, 1e4):
        zr = impedance_bus(_inductif(NOEUDS="resistance", R_NOEUDS=r), 2, f)
        ecart = np.abs(z - zr).max() / np.abs(z).max()
        assert ecart < 2 * np.abs(z).max() / r
        assert ecart > 0.2 * np.abs(z).max() / r


def test_linearisation_le_long_refusee(ideal):
    tr, _ = _simuler(ideal, Evenement(0.01, "SM1.p_ref", 0.05))
    with pytest.raises(ValueError, match="ideaux"):
        linearisation(ideal, tr, 0.05)


def test_l_atelier_prend_bdf2():
    from g2egrid_api import analyse

    t = charger("wscc9").copie(NOEUDS="ideal")
    r = analyse.temporel(
        t, "SM1.p_ref", 0.05, duree=0.05, n=11, voie="non_lineaire", pas=1e-4,
        mesures=["V_bus1"],
    )  # fmt: skip
    assert r["integrateur"] == "BDF2"
    assert np.isfinite(r["mesures"]["V_bus1"]).all()
    a = analyse.apres_evenement(
        t, [{"genre": "ouverture", "instant": 0.02, "cible": "L3_5_7"}],
        duree=0.06, n=13, pas=1e-4,
    )  # fmt: skip
    assert {b["bus"] for b in a["bus"]} == set(range(1, 10))
    assert all(b["V0"] is not None for b in a["bus"])
