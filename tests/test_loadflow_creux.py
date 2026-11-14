"""Load flow AC creux (FR-027, etape 2) : equivalence avec la boucle d'origine.

Le Newton polaire de ``run_loadflow`` travaille sur un Y-bus CSC et un jacobien
creux (formules ``dSbus_dV`` de MATPOWER), factorise par SuperLU. L'ORACLE est
la double boucle Python d'origine, copiee ci-dessous TELLE QUELLE (elle a ete
retiree du code par FR-027 etape 1 puis etape 2, parce qu'elle coutait 13,6 s
par load flow sur ieee39_mixte). Elle sert de verite de reference a trois
niveaux :

* le jacobien, terme a terme, sur les 23 reseaux de
  ``reference_loadflow.json``, a la solution ET en un point quelconque ;
* le Newton complet (mêmes iterations, tensions a 1e-12 pu) ;
* le Y-bus.

Les filets de non-regression de ``test_loadflow_acdc.py`` (reference stockee a
1e-11 pu, memes iterations, temoins de mutation) restent le garde-fou du
resultat ; ce fichier-ci garde la METHODE.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from test_loadflow_acdc import ALTERNATIFS

from g2egrid.catalogue import charger
from g2egrid.loadflow import (
    PQ,
    PV,
    SLACK,
    build_ybus,
    jacobien_polaire,
    jacobien_polaire_dense,
    run_loadflow,
)
from g2egrid.tables import assembler

#: Ecart admis entre le jacobien creux et celui de la boucle d'origine. Les
#: deux ecrivent les memes derivees avec un ordre d'operations different : la
#: difference est de l'arrondi, relatif a la plus grande entree (|J| atteint
#: 1e3 a 1e4 pu sur les reseaux a branches courtes).
TOL_J = 1e-13


# ----------------------------------------------------------------------
# l'oracle : la boucle d'origine
# ----------------------------------------------------------------------
def jacobien_boucle(Y, V, th, P, Q, idx_th, idx_v):
    """La double boucle de ``run_loadflow`` d'avant FR-027, sans changement."""
    G, B = Y.real, Y.imag
    nt, nv = len(idx_th), len(idx_v)
    J = np.zeros((nt + nv, nt + nv))
    cos = np.cos(th[:, None] - th[None, :])
    sin = np.sin(th[:, None] - th[None, :])
    for a, i in enumerate(idx_th):
        for b, j in enumerate(idx_th):
            if i == j:
                J[a, b] = -Q[i] - B[i, i] * V[i] ** 2
            else:
                J[a, b] = V[i] * V[j] * (G[i, j] * sin[i, j] - B[i, j] * cos[i, j])
        for b, j in enumerate(idx_v):
            if i == j:
                J[a, nt + b] = P[i] / V[i] + G[i, i] * V[i]
            else:
                J[a, nt + b] = V[i] * (G[i, j] * cos[i, j] + B[i, j] * sin[i, j])
    for a, i in enumerate(idx_v):
        for b, j in enumerate(idx_th):
            if i == j:
                J[nt + a, b] = P[i] - G[i, i] * V[i] ** 2
            else:
                J[nt + a, b] = (
                    -V[i] * V[j] * (G[i, j] * cos[i, j] + B[i, j] * sin[i, j])
                )
        for b, j in enumerate(idx_v):
            if i == j:
                J[nt + a, nt + b] = Q[i] / V[i] - B[i, i] * V[i]
            else:
                J[nt + a, nt + b] = V[i] * (G[i, j] * sin[i, j] - B[i, j] * cos[i, j])
    return J


def newton_boucle(Y, btype, Vsp, theta_sp, Psp, Qsp, tol=1e-10, max_iter=40):
    """Le Newton d'origine (jacobien en boucle, ``np.linalg.solve`` dense)."""
    V = np.array(Vsp, dtype=float).copy()
    th = np.array(theta_sp, dtype=float).copy()
    pv = np.where(btype == PV)[0]
    pq = np.where(btype == PQ)[0]
    idx_th = np.concatenate((pv, pq))
    idx_v = pq
    nt = len(idx_th)
    it = 0
    for it in range(1, max_iter + 1):
        Vc = V * np.exp(1j * th)
        S = Vc * np.conj(Y @ Vc)
        P, Q = S.real, S.imag
        mis = np.concatenate((Psp[idx_th] - P[idx_th], Qsp[idx_v] - Q[idx_v]))
        if np.max(np.abs(mis)) < tol:
            break
        J = jacobien_boucle(Y, V, th, P, Q, idx_th, idx_v)
        dx = np.linalg.solve(J, mis)
        th[idx_th] += dx[:nt]
        V[idx_v] += dx[nt:]
    return V * np.exp(1j * th), it


def donnees(nom: str):
    """``(Y dense, btype, Vsp, theta_sp, Psp, Qsp)`` d'un reseau du catalogue."""
    net = assembler(charger(nom))
    idx = net._bus_index()
    n = len(idx)
    btype = np.zeros(n, dtype=int)
    Vsp, th = np.ones(n), np.zeros(n)
    Psp, Qsp = np.zeros(n), np.zeros(n)
    for b, i in idx.items():
        s = net.bus_specs[b]
        btype[i], Vsp[i], th[i] = s.btype, s.V, s.theta
        Psp[i], Qsp[i] = s.Pg - s.Pl, s.Qg - s.Ql
    return net.build_ybus(), btype, Vsp, th, Psp, Qsp


# ----------------------------------------------------------------------
# le jacobien
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nom", ALTERNATIFS)
def test_jacobien_creux_egale_la_boucle(nom):
    """``max|J_creux - J_boucle| <= 1e-13`` (relatif a max|J|), a la solution
    et en un point perturbe (angles et modules quelconques : la diagonale ne
    se reduit plus a la loi des noeuds)."""
    Y, btype, Vsp, th0, Psp, Qsp = donnees(nom)
    pv = np.where(btype == PV)[0]
    pq = np.where(btype == PQ)[0]
    idx_th, idx_v = np.concatenate((pv, pq)), pq
    lf = run_loadflow(Y, btype, Vsp, th0, Psp, Qsp)
    rng = np.random.default_rng(7)
    points = [
        (np.abs(lf.V), np.angle(lf.V)),
        (
            np.abs(lf.V) * (1 + 0.05 * rng.standard_normal(Y.shape[0])),
            np.angle(lf.V) + 0.1 * rng.standard_normal(Y.shape[0]),
        ),
    ]
    for V, th in points:
        Vc = V * np.exp(1j * th)
        S = Vc * np.conj(Y @ Vc)
        Jb = jacobien_boucle(Y, V, th, S.real, S.imag, idx_th, idx_v)
        Jc = jacobien_polaire(sp.csc_array(Y), Vc, idx_th, idx_v)
        assert sp.issparse(Jc) and Jc.shape == Jb.shape
        Jd = jacobien_polaire_dense(Y, V, th, S.real, S.imag, idx_th, idx_v)
        tol = TOL_J * max(1.0, np.abs(Jb).max())
        assert np.abs(Jc.toarray() - Jb).max() <= tol, nom
        assert np.abs(Jd - Jb).max() <= tol, nom


def test_oracle_voit_un_jacobien_faux():
    """TEMOIN : le test d'equivalence echoue si une derivee est fausse.

    Le signe de ``dQ/dV`` est inverse sur un coefficient : l'ecart depasse
    d'un facteur 1e10 la tolerance.
    """
    Y, btype, Vsp, th0, Psp, Qsp = donnees("ieee14")
    pq = np.where(btype == PQ)[0]
    pv = np.where(btype == PV)[0]
    idx_th, idx_v = np.concatenate((pv, pq)), pq
    lf = run_loadflow(Y, btype, Vsp, th0, Psp, Qsp)
    V, th = np.abs(lf.V), np.angle(lf.V)
    S = lf.V * np.conj(Y @ lf.V)
    Jb = jacobien_boucle(Y, V, th, S.real, S.imag, idx_th, idx_v)
    Jc = jacobien_polaire(sp.csc_array(Y), lf.V, idx_th, idx_v).toarray()
    Jc[-1, -1] *= -1.0
    assert np.abs(Jc - Jb).max() > 1e10 * TOL_J


# ----------------------------------------------------------------------
# le Newton complet
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nom", ALTERNATIFS)
def test_newton_creux_egale_newton_boucle(nom):
    """Memes iterations, tensions a 1e-12 pu de la boucle d'origine."""
    Y, btype, Vsp, th0, Psp, Qsp = donnees(nom)
    V_ref, it_ref = newton_boucle(Y, btype, Vsp, th0, Psp, Qsp)
    for methode in ("creux", "dense", "auto"):
        lf = run_loadflow(Y, btype, Vsp, th0, Psp, Qsp, methode=methode)
        assert lf.converged and lf.iterations == it_ref, (nom, methode)
        assert np.abs(lf.V - V_ref).max() <= 1e-12, (nom, methode)
        # Y creux ou dense : meme chemin, meme resultat
        lf2 = run_loadflow(sp.csc_array(Y), btype, Vsp, th0, Psp, Qsp, methode=methode)
        assert np.array_equal(lf.V, lf2.V), (nom, methode)


@pytest.mark.parametrize("nom", ["ieee57", "wscc9_cigre"])
def test_converge_au_dernier_pas_de_la_boucle(nom):
    """Si le dernier pas permis tombe sur la solution, ``converged`` est vrai.

    Avant FR-027 la boucle sortait sans reevaluer le mismatch : un Newton qui
    atteignait la solution a son ``max_iter``-ieme pas restait declare non
    converge (mismatch perime).
    """
    Y, btype, Vsp, th0, Psp, Qsp = donnees(nom)
    ref = run_loadflow(Y, btype, Vsp, th0, Psp, Qsp)
    pas = ref.iterations - 1  # nombre de pas pour descendre sous tol
    juste = run_loadflow(Y, btype, Vsp, th0, Psp, Qsp, max_iter=pas)
    assert juste.converged and juste.mismatch < 1e-10
    assert np.abs(juste.V - ref.V).max() <= 1e-13
    trop_court = run_loadflow(Y, btype, Vsp, th0, Psp, Qsp, max_iter=pas - 1)
    assert not trop_court.converged and trop_court.iterations == pas - 1


def test_jacobien_singulier_leve_linalgerror():
    """Un bus PQ sans branche : jacobien singulier, ``LinAlgError`` comme avant
    (``np.linalg.solve``), pas un ``RuntimeError`` de SuperLU."""
    Y = build_ybus(3, [(0, 1, 0.01, 0.1, 1.0)])  # le bus 2 est isole
    btype = np.array([1, 3, 3])
    with pytest.raises(np.linalg.LinAlgError):
        run_loadflow(
            Y, btype, np.ones(3), np.zeros(3), np.array([0, -0.1, -0.1]), np.zeros(3)
        )


# ----------------------------------------------------------------------
# polissage et pas amorti
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nom", ALTERNATIFS)
def test_polissage_ne_change_pas_les_iterations(nom):
    """Le polissage est un pas de plus, hors du compte : memes iterations. V
    bouge de moins de 5e-11 pu (1,7e-11 au pire, wscc9_cigre : c'est pourquoi
    il n'est pas le defaut, le filet de FR-001 tolere 1e-11) et le mismatch ne
    se degrade pas au-dela du plancher d'arrondi."""
    Y, btype, Vsp, th0, Psp, Qsp = donnees(nom)
    a = run_loadflow(Y, btype, Vsp, th0, Psp, Qsp)
    b = run_loadflow(Y, btype, Vsp, th0, Psp, Qsp, polissage=True)
    assert a.iterations == b.iterations and b.converged
    assert np.abs(a.V - b.V).max() < 5e-11
    assert b.mismatch <= max(a.mismatch, 2e-12)


@pytest.mark.parametrize("nom", ALTERNATIFS)
def test_pas_amorti_sans_effet_quand_le_pas_plein_descend(nom):
    """``amorti=True`` ne bouge ni les iterations ni les tensions sur le
    catalogue : le pas plein y reduit toujours le mismatch en norme 2."""
    Y, btype, Vsp, th0, Psp, Qsp = donnees(nom)
    a = run_loadflow(Y, btype, Vsp, th0, Psp, Qsp)
    b = run_loadflow(Y, btype, Vsp, th0, Psp, Qsp, amorti=True)
    assert a.iterations == b.iterations
    assert np.array_equal(a.V, b.V)


#: (reseau, ecart d'angle, voie). Les departs lointains sont chaotiques : un
#: cas ou le Newton plein diverge sur une voie peut converger sur l'autre (nordic44
#: 1,5 rad converge en dense, a la 33e iteration) ; chaque voie a donc sa liste.
CAS_DEPART_LOINTAIN = [
    (nom, a, m)
    for m in ("creux", "dense")
    for nom, a in [
        ("ieee14", 1.5),
        ("ieee39", 1.0),
        ("ieee33", 0.3),
        ("kundur2area", 1.0),
    ]
] + [("nordic44", 1.5, "creux")]


@pytest.mark.parametrize(("nom", "ecart_angle", "methode"), CAS_DEPART_LOINTAIN)
def test_pas_amorti_rattrape_un_depart_lointain(nom, ecart_angle, methode):
    """Un depart a angles etales (rampe de 0,3 a 1,5 rad) fait diverger le
    Newton plein (non converge en 40 iterations) ; le pas amorti converge, vers
    la meme solution que le depart plat."""
    Y, btype, Vsp, th0, Psp, Qsp = donnees(nom)
    depart = th0 + ecart_angle * np.linspace(0.0, 1.0, th0.size)
    depart[btype == SLACK] = th0[btype == SLACK]  # l'angle du slack est une donnee
    assert not run_loadflow(Y, btype, Vsp, depart, Psp, Qsp, methode=methode).converged
    lf = run_loadflow(Y, btype, Vsp, depart, Psp, Qsp, amorti=True, methode=methode)
    plat = run_loadflow(Y, btype, Vsp, th0, Psp, Qsp)
    assert lf.converged
    assert np.abs(lf.V - plat.V).max() < 1e-9


# ----------------------------------------------------------------------
# le Y-bus : un seul constructeur, creux
# ----------------------------------------------------------------------
def ybus_boucle(n_bus, branches):
    """Le ``loadflow.build_ybus`` dense d'origine (code mort jusqu'ici)."""
    Y = np.zeros((n_bus, n_bus), dtype=complex)
    for f, t, R, X, n in branches:
        y = 1.0 / complex(R, X)
        Y[f, f] += y / (n * n)
        Y[f, t] -= y / n
        Y[t, f] -= y / n
        Y[t, t] += y
    return Y


def test_build_ybus_creux_egale_la_boucle():
    """Lignes en parallele, transformateur a prise, shunts doubles."""
    branches = [
        (0, 1, 0.01, 0.1, 1.0),
        (0, 1, 0.02, 0.2, 1.0),  # parallele
        (1, 2, 0.0, 0.05, 1.05),  # transformateur
        (2, 3, 0.03, 0.12, 1.0),
    ]
    Y = build_ybus(4, branches, shunts=[(1, 0.2j), (1, 0.01 + 0.0j), (3, -0.5j)])
    assert sp.issparse(Y) and Y.format == "csc"
    attendu = ybus_boucle(4, branches)
    attendu[1, 1] += 0.2j + 0.01
    attendu[3, 3] += -0.5j
    assert np.abs(Y.toarray() - attendu).max() <= 1e-15
    with pytest.raises(ValueError):
        build_ybus(2, [(0, 1, 0.0, 0.0, 1.0)])


@pytest.mark.parametrize("nom", ["smib", "wscc9", "ieee39", "nordic44"])
def test_network_build_ybus_est_la_version_dense_du_creux(nom):
    net = assembler(charger(nom))
    Yc = net.build_ybus(creux=True)
    Yd = net.build_ybus()
    assert sp.issparse(Yc) and isinstance(Yd, np.ndarray)
    assert np.array_equal(Yc.toarray(), Yd)
    # les branches seules, contre l'oracle dense
    idx = net._bus_index()
    branches = [
        (
            idx[e.buses[0]],
            idx[e.buses[1]],
            e.pdict["R"],
            e.pdict["X"],
            e.pdict.get("n", 1.0),
        )
        for e in net.elements
        if len(e.ports) == 2
    ]
    sans_shunt = build_ybus(len(idx), branches).toarray()
    assert np.abs(sans_shunt - ybus_boucle(len(idx), branches)).max() <= 1e-12


# ----------------------------------------------------------------------
# temoins de mutation de la voie creuse
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nom", ["wscc9", "ieee14_ieee33"])
def test_le_filet_voit_un_jacobien_creux_degrade(nom, monkeypatch):
    """TEMOIN, voie creuse : un jacobien faux de 1 % ralentit le Newton, et le
    filet de ``test_loadflow_acdc.py`` le voit par les iterations.

    Le temoin d'origine (``test_le_filet_voit_un_jacobien_degrade``) degrade
    ``np.linalg.solve``, donc la voie DENSE (< ``SEUIL_CREUX`` bus). Celui-ci
    degrade ``splu``, la voie creuse, forcee sur ``wscc9``.
    """
    from test_loadflow_acdc import _reference, ecarts_alternatifs

    import g2egrid.loadflow as lfmod

    factoriser = lfmod.splu
    monkeypatch.setattr(lfmod, "splu", lambda J: factoriser(1.01 * J))
    net = assembler(charger(nom))
    attendu = _reference(nom)
    lf = net.run_loadflow(methode="creux")
    ecarts = ecarts_alternatifs(lf, attendu)
    assert ecarts == [
        f"{lf.iterations} iterations au lieu de {attendu['iterations']}"
    ], ecarts
    assert lf.iterations > attendu["iterations"]


@pytest.mark.parametrize("nom", ALTERNATIFS)
def test_les_deux_voies_contre_reference_stockee(nom):
    """Chaque voie, FORCEE, donne la reference de FR-001 sur les 23 reseaux :
    ``max|dV| <= 1e-11`` pu et memes iterations. Le filet de
    ``test_loadflow_acdc.py`` n'emprunte, lui, que la voie du seuil
    (``SEUIL_CREUX``)."""
    from test_loadflow_acdc import _reference, ecarts_alternatifs

    net = assembler(charger(nom))
    for methode in ("creux", "dense"):
        lf = net.run_loadflow(methode=methode)
        assert not ecarts_alternatifs(lf, _reference(nom)), (nom, methode)


def test_methode_inconnue():
    Y, btype, Vsp, th0, Psp, Qsp = donnees("smib")
    with pytest.raises(ValueError, match="methode"):
        run_loadflow(Y, btype, Vsp, th0, Psp, Qsp, methode="gradient")
