"""Validation des elements de convertisseur et de machine qui n'avaient pas de
test : SMTR, GFL2, GFM2 et GFM4 (equilibre exact, jacobiens contre differences
finies), equivalence SMTR = SM + transformateur, et comparaison a la
REFERENCE MATLAB (l'admittance dq de ``bode_equipements_matlab.csv``).

Le CSV est la seule sortie MATLAB exploitable ici (MATLAB n'est pas
disponible) : l'admittance ``Y(jw)`` 2x2 de chaque equipement face a un reseau
infini, 0,01 a 1000 Hz. Les seuils sont justifies a cote de chaque test.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import brentq

from g2egrid import Network
from g2egrid.elements import (
    GFL1,
    GFL2,
    GFM1,
    GFM2,
    GFM3,
    GFM4,
    IB,
    Line,
    SynchronousMachine,
    SynchronousMachineTR,
    Transformer,
)
from g2egrid.frequence import admittance_dq
from test_limitation_impedance_virtuelle import _jacobiens_fd

SB, UN, F0 = 1e6, 20e3, 50.0
CSV = (
    pathlib.Path(__file__).resolve().parents[2]
    / "matlab/G2Elin_V2/etudes/bode_equipements_matlab.csv"
)


def banc(cls, P=0.5, Q=0.0, **kw):
    """L'equipement seul face a un reseau infini, V = 1 pu : le banc du CSV."""
    net = Network(f0=F0, Sb=SB)
    net.add_bus(1, UN, btype=3, Pg=P, Qg=Q, C=1e-6)
    net.add_bus(2, UN, btype=1, V=1.0, C=1e-6)
    net.add(cls("EQ", 1, Sn=SB, Sb=SB, wb=net.wb, **kw))
    net.add(Line("L", 1, 2, R=1e-6, X=1e-5, wb=net.wb))
    net.add(IB("GRID", 2, R=1e-6, X=1e-5, wb=net.wb, Sn=1e4 * SB, Sb=SB))
    net.run_loadflow()
    net.solve_equilibrium(tol=1e-13)
    return next(e for e in net.elements if e.name == "EQ")


SMTR_KW = dict(Rt=0.0002, Lt=0.024)
# (nom, classe, reglages). SMTR en init physique (equilibre exact) ; l'init de
# V1 (``symSG``) n'est PAS un equilibre (|f| = 114 sur la machine, algebriques
# inconsistantes a 4e-6) : elle ne sert qu'a la comparaison au CSV MATLAB.
CAS = [
    ("SMTR", SynchronousMachineTR, SMTR_KW),
    ("GFL2", GFL2, {}),
    ("GFM2", GFM2, {}),
    ("GFM4", GFM4, {}),
]


# ----------------------------------------------------------------------
# equilibre et jacobiens
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nom,cls,kw", CAS, ids=[c[0] for c in CAS])
def test_equilibre_exact(nom, cls, kw):
    """|f| petit a l'equilibre. Les equations a constante de temps courte
    (didc/dt a Tdc = 1e-4 s) amplifient l'arrondi de 1e4 : seuil 1e-8, soit un
    arrondi de 1e-12 sur la grandeur."""
    e = banc(cls, **kw)
    f = np.asarray(e.compiled().f(e.x0, e.z0, e.u0, e.p))
    assert np.abs(f).max() < 1e-8


@pytest.mark.parametrize("nom,cls,kw", CAS, ids=[c[0] for c in CAS])
def test_jacobiens_contre_differences_finies(nom, cls, kw):
    """A, B, C, D symboliques contre differences finies centrees, erreur
    relative a la plus grande entree : 1e-8. h = 1e-7 : l'erreur de troncature
    est en h^2 (SMTR : 3,5e-8 a h = 1e-6, 1,7e-9 a h = 1e-7), l'arrondi ne
    domine qu'a h = 1e-8."""
    assert _jacobiens_fd(banc(cls, **kw), h=1e-7) < 1e-8


# ----------------------------------------------------------------------
# SMTR == SM + transformateur
# ----------------------------------------------------------------------
def _reseau_machine(avec_tr: bool, V3: float = 1.0):
    net = Network(f0=F0, Sb=SB)
    if avec_tr:
        net.add_bus(1, UN, btype=2, V=1.0, Pg=0.8, C=1e-3)
    else:
        net.add_bus(1, UN, btype=3, C=1e-3)
        net.add_bus(3, UN, btype=2, V=V3, Pg=0.8, C=1e-9)
    net.add_bus(2, UN, btype=1, V=1.0, C=1e-3)
    if avec_tr:
        net.add(SynchronousMachineTR("SM1", 1, Sn=SB, Sb=SB, wb=net.wb, **SMTR_KW))
    else:
        net.add(SynchronousMachine("SM1", 3, Sn=SB, Sb=SB, wb=net.wb))
        net.add(Transformer("TR1", 3, 1, R=0.0002, X=0.024, wb=net.wb))
    net.add(Line("L1", 1, 2, R=0.01, X=0.2, wb=net.wb))
    net.add(IB("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
    net.run_loadflow()
    return net


def test_smtr_equivaut_a_sm_plus_transformateur():
    """SMTR est la machine et son transformateur dans une seule classe : meme
    modele, a la tension de la machine pres. On cale la tension V3 de la machine
    nue pour que le bus reseau soit a 1 pu dans les deux cas, puis on compare
    les modes electromecaniques (0,1 a 20 Hz) : ils coincident a 1e-3."""
    V3 = brentq(
        lambda v: abs(_reseau_machine(False, v).lf.V[0]) - 1.0, 0.9, 1.2, xtol=1e-12
    )
    a, b = _reseau_machine(True), _reseau_machine(False, V3)
    for n in (a, b):
        n.solve_equilibrium()

    def modes(net):
        ev = net.build_state_space().eig()
        ev = ev[
            (np.abs(ev.imag) > 2 * np.pi * 0.1) & (np.abs(ev.imag) < 2 * np.pi * 20)
        ]
        return np.sort_complex(ev)

    ma, mb = modes(a), modes(b)
    assert len(ma) == len(mb) >= 4
    assert np.abs(ma - mb).max() < 1e-3
    assert a.build_state_space().is_stable() and b.build_state_space().is_stable()


# ----------------------------------------------------------------------
# reference MATLAB : bode_equipements_matlab.csv
# ----------------------------------------------------------------------
pytestmark_csv = pytest.mark.skipif(not CSV.exists(), reason="CSV MATLAB absent")


def _y_matlab(T, nom):
    Y = np.empty((len(T), 2, 2), complex)
    for i in range(2):
        for j in range(2):
            Y[:, i, j] = (
                T[f"{nom}_Y{i + 1}{j + 1}_re"] + 1j * T[f"{nom}_Y{i + 1}{j + 1}_im"]
            )
    return Y


def _ecart(el, T, nom):
    """Ecart relatif max, frequence par frequence, a la plus grande entree de
    Y MATLAB."""
    Ym = _y_matlab(T, nom)
    Yp = admittance_dq(el, T["f_Hz"].to_numpy())
    return float((np.abs(Yp - Ym).max(axis=(1, 2)) / np.abs(Ym).max(axis=(1, 2))).max())


#: Seuil sur l'ecart relatif max. Le CSV est ecrit avec ~6 chiffres
#: significatifs et les entrees de Y s'etendent sur plusieurs decades : les
#: ecarts constates sont de 6e-6 a 1e-5 (SM, GFM1-4, GFL1-2). 5e-5 laisse un
#: facteur 5 sans cacher un vrai ecart de modele (un terme oublie donne 1e-2).
SEUIL = 5e-5
#: SMTR: l'equilibre de V1 (``symSG``) est reproduit a 2,4e-4 ; l'init
#: PHYSIQUE donne 0,36 en tres basse frequence : ECART VOULU, l'equilibre de V1
#: ne satisfait pas les equations (il est faux de quelques 1e-3 pu) et le port
#: Python le corrige. Seuil 1e-3 pour symSG.
SEUIL_SMTR_SYMSG = 1e-3

EQUIPEMENTS = {
    "SM": (SynchronousMachine, {}),
    "GFM1": (GFM1, {}),
    "GFM2": (GFM2, {}),
    "GFM3": (GFM3, {}),
    "GFM4": (GFM4, {}),
    "GFL1": (GFL1, {}),
    "GFL2": (GFL2, {}),
}


@pytestmark_csv
@pytest.mark.parametrize("nom", list(EQUIPEMENTS))
def test_admittance_contre_matlab(nom):
    cls, kw = EQUIPEMENTS[nom]
    T = pd.read_csv(CSV)
    assert _ecart(banc(cls, **kw), T, nom) < SEUIL


@pytestmark_csv
def test_admittance_smtr_contre_matlab():
    T = pd.read_csv(CSV)
    sym = banc(SynchronousMachineTR, init="symSG", **SMTR_KW)
    assert _ecart(sym, T, "SMTR") < SEUIL_SMTR_SYMSG
    # init physique : l'ecart est reel et dit ce qu'il est (equilibre V1 faux)
    phys = banc(SynchronousMachineTR, **SMTR_KW)
    assert _ecart(phys, T, "SMTR") > 0.05


@pytestmark_csv
def test_profil_harmonise_contre_matlab():
    """Profil a parametres communs (statisme 4 %, H = 5 s, X de couplage 0,15) :
    les GFM, GFL et la SMTR donnent la meme admittance que MATLAB."""
    T = pd.read_csv(CSV)
    XT, RT, STAT, H, KA, TAU_CL = 0.15, 0.005, 0.04, 5.0, 300.0, 1e-3
    wf = 1 / (2 * STAT * H)
    conv = dict(R2_pu=RT, L2_pu=XT, tau_CL=TAU_CL)
    gfm = dict(conv, mp=STAT, wf_p=wf, wf_q=wf)
    cas = {
        "H_GFM1_nq4": (GFM1, dict(gfm, nq=0.04)),
        "H_GFM1_nqAVR": (GFM1, dict(gfm, nq=1 / KA)),
        "H_GFM3_nq4": (GFM3, dict(gfm, nq=0.04)),
        "H_GFL1": (GFL1, conv),
        "H_GFL2": (GFL2, conv),
    }
    for nom, (cls, kw) in cas.items():
        assert _ecart(banc(cls, **kw), T, nom) < SEUIL, nom
    sm = dict(Rt=RT, Lt=XT, H=H, KD=1.0, mp=100 * STAT, Ka=KA)
    el = banc(SynchronousMachineTR, init="symSG", **sm)
    assert _ecart(el, T, "H_SMTR") < SEUIL_SMTR_SYMSG
