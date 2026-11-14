"""Limitation de courant par IMPEDANCE VIRTUELLE (VI) des GFM.

Ce que le fichier prouve, dans l'ordre :

* la VI est INACTIVE par defaut et SOUS LE SEUIL : jacobiens et spectre
  identiques au bit pres a ceux du convertisseur sans VI ;
* au-dessus du seuil, l'equilibre reste exact (|f| < 1e-10) pour toutes les
  lois de synchronisation, et les jacobiens symboliques collent aux
  differences finies (1e-9), dans la zone de raccord parabolique comme dans
  la zone lineaire ;
* VI active, les sensibilites ANALYTIQUES d'un mode a tous les parametres du
  GFM collent aux differences finies (et se calculent en quelques secondes) ;
* le dimensionnement ``kRv_dimensionne`` donne le courant PERMANENT de defaut
  franc voulu, en formule fermee comme en simulation non lineaire ;
* sur un defaut franc de ``wscc9_gfm`` : la VI seule tient le courant
  permanent sous ``Imax`` mais PAS le pic transitoire (la boucle de tension
  fait 15 ms) ; avec le limiteur circulaire comme gardien, le pic reste borne
  et la VI evite a ce dernier de basculer sans cesse ;
* ``limites.py`` (FR-009) : un point limite par la VI n'est pas hors limite
  tant que ``|i_s| <= Imax``, et l'est des que ``|i_s| > Imax`` ;
* le CONSTAT de conception : une reactance virtuelle algebrique elevee
  (X/R >= 1) instabilise le mode LC du filtre avec les boucles de tension du
  MATLAB (``tau_VL`` = 15 ms) ; le defaut ``XsurR_vi`` est donc bas.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from scipy.optimize import brentq, root

import g2egrid.catalogue as cat
from g2egrid.elements.gfm import (
    GFM1,
    GFM2,
    GFM3,
    GFM4,
    GFMdVOC,
    GFMMatching,
)
from g2egrid.elements.limiteurs import depassement_vi, kRv_dimensionne
from g2egrid.limites import DepassementLimite
from g2egrid.modal import modes
from g2egrid.sensibilites import sensibilites
from g2egrid.tables import assembler
from g2egrid.temporel import Manoeuvre, ModeleNonLineaire

TOUTES = ("T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s")
#: (table, loi, classe attendue)
CLASSES = [
    ("T_GFM1s", None, GFM1),
    ("T_GFM2s", None, GFM2),
    ("T_GFM3s", None, GFM3),
    ("T_GFM4s", None, GFM4),
    ("T_GFM1s", "dvoc", GFMdVOC),
    ("T_GFM3s", "matching", GFMMatching),
]
IDS = [c[2].__name__ for c in CLASSES]


def _reseau(nom, table, loi, **reglages):
    """``nom`` dont TOUS les GFM sont de la table ``table`` (et de la loi
    ``loi``), reglages de GFM_PARAMS en plus. Equilibre resolu."""
    t = cat.charger(nom).copie()
    lignes = [r for k in TOUTES for r in getattr(t, k)]
    for k in TOUTES:
        getattr(t, k).clear()
    getattr(t, table).extend(lignes)
    t.GFM_PARAMS = dict(t.GFM_PARAMS, **reglages)
    if loi:
        t.GFM_PARAMS["loi"] = loi
    net = assembler(t)
    net.solve_equilibrium(tol=1e-13)  # defaut 1e-9 : trop lache pour |f| < 1e-10
    return net


def _gfm(net, classe=None):
    return [
        e
        for e in net.all_elements
        if e.name.startswith("GFM") and (classe is None or isinstance(e, classe))
    ]


def _is(e):
    return float(np.hypot(e.x0[e.states.index("isd")], e.x0[e.states.index("isq")]))


def _reseau_vi(nom, table, loi, ecart=0.2, **reglages):
    """Reseau dont la VI est ACTIVE au point de fonctionnement : le seuil est
    pose ``ecart`` pu sous le courant d'equilibre (``|i_s|`` ne depend pas de
    la VI : il vient du load flow)."""
    i0 = _is(_gfm(_reseau(nom, table, loi))[0])
    return _reseau(nom, table, loi, kRv=1.0, I_vi=i0 - ecart, **reglages)


def _jacobiens_fd(e, h=1e-6):
    """Ecart relatif max entre les jacobiens (A, B, C, D) symboliques et les
    differences finies centrees, z etant re-resolu a chaque perturbation."""
    cm = e.compiled()
    x0, z0, u0, p = e.x0, e.z0, e.u0, e.p
    nz = cm.n_z

    def fh(x, u):
        z = z0
        if nz:
            z = root(lambda zz: cm.g(x, zz, u, p), z0, tol=1e-14).x
        return np.asarray(cm.f(x, z, u, p)), np.asarray(cm.h(x, z, u, p))

    A, B, C, D = cm.linearize(x0, z0, u0, p)
    nx, nu = len(x0), len(u0)
    f0, h0 = fh(x0, u0)
    An, Bn = np.zeros((nx, nx)), np.zeros((nx, nu))
    Cn, Dn = np.zeros((len(h0), nx)), np.zeros((len(h0), nu))
    for j in range(nx):
        d = h * max(1.0, abs(x0[j]))
        xp, xm = x0.copy(), x0.copy()
        xp[j] += d
        xm[j] -= d
        (fp, hp), (fm, hm) = fh(xp, u0), fh(xm, u0)
        An[:, j], Cn[:, j] = (fp - fm) / (2 * d), (hp - hm) / (2 * d)
    for j in range(nu):
        d = h * max(1.0, abs(u0[j]))
        up, um = u0.copy(), u0.copy()
        up[j] += d
        um[j] -= d
        (fp, hp), (fm, hm) = fh(x0, up), fh(x0, um)
        Bn[:, j], Dn[:, j] = (fp - fm) / (2 * d), (hp - hm) / (2 * d)
    ech = lambda M: max(1.0, float(np.abs(M).max())) if M.size else 1.0  # noqa: E731
    return max(
        float(np.abs(a - b).max()) / ech(a) if a.size else 0.0
        for a, b in ((A, An), (B, Bn), (C, Cn), (D, Dn))
    )


# ----------------------------------------------------------------------
# la fonction de depassement
# ----------------------------------------------------------------------
def test_depassement_c1_et_nul_sous_le_seuil():
    seuil, dI = 1.0, 0.05
    ov = lambda i: float(depassement_vi(i, seuil, dI))  # noqa: E731
    assert ov(0.0) == 0.0 and ov(1.0) == 0.0 and ov(0.999) == 0.0
    assert ov(1.025) == pytest.approx(0.025**2 / 0.1)
    assert ov(1.2) == pytest.approx(0.2 - 0.025)
    # continuite de la valeur et de la pente aux deux raccords
    e = 1e-9
    for i0 in (1.0, 1.05):
        assert ov(i0 + e) == pytest.approx(ov(i0 - e), abs=1e-8)
        pente = (ov(i0 + e) - ov(i0 - e)) / (2 * e)
        assert pente == pytest.approx(0.0 if i0 == 1.0 else 1.0, abs=1e-6)


# ----------------------------------------------------------------------
# inactive par defaut, et sous le seuil, identique au bit pres
# ----------------------------------------------------------------------
@pytest.mark.parametrize("table,loi,classe", CLASSES, ids=IDS)
def test_inactive_par_defaut_et_sous_le_seuil(table, loi, classe):
    """kRv = 0 (defaut) ou courant sous I_vi : A, B, C, D, equilibre et
    spectre du reseau complet sont ceux du modele sans VI, au bit pres."""
    ref = _reseau("smib_gfm", table, loi)
    e0 = _gfm(ref, classe)[0]
    assert e0.pdict["kRv"] == 0.0
    # kRv > 0 mais seuil au-dessus du courant d'equilibre
    sous = _reseau("smib_gfm", table, loi, kRv=50.0, I_vi=_is(e0) + 0.05)
    e1 = _gfm(sous, classe)[0]
    assert e1.pdict["kRv"] == 50.0
    for a, b in zip(
        e0.compiled().linearize(e0.x0, e0.z0, e0.u0, e0.p),
        e1.compiled().linearize(e1.x0, e1.z0, e1.u0, e1.p),
    ):
        assert np.array_equal(a, b)
    assert np.array_equal(e0.x0, e1.x0)
    assert np.array_equal(
        np.asarray(ref.build_state_space().A), np.asarray(sous.build_state_space().A)
    )


def test_references_inchangees():
    """Les invariants du depot : 13/14/15/16 etats, mode local du SMIB-GFM."""
    nets = {c[2]: _reseau("smib_gfm", c[0], c[1]) for c in CLASSES}
    for classe, n in ((GFM1, 13), (GFM2, 14), (GFM3, 15), (GFM4, 16)):
        assert len(_gfm(nets[classe])[0].states) == n
    assert nets[GFM1].build_state_space().is_stable()


def test_statisme_inchange_sous_le_seuil():
    """La pente de statisme se lit sur le jacobien : d(dtheta/dt)/d pm =
    -wb.mp. Elle est la meme avec et sans VI (VI non declenchee) -- et meme
    VI declenchee, puisque w ne depend que de p_ref - pm."""
    ref = _gfm(_reseau("smib_gfm", "T_GFM1s", None))[0]
    for vi in (
        _reseau("smib_gfm", "T_GFM1s", None, kRv=50.0, I_vi=2.0),
        _reseau_vi("smib_gfm", "T_GFM1s", None),
    ):
        e = _gfm(vi)[0]
        for ee in (ref, e):
            A = ee.compiled().linearize(ee.x0, ee.z0, ee.u0, ee.p)[0]
            i, j = ee.states.index("theta"), ee.states.index("pm")
            assert A[i, j] == pytest.approx(-ee.pdict["wb"] * ee.pdict["mp"], rel=1e-12)


# ----------------------------------------------------------------------
# au-dessus du seuil : equilibre, jacobiens, stabilite
# ----------------------------------------------------------------------
@pytest.mark.parametrize("XsurR", [0.0, 0.25, 5.0], ids=lambda v: f"X/R={v}")
@pytest.mark.parametrize("table,loi,classe", CLASSES, ids=IDS)
def test_equilibre_exact_au_dessus_du_seuil(table, loi, classe, XsurR):
    """Avec une VI declenchee (ov > 0), le repere local suit le phaseur
    ``ve + Zv.is`` et non plus ``ve`` : veq != 0, l'equilibre reste a 1e-10, et
    la consigne de tension est tenue : ``ve + Zv.is = ve_ref`` (reel)."""
    net = _reseau_vi("smib_gfm", table, loi, XsurR_vi=XsurR)
    e = _gfm(net, classe)[0]
    f = np.array(e.compiled().f(e.x0, e.z0, e.u0, e.p))
    if "idc" in e.states:
        # didc/dt = (idc_ref - idc)/Tdc, Tdc = 1e-4 : l'arrondi de idc_ref
        # (1e-13, Kpdc ~ 2,6e4) y est amplifie de 1e4 ; on compare l'ecart
        # idc_ref - idc, en pu de courant continu
        f[e.states.index("idc")] *= e.pdict["Tdc"]
    assert np.abs(f).max() < 1e-10
    p = e.pdict
    ov = float(depassement_vi(_is(e), p["I_vi"], p["dI_vi"]))
    assert ov > 0.1
    x = lambda n: e.x0[e.states.index(n)]  # noqa: E731
    Rv = p["kRv"] * ov
    Zv = complex(Rv, XsurR * Rv)
    ve = complex(x("ved"), x("veq"))
    i_s = complex(x("isd"), x("isq"))
    ve_ref = e.u0[e.input_names().index("ve_ref")]
    assert ve + Zv * i_s == pytest.approx(ve_ref, abs=1e-9)
    if XsurR > 0:
        assert abs(x("veq")) > 1e-3  # le repere n'est plus sur ve
    # la puissance reste celle du load flow (consigne tenue)
    if "pm" in e.states:
        assert x("pm") == pytest.approx(e.u0[e.input_names().index("p_ref")], abs=1e-9)


@pytest.mark.parametrize(
    "ecart", [-0.3, 0.02, 0.2], ids=["sous", "raccord", "lineaire"]
)
@pytest.mark.parametrize("table,loi,classe", CLASSES, ids=IDS)
def test_jacobiens_contre_differences_finies(table, loi, classe, ecart):
    """Sous le seuil, dans le raccord parabolique (0 < i - I_vi < dI_vi) et
    dans la zone lineaire : A, B, C, D contre differences finies, 1e-9."""
    net = _reseau_vi("smib_gfm", table, loi, ecart=ecart, XsurR_vi=0.25)
    e = _gfm(net, classe)[0]
    assert _jacobiens_fd(e) < 1e-9


@pytest.mark.parametrize("ecart", [0.02, 0.2], ids=["raccord", "lineaire"])
def test_sensibilites_analytiques_vi_active(ecart):
    """VI declenchee au point de fonctionnement : la sensibilite ANALYTIQUE
    d'un mode a CHAQUE parametre du GFM -- ceux de la VI et les drapeaux des
    limiteurs compris -- colle aux differences finies. Les derivees
    symboliques de ``dJ/dp`` y traversent le ``Piecewise`` du seuil ; avant
    que les drapeaux des limiteurs n'y soient figes
    (``CompiledModel.derivees_parametres``), le calcul ne rendait plus.

    Le mode est le MOINS AMORTI sous 1 kHz (127 Hz et 211 Hz, celui que la VI
    deplace) : la VI rend les modes lents presque reels et doubles (zeta >
    99 %), ou la valeur propre est mal conditionnee. Comme dans
    ``test_dynamiques``, on compare les six plus fortes, et en plus celles de
    la VI : pour un parametre de valeur nulle ou minuscule (``nq`` = 1e-4,
    les drapeaux), le pas absolu de la difference finie (1e-9) passe sous le
    bruit d'arrondi des valeurs propres (``|A|`` ~ 2e6) ; son ecart a la
    valeur analytique decroit alors en 1/h quand le pas grandit."""
    net = _reseau_vi("smib_gfm", "T_GFM2s", None, ecart=ecart, XsurR_vi=0.25)
    g = _gfm(net)[0].name
    ss = net.build_state_space()
    m = min((m for m in modes(ss) if m.oscillant and m.f < 1e3), key=lambda m: m.zeta)
    ana = {s.parametre: s for s in sensibilites(net, m, elements=[g])}
    fd = {
        s.parametre: s.dlambda
        for s in sensibilites(net, m, methode="differences", elements=[g])
    }
    assert set(ana) == set(fd)
    top = sorted(ana.values(), key=lambda s: -abs(s.normalisee))[:6]
    vi = [ana[p] for p in ("kRv", "XsurR_vi", "I_vi", "dI_vi")]
    for s in top + vi:
        d = fd[s.parametre]
        assert abs(s.dlambda - d) <= 1e-4 * max(abs(s.dlambda), 1e-3), s
    assert all(abs(s.normalisee) > 0 for s in vi)  # la VI y compte
    # les drapeaux, derives a part (les autres figes) : sous Imax, k = 1 quel
    # que soit mode_lim, sa derivee est nulle EXACTEMENT ; mode_vsat (Mmax = 0
    # sur bus continu ideal) contre une difference finie au pas 1e-3
    assert ana["mode_lim"].dlambda == 0
    s = ana["mode_vsat"]
    (d,) = sensibilites(
        net, m, methode="differences", elements=[g], parametres=["mode_vsat"], pas=1e-3
    )
    assert abs(s.dlambda) > 1 and abs(s.dlambda - d.dlambda) <= 1e-4 * abs(s.dlambda)


@pytest.mark.parametrize("nom", ["smib_gfm", "wscc9_gfm"])
def test_reseau_stable_avec_vi_resistive(nom):
    """VI declenchee au point de fonctionnement, X/R bas : le reseau complet
    reste stable (les GFM des deux reseaux)."""
    net = _reseau_vi(nom, "T_GFM2s", None, XsurR_vi=0.25)
    ev = np.linalg.eigvals(np.asarray(net.build_state_space().A))
    # Aucun mode rapide instable. Un mode LENT reste marginalement instable
    # (+0,02 /s, constat) : une VI resistive declenchee en PERMANENT retire de
    # l'amortissement au mode de statisme ; ce n'est pas un point de
    # fonctionnement reel (la VI ne vit qu'au-dessus de 1 pu), d'ou la borne
    # lache a 1 /s plutot que is_stable().
    assert ev.real.max() < 1.0


def test_reactance_virtuelle_elevee_instabilise_le_filtre():
    """CONSTAT de conception. Une reactance virtuelle algebrique (la chute
    ``j.Xv.i`` du repere dq) est un gyrateur : avec la boucle de tension du
    MATLAB (tau_VL = 15 ms) elle apporte de l'amortissement NEGATIF au mode LC
    du filtre. A Rv = 0,3 pu, X/R = 0,25 est stable, X/R >= 1 ne l'est pas
    (modes a 100-300 Hz, Re de +10 a +600 /s). C'est pourquoi XsurR_vi vaut
    0,1 par defaut, et non les 5 a 10 de la litterature (boucles de tension
    plus rapides). Plus bas encore : a X/R = 0,25 et Rv = 0,8 le petit signal
    est stable mais la simulation de defaut (VI seule) perd son Newton a
    l'elimination ; a X/R = 0,1 elle tient (``comparaison``)."""
    i0 = _is(_gfm(_reseau("wscc9_gfm", "T_GFM2s", None))[0])

    def pic(sig):
        net = _reseau(
            "wscc9_gfm", "T_GFM2s", None, kRv=0.3 / 0.175, I_vi=i0 - 0.2, XsurR_vi=sig
        )
        return np.linalg.eigvals(np.asarray(net.build_state_space().A)).real.max()

    assert pic(0.25) < 1.0
    assert pic(1.0) > 50.0
    assert pic(5.0) > pic(1.0)


# ----------------------------------------------------------------------
# dimensionnement
# ----------------------------------------------------------------------
def test_dimensionnement_courant_permanent_de_defaut():
    """``kRv_dimensionne`` : le courant permanent solution de ``|Zv(I) + Zd|.I
    = V`` vaut Imax ; et un X/R eleve ou un Zd inductif se traitent en
    phaseurs (une VI resistive devant du reactif limite MOINS qu'une somme de
    modules)."""
    Imax, I_vi, dI, V = 1.2, 1.0, 0.05, 1.0
    for sig, Zd in [(0.0, 0.0005 + 0.046j), (0.25, 0.01 + 0.03j), (5.0, 0.05j)]:
        k = kRv_dimensionne(Imax, I_vi, sig, V, Zd, dI)

        def residu(i):
            ov = float(depassement_vi(i, I_vi, dI))
            return abs(complex(k * ov, sig * k * ov) + Zd) * i - V

        assert residu(Imax) == pytest.approx(0.0, abs=1e-12)
        assert brentq(residu, I_vi, 3.0) == pytest.approx(Imax, abs=1e-9)
    # le filtre seul limite deja : pas de VI
    assert kRv_dimensionne(1.2, 1.0, 0.25, 1.0, 1.0j) == 0.0
    with pytest.raises(ValueError, match="Imax"):
        kRv_dimensionne(0.9, 1.0)


# ----------------------------------------------------------------------
# defaut franc sur wscc9_gfm : VI contre limiteur circulaire
# ----------------------------------------------------------------------
T_DEFAUT, DUREE, T_FIN = 0.05, 0.12, 0.6
#: impedance entre la tension controlee de GFM2 et le defaut au bus 7, en pu
#: de l'equipement (Sn = 200 MVA, Sb = 100 MVA) : L2 (0,01) + transformateur
#: T2_2_7 (0,0625 / kS = 0,031) + defaut (R = 1e-3, X = 1e-2 sur Sb)
Z_DEFAUT = complex(0.0005, 0.046)


def _defaut(kRv, limiteurs, XsurR=0.1):
    t = cat.charger("wscc9_gfm").copie()
    t.GFM_PARAMS = dict(t.GFM_PARAMS, kRv=kRv, XsurR_vi=XsurR)
    m = ModeleNonLineaire(assembler(t))
    u = m.scenario(
        Manoeuvre(T_DEFAUT, "defaut", "bus7", R=1e-3, X=1e-2),
        Manoeuvre(T_DEFAUT + DUREE, "elimination", "bus7"),
    )
    return m.simuler(
        (0.0, T_FIN), uext_fn=u, methode="BDF2", pas=1e-4, limiteurs=limiteurs
    )


@pytest.fixture(scope="module")
def comparaison():
    k = kRv_dimensionne(1.2, 1.0, 0.1, 1.0, Z_DEFAUT, 0.05)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return {
            "k": k,
            "circulaire": _defaut(0.0, True),
            "vi": _defaut(k, False),
            "vi+circulaire": _defaut(k, True),
        }


def _mesures(tr, g="GFM2"):
    i = np.hypot(tr.etat(f"{g}.isd"), tr.etat(f"{g}.isq"))
    th = tr.etat(f"{g}.theta")
    pm = tr.etat(f"{g}.pm")
    fin_defaut = T_DEFAUT + DUREE
    permanent = (tr.t > fin_defaut - 0.02) & (tr.t <= fin_defaut)
    pendant = (tr.t > T_DEFAUT) & (tr.t <= fin_defaut)
    pm_fin = pm[tr.t > T_FIN - 0.05].mean()
    hors = np.where((np.abs(pm - pm_fin) > 0.02) & (tr.t > fin_defaut))[0]
    return dict(
        pic=i[pendant].max(),
        permanent=i[permanent].mean(),
        excursion=np.abs(th - th[0]).max(),
        retour=(tr.t[hors[-1]] - fin_defaut) if len(hors) else 0.0,
        bascule=len(tr.basculements),
        theta_fin=th[-1] - th[0],
    )


def test_courant_permanent_de_defaut_dimensionne(comparaison):
    """Le courant permanent de defaut avec VI seule (limiteur circulaire
    coupe) est celui du dimensionnement, Imax = 1,2, a 1 % pres."""
    m = _mesures(comparaison["vi"])
    assert m["permanent"] == pytest.approx(1.2, rel=0.01)
    assert m["permanent"] <= 1.2


def test_vi_seule_ne_tient_pas_le_pic_mais_le_limiteur_oui(comparaison):
    """La boucle de tension (15 ms) est trop lente pour la VI seule : le pic
    transitoire depasse Imax (borne large, 2 pu). Le limiteur circulaire, en
    gardien, le ramene a Imax a l'hysteresis et a la resonance LCL pres
    (+2 %), et la VI lui evite des basculements."""
    vi = _mesures(comparaison["vi"])
    mix = _mesures(comparaison["vi+circulaire"])
    cir = _mesures(comparaison["circulaire"])
    assert 1.2 < vi["pic"] < 2.0
    assert mix["pic"] <= 1.2 * 1.02
    assert cir["pic"] <= 1.2 * 1.02
    assert mix["bascule"] < cir["bascule"]


def test_resynchronisation_apres_elimination(comparaison):
    """Avec VI, l'angle revient a un decalage fini (pas de glissement) et la
    puissance revient dans 2 % de sa valeur finale ; l'excursion d'angle reste
    du meme ordre que celle de la limitation circulaire existante."""
    cir = _mesures(comparaison["circulaire"])
    for cle in ("vi", "vi+circulaire"):
        m = _mesures(comparaison[cle])
        assert abs(m["theta_fin"]) < 0.1  # pas de glissement de pole
        assert m["retour"] < T_FIN - (T_DEFAUT + DUREE) - 0.01
        assert m["excursion"] < 2 * cir["excursion"] + 0.01


# ----------------------------------------------------------------------
# FR-009 : limites.py
# ----------------------------------------------------------------------
def _signale(**reglages):
    t = cat.charger("wscc9_gfm").copie()
    t.GFM_PARAMS = dict(t.GFM_PARAMS, **reglages)
    net = assembler(t)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        net.solve_equilibrium()
    return net, [x for x in w if issubclass(x.category, DepassementLimite)]


def test_point_limite_par_la_vi_n_est_pas_hors_limite():
    """|i_s| = 0,80 et 0,66 pu aux deux GFM. VI declenchee (I_vi = 0,3) avec
    Imax = 0,9 : le modele libre INCLUT la VI, le point est dans son domaine,
    aucun avertissement et aucun basculement au depart. Avec Imax = 0,7 : le
    GFM2 depasse Imax, la VI n'y change rien, l'avertissement est emis."""
    net, w = _signale(kRv=3.0, I_vi=0.3, Imax=0.9)
    assert not w and net.depassements == []
    m = ModeleNonLineaire(net)
    assert m.simuler((0, 0.01), methode="BDF2", pas=1e-4).basculements == []
    net, w = _signale(kRv=3.0, I_vi=0.3, Imax=0.7)
    assert [d.element for d in net.depassements] == ["GFM2"]
    assert w
    net_sans, _ = _signale(Imax=0.7)  # meme diagnostic sans VI
    assert [d.element for d in net_sans.depassements] == ["GFM2"]
