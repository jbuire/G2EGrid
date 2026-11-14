"""PSS et AVR de la machine synchrone : conventions, deuxieme voie, ANDES (FR-007).

Avec ses constantes par defaut (``K_PSS = 20``, passe-bas 0,03 s, washout
2 s, avance 0,05/0,02 puis RETARD 3/5,4), le PSS DESAMORTIT le mode local du
SMIB, de wscc9 et de kundur2area. Erreur de signe, ou reglage inadapte ? Ce
fichier tranche par des voies independantes du jumeau MATLAB, qui partage les
equations de la machine et ne prouve donc rien sur ce point.

1. Conventions, lues sur le modele linearise de l'element : l'entree du PSS
   est ``w - 1``, sa fonction de transfert est celle du schema IEEE PSS1A a
   gain positif, sa sortie s'ajoute a la consigne de l'AVR, et le retour de
   stabilisation de l'AVR est retranche. Ce sont les conventions de
   l'IEEEST et de l'EXST1 d'ANDES.
2. Deuxieme voie, analytique : GEP(s), de la consigne de l'AVR au couple
   electrique, vitesse(s) et angle(s) des machines figes. A la frequence du
   mode local, GEP est en retard de PLUS de 90 degres, et le PSS n'apporte
   que 0 degre : le couple du PSS y est en opposition partielle avec l'ecart
   de vitesse, il desamortit. Le retard vient du retour de stabilisation de
   l'AVR (``Ka.Kfd/Tfd = 3``) : un modele de Heffron-Phillips ecrit ici a la
   main le retrouve, et sans ce retour le meme PSS amortit. La sensibilite
   de la valeur propre au gain (differences finies, residu exact) a le
   signe que GEP predit.
3. ANDES, AVR et PSS actifs (EXST1 + IEEEST, memes structures), sur le
   reseau a deux zones de Kundur : meme desamortissement des modes locaux,
   a moins de 1 point d'amortissement ; l'ecart du mode inter-zones (jusqu'a
   1,7 point) vient de la f.e.m. de vitesse du stator, ce que le dernier test
   demontre en la retirant.

Verdict : pas d'erreur de signe ni de convention ; un PSS sans avance de
phase associe a un AVR qui retarde de plus de 90 degres au-dela de ~1,5 Hz,
et des modes locaux que le regulateur de vitesse du catalogue (statisme
0,5 %, 0,2 s) pousse vers 2 Hz. Voir ``docs/equivalence_matlab.md`` §G.

Les tolerances des valeurs de reference couvrent l'arrondi de la valeur
stockee (5 a 6 chiffres) : le calcul est deterministe (resolutions lineaires
denses), et une derive au-dela signale un changement du modele, a expliquer.
"""

from __future__ import annotations

import importlib.util
import math

import numpy as np
import pytest
import scipy.linalg as sl
from scipy.optimize import brentq

from g2egrid.catalogue import charger
from g2egrid.elements.sync_machine import SM_DEFAULT, SynchronousMachine
from g2egrid.tables import assembler

#: bande des modes electromecaniques, comme la fiche FR-007 [Hz]
BANDE_EM = (0.1, 3.0)
RESEAUX = ("smib", "wscc9", "kundur2area")
ANDES = importlib.util.find_spec("andes") is not None
sans_andes = pytest.mark.skipif(
    not ANDES, reason="ANDES absent : uv sync --group validation"
)


# ----------------------------------------------------------------------
# outils
# ----------------------------------------------------------------------
def _reseau(nom: str, **sm):
    """Reseau du catalogue a l'equilibre, reglages des machines surcharges."""
    t = charger(nom)
    t.SM_PARAMS = {**(t.SM_PARAMS or {}), **sm}
    net = assembler(t)
    net.solve_equilibrium()
    return net, net.build_state_space()


def _machines(net) -> list[str]:
    return [e.name for e in net.elements if isinstance(e, SynchronousMachine)]


def _moins_amorti(ss) -> complex:
    """Mode le moins amorti de la bande electromecanique (definition FR-007)."""
    lam = np.linalg.eigvals(ss.A)
    f = lam.imag / (2 * np.pi)
    lam = lam[(f > BANDE_EM[0]) & (f < BANDE_EM[1])]
    return lam[np.argmin(-lam.real / np.abs(lam))]


def _f_zeta(lam: complex) -> tuple[float, float]:
    return lam.imag / (2 * np.pi), -100 * lam.real / abs(lam)


def pss_analytique(p: dict, s):
    """Le PSS d'apres le schema PSS1A : transducteur, gain, washout, deux
    avances-retards. Ecart de vitesse (pu) -> signal ajoute a la consigne."""
    return (
        p["K_PSS"]
        * s * p["T_HP"] / (1 + s * p["T_HP"])
        / (1 + s * p["T_LP"])
        * (1 + s * p["T1n"]) / (1 + s * p["T1d"])
        * (1 + s * p["T2n"]) / (1 + s * p["T2d"])
    )  # fmt: skip


def gep(ss, s, machine: str, figees: list[str]):
    """GEP(s) : consigne ``v_ref`` de ``machine`` -> son couple ``Ce``.

    Les etats de vitesse et d'angle des machines ``figees`` sont retires
    (ecart nul) : c'est la definition de Heffron-Phillips, ``Dd = Dw = 0``.
    Pour un reseau multi-machines, on les fige toutes, sans quoi GEP porte
    les resonances des autres modes electromecaniques.
    """
    fig = {f"{m}.{n}" for m in figees for n in ("Dwr", "theta")}
    k = [i for i, n in enumerate(ss.states) if n not in fig]
    A = ss.A[np.ix_(k, k)]
    B = ss.B[k, ss.iu(f"{machine}.v_ref")]
    C = ss.C[ss.iy(f"{machine}.Ce"), k]
    D = ss.D[ss.iy(f"{machine}.Ce"), ss.iu(f"{machine}.v_ref")]
    return C @ np.linalg.solve(s * np.eye(len(k)) - A, B) + D


# ----------------------------------------------------------------------
# 1. conventions de la chaine dw -> PSS -> AVR
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def smib():
    net, ss = _reseau("smib")
    sm = next(e for e in net.elements if isinstance(e, SynchronousMachine))
    A, B, C, _ = sm.compiled().linearize(sm.x0, sm.z0, sm.u0, sm.p)
    return net, ss, sm, A, B, C


def test_entree_du_pss_est_w_moins_1(smib):
    """``Dwr`` est l'ecart a la vitesse de reference ``w_ref`` (entree
    externe, 1 pu) : ``Omega = w_ref + Dwr``, et ``Dwr`` est nul a
    l'equilibre (au raffinement de Newton pres, 1e-30 mesure)."""
    _net, _ss, sm, _A, _B, C = smib
    cm = sm.compiled()
    y = cm.h(sm.x0, sm.z0, sm.u0, sm.p)
    om = cm.iy("Omega")
    assert y[om] == pytest.approx(1.0, abs=1e-12)
    assert abs(sm.x0[sm.states.index("Dwr")]) < 1e-12
    assert C[om, sm.states.index("Dwr")] == 1.0


@pytest.mark.parametrize("f", [0.05, 0.5, 1.7, 5.0])
def test_pss_est_le_schema_pss1a_a_gain_positif(smib, f):
    """La fonction de transfert ``Dwr -> vpss`` du modele linearise est celle
    du schema PSS1A ecrite a la main, gain ``K_PSS = +20``. Meme fraction
    rationnelle calculee de deux facons : 1e-10 en relatif."""
    _net, _ss, sm, A, _B, _C = smib
    i = [sm.states.index(n) for n in ("Dw1", "v1", "v2", "vpss")]
    j = sm.states.index("Dwr")
    # le bloc PSS n'est attaque que par la vitesse
    autres = [k for k in range(len(sm.states)) if k not in i and k != j]
    assert not np.any(A[np.ix_(i, autres)])
    s = 2j * np.pi * f
    h = np.linalg.solve(s * np.eye(4) - A[np.ix_(i, i)], A[i, j])[3]
    assert h == pytest.approx(pss_analytique(sm.pdict, s), rel=1e-10)
    assert sm.pdict["K_PSS"] == 20.0


def test_sortie_du_pss_ajoutee_a_la_consigne_de_l_avr(smib):
    """``de2/dt = Ka.(v_ref - e1 - e3 + vpss)/Ta`` : le PSS entre avec le
    signe de la consigne, la mesure et le retour de stabilisation sont
    retranches -- les conventions de l'EXST1."""
    _net, _ss, sm, A, B, _C = smib
    p = sm.pdict
    k = sm.states.index("e2")
    gain = p["Ka"] / p["Ta"]
    assert B[k, sm.input_names().index("v_ref")] == pytest.approx(gain, rel=1e-14)
    assert A[k, sm.states.index("vpss")] == pytest.approx(gain, rel=1e-14)
    assert A[k, sm.states.index("e1")] == pytest.approx(-gain, rel=1e-14)
    assert A[k, sm.states.index("e3")] == pytest.approx(-gain, rel=1e-14)


# ----------------------------------------------------------------------
# 2. deuxieme voie : GEP, Heffron-Phillips, sensibilite
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def smib_sans_pss():
    net, ss = _reseau("smib", K_PSS=0.0)
    return net, ss, _moins_amorti(ss)


def test_gep_du_smib_en_retard_de_plus_de_90_degres(smib_sans_pss):
    """SMIB du catalogue, PSS coupe : mode local 1,6816 Hz, 19,054 %.
    A cette frequence, GEP = 1,3878 a -102,646 degres, le PSS par defaut
    +0,234 degre : ``Re(GEP.PSS) = -3,493`` pu de couple par pu de vitesse,
    un amortissement NEGATIF."""
    net, ss, lam = smib_sans_pss
    f, z = _f_zeta(lam)
    assert f == pytest.approx(1.6816, abs=5e-5)
    assert z == pytest.approx(19.054, abs=5e-4)
    s = 2j * np.pi * f
    g = gep(ss, s, "SM1", _machines(net))
    h = pss_analytique(SM_DEFAULT, s)
    assert abs(g) == pytest.approx(1.38780, rel=1e-5)
    assert math.degrees(np.angle(g)) == pytest.approx(-102.6462, abs=2e-4)
    assert math.degrees(np.angle(h)) == pytest.approx(0.2342, abs=2e-4)
    assert (g * h).real == pytest.approx(-3.49332, rel=1e-5)


@pytest.mark.parametrize("nom", RESEAUX)
def test_couple_du_pss_desamortissant_pour_chaque_machine(nom):
    """Sur chaque reseau, a la frequence du mode local (PSS coupe), la phase
    de ``GEP.PSS`` est au-dela de -90 degres pour CHAQUE machine : le PSS par
    defaut ne peut que desamortir ce mode. Mesure : -97,6 a -106,7 degres."""
    net, ss = _reseau(nom, K_PSS=0.0)
    s = 1j * _moins_amorti(ss).imag
    for m in _machines(net):
        phase = math.degrees(
            np.angle(gep(ss, s, m, _machines(net)) * pss_analytique(SM_DEFAULT, s))
        )
        assert -180.0 < phase < -95.0, (m, phase)


def _gep_heffron_phillips(f: float, Kfd: float) -> complex:
    """GEP du SMIB du catalogue par le modele a un axe de Heffron-Phillips.

    Voie independante : rien de g2egrid. Machine a un axe (``Xd``, ``Xq``,
    ``X'd``, ``T'd0`` deduits des inductances de ``SM_DEFAULT``), reseau
    reduit a ``Re + jXe`` (transfo 0,15 + ligne 0,02 + j0,5 + source) vers un
    bus infini a 1 pu, point de fonctionnement P = 0,9 pu sous 1,03 pu,
    charges de ligne ignorees. K1..K6 par differences finies sur les
    equations algebriques du stator ; AVR complet (mesure, amplificateur,
    excitatrice, retour de stabilisation).
    """
    p = SM_DEFAULT
    wb = 2 * math.pi * 50.0
    Xd, Xq = p["Ll"] + p["Lad"], p["Ll"] + p["Laq"]
    Xdp = p["Ll"] + p["Lad"] * p["Lfd"] / (p["Lad"] + p["Lfd"])
    Td0p = (p["Lad"] + p["Lfd"]) / (wb * p["Rfd"])
    Ra = p["Ra"]
    Re, Xe = 0.02 + 1e-4, 0.15 + 0.5 + 1e-3
    P, Vt, Vinf = 0.9, 1.03, 1.0
    Ze = complex(Re, Xe)

    def p_borne(th):
        V = Vt * np.exp(1j * th)
        return (V * np.conj((V - Vinf) / Ze)).real

    th = brentq(lambda t: p_borne(t) - P, 0.0, 1.5)
    V = Vt * np.exp(1j * th)
    I = (V - Vinf) / Ze
    delta0 = np.angle(V + complex(Ra, Xq) * I)
    idq = I * np.exp(-1j * (delta0 - math.pi / 2))
    vdq = V * np.exp(-1j * (delta0 - math.pi / 2))
    Eqp0 = vdq.imag + Ra * idq.imag + Xdp * idq.real

    def stator(delta, Eqp):
        a = Vinf * np.exp(-1j * (delta - math.pi / 2))
        M = np.array([[Re + Ra, -Xe - Xq], [Xe + Xdp, Re + Ra]])
        i_d, i_q = np.linalg.solve(M, [-a.real, Eqp - a.imag])
        vd = a.real + Re * i_d - Xe * i_q
        vq = a.imag + Xe * i_d + Re * i_q
        Te = vd * i_d + vq * i_q + Ra * (i_d**2 + i_q**2)
        return np.array([Te, i_d, math.hypot(vd, vq)])

    h = 1e-6
    dd = (stator(delta0 + h, Eqp0) - stator(delta0 - h, Eqp0)) / (2 * h)
    dE = (stator(delta0, Eqp0 + h) - stator(delta0, Eqp0 - h)) / (2 * h)
    K2, K3, K6 = dE[0], 1 / (1 + (Xd - Xdp) * dE[1]), dE[2]
    assert dd[0] > 0 and K2 > 0 and K6 > 0  # point de fonctionnement normal
    s = 2j * math.pi * f
    exc = (p["Ka"] * p["Ke"]) / (
        (1 + s * p["Ta"]) * (1 + s * p["Te"])
        + p["Ka"] * p["Ke"] * Kfd * s / (1 + s * p["Tfd"])
    )
    return K2 * K3 * exc / (1 + s * K3 * Td0p + K3 * K6 * exc / (1 + s * p["Tr"]))


def test_heffron_phillips_retrouve_le_retard_et_sa_cause(smib_sans_pss):
    """Deux modeles, un resultat : a 1,6816 Hz, GEP est en retard de plus de
    90 degres AVEC le retour de stabilisation de l'AVR (-91,4 degres a un
    axe, -102,6 complet) et de moins de 40 degres SANS (-26,8 et -32,7). Les
    deux modeles different de 11 degres au plus, et de 2,5 % en module : le
    modele a un axe n'a ni amortisseurs, ni dynamique du stator et du
    reseau, ni charges de ligne -- d'ou une marge de 15 degres et 5 %."""
    net, ss, lam = smib_sans_pss
    f = lam.imag / (2 * np.pi)
    s = 2j * np.pi * f
    g_complet = gep(ss, s, "SM1", _machines(net))
    net0, ss0 = _reseau("smib", K_PSS=0.0, Kfd=0.0)
    g_complet0 = gep(ss0, s, "SM1", _machines(net0))
    g_hp = _gep_heffron_phillips(f, SM_DEFAULT["Kfd"])
    g_hp0 = _gep_heffron_phillips(f, 0.0)

    def deg(x):
        return math.degrees(np.angle(x))

    assert deg(g_hp) == pytest.approx(-91.401, abs=1e-3)
    assert deg(g_hp0) == pytest.approx(-26.786, abs=1e-3)
    assert deg(g_complet0) == pytest.approx(-32.743, abs=1e-3)
    assert abs(g_hp) == pytest.approx(abs(g_complet), rel=0.05)
    for complet, hp in ((g_complet, g_hp), (g_complet0, g_hp0)):
        assert abs(deg(complet) - deg(hp)) < 15.0
    assert deg(g_complet) < -90.0 and deg(g_hp) < -90.0
    assert deg(g_complet0) > -40.0 and deg(g_hp0) > -40.0


def test_sensibilite_au_gain_signe_de_gep(smib_sans_pss):
    """dlambda/dK_PSS en K = 0 sur le mode local du SMIB, trois facons :
    differences finies (pas 1e-3), residu exact (vecteurs propres a gauche et
    a droite, boucle v_ref -> Dwr fermee par le PSS), et approximation de
    l'equation de mouvement ``-GEP.PSS/(4H)``. Les deux premieres coincident
    a 1e-3 (erreur de troncature et bruit de ``eig``) ; la troisieme donne
    le meme signe aux deux parties, a 30 % pres en module (elle ignore le
    couplage du mode aux autres etats). Partie reelle POSITIVE : un gain
    positif pousse le mode vers l'instabilite."""
    net, ss, lam = smib_sans_pss
    dK = 1e-3
    proches = []
    for K in (dK, -dK):
        _n, s_k = _reseau("smib", K_PSS=K)
        ev = np.linalg.eigvals(s_k.A)
        proches.append(ev[np.argmin(np.abs(ev - lam))])
    d_fd = (proches[0] - proches[1]) / (2 * dK)

    ev, VL, VR = sl.eig(ss.A, left=True, right=True)
    k = np.argmin(np.abs(ev - lam))
    v, u = VR[:, k], VL[:, k].conj()
    b = ss.B[:, ss.iu("SM1.v_ref")]
    residu = v[ss.states.index("SM1.Dwr")] * (u @ b) / (u @ v)
    pss1 = pss_analytique({**SM_DEFAULT, "K_PSS": 1.0}, lam)
    d_res = residu * pss1

    sm = next(e for e in net.elements if isinstance(e, SynchronousMachine))
    d_gep = -gep(ss, lam, "SM1", _machines(net)) * pss1 / (4 * sm.pdict["H"])

    assert d_fd == pytest.approx(d_res, rel=1e-3)
    assert d_fd == pytest.approx(0.03519 + 0.03773j, rel=1e-3)
    assert d_fd.real > 0 and d_gep.real > 0
    assert d_fd.imag > 0 and d_gep.imag > 0
    assert abs(d_gep - d_fd) / abs(d_fd) < 0.30


#: mode le moins amorti (f [Hz], zeta [%]) a K_PSS = 0 puis 20 (fiche FR-007)
REFERENCES_FR007 = {
    "smib": (1.6816, 19.054, 1.7952, 13.106),
    "wscc9": (2.0599, 14.362, 2.2000, 8.696),
    "kundur2area": (1.9838, 14.623, 2.1172, 9.496),
}


@pytest.mark.parametrize("nom", RESEAUX)
def test_pss_par_defaut_desamortit_le_mode_local(nom):
    """Comportement ACTUEL, chiffre (fiche FR-007) : le PSS par defaut fait
    perdre 5 a 6 points d'amortissement au mode le moins amorti."""
    f0, z0, f20, z20 = REFERENCES_FR007[nom]
    for K, f_ref, z_ref in ((0.0, f0, z0), (20.0, f20, z20)):
        f, z = _f_zeta(_moins_amorti(_reseau(nom, K_PSS=K)[1]))
        assert f == pytest.approx(f_ref, abs=5e-5)
        assert z == pytest.approx(z_ref, abs=5e-4)


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="FR-007 : reglage par defaut inadapte (PSS sans avance face a un "
    "GEP en retard de plus de 90 degres), pas une erreur de signe ; "
    "correction a decider",
)
@pytest.mark.parametrize("nom", RESEAUX)
def test_pss_par_defaut_devrait_amortir_le_mode_local(nom):
    """Comportement ATTENDU d'un stabilisateur : il ne fait pas perdre
    d'amortissement au mode le moins amorti du reseau."""
    z0 = _f_zeta(_moins_amorti(_reseau(nom, K_PSS=0.0)[1]))[1]
    z20 = _f_zeta(_moins_amorti(_reseau(nom)[1]))[1]
    assert z20 > z0


@pytest.mark.parametrize("nom", RESEAUX)
def test_sans_retour_de_stabilisation_le_pss_par_defaut_amortit(nom):
    """Recoupement de la cause : avec ``Kfd = 0``, le retard de GEP retombe
    sous 60 degres a la frequence du mode, et le meme PSS amortit (SMIB
    13,9 -> 21,5 %, wscc9 13,6 -> 14,0 %, kundur2area 13,3 -> 19,2 %)."""
    net, ss = _reseau(nom, K_PSS=0.0, Kfd=0.0)
    lam0 = _moins_amorti(ss)
    s = 1j * lam0.imag
    for m in _machines(net):
        phase = math.degrees(
            np.angle(gep(ss, s, m, _machines(net)) * pss_analytique(SM_DEFAULT, s))
        )
        assert -60.0 < phase < 0.0, (m, phase)
    z20 = _f_zeta(_moins_amorti(_reseau(nom, Kfd=0.0)[1]))[1]
    assert z20 > _f_zeta(lam0)[1]


# ----------------------------------------------------------------------
# 3. ANDES, regulateurs actifs, Kundur 12.6
# ----------------------------------------------------------------------
class MachineSansFemVitesse(SynchronousMachine):
    """La machine du paquet, stator ecrit sans la vitesse : ``w_ref.phi`` au
    lieu de ``wr.phi``, comme le ``Flux0`` du GENROU d'ANDES. Sert a
    demontrer d'ou vient l'ecart du mode inter-zones, rien d'autre."""

    @staticmethod
    def equations(s):
        f, g, h = SynchronousMachine.equations(s)
        f = list(f)
        # f[0] = wb.(ved + Ra.igd + wr.phi_q), f[1] = wb.(veq + Ra.igq - wr.phi_d)
        f[0] = f[0] - s.wb * s.Dwr * s.phi_q
        f[1] = f[1] + s.wb * s.Dwr * s.phi_d
        return f, g, h


def _cas_andes():
    from g2egrid.andes_croise import AVR_PSS, SANS_REGULATEURS

    return {
        "sans": (SANS_REGULATEURS, False),
        "avr": ({**AVR_PSS, "K_PSS": 0.0}, True),
        "avr_pss": (AVR_PSS, True),
        "avr_pss_inverse": ({**AVR_PSS, "K_PSS": -20.0}, True),
    }


CAS_ANDES = ("sans", "avr", "avr_pss", "avr_pss_inverse")


@pytest.fixture(scope="module")
def croisement():
    from g2egrid.andes_croise import comparer, kundur_deux_zones

    return {
        nom: comparer(kundur_deux_zones(reglage=reglage), regulateurs=reg)
        for nom, (reglage, reg) in _cas_andes().items()
    }


@sans_andes
@pytest.mark.parametrize("cas", CAS_ANDES)
def test_andes_memes_modes(croisement, cas):
    """Trois modes electromecaniques des deux cotes (inter-zones, deux
    locaux), frequences a 1,5 % (mesure : 1,06 % au plus, AVR seul -- GENROU
    contre modele de Park)."""
    r = croisement[cas]
    assert r["n_andes"] == r["n_g2egrid"] == 3
    for p in r["modes"]:
        assert p["g2egrid"][0] == pytest.approx(p["andes"][0], rel=1.5e-2)


@sans_andes
@pytest.mark.parametrize("cas", CAS_ANDES)
def test_andes_modes_locaux_a_moins_d_un_point(croisement, cas):
    """Critere FR-007 : amortissement des modes locaux a moins de 1 point
    (mesure : 0,01 / 0,10 / 0,27 / 0,47 point selon le cas). Le mode
    inter-zones, a moins de 2 points (mesure : 0,35 / 1,48 / 1,67 / 1,12),
    ecart explique par le test suivant."""
    inter, *locaux = croisement[cas]["modes"]
    for p in locaux:
        assert abs(p["g2egrid"][1] - p["andes"][1]) < 1.0, p
    assert abs(inter["g2egrid"][1] - inter["andes"][1]) < 2.0, inter


@sans_andes
def test_andes_confirme_le_desamortissement(croisement):
    """Les deux outils voient la meme chose. Le PSS par defaut fait passer
    les modes locaux d'environ 4,5 % a moins de 2 %, et amortit fortement le
    mode inter-zones (vers 24 %) ; inverse, il destabilise les trois. Si le
    signe de la chaine differait de la convention IEEE (celle d'ANDES), le
    cas ``avr_pss`` de ce paquet ressemblerait au cas inverse d'ANDES."""
    for outil in ("andes", "g2egrid"):

        def zeta(cas, outil=outil):
            return [p[outil][1] for p in croisement[cas]["modes"]]

        inter_avr, *loc_avr = zeta("avr")
        inter_pss, *loc_pss = zeta("avr_pss")
        assert all(z < 2.0 for z in loc_pss) and all(z > 4.0 for z in loc_avr), outil
        assert inter_pss > 20.0 > inter_avr, outil
        assert all(z < 0.0 for z in zeta("avr_pss_inverse")), outil


@sans_andes
@pytest.mark.parametrize("cas", ["sans", "avr", "avr_pss", "avr_pss_inverse"])
def test_ecart_inter_zones_vient_de_la_fem_de_vitesse(monkeypatch, cas):
    """Le GENROU d'ANDES ecrit son stator sans la vitesse. Retirer ici la
    f.e.m. de vitesse ``Dwr.phi`` ramene TOUS les modes, inter-zones compris,
    a moins de 0,5 point d'ANDES (mesure : 0,13 / 0,14 / 0,20 / 0,31). L'AVR
    amplifie ce terme, puisqu'il agit sur la tension qu'il regule : c'est
    pourquoi l'ecart inter-zones passe de 0,35 point sans regulateur a
    1,5 point avec."""
    from g2egrid import andes_croise

    def assembler_sans_fem(t):
        net = assembler(t)
        for e in net.elements:
            if type(e) is SynchronousMachine:
                e.__class__ = MachineSansFemVitesse
        return net

    monkeypatch.setattr(andes_croise, "assembler", assembler_sans_fem)
    reglage, reg = _cas_andes()[cas]
    r = andes_croise.comparer(
        andes_croise.kundur_deux_zones(reglage=reglage), regulateurs=reg
    )
    assert r["n_andes"] == r["n_g2egrid"] == 3
    for p in r["modes"]:
        assert abs(p["g2egrid"][1] - p["andes"][1]) < 0.5, p


@sans_andes
def test_andes_refuse_ce_qu_il_ne_traduit_pas():
    """Deux garde-fous de la traduction : le GENROU impose X''d = X''q (la
    machine par defaut a 0,230 et 0,250), et le regulateur de vitesse n'a pas
    d'homologue."""
    from g2egrid.andes_croise import kundur_deux_zones, vers_andes

    with pytest.raises(ValueError, match="X''d"):
        vers_andes(charger("kundur2area"), regulateurs=True)
    with pytest.raises(ValueError, match="regulateur de vitesse"):
        vers_andes(kundur_deux_zones(reglage={}), regulateurs=True)
