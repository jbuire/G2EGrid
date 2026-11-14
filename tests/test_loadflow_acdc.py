"""Tests du LOAD FLOW AC/DC UNIFIE et des tables continues (lot 4).

L'ordre de ce fichier est celui dans lequel les tests ont ete ecrits, et ce
n'est pas un detail de presentation :

1. **le jacobien contre les differences finies**, avant tout le reste. Un
   jacobien faux ne fait pas diverger un Newton, il le fait converger plus
   lentement : sur les reseaux simples on ne voit rien, et l'erreur ressort
   sur le cas qui compte, six heures plus tard ;
2. **la non-regression** du load flow alternatif sur tous les reseaux
   alternatifs du catalogue. Le lot 3 pinait le TEXTE de ``loadflow.py`` ;
   le lot 4 doit y ecrire, donc la garantie change de nature : ce qui est
   preserve n'est plus le texte du load flow alternatif mais son RESULTAT.
   Ce resultat a d'abord ete compare octet pour octet. Depuis FR-001, la
   comparaison se fait a ``TOL_V = 1e-11`` pu, au meme nombre d'iterations,
   parce que l'octet ne survit pas d'une machine a l'autre. Un second filet,
   independant, recalcule les equations du load flow sans aucune reference ;
3. le bilan de puissance a deux puis trois terminaux, PERTES RECALCULEES A
   PART -- sinon le test est une tautologie ;
4. ``residu_equilibre``, qui dit si le point de fonctionnement mixte est
   vraiment un point fixe du modele non lineaire ;
5. les erreurs explicites : ile continue sans reference, ile coupee, mode
   impossible pour la classe.

Ajoutees ensuite, pour FR-001 :

10. **la non-regression du load flow MIXTE** sur les six reseaux continus du
    catalogue, avec les deux memes filets que l'alternatif : une reference
    (``reference_loadflow_acdc.json``, V et Vdc a ``TOL_V``) et les
    equations recalculees, bilan de chaque station compris ;
11. **une limite connue**, en ``xfail`` strict : la tolerance d'arret du
    Newton mixte (1e-12) est sous le plancher d'arrondi d'un cable continu
    court, et le load flow y leve « non converge » au point fixe.
"""

import base64
import functools
import json
import pathlib

import numpy as np
import pytest

from g2egrid import Network
from g2egrid.catalogue import RESEAUX, charger
from g2egrid.element import InitContext
from g2egrid.elements import (
    GFL2DC,
    GFM3DC,
    GFM4DC,
    IB,
    BusDC,
    ChargeDC_P,
    LineDC,
    SourceDC,
)
from g2egrid.elements.mmc import MMC
from g2egrid.loadflow import (
    DC_P,
    DC_SLACK,
    MODES_VSC,
    PQ,
    PV,
    SLACK,
    ProblemeACDC,
    run_loadflow,
)
from g2egrid.tables import (
    TABLES,
    Tables,
    assembler,
    tables_depuis_json,
    tables_vers_json,
)
from g2egrid.temporel import ModeleNonLineaire

REFERENCE = pathlib.Path(__file__).resolve().parent / "reference_loadflow.json"
REFERENCE_ACDC = (
    pathlib.Path(__file__).resolve().parent / "reference_loadflow_acdc.json"
)


# ----------------------------------------------------------------------
# fixtures : les deux reseaux du lot, charges depuis leurs notebooks
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def hvdc2t() -> Network:
    """Liaison a deux terminaux, load flow resolu (pas d'equilibre)."""
    net = assembler(charger("hvdc2t"))
    net.run_loadflow()
    return net


@pytest.fixture(scope="module")
def mtdc3t() -> Network:
    """Reseau continu a trois terminaux, load flow resolu."""
    net = assembler(charger("mtdc3t"))
    net.run_loadflow()
    return net


@pytest.fixture(scope="module")
def hvdc2t_resolu() -> Network:
    """Liaison a deux terminaux, equilibre par element compris."""
    net = assembler(charger("hvdc2t"))
    net.solve_equilibrium()
    return net


@pytest.fixture(scope="module")
def mtdc3t_resolu() -> Network:
    net = assembler(charger("mtdc3t"))
    net.solve_equilibrium()
    return net


def pertes_convertisseur(el, etats=None, vdc=None) -> float:
    """``kS.(R1.|is|^2 + R2.|ig|^2 + Gdc.vdc^2)``, en pu base systeme.

    LES PERTES SONT RECALCULEES DEPUIS LES ETATS DU MODELE, et n'empruntent
    rien au load flow. C'est ce qui fait du bilan de puissance un test et non
    une identite : si le couplage avait pose un rendement au lieu du vrai
    courant de port, les deux comptes ne tomberaient pas ensemble.

    Le ``R2.|ig|^2`` n'est pas un ajout de confort : entre le condensateur du
    filtre, ou le modele mesure ``p``, et le noeud du reseau, il y a la branche
    R2-L2. La puissance qu'un VSC pose sur son BUS est donc celle du
    condensateur moins cette perte-la, et le bilan noeud-a-noeud la reclame.

    Un :class:`~g2egrid.elements.mmc.MMC` n'a ni filtre LCL ni condensateur
    continu : ses pertes sont celles de ``R1`` cote alternatif et de ``Rdc``
    (bras) cote continu, ``kS.(R1.|is|^2 + Rdc.idc^2)``.

    Par defaut, les etats et la tension continue sont ceux de l'equilibre
    (``el.x0``, ``el.u0``) ; ``etats`` et ``vdc`` permettent de les prendre
    ailleurs, par exemple au point du load flow (:func:`etats_station`).
    """
    p = el.pdict
    x = dict(zip(el.states, el.x0)) if etats is None else etats
    vdc = el.u0[el.port_dc_input_idx(0)] if vdc is None else vdc
    if isinstance(el, MMC):
        return p["kS"] * (
            p["R1_pu"] * (x["isd"] ** 2 + x["isq"] ** 2) + p["Rdc_pu"] * x["idc"] ** 2
        )
    if not isinstance(el, (GFM3DC, GFM4DC, GFL2DC)):
        pytest.fail(
            f"{el.name}: pertes de {type(el).__name__} inconnues ; les ecrire "
            "dans pertes_convertisseur avant de la mettre sous le filet"
        )
    return p["kS"] * (
        p["R1_pu"] * (x["isd"] ** 2 + x["isq"] ** 2)
        + p["R2_pu"] * (x["igd"] ** 2 + x["igq"] ** 2)
        + p["Gdc_pu"] * vdc**2
    )


def bilan(net: Network) -> dict[str, float]:
    """``somme(P_ac) + pertes convertisseurs + pertes cables``, et ses termes."""
    lf = net.lf
    idx = net._bus_index()
    prog = net._injections_programmees()
    p_ac = {}
    pertes = {}
    for el in net._elements_vsc():
        p_ac[el.name] = float(lf.S[idx[el.buses[0]]].real - prog[el.buses[0]][0])
        pertes[el.name] = pertes_convertisseur(el)
    pos = {b: j for j, b in enumerate(sorted(net.bus_specs_dc))}
    G = net._reseau_dc(pos)[0]
    cables = float(lf.Vdc @ (G @ lf.Vdc))
    return {
        "p_ac": p_ac,
        "pertes": pertes,
        "cables": cables,
        "bilan": sum(p_ac.values()) + sum(pertes.values()) + cables,
    }


# ======================================================================
# 1) LE JACOBIEN CONTRE LES DIFFERENCES FINIES -- en premier
# ======================================================================
def jacobien_numerique(pb, x, h_rel=1e-7) -> np.ndarray:
    """``-dF/dx`` par differences CENTREES, la convention de ``pb.jacobien``.

    Centrees et non avancees : l'erreur de troncature tombe en ``h^2``, ce qui
    laisse assez de marge entre elle et l'arrondi pour que l'ecart mesure soit
    celui du jacobien, et non celui de la methode qui le mesure.
    """
    J = np.zeros((pb.n_f, pb.n_x))
    for k in range(pb.n_x):
        h = h_rel * max(1.0, abs(x[k]))
        xp, xm = x.copy(), x.copy()
        xp[k] += h
        xm[k] -= h
        J[:, k] = -(pb.residu(xp) - pb.residu(xm)) / (2.0 * h)
    return J


@pytest.mark.parametrize("nom", ["hvdc2t", "mtdc3t"])
def test_jacobien_contre_differences_finies(nom):
    """Terme a terme, le jacobien analytique EST la derivee du residu.

    Eprouve sur le vrai reseau, et pas seulement au point de depart : un
    jacobien peut etre juste a plat (ou les angles sont nuls et les sinus
    aussi) et faux partout ailleurs. On le mesure donc en trois points --
    l'amorce, un point perturbe au hasard, et la solution.

    Ce test couvre les quatre blocs a la fois : ``J1..J4`` du cote alternatif,
    la colonne ``-1`` du P des VSC en mode ``Vdc...``, les derivees du
    couplage, et la regle de chaine ``dP_dc/dq . dQ/dx`` par laquelle un
    grid-forming sur bus PV fait dependre son bilan continu des angles
    alternatifs. C'est cette derniere qui manquerait le plus discretement.
    """
    net = assembler(charger(nom))
    pb = net._probleme_acdc()
    rng = np.random.default_rng(12345)
    depart = pb.depart()
    solution = pb.resoudre().V  # pour verifier que la resolution passe
    assert solution is not None

    pires = []
    for etiquette, x in (
        ("amorce", depart),
        ("perturbe", depart + 0.02 * rng.standard_normal(pb.n_x)),
        ("solution", _x_solution(net, pb)),
    ):
        Ja = pb.jacobien(x)
        Jn = jacobien_numerique(pb, x)
        echelle = max(float(np.abs(Jn).max()), 1e-12)
        err = float(np.abs(Ja - Jn).max()) / echelle
        pires.append((etiquette, err))
        assert err < 1e-6, f"{nom} au point '{etiquette}' : ecart relatif {err:.3e}"
    print(f"\n{nom} : ecart jacobien/DF " + ", ".join(f"{e}={v:.2e}" for e, v in pires))


def _x_solution(net: Network, pb) -> np.ndarray:
    """Le vecteur d'inconnues a la solution du load flow."""
    lf = net.run_loadflow()
    th = np.angle(lf.V)
    V = np.abs(lf.V)
    return np.concatenate((th[pb.idx_th], V[pb.idx_v], lf.Vdc[pb.idx_dcv], lf.P_vsc))


def test_jacobien_est_bien_conditionne(hvdc2t):
    """Le systeme est carre et inversible : ni inconnue en trop, ni equation."""
    pb = hvdc2t._probleme_acdc()
    assert pb.n_x == pb.n_f
    J = pb.jacobien(pb.depart())
    assert np.linalg.cond(J) < 1e10, f"conditionnement {np.linalg.cond(J):.2e}"


# ======================================================================
# 2) LA NON REGRESSION DU LOAD FLOW ALTERNATIF : deux filets independants
# ======================================================================
# Le PREMIER compare le resultat a une reference stockee : il voit tout ce qui
# deplace le point de fonctionnement, que le deplacement soit juste ou faux,
# et, par le nombre d'iterations, un jacobien faux qui ne fait que ralentir le
# Newton. Le SECOND n'a besoin d'aucune reference : il recalcule les equations
# du load flow depuis le Y-bus et la solution. Il voit un point qui n'en est
# pas un, meme si l'on a regenere la reference dessus. Mais il recalcule avec
# le MEME ``build_ybus`` que le load flow : il garantit le Newton, pas la
# matrice. Une erreur dans le Y-bus (signe d'une resistance, admittance
# serie) n'est vue que par le premier filet -- tant que la reference n'est pas
# regeneree dessus -- et par ``tests/test_interop.py``, qui compare a
# MATPOWER (case9, case57). Revue de FR-001 : 166 mutants sur 166 a R
# inverse, 155 sur 159 a admittance serie inversee passaient le second filet.
#
# CE QUI ETAIT GARANTI ET NE L'EST PLUS. Jusqu'a FR-001, le premier filet
# comparait ``lf.V.tobytes()`` a la reference, octet pour octet. Il attrapait
# une expression deplacee qui ne change que l'arrondi, mais seulement sur la
# machine qui avait produit la reference. Sur une autre (Linux, NumPy 2.5.3,
# OpenBLAS), 18 des 23 reseaux s'en ecartaient de 1e-17 a 1,8e-13 pu, au meme
# nombre d'iterations : c'est l'arrondi de BLAS/LAPACK et des sin/cos
# vectorises de NumPy, pas un changement de modele. Une tolerance ne voit plus
# une expression deplacee de 1 ulp. C'est le prix, assume, d'un filet qui
# passe sur toutes les machines.


#: Les reseaux du catalogue qui sont PUREMENT ALTERNATIFS, donc compares a
#: ``reference_loadflow.json``. La liste se deduit des tables plutot que de
#: s'ecrire a la main : ajouter un reseau alternatif au catalogue le fait
#: entrer d'office sous le filet, et c'est exactement ce qu'on veut.
def _alternatifs() -> list[str]:
    from g2egrid.catalogue import charger

    out = []
    for nom in RESEAUX:
        t = charger(nom)
        if not (t.T_BusDC or t.T_VSCs or t.T_LinesDC):
            out.append(nom)
    return out


ALTERNATIFS = _alternatifs()
#: Les autres : ils ont leur propre reference (section 10).
CONTINUS = [nom for nom in RESEAUX if nom not in ALTERNATIFS]

#: ECART ADMIS ENTRE UN LOAD FLOW ET SA REFERENCE, en pu : le module de la
#: difference des phaseurs, au pire bus (et, pour les reseaux continus, au
#: pire noeud continu). Cet ecart vaut :
#:
#: - 55 fois le bruit MESURE d'une machine a l'autre (1,8e-13, ieee39_mixte) ;
#: - 28 fois le bruit EMULE au pire, en perturbant de 4 ulp chaque coefficient
#:   du jacobien (3,6e-13, cigre_mv_dc). Sur les reseaux alternatifs, ce bruit
#:   emule vaut 1 a 2 fois le bruit mesure, et il sert de mesure la ou l'on
#:   n'a qu'une machine (les reseaux continus) ;
#: - 400 fois moins que le plus petit effet d'une consigne PV deplacee de
#:   1e-6 (4,1e-9 pu, consigne de puissance de nordic44) : voir
#:   :func:`test_le_filet_voit_une_consigne_pv_deplacee`.
#:
#: Ce qui passe en dessous : une consigne de puissance deplacee de moins
#: d'environ 3e-9 pu (sur nordic44, 1e-9 pu ne deplace V que de 4e-12 pu).
TOL_V = 1e-11

#: Le seuil auquel on COMPTE les iterations : la tolerance d'arret par defaut
#: de :func:`~g2egrid.loadflow.run_loadflow`. Le compte est portable : de part
#: et d'autre de 1e-10, le mismatch de chaque reseau alternatif reste a plus
#: de 6e-11 du seuil (wscc3sm_* : 1,6e-10 a l'avant-derniere iteration), pour
#: un plancher d'arrondi d'au plus 1e-13 (aucun compte n'a bouge sous le bruit
#: emule).
TOL_ITERATIONS = 1e-10


@functools.cache
def _resolu(nom: str):
    """Un reseau du catalogue et son load flow, calcules une fois par session.

    Les deux filets les partagent et ne les modifient pas. Un temoin de
    mutation construit TOUJOURS son propre reseau.
    """
    net = assembler(charger(nom))
    return net, net.run_loadflow()


def _reference(nom: str) -> dict:
    return json.loads(REFERENCE.read_text(encoding="utf-8"))[nom]


def ecarts_alternatifs(lf, attendu: dict) -> list[str]:
    """Ce qui separe un load flow de sa reference : une liste vide si rien.

    Une liste et non des assertions : les temoins de mutation verifient que
    cette fonction PARLE quand il le faut, et le filet emploie la meme, pas une
    copie. On ecrit ``not dV <= TOL_V`` et non ``dV > TOL_V``, pour qu'un NaN
    echoue.
    """
    if lf.V.size != attendu["n_bus"]:
        return [f"{lf.V.size} bus au lieu de {attendu['n_bus']}"]
    V_ref = np.frombuffer(base64.b64decode(attendu["V"]), dtype=complex)
    out = []
    dV = float(np.abs(lf.V - V_ref).max())
    if not dV <= TOL_V:
        out.append(f"max|dV| = {dV:.3e} pu > {TOL_V:.0e}")
    if lf.iterations != attendu["iterations"]:
        out.append(f"{lf.iterations} iterations au lieu de {attendu['iterations']}")
    return out


def injections_polaires(Y: np.ndarray, V: np.ndarray):
    """``(P, Q)`` injectes a chaque bus, en ecriture POLAIRE.

    ``P_i = V_i.somme_k V_k.(G_ik.cos th_ik + B_ik.sin th_ik)``,
    ``Q_i = V_i.somme_k V_k.(G_ik.sin th_ik - B_ik.cos th_ik)``.

    Le load flow, lui, calcule ``V.conj(Y.V)`` en complexe. Ce sont deux
    chemins d'arrondi distincts : une erreur de signe ou d'indice dans l'un ne
    se retrouverait pas dans l'autre.
    """
    Vm, th = np.abs(V), np.angle(V)
    G, B = Y.real, Y.imag
    d = th[:, None] - th[None, :]
    P = Vm * ((G * np.cos(d) + B * np.sin(d)) @ Vm)
    Q = Vm * ((G * np.sin(d) - B * np.cos(d)) @ Vm)
    return P, Q


def residus_alternatifs(net: Network, lf) -> dict[str, float]:
    """Les equations du load flow, recalculees ICI, bus par bus.

    Ni ``lf.mismatch`` ni ``lf.S`` : seulement le Y-bus et ``lf.V``, contre
    les consignes de ``net.bus_specs``. La fonction rend le pire ecart de
    chaque famille d'equations, plus les pertes du reseau (``somme(P)``).
    """
    P, Q = injections_polaires(net.build_ybus(), lf.V)
    Vm, th = np.abs(lf.V), np.angle(lf.V)
    r = dict.fromkeys(("P", "Q", "V impose", "theta impose"), 0.0)
    for b, i in net._bus_index().items():
        s = net.bus_specs[b]
        if s.btype in (PV, PQ):
            r["P"] = max(r["P"], abs(P[i] - (s.Pg - s.Pl)))
        if s.btype == PQ:
            r["Q"] = max(r["Q"], abs(Q[i] - (s.Qg - s.Ql)))
        if s.btype in (SLACK, PV):
            r["V impose"] = max(r["V impose"], abs(Vm[i] - s.V))
        if s.btype == SLACK:
            ecart = abs(np.angle(np.exp(1j * (th[i] - s.theta))))
            r["theta impose"] = max(r["theta impose"], ecart)
    r["pertes"] = float(P.sum())
    return r


def test_reference_couvre_tous_les_reseaux_alternatifs():
    """La reference n'a pas silencieusement perdu un reseau.

    Le controle se fait dans LES DEUX SENS, et le second est le vrai garde-fou :
    un reseau alternatif ajoute au catalogue sans entree de reference
    echapperait a la non-regression sans que rien ne le signale. Un reseau
    CONTINU, lui, n'a pas a y figurer : il a sa propre reference (section 10).
    """
    ref = json.loads(REFERENCE.read_text(encoding="utf-8"))
    manquants = set(ALTERNATIFS) - set(ref)
    assert not manquants, f"reseaux alternatifs hors reference : {sorted(manquants)}"
    disparus = set(ref) - set(RESEAUX)
    assert not disparus, (
        f"reseaux de la reference sortis du catalogue : {sorted(disparus)}"
    )


def test_les_reseaux_continus_ont_leur_propre_reference():
    """Et on sait lesquels : la distinction doit rester explicite.

    Ils ne figurent pas dans ``reference_loadflow.json``, qui date d'avant le
    load flow mixte. Ils sont dans ``reference_loadflow_acdc.json``, qui doit
    les couvrir exactement, dans les deux sens, pour la meme raison que
    ci-dessus.
    """
    assert set(CONTINUS) == {
        "hvdc2t",
        "mtdc3t",
        "cigre_b4",
        "cigre_b4_dcs1",
        "cigre_b4_mmc",
        "cigre_mv_dc",
    }, CONTINUS
    ref = json.loads(REFERENCE_ACDC.read_text(encoding="utf-8"))["reseaux"]
    assert set(ref) == set(CONTINUS), set(ref) ^ set(CONTINUS)


@pytest.mark.parametrize("nom", ALTERNATIFS)
def test_load_flow_alternatif_contre_reference(nom):
    """LE PREMIER FILET : ``max|V - V_ref| <= TOL_V``, et le meme nombre
    d'iterations.

    Le lot 3 gardait le TEXTE de ``loadflow.py`` et de ``Network.run_loadflow``
    fige. Le lot 4 a du y ecrire, ne serait-ce que la garde
    ``if self.bus_specs_dc:``, et la garantie s'est deplacee du texte vers le
    resultat. FR-001 l'a ensuite fait passer de l'octet a ``TOL_V``, pour
    qu'elle tienne d'une machine a l'autre (voir l'en-tete de la section).

    Les ITERATIONS se comparent a part, et pas apres V : l'ancien test les
    verifiait apres ``tobytes()``. Un ecart de V masquait donc celui des
    iterations, et ``wscc9_gfm`` (5 iterations au lieu de 6) n'a jamais ete
    signale.

    La reference a ete produite depuis un export ``git archive HEAD`` extrait a
    part, pas depuis l'arbre de travail. Seule celle de ``wscc9_gfm`` a ete
    regeneree pour FR-001, depuis un arbre commite : ses bus GFM sont PV depuis
    ``d263b99`` (convention GFM a batterie). L'outil est
    ``developpement/regenerer_reference_loadflow.py``.
    """
    net, lf = _resolu(nom)
    assert not net.bus_specs_dc, f"{nom} n'est pas un reseau purement alternatif"
    ecarts = ecarts_alternatifs(lf, _reference(nom))
    assert not ecarts, f"{nom} : le load flow alternatif a bouge : " + " ; ".join(
        ecarts
    )


def test_run_loadflow_module_inchangee_sur_un_cas_direct():
    """La fonction de ``loadflow.py`` elle-meme, appelee sans passer par Network.

    Le test precedent passe par ``Network.run_loadflow`` ; celui-ci attaque la
    fonction de module directement, pour que la garantie ne depende pas du
    chemin d'appel.
    """
    net = assembler(charger("wscc9"))
    idx = net._bus_index()
    n = len(idx)
    Y = net.build_ybus()
    btype = np.zeros(n, dtype=int)
    Vsp, th = np.ones(n), np.zeros(n)
    Psp, Qsp = np.zeros(n), np.zeros(n)
    for b, i in idx.items():
        s = net.bus_specs[b]
        btype[i], Vsp[i], th[i] = s.btype, s.V, s.theta
        Psp[i], Qsp[i] = s.Pg - s.Pl, s.Qg - s.Ql
    lf = run_loadflow(Y, btype, Vsp, th, Psp, Qsp)
    ecarts = ecarts_alternatifs(lf, _reference("wscc9"))
    assert not ecarts, " ; ".join(ecarts)


@pytest.mark.parametrize("nom", ALTERNATIFS)
def test_load_flow_alternatif_resout_ses_equations(nom):
    """LE SECOND FILET : la solution satisfait ses equations, sans reference.

    - ``P`` aux bus PV et PQ, ``Q`` aux bus PQ, sous 1e-10 pu. C'est la
      promesse de :func:`run_loadflow` (``tol = 1e-10`` sur le mismatch),
      tenue ici par une autre ecriture. Mesure : 3,2e-11 au pire
      (wscc9_cigre6). L'ecriture polaire n'y ajoute que l'arrondi, quelques
      1e-13 au plus (|Y| <= 1,9e3 pu dans le catalogue).
    - ``|V|`` des bus PV et slack, et ``theta`` du slack, a 1e-14 pres. Ce ne
      sont pas des inconnues : le Newton ne les touche pas, et
      ``abs(V.exp(j.theta))`` ne les arrondit que d'un ulp (2,2e-16 mesure).
    - Les pertes ``somme(P)`` sont strictement positives. Le Y-bus est
      passif (resistances et conductances shunt positives), et tous les
      reseaux du catalogue ont des resistances : au moins 6,2e-5 pu mesure
      (smib_sm). Ce controle ne protege PAS ``build_ybus`` : une seule
      branche a resistance inversee laisse le total positif (voir l'en-tete
      de la section).
    """
    net, lf = _resolu(nom)
    r = residus_alternatifs(net, lf)
    assert r["P"] < 1e-10 and r["Q"] < 1e-10, (nom, r)
    assert r["V impose"] < 1e-14 and r["theta impose"] < 1e-14, (nom, r)
    assert r["pertes"] > 0.0, (nom, r)


@pytest.mark.parametrize("consigne", ["V", "Pg"])
@pytest.mark.parametrize("nom", ["wscc9", "nordic44"])
def test_le_filet_voit_une_consigne_pv_deplacee(nom, consigne):
    """TEMOIN DE MUTATION : une consigne PV deplacee de 1e-6 fait echouer le
    premier filet.

    Un filet qui ne voit rien n'echoue jamais : celui-ci doit prouver qu'il
    voit ce qu'on attend de lui. La consigne de tension se retrouve telle
    quelle dans V (1e-6 au bus, cent mille fois ``TOL_V``). Celle de puissance
    passe par le reseau, et ``nordic44`` est le reseau du catalogue ou elle
    deplace le moins V (4,1e-9 pu, mesure sur le premier bus PV de chaque
    reseau).

    Le second filet, lui, doit RESTER vert : le load flow mute est juste, il
    resout d'autres equations. C'est ce qui rend les deux filets independants.
    """
    net = assembler(charger(nom))
    b = min(b for b, s in net.bus_specs.items() if s.btype == PV)
    s = net.bus_specs[b]
    setattr(s, consigne, getattr(s, consigne) + 1e-6)
    lf = net.run_loadflow()
    ecarts = ecarts_alternatifs(lf, _reference(nom))
    assert any(e.startswith("max|dV|") for e in ecarts), (nom, b, consigne, ecarts)
    r = residus_alternatifs(net, lf)
    assert r["P"] < 1e-10 and r["Q"] < 1e-10, r


def test_le_filet_voit_un_jacobien_degrade(monkeypatch):
    """TEMOIN : un jacobien faux de 1 % ne deplace pas la solution, il ralentit
    le Newton, et le filet le voit par le nombre d'iterations SEULEMENT.

    Mesure sur ``wscc9`` : 7 iterations au lieu de 5, et V a 3,1e-12 pu de la
    reference, sous ``TOL_V``. Le compte des iterations n'est donc pas un
    detail du filet : sans lui, une erreur de derivee passerait.
    """
    net = assembler(charger("wscc9"))
    attendu = _reference("wscc9")
    resoudre = np.linalg.solve
    monkeypatch.setattr(np.linalg, "solve", lambda J, F: resoudre(1.01 * J, F))
    lf = net.run_loadflow()
    ecarts = ecarts_alternatifs(lf, attendu)
    assert ecarts == [
        f"{lf.iterations} iterations au lieu de {attendu['iterations']}"
    ], ecarts
    assert lf.iterations > attendu["iterations"]


# ======================================================================
# 3) LE BILAN A DEUX TERMINAUX
# ======================================================================
def test_bilan_deux_terminaux(hvdc2t_resolu):
    """``P_ac1 + P_ac2 + pertes = 0`` -- le coeur du lot, en un nombre.

    Les deux zones alternatives n'ont AUCUN lien a courant alternatif : tout ce
    qui sort de l'une par son convertisseur entre dans l'autre par le sien,
    moins ce que les deux convertisseurs et le cable dissipent. Si le load flow
    avait resolu le continu a part, ou avait pose un rendement, ce nombre ne
    serait pas nul -- il serait petit, ce qui est pire.
    """
    b = bilan(hvdc2t_resolu)
    assert len(b["p_ac"]) == 2
    assert abs(b["bilan"]) < 1e-10, (
        f"bilan {b['bilan']:.3e} ; P_ac {b['p_ac']}, "
        f"pertes {b['pertes']}, cables {b['cables']:.6e}"
    )
    # les deux convertisseurs travaillent bien en sens opposes
    valeurs = sorted(b["p_ac"].values())
    assert valeurs[0] < 0 < valeurs[1]
    # et les pertes sont strictement positives : elles ne peuvent pas etre nulles
    assert all(p > 0 for p in b["pertes"].values())
    assert b["cables"] > 0


def test_loi_des_noeuds_continue_a_deux_terminaux(hvdc2t_resolu):
    """A chaque noeud continu, la somme des courants de port est nulle.

    Le pendant, au niveau du CABLAGE, de ce que le bilan verifie au niveau des
    puissances : le load flow et le modele dynamique doivent tomber sur le meme
    point, pas seulement sur le meme total.
    """
    somme = {b: 0.0 for b in hvdc2t_resolu.bus_specs_dc}
    for el in hvdc2t_resolu.elements:
        if not el.ports_dc:
            continue
        y = el.compiled().h(el.x0, el.z0, el.u0, el.p)
        for pk in range(len(el.ports_dc)):
            somme[el.buses[len(el.ports) + pk]] += float(y[el.port_dc_output_idx(pk)])
    for b, i in somme.items():
        assert abs(i) < 1e-9, f"bus continu {b}: somme des courants {i:.3e}"


def test_tension_continue_tenue_par_le_vsc(hvdc2t_resolu):
    """Le noeud de type 1 est EXACTEMENT a sa consigne, l'autre s'en ecarte."""
    net = hvdc2t_resolu
    slack = [b for b, s in net.bus_specs_dc.items() if s.btype == DC_SLACK]
    libre = [b for b, s in net.bus_specs_dc.items() if s.btype == DC_P]
    assert len(slack) == 1 and len(libre) == 1
    assert net.v_dc[slack[0]] == pytest.approx(net.bus_specs_dc[slack[0]].V, abs=1e-12)
    # le courant circule du noeud libre vers le noeud tenu : il est plus haut
    assert net.v_dc[libre[0]] > net.v_dc[slack[0]]


# ======================================================================
# 4) LE MTDC A TROIS TERMINAUX
# ======================================================================
def test_bilan_trois_terminaux(mtdc3t_resolu):
    """Meme bilan, mais le reseau continu est un graphe et non un fil.

    Trois terminaux, trois modes differents (``PV``, ``VdcQ``, ``PQ``) et deux
    cables : si le laplacien des conductances etait mal monte -- un signe, un
    indice -- le bilan le dirait ici et pas a deux terminaux, ou toute matrice
    antisymetrique 2x2 a l'air juste.
    """
    b = bilan(mtdc3t_resolu)
    assert len(b["p_ac"]) == 3
    assert abs(b["bilan"]) < 1e-10, (
        f"bilan {b['bilan']:.3e} ; P_ac {b['p_ac']}, "
        f"pertes {b['pertes']}, cables {b['cables']:.6e}"
    )


def test_trois_modes_exerces(mtdc3t_resolu):
    """Le reseau a trois terminaux exerce bien les trois modes du lot."""
    assert sorted(mtdc3t_resolu.vsc_modes.values()) == ["PQ", "PV", "VdcQ"]


def test_puissance_programmee_respectee(mtdc3t_resolu):
    """Un VSC qui n'est pas en mode ``Vdc...`` tient sa puissance au watt pres.

    C'est la contrepartie du test precedent : le Newton a le droit de bouger la
    puissance du SEUL convertisseur qui tient la tension continue. S'il touchait
    aux autres, le bilan boucleraient quand meme -- et le reseau ne ferait plus
    ce qu'on lui a demande.
    """
    net = mtdc3t_resolu
    idx = net._bus_index()
    prog = net._injections_programmees()
    for el in net._elements_vsc():
        bus = el.buses[0]
        if net.vsc_modes[bus].startswith("Vdc"):
            continue
        s = net.bus_specs[bus]
        programme = (s.Pg - s.Pl) - prog[bus][0]
        obtenu = float(net.lf.S[idx[bus]].real - prog[bus][0])
        assert obtenu == pytest.approx(programme, abs=1e-10), el.name


# ======================================================================
# 5) LE POINT DE FONCTIONNEMENT EST UN POINT FIXE DU MODELE NON LINEAIRE
# ======================================================================
@pytest.mark.parametrize("fixture", ["hvdc2t_resolu", "mtdc3t_resolu"])
def test_residu_equilibre(fixture, request):
    """``residu_equilibre() < 1e-8`` : l'epreuve du lot.

    Elle ne mesure pas le load flow, elle mesure l'ACCORD entre le load flow et
    le modele dynamique complet. Un point de fonctionnement peut satisfaire un
    bilan de puissance a 1e-15 et n'etre pas un point fixe du modele : il
    suffit que la puissance continue attribuee a un convertisseur ne soit pas
    celle que ses treize equations rendent a cette tension-la. C'est
    exactement ce que la version Kelada ratait de 0,08 pu.
    """
    net = request.getfixturevalue(fixture)
    net.build_state_space()
    dx, dzu = ModeleNonLineaire(net).residu_equilibre()
    assert dx < 1e-8, f"{fixture} : |dx/dt| = {dx:.3e}"
    assert dzu < 1e-8, f"{fixture} : |z - z0| + |u - u0| = {dzu:.3e}"
    print(f"\n{fixture} : |dx/dt| = {dx:.3e}, dzu = {dzu:.3e}")


def test_liaison_deux_terminaux_est_stable(hvdc2t_resolu):
    """La liaison a deux terminaux n'a aucun mode a droite.

    Pas une exigence du load flow, mais la preuve qu'on livre un reseau de test
    utilisable : un reseau dont le point de fonctionnement est juste mais qui
    diverge en temporel ne sert a rien pour la suite.
    """
    ss = hvdc2t_resolu.build_state_space()
    assert ss.is_stable(), f"max Re(lambda) = {ss.eig().real.max():+.4e}"


def test_mtdc_trois_terminaux_est_stable(mtdc3t_resolu):
    """Idem a trois terminaux.

    C'est ici que la resistance des cables continus gagne sa valeur : le reseau
    continu a une resonance LC entre l'inductance des cables et la capacite des
    noeuds, que les terminaux a puissance constante DESAMORTISSENT (une charge
    a puissance constante a une conductance incrementale negative). Mesure : a
    0,02 ohm/km ce mode sort a +66 s^-1 et 426 Hz.
    """
    ss = mtdc3t_resolu.build_state_space()
    assert ss.is_stable(), f"max Re(lambda) = {ss.eig().real.max():+.4e}"


# ======================================================================
# 6) LES ILES CONTINUES : sans reference, ou coupees
# ======================================================================
def _tables_hvdc2t_modifiees(**modifs) -> Tables:
    return charger("hvdc2t").copie(**modifs)


def test_ile_continue_sans_reference():
    """Deux convertisseurs a puissance imposee : rien ne fixe le niveau.

    Vu du continu, un VSC est un terminal a PUISSANCE constante. Une ile qui
    n'en contient que de ceux-la n'a aucune equation qui fixe son niveau de
    tension, et les pertes des cables ne sont a la charge de personne. L'erreur
    doit le dire, et NOMMER l'ile.
    """
    t = _tables_hvdc2t_modifiees(
        T_BusDC=[(10, 2, 200e3, 1.00, 0.0), (11, 2, 200e3, 1.00, 0.0)],
        T_VSCs=[(1, 2, 10, "PV"), (2, 3, 11, "PQ")],
    )
    net = assembler(t)
    with pytest.raises(ValueError, match="aucune reference de tension") as exc:
        net.run_loadflow()
    assert "10" in str(exc.value) and "11" in str(exc.value)


def test_ile_continue_non_connexe():
    """Un reseau continu coupe en deux : la moitie orpheline est nommee.

    On ajoute un troisieme noeud continu et un troisieme convertisseur, SANS
    cable pour le relier au reste. Son ile n'a pas de reference, et le message
    doit designer cette ile-la -- pas l'autre, qui va tres bien.
    """
    t = _tables_hvdc2t_modifiees(
        Y_DER=list(charger("hvdc2t").Y_DER)
        + [
            (5, 3, 230e3, 1.00, 0.0, 25e6, 0, 0, 0),
            (6, 1, 230e3, 1.00, 0.0, 0, 0, 25e6, 8e6),
        ],
        T_BusDC=list(charger("hvdc2t").T_BusDC) + [(12, 2, 200e3, 1.00, 0.0)],
        T_GFL2s=list(charger("hvdc2t").T_GFL2s) + [(2, 5, 200e6, 200e3)],
        T_VSCs=list(charger("hvdc2t").T_VSCs) + [(3, 5, 12, "PQ")],
    )
    t.T_Lines_SI = list(t.T_Lines_SI) + [(3, 5, 6, 150.0, 0.03, 0.30, 4.0e-6)]
    t.T_IBs = list(t.T_IBs) + [(3, 6, 1000e6, 1.00, 0, 0.005, 0.05)]
    net = assembler(t)
    with pytest.raises(ValueError, match="ile continue") as exc:
        net.run_loadflow()
    message = str(exc.value)
    assert "12" in message, message
    assert "aucune reference" in message


def test_deux_references_dans_la_meme_ile():
    """Un NOEUD continu ne tient qu'une tension : deux references au meme
    noeud, c'est surdetermine. (Deux references a des noeuds distincts d'une
    meme ile sont admises : voir ``test_references_a_des_noeuds_distincts``.)"""
    net = Network(f0=50.0, Sb=100e6)
    net.add_bus(1, 230e3, btype=3, V=1.0, C=0.1)
    net.add_bus(2, 230e3, btype=1, V=1.0, C=0.1)
    net.add(BusDC(10, 200e3, C=0.3))
    net.bus_specs_dc[10].btype = DC_SLACK
    net.add(GFL2DC("CONV", 1, 10, Sn=100e6, Sb=100e6, wb=net.wb, Cdc_pu=0.3))
    net.vsc_modes[1] = "VdcQ"
    net.add(IB("GRID", 2, R=1e-3, X=1e-2, wb=net.wb))
    net.add(SourceDC("SDC", 10, R=0.05, X=0.5, wb=net.wb, E=1.0))
    with pytest.raises(ValueError, match="2 references de tension"):
        net.run_loadflow()


def test_ile_purement_passive_ne_reclame_pas_de_reference():
    """Un reseau continu SANS convertisseur n'a rien a prouver.

    C'est le cas des reseaux du lot 3 : une source derriere R, une charge
    resistive, un cable. Leur injection depend de leur tension, donc ils fixent
    le niveau d'eux-memes, et la garde ne doit pas se declencher.
    """
    net = Network(f0=50.0, Sb=1e6)
    net.add(BusDC(10, 640e3, C=0.05))
    net.add(BusDC(11, 640e3, C=0.05))
    net.add(SourceDC("SDC", 10, R=0.01, X=0.05, wb=net.wb, E=1.0))
    net.add(LineDC("CAB", 10, 11, R=0.01, X=0.3, wb=net.wb))
    net.add(ChargeDC_P("PDC", 11, P=0.4))
    net.solve_equilibrium()
    dx, dzu = ModeleNonLineaire(net).residu_equilibre()
    assert dx < 1e-8 and dzu < 1e-8


# ======================================================================
# 7) LES MODES : valides contre l'attribut de classe modes_dc
# ======================================================================
def test_modes_dc_declares():
    """Les modes retenus, ecrits en dur : ce test est la pour CASSER."""
    assert GFM3DC.modes_dc == ("PV", "VdcV", "Vtheta")
    assert GFM4DC.modes_dc == ("PV", "VdcV", "Vtheta")
    assert GFL2DC.modes_dc == ("VdcQ", "PQ", "droop", "droopI")
    for classe in (GFM3DC, GFM4DC, GFL2DC):
        assert set(classe.modes_dc) <= set(MODES_VSC)


@pytest.mark.parametrize(
    "mode, classe, admis",
    [
        ("VdcQ", "GFM3DC", "PV, VdcV, Vtheta"),
        ("PQ", "GFM3DC", "PV, VdcV, Vtheta"),
        ("PV", "GFL2DC", "VdcQ, PQ, droop, droopI"),
        ("VdcV", "GFL2DC", "VdcQ, PQ, droop, droopI"),
    ],
)
def test_mode_invalide_pour_la_classe(mode, classe, admis):
    """Un mode impossible est refuse A L'ASSEMBLAGE, en nommant la classe.

    Pas au Newton, trois etapes plus loin, sous la forme d'un systeme mal pose :
    un grid-forming tient la tension de son bus alternatif, il n'a pas de
    ``q_ref`` a suivre ; un grid-following ne la tient pas, il n'a pas de
    ``ve_ref``. La question se tranche en lisant la classe, donc elle se tranche
    ici.
    """
    bus_ac = 2 if classe == "GFM3DC" else 3
    t = _tables_hvdc2t_modifiees()
    t.T_VSCs = [
        (i, b, d, mode if b == bus_ac else m) for i, b, d, m in charger("hvdc2t").T_VSCs
    ]
    with pytest.raises(ValueError) as exc:
        assembler(t)
    message = str(exc.value)
    assert classe in message, message
    assert admis in message, message
    assert mode in message


def test_mode_inconnu_tout_court():
    """Un mode qui n'est dans aucune liste : le message donne les quatre."""
    t = _tables_hvdc2t_modifiees(T_VSCs=[(1, 2, 10, "PVdc"), (2, 3, 11, "VdcQ")])
    with pytest.raises(ValueError, match="PV, VdcV"):
        assembler(t)


def test_vsc_sur_un_bus_sans_convertisseur():
    """``T_VSCs`` ne pose pas de machine, elle en DESIGNE une."""
    t = _tables_hvdc2t_modifiees(
        T_VSCs=[(1, 2, 10, "PV"), (2, 3, 11, "VdcQ"), (3, 1, 10, "PQ")]
    )
    with pytest.raises(ValueError, match="aucun convertisseur declare"):
        assembler(t)


def test_vsc_sur_une_famille_sans_port_continu():
    """``GFM1`` a un bus continu IDEAL : il n'a rien a sortir dans le reseau."""
    t = _tables_hvdc2t_modifiees(
        T_GFM3s=[],
        T_GFM1s=[(1, 2, 200e6, 200e3)],
        T_VSCs=[(1, 2, 10, "PV"), (2, 3, 11, "VdcQ")],
    )
    with pytest.raises(ValueError, match="pas de variante a port"):
        assembler(t)


def test_bus_ac_est_la_cle_de_t_vscs():
    """Deux lignes de ``T_VSCs`` sur le meme bus alternatif : refuse."""
    t = _tables_hvdc2t_modifiees(
        T_VSCs=[(1, 2, 10, "PV"), (2, 3, 11, "VdcQ"), (3, 2, 11, "PV")]
    )
    with pytest.raises(ValueError, match="La cle de cette table est BusAC"):
        assembler(t)


def test_mode_vdc_exige_un_noeud_de_type_1():
    """Tenir la tension d'un noeud, c'est en faire le slack de son ile."""
    t = _tables_hvdc2t_modifiees(
        T_BusDC=[(10, 2, 200e3, 1.00, 0.0), (11, 2, 200e3, 1.00, 0.0)]
    )
    net = assembler(t)
    with pytest.raises(ValueError, match="declarer le bus 11 de type 1"):
        net.run_loadflow()


def test_bases_continues_incoherentes():
    """``Vdc_V`` du convertisseur et ``Un_V`` de son noeud doivent coincider.

    Reserve du lot 3, tranchee : ``kVdc`` ne pondere que le courant rendu au
    noeud, pas la tension que lisent les boucles internes. Plutot que de
    fabriquer silencieusement un modele faux, on refuse et on dit pourquoi.
    """
    t = _tables_hvdc2t_modifiees(T_GFM3s=[(1, 2, 200e6, 150e3)])
    with pytest.raises(ValueError, match="kVdc = 1"):
        assembler(t)


# ======================================================================
# 8) LES TABLES CONTINUES : capacites, contrat de colonnes, codecs
# ======================================================================
def test_capacite_du_noeud_continu_est_calculee_par_les_tables():
    """``C = somme(C/2 des cables) + somme(Cdc des VSC) + shunts``.

    Le lot 3 laissait l'appelant renseigner ``BusDC.C`` a la main et l'a signale
    comme une dette. C'est ``g2egrid.tables`` qui doit le faire, exactement comme
    ``Cbus`` du cote alternatif.
    """
    t = charger("hvdc2t")
    net = assembler(t)
    wb = net.wb
    Zb = 200e3**2 / t.Sb
    _id, i, j, km, _r, _l, c_km = t.T_LinesDC[0]
    demi = wb * (c_km * 1e-6 * km) * Zb / 2.0
    shunt = {int(r[1]): r[2] for r in t.T_ShuntsDC}
    # bus 10 : GFM3DC (Cdc_pu = 0.1 base convertisseur, Sn/Sb = 2)
    attendu_10 = demi + shunt[10] + 0.1 * (200e6 / t.Sb)
    # bus 11 : GFL2DC (Cdc_pu = 0.2 base convertisseur)
    attendu_11 = demi + shunt[11] + 0.2 * (200e6 / t.Sb)
    assert net.bus_specs_dc[10].C == pytest.approx(attendu_10, rel=1e-12)
    assert net.bus_specs_dc[11].C == pytest.approx(attendu_11, rel=1e-12)
    assert i == 10 and j == 11


def test_capacite_du_cable_va_aux_deux_noeuds():
    """Rallonger le cable augmente la capacite des DEUX noeuds, par moities."""
    t = charger("hvdc2t")
    avant = {b: s.C for b, s in assembler(t).bus_specs_dc.items()}
    t2 = t.copie(T_LinesDC=[(1, 10, 11, 200.0, 0.050, 0.35, 0.012)])
    apres = {b: s.C for b, s in assembler(t2).bus_specs_dc.items()}
    d10, d11 = apres[10] - avant[10], apres[11] - avant[11]
    assert d10 == pytest.approx(d11, rel=1e-12)
    assert d10 > 0


def test_gfl2dc_recoit_la_capacite_reelle_de_son_noeud():
    """``Cdc_pu`` du GFL2DC n'est plus un condensateur, c'est un REGLAGE.

    Le condensateur est dans le noeud ; mais le calage de ``Kpdc`` en a
    toujours besoin, et ce qu'il lui faut est la capacite REELLE du noeud,
    ramenee a la base du convertisseur.
    """
    t = charger("hvdc2t")
    net = assembler(t)
    gfl = next(e for e in net.elements if isinstance(e, GFL2DC))
    bus_dc = gfl.buses[1]
    attendu = net.bus_specs_dc[bus_dc].C * (t.Sb / gfl.Sn)
    assert gfl.pdict["Cdc_pu"] == pytest.approx(attendu, rel=1e-12)


def test_contrat_de_colonnes_complet():
    """``TABLES`` decrit exactement les champs-tables de ``Tables``."""
    champs = {n for n in Tables.__dataclass_fields__ if n.startswith(("Y_", "T_"))}
    assert set(TABLES) == champs, set(TABLES) ^ champs
    for nom in ("T_BusDC", "T_LinesDC", "T_ShuntsDC", "T_VSCs"):
        assert nom in TABLES
    assert TABLES["T_BusDC"] == ["BusID", "Type", "Un_V", "V_pu", "Pinj_W"]
    assert TABLES["T_VSCs"] == ["ID", "BusAC", "BusDC", "Mode"]


def test_largeur_des_lignes_conforme_au_contrat():
    """Aucune colonne cachee : chaque ligne a la largeur que TABLES annonce."""
    for nom in ("hvdc2t", "mtdc3t"):
        t = charger(nom)
        for table in ("T_BusDC", "T_LinesDC", "T_ShuntsDC", "T_VSCs"):
            for ligne in getattr(t, table):
                assert len(ligne) == len(TABLES[table]), (nom, table, ligne)


@pytest.mark.parametrize("nom", ["hvdc2t", "mtdc3t", "wscc9"])
def test_codecs_aller_retour(nom):
    """``tables_depuis_json(tables_vers_json(t))`` rend les memes tables.

    Les codecs ont quitte ``g2egrid_api.analyse`` pour vivre a cote de
    l'assembleur qui les lit : deux copies du contrat auraient derive des que
    les tables continues sont apparues.
    """
    t = charger(nom)
    t2 = tables_depuis_json(json.loads(json.dumps(tables_vers_json(t))))
    for table in TABLES:
        assert [tuple(r) for r in getattr(t, table)] == list(getattr(t2, table)), table
    assert (t2.f0, t2.Sb, t2.nom) == (t.f0, t.Sb, t.nom)
    assert assembler(t2).bus_specs_dc.keys() == assembler(t).bus_specs_dc.keys()


def test_codec_json_porte_les_tables_continues():
    """Le JSON contient les quatre tables continues, meme vides."""
    j = tables_vers_json(charger("wscc9"))
    for nom in ("T_BusDC", "T_LinesDC", "T_ShuntsDC", "T_VSCs"):
        assert j[nom] == []
    j2 = tables_vers_json(charger("hvdc2t"))
    assert j2["T_VSCs"] == [[1, 2, 10, "PV"], [2, 3, 11, "VdcQ"]]


# ======================================================================
# 9) LE SLACK CONTINU ALIMENTE PAR UNE SOURCE
# ======================================================================
def test_slack_continu_sur_source_deduit_sa_fem():
    """Un noeud continu de type 1 tenu par une source : ``E`` se DEDUIT.

    Le pendant exact de ce que fait ``IB.initialize`` du cote alternatif -- au
    slack, la tension est une donnee et la f.e.m. s'en deduit, pas l'inverse.
    Cela n'arrive QUE sur un noeud de type 1, qui se declare explicitement : un
    noeud de type 2 laisse la source avec le ``E`` qu'on lui a donne, et c'est
    sa tension qui s'ajuste (le comportement du lot 3).
    """
    net = Network(f0=50.0, Sb=1e6)
    net.add(BusDC(10, 640e3, C=0.05, V=1.02))
    net.bus_specs_dc[10].btype = DC_SLACK
    net.add(SourceDC("SDC", 10, R=0.05, X=0.5, wb=net.wb, E=1.0))
    net.add(ChargeDC_P("PDC", 10, P=0.3))
    net.solve_equilibrium()

    assert net.v_dc[10] == pytest.approx(1.02, abs=1e-12)
    src = next(e for e in net.elements if isinstance(e, SourceDC))
    # la source doit fournir exactement la charge : i = P/v
    assert src.x0[0] == pytest.approx(0.3 / 1.02, rel=1e-9)
    assert src.E == pytest.approx(1.02 + 0.05 * 0.3 / 1.02, rel=1e-9)
    dx, dzu = ModeleNonLineaire(net).residu_equilibre()
    assert dx < 1e-8 and dzu < 1e-8


def test_noeud_de_type_2_laisse_la_fem_de_la_source():
    """Sur un noeud de type 2, rien n'est deduit : c'est le chemin du lot 3."""
    net = Network(f0=50.0, Sb=1e6)
    net.add(BusDC(10, 640e3, C=0.05))
    net.add(SourceDC("SDC", 10, R=0.05, X=0.5, wb=net.wb, E=1.0))
    net.add(ChargeDC_P("PDC", 10, P=0.3))
    net.solve_equilibrium()
    src = next(e for e in net.elements if isinstance(e, SourceDC))
    assert src.E == 1.0
    # v(E - v)/R = P  =>  v^2 - E.v + R.P = 0
    v = net.v_dc[10]
    assert v * (1.0 - v) / 0.05 == pytest.approx(0.3, abs=1e-10)


def test_bus_dc_par_defaut_est_de_type_p_impose():
    """Le defaut est le comportement du lot 3, pas l'inverse."""
    net = Network(f0=50.0, Sb=1e6)
    net.add(BusDC(10, 640e3, C=0.05))
    assert net.bus_specs_dc[10].btype == DC_P
    assert net.bus_specs_dc[10].Pinj == 0.0


def test_pinj_du_noeud_continu_est_prise_en_compte():
    """``Pinj_W`` de ``T_BusDC`` entre bien dans le bilan du noeud."""
    t = charger("hvdc2t")
    base = assembler(t)
    base.run_loadflow()
    charge = t.copie(T_BusDC=[(10, 2, 200e3, 1.00, -10e6), (11, 1, 200e3, 1.00, 0.0)])
    net = assembler(charge)
    lf = net.run_loadflow()
    assert net.bus_specs_dc[10].Pinj == pytest.approx(-0.1, rel=1e-12)
    # le noeud 10 doit fournir 0,1 pu de plus : sa tension descend
    assert lf.Vdc[0] < base.lf.Vdc[0]


# ======================================================================
# 10) LA NON REGRESSION DU LOAD FLOW MIXTE : les memes deux filets
# ======================================================================
# Les reseaux continus du catalogue n'avaient aucun filet de non-regression
# du load flow. Les tests des sections 3 a 5 y verifient des proprietes (bilan
# global, loi des noeuds, point fixe du modele) sur hvdc2t et mtdc3t
# seulement, et test_cigre_b4.py compare cigre_b4 a PowerFactory a 1e-3 pu.
# Rien ne disait si le point de fonctionnement des six reseaux bougeait.
#
# LE PREMIER FILET compare V (alternatif) et Vdc (continu) a
# ``reference_loadflow_acdc.json``, a la meme ``TOL_V`` que l'alternatif. Le
# bruit emule (4 ulp sur chaque coefficient du jacobien et sur chaque p_dc
# de station) y vaut au plus 3,6e-13 pu sur V (cigre_mv_dc) et 9e-16 sur Vdc.
#
# ICI, LES ITERATIONS NE SE COMPTENT PAS A LA TOLERANCE D'ARRET. Le Newton
# mixte s'arrete a 1e-12 (``run_loadflow_acdc``). Or, sur cigre_b4,
# cigre_b4_mmc et cigre_mv_dc, le plancher d'arrondi du residu vaut de 3e-13
# a 1,3e-12, mesure en poursuivant le Newton apres convergence. Le nombre
# d'iterations pour descendre sous 1e-12 depend donc de la machine : sous
# bruit emule, cigre_b4 passe de 6 a 7 iterations et cigre_b4_mmc de 6 a 8,
# et une consigne deplacee de 1e-9 suffit aussi. On compte donc les
# iterations a ``TOL_ITERATIONS = 1e-10``, celle de l'alternatif : chaque
# mismatch y est a plus de cent fois le plancher, de part et d'autre du
# seuil. Les tensions, elles, se comparent au resultat publie (tolerance par
# defaut) : sous le meme bruit, elles ne bougent que de 4e-14.
#
# LE SECOND FILET recalcule les equations : celles de chaque bus alternatif
# selon le mode de sa station, celles de chaque noeud continu, et le bilan de
# chaque station, ``P_ac + P_dc + pertes = 0``. Les pertes y sont recalculees
# depuis les etats du convertisseur au point du load flow.
#
# Il etend ``test_bilan_*`` (bilan GLOBAL, apres ``solve_equilibrium``, sur
# deux reseaux) aux six reseaux, station par station, au point du load flow
# lui-meme. Il ne remplace ni ``test_residu_equilibre`` (accord avec le modele
# dynamique complet) ni la comparaison a PowerFactory.


@functools.cache
def _resolu_mixte(nom: str):
    """Un reseau continu du catalogue : son load flow publie (tolerance par
    defaut), et le nombre d'iterations qu'il faut pour descendre sous
    ``TOL_ITERATIONS``. Les filets partagent ce resultat sans le modifier."""
    net = assembler(charger(nom))
    iterations = net.run_loadflow(tol=TOL_ITERATIONS).iterations
    return net, net.run_loadflow(), iterations


def _reference_acdc(nom: str) -> dict:
    return json.loads(REFERENCE_ACDC.read_text(encoding="utf-8"))["reseaux"][nom]


def ecarts_mixtes(net: Network, lf, iterations: int, attendu: dict) -> list[str]:
    """Ce qui separe un load flow mixte de sa reference : une liste vide si rien.

    La reference est indexee par NUMERO de bus, et pas par position : un bus
    ajoute ou retire se signale par son nom, au lieu de decaler en silence
    toutes les comparaisons qui suivent.
    """
    bus, bus_dc = sorted(net.bus_specs), sorted(net.bus_specs_dc)
    if {str(b) for b in bus} != set(attendu["V"]):
        return [f"bus alternatifs {bus} au lieu de {list(attendu['V'])}"]
    if {str(b) for b in bus_dc} != set(attendu["Vdc"]):
        return [f"noeuds continus {bus_dc} au lieu de {list(attendu['Vdc'])}"]
    V_ref = np.array([complex(*attendu["V"][str(b)]) for b in bus])
    Vdc_ref = np.array([attendu["Vdc"][str(b)] for b in bus_dc])
    out = []
    dV = float(np.abs(lf.V - V_ref).max())
    if not dV <= TOL_V:
        out.append(f"max|dV| = {dV:.3e} pu > {TOL_V:.0e}")
    dVdc = float(np.abs(lf.Vdc - Vdc_ref).max())
    if not dVdc <= TOL_V:
        out.append(f"max|dVdc| = {dVdc:.3e} pu > {TOL_V:.0e}")
    if iterations != attendu["iterations_a_1e-10"]:
        out.append(
            f"{iterations} iterations sous {TOL_ITERATIONS:.0e} au lieu de "
            f"{attendu['iterations_a_1e-10']}"
        )
    return out


def etats_station(net: Network, el, vm, p_ac, q_ac, vdc) -> dict[str, float]:
    """Les etats d'un convertisseur au point ``(vm, p_ac, q_ac, vdc)``.

    Ce sont ses formules d'initialisation, sous le contexte que pose le
    couplage du load flow : tension alternative sur l'axe reel et
    ``theta_g = 0``. Le courant continu ne depend pas de la phase.
    """
    ctx = InitContext(
        v={el.ports[0]: complex(vm, 0.0), el.ports_dc[0]: float(vdc)},
        s={el.ports[0]: complex(p_ac, q_ac)},
        w=net.w_sys,
        f0=net.f0,
        theta_g=0.0,
    )
    return el.initialize(ctx)


def residus_mixtes(net: Network, lf) -> dict[str, float]:
    """Les equations du load flow MIXTE, recalculees ICI, et le bilan des
    stations.

    De ``net._probleme_acdc()``, on ne prend que les DONNEES des stations :
    mode, consignes, ``p_dc`` (le modele de :class:`~g2egrid.loadflow.VSCSpec`).
    Les equations sont reecrites ici :

    - aux bus alternatifs, ``P`` et ``Q`` en ecriture polaire depuis le Y-bus,
      contre la loi du mode de la station (voir ``MODES_VSC``). Une station en
      ``Vdc...`` y prend ``lf.P_vsc``, son inconnue de puissance : l'equation
      verifie que la valeur rendue est bien celle que le bus injecte ;
    - aux noeuds continus de type 2, et a ceux que tient une station, le
      bilan ``Pinj + P(elements) + P_dc(stations) - V_j.somme G_jk.(V_j -
      V_k)``. Les conductances sont reprises des cables, et le flux est ecrit
      en DIFFERENCES de tension, la ou le load flow forme ``V.(G.V)`` ;
    - a chaque station, ``P_ac + P_dc + pertes = 0``. ``P_ac`` est ce que la
      station pose sur son bus (le load flow, moins les autres elements),
      ``P_dc`` sort de son ``p_dc`` a ce point, et les pertes de ses etats
      (:func:`pertes_convertisseur`).
    """
    pb = net._probleme_acdc()
    P, Q = injections_polaires(net.build_ybus(), lf.V)
    Vm, th = np.abs(lf.V), np.angle(lf.V)
    vdc = lf.Vdc
    elements = net._elements_vsc()
    assert [el.name for el in elements] == [c.nom for c in pb.vscs]
    par_bus = {c.iac: c for c in pb.vscs}
    libres = [c.nom for c in pb.vscs if c.tient_vdc]  # l'ordre de lf.P_vsc
    noms = ("P", "Q", "V impose", "theta impose", "Vdc impose", "bilan DC")
    r = dict.fromkeys(noms + ("Pdc rendu", "station"), 0.0)

    # --- les bus alternatifs ----------------------------------------------
    for b, i in net._bus_index().items():
        s = net.bus_specs[b]
        if s.btype == SLACK:  # y compris l'ile qu'une station Vtheta forme
            r["V impose"] = max(r["V impose"], abs(Vm[i] - s.V))
            ecart = abs(np.angle(np.exp(1j * (th[i] - s.theta))))
            r["theta impose"] = max(r["theta impose"], ecart)
            continue
        p_sp, q_sp = s.Pg - s.Pl, s.Qg - s.Ql
        c = par_bus.get(i)
        if c is not None:
            v = vdc[c.idc]
            if c.tient_vdc:
                p_sp += lf.P_vsc[libres.index(c.nom)]
            elif c.mode == "droop":  # p = p_ref + s.(v_dc - v_ref)
                p_sp += c.pente * (v - c.vdc_ref)
            elif c.mode == "droopI":  # p = v_dc.(p_ref/v_ref + s.(v_dc - v_ref))
                p_sp += v * (c.p_ac / c.vdc_ref + c.pente * (v - c.vdc_ref)) - c.p_ac
            if s.btype == PQ:  # statisme Q-V : q = q_ref + kqv.(v_ref - |V|)
                q_sp += c.kqv * (c.v_ref_ac - Vm[i])
        r["P"] = max(r["P"], abs(P[i] - p_sp))
        if s.btype == PQ:
            r["Q"] = max(r["Q"], abs(Q[i] - q_sp))
        else:
            r["V impose"] = max(r["V impose"], abs(Vm[i] - s.V))

    # --- les noeuds continus ----------------------------------------------
    pos = {b: j for j, b in enumerate(sorted(net.bus_specs_dc))}
    n = len(pos)
    g = np.zeros((n, n))  # conductance de cable entre deux noeuds
    p_propre = np.zeros(n)
    for el in net.elements:
        if el.ports or not el.ports_dc:
            continue
        if isinstance(el, LineDC):
            f, t = pos[el.buses[0]], pos[el.buses[1]]
            g[f, t] += 1.0 / el.pdict["R"]
            g[t, f] += 1.0 / el.pdict["R"]
        else:
            j = pos[el.buses[0]]
            p_propre[j] += net._injection_dc_propre(el)(vdc[j])
    p_dc = np.zeros(n)
    pertes_min = np.inf
    for el, c in zip(elements, pb.vscs):
        i, j = c.iac, c.idc
        p_ac, q_ac = P[i] - c.p_autres, Q[i] - c.q_autres
        p = c.p_dc(Vm[i], p_ac, q_ac, vdc[j])
        p_dc[j] += p
        etats = etats_station(net, el, Vm[i], p_ac, q_ac, vdc[j])
        pertes = pertes_convertisseur(el, etats, vdc[j])
        pertes_min = min(pertes_min, pertes)
        r["station"] = max(r["station"], abs(p_ac + p + pertes))
    flux = vdc * ((g * (vdc[:, None] - vdc[None, :])).sum(axis=1))
    pinj = np.array([net.bus_specs_dc[b].Pinj for b in sorted(pos)])
    bilan = pinj + p_propre + p_dc - flux
    tenus = {c.idc for c in pb.vscs if c.tient_vdc}
    for b, j in pos.items():
        s = net.bus_specs_dc[b]
        if s.btype == DC_P or j in tenus:
            r["bilan DC"] = max(r["bilan DC"], abs(bilan[j]))
        if s.btype == DC_SLACK:
            r["Vdc impose"] = max(r["Vdc impose"], abs(vdc[j] - s.V))
    r["Pdc rendu"] = float(np.abs(lf.Pdc - p_dc).max())
    r["pertes station min"] = float(pertes_min)
    r["pertes cables"] = float(0.5 * (g * (vdc[:, None] - vdc[None, :]) ** 2).sum())
    return r


@pytest.mark.parametrize("nom", CONTINUS)
def test_load_flow_mixte_contre_reference(nom):
    """LE PREMIER FILET, cote continu : V et Vdc a ``TOL_V``, et le meme
    nombre d'iterations sous 1e-10.

    La reference est ``reference_loadflow_acdc.json``, produite par
    ``developpement/regenerer_reference_loadflow.py continu`` depuis un arbre
    commite. Sa provenance (commit, versions, plateforme) est dans le fichier.
    """
    net, lf, iterations = _resolu_mixte(nom)
    ecarts = ecarts_mixtes(net, lf, iterations, _reference_acdc(nom))
    assert not ecarts, f"{nom} : le load flow mixte a bouge : " + " ; ".join(ecarts)


@pytest.mark.parametrize("nom", CONTINUS)
def test_load_flow_mixte_resout_ses_equations(nom):
    """LE SECOND FILET, cote continu : chaque equation sous 1e-10 pu.

    1e-10 est le seuil du second filet alternatif. Le Newton mixte s'arrete
    a 1e-12, et le plancher d'arrondi de ces reseaux vaut au plus 1,3e-12 :
    la marge est de deux ordres de grandeur. Le bilan des stations n'est pas
    une equation du Newton, et il tombe au niveau de l'arrondi des puissances
    (jusqu'a 24 pu sur cigre_b4).

    Les grandeurs imposees (``|V|`` des bus PV et slack, ``theta`` du slack,
    ``Vdc`` des noeuds de type 1) sont tenues a 1e-14 pres : ce ne sont pas
    des inconnues. Les pertes de chaque station, recalculees depuis ses
    etats, sont strictement positives.
    """
    net, lf, _ = _resolu_mixte(nom)
    r = residus_mixtes(net, lf)
    for cle in ("P", "Q", "bilan DC", "Pdc rendu", "station"):
        assert r[cle] < 1e-10, (nom, cle, r)
    for cle in ("V impose", "theta impose", "Vdc impose"):
        assert r[cle] < 1e-14, (nom, cle, r)
    assert r["pertes station min"] > 0.0, (nom, r)


def _station_a_puissance_programmee(net: Network) -> int:
    """Le bus alternatif de la premiere station dont la puissance est une
    CONSIGNE : ni ``Vdc...`` (sa puissance est l'inconnue) ni ``Vtheta``
    (elle est le slack de son ile)."""
    return min(
        el.buses[0]
        for el in net._elements_vsc()
        if net.bus_specs[el.buses[0]].btype != SLACK
        and not net.vsc_modes.get(el.buses[0], "").startswith("Vdc")
    )


@pytest.mark.parametrize("mutation", ["Vdc", "P"])
@pytest.mark.parametrize("nom", ["hvdc2t", "cigre_b4"])
def test_le_filet_continu_voit_une_consigne_deplacee(nom, mutation):
    """TEMOIN DE MUTATION, cote continu : la consigne de tension d'un noeud de
    type 1, ou celle de puissance d'une station, deplacee de 1e-6.

    La consigne ``Vdc`` se retrouve telle quelle dans ``Vdc``. Celle de
    puissance deplace V d'au moins 2,9e-8 pu sur les deux reseaux (mesure).
    ``hvdc2t`` est le cas d'ecole. ``cigre_b4`` exerce les modes ``droop`` et
    ``Vtheta``, et le reseau continu maille.

    Comme cote alternatif, le second filet RESTE vert sur le reseau mute.
    """
    net = assembler(charger(nom))
    if mutation == "Vdc":
        b = min(b for b, s in net.bus_specs_dc.items() if s.btype == DC_SLACK)
        net.bus_specs_dc[b].V += 1e-6
        attendu = "max|dVdc|"
    else:
        b = _station_a_puissance_programmee(net)
        net.bus_specs[b].Pg += 1e-6
        attendu = "max|dV|"
    iterations = net.run_loadflow(tol=TOL_ITERATIONS).iterations
    lf = net.run_loadflow()
    ecarts = ecarts_mixtes(net, lf, iterations, _reference_acdc(nom))
    assert any(e.startswith(attendu) for e in ecarts), (nom, b, mutation, ecarts)
    r = residus_mixtes(net, lf)
    assert max(r[c] for c in ("P", "Q", "bilan DC", "station")) < 1e-10, r


def test_le_filet_continu_voit_un_jacobien_degrade(monkeypatch):
    """TEMOIN : un jacobien mixte faux de 1 % ralentit le Newton, sans
    deplacer la solution. Le filet le voit par les iterations SEULEMENT.

    Mesure sur ``hvdc2t`` : 7 iterations sous 1e-10 au lieu de 4, et V a
    1,2e-15 pu de la reference.
    """
    net = assembler(charger("hvdc2t"))
    attendu = _reference_acdc("hvdc2t")
    jacobien = ProblemeACDC.jacobien
    monkeypatch.setattr(
        ProblemeACDC, "jacobien", lambda self, x: 1.01 * jacobien(self, x)
    )
    iterations = net.run_loadflow(tol=TOL_ITERATIONS).iterations
    ecarts = ecarts_mixtes(net, net.run_loadflow(), iterations, attendu)
    assert len(ecarts) == 1 and "iterations" in ecarts[0], ecarts
    assert iterations > attendu["iterations_a_1e-10"]


# ======================================================================
# 11) UNE LIMITE CONNUE : la tolerance d'arret du Newton mixte
# ======================================================================
@pytest.mark.xfail(
    strict=True,
    raises=RuntimeError,
    reason=(
        "FR-001 : tol = 1e-12 (run_loadflow_acdc) est sous le plancher "
        "d'arrondi eps.|G_dc| d'un cable continu court ; le Newton atteint le "
        "point fixe et leve 'non converge' apres 40 iterations"
    ),
)
def test_cable_continu_court_converge():
    """``hvdc2t`` avec un cable de 30 m : le load flow doit converger.

    Il ne converge pas. ``run_loadflow_acdc`` s'arrete quand ``max|F| < 1e-12``
    pu. Mais une tension continue ne se represente qu'a un ulp pres (2,2e-16
    autour de 1 pu), et une erreur d'un ulp sur ``V_j`` deplace le bilan du
    noeud de ``G_jk.ulp``. A 30 m (``G = 2,7e5`` pu), le residu ne descend pas
    sous 2,3e-11 : le Newton est au point fixe des la cinquieme iteration, puis
    tourne jusqu'a la quarantieme et leve ``RuntimeError: load flow AC/DC non
    converge``. Mesure : 300 m (``G = 2,7e4``) echoue aussi, a 2,1e-12 ; 100 m
    converge par chance (2,1e-13) ; 1 km passe a 8,6e-13. Une station
    dos-a-dos, modelisee par un cable de quelques dizaines de metres, ne
    passe donc pas.

    Le catalogue frole deja ce plancher : 3e-13 a 1,3e-12 sur cigre_b4 et
    cigre_b4_mmc (|G_dc| = 4,3e3 pu), 2e-13 a 6e-13 sur cigre_mv_dc (|Y| =
    1,9e3 pu, cote alternatif). Il converge en ce sens que l'arrondi y tombe
    une fois sous 1e-12.

    Ce que le test exigera, une fois le critere d'arret corrige : la
    convergence, les equations sous 1e-10, et une chute de tension continue
    proportionnelle a la longueur. Le cable de 10 km sert de comparaison :
    les pertes de son cable deplacent le courant de 6e-4 en relatif, d'ou
    ``rel = 2e-3``.
    """
    t = charger("hvdc2t")
    ident, i, j, _km, r_km, l_km, c_km = t.T_LinesDC[0]
    long_ = assembler(t.copie(T_LinesDC=[(ident, i, j, 10.0, r_km, l_km, c_km)]))
    long_.run_loadflow()
    court = assembler(t.copie(T_LinesDC=[(ident, i, j, 0.03, r_km, l_km, c_km)]))
    lf = court.run_loadflow()  # RuntimeError aujourd'hui

    def chute_par_km(net, km):
        return (net.lf.Vdc[0] - net.lf.Vdc[1]) / km

    assert chute_par_km(court, 0.03) == pytest.approx(
        chute_par_km(long_, 10.0), rel=2e-3
    )
    r = residus_mixtes(court, lf)
    assert max(r[c] for c in ("P", "Q", "bilan DC", "station")) < 1e-10, r


# ======================================================================
# 12) NOEUD DE TYPE 1 SANS TENANT, JACOBIEN SINGULIER, PLANCHER DU RESIDU
# ======================================================================
def test_noeud_type_un_sans_tenant_refuse():
    """Un noeud de type 1 que rien ne tient serait un slack IDEAL du load flow
    (tension imposee, bilan jamais ecrit) alors que le modele dynamique n'y met
    aucune source : le point calcule n'y serait pas un equilibre. Sur mtdc3t,
    le noeud 12 passe de type 2 a type 1 alors que sa station est en PQ : l'ile
    a ailleurs une reference (GFL en VdcQ au noeud 11), et c'est ce cas, plus
    sournois que l'ile sans reference, que la garde doit attraper."""
    t = charger("mtdc3t")
    bus = [(*r[:1], 1, *r[2:]) if int(r[0]) == 12 else r for r in t.T_BusDC]
    net = assembler(t.copie(T_BusDC=bus))
    with pytest.raises(ValueError, match=r"noeud continu 12 : declare de type 1"):
        net.run_loadflow()
    # le reseau d'origine, lui, passe
    assembler(t).run_loadflow()


def test_jacobien_singulier_rend_un_resultat_non_converge(hvdc2t):
    """Une ``LinAlgError`` du Newton ne sort plus en exception brute : le
    resultat est non converge et porte un message qui dit pourquoi."""
    pb = hvdc2t._probleme_acdc()
    pb.jacobien = lambda x: np.zeros((pb.n_x, pb.n_x))
    res = pb.resoudre()
    assert not res.converged
    assert "jacobien singulier" in res.message
    assert res.iterations == 1


def test_jacobien_singulier_message_du_reseau(monkeypatch):
    net = assembler(charger("hvdc2t"))
    vrai = ProblemeACDC.jacobien

    def singulier(self, x):
        return np.zeros_like(vrai(self, x))

    monkeypatch.setattr(ProblemeACDC, "jacobien", singulier)
    with pytest.raises(RuntimeError, match="jacobien singulier"):
        net.run_loadflow()


@pytest.mark.parametrize(
    "nom, plafond",
    [
        ("hvdc2t", 1e-10),
        ("mtdc3t", 1e-9),
        ("cigre_b4_dcs1", 1e-7),
        ("cigre_b4_mmc", 5e-7),
        ("cigre_b4", 1e-6),
        ("cigre_mv_dc", 1e-4),
    ],
)
def test_residu_equilibre_du_catalogue_continu(nom, plafond):
    """Le residu ``|dx/dt|`` au point de fonctionnement mixte, tel que le load
    flow le laisse. Il n'est pas nul pour une raison arithmetique et pas
    physique : le Newton s'arrete a un mismatch de puissance de l'ordre de
    4e-13 (plancher d'arrondi, cf. ``run_loadflow_acdc``), et la loi des noeuds
    s'ecrit ``dv/dt = (wb/C).isum``. Sur ``cigre_mv_dc`` la capacite d'un bus
    alternatif n'est que de 3,8e-6 pu : le gain wb/C vaut 8e7 et un courant
    residuel de 3e-13 ressort a 2e-5. Une iteration de polissage n'y change
    rien (le mismatch plafonne a 4e-13 au mieux) ; on borne donc le residu et,
    plus bas, le courant residuel qui, lui, est a l'arrondi."""
    net = assembler(charger(nom))
    net.solve_equilibrium()
    dx, dzu = ModeleNonLineaire(net).residu_equilibre()
    assert dx < plafond, f"{nom} : |dx/dt| = {dx:.3e}"
    assert dzu < 1e-9, f"{nom} : dzu = {dzu:.3e}"


def test_residu_cigre_mv_dc_est_un_courant_a_l_arrondi():
    """Le residu de ``cigre_mv_dc`` rapporte au gain wb/C du noeud le plus
    leger est un courant de l'ordre de l'arrondi (< 1e-12 pu)."""
    net = assembler(charger("cigre_mv_dc"))
    net.solve_equilibrium()
    dx, _ = ModeleNonLineaire(net).residu_equilibre()
    c_min = min(s.C for s in net.bus_specs.values() if s.C > 0)
    assert dx / (net.wb / c_min) < 1e-12
