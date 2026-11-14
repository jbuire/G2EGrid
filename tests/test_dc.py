"""Tests des PORTS CONTINUS : liaisons VSC-HVDC / MTDC.

Un port continu n'apporte qu'UN signal de chaque cote (``v_<p>`` en entree,
``i_<p>`` en sortie) la ou un port alternatif en apporte deux. Tout ce qui
suit verifie que cette asymetrie est tenue partout : dans les indices de
signaux, dans le cablage, dans la loi des noeuds, et dans les sorties
exposees.

Le test central est :func:`test_equivalence_gfm3_decoupe` : ``GFM3`` porte son
condensateur et sa boucle de courant DC en interne ; ``GFM3DC + BusDC +
SourceDC`` les sort dans le reseau. Les deux doivent avoir le MEME spectre,
sinon le decoupage a change la physique.
"""

import ast
import pathlib
import subprocess

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

from g2egrid import Network
from g2egrid.elements import (
    GFL2DC,
    GFM1,
    GFM3,
    GFM3DC,
    GFM4,
    GFM4DC,
    IB,
    Branch,
    Bus,
    BusDC,
    ChargeDC_P,
    ChargeDC_R,
    Line,
    LineDC,
    SourceDC,
)
from g2egrid.elements.gfm import GFM_DC_DEFAULT
from g2egrid.temporel import ModeleNonLineaire, _plan_mesures, mesures

RACINE = pathlib.Path(__file__).resolve().parents[2]


# ----------------------------------------------------------------------
# fixtures
# ----------------------------------------------------------------------
def reseau_rc(R=2.0, C=0.05):
    """Le plus petit reseau continu qui existe : une capacite et une resistance."""
    net = Network(f0=50.0, Sb=1e6)
    net.add(BusDC(10, 640e3, C=C))
    net.add(ChargeDC_R("RL", 10, R=R))
    return net


def reseau_cable(C1=0.05, C2=0.08, R=0.01, X=0.3):
    """Deux noeuds continus relies par un cable : un mode de jauge et une paire."""
    net = Network(f0=50.0, Sb=1e6)
    net.add(BusDC(10, 640e3, C=C1))
    net.add(BusDC(11, 640e3, C=C2))
    net.add(LineDC("CAB", 10, 11, R=R, X=X, wb=net.wb))
    return net


def reseau_gfm3(Pg=0.5, C=1e-3, classe=GFM3):
    """GFM monolithique (bus 1) -- ligne -- reseau amont (bus 2, slack)."""
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=3, V=1.0, Pg=Pg, C=C)
    net.add_bus(2, 20e3, btype=1, V=1.0, C=C)
    net.add(classe("CONV", 1, Sn=1e6, Sb=1e6, wb=net.wb))
    net.add(Line("L1", 1, 2, R=0.02, X=0.20, wb=net.wb))
    net.add(IB("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
    return net


def equivalent_dc(ref: Network, classe=GFM3DC, Pg=0.5, C=1e-3):
    """Le meme reseau, mais le cote continu du GFM sorti dans le reseau.

    ``ref`` est le reseau :func:`reseau_gfm3` DEJA resolu : c'est lui qui
    fournit le ``Kpdc`` auto-regle et les consignes ``p_ref``/``vdc_ref`` du
    point de fonctionnement, dont se deduisent les trois parametres de la
    source continue (cf. la docstring de ``SourceDC``).
    """
    gfm = next(e for e in ref.elements if "Kpdc" in e.params)
    Kpdc = gfm.pdict["Kpdc"]
    Tdc = gfm.pdict["Tdc"]
    Cdc = gfm.pdict["Cdc_pu"]
    p_ref = gfm.u0[gfm.svc_inputs.index("p_ref")]
    vdc_ref = gfm.u0[gfm.svc_inputs.index("vdc_ref")]

    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=3, V=1.0, Pg=Pg, C=C)
    net.add_bus(2, 20e3, btype=1, V=1.0, C=C)
    net.add(BusDC(10, 640e3, C=Cdc, V=vdc_ref))
    net.add(classe("CONV", 1, 10, Sn=1e6, Sb=1e6, wb=net.wb))
    net.add(Line("L1", 1, 2, R=0.02, X=0.20, wb=net.wb))
    net.add(IB("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
    net.add(
        SourceDC(
            "SDC",
            10,
            R=1.0 / Kpdc,
            X=net.wb * Tdc / Kpdc,
            wb=net.wb,
            E=vdc_ref + p_ref / (Kpdc * vdc_ref),
        )
    )
    return net


def apparier(la, lb):
    """Ecart maximal apparie entre deux spectres (appariement hongrois)."""
    C = np.abs(la[:, None] - lb[None, :]) / np.maximum(1.0, np.abs(la[:, None]))
    r, c = linear_sum_assignment(C)
    return float(C[r, c].max())


def courants_dc(net: Network) -> dict[int, float]:
    """Somme des courants de port continus, par noeud continu, a l'equilibre."""
    somme = {b: 0.0 for b in net.bus_specs_dc}
    for el in net.elements:
        if not el.ports_dc:
            continue
        y = el.compiled().h(el.x0, el.z0, el.u0, el.p)
        for pk in range(len(el.ports_dc)):
            b = el.buses[len(el.ports) + pk]
            somme[b] += float(y[el.port_dc_output_idx(pk)])
    return somme


# ----------------------------------------------------------------------
# 1) le circuit RC continu
# ----------------------------------------------------------------------
def test_rc_continu_pole_exact():
    """Un noeud capacitif et une charge resistive : ``lambda = -wb/(R.C)``.

    C'est ``-1/(R.C)`` en temps pu ; le modele integre en secondes, d'ou le
    ``wb`` de ``dv/dt = (wb/C).isum``. Le pole doit sortir EXACTEMENT, a la
    precision de ``eigvals`` : il n'y a pas d'approximation dans la chaine,
    seulement un element, un cablage et une loi des noeuds.
    """
    R, C = 2.0, 0.05
    net = reseau_rc(R, C)
    ss = net.build_state_space()
    assert len(ss.states) == 1
    lam = ss.eig()
    attendu = -net.wb / (R * C)
    assert abs(lam[0] - attendu) / abs(attendu) < 1e-12


# ----------------------------------------------------------------------
# 2) deux noeuds et un cable
# ----------------------------------------------------------------------
def test_cable_continu_spectre_analytique():
    """Spectre de deux capacites reliees par un R-L continu.

    Sur l'ecart ``u = v1 - v2`` le systeme se reduit a un second ordre ::

        lambda^2 + (wb.R/X).lambda + wb^2.(1/C1 + 1/C2)/X = 0

    et la translation en bloc des deux tensions laisse tout inchange : c'est
    le troisieme mode, exactement nul.
    """
    C1, C2, R, X = 0.05, 0.08, 0.01, 0.3
    net = reseau_cable(C1, C2, R, X)
    ss = net.build_state_space()
    assert len(ss.states) == 3

    wb = net.wb
    racines = np.roots([1.0, wb * R / X, wb**2 * (1 / C1 + 1 / C2) / X])
    attendu = np.concatenate((racines, [0.0]))
    assert apparier(ss.eig(), attendu) < 1e-12


# ----------------------------------------------------------------------
# 3) LE test : le decoupage ne change pas la physique
# ----------------------------------------------------------------------
def test_equivalence_gfm3_decoupe():
    """``GFM3DC + BusDC + SourceDC`` a EXACTEMENT le spectre de ``GFM3``.

    Les deux etats continus de ``GFM3`` sont sortis dans le reseau, pas
    supprimes : le condensateur devient la capacite du noeud continu, la
    boucle de courant DC devient la source. Le nombre d'etats est donc le
    meme des deux cotes, et les valeurs propres aussi -- sinon le decoupage
    a introduit ou perdu de la dynamique.
    """
    ref = reseau_gfm3()
    ref.solve_equilibrium()
    ssa = ref.build_state_space()

    net = equivalent_dc(ref)
    ssb = net.build_state_space()

    assert len(ssa.states) == len(ssb.states)
    assert apparier(ssa.eig(), ssb.eig()) < 1e-8


def test_equivalence_tension_continue():
    """Le noeud continu se cale sur le ``vdc`` que ``GFM3`` trouvait en interne.

    Corollaire du precedent, mais il isole la cause la plus probable d'un
    ecart de spectre : deux points de LINEARISATION differents.
    """
    ref = reseau_gfm3()
    ref.solve_equilibrium()
    gfm = next(e for e in ref.elements if isinstance(e, GFM3))
    vdc_ref = gfm.x0[gfm.states.index("vdc")]
    idc_ref = gfm.x0[gfm.states.index("idc")]

    net = equivalent_dc(ref)
    net.solve_equilibrium()
    bus = net.buses_dc[10]
    src = next(e for e in net.elements if isinstance(e, SourceDC))
    assert abs(bus.x0[0] - vdc_ref) < 1e-10
    assert abs(src.x0[0] - idc_ref) < 1e-10


def test_gfm4dc_meme_nombre_d_etats():
    """``GFM4DC`` perd les deux memes etats que ``GFM3DC``, pas d'autres."""
    assert len(GFM3DC.states) == 13
    assert len(GFM4DC.states) == 14


def test_equivalence_gfm4_decoupe():
    """Meme equivalence pour le statisme du 2e ordre : ``GFM4`` = ``GFM4DC`` decoupe.

    Le statisme n'a rien a voir avec le cote continu ; si l'equivalence ne
    tenait que pour ``GFM3``, c'est que le decoupage aurait attrape au passage
    quelque chose qui ne lui appartient pas.
    """
    ref = reseau_gfm3(classe=GFM4)
    ref.solve_equilibrium()
    net = equivalent_dc(ref, classe=GFM4DC)
    ssa, ssb = ref.build_state_space(), net.build_state_space()
    assert len(ssa.states) == len(ssb.states)
    assert apparier(ssa.eig(), ssb.eig()) < 1e-8


def test_gfl2dc_s_assemble():
    """``GFL2DC`` : 12 etats, le noeud continu en porte un, la source un autre.

    PAS d'equivalence spectrale avec ``GFL2`` ici, et ce n'est pas un defaut :
    ``GFL2`` alimente son condensateur par un courant de consigne filtre
    (``idc_ref`` constant, donc aucune reaction a ``vdc``) la ou
    :class:`SourceDC` est une f.e.m. derriere R-L. Les deux cotes continus ne
    sont pas le meme circuit -- seul le grid-forming avait, lui, une boucle
    exactement transposable.
    """
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=3, V=1.0, Pg=0.5, C=1e-3)
    net.add_bus(2, 20e3, btype=1, V=1.0, C=1e-3)
    net.add(BusDC(10, 640e3, C=0.2))
    net.add(GFL2DC("CONV", 1, 10, Sn=1e6, Sb=1e6, wb=net.wb, Cdc_pu=0.2))
    net.add(Line("L1", 1, 2, R=0.02, X=0.20, wb=net.wb))
    net.add(IB("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
    net.add(SourceDC("SDC", 10, R=1e-3, X=1e-2, wb=net.wb, E=1.0))

    net.solve_equilibrium()
    assert len(GFL2DC.states) == 12
    for b, i in courants_dc(net).items():
        assert abs(i) < 1e-9, f"bus {b}: somme des courants {i:.3e}"
    ss = net.build_state_space()
    assert len(ss.states) == 12 + 1 + 1 + 2 + 2 + 2 + 3
    assert not [o for o in ss.outputs if o.endswith(".i_dc")]
    assert "CONV.vdc_out" in ss.outputs


# ----------------------------------------------------------------------
# 4) la loi des noeuds continue
# ----------------------------------------------------------------------
def test_loi_des_noeuds_continue():
    """La somme des courants a chaque noeud continu est nulle a l'equilibre."""
    ref = reseau_gfm3()
    ref.solve_equilibrium()
    net = equivalent_dc(ref)
    net.solve_equilibrium()
    for b, i in courants_dc(net).items():
        assert abs(i) < 1e-9, f"bus {b}: somme des courants {i:.3e}"


def test_loi_des_noeuds_continue_mtdc():
    """Meme chose sur une liaison a deux noeuds, cable et charge a P constante."""
    net = Network(f0=50.0, Sb=1e6)
    net.add(BusDC(10, 640e3, C=0.05))
    net.add(BusDC(11, 640e3, C=0.05))
    net.add(SourceDC("SDC", 10, R=0.01, X=0.05, wb=net.wb, E=1.0))
    net.add(LineDC("CAB", 10, 11, R=0.01, X=0.3, wb=net.wb))
    net.add(ChargeDC_P("PDC", 11, P=0.4))
    net.solve_equilibrium()
    for b, i in courants_dc(net).items():
        assert abs(i) < 1e-9, f"bus {b}: somme des courants {i:.3e}"
    ss = net.build_state_space()
    assert len(ss.states) == 5  # deux tensions, trois courants
    assert ss.eig().real.max() < 1e-6


def test_vsc_alimente_par_un_cable_continu():
    """Deux noeuds continus, un cable, un VSC a un bout et une source a l'autre.

    La configuration que le lot vise vraiment : le convertisseur ne partage
    plus son condensateur avec sa source, un cable les separe. C'est aussi le
    seul cas ou le Newton continu a plus d'une inconnue et un vrai couplage
    entre elles.
    """
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=3, V=1.0, Pg=0.5, C=1e-3)
    net.add_bus(2, 20e3, btype=1, V=1.0, C=1e-3)
    net.add(BusDC(10, 640e3, C=0.08))
    net.add(BusDC(11, 640e3, C=0.02))
    net.add(GFM3DC("CONV", 1, 10, Sn=1e6, Sb=1e6, wb=net.wb))
    net.add(Line("L1", 1, 2, R=0.02, X=0.20, wb=net.wb))
    net.add(IB("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
    net.add(LineDC("CAB", 11, 10, R=0.01, X=0.1, wb=net.wb))
    net.add(SourceDC("SDC", 11, R=0.05, X=0.5, wb=net.wb, E=1.05))

    net.solve_equilibrium()
    for b, i in courants_dc(net).items():
        assert abs(i) < 1e-9, f"bus {b}: somme des courants {i:.3e}"
    # le cable transporte : la tension chute de la source vers le convertisseur
    assert net.v_dc[11] > net.v_dc[10]
    ss = net.build_state_space()
    assert len(ss.states) == 13 + 2 + 2 + 1 + 1 + 2 + 3 + 1 + 1
    assert not [o for o in ss.outputs if o.endswith(".i_dc")]


def test_source_continue_trop_faible():
    """Une source continue qui ne peut pas fournir la charge : message explicite.

    Vu du continu, un VSC est une charge a puissance constante ; une f.e.m.
    derriere ``R`` plafonne a ``E^2/(4R)``. Au-dela il n'existe aucune tension
    d'equilibre, et l'echec doit le dire plutot que de rendre un ``sol.message``
    de MINPACK tout seul.
    """
    net = Network(f0=50.0, Sb=1e6)
    net.add(BusDC(10, 640e3, C=0.05))
    net.add(SourceDC("SDC", 10, R=1.0, X=0.5, wb=net.wb, E=1.0))
    net.add(ChargeDC_P("PDC", 10, P=0.9))  # au-dela de E^2/(4R) = 0.25
    with pytest.raises(RuntimeError, match="puissance constante"):
        net.solve_equilibrium()


# ----------------------------------------------------------------------
# 5) le point de fonctionnement est un point fixe du modele non lineaire
# ----------------------------------------------------------------------
def test_residu_equilibre_reseau_continu():
    """``residu_equilibre`` nul sur un reseau PUREMENT continu.

    C'est la verification croisee entre les equilibres resolus PAR ELEMENT et
    le reseau entier : si le Newton sur les tensions continues s'etait arrete
    trop tot, c'est ici que cela se verrait, et ce reseau n'a aucun element
    alternatif pour diluer la mesure.
    """
    net = Network(f0=50.0, Sb=1e6)
    net.add(BusDC(10, 640e3, C=0.05))
    net.add(BusDC(11, 640e3, C=0.05))
    net.add(SourceDC("SDC", 10, R=0.01, X=0.05, wb=net.wb, E=1.0))
    net.add(LineDC("CAB", 10, 11, R=0.01, X=0.3, wb=net.wb))
    net.add(ChargeDC_P("PDC", 11, P=0.4))
    m = ModeleNonLineaire(net)
    dx, dzu = m.residu_equilibre()
    assert dx < 1e-8, f"|dx/dt| = {dx:.3e}"
    assert dzu < 1e-8, f"|z - z0| + |u - u0| = {dzu:.3e}"


def test_residu_equilibre_liaison_mixte():
    """Le decoupage continu n'AJOUTE aucun residu au reseau alternatif.

    ``|dx/dt|`` ne descend pas a 1e-8 sur ce reseau-ci, et ce n'est pas le
    fait du continu : les deux noeuds alternatifs y ont ``C = 1e-3``, donc
    ``dv/dt = (wb/C).isum`` amplifie de 3,1e5 un residu de courant de 1e-13.
    C'est la raison pour laquelle le test alternatif equivalent
    (``test_equilibre_est_un_point_fixe_du_reseau_entier``) borne ce terme a
    1e-5 et pas a 1e-8.

    Ce qui SE mesure a 1e-8, et qui est ce que le continu doit garantir :
    la fermeture du cablage (``dzu``), la loi des noeuds continue, et le fait
    que le reseau decoupe rende le MEME residu que le ``GFM3`` monolithique,
    au bit pres.
    """
    ref = reseau_gfm3()
    ref.solve_equilibrium()
    net = equivalent_dc(ref)
    net.solve_equilibrium()

    dx_ref = ModeleNonLineaire(ref).residu_equilibre()[0]
    m = ModeleNonLineaire(net)
    dx, dzu = m.residu_equilibre()

    assert dzu < 1e-8, f"|z - z0| + |u - u0| = {dzu:.3e}"
    assert dx == pytest.approx(dx_ref, rel=1e-12), (dx, dx_ref)
    assert dx < 1e-5, f"|dx/dt| = {dx:.3e}"

    # la part continue du residu, elle, tient largement 1e-8
    c = m.c
    xdot, _, _ = m.derivee(c.x0(), m.uext0(), c.z0(), c.u0())
    continus = [
        abs(v)
        for n, v in zip(c.state_names, xdot)
        if n in ("bus10.v", "SDC.i") or n.startswith("CONV.")
    ]
    assert max(continus) < 1e-8, f"pire derivee continue {max(continus):.3e}"


# ----------------------------------------------------------------------
# 6) les indices de port, figes en litteraux
# ----------------------------------------------------------------------
@pytest.mark.parametrize(
    "cls, entrees, sorties",
    [
        (
            Bus,
            ("isum_d", "isum_q", "w_ref", "theta_g"),
            ("vd_o", "vq_o", "vmag"),
        ),
        (
            Branch,
            ("w_ref", "theta_g", "vd_f", "vq_f", "vd_t", "vq_t"),
            ("id_f", "iq_f", "id_t", "iq_t", "p_f", "q_f"),
        ),
        (
            IB,
            ("E_ref", "delta_ref", "w_set", "w_ref", "theta_g", "vd_t", "vq_t"),
            ("id_t", "iq_t", "P", "Q", "Vt", "theta_out", "omega"),
        ),
        (
            GFM3,
            ("p_ref", "q_ref", "ve_ref", "vdc_ref", "w_ref", "theta_g", "vd_t", "vq_t"),
            ("id_t", "iq_t", "P", "Q", "w", "Vt", "theta_out", "vdc_out"),
        ),
        (
            GFM3DC,
            ("p_ref", "q_ref", "ve_ref", "w_ref", "theta_g", "vd_t", "vq_t", "v_dc"),
            ("id_t", "iq_t", "i_dc", "P", "Q", "w", "Vt", "theta_out", "vdc_out"),
        ),
        (
            BusDC,
            ("isum", "w_ref", "theta_g"),
            ("v_o",),
        ),
        (
            LineDC,
            ("w_ref", "theta_g", "v_f", "v_t"),
            ("i_f", "i_t", "p_f"),
        ),
        (
            SourceDC,
            ("E_ref", "w_ref", "theta_g", "v_t"),
            ("i_t", "P", "Vt"),
        ),
    ],
)
def test_signaux_figes(cls, entrees, sorties):
    """Les noms ET LES POSITIONS des signaux, ecrits en dur.

    Ce test est la pour CASSER. Les ports continus s'ajoutent en queue
    justement pour qu'aucun signal existant ne bouge ; le jour ou quelqu'un
    en insere un au milieu, c'est ici que cela doit se voir, et pas dans
    trois modeles silencieusement decales de deux lignes.
    """
    assert cls.input_names() == entrees
    assert cls.output_names() == sorties


def test_n_port_signaux():
    """2 par port alternatif, 1 par port continu -- la frontiere des sorties."""
    assert GFM3.n_port_signaux() == 2
    assert GFM3DC.n_port_signaux() == 3
    assert LineDC.n_port_signaux() == 2
    assert BusDC.n_port_signaux() == 0
    # la frontiere tombe bien la ou commencent les sorties de service
    assert GFM3DC.output_names()[GFM3DC.n_port_signaux()] == "P"
    assert GFM3DC.input_names()[GFM3DC.port_dc_input_idx(0)] == "v_dc"
    assert GFM3DC.output_names()[GFM3DC.port_dc_output_idx(0)] == "i_dc"
    # Ajouter un port continu ne deplace RIEN du cote alternatif : a entrees
    # de service egales (GFM1 en a les memes trois que GFM3DC), les indices
    # du port alternatif, de w_ref et de theta_g sont identiques.
    assert GFM3DC.port_input_idx(0) == GFM1.port_input_idx(0) == (5, 6)
    assert GFM3DC.port_output_idx(0) == GFM1.port_output_idx(0) == (0, 1)
    assert GFM3DC.wref_idx() == GFM1.wref_idx() == 3
    assert GFM3DC.thetaref_idx() == GFM1.thetaref_idx() == 4


# ----------------------------------------------------------------------
# 7) un port sur le mauvais genre de bus
# ----------------------------------------------------------------------
def test_port_alternatif_sur_bus_continu():
    net = Network(f0=50.0, Sb=1e6)
    net.add(BusDC(10, 640e3, C=0.05))
    with pytest.raises(ValueError, match=r"L1.*\b10\b.*continu"):
        net.add(Line("L1", 10, 10, R=0.01, X=0.1, wb=net.wb))


def test_port_continu_sur_bus_alternatif():
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=1, V=1.0, C=1e-3)
    with pytest.raises(ValueError, match=r"SDC.*\b1\b.*alternatif"):
        net.add(SourceDC("SDC", 1, R=0.01, X=0.05, wb=net.wb))


def test_numero_de_bus_partage_refuse():
    """Le cablage repere un noeud par son seul numero : pas de doublon de genre."""
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=1, V=1.0, C=1e-3)
    with pytest.raises(ValueError, match="deja declare comme alternatif"):
        net.add(BusDC(1, 640e3, C=0.05))


def test_bus_continu_sans_capacite_refuse():
    net = Network(f0=50.0, Sb=1e6)
    with pytest.raises(ValueError, match="sans capacite"):
        net.add(BusDC(10, 640e3, C=0.0))


def test_genre_verifie_aussi_au_cablage():
    """Le cablage refait la verification, pour un reseau bati autrement.

    ``add`` la fait deja, mais c'est ``cablage`` qui ecrirait la matrice
    fausse : la garde doit etre la ou est le danger.
    """
    net = reseau_rc()
    charge = net.elements[0]
    net.add_bus(1, 20e3, btype=1, V=1.0, C=1e-3)
    charge.buses[0] = 1  # on force le port continu sur un noeud alternatif
    with pytest.raises(ValueError, match=r"RL.*\b1\b.*alternatif"):
        net.cablage()


# ----------------------------------------------------------------------
# 8) i_dc n'est pas une sortie exposee
# ----------------------------------------------------------------------
def test_courant_continu_absent_des_sorties():
    """``i_dc`` est un signal INTERNE du cablage, comme ``id_t``/``iq_t``.

    S'il fuyait dans les sorties exposees, ce ne serait pas seulement une
    sortie de trop : toutes les lignes suivantes de la matrice C seraient
    decalees.
    """
    ref = reseau_gfm3()
    ref.solve_equilibrium()
    net = equivalent_dc(ref)
    ss = net.build_state_space()
    for nom in ss.outputs:
        assert not nom.endswith(".i_dc"), nom
        assert not nom.endswith(".id_t"), nom
        assert not nom.endswith(".i_t"), nom
    # ... mais la tension continue, elle, est bien exposee
    assert "CONV.vdc_out" in ss.outputs
    assert "bus10.v_o" in ss.outputs
    assert "SDC.P" in ss.outputs


def test_sorties_exposees_du_reseau_continu():
    """Un reseau purement continu expose ses tensions et ses puissances."""
    net = reseau_cable()
    ss = net.build_state_space()
    assert ss.outputs == ["bus10.v_o", "bus11.v_o", "CAB.p_f"]


# ----------------------------------------------------------------------
# 9) les fichiers interdits n'ont pas bouge
# ----------------------------------------------------------------------
def _git_show(chemin: str) -> str:
    return subprocess.run(
        ["git", "show", f"HEAD:{chemin}"],
        cwd=RACINE,
        capture_output=True,
        text=True,
        check=True,
    ).stdout


def _source_fonction(texte: str, nom: str) -> str:
    """Le code source d'une fonction nommee, extrait d'un fichier Python."""
    arbre = ast.parse(texte)
    for noeud in ast.walk(arbre):
        if isinstance(noeud, ast.FunctionDef) and noeud.name == nom:
            return ast.get_source_segment(texte, noeud)
    raise AssertionError(f"fonction '{nom}' introuvable")


# Deux garde-fous pinaient ici le CONTENU de ``loadflow.py``, de ``tables.py``
# et de ``run_loadflow`` sur celui de HEAD. Ils ont fait leur office -- ils ont
# tenu le lot des ports continus a distance du load flow -- mais ils mesuraient
# le texte, pas le comportement, et se comparaient a HEAD : ils tombaient donc
# au rouge devant toute modification non commitee de ces fichiers, y compris
# celle que le load flow AC/DC doit legitimement y faire (envelopper le corps
# alternatif de ``run_loadflow`` dans un ``if not self.bus_specs_dc``).
#
# La garantie ne disparait pas, elle CHANGE DE NATURE : ce qui doit etre
# preserve n'est pas le texte du load flow alternatif, c'est son resultat. Le
# controle vit desormais dans ``tests/test_loadflow_acdc.py``, qui rejoue les
# reseaux du catalogue et compare les tensions a une reference, a 1e-11 pu
# (FR-001 : l'octet ne survivait pas d'une machine a l'autre).
#
# ``build_ybus`` a ete figee ici de la meme facon, par son texte, tant qu'aucun
# lot n'avait de raison d'y toucher. Les reactances shunt en ont une (elles
# inscrivent leur admittance au Y-bus, ``Element.admittance_ybus``) : la garde
# passe donc, elle aussi, du texte au comportement.
def test_build_ybus_intacte():
    """``build_ybus`` ne voit rien du continu, et ne doit rien en voir.

    ``bus_specs_dc`` est volontairement SEPARE de ``bus_specs`` pour que la
    construction de la matrice d'admittance alternative ne puisse pas, meme
    par accident, rencontrer un bus continu -- elle y leverait un ``KeyError``,
    ou pire, l'integrerait silencieusement. On le verifie sur un reseau mixte :
    le Y-bus a exactement la taille du reseau alternatif, et aucun element
    continu ne declare d'admittance a y inscrire.
    """
    import g2egrid.catalogue as cat
    from g2egrid.tables import assembler

    net = assembler(cat.charger("hvdc2t"))
    Y = net.build_ybus()
    assert Y.shape == (len(net.bus_specs),) * 2
    assert net.bus_specs_dc and not set(net.bus_specs) & set(net.bus_specs_dc)
    for el in net.elements:
        if el.ports_dc and not el.ports:
            assert el.shunt_admittance() is None, el.name
            assert el.admittance_ybus() is None, el.name


# ----------------------------------------------------------------------
# garde-fous du lot 4
# ----------------------------------------------------------------------
def test_modeles_compiles_distincts():
    """Chaque variante a port DC compile SOUS SON PROPRE NOM.

    ``compile_spec`` met en cache par ``spec.name`` : deux classes qui
    partageraient ce nom se rendraient mutuellement le modele de l'autre --
    silencieusement, avec les bons parametres et les mauvaises equations.
    """
    noms = [c.spec().name for c in (GFM3, GFM3DC, GFM4, GFM4DC, GFL2DC)]
    assert len(set(noms)) == len(noms), noms
    assert GFM3DC.compiled() is not GFM3.compiled()
    assert GFM3DC.compiled().n_x == 13 and GFM3.compiled().n_x == 15
    # les noms de signaux ne se marchent pas dessus non plus : ModelSpec
    # refuse un meme nom dans deux groupes, et v_dc / i_dc ne sont pas vdc
    spec = GFM3DC.spec()
    tous = spec.params + spec.states + spec.algebs + spec.inputs + spec.outputs
    assert len(set(tous)) == len(tous)


def test_elements_continus_sans_shunt_admittance():
    """Aucun element continu ne declare de shunt : le Y-bus le lirait.

    ``build_ybus`` indexe par ``self._bus_index()``, qui ne connait que les
    noeuds ALTERNATIFS : un shunt declare par un element continu y leverait un
    ``KeyError`` sur son propre numero de bus.
    """
    wb = 2 * np.pi * 50
    elements = [
        BusDC(10, 640e3, C=0.05, wb=wb),
        LineDC("CAB", 10, 11, R=0.01, X=0.3, wb=wb),
        SourceDC("SDC", 10, R=0.01, X=0.05, wb=wb),
        ChargeDC_R("RL", 10, R=2.0),
        ChargeDC_P("PDC", 10, P=0.3),
    ]
    for el in elements:
        assert el.shunt_admittance() is None, el.name


def test_capacite_du_bus_continu_vaut_celle_du_convertisseur():
    """Le noeud continu EST le condensateur du VSC : c'est la regle d'assemblage."""
    ref = reseau_gfm3()
    ref.solve_equilibrium()
    net = equivalent_dc(ref)
    assert net.bus_specs_dc[10].C == GFM_DC_DEFAULT["Cdc_pu"]
    assert 10 not in net.bus_specs  # invisible du load flow alternatif
    assert set(net._bus_index()) == {1, 2}


# ----------------------------------------------------------------------
# les mesures physiques sur un reseau continu
# ----------------------------------------------------------------------
def test_mesures_reconnaissent_un_noeud_continu():
    """``BusDC`` derive de ``Bus`` : le lire comme un noeud alternatif serait
    une erreur SILENCIEUSE.

    Un noeud alternatif a deux sorties (``vd_o``, ``vq_o``), un noeud continu
    une seule. Prendre ``y[sy], y[sy + 1]`` sur un noeud continu lit donc la
    premiere sortie de l'element SUIVANT comme son ``vq`` : pas une exception,
    un module de tension faux. D'ou ``Vdc_busN``, et ni angle ni frequence --
    un continu n'a pas de phase.
    """
    ref = reseau_gfm3()
    ref.solve_equilibrium()
    net = equivalent_dc(ref)
    net.solve_equilibrium()
    net.build_state_space()
    m = ModeleNonLineaire(net)

    plan = _plan_mesures(m)
    assert "Vdc_bus10" in plan
    for interdit in ("V_bus10", "angle_bus10", "f_bus10", "fm_bus10"):
        assert interdit not in plan, f"{interdit} n'a pas de sens sur un noeud continu"

    tr = m.simuler((0.0, 0.02), t_eval=np.linspace(0.0, 0.02, 11))
    v = mesures(m, tr, quoi=["Vdc_bus10", "P_CONV_dc", "P_SDC_dc", "P_CONV"])

    # la tension continue mesuree EST celle du point de linearisation
    bus = next(e for e in m.c.els if getattr(e, "dc", False))
    assert v["Vdc_bus10"][0] == pytest.approx(bus.x0[0], rel=1e-12)

    # loi des noeuds continue : ce que la source injecte, le convertisseur
    # le prend -- au signe pres et a la precision du point de fonctionnement
    assert v["P_SDC_dc"][0] == pytest.approx(-v["P_CONV_dc"][0], rel=1e-9)

    # et les pertes du convertisseur sont du bon cote : il tire du continu
    # un peu plus qu'il ne rend en alternatif
    assert -v["P_CONV_dc"][0] > v["P_CONV"][0] > 0


# ----------------------------------------------------------------------
# convertisseur DC/DC : un PV tient sa puissance, une batterie sa tension
# ----------------------------------------------------------------------
def _lf_mtdc(**tables):
    from g2egrid.catalogue import charger
    from g2egrid.tables import assembler

    net = assembler(charger("mtdc3t").copie(**tables))
    return net.run_loadflow()


def test_dcdc_pv_equivaut_a_une_charge_negative():
    pv = _lf_mtdc(T_DCDC=[(1, 12, 5e4, "pv", 2e4)])
    charge = _lf_mtdc(T_LoadsDC=[(1, 12, -2e4, "P")])
    assert pv.converged and charge.converged
    assert np.allclose(pv.Vdc, charge.Vdc, atol=1e-9)
    assert np.allclose(pv.V, charge.V, atol=1e-9)


def test_dcdc_batterie_a_vide_ne_deplace_pas_le_point():
    # la batterie est une reference secondaire : le VSC en Vdc reste la
    # reference de l'ile, et une batterie sans consigne n'injecte rien
    base = _lf_mtdc()
    bat = _lf_mtdc(T_DCDC=[(1, 12, 5e4, "batterie", 0.0)])
    assert bat.converged
    assert np.allclose(bat.Vdc, base.Vdc, atol=1e-6)


def test_charge_continue_resistive_suit_la_tension():
    lf = _lf_mtdc(T_LoadsDC=[(1, 12, 1e4, "R")])
    assert lf.converged


def test_dcdc_type_inconnu_refuse():
    with pytest.raises(ValueError, match="inconnu"):
        _lf_mtdc(T_DCDC=[(1, 12, 5e4, "eolien", 0.0)])
