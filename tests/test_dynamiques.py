"""Selection des dynamiques, resolution globale, repere flottant, sensibilites.

Quatre idees reprises du travail de F. Kelada, chacune verifiee par une
EQUIVALENCE exacte plutot que par une valeur de reference :

1. supprimer une dynamique (``dx/dt = 0`` ou ``x = x0``) fait des variantes
   GFM1..3 des SELECTIONS du modele complet GFM4 -- memes spectres ;
2. la resolution globale reduit l'indice a un noeud ideal : deux branches en
   serie ont le spectre d'une branche unique ;
3. un repere flottant (centre d'inertie) donne le spectre du repere lie au
   slack, a un mode de jauge pres -- et survit au declenchement du slack ;
4. la sensibilite analytique d'un mode a un parametre egale la difference
   finie, point de fonctionnement fixe.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2egrid import Network, catalogue
from g2egrid.elements import (
    GFM1,
    GFM2,
    GFM3,
    GFM4,
    IdealBus,
    Line,
    SynchronousMachine,
    Transformer,
    VoltageSource,
)
from g2egrid.modal import modes
from g2egrid.sensibilites import sensibilites
from g2egrid.tables import assembler, tables_depuis_json, tables_vers_json


def _ecart(a, b) -> float:
    la, lb = np.sort_complex(a.eig()), np.sort_complex(b.eig())
    assert len(la) == len(lb)
    return float(np.max(np.abs(la - lb) / np.maximum(1.0, np.abs(la))))


def _smib_convertisseur(fabrique):
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=2, V=1.0, Pg=0.5, C=1e-3)
    net.add_bus(2, 20e3, btype=1, V=1.0, C=1e-3)
    net.add(fabrique(net))
    net.add(Line("L1", 1, 2, R=0.01, X=0.2, wb=net.wb))
    net.add(VoltageSource("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
    net.run_loadflow()
    return net


def _gfm(cls, **red):
    def fab(net):
        el = cls("G", 1, Sn=1e6, Sb=1e6, wb=net.wb)
        return el.reduire(**red) if red else el

    return fab


# ------------------------------------------------------ 1. dynamiques
@pytest.mark.parametrize(
    "variante, reduction",
    [
        (GFM1, dict(statique=["Dw"], fige=["vdc", "idc"])),
        (GFM2, dict(fige=["vdc", "idc"])),
        (GFM3, dict(statique=["Dw"])),
    ],
)
def test_les_variantes_gfm_sont_des_selections_de_gfm4(variante, reduction):
    ref = _smib_convertisseur(_gfm(variante)).build_state_space()
    red = _smib_convertisseur(_gfm(GFM4, **reduction)).build_state_space()
    assert len(ref.states) == len(red.states)
    assert _ecart(ref, red) < 1e-8


def test_reduire_change_la_classe_sans_perdre_la_famille():
    el = GFM4("G", 1, Sn=1e6, Sb=1e6, wb=314.0).reduire(
        statique="Dw", fige={"vdc": 1.0}
    )
    assert isinstance(el, GFM4)
    assert "Dw" in el.algebs and "Dw" not in el.states
    assert "vdc" in el.params and el.pdict["vdc"] == 1.0
    assert el.dynamiques()["idc"] == "dynamique"
    # les reductions se cumulent, et deux elements reduits pareil partagent
    # leur classe (donc leur modele compile)
    el.reduire(fige="idc")
    autre = GFM4("H", 1, Sn=1e6, Sb=1e6, wb=314.0).reduire(
        fige=["idc", "vdc"], statique=["Dw"]
    )
    assert type(el) is type(autre)


def test_reduire_refuse_un_etat_inconnu():
    el = GFM1("G", 1, Sn=1e6, Sb=1e6, wb=314.0)
    with pytest.raises(ValueError, match="Dw"):
        el.reduire(statique=["Dw"])
    assert el.reduire(statique=["Dw"], strict=False) is el  # rien a faire


def test_dynamiques_par_famille_et_par_bus_dans_les_tables():
    t = catalogue.charger("wscc9_gfm")
    ref = assembler(t)
    ref.run_loadflow()
    t4 = t.copie(
        T_GFM1s=[],
        T_GFM4s=list(t.T_GFM1s),
        GFM_PARAMS={"dynamiques": {"Dw": "statique", "vdc": "fige", "idc": "fige"}},
    )
    net = assembler(t4)
    net.run_loadflow()
    assert _ecart(ref.build_state_space(), net.build_state_space()) < 1e-8

    # un bus rend a son GFM la dynamique que la famille supprime
    b = int(t.T_GFM1s[-1][1])
    t5 = t4.copie(PARAMS_PAR_BUS={b: {"dynamiques": {"Dw": "dynamique"}}})
    net5 = assembler(t5)
    sel = {e.name: e.dynamiques()["Dw"] for e in net5.elements if isinstance(e, GFM4)}
    assert sel[f"GFM{b}"] == "dynamique"
    assert sorted(sel.values()) == ["dynamique", "statique"]

    with pytest.raises(ValueError, match="Dww"):
        assembler(t.copie(GFM_PARAMS={"dynamiques": {"Dww": "statique"}}))


# ------------------------------------------------ 2. resolution globale
def _serie(noeuds="ideal", en_serie=True, C3=0.0, **kw):
    net = Network(f0=50.0, Sb=1e6, noeuds_sans_capacite=noeuds, **kw)
    net.add_bus(1, 20e3, btype=2, V=1.0, Pg=0.8, C=1e-3)
    net.add_bus(2, 20e3, btype=1, V=1.0, C=1e-3)
    net.add(SynchronousMachine("SM1", 1, Sn=1e6, Sb=1e6, wb=net.wb))
    if en_serie:
        net.add_bus(3, 20e3, btype=3, C=C3)
        net.add(Line("L1", 1, 3, R=0.004, X=0.08, wb=net.wb))
        net.add(Transformer("T1", 3, 2, R=0.006, X=0.12, wb=net.wb))
    else:
        net.add(Line("L1", 1, 2, R=0.01, X=0.20, wb=net.wb))
    net.add(VoltageSource("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
    net.run_loadflow()
    return net


def test_deux_branches_en_serie_a_un_noeud_ideal_font_une_branche():
    unique = _serie(en_serie=False).build_state_space()
    net = _serie()
    assert isinstance(net.buses[3], IdealBus)
    serie = net.build_state_space()
    # deux courants de branche de plus, et deux contraintes i1 = i2 : le
    # modele reduit a exactement autant d'etats que la branche unique
    assert len(serie.elimines) == 2 and len(serie.states) == len(unique.states)
    assert _ecart(unique, serie) < 1e-9
    # x = T.x_r reconstruit les etats elimines : T est de rang plein
    assert serie.T.shape == (len(serie.states) + 2, len(serie.states))


def test_le_noeud_ideal_est_la_limite_d_une_capacite_qui_s_annule():
    ideal = _serie().build_state_space().eig()
    petit = _serie(noeuds="resistance", C3=1e-7).build_state_space().eig()
    lents = np.sort_complex(ideal[np.abs(ideal) < 1e3])
    petits = np.sort_complex(petit[np.abs(petit) < 1e3])
    assert np.max(np.abs(lents - petits)) < 1e-4


def test_le_dae_global_rassemble_toutes_les_equations():
    dae = _serie().dae_lineaire()
    assert "bus3.vd" in dae.algebriques and "SM1.i_d" in dae.algebriques
    assert dae.BZ.shape == (len(dae.algebriques), len(dae.algebriques))
    assert np.linalg.matrix_rank(dae.BZ) < len(dae.algebriques)  # indice 2


def test_participations_dans_la_base_complete():
    """Deux branches en serie a un noeud ideal : dans la base complete, les
    participations sont la limite R -> infini de celles du modele a
    resistance de causalite (memes etats, memes modes), et somment a 1."""
    from g2egrid import modal

    net = _serie()
    net.solve_equilibrium()
    ss = net.build_state_space()
    assert ss.elimines and ss.T_gauche is not None
    assert np.allclose(ss.T_gauche @ ss.T, np.eye(ss.T.shape[1]))
    lam, _V, p, noms = modal.spectre_complet(ss)
    assert noms == ss.etats_complets and p.shape == (lam.size, len(noms))
    assert np.allclose(p.sum(axis=1), 1.0)

    lim = _serie(noeuds="resistance", r_causalite=1e4)
    lim.solve_equilibrium()
    lam_r, _, p_r, noms_r = modal.spectre_complet(lim.build_state_space())
    assert sorted(noms_r) == sorted(noms)
    ordre = [noms_r.index(n) for n in noms]
    for k in np.flatnonzero((np.abs(lam) > 1e-6) & (lam.imag >= 0)):
        j = int(np.argmin(np.abs(lam_r - lam[k])))
        assert abs(lam_r[j] - lam[k]) < 1e-3 * abs(lam[k])
        assert np.abs(p_r[j][ordre] - p[k]).max() < 1e-3


def test_le_pivot_de_la_reduction_ne_depend_pas_de_l_arrondi():
    """Deux courants de meme poids dans une contrainte : le premier dans
    l'ordre des etats est elimine, que l'arrondi penche d'un cote ou de
    l'autre ; une colonne nettement plus faible n'est jamais prise."""
    from g2egrid.network import _pivots_stables

    for eps in (1e-15, -1e-15):
        assert _pivots_stables(np.array([[1.0, 1.0 + eps, 0.0]])) == [0]
    assert _pivots_stables(np.array([[0.01, 1.0, 0.5]])) == [1]
    K = np.array([[1.0, -1.0, 0.0, 0.0], [0.0, 0.0, 1.0, -1.0]])
    assert _pivots_stables(K) == [0, 2]


def test_le_noeud_ideal_s_integre_en_bdf2_comme_la_branche_unique():
    """Radau resout l'algebrique a ``x`` fige, ce qu'un noeud ideal d'indice
    2 n'a pas : il est refuse. BDF2, choisi d'office, integre les deux
    branches en serie comme la branche unique equivalente."""
    from g2egrid.temporel import ModeleNonLineaire

    trajectoires = []
    for en_serie in (True, False):
        net = _serie(en_serie=en_serie)
        net.solve_equilibrium()
        m = ModeleNonLineaire(net)
        assert m.ideaux == (["bus3"] if en_serie else [])
        if en_serie:
            with pytest.raises(ValueError, match="BDF2"):
                m.simuler((0.0, 0.1), methode="Radau")
        u = m.echelon("SM1.p_ref", 0.05, 0.01)
        trajectoires.append(m.simuler((0.0, 0.3), uext_fn=u, methode="BDF2", pas=1e-4))
    serie, unique = trajectoires
    w_s, w_u = serie.etat("SM1.Dwr"), unique.etat("SM1.Dwr")
    assert np.abs(w_u).max() > 1e-4  # l'echelon fait bouger la machine
    assert np.abs(w_s - w_u).max() < 1e-6 * np.abs(w_u).max() + 1e-12


# ------------------------------------------------ 3. repere flottant
@pytest.mark.parametrize("nom", ["wscc9", "wscc9_gfm"])
def test_le_repere_flottant_ne_change_pas_la_physique(nom):
    lie = catalogue.construire(nom).build_state_space()
    net = catalogue.construire(nom)
    net.reference_flottante("coi")
    flottant = net.build_state_space()
    assert len(flottant.states) == len(lie.states) + 1  # l'angle du repere
    la, lb = lie.eig(), flottant.eig()
    la, lb = (
        np.sort_complex(la[np.abs(la) > 1e-6]),
        np.sort_complex(lb[np.abs(lb) > 1e-6]),
    )
    assert len(la) == len(lb)
    assert np.max(np.abs(la - lb) / np.maximum(1.0, np.abs(la))) < 1e-9
    # un mode de jauge de plus : le repere a son propre angle absolu
    assert (np.abs(flottant.eig()) < 1e-6).sum() == (np.abs(lie.eig()) < 1e-6).sum() + 1


def test_le_centre_d_inertie_pondere_les_machines_par_leur_inertie():
    net = catalogue.construire("wscc9")
    net.reference_flottante("coi")
    poids = net._poids_suivis()
    assert abs(sum(poids.values()) - 1.0) < 1e-12
    sm = {e.name: e for e in net.elements if e.name in poids}
    brut = {n: 2 * e.pdict["H"] * e.pdict["kS"] for n, e in sm.items()}
    tot = sum(brut.values())
    for n in poids:
        assert abs(poids[n] - brut[n] / tot) < 1e-12


@pytest.mark.parametrize("reference", ["coi", "nominale"])
def test_le_slack_peut_etre_declenche_sous_repere_flottant(reference):
    t = catalogue.charger("wscc9")
    slack = next(r for r in t.T_SMs if int(t.bus(int(r[1]))[1]) == 1)
    t2 = t.copie(REFERENCE=reference, OUVERTS=[("T_SMs", slack[0])])
    net = assembler(t2)
    net.run_loadflow()
    assert net.reference == "Repere"
    ss = net.build_state_space()
    lam = ss.eig()
    assert np.max(lam.real[np.abs(lam) > 1e-6]) < 0
    assert f"SM{int(slack[1])}.theta" not in ss.states


def test_reference_et_noeuds_passent_par_le_json():
    # ieee14 n'est pas un reseau de parite MATLAB : il suit les defauts
    t = catalogue.charger("ieee14").copie(REFERENCE="coi", NOEUDS="resistance")
    t2 = tables_depuis_json(tables_vers_json(t))
    assert (t2.REFERENCE, t2.NOEUDS) == ("coi", "resistance")
    d = tables_vers_json(catalogue.charger("ieee14"))
    assert "REFERENCE" not in d and "NOEUDS" not in d  # defauts : rien d'ecrit


# ------------------------------------------------ 4. sensibilites
@pytest.mark.parametrize("nom", ["smib_gfm", "wscc9"])
def test_sensibilite_analytique_egale_difference_finie(nom):
    net = catalogue.construire(nom)
    ss = net.build_state_space()
    m = next(m for m in modes(ss) if m.oscillant and m.f < 5)
    ana = sensibilites(net, m)
    top = ana[:6]
    fd = sensibilites(
        net, m, methode="differences", elements=sorted({s.element for s in top})
    )
    par = {(s.element, s.parametre): s.dlambda for s in fd}
    for s in top:
        d = par[(s.element, s.parametre)]
        assert abs(s.dlambda - d) <= 1e-4 * max(abs(s.dlambda), 1e-3), s


def test_sensibilite_d_un_etat_fige():
    """Un etat fige devient un parametre : sa sensibilite se lit comme une autre."""
    net = _smib_convertisseur(_gfm(GFM4, fige=["vdc", "idc"]))
    ss = net.build_state_space()
    k = int(np.argmax(ss.eig().real + 1e3 * (np.abs(ss.eig().imag) > 1)))
    s = sensibilites(net, k, elements=["G"], parametres=["vdc"])
    assert len(s) == 1 and s[0].valeur == pytest.approx(1.0)
    fd = sensibilites(net, k, methode="differences", elements=["G"], parametres=["vdc"])
    assert abs(s[0].dlambda - fd[0].dlambda) <= 1e-4 * max(abs(s[0].dlambda), 1e-3)


def test_sensibilite_refusee_sur_un_reseau_reduit_en_indice():
    net = _serie()
    with pytest.raises(NotImplementedError, match="differences"):
        sensibilites(net, 0)


def test_l_api_rend_les_parametres_qui_deplacent_un_mode():
    from g2egrid_api import analyse

    t = catalogue.charger("wscc9")
    m = next(m for m in analyse.modes(t)["modes"] if abs(m["im"]) > 1 and m["f"] < 5)
    d = analyse.detail_mode(t, m["indice"], n_parametres=5)
    assert len(d["parametres"]) == 5
    assert {"element", "parametre", "valeur", "re", "im", "df", "dzeta"} <= set(
        d["parametres"][0]
    )
    assert analyse.detail_mode(t, m["indice"], n_parametres=0)["parametres"] == []
