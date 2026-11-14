"""Tests de la simulation temporelle non lineaire (``g2egrid.temporel``).

Le fil conducteur est le meme que pour le reste du paquet : confronter la voie
non lineaire a une reference qui ne partage pas son code. Le point de
fonctionnement du load flow sert de reference aux mesures, le modele d'etat
lineaire sert de reference a la linearisation instantanee, et la superposition
au premier ordre sert de reference aux scenarios composes.
"""

import types

import numpy as np
import pytest

from g2egrid import Network
from g2egrid.catalogue import construire
from g2egrid.elements import IB, Line, SynchronousMachine
from g2egrid.temporel import (
    Evenement,
    ModeleNonLineaire,
    Trajectoire,
    _premier_ordre,
    linearisation,
    mesures,
)


# ----------------------------------------------------------------------
# fixtures
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def net():
    """SMIB du catalogue : 28 etats, le plus rapide a simuler."""
    return construire("smib")


@pytest.fixture(scope="module")
def ss(net):
    return net.build_state_space()


@pytest.fixture(scope="module")
def m(net):
    return ModeleNonLineaire(net)


@pytest.fixture(scope="module")
def plate(m):
    """Trajectoire NON PERTURBEE : le reseau doit y rester ou il est.

    C'est la reference des mesures : sans perturbation, tout ce qui en sort
    doit redonner le load flow, coefficient par coefficient.
    """
    return m.simuler((0.0, 0.02), t_eval=np.array([0.0, 0.01, 0.02]))


@pytest.fixture(scope="module")
def perturbee(m):
    """Echelon de 5 % sur ``p_ref`` a 50 ms, une demi-seconde de transitoire."""
    uext_fn = m.scenario(Evenement(0.05, "SM1.p_ref", 0.05))
    tr = m.simuler((0.0, 0.5), uext_fn=uext_fn, t_eval=np.linspace(0.0, 0.5, 51))
    return tr, uext_fn


# ======================================================================
# B. evenements et scenarios
# ======================================================================
def test_evenement_echelon():
    ev = Evenement(0.1, "SM1.p_ref", 0.05)
    assert ev.valeur(0.0) == 0.0
    assert ev.valeur(0.099) == 0.0
    assert ev.valeur(0.1) == 0.05  # l'echelon appartient a l'instant qui le porte
    assert ev.valeur(10.0) == 0.05
    assert ev.discontinuites() == (0.1,)


def test_evenement_rampe():
    """La rampe monte lineairement, et ses DEUX bouts sont des cassures."""
    ev = Evenement(1.0, "SM1.p_ref", 0.2, genre="rampe", duree=0.5)
    assert ev.valeur(0.99) == 0.0
    assert ev.valeur(1.0) == 0.0
    assert ev.valeur(1.25) == pytest.approx(0.1)
    assert ev.valeur(1.5) == pytest.approx(0.2)
    assert ev.valeur(9.0) == pytest.approx(0.2)
    assert ev.discontinuites() == (1.0, 1.5)


def test_evenement_refuse_ce_qui_n_a_pas_de_sens():
    with pytest.raises(ValueError, match="genre"):
        Evenement(0.0, "SM1.p_ref", 1.0, genre="creneau")
    with pytest.raises(ValueError, match="duree"):
        Evenement(0.0, "SM1.p_ref", 1.0, genre="rampe")


def test_scenario_compose_plusieurs_entrees(m):
    """Les evenements s'ajoutent, y compris sur des entrees differentes."""
    a = Evenement(0.1, "SM1.p_ref", 0.01)
    b = Evenement(0.5, "SM1.v_ref", 0.005, genre="rampe", duree=0.2)
    uext_fn = m.scenario(a, b)
    assert uext_fn.instants == (0.1, 0.5, 0.7)

    ip = m.c.iu("SM1.p_ref")
    iv = m.c.iu("SM1.v_ref")
    u0 = m.uext0()
    assert uext_fn(0.0) == pytest.approx(u0)
    assert uext_fn(0.2)[ip] == pytest.approx(u0[ip] + 0.01)
    assert uext_fn(0.2)[iv] == pytest.approx(u0[iv])
    assert uext_fn(0.6)[iv] == pytest.approx(u0[iv] + 0.0025)
    assert uext_fn(1.0)[iv] == pytest.approx(u0[iv] + 0.005)
    assert uext_fn(1.0)[ip] == pytest.approx(u0[ip] + 0.01)


def test_deux_evenements_sur_la_meme_entree_se_cumulent(m):
    """Un creneau, c'est un echelon et son retrait."""
    uext_fn = m.scenario(
        Evenement(0.1, "SM1.p_ref", 0.02), Evenement(0.3, "SM1.p_ref", -0.02)
    )
    i = m.c.iu("SM1.p_ref")
    u0 = m.uext0()[i]
    assert uext_fn(0.0)[i] == pytest.approx(u0)
    assert uext_fn(0.2)[i] == pytest.approx(u0 + 0.02)
    assert uext_fn(0.4)[i] == pytest.approx(u0)
    assert uext_fn.instants == (0.1, 0.3)


def test_echelon_reste_un_cas_particulier_de_scenario(m):
    """La signature historique ne bouge pas, et gagne ses instants."""
    ancien = m.echelon("SM1.p_ref", 0.05, t0=0.1)
    nouveau = m.scenario(Evenement(0.1, "SM1.p_ref", 0.05))
    assert ancien.instants == (0.1,)
    for t in (0.0, 0.05, 0.1, 0.5):
        assert ancien(t) == pytest.approx(nouveau(t))


def test_scenario_refuse_une_entree_inconnue(m):
    with pytest.raises(ValueError):
        m.scenario(Evenement(0.1, "SM1.pas_une_entree", 1.0))


# ======================================================================
# B. decoupage aux discontinuites
# ======================================================================
def test_decoupage_aux_discontinuites_coute_moins_cher(m):
    """Meme trajectoire, moins de resolutions de Newton.

    Sans decoupage, Radau ne decouvre l'echelon qu'apres l'avoir franchi et
    ramene son pas par bissections jusqu'a le border ; chaque bissection coute
    une resolution du systeme couple, et ne sert a rien. Le test verifie les
    deux choses qui comptent : que le gain existe, et qu'il ne s'obtient pas
    en integrant autre chose.
    """
    t_eval = np.linspace(0.0, 0.5, 51)
    uext_fn = m.echelon("SM1.p_ref", 0.05, t0=0.1)
    coupe = m.simuler((0.0, 0.5), uext_fn=uext_fn, t_eval=t_eval)
    bloc = m.simuler((0.0, 0.5), uext_fn=uext_fn, t_eval=t_eval, instants=[])

    assert coupe.t == pytest.approx(bloc.t)
    ecart = np.abs(coupe.x - bloc.x).max() / np.abs(bloc.x - m.c.x0()).max()
    assert ecart < 1e-3, f"les deux integrations divergent ({ecart:.2e})"
    assert coupe.n_solves < bloc.n_solves / 1.3, (coupe.n_solves, bloc.n_solves)


def test_instants_explicites_valent_ceux_du_scenario(m):
    """``instants=`` doit marcher sur un ``uext`` quelconque, sans attribut."""
    t_eval = np.linspace(0.0, 0.3, 31)
    u0 = m.uext0()
    i = m.c.iu("SM1.p_ref")

    def brut(t):  # une fonction a la main, sans .instants
        u = u0.copy()
        if t >= 0.1:
            u[i] += 0.05
        return u

    nu = m.simuler((0.0, 0.3), uext_fn=brut, t_eval=t_eval)
    dit = m.simuler((0.0, 0.3), uext_fn=brut, t_eval=t_eval, instants=[0.1])
    assert dit.n_solves < nu.n_solves
    ampli = np.abs(nu.x - m.c.x0()).max()
    assert np.abs(nu.x - dit.x).max() / ampli < 1e-3


def test_decoupage_ne_pollue_pas_les_instants_rendus(m):
    """Les bornes d'intervalle ne doivent pas se glisser dans la trajectoire.

    Un point de jonction appartient aux deux intervalles : sans precaution il
    sortirait deux fois, et ``mesures`` rendrait des tableaux plus longs que
    ``t_eval``.
    """
    t_eval = np.linspace(0.0, 0.3, 16)
    uext_fn = m.scenario(
        Evenement(0.1, "SM1.p_ref", 0.01),
        Evenement(0.15, "SM1.v_ref", 0.002, genre="rampe", duree=0.05),
    )
    tr = m.simuler((0.0, 0.3), uext_fn=uext_fn, t_eval=t_eval)
    assert tr.t == pytest.approx(t_eval)

    # sans t_eval : les pas de l'integrateur, mais sans doublon a la jonction
    libre = m.simuler((0.0, 0.3), uext_fn=uext_fn)
    assert np.all(np.diff(libre.t) > 0.0)
    assert libre.t[0] == 0.0 and libre.t[-1] == pytest.approx(0.3)


# ======================================================================
# F. superposition
# ======================================================================
def test_superposition_de_deux_echelons_sur_deux_entrees(m):
    """A petite amplitude, deux perturbations s'ajoutent.

    C'est la definition meme du premier ordre : le residu de superposition est
    un terme en amplitude au CARRE. Le test le verifie deux fois -- l'ecart est
    petit, et il est multiplie par cent quand l'amplitude l'est par dix. La
    premiere assertion seule se satisferait d'un modele trop lineaire (un
    modele qui aurait perdu ses non-linearites), la seconde ne le tolererait
    pas.
    """
    t_eval = np.linspace(0.0, 0.6, 61)
    x0 = m.c.x0()

    def ecart_de_superposition(amp):
        a = Evenement(0.05, "SM1.p_ref", amp)
        b = Evenement(0.25, "SM1.v_ref", amp)
        xa = m.simuler((0.0, 0.6), uext_fn=m.scenario(a), t_eval=t_eval).x - x0
        xb = m.simuler((0.0, 0.6), uext_fn=m.scenario(b), t_eval=t_eval).x - x0
        xab = m.simuler((0.0, 0.6), uext_fn=m.scenario(a, b), t_eval=t_eval).x - x0
        return np.abs(xab - (xa + xb)).max(), np.abs(xab).max()

    petit, echelle = ecart_de_superposition(1e-3)
    assert petit / echelle < 1e-3, f"superposition violee des 1e-3 ({petit:.2e})"

    grand, _ = ecart_de_superposition(1e-2)
    assert 30.0 < grand / petit < 300.0, (
        f"l'ecart de superposition ne croit pas comme le carre de l'amplitude "
        f"(rapport {grand / petit:.1f} pour une amplitude x10)"
    )


# ======================================================================
# C. mesures
# ======================================================================
def test_mesures_reproduisent_le_load_flow(m, net, plate):
    """Au point de fonctionnement, les mesures DOIVENT redonner le load flow.

    Verification croisee decisive : le load flow resout un systeme de phaseurs
    en regime permanent, les mesures reconstruisent des grandeurs a partir des
    etats dq et des courants de port du modele dynamique. Rien de commun entre
    les deux chemins sinon le point de fonctionnement lui-meme. Le premier
    instant rendu porte exactement ``x0``, d'ou la tolerance a 1e-8.
    """
    mes = mesures(m, plate)
    idx = net._bus_index()
    for b, i in idx.items():
        V = net.lf.V[i]
        assert mes[f"V_bus{b}"][0] == pytest.approx(abs(V), abs=1e-8)
        assert mes[f"angle_bus{b}"][0] == pytest.approx(
            np.degrees(np.angle(V)), abs=1e-8
        )


def test_mesures_p_q_valent_les_injections_du_load_flow(m, net, plate):
    """``P_<element>`` et ``Q_<element>`` sont les injections attribuees par le
    load flow, convention du paquet : positif = fourni au reseau.

    ``Network._element_injections`` est justement l'attribution que le load
    flow fait a chaque element ; c'est la reference exacte, pas une
    approximation.
    """
    mes = mesures(m, plate)
    injections = net._element_injections()
    vus = 0
    for k, el in enumerate(net.elements):
        if len(el.ports) != 1 or k not in injections:
            continue
        S = injections[k]
        assert mes[f"P_{el.name}"][0] == pytest.approx(S.real, abs=1e-8), el.name
        assert mes[f"Q_{el.name}"][0] == pytest.approx(S.imag, abs=1e-8), el.name
        vus += 1
    assert vus >= 2, "aucun element a un port compare au load flow"


def test_mesures_transits_de_branche_bouclent_sur_les_pertes(m, net, plate):
    """``P_from + P_to = -R.|i|^2`` : la branche n'absorbe que ses pertes joule.

    C'est la verification interne qui remplace la comparaison au transit de
    load flow : celle-ci est hors de portee, puisque la susceptance des lignes
    n'est pas dans la branche mais repartie dans la capacite des bus, d'ou elle
    n'est plus separable ligne par ligne (cf. la docstring de ``mesures``).
    """
    mes = mesures(m, plate)
    vues = 0
    for el in net.elements:
        if len(el.ports) != 2:
            continue
        i = complex(*el.x0)
        pertes = -el.pdict["R"] * abs(i) ** 2
        somme = mes[f"P_{el.name}_from"][0] + mes[f"P_{el.name}_to"][0]
        assert somme == pytest.approx(pertes, abs=1e-9), el.name
        vues += 1
    assert vues >= 1


def test_mesures_frequence_vaut_f0_a_l_equilibre(m, net, plate):
    """Sans perturbation, ``dphi/dt`` est nul et la mesure vaut exactement f0.

    Le piege serait de confondre la pulsation du repere dq avec celle du
    reseau : la tension d'un bus est FIXE dans le repere a l'equilibre, donc sa
    derivee d'angle y est nulle, et toute la frequence vient de ``w_repere``.
    Une mesure qui rendrait 0 Hz ou 100 Hz signalerait ce contresens.
    """
    mes = mesures(m, plate)
    freqs = [n for n in mes if n.startswith("f_bus")]
    assert freqs, "aucun bus capacitif : le test ne verifie rien"
    for nom in freqs:
        assert mes[nom] == pytest.approx(net.f0, abs=1e-6)
        assert mes["fm" + nom[1:]] == pytest.approx(net.f0, abs=1e-6)


def test_mesures_suivent_le_transitoire(m, perturbee):
    """Apres un echelon de puissance, les mesures doivent BOUGER.

    Un test qui ne verifierait que le point de fonctionnement passerait aussi
    avec des mesures figees ; celui-ci verifie qu'elles vivent, et dans le bon
    sens : plus de puissance active demandee, plus de puissance injectee, et
    une frequence qui s'ecarte de 50 Hz le temps du transitoire.
    """
    tr, uext_fn = perturbee
    mes = mesures(m, tr, uext_fn, quoi=["P_SM1", "V_bus2", "f_bus2", "fm_bus2"])
    assert mes["P_SM1"][-1] > mes["P_SM1"][0] + 0.01
    assert np.abs(mes["f_bus2"] - 50.0).max() > 1e-4
    # un filtre passe-bas ne sort pas de la plage de son entree
    f, fm = mes["f_bus2"], mes["fm_bus2"]
    assert fm[0] == pytest.approx(f[0])
    assert f.min() - 1e-12 <= fm.min() and fm.max() <= f.max() + 1e-12
    assert np.abs(mes["V_bus2"] - mes["V_bus2"][0]).max() > 1e-5


def test_mesures_quoi_restreint_et_valide(m, plate):
    mes = mesures(m, plate, quoi=["V_bus2", "P_SM1"])
    assert list(mes) == ["V_bus2", "P_SM1"]
    with pytest.raises(KeyError, match="mesure"):
        mesures(m, plate, quoi=["V_bus2", "Z_bus42"])


def test_noeud_resistif_n_offre_pas_de_frequence(m, net, plate):
    """Un noeud sans capacite n'a pas d'etat de tension, donc pas de dv/dt.

    Le catalogue de mesures doit le dire en n'exposant pas la grandeur, plutot
    qu'en rendant un zero ou un NaN que personne ne remarquerait.
    """
    from g2egrid.elements import ResistiveBus

    resistifs = [b for b, el in net.buses.items() if isinstance(el, ResistiveBus)]
    assert resistifs, "ce reseau n'a aucun noeud resistif : le test ne verifie rien"
    mes = mesures(m, plate)
    for b in resistifs:
        assert f"V_bus{b}" in mes and f"angle_bus{b}" in mes
        assert f"f_bus{b}" not in mes and f"fm_bus{b}" not in mes


def test_filtre_premier_ordre_exact_a_pas_irregulier():
    """Le filtre doit rendre la reponse CONTINUE, quel que soit le maillage.

    Sur un echelon, la reponse exacte du premier ordre est
    ``1 - exp(-(t - t0)/tau)`` ; avec ``a = 1 - exp(-dt/tau)`` la recurrence la
    redonne terme a terme, sur un pas quelconque. Avec le ``a = dt/tau`` d'un
    filtre numerique ordinaire, elle en differerait d'autant plus que le pas
    est grand -- c'est-a-dire hors des transitoires, la ou l'integrateur se
    relache.
    """
    rng = np.random.default_rng(3)
    t = np.cumsum(rng.uniform(1e-4, 2e-2, 60))
    t = np.concatenate([[0.0], t])
    v = np.ones(len(t))
    v[0] = 0.0  # echelon applique juste apres t[0]
    tau = 1.0 / 50.0
    y = _premier_ordre(v, t, tau)
    exact = 1.0 - np.exp(-(t - t[0]) / tau)
    assert y[0] == 0.0
    assert np.abs(y[1:] - exact[1:]).max() < 1e-14

    naif = np.zeros(len(t))
    for k in range(1, len(t)):
        naif[k] = naif[k - 1] + (t[k] - t[k - 1]) / tau * (v[k] - naif[k - 1])
    assert np.abs(naif - exact).max() > 1e-3, "le filtre naif devrait, lui, deriver"


# ======================================================================
# D. linearisation le long d'une trajectoire
# ======================================================================
def test_linearisation_a_l_origine_redonne_la_matrice_d_etat(m, ss, plate):
    """A ``t = 0`` sur une trajectoire non perturbee, ``A`` est celle du modele.

    Les deux matrices sortent de chemins differents -- complement de Schur sur
    les jacobiens assembles d'un cote, theoreme des fonctions implicites sur le
    systeme couple de l'autre -- et doivent coincider a l'arrondi pres.
    """
    A = linearisation(m, plate, 0.0)
    assert A.shape == ss.A.shape
    assert np.abs(A - ss.A).max() / np.abs(ss.A).max() < 1e-10


def test_linearisation_prend_l_instant_rendu_le_plus_proche(m, perturbee):
    """Pas d'interpolation : un ``t`` entre deux points retombe sur le voisin."""
    tr, uext_fn = perturbee
    exact = linearisation(m, tr, 0.20, uext_fn)
    entre = linearisation(m, tr, 0.2049, uext_fn)
    assert np.abs(exact - entre).max() == 0.0
    autre = linearisation(m, tr, 0.2051, uext_fn)
    assert np.abs(exact - autre).max() > 0.0


def test_linearisation_bouge_pendant_le_transitoire(m, ss, perturbee):
    """La matrice instantanee doit s'ecarter de celle de l'equilibre.

    Si elle ne bougeait pas, la fonction ne repondrait a rien : c'est justement
    parce que ``A`` se deplace avec l'etat qu'un mode peut n'etre dangereux que
    pendant le transitoire.
    """
    tr, uext_fn = perturbee
    ref = np.abs(ss.A).max()
    ecarts = [
        np.abs(linearisation(m, tr, t, uext_fn) - ss.A).max() / ref
        for t in (0.0, 0.1, 0.5)
    ]
    assert ecarts[0] < 1e-10
    assert ecarts[1] > 1e-5
    assert ecarts[2] > ecarts[1]


def test_linearisation_refuse_une_trajectoire_vide(m):
    vide = Trajectoire(t=np.zeros(0), x=np.zeros((0, m.c.nx)), etats=m.etats)
    with pytest.raises(ValueError, match="vide"):
        linearisation(m, vide, 0.0)


# ======================================================================
# E. robustesse
# ======================================================================
def test_diagnostic_nomme_l_equation_algebrique_fautive(m):
    """Le message d'echec doit remonter du numero de ligne au nom d'element.

    Un ``sol.message`` de MINPACK ne dit rien du reseau. On fabrique ici un
    residu dont on SAIT ou se trouve le pire terme -- le courant d'excitation
    de la machine -- et on verifie que le diagnostic le retrouve.
    """
    c = m.c
    k = [e.name for e in c.els].index("SM1")
    zu = np.concatenate([c.z0(), c.u0()])
    zu[c.oz[k] + c.cms[k].spec.algebs.index("i_fd")] += 500.0
    faux = types.SimpleNamespace(x=zu, message="not making good progress")

    msg = m._diagnostic(faux, m._x(c.x0()), m.uext0(), 0.42)
    assert "t = 0.42 s" in msg
    assert "equation algebrique 'i_fd' de SM1" in msg
    assert "equations g" in msg and "loi de cablage" in msg
    assert "not making good progress" in msg


def test_diagnostic_nomme_la_ligne_de_cablage_fautive(m):
    """Meme chose du cote de la loi de cablage : l'entree et son element."""
    c = m.c
    k = [e.name for e in c.els].index("L1_2_3")
    el = c.els[k]
    zu = np.concatenate([c.z0(), c.u0()])
    zu[c.nz + c.ou[k] + el.port_input_idx(0)[0]] += 7.0
    faux = types.SimpleNamespace(x=zu, message="bidon")

    msg = m._diagnostic(faux, m._x(c.x0()), m.uext0(), 1.5)
    assert "t = 1.5 s" in msg
    assert "loi de cablage sur l'entree 'vd_f' de L1_2_3" in msg


def test_newton_en_echec_leve_un_message_exploitable(m, monkeypatch):
    """Bout en bout : un systeme sans solution doit produire ce message.

    On rend une equation de SM1 impossible (un residu constant, que rien ne
    peut annuler). Un etat seulement "absurde" -- des flux a 1e6 -- ne
    suffit plus : le Newton creux y trouve la solution, qui existe (des
    courants de 1e8), la ou le dogleg de MINPACK echouait.
    """
    c = m.c
    k = [e.name for e in c.els].index("SM1")
    ligne = c.oz[k]
    vrai = m.residu

    def impossible(zu, xs, uext):
        r = vrai(zu, xs, uext)
        r[ligne] = 1.0
        return r

    monkeypatch.setattr(m, "residu", impossible)
    with pytest.raises(RuntimeError) as capture:
        m.resoudre_algebrique(c.x0(), m.uext0(), t=0.7)
    msg = str(capture.value)
    assert "t = 0.7 s" in msg and "|residu|" in msg
    assert "de SM1" in msg


def test_trajectoire_non_finie_est_refusee(m):
    """Un NaN dans la trajectoire doit lever, en nommant l'etat et l'instant.

    Newton converge parfaitement sur des ``NaN`` -- aucune comparaison a une
    tolerance n'echoue -- et ``solve_ivp`` peut donc rendre ``success = True``
    avec une trajectoire entierement polluee. Sans ce garde-fou, les mesures
    rendraient des tableaux de ``NaN`` sans que rien ne l'ait signale.
    """
    y = np.zeros((m.c.nx, 4))
    y[3, 2] = np.nan
    faux = types.SimpleNamespace(y=y, t=np.array([0.0, 0.1, 0.2, 0.3]))
    with pytest.raises(RuntimeError, match="divergente"):
        m._verifier_fini(faux)
    msg = str(pytest.raises(RuntimeError, m._verifier_fini, faux).value)
    assert m.c.state_names[3] in msg and "t = 0.2 s" in msg


@pytest.fixture(scope="module")
def instable():
    """SMIB a amortissement NEGATIF : son mode local diverge lentement.

    Une instabilite douce (Re(lambda) de l'ordre de +2 s^-1) est ce qu'il faut
    ici : elle diverge franchement en deux secondes tout en restant integrable
    en une seconde de calcul. Les instabilites raides d'un convertisseur -- un
    GFL a PLL trop rapide sur reseau faible part a +300 s^-1 -- font ramer
    l'integrateur pendant des minutes avant de deborder, ce qui n'est pas
    testable.
    """
    net = Network(f0=50.0, Sb=1e6)
    net.add_bus(1, 20e3, btype=2, V=1.0, Pg=0.8, C=1e-3)
    net.add_bus(2, 20e3, btype=1, V=1.0, C=1e-3)
    net.add(SynchronousMachine("SM1", 1, Sn=1e6, Sb=1e6, wb=net.wb, KD=-60.0))
    net.add(Line("L1", 1, 2, R=0.01, X=0.20, wb=net.wb))
    net.add(IB("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
    net.run_loadflow()
    net.solve_equilibrium()
    return net


def test_reseau_instable_diverge_proprement(instable):
    """Une divergence doit se VOIR : elle croit, elle reste finie, ou elle leve.

    Le piege serait une trajectoire de ``NaN`` rendue sans erreur. On verifie
    donc que tout est fini, que l'ecart au point de fonctionnement croit, et
    que son taux de croissance est bien celui que le modele lineaire annonce.
    """
    ev = instable.build_state_space().eig()
    taux = ev.real[np.abs(ev) > 1e-6].max()
    assert taux > 0.5, "ce reseau n'est pas instable : le test ne verifie rien"

    m = ModeleNonLineaire(instable)
    t_eval = np.linspace(0.0, 2.0, 41)
    tr = m.simuler(
        (0.0, 2.0), uext_fn=m.echelon("SM1.p_ref", 1e-3, t0=0.02), t_eval=t_eval
    )
    assert np.isfinite(tr.x).all(), "des NaN sont passes sans etre signales"

    ecart = np.abs(tr.x - m.c.x0()).max(axis=1)
    assert ecart[-1] > 10.0 * ecart[4], f"pas de divergence visible ({ecart[-1]:.2e})"
    pente = np.polyfit(t_eval[20:], np.log(ecart[20:]), 1)[0]
    assert pente == pytest.approx(taux, rel=0.25), (
        f"croissance mesuree {pente:.3f} s^-1 contre {taux:.3f} annonce par A"
    )
