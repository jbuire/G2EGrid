"""Tests de l'identification modale (``g2egrid.identification``).

Deux references independantes, dans cet ordre de cout :

* un signal SYNTHETIQUE dont on connait les poles a la machine pres. C'est
  lui qui dit si l'algebre est juste, et rien d'autre ne le dit : sur un vrai
  reseau, un pole faux de 1 % reste parfaitement credible. Les tolerances y
  sont donc serrees a l'os -- 1e-8 en relatif sans bruit -- et toute
  degradation ulterieure se lira ici avant de se lire ailleurs ;

* le spectre de la matrice ``A`` de ``wscc9``, qui ne partage aucune ligne de
  code de resolution avec la simulation non lineaire sur laquelle on
  identifie. C'est l'epreuve qui a un sens physique : elle dit que ``A``
  decrit vraiment ce que fait le modele non lineaire quand on le perturbe.

Les valeurs attendues sur ``wscc9`` ne sont pas des voeux : le spectre de
``A`` donne les deux modes electromecaniques a 2,2000 Hz / 8,696 % et
2,0115 Hz / 11,251 %, et l'identification sur les trois vitesses rotor les
rend a mieux que 0,1 % en frequence et 1 % en amortissement. Les tolerances
du contrat -- 2 % et 15 % -- laissent donc un facteur 20 de marge, ce qui est
voulu : le test doit se casser sur une regression du modele, pas sur un
changement de version de LAPACK.
"""

import numpy as np
import pytest

from g2egrid import modal
from g2egrid.catalogue import construire
from g2egrid.identification import (
    Identification,
    ModeIdentifie,
    apparier,
    identifier,
    reconstruire,
    reechantillonner,
    sur_trajectoire,
)
from g2egrid.modal import Mode
from g2egrid.temporel import Evenement, ModeleNonLineaire

METHODES = ("pencil", "era")

#: Les poles du signal synthetique : 1,3 Hz peu amorti et 3,7 Hz plus amorti,
#: soit l'ordre de grandeur d'un mode electromecanique et d'un mode de
#: regulation. Les deux sont assez separes pour qu'aucune methode n'ait
#: d'excuse, et c'est le but : ce test mesure l'algebre, pas la resolution.
POLES = np.array([-0.4 + 2j * np.pi * 1.3, -1.1 + 2j * np.pi * 3.7])

#: Residus (mode, voie) du signal synthetique. La seconde voie porte le
#: second mode en opposition, comme deux machines d'un mode inter-machines.
RESIDUS = np.array([[1.0 + 0.0j, 0.5 - 0.2j], [0.3 + 0.4j, -0.8 + 0.1j]])

T_ECHELON = 0.1  # instant de l'echelon de consigne sur wscc9 [s]


# ----------------------------------------------------------------------
# fixtures
# ----------------------------------------------------------------------
def _signal(t: np.ndarray, poles=POLES, residus=RESIDUS) -> np.ndarray:
    """Somme d'exponentielles amorties REELLE : ``2.Re(somme r_i.exp(lam_i.t))``.

    Ecrire le signal comme la partie reelle d'une somme sur les seuls poles
    de partie imaginaire POSITIVE, et non comme une somme sur les paires
    conjuguees, evite l'erreur classique du facteur deux : l'identification
    rend les deux membres de chaque paire, chacun avec le residu ``r_i``, et
    c'est bien ``2.|r_i|`` qui est l'amplitude physique de l'oscillation.
    """
    return 2.0 * np.real(np.exp(np.outer(t, poles)) @ residus)


@pytest.fixture(scope="module")
def synthetique():
    """6 s a 100 Hz : 601 points, une dizaine de periodes du mode lent."""
    t = np.arange(0.0, 6.0 + 1e-9, 0.01)
    return t, _signal(t)


@pytest.fixture(scope="module")
def wscc9():
    """Le reseau a trois machines et son modele d'etat lineaire."""
    net = construire("wscc9")
    return net, net.build_state_space()


@pytest.fixture(scope="module")
def trajectoire(wscc9):
    """Reponse NON LINEAIRE a un echelon de 0,01 pu sur ``SM1.p_ref``.

    L'amplitude est petite devant la puissance de la machine : on cherche la
    reponse du systeme LINEARISE, et un echelon de 0,1 pu ferait sortir les
    ecarts du domaine ou la comparaison au spectre de ``A`` a un sens.

    QUATRE SECONDES -- le moins amorti des deux modes (8,7 % a 2,2 Hz) met
    environ 1,6 s a perdre 95 % de son amplitude, mais c'est le NOMBRE de
    periodes vues qui fait la precision de l'amortissement identifie. Mesure
    sur ce cas, a rang impose 16, l'erreur d'amortissement en fonction de la
    longueur de fenetre vaut 56 % sur 0,6 s, 37 % sur 1 s, 37 % sur 1,5 s,
    7,7 % sur 2 s, 0,8 % sur 3 s et 0,1 % sur 4 s. La tolerance de 15 % du
    contrat tombe donc en dessous de deux secondes, et deux secondes la
    tiendraient de justesse ; quatre la tiennent avec deux ordres de
    grandeur de marge, pour cinq secondes de calcul.

    Le residu de reconstruction ne previent PAS de cette erreur-la : sur
    0,6 s il vaut 0,0002 %, soit dix fois mieux que sur 4 s (cf.
    ``Identification``). C'est la seule facon de se tromper que
    l'identification ne sache pas signaler elle-meme.

    ``t_eval`` uniforme : c'est la condition d'entree de :func:`identifier`,
    et la demander a ``simuler`` ne coute rien (Radau interpole sur ses
    propres pas, il n'en rajoute pas).
    """
    net, _ = wscc9
    m = ModeleNonLineaire(net)
    dt, duree = 0.005, 4.0
    t_eval = np.arange(0.0, duree + 0.5 * dt, dt)
    uext_fn = m.scenario(Evenement(T_ECHELON, "SM1.p_ref", 0.01))
    tr = m.simuler((0.0, duree), uext_fn=uext_fn, t_eval=t_eval)
    # On ne donne a identifier que ce qui SUIT l'echelon : avant, le signal
    # est plat, et une fenetre qui commencerait avant l'evenement demanderait
    # au modele en somme d'exponentielles de representer aussi le saut.
    apres = tr.t >= T_ECHELON - 1e-12
    noms = [f"SM{i}.Dwr" for i in (1, 2, 3)]
    y = np.column_stack([tr.etat(n)[apres] for n in noms])
    return tr.t[apres], y, noms


@pytest.fixture(scope="module")
def electromecaniques(wscc9):
    """Les deux modes electromecaniques predits, du plus "rotor" au moins.

    Meme tri que ``tests/test_modal.py`` : la participation cumulee des
    vitesses rotor est ce qui separe un mode electromecanique d'un mode
    d'excitation tombe dans la meme bande de frequence -- et il y en a un sur
    ``wscc9``, a 1,84 Hz.
    """
    _, ss = wscc9
    lam = ss.eig()
    f = modal.frequence(lam)
    zeta = modal.amortissement(lam)
    p = modal.participations(ss)
    iv = [i for i, n in enumerate(ss.states) if n.endswith(".Dwr")]
    garde = np.flatnonzero((lam.imag > 1e-3) & (f > 0.5) & (f < 5.0))
    ordre = sorted((int(k) for k in garde), key=lambda k: -p[k, iv].sum())
    return [
        Mode(
            indice=k,
            valeur=complex(lam[k]),
            f=float(f[k]),
            zeta=float(zeta[k]),
            etats=[],
        )
        for k in ordre[:2]
    ]


@pytest.fixture(scope="module")
def predits(wscc9):
    """Les modes oscillants lents du modele lineaire, candidats a l'appariement."""
    _, ss = wscc9
    return modal.modes(ss, f_max=6.0, zeta_max=98.0, oscillants=True)


def _poles_positifs(ident: Identification) -> np.ndarray:
    """Les valeurs propres identifiees de partie imaginaire > 0, par frequence."""
    z = np.array([m.valeur for m in ident.modes if m.valeur.imag > 1e-9])
    return z[np.argsort(z.imag)]


# ======================================================================
# A. l'algebre, sur un signal dont on connait la reponse
# ======================================================================
@pytest.mark.parametrize("methode", METHODES)
def test_poles_connus_retrouves_a_la_machine_pres(synthetique, methode):
    """Sans bruit, les deux estimateurs doivent etre EXACTS.

    Une somme finie d'exponentielles echantillonnee uniformement est
    exactement de rang fini : la matrice de Hankel a exactement 4 valeurs
    singulieres non nulles, et les poles sont exactement les valeurs propres
    du decalage. Il n'y a donc aucune erreur de modele a tolerer, seulement
    l'arrondi flottant. Exiger 1e-8 quand la mesure donne 5e-15 laisse six
    decades de marge a l'arithmetique, et pas une a une erreur de convention.

    ``offset=False`` : le signal tend vers zero de lui-meme, et retrancher la
    moyenne de sa queue -- qui vaut encore 9 % de l'amplitude initiale a 6 s
    -- introduirait justement le biais que ce test doit exclure.
    """
    t, y = synthetique
    ident = identifier(t, y, methode=methode, offset=False)
    assert ident.rang == 4  # deux paires conjuguees, rien de plus
    err = np.abs(_poles_positifs(ident) - POLES) / np.abs(POLES)
    assert err.max() < 1e-8, f"{methode} : erreur relative {err}"


@pytest.mark.parametrize("methode", METHODES)
def test_poles_retrouves_sous_bruit(synthetique, methode):
    """Avec 0,1 % de bruit blanc, l'erreur doit rester de l'ordre du bruit.

    La troncature au rang 4 est ce qui debruite : les 4 premieres valeurs
    singulieres portent le signal, les 597 autres le bruit. On attend donc
    une erreur bien INFERIEURE au bruit d'entree -- la mesure donne 1,5e-5
    pour 1e-3 de bruit, soit un facteur 60 -- et on tolere 1e-3, c'est-a-dire
    "pas pire que le bruit". Un estimateur qui se contenterait de suivre le
    bruit passerait tout juste, et c'est voulu : ce test ne certifie pas une
    performance, il interdit une divergence.
    """
    t, y = synthetique
    rng = np.random.default_rng(0)
    bruit = 1e-3 * np.abs(y).max() * rng.standard_normal(y.shape)
    ident = identifier(t, y + bruit, methode=methode, rang=4, offset=False)
    err = np.abs(_poles_positifs(ident) - POLES) / np.abs(POLES)
    assert err.max() < 1e-3, f"{methode} : erreur relative {err}"


@pytest.mark.parametrize("methode", METHODES)
def test_offset_retranche_ne_deplace_pas_les_poles(synthetique, methode):
    """Un signal decale d'une constante doit rendre les MEMES poles.

    C'est la situation d'apres echelon : le systeme s'arrete sur un autre
    equilibre, et le signal garde une composante continue qu'aucune
    exponentielle amortie ne sait representer. Si le retrait de l'offset
    fonctionne, la constante n'a aucun effet sur les poles ; s'il est mal
    fait, elle les deplace tous.
    """
    t, y = synthetique
    ident = identifier(t, y + 0.37, methode=methode)
    err = np.abs(_poles_positifs(ident)[:2] - POLES) / np.abs(POLES)
    assert err.max() < 1e-6, f"{methode} : erreur relative {err}"


@pytest.mark.parametrize("methode", METHODES)
def test_reconstruction_rend_le_signal(synthetique, methode):
    """Les modes identifies doivent RESYNTHETISER le signal, pas lui ressembler.

    Le residu de reconstruction est le seul garde-fou contre une
    identification qui rend des poles plausibles sans rapport avec le signal.
    On verifie les deux bouts : que le residu annonce est minuscule, et qu'il
    dit bien ce qu'il pretend -- la resynthese explicite par
    :func:`reconstruire` doit redonner le signal fourni.
    """
    t, y = synthetique
    ident = identifier(t, y, methode=methode, offset=False)
    assert ident.residu < 1e-10
    rec = reconstruire(ident.modes, ident.t - ident.t0)
    assert rec.shape == ident.y.shape
    assert np.abs(rec - ident.y).max() < 1e-10 * np.abs(y).max()
    assert np.abs(ident.reconstruction() - rec).max() == 0.0


def test_residu_denonce_un_rang_trop_bas(synthetique):
    """Un rang insuffisant rend des modes credibles qui n'expliquent rien.

    C'est le piege que le residu existe pour signaler : a rang 2, la methode
    rend une paire conjuguee de frequence et d'amortissement parfaitement
    presentables -- rien dans les valeurs propres ne trahit quoi que ce soit
    -- alors qu'elle ignore la moitie du signal. Seul le residu le dit.
    """
    t, y = synthetique
    pauvre = identifier(t, y, rang=2, offset=False)
    riche = identifier(t, y, rang=4, offset=False)
    (pole,) = _poles_positifs(pauvre)
    assert 1.0 < modal.frequence(pole)[0] < 4.0  # credible, et pourtant faux
    assert pauvre.residu > 0.1
    assert riche.residu < 1e-10


def test_multi_voies_voit_un_mode_absent_d_une_voie():
    """Un mode dont une voie est un noeud n'est identifiable que sur les autres.

    C'est l'argument central du multi-voies, et il n'a rien de theorique :
    sur ``wscc9``, la vitesse de SM1 est presque un noeud du mode a 2,2 Hz
    (cf. :func:`test_une_seule_voie_ne_suffit_pas`). On le reproduit ici a
    l'etat pur -- residu exactement nul sur la premiere voie -- pour separer
    la cause de ses consequences : sur une voie, le mode n'existe pas ; sur
    les deux, il ressort exactement.
    """
    t = np.arange(0.0, 6.0 + 1e-9, 0.01)
    residus = np.array([[1.0 + 0.0j, 0.4 + 0.1j], [0.0 + 0.0j, 0.9 - 0.3j]])
    y = _signal(t, POLES, residus)

    seule = identifier(t, y[:, [0]], offset=False)
    assert seule.rang == 2  # la voie 0 ne porte qu'un seul mode
    assert np.abs(_poles_positifs(seule) - POLES[0]).max() < 1e-8

    deux = identifier(t, y, offset=False)
    assert deux.rang == 4
    err = np.abs(_poles_positifs(deux) - POLES) / np.abs(POLES)
    assert err.max() < 1e-8


def test_energies_trient_et_totalisent_cent(synthetique):
    """L'energie est une part : elle se somme a 100 % et classe les modes."""
    t, y = synthetique
    ident = identifier(t, y, rang=6, offset=False)
    energies = [m.energie for m in ident.modes]
    assert sum(energies) == pytest.approx(100.0)
    assert energies == sorted(energies, reverse=True)
    # les modes de troncature ne portent rien : c'est ce qui les designe
    assert energies[4] < 1e-6 * energies[0]


def test_conjugues_rendus_par_paire(synthetique):
    """L'identification rend les DEUX membres de chaque paire conjuguee.

    C'est la convention du module, et elle est verifiable : la matrice de
    Hankel d'un signal reel est reelle, donc ses valeurs propres sont
    conjuguees au bit pres. ``Identification.oscillants`` est ce qui ramene
    la liste a la convention de ``modal.modes(..., oscillants=True)``.
    """
    t, y = synthetique
    ident = identifier(t, y, offset=False)
    hauts = ident.oscillants()
    assert len(hauts) == 2
    for m in hauts:
        assert any(n.valeur == m.valeur.conjugate() for n in ident.modes), (
            f"{m.valeur} sans conjugue exact"
        )


@pytest.mark.parametrize("methode", METHODES)
def test_amplitude_et_phase_referees_au_debut_de_la_fenetre(synthetique, methode):
    """Les residus comptent a partir de ``t0``, et c'est verifiable au chiffre pres.

    La convention n'est pas gratuite (cf. :func:`reconstruire`), et elle se
    voit : la fenetre par defaut demarre 10 ms apres l'origine du signal, et
    le residu rendu vaut alors exactement ``r . exp(lam.t0)`` -- ici
    ``1 + 0j`` devient ``0,9927 + 0,0813j``, soit 4,7 degres de phase de plus
    a 1,3 Hz. Un test qui comparerait aux residus injectes SANS ce facteur
    echouerait de 8 %, et la tentation serait de relacher la tolerance : ce
    serait effacer la convention au lieu de la verifier.
    """
    t, y = synthetique
    ident = identifier(t, y, methode=methode, offset=False)
    attendus = RESIDUS * np.exp(POLES * (ident.t0 - t[0]))[:, None]
    for k, m in enumerate(sorted(ident.oscillants(), key=lambda m: m.f)):
        assert np.abs(m.residu - attendus[k]).max() < 1e-8 * np.abs(RESIDUS).max()
        assert m.voies == ["y0", "y1"]
    assert ident.modes[0].voie_dominante in ("y0", "y1")


# ----------------------------------------------------------------------
# A bis. les conditions d'entree
# ----------------------------------------------------------------------
def test_pas_non_uniforme_refuse_avec_le_remede(synthetique):
    """Un pas irregulier doit lever, en disant quoi faire.

    Silencieusement degrade serait le pire des cas : les deux methodes
    rendraient des poles d'allure normale, faux d'un facteur qui depend de la
    facon dont l'integrateur a choisi ses pas -- donc irreproductible.
    """
    t, y = synthetique
    tordu = t.copy()
    tordu[100:] += 0.003
    with pytest.raises(ValueError, match="pas non uniforme"):
        identifier(tordu, y, offset=False)


def test_reechantillonner_rattrape_une_grille_irreguliere(synthetique):
    """Sur une grille tiree au hasard, la spline cubique suffit a identifier.

    Le cas d'usage est l'enregistrement dont on ne maitrise pas la grille, ou
    la trajectoire Radau deja simulee sans ``t_eval``. Le prix a payer pour
    ne pas avoir demande les bons instants se lit dans le chiffre : 3e-6
    d'erreur relative sur les poles, contre 5e-15 sur grille propre, soit
    neuf decades perdues -- mais 3e-6 reste negligeable devant tout ce qui
    limite une identification reelle. La tolerance est mise a 1e-4 pour
    laisser passer une autre version de ``CubicSpline``, pas pour couvrir
    une regression.
    """
    t, _ = synthetique
    rng = np.random.default_rng(1)
    irregulier = np.sort(rng.uniform(t[0], t[-1], 900))
    tu, yu = reechantillonner(irregulier, _signal(irregulier), dt=0.01)
    assert np.abs(np.diff(tu) - 0.01).max() < 1e-12
    ident = identifier(tu, yu, rang=4, offset=False)
    err = np.abs(_poles_positifs(ident) - POLES) / np.abs(POLES)
    assert err.max() < 1e-4, f"erreur relative {err}"


def test_fenetre_par_defaut_saute_le_retard(synthetique):
    """Par defaut la fenetre commence ``RETARD_DEFAUT`` apres le premier point."""
    t, y = synthetique
    ident = identifier(t, y, offset=False)
    assert ident.t0 == pytest.approx(t[0] + 0.01)  # premier point apres 5 ms
    tardif = identifier(t, y, debut=1.0, fin=4.0, offset=False)
    assert tardif.t0 == pytest.approx(1.0)
    assert tardif.t[-1] == pytest.approx(4.0)
    assert np.abs(_poles_positifs(tardif) - POLES).max() < 1e-6


def test_arguments_refuses(synthetique):
    """Ce qui n'a pas de sens doit lever, pas produire un resultat."""
    t, y = synthetique
    with pytest.raises(ValueError, match="methode"):
        identifier(t, y, methode="prony")
    with pytest.raises(ValueError, match="noms de voies"):
        identifier(t, y, voies=["a", "b", "c"])
    with pytest.raises(ValueError, match="fenetre"):
        identifier(t, y, debut=5.999)
    with pytest.raises(ValueError, match="rang"):
        identifier(t, y, rang=10_000)
    with pytest.raises(ValueError, match="n OU dt"):
        reechantillonner(t, y, n=100, dt=0.01)


# ======================================================================
# B. appariement
# ======================================================================
def _mode(valeur: complex, indice: int = 0) -> Mode:
    """Un ``Mode`` lineaire construit a la main."""
    return Mode(
        indice=indice,
        valeur=valeur,
        f=float(modal.frequence(valeur)[0]),
        zeta=float(modal.amortissement(valeur)[0]),
        etats=[],
    )


def _identifie(valeur: complex, energie: float = 50.0) -> ModeIdentifie:
    """Un ``ModeIdentifie`` construit a la main, sans passer par un signal."""
    return ModeIdentifie(
        valeur=valeur,
        f=float(modal.frequence(valeur)[0]),
        zeta=float(modal.amortissement(valeur)[0]),
        amplitude=np.array([1.0]),
        phase=np.array([0.0]),
        energie=energie,
        voies=["y0"],
    )


def _glouton(d: np.ndarray) -> list[tuple[int, int]]:
    """Appariement glouton : la plus petite distance d'abord, et on verrouille.

    Cette fonction n'est PAS dans le module : elle n'y sert a rien. Elle est
    ici parce qu'un test qui affirme "le global fait mieux que le glouton"
    sans jamais faire tourner le glouton n'affirme rien du tout.
    """
    d = d.copy()
    paires = []
    for _ in range(min(d.shape)):
        i, j = np.unravel_index(np.argmin(d), d.shape)
        paires.append((int(i), int(j)))
        d[i, :] = np.inf
        d[:, j] = np.inf
    return paires


def test_appariement_global_la_ou_le_glouton_se_trompe():
    """Deux modes voisins : le glouton verrouille la mauvaise paire.

    Le modele lineaire predit 1,59 Hz et 1,91 Hz ; l'identification rend
    1,46 Hz et 1,72 Hz. La paire la plus serree du tableau est
    (1,72 identifie <-> 1,59 predit, a 0,0796) : le glouton la prend
    d'abord, et il ne lui reste plus qu'a marier 1,46 a 1,91, a 0,2367.
    L'appariement global, lui, lit la configuration pour ce qu'elle est --
    les deux modes identifies sont decales VERS LE BAS -- et apparie 1,46 a
    1,59 et 1,72 a 1,91 : somme 0,184 contre 0,316, soit 1,7 fois moins.

    Ce n'est pas une subtilite d'algorithmique : les deux lectures donnent
    des diagnostics opposes. Le glouton conclut que le mode a 1,91 Hz est
    identifie 24 % trop bas -- donc qu'il y a une erreur quelque part dans
    le modele ; le global conclut que les deux modes sont identifies 8,5 et
    10 % trop bas -- donc que c'est l'identification qui derive, ce qui
    n'appelle pas du tout la meme suite.
    """
    lin = [_mode(-1 + 10.0j, 0), _mode(-1 + 12.0j, 1)]
    ident = [_identifie(-1 + 9.15j), _identifie(-1 + 10.8j)]

    vi = np.array([m.valeur for m in ident])
    vl = np.array([m.valeur for m in lin])
    d = np.abs(vi[:, None] - vl[None, :]) / np.abs(vl)[None, :]
    assert _glouton(d) == [(1, 0), (0, 1)], "le glouton doit se tromper ici"
    somme_gloutonne = d[1, 0] + d[0, 1]

    ap = apparier(ident, lin, seuil=0.3)
    assert len(ap.paires) == 2
    obtenu = {(p.identifie.valeur, p.lineaire.valeur) for p in ap.paires}
    assert obtenu == {
        (ident[0].valeur, lin[0].valeur),
        (ident[1].valeur, lin[1].valeur),
    }
    assert sum(p.distance for p in ap.paires) < somme_gloutonne


def test_apparier_indifferent_a_la_convention_des_conjugues():
    """Presents ou non, les conjugues ne doivent rien changer au resultat.

    ``modal.modes(..., oscillants=True)`` ne rend que la moitie superieure du
    spectre, l'identification rend les deux moities, et rien ne garantit que
    l'appelant y pense. Le repliement est donc fait des deux cotes, et ce
    test le verifie dans les quatre combinaisons : le meme appariement doit
    sortir, avec les memes ecarts.
    """
    lin = [_mode(-1 + 10.0j, 0), _mode(-2 + 20.0j, 1)]
    ident = [_identifie(-1.05 + 10.1j), _identifie(-2.1 + 20.2j)]
    lin_double = lin + [_mode(m.valeur.conjugate(), 9) for m in lin]
    ident_double = ident + [_identifie(m.valeur.conjugate()) for m in ident]

    reference = None
    for i in (ident, ident_double):
        for lo in (lin, lin_double):
            ap = apparier(i, lo)
            obtenu = sorted(
                (p.lineaire.f, round(p.ecart_f, 9), round(p.ecart_zeta, 9))
                for p in ap.paires
            )
            assert len(ap.paires) == 2
            assert not ap.identifies_seuls and not ap.lineaires_seuls
            if reference is None:
                reference = obtenu
            assert obtenu == reference


def test_apparier_refuse_au_dela_du_seuil_et_rend_les_orphelins():
    """Un mode predit que rien n'approche doit ressortir NON APPARIE.

    Un appariement qui marie tout le monde ne dit rien : c'est la liste des
    laisses-pour-compte qui porte le diagnostic -- mode inobservable sur les
    voies choisies, mode inexcite, ou mode qui n'existe pas dans le systeme.
    """
    lin = [_mode(-1 + 10.0j, 0), _mode(-3 + 60.0j, 1)]
    ident = [_identifie(-1.02 + 10.05j), _identifie(-0.5 + 2.0j, energie=1.0)]
    ap = apparier(ident, lin, seuil=0.1)
    assert len(ap.paires) == 1
    assert ap.paires[0].lineaire.valeur == lin[0].valeur
    assert [m.valeur for m in ap.lineaires_seuls] == [lin[1].valeur]
    assert [m.valeur for m in ap.identifies_seuls] == [ident[1].valeur]
    assert "non identifie" in ap.tableau()
    assert "non predit" in ap.tableau()


def test_energie_min_ecarte_les_modes_de_troncature():
    """Un mode parasite bien place ne doit pas voler la place du vrai.

    A energie negligeable, un mode de troncature tombe parfois plus pres d'un
    mode predit que le mode physique lui-meme. ``energie_min`` est ce qui
    l'ecarte AVANT l'appariement, et non apres : une fois entre dans le
    probleme d'affectation, il deplace tout le reste.
    """
    lin = [_mode(-1 + 10.0j, 0)]
    vrai = _identifie(-1.05 + 10.2j, energie=30.0)
    parasite = _identifie(-1.0 + 10.02j, energie=1e-4)
    ap = apparier([parasite, vrai], lin, energie_min=0.01)
    assert len(ap.paires) == 1
    assert ap.paires[0].identifie.valeur == vrai.valeur


def test_ecarts_f_et_zeta_dans_les_bonnes_unites():
    """``ecart_zeta`` est en POINTS, ``ecart_zeta_rel`` en %, et ils different.

    Confondre les deux est l'erreur qui fait passer un test qui devrait
    echouer : sur un mode a 8,7 % d'amortissement, un point vaut 11,5 %
    relatifs. Le contrat de validation parle de 15 % RELATIFS.
    """
    lin = _mode(-1.0 + 10.0j, 0)  # zeta = 9,95 %
    ident = _identifie(-1.1 + 10.0j)  # zeta = 10,93 %
    (p,) = apparier([ident], [lin]).paires
    assert p.ecart_zeta == pytest.approx(ident.zeta - lin.zeta)
    assert p.ecart_zeta == pytest.approx(0.98, abs=0.02)  # points
    assert p.ecart_zeta_rel == pytest.approx(100.0 * p.ecart_zeta / lin.zeta)
    assert p.ecart_zeta_rel == pytest.approx(9.85, abs=0.1)  # %
    assert p.ecart_f == pytest.approx(0.0, abs=1e-12)  # meme partie imaginaire


# ======================================================================
# C. le test decisif : wscc9, identification contre spectre de A
# ======================================================================
@pytest.mark.parametrize("methode", METHODES)
def test_modes_electromecaniques_identifies_sur_la_trajectoire(
    trajectoire, predits, electromecaniques, methode
):
    """LE test du lot : les modes de ``A`` se lisent sur la trajectoire EMT.

    On simule le modele NON LINEAIRE, on identifie sur les trois vitesses
    rotor, et on apparie au spectre de la matrice ``A`` -- qui n'a pas servi
    a la simulation. Les deux modes electromecaniques doivent se retrouver a
    moins de 2 % en frequence et 15 % en amortissement.

    Ce que ce test attrape et qu'aucun autre n'attrape : les deux premieres
    epreuves du contrat (point fixe, et ecart en carre de l'amplitude) disent
    que le modele lineaire est coherent avec LUI-MEME au voisinage de son
    point de linearisation. Celle-ci est la seule ou la reference -- le
    spectre de ``A`` -- et la mesure -- la trajectoire de Radau -- ne
    partagent aucun code de resolution.
    """
    t, y, noms = trajectoire
    ident = identifier(t, y, methode=methode, voies=noms)
    ap = apparier(ident.modes, predits, energie_min=0.01)
    par_f = {round(p.lineaire.f, 6): p for p in ap.paires}

    for attendu in electromecaniques:
        p = par_f.get(round(attendu.f, 6))
        assert p is not None, (
            f"{methode} : mode a {attendu.f:.4f} Hz non apparie.\n{ap.tableau()}"
        )
        assert abs(p.ecart_f) < 2.0, (
            f"{methode} : {attendu.f:.4f} Hz identifie a {p.identifie.f:.4f} Hz "
            f"({p.ecart_f:+.2f} %)\n{ap.tableau()}"
        )
        assert abs(p.ecart_zeta_rel) < 15.0, (
            f"{methode} : zeta de {attendu.zeta:.3f} % identifie a "
            f"{p.identifie.zeta:.3f} % ({p.ecart_zeta_rel:+.1f} %)\n{ap.tableau()}"
        )


@pytest.mark.parametrize("methode", METHODES)
def test_la_reconstruction_explique_la_trajectoire(trajectoire, methode):
    """Les modes identifies doivent resynthetiser la trajectoire simulee.

    Sans ce controle, le test precedent se contenterait de modes bien places.
    Un residu de 1e-4 dit que les 14 modes retenus rendent compte de ce que
    93 etats ont fait : c'est ce qui autorise a lire le tableau des ecarts
    comme une comparaison, et non comme une coincidence.
    """
    t, y, noms = trajectoire
    ident = identifier(t, y, methode=methode, voies=noms)
    assert ident.residu < 1e-3, f"{methode} : residu {100 * ident.residu:.3f} %"
    rec = reconstruire(ident.modes, ident.t - ident.t0)
    erreur = np.abs(rec - ident.y).max() / np.abs(ident.y).max()
    assert erreur < 1e-2


def test_robustesse_au_rang(trajectoire, predits, electromecaniques):
    """Passe un rang suffisant, les modes dominants ne doivent plus bouger.

    C'est la question qu'on pose toujours a une identification : "et si
    j'avais pris un autre rang ?". La reponse mesuree ici est nette -- de 14
    a 30, les deux modes electromecaniques bougent de moins de 0,2 % en
    frequence -- et elle a une raison : les modes surnumeraires n'ont pas
    d'energie a prendre, ils se partagent le bruit de l'integration. Le rang
    14 choisi automatiquement n'est pas un reglage heureux, c'est un
    plancher.

    Le rang 8 est la pour montrer le contraire, et le montrer FRANCHEMENT :
    a ce rang, le mode a 2,2 Hz n'est pas mal identifie, il est absent -- et
    le residu de reconstruction, lui, le signale (2,6 % contre 0,01 %). Un
    rang trop bas ne degrade pas, il supprime.
    """
    t, y, noms = trajectoire
    cibles = [round(m.f, 6) for m in electromecaniques]

    def paires_et_residu(rang):
        ident = identifier(t, y, voies=noms, rang=rang)
        ap = apparier(ident.modes, predits, energie_min=0.01)
        par_f = {round(p.lineaire.f, 6): p for p in ap.paires}
        return {c: par_f.get(c) for c in cibles}, ident.residu

    pauvre, residu_pauvre = paires_et_residu(8)
    assert any(p is None for p in pauvre.values())
    assert residu_pauvre > 1e-2  # et le residu le dit

    reference = None
    for rang in (14, 16, 20, 30):
        obtenu, residu = paires_et_residu(rang)
        assert residu < 1e-3
        valeurs = {}
        for c, p in obtenu.items():
            assert p is not None, f"rang {rang} : mode a {c} Hz perdu"
            assert abs(p.ecart_f) < 2.0
            assert abs(p.ecart_zeta_rel) < 15.0
            valeurs[c] = (p.identifie.f, p.identifie.zeta)
        if reference is None:
            reference = valeurs
        for c in cibles:
            f_ref, z_ref = reference[c]
            f_k, z_k = valeurs[c]
            assert abs(f_k - f_ref) / f_ref < 5e-3, f"rang {rang} : f a {c} Hz"
            assert abs(z_k - z_ref) / z_ref < 5e-2, f"rang {rang} : zeta a {c} Hz"


def test_une_seule_voie_ne_suffit_pas(trajectoire, predits, electromecaniques):
    """Le multi-voies n'est pas un confort : sur une voie, le critere tombe.

    La vitesse de SM1 est presque un noeud du mode a 2,2 Hz -- son residu y
    vaut 2,0e-7 contre 1,5e-6 sur SM3, sept fois moins. Identifier sur elle
    seule rend ce mode a -7,7 % en frequence et +100 % en amortissement
    (8,7 % predit, 17,4 % identifie) : le contrat de 2 % et 15 % n'est pas
    tenu, et il ne le serait pas davantage en simulant plus longtemps. Les
    trois voies ensemble le tiennent avec un facteur 20 de marge.

    C'est la justification du Hankel empile par blocs, et elle est physique
    avant d'etre numerique : un mode inter-machines n'est pas observable de
    la meme facon partout.
    """
    t, y, noms = trajectoire
    cibles = [round(m.f, 6) for m in electromecaniques]

    def pires_ecarts(colonnes, noms_colonnes):
        ident = identifier(t, y[:, colonnes], voies=noms_colonnes)
        ap = apparier(ident.modes, predits, energie_min=0.01)
        par_f = {round(p.lineaire.f, 6): p for p in ap.paires}
        if any(par_f.get(c) is None for c in cibles):
            return np.inf, np.inf  # mode perdu : pire que n'importe quel ecart
        return (
            max(abs(par_f[c].ecart_f) for c in cibles),
            max(abs(par_f[c].ecart_zeta_rel) for c in cibles),
        )

    seules = [pires_ecarts([i], [noms[i]]) for i in range(3)]
    ensemble_f, ensemble_z = pires_ecarts([0, 1, 2], noms)

    assert ensemble_f < 2.0 and ensemble_z < 15.0
    assert max(f for f, _ in seules) > 2.0, (
        f"aucune voie seule ne rate le critere : {seules}"
    )
    assert ensemble_f < min(f for f, _ in seules)


def test_pencil_et_era_sont_d_accord(trajectoire):
    """Deux algebres differentes sur le meme signal doivent converger.

    C'est le controle croise, et c'est lui qui distingue un mode du signal
    d'un artefact de troncature. Avec chacune sa forme de Hankel optimale
    (cf. ``_forme_hankel``), les deux methodes rendent ici les memes modes a
    mieux que 1 % : ce qu'elles disent est dans le signal.
    """
    t, y, noms = trajectoire
    a = identifier(t, y, methode="pencil", voies=noms)
    b = identifier(t, y, methode="era", voies=noms)
    assert a.rang == b.rang
    fa = np.array(sorted(m.f for m in a.oscillants(energie_min=0.05)))
    fb = np.array(sorted(m.f for m in b.oscillants(energie_min=0.05)))
    assert len(fa) == len(fb), f"{fa} vs {fb}"
    assert np.max(np.abs(fa - fb) / fa) < 1e-2


# ----------------------------------------------------------------------
# l'entree "trajectoire"
# ----------------------------------------------------------------------
def test_sur_trajectoire_vaut_la_decoupe_a_la_main(wscc9, trajectoire):
    """L'adaptateur doit rendre EXACTEMENT ce que fait la plomberie manuelle.

    C'est sa seule raison d'etre : laisser chaque appelant decouper la
    fenetre et empiler les voies lui-meme est le meilleur moyen que deux
    appelants ne le fassent pas pareil, et que l'un des deux se trompe d'un
    point de part et d'autre de l'evenement.
    """
    net, _ = wscc9
    m = ModeleNonLineaire(net)
    dt, duree = 0.005, 4.0
    t_eval = np.arange(0.0, duree + 0.5 * dt, dt)
    tr = m.simuler(
        (0.0, duree),
        uext_fn=m.scenario(Evenement(T_ECHELON, "SM1.p_ref", 0.01)),
        t_eval=t_eval,
    )
    noms = [f"SM{i}.Dwr" for i in (1, 2, 3)]
    par_adaptateur = sur_trajectoire(tr, noms, evenement=T_ECHELON)

    t, y, noms_ref = trajectoire
    a_la_main = identifier(t, y, voies=noms_ref)

    assert par_adaptateur.rang == a_la_main.rang
    assert par_adaptateur.residu == pytest.approx(a_la_main.residu, rel=1e-12)
    fa = sorted(m_.f for m_ in par_adaptateur.modes)
    fb = sorted(m_.f for m_ in a_la_main.modes)
    assert np.allclose(fa, fb, rtol=1e-12)


def test_sur_trajectoire_jette_ce_qui_precede_l_evenement(wscc9):
    """Avant l'evenement le signal est plat : l'inclure serait demander au
    modele en somme d'exponentielles de representer aussi le saut."""
    net, _ = wscc9
    m = ModeleNonLineaire(net)
    dt, duree = 0.005, 2.0
    t_eval = np.arange(0.0, duree + 0.5 * dt, dt)
    tr = m.simuler(
        (0.0, duree),
        uext_fn=m.scenario(Evenement(T_ECHELON, "SM1.p_ref", 0.01)),
        t_eval=t_eval,
    )
    noms = ["SM2.Dwr"]
    apres = sur_trajectoire(tr, noms, evenement=T_ECHELON)
    tout = sur_trajectoire(tr, noms)
    # la fenetre gardee est plus courte d'exactement l'avant-evenement
    assert len(apres.t) == len(tout.t) - int(round(T_ECHELON / dt))


def test_sur_trajectoire_refuse_un_etat_absent(wscc9):
    net, _ = wscc9
    m = ModeleNonLineaire(net)
    tr = m.simuler((0.0, 0.2), t_eval=np.linspace(0.0, 0.2, 41))
    with pytest.raises(ValueError, match="SM9.Dwr"):
        sur_trajectoire(tr, ["SM1.Dwr", "SM9.Dwr"])


def test_sur_trajectoire_refuse_une_fenetre_vide(wscc9):
    """Un evenement apres la fin de la simulation doit le DIRE."""
    net, _ = wscc9
    m = ModeleNonLineaire(net)
    tr = m.simuler((0.0, 0.2), t_eval=np.linspace(0.0, 0.2, 41))
    with pytest.raises(ValueError, match="points"):
        sur_trajectoire(tr, ["SM1.Dwr"], evenement=0.5)
