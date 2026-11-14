"""Tests de l'analyse modale (``g2egrid.modal``).

Chaque test confronte le module a une reference INDEPENDANTE, pas a lui-meme :
differences finies sur la matrice ``A`` pour la sensibilite, ``scipy.linalg.
expm`` pour la reponse libre, la discretisation exacte de
``g2egrid.separation`` pour la reponse indicielle, et ``numpy.linalg.eig`` pour
le spectre partiel. C'est la seule facon de detecter une erreur de convention
-- vecteurs a gauche mal apparies, transposition, facteur de normalisation --
qui laisserait le module coherent avec lui-meme tout en etant faux.
"""

import numpy as np
import pytest
from scipy.linalg import expm

from g2egrid import modal
from g2egrid.catalogue import construire
from g2egrid.separation import reponse_indicielle


# ----------------------------------------------------------------------
# fixtures
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def smib():
    """28 etats : le plus petit reseau du catalogue qui ait un mode local."""
    return construire("smib").build_state_space()


@pytest.fixture(scope="module")
def wscc9():
    """93 etats, trois machines : il y faut des modes INTER-machines."""
    return construire("wscc9").build_state_space()


def _electromecaniques(ss):
    """Les modes oscillants de 0,5 a 3 Hz, du plus "rotor" au moins.

    Le tri se fait sur la participation cumulee des vitesses rotor : c'est ce
    qui separe un mode electromecanique d'un mode d'excitation qui se trouve
    dans la meme bande de frequence.
    """
    lam = ss.eig()
    f = modal.frequence(lam)
    p = modal.participations(ss)
    iv = [i for i, n in enumerate(ss.states) if n.endswith(".Dwr")]
    garde = np.flatnonzero((lam.imag > 1e-3) & (f > 0.5) & (f < 3.0))
    return sorted((int(k) for k in garde), key=lambda k: -p[k, iv].sum())


# ----------------------------------------------------------------------
# sensibilite aux coefficients de A
# ----------------------------------------------------------------------
def test_sensibilite_vs_differences_finies(smib):
    """``dlambda/dA[i,j]`` doit coller au deplacement reel de la valeur propre.

    On perturbe ``A[i, j]`` de 1e-6 en RELATIF -- les coefficients de ``A``
    couvrent ici huit decades, de la constante de temps d'un PSS aux 1e5 rad/s
    d'une capacite de bus, et un pas absolu serait tantot du bruit tantot une
    perturbation d'ordre un. La valeur propre perturbee est reappariee par
    proximite : a ce pas, aucune autre n'est a portee.
    """
    A = smib.A
    k = _electromecaniques(smib)[0]
    lam = smib.eig()[k]
    S = modal.sensibilite(smib, k)

    # les coefficients les plus sensibles : ceux qu'un ingenieur irait regler
    for plat in np.argsort(-np.abs(S), axis=None)[:6]:
        i, j = np.unravel_index(plat, S.shape)
        assert A[i, j] != 0.0, "le masque doit ecarter les zeros structurels"
        h = 1e-6 * abs(A[i, j])
        Ap = A.copy()
        Ap[i, j] += h
        lp = np.linalg.eigvals(Ap)
        deplace = lp[np.argmin(np.abs(lp - lam))]
        num = (deplace - lam) / h
        err = abs(num - S[i, j]) / abs(S[i, j])
        assert err < 2e-3, (
            f"sensibilite fausse en ({smib.states[i]}, {smib.states[j]}) : "
            f"analytique {S[i, j]:.4e}, differences finies {num:.4e}"
        )


def test_sensibilite_masque_les_zeros_structurels(smib):
    """Un coefficient nul de ``A`` doit avoir une sensibilite nulle.

    Sans ce masque, les plus grandes valeurs du tableau designent des couples
    d'etats que le modele ne relie par aucune equation : la sensibilite y est
    parfaitement definie et parfaitement inutilisable, puisque aucun parametre
    physique ne peut faire bouger ce coefficient.
    """
    S = modal.sensibilite(smib, _electromecaniques(smib)[0])
    assert np.all(S[smib.A == 0.0] == 0.0)
    assert np.any(np.abs(S[smib.A != 0.0]) > 0.0)


def test_sensibilite_hors_spectre_leve(smib):
    with pytest.raises(IndexError):
        modal.sensibilite(smib, len(smib.states))


def test_sensibilite_top_nomme_les_etats(smib):
    """``sensibilite_top`` doit rendre des couples (ligne, colonne) ORDONNES.

    ``A[i, j]`` est l'effet de l'etat ``j`` sur la derivee de l'etat ``i`` : le
    couple rendu designe un coefficient d'equation, pas une paire symetrique.
    """
    k = _electromecaniques(smib)[0]
    top = modal.sensibilite_top(smib, k, n=5)
    S = np.abs(modal.sensibilite(smib, k))
    assert len(top) == 5
    assert [v for _, _, v in top] == sorted((v for _, _, v in top), reverse=True)
    for ligne, colonne, v in top:
        i, j = smib.states.index(ligne), smib.states.index(colonne)
        assert S[i, j] == pytest.approx(v)


# ----------------------------------------------------------------------
# reponses en forme close
# ----------------------------------------------------------------------
def test_reponse_libre_vs_expm(smib):
    """``Vr.diag(exp(lam.t)).Wl.x0`` doit valoir ``expm(A.t).x0``.

    ``expm`` (Pade avec mise a l'echelle et elevations au carre) ne partage
    aucun code avec la decomposition modale : si les vecteurs a gauche etaient
    mal apparies a ceux de droite, ou ``Wl`` transposee, les deux differeraient
    des le premier instant non nul.
    """
    x0 = np.random.default_rng(0).standard_normal(len(smib.states))
    t = np.concatenate([np.linspace(0.0, 1e-3, 4), np.linspace(1e-2, 1.0, 6)])
    a = modal.reponse_libre(smib, x0, t)
    b = np.array([expm(smib.A * tk) @ x0 for tk in t])
    assert a.shape == (len(t), len(smib.states))
    assert np.abs(a - b).max() / np.abs(b).max() < 1e-8


def test_reponse_libre_part_de_x0(smib):
    """A ``t = 0`` la forme close doit redonner ``x0`` exactement."""
    x0 = np.random.default_rng(1).standard_normal(len(smib.states))
    a = modal.reponse_libre(smib, x0, np.array([0.0]))
    assert np.abs(a[0] - x0).max() < 1e-10


def test_reponse_libre_refuse_une_taille_fausse(smib):
    with pytest.raises(ValueError):
        modal.reponse_libre(smib, np.zeros(3), np.array([0.0, 1.0]))


def test_reponse_indicielle_modale_vs_discretisation_exacte(smib):
    """La somme modale doit valoir la reponse indicielle par bloqueur d'ordre zero.

    ``separation.reponse_indicielle`` passe par l'exponentielle de la matrice
    augmentee : aucune valeur propre n'y intervient. Les deux voies partent
    donc d'objets differents pour arriver au meme ``x(t)``.
    """
    u = np.zeros(len(smib.inputs))
    u[smib.iu("SM1.p_ref")] = 1.0
    t = np.linspace(0.0, 2.0, 401)
    xm = modal.reponse_indicielle_modale(smib, u, t)
    xs = reponse_indicielle(smib.A, smib.B, t, u)
    assert xm.shape == xs.shape == (len(t), len(smib.states))
    assert np.abs(xm - xs).max() / np.abs(xs).max() < 1e-9


def test_contributions_par_mode_somment_a_la_reponse(smib):
    """``par_mode`` doit DECOMPOSER la reponse, pas en donner une autre.

    C'est tout l'interet de la fonction : dire quel mode porte le depassement.
    Si la somme des parts ne redonnait pas la trajectoire, la lecture serait
    fausse -- et le mode de jauge, dont la contribution est le passage a la
    limite ``t`` et non ``(exp(lam.t)-1)/lam``, est justement l'endroit ou une
    erreur de cas limite passerait inapercue.
    """
    u = np.zeros(len(smib.inputs))
    u[smib.iu("SM1.p_ref")] = 1.0
    t = np.linspace(0.0, 1.0, 101)
    etats = ["SM1.Dwr", "SM1.Pm", "SM1.theta"]
    x, contributions = modal.reponse_indicielle_modale(
        smib, u, t, par_mode=True, etats=etats
    )
    assert contributions.shape == (len(smib.states), len(t), len(etats))
    assert np.abs(contributions.sum(axis=0) - x).max() < 1e-12
    seul = modal.reponse_indicielle_modale(smib, u, t, etats=etats)
    assert np.abs(seul - x).max() < 1e-12


def test_reponse_indicielle_modale_refuse_une_matrice_nue(smib):
    """Sans ``B``, la question n'a pas de sens : il faut un ``StateSpace``."""
    with pytest.raises(TypeError):
        modal.reponse_indicielle_modale(smib.A, np.zeros(2), np.array([0.0]))


# ----------------------------------------------------------------------
# spectre partiel
# ----------------------------------------------------------------------
def test_valeurs_propres_partielles_vs_spectre_dense(smib):
    """Le shift-and-invert doit rendre les k valeurs propres les plus proches.

    On compare les DISTANCES a sigma, pas les valeurs propres appariees une a
    une : quand la k-ieme et la (k+1)-ieme sont a egale distance -- et c'est
    le cas des que la coupure tombe au milieu d'une paire conjuguee, ce qui
    arrive sur SMIB a sigma = -0,1 avec k = 8, ou les deux valent 11,3649 --
    prendre l'une ou l'autre sont deux reponses egalement justes. Exiger un
    appariement ferait echouer le test sur un comportement correct.

    On verifie en revanche que chaque valeur rendue appartient bien au spectre
    et que son vecteur propre va avec : c'est la seule chose qui ne se negocie
    pas.
    """
    lam = np.linalg.eigvals(smib.A)
    for sigma in (-0.1, -1.0, -20.0):
        for k in (4, 8):
            lp, Vr = modal.valeurs_propres_partielles(smib, k, sigma)
            d_lu = np.sort(np.abs(lp - sigma))
            d_ref = np.sort(np.abs(lam - sigma))[:k]
            assert np.abs(d_lu - d_ref).max() < 1e-8, f"sigma={sigma}, k={k}"
            # chaque valeur rendue est bien une valeur propre de A
            for v in lp:
                assert np.abs(lam - v).min() < 1e-8, f"{v} hors du spectre"
            # et les vecteurs propres vont avec les valeurs propres
            assert np.abs(smib.A @ Vr - Vr * lp).max() < 1e-6


def test_valeurs_propres_partielles_borne_arpack(smib):
    n = len(smib.states)
    for k in (0, n - 1, n):
        with pytest.raises(ValueError):
            modal.valeurs_propres_partielles(smib, k, -1.0)


def test_spectre_partiel_a_sigma_nul_bute_sur_le_mode_de_jauge(smib):
    """``sigma = 0`` est impraticable des que le reseau a un mode nul.

    Un reseau sans bus infini -- et meme ce SMIB, dont l'angle de la source est
    un integrateur pur -- a une valeur propre exactement nulle. ``A - 0.I`` est
    alors exactement singuliere et la factorisation de SuperLU echoue. Ce test
    fige le comportement pour qu'on ne le decouvre pas en production : la
    docstring de ``valeurs_propres_partielles`` recommande ``sigma = 0`` "pour
    les modes lents", ce qui est precisement le cas ou il casse. Prendre un
    ``sigma`` legerement decale (-0,1) donne le meme resultat sans le piege.
    """
    assert np.min(np.abs(np.linalg.eigvals(smib.A))) < 1e-8
    lam, _ = modal.valeurs_propres_partielles(smib, 4, 0.0)
    ref = np.linalg.eigvals(smib.A)
    ref = ref[np.argsort(np.abs(ref))][:4]
    assert np.abs(np.sort_complex(lam) - np.sort_complex(ref)).max() < 1e-8


def test_spectre_partiel_shift_complexe(smib):
    """Le cas d'usage documente : ``sigma = -1 + 2j.pi.f`` autour de 1,5 Hz.

    C'est ici que se voit le piege du shift-invert REEL d'ARPACK : devant une
    matrice reelle il classe par resonance et non par distance a ``sigma``, et
    rendait des modes a 15 et 27 de distance au lieu de 9,4. Le correctif est
    dans ``modal.valeurs_propres_partielles`` : passer A en complexe.
    """
    lam = np.linalg.eigvals(smib.A)
    sigma = -1.0 + 2j * np.pi * 1.5
    lp, _ = modal.valeurs_propres_partielles(smib, 6, sigma)
    ref = lam[np.argsort(np.abs(lam - sigma))[:6]]
    assert np.abs(np.sort_complex(lp) - np.sort_complex(ref)).max() < 1e-8


# ----------------------------------------------------------------------
# forme modale
# ----------------------------------------------------------------------
def test_mode_inter_machines_en_opposition_de_phase(wscc9):
    """Deux vitesses rotor a 180 degres l'une de l'autre : c'est la definition.

    Le vecteur propre a droite n'est defini qu'a un facteur complexe pres, mais
    la phase RELATIVE de deux composantes, elle, ne depend d'aucune
    normalisation : c'est ce que ``forme_modale`` doit rendre lisible. Sur les
    deux modes du WSCC 9 bus ou les vitesses rotor participent le plus -- 2,20
    et 2,01 Hz -- deux machines doivent osciller l'une contre l'autre.
    """
    lam = wscc9.eig()
    for k in _electromecaniques(wscc9)[:2]:
        forme = modal.forme_modale(wscc9, k, n=12)
        vitesses = [(n, mod, ph) for n, mod, ph in forme if n.endswith(".Dwr")]
        assert len(vitesses) >= 2, f"mode {k} : une seule vitesse rotor dans la forme"
        (na, _, pha), (nb, _, phb) = vitesses[:2]
        ecart = abs((phb - pha) % 360.0 - 180.0)
        assert ecart < 30.0, (
            f"mode {k} a {modal.frequence(lam[k])[0]:.2f} Hz : {na} et {nb} "
            f"a {phb - pha:+.1f} degres, ce n'est pas un mode inter-machines"
        )


def test_forme_modale_normalisee_sur_l_etat_dominant(wscc9):
    """L'etat le plus participant sert de reference : module 1, phase 0."""
    forme = modal.forme_modale(wscc9, _electromecaniques(wscc9)[0], n=5)
    nom, module, phase = forme[0]
    assert module == pytest.approx(1.0)
    assert phase == pytest.approx(0.0, abs=1e-9)
    assert all(-180.0 < ph <= 180.0 for _, _, ph in forme)


# ----------------------------------------------------------------------
# amortissement, frequence, modes
# ----------------------------------------------------------------------
def test_amortissement_et_frequence_coherents_avec_modes(smib):
    """``modes()`` ne doit pas recalculer autrement ce que les deux fonctions
    vectorielles donnent, et doit trier du moins au plus amorti."""
    lam = smib.eig()
    zeta, f = modal.amortissement(lam), modal.frequence(lam)
    liste = modal.modes(smib)
    assert len(liste) == len(lam)
    for mode in liste:
        assert mode.valeur == lam[mode.indice]
        assert mode.zeta == pytest.approx(zeta[mode.indice])
        assert mode.f == pytest.approx(f[mode.indice])
        assert mode.stable == (mode.valeur.real < 0)
        assert mode.dominant == mode.etats[0][0]
    assert [m.zeta for m in liste] == sorted(m.zeta for m in liste)


def test_amortissement_est_bien_zeta(smib):
    """``zeta = -Re(lam)/|lam|``, et 100 % pour la valeur propre nulle.

    Le mode de jauge n'a pas d'amortissement defini ; le compter a 100 % le
    range avec les modes non oscillants plutot qu'en tete de la liste des
    modes critiques, ou il masquerait le vrai mode le moins amorti.
    """
    lam = np.array([-1 + 1j, -1 - 1j, -3 + 0j, 0j])
    attendu = 100.0 * np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 1.0, 1.0])
    assert modal.amortissement(lam) == pytest.approx(attendu)
    assert modal.frequence(lam) == pytest.approx(
        np.array([1.0, 1.0, 0.0, 0.0]) / (2 * np.pi)
    )


def test_participations_normalisees_par_mode(smib):
    """Chaque ligne de la matrice de participation somme a 1."""
    lam, Vr, p = modal.spectre(smib)
    assert p.shape == (len(lam), len(lam))
    assert np.all(p >= 0.0)
    assert p.sum(axis=1) == pytest.approx(np.ones(len(lam)))
    assert np.abs(modal.participations(smib) - p).max() == 0.0


def test_modes_filtre_sur_frequence_et_amortissement(smib):
    """Les filtres de ``modes()`` doivent se composer, pas s'ignorer."""
    tous = modal.modes(smib)
    lents = modal.modes(smib, f_max=3.0)
    peu_amortis = modal.modes(smib, zeta_max=50.0)
    oscillants = modal.modes(smib, oscillants=True)
    assert 0 < len(lents) < len(tous)
    assert 0 < len(peu_amortis) < len(tous)
    assert all(m.f <= 3.0 for m in lents)
    assert all(m.zeta <= 50.0 for m in peu_amortis)
    assert all(m.valeur.imag > 0 and m.oscillant for m in oscillants)
