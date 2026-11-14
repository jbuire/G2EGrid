"""Le balayage parametrique : appliquer une cible, suivre les modes.

Ce qui peut casser ici n'est pas le calcul modal -- il est teste ailleurs --
mais les trois choses que le balayage ajoute par-dessus : designer un
parametre sans toucher aux tables d'origine, suivre un mode d'un pas au
suivant, et ne pas confondre le mode de jauge avec un mode physique.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2egrid.catalogue import charger
from g2egrid_api import analyse


@pytest.fixture(scope="module")
def smib_gfm():
    return charger("smib_gfm")


# ----------------------------------------------------------------------
# designer un parametre
# ----------------------------------------------------------------------
def test_appliquer_ne_touche_pas_l_original(smib_gfm):
    """La cible s'applique a une COPIE.

    Le cache des reseaux assembles est indexe par une empreinte des tables :
    muter l'original ferait repondre le cache pour un reseau qu'il n'a pas
    assemble, et le balayage rendrait le meme spectre a tous les pas.
    """
    avant = dict(smib_gfm.GFM_PARAMS)
    t2 = analyse._appliquer(
        smib_gfm, {"genre": "famille", "famille": "GFM_PARAMS", "cle": "mp"}, 0.123
    )
    assert t2.GFM_PARAMS["mp"] == 0.123
    assert smib_gfm.GFM_PARAMS == avant, "les tables d'origine ont ete modifiees"


def test_appliquer_par_bus(smib_gfm):
    bus = int(smib_gfm.Y_DER[0][0])
    t2 = analyse._appliquer(smib_gfm, {"genre": "bus", "bus": bus, "cle": "mp"}, 0.077)
    assert t2.PARAMS_PAR_BUS[bus]["mp"] == 0.077
    assert (
        bus not in smib_gfm.PARAMS_PAR_BUS or "mp" not in smib_gfm.PARAMS_PAR_BUS[bus]
    )


def test_appliquer_cellule_de_table(smib_gfm):
    """Une cellule de table, designee par son NOM de colonne.

    Un indice numerique se decalerait a la premiere colonne ajoutee ; le nom
    est ce que le contrat ``TABLES`` garantit.
    """
    table = "T_Lines_pu" if smib_gfm.T_Lines_pu else "T_Lines_SI"
    col = analyse.TABLES[table][3]
    t2 = analyse._appliquer(
        smib_gfm,
        {"genre": "table", "table": table, "ligne": 0, "colonne": col},
        0.42,
    )
    j = analyse.TABLES[table].index(col)
    assert getattr(t2, table)[0][j] == 0.42
    assert getattr(smib_gfm, table)[0][j] != 0.42
    # et le type de la ligne est conserve
    assert type(getattr(t2, table)[0]) is type(getattr(smib_gfm, table)[0])


@pytest.mark.parametrize(
    "cible, motif",
    [
        ({"genre": "chose"}, "genre"),
        ({"genre": "famille", "famille": "XX_PARAMS", "cle": "mp"}, "famille"),
        (
            {"genre": "table", "table": "T_Inexistante", "ligne": 0, "colonne": "a"},
            "table",
        ),
        (
            {"genre": "table", "table": "Y_DER", "ligne": 999, "colonne": "Type"},
            "ligne",
        ),
        ({"genre": "table", "table": "Y_DER", "ligne": 0, "colonne": "Zzz"}, "colonne"),
    ],
)
def test_cible_invalide_leve_un_message_utile(smib_gfm, cible, motif):
    with pytest.raises(ValueError, match=motif):
        analyse._appliquer(smib_gfm, cible, 1.0)


# ----------------------------------------------------------------------
# suivre un mode
# ----------------------------------------------------------------------
def test_suivre_retrouve_une_permutation():
    """L'appariement rend l'ordre qui remet les modes en face des bons."""
    prec = np.array([-1 + 10j, -2 + 20j, -3 + 30j])
    melange = prec[[2, 0, 1]] + 0.01  # meme spectre, ordre change, bruit
    ordre = analyse._suivre(prec, melange)
    assert np.abs(melange[ordre] - prec).max() < 0.02


def _glouton(prec: np.ndarray, cour: np.ndarray) -> np.ndarray:
    """L'appariement naif : servir d'abord le couple le plus proche.

    Present ici pour que le test contraste deux comportements calcules, et
    non un comportement calcule contre un resultat que j'aurais recopie a la
    main -- lequel peut etre faux sans que rien ne le signale.
    """
    cout = np.abs(prec[:, None] - cour[None, :])
    ordre = np.full(len(cour), -1)
    libres_l, libres_c = set(range(len(prec))), set(range(len(cour)))
    while libres_l:
        i, j = min(
            ((i, j) for i in libres_l for j in libres_c), key=lambda ij: cout[ij]
        )
        ordre[i] = j
        libres_l.discard(i)
        libres_c.discard(j)
    return ordre


def test_suivre_bat_le_glouton_sur_un_croisement():
    """Deux modes qui se croisent : le glouton se trompe, l'appariement non.

    C'est exactement le cas interessant d'un lieu des racines -- deux modes
    se rencontrent puis repartent -- et le seul ou le choix d'un appariement
    GLOBAL se voit. Ici le glouton sert d'abord (-1,2 -> -1,15), le couple le
    moins cher, et se retrouve force de relier -1 a -1,6 : cout total 0,65,
    contre 0,55 pour l'appariement global.
    """
    prec = np.array([-1 + 0j, -1.2 + 0j])
    cour = np.array([-1.6 + 0j, -1.15 + 0j])
    cout = np.abs(prec[:, None] - cour[None, :])

    o_global = analyse._suivre(prec, cour)
    o_glouton = _glouton(prec, cour)
    total = lambda o: cout[np.arange(len(prec)), o].sum()  # noqa: E731

    assert not np.array_equal(o_global, o_glouton), "le cas ne distingue rien"
    assert total(o_global) < total(o_glouton)
    # et c'est bien l'optimum : avec deux modes, il n'y a que deux choix
    assert total(o_global) == pytest.approx(min(0.6 + 0.05, 0.15 + 0.4))


# ----------------------------------------------------------------------
# le balayage lui-meme
# ----------------------------------------------------------------------
def test_balayage_lieux_continus_et_paire_conjuguee(smib_gfm):
    """Un lieu est une courbe, et une paire conjuguee traverse ensemble.

    Deux controles qui tiennent ensemble : si le suivi etait faux, les lieux
    sauteraient d'un pas a l'autre, et les deux moities d'une paire
    conjuguee ne croiseraient pas l'axe a la meme valeur.
    """
    v = list(np.linspace(0.005, 0.10, 30))
    d = analyse.balayage(
        smib_gfm,
        {"genre": "famille", "famille": "GFM_PARAMS", "cle": "mp"},
        v,
        f_max=60.0,
        n_modes=6,
    )
    assert len(d["pas"]) == len(v)
    assert all("modes" in p for p in d["pas"]), "un pas n'a pas abouti"

    L = np.array([[complex(m["re"], m["im"]) for m in p["modes"]] for p in d["pas"]])
    # continuite : un saut d'un pas au suivant reste petit devant l'etendue
    saut = np.abs(np.diff(L, axis=0))
    etendue = np.abs(L.max() - L.min())
    assert saut.max() < 0.25 * etendue, (
        f"saut {saut.max():.2f} pour une etendue {etendue:.2f}"
    )

    # la paire conjuguee qui devient instable croise a la meme valeur
    croisements = []
    for j in range(L.shape[1]):
        re = L[:, j].real
        if re.max() > 0 > re.min():
            k = int(np.argmax(re > 0))
            croisements.append(
                v[k - 1] + (v[k] - v[k - 1]) * (-re[k - 1]) / (re[k] - re[k - 1])
            )
    assert len(croisements) == 2, f"{len(croisements)} traversees, attendu une paire"
    assert abs(croisements[0] - croisements[1]) < 1e-12


def test_balayage_ecarte_le_mode_de_jauge():
    """Un reseau sans bus infini porte un zero exact : ce n'est pas un mode.

    Son amortissement vaut -Re/|lambda|, soit du bruit d'arrondi divise par
    du bruit d'arrondi : il arrive en tete du classement par amortissement,
    occupe un lieu a lui seul et fait afficher "max Re = +0.000" sur un
    reseau parfaitement stable.
    """
    t = charger("wscc9_gfm")
    d = analyse.balayage(
        t,
        {"genre": "famille", "famille": "GFM_PARAMS", "cle": "mp"},
        [0.02, 0.04],
        f_max=60.0,
        n_modes=6,
    )
    for p in d["pas"]:
        assert p["jauge"] == 1, "le mode de jauge n'a pas ete reconnu"
        assert all(abs(complex(m["re"], m["im"])) > analyse.JAUGE for m in p["modes"])
        # sans le filtre, max_re valait le bruit de la jauge, autour de 0
        assert p["max_re"] < -1e-3, f"max_re = {p['max_re']}, la jauge domine encore"


def test_balayage_borne_le_nombre_de_pas(smib_gfm):
    """Au-dela d'une soixantaine de pas, c'est un script qu'il faut."""
    cible = {"genre": "famille", "famille": "GFM_PARAMS", "cle": "mp"}
    with pytest.raises(ValueError, match="valeurs"):
        analyse.balayage(smib_gfm, cible, list(np.linspace(0.01, 0.1, 61)))
    with pytest.raises(ValueError, match="valeurs"):
        analyse.balayage(smib_gfm, cible, [])


def test_balayage_garde_les_pas_sans_point_de_fonctionnement(smib_gfm):
    """Un pas qui ne converge pas est un RESULTAT, pas un trou.

    C'est souvent la fin du domaine de validite ; l'effacer ferait croire
    que le balayage s'arrete sans raison.
    """
    # une reactance de ligne absurde : le load flow n'a plus de solution
    table = "T_Lines_pu" if smib_gfm.T_Lines_pu else "T_Lines_SI"
    col = next(c for c in analyse.TABLES[table] if c.startswith("X"))
    d = analyse.balayage(
        smib_gfm,
        {"genre": "table", "table": table, "ligne": 0, "colonne": col},
        [0.1, 1e6],
        f_max=60.0,
    )
    assert len(d["pas"]) == 2
    rates = [p for p in d["pas"] if "erreur" in p]
    assert rates, "aucun pas en echec : le cas limite n'est pas atteint"
    assert rates[0]["erreur"], "un echec sans message n'aide personne"
