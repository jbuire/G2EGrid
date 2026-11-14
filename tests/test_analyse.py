"""Les fonctions d'analyse de l'interface, appelees SANS HTTP.

C'est la regle d'architecture du module : ``g2egrid_api/analyse.py`` ne sait
pas ce qu'est HTTP, ``app.py`` n'est qu'une couche de transport. Ces tests en
sont la consequence -- et la preuve : tout ce que fait l'interface se teste
et s'appelle depuis un script, par les memes fonctions.

Ce qui est verifie ici n'est pas le calcul modal ni la simulation, testes
ailleurs, mais ce que cette couche ajoute par-dessus : le choix des defauts,
les conversions, et surtout les endroits ou elle refuse de rendre un resultat
qui aurait l'air juste sans l'etre.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2egrid.catalogue import charger
from g2egrid.tables import TABLES, Tables
from g2egrid_api import analyse


@pytest.fixture(scope="module")
def smib():
    return charger("smib")


# ----------------------------------------------------------------------
# le detail d'un mode
# ----------------------------------------------------------------------
def test_detail_mode_phase_de_reference(smib):
    """L'etat le plus participant sert de reference de phase.

    Le vecteur propre a droite n'est defini qu'a un facteur complexe pres :
    ``eig`` le rend de norme 1 mais de phase arbitraire, qui change d'une
    version de LAPACK a l'autre. Sans cette normalisation, le diagramme
    polaire tournerait au gre de la bibliotheque d'algebre.
    """
    d = analyse.modes(smib, f_max=5.0)
    k = d["modes"][0]["indice"]
    det = analyse.detail_mode(smib, k)
    assert det["forme"][0]["module"] == pytest.approx(1.0)
    assert det["forme"][0]["phase"] == pytest.approx(0.0, abs=1e-9)
    assert all(-180.0 < e["phase"] <= 180.0 for e in det["forme"])


def test_detail_mode_est_le_meme_mode_que_la_table(smib):
    """L'indice est celui du SPECTRE, pas le rang dans la table affichee.

    La table est triee par amortissement et filtree en frequence : un rang
    dedans ne designe rien de stable. Confondre les deux ouvrirait la forme
    modale d'un autre mode que celui qu'on a clique.
    """
    d = analyse.modes(smib, f_max=5.0)
    for m in d["modes"][:3]:
        det = analyse.detail_mode(smib, m["indice"])
        assert det["f"] == pytest.approx(m["f"], rel=1e-12)
        assert det["zeta"] == pytest.approx(m["zeta"], rel=1e-12)
        assert det["re"] == pytest.approx(m["re"], rel=1e-12)


def test_detail_mode_hors_spectre(smib):
    with pytest.raises(ValueError, match="hors du spectre"):
        analyse.detail_mode(smib, 10_000)


def test_sensibilites_designent_un_coefficient_du_modele(smib):
    """``A[i, j]`` : l'etat COLONNE agit sur la derivee de l'etat LIGNE.

    Ce n'est pas une paire symetrique, et le dire compte : la sensibilite
    designe une equation ou agir, pas deux etats qui se ressemblent.
    """
    d = analyse.modes(smib, f_max=5.0)
    det = analyse.detail_mode(smib, d["modes"][0]["indice"], n_sensibilites=6)
    assert len(det["sensibilites"]) == 6
    v = [e["valeur"] for e in det["sensibilites"]]
    assert v == sorted(v, reverse=True), "les sensibilites ne sont pas triees"
    etats = set(analyse.signaux(smib)["etats"])
    for e in det["sensibilites"]:
        assert e["ligne"] in etats and e["colonne"] in etats


# ----------------------------------------------------------------------
# mesures et scenarios
# ----------------------------------------------------------------------
def test_noms_mesures_suit_le_reseau(smib):
    """Le catalogue des mesures se deduit du reseau ASSEMBLE.

    L'ecrire a la main dans l'interface en ferait un doublon qui derive : un
    bus de plus, et il manque quatre mesures.
    """
    noms = analyse.noms_mesures(smib)
    bus = [int(r[0]) for r in smib.Y_DER]
    for b in bus:
        assert f"V_bus{b}" in noms
        assert f"angle_bus{b}" in noms
    assert noms == sorted(noms)
    assert analyse.signaux(smib)["mesures"] == noms


def test_mesures_a_l_equilibre_reproduisent_le_load_flow(smib):
    """Sans perturbation, les mesures doivent redonner le point de depart.

    C'est le controle croise entre deux chemins qui ne partagent rien : le
    load flow resout un systeme algebrique en complexes, les mesures
    relisent les etats dq du modele non lineaire.
    """
    lf = analyse.loadflow(smib)
    quoi = [f"V_bus{b['id']}" for b in lf["bus"]]
    d = analyse.temporel(
        smib, "SM1.p_ref", 0.0, duree=0.1, n=11, voie="non_lineaire", mesures=quoi
    )
    for b in lf["bus"]:
        mesure = np.asarray(d["mesures"][f"V_bus{b['id']}"])
        assert mesure[0] == pytest.approx(b["V"], rel=1e-6), f"bus {b['id']}"
        assert mesure.std() < 1e-9, "la trajectoire non perturbee derive"


def test_mesures_seulement_sur_la_voie_non_lineaire(smib):
    """Les mesures se reconstruisent a partir des (x, z, u, y) re-resolus.

    Une exponentielle de matrice sur les etats ne les fournit pas : la voie
    lineaire ne peut pas les rendre, et pretendre le contraire vaudrait
    mieux qu'un silence -- d'ou le fait qu'elles ne sortent QUE du non
    lineaire.
    """
    d = analyse.temporel(
        smib, "SM1.p_ref", 0.01, duree=0.2, n=21, voie="lineaire", mesures=["V_bus2"]
    )
    assert d["mesures"] == {}
    assert d["voies"] == ["lineaire"]


def test_mesure_inconnue_est_refusee(smib):
    with pytest.raises(ValueError, match="mesure"):
        analyse.temporel(smib, "SM1.p_ref", 0.01, duree=0.1, n=11, mesures=["V_bus999"])


def test_scenario_coupe_l_integration_a_la_date(smib):
    """Avant le premier evenement, RIEN ne doit bouger.

    ``simuler`` decoupe l'integration aux discontinuites et repart de l'etat
    atteint. Sans ce decoupage, Radau franchirait la marche avec un pas trop
    grand et laisserait fuir du transitoire en amont de la date.
    """
    ev = [
        {"instant": 0.5, "entree": "SM1.p_ref", "amplitude": 0.02},
        {"instant": 1.0, "entree": "SM1.p_ref", "amplitude": -0.02},
    ]
    d = analyse.temporel(
        smib, "SM1.p_ref", duree=1.5, n=151, evenements=ev, mesures=["f_bus2"]
    )
    t = np.asarray(d["t"])
    f = np.asarray(d["mesures"]["f_bus2"])
    avant = f[t < 0.5 - 1e-12]
    apres = f[t >= 0.5]
    # On compare les deux plutot que de fixer un seuil : ce qui subsiste
    # avant la date n'est pas nul -- c'est le bruit de l'integrateur, 3e-8 Hz
    # sur 50, reconstruit par une division par |v|^2 -- mais il doit etre
    # sans commune mesure avec l'excursion que l'evenement provoque, de
    # l'ordre du millihertz. Deux ordres de grandeur separent les deux, et
    # c'est cela qui a un sens, pas une tolerance ecrite a la main.
    bruit = float(np.abs(avant - avant[0]).max())
    excursion = float(np.abs(apres - avant[0]).max())
    assert excursion > 100 * bruit, (
        f"bruit avant la date {bruit:.2e} Hz, excursion apres {excursion:.2e} Hz"
    )
    assert len(d["evenements"]) == 2


def test_scenario_a_plusieurs_evenements_omet_la_voie_lineaire(smib):
    """La reponse lineaire est calculee en forme close : elle ne suit qu'un
    echelon unique. Au-dela, l'omettre vaut mieux que la rendre fausse."""
    ev = [
        {"instant": 0.1, "entree": "SM1.p_ref", "amplitude": 0.01},
        {"instant": 0.3, "entree": "SM1.p_ref", "amplitude": 0.01},
    ]
    d = analyse.temporel(smib, "SM1.p_ref", duree=0.5, n=51, evenements=ev)
    assert d["voies"] == ["non_lineaire"]
    assert d["ecart"] is None


def test_echelon_simple_garde_les_deux_voies(smib):
    """Le cas courant ne doit pas payer l'existence du cas general."""
    d = analyse.temporel(smib, "SM1.p_ref", 0.01, duree=0.5, n=51)
    assert d["voies"] == ["lineaire", "non_lineaire"]
    assert d["ecart"]["relatif"] < 0.05


# ----------------------------------------------------------------------
# validation
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def validation_smib(smib):
    return analyse.validation(smib, duree=3.0)


def test_validation_point_fixe(validation_smib):
    """Le point de fonctionnement EST un point fixe du modele non lineaire.

    Sinon on linearise autour d'un point qui n'existe pas : invisible en
    modal, fatal en temporel, ou la trajectoire derive sans perturbation.
    """
    pf = validation_smib["point_fixe"]
    assert pf["residu_x"] < 1e-6
    assert pf["residu_algebrique"] < 1e-6
    assert pf["verdict"]
    # aucun convertisseur hors limite (FR-009, tests/test_limite_equilibre.py)
    assert validation_smib["limites"] == {"depassements": [], "verdict": True}


def test_validation_ecart_en_carre_de_l_amplitude(validation_smib):
    """Doubler l'amplitude doit QUADRUPLER l'ecart entre les deux voies.

    C'est la definition d'une linearisation au premier ordre, et le test le
    plus discriminant du lot : un rapport proche de 2 dirait que l'ecart est
    du premier ordre, c'est-a-dire que ``A`` n'est pas la bonne matrice.
    """
    c = validation_smib["carre"]
    assert c["rapport"] == pytest.approx(4.0, rel=0.1)
    assert c["verdict"]


def test_validation_identifie_les_modes_du_spectre(validation_smib):
    """Les modes lus sur la trajectoire non lineaire sont ceux de ``A``.

    Les deux chemins ne partagent aucune ligne de code de resolution : leur
    accord ne peut pas venir d'une erreur commune. C'est la seule des trois
    epreuves qui ait un sens physique.
    """
    i = validation_smib["identification"]
    assert i["paires"], "aucun mode apparie"
    assert i["residu"] < 1e-3
    for p in i["paires"]:
        assert abs(p["ecart_f"]) < 2.0, f"{p['f_predit']} Hz : {p['ecart_f']} % en f"
        assert abs(p["ecart_zeta_rel"]) < 15.0, f"{p['f_predit']} Hz en zeta"


def test_validation_ecarte_ce_qui_n_est_pas_une_oscillation(smib):
    """Un mode amorti a plus de 95 % decroit avant d'avoir boucle une periode.

    Le filtre ``oscillants`` ne regarde que la partie imaginaire et garde de
    telles paires ; les exiger de l'identification ferait echouer la
    validation sur une propriete du SIGNAL, pas du modele.
    """
    d = analyse.validation(smib, duree=2.0)
    i = d["identification"]
    vus = [p["zeta_predit"] for p in i["paires"]]
    vus += [x["zeta"] for x in i["predits_seuls"]]
    assert all(z <= analyse.ZETA_OSCILLANT for z in vus)


def test_validation_refuse_une_entree_inconnue(smib):
    with pytest.raises(ValueError, match="entree inconnue"):
        analyse.validation(smib, entree="pas_une_entree")


def test_validation_choisit_une_vraie_consigne(validation_smib):
    """Pas ``w_ref``, qui est l'entree GLOBALE.

    Un echelon de 0,01 pu sur la pulsation du repere deplace le repere de
    tout le reseau d'un demi-hertz : les deux voies s'ecartent alors de
    386 % sur ce reseau, contre 0,15 % sur une vraie consigne. Le defaut
    aurait fait conclure que le modele lineaire est faux.
    """
    assert validation_smib["entree"] != "w_ref"
    assert validation_smib["entree"].endswith(".p_ref")


# ----------------------------------------------------------------------
# validation d'un reseau mixte alternatif / continu
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nom", ["hvdc2t", "mtdc3t", "wscc9", "smib"])
def test_valider_ne_denonce_pas_un_reseau_correct(nom):
    """Un reseau du catalogue ne doit produire AUCUNE erreur.

    Le piege etait le continu : une liaison a courant continu separe les
    iles alternatives, et le validateur y voyait deux fautes qui n'en sont
    pas -- "deux bus slack" la ou chaque ile garde legitimement sa propre
    reference d'angle, et "bus isole" sur des terminaux relies par le cote
    continu, que les tables de branches alternatives ne montrent pas.
    """
    pb = [x for x in analyse.valider(charger(nom)) if x["gravite"] == "erreur"]
    assert pb == [], [x["message"] for x in pb]


def test_valider_exige_une_reference_par_ile():
    """Une ile alternative sans slack n'a pas de repere d'angle.

    C'est l'envers du test precedent : accepter plusieurs slacks ne doit pas
    revenir a n'en exiger aucun.
    """
    t = charger("hvdc2t").copie()
    lignes = [list(r) for r in t.Y_DER]
    for r in lignes:
        if int(r[1]) == 1:
            r[1] = 2
            break
    t.Y_DER = [tuple(r) for r in lignes]
    msgs = [x["message"] for x in analyse.valider(t) if x["gravite"] == "erreur"]
    assert any("aucun bus slack" in m for m in msgs), msgs
    assert any("ile" in m for m in msgs), "le message doit nommer l'ile fautive"


def test_valider_refuse_deux_slacks_dans_la_meme_ile():
    t = charger("wscc9").copie()
    lignes = [list(r) for r in t.Y_DER]
    lignes[1][1] = 1
    t.Y_DER = [tuple(r) for r in lignes]
    msgs = [x["message"] for x in analyse.valider(t) if x["gravite"] == "erreur"]
    assert any("2 bus slack" in m for m in msgs), msgs


# ----------------------------------------------------------------------
# la marge de stabilite se lit sur les modes physiques
# ----------------------------------------------------------------------
@pytest.mark.parametrize(
    "nom, jauge_attendue", [("smib", 1), ("wscc9", 1), ("hvdc2t", 3), ("mtdc3t", 5)]
)
def test_modes_ecarte_la_jauge_de_la_marge(nom, jauge_attendue):
    """``max Re`` ne doit pas etre le bruit d'arrondi d'un zero exact.

    Un reseau dont la reference d'angle n'est tenue par aucun bus infini
    porte un zero exact ; une liaison continue en fait porter plusieurs. Sur
    ``mtdc3t``, prendre le maximum sur tout le spectre affichait
    "max Re = +0,00e+00" pour une marge reelle de -0,880 : le nombre qu'on
    lit en premier ne distinguait plus un reseau juste stable d'un reseau
    confortable.
    """
    d = analyse.modes(charger(nom), f_max=60.0)
    assert d["jauge"] == jauge_attendue
    assert d["max_re"] < -1e-3, f"max Re = {d['max_re']}, la jauge domine encore"
    assert d["stable"]
    assert all(abs(complex(m["re"], m["im"])) > analyse.JAUGE for m in d["modes"]), (
        "un mode de jauge figure encore dans la table"
    )


def test_modes_voit_toujours_une_vraie_instabilite():
    """Ecarter la jauge ne doit pas ecarter le danger.

    Le meme reseau a 2 % et a 12 % de statisme : l'un est sain, l'autre
    franchement instable, et la difference doit rester visible.
    """
    t = charger("smib_gfm").copie()
    t.GFM_PARAMS["mp"] = 0.02
    sain = analyse.modes(t, f_max=60.0)
    t.GFM_PARAMS["mp"] = 0.12
    fatal = analyse.modes(t, f_max=60.0)
    assert sain["stable"] and sain["max_re"] < 0
    assert not fatal["stable"] and fatal["max_re"] > 1.0


# ----------------------------------------------------------------------
# les raccordements continus, vus par le validateur
# ----------------------------------------------------------------------
def _hvdc():
    """Une copie de ``hvdc2t``, a abimer de differentes facons."""
    return charger("hvdc2t").copie()


def test_valider_refuse_un_vsc_sur_le_slack():
    """Un convertisseur tient exactement DEUX grandeurs.

    Sur la reference d'angle, les deux puissances sont libres : il n'y a
    rien a tenir. La regle surprend en dessinant, parce que le premier
    equipement pose devient justement la reference -- d'ou l'interet de la
    dire pendant qu'on dessine et non au premier load flow.
    """
    t = _hvdc()
    ac = int(t.T_VSCs[0][1])
    lignes = [list(r) for r in t.Y_DER]
    for r in lignes:
        r[1] = 1 if int(r[0]) == ac else (2 if int(r[1]) == 1 else r[1])
    t.Y_DER = [tuple(r) for r in lignes]
    msgs = [x["message"] for x in analyse.valider(t) if x["gravite"] == "erreur"]
    assert any("slack" in m and "deux grandeurs" in m for m in msgs), msgs


def test_valider_refuse_un_vsc_sur_un_bus_sans_convertisseur():
    t = _hvdc()
    lignes = [list(r) for r in t.T_VSCs]
    lignes[0][1] = 9999
    t.T_VSCs = [tuple(r) for r in lignes]
    msgs = [x["message"] for x in analyse.valider(t) if x["gravite"] == "erreur"]
    assert any("ni GFM3" in m for m in msgs), msgs


def test_valider_refuse_un_bus_continu_en_l_air():
    t = _hvdc()
    t.T_BusDC = list(t.T_BusDC) + [(999, 2, 40e3, 1.0, 0.0)]
    msgs = [x["message"] for x in analyse.valider(t) if x["gravite"] == "erreur"]
    assert any("999" in m and "isole" in m for m in msgs), msgs


def test_positions_placent_aussi_les_noeuds_continus():
    """Un noeud sans position atterrit au centre du canevas.

    Sur une liaison a deux terminaux, les deux s'y empilent et le schema
    devient illisible. Les cables continus et les raccordements de T_VSCs
    comptent donc comme aretes du placement -- les seconds ne sont pas des
    branches du modele, mais ce sont eux qui disent quel noeud continu va
    pres de quel poste.
    """
    for nom in ("hvdc2t", "mtdc3t"):
        t = charger(nom)
        pos = analyse.depuis_catalogue(nom)["positions"]
        attendus = {str(int(r[0])) for r in t.Y_DER} | {
            str(int(r[0])) for r in t.T_BusDC
        }
        assert set(pos) == attendus, f"{nom} : {attendus ^ set(pos)}"
        distinctes = {tuple(v) for v in pos.values()}
        assert len(distinctes) == len(pos), f"{nom} : des noeuds superposes"


# ----------------------------------------------------------------------
# la cellule de donnees : ce qu'on dessine se relit en code
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nom", ["smib", "wscc9", "hvdc2t", "mtdc3t"])
def test_la_cellule_exportee_redonne_le_meme_reseau(nom):
    """Exporter puis relire doit rendre le reseau A L'IDENTIQUE.

    C'est le contrat de l'atelier : ce qu'on dessine EST les tables, et la
    cellule exportee est collable telle quelle dans un notebook. Si
    l'aller-retour perdait une colonne -- ou, pire, arrondissait une valeur
    -- l'interface et les notebooks decriraient deux reseaux differents en
    pretendant decrire le meme.

    On compare les tables ET la matrice ``A`` : la premiere comparaison dit
    ou ca casse, la seconde dit que ca n'a pas casse.
    """
    t0 = charger(nom)
    g = {"np": np, "math": __import__("math")}
    exec(compile(analyse.exporter(t0), f"{nom}:export", "exec"), g)  # noqa: S102
    t1 = Tables(
        **{k: list(g.get(k, [])) for k in TABLES},
        **{
            k: dict(g.get(k, {}))
            for k in ("SM_PARAMS", "GFM_PARAMS", "GFL_PARAMS", "PARAMS_PAR_BUS")
        },
        f0=float(g.get("f0", 50.0)),
        Sb=float(g.get("Sb", 100e6)),
        nom=nom,
        **({"NOEUDS": g["NOEUDS"]} if "NOEUDS" in g else {}),
        **({"R_NOEUDS": float(g["R_NOEUDS"])} if "R_NOEUDS" in g else {}),
    )
    assert (t1.NOEUDS, t1.R_NOEUDS) == (t0.NOEUDS, t0.R_NOEUDS)
    for k in TABLES:
        a = [tuple(r) for r in getattr(t0, k)]
        b = [tuple(r) for r in getattr(t1, k)]
        assert a == b, f"{nom} : la table {k} differe apres aller-retour"

    A0 = analyse.construire(t0).etat_stationnaire().A
    A1 = analyse.construire(t1).etat_stationnaire().A
    assert A0.shape == A1.shape
    assert np.abs(A0 - A1).max() == 0.0, f"{nom} : la matrice A a bouge"
