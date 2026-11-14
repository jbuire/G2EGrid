"""Les routes des zones, des fichiers et des exports.

Deux niveaux, parce qu'ils ne cassent pas pour les memes raisons.

La LOGIQUE vit dans ``g2egrid_api.analyse``, qui ne sait pas ce qu'est HTTP :
c'est la que se verifie qu'un reseau ouvert peut etre renvoye tel quel pour
etre aplati ou enregistre, et qu'il revient identique -- tables et matrice A.

Le TRANSPORT vit dans ``g2egrid_api.app``. Ce qu'il ajoute, c'est une regle :
un reseau mal construit est une reponse a rendre a l'utilisateur (422, avec
le message du modele), jamais une panne du serveur (500). Elle se verifie en
appelant les fonctions de route, qui restent des fonctions Python ordinaires
sous leur decorateur.

Aucun client HTTP n'est monte ici : ``fastapi.testclient`` demande le paquet
``httpx``, que l'environnement n'installe pas. Ce que ces tests ne couvrent
donc pas, c'est la serialisation des corps de requete par pydantic.
"""

from __future__ import annotations

import json
import os
import pathlib

import numpy as np
import pytest

fastapi = pytest.importorskip("fastapi")

import g2egrid.catalogue as cat  # noqa: E402
from g2egrid_api import analyse  # noqa: E402
from g2egrid_api import app as routes  # noqa: E402
from g2egrid_api.schemas import CorpsFichier, CorpsReseau  # noqa: E402


@pytest.fixture
def dossier_temporaire(tmp_path, monkeypatch):
    """Les enregistrements vont dans un dossier jetable.

    ``g2egrid/reseaux/`` porte les exemples natifs, suivis par git : un test
    qui y ecrirait laisserait des fichiers derriere lui, ou en ecraserait.
    """
    monkeypatch.setattr(cat, "DOSSIER_RESEAUX", tmp_path)
    return tmp_path


@pytest.fixture(scope="module")
def _lecture_zones_une_fois():
    return analyse.depuis_catalogue("wscc3sm_zones")


@pytest.fixture
def lecture_zones(_lecture_zones_une_fois):
    """La lecture de ``wscc3sm_zones``, faite AVANT toute redirection.

    Une fixture de module est instanciee avant les fixtures de fonction, donc
    avant que ``dossier_temporaire`` ne detourne ``DOSSIER_RESEAUX``. Lire
    apres le detournement cherchait le fichier natif dans le dossier jetable,
    ne l'y trouvait pas, et retombait sur le notebook -- dont la cellule de
    donnees ne s'execute pas seule. Une copie par test : certains la modifient.
    """
    import copy

    return copy.deepcopy(_lecture_zones_une_fois)


def _A(t):
    return analyse.construire(t).etat_stationnaire().A


# ----------------------------------------------------------------------
# lecture : les deux formes du reseau
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nom", ["smib", "hvdc2t", "wscc3sm_zones"])
def test_la_lecture_garde_le_contrat_a_plat(nom):
    """``tables`` et ``positions`` restent le reseau APLATI, partout.

    C'est ce que lisent le schema, l'onglet Tables et toutes les routes de
    calcul : les zones s'ajoutent a la reponse, elles ne la remplacent pas.
    """
    d = analyse.depuis_catalogue(nom)
    assert {"tables", "positions", "transport", "zones", "instances"} <= set(d)
    t = cat.charger(nom)
    assert len(d["tables"]["Y_DER"]) == len(t.Y_DER)
    attendus = {str(int(r[0])) for r in t.Y_DER} | {str(int(r[0])) for r in t.T_BusDC}
    assert set(d["positions"]) == attendus


def test_une_zone_se_lit_une_fois():
    """Trois instances d'une meme distribution, UNE description.

    Rendre seulement le reseau a plat revenait a montrer trois copies d'une
    zone -- c'est-a-dire a perdre exactement ce que les zones apportent.
    """
    d = analyse.depuis_catalogue("wscc3sm_zones")
    assert len(d["zones"]) == 1
    assert len(d["instances"]) == 3
    assert len({i["zone"] for i in d["instances"]}) == 1
    # chaque instance dit quels numeros elle occupe a plat, pour que
    # l'interface relie un bus du schema a sa zone
    for i in d["instances"]:
        assert i["bus"], "une instance sans correspondance vers le reseau a plat"


# ----------------------------------------------------------------------
# aller-retour : ce qu'on ouvre, on doit pouvoir le renvoyer
# ----------------------------------------------------------------------
def test_la_lecture_se_renvoie_a_aplatir():
    """L'interface ouvre un reseau, l'edite, le renvoie : ca doit marcher.

    La reponse de lecture et le fichier n'ont pas la meme forme, et renvoyer
    la premiere a ``aplatir`` echouait sur "format 'None' inconnu" -- on
    pouvait ouvrir un reseau a zones, pas l'aplatir ni le sauver.
    """
    d = analyse.depuis_catalogue("wscc3sm_zones")
    plat = analyse.aplatir_reseau(d)
    tp = plat.get("tables", plat)
    ref = cat.charger("wscc3sm_cigre3")
    for k in analyse.TABLES:
        assert [list(r) for r in tp.get(k, [])] == [list(r) for r in getattr(ref, k)], k


def test_enregistrer_puis_relire_redonne_la_meme_matrice(
    dossier_temporaire, lecture_zones
):
    """Le fichier est la verite : relu, il doit donner EXACTEMENT le meme A."""
    analyse.enregistrer_fichier(lecture_zones, "essai_zones")
    assert (dossier_temporaire / "sauvegardes" / "essai_zones.g2egrid.json").exists()

    relu = cat.charger("essai_zones")
    ref = cat.charger("wscc3sm_cigre3")  # notebook : non touche par la redirection
    A0, A1 = _A(ref), _A(relu)
    assert A0.shape == A1.shape
    assert np.abs(A0 - A1).max() == 0.0

    # et la structure a survecu : une zone, trois instances
    d2 = analyse.depuis_catalogue("essai_zones")
    assert len(d2["zones"]) == 1 and len(d2["instances"]) == 3


def test_un_reseau_dessine_a_la_main_s_enregistre(dossier_temporaire):
    """Un reseau du constructeur n'a que ``tables`` et ``positions``.

    Il n'a ni transport ni zone : son transport est le reseau entier. Il doit
    pouvoir s'enregistrer quand meme -- c'est le cas le plus frequent.
    """
    t = cat.charger("smib")
    d = {
        "nom": "dessine",
        "tables": analyse.tables_vers_json(t),
        "positions": analyse.depuis_catalogue("smib")["positions"],
    }
    analyse.enregistrer_fichier(d, "dessine")
    relu = cat.charger("dessine")
    assert [list(r) for r in relu.Y_DER] == [list(r) for r in t.Y_DER]
    assert np.abs(_A(t) - _A(relu)).max() == 0.0


def test_fichiers_liste_ce_qui_est_enregistre(dossier_temporaire):
    assert analyse.fichiers() == []
    analyse.enregistrer_fichier(
        {"tables": analyse.tables_vers_json(cat.charger("smib"))}, "un"
    )
    noms = [f["nom"] for f in analyse.fichiers()]
    assert noms == ["un"]


def test_les_exemples_natifs_ne_bougent_pas():
    """Aucun test de ce fichier n'ecrit dans ``g2egrid/reseaux/``."""
    natifs = sorted(p.name for p in pathlib.Path(cat.DOSSIER_RESEAUX).iterdir())
    assert "smib.g2egrid.json" in natifs
    assert not any(n.startswith(("essai", "dessine", "un.")) for n in natifs)


@pytest.mark.parametrize(
    "nom", ["../evasion", "a/b", "x" * 80, "con fichier", ".cache"]
)
def test_un_nom_de_fichier_dangereux_est_refuse(dossier_temporaire, nom):
    """Le nom vient du navigateur : il ne doit pas pouvoir sortir du dossier."""
    with pytest.raises(ValueError):
        analyse.enregistrer_fichier(
            {"tables": analyse.tables_vers_json(cat.charger("smib"))}, nom
        )
    assert list(dossier_temporaire.iterdir()) == [], "un fichier a ete ecrit quand meme"


def test_un_nom_vide_prend_celui_du_reseau(dossier_temporaire):
    """Ne pas nommer n'est pas une faute : le reseau porte deja un nom."""
    t = cat.charger("smib")
    analyse.enregistrer_fichier(
        {"nom": "mon_smib", "tables": analyse.tables_vers_json(t)}, ""
    )
    assert [f.name for f in (dossier_temporaire / "sauvegardes").iterdir()] == [
        "mon_smib.g2egrid.json"
    ]


# ----------------------------------------------------------------------
# exports
# ----------------------------------------------------------------------
@pytest.mark.parametrize("quoi", ["donnees", "notebook", "matlab"])
def test_les_trois_exports_rendent_un_contenu(quoi):
    d = analyse.exporter_vers(cat.charger("wscc9"), quoi)
    assert d["quoi"] == quoi
    assert d["contenu"], "export vide"
    assert d["fichier"]
    assert isinstance(d["avertissements"], list)


def test_l_export_matlab_dit_ce_qu_il_perd():
    """Les shunts partent en commentaire : l'export doit le DIRE.

    Un ``.m`` ampute reste executable et decrit pourtant un autre reseau :
    avaler l'avertissement donnerait un export qui a l'air complet.
    """
    t = cat.charger("wscc3sm_cigre3")
    assert t.T_Shunts, "le cas choisi n'a plus de shunt a perdre"
    d = analyse.exporter_vers(t, "matlab")
    assert d["avertissements"], "des shunts ont ete perdus sans le dire"


def test_l_export_matlab_refuse_le_continu():
    with pytest.raises(ValueError, match="continu"):
        analyse.exporter_vers(cat.charger("hvdc2t"), "matlab")


def test_un_format_d_export_inconnu_est_refuse():
    with pytest.raises(ValueError, match="inconnu"):
        analyse.exporter_vers(cat.charger("smib"), "pdf")


# ----------------------------------------------------------------------
# transport : une erreur de modele est une 422, jamais une 500
# ----------------------------------------------------------------------
def _statut(fn, *a):
    """Le code HTTP que rendrait la route, sans monter de serveur."""
    try:
        fn(*a)
    except fastapi.HTTPException as err:
        return err.status_code, str(err.detail)
    return 200, ""


def test_les_routes_rendent_200_sur_un_reseau_correct(
    dossier_temporaire, lecture_zones
):
    corps = CorpsReseau(tables=analyse.tables_vers_json(cat.charger("wscc9")))
    assert _statut(routes.fichiers)[0] == 200
    assert _statut(routes.exporter_vers, "donnees", corps)[0] == 200
    assert _statut(routes.zones_aplatir, CorpsFichier(reseau=lecture_zones))[0] == 200
    assert (
        _statut(
            routes.enregistrer_fichier, CorpsFichier(reseau=lecture_zones, nom="ok")
        )[0]
        == 200
    )


def test_le_continu_vers_matlab_est_une_422():
    corps = CorpsReseau(tables=analyse.tables_vers_json(cat.charger("hvdc2t")))
    code, detail = _statut(routes.exporter_vers, "matlab", corps)
    assert code == 422
    assert "continu" in detail


def test_un_reseau_a_zones_incoherent_est_une_422(lecture_zones):
    """Les erreurs d'aplatissement sont redigees pour un utilisateur.

    Un prefixe qui chevauche le transport, un bus amont absent : ce sont des
    fautes de saisie, et elles doivent revenir avec leur message.
    """
    d = lecture_zones
    d["instances"][0]["bus_amont"] = 99999
    code, detail = _statut(routes.zones_aplatir, CorpsFichier(reseau=d))
    assert code == 422, detail
    assert "99999" in detail


def test_un_corps_malforme_est_une_422_et_non_une_500(lecture_zones):
    """Une instance sans ``zone`` ne doit pas faire tomber le serveur."""
    d = lecture_zones
    del d["instances"][0]["zone"]
    code, _ = _statut(routes.zones_aplatir, CorpsFichier(reseau=d))
    assert code == 422


# ----------------------------------------------------------------------
# ouvrir et telecharger depuis le navigateur
# ----------------------------------------------------------------------
def test_le_telechargement_est_un_vrai_fichier(dossier_temporaire, lecture_zones):
    """Ce qu'on telecharge doit pouvoir etre depose dans ``g2egrid/reseaux/``.

    Telecharger la forme de lecture donnait un fichier que le catalogue ne
    relisait pas. On ecrit ici le contenu converti a la main dans le dossier,
    et le catalogue doit le relire -- et redonner la meme matrice A.
    """
    import json

    contenu = analyse.convertir_fichier(lecture_zones)
    assert contenu["format"] == "g2egrid"
    # depose a la main dans le dossier des exemples, comme on le ferait
    (dossier_temporaire / "depose.g2egrid.json").write_text(
        json.dumps(contenu), encoding="utf-8"
    )
    relu = cat.charger("depose")
    assert np.abs(_A(cat.charger("wscc3sm_cigre3")) - _A(relu)).max() == 0.0


@pytest.mark.parametrize("forme", ["natif", "lecture", "a_plat"])
def test_ouvrir_accepte_les_trois_formes(forme, lecture_zones):
    """Le navigateur ne connait pas le format : le serveur les accepte tous."""
    if forme == "natif":
        d = analyse.convertir_fichier(lecture_zones)
    elif forme == "lecture":
        d = lecture_zones
    else:
        t = cat.charger("wscc3sm_cigre3")
        d = {"tables": analyse.tables_vers_json(t), "positions": {}}
    lu = analyse.lire_fichier(d)
    assert {"tables", "positions", "transport", "zones", "instances"} <= set(lu)
    assert len(lu["tables"]["Y_DER"]) == 54
    # un reseau a plat n'a pas de zone, et c'est une reponse, pas une erreur
    assert (len(lu["zones"]), len(lu["instances"])) == (
        (0, 0) if forme == "a_plat" else (1, 3)
    )


def test_un_telechargement_incoherent_est_refuse_avant_d_ecrire(lecture_zones):
    """Les garde-fous de l'aplatissement jouent avant qu'un fichier parte."""
    lecture_zones["instances"][0]["bus_amont"] = 99999
    code, detail = _statut(routes.convertir_fichier, CorpsFichier(reseau=lecture_zones))
    assert code == 422 and "99999" in detail


def test_l_ancien_fichier_du_constructeur_s_ouvre():
    """L'ancien bouton Enregistrer ecrivait un faux fichier natif.

    Il portait {"format": "g2egrid", "version": 1} mais rangeait ses tables a
    plat, sans ``transport``. Pris pour natif, il s'ouvrait en un reseau VIDE
    -- 0 bus, et aucune erreur. Les utilisateurs ont de ces fichiers : ils
    doivent se rouvrir, et a l'identique.
    """
    t = cat.charger("wscc9")
    vieux = {
        "format": "g2egrid",
        "version": 1,
        "nom": "vieux",
        "tables": analyse.tables_vers_json(t),
        "positions": {},
    }
    lu = analyse.lire_fichier(vieux)
    assert len(lu["tables"]["Y_DER"]) == len(t.Y_DER) == 9
    assert np.abs(_A(t) - _A(analyse.tables_depuis_json(lu["tables"]))).max() == 0.0


def test_une_sauvegarde_ne_rejoint_pas_les_exemples(dossier_temporaire):
    """Enregistrer depuis l'interface ne doit pas casser la suite de tests.

    Les sauvegardes allaient dans le dossier des exemples ; la
    non-regression bit a bit, qui enumere les reseaux alternatifs du
    catalogue et exige une reference pour chacun, tombait alors. Une
    sauvegarde va dans son sous-dossier, apparait tout de suite dans la
    liste du catalogue -- sans redemarrer -- mais n'entre pas dans RESEAUX,
    la liste figee des exemples que les tests couvrent.
    """
    analyse.enregistrer_fichier(
        {"nom": "perso", "tables": analyse.tables_vers_json(cat.charger("smib"))},
        "perso",
    )
    assert (dossier_temporaire / "sauvegardes" / "perso.g2egrid.json").exists()
    assert not (dossier_temporaire / "perso.g2egrid.json").exists()
    assert "perso" not in cat.RESEAUX
    assert "perso" in [f["nom"] for f in analyse.catalogue()]
    assert [f["origine"] for f in analyse.fichiers() if f["nom"] == "perso"] == [
        "sauvegarde"
    ]


def test_une_sauvegarde_ne_masque_pas_un_exemple(dossier_temporaire):
    """Un exemple et une sauvegarde de meme nom : l'exemple gagne."""
    import json

    exemple = analyse.convertir_fichier(
        {"nom": "double", "tables": analyse.tables_vers_json(cat.charger("smib"))}
    )
    (dossier_temporaire / "double.g2egrid.json").write_text(
        json.dumps(exemple), encoding="utf-8"
    )
    analyse.enregistrer_fichier(
        {"nom": "double", "tables": analyse.tables_vers_json(cat.charger("wscc9"))},
        "double",
    )
    assert len(cat.charger("double").Y_DER) == 3  # smib, pas wscc9


def test_definir_zone_rend_un_reseau_en_zones():
    d = analyse.depuis_catalogue("ieee14")
    r = analyse.definir_zone(d, [8], "z")
    assert list(r["zones"]) == ["z"]
    assert len(r["instances"]) == 1


def test_robustesse_rend_des_modes_classes():
    t = cat.charger("wscc9")
    d = analyse.robustesse(t, zeta_max=30, puissance=True, deconnexion=False)
    assert d["n_variantes"] > 0 and d["modes"]
    assert all(m["score_deconnexion"] == 0 for m in d["modes"])
    scores = [m["score_puissance"] for m in d["modes"]]
    assert scores == sorted(scores, reverse=True)


def test_convention_gfm_batterie_gfl_pv_signalee():
    t = cat.charger("wscc9_cigre")
    assert not [p for p in analyse.valider(t) if "GF" in p["message"]]
    y = [list(r) for r in t.Y_DER]
    for r in y:
        if int(r[0]) == 2:  # le GFM
            r[1] = 3
        if int(r[0]) == 3:  # un GFL
            r[1] = 2
    msgs = [p["message"] for p in analyse.valider(t.copie(Y_DER=y))]
    assert any(m.startswith("GFM du bus 2 en PQ") for m in msgs)
    assert any(m.startswith("GFL du bus 3 en PV") for m in msgs)


def test_importer_refuse_un_nom_qui_sort_du_catalogue(tmp_path):
    """Un nom de catalogue est un nom de fichier nu, jamais un chemin.

    Non-regression : ``{"catalogue": "../../x"}`` faisait executer la
    cellule de donnees de n'importe quel ``.ipynb`` du disque. Le notebook
    piege ecrirait un temoin : il ne doit pas exister apres l'appel.
    """
    temoin = tmp_path / "temoin"
    cellule = f"open({str(temoin)!r}, 'w').close()\nY_DER = []\n"
    nb = {"cells": [{"cell_type": "code", "source": [cellule]}]}
    (tmp_path / "piege.ipynb").write_text(json.dumps(nb), encoding="utf-8")
    nom = os.path.relpath(tmp_path / "piege", cat.DOSSIER)
    d = analyse.depuis_catalogue("wscc9")
    code, msg = _statut(
        routes.reseau_importer,
        routes.CorpsImport(reseau={"tables": d["tables"]}, source={"catalogue": nom}),
    )
    assert code == 422 and "nom de reseau invalide" in msg
    assert not temoin.exists()


def test_importer_un_reseau_du_catalogue():
    d = analyse.depuis_catalogue("wscc9")
    o = analyse.importer(
        {"tables": d["tables"], "positions": d["positions"]},
        {"catalogue": "ieee14"},
        {"mode": "transfo", "bus_importe": 1, "bus_existant": 5},
    )
    assert len(o["tables"]["Y_DER"]) == 23
    assert o["import"]["bus"]["1"] == 101
    assert set(o["positions"]) == {str(int(r[0])) for r in o["tables"]["Y_DER"]}
    assert not [
        p for p in analyse.valider(analyse._tables_de(o)) if p["gravite"] == "erreur"
    ]
    code, msg = _statut(
        routes.reseau_importer,
        routes.CorpsImport(
            reseau={"tables": d["tables"]},
            source={"catalogue": "wscc9"},
            raccord={"mode": "fusion", "bus_importe": 1, "bus_existant": 1},
        ),
    )
    assert code == 422 and "deux unites" in msg
