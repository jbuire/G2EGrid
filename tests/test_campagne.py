"""La campagne (``etudes/campagne``) : le spectre reduit que garde la
comparaison, le classement des modes instables d'un modele reduit, et les
trois niveaux de fidelite de l'analyse (detection, categorie, mecanisme)."""

from __future__ import annotations

import numpy as np
import pytest
from _etudes import dossier_etude
from g2egrid.catalogue import charger
from g2egrid.robustesse import _spectre
from g2egrid.simplification import comparer, preregle

ICI = dossier_etude("campagne")

import analyse
import campagne


@pytest.fixture(scope="module")
def wscc():
    return charger("wscc9")


def test_comparer_garde_le_spectre_reduit(wscc):
    """Chaque mode retrouve ou nouveau designe SON mode du spectre reduit,
    vecteurs propres compris quand on les demande."""
    c = comparer(wscc, preregle(wscc, "flux_stator"), vecteurs=True)
    red = c.reduit
    assert red.V is not None and red.V.shape == (len(red.etats), red.lam.size)
    vus = [m for m in c.modes if m.statut in ("retrouve", "nouveau")]
    assert vus
    for m in vus:
        assert red.lam[m.j_reduit] == pytest.approx(m.lam_reduit)
    # sans la demande, pas de vecteurs : le spectre reste leger
    assert comparer(wscc, preregle(wscc, "flux_stator")).reduit.V is None


def test_modes_instables_d_un_modele_reduit_classes():
    """Un modele reduit instable (tirage 3 de WSCC9, stator et amortisseurs
    algebriques : il invente des resonances) : ses modes instables sont
    classes, du plus au moins instable, et rattaches au mode du complet dont
    ils sont l'image."""
    t = campagne.construire(campagne.tirer(3, "wscc9"))
    _, complet, inerties = campagne.modele_complet(t)
    retenus = campagne._modes_retenus(complet)
    simp = preregle(t, "flux_stator", "amortisseurs")
    c = comparer(t, simp, complet=complet, vecteurs=True)
    inst = campagne.instables_reduits(c, complet, retenus, inerties)
    assert len(inst) == int(np.sum(c.reduit.lam.real > 0)) > 0
    assert [e["re"] for e in inst] == sorted((e["re"] for e in inst), reverse=True)
    for e in inst:
        assert e["categorie"] in analyse.CATEGORIES
        assert e["categorie_complet"] in (None, *analyse.CATEGORIES)
        assert e["complet"] is None or 0 <= e["complet"] < len(retenus)
        assert abs(sum(e["groupes"].values()) - 1) < 0.1


def _tirage(modes, reduit, max_re=None):
    """Un tirage synthetique : ``modes`` du complet (categorie, groupe
    dominant, instable), ``reduit`` les modes instables du modele reduit
    (categorie, groupe dominant, indice du mode du complet ou None)."""
    ms = [
        {"categorie": c, "groupes": {g: 0.8}, "instable": i, "re": 1.0 if i else -1.0}
        for c, g, i in modes
    ]
    inst = [
        {"categorie": c, "groupes": {g: 0.7}, "complet": k, "re": 1.0}
        for c, g, k in reduit
    ]
    return {
        "systeme": "wscc9",
        "i": 0,
        "statut": "ok",
        "max_re": max_re
        if max_re is not None
        else (1.0 if any(m["instable"] for m in ms) else -1.0),
        "modes": ms,
        "reduits": {"stator": {"max_re": 1.0 if inst else -1.0, "instables": inst}},
    }


def test_trois_niveaux_de_fidelite():
    R = [
        # angle garde, meme mecanisme : vrai positif a tous les niveaux
        _tirage([("angle", "SM_rotor", True)], [("angle", "SM_rotor", 0)]),
        # angle garde instable mais vu en tension : detecte, mal classe
        _tirage([("angle", "SM_rotor", True)], [("tension", "SM_excitation", 0)]),
        # interaction rapide perdue
        _tirage([("convertisseur_rapide", "GFM_courant", True)], []),
        # tirage stable ou le reduit invente une resonance
        _tirage([("angle", "SM_rotor", False)], [("resonance", "reseau", None)]),
        # tirage stable, reduit stable
        _tirage([("angle", "SM_rotor", False)], []),
        # meme categorie, autre groupe dominant : categorie gardee, pas le mecanisme
        _tirage([("angle", "SM_rotor", True)], [("angle", "SM_pss", 0)]),
    ]
    f = analyse.fidelite(R)["stator"]
    d = f["detection"]
    assert (d["tp"], d["fn"], d["fp"], d["tn"]) == (3, 1, 1, 1)
    a = f["categorie"]["angle"]
    assert (a["tp"], a["fn"], a["fp"]) == (2, 1, 0)
    assert f["categorie"]["tension"]["fp"] == 1
    assert f["categorie"]["resonance"]["fp"] == 1
    assert f["categorie"]["convertisseur_rapide"]["rappel"] == 0.0
    m = f["mecanisme"]["angle"]
    assert (m["tp"], m["fn"], m["fp"]) == (1, 2, 1)
    assert f["mecanisme"]["angle"]["rappel"] == pytest.approx(1 / 3)
    conf = f["confusion"]
    assert conf["angle>angle"] == 2 and conf["angle>tension"] == 1
    assert conf["convertisseur_rapide>stable"] == 1
    assert conf["stable>resonance"] == 1 and conf["stable>stable"] == 1
    # des trois modes gardes instables, deux gardent leur categorie, un son
    # mecanisme
    g = f["gardes"]
    assert g["n"] == 3
    assert g["categorie"] == pytest.approx(2 / 3)
    assert g["mecanisme"] == pytest.approx(1 / 3)
    # tirages dits instables a raison : le deuxieme annonce une tension que
    # le complet n'a pas
    assert g["tirages"] == 3
    assert g["tirages_memes"] == pytest.approx(2 / 3)
    assert g["tirages_inclus"] == pytest.approx(2 / 3)


def test_un_tirage_non_couvert_n_est_pas_juge():
    """Un modele reduit instable dont les modes ne sont pas classes (campagne
    anterieure au complement) : le tirage est hors de l'analyse des
    categories."""
    r = _tirage([("angle", "SM_rotor", True)], [("angle", "SM_rotor", 0)])
    del r["reduits"]["stator"]["instables"]
    assert not analyse.classes(r)
    assert analyse.fidelite([r]) == {}
    r["complement"] = True
    assert analyse.classes(r)


def test_lent_ou_rapide_par_les_participations():
    modes = [
        {"categorie": "convertisseur_lent", "groupes": {"GFL_pll": 0.6, "reseau": 0.3}},
        {
            "categorie": "convertisseur_lent",
            "groupes": {"GFL_courant": 0.5, "reseau": 0.4},
        },
        {
            "categorie": "convertisseur_rapide",
            "groupes": {"GFM_courant": 0.5, "reseau": 0.5},
        },
        {"categorie": "angle", "groupes": {"SM_rotor": 1.0}},
    ]
    x = analyse.lent_ou_rapide(modes)
    assert x["n"] == 3
    assert x["accord"] == pytest.approx(2 / 3)
    assert x["table"]["convertisseur_lent>rapide"] == 1


def test_spectre_sans_vecteurs_par_defaut(wscc):
    assert _spectre(wscc).V is None
    s = _spectre(wscc, vecteurs=True)
    assert s.V.shape == (len(s.etats), s.lam.size)


def test_parts_retirees_egales_a_celles_de_la_comparaison():
    """La participation negligee calculee sur le seul modele complet est celle
    que mesure la comparaison au modele reduit, etats elimines par reduction
    d'indice compris (PI de PLL instantane : v_q du condensateur du GFL)."""
    for i in range(10):
        t = campagne.construire(campagne.tirer(i, "three_machine_ib"))
        try:
            net, complet, _ = campagne.modele_complet(t)
        except Exception:
            continue
        retenus = campagne._modes_retenus(complet)
        gfl = any(type(el).__name__.startswith("GFL") for el in net.all_elements)
        if retenus and gfl:
            break
    else:
        pytest.fail("aucun tirage a modes retenus et a GFL")
    parts = campagne.parts_retirees(t, net, complet, retenus)
    from g2egrid.simplification import fusionner

    for nom, pre in campagne.REDUITS.items():
        simp = fusionner(*(preregle(t, x) for x in pre))
        c = comparer(t, simp, complet=complet)
        for a, k in enumerate(retenus):
            assert parts[a][nom] == pytest.approx(c.modes[k].part_retiree, abs=1e-3)


def test_melange_et_localisation():
    assert analyse.melange({"angle": 1.0}) == pytest.approx(0.0)
    # deux familles a parts egales : ln 2 / ln 5
    assert analyse.melange({"angle": 0.5, "tension": 0.5}) == pytest.approx(
        np.log(2) / np.log(5)
    )
    assert analyse.melange({k: 0.2 for k in "abcde"}) == pytest.approx(1.0)
    assert analyse.localisation([("SM2", 0.7), ("SM1", 0.2)]) == "local"
    assert analyse.localisation([("SM2", 0.3), ("SM1", 0.2), ("L1", 0.1)]) == "regional"
    assert analyse.localisation([("L1", 0.2), ("L2", 0.1), ("L3", 0.1)]) == "etendu"


def test_retenue_selon_la_part_negligee():
    """Un mode garde par un modele reduit qui en neglige peu, perdu par un
    autre qui en neglige beaucoup : chacun dans sa classe de part."""
    r = {
        "statut": "ok",
        "modes": [
            {
                "categorie": "angle",
                "instable": True,
                "parts_retirees": {"stator": 0.002, "excitation": 0.6, "bus_dc": 0.0},
            }
        ],
        "reduits": {
            "stator": {"modes": [{"statut": "retrouve", "re": 0.5}]},
            "excitation": {"modes": [{"statut": "retrouve", "re": -0.2}]},
            "bus_dc": {"modes": [{"statut": "retrouve", "re": 0.5}]},
        },
    }
    x = analyse.retenue_selon_part([r])
    assert x["toutes"][0] == {"n": 1, "gardes": 1}  # < 1 % : garde
    assert x["toutes"][-1] == {"n": 1, "gardes": 0}  # >= 50 % : perdu
    # le bus continu ne neglige rien : hors du compte
    assert sum(v["n"] for v in x["toutes"]) == 2


def test_sorties_dans_les_deux_langues(tmp_path):
    """Les tables du papier s'ecrivent en francais et en anglais a partir du
    seul bilan, avec les memes macros."""
    import sorties

    R = [
        _tirage([("angle", "SM_rotor", True)], [("angle", "SM_rotor", 0)]),
        _tirage([("angle", "SM_rotor", False)], [("resonance", "reseau", None)]),
    ]
    for r in R:
        for m in r["modes"]:
            m.update(f=1.2, nature="oscillant", parts={"angle": 0.8, "tension": 0.2})
            m["elements"] = [("SM1", 0.7)]
            m["seconde"] = "tension"
    cat = {
        "n": 2,
        "ok": 2,
        "stables": 1,
        "instables": 1,
        "sans_point": 0,
        "n_modes": 1,
        "types": {
            t: {"tirages": 0, "modes": 0, "f_min": None, "f_max": None}
            for t in sorties.TYPES
        },
        "causes": {},
    }
    b = {
        "catalogue": cat,
        "detection": {},
        "par_reseau": {},
        "fidelite": analyse.fidelite(R),
        "lent_rapide": analyse.lent_ou_rapide([]),
        "couverture": analyse.couverture(R),
        "multiaxe": analyse.multiaxe(
            [m for r in R for m in r["modes"] if m["instable"]]
        ),
        "retenue": analyse.retenue_selon_part(R),
    }
    for langue in sorties.LANGUES:
        f = tmp_path / f"resultats_{langue}.tex"
        sorties.latex(b, f, langue)
        texte = f.read_text(encoding="utf-8")
        for macro in (
            "\\TableFidelite",
            "\\TableSynthese",
            "\\TableMultiaxe",
            "fid:stator:det_f1",
        ):
            assert macro in texte
    assert "rappel" in (tmp_path / "resultats_fr.tex").read_text(encoding="utf-8")
    assert "Reduced model" in (tmp_path / "resultats_en.tex").read_text(
        encoding="utf-8"
    )


def test_le_classement_ne_depend_pas_du_traitement_des_noeuds():
    """FR-044 : tirage 33 d'ieee14, une resonance a 60,5 Hz. Aux noeuds
    ideaux, la reduction d'indice elimine 14 courants de ligne ; lue sur les
    seuls etats gardes, la part du reseau tombait de 0,80 a 0,35 et le mode
    passait « convertisseur rapide » a valeur propre identique. Dans la base
    complete, c'est la limite R -> infini du modele a resistance de
    causalite : meme categorie, meme part des commandes."""
    p = campagne.tirer(33, "ieee14")
    vus = []
    for reglage in (
        {"NOEUDS": "ideal"},
        {"NOEUDS": "resistance", "R_NOEUDS": 1e4},
    ):
        t = campagne.construire(p, tables=charger("ieee14").copie(**reglage))
        _, complet, inerties = campagne.modele_complet(t)
        k = int(np.argmax(complet.lam.real))
        assert abs(complet.lam[k].imag) / (2 * np.pi) == pytest.approx(60.51, abs=0.01)
        vus.append(campagne.classer_mode(complet, k, inerties))
    ideal, resistif = vus
    assert ideal["categorie"] == resistif["categorie"] == "resonance"
    assert ideal["part_commande"] == pytest.approx(resistif["part_commande"], abs=1e-3)
