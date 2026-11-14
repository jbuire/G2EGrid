"""FR-051, criteres 4, 5 et 6 : equations -> blocs (``depuis_element``).

Tests ecrits A L'AVEUGLE depuis le seul cahier
``.claude/equipe/cahiers/FR-051_editeur_blocs.md`` : ils visent le contrat de
``depuis_element(classe_ou_nom, parametres=None, detail="blocs") -> schema``.

CRITERE 4 (aller-retour DAE -> blocs -> DAE). Pour GFL1, GFL2, GFM1, GFM2,
GFMdVOC, GFMMatching, IB, la machine synchrone (AVR et PSS compris) et une
charge R-L, l'element du schema ``depuis_element(C, parametres)`` (les deux
niveaux de detail) est pose a la place du natif sur un reseau du catalogue qui
le contient, et on exige : meme load flow, initialisation reussie, memes
valeurs propres que le natif a 1e-8 (relatif), mode de jauge exclu.

- Reseaux du catalogue : ``smib_gfl`` (GFL2, et son GFL1 en deplacant la ligne
  de ``T_GFL2s`` vers ``T_GFL1s``), ``wscc9_gfm`` (GFM1 du bus 2), ``smib_gfm``
  (GFM2 par deplacement vers ``T_GFM2s``, GFMdVOC avec ``loi = "dvoc"``,
  GFMMatching avec ``loi = "matching"``), ``smib_sm`` (IB du bus 1), ``smib``
  (machine synchrone), ``wscc9`` (charge R-L du bus 5).
- Les ``parametres`` donnes sont ceux de l'element natif ``el.pdict`` (apres
  initialisation : pour la charge, R et X recalcules a la tension du load
  flow), sans les reserves ``kS`` et ``wb``. Lecture du cahier : ``parametres``
  nomme les parametres natifs de la classe (les cles de ``Element.params``).
- La charge est remplacee en mettant son ``Pl, Ql`` a zero dans ``Y_DER`` et
  en declarant la meme puissance en generation NEGATIVE au meme bus (meme
  injection nette pour le load flow) avant d'y poser l'element a blocs,
  ``Sn_VA = Sb``.
- 1e-8 : seuil du cahier. L'ecart relatif est ``|a - b| / max(1, |a|)`` apres
  appariement optimal (``test_fr051_aide.ecart_spectres``) ; le bruit de
  ``eig`` (~1e-13) est quatre ordres de grandeur sous le seuil.
- Load flow : meme tension de bus a 1e-10 pu pres (memes entrees, meme
  Newton ; l'element a blocs n'y ajoute rien).

CRITERE 5 (edition). ``Kipll`` x 2 dans les parametres du schema genere d'un
GFL donne le spectre du natif avec ce parametre double, a 1e-8 ; meme chose
pour ``mp`` x 2 sur un GFM2. On verifie aussi que l'edition DEPLACE le spectre
(sinon le test serait vide).

CRITERE 6 (disposition). Sur chaque schema genere (tous niveaux, deux niveaux
de detail), chaque bloc a ``x`` et ``y`` et les boites 92 x 50 sont disjointes.

AUSSI, depuis le cahier : un integrateur (K = 1) par etat nomme d'apres l'etat ;
un sous-systeme ``d<etat>/dt`` par equation d'etat (detail ``blocs``) ; entrees
``vd, vq, w, theta_g`` ; sorties ``id, iq`` (+ ``theta, w`` pour un formeur) et
jamais P, Q, Vt ; deux consignes calees dont ``p_ref`` ; parametres natifs et
drapeaux de mode a 0, ``kS`` et ``wb`` absents ; ``amorce`` ; refus des
elements a port continu.

Lectures du cahier retenues : « nomme d'apres l'etat » = l'identifiant de
l'integrateur contient le nom de l'etat ; le refus d'un element a port continu
est une ``SchemaInvalide`` (donc une ``ValueError``) dont le message contient
« continu » ; ``amorce.correspondance`` va de variables du schema (etats ou
algebriques compilees) vers des variables natives (``states`` ou ``algebs``).

Lent : les 18 comparaisons de spectre assemblent chacune deux reseaux (de 3
bus a 9 bus avec 3 machines synchrones) : de l'ordre de la seconde chacune,
la charge et les GFM1 de ``wscc9`` quelques secondes.
"""

from __future__ import annotations

import copy

import g2egrid.catalogue as cat
import numpy as np
import pytest
import sympy as sp
from test_fr051_aide import (
    copie,
    ecart_spectres,
    recouvrements,
    sans_jauge,
    spectre,
    tous_les_niveaux,
)

from g2egrid.elements import (
    GFL1,
    GFL2,
    GFL2DC,
    GFM1,
    GFM2,
    GFM3DC,
    GFM4DC,
    IB,
    MMC,
    GFMdVOC,
    GFMMatching,
    RLLoad,
    SynchronousMachine,
)
from g2egrid.elements import blocs as B
from g2egrid.elements.blocs import SchemaInvalide, compiler


def depuis_element(*args, **kwargs):
    """Appel tardif : le module se collecte meme si la fonction manque."""
    return B.depuis_element(*args, **kwargs)


TOL_SPECTRE = 1e-8  # cahier, criteres 4 et 5
TOL_LF = 1e-10
DETAILS = ("blocs", "equations")
RESERVES = ("kS", "wb")


# ------------------------------------------------------------ les reseaux
def _charger(nom: str):
    return copy.deepcopy(cat.charger(nom))


def _gfl1():
    t = _charger("smib_gfl")
    t.T_GFL1s, t.T_GFL2s = t.T_GFL2s, []
    return t, "T_GFL1s", 0, "GFL2"


def _gfl2():
    return _charger("smib_gfl"), "T_GFL2s", 0, "GFL2"


def _gfm1():
    return _charger("wscc9_gfm"), "T_GFM1s", 0, "GFM2"


def _gfm2():
    t = _charger("smib_gfm")
    t.T_GFM2s, t.T_GFM3s = t.T_GFM3s, []
    return t, "T_GFM2s", 0, "GFM2"


def _gfm_dvoc():
    t, table, idx, nom = _gfm2()
    t.GFM_PARAMS = dict(t.GFM_PARAMS, loi="dvoc")
    return t, table, idx, nom


def _gfm_matching():
    t = _charger("smib_gfm")
    t.GFM_PARAMS = dict(t.GFM_PARAMS, loi="matching")
    return t, "T_GFM3s", 0, "GFM2"


def _ib():
    return _charger("smib_sm"), "T_IBs", 0, "IB1"


def _sm():
    return _charger("smib"), "T_SMs", 0, "SM1"


def _charge():
    return _charger("wscc9"), "Y_DER", 4, "Load5"


CAS = {
    "GFL1": (GFL1, _gfl1),
    "GFL2": (GFL2, _gfl2),
    "GFM1": (GFM1, _gfm1),
    "GFM2": (GFM2, _gfm2),
    "GFMdVOC": (GFMdVOC, _gfm_dvoc),
    "GFMMatching": (GFMMatching, _gfm_matching),
    "IB": (IB, _ib),
    "SM": (SynchronousMachine, _sm),
    "RLLoad": (RLLoad, _charge),
}


def _vers_blocs(t, table: str, index: int, schema: dict):
    """Pose ``schema`` a la place de la ligne ``table[index]`` : meme bus, meme
    ``Sn_VA`` (cahier, U4). Pour une charge (ligne de ``Y_DER``), la meme
    puissance est declaree en generation negative."""
    if table == "Y_DER":
        ligne = list(t.Y_DER[index])
        bus, Pl, Ql = int(ligne[0]), ligne[7], ligne[8]
        ligne[5], ligne[6], ligne[7], ligne[8] = -Pl, -Ql, 0.0, 0.0
        t.Y_DER[index] = tuple(ligne)
        Sn = t.Sb
    else:
        ligne = getattr(t, table).pop(index)
        bus, Sn = int(ligne[1]), ligne[2]
    t.T_Blocs = [*t.T_Blocs, (len(t.T_Blocs) + 1, bus, "modele", Sn)]
    t.MODELES_BLOCS = {**t.MODELES_BLOCS, "modele": schema}
    return t


_NATIFS: dict[str, tuple] = {}


def _natif(cas: str):
    """(tables, table, index, load flow, spectre, parametres natifs), calcule
    une fois par cas."""
    if cas not in _NATIFS:
        cls, fab = CAS[cas]
        t, table, index, nom = fab()
        net, lf, lam = spectre(copy.deepcopy(t))
        el = next(e for e in net.all_elements if e.name == nom)
        assert type(el) is cls, (type(el), cls)
        params = {k: v for k, v in el.pdict.items() if k not in RESERVES}
        _NATIFS[cas] = (t, table, index, lf, lam, params)
    return _NATIFS[cas]


# ================================================= critere 4 : aller-retour
@pytest.mark.parametrize("detail", DETAILS)
@pytest.mark.parametrize("cas", list(CAS))
def test_aller_retour_meme_load_flow_et_memes_valeurs_propres(cas, detail):
    cls, _ = CAS[cas]
    t, table, index, lf0, lam0, params = _natif(cas)
    schema = depuis_element(cls, parametres=params, detail=detail)
    tb = _vers_blocs(copy.deepcopy(t), table, index, schema)
    _net, lf1, lam1 = spectre(tb)  # le load flow converge, l'element s'initialise
    assert np.max(np.abs(lf0.V - lf1.V)) < TOL_LF
    assert ecart_spectres(sans_jauge(lam0), sans_jauge(lam1)) < TOL_SPECTRE


def test_classe_ou_nom_donne_le_meme_schema():
    assert depuis_element("GFL1") == depuis_element(GFL1)
    assert depuis_element("GFM2", detail="equations") == depuis_element(
        GFM2, detail="equations"
    )


# ===================================== contenu du schema genere (cahier)
@pytest.mark.parametrize("detail", DETAILS)
@pytest.mark.parametrize("cas", list(CAS))
def test_structure_du_schema_genere(cas, detail):
    cls, _ = CAS[cas]
    _t, _tab, _i, _lf, _lam, params = _natif(cas)
    s = depuis_element(cls, parametres=params, detail=detail)
    blocs = s["blocs"]
    c = compiler(copie(s))
    tous = [b for _, bs in tous_les_niveaux(s) for b in bs]

    # un integrateur par etat, nomme d'apres l'etat. Cahier amende a la
    # reunion (2026-11-10) : le K de l'integrateur porte le facteur
    # PARAMETRIQUE de la derivee (wb/L1_pu, KiCL...) au lieu de 1 ; avec K = 1
    # ces facteurs (6e5) passaient dans les equations algebriques et
    # l'assembleur refusait le modele, mal conditionne. K ne depend donc que
    # des parametres du modele (et des reserves wb, kS), jamais d'un signal.
    integ = [b for b in blocs if b["type"] == "integrateur"]
    assert len(integ) == len(cls.states) == len(c.etats)
    admis = set(s["parametres"]) | set(RESERVES)
    for b in integ:
        libres = {str(x) for x in sp.sympify(str(b.get("K", 1))).free_symbols}
        assert libres <= admis, (b["id"], b.get("K"))
    for etat in cls.states:
        assert any(etat in b["id"] for b in integ), etat

    # equations d'etat : un sous-systeme d<etat>/dt (blocs), ou des blocs
    # fonction (equations)
    noms_ss = {b.get("nom") for b in tous if b["type"] == "sous_systeme"}
    if detail == "blocs":
        # cahier amende a la reunion : une derivee qui n'est qu'un signal (une
        # fois son facteur parametrique passe dans le K de l'integrateur, cf.
        # plus haut) se relie directement a l'integrateur, sans sous-systeme
        # vide ; toute autre derivee a son sous-systeme d<etat>/dt
        for b in integ:
            etat = next(e for e in cls.states if e in b["id"])
            source = next(
                (lien["de"] for lien in s["liens"] if lien["vers"] == b["id"]), None
            )
            assert source is not None, etat
            type_source = next(x["type"] for x in blocs if x["id"] == source)
            if type_source == "sous_systeme":
                assert f"d{etat}/dt" in noms_ss, etat
    else:
        assert not any(n and n.endswith("/dt") for n in noms_ss)
        assert sum(b["type"] == "fonction" for b in tous) >= len(cls.states)

    # entrees et sorties du reseau
    entrees = {b["signal"] for b in blocs if b["type"] == "entree"}
    assert {"vd", "vq"} <= entrees <= {"vd", "vq", "w", "theta_g"}
    sorties = {b["signal"] for b in tous if b["type"] == "sortie"}
    assert {"id", "iq"} <= sorties <= {"id", "iq", "theta", "w"}

    # consignes : deux calees, dont p_ref
    consignes = [b for b in tous if b["type"] == "consigne"]
    assert sum(b.get("calee", True) for b in consignes) == 2
    # cahier amende a la reunion : p_ref est calee quand la classe en a une
    # (IB, la charge et GFL2, qui tient la tension continue, n'en ont pas qui
    # commande P : voir blocs_inverse.CONSIGNES)
    if "p_ref" in cls.svc_inputs:
        assert any(b["nom"] == "p_ref" and b.get("calee", True) for b in consignes)
    assert c.etats and len(c.consignes) >= 2

    # parametres natifs, drapeaux de mode a 0, reserves absentes
    p = s["parametres"]
    assert not set(RESERVES) & set(p)
    for k, v in p.items():
        if k in params and not k.startswith("mode_"):
            assert v == pytest.approx(params[k], rel=1e-12), k
    for drapeau in getattr(cls, "modes", ()):
        assert p.get(drapeau) == 0, drapeau

    # amorce de l'initialisation
    amorce = s["amorce"]
    assert amorce["classe"] == cls.__name__
    # les consignes sont aussi des inconnues de l'initialisation : la
    # correspondance peut les apparier aux entrees de service natives
    noms_schema = set(c.etats) | set(c.algebriques) | set(c.consignes)
    natives = set(cls.states) | set(cls.algebs) | set(cls.svc_inputs)
    corr = amorce["correspondance"]
    assert corr
    assert set(corr) <= noms_schema
    assert set(corr.values()) <= natives


def test_les_courants_non_formeurs_n_ont_pas_de_sortie_de_service():
    """GFL et charge ne forment pas le repere : ni ``theta`` ni ``w`` en
    sortie. Un GFM ou un IB (formeurs) les portent."""
    for cas, attendu in (("GFL1", False), ("RLLoad", False), ("GFM2", True)):
        cls, _ = CAS[cas]
        _t, _tab, _i, _lf, _lam, params = _natif(cas)
        s = depuis_element(cls, parametres=params)
        sorties = {
            b["signal"]
            for _, bs in tous_les_niveaux(s)
            for b in bs
            if b["type"] == "sortie"
        }
        assert ({"theta", "w"} <= sorties) is attendu, cas


def test_parametres_par_defaut_et_surcharge():
    from g2egrid.elements.gfl import GFL_DEFAULT

    s = depuis_element("GFL1")
    compiler(copie(s))
    assert s["parametres"]["R1_pu"] == pytest.approx(GFL_DEFAULT["R1_pu"])
    s2 = depuis_element("GFL1", parametres={"R1_pu": 0.021})
    assert s2["parametres"]["R1_pu"] == pytest.approx(0.021)


@pytest.mark.parametrize("classe", [GFM3DC, GFM4DC, GFL2DC, MMC, "GFM3DC"])
def test_element_a_port_continu_refuse(classe):
    """Port alternatif ``t`` seulement en v1 : refus clair."""
    with pytest.raises(SchemaInvalide, match="continu"):
        depuis_element(classe)


def test_detail_inconnu_refuse():
    with pytest.raises(SchemaInvalide):
        depuis_element(GFL1, detail="n_importe_quoi")


# ================================================================ critere 5
@pytest.mark.parametrize(
    "cas, cle",
    [("GFL1", "Kipll"), ("GFL2", "Kipll"), ("GFM2", "mp")],
)
def test_multiplier_un_gain_du_schema_donne_le_spectre_du_natif(cas, cle):
    cls, _ = CAS[cas]
    t, table, index, _lf, lam0, params = _natif(cas)
    # natif avec le parametre double
    tn = copy.deepcopy(t)
    if cas.startswith("GFL"):
        tn.GFL_PARAMS = dict(tn.GFL_PARAMS, **{cle: 2 * params[cle]})
    else:
        tn.GFM_PARAMS = dict(tn.GFM_PARAMS, **{cle: 2 * params[cle]})
    _n, _lf, lam_natif = spectre(tn)
    # le schema genere, parametre double
    s = depuis_element(cls, parametres=params)
    s["parametres"][cle] = 2 * s["parametres"][cle]
    _n, _lf, lam_blocs = spectre(_vers_blocs(copy.deepcopy(t), table, index, s))
    assert ecart_spectres(sans_jauge(lam_natif), sans_jauge(lam_blocs)) < TOL_SPECTRE
    # l'edition compte : le spectre a bien bouge (sinon le test est vide)
    assert ecart_spectres(sans_jauge(lam0), sans_jauge(lam_natif)) > 1e-3


# ================================================================ critere 6
@pytest.mark.parametrize("detail", DETAILS)
@pytest.mark.parametrize("cas", list(CAS))
def test_aucun_recouvrement_a_aucun_niveau(cas, detail):
    cls, _ = CAS[cas]
    s = depuis_element(cls, detail=detail)
    n_niveaux = 0
    for chemin, blocs in tous_les_niveaux(s):
        n_niveaux += 1
        for b in blocs:
            assert isinstance(b.get("x"), int | float), (chemin, b["id"])
            assert isinstance(b.get("y"), int | float), (chemin, b["id"])
        assert recouvrements(blocs) == [], chemin
    if detail == "blocs":
        assert n_niveaux > 1  # des sous-systemes ont bien ete poses


@pytest.mark.parametrize("cas", ["GFL2", "GFM2"])
def test_flux_de_gauche_a_droite(cas):
    """Disposition en flux de gauche a droite : une entree du reseau est a
    gauche d'une sortie vers le reseau (premier niveau)."""
    cls, _ = CAS[cas]
    s = depuis_element(cls)
    xe = [b["x"] for b in s["blocs"] if b["type"] == "entree"]
    xs = [b["x"] for b in s["blocs"] if b["type"] == "sortie"]
    assert max(xe) < max(xs) and min(xe) <= min(xs)


def test_depuis_element_est_exporte_avec_ses_voisins():
    for nom in (
        "TYPES",
        "compiler",
        "classe_element",
        "aplatir",
        "regrouper",
        "degrouper",
        "depuis_element",
        "SchemaInvalide",
    ):
        assert hasattr(B, nom), nom
