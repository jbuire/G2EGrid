"""FR-051, suites de la relecture : ce que la relecture a trouve, verrouille.

* toute classe native a un seul port alternatif ``t`` donne, par
  ``depuis_element``, un schema qui COMPILE, ou un refus explicite
  (``SchemaInvalide``, 422 par l'API) -- jamais un schema inutilisable
  (``ZLoad``, sans entree de service, en rendait un) ni une erreur interne
  (``ReactanceShunt``, a argument obligatoire, donnait un 500) ;
* un parametre inconnu ou non numerique est refuse, pas ignore ;
* NaN, l'infini et les booleens sont refuses dans un champ numerique de bloc,
  et la bibliotheque refuse un schema qui en contient (ce n'est pas du JSON) ;
* un parametre recale par l'amorce (``R``, ``X`` d'une charge) se recale a
  CHAQUE initialisation, comme dans le natif : une seconde initialisation au
  point B donne ce que donne un element neuf initialise au point B.
"""

from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from g2egrid import elements
from g2egrid.elements.blocs import SchemaInvalide, classe_element, compiler
from g2egrid.elements.blocs_inverse import depuis_element
from g2egrid.noyau.element import Element, InitContext


def _classes_un_port_t() -> list[type]:
    out = []
    for nom in dir(elements):
        c = getattr(elements, nom)
        if (
            inspect.isclass(c)
            and issubclass(c, Element)
            and c is not Element
            and tuple(c.ports) == ("t",)
            and not c.ports_dc
        ):
            out.append(c)
    return sorted(set(out), key=lambda c: c.__name__)


@pytest.mark.parametrize("cls", _classes_un_port_t(), ids=lambda c: c.__name__)
def test_toute_classe_a_un_port_compile_ou_est_refusee(cls):
    try:
        s = depuis_element(cls)
    except SchemaInvalide:
        return  # refus explicite : c'est permis
    c = compiler(s)  # ne doit pas lever : le schema rendu est utilisable
    assert len(c.consignes) >= 2


def test_classe_sans_consigne_refusee_clairement():
    with pytest.raises(SchemaInvalide, match="entree"):
        depuis_element("ZLoad")


def test_argument_obligatoire_manquant_refuse_sans_erreur_interne():
    with pytest.raises(SchemaInvalide):
        depuis_element("ReactanceShunt")


@pytest.mark.parametrize(
    "parametres, mot",
    [
        ({"Kpll_inexistant": 1.0}, "inconnu"),
        ({"Kipll": "abc"}, "nombre"),
        ({"Kipll": None}, "nombre"),
        ({"Kipll": float("nan")}, "nombre"),
        ({"Kipll": True}, "nombre"),
    ],
)
def test_parametre_invalide_refuse(parametres, mot):
    with pytest.raises(SchemaInvalide, match=mot):
        depuis_element("GFL1", parametres=parametres)


def _gain(K) -> dict:
    return {
        "nom": "g",
        "parametres": {},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref"},
            {"id": "qr", "type": "consigne", "nom": "Q_ref"},
            {"id": "k", "type": "gain", "K": K},
            {"id": "od", "type": "sortie", "signal": "id"},
            {"id": "oq", "type": "sortie", "signal": "iq"},
        ],
        "liens": [
            {"de": "pr", "vers": "k"},
            {"de": "k", "vers": "od"},
            {"de": "qr", "vers": "oq"},
        ],
    }


@pytest.mark.parametrize("K", [float("nan"), float("inf"), -math.inf, True])
def test_valeur_non_finie_ou_booleenne_refusee(K):
    compiler(_gain(2.0))  # temoin : le meme schema sain compile
    with pytest.raises(SchemaInvalide, match="illisible"):
        compiler(_gain(K))


def test_bibliotheque_refuse_nan(tmp_path, monkeypatch):
    import g2egrid.reseau.catalogue as cat
    from g2egrid_api.analyse import blocs_edition

    monkeypatch.setattr(cat, "DOSSIER_RESEAUX", tmp_path)
    blocs_edition.bibliotheque_enregistrer(_gain(2.0))  # temoin
    with pytest.raises(SchemaInvalide, match="non finie"):
        blocs_edition.bibliotheque_enregistrer(_gain(float("nan")))


def test_seconde_initialisation_recale_comme_un_element_neuf():
    """Charge R-L : ``R`` et ``X`` sont recales par l'amorce sur le point de
    fonctionnement. Initialisee au point A puis au point B, l'instance doit
    rendre ce que rend une instance neuve initialisee au point B."""
    s = depuis_element("RLLoad")
    if not (s.get("amorce") or {}).get("recales"):
        pytest.skip("aucun parametre recale pour cette classe")
    cls = classe_element(s)
    reg = dict(cls.compilation.parametres)
    a = InitContext(v={"t": 1.0 + 0j}, s={"t": -(0.5 + 0.2j)})
    b = InitContext(v={"t": 0.98 + 0.05j}, s={"t": -(0.8 + 0.3j)})

    el = cls("L", 1, **reg)
    el.initialize(a)
    deux_fois = el.initialize(b)
    p_deux = dict(el.pdict)

    neuf = cls("L", 1, **reg)
    une_fois = neuf.initialize(b)
    for k in s["amorce"]["recales"]:
        if k in neuf.pdict:
            assert p_deux[k] == pytest.approx(neuf.pdict[k], rel=1e-12), k
    for k, v in une_fois.items():
        assert deux_fois[k] == pytest.approx(v, rel=1e-9, abs=1e-12), k
    assert np.isfinite(list(deux_fois.values())).all()
