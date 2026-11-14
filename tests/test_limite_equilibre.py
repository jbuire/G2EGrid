"""Un point d'equilibre au-dela d'une limite de convertisseur est signale (FR-009).

Le modele linearise est calcule en mode libre ; au-dela de ``Imax`` (ou de
``Mmax``), il ne represente plus le convertisseur, qui serait en limitation.
On verifie :

* le cas du constat : ``wscc9_gfm``, 200 MW demandes au GFM du bus 3
  (Sn = 130 MVA). ``solve_equilibrium`` range dans ``net.depassements`` le
  seul ``GFM3``, a ``|is| = 1,4998 pu`` pour ``Imax = 1,2 pu``, et emet un
  avertissement ``DepassementLimite`` ; ``GFM2`` (0,79 pu) n'y est pas ;
* le silence au nominal, sur tout le catalogue (l'avertissement y est rendu
  bloquant) ;
* que la grandeur comparee est EXACTEMENT celle de la garde du limiteur, pour
  toutes les familles qui en portent un (GFM, GFL, MMC, ports continus), et
  qu'elle vaut a l'equilibre le courant COTE CONVERTISSEUR ``|is|`` -- pas
  le courant cote reseau ``|ig|`` ;
* l'accord avec la simulation : un depassement est signale si et seulement
  si la simulation, limiteurs surveilles, demarre en limitation ;
* la modulation face a ``Mmax``, recoupee par l'equation de l'inductance ;
* le seuil, un drapeau deja leve, un convertisseur bloque, une variante
  reduite (integrateurs de la boucle de courant figes) ;
* la continuation et la robustesse, qui marquent sans s'arreter ;
* l'avis de la validation de l'atelier (``analyse.valider``).
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import g2egrid.catalogue as cat
from g2egrid import continuation, limites, robustesse
from g2egrid.elements.limiteurs import LimiteurCourant, LimiteurCourantTension
from g2egrid.limites import TOLERANCE, DepassementLimite
from g2egrid.tables import assembler
from g2egrid.temporel import ModeleNonLineaire

#: |is| du GFM3 de wscc9_gfm a 200 MW (pu de 130 MVA) : le chiffre du constat
#: FR-009 (script de la fiche), a la tolerance du Newton d'equilibre pres
IS_200MW = 1.4998447784797937


def _wscc9_gfm(p3_mw: float | None = None):
    """Les tables de wscc9_gfm, la consigne P du GFM du bus 3 changee."""
    t = cat.charger("wscc9_gfm")
    if p3_mw is not None:
        t.Y_DER = [
            tuple(r) if r[0] != 3 else (*r[:5], p3_mw * 1e6, *r[6:]) for r in t.Y_DER
        ]
    return t


def _resolu(t):
    net = assembler(t)
    net.solve_equilibrium()
    return net


def _element(net, nom):
    return next(e for e in net.all_elements if e.name == nom)


def _module(el, a, b):
    x = dict(zip(el.states, el.x0))
    return float(np.hypot(x[a], x[b]))


@pytest.fixture(scope="module")
def surcharge():
    """wscc9_gfm, 200 MW au GFM du bus 3 : l'avertissement est attendu."""
    net = assembler(_wscc9_gfm(200.0))
    with pytest.warns(DepassementLimite, match="GFM3") as w:
        net.solve_equilibrium()
    net.avertissements = [str(x.message) for x in w]
    return net


# ------------------------------------------------------------ le constat
def test_200_mw_au_gfm3_signale(surcharge):
    (d,) = surcharge.depassements
    assert (d.element, d.bus, d.drapeau, d.grandeur) == (
        "GFM3",
        3,
        "mode_lim",
        "courant",
    )
    assert d.valeur == pytest.approx(IS_200MW, abs=1e-8)
    assert d.limite == 1.2 and d.Sn == 130e6
    assert d.rapport == pytest.approx(IS_200MW / 1.2, abs=1e-8)
    # un seul avertissement, qui nomme l'element, le courant et la limite
    (msg,) = surcharge.avertissements
    assert "GFM3 (bus 3) : courant 1.4998 pu > Imax = 1.2 pu" in msg
    assert "GFM2" not in msg
    # GFM2 reste loin de sa limite (0,79 pu)
    assert _module(_element(surcharge, "GFM2"), "isd", "isq") < 0.8


def test_courant_cote_convertisseur(surcharge):
    """La grandeur comparee est |is| (cote convertisseur, celle que sature le
    limiteur), pas |ig| (cote reseau), qui en differe du courant du
    condensateur du filtre : 2,7e-3 pu ici."""
    el = _element(surcharge, "GFM3")
    (d,) = surcharge.depassements
    i_s, i_g = _module(el, "isd", "isq"), _module(el, "igd", "igq")
    assert d.valeur == pytest.approx(i_s, abs=1e-12)
    assert abs(i_g - i_s) > 2e-3


def test_le_modele_libre_reste_calculable(surcharge):
    """Pas d'exception : on peut explorer au-dela de la limite, en le sachant."""
    ss = surcharge.build_state_space()
    assert np.linalg.eigvals(ss.A).real.max() < 1e-6


def test_wscc9_gfm_nominal_silencieux():
    with warnings.catch_warnings():
        warnings.simplefilter("error", DepassementLimite)
        net = _resolu(_wscc9_gfm())
    assert net.depassements == []


@pytest.mark.parametrize("nom", list(cat.RESEAUX))
def test_catalogue_nominal_silencieux(nom):
    """Aucun reseau du catalogue n'a de convertisseur hors limite au point
    nominal (marge mesuree : 81 % de Imax au plus, GFL113 de cigre_b4 ; 87 %
    de Mmax, GFM141)."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", DepassementLimite)
        net = _resolu(cat.charger(nom))
    assert net.depassements == []


# ------------------------------------------------- la grandeur comparee
@pytest.mark.parametrize(
    "nom", ["wscc9_gfm", "smib_gfm", "smib_gfl", "cigre_mv_dc", "cigre_b4_mmc"]
)
def test_meme_grandeur_que_la_garde(nom):
    """``grandeurs_limitees`` rend la grandeur de la garde au bit pres
    (``g = valeur - limite``), et la consigne de courant libre vaut a
    l'equilibre ``|is|`` : ``dM_CL/dt = KiCL.(i_ref - is) = 0``. L'ecart
    tolere, 1e-9 pu, borne ``|dM_CL/dt| / KiCL`` au residu d'equilibre (1e-8)
    et au plus petit KiCL du paquet (58, GFL)."""
    net = _resolu(cat.charger(nom))
    vus = set()
    for el in net.all_elements:
        if not isinstance(el, LimiteurCourant):
            continue
        vus.add(type(el).__name__)
        cm = el.compiled()
        gardes = {d: g for d, g, _v, _s in el.gardes(cm, el.x0, el.z0, el.u0, el.p)}
        lues = el.grandeurs_limitees(cm, el.x0, el.z0, el.u0, el.p)
        assert {d for d, _v, _l in lues} == set(gardes)
        for d, valeur, limite in lues:
            assert valeur - limite == gardes[d]
            if d == "mode_lim":
                assert valeur == pytest.approx(_module(el, "isd", "isq"), abs=1e-9)
    attendus = {
        "wscc9_gfm": {"GFM1"},
        "smib_gfm": {"GFM3"},
        "smib_gfl": {"GFL2"},
        "cigre_mv_dc": {"GFL1", "GFL2DC", "GFM3DC"},
        "cigre_b4_mmc": {"GFL1", "GFM3DC", "MMC"},
    }[nom]
    assert attendus <= vus


def test_accord_avec_la_simulation(surcharge):
    """Le reseau signale un depassement si et seulement si la simulation,
    limiteurs surveilles, bascule ce limiteur des t = 0."""
    for net in (surcharge, _resolu(_wscc9_gfm())):
        tr = ModeleNonLineaire(net).simuler((0, 1e-3), methode="BDF2", pas=1e-4)
        au_depart = {(n, d) for t, n, d, v in tr.basculements if t == 0.0 and v}
        assert au_depart == {(d.element, d.drapeau) for d in net.depassements}
    assert au_depart == set()  # le nominal


# ------------------------------------------------------------ modulation
def _modulation_lcl(el) -> float:
    """``|e|/v_dc`` par l'equation de l'inductance cote onduleur a l'equilibre,
    en repere local : ``e = ve + (R1 + j.w.L1).is``, ``w = 1``."""
    x, p = dict(zip(el.states, el.x0)), el.pdict
    i_s = complex(x["isd"], x["isq"])
    e = complex(x["ved"], x["veq"]) + complex(p["R1_pu"], p["L1_pu"]) * i_s
    return abs(e) / x.get("vdc", 1.0)


def test_modulation_au_dela_de_mmax():
    """Saturation de tension reglee sous la modulation d'equilibre (pres de
    1 pu) : chaque GFM est signale en modulation, la valeur recoupee par
    l'equation de l'inductance. Reglee a MMAX = 1,15, rien."""
    t = _wscc9_gfm()
    t.GFM_PARAMS = dict(t.GFM_PARAMS, Mmax=0.9)
    with pytest.warns(DepassementLimite, match="modulation"):
        net = _resolu(t)
    vsat = {d.element: d for d in net.depassements if d.drapeau == "mode_vsat"}
    assert set(vsat) == {"GFM2", "GFM3"}
    for nom, d in vsat.items():
        assert d.grandeur == "modulation" and d.limite == 0.9
        assert d.valeur == pytest.approx(_modulation_lcl(_element(net, nom)), abs=1e-9)
    t.GFM_PARAMS = dict(t.GFM_PARAMS, Mmax=1.15)
    with warnings.catch_warnings():
        warnings.simplefilter("error", DepassementLimite)
        assert _resolu(t).depassements == []


def test_modulation_du_gfl():
    """GFL : la modulation libre est algebrique (``md``, ``mq``)."""
    t = cat.charger("wscc9_faible")
    t.GFL_PARAMS = dict(t.GFL_PARAMS, Mmax=0.9)
    with pytest.warns(DepassementLimite):
        net = _resolu(t)
    vsat = [d for d in net.depassements if d.drapeau == "mode_vsat"]
    assert {d.element for d in vsat} == {"GFL2", "GFL3"}
    for d in vsat:
        el = _element(net, d.element)
        z = dict(zip(el.algebs, el.z0))
        assert d.valeur == pytest.approx(np.hypot(z["md"], z["mq"]), abs=1e-12)
        assert isinstance(el, LimiteurCourantTension)


# ------------------------------------------- seuil, drapeau, blocage, reduction
def test_seuil_de_la_simulation():
    """Le seuil est ``TOLERANCE`` (1e-9), celui de la simulation au depart :
    une limite posee exactement sur le courant, ou a moins de 1e-9 dessous,
    n'est pas un depassement ; 2e-9 dessous en est un."""
    net = _resolu(_wscc9_gfm())
    el = _element(net, "GFM3")
    i0 = _module(el, "isd", "isq")
    for imax, signale in (
        (i0, False),
        (i0 - 0.5 * TOLERANCE, False),
        (i0 - 2 * TOLERANCE, True),
    ):
        el.set_param("Imax", imax)
        assert bool(limites.depassements(net)) is signale


def test_drapeau_leve_ou_bloque_non_signale():
    """Drapeau deja leve : le modele est celui du mode limite, coherent avec
    le point. Convertisseur bloque : son limiteur n'est plus surveille."""
    with pytest.warns(DepassementLimite):
        net = _resolu(_wscc9_gfm(200.0))
    el = _element(net, "GFM3")
    el.set_param("mode_lim", 1.0)
    assert limites.depassements(net) == []
    el.set_param("mode_lim", 0.0)
    assert len(limites.depassements(net)) == 1
    el.bloque = True
    assert limites.depassements(net) == []


def test_variante_a_boucle_de_courant_figee():
    """Integrateurs de la boucle de courant figes (groupe ``boucle_i`` de
    l'atelier) : la garde ne lit plus ``dM_CL/dt``, prend ``is`` -- exact a
    l'equilibre -- et le depassement est le meme que sur le modele complet.
    Avant, la garde echouait (``ValueError``) faute d'etat ``M_CLd``."""
    net = assembler(_wscc9_gfm(200.0))
    for el in net.all_elements:
        if isinstance(el, LimiteurCourant):
            el.reduire(fige=("M_CLd", "M_CLq"))
    with pytest.warns(DepassementLimite):
        net.solve_equilibrium()
    (d,) = net.depassements
    assert d.element == "GFM3"
    assert d.valeur == pytest.approx(IS_200MW, abs=1e-8)
    el = _element(net, "GFM3")
    (g,) = el.gardes(el.compiled(), el.x0, el.z0, el.u0, el.p)
    assert g[1] == pytest.approx(IS_200MW - 1.2, abs=1e-8)


# --------------------------------------------- continuation et robustesse
def test_continuation_marque_sans_s_arreter():
    """P du GFM3 de 120 a 210 MW : les deux derniers points sont hors limite,
    la suite va jusqu'au bout, la frontiere est raffinee la ou |is| = Imax
    (a la tolerance de la dichotomie : 1e-4 x 150 MW, soit 1,1e-4 pu), et un
    seul avertissement resume la suite."""

    def fabrique(p_mw):
        return _resolu(_wscc9_gfm(p_mw))

    with pytest.warns(DepassementLimite) as w:
        r = continuation.continuer(fabrique, [120.0, 150.0, 180.0, 210.0])
    (avert,) = [x for x in w if issubclass(x.category, DepassementLimite)]
    assert "continuation, 2 point(s) sur 4 hors limite" in str(avert.message)
    assert list(r.mu) == [120.0, 150.0, 180.0, 210.0]
    assert list(r.hors_limite) == [False, False, True, True]
    assert [d.element for d in r.depassements[2]] == ["GFM3"]
    assert 150.0 < r.mu_hors_limite < 180.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DepassementLimite)
        dehors = _element(fabrique(r.mu_hors_limite), "GFM3")
        dedans = _element(fabrique(r.mu_hors_limite - 0.02), "GFM3")
    assert 0.0 < _module(dehors, "isd", "isq") - 1.2 < 2e-4
    assert _module(dedans, "isd", "isq") < 1.2
    assert "hors limite des convertisseurs a partir de mu" in r.resume()


def test_continuation_nominale_silencieuse():
    with warnings.catch_warnings():
        warnings.simplefilter("error", DepassementLimite)
        r = continuation.continuer(lambda p: _resolu(_wscc9_gfm(p)), [85.0, 120.0])
    assert r.mu_hors_limite is None and not r.hors_limite.any()


def test_robustesse_marque_la_variante():
    """Base a 150 MW (|is| = 1,13 pu) : seule la variante « GFM3 P +10% »
    (165 MW, 1,24 pu) sort de la limite. Elle est marquee, pas retiree."""
    with pytest.warns(DepassementLimite) as w:
        e = robustesse.etudier(_wscc9_gfm(150.0), deconnexion=False)
    assert len([x for x in w if issubclass(x.category, DepassementLimite)]) == 1
    assert e.depassements == []
    ((p, deps),) = e.hors_limite
    assert p.libelle == "GFM3 P +10%"
    assert [d.element for d in deps] == ["GFM3"]
    assert deps[0].valeur == pytest.approx(1.2414568657117073, abs=1e-8)
    for m in e.modes:
        for x in m.effets:
            assert x.hors_limite is (x.perturbation.libelle == "GFM3 P +10%")
        assert any(x.hors_limite for x in m.suivis())  # garde dans les scores
    assert "HORS LIMITE apres GFM3 P +10%" in robustesse.resume(e)


# ------------------------------------------------------------- l'atelier
def test_avis_de_la_validation():
    pytest.importorskip("fastapi")
    from g2egrid_api import analyse

    analyse.vider_cache()  # un equilibre deja resolu n'avertirait plus
    with pytest.warns(DepassementLimite):
        avis = [x for x in analyse.valider(_wscc9_gfm(200.0)) if "depassement" in x]
    (a,) = avis
    assert a["gravite"] == "avis" and a["message"].startswith("GFM3 (bus 3)")
    assert a["depassement"]["valeur"] == pytest.approx(IS_200MW, abs=1e-8)
    assert a["depassement"]["limite"] == 1.2
    with warnings.catch_warnings():
        warnings.simplefilter("error", DepassementLimite)
        assert not [x for x in analyse.valider(_wscc9_gfm()) if "depassement" in x]


def test_atelier_recoit_les_marques_hors_limite():
    """Revue de FR-009 : le JSON que l'atelier recoit porte les marques.

    Robustesse : la variante hors limite est nommee avec son detail, ses
    effets sont marques, le reste non. Continuation : masque par ``mu``,
    premier ``mu`` hors limite et detail, par le meme assembleur que les deux
    routes de suite (``continuation_parametre``, ``balayage_statismes``). Le
    tout passe par ``json.dumps``.
    """
    pytest.importorskip("fastapi")
    import json

    from g2egrid_api import analyse

    with pytest.warns(DepassementLimite):
        d = analyse.robustesse(_wscc9_gfm(150.0), zeta_max=30, deconnexion=False)
    json.dumps(d)
    assert d["depassements"] == []
    ((h,),) = [d["hors_limite"]]
    assert h["libelle"] == "GFM3 P +10%"
    assert [x["element"] for x in h["depassements"]] == ["GFM3"]
    effets = [x for m in d["modes"] for x in m["effets"]]
    assert effets
    for x in effets:
        assert x["hors_limite"] is (x["libelle"] == "GFM3 P +10%")

    with pytest.warns(DepassementLimite):
        r = continuation.continuer(lambda p: _resolu(_wscc9_gfm(p)), [120.0, 180.0])
    s = analyse._limites_suite(r)
    json.dumps(s)
    assert s["hors_limite"] == [False, True]
    assert 120.0 < s["mu_hors_limite"] < 180.0
    assert s["depassements"][0] == []
    assert [x["element"] for x in s["depassements"][1]] == ["GFM3"]
    assert all(isinstance(b["hors_limite"], bool) for b in analyse._bifurcations(r))
