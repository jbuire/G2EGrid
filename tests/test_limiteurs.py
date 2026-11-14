"""Limitation de courant des convertisseurs, butees d'excitation, defauts."""

from __future__ import annotations

import numpy as np
import pytest

import g2egrid.catalogue as cat
from g2egrid.tables import assembler
from g2egrid.temporel import Manoeuvre, ModeleNonLineaire


def _courant(tr, nom):
    return np.hypot(tr.etat(f"{nom}.isd"), tr.etat(f"{nom}.isq"))


@pytest.fixture(scope="module")
def defaut():
    """Defaut franc au bus 5 de wscc9_gfm, de 0,05 a 0,10 s, avec et sans
    limiteurs (BDF2 : Radau suivrait les resonances du defaut a 10 us)."""
    m = ModeleNonLineaire(cat.construire("wscc9_gfm"))
    u = m.scenario(
        Manoeuvre(0.05, "defaut", "bus5", R=1e-3, X=1e-2),
        Manoeuvre(0.10, "elimination", "bus5"),
    )
    kw = dict(uext_fn=u, methode="BDF2", pas=1e-4)
    return (
        m.simuler((0.0, 0.14), limiteurs=True, **kw),
        m.simuler((0.0, 0.14), limiteurs=False, **kw),
    )


def test_la_limitation_tient_le_courant(defaut):
    avec, sans = defaut
    pendant = (avec.t > 0.051) & (avec.t < 0.10)
    for g in ("GFM2", "GFM3"):
        assert _courant(avec, g)[pendant].max() < 1.25
        assert _courant(sans, g)[(sans.t > 0.051) & (sans.t < 0.10)].max() > 2.0
    entrees = [b for b in avec.basculements if b[2] == "mode_lim" and b[3] == 1.0]
    assert {b[1] for b in entrees} == {"GFM2", "GFM3"}
    assert min(b[0] for b in entrees) == pytest.approx(0.05, abs=2e-3)


def test_la_butee_d_excitation(defaut):
    """La tension s'effondre, l'AVR pousse : l'amplificateur bute a VRmax."""
    avec, sans = defaut
    assert any(
        b[1] == "SM1" and b[2] == "mode_e2" and b[3] == 1.0 for b in avec.basculements
    )
    assert avec.etat("SM1.e2").max() < 5.0 + 1e-3
    assert sans.etat("SM1.e2").max() > 5.0


def test_sans_franchissement_rien_ne_change():
    """Un petit echelon ne touche aucune limite : les deux simulations sont
    identiques, et le modele linearise (mode libre) aussi."""
    m = ModeleNonLineaire(cat.construire("wscc9_gfm"))
    u = m.echelon("GFM2.p_ref", 0.005, 0.02)
    a = m.simuler((0, 0.1), uext_fn=u, limiteurs=True)
    b = m.simuler((0, 0.1), uext_fn=u, limiteurs=False)
    assert a.basculements == []
    assert np.array_equal(a.x, b.x)


def test_une_unite_surchargee_part_en_limitation():
    """Imax sous le courant d'equilibre : la simulation commence limitee, et
    l'equilibre l'avait signale (FR-009, ``tests/test_limite_equilibre.py``)."""
    from g2egrid.limites import DepassementLimite

    t = cat.charger("wscc9_gfm")
    bus = int(t.T_GFM1s[0][1])
    t.PARAMS_PAR_BUS = {
        **t.PARAMS_PAR_BUS,
        bus: {**t.PARAMS_PAR_BUS.get(bus, {}), "Imax": 0.3},
    }
    net = assembler(t)
    with pytest.warns(DepassementLimite, match=f"GFM{bus} "):
        net.solve_equilibrium()
    m = ModeleNonLineaire(net)
    tr = m.simuler((0, 0.01), methode="BDF2", pas=1e-4)
    assert (0.0, f"GFM{bus}", "mode_lim", 1.0) in tr.basculements


def test_saturation_de_tension_sur_defaut():
    """Saturation de tension reglee (Mmax = 1,15) : sur un defaut franc, la
    modulation sature AVANT que la limitation de courant ne prenne le relais,
    et la modulation reste dans |m| <= Mmax."""
    from g2egrid.elements.limiteurs import MMAX

    t = cat.charger("wscc9_gfm")
    t.GFM_PARAMS = dict(t.GFM_PARAMS, Mmax=MMAX)
    m = ModeleNonLineaire(assembler(t))
    u = m.scenario(Manoeuvre(0.05, "defaut", "bus5", R=1e-3, X=1e-2))
    tr = m.simuler((0.0, 0.07), uext_fn=u, methode="BDF2", pas=1e-5)
    genres = [(b[2], b[3]) for b in tr.basculements]
    assert ("mode_vsat", 1.0) in genres and ("mode_lim", 1.0) in genres
    premier = {}
    for b in tr.basculements:
        premier.setdefault(b[2], b[0])
    assert premier["mode_vsat"] < premier["mode_lim"]


def test_saturation_de_tension_inactive_par_defaut():
    """Mmax = 0 (defaut) : pas de garde de tension, le modele est celui
    d'avant au bit pres."""
    net = cat.construire("wscc9_gfm")
    el = next(e for e in net.elements if e.name.startswith("GFM"))
    assert el.pdict["Mmax"] == 0.0 and el.pdict["mode_vsat"] == 0.0


def test_garde_de_courant_independante_de_la_saturation_de_tension():
    """La garde de courant lit la consigne LIBRE sur l'integrateur de la
    boucle de courant ; la saturation de tension y ajoute son retro-calcul.
    Ce terme ne doit pas fausser la lecture : meme garde de courant, que la
    saturation soit active ou non, a etat donne."""
    from g2egrid.elements.limiteurs import MMAX
    from g2egrid.temporel import ModeleNonLineaire

    t = cat.charger("wscc9_gfm")
    t.GFM_PARAMS = dict(t.GFM_PARAMS, Mmax=MMAX, Taw_v=5e-5)
    m = ModeleNonLineaire(assembler(t))
    k = next(i for i, e in enumerate(m.c.els) if e.name.startswith("GFM"))
    el, cm = m.c.els[k], m.c.cms[k]
    x, z, u = (np.array(v, dtype=float) for v in (el.x0, el.z0, el.u0))
    # un etat hors equilibre, ou la modulation libre depasse Mmax
    x[el.states.index("M_CLd")] += 3.0
    p = np.array(el.p, dtype=float)
    gardes = {}
    for v in (0.0, 1.0):
        p[el.params.index("mode_vsat")] = v
        gardes[v] = dict((d, g) for d, g, _v, _s in el.gardes(cm, x, z, u, p))
    assert gardes[0.0]["mode_lim"] == pytest.approx(gardes[1.0]["mode_lim"], abs=1e-12)


def test_saturation_de_tension_relache_avec_le_retro_calcul_par_defaut():
    """Saturation reglee (Mmax = 1,15), retro-calcul PAR DEFAUT (tau_CL/2) :
    sur un defaut a la charge du bus 5, la modulation sature quelques
    millisecondes puis relache, et la limitation de courant tient le courant
    pendant et apres le defaut. (Avec l'ancien defaut de 20 ms, le Newton
    lachait des 51 ms ; voir GFM_DEFAULT.)"""
    from g2egrid.elements.limiteurs import MMAX

    t = cat.charger("wscc9_gfm")
    t.GFM_PARAMS = dict(t.GFM_PARAMS, Mmax=MMAX)
    net = assembler(t)
    gfms = [e for e in net.elements if e.name.startswith("GFM")]
    for e in gfms:
        assert e.pdict["Taw_v"] == pytest.approx(e.tuning["tau_CL"] / 2)
    m = ModeleNonLineaire(net)
    u = m.scenario(
        Manoeuvre(0.05, "defaut", "bus5", R=1e-3, X=1e-2),
        Manoeuvre(0.10, "elimination", "bus5"),
    )
    tr = m.simuler((0.0, 0.16), uext_fn=u, methode="BDF2", pas=2e-5)
    for e in gfms:
        vsat = [b[3] for b in tr.basculements if b[1] == e.name and b[2] == "mode_vsat"]
        assert vsat[:2] == [1.0, 0.0]  # elle s'enclenche, puis relache
        i = _courant(tr, e.name)
        assert i[(tr.t > 0.051) & (tr.t < 0.10)].max() < 1.3
        assert i[tr.t > 0.13].max() < 1.25


def test_saturation_de_tension_du_gfl():
    """GFL (WSCC 9 faible, SCR 1,3 et 1,95) : saturation inactive par defaut,
    retro-calcul par defaut tau_CL/2. Reglee (Mmax = 1,15), sur un defaut au
    bus 7 elimine a 0,10 s, chaque GFL part en limitation de courant pendant
    le defaut, sature a l'elimination (la tension revient, la boucle de
    courant demande plus que l'onduleur ne sait faire), puis relache.

    Aux noeuds ideaux (le defaut, FR-005), la modulation de GFL2 sature aussi
    5 ms a l'apparition du defaut, vers 0,053 s ; une resistance de causalite
    de 1e5 pu donne les memes basculements a 1e-5 s pres, celle de 100 pu de
    Node.m les masquait en chargeant les noeuds."""
    from g2egrid.elements.gfl import GFL_DEFAULT
    from g2egrid.elements.limiteurs import MMAX

    t = cat.charger("wscc9_faible")
    for e in assembler(t).elements:
        if e.name.startswith("GFL"):
            assert e.pdict["Mmax"] == 0.0 and e.pdict["mode_vsat"] == 0.0
            assert e.pdict["Taw_v"] == pytest.approx(GFL_DEFAULT["tau_CL"] / 2)
    t.GFL_PARAMS = dict(t.GFL_PARAMS, Mmax=MMAX)
    net = assembler(t)
    gfls = [e.name for e in net.elements if e.name.startswith("GFL")]
    m = ModeleNonLineaire(net)
    u = m.scenario(
        Manoeuvre(0.05, "defaut", "bus7", R=1e-3, X=1e-2),
        Manoeuvre(0.10, "elimination", "bus7"),
    )
    tr = m.simuler((0.0, 0.16), uext_fn=u, methode="BDF2", pas=1e-4)
    for nom in gfls:
        b = [(s, d, v) for s, n, d, v in tr.basculements if n == nom]
        assert ("mode_lim", 1.0) in [x[1:] for x in b if x[0] < 0.10]
        apres = [x for x in b if x[1] == "mode_vsat" and x[0] > 0.10]
        assert apres[0][2] == 1.0  # elle sature a l'elimination
        assert 0.0 in [x[2] for x in apres]  # puis relache
        i = _courant(tr, nom)
        assert i[(tr.t > 0.051) & (tr.t < 0.10)].max() < 1.3


# ------------------------------------------- bus continu dynamique (MVDC)
def _defaut_bus15(reglages_station: dict | None = None):
    """cigre_mv_dc, defaut au bus 15 (station GFL2DC qui tient la tension de
    la liaison MVDC) de 0,05 a 0,10 s, simule jusqu'a 0,6 s ;
    ``reglages_station`` s'applique a cette station."""
    t = cat.charger("cigre_mv_dc")
    if reglages_station:
        t.PARAMS_PAR_BUS[15] = dict(reglages_station)
    m = ModeleNonLineaire(assembler(t))
    u = m.scenario(
        Manoeuvre(0.05, "defaut", "bus15", R=1e-3, X=1e-2),
        Manoeuvre(0.10, "elimination", "bus15"),
    )
    return m.simuler((0.0, 0.6), uext_fn=u, methode="BDF2", pas=1e-4)


def test_saturation_active_par_defaut_sur_bus_continu_dynamique():
    """Les GFM dont la tension continue est une variable (ici la station
    GFM3DC de la liaison MVDC de cigre_mv_dc) ont la saturation de tension
    active par defaut ; les GFL, et les GFM sur bus continu ideal, non. Un
    reglage de famille l'emporte, comme pour tout parametre."""
    from g2egrid.elements.limiteurs import MMAX

    net = cat.construire("cigre_mv_dc")
    mmax = {
        e.name: e.pdict["Mmax"] for e in net.elements if e.name[:3] in ("GFM", "GFL")
    }
    assert mmax == {
        "GFM16": MMAX,
        "GFL15": 0.0,
        "GFL6": 0.0,
        "GFL11": 0.0,
        "GFL14": 0.0,
    }
    # GFM2 : bus continu ideal
    gfm = next(e for e in cat.construire("wscc9_gfm").elements if e.name[:3] == "GFM")
    assert gfm.pdict["Mmax"] == 0.0
    t = cat.charger("cigre_mv_dc")
    t.GFM_PARAMS = dict(t.GFM_PARAMS, Mmax=0.0)
    gfm = next(e for e in assembler(t).elements if e.name == "GFM16")
    assert gfm.pdict["Mmax"] == 0.0


def test_defaut_a_la_station_qui_tient_la_tension_continue():
    """Pendant le defaut au bus 15, la station GFM du bus 16 soutient son
    reseau en puisant dans la liaison continue, que plus rien ne recharge.
    Modulation non bornee, v_dc s'effondrait en un temps fini et le Newton
    lachait a 83 ms ; saturee, la modulation ne sait plus fabriquer la
    tension alternative quand v_dc tombe : le transit s'arrete, v_dc se tient
    vers 0,4 pu, puis revient a 1 pu apres l'elimination."""
    tr = _defaut_bus15()
    assert tr.t[-1] == pytest.approx(0.6)
    vdc = tr.etat("bus101.v")
    assert 0.3 < vdc.min() < 0.5
    assert abs(vdc[-1] - 1.0) < 0.01
    assert ("GFM16", "mode_vsat", 1.0) in [b[1:] for b in tr.basculements]


def test_boucle_continue_proportionnelle_garde_son_biais():
    """La boucle de tension continue des GFL2 est proportionnelle (Kidc = 0) :
    son etat M_d, le courant du point de fonctionnement, ne doit pas bouger
    quand la station passe en limitation de courant. Sinon, aucune integrale
    ne le ramenant, v_dc garde l'ecart : 4 % apres ce defaut, la station
    GFL2DC (et elle seule) reglee avec la saturation de tension, qui la fait
    limiter a l'elimination."""
    from g2egrid.elements.limiteurs import MMAX

    tr = _defaut_bus15({"Mmax": MMAX})
    assert ("GFL15", "mode_lim", 1.0) in [b[1:] for b in tr.basculements]
    M_d = tr.etat("GFL15.M_d")
    assert np.abs(M_d - M_d[0]).max() < 1e-12
    assert abs(tr.etat("bus101.v")[-1] - 1.0) < 0.01
