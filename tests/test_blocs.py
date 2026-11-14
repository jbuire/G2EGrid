"""Elements decrits par schema-bloc (:mod:`g2egrid.blocs`) : compilation en
equations, assemblage dans un reseau, equilibre au point du load flow et
valeurs propres attendues."""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from g2egrid.blocs import SchemaInvalide, classe_element, compiler
from g2egrid.tables import Tables, assembler, tables_depuis_json, tables_vers_json

TAU = 0.05


def source_filtree(tau: float = TAU) -> dict:
    """P et Q suivent leurs consignes par un premier ordre."""
    return {
        "nom": "source_filtree",
        "parametres": {"tau": tau},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref"},
            {"id": "qr", "type": "consigne", "nom": "Q_ref"},
            {"id": "fp", "type": "premier_ordre", "K": 1, "tau": "tau"},
            {"id": "fq", "type": "premier_ordre", "K": 1, "tau": "2*tau"},
            {"id": "op", "type": "sortie", "signal": "P"},
            {"id": "oq", "type": "sortie", "signal": "Q"},
        ],
        "liens": [
            {"de": "pr", "vers": "fp"},
            {"de": "qr", "vers": "fq"},
            {"de": "fp", "vers": "op"},
            {"de": "fq", "vers": "oq"},
        ],
    }


def source_asservie() -> dict:
    """Asservissement de P en boucle fermee : ``P = x``, PI sur
    ``P_ref - P_mes``, avec une boucle algebrique (gain sur la somme)."""
    return {
        "nom": "source_pi",
        "parametres": {"Kp": 0.5, "Ki": 20.0, "tm": 0.02},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref"},
            {"id": "qr", "type": "consigne", "nom": "Q_ref"},
            {"id": "e", "type": "somme", "signes": "+-"},
            {"id": "c", "type": "pi", "Kp": "Kp", "Ki": "Ki"},
            {"id": "m", "type": "premier_ordre", "K": 1, "tau": "tm"},
            {"id": "op", "type": "sortie", "signal": "P"},
            {"id": "oq", "type": "sortie", "signal": "Q"},
        ],
        "liens": [
            {"de": "pr", "vers": "e", "entree": 0},
            {"de": "m", "vers": "e", "entree": 1},
            {"de": "e", "vers": "c"},
            {"de": "c", "vers": "m"},
            {"de": "c", "vers": "op"},
            {"de": "qr", "vers": "oq"},
        ],
    }


def reseau(schema: dict, P: float = 20e6, Q: float = 5e6) -> Tables:
    """Un bus infini (bus 1) relie par une ligne au bus 2, ou l'element a
    schema-bloc injecte ``P``, ``Q`` (bus PQ, type 3 ; le slack est le
    type 1)."""
    return Tables(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3, 20e3, 1.0, 0.0, P, Q, 0.0, 0.0),
        ],
        T_Lines_pu=[(1, 1, 2, 1.0, 0.01, 0.05, 0.0, 20e3, 100e6)],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
        # un condensateur au bus 2 : sans lui, le noeud n'a que sa resistance
        # de causalite, et une injection a puissance constante y est instable
        T_Shunts=[(1, 2, 0.01)],
        T_Blocs=[(1, 2, "src")],
        MODELES_BLOCS={"src": schema},
    )


def _spectre(t: Tables):
    net = assembler(t)
    assert net.run_loadflow().converged
    net.solve_equilibrium()
    ss = net.build_state_space()
    return net, ss, np.linalg.eigvals(ss.A)


def test_compilation_premier_ordre():
    c = compiler(source_filtree())
    # + le filtre de mesure de tension d'une sortie en puissance
    assert c.etats == ["x_fp", "x_fq", "vd_mes", "vq_mes"]
    assert c.parametres["tau_mes"] > 0
    assert c.consignes == ["P_ref", "Q_ref"]
    assert set(c.algebriques) == {"y_fp", "y_fq"}
    tex = [ligne["latex"] for ligne in c.latex()]
    assert any(r"\dot{x_{fp}}" in s for s in tex)


def test_source_filtree_equilibre_et_valeurs_propres():
    _net, _ss, lam = _spectre(reseau(source_filtree()))
    # les deux poles du schema, -1/tau et -1/(2 tau), sont dans le spectre
    for p in (-1 / TAU, -1 / (2 * TAU)):
        assert np.min(np.abs(lam - p)) < 1e-6 * abs(p)
    assert lam.real.max() < 1e-9


def test_source_de_courant():
    """Sorties ``id``/``iq`` directes : pas de filtre de mesure."""
    s = source_filtree()
    s["blocs"][4]["signal"], s["blocs"][5]["signal"] = "id", "iq"
    c = compiler(s)
    assert c.etats == ["x_fp", "x_fq"] and "tau_mes" not in c.parametres
    _net, _ss, lam = _spectre(reseau(s))
    for p in (-1 / TAU, -1 / (2 * TAU)):
        assert np.min(np.abs(lam - p)) < 1e-6 * abs(p)


def test_fonction_et_expression():
    s = source_filtree()
    s["blocs"].append({"id": "g", "type": "fonction", "expr": "tau*u^2", "n": 1})
    s["blocs"].append({"id": "k", "type": "gain", "K": "1/tau"})
    s["liens"] += [{"de": "fp", "vers": "g"}, {"de": "g", "vers": "k"}]
    c = compiler(s)
    i = c.algebriques.index("y_k")
    y_k, y_g, y_fp, tau = sp.symbols("y_k y_g y_fp tau")
    assert sp.simplify(c.g[i] - (y_k - y_g / tau)) == 0
    j = c.algebriques.index("y_g")
    assert sp.simplify(c.g[j] - (y_g - tau * y_fp**2)) == 0


def test_source_asservie_boucle_fermee():
    s = source_asservie()
    _net, _ss, lam = _spectre(reseau(s))
    # boucle fermee : tm.s^2 + (1 + Kp).s + Ki = 0
    Kp, Ki, tm = (s["parametres"][k] for k in ("Kp", "Ki", "tm"))
    attendus = np.roots([tm, 1 + Kp, Ki])
    for p in attendus:
        assert np.min(np.abs(lam - p)) < 1e-6 * abs(p)


def test_parametre_surcharge_par_bus():
    t = reseau(source_filtree())
    t.PARAMS_PAR_BUS[2] = {"tau": 0.2}
    _net, _ss, lam = _spectre(t)
    assert np.min(np.abs(lam + 5.0)) < 1e-6


def test_aller_retour_json():
    t = reseau(source_filtree())
    t2 = tables_depuis_json(tables_vers_json(t))
    assert t2.T_Blocs == t.T_Blocs
    assert t2.MODELES_BLOCS == t.MODELES_BLOCS


def test_classe_mise_en_cache():
    assert classe_element(source_filtree()) is classe_element(source_filtree())
    assert classe_element(source_filtree(0.1)) is not classe_element(source_filtree())


@pytest.mark.parametrize(
    "modif, message",
    [
        (lambda s: s["liens"].pop(0), "non reliee"),
        (lambda s: s["blocs"].pop(1), "bloc inconnu"),
        (lambda s: s["blocs"][2].update(tau="beta"), "parametres du schema"),
        (lambda s: s["blocs"][4].update(signal="id"), "deux sorties"),
        (lambda s: s["blocs"][0].update(type="truc"), "inconnu"),
    ],
)
def test_erreurs_de_compilation(modif, message):
    s = source_filtree()
    modif(s)
    with pytest.raises(SchemaInvalide, match=message):
        compiler(s)


def test_modele_absent():
    t = reseau(source_filtree())
    t.MODELES_BLOCS = {}
    with pytest.raises(ValueError, match="MODELES_BLOCS"):
        assembler(t)


def test_api_compiler_et_erreur_422():
    pytest.importorskip("fastapi")
    from fastapi import HTTPException

    from g2egrid_api import app as routes
    from g2egrid_api.schemas import CorpsBlocs

    d = routes.blocs_compiler(CorpsBlocs(schema_bloc=source_filtree()))
    assert d["consignes"] == ["P_ref", "Q_ref"]
    assert any(r"\dot" in e["latex"] for e in d["equations"])
    assert d["entrees_reseau"] == ["vd_t", "vq_t"]
    assert "premier_ordre" in routes.blocs_types()["types"]
    s = source_filtree()
    s["liens"].pop()
    with pytest.raises(HTTPException) as err:
        routes.blocs_compiler(CorpsBlocs(schema_bloc=s))
    assert err.value.status_code == 422


def test_topologie_et_validation_voient_l_element():
    from g2egrid_api import analyse

    t = reseau(source_filtree())
    bus2 = next(b for b in analyse.topologie(t)["bus"] if b["id"] == 2)
    assert "Blocs" in [e["famille"] for e in bus2["equipements"]]
    assert not [p for p in analyse.valider(t) if p["gravite"] == "erreur"]


# ------------------------------------------------- GFL complet, bloc a bloc
def _reseau_gfl(element: str, X: float = 0.15) -> Tables:
    """Bus infini, ligne de reactance ``X``, GFL de 10 MVA a 8 MW au bus 2 :
    ``element`` vaut ``"GFL1"`` (le modele du moteur) ou ``"blocs"`` (le meme
    en schema-bloc)."""
    from g2egrid.blocs import schema_gfl

    t = Tables(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3, 20e3, 1.0, 0.0, 8e6, 1e6, 0.0, 0.0),
        ],
        T_Lines_pu=[(1, 1, 2, 1.0, 0.01, X, 0.0, 20e3, 100e6)],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
    )
    if element == "blocs":
        t.T_Blocs = [(1, 2, "gfl", 10e6)]
        t.MODELES_BLOCS = {"gfl": schema_gfl()}
    else:
        t.T_GFL1s = [(1, 2, 10e6, 40e3)]
    return t


@pytest.mark.parametrize("X", [0.15, 0.6])
def test_gfl_en_blocs_egale_gfl1(X):
    """Le GFL dessine bloc a bloc (LCL, PLL, boucles de courant et de
    puissance) a EXACTEMENT les valeurs propres du GFL1 du moteur, sur un
    reseau fort comme sur un reseau faible."""
    lams = []
    for element in ("GFL1", "blocs"):
        _net, _ss, lam = _spectre(_reseau_gfl(element, X))
        lams.append(np.sort_complex(lam[np.abs(lam) > 1e-6]))
    a, b = lams
    assert a.size == b.size
    assert np.max(np.abs(a - b) / np.maximum(1.0, np.abs(a))) < 1e-7


def test_gfl_en_blocs_structure():
    from g2egrid.blocs import schema_gfl

    c = compiler(schema_gfl())
    assert len(c.etats) == 12  # comme GFL1
    assert c.consignes == ["P_ref", "Q_ref"]
    assert c.reserves == ["wb", "kS"]
    assert set(c.angles) == {"x_th", "y_th"}
    assert "y_pk_d" in c.algebriques and "x_rl1_q" in c.etats


def test_park_puis_park_inverse_rend_l_entree():
    s = {
        "nom": "aller_retour",
        "blocs": [
            {"id": "a", "type": "consigne", "nom": "id0"},
            {"id": "b", "type": "consigne", "nom": "iq0"},
            {"id": "t", "type": "entree", "signal": "theta_g"},
            {"id": "p", "type": "park"},
            {"id": "pi_", "type": "park_inverse"},
            {"id": "od", "type": "sortie", "signal": "id"},
            {"id": "oq", "type": "sortie", "signal": "iq"},
        ],
        "liens": [
            {"de": "a", "vers": "p", "entree": 0},
            {"de": "b", "vers": "p", "entree": 1},
            {"de": "t", "vers": "p", "entree": 2},
            {"de": "p", "vers": "pi_", "entree": 0},
            {"de": "p", "vers": "pi_", "entree": 1, "sortie": 1},
            {"de": "t", "vers": "pi_", "entree": 2},
            {"de": "pi_", "vers": "od"},
            {"de": "pi_", "vers": "oq", "sortie": 1},
        ],
    }
    c = compiler(s)
    sol = sp.solve(c.g, [sp.Symbol(n) for n in c.algebriques], dict=True)[0]
    id0, iq0 = sp.symbols("id0 iq0")
    assert sp.simplify(c.h[0].xreplace(sol) - id0) == 0
    assert sp.simplify(c.h[1].xreplace(sol) - iq0) == 0


def test_erreurs_des_blocs_vectoriels():
    from g2egrid.blocs import schema_gfl

    s = schema_gfl()
    s["liens"][0]["sortie"] = 1  # vd n'a qu'une sortie
    with pytest.raises(SchemaInvalide, match="sortie 1 inexistante"):
        compiler(s)
    s = schema_gfl()
    s["parametres"]["wb"] = 314.0
    with pytest.raises(SchemaInvalide, match="reserve"):
        compiler(s)


def test_tous_les_exemples_compilent():
    from g2egrid.blocs import exemples

    ex = exemples()
    assert {
        "source_filtree",
        "asservissement",
        "statisme_qv",
        "gfl_complet",
        "gfm_complet",
    } <= set(ex)
    for v in ex.values():
        compiler(v["schema"])


# ------------------------------------------- GFM complet (VSM), bloc a bloc
def _deux_spectres(construire) -> tuple[np.ndarray, np.ndarray]:
    out = []
    for element in ("moteur", "blocs"):
        net, _ss, lam = _spectre(construire(element))
        out.append((net, np.sort_complex(lam[np.abs(lam) > 1e-6])))
    return out


def _ecart(a: np.ndarray, b: np.ndarray) -> float:
    assert a.size == b.size
    return float(np.max(np.abs(a - b) / np.maximum(1.0, np.abs(a))))


def _reseau_gfm(element: str, X: float = 0.15, typ: int = 2) -> Tables:
    """Bus infini, ligne de reactance ``X``, GFM de 10 MVA a 8 MW au bus 2
    (PV si ``typ = 2``) : le GFM2 du moteur, ou le meme en schema-bloc."""
    from g2egrid.blocs import schema_gfm

    t = Tables(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, typ, 20e3, 1.01, 0.0, 8e6, 1e6, 0.0, 0.0),
        ],
        T_Lines_pu=[(1, 1, 2, 1.0, 0.01, X, 0.0, 20e3, 100e6)],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
    )
    if element == "blocs":
        t.T_Blocs = [(1, 2, "gfm", 10e6)]
        t.MODELES_BLOCS = {"gfm": schema_gfm()}
    else:
        t.T_GFM2s = [(1, 2, 10e6, 40e3)]
    return t


@pytest.mark.parametrize("X, typ", [(0.15, 2), (0.6, 2), (0.15, 3)])
def test_gfm_en_blocs_egale_gfm2(X, typ):
    """Le GFM dessine bloc a bloc (VSM, statisme Q-V, boucles de tension et
    de courant, LCL) a EXACTEMENT les valeurs propres du GFM2 du moteur."""
    (_n1, a), (_n2, b) = _deux_spectres(lambda e: _reseau_gfm(e, X, typ))
    assert _ecart(a, b) < 1e-7


def test_gfm_en_blocs_forme_un_ilot_seul():
    """Pose au bus slack d'un ilot sans autre source, le GFM a schema-bloc
    porte le repere du reseau (sorties theta et w) et tient l'ilot : memes
    valeurs propres qu'avec le GFM2."""
    from g2egrid.blocs import schema_gfm

    def ilot(element):
        t = Tables(
            Y_DER=[
                (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                (2, 3, 20e3, 1.0, 0.0, 0.0, 0.0, 6e6, 2e6),
            ],
            T_Lines_pu=[(1, 1, 2, 1.0, 0.01, 0.05, 0.0, 20e3, 100e6)],
        )
        if element == "blocs":
            t.T_Blocs = [(1, 1, "gfm", 10e6)]
            t.MODELES_BLOCS = {"gfm": schema_gfm()}
        else:
            t.T_GFM2s = [(1, 1, 10e6, 40e3)]
        return t

    (_n1, a), (net, b) = _deux_spectres(ilot)
    assert _ecart(a, b) < 1e-7
    assert b.real.max() < 0
    k, _ty = net._element_reference()
    assert net.all_elements[k].name.startswith("Blocs")
    assert net.all_elements[k] in net.formeurs()


def test_gfm_en_blocs_structure():
    from g2egrid.blocs import schema_gfm

    c = compiler(schema_gfm())
    assert len(c.etats) == 14  # comme GFM2
    assert c.consignes == ["P_ref", "E_ref"]
    assert c.services == ["theta_out", "w"]
    p = c.parametres
    assert p["H"] == pytest.approx(1 / (2 * 0.001 * 0.1 * 2 * np.pi * 50))
    assert p["D"] == pytest.approx(1000.0)


# ------------------------------------------------ limitation de courant
def _reseau_defaut(genre: str, element: str) -> Tables:
    """Bus infini, deux lignes, convertisseur de 10 MVA au bus 2 : le GFL1 /
    GFM2 du moteur, ou le meme en schema-bloc (limiteur compris)."""
    from g2egrid.blocs import schema_gfl, schema_gfm

    t = Tables(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 2 if genre == "gfm" else 3, 20e3, 1.01, 0.0, 8e6, 1e6, 0.0, 0.0),
            (3, 3, 20e3, 1.0, 0.0, 0.0, 0.0, 4e6, 1e6),
        ],
        T_Lines_pu=[
            (1, 1, 3, 1.0, 0.01, 0.15, 0.0, 20e3, 100e6),
            (2, 3, 2, 1.0, 0.005, 0.05, 0.0, 20e3, 100e6),
        ],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
    )
    if element == "blocs":
        t.T_Blocs = [(1, 2, genre, 10e6)]
        t.MODELES_BLOCS = {genre: (schema_gfl if genre == "gfl" else schema_gfm)()}
    elif genre == "gfl":
        t.T_GFL1s = [(1, 2, 10e6, 40e3)]
    else:
        t.T_GFM2s = [(1, 2, 10e6, 40e3)]
    return t


def _simuler_defaut(
    genre: str,
    element: str,
    limiteurs: bool = True,
    man=(),
    reglages: dict | None = None,
    pas: float = 1e-4,
    bus: str = "bus3",
):
    """Defaut au ``bus`` de 0,05 a 0,10 s (BDF2, pas ``pas``) : (temps, |i|
    du convertisseur, basculements). ``reglages`` s'applique au convertisseur
    (PARAMS_PAR_BUS du schema-bloc, famille du moteur)."""
    from g2egrid.temporel import Manoeuvre, ModeleNonLineaire

    t = _reseau_defaut(genre, element)
    if reglages:
        if element == "blocs":
            t.PARAMS_PAR_BUS[2] = dict(reglages)
        elif genre == "gfm":
            t.GFM_PARAMS = dict(t.GFM_PARAMS, **reglages)
        else:
            t.GFL_PARAMS = dict(t.GFL_PARAMS, **reglages)
    net = assembler(t)
    net.run_loadflow()
    net.solve_equilibrium()
    el = next(e for e in net.elements if e.name.startswith(("Blocs", "GFL", "GFM")))
    m = ModeleNonLineaire(net)
    u = m.scenario(
        Manoeuvre(0.05, "defaut", bus, R=1e-3, X=1e-2),
        Manoeuvre(0.10, "elimination", bus),
        *(Manoeuvre(t, g, el.name) for t, g in man),
    )
    tr = m.simuler((0.0, 0.16), uext_fn=u, methode="BDF2", pas=pas, limiteurs=limiteurs)
    d, q = ("x_rl1_d", "x_rl1_q") if element == "blocs" else ("isd", "isq")
    i = np.hypot(tr.etat(f"{el.name}.{d}"), tr.etat(f"{el.name}.{q}"))
    bascule = [(round(b[0], 6), b[2], b[3]) for b in tr.basculements]
    return tr.t, i, bascule


@pytest.mark.parametrize("genre, tol", [("gfl", 1e-6), ("gfm", 1e-4)])
def test_limitation_des_blocs_comme_le_moteur(genre, tol):
    """Sur un defaut, le limiteur du schema-bloc bascule aux memes instants
    que celui du GFL1 / GFM2 du moteur, et le courant du convertisseur suit
    la meme trajectoire, limitee a Imax = 1,2 pu."""
    # la saturation de tension des schemas-blocs est reglee (Mmax = 1,15) ;
    # celle du moteur ne l'est pas par defaut : on compare sans elle
    t_m, i_m, b_m = _simuler_defaut(genre, "moteur")
    t_b, i_b, b_b = _simuler_defaut(genre, "blocs", reglages={"Mmax": 0.0})
    assert b_b == b_m and ("mode_lim", 1.0) in [b[1:] for b in b_b]
    assert t_b.size == t_m.size and np.abs(i_b - i_m).max() < tol
    pendant = (t_b > 0.051) & (t_b < 0.10)
    assert i_b[pendant].max() < 1.3


def test_sans_limiteur_le_gfm_en_blocs_depasse():
    t, i, bascule = _simuler_defaut(
        "gfm", "blocs", limiteurs=False, reglages={"Mmax": 0.0}
    )
    assert bascule == []
    assert i[(t > 0.051) & (t < 0.10)].max() > 2.0


def test_blocage_d_un_convertisseur_en_blocs():
    """La manoeuvre ``blocage`` du moteur trouve le limiteur du schema-bloc
    (drapeau ``mode_lim``, parametre ``Imax``) : le courant s'eteint."""
    t, i, _b = _simuler_defaut(
        "gfm", "blocs", man=[(0.12, "blocage")], reglages={"Mmax": 0.0}
    )
    assert i[t > 0.15].max() < 0.02


def test_limiteur_et_pi_aw_compiles():
    from g2egrid.blocs import classe_element, schema_gfl

    c = compiler(schema_gfl())
    assert [d for d, _m, _i in c.limiteurs] == ["mode_lim", "mode_vsat"]
    assert c.parametres["mode_lim"] == 0.0 and c.parametres["Imax"] == 1.2
    assert c.parametres["mode_vsat"] == 0.0 and c.parametres["Mmax"] == 1.15
    assert {"y_lim_aw_d", "y_vsat_aw_q"} <= set(c.algebriques)
    cls = classe_element(schema_gfl())
    assert cls.modes == ("mode_lim", "mode_vsat")


# ------------------------------------------ saturation de tension (GFM, GFL)
VSAT = {"Mmax": 1.15}  # le retro-calcul rapide est le defaut du moteur


def test_saturation_de_tension_du_gfl_en_blocs():
    """Defaut au bus 3, elimine a 0,10 s : le GFL a schema-bloc part en
    limitation de courant pendant le defaut ; a l'elimination, la tension
    revient et la boucle de courant demande plus que |m| <= 1,15 : la
    modulation sature quelques millisecondes, puis relache. Le GFL1 du moteur
    regle de meme bascule aux MEMES instants et son courant suit la meme
    trajectoire."""
    t_b, i_b, b_b = _simuler_defaut("gfl", "blocs")
    t_m, i_m, b_m = _simuler_defaut("gfl", "moteur", reglages=VSAT)
    vsat = [b for b in b_b if b[1] == "mode_vsat"]
    assert vsat[0][1:] == ("mode_vsat", 1.0) and 0.10 <= vsat[0][0] < 0.101
    assert vsat[-1][2] == 0.0  # elle relache
    assert ("mode_lim", 1.0) in [b[1:] for b in b_b if b[0] < 0.10]
    assert [b[1:] for b in b_b] == [b[1:] for b in b_m]
    assert max(abs(a[0] - b[0]) for a, b in zip(b_b, b_m)) < 1e-5
    grille = np.linspace(0.0, 0.16, 8000)
    ecart = np.interp(grille, t_b, i_b) - np.interp(grille, t_m, i_m)
    assert np.abs(ecart).max() < 1e-6
    assert i_b[(t_b > 0.051) & (t_b < 0.10)].max() < 1.3
    assert i_b[-1] < 1.0  # revenu au point de fonctionnement


def test_saturation_de_tension_du_gfm_en_blocs():
    """Defaut franc aux bornes du GFM : la modulation sature d'abord (|m| <=
    1,15), la limitation de courant prend le relais, la saturation relache ;
    le courant reste tenu pendant et apres le defaut. Le GFM2 du moteur regle
    de meme bascule aux MEMES instants (a 1 us pres, 33 basculements) et son
    courant suit la meme trajectoire."""
    t_b, i_b, b_b = _simuler_defaut("gfm", "blocs", pas=2e-5, bus="bus2")
    t_m, i_m, b_m = _simuler_defaut(
        "gfm", "moteur", reglages=VSAT, pas=2e-5, bus="bus2"
    )
    vsat = [b for b in b_b if b[1] == "mode_vsat"]
    assert vsat[0] == (0.05, "mode_vsat", 1.0)
    assert ("mode_vsat", 0.0) in [b[1:] for b in vsat]
    assert [b[1:] for b in b_b] == [b[1:] for b in b_m]
    assert max(abs(a[0] - b[0]) for a, b in zip(b_b, b_m)) < 1e-5
    grille = np.linspace(0.0, 0.16, 8000)
    ecart = np.interp(grille, t_b, i_b) - np.interp(grille, t_m, i_m)
    assert np.abs(ecart).max() < 5e-3
    assert i_b[(t_b > 0.051) & (t_b < 0.10)].max() < 1.3
    assert i_b[t_b > 0.13].max() < 1.25


def test_retro_calcul_lent_colle_la_saturation():
    """Avec le retro-calcul de 20 ms du moteur, l'integrateur de courant
    s'emballe : la saturation ne relache plus et le courant monte bien
    au-dela de Imax. D'ou Taw_v = tau_CL/2 dans le schema-bloc."""
    t, i, bascule = _simuler_defaut("gfm", "blocs", pas=2e-5, reglages={"Taw_v": 20e-3})
    assert [b for b in bascule if b[1] == "mode_vsat"] == [(0.05, "mode_vsat", 1.0)]
    assert i[t > 0.13].max() > 2.0
