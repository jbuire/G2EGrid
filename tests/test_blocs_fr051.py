"""FR-051, partie serveur : nouveaux blocs elementaires (fonctions de transfert
contre les formules analytiques, saturation, implicite), sous-systemes
(aplatir, regrouper, degrouper, imbrication, spectre identique sur un banc
SMIB) et regroupement a l'interieur d'un sous-systeme."""

from __future__ import annotations

import copy

import numpy as np
import pytest
import sympy as sp

from g2egrid.elements.blocs import (
    SchemaInvalide,
    aplatir,
    compiler,
    degrouper,
    exemples,
    regrouper,
    schema_gfl,
    schema_gfm,
)
from g2egrid.reseau.tables import Tables, assembler

OMEGAS = (0.1, 1.0, 10.0, 100.0, 1000.0)


def _banc(bloc: dict, params: dict | None = None) -> dict:
    """Un schema dont la consigne P_ref traverse ``bloc`` (id ``b``)."""
    return {
        "nom": "banc",
        "parametres": params or {},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref"},
            {"id": "qr", "type": "consigne", "nom": "Q_ref"},
            {"id": "b", **bloc},
            {"id": "op", "type": "sortie", "signal": "P"},
            {"id": "oq", "type": "sortie", "signal": "Q"},
        ],
        "liens": [
            {"de": "pr", "vers": "b"},
            {"de": "b", "vers": "op"},
            {"de": "qr", "vers": "oq"},
        ],
    }


def _transfert(schema: dict, w: float, bloc_id: str = "b") -> complex:
    """``y_b / P_ref`` a ``j.w`` : lineaire, par resolution du systeme
    algebro-differentiel compile (x' = f, 0 = g)."""
    c = compiler(schema)
    xs = [sp.Symbol(n) for n in c.etats]
    ys = [sp.Symbol(n) for n in c.algebriques]
    u = sp.Symbol("P_ref")
    pts = {sp.Symbol(k): v for k, v in c.parametres.items()}
    zero = {s: 0.0 for s in xs + ys + [u]} | pts
    f, g = sp.Matrix(c.f), sp.Matrix(c.g)
    A = np.array(f.jacobian(xs).subs(zero).tolist(), dtype=complex).reshape(
        len(xs), len(xs)
    )
    Fy = np.array(f.jacobian(ys).subs(zero).tolist(), dtype=complex).reshape(
        len(xs), len(ys)
    )
    Gx = np.array(g.jacobian(xs).subs(zero).tolist(), dtype=complex).reshape(
        len(ys), len(xs)
    )
    Gy = np.array(g.jacobian(ys).subs(zero).tolist(), dtype=complex).reshape(
        len(ys), len(ys)
    )
    Fu = np.array(f.jacobian([u]).subs(zero).tolist(), dtype=complex).reshape(
        len(xs), 1
    )
    Gu = np.array(g.jacobian([u]).subs(zero).tolist(), dtype=complex).reshape(
        len(ys), 1
    )
    n = len(xs)
    M = np.block([[1j * w * np.eye(n) - A, -Fy], [-Gx, -Gy]])
    sol = np.linalg.solve(M, np.vstack([Fu, Gu]))
    return complex(sol[n + c.algebriques.index(f"y_{bloc_id}")][0])


CAS = [
    (
        {"type": "pid", "Kp": 2.0, "Ki": 30.0, "Kd": 0.05, "Tf": 0.01},
        lambda s: 2.0 + 30.0 / s + 0.05 * s / (1 + 0.01 * s),
    ),
    (
        {"type": "avance_retard", "K": 1.5, "T1": 0.2, "T2": 0.05},
        lambda s: 1.5 * (1 + 0.2 * s) / (1 + 0.05 * s),
    ),
    (
        {"type": "second_ordre", "K": 2.0, "wn": 30.0, "zeta": 0.4},
        lambda s: 2.0 * 900.0 / (s**2 + 2 * 0.4 * 30.0 * s + 900.0),
    ),
    (
        {"type": "fonction_transfert", "num": [3.0], "den": [1.0, 2.0]},
        lambda s: 3.0 / (s + 2.0),
    ),
    (
        {"type": "fonction_transfert", "num": [1.0, 4.0], "den": [2.0, 3.0, 10.0]},
        lambda s: (s + 4.0) / (2 * s**2 + 3 * s + 10.0),
    ),
    (  # ordre 3, avec transmission directe
        {
            "type": "fonction_transfert",
            "num": [0.5, 2.0, 3.0, 40.0],
            "den": [1.0, 6.0, 11.0, 6.0],
        },
        lambda s: (
            (0.5 * s**3 + 2 * s**2 + 3 * s + 40.0) / (s**3 + 6 * s**2 + 11 * s + 6.0)
        ),
    ),
    (  # coefficients en expressions des parametres
        {"type": "fonction_transfert", "num": ["k"], "den": [1, "2*a", "a^2"]},
        lambda s: 5.0 / (s + 3.0) ** 2,
    ),
    ({"type": "saturation", "min": -1.0, "max": 2.0}, lambda s: 1.0),
]


@pytest.mark.parametrize("bloc, formule", CAS, ids=[c[0]["type"] for c in CAS])
def test_fonction_de_transfert_egale_la_formule(bloc, formule):
    """La reponse en frequence de chaque bloc (resolue sur ses equations
    compilees) egale la formule analytique a 1e-12 relatif, de 0,1 a
    1000 rad/s ; la saturation, libre, vaut 1."""
    schema = _banc(bloc, {"k": 5.0, "a": 3.0})
    for w in OMEGAS:
        attendu = formule(1j * w)
        assert abs(_transfert(schema, w) - attendu) <= 1e-12 * abs(attendu)


def test_fonction_de_transfert_erreurs():
    with pytest.raises(SchemaInvalide, match="impropre"):
        compiler(_banc({"type": "fonction_transfert", "num": [1, 0, 0], "den": [1, 1]}))
    with pytest.raises(SchemaInvalide, match="liste"):
        compiler(_banc({"type": "fonction_transfert", "num": 1, "den": [1, 1]}))
    with pytest.raises(SchemaInvalide, match="strictement positif"):
        compiler(_banc({"type": "pid", "Kp": 1, "Ki": 1, "Kd": 1, "Tf": 0}))


def test_noms_d_etats_stables():
    c = compiler(_banc({"type": "fonction_transfert", "num": [1], "den": [1, 2, 3]}))
    assert c.etats[:2] == ["x_b_1", "x_b_2"]
    c = compiler(_banc({"type": "pid", "Kp": 1, "Ki": 1, "Kd": 1, "Tf": 1}))
    assert c.etats[:2] == ["x_b_1", "x_b_2"]
    c = compiler(_banc({"type": "avance_retard", "K": 1, "T1": 1, "T2": 1}))
    assert c.etats[0] == "x_b"


def test_categories_de_tous_les_types():
    from g2egrid.elements.blocs import CATEGORIES, TYPES

    cles = {k for k, _t in CATEGORIES}
    assert all(t["categorie"] in cles for t in TYPES.values())
    assert TYPES["pid"]["categorie"] == "dynamique"
    assert TYPES["port_entree"]["categorie"] == "structure"


def test_implicite_egale_le_calcul_explicite():
    """``y`` tel que ``y^3 + y - u = 0`` : meme solution que le calcul direct
    (Newton sur le residu compile), et le champ ``y`` est exige."""
    s = _banc({"type": "implicite", "expr": "y^3 + y - u", "n": 1})
    c = compiler(s)
    y = sp.Symbol("y_b")
    res = [
        e
        for e in c.g
        if y in e.free_symbols and sp.Symbol("y_pr") not in e.free_symbols
    ]
    u = 2.0
    sol = float(sp.nsolve(res[0].subs(sp.Symbol("P_ref"), u), y, 1.0))
    assert abs(sol**3 + sol - u) < 1e-12
    assert any("0 =" in e["texte"] for e in c.latex())
    with pytest.raises(SchemaInvalide, match="ne contient pas y"):
        compiler(_banc({"type": "implicite", "expr": "u - 1", "n": 1}))


def test_consigne_non_calee():
    s = _banc({"type": "gain", "K": 1})
    s["blocs"].append(
        {"id": "c3", "type": "consigne", "nom": "X", "calee": False, "valeur": 0.3}
    )
    c = compiler(s)
    assert c.consignes == ["P_ref", "Q_ref", "X"] and c.calees == ["P_ref", "Q_ref"]
    assert float(c.valeurs["X"]) == 0.3
    s["blocs"][1]["calee"] = False
    with pytest.raises(SchemaInvalide, match="calee"):
        compiler(s)


def _reseau(schema: dict) -> Tables:
    return Tables(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3, 20e3, 1.0, 0.0, 20e6, 5e6, 0.0, 0.0),
        ],
        T_Lines_pu=[(1, 1, 2, 1.0, 0.01, 0.05, 0.0, 20e3, 100e6)],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
        T_Shunts=[(1, 2, 0.01)],
        T_Blocs=[(1, 2, "src")],
        MODELES_BLOCS={"src": schema},
    )


def _source_saturee() -> dict:
    """``id = P_ref + 0,1.sat(V)`` : la saturation (bornes 0,6 et 1,5) lit le
    module de la tension du bus, que filtre ensuite un premier ordre (etat
    ``x_f``) ; ``iq = Q_ref``."""
    return {
        "nom": "sat",
        "parametres": {"tau": 0.002, "lo": 0.6, "hi": 1.5},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref"},
            {"id": "qr", "type": "consigne", "nom": "Q_ref"},
            {"id": "v", "type": "entree", "signal": "V"},
            {"id": "sat", "type": "saturation", "min": "lo", "max": "hi"},
            {"id": "f", "type": "premier_ordre", "K": 0.1, "tau": "tau"},
            {"id": "s", "type": "somme", "signes": "++"},
            {"id": "op", "type": "sortie", "signal": "id"},
            {"id": "oq", "type": "sortie", "signal": "iq"},
        ],
        "liens": [
            {"de": "v", "vers": "sat"},
            {"de": "sat", "vers": "f"},
            {"de": "pr", "vers": "s", "entree": 0},
            {"de": "f", "vers": "s", "entree": 1},
            {"de": "s", "vers": "op"},
            {"de": "qr", "vers": "oq"},
        ],
    }


def _simuler_saturation(limiteurs: bool):
    """Defaut franc au bus 3 de 0,05 a 0,10 s : (temps, x_f, basculements)."""
    from g2egrid.temporel import Manoeuvre, ModeleNonLineaire

    t = Tables(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3, 20e3, 1.0, 0.0, 4e6, 1e6, 0.0, 0.0),
            (3, 3, 20e3, 1.0, 0.0, 0.0, 0.0, 4e6, 1e6),
        ],
        T_Lines_pu=[
            (1, 1, 3, 1.0, 0.01, 0.15, 0.0, 20e3, 100e6),
            (2, 3, 2, 1.0, 0.005, 0.05, 0.0, 20e3, 100e6),
        ],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
        T_Blocs=[(1, 2, "sat", 10e6)],
        MODELES_BLOCS={"sat": _source_saturee()},
    )
    net = assembler(t)
    net.run_loadflow()
    net.solve_equilibrium()
    el = next(e for e in net.elements if e.name.startswith("Blocs"))
    m = ModeleNonLineaire(net)
    u = m.scenario(
        Manoeuvre(0.05, "defaut", "bus3", R=1e-3, X=1e-2),
        Manoeuvre(0.10, "elimination", "bus3"),
    )
    tr = m.simuler(
        (0.0, 0.16), uext_fn=u, methode="BDF2", pas=1e-4, limiteurs=limiteurs
    )
    bascule = [(round(b[0], 4), b[2], b[3]) for b in tr.basculements]
    return tr.t, tr.etat(f"{el.name}.x_f"), bascule


def test_saturation_borne_la_sortie_en_simulation():
    """Libre au point de fonctionnement (sortie = entree) ; sur un defaut, la
    tension passe sous la borne basse : la garde bascule le mode et la sortie
    reste a 0,6 (``x_f`` suit 0,1 x 0,6) ; sans les gardes, elle suit la
    tension jusqu'a presque zero ; la tension revenue, le mode relache."""
    t, x, bascule = _simuler_saturation(True)
    assert 0.0095 < x[t < 0.04].mean() < 0.11  # libre : 0,1 x V
    pendant = (t > 0.06) & (t < 0.10)
    assert x[pendant].min() >= 0.1 * 0.6 - 1e-4
    assert ("mode_sat", 1.0) in [b[1:] for b in bascule]
    assert ("mode_sat", 0.0) in [b[1:] for b in bascule]
    t2, x2, bascule2 = _simuler_saturation(False)
    assert bascule2 == []
    assert x2[(t2 > 0.06) & (t2 < 0.10)].min() < 0.1 * 0.5


# ------------------------------------------------------------ sous-systemes
def _sans_positions(s):
    if isinstance(s, dict):
        return {
            k: _sans_positions(v)
            for k, v in s.items()
            if not (k in ("x", "y") and "type" in s)
        }
    if isinstance(s, list):
        return [_sans_positions(v) for v in s]
    return s


def _canon(schema: dict) -> dict:
    """Le schema sans positions, blocs et liens tries (l'ordre des listes
    n'est pas du contenu)."""
    s = _sans_positions(copy.deepcopy(schema))

    def trie(niv):
        niv["blocs"] = sorted(niv["blocs"], key=lambda b: b["id"])
        for b in niv["blocs"]:
            if b["type"] == "sous_systeme":
                trie(b["schema"])
        niv["liens"] = sorted(
            niv["liens"],
            key=lambda lien: (
                lien["de"],
                lien["vers"],
                lien.get("entree", 0),
                lien.get("sortie", 0),
            ),
        )

    trie(s)
    return s


def _eq_texte(c) -> set[str]:
    return {sp.sstr(e) for e in c.f + c.g + c.h}


def _renomme(c, correspondance: dict[str, str]) -> set[str]:
    sub = {sp.Symbol(a): sp.Symbol(b) for a, b in correspondance.items()}
    return {sp.sstr(e.xreplace(sub)) for e in c.f + c.g + c.h}


def _mini() -> dict:
    return {
        "nom": "mini",
        "parametres": {"k": 2.0, "tau": 0.1},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref", "x": 10, "y": 10},
            {"id": "qr", "type": "consigne", "nom": "Q_ref", "x": 10, "y": 90},
            {"id": "g", "type": "gain", "K": "k", "x": 150, "y": 10},
            {
                "id": "f",
                "type": "premier_ordre",
                "K": 1,
                "tau": "tau",
                "x": 300,
                "y": 10,
            },
            {"id": "s", "type": "somme", "signes": "+-", "x": 450, "y": 10},
            {"id": "op", "type": "sortie", "signal": "id", "x": 600, "y": 10},
            {"id": "oq", "type": "sortie", "signal": "iq", "x": 600, "y": 90},
        ],
        "liens": [
            {"de": "pr", "vers": "g", "entree": 0},
            {"de": "g", "vers": "f", "entree": 0},
            {"de": "f", "vers": "s", "entree": 0},
            {"de": "pr", "vers": "s", "entree": 1},
            {"de": "s", "vers": "op", "entree": 0},
            {"de": "qr", "vers": "oq", "entree": 0},
        ],
    }


def test_regrouper_degrouper_aller_retour():
    s = _mini()
    r = regrouper(s, ["g", "f"], "filtre")
    ss = next(b for b in r["blocs"] if b["type"] == "sous_systeme")
    assert ss["nom"] == "filtre" and "__" not in ss["id"]
    types = [b["type"] for b in ss["schema"]["blocs"]]
    assert types.count("port_entree") == 1 and types.count("port_sortie") == 1
    assert _canon(degrouper(r, ss["id"])) == _canon(s)
    # le schema recu n'est pas modifie
    assert _canon(s) == _canon(_mini())
    # blocs contigus : l'ordre des listes revient aussi
    assert _sans_positions(degrouper(r, ss["id"])) == _sans_positions(s)


def test_un_port_par_signal_source():
    """``pr`` alimente deux blocs du groupe : un seul port d'entree."""
    s = _mini()
    r = regrouper(s, ["g", "s"], "x")
    ss = next(b for b in r["blocs"] if b["type"] == "sous_systeme")
    ports = [b for b in ss["schema"]["blocs"] if b["type"] == "port_entree"]
    assert len(ports) == 2  # pr (une fois) et f
    outer = [lien for lien in r["liens"] if lien["vers"] == ss["id"]]
    assert sorted((lien["de"], lien["entree"]) for lien in outer) == [
        ("f", 1),
        ("pr", 0),
    ]


def test_compiler_egale_compiler_aplati_et_regroupement_neutre():
    s = _mini()
    r = regrouper(s, ["g", "f"], "filtre", id="flt")
    c0, c1 = compiler(s), compiler(r)
    assert _eq_texte(compiler(aplatir(r))) == _eq_texte(c1)
    ren = {"y_g": "y_flt__g", "y_f": "y_flt__f", "x_f": "x_flt__f"}
    assert _renomme(c0, ren) == _eq_texte(c1)


def test_regrouper_erreurs():
    s = _mini()
    for ids, nom, msg in [
        ([], "x", "aucun bloc"),
        (["zz"], "x", "inconnu"),
        (["g", "g"], "x", "deux fois"),
        (["g"], "  ", "nom"),
    ]:
        with pytest.raises(SchemaInvalide, match=msg):
            regrouper(s, ids, nom)
    with pytest.raises(SchemaInvalide, match="existe deja"):
        regrouper(s, ["g"], "x", id="pr")
    r = regrouper(s, ["g", "f"], "filtre")
    ss = next(b for b in r["blocs"] if b["type"] == "sous_systeme")
    port = next(b for b in ss["schema"]["blocs"] if b["type"] == "port_entree")
    with pytest.raises(SchemaInvalide, match="port"):
        regrouper(ss["schema"], [port["id"]], "p")
    with pytest.raises(SchemaInvalide, match="pas un sous-systeme"):
        degrouper(s, "g")


def test_ports_hors_sous_systeme_et_double_souligne():
    s = _mini()
    s["blocs"].append({"id": "p", "type": "port_entree", "index": 0, "nom": "a"})
    with pytest.raises(SchemaInvalide, match="ne sont admis que dans un sous-systeme"):
        compiler(s)
    r = regrouper(_mini(), ["g"], "gg")
    ss = next(b for b in r["blocs"] if b["type"] == "sous_systeme")
    ss["schema"]["blocs"][1]["id"] = "a__b"
    with pytest.raises(SchemaInvalide, match="double souligne"):
        compiler(r)


def test_sortie_de_sous_systeme_non_reliee():
    r = regrouper(_mini(), ["g"], "gg")
    ss = next(b for b in r["blocs"] if b["type"] == "sous_systeme")
    ss["schema"]["liens"] = [
        lien for lien in ss["schema"]["liens"] if lien["vers"] != "s1"
    ]
    with pytest.raises(SchemaInvalide, match="n'est reliee a rien"):
        compiler(r)


def test_trois_niveaux_d_imbrication():
    s = _mini()
    r = regrouper(s, ["g", "f"], "niv1", id="n1")
    r = regrouper(r, ["g"], "niv2", id="n2", chemin=["n1"])
    r = regrouper(r, ["g"], "niv3", id="n3", chemin=["n1", "n2"])
    plat = aplatir(r)
    assert {b["id"] for b in plat["blocs"]} >= {"n1__n2__n3__g", "n1__f"}
    assert not any(b["type"] in ("sous_systeme", "port_entree") for b in plat["blocs"])
    ren = {
        "y_g": "y_n1__n2__n3__g",
        "y_f": "y_n1__f",
        "x_f": "x_n1__f",
    }
    assert _renomme(compiler(s), ren) == _eq_texte(compiler(r))
    # on remonte : degrouper de l'interieur vers l'exterieur redonne l'etat de depart
    d = degrouper(r, "n3", chemin=["n1", "n2"])
    d = degrouper(d, "n2", chemin=["n1"])
    d = degrouper(d, "n1")
    assert _canon(d) == _canon(s)


def test_profondeur_maximale():
    s = _mini()
    for k in range(9):
        s = regrouper(
            s,
            ["g"] if k == 0 else [f"n{k - 1}"],
            f"n{k}",
            id=f"n{k}",
            chemin=[f"n{j}" for j in range(k - 1, -1, -1)][::-1][:0],
        )
    # neuf niveaux : refuse
    with pytest.raises(SchemaInvalide, match="niveaux"):
        compiler(s)


def test_regrouper_dans_un_sous_systeme_par_le_niveau_recu():
    """L'interface envoie ``{nom, parametres, blocs, liens}`` du niveau
    affiche, ports compris, et remplace blocs/liens par ceux du resultat."""
    r = regrouper(_mini(), ["g", "f"], "filtre", id="flt")
    ss = next(b for b in r["blocs"] if b["type"] == "sous_systeme")
    niveau = {"nom": "mini", "parametres": {"k": 2.0, "tau": 0.1}} | copy.deepcopy(
        ss["schema"]
    )
    avant = {b["id"] for b in niveau["blocs"]}
    out = regrouper(niveau, ["g"], "gain")
    nouveaux = {b["id"] for b in out["blocs"]} - avant
    assert len(nouveaux) == 1
    assert avant - {b["id"] for b in out["blocs"]} == {"g"}
    ss["schema"] = {"blocs": out["blocs"], "liens": out["liens"]}
    ren = {"y_g": "y_flt__gain__g", "y_f": "y_flt__f", "x_f": "x_flt__f"}
    assert _renomme(compiler(_mini()), ren) == _eq_texte(compiler(r))
    # et degrouper du meme niveau
    nid = nouveaux.pop()
    back = degrouper(out, nid)
    assert _canon({"blocs": back["blocs"], "liens": back["liens"]}) == _canon(
        {"blocs": niveau["blocs"], "liens": niveau["liens"]}
    )


def _mini_flt() -> dict:
    return regrouper(_mini(), ["g", "f"], "filtre", id="flt")


# ------------------------------------------- gfl / gfm complets, banc SMIB
def _banc_smib(schema: dict, nom: str) -> Tables:
    t = Tables(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 3 if nom == "gfl" else 2, 20e3, 1.0, 0.0, 8e6, 1e6, 0.0, 0.0),
        ],
        T_Lines_pu=[(1, 1, 2, 1.0, 0.01, 0.15, 0.0, 20e3, 100e6)],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
        T_Blocs=[(1, 2, "m", 10e6)],
        MODELES_BLOCS={"m": schema},
    )
    return t


def _spectre(t: Tables) -> np.ndarray:
    net = assembler(t)
    assert net.run_loadflow().converged
    net.solve_equilibrium()
    lam = np.linalg.eigvals(net.build_state_space().A)
    return np.sort_complex(lam[np.abs(lam) > 1e-6])


@pytest.mark.parametrize(
    "cle, ids",
    [
        ("gfl_complet", ["pk", "pll", "wp", "dw", "th", "wf"]),
        ("gfl_complet", ["rl1", "cf", "rl2"]),
        ("gfl_complet", ["pe", "qe", "ep", "cp", "eqp", "cq"]),
        ("gfm_complet", ["pe", "qe", "pm", "qm", "ep", "inertie", "dw", "amort"]),
        ("gfm_complet", ["rl1", "cf", "rl2", "pinv", "kd", "kq"]),
    ],
)
def test_regrouper_gfl_gfm_garde_equations_et_spectre(cle, ids):
    s = exemples()[cle]["schema"]
    r = regrouper(s, ids, "groupe", id="grp")
    ren = {}
    for b in s["blocs"]:
        if b["id"] in ids:
            ren |= {
                f"y_{b['id']}": f"y_grp__{b['id']}",
                f"x_{b['id']}": f"x_grp__{b['id']}",
            }
            for suf in ("d", "q", "aw_d", "aw_q"):
                ren |= {
                    f"y_{b['id']}_{suf}": f"y_grp__{b['id']}_{suf}",
                    f"x_{b['id']}_{suf}": f"x_grp__{b['id']}_{suf}",
                }
    c0, c1 = compiler(s), compiler(r)
    assert _renomme(c0, ren) == _eq_texte(c1)
    nom = "gfl" if "gfl" in cle else "gfm"
    a = _spectre(_banc_smib(s, nom))
    b = _spectre(_banc_smib(r, nom))
    assert a.size == b.size
    assert np.max(np.abs(a - b) / np.maximum(1.0, np.abs(a))) < 1e-12


def test_exemples_v1_compilent_a_l_identique():
    """Un schema de la version 1 compile comme avant : memes noms, memes
    equations (les exemples, ``schema_gfl``, ``schema_gfm``)."""
    for v in exemples().values():
        c = compiler(v["schema"])
        assert _eq_texte(c) == _eq_texte(compiler(aplatir(v["schema"])))
    assert "x_rl1_d" in compiler(schema_gfl()).etats
    assert len(compiler(schema_gfm()).etats) == 14
