# ruff: noqa: C408
"""FR-051, criteres 1 et 2 : compatibilite de la version 1 et blocs elementaires.

Tests ecrits A L'AVEUGLE, depuis le seul cahier
``.claude/equipe/cahiers/FR-051_editeur_blocs.md`` : ils visent le contrat
(``TYPES``, ``compiler``, ``classe_element``, noms de champs), pas
l'implementation.

CRITERE 1 (compatibilite). Les cinq schemas de ``exemples()`` tels que la
version 1 les compile (figes dans ``test_fr051_donnees_v1``, avec l'empreinte
de leurs equations f, g, h et de leurs noms) doivent compiler a l'identique :
meme empreinte a 1e-11 (l'empreinte somme des valeurs a ~1e-16 pres ; 1e-11
laisse passer un reordonnancement de sommes, pas une equation modifiee).
Chaque type porte une ``categorie`` parmi les sept du cahier.

CRITERE 2 (blocs elementaires). La fonction de transfert de ``pid``,
``avance_retard``, ``second_ordre``, ``fonction_transfert`` (quatre jeux de
coefficients dont deux d'ordre 3) et ``saturation`` (libre) est calculee sur
le modele LINEARISE d'un banc a un seul bloc (jacobiens des equations
compilees, resolution en 40 chiffres) et comparee a la formule ANALYTIQUE
ecrite a la main, en 5 pulsations de 0,1 a 1000 rad/s, a 1e-12 relatif (le
seuil du cahier ; la formule en complexe double est exacte a ~1e-16, la
resolution est en precision etendue : la marge est de quatre ordres de
grandeur). ``saturation`` borne la sortie en simulation temporelle ;
``implicite`` donne la solution d'un calcul explicite.

Lectures du cahier retenues :

- le nombre d'entrees de ``implicite`` se donne dans le champ ``n`` (comme
  ``fonction``), le residu s'ecrit en ``u1..un`` et ``y`` ;
- ``Tf <= 0`` (pid) et ``T2 <= 0`` (avance_retard), NUMERIQUES, sont refuses
  (``SchemaInvalide``) ; un degre du numerateur superieur a celui du
  denominateur aussi ;
- ``port_entree``, ``port_sortie`` et ``sous_systeme`` sont de categorie
  ``structure`` ; ``sortie`` de categorie ``sorties``.
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp
from g2egrid.tables import Tables, assembler
from test_fr051_aide import (
    banc_bloc,
    copie,
    empreinte,
    erreur_relative,
    fdt,
    memes_empreintes,
    valeurs_statiques,
)
from test_fr051_donnees_v1 import EMPREINTES_V1, SCHEMAS_V1

from g2egrid.elements import blocs as B
from g2egrid.elements.blocs import SchemaInvalide, classe_element, compiler

OMEGAS = (0.1, 1.0, 10.0, 100.0, 1000.0)
TOL = 1e-12  # cahier, critere 2
CATEGORIES = {
    "sources",
    "operations",
    "dynamique",
    "repere_dq",
    "limiteurs",
    "structure",
    "sorties",
}
TYPES_V1 = {
    "entree": 0,
    "constante": 0,
    "consigne": 0,
    "gain": 1,
    "produit": 2,
    "division": 2,
    "integrateur": 1,
    "premier_ordre": 1,
    "pi": 1,
    "derivee": 1,
    "park": 3,
    "park_inverse": 3,
    "rl_dq": 5,
    "c_dq": 5,
    "limiteur_courant": 2,
    "saturation_tension": 2,
    "pi_aw": 2,
    "sortie": 1,
}  # + somme et fonction (entrees lues dans le bloc) : 20 types


# ============================================================ critere 1
@pytest.mark.parametrize("nom", sorted(SCHEMAS_V1))
def test_schema_v1_compile_a_l_identique(nom):
    """Un schema de la version 1 (fige) compile aux memes equations et aux
    memes noms de variables qu'en version 1."""
    c = compiler(copie(SCHEMAS_V1[nom]))
    assert memes_empreintes(empreinte(c), EMPREINTES_V1[nom])


def test_les_types_de_la_version_1_restent():
    for t, n in TYPES_V1.items():
        assert t in B.TYPES, t
        assert B.TYPES[t]["entrees"] == n, t
    assert {"somme", "fonction"} <= set(B.TYPES)


def test_chaque_type_a_une_categorie_du_cahier():
    for nom, t in B.TYPES.items():
        assert t.get("categorie") in CATEGORIES, f"{nom} : {t.get('categorie')!r}"
    # les sept categories sont toutes employees (la palette est construite
    # depuis ces donnees : une categorie vide serait un trou dans l'interface)
    assert {t["categorie"] for t in B.TYPES.values()} == CATEGORIES


def test_nouveaux_types_declares():
    for nom in (
        "pid",
        "avance_retard",
        "second_ordre",
        "fonction_transfert",
        "saturation",
        "implicite",
        "sous_systeme",
        "port_entree",
        "port_sortie",
    ):
        assert nom in B.TYPES, nom
    for nom in ("sous_systeme", "port_entree", "port_sortie"):
        assert B.TYPES[nom]["categorie"] == "structure", nom  # lecture du cahier
    assert B.TYPES["sortie"]["categorie"] == "sorties"
    assert B.TYPES["pid"]["entrees"] == 1
    assert B.TYPES["saturation"]["entrees"] == 1


def test_les_exemples_v2_compilent_et_gardent_leurs_tailles():
    ex = B.exemples()
    assert {
        "source_filtree",
        "asservissement",
        "statisme_qv",
        "gfl_complet",
        "gfm_complet",
    } <= set(ex)
    taille = {"gfl_complet": (12, 2), "gfm_complet": (14, 2)}
    for cle, v in ex.items():
        c = compiler(v["schema"])
        if cle in taille:
            assert (len(c.etats), len(c.consignes)) == taille[cle], cle


# ============================================================ critere 2
def _pid(Kp, Ki, Kd, Tf):
    return lambda s: Kp + Ki / s + Kd * s / (1 + Tf * s)


def _avance_retard(K, T1, T2):
    return lambda s: K * (1 + T1 * s) / (1 + T2 * s)


def _second_ordre(K, wn, zeta):
    return lambda s: K * wn**2 / (s**2 + 2 * zeta * wn * s + wn**2)


def _poly(coef):
    return lambda s: sum(c * s ** (len(coef) - 1 - k) for k, c in enumerate(coef))


def _ft(num, den):
    return lambda s: _poly(num)(s) / _poly(den)(s)


CAS_FDT = [
    # (id, bloc, parametres du schema, formule analytique, nombre d'etats)
    (
        "pid_nombres",
        dict(type="pid", Kp=2.0, Ki=3.0, Kd=0.5, Tf=0.01),
        {},
        _pid(2.0, 3.0, 0.5, 0.01),
        2,
    ),
    (
        "pid_parametres",
        dict(type="pid", Kp="kp", Ki="ki", Kd="kd", Tf="tf"),
        dict(kp=0.7, ki=40.0, kd=0.02, tf=0.05),
        _pid(0.7, 40.0, 0.02, 0.05),
        2,
    ),
    (
        "pid_sans_integrale",
        dict(type="pid", Kp=1.5, Ki=0.0, Kd=0.3, Tf=0.02),
        {},
        _pid(1.5, 0.0, 0.3, 0.02),
        2,
    ),
    (
        "avance_retard_a",
        dict(type="avance_retard", K=2.0, T1=0.5, T2=0.1),
        {},
        _avance_retard(2.0, 0.5, 0.1),
        1,
    ),
    (
        "avance_retard_b",
        dict(type="avance_retard", K="k", T1="t1", T2="t2"),
        dict(k=1.0, t1=0.01, t2=0.2),
        _avance_retard(1.0, 0.01, 0.2),
        1,
    ),
    (
        "avance_retard_retard_pur",
        dict(type="avance_retard", K=0.5, T1=0.0, T2=1.0),
        {},
        _avance_retard(0.5, 0.0, 1.0),
        1,
    ),
    (
        "second_ordre_sous_amorti",
        dict(type="second_ordre", K=1.0, wn=10.0, zeta=0.3),
        {},
        _second_ordre(1.0, 10.0, 0.3),
        2,
    ),
    (
        "second_ordre_amorti",
        dict(type="second_ordre", K=2.5, wn=50.0, zeta=0.7),
        {},
        _second_ordre(2.5, 50.0, 0.7),
        2,
    ),
    (
        "second_ordre_sur_amorti",
        dict(type="second_ordre", K="k", wn="wn", zeta="z"),
        dict(k=0.8, wn=3.0, z=1.5),
        _second_ordre(0.8, 3.0, 1.5),
        2,
    ),
    (
        "ft_ordre_2",
        dict(type="fonction_transfert", num=[1.0], den=[1.0, 2.0, 1.0]),
        {},
        _ft([1.0], [1.0, 2.0, 1.0]),
        2,
    ),
    (
        "ft_ordre_2_biproper",
        dict(type="fonction_transfert", num=[2.0, 1.0, 0.5], den=[0.5, 3.0, 4.0]),
        {},
        _ft([2.0, 1.0, 0.5], [0.5, 3.0, 4.0]),
        2,
    ),
    (
        "ft_ordre_3",
        dict(type="fonction_transfert", num=[1.0, 2.0, 3.0], den=[1.0, 6.0, 11.0, 6.0]),
        {},
        _ft([1.0, 2.0, 3.0], [1.0, 6.0, 11.0, 6.0]),
        3,
    ),
    (
        "ft_ordre_3_expressions",
        dict(
            type="fonction_transfert",
            num=["K*T", "K"],
            den=["T^3", "3*T^2", "3*T", "1"],
        ),
        dict(K=4.0, T=0.2),
        _ft([0.8, 4.0], [0.008, 0.12, 0.6, 1.0]),
        3,
    ),
    (
        "ft_ordre_3_non_unitaire",
        dict(
            type="fonction_transfert",
            num=[0.0, 5.0, 1.0, 20.0],
            den=[2.0, 1.0, 10.0, 100.0],
        ),
        {},
        _ft([0.0, 5.0, 1.0, 20.0], [2.0, 1.0, 10.0, 100.0]),
        3,
    ),
    (
        "saturation_libre",
        dict(type="saturation", min=-0.5, max=0.5),
        {},
        lambda s: 1.0 + 0.0 * s,
        0,
    ),
    (
        "saturation_libre_parametres",
        dict(type="saturation", min="lo", max="hi"),
        dict(lo=-3.0, hi=7.0),
        lambda s: 1.0 + 0.0 * s,
        0,
    ),
]


@pytest.mark.parametrize(
    "bloc, params, formule, n_etats",
    [pytest.param(*c[1:], id=c[0]) for c in CAS_FDT],
)
def test_fonction_de_transfert_analytique(bloc, params, formule, n_etats):
    """G(jw) du modele linearise = formule analytique, a 1e-12 relatif, en
    cinq pulsations de 0,1 a 1000 rad/s ; et le bloc a le nombre d'etats
    annonce (pid 2, avance_retard 1, second_ordre 2, fonction_transfert deg
    den, saturation 0)."""
    schema = banc_bloc(bloc, params)
    assert len(compiler(copie(schema)).etats) == n_etats
    for w in OMEGAS:
        ref = complex(formule(1j * w))
        assert erreur_relative(fdt(schema, w), ref) < TOL, f"w = {w}"


@pytest.mark.parametrize("cas", ["pid_nombres", "second_ordre_amorti", "ft_ordre_3"])
def test_noms_des_etats_stables_et_distincts(cas):
    """Les noms des etats d'un bloc a plusieurs etats suivent une convention
    stable (lecture du cahier : ils commencent par ``x_<id>``, sont distincts,
    et sont les memes d'une compilation a l'autre)."""
    bloc, params = next((c[1], c[2]) for c in CAS_FDT if c[0] == cas)
    a = compiler(banc_bloc(bloc, params)).etats
    b = compiler(banc_bloc(bloc, params)).etats
    assert a == b and len(set(a)) == len(a) > 1
    assert all(n.startswith("x_b") for n in a), a


@pytest.mark.parametrize(
    "bloc",
    [
        dict(type="pid", Kp=1, Ki=1, Kd=1, Tf=0.0),
        dict(type="pid", Kp=1, Ki=1, Kd=1, Tf=-0.1),
        dict(type="avance_retard", K=1, T1=0.1, T2=0.0),
        dict(type="avance_retard", K=1, T1=0.1, T2=-1.0),
        # impropre : degre du numerateur (2) > degre du denominateur (1)
        dict(type="fonction_transfert", num=[1.0, 0.0, 0.0], den=[1.0, 1.0]),
    ],
    ids=["pid_Tf_nul", "pid_Tf_neg", "ar_T2_nul", "ar_T2_neg", "ft_impropre"],
)
def test_parametres_interdits_sont_refuses(bloc):
    """Temoin : le meme bloc avec des valeurs admises compile ; refus : la
    valeur interdite seule fait lever ``SchemaInvalide``."""
    ok = dict(bloc)
    if ok["type"] == "pid":
        ok["Tf"] = 0.1
    elif ok["type"] == "avance_retard":
        ok["T2"] = 0.1
    else:
        ok["num"] = [1.0, 0.0]
    compiler(banc_bloc(ok))
    with pytest.raises(SchemaInvalide):
        compiler(banc_bloc(bloc))


def test_saturation_libre_declare_son_drapeau_de_mode():
    """Comme ``limiteur_courant`` : un drapeau ``mode_<id>`` (parametre nul :
    fonctionnement libre) que la simulation temporelle bascule."""
    c = compiler(banc_bloc(dict(type="saturation", min=-1.0, max=1.0)))
    assert c.parametres.get("mode_b") == 0.0, sorted(c.parametres)


# ---- saturation : elle borne la sortie en simulation temporelle
def _reseau_source(schema: dict) -> Tables:
    """Bus infini, ligne, bus 2 ou l'element a blocs injecte 20 MW / 5 Mvar
    (base systeme 100 MVA : ``id`` de 0,2 pu, ``kS = 1``)."""
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


def _source_saturee(lo: float, hi: float) -> dict:
    """``id`` = filtre(saturation(P_ref)), ``iq`` = Q_ref."""
    return {
        "nom": "source_saturee",
        "parametres": {"tau": 0.005},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref"},
            {"id": "qr", "type": "consigne", "nom": "Q_ref"},
            {"id": "sat", "type": "saturation", "min": lo, "max": hi},
            {"id": "fp", "type": "premier_ordre", "K": 1, "tau": "tau"},
            {"id": "oi", "type": "sortie", "signal": "id"},
            {"id": "oq", "type": "sortie", "signal": "iq"},
        ],
        "liens": [
            {"de": "pr", "vers": "sat"},
            {"de": "sat", "vers": "fp"},
            {"de": "fp", "vers": "oi"},
            {"de": "qr", "vers": "oq"},
        ],
    }


def _simuler_echelon(amplitude: float, lo=-0.3, hi=0.3):
    from g2egrid.temporel import Evenement, ModeleNonLineaire

    net = assembler(_reseau_source(_source_saturee(lo, hi)))
    assert net.run_loadflow().converged
    net.solve_equilibrium()
    el = next(e for e in net.elements if e.name.startswith("Blocs"))
    m = ModeleNonLineaire(net)
    u = m.scenario(Evenement(0.05, f"{el.name}.P_ref", amplitude))
    tr = m.simuler((0.0, 0.3), uext_fn=u, methode="BDF2", pas=1e-3, limiteurs=True)
    return tr.t, tr.etat(f"{el.name}.x_fp")


@pytest.mark.parametrize(
    "amplitude, attendu",
    [
        (+0.5, 0.3),  # P_ref = 0,70 : sature a max
        (-0.9, -0.3),  # P_ref = -0,70 : sature a min
        (+0.05, None),  # P_ref = 0,25 : reste libre, y = u
    ],
)
def test_saturation_borne_la_sortie_en_simulation_temporelle(amplitude, attendu):
    """Le point de fonctionnement (id = 0,1995 pu, constate sur le banc avec un
    gain a la place de la saturation) est libre ; apres l'echelon de consigne,
    le filtre aval (tau = 5 ms, regime etabli atteint bien avant 0,3 s) suit la
    sortie de la saturation : elle vaut ``max`` ou ``min`` si la consigne les
    depasse, la consigne sinon (``attendu`` None : x0 + amplitude). Le filtre
    ne depasse jamais la borne (1e-9) : la garde a bascule, la saturation ne
    s'est pas contentee d'etre lineaire."""
    _t, x = _simuler_echelon(amplitude)
    assert x[0] == pytest.approx(0.2, abs=5e-3)
    if attendu is None:
        attendu = x[0] + amplitude
    assert x[-1] == pytest.approx(attendu, abs=1e-6)
    assert x.max() <= max(0.3, attendu) + 1e-9
    assert x.min() >= min(-0.3, attendu) - 1e-9


# ---- implicite
def test_implicite_donne_la_solution_d_un_calcul_explicite():
    """``y`` tel que ``u1.y - u2.(1 + y) = 0`` vaut ``u2/(u1 - u2)`` : le meme
    nombre que la chaine explicite somme(+-) puis division (a la main :
    u1 = 3, u2 = 4 donne y = -4)."""
    s = {
        "nom": "implicite_vs_explicite",
        "parametres": {},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref"},
            {"id": "qr", "type": "consigne", "nom": "Q_ref"},
            {"id": "a", "type": "constante", "valeur": 3.0},
            {"id": "b", "type": "constante", "valeur": 4.0},
            {"id": "im", "type": "implicite", "n": 2, "expr": "u1*y - u2*(1 + y)"},
            {"id": "d", "type": "somme", "signes": "+-"},
            {"id": "q", "type": "division"},
            {"id": "oi", "type": "sortie", "signal": "id"},
            {"id": "oq", "type": "sortie", "signal": "iq"},
        ],
        "liens": [
            {"de": "a", "vers": "im", "entree": 0},
            {"de": "b", "vers": "im", "entree": 1},
            {"de": "a", "vers": "d", "entree": 0},
            {"de": "b", "vers": "d", "entree": 1},
            {"de": "b", "vers": "q", "entree": 0},
            {"de": "d", "vers": "q", "entree": 1},
            {"de": "im", "vers": "oi"},
            {"de": "q", "vers": "oq"},
        ],
    }
    implicite, explicite = valeurs_statiques(s)
    assert explicite == pytest.approx(-4.0, abs=1e-12)
    assert implicite == pytest.approx(explicite, abs=1e-12)
    assert not compiler(s).etats  # contrainte algebrique : aucun etat


def test_implicite_non_lineaire_linearise_comme_la_derivee_implicite():
    """``y + 0,1.y^3 = u``  =>  dy/du = 1/(1 + 0,3.y0^2) au point (u0, y0),
    theoreme des fonctions implicites (formule a la main)."""
    s = banc_bloc(dict(type="implicite", n=1, expr="y + 0.1*y^3 - u1"))
    u0 = 0.7
    racines = np.roots([0.1, 0.0, 1.0, -u0])  # une seule racine reelle (derivee > 0)
    y0 = float(racines[np.isclose(racines.imag, 0.0)][0].real)
    point = {sp.Symbol("vd_t"): sp.Float(u0), sp.Symbol("y_b"): sp.Float(y0)}
    ref = 1.0 / (1.0 + 0.3 * y0**2)
    g = fdt(s, 1.0, point=point)
    assert abs(g.imag) < 1e-12
    assert abs(g.real - ref) < TOL * ref


# ---- consignes calees ou non
def _source_a_biais(valeur: float, calee: bool = False) -> dict:
    """``id`` = filtre(P_ref) + biais ; le biais est une consigne NON calee
    fixee a ``valeur`` (cahier : entree de service fixee a ``valeur``)."""
    return {
        "nom": "source_biais",
        "parametres": {"tau": 0.005},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref"},
            {"id": "qr", "type": "consigne", "nom": "Q_ref"},
            {
                "id": "bi",
                "type": "consigne",
                "nom": "biais",
                "calee": calee,
                "valeur": valeur,
            },
            {"id": "fp", "type": "premier_ordre", "K": 1, "tau": "tau"},
            {"id": "s", "type": "somme", "signes": "++"},
            {"id": "oi", "type": "sortie", "signal": "id"},
            {"id": "oq", "type": "sortie", "signal": "iq"},
        ],
        "liens": [
            {"de": "pr", "vers": "fp"},
            {"de": "fp", "vers": "s", "entree": 0},
            {"de": "bi", "vers": "s", "entree": 1},
            {"de": "s", "vers": "oi"},
            {"de": "qr", "vers": "oq"},
        ],
    }


def _x_equilibre(schema: dict, nom: str = "x_fp") -> float:
    net = assembler(_reseau_source(schema))
    assert net.run_loadflow().converged
    net.solve_equilibrium()
    el = next(e for e in net.elements if e.name.startswith("Blocs"))
    return float(el.x0[el.states.index(nom)])


def test_consigne_non_calee_est_fixee_a_sa_valeur():
    """Le courant injecte est fixe par le load flow : si le biais (non cale)
    passe de 0 a 0,25, la consigne calee P_ref baisse de 0,25 et l'etat du
    filtre (= P_ref a l'equilibre) aussi : exactement 0,25 d'ecart."""
    x0 = _x_equilibre(_source_a_biais(0.0))
    x1 = _x_equilibre(_source_a_biais(0.25))
    assert x0 - x1 == pytest.approx(0.25, abs=1e-9)


def test_trois_consignes_calees_sont_refusees_a_l_initialisation():
    """L'initialisation exige TOUJOURS exactement deux consignes calees : une
    troisieme (calee par defaut) est une erreur claire, jamais un Newton
    sous-determine. Le message peut venir de la compilation ou de
    l'initialisation, d'ou le bloc ``try`` unique."""
    classe_element(_source_a_biais(0.25))  # temoin : trois dont une non calee, admis
    s = _source_a_biais(0.25, calee=True)
    with pytest.raises((SchemaInvalide, ValueError)):
        classe_element(s)
        _x_equilibre(s)


def test_consigne_calee_par_defaut():
    """Sans ``calee`` ni ``valeur``, une consigne se comporte comme en
    version 1 : calee (un schema de la version 1 a deux consignes sans autre
    champ, et s'initialise)."""
    s = copie(SCHEMAS_V1["source_filtree"])
    assert all("calee" not in b for b in s["blocs"] if b["type"] == "consigne")
    assert len(compiler(s).consignes) == 2
    assert _x_equilibre(s, "x_fp") == pytest.approx(0.2, abs=1e-6)
