"""Imports et export (``g2egrid.interop``), contre les load flows des outils
d'origine : PYPOWER (MATPOWER) et pandapower, solutions gardees dans
``tests/donnees/interop/references.json`` (aucun des deux outils n'est
necessaire pour les tests)."""

import json
import pathlib

import numpy as np
import pytest

from g2egrid import interop
from g2egrid.catalogue import charger
from g2egrid.tables import assembler

D = pathlib.Path(__file__).resolve().parent / "donnees" / "interop"
REF = json.loads((D / "references.json").read_text())


def _lf(t, ideal=True):
    if ideal:  # sans resistance de causalite : la loi des noeuds exacte
        t.NOEUDS = "ideal"
    net = assembler(t)
    return net.run_loadflow(), net._bus_index()


def _ecarts(t, nom, decalage=0):
    """Ecart max des modules, et des angles a un decalage commun pres."""
    lf, ix = _lf(t)
    r = REF[nom]
    dv, da = [], []
    for b, vm, va in zip(r["bus"], r["vm"], r["va"]):
        if b + decalage not in ix:
            continue  # bus fusionne a un autre par un interrupteur ferme
        V = lf.V[ix[b + decalage]]
        dv.append(abs(abs(V) - vm))
        da.append((np.degrees(np.angle(V)) - va + 180) % 360 - 180)
    return max(dv), np.ptp(da)


@pytest.mark.parametrize("nom", ["case9", "case57"])
def test_matpower_load_flow_exact(nom):
    t = interop.depuis_matpower(D / f"{nom}.m", x_sur_r_max=None)
    dv, da = _ecarts(t, nom)
    assert dv < 1e-10 and da < 1e-8


def test_matpower_redonne_ieee57_du_catalogue():
    t = interop.depuis_matpower(D / "case57.m", nom="IEEE57", f0=60.0)
    c = charger("ieee57")
    for k in ("Y_DER", "T_Lines_pu", "T_Transfos", "T_SMs", "T_Shunts"):
        assert getattr(t, k) == getattr(c, k), k
    assert any("BASE_KV" in a for a in t.avertissements)


@pytest.mark.parametrize("rev", [33, 34])
def test_psse_comme_matpower(rev):
    """Le meme case9 ecrit en RAW (une branche en transformateur, plus un
    trois-enroulements a ignorer) donne le meme load flow."""
    a = interop.depuis_matpower(D / "case9.m", x_sur_r_max=None)
    b = interop.depuis_psse(D / f"case9_v{rev}.raw", x_sur_r_max=None)
    assert len(b.T_Transfos) == 1 and b.f0 == 60.0
    assert any("trois enroulements" in x for x in b.avertissements)
    np.testing.assert_allclose(_lf(b)[0].V, _lf(a)[0].V, atol=1e-12)


@pytest.mark.parametrize("nom", ["pp_simple", "pp_mv_oberrhein"])
def test_pandapower(nom):
    """Modules a 1e-5 (la magnetisation des transformateurs, ignoree, pese
    5e-6) ; les angles a 150 deg pres de chaque cote d'un couplage Yd."""
    t = interop.depuis_pandapower(D / f"{nom}.json", x_sur_r_max=None)
    dv, _ = _ecarts(t, nom, decalage=1)
    assert dv < 1e-5
    assert any("couplage 150" in a for a in t.avertissements)


def test_pandapower_interrupteurs_et_sgen():
    t = interop.depuis_pandapower(D / "pp_simple.json")
    # 7 bus pandapower, deux paires fusionnees par des interrupteurs fermes
    assert len(t.Y_DER) == 5
    assert len(t.T_GFL1s) == 1 and len(t.T_SMs) == 1 and len(t.T_IBs) == 1
    assert any(r[1] == 1 for r in t.Y_DER)  # un slack


def test_aller_retour_matpower(tmp_path):
    for t in (interop.depuis_matpower(D / "case57.m"), charger("ieee39")):
        texte = interop.vers_matpower(t, tmp_path / "cas.m")
        t2 = interop.depuis_matpower(tmp_path / "cas.m", f0=t.f0, x_sur_r_max=None)
        assert "mpc.branch" in texte
        V1, V2 = _lf(t)[0].V, _lf(t2)[0].V
        np.testing.assert_allclose(np.abs(V2), np.abs(V1), atol=1e-8)


def test_charges():
    z = interop.depuis_matpower(D / "case9.m")
    p = interop.depuis_matpower(D / "case9.m", charges="P")
    assert sum(r[7] for r in z.Y_DER) > 0 and not z.T_Loads
    assert sum(r[7] for r in p.Y_DER) == 0 and len(p.T_Loads) == 3
    with pytest.raises(ValueError):
        interop.depuis_matpower(D / "case9.m", charges="X")


def test_importer_par_extension(tmp_path):
    assert len(interop.importer(D / "case9.m").Y_DER) == 9
    assert len(interop.importer(D / "case9_v33.raw").Y_DER) == 9
    assert len(interop.importer(D / "pp_simple.json").Y_DER) == 5
    (tmp_path / "x.xyz").write_text("")
    with pytest.raises(ValueError):
        interop.importer(tmp_path / "x.xyz")


def test_atelier_importe_et_exporte():
    analyse = pytest.importorskip("g2egrid_api.analyse")
    for nom in ("case9.m", "case9_v33.raw", "pp_simple.json"):
        d = analyse.importer_externe(nom, (D / nom).read_text())
        assert d["tables"]["Y_DER"] and len(d["positions"]) == len(d["tables"]["Y_DER"])
    assert analyse.importer_externe(
        "pp_simple.json", (D / "pp_simple.json").read_text()
    )["avertissements"]
    with pytest.raises(ValueError):
        analyse.importer_externe("x.json", '{"pas": "pandapower"}')
    e = analyse.exporter_vers(charger("wscc9"), "matpower")
    t = interop.depuis_matpower(e["contenu"], x_sur_r_max=None)
    assert len(t.Y_DER) == 9 and e["avertissements"]


def test_conductance_shunt_signalee_et_en_puissance_constante():
    """``Gs`` (MW a 1 pu) devient une charge a puissance constante, pas une
    impedance constante comme dans MATPOWER : l'import le dit, et la charge tire
    bien ``Gs`` MW meme a V < 1 pu (MATPOWER : ``Gs.V^2``, moins que cela)."""
    mpc = {
        "baseMVA": 100.0,
        "bus": [
            [1, 3, 0, 0, 0, 0, 1, 1.0, 0, 20, 1, 1.1, 0.9],
            [2, 1, 0, 0, 5.0, 0, 1, 1.0, 0, 20, 1, 1.1, 0.9],
        ],
        "gen": [[1, 0, 0, 100, -100, 1.0, 100, 1, 100, 0]],
        "branch": [[1, 2, 0.05, 0.2, 0, 0, 0, 0, 0, 0, 1, -360, 360]],
    }
    t = interop.depuis_matpower(mpc, x_sur_r_max=None)
    assert any(
        "conductance shunt" in a and "puissance constante" in a
        for a in t.avertissements
    )
    lf, ix = _lf(t)
    assert lf.S[ix[2]].real * 100 == pytest.approx(-5.0, abs=1e-8)  # MW
    assert abs(lf.V[ix[2]]) < 0.998  # a V < 1 pu, Gs.V^2 serait plus petit
