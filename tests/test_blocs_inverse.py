"""FR-051, compilateur inverse (equations -> blocs) : ``depuis_element`` ecrit
le schema-bloc d'une classe NATIVE ; pose a la place du natif dans un reseau
du catalogue (``T_Blocs``), il donne le meme load flow, s'initialise et a les
memes valeurs propres (a 1e-8 relatif, mode de jauge exclu) -- pour les deux
niveaux de detail. Edition d'un gain, disposition sans recouvrement, refus
des elements a port continu."""

from __future__ import annotations

import copy

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

from g2egrid.elements.blocs import (
    BLOC_H,
    BLOC_L,
    SchemaInvalide,
    aplatir,
    compiler,
    depuis_element,
)
from g2egrid.reseau.catalogue import charger
from g2egrid.reseau.tables import assembler


def _resoudre(t):
    net = assembler(t)
    assert net.run_loadflow().converged
    net.solve_equilibrium()
    lam = np.linalg.eigvals(net.build_state_space().A)
    return net, lam[np.abs(lam) > 1e-6]


def _ecart_spectres(a: np.ndarray, b: np.ndarray) -> float:
    """Ecart relatif maximal entre deux spectres APPARIES (affectation
    optimale, pas un simple tri)."""
    assert a.size == b.size, (a.size, b.size)
    cout = np.abs(a[:, None] - b[None, :])
    i, j = linear_sum_assignment(cout)
    return float(np.max(np.abs(a[i] - b[j]) / np.maximum(1.0, np.abs(a[i]))))


def _poser(t, table: str, k: int, schema: dict, load: bool = False):
    """``t`` ou la ligne ``k`` de ``table`` est remplacee par un element a
    blocs de ``schema`` (meme bus, meme ``Sn_VA``)."""
    t = copy.deepcopy(t)
    if load:
        r = list(t.Y_DER[k])
        bus, P, Q = int(r[0]), r[7], r[8]
        r[5], r[6], r[7], r[8] = r[5] - P, r[6] - Q, 0.0, 0.0
        t.Y_DER[k] = tuple(r)
        sn = 1e6
    else:
        rows = list(getattr(t, table))
        row = rows.pop(k)
        setattr(t, table, rows)
        bus, sn = int(row[1]), row[2]
    t.T_Blocs = [*t.T_Blocs, (900, bus, "modele_fr051", sn)]
    t.MODELES_BLOCS = {**t.MODELES_BLOCS, "modele_fr051": schema}
    return t


def _gfl_vers(t, table: str, extra: dict | None = None):
    """Un reseau ou le premier GFL1 passe dans ``table`` (autre classe)."""
    t = copy.deepcopy(t)
    row = t.T_GFL1s.pop(0)
    setattr(t, table, [*getattr(t, table), row])
    if extra:
        t.GFM_PARAMS = {**t.GFM_PARAMS, **extra}
    return t


def _cas(nom: str):
    """``(tables natives, table, indice, charge ?, classe, parametres)``."""
    if nom == "GFL1":
        return charger("wscc9_cigre6"), "T_GFL1s", 0, False, "GFL1"
    if nom == "GFL2":
        return (
            _gfl_vers(charger("wscc9_cigre6"), "T_GFL2s"),
            "T_GFL2s",
            0,
            False,
            "GFL2",
        )
    if nom == "GFM1":
        return (
            _gfl_vers(charger("wscc9_cigre6"), "T_GFM1s"),
            "T_GFM1s",
            0,
            False,
            "GFM1",
        )
    if nom == "GFM2":
        return (
            _gfl_vers(charger("wscc9_cigre6"), "T_GFM2s"),
            "T_GFM2s",
            0,
            False,
            "GFM2",
        )
    if nom == "GFMdVOC":
        t = _gfl_vers(charger("wscc9_cigre6"), "T_GFM1s", {"loi": "dvoc"})
        return t, "T_GFM1s", 0, False, "GFMdVOC"
    if nom == "GFMMatching":
        t = _gfl_vers(charger("wscc9_cigre6"), "T_GFM3s", {"loi": "matching"})
        return t, "T_GFM3s", 0, False, "GFMMatching"
    if nom == "GFM3":
        return (
            _gfl_vers(charger("wscc9_cigre6"), "T_GFM3s"),
            "T_GFM3s",
            0,
            False,
            "GFM3",
        )
    if nom == "GFM4":
        return (
            _gfl_vers(charger("wscc9_cigre6"), "T_GFM4s"),
            "T_GFM4s",
            0,
            False,
            "GFM4",
        )
    if nom == "IB":
        return charger("smib"), "T_IBs", 0, False, "IB"
    if nom == "SynchronousMachine":
        return charger("wscc9"), "T_SMs", 1, False, "SynchronousMachine"
    if nom == "RLLoad":
        t = charger("wscc9")
        k = next(i for i, r in enumerate(t.Y_DER) if r[7])
        return t, "Y_DER", k, True, "RLLoad"
    raise KeyError(nom)


CLASSES = [
    "GFL1",
    "GFL2",
    "GFM1",
    "GFM2",
    "GFMdVOC",
    "GFMMatching",
    "GFM3",
    "GFM4",
    "IB",
    "SynchronousMachine",
    "RLLoad",
]


@pytest.fixture(scope="module")
def natifs():
    """Load flow et spectre du reseau natif de chaque classe (calcules une
    fois)."""
    out = {}
    for nom in CLASSES:
        t, table, k, load, classe = _cas(nom)
        net, lam = _resoudre(t)
        out[nom] = (t, table, k, load, classe, net.lf.V.copy(), lam)
    return out


@pytest.mark.parametrize("detail", ["equations", "blocs"])
@pytest.mark.parametrize("nom", CLASSES)
def test_aller_retour_dae_blocs_dae(natifs, nom, detail):
    """Critere 4 : meme load flow, memes valeurs propres a 1e-8."""
    t, table, k, load, classe, V0, lam0 = natifs[nom]
    schema = depuis_element(classe, detail=detail, f0=t.f0)
    assert schema["amorce"]["classe"] == classe
    t_blocs = _poser(t, table, k, schema, load)
    net, lam = _resoudre(t_blocs)
    assert np.max(np.abs(net.lf.V - V0)) < 1e-8
    assert _ecart_spectres(lam0, lam) < 1e-8


@pytest.mark.parametrize("nom", CLASSES)
def test_disposition_sans_recouvrement(nom):
    """Critere 6 : a tous les niveaux, les boites BLOC_L x BLOC_H des blocs
    sont disjointes."""
    for detail in ("equations", "blocs"):
        schema = depuis_element(nom, detail=detail)

        def verifier(niveau, ou):
            bs = niveau["blocs"]
            assert all("x" in b and "y" in b for b in bs), ou
            for i, a in enumerate(bs):
                for b in bs[i + 1 :]:
                    assert (
                        abs(a["x"] - b["x"]) >= BLOC_L or abs(a["y"] - b["y"]) >= BLOC_H
                    ), (ou, a["id"], b["id"])
                if a["type"] == "sous_systeme":
                    verifier(a["schema"], f"{ou}>{a['id']}")

        verifier(schema, nom)


def test_structure_du_schema_genere():
    s = depuis_element("GFL1")
    types = [b["type"] for b in s["blocs"]]
    assert types.count("integrateur") == 12
    assert (
        types.count("sous_systeme") >= 15
    )  # f (sauf les derivees reduites a une variable), g, h
    assert {b["nom"] for b in s["blocs"] if b["type"] == "consigne"} == {
        "p_ref",
        "q_ref",
    }
    assert {b["signal"] for b in s["blocs"] if b["type"] == "sortie"} == {"id", "iq"}
    assert "kS" not in s["parametres"] and "wb" not in s["parametres"]
    assert s["parametres"]["mode_lim"] == 0.0
    ss = next(b for b in s["blocs"] if b["type"] == "sous_systeme")
    assert ss["nom"].startswith("d") or ss["nom"].startswith("m") or ss["nom"]
    noms_ports = [b["nom"] for b in ss["schema"]["blocs"] if b["type"] == "port_entree"]
    assert noms_ports  # les ports portent le nom de la variable
    eq = depuis_element("GFL1", detail="equations")
    assert not any(b["type"] == "sous_systeme" for b in eq["blocs"])
    assert sum(b["type"] == "fonction" for b in eq["blocs"]) >= 12


def test_formeurs_exposent_angle_et_vitesse():
    for nom in ("GFM1", "GFMdVOC", "IB", "SynchronousMachine"):
        s = depuis_element(nom)
        assert {"theta", "w"} <= {
            b["signal"] for b in s["blocs"] if b["type"] == "sortie"
        }
    s = depuis_element("GFL1")
    assert {b["signal"] for b in s["blocs"] if b["type"] == "sortie"} == {"id", "iq"}


def test_consignes_calees_et_fixees():
    s = depuis_element("GFM2")
    c = {b["nom"]: b for b in s["blocs"] if b["type"] == "consigne"}
    assert "calee" not in c["p_ref"] and "calee" not in c["ve_ref"]
    assert c["q_ref"]["calee"] is False and c["q_ref"]["valeur"] == 0.0
    cp = compiler(s)
    assert cp.calees == ["p_ref", "ve_ref"]


def test_edition_d_un_gain_egale_le_natif(natifs):
    """Critere 5 : Kipll x 2 dans le schema genere = Kipll x 2 sur le natif."""
    t, table, k, load, _classe, _V0, _lam0 = natifs["GFL1"]
    Kipll = depuis_element("GFL1", f0=t.f0)["parametres"]["Kipll"]
    # natif : le reglage par bus
    tn = copy.deepcopy(t)
    bus = int(getattr(tn, table)[k][1])
    tn.PARAMS_PAR_BUS = {**tn.PARAMS_PAR_BUS, bus: {"Kipll": 2 * Kipll}}
    _n, lam_nat = _resoudre(tn)
    # blocs : le meme gain edite dans le schema
    s = depuis_element("GFL1", f0=t.f0)
    s["parametres"]["Kipll"] = 2 * Kipll
    _n, lam_blocs = _resoudre(_poser(t, table, k, s, load))
    assert _ecart_spectres(lam_nat, lam_blocs) < 1e-8
    # et le gain a bien bouge le spectre
    _n, lam_base = _resoudre(t)
    assert _ecart_spectres(lam_base, lam_nat) > 1e-3


def test_parametres_donnes():
    s = depuis_element("GFL1", {"Kipll": 1234.5, "R1_pu": 0.02})
    assert s["parametres"]["Kipll"] == 1234.5 and s["parametres"]["R1_pu"] == 0.02


def test_amorce_recales():
    assert "Kpdc" in depuis_element("GFL2")["amorce"]["recales"]
    assert "theta_g0" in depuis_element("SynchronousMachine")["amorce"]["recales"]
    assert {"R", "X"} <= set(depuis_element("RLLoad")["amorce"]["recales"])
    assert "ve_ref" in depuis_element("GFMdVOC")["amorce"]["recales"]


def test_sans_amorce_l_initialisation_repart_de_zero(natifs):
    """L'amorce n'est qu'un point de depart : sans elle, le Newton d'avant
    atteint le meme point (GFL1)."""
    t, table, k, load, _classe, _V0, lam0 = natifs["GFL1"]
    s = depuis_element("GFL1", f0=t.f0)
    s.pop("amorce")
    _n, lam = _resoudre(_poser(t, table, k, s, load))
    assert _ecart_spectres(lam0, lam) < 1e-8


@pytest.mark.parametrize("nom", ["GFM3DC", "GFL2DC", "MMC", "Line", "Transformer"])
def test_elements_a_port_continu_ou_multiple_refuses(nom):
    with pytest.raises(SchemaInvalide, match="port"):
        depuis_element(nom)


def test_classe_inconnue_et_detail_invalide():
    with pytest.raises(SchemaInvalide, match="inconnue"):
        depuis_element("Nexiste")
    with pytest.raises(SchemaInvalide, match="detail"):
        depuis_element("GFL1", detail="xx")


def test_aplatir_du_schema_genere_compile():
    s = depuis_element("GFM1")
    assert compiler(aplatir(s)).etats == compiler(s).etats
