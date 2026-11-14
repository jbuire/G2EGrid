"""Impedance virtuelle des GFM (``Rv_pu``, ``Lv_pu``, ``zv_dyn``) : equilibre
analytique exact pour toutes les variantes, egalite moteur / schema-bloc, et
le cas d'usage -- un GFM a statisme du 2e ordre, instable face a un reseau
fort, stabilise par une inductance virtuelle dynamique ; la forme quasi
statique, elle, destabilise la boucle de tension."""

import numpy as np
import pytest

from g2egrid.blocs import schema_gfm
from g2egrid.catalogue import charger
from g2egrid.elements.gfm import GFM_DEFAULT, _GFMBase
from g2egrid.tables import Tables, assembler

TOUTES = ("T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s")
ZV = dict(Lv_pu=0.1, Rv_pu=0.02)


def _gfm(net) -> list:
    return [e for e in net.all_elements if isinstance(e, _GFMBase)]


def _equilibre_analytique(t: Tables):
    """Le reseau a l'equilibre, celui des GFM SANS le Newton de raffinement :
    leur point est celui de ``initialize``."""
    net = assembler(t)
    for e in _gfm(net):
        e.raffiner_equilibre = False
    net.run_loadflow()
    net.solve_equilibrium()
    return net


def _residu(net) -> float:
    # 1e-8 comme tests/test_gfm_lois.py : la boucle de courant continue des
    # GFM3/GFM4 du SMIB (Tdc = 1 ms, Kpdc ~ 2,6e4) laisse ~1e-9 sur idc, avec
    # ou sans impedance virtuelle
    return max(
        float(np.abs(e.compiled().residual(e.x0, e.z0, e.u0, e.p)).max())
        for e in _gfm(net)
    )


def _pire(net) -> complex:
    lam = np.linalg.eigvals(net.build_state_space().A)
    lam = lam[np.abs(lam) > 1e-6]
    return complex(lam[np.argmax(lam.real)])


def _smib(table: str = "T_GFM2s", x_ib: float = 0.024, **reglages) -> Tables:
    """Le SMIB du catalogue, le GFM dans ``table``, bus infini de reactance
    ``x_ib`` (pu, R = X/10)."""
    t = charger("smib_gfm").copie()
    lignes = [r for k in TOUTES for r in getattr(t, k)]
    for k in TOUTES:
        getattr(t, k).clear()
    getattr(t, table).extend(lignes)
    ib = list(t.T_IBs[0])
    ib[5], ib[6] = x_ib / 10, x_ib
    t.T_IBs = [tuple(ib)]
    t.GFM_PARAMS = dict(t.GFM_PARAMS) | reglages
    return t


def test_par_defaut_pas_d_impedance_virtuelle():
    assert GFM_DEFAULT["Rv_pu"] == 0.0 and GFM_DEFAULT["Lv_pu"] == 0.0


@pytest.mark.parametrize("zv_dyn", [1.0, 0.0])
@pytest.mark.parametrize(
    "table, loi",
    [
        ("T_GFM1s", None),
        ("T_GFM2s", None),
        ("T_GFM3s", None),
        ("T_GFM4s", None),
        ("T_GFM1s", "dvoc"),
        ("T_GFM3s", "matching"),
    ],
)
def test_equilibre_analytique_exact(table, loi, zv_dyn):
    """Le repere local s'aligne sur la f.e.m. virtuelle : l'initialisation
    est un point fixe exact, sans Newton, sous les deux formes."""
    reg = dict(ZV, zv_dyn=zv_dyn) | ({"loi": loi} if loi else {})
    net = _equilibre_analytique(_smib(table, **reg))
    assert _residu(net) < 1e-8
    e = _gfm(net)[0]
    x = dict(zip(e.states, e.x0, strict=True))
    # la tension du condensateur n'est plus sur l'axe d : elle suit la f.e.m.
    # de la chute virtuelle
    assert abs(x["veq"]) > 1e-3


def test_equilibre_analytique_exact_port_continu():
    t = charger("hvdc2t").copie()
    t.GFM_PARAMS = dict(t.GFM_PARAMS) | ZV
    net = _equilibre_analytique(t)
    assert any(type(e).__name__ == "GFM3DC" for e in _gfm(net))
    assert _residu(net) < 1e-8


def test_meme_equilibre_dynamique_et_quasi_statique():
    """A l'equilibre, d(i_g)/dt = 0 : les deux formes donnent le meme point."""
    a = _gfm(_equilibre_analytique(_smib(**ZV, zv_dyn=1.0)))[0]
    b = _gfm(_equilibre_analytique(_smib(**ZV, zv_dyn=0.0)))[0]
    assert np.allclose(a.x0, b.x0, atol=1e-12)


def test_reseau_fort_stabilise_par_l_inductance_virtuelle():
    """GFM a statisme du 2e ordre (inertie virtuelle) face a un bus infini
    tres fort (x = 0,01 pu) : la boucle de puissance oscille et diverge (mode
    de quelques Hz). L'inductance virtuelle dynamique la ramene dans le demi-plan
    gauche ; la forme quasi statique, a la meme valeur, destabilise la boucle
    de tension (mode rapide, au-dela de 60 Hz)."""
    reg = dict(mp=0.02)
    sans = _pire(_equilibre_analytique(_smib(x_ib=0.01, **reg)))
    assert sans.real > 1.0 and abs(sans.imag) / (2 * np.pi) < 20
    dyn = _pire(_equilibre_analytique(_smib(x_ib=0.01, **reg, **ZV, zv_dyn=1.0)))
    assert dyn.real < -1.0
    qs = _pire(_equilibre_analytique(_smib(x_ib=0.01, **reg, **ZV, zv_dyn=0.0)))
    assert qs.real > 50.0 and abs(qs.imag) / (2 * np.pi) > 60


def test_sur_reseau_faible_l_inductance_virtuelle_ne_gene_pas():
    reg = dict(mp=0.02)
    sans = _pire(_equilibre_analytique(_smib(x_ib=0.5, **reg)))
    avec = _pire(_equilibre_analytique(_smib(x_ib=0.5, **reg, **ZV)))
    assert sans.real < 0 and avec.real < 0


# ------------------------------------------------- moteur = schema-bloc
def _reseau_gfm(element: str, reglages: dict) -> Tables:
    """Bus infini, ligne, GFM de 10 MVA a 8 MW au bus 2 : le GFM2 du moteur,
    ou le meme en schema-bloc (cf. tests/test_blocs.py)."""
    t = Tables(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (2, 2, 20e3, 1.01, 0.0, 8e6, 1e6, 0.0, 0.0),
        ],
        T_Lines_pu=[(1, 1, 2, 1.0, 0.01, 0.15, 0.0, 20e3, 100e6)],
        T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
    )
    noms = {"Rv_pu": "Rv", "Lv_pu": "Lv", "zv_dyn": "zv_dyn"}
    if element == "blocs":
        t.T_Blocs = [(1, 2, "gfm", 10e6)]
        t.MODELES_BLOCS = {
            "gfm": schema_gfm(**{noms[k]: v for k, v in reglages.items()})
        }
    else:
        t.T_GFM2s = [(1, 2, 10e6, 40e3)]
        t.GFM_PARAMS = dict(reglages)
    return t


@pytest.mark.parametrize("zv_dyn", [1.0, 0.0])
def test_gfm_en_blocs_egale_gfm2_avec_impedance_virtuelle(zv_dyn):
    spectres = []
    for element in ("moteur", "blocs"):
        net = assembler(_reseau_gfm(element, dict(ZV, zv_dyn=zv_dyn)))
        net.run_loadflow()
        net.solve_equilibrium()
        lam = np.linalg.eigvals(net.build_state_space().A)
        spectres.append(np.sort_complex(lam[np.abs(lam) > 1e-6]))
    a, b = spectres
    assert a.size == b.size
    assert np.max(np.abs(a - b) / np.maximum(1.0, np.abs(a))) < 1e-7
