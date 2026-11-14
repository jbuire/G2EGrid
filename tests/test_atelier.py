"""L'Atelier : les onglets de l'interface, en code, et au bit pres."""

from __future__ import annotations

import numpy as np
import pytest

import g2egrid.catalogue as cat
from g2egrid_api import analyse
from g2egrid_api.atelier import Atelier


def test_construire_en_code_puis_calculer():
    a = Atelier.vide("demo")
    b1 = a.bus(230e3)
    b2 = a.bus(230e3, Pl_W=80e6, Ql_var=20e6)
    a.ligne(b1, b2, R_pu=0.01, X_pu=0.1)
    a.equipement("IB", b1, slack=True)
    a.equipement("GFM1", b2, Sn_VA=100e6, P_W=50e6)
    lf = a.loadflow()
    assert lf.attrs["mismatch"] < 1e-8
    assert len(a.transits()) == 3  # la ligne et les deux transfos de raccordement
    T = a.T_eig()
    assert {"Re", "Im", "f_und [Hz]", "f_dam [Hz]", "zeta [%]", "Etat 1", "%1"} <= set(
        T.columns
    )


def test_meme_resultat_que_l_interface():
    """Un calcul en code et le meme calcul par l'API rendent LA MEME chose."""
    a = Atelier.catalogue("wscc9")
    t = cat.charger("wscc9")
    assert a.T_eig()["Re"].tolist() == [m["re"] for m in analyse.t_eig(t)["modes"]]
    assert np.array_equal(a.etat().A, analyse.construire(t).etat_stationnaire().A)


def test_un_slack_tout_pq_en_code():
    a = Atelier.catalogue("wscc9")
    v0 = a.loadflow()["V"].to_numpy()
    a.un_slack_tout_pq(slack_bus=1)
    assert (a.table("Y_DER")["Type"] == 3).sum() == 8
    assert np.abs(a.loadflow()["V"].to_numpy() - v0).max() < 1e-9


def test_un_reseau_a_zones_en_code():
    a = Atelier.catalogue("wscc3sm_zones")
    assert a.reseau is not None and len(a.tables.Y_DER) == 54
    im = a.impact(par_poste=False)
    assert im["ecart_tension"] < 1e-10
    with pytest.raises(ValueError, match="zones"):
        a.bus(20e3)


def test_enregistrer_en_code_ouvrir_dans_l_interface(tmp_path, monkeypatch):
    """Ecrit en code, le reseau apparait dans la liste de l'interface."""
    monkeypatch.setattr(cat, "DOSSIER_RESEAUX", tmp_path)
    a = Atelier.vide("demo")
    b1 = a.bus(20e3)
    a.equipement("IB", b1, slack=True)
    a.charge(b1, 1e6, 0.1e6)
    a.enregistrer("demo_code")
    assert "demo_code" in [f["nom"] for f in analyse.catalogue()]
    relu = Atelier.catalogue("demo_code")
    assert relu.table("Y_DER").shape[0] == a.table("Y_DER").shape[0]


def test_le_code_genere_s_execute(tmp_path, monkeypatch):
    """Le script de l'onglet Code tourne tel quel et redonne le meme T_eig."""
    d = analyse.code_python(analyse.depuis_catalogue("smib"))
    g: dict = {}
    exec(compile(d["python"], "code", "exec"), g)  # noqa: S102
    assert g["T_eig"]["Re"].tolist() == Atelier.catalogue("smib").T_eig()["Re"].tolist()
