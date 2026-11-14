"""Les fiches d'equipement (:mod:`g2egrid.fiches`) : equations lisibles,
ordonnees par sous-systeme, et EXACTEMENT celles du modele compile."""

import pytest

from g2egrid.catalogue import charger
from g2egrid.equations import classes
from g2egrid.fiches import FICHES, SCHEMAS, fiche, symbole_affiche, verifier


@pytest.mark.parametrize("nom", sorted(FICHES))
def test_la_fiche_est_le_modele(nom):
    """Chaque signal intermediaire substitue par sa definition, chaque etat,
    variable algebrique et sortie redigee coincide avec le modele compile
    (fonctionnement libre) ; aucun etat n'est oublie."""
    ecarts = verifier(nom)
    assert ecarts
    assert all(e == 0 for e in ecarts.values()), {
        k: v for k, v in ecarts.items() if v != 0
    }


@pytest.mark.parametrize("nom", sorted(FICHES))
def test_chaque_etat_une_fois(nom):
    f = fiche(nom)
    etats = [
        li["nom"] for b in f["blocs"] for li in b["lignes"] if li["genre"] == "etat"
    ]
    assert sorted(etats) == sorted(f["etats"])
    assert f["schema"] == SCHEMAS[nom]
    assert all(
        li["mathml"].startswith("<math") for b in f["blocs"] for li in b["lignes"]
    )


def test_les_convertisseurs_ont_leur_fiche():
    """Tous les GFM, GFL et MMC exportes ont une fiche redigee."""
    convertisseurs = {n for n in classes() if n.startswith(("GFM", "GFL", "MMC"))}
    assert convertisseurs <= set(FICHES)


def test_ordre_des_sous_systemes():
    titres = [b["titre"] for b in fiche("GFL1")["blocs"]]
    assert titres.index("PLL") < titres.index("Boucles externes P et Q")
    assert titres.index("Boucle de courant") < titres.index("Onduleur et filtre LCL")
    titres = [b["titre"] for b in fiche("GFM4")["blocs"]]
    assert any("2e ordre" in t for t in titres)
    assert titres[-2] == "Bus continu"


def test_ordre_redige_des_termes():
    """L'affichage garde l'ordre ecrit : i* - i, pas -i + i*."""
    f = fiche("GFL1")
    cl = next(b for b in f["blocs"] if b["titre"] == "Boucle de courant")
    md = next(li for li in cl["lignes"] if li["nom"] == "md")
    assert md["latex"].startswith(r"m_{d} = K_{pCL} \left(i^{*}_{sd} - i_{sd}\right)")


def test_fiche_generique():
    f = fiche("SynchronousMachine")
    assert not f["ecrite"] and f["schema"] is None
    assert [b["titre"] for b in f["blocs"]][0] == "Equations d'etat"


@pytest.mark.parametrize(
    "nom, attendu",
    [
        ("isd", "i_sd"),
        ("KpCL", "K_pCL"),
        ("R1_pu", "R_1"),
        ("p_ref", "p^*"),
        ("vdc_ref", "v_dc^*"),
        ("wb", "ω_b"),
    ],
)
def test_symboles(nom, attendu):
    assert symbole_affiche(nom).name == attendu


def test_api_fiche():
    fastapi = pytest.importorskip("fastapi")  # noqa: F841
    from fastapi.testclient import TestClient

    from g2egrid.tables import tables_vers_json
    from g2egrid_api.app import app

    client = TestClient(app, headers={"Host": "127.0.0.1"})
    t = tables_vers_json(charger("cigre_b4_dcs1"))
    r = client.post("/api/fiche", json={"tables": t, "table": "T_GFL2s", "index": 1})
    if r.status_code == 403:
        pytest.skip("atelier protege par jeton dans cet environnement")
    assert r.status_code == 200, r.text
    d = r.json()
    assert d["element"]["nom"] == "GFL121" and d["classe"] == "GFL2DC"
    assert d["equilibre"]["ved"] == pytest.approx(1.0, abs=0.05)
    assert any(p["nom"] == "Kpdc" for p in d["parametres"])
