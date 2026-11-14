"""FR-051 et FR-052, audit du 2026-11-12 : ce que la relecture d'ensemble a
trouve, verrouille.

* FORME : un schema-bloc ou un modele a equations mal forme (``blocs`` qui
  n'est pas une liste, bloc qui n'est pas un objet, ``parametres`` en
  liste...) est REFUSE avec la raison (``SchemaInvalide``, 422 par l'API) --
  ces cas donnaient une erreur interne 500 ;
* VALEURS : un parametre NaN ou infini, une constante ``1e999``, un diviseur
  nul (``tau = 0`` d'un premier ordre), une saturation ``min > max`` sont
  refuses a la compilation au lieu de donner un modele qui divise par zero ;
* un schema APLATI ne se regroupe plus en silence en un schema illisible ;
* REMPLACEMENT par le serveur : un element natif remplace par son
  schema-bloc garde son load flow et son spectre (1e-8), une charge de
  ``T_Loads`` comprise (sa puissance etait perdue par le remplacement fait
  dans le navigateur), et REVIENT a son etat natif exactement ;
* RETOUR AU NATIF robuste : apres des editions faites entre-temps (consigne
  de tension du bus de la charge, bus ajoute, identifiant repris), il rend
  les puissances de l'element sans ecraser ces editions ;
* un SCHEMA-BLOC dont on modifie une equation n'est pas ecrase : le dessin
  reste dans ``MODELES_BLOCS`` ;
* le PARAMETRE d'un modele se modifie par l'API pour un element converti, et
  est refuse (avec la marche a suivre) pour un element natif ;
* la BIBLIOTHEQUE refuse un modele mal forme ou trop gros, accepte un modele
  a equations, et rend 404 pour un modele absent.
"""

from __future__ import annotations

import copy
import json

import numpy as np
import pytest

import g2egrid.reseau.catalogue as cat
from g2egrid.elements import blocs
from g2egrid.elements.blocs import SchemaInvalide, compiler, exemples
from g2egrid.elements.blocs_equations import equations_depuis_schema
from g2egrid.reseau import equations_reseau as er
from g2egrid.reseau.catalogue import charger
from g2egrid.reseau.tables import assembler, tables_vers_json

TOL_SPECTRE = 1e-8


def _base() -> dict:
    return copy.deepcopy(exemples()["source_filtree"]["schema"])


def _spectre(t):
    net = assembler(t)
    lf = net.run_loadflow()
    assert lf.converged
    net.solve_equilibrium()
    lam = np.linalg.eigvals(net.build_state_space().A)
    return lf, lam[np.abs(lam) > 1e-6]


def _ecart(a, b) -> float:
    """Ecart relatif entre deux spectres, apparies au plus proche."""
    from scipy.optimize import linear_sum_assignment

    assert len(a) == len(b), (len(a), len(b))
    d = np.abs(a[:, None] - b[None, :])
    i, j = linear_sum_assignment(d)
    return float(np.max(d[i, j] / np.maximum(np.abs(a[i]), 1.0)))


# ------------------------------------------------------------------ forme
def _mal_formes():
    def v(f):
        s = _base()
        f(s)
        return s

    return {
        "blocs_texte": v(lambda s: s.update(blocs="abc")),
        "bloc_nombre": v(lambda s: s.update(blocs=[3])),
        "liens_texte": v(lambda s: s.update(liens="x")),
        "lien_nombre": v(lambda s: s.update(liens=[3])),
        "lien_entree_texte": v(lambda s: s["liens"][0].update(entree="a")),
        "parametres_texte": v(lambda s: s.update(parametres="x")),
        "parametres_liste": v(lambda s: s.update(parametres=[1])),
        "parametre_texte": v(lambda s: s["parametres"].update(tau="abc")),
        "parametre_nan": v(lambda s: s["parametres"].update(tau=float("nan"))),
        "parametre_infini": v(lambda s: s["parametres"].update(tau=float("inf"))),
        "parametre_booleen": v(lambda s: s["parametres"].update(tau=True)),
        "sous_systeme_schema_nombre": v(
            lambda s: s["blocs"].append(
                {"id": "ss", "type": "sous_systeme", "schema": 3}
            )
        ),
    }


@pytest.mark.parametrize("cas", list(_mal_formes()))
def test_schema_mal_forme_refuse_avec_sa_raison(cas):
    compiler(_base())  # temoin : le schema sain compile
    with pytest.raises(SchemaInvalide):
        compiler(_mal_formes()[cas])


@pytest.mark.parametrize(
    "cas", ["blocs_texte", "bloc_nombre", "liens_texte", "lien_nombre"]
)
def test_structure_refuse_un_schema_mal_forme(cas):
    s = _mal_formes()[cas]
    with pytest.raises(SchemaInvalide):
        blocs.aplatir(s)
    with pytest.raises(SchemaInvalide):
        blocs.regrouper(s, ["fp"], "g")
    with pytest.raises(SchemaInvalide):
        blocs.degrouper(s, "fp")


def _modele_eq() -> dict:
    return equations_depuis_schema(_base())


@pytest.mark.parametrize(
    "modif",
    [
        lambda m: m.update(etats="abc"),
        lambda m: m.update(etats=[3]),
        lambda m: m.update(algebriques={"a": 1}),
        lambda m: m.update(consignes=[3]),
        lambda m: m.update(sorties=["id"]),
        lambda m: m.update(parametres=[1]),
        lambda m: m["parametres"].update(tau=float("nan")),
        lambda m: m.update(angles="theta"),
    ],
)
def test_modele_a_equations_mal_forme_refuse(modif):
    compiler(_modele_eq())  # temoin
    m = _modele_eq()
    modif(m)
    with pytest.raises(SchemaInvalide):
        compiler(m)


# ---------------------------------------------------------------- valeurs
def test_constante_hors_de_portee_refusee():
    """``1e999`` se lit comme l'infini : SymPy en faisait un flottant fini,
    NumPy un ``inf`` a l'execution."""
    s = _base()
    s["blocs"][2]["K"] = "1e999"
    with pytest.raises(SchemaInvalide, match="hors de portee"):
        compiler(s)
    m = _modele_eq()
    m["etats"][0]["derivee"] = "1e999"
    with pytest.raises(SchemaInvalide, match="hors de portee"):
        compiler(m)


def test_diviseur_nul_refuse_et_negatif_admis():
    s = _base()
    s["parametres"]["tau"] = 0.0
    with pytest.raises(SchemaInvalide, match="tau = 0"):
        compiler(s)
    s["parametres"]["tau"] = -0.05  # instable mais bien pose : admis
    compiler(s)


def _avec_saturation(lo, hi) -> dict:
    s = _base()
    s["blocs"].append({"id": "sat", "type": "saturation", "min": lo, "max": hi})
    for lien in s["liens"]:
        if lien["de"] == "fp" and lien["vers"] == "op":
            lien["vers"] = "sat"
    s["liens"].append({"de": "sat", "vers": "op", "entree": 0})
    return s


def test_saturation_min_superieur_a_max_refusee():
    compiler(_avec_saturation(-1, 1))  # temoin
    with pytest.raises(SchemaInvalide, match="depasse max"):
        compiler(_avec_saturation(1, -1))


def test_regrouper_un_schema_aplati_dit_pourquoi():
    s = _base()
    s = blocs.regrouper(s, ["fp"], "filtre")
    plat = blocs.aplatir(s)
    ids = [b["id"] for b in plat["blocs"] if "__" in b["id"]]
    assert ids, "l'aplatissement renomme en __"
    with pytest.raises(SchemaInvalide, match="aplati"):
        blocs.regrouper(plat, ids[:1], "encore")


# ------------------------------------------------- remplacement par le serveur
def test_remplacer_par_blocs_puis_revenir_gfl():
    t = charger("smib_gfl")
    lf0, lam0 = _spectre(t)
    avant = tables_vers_json(t)
    n, cle = er.remplacer_par_modele(t, ("T_GFL2s", 0), "blocs")
    assert tables_vers_json(t) == avant, "les tables recues ne bougent pas"
    assert not n.T_GFL2s and [r[2] for r in n.T_Blocs] == [cle]
    assert n.MODELES_BLOCS[cle]["origine"]["element"].startswith("GFL")
    lf1, lam1 = _spectre(n)
    assert np.max(np.abs(lf0.V - lf1.V)) < 1e-10
    assert _ecart(lam0, lam1) < TOL_SPECTRE
    retour = er.revenir_au_natif(n, ("T_Blocs", 0))
    assert tables_vers_json(retour) == avant


def test_remplacer_deja_modele_ou_non_remplacable_refuse():
    t = charger("smib_gfl")
    n, _cle = er.remplacer_par_modele(t, ("T_GFL2s", 0), "blocs")
    with pytest.raises(SchemaInvalide, match="deja"):
        er.remplacer_par_modele(n, ("T_Blocs", 0), "blocs")
    with pytest.raises(SchemaInvalide, match="ne se remplace pas"):
        er.remplacer_par_modele(t, ("T_Lines_SI", 0), "blocs")
    with pytest.raises(SchemaInvalide, match="genre"):
        er.remplacer_par_modele(t, ("T_GFL2s", 0), "autre")


def _reseau_a_charge_speciale():
    """Un reseau du catalogue qui porte une charge de ``T_Loads`` seule a son
    bus, et l'indice de sa ligne."""
    for nom in ("cigre_b4", "nordic44", "ieee118"):
        t = charger(nom)
        for i in range(len(t.T_Loads)):
            try:
                er.remplacer_par_modele(t, ("T_Loads", i), "equations")
            except SchemaInvalide:
                continue
            return nom, t, i
    pytest.skip("aucune charge de T_Loads remplacable dans le catalogue")


def test_charge_de_t_loads_remplacee_garde_sa_puissance():
    """Le remplacement fait dans le navigateur retirait la ligne de
    ``T_Loads`` sans reporter sa puissance au bus : la charge disparaissait du
    load flow. Le serveur la declare en generation negative."""
    _nom, t, i = _reseau_a_charge_speciale()
    lf0 = assembler(t).run_loadflow()
    avant = tables_vers_json(t)
    n, _cle = er.remplacer_par_modele(t, ("T_Loads", i), "equations")
    lf1 = assembler(n).run_loadflow()
    assert lf1.converged
    assert np.max(np.abs(lf0.V - lf1.V)) < 1e-10
    k = next(
        j
        for j, r in enumerate(n.T_Blocs)
        if r[2].startswith(("PLoad", "ILoad", "ZLoad", "RLLoad"))
    )
    assert tables_vers_json(er.revenir_au_natif(n, ("T_Blocs", k))) == avant


# ----------------------------------------------------- retour au natif robuste
def test_revenir_au_natif_garde_les_editions_faites_entre_temps():
    """Charge de ``Y_DER`` (wscc9, bus 5) convertie ; l'utilisateur change
    ensuite la consigne de tension de CE bus et ajoute un bus en tete de
    table. Le retour rend Pl, Ql et Pg, Qg de la charge, a la bonne ligne,
    sans ecraser la consigne."""
    t = charger("wscc9")
    i5 = next(i for i, r in enumerate(t.Y_DER) if int(r[0]) == 5)
    origine = list(t.Y_DER[i5])
    assert origine[7] or origine[8], "le bus 5 porte une charge"
    n = er.modifier_equation(
        t, "Load5", "i_d", "0*i_d + " + _derivee(t, "Load5", "i_d")
    )
    if not n.T_Blocs:  # equation inchangee -> pas de conversion : la forcer
        n, _ = er.remplacer_par_modele(t, "Load5", "equations")
    k = next(i for i, r in enumerate(n.Y_DER) if int(r[0]) == 5)
    ligne = list(n.Y_DER[k])
    assert ligne[7] == 0 and ligne[8] == 0
    ligne[3] = 1.02  # l'utilisateur retouche la consigne de tension du bus
    n.Y_DER[k] = tuple(ligne)
    n.Y_DER.insert(
        0, (99, 3, origine[2], 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    )  # et ajoute un bus
    retour = er.revenir_au_natif(n, ("T_Blocs", 0))
    assert not retour.T_Blocs
    r5 = next(r for r in retour.Y_DER if int(r[0]) == 5)
    assert r5[3] == 1.02, "la retouche de l'utilisateur reste"
    for c in (5, 6, 7, 8):
        assert r5[c] == pytest.approx(origine[c], abs=1e-6), c
    assert int(retour.Y_DER[0][0]) == 99, "le bus ajoute n'est pas ecrase"


def _derivee(t, element, variable) -> str:
    m = er.modele_element(t, element)
    return next(e["derivee"] for e in m["etats"] if e["nom"] == variable)


def test_revenir_au_natif_identifiant_repris_entre_temps():
    t = charger("smib_gfl")
    ligne = t.T_GFL2s[0]
    n, _cle = er.remplacer_par_modele(t, ("T_GFL2s", 0), "equations")
    n.T_GFL2s.append((ligne[0], 99, *ligne[2:]))  # un autre prend le meme ID
    retour = er.revenir_au_natif(n, ("T_Blocs", 0))
    ids = [int(r[0]) for r in retour.T_GFL2s]
    assert len(ids) == 2 and len(set(ids)) == 2, ids


def test_revenir_au_natif_bus_supprime_dit_pourquoi():
    t = charger("wscc9")
    n, _ = er.remplacer_par_modele(t, "Load5", "equations")
    n.Y_DER[:] = [r for r in n.Y_DER if int(r[0]) != 5]
    with pytest.raises(SchemaInvalide, match="n'est plus dans"):
        er.revenir_au_natif(n, ("T_Blocs", 0))


# ------------------------------------------------ un dessin n'est pas ecrase
def test_modifier_une_equation_d_un_schema_bloc_garde_le_dessin():
    t = charger("smib_gfl")
    n, cle = er.remplacer_par_modele(t, ("T_GFL2s", 0), "blocs")
    dessin = copy.deepcopy(n.MODELES_BLOCS[cle])
    assert "blocs" in dessin
    m = er.modele_element(n, ("T_Blocs", 0))
    e = m["etats"][0]
    n2 = er.modifier_equation(
        n, ("T_Blocs", 0), e["nom"], f"2*({e['derivee']})", "etat"
    )
    nouveau = n2.T_Blocs[0][2]
    assert nouveau != cle and n2.MODELES_BLOCS[nouveau]["genre"] == "equations"
    sans_origine = {k: v for k, v in dessin.items() if k != "origine"}
    reste = {k: v for k, v in n2.MODELES_BLOCS[cle].items() if k != "origine"}
    assert reste == sans_origine, "le schema-bloc d'origine est intact"
    # et l'origine a suivi l'element : le retour au natif marche encore
    retour = er.revenir_au_natif(n2, ("T_Blocs", 0))
    assert retour.T_GFL2s == t.T_GFL2s and not retour.T_Blocs


# -------------------------------------------------------------------- API
@pytest.fixture
def client(tmp_path, monkeypatch):
    pytest.importorskip("fastapi")
    from fastapi.testclient import TestClient

    from g2egrid_api.app import app

    monkeypatch.setattr(cat, "DOSSIER_RESEAUX", tmp_path)
    c = TestClient(app, headers={"Host": "127.0.0.1"})
    if c.get("/api/blocs/types").status_code == 403:
        pytest.skip("atelier protege par jeton dans cet environnement")
    return c


@pytest.fixture(scope="module")
def smib_gfl():
    return tables_vers_json(charger("smib_gfl"))


def _post_brut(client, route, corps):
    """Envoie le JSON tel quel, NaN compris (``json.loads`` du serveur l'accepte)."""
    return client.post(
        route,
        content=json.dumps(corps, allow_nan=True),
        headers={"Content-Type": "application/json"},
    )


@pytest.mark.parametrize("cas", list(_mal_formes()))
def test_api_compiler_mal_forme_422_jamais_500(client, cas):
    r = _post_brut(client, "/api/blocs/compiler", {"schema_bloc": _mal_formes()[cas]})
    assert r.status_code == 422, (r.status_code, r.text[:200])
    assert "erreur interne" not in r.text


@pytest.mark.parametrize("route", ["/api/blocs/regrouper", "/api/blocs/aplatir"])
@pytest.mark.parametrize("cas", ["blocs_texte", "bloc_nombre", "lien_nombre"])
def test_api_structure_mal_formee_422(client, route, cas):
    corps = {"schema_bloc": _mal_formes()[cas], "ids": ["fp"], "nom": "g"}
    r = _post_brut(client, route, corps)
    assert r.status_code == 422, (r.status_code, r.text[:200])


def test_api_verifier_modele_mal_forme_sans_500(client):
    m = _modele_eq()
    m["etats"] = "abc"
    r = client.post("/api/equations/verifier", json={"modele": m})
    assert r.status_code == 200 and r.json()["ok"] is False
    assert "liste" in r.json()["message"]


def test_api_remplacer_puis_revenir(client, smib_gfl):
    r = client.post(
        "/api/blocs/remplacer",
        json={"tables": smib_gfl, "table": "T_GFL2s", "index": 0},
    )
    assert r.status_code == 200, r.text[:300]
    d = r.json()
    assert d["modele"] in d["tables"]["MODELES_BLOCS"]
    assert d["schema"]["blocs"] and not d["tables"]["T_GFL2s"]
    r2 = client.post(
        "/api/equations/revenir",
        json={"tables": d["tables"], "table": "T_Blocs", "index": 0},
    )
    assert r2.status_code == 200, r2.text[:300]
    assert r2.json()["tables"] == smib_gfl
    r3 = client.post(
        "/api/blocs/remplacer",
        json={"tables": smib_gfl, "table": "T_Lines_SI", "index": 0},
    )
    assert r3.status_code == 422 and "ne se remplace pas" in r3.json()["detail"]


def test_api_parametre_d_un_modele(client, smib_gfl):
    r = client.post(
        "/api/equations/parametre",
        json={
            "tables": smib_gfl,
            "table": "T_GFL2s",
            "index": 0,
            "nom": "Kipll",
            "valeur": 1.0,
        },
    )
    assert r.status_code == 422 and "PARAMS_PAR_BUS" in r.json()["detail"]
    d = client.post(
        "/api/blocs/remplacer",
        json={"tables": smib_gfl, "table": "T_GFL2s", "index": 0},
    ).json()
    cle = d["modele"]
    p0 = d["tables"]["MODELES_BLOCS"][cle]["parametres"]["Kipll"]
    corps = {
        "tables": d["tables"],
        "table": "T_Blocs",
        "index": 0,
        "nom": "Kipll",
        "valeur": 2 * p0,
    }
    r = client.post("/api/equations/parametre", json=corps)
    assert r.status_code == 200, r.text[:300]
    t2 = r.json()["tables"]
    cle2 = t2["T_Blocs"][0][2]
    assert t2["MODELES_BLOCS"][cle2]["parametres"]["Kipll"] == pytest.approx(2 * p0)
    r = client.post("/api/equations/parametre", json={**corps, "nom": "nexiste_pas"})
    assert r.status_code == 422
    r = _post_brut(
        client, "/api/equations/parametre", {**corps, "valeur": float("nan")}
    )
    assert r.status_code == 422


def test_bibliotheque_forme_taille_equations_et_404(client):
    mal = _mal_formes()["bloc_nombre"]
    r = _post_brut(client, "/api/blocs/bibliotheque", {"schema_bloc": mal})
    assert r.status_code == 422
    gros = _base()
    gros["blocs"] += [
        {"id": f"b{k}", "type": "constante", "valeur": 1, "note": "z" * 2000}
        for k in range(1500)
    ]
    r = client.post("/api/blocs/bibliotheque", json={"schema_bloc": gros})
    assert r.status_code == 422 and "trop gros" in r.json()["detail"]
    m = _modele_eq()
    m["nom"] = "modele_eq"
    r = client.post("/api/blocs/bibliotheque", json={"schema_bloc": m})
    assert r.status_code == 200, r.text[:200]
    assert client.get("/api/blocs/bibliotheque/modele_eq").json()["schema"] == m
    assert client.get("/api/blocs/bibliotheque/absent").status_code == 404
    assert client.delete("/api/blocs/bibliotheque/absent").status_code == 404
    assert client.delete("/api/blocs/bibliotheque/modele_eq").status_code == 200


# ------------------------------------------- le lieu de l'element, rendu
def test_retour_au_natif_rend_le_lieu_de_l_element():
    """L'interface suit l'element sans le deviner : la ligne reinseree pour
    un convertisseur, celle de ``Y_DER`` pour la charge d'un bus."""
    t = charger("smib_gfl")
    n, _cle = er.remplacer_par_modele(t, ("T_GFL2s", 0), "blocs")
    retour, lieu = er.retour_au_natif(n, ("T_Blocs", 0))
    assert lieu == ("T_GFL2s", 0)
    assert tables_vers_json(retour) == tables_vers_json(t)
    w = charger("wscc9")
    n, _cle = er.remplacer_par_modele(w, "Load5", "equations")
    _retour, (table, k) = er.retour_au_natif(n, ("T_Blocs", 0))
    assert table == "Y_DER" and int(w.Y_DER[k][0]) == 5


def test_api_reponses_donnent_le_lieu(client, smib_gfl):
    corps = {"tables": smib_gfl, "table": "T_GFL2s", "index": 0}
    d = client.post("/api/blocs/remplacer", json=corps).json()
    assert (d["table"], d["index"]) == ("T_Blocs", 0)
    e = client.post(
        "/api/equations/element",
        json={"tables": d["tables"], "table": "T_Blocs", "index": 0},
    ).json()
    assert e["natif"] == "GFL2" and e["schema_bloc"] is True
    natif = client.post("/api/equations/element", json=corps).json()
    assert "natif" not in natif and "schema_bloc" not in natif
    r = client.post(
        "/api/equations/revenir",
        json={"tables": d["tables"], "table": "T_Blocs", "index": 0},
    ).json()
    assert (r["table"], r["index"]) == ("T_GFL2s", 0)
    # une equation modifiee sur le natif : l'element est passe dans T_Blocs
    m = natif["modele"]["etats"][0]
    r = client.post(
        "/api/equations/modifier",
        json={**corps, "variable": m["nom"], "expression": f"2*({m['derivee']})"},
    ).json()
    assert (r["table"], r["index"]) == ("T_Blocs", 0)
    e = client.post(
        "/api/equations/element",
        json={"tables": r["tables"], "table": "T_Blocs", "index": 0},
    ).json()
    assert e["natif"] == "GFL2" and e["schema_bloc"] is False


def test_api_element_modele_ecrit_a_la_main_sans_natif(client):
    """Un modele pose a la main dans ``T_Blocs`` ne remplace aucun element
    natif : ``natif`` est ``None`` (l'interface ne propose pas d'y revenir),
    et le retour est refuse avec la raison."""
    t = charger("smib_gfl")
    n, cle = er.remplacer_par_modele(t, ("T_GFL2s", 0), "blocs")
    n.MODELES_BLOCS[cle].pop("origine")
    tables = tables_vers_json(n)
    corps = {"tables": tables, "table": "T_Blocs", "index": 0}
    e = client.post("/api/equations/element", json=corps).json()
    assert e["natif"] is None
    r = client.post("/api/equations/revenir", json=corps)
    assert r.status_code == 422 and "origine native" in r.json()["detail"]


# --------------------------- le parametre d'un schema-bloc garde le dessin
def test_parametre_d_un_schema_bloc_garde_le_dessin():
    """Changer une VALEUR ne doit pas faire d'un schema-bloc un modele a
    equations : le dessin reste, sous le meme nom, avec son origine ; le
    spectre suit la valeur, et le retour au natif marche encore."""
    t = charger("smib_gfl")
    n, cle = er.remplacer_par_modele(t, ("T_GFL2s", 0), "blocs")
    _lf, lam0 = _spectre(n)
    k0 = n.MODELES_BLOCS[cle]["parametres"]["Kipll"]
    n2 = er.modifier_parametre(n, ("T_Blocs", 0), "Kipll", 3 * k0)
    assert [r[2] for r in n2.T_Blocs] == [cle] and list(n2.MODELES_BLOCS) == [cle]
    m = n2.MODELES_BLOCS[cle]
    assert "blocs" in m and m.get("genre") != "equations"
    assert m["parametres"]["Kipll"] == pytest.approx(3 * k0)
    assert m["origine"] == n.MODELES_BLOCS[cle]["origine"]
    assert n.MODELES_BLOCS[cle]["parametres"]["Kipll"] == k0, "les tables recues"
    _lf, lam2 = _spectre(n2)
    assert _ecart(lam0, lam2) > 1e-3, "tripler le gain de la PLL deplace le spectre"
    with pytest.raises(SchemaInvalide, match="pas de parametre"):
        er.modifier_parametre(n, ("T_Blocs", 0), "nexiste_pas", 1.0)
    retour = er.revenir_au_natif(n2, ("T_Blocs", 0))
    assert tables_vers_json(retour) == tables_vers_json(t)


def test_parametre_refuse_ne_change_rien():
    """Une valeur que la compilation refuse (diviseur nul) laisse les tables
    intactes."""
    t = charger("smib_gfl")
    n, cle = er.remplacer_par_modele(t, ("T_GFL2s", 0), "blocs")
    avant = tables_vers_json(n)
    diviseurs = [p for p in ("L1_pu", "C_pu", "L2_pu", "Cdc_pu", "Tdc")
                 if p in n.MODELES_BLOCS[cle]["parametres"]]  # fmt: skip
    refus = 0
    for p in diviseurs:
        try:
            er.modifier_parametre(n, ("T_Blocs", 0), p, 0.0)
        except SchemaInvalide:
            refus += 1
    assert refus, f"aucun diviseur nul refuse parmi {diviseurs}"
    assert tables_vers_json(n) == avant


def test_loadflow_converge_dit_ce_que_dit_le_noyau():
    """``converge`` de l'API etait ``iterations > 0`` : c'est maintenant le
    drapeau du load flow lui-meme."""
    from g2egrid_api.analyse import loadflow

    d = loadflow(charger("wscc9"))
    assert d["converge"] is True and d["mismatch"] < 1e-8
