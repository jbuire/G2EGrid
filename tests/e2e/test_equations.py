"""FR-052, partie interface : les equations dynamiques modifiables, dans Chromium.

Deux familles de tests, comme pour FR-051 :

* ceux qui tournent contre un serveur SIMULE (``page.route`` sur
  ``/api/equations/*``, qui suit le contrat du cahier des charges) : ils
  eprouvent l'interface seule -- liste, filtres, edition sur place, section de
  la fiche et des proprietes, passage entre equations et blocs, annuler /
  retablir, messages d'erreur -- et passent contre n'importe quel serveur ;
* ceux du critere 8 contre le VRAI serveur (modifier dans l'onglet, relancer
  l'analyse modale, annuler ; fiche d'un GFM ; ajouter un etat ; revenir au
  natif) : ils se SAUTENT tant que le serveur n'a pas les routes.

Les fixtures (serveur, navigateur, page) sont celles de ``test_navigateur``.
"""

from __future__ import annotations

import copy
import json

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("playwright.sync_api")

import test_navigateur as _nav

serveur = _nav.serveur
navigateur = _nav.navigateur
page = _nav.page
_charger = _nav._charger


# --------------------------------------------------------- serveur simule
def _modele_gfl1() -> dict:
    return {
        "genre": "equations",
        "nom": "GFL1_mod",
        "parametres": {"Kipll": 1200.0, "Kppll": 60.0},
        "consignes": [
            {"nom": "p_ref", "calee": True},
            {"nom": "q_ref", "calee": True},
        ],
        "etats": [
            {"nom": "theta_pll", "derivee": "wb*(w_pll - w)"},
            {"nom": "M_pll", "derivee": "Kipll*vq"},
        ],
        "algebriques": [{"nom": "w_pll", "expr": "1 + Kppll*vq + M_pll"}],
        "sorties": {"id": "p_ref/V", "iq": "-q_ref/V"},
        "amorce": {"classe": "GFL1"},
    }


class FauxServeur:
    """Les routes /api/equations/* du contrat, a la main, pour un GFL1."""

    def __init__(self, page) -> None:
        self.page = page
        self.journal: list[tuple[str, dict]] = []
        self.compteur = 0
        self.refus: str | None = (
            None  # un message : la prochaine ecriture est refusee (422)
        )
        page.route("**/api/equations/*", self.gerer)

    # -- les lignes
    def _lignes(self, t: dict) -> list[dict]:
        out = []
        for i, r in enumerate(t.get("T_GFL1s", [])):
            out += self._lignes_modele(
                _modele_gfl1(), f"GFL1_{r[0]}", "T_GFL1s", i, int(r[1]), "GFL1", False
            )
        for i, r in enumerate(t.get("T_Blocs", [])):
            m = t["MODELES_BLOCS"][r[2]]
            out += self._lignes_modele(
                m, f"Blocs_{r[0]}", "T_Blocs", i, int(r[1]), m["amorce"]["classe"], True
            )
        for i, r in enumerate(t.get("T_Lines_pu", [])):
            out.append(
                {
                    "element": f"Ligne_{r[0]}", "table": "T_Lines_pu", "index": i,
                    "bus": int(r[1]), "classe": "Ligne", "genre": "etat",
                    "variable": "iLd", "expression": "wb*(vd - R*iLd)/L",
                    "modifiable": False, "modifiee": False,
                    "pourquoi": "branche : ses equations sont celles du reseau, non modifiables",
                }
            )  # fmt: skip
        return out

    @staticmethod
    def _lignes_modele(m, nom, table, i, bus, classe, modifiee):
        base = {
            "element": nom, "table": table, "index": i, "bus": bus, "classe": classe,
            "modifiable": True, "modifiee": modifiee,
        }  # fmt: skip
        out = [
            {**base, "genre": "etat", "variable": e["nom"], "expression": e["derivee"]}
            for e in m["etats"]
        ]
        out += [
            {
                **base,
                "genre": "algebrique",
                "variable": a["nom"],
                "expression": a["expr"],
            }
            for a in m["algebriques"]
        ]
        out += [
            {**base, "genre": "sortie", "variable": k, "expression": v}
            for k, v in m["sorties"].items()
        ]
        return out

    # -- ecriture
    def _modele_de(self, t: dict, table: str, index: int) -> dict:
        if table == "T_Blocs":
            return t["MODELES_BLOCS"][t["T_Blocs"][index][2]]
        return _modele_gfl1()

    def _poser_modele(self, t: dict, table: str, index: int, m: dict) -> dict:
        """Le modele ``m`` pour l'element : natif converti (le modele garde
        son ``origine``, comme le vrai serveur), ou modele propre. Rend la
        reponse du contrat : ``{tables, modele, table, index}``."""
        if table != "T_Blocs":
            ligne = t[table].pop(index)
            self.compteur += 1
            nom = f"GFL1_{ligne[0]}"
            m["nom"] = nom
            m["origine"] = {
                "element": nom,
                "id_bloc": self.compteur + 100,
                "retraits": [{"table": table, "index": index, "ligne": ligne}],
                "remplacements": [],
            }
            t.setdefault("MODELES_BLOCS", {})[nom] = m
            t.setdefault("T_Blocs", []).append(
                [self.compteur + 100, ligne[1], nom, ligne[2]]
            )
            index = len(t["T_Blocs"]) - 1
        else:
            t["MODELES_BLOCS"][t["T_Blocs"][index][2]] = m
        return {"tables": t, "modele": m, "table": "T_Blocs", "index": index}

    def _revenir(self, t: dict, index: int) -> dict:
        """Le retour au natif du contrat : la ligne gardee dans ``origine``
        retrouve sa table, le modele s'en va (et ``MODELES_BLOCS`` vide n'est
        pas ecrit, comme ``tables_vers_json``)."""
        ligne = t["T_Blocs"][index]
        m = t["MODELES_BLOCS"][ligne[2]]
        origine = m.get("origine")
        if not origine or origine["id_bloc"] != ligne[0]:
            raise ValueError(f"{ligne[2]} n'a pas d'origine native connue")
        r = origine["retraits"][0]
        t["T_Blocs"].pop(index)
        k = min(r["index"], len(t[r["table"]]))
        t[r["table"]].insert(k, r["ligne"])
        del t["MODELES_BLOCS"][ligne[2]]
        if not t["MODELES_BLOCS"]:
            del t["MODELES_BLOCS"]
        return {"tables": t, "table": r["table"], "index": k}

    def gerer(self, route) -> None:
        req = route.request
        if req.method != "POST":
            return route.fallback()
        chemin = req.url.split("/api/equations/")[1].split("?")[0]
        corps = json.loads(req.post_data or "{}")
        self.journal.append((chemin, corps))
        if chemin in ("modifier", "ajouter", "supprimer", "parametre") and self.refus:
            msg, self.refus = self.refus, None
            return route.fulfill(status=422, json={"detail": msg})
        t = copy.deepcopy(corps.get("tables") or {})
        try:
            rep = self._repondre(chemin, corps, t)
        except ValueError as err:
            return route.fulfill(status=422, json={"detail": str(err)})
        route.fulfill(json=rep)

    def _repondre(self, chemin: str, corps: dict, t: dict) -> dict:
        if chemin == "reseau":
            return {"lignes": self._lignes(t)}
        if chemin == "verifier":
            return {"ok": True, "equations": []}
        if chemin == "depuis_schema":
            return {"modele": _modele_gfl1()}
        if chemin == "vers_blocs":
            return {"schema": {"nom": "x", "parametres": {}, "blocs": [], "liens": []}}
        table, index = corps["table"], corps["index"]
        if table not in ("T_GFL1s", "T_Blocs"):
            raise ValueError("cet element n'a pas d'equations modifiables")
        if chemin == "revenir":
            if table != "T_Blocs":
                raise ValueError("cet element est deja natif")
            return self._revenir(t, index)
        m = copy.deepcopy(self._modele_de(t, table, index))
        if chemin == "element":
            rep = {"modele": m}
            if table == "T_Blocs":
                origine = m.get("origine") or {}
                rep["natif"] = origine.get("element")
                rep["schema_bloc"] = False
            return rep
        if chemin == "parametre":
            if table != "T_Blocs":
                raise ValueError(
                    "element natif : ses parametres se reglent dans ses proprietes"
                )
            if corps["nom"] not in m["parametres"]:
                raise ValueError(f"pas de parametre {corps['nom']!r}")
            m["parametres"][corps["nom"]] = corps["valeur"]
            return self._poser_modele(t, table, index, m)
        v = corps["variable"]
        if chemin == "modifier":
            e = corps["expression"]
            if "inconnu" in e:
                raise ValueError(f"expression refusee : nom inconnu dans {e!r}")
            for x in m["etats"]:
                if x["nom"] == v:
                    x["derivee"] = e
            for x in m["algebriques"]:
                if x["nom"] == v:
                    x["expr"] = e
            if v in m["sorties"]:
                m["sorties"][v] = e
        elif chemin == "ajouter":
            if corps["genre"] == "etat":
                m["etats"].append({"nom": v, "derivee": corps["expression"]})
            else:
                m["algebriques"].append({"nom": v, "expr": corps["expression"]})
        elif chemin == "supprimer":
            m["etats"] = [x for x in m["etats"] if x["nom"] != v]
            m["algebriques"] = [x for x in m["algebriques"] if x["nom"] != v]
        return self._poser_modele(t, table, index, m)


def _poser_gfl1(page):
    """Un GFL1 sur le bus 2 de smib. smib est aussi le reseau de demarrage :
    on laisse finir son chargement, qui ecraserait sinon les tables modifiees."""
    page.wait_for_timeout(1500)
    page.evaluate(
        """() => {
            const t = E.tables, Sn = t.Sb;
            t.T_GFL1s = [[1, 2, Sn, 40e3]];
            t.Y_DER = t.Y_DER.map(r => +r[0] === 2 ? [...r.slice(0, 5), 0.3 * Sn, 0, ...r.slice(7)] : r);
            modifie();
        }"""
    )
    page.wait_for_timeout(700)  # l'instantane suivant ne se fond pas avec celui-ci


def _onglet_equations(page):
    page.click('nav button[data-onglet="equations"]')
    page.wait_for_selector("#eq-liste .eq-ligne", timeout=15000)


@pytest.fixture
def faux(page):
    _charger(page, "smib")
    _poser_gfl1(page)
    return FauxServeur(page)


def _expr(page, element, variable):
    return page.locator(
        f'#eq-liste tr.eq-ligne[data-element="{element}"][data-variable="{variable}"] td.eq-expr'
    )


def _confirmer(page):
    """La confirmation « fonctionnement libre » d'une conversion native -> modele."""
    page.locator('[data-confirmation="gardes"]').wait_for()
    page.click('[data-action="confirmer-gardes"]')


def _ecrire_et_valider(page, cellule, texte, confirmer=False):
    cellule.click()
    champ = cellule.locator("input.eq-champ")
    champ.fill(texte)
    champ.press("Enter")
    if confirmer:
        _confirmer(page)


# ------------------------------------------------------------ sans serveur
def test_onglet_liste_etats_seuls_puis_filtres(page, faux):
    """Etats seuls par defaut ; les genres se cochent ; filtre par
    equipement ; recherche ; ligne non modifiable grisee avec sa raison."""
    _onglet_equations(page)
    lignes = page.locator("#eq-liste tr.eq-ligne")
    # par defaut : les deux etats du GFL1 (dx/dt = f) et l'etat de la ligne
    genres = {g for g in page.eval_on_selector_all(
        "#eq-liste tr.eq-ligne", "rs => rs.map(r => r.dataset.genre)")}  # fmt: skip
    assert genres == {"etat"}
    assert page.locator('#eq-liste tr[data-element="GFL1_1"]').count() == 2
    # variable a gauche sous la forme dx/dt, expression a droite
    assert (
        page.locator('#eq-liste tr[data-variable="M_pll"] td.eq-var').inner_text()
        == "dM_pll/dt"
    )
    assert _expr(page, "GFL1_1", "M_pll").inner_text() == "Kipll*vq"
    # cocher les algebriques et les sorties
    page.check("#eq-g-algebrique")
    page.check("#eq-g-sortie")
    assert page.locator('#eq-liste tr[data-element="GFL1_1"]').count() == 5
    # filtre par equipement
    page.select_option("#eq-element", "GFL1_1")
    assert {r for r in page.eval_on_selector_all(
        "#eq-liste tr.eq-ligne", "rs => rs.map(r => r.dataset.element)")} == {"GFL1_1"}  # fmt: skip
    page.select_option("#eq-element", "")
    # recherche dans l'expression
    page.fill("#eq-recherche", "Kppll")
    page.wait_for_function(
        "document.querySelectorAll('#eq-liste tr.eq-ligne').length === 1"
    )
    assert lignes.first.get_attribute("data-variable") == "w_pll"
    page.fill("#eq-recherche", "")
    page.wait_for_function(
        "document.querySelectorAll('#eq-liste tr.eq-ligne').length > 1"
    )
    # la ligne de branche : grisee, raison en info-bulle, pas modifiable
    ligne = page.locator('#eq-liste tr[data-modifiable="0"]').first
    assert "eq-grise" in ligne.get_attribute("class")
    assert "non modifiables" in ligne.get_attribute("title")
    ligne.locator("td.eq-expr").click()
    assert ligne.locator("input").count() == 0
    # le compte est affiche
    assert "equation(s) sur" in page.inner_text("#eq-compte")


def test_edition_sur_place_echap_annule_entree_ecrit_annuler_rend(page, faux):
    _onglet_equations(page)
    cellule = _expr(page, "GFL1_1", "M_pll")
    # Echap : rien n'est envoye
    cellule.click()
    champ = cellule.locator("input.eq-champ")
    champ.fill("2*Kipll*vq")
    champ.press("Escape")
    assert cellule.locator("input").count() == 0
    assert cellule.inner_text() == "Kipll*vq"
    assert not [j for j in faux.journal if j[0] == "modifier"]
    # Entree : le corps suit le contrat
    _ecrire_et_valider(page, cellule, "2*Kipll*vq", confirmer=True)
    page.wait_for_function("(E.tables.T_Blocs || []).length === 1")
    envoi = next(c for n, c in faux.journal if n == "modifier")
    assert set(envoi) == {"tables", "table", "index", "variable", "expression"}
    assert (envoi["table"], envoi["index"], envoi["variable"]) == (
        "T_GFL1s",
        0,
        "M_pll",
    )
    assert envoi["expression"] == "2*Kipll*vq"
    # le natif a quitte sa table, le modele est dans le reseau
    assert page.evaluate("E.tables.T_GFL1s.length") == 0
    assert page.evaluate("Object.keys(E.tables.MODELES_BLOCS)") == ["GFL1_1"]
    # la liste se redessine : equipement marque, expression changee
    page.wait_for_selector('#eq-liste tr[data-element="Blocs_101"][data-modifiee="1"]')
    assert _expr(page, "Blocs_101", "M_pll").inner_text() == "2*Kipll*vq"
    assert page.locator("#eq-liste .eq-badge").count() >= 1
    # la validation du reseau repart
    page.wait_for_selector(
        '.eq-validation[data-etat="ok"], .eq-validation[data-etat="erreur"]'
    )
    # annuler rend le natif, retablir le re-ecrit
    page.wait_for_timeout(700)
    page.keyboard.press("Control+z")
    page.wait_for_function(
        "E.tables.T_GFL1s.length === 1 && !(E.tables.T_Blocs || []).length"
    )
    page.click('nav button[data-onglet="equations"]')
    page.wait_for_selector('#eq-liste tr[data-element="GFL1_1"]')
    assert _expr(page, "GFL1_1", "M_pll").inner_text() == "Kipll*vq"
    page.keyboard.press("Control+y")
    page.wait_for_function("(E.tables.T_Blocs || []).length === 1")


def test_conversion_native_avertit_fonctionnement_libre(page, faux):
    """Remplacer un natif par un modele (equations ou blocs) avertit que les
    gardes (limitation de courant, saturation de tension) ne sont pas reprises ;
    Annuler, Echap et le fond laissent l'element natif et n'appellent pas le serveur."""
    _onglet_equations(page)
    cellule = _expr(page, "GFL1_1", "M_pll")
    cellule.click()
    cellule.locator("input.eq-champ").fill("2*Kipll*vq")
    cellule.locator("input.eq-champ").press("Enter")
    boite = page.locator('[data-confirmation="gardes"]')
    boite.wait_for()
    texte = boite.inner_text()
    assert "fonctionnement libre" in texte
    assert "limitation de courant" in texte and "saturation de tension" in texte
    assert "ne se limitera plus en simulation temporelle" in texte
    # Annuler : rien n'est envoye, l'element reste natif, le champ montre pourquoi
    page.click('[data-action="annuler-gardes"]')
    assert boite.count() == 0
    assert not [j for j in faux.journal if j[0] == "modifier"]
    assert page.evaluate("E.tables.T_GFL1s.length") == 1
    page.wait_for_selector("#eq-liste .eq-erreur:not([hidden])")
    assert "reste natif" in page.inner_text("#eq-liste .eq-erreur")
    # Echap annule aussi, sans fermer autre chose
    cellule.locator("input.eq-champ").press("Enter")
    boite.wait_for()
    page.keyboard.press("Escape")
    assert boite.count() == 0 and not [j for j in faux.journal if j[0] == "modifier"]
    # « Remplacer par son modele a blocs » : meme avertissement, avant tout appel
    page.evaluate("ouvrirFiche('T_GFL1s', 0)")
    page.wait_for_selector('#fiche [data-action="remplacer-par-blocs"]', timeout=30000)
    page.click('#fiche [data-action="remplacer-par-blocs"]')
    boite.wait_for()
    assert (
        "modele a blocs" in boite.inner_text()
        and "limitation de courant" in boite.inner_text()
    )
    page.click('[data-action="annuler-gardes"]')
    assert page.evaluate("E.tables.T_GFL1s.length") == 1
    assert not page.evaluate("(E.tables.T_Blocs || []).length")


def test_expression_refusee_reste_dans_le_champ_tables_intactes(page, faux):
    _onglet_equations(page)
    avant = page.evaluate("JSON.stringify(E.tables)")
    cellule = _expr(page, "GFL1_1", "M_pll")
    _ecrire_et_valider(page, cellule, "inconnu*2", confirmer=True)
    page.wait_for_selector("#eq-liste .eq-erreur:not([hidden])")
    msg = page.inner_text("#eq-liste .eq-erreur")
    assert "nom inconnu" in msg and "M_pll" in msg
    # le champ garde ce qu'on a tape, les tables ne bougent pas
    assert cellule.locator("input.eq-champ").input_value() == "inconnu*2"
    assert page.evaluate("JSON.stringify(E.tables)") == avant
    # une expression vide est refusee avant tout appel
    n = len(faux.journal)
    cellule.locator("input.eq-champ").fill("  ")
    cellule.locator("input.eq-champ").press("Enter")
    assert "vide" in page.inner_text("#eq-liste .eq-erreur")
    assert len(faux.journal) == n
    cellule.locator("input.eq-champ").press("Escape")


def test_noms_et_expressions_jamais_interpretes_comme_html(page, faux):
    """Un nom ou une expression venus du serveur s'affichent comme du texte."""
    _onglet_equations(page)
    page.unroute("**/api/equations/*")
    piege = '<img src=x onerror="window.__pwned=1">'

    def repondre(route):
        if route.request.method != "POST":
            return route.fallback()
        route.fulfill(
            json={
                "lignes": [
                    {
                        "element": piege, "table": "T_GFL1s", "index": 0, "bus": 2,
                        "classe": "GFL1", "genre": "etat", "variable": piege,
                        "expression": piege + " + 1", "modifiable": True, "modifiee": True,
                    }
                ]
            }
        )  # fmt: skip

    page.route("**/api/equations/*", repondre)
    page.click("#eq-actualiser")
    # la liste precedente est encore la : attendre la NOUVELLE avant de regarder
    page.wait_for_function(
        "document.querySelector('#eq-liste').innerText.includes('onerror')"
    )
    assert page.locator("#eq-liste img").count() == 0
    assert page.evaluate("window.__pwned") is None
    assert piege in page.inner_text("#eq-liste")
    # idem dans la liste deroulante des equipements et en mode LaTeX
    assert page.locator("#eq-element img").count() == 0
    page.check("#eq-latex")
    assert page.locator("#eq-liste img").count() == 0
    # le champ d'edition porte le texte brut
    page.uncheck("#eq-latex")
    page.locator("#eq-liste td.eq-expr").click()
    assert piege in page.locator("#eq-liste input.eq-champ").input_value()
    page.locator("#eq-liste input.eq-champ").press("Escape")


def test_affichage_latex(page, faux):
    _onglet_equations(page)
    page.check("#eq-latex")
    assert (
        "\\frac{d M_{pll}}{dt}"
        in page.locator('#eq-liste tr[data-variable="M_pll"] td.eq-var').inner_text()
    )
    assert "\\cdot" in _expr(page, "GFL1_1", "M_pll").inner_text()
    page.click("#eq-copier")
    assert "\\frac{d M_{pll}}{dt} = " in page.inner_text("#eq-latex-sortie")
    page.uncheck("#eq-latex")
    assert (
        page.locator('#eq-liste tr[data-variable="M_pll"] td.eq-var').inner_text()
        == "dM_pll/dt"
    )


def test_route_absente_message_clair(page):
    _charger(page, "smib")
    page.route(
        "**/api/equations/reseau",
        lambda r: (
            r.fulfill(status=404, json={"detail": "Not Found"})
            if r.request.method == "POST"
            else r.fallback()
        ),
    )
    page.click('nav button[data-onglet="equations"]')
    # « calcul des equations... » s'affiche d'abord : attendre le message final
    page.wait_for_function(
        "document.querySelector('#eq-etat').textContent.includes('version du serveur')"
    )
    msg = page.inner_text("#eq-etat")
    assert "version du serveur" in msg and "/api/equations/reseau" in msg
    assert page.locator("#eq-liste .eq-ligne").count() == 0


def test_section_equations_fiche_et_proprietes(page, faux):
    """La section Equations de la fiche et des proprietes : modifier, ajouter
    un etat, ajouter une variable, supprimer, revenir au modele natif."""
    page.click('nav button[data-onglet="constructeur"]')
    page.evaluate("selectionner('T_GFL1s', 0)")
    sec = page.locator('#inspecteur [data-section-equations="inspecteur"]')
    sec.wait_for()
    sec.locator("summary").click()  # repliee par defaut dans les proprietes
    sec.locator("tr.eq-ligne").first.wait_for()
    assert sec.locator("tr.eq-ligne").count() == 5
    assert (
        sec.locator('[data-action="revenir-natif"]').count() == 0
    )  # natif : rien a rendre
    # les parametres d'un natif : lus, pas modifiables ici
    par = sec.locator("[data-section-parametres]")
    assert "(2)" in par.locator("summary").inner_text()
    par.locator("summary").click()
    assert par.locator("input").count() == 0
    assert "Kipll = 1200" in par.locator('[data-parametre="Kipll"]').inner_text()
    assert "Reglages" in par.inner_text()
    # modifier sur place
    cellule = sec.locator('tr[data-variable="M_pll"] td.eq-expr')
    cellule.click()
    cellule.locator("input.eq-champ").fill("Kipll*vq/2")
    cellule.locator("input.eq-champ").press("Enter")
    _confirmer(page)
    page.wait_for_function("(E.tables.T_Blocs || []).length === 1")
    # l'equipement a change de table : la selection le suit, la section reste ouverte
    assert page.evaluate("E.sel.table") == "T_Blocs"
    sec = page.locator('#inspecteur [data-section-equations="inspecteur"]')
    sec.locator('tr[data-variable="M_pll"]').wait_for()
    assert (
        sec.locator('tr[data-variable="M_pll"] td.eq-expr').inner_text() == "Kipll*vq/2"
    )
    # ajouter un etat
    sec.locator('[data-action="ajouter-etat"]').click()
    sec.locator('[data-champ="nom"]').fill("x_test")
    sec.locator('.eq-ajout [data-champ="expression"]').fill("wb*(w - 1)")
    sec.locator('[data-champ="depart"]').fill("0.5")
    sec.locator('[data-action="valider-ajout"]').click()
    sec.locator('tr[data-variable="x_test"]').wait_for()
    envoi = [c for n, c in faux.journal if n == "ajouter"][-1]
    assert (envoi["genre"], envoi["variable"], envoi["expression"], envoi["depart"]) == (
        "etat", "x_test", "wb*(w - 1)", 0.5,
    )  # fmt: skip
    assert (
        sec.locator('tr[data-variable="x_test"]').get_attribute("data-genre") == "etat"
    )
    # ajouter une variable algebrique
    sec.locator('[data-action="ajouter-variable"]').click()
    sec.locator('[data-champ="nom"]').fill("z_test")
    sec.locator('.eq-ajout [data-champ="expression"]').fill("2*V")
    sec.locator('[data-action="valider-ajout"]').click()
    sec.locator('tr[data-variable="z_test"]').wait_for()
    assert [c for n, c in faux.journal if n == "ajouter"][-1]["genre"] == "algebrique"
    # un nom invalide est refuse sur place, sans appel
    n = len(faux.journal)
    sec.locator('[data-action="ajouter-etat"]').click()
    sec.locator('[data-champ="nom"]').fill("2mauvais")
    sec.locator('[data-action="valider-ajout"]').click()
    assert "nom invalide" in sec.locator(".eq-ajout .pb").inner_text()
    assert len(faux.journal) == n
    sec.locator(
        '[data-action="ajouter-etat"]'
    ).click()  # referme et rouvre : formulaire vierge
    # supprimer
    sec.locator('[data-supprimer="x_test"]').click()
    page.wait_for_function(
        "![...document.querySelectorAll('#inspecteur tr.eq-ligne')].some(r => r.dataset.variable === 'x_test')"
    )
    # un parametre du modele converti se modifie sur place : le serveur l'ecrit
    sec = page.locator('#inspecteur [data-section-equations="inspecteur"]')
    champ = sec.locator('[data-section-parametres] input[data-parametre="Kipll"]')
    champ.wait_for()  # la section des parametres est restee ouverte
    assert champ.input_value() == "1200"
    champ.fill("2400")
    champ.press("Enter")
    page.wait_for_function("E.tables.MODELES_BLOCS.GFL1_1.parametres.Kipll === 2400")
    envoi = [c for n, c in faux.journal if n == "parametre"][-1]
    assert set(envoi) == {"tables", "table", "index", "nom", "valeur"}
    assert (envoi["table"], envoi["nom"], envoi["valeur"]) == ("T_Blocs", "Kipll", 2400)
    # un refus du serveur rend la valeur et dit pourquoi
    sec = page.locator('#inspecteur [data-section-equations="inspecteur"]')
    champ = sec.locator('[data-section-parametres] input[data-parametre="Kppll"]')
    champ.wait_for()
    faux.refus = "Kppll = 0 : division par zero"
    champ.fill("0")
    champ.press("Enter")
    sec.locator(".eq-param-erreur").wait_for()
    assert "division par zero" in sec.locator(".eq-param-erreur").inner_text()
    assert champ.input_value() == "60"
    assert page.evaluate("E.tables.MODELES_BLOCS.GFL1_1.parametres.Kppll") == 60
    # revenir au modele natif : c'est le SERVEUR qui rend la ligne d'origine
    page.wait_for_timeout(700)
    sec.locator('[data-action="revenir-natif"]').click()
    page.wait_for_function(
        "E.tables.T_GFL1s.length === 1 && !(E.tables.T_Blocs || []).length"
    )
    envoi = [c for n, c in faux.journal if n == "revenir"][-1]
    assert set(envoi) == {"tables", "table", "index"}
    assert (envoi["table"], envoi["index"]) == ("T_Blocs", 0)
    assert page.evaluate("E.tables.T_GFL1s[0]") == [
        1,
        2,
        page.evaluate("E.tables.Sb"),
        40e3,
    ]
    assert page.evaluate("Object.keys(E.tables.MODELES_BLOCS || {}).length") == 0
    assert page.evaluate("E.sel.table") == "T_GFL1s"
    # annuler le retour : l'element est de nouveau a equations
    page.wait_for_timeout(700)
    page.keyboard.press("Control+z")
    page.wait_for_function("(E.tables.T_Blocs || []).length === 1")


def _poser_modele_converti(page, origine: bool) -> None:
    """Le reseau tel qu'un fichier le rend : le GFL1 deja remplace par son
    modele a equations (rien n'en reste dans la memoire du navigateur)."""
    m = _modele_gfl1()
    m["nom"] = "GFL1_1"
    if origine:
        m["origine"] = {
            "element": "GFL1_1",
            "id_bloc": 7,
            "retraits": [
                {"table": "T_GFL1s", "index": 0, "ligne": [1, 2, 100e6, 33e3]}
            ],
            "remplacements": [],
        }
    page.evaluate(
        """m => {
            E.tables.T_GFL1s = [];
            E.tables.T_Blocs = [[7, 2, 'GFL1_1', 100e6]];
            E.tables.MODELES_BLOCS = {GFL1_1: m};
            modifie();
            selectionner('T_Blocs', 0);
        }""",
        m,
    )


def test_retour_au_natif_apres_rechargement_rend_la_vraie_ligne(page, faux):
    """Le retour au natif ne depend pas de la memoire du navigateur : apres un
    rechargement, c'est l'origine gardee DANS LE MODELE qui rend la ligne --
    avec ses valeurs (33 kV ici), pas une ligne par defaut de la classe."""
    page.click('nav button[data-onglet="constructeur"]')
    _poser_modele_converti(page, origine=True)
    sec = page.locator('#inspecteur [data-section-equations="inspecteur"]')
    sec.wait_for()
    sec.locator("summary").first.click()
    sec.locator('[data-action="revenir-natif"]').wait_for()
    assert (
        "remplace l'equipement natif GFL1_1" in sec.locator(".eq-resume").inner_text()
    )
    page.wait_for_timeout(700)
    sec.locator('[data-action="revenir-natif"]').click()
    page.wait_for_function("E.tables.T_GFL1s.length === 1 && !E.tables.T_Blocs.length")
    assert page.evaluate("E.tables.T_GFL1s[0]") == [1, 2, 100e6, 33e3]
    assert page.evaluate("Object.keys(E.tables.MODELES_BLOCS || {})") == []
    assert page.evaluate("[E.sel.table, E.sel.index]") == ["T_GFL1s", 0]


def test_modele_sans_origine_native_pas_de_retour(page, faux):
    """Un modele ecrit a la main ne remplace aucun equipement natif : la
    section le dit, et ne propose pas d'y « revenir »."""
    page.click('nav button[data-onglet="constructeur"]')
    _poser_modele_converti(page, origine=False)
    sec = page.locator('#inspecteur [data-section-equations="inspecteur"]')
    sec.wait_for()
    sec.locator("summary").first.click()
    sec.locator("tr.eq-ligne").first.wait_for()
    assert (
        "ne remplace aucun equipement natif" in sec.locator(".eq-resume").inner_text()
    )
    assert sec.locator('[data-action="revenir-natif"]').count() == 0


def test_section_equations_de_la_fiche_modifie_et_valide(page, faux):
    page.evaluate("ouvrirFiche('T_GFL1s', 0)")
    sec = page.locator('#fiche [data-section-equations="fiche"]')
    sec.locator("tr.eq-ligne").first.wait_for(timeout=30000)  # ouverte par defaut
    cellule = sec.locator('tr[data-variable="theta_pll"] td.eq-expr')
    cellule.click()
    cellule.locator("input.eq-champ").fill("wb*(w_pll - w)*1.1")
    cellule.locator("input.eq-champ").press("Enter")
    _confirmer(page)
    page.wait_for_function("E.ficheRef.table === 'T_Blocs'")
    sec = page.locator('#fiche [data-section-equations="fiche"]')
    sec.locator('tr[data-variable="theta_pll"]').wait_for()
    assert (
        sec.locator('tr[data-variable="theta_pll"] td.eq-expr').inner_text()
        == "wb*(w_pll - w)*1.1"
    )
    # un refus du serveur reste dans la fiche, sans rien changer
    faux.refus = "expression refusee : division par zero dans 1/0"
    c2 = sec.locator('tr[data-variable="M_pll"] td.eq-expr')
    c2.click()
    c2.locator("input.eq-champ").fill("1/0")
    c2.locator("input.eq-champ").press("Enter")
    sec.locator(".eq-erreur:not([hidden])").wait_for()
    assert "division par zero" in sec.locator(".eq-erreur").inner_text()
    page.locator('#fiche [data-section-equations="fiche"] input.eq-champ').press(
        "Escape"
    )
    assert page.locator("#fiche").is_visible()  # Echap n'a ferme que le champ


def test_modele_a_equations_dans_l_onglet_blocs(page, faux):
    """Un modele a equations ouvert dans l'onglet Blocs : panneau d'equations
    (au lieu du dessin), correction verifiee par le serveur, passage en blocs."""
    page.evaluate(
        """m => { E.tables.MODELES_BLOCS = {GFL1_mod: m}; ouvrirModeleBlocs('GFL1_mod'); modifie(); }""",
        _modele_gfl1(),
    )
    page.wait_for_timeout(700)
    page.click('nav button[data-onglet="blocs"]')
    page.wait_for_selector("#blocs-eq-modele tr.eq-ligne")
    assert page.locator("#blocs-espace").is_hidden()
    assert page.inner_text("#btn-blocs-convertir") == "Convertir en blocs"
    # correction verifiee a blanc par /verifier
    cellule = page.locator('#blocs-eq-modele tr[data-variable="M_pll"] td.eq-expr')
    cellule.click()
    cellule.locator("input.eq-champ").fill("3*Kipll*vq")
    cellule.locator("input.eq-champ").press("Enter")
    page.wait_for_function(
        "E.tables.MODELES_BLOCS.GFL1_mod.etats[1].derivee === '3*Kipll*vq'"
    )
    assert any(n == "verifier" and "modele" in c for n, c in faux.journal)
    # un parametre du modele : verifie a blanc lui aussi, puis garde
    par = page.locator("#blocs-eq-modele [data-section-parametres]")
    par.locator("summary").click()
    champ = par.locator('input[data-parametre="Kppll"]')
    champ.fill("75")
    champ.press("Enter")
    page.wait_for_function("E.tables.MODELES_BLOCS.GFL1_mod.parametres.Kppll === 75")
    verif = [c for n, c in faux.journal if n == "verifier"][-1]
    assert verif["modele"]["parametres"]["Kppll"] == 75
    # compiler : la compilation a blanc
    page.click("#btn-blocs-compiler")
    page.wait_for_selector("#blocs-eq-modele .eq-verif.ok")
    # le passage en blocs cree un nouveau modele, ouvert comme schema-bloc
    page.click("#btn-blocs-convertir")
    page.wait_for_function("B.modele === 'GFL1_mod_blocs'")
    assert page.locator("#blocs-eq-modele").is_hidden()
    assert page.locator("#blocs-espace").is_visible()
    assert page.inner_text("#btn-blocs-convertir") == "Voir en equations"
    assert "GFL1_mod" in page.evaluate(
        "Object.keys(E.tables.MODELES_BLOCS)"
    )  # l'original reste


# ------------------------------------------------- contre le vrai serveur
def _serveur_fr052(page) -> bool:
    return page.evaluate(
        """async () => {
            const r = await fetch('/api/equations/reseau', {method: 'POST',
                headers: {'Content-Type': 'application/json'}, body: JSON.stringify({tables: E.tables})});
            return r.status !== 404 && r.status !== 405;
        }"""
    )


@pytest.fixture
def page_fr052(page):
    _charger(page, "wscc9_gfm")
    if not _serveur_fr052(page):
        pytest.skip("le serveur n'a pas encore les routes de FR-052")
    return page


def _spectre(page):
    """Le spectre d'un calcul modal, tel que l'atelier le recoit."""
    d = page.evaluate("api('/api/modal', corpsReseau({f_max: null, n_etats: 0}))")
    cles = [k for k in d if k in ("valeurs_propres", "modes", "eigs", "lambda")]
    assert cles, list(d)
    return json.dumps(d[cles[0]], sort_keys=True)


def test_critere8_onglet_modifier_modal_annuler(page_fr052):
    """Critere 8 (1/3) : multiplier par 2 un facteur d'une derivee dans la
    liste globale, relancer l'analyse modale (le spectre change), annuler (il
    revient)."""
    page = page_fr052
    page.wait_for_timeout(700)
    avant = _spectre(page)
    page.click('nav button[data-onglet="equations"]')
    page.wait_for_selector("#eq-liste tr.eq-ligne", timeout=60000)
    ligne = page.locator('#eq-liste tr[data-modifiable="1"][data-genre="etat"]').first
    element = ligne.get_attribute("data-element")
    variable = ligne.get_attribute("data-variable")
    expression = ligne.locator("td.eq-expr").inner_text()
    ligne.locator("td.eq-expr").click()
    champ = ligne.locator("input.eq-champ")
    champ.fill(f"2*({expression})")
    champ.press("Enter")
    _confirmer(page)
    page.wait_for_selector(
        f'#eq-liste tr[data-element="{element}"][data-variable="{variable}"][data-modifiee="1"]',
        timeout=60000,
    )
    page.wait_for_selector(
        '.eq-validation[data-etat="ok"], .eq-validation[data-etat="erreur"]',
        timeout=60000,
    )
    apres = _spectre(page)
    assert apres != avant, "doubler une derivee doit changer le spectre"
    page.wait_for_timeout(700)
    page.keyboard.press("Control+z")
    page.wait_for_function("!(E.tables.T_Blocs || []).length", timeout=30000)
    assert _spectre(page) == avant


def _un_gfm(page) -> tuple[str, int]:
    return page.evaluate(
        """() => {
            for (const t of Object.keys(E.tables)) if (/^T_GFM\\ds$/.test(t) && E.tables[t].length) return [t, 0];
            return null;
        }"""
    )


def test_critere8_fiche_gfm_ajouter_un_etat_revenir_au_natif(page_fr052):
    """Critere 8 (2/3) : modifier depuis la fiche d'un GFM, ajouter un etat,
    revenir au natif (spectre de depart rendu)."""
    page = page_fr052
    table, i = _un_gfm(page)
    page.wait_for_timeout(700)
    avant = _spectre(page)
    page.evaluate(f"ouvrirFiche('{table}', {i})")
    sec = page.locator('#fiche [data-section-equations="fiche"]')
    sec.locator("tr.eq-ligne").first.wait_for(timeout=60000)
    cellule = sec.locator('tr[data-genre="etat"] td.eq-expr').first
    exp = cellule.inner_text()
    cellule.click()
    cellule.locator("input.eq-champ").fill(f"({exp})*1.0")
    cellule.locator("input.eq-champ").press("Enter")
    _confirmer(page)
    page.wait_for_function("E.ficheRef.table === 'T_Blocs'", timeout=60000)
    sec = page.locator('#fiche [data-section-equations="fiche"]')
    # la section recharge sa liste apres le changement de table : compter
    # avant qu'elle soit la donnait 0
    sec.locator('tr[data-genre="etat"]').first.wait_for(timeout=60000)
    n0 = sec.locator('tr[data-genre="etat"]').count()
    sec.locator('[data-action="ajouter-etat"]').click()
    sec.locator('[data-champ="nom"]').fill("x_ajoute")
    sec.locator('.eq-ajout [data-champ="expression"]').fill("-x_ajoute")
    sec.locator('[data-action="valider-ajout"]').click()
    sec.locator('tr[data-variable="x_ajoute"]').wait_for(timeout=60000)
    assert sec.locator('tr[data-genre="etat"]').count() == n0 + 1
    sec.locator('[data-action="revenir-natif"]').click()
    page.wait_for_function(
        f"E.tables.{table}.length > 0 && !(E.tables.T_Blocs || []).length"
    )
    assert _spectre(page) == avant


def _tensions(page) -> list[float]:
    d = page.evaluate("api('/api/loadflow', corpsReseau())")
    return [b["V"] for b in d["bus"]]


def test_charge_remplacee_par_blocs_garde_le_load_flow_et_revient(page):
    """Le defaut trouve par l'audit, de bout en bout contre le vrai serveur :
    « Remplacer par son modele a blocs » sur une charge de ``T_Loads`` retirait
    sa ligne dans le navigateur SANS reporter sa puissance a son bus (la charge
    sortait du load flow), et « Revenir au modele natif » posait une ligne par
    defaut. Le serveur fait les deux : le load flow ne bouge pas, et le retour
    rend les tables de depart, a l'identique."""
    _charger(page, "cigre_b4", _nav.B4)
    if not _serveur_fr052(page):
        pytest.skip("le serveur n'a pas encore les routes de FR-052")
    page.wait_for_timeout(700)
    avant = page.evaluate("JSON.stringify(E.tables)")
    v0 = _tensions(page)
    i = page.evaluate("E.tables.T_Loads.length") - 1  # la charge seule a son bus
    ligne = page.evaluate(f"E.tables.T_Loads[{i}]")
    page.evaluate(f"ouvrirFiche('T_Loads', {i})")
    page.wait_for_selector('#fiche [data-action="remplacer-par-blocs"]', timeout=60000)
    page.click('#fiche [data-action="remplacer-par-blocs"]')
    _confirmer(page)
    page.wait_for_function("(E.tables.T_Blocs || []).length === 1", timeout=60000)
    assert page.evaluate("E.tables.T_Loads.length") == i
    modele = page.evaluate("E.tables.T_Blocs[0][2]")
    origine = page.evaluate(f"E.tables.MODELES_BLOCS['{modele}'].origine")
    assert origine["retraits"][0]["ligne"] == ligne
    assert origine["remplacements"], "la puissance de la charge est reportee a son bus"
    v1 = _tensions(page)
    assert max(abs(a - b) for a, b in zip(v0, v1)) < 1e-9
    # retour au natif depuis les proprietes de l'element
    page.click('nav button[data-onglet="constructeur"]')
    page.evaluate("selectionner('T_Blocs', 0)")
    sec = page.locator('#inspecteur [data-section-equations="inspecteur"]')
    sec.wait_for()
    sec.locator("summary").first.click()
    sec.locator('[data-action="revenir-natif"]').wait_for(timeout=60000)
    page.wait_for_timeout(700)
    sec.locator('[data-action="revenir-natif"]').click()
    page.wait_for_function("!(E.tables.T_Blocs || []).length", timeout=60000)
    apres = json.loads(page.evaluate("JSON.stringify(E.tables)"))
    if not apres.get("MODELES_BLOCS"):
        apres.pop("MODELES_BLOCS", None)  # vide ou absent : la meme chose
    assert apres == json.loads(avant)
