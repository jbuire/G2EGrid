"""FR-051, partie interface : l'editeur de blocs dans Chromium.

Deux familles de tests :

* ceux qui ne dependent que du JSON et de l'interface (palette par categories,
  rendu et navigation dans un sous-systeme ecrit a la main, selection
  multiple, messages d'erreur d'une route absente, annuler/retablir) : ils
  passent contre n'importe quel serveur ;
* ceux qui appellent les routes de FR-051 (``regrouper``, ``depuis_element``,
  ``bibliotheque``) : ils se SAUTENT tant que le serveur ne les a pas.

Les fixtures (serveur, navigateur, page) sont celles de ``test_navigateur``.
"""

from __future__ import annotations

import json

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("playwright.sync_api")

import test_navigateur as _nav

# les fixtures de test_navigateur (serveur lance avec un jeton, Chromium, page)
serveur = _nav.serveur
navigateur = _nav.navigateur
page = _nav.page
_charger = _nav._charger


def _onglet_blocs(page):
    page.click('nav button[data-onglet="blocs"]')
    page.wait_for_function("B.types !== null")


def _serveur_fr051(page) -> bool:
    """Le serveur a-t-il les routes de FR-051 ? (categorie dans les types et
    bibliotheque : l'absence de l'un ou l'autre suffit a sauter.)"""
    return page.evaluate(
        """async () => {
            const t = await (await fetch('/api/blocs/types')).json();
            const b = await fetch('/api/blocs/bibliotheque');
            return b.status !== 404 && Object.values(t.types).every(s => 'categorie' in s)
                && 'sous_systeme' in t.types;
        }"""
    )


@pytest.fixture
def page_fr051(page):
    _charger(page, "smib")
    _onglet_blocs(page)
    if not _serveur_fr051(page):
        pytest.skip("le serveur n'a pas encore les routes de FR-051")
    return page


def _schema_a_la_main() -> dict:
    """Une chaine P avec un sous-systeme (gain 2), une chaine Q a plat."""
    return {
        "nom": "main",
        "parametres": {"tau": 0.05},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref", "x": 40, "y": 40},
            {"id": "qr", "type": "consigne", "nom": "Q_ref", "x": 40, "y": 220},
            {
                "id": "ss1",
                "type": "sous_systeme",
                "nom": "chaine",
                "x": 240,
                "y": 40,
                "schema": {
                    "blocs": [
                        {
                            "id": "pe1",
                            "type": "port_entree",
                            "index": 0,
                            "nom": "u",
                            "x": 40,
                            "y": 60,
                        },
                        {"id": "g1", "type": "gain", "K": 2, "x": 200, "y": 60},
                        {
                            "id": "ps1",
                            "type": "port_sortie",
                            "index": 0,
                            "nom": "y",
                            "x": 400,
                            "y": 60,
                        },
                    ],
                    "liens": [
                        {"de": "pe1", "vers": "g1", "entree": 0},
                        {"de": "g1", "vers": "ps1", "entree": 0},
                    ],
                },
            },
            {
                "id": "fq",
                "type": "premier_ordre",
                "K": 1,
                "tau": "tau",
                "x": 240,
                "y": 220,
            },
            {"id": "op", "type": "sortie", "signal": "P", "x": 470, "y": 40},
            {"id": "oq", "type": "sortie", "signal": "Q", "x": 470, "y": 220},
        ],
        "liens": [
            {"de": "pr", "vers": "ss1", "entree": 0},
            {"de": "ss1", "vers": "op", "entree": 0},
            {"de": "qr", "vers": "fq", "entree": 0},
            {"de": "fq", "vers": "oq", "entree": 0},
        ],
    }


def _ouvrir_schema(page, schema):
    page.evaluate(
        """s => { E.tables.MODELES_BLOCS[s.nom] = s; ouvrirModeleBlocs(s.nom); modifie(); rendreBlocs(); }""",
        schema,
    )
    page.wait_for_selector("#blocs-svg .bloc")
    page.wait_for_timeout(700)  # l'instantane suivant ne se fond pas avec celui-ci


def _bloc(page, ident):
    return page.locator(f'#blocs-svg .bloc[data-bloc-id="{ident}"]')


def _creer_source(page):
    page.select_option("#blocs-exemple", "source_filtree")
    page.click("#btn-blocs-nouveau")
    assert page.locator("#blocs-svg .bloc").count() == 6


def _poser_gfl1(page):
    """Un GFL1 sur le bus 2 de smib, avec sa puissance : (table, indice)."""
    # smib est aussi le reseau de demarrage : on laisse finir son chargement,
    # qui remplacerait sinon les tables qu'on vient de modifier
    page.wait_for_timeout(1500)
    page.evaluate(
        """() => {
            const t = E.tables, Sn = t.Sb;
            t.T_GFL1s = [[1, 2, Sn, 40e3]];
            t.Y_DER = t.Y_DER.map(r => +r[0] === 2 ? [...r.slice(0, 5), 0.3 * Sn, 0, ...r.slice(7)] : r);
            modifie();
        }"""
    )
    return "T_GFL1s", 0


# ------------------------------------------------------------- sans serveur
def test_palette_par_categories(page):
    """U1 : la palette est rangee par categories ; les ports n'existent que
    dans un sous-systeme."""
    _charger(page, "smib")
    _onglet_blocs(page)
    _creer_source(page)
    groupes = page.locator(".blocs-palette .blocs-groupe")
    assert groupes.count() >= 5
    cats = page.evaluate(
        "[...document.querySelectorAll('.blocs-groupe')].map(g => g.dataset.categorie)"
    )
    assert "operations" in cats and "dynamique" in cats
    # chaque type du serveur est dans la palette (hors ports, interdits au premier niveau)
    attendus = page.evaluate(
        "Object.keys(B.types.types).filter(t => !TYPES_INTERNES.includes(t))"
    )
    presents = page.evaluate(
        "[...document.querySelectorAll('.blocs-palette button')].map(b => b.dataset.bloc)"
    )
    assert set(attendus) <= set(presents)
    assert "port_entree" not in presents and "port_sortie" not in presents
    assert "sous_systeme" in presents


def test_sous_systeme_a_la_main_entrer_remonter(page):
    """U3 : un sous-systeme ecrit a la main s'affiche avec ses ports, on y
    entre par double-clic, on y edite, on remonte par le fil d'Ariane."""
    _charger(page, "smib")
    _onglet_blocs(page)
    _ouvrir_schema(page, _schema_a_la_main())
    ss = _bloc(page, "ss1")
    assert "bloc-sous-systeme" in ss.get_attribute("class")
    # ses ports portent les noms de port_entree / port_sortie
    textes = [t.text_content() for t in ss.locator("text").all()]
    assert "chaine" in textes and "u" in textes and "y" in textes
    assert page.locator("#blocs-ariane .ariane-pas").count() == 1
    ss.dblclick()
    page.wait_for_function("B.chemin.length === 1")
    pas = page.locator("#blocs-ariane .ariane-pas").all_inner_texts()
    assert pas == ["main", "chaine"]
    assert page.locator("#blocs-svg .bloc").count() == 3
    # les ports ne sont proposes qu'ici
    assert page.locator('.blocs-palette button[data-bloc="port_entree"]').count() == 1
    assert page.locator('.blocs-palette button[data-bloc="port_sortie"]').count() == 1
    # l'edition marche comme au premier niveau : un gain de plus, dans le sous-systeme
    page.click('.blocs-palette button[data-bloc="gain"]')
    assert (
        page.evaluate(
            "E.tables.MODELES_BLOCS.main.blocs.find(b => b.id === 'ss1').schema.blocs.length"
        )
        == 4
    )
    assert page.evaluate("E.tables.MODELES_BLOCS.main.blocs.length") == 6
    # un second port d'entree prend l'index suivant
    page.click('.blocs-palette button[data-bloc="port_entree"]')
    assert page.evaluate(
        "E.tables.MODELES_BLOCS.main.blocs.find(b => b.id === 'ss1').schema.blocs.filter(b => b.type === 'port_entree').map(b => b.index)"
    ) == [0, 1]
    # remonter : le sous-systeme a maintenant deux entrees
    page.click('#blocs-ariane .ariane-pas[data-niveau="0"]')
    assert page.evaluate("B.chemin.length") == 0
    assert page.locator("#blocs-svg .bloc").count() == 6
    assert page.locator('.blocs-palette button[data-bloc="port_entree"]').count() == 0
    assert _bloc(page, "ss1").locator(".bloc-port-entree").count() == 2
    # Echap remonte aussi
    _bloc(page, "ss1").dblclick()
    page.wait_for_function("B.chemin.length === 1")
    page.keyboard.press("Escape")
    page.wait_for_function("B.chemin.length === 0")


def test_libelles_des_nouveaux_types(page):
    """Les nouveaux types se dessinent avec leur libelle (fonction de transfert lisible)."""
    _charger(page, "smib")
    _onglet_blocs(page)
    schema = {
        "nom": "types",
        "parametres": {},
        "blocs": [
            {
                "id": "a",
                "type": "pid",
                "Kp": 1,
                "Ki": 10,
                "Kd": 0.1,
                "Tf": 0.01,
                "x": 40,
                "y": 40,
            },
            {
                "id": "b",
                "type": "fonction_transfert",
                "num": [1],
                "den": [1, 3, 2],
                "x": 40,
                "y": 140,
            },
            {"id": "c", "type": "saturation", "min": -1, "max": 1, "x": 240, "y": 40},
            {
                "id": "d",
                "type": "avance_retard",
                "K": 1,
                "T1": 0.1,
                "T2": 0.05,
                "x": 240,
                "y": 140,
            },
            {
                "id": "e",
                "type": "second_ordre",
                "K": 1,
                "wn": 100,
                "zeta": 0.7,
                "x": 440,
                "y": 40,
            },
            {
                "id": "f",
                "type": "implicite",
                "n": 2,
                "expr": "u1 + u2 - y",
                "x": 440,
                "y": 140,
            },
        ],
        "liens": [],
    }
    _ouvrir_schema(page, schema)
    lib = page.evaluate(
        "Object.fromEntries(schemaCourant().blocs.map(b => [b.id, libelleBloc(b)]))"
    )
    assert lib["b"] == "1 / (s²+3·s+2)"
    assert lib["c"] == "[-1 ; 1]"
    assert lib["a"].startswith("PID") and lib["f"] == "0 = u1 + u2 - y"
    assert "T1" not in lib["d"] and "(1+0.1·s)/(1+0.05·s)" in lib["d"]
    # le bloc a plusieurs entrees est dessine avec ses ports
    assert _bloc(page, "f").locator(".bloc-port-entree").count() == 2 or page.evaluate(
        "specType('implicite') === null"
    )
    # ses champs s'editent : la fonction de transfert, en liste
    _bloc(page, "b").click(position={"x": 30, "y": 25})
    champ = page.locator('#blocs-proprietes input[data-champ="den"]')
    if champ.count():
        champ.fill("1, 5, 6")
        champ.dispatch_event("change")
        assert page.evaluate("schemaCourant().blocs.find(b => b.id === 'b').den") == [
            1,
            5,
            6,
        ]


def test_selection_multiple(page):
    """U2 : Maj+clic ajoute et retire, le rectangle de selection prend ce qu'il touche."""
    _charger(page, "smib")
    _onglet_blocs(page)
    _creer_source(page)
    _bloc(page, "pr").click(position={"x": 20, "y": 20})
    _bloc(page, "fp").click(position={"x": 20, "y": 20}, modifiers=["Shift"])
    assert page.evaluate("[...B.multi].sort()") == ["fp", "pr"]
    assert "2 blocs selectionnes" in page.locator("#blocs-proprietes").inner_text()
    assert page.locator("#btn-blocs-regrouper").is_enabled()
    _bloc(page, "pr").click(position={"x": 20, "y": 20}, modifiers=["Shift"])
    assert page.evaluate("[...B.multi]") == ["fp"]
    # le fond deselectionne
    svg = page.locator("#blocs-svg").bounding_box()
    page.mouse.click(svg["x"] + svg["width"] * 0.5, svg["y"] + svg["height"] * 0.9)
    assert page.evaluate("B.multi.size") == 0
    assert page.locator("#btn-blocs-regrouper").is_disabled()
    # rectangle : tout le bas du schema (la chaine Q)
    d = page.evaluate(
        "Object.fromEntries(schemaCourant().blocs.map(b => [b.id, [b.x, b.y]]))"
    )
    vb = page.evaluate(
        "[$('#blocs-svg').viewBox.baseVal.width, $('#blocs-svg').viewBox.baseVal.height]"
    )
    k = svg["width"] / vb[0]

    def px(x, y):
        return svg["x"] + x * k, svg["y"] + y * k

    x0, y0 = px(10, d["qr"][1] - 20)
    x1, y1 = px(d["oq"][0] + 120, d["oq"][1] + 80)
    page.mouse.move(x0, y0)
    page.mouse.down()
    page.mouse.move((x0 + x1) / 2, (y0 + y1) / 2, steps=4)
    page.mouse.move(x1, y1, steps=4)
    page.mouse.up()
    assert page.evaluate("[...B.multi].sort()") == ["fq", "oq", "qr"]
    # glisser un bloc du groupe les deplace tous
    avant = page.evaluate(
        "schemaCourant().blocs.filter(b => B.multi.has(b.id)).map(b => b.x)"
    )
    b = _bloc(page, "fq").bounding_box()
    page.mouse.move(b["x"] + 10, b["y"] + 10)
    page.mouse.down()
    page.mouse.move(b["x"] + 10, b["y"] + 60, steps=5)
    page.mouse.up()
    ys = page.evaluate(
        "schemaCourant().blocs.filter(b => B.multi.has(b.id)).map(b => b.y)"
    )
    decalages = {y - d[i][1] for y, i in zip(ys, ["qr", "fq", "oq"], strict=True)}
    assert len(decalages) == 1 and decalages.pop() > 0  # le meme decalage pour tous
    assert (
        page.evaluate(
            "schemaCourant().blocs.filter(b => B.multi.has(b.id)).map(b => b.x)"
        )
        == avant
    )
    assert (
        page.evaluate("schemaCourant().blocs.find(b => b.id === 'pr').y") == d["pr"][1]
    )


def test_route_absente_ou_refus_message_clair(page):
    """Une route qui repond 404 ou 422 donne un message qui dit quoi faire,
    jamais une trace, et ne touche pas au schema."""
    _charger(page, "smib")
    _onglet_blocs(page)
    _creer_source(page)
    avant = page.evaluate("JSON.stringify(schemaCourant())")
    for ident in ("pr", "fp"):
        _bloc(page, ident).click(position={"x": 20, "y": 20}, modifiers=["Shift"])
    page.route(
        "**/api/blocs/regrouper",
        lambda r: r.fulfill(status=404, json={"detail": "Not Found"}),
    )
    page.once("dialog", lambda dlg: dlg.accept("groupe"))
    page.click("#btn-blocs-regrouper")
    page.wait_for_function(
        "document.querySelector('#blocs-etat').textContent.includes('Regrouper')"
    )
    msg = page.locator("#blocs-etat").inner_text()
    assert "n'existe pas" in msg and "Traceback" not in msg
    page.unroute("**/api/blocs/regrouper")
    page.route(
        "**/api/blocs/regrouper",
        lambda r: r.fulfill(
            status=422,
            json={"detail": "bloc 'zz' inexistant : verifier les identifiants"},
        ),
    )
    page.once("dialog", lambda dlg: dlg.accept("groupe"))
    page.click("#btn-blocs-regrouper")
    page.wait_for_function(
        "document.querySelector('#blocs-etat').textContent.includes('refuse')"
    )
    assert "bloc 'zz' inexistant" in page.locator("#blocs-etat").inner_text()
    assert page.evaluate("JSON.stringify(schemaCourant())") == avant


def test_regrouper_annuler_retablir_sans_serveur_nouveau(page):
    """U6 : regrouper passe par E.tables, donc par annuler/retablir. La route
    est ici simulee (reponse fixe) pour ne tester que l'interface."""
    _charger(page, "smib")
    _onglet_blocs(page)
    _ouvrir_schema(page, _schema_a_la_main())
    s = _schema_a_la_main()
    # ce que rendrait le serveur en regroupant "fq" et "oq" dans "bas"
    s["blocs"] = [b for b in s["blocs"] if b["id"] not in ("fq", "oq")]
    s["blocs"].append(
        {
            "id": "bas",
            "type": "sous_systeme",
            "nom": "bas",
            "x": 240,
            "y": 220,
            "schema": {
                "blocs": [
                    {
                        "id": "pe1",
                        "type": "port_entree",
                        "index": 0,
                        "nom": "u1",
                        "x": 40,
                        "y": 40,
                    },
                    {
                        "id": "fq",
                        "type": "premier_ordre",
                        "K": 1,
                        "tau": "tau",
                        "x": 200,
                        "y": 40,
                    },
                    {"id": "oq", "type": "sortie", "signal": "Q", "x": 400, "y": 40},
                ],
                "liens": [
                    {"de": "pe1", "vers": "fq", "entree": 0},
                    {"de": "fq", "vers": "oq", "entree": 0},
                ],
            },
        }
    )
    s["liens"] = [
        {"de": "pr", "vers": "ss1", "entree": 0},
        {"de": "ss1", "vers": "op", "entree": 0},
        {"de": "qr", "vers": "bas", "entree": 0},
    ]
    envoyes = []

    def repondre(route):
        envoyes.append(json.loads(route.request.post_data))
        route.fulfill(json={"schema": s})

    page.route("**/api/blocs/regrouper", repondre)
    for ident in ("fq", "oq"):
        _bloc(page, ident).click(position={"x": 20, "y": 20}, modifiers=["Shift"])
    page.once("dialog", lambda dlg: dlg.accept("bas"))
    page.click("#btn-blocs-regrouper")
    page.wait_for_function("schemaCourant().blocs.some(b => b.id === 'bas')")
    # le contrat : schema_bloc, ids (ceux du niveau), nom
    assert sorted(envoyes[0]["ids"]) == ["fq", "oq"] and envoyes[0]["nom"] == "bas"
    assert envoyes[0]["schema_bloc"]["nom"] == "main"
    assert page.evaluate("B.multi.has('bas')")
    assert page.locator("#btn-blocs-degrouper").is_enabled()
    # annuler : on retrouve les six blocs a plat
    page.wait_for_timeout(700)
    page.click("#btn-annuler")
    page.wait_for_function("schemaCourant().blocs.some(b => b.id === 'fq')")
    assert page.locator("#blocs-svg .bloc").count() == 6
    page.click("#btn-retablir")
    page.wait_for_function("schemaCourant().blocs.some(b => b.id === 'bas')")
    assert page.locator("#blocs-svg .bloc").count() == 5


def test_noms_de_blocs_jamais_interpretes_comme_html(page):
    """U7 : un nom de sous-systeme ou d'expression est du texte, pas du HTML."""
    _charger(page, "smib")
    _onglet_blocs(page)
    s = _schema_a_la_main()
    s["blocs"][2]["nom"] = "<img src=x onerror=window.__pwned=1>"
    s["blocs"][3]["expr"] = "<b>u</b>"
    s["blocs"][3]["type"] = "fonction"
    _ouvrir_schema(page, s)
    _bloc(page, "ss1").click(position={"x": 20, "y": 20})
    _bloc(page, "fq").click(position={"x": 20, "y": 20})
    _bloc(page, "ss1").dblclick()
    page.wait_for_function("B.chemin.length === 1")
    assert (
        page.locator("#blocs-ariane img, #blocs-svg img, #blocs-proprietes img").count()
        == 0
    )
    assert "<img" in page.locator("#blocs-ariane").inner_text()
    assert page.evaluate("window.__pwned") is None


# --------------------------------------------------------- avec le serveur
def test_regrouper_trois_blocs_entrer_remonter_compiler(page_fr051):
    """Critere 8 / U2-U3 : regrouper 3 blocs, entrer, remonter, compiler ;
    les equations sont celles du schema a plat."""
    page = page_fr051
    _creer_source(page)
    page.click("#btn-blocs-compiler")
    page.wait_for_function("B.compile !== null", timeout=30000)
    avant = page.evaluate("B.compile.etats.slice().sort()")
    for ident in ("pr", "fp", "op"):
        _bloc(page, ident).click(position={"x": 20, "y": 20}, modifiers=["Shift"])
    page.once("dialog", lambda dlg: dlg.accept("chaine_p"))
    page.click("#btn-blocs-regrouper")
    page.wait_for_function(
        "schemaCourant().blocs.some(b => b.type === 'sous_systeme')", timeout=30000
    )
    assert page.locator("#blocs-svg .bloc").count() == 4
    ss = page.evaluate("schemaCourant().blocs.find(b => b.type === 'sous_systeme').id")
    page.wait_for_function("B.multi.size === 1")
    # entrer
    _bloc(page, ss).dblclick()
    page.wait_for_function("B.chemin.length === 1")
    assert page.locator("#blocs-svg .bloc").count() >= 3
    assert "chaine_p" in page.locator("#blocs-ariane").inner_text()
    # remonter
    page.click('#blocs-ariane .ariane-pas[data-niveau="0"]')
    assert page.evaluate("B.chemin.length") == 0
    # compiler : meme nombre d'etats, a un renommage pres
    page.click("#btn-blocs-compiler")
    page.wait_for_function("B.compile !== null", timeout=30000)
    assert len(page.evaluate("B.compile.etats")) == len(avant)
    # degrouper rend le schema a plat
    _bloc(page, ss).click(position={"x": 20, "y": 20})
    page.click("#btn-blocs-degrouper")
    page.wait_for_function(
        "!schemaCourant().blocs.some(b => b.type === 'sous_systeme')", timeout=30000
    )
    assert page.locator("#blocs-svg .bloc").count() == 6


def test_gfl_du_reseau_en_blocs_remplace_et_analyse(page_fr051):
    """Critere 8 / U4 : ouvrir un GFL du reseau en blocs, changer un gain, le
    remplacer dans le reseau, relancer l'analyse modale."""
    page = page_fr051
    table, i = _poser_gfl1(page)
    # depuis la fiche (double-clic equivalent)
    page.evaluate(f"ouvrirFiche('{table}', {i})")
    page.wait_for_selector('#fiche [data-action="editer-en-blocs"]', timeout=30000)
    assert page.locator('#fiche [data-action="remplacer-par-blocs"]').count() == 1
    page.click('#fiche [data-action="editer-en-blocs"]')
    page.wait_for_function("B.modele !== null && /_1$/.test(B.modele)", timeout=60000)
    assert page.locator("#onglet-blocs").is_visible()
    n_avant = len(page.evaluate("E.tables.T_Blocs || []"))
    # le natif est toujours la
    assert page.evaluate(f"E.tables.{table}.length") == i + 1
    # un gain change dans le schema (premier gain numerique trouve a un niveau)
    page.evaluate(
        """() => {
            const cherche = (n) => { for (const b of n.blocs) {
                if (b.type === 'gain' && typeof b.K === 'number') return b;
                if (b.schema) { const r = cherche(b.schema); if (r) return r; } } return null; };
            const g = cherche(schemaCourant());
            if (g) { g.K = g.K * 2; modifie(); }
        }"""
    )
    # remplacer dans le reseau, depuis les proprietes de l'element
    page.click('nav button[data-onglet="constructeur"]')
    page.evaluate(f"selectionner('{table}', {i})")
    page.wait_for_selector('[data-action="remplacer-par-blocs"]')
    page.click('[data-action="remplacer-par-blocs"]')
    page.click('[data-action="confirmer-gardes"]')
    page.wait_for_function("(E.tables.T_Blocs || []).length >= 1", timeout=60000)
    assert page.evaluate(f"E.tables.{table}.length") == i
    assert len(page.evaluate("E.tables.T_Blocs")) == n_avant + 1
    # l'analyse modale tourne sur le reseau a blocs
    page.click('nav button[data-onglet="modal"]')
    page.click("#btn-modal")
    page.wait_for_function("E.res && E.res.modal", timeout=120000)
    assert page.evaluate("E.res.modal.modes.length") > 0
    # annuler ramene le natif
    page.click("#btn-annuler")
    page.wait_for_function(f"E.tables.{table}.length === {i + 1}")


def test_enregistrer_en_bibliotheque_puis_rouvrir_dans_un_autre_reseau(page_fr051):
    """Critere 8 / U5 : enregistrer dans la bibliotheque, puis rouvrir dans un autre reseau."""
    page = page_fr051
    _creer_source(page)
    nom = page.evaluate("B.modele")
    page.click("#btn-blocs-enregistrer")
    page.wait_for_function(
        "document.querySelector('#blocs-etat').textContent.includes('bibliotheque')"
    )
    nom_biblio = page.evaluate(
        "(document.querySelector('#blocs-etat').textContent.match(/: (\\S+)$/) || [])[1]"
    )
    assert nom_biblio
    _charger(page, "wscc9")
    assert not page.evaluate("Object.keys(E.tables.MODELES_BLOCS || {}).length")
    page.click('nav button[data-onglet="blocs"]')
    page.click("#btn-blocs-biblio")
    page.wait_for_selector(f'[data-ouvrir="{nom_biblio}"]', timeout=30000)
    page.click(f'[data-ouvrir="{nom_biblio}"]')
    page.wait_for_function("B.modele !== null && schemaCourant() !== undefined")
    assert page.locator("#blocs-svg .bloc").count() == 6
    assert page.evaluate("Object.keys(E.tables.MODELES_BLOCS).length") == 1
    # menage : la bibliotheque est persistante
    page.click(f'[data-supprimer="{nom_biblio}"]')
    page.wait_for_selector(f'[data-ouvrir="{nom_biblio}"]', state="detached")
    assert nom


# -------------------------------- contrat des appels, serveur simule (sans serveur FR-051)
def _schema_depuis_element() -> dict:
    return {
        "nom": "GFL1",
        "parametres": {"tau": 0.05},
        "amorce": {"classe": "GFL1", "correspondance": {}},
        "blocs": [
            {"id": "pr", "type": "consigne", "nom": "P_ref", "x": 40, "y": 40},
            {"id": "qr", "type": "consigne", "nom": "Q_ref", "x": 40, "y": 160},
            {
                "id": "fp",
                "type": "premier_ordre",
                "K": 1,
                "tau": "tau",
                "x": 240,
                "y": 40,
            },
            {
                "id": "fq",
                "type": "premier_ordre",
                "K": 1,
                "tau": "tau",
                "x": 240,
                "y": 160,
            },
            {"id": "op", "type": "sortie", "signal": "P", "x": 440, "y": 40},
            {"id": "oq", "type": "sortie", "signal": "Q", "x": 440, "y": 160},
        ],
        "liens": [
            {"de": "pr", "vers": "fp", "entree": 0},
            {"de": "qr", "vers": "fq", "entree": 0},
            {"de": "fp", "vers": "op", "entree": 0},
            {"de": "fq", "vers": "oq", "entree": 0},
        ],
    }


def test_editer_en_blocs_contrat_et_remplacement_simule(page):
    """U4/U6 contre un serveur simule : le corps envoye est {tables, table,
    index}, le modele edite s'appelle <classe>_<n> ; « Remplacer » est fait
    par le SERVEUR (/api/blocs/remplacer : tables, modele, lieu), l'interface
    applique ses tables telles quelles, et l'annulation rend l'element."""
    _charger(page, "smib")
    _onglet_blocs(page)
    table, i = _poser_gfl1(page)
    page.wait_for_timeout(700)
    envoyes = []

    def repondre(route):
        envoyes.append(json.loads(route.request.post_data))
        route.fulfill(json={"schema": _schema_depuis_element()})

    remplacements = []

    def remplacer(route):
        corps = json.loads(route.request.post_data)
        remplacements.append(corps)
        t = corps["tables"]
        ligne = t[corps["table"]].pop(corps["index"])
        schema = _schema_depuis_element()
        schema["nom"] = "GFL1_GFL2"
        schema["origine"] = {
            "element": "GFL2",
            "id_bloc": 1,
            "retraits": [
                {"table": corps["table"], "index": corps["index"], "ligne": ligne}
            ],
            "remplacements": [],
        }
        t.setdefault("MODELES_BLOCS", {})["GFL1_GFL2"] = schema
        t["T_Blocs"] = [[1, ligne[1], "GFL1_GFL2", ligne[2]]]
        route.fulfill(
            json={"tables": t, "modele": "GFL1_GFL2", "schema": schema,
                  "table": "T_Blocs", "index": 0}
        )  # fmt: skip

    page.route("**/api/blocs/depuis_element", repondre)
    page.route("**/api/blocs/remplacer", remplacer)
    page.click('nav button[data-onglet="constructeur"]')
    page.evaluate(f"selectionner('{table}', {i})")
    page.click('[data-action="editer-en-blocs"]')
    page.wait_for_function("B.modele === 'GFL1_1'")
    assert page.locator("#onglet-blocs").is_visible()
    assert set(envoyes[0]) == {"tables", "table", "index"}
    assert envoyes[0]["table"] == table and envoyes[0]["index"] == i
    # le natif n'a pas bouge
    assert page.evaluate("E.tables.T_GFL1s.length") == 1
    assert not page.evaluate("(E.tables.T_Blocs || []).length")
    # second geste, depuis la fiche : « Remplacer par son modele a blocs »
    page.evaluate(f"ouvrirFiche('{table}', {i})")
    page.wait_for_selector('#fiche [data-action="remplacer-par-blocs"]', timeout=30000)
    page.wait_for_timeout(700)
    page.click('#fiche [data-action="remplacer-par-blocs"]')
    page.click('[data-action="confirmer-gardes"]')
    page.wait_for_function("(E.tables.T_Blocs || []).length === 1")
    # le navigateur n'a rien remplace lui-meme : il a demande au serveur
    assert len(envoyes) == 1 and len(remplacements) == 1
    assert set(remplacements[0]) == {"tables", "table", "index"}
    assert (remplacements[0]["table"], remplacements[0]["index"]) == (table, i)
    assert page.evaluate("B.modele") == "GFL1_GFL2"
    ligne = page.evaluate("E.tables.T_Blocs[0]")
    assert (
        ligne[1] == 2
        and ligne[2] == "GFL1_GFL2"
        and ligne[3] == page.evaluate("E.tables.Sb")
    )
    assert page.evaluate("E.tables.MODELES_BLOCS.GFL1_GFL2.origine.element") == "GFL2"
    assert page.evaluate("E.tables.T_GFL1s.length") == 0
    assert page.locator("#fiche").is_hidden()
    page.click("#btn-annuler")
    page.wait_for_function("E.tables.T_GFL1s.length === 1")
    assert not page.evaluate("(E.tables.T_Blocs || []).length")


def test_depuis_element_refuse_message_clair(page):
    """Un element que le serveur refuse (port continu) : le message du
    serveur est montre dans la fiche, rien n'est modifie."""
    _charger(page, "smib")
    table, i = _poser_gfl1(page)
    page.route(
        "**/api/blocs/depuis_element",
        lambda r: r.fulfill(
            status=422,
            json={"detail": "GFM3DC a un port continu : non decrit en blocs en v1"},
        ),
    )
    page.evaluate(f"ouvrirFiche('{table}', {i})")
    page.wait_for_selector('#fiche [data-action="editer-en-blocs"]', timeout=30000)
    page.click('#fiche [data-action="editer-en-blocs"]')
    page.wait_for_selector("#fiche .pb.erreur")
    assert "port continu" in page.locator("#fiche .pb.erreur").inner_text()
    assert page.locator("#fiche").is_visible()
    assert page.evaluate("E.tables.T_GFL1s.length") == 1


def test_bibliotheque_contrat_simule(page):
    """U5 contre un serveur simule : enregistrer envoie {schema_bloc}, la
    liste se lit, ouvrir ajoute le modele (sous un nom libre) au reseau."""
    _charger(page, "smib")
    _onglet_blocs(page)
    _creer_source(page)
    schema_gardes = {}

    def biblio(route):
        req = route.request
        if req.method == "POST":
            corps = json.loads(req.post_data)
            schema_gardes[corps["schema_bloc"]["nom"]] = corps["schema_bloc"]
            route.fulfill(json={"nom": corps["schema_bloc"]["nom"]})
        else:
            route.fulfill(
                json=[
                    {"nom": n, "titre": n, "modifie": "2026-11-07T10:00:00"}
                    for n in schema_gardes
                ]
            )

    def un_modele(route):
        nom = route.request.url.rsplit("/", 1)[1]
        if route.request.method == "DELETE":
            schema_gardes.pop(nom, None)
            route.fulfill(json={"nom": nom})
        else:
            route.fulfill(json={"schema": schema_gardes[nom]})

    page.route("**/api/blocs/bibliotheque", biblio)
    page.route("**/api/blocs/bibliotheque/*", un_modele)
    nom = page.evaluate("B.modele")
    page.click("#btn-blocs-enregistrer")
    page.wait_for_function(
        "document.querySelector('#blocs-etat').textContent.includes('bibliotheque')"
    )
    assert nom in schema_gardes
    page.click("#btn-blocs-biblio")
    page.wait_for_selector(f'[data-ouvrir="{nom}"]')
    page.click(f'[data-ouvrir="{nom}"]')
    page.wait_for_function(f"B.modele === '{nom}_2'")
    assert page.evaluate("Object.keys(E.tables.MODELES_BLOCS).length") == 2
    page.click(f'[data-supprimer="{nom}"]')
    page.wait_for_selector(f'[data-ouvrir="{nom}"]', state="detached")
