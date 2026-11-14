"""L'atelier dans un vrai navigateur : les gestes de l'utilisateur, de bout
en bout, contre un serveur lance pour l'occasion -- jeton compris.

Ces tests demandent Playwright et un Chromium ::

    uv sync --all-extras --group e2e
    uv run playwright install chromium        # sauf s'il est deja la
    uv run pytest tests/e2e

Sans Playwright, ils sont sautes ; ``G2EGRID_CHROMIUM`` designe un Chromium
deja installe si Playwright ne trouve pas le sien.
"""

from __future__ import annotations

import glob
import json
import os
import pathlib
import socket
import subprocess
import sys
import time
import urllib.request

import pytest

pytest.importorskip("fastapi")
sync_api = pytest.importorskip("playwright.sync_api")

ICI = pathlib.Path(__file__).resolve().parents[2]
DONNEES = ICI / "tests" / "donnees" / "interop"
JETON = "jeton-e2e-0123456789abcdef"


def _port_libre() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture(scope="module")
def serveur(tmp_path_factory):
    port = _port_libre()
    journal = open(tmp_path_factory.mktemp("serveur") / "uvicorn.log", "w")
    env = dict(os.environ, G2EGRID_JETON=JETON)
    p = subprocess.Popen(
        [sys.executable, "-m", "uvicorn", "g2egrid_api.app:app", "--port", str(port),
         "--log-level", "warning"],
        cwd=ICI, env=env, stdout=journal, stderr=journal,
    )  # fmt: skip
    base = f"http://127.0.0.1:{port}"
    fin = time.time() + 60
    while time.time() < fin:
        try:
            urllib.request.urlopen(base + "/api/defauts", timeout=1)
        except urllib.error.HTTPError:
            break  # 403 : il repond, et il garde la porte
        except OSError:
            time.sleep(0.3)
    yield base
    p.terminate()
    p.wait(10)
    journal.close()


@pytest.fixture(scope="module")
def navigateur():
    with sync_api.sync_playwright() as p:
        try:
            b = p.chromium.launch()
        except Exception:  # noqa: BLE001
            chemins = [os.environ.get("G2EGRID_CHROMIUM", "")] + glob.glob(
                "/opt/pw-browsers/chromium-*/chrome-linux/chrome"
            )
            chemin = next((c for c in chemins if c and os.path.exists(c)), None)
            if chemin is None:
                pytest.skip("pas de Chromium pour Playwright")
            b = p.chromium.launch(executable_path=chemin)
        yield b
        b.close()


@pytest.fixture
def page(serveur, navigateur):
    # l'interface prend la langue du navigateur : ces tests lisent l'atelier
    # en francais (tests/e2e/test_langue.py eprouve l'anglais)
    ctx = navigateur.new_context(
        viewport={"width": 1500, "height": 1100}, accept_downloads=True, locale="fr-FR"
    )
    pg = ctx.new_page()
    erreurs: list[str] = []
    pg.on("pageerror", lambda e: erreurs.append(str(e)))
    pg.goto(f"{serveur}/?jeton={JETON}")
    pg.wait_for_function("typeof E !== 'undefined' && E.defauts", timeout=30000)
    pg.evaluate("localStorage.clear()")
    yield pg
    assert not erreurs, erreurs
    ctx.close()


def _charger(pg, nom, attente=None):
    pg.select_option("#catalogue", nom)
    pg.click("#btn-charger")
    pg.wait_for_function(
        attente
        or f"E.tables && (E.tables.nom || '').toLowerCase() === '{nom}'.toLowerCase()",
        timeout=30000,
    )


#: le CIGRE B4 ne porte pas son nom de catalogue : on attend ses onze stations
B4 = "E.tables && (E.tables.T_VSCs || []).length === 11"


def test_sans_jeton_la_porte_est_fermee(serveur, navigateur):
    ctx = navigateur.new_context(locale="fr-FR")
    pg = ctx.new_page()
    r = pg.goto(serveur + "/")
    assert r.status == 403 and "jeton" in pg.content()
    ctx.close()


def test_nouveau_pendant_le_demarrage_n_est_pas_ecrase(serveur, navigateur):
    """Le reseau de demarrage arrive APRES un clic sur Nouveau : il ne doit
    pas l'ecraser. On retient sa reponse le temps de cliquer."""
    ctx = navigateur.new_context(locale="fr-FR")
    pg = ctx.new_page()
    erreurs: list[str] = []
    pg.on("pageerror", lambda e: erreurs.append(str(e)))
    retenues = []
    pg.route("**/api/reseaux/*", lambda route: retenues.append(route))
    pg.goto(f"{serveur}/?jeton={JETON}")
    pg.wait_for_function("typeof E !== 'undefined' && E.defauts", timeout=30000)
    for _ in range(100):
        if retenues:
            break
        pg.wait_for_timeout(100)
    assert retenues, "le chargement du demarrage n'est pas parti"
    pg.click("#btn-nouveau")
    pg.wait_for_function("E.tables && E.tables.nom === 'nouveau'")
    for r in retenues:
        r.continue_()
    pg.wait_for_function(
        "!document.querySelector('#etat').textContent.startsWith('chargement')"
    )
    assert pg.evaluate("[E.tables.nom, E.tables.Y_DER.length]") == ["nouveau", 0]
    assert not erreurs, erreurs
    ctx.close()


def test_annuler_retablir(page):
    _charger(page, "wscc9")
    page.evaluate("E.tables.Y_DER.pop(); modifie()")
    assert page.evaluate("E.tables.Y_DER.length") == 8
    page.keyboard.press("Control+z")
    assert page.evaluate("E.tables.Y_DER.length") == 9
    page.keyboard.press("Control+y")
    assert page.evaluate("E.tables.Y_DER.length") == 8


def test_variations_liste_les_parametres(page):
    """L'onglet Variations propose les parametres de l'element choisi.

    Regression : une seconde ``remplirParametres`` (celle de la continuation)
    ecrasait la premiere, et la liste des parametres restait vide."""
    _charger(page, "wscc9")
    page.click("button[data-onglet=balayage]")
    page.wait_for_function(
        "document.querySelectorAll('#bal-element option').length > 0"
    )
    assert page.eval_on_selector_all("#bal-param option", "os => os.length") > 0
    page.select_option("#bal-element", index=1)
    noms = page.eval_on_selector_all(
        "#bal-param option", "os => os.map(o => o.textContent)"
    )
    assert noms, "aucun parametre pour le deuxieme element"
    assert page.inner_text("#bal-actuelle").strip()


def test_un_bouton_par_identifiant():
    """Aucun identifiant en double dans la page : ``$("#x")`` ne relie que le
    premier, et le second bouton restait inerte (« Mesurer l'impact »)."""
    import re
    from collections import Counter

    html = (ICI / "g2egrid_api" / "web" / "index.html").read_text(encoding="utf-8")
    doubles = [
        k for k, n in Counter(re.findall(r'id="([^"]+)"', html)).items() if n > 1
    ]
    assert not doubles, doubles


def test_reparer(page):
    _charger(page, "wscc9")
    page.evaluate(
        "E.tables.Y_DER = E.tables.Y_DER.map(r => r[1] === 1 ? [r[0], 2, ...r.slice(2)] : r); modifie()"
    )
    page.wait_for_selector("button.reparer", timeout=15000)
    page.locator("button.reparer").first.click()
    page.wait_for_function(
        "document.querySelector('#etat').textContent.startsWith('repare')",
        timeout=15000,
    )
    assert page.evaluate("E.tables.Y_DER.filter(r => r[1] === 1).length") == 1


def test_mode_sur_le_schema(page):
    _charger(page, "wscc9")
    page.click('nav button[data-onglet="modal"]')
    page.click("#btn-modal")
    page.wait_for_function("E.res.modal", timeout=60000)
    k = page.evaluate(
        "(() => { const m = E.res.modal.modes.filter(x => x.f > 0.3 && x.f < 3);"
        " m.sort((a, b) => a.zeta - b.zeta); return m[0].indice; })()"
    )
    page.evaluate(f"ouvrirMode({k})")
    page.wait_for_function(
        "!document.querySelector('#rangee-mode-schema').hidden", timeout=60000
    )
    assert page.locator("#schema-mode .phaseurs line").count() >= 2


def test_ouvrir_un_cas_matpower(page):
    page.set_input_files("#fichier", str(DONNEES / "case9.m"))
    page.wait_for_function(
        "document.querySelector('#etat').textContent.includes('case9.m')", timeout=30000
    )
    assert page.evaluate("E.tables.Y_DER.length") == 9


def test_rapport_puis_reouverture(page, tmp_path):
    _charger(page, "smib")
    with page.expect_download(timeout=120000) as dl:
        page.click("#btn-rapport")
    chemin = tmp_path / "rapport.html"
    dl.value.save_as(chemin)
    texte = chemin.read_text()
    d = json.loads(texte.split('id="g2egrid-tables">', 1)[1].split("</script>", 1)[0])
    assert d["empreinte"] in texte
    _charger(page, "wscc9")
    page.set_input_files("#fichier", str(chemin))
    page.wait_for_function(
        "document.querySelector('#etat').textContent.includes('rapport rouvert')",
        timeout=30000,
    )
    assert page.evaluate("E.tables.Y_DER.length") == len(d["tables"]["Y_DER"])


def test_frise_et_rejeu(page):
    _charger(page, "wscc9")
    page.click('nav button[data-onglet="temporel"]')
    page.wait_for_function("E.res.signaux", timeout=60000)
    page.fill("#duree", "1")
    page.click("#bloc-scenario summary")
    page.click("#btn-ajouter-evt")
    page.evaluate(
        "E.scenario[0] = {...E.scenario[0], genre: 'charge', cible: 'Load5', instant: 0.2,"
        " dp: 0.2, dq: 0.1}; rendreScenario()"
    )
    assert page.locator("#frise polygon").count() == 1
    page.click("#btn-temporel")
    page.wait_for_function("E.res.temporel", timeout=300000)
    page.wait_for_selector("#rangee-rejeu:not([hidden])")
    page.fill("#rejeu-duree", "1")
    page.click("#btn-rejeu")
    page.wait_for_function(
        "document.querySelector('#btn-rejeu').textContent === 'Lecture'", timeout=20000
    )
    assert "charge" in page.text_content("#rejeu-evts")


# l'outil de la palette, la source attendue sur le symbole et la commande
# que l'editeur accorde au convertisseur
SOURCES = [
    ("gfm_batterie", "batterie", "PV"),
    ("gfm_pv", "PV", "VdcV"),
    ("gfm_dc", "DC", "VdcV"),
    ("gfl_pv", "PV", "VdcQ"),
    ("gfl_batterie", "batterie", "PQ"),
    ("gfl_dc", "DC", "VdcQ"),
]


@pytest.mark.parametrize("outil,source,mode", SOURCES)
def test_palette_convertisseur_avec_sa_source(page, outil, source, mode):
    """Un clic sur l'outil, un clic sur un bus : le convertisseur, son
    transformateur, son bus continu, son condensateur et sa source sont poses,
    le symbole porte le pictogramme de la source, et le reseau obtenu a un
    point de fonctionnement stable."""
    from g2egrid.robustesse import _spectre
    from g2egrid.tables import tables_depuis_json

    _charger(page, "smib")
    assert page.locator("#palette-outils .sep").count() >= 3  # une seule palette
    page.click(f'button[data-outil="{outil}"]')
    page.locator('#canevas circle[data-bus="2"]').first.click()
    t = page.evaluate("E.tables")
    assert len(t["T_VSCs"]) == 1 and t["T_VSCs"][0][3] == mode
    assert len(t["T_ShuntsDC"]) == 1
    n_dcdc = 0 if source == "DC" else 1
    assert len(t["T_DCDC"]) == n_dcdc
    if n_dcdc:
        assert t["T_DCDC"][0][3].lower() == source.lower()
    # la source est dessinee du cote continu : en entier (batterie, PV), ou
    # c'est le bus continu, relie a la face continue du convertisseur
    if source == "DC":
        assert page.locator('#canevas [data-sel-table="T_VSCs"]').count() >= 1
    else:
        assert (
            page.locator(f'#canevas .symbole-source[data-source="{source}"]').count()
            >= 1
        )
    if source != "DC":  # un convertisseur seul sur son bus DC n'a pas de source
        assert float(_spectre(tables_depuis_json(t)).lam.real.max()) < 1e-3


def test_palette_gfm_batterie_premier_formeur(page):
    """Dans une ile sans reference d'angle, le GFM batterie se pose a bus
    continu ideal : il peut alors etre le slack."""
    page.click("#btn-nouveau")
    page.wait_for_function("E.tables && E.tables.Y_DER.length === 0")
    page.click('button[data-outil="bus"]')
    page.locator("#canevas").click(position={"x": 300, "y": 200})
    page.wait_for_function("E.tables.Y_DER.length === 1")
    page.click('button[data-outil="gfm_batterie"]')
    page.locator("#canevas circle[data-bus]").first.click()
    page.wait_for_function("E.tables.T_GFM1s.length + E.tables.T_GFM3s.length === 1")
    t = page.evaluate("E.tables")
    assert len(t["T_GFM1s"]) == 1 and not t["T_VSCs"]
    assert any(r[1] == 1 for r in t["Y_DER"])
    assert page.locator('#canevas .symbole-source[data-source="batterie"]').count() >= 1


def test_statisme_dans_l_inspecteur(page):
    """Une station en statisme montre k, P_ref, Vdc_ref et sa caracteristique ;
    apres un load flow, le point de fonctionnement y est marque."""
    _charger(page, "cigre_b4", B4)
    page.evaluate("calculerLoadflow()")
    page.wait_for_function("E.res.lf && E.res.lf.bus_dc.length === 13", timeout=60000)
    i = page.evaluate("E.tables.T_VSCs.findIndex(r => r[1] === 111)")
    page.evaluate(f"selectionner('T_VSCs', {i})")
    bloc = page.locator("#inspecteur .bloc-statisme, .bloc-statisme").first
    assert bloc.is_visible()
    k = bloc.locator("input").first.input_value()
    assert abs(float(k) - 19.9992) < 1e-3
    assert page.locator(".carac-statisme circle[stroke]").count() == 1


def test_passer_une_station_en_statisme(page):
    """Choisir ``droop`` met les tables en coherence (noeud continu de type 2,
    statisme par defaut) et le reseau garde un point de fonctionnement."""
    from g2egrid.tables import assembler, tables_depuis_json

    _charger(page, "mtdc3t")
    page.evaluate("selectionner('T_VSCs', 2)")  # GFL du bus 5, en PQ
    page.locator("select").filter(
        has=page.locator('option[value="droop"]')
    ).first.select_option("droop")
    t = page.evaluate("E.tables")
    assert t["T_VSCs"][2][3] == "droop"
    assert next(r for r in t["T_BusDC"] if r[0] == 12)[1] == 2
    assert t["PARAMS_PAR_BUS"]["5"]["kdroop"] == 0.05
    assert page.locator(".bloc-statisme").count() == 1
    lf = assembler(tables_depuis_json(t)).run_loadflow()
    assert lf.converged


def test_balayage_des_statismes(page):
    """Le balayage, de l'onglet Limites et impedances : partage de puissance,
    tensions continues et marge, un point par facteur."""
    _charger(page, "cigre_b4", B4)
    page.click('nav button[data-onglet="stabilite"]')
    page.fill("#stat-pas", "3")
    page.click("#btn-stat")
    page.wait_for_function("E.res.statismes", timeout=180000)
    d = page.evaluate("E.res.statismes")
    assert d["stations"] == [111, 113, 114] and len(d["mu"]) == 3
    # plus raide (mu < 1), la station de Cb-B2 s'ecarte davantage de P_ref
    p113 = d["p_mw"]["113"]
    assert abs(p113[0] - 1700) > abs(p113[-1] - 1700)
    assert page.locator("#carte-statp path, #carte-statp polyline").count() > 0


def test_modes_par_domaine(page):
    """Onglet Modes d'un reseau AC/DC : filtre par domaine, colonne du
    domaine, et comparaison au continu rigide."""
    _charger(page, "cigre_b4", B4)
    page.click('nav button[data-onglet="modal"]')
    page.click("#btn-modal")
    page.wait_for_function("E.res.modal", timeout=120000)
    assert page.locator("#bloc-filtre-domaine").is_visible()
    tous = page.locator("#table-modes tbody tr").count()
    page.select_option("#filtre-domaine", "DC")
    dc = page.locator("#table-modes tbody tr").count()
    assert 0 < dc < tous
    assert page.locator("#carte-modes rect[fill]").count() >= 1  # carres : DC
    page.click("#btn-impact-dc")
    page.wait_for_function("E.res.impact_dc", timeout=120000)
    assert page.locator("#impact-resultat table").count() >= 2


def test_defaut_dc_et_protections(page):
    """Onglet Temporel de DCS1 : un defaut continu dans la frise, les
    protections DC cochees ; les blocages sont listes sous les courbes."""
    _charger(page, "cigre_b4_dcs1", "E.tables && (E.tables.T_VSCs || []).length === 2")
    page.click('nav button[data-onglet="temporel"]')
    page.wait_for_function("E.res.signaux", timeout=60000)
    assert page.evaluate("E.res.signaux.manoeuvres.bus_dc") == [201, 221]
    page.fill("#duree", "0.05")
    page.click("#bloc-scenario summary")
    page.click("#btn-ajouter-evt")
    page.evaluate(
        "E.scenario[0] = {...E.scenario[0], genre: 'defaut_dc', cible: 'bus221', instant: 0.02};"
        " rendreScenario()"
    )
    page.check("#scenario-actif")
    page.check("#protections-dc")
    page.fill("#pas", "25")
    page.click("#btn-temporel")
    page.wait_for_function("E.res.temporel", timeout=300000)
    assert page.evaluate("E.res.temporel.declenchements.length") >= 2
    assert "sous tension dc bus221" in page.text_content(".declenchements").replace(
        "_", " "
    )


def test_station_mmc_et_cable_par_pole(page):
    """L'inspecteur d'une station a port continu propose son modele (VSC ou
    MMC) ; le choix s'ecrit dans PARAMS_PAR_BUS, les parametres affiches
    deviennent ceux du MMC, et le reseau garde un point de fonctionnement.
    Celui d'un cable continu propose sa configuration et ses circuits."""
    _charger(page, "cigre_b4_dcs1", "E.tables && (E.tables.T_VSCs || []).length === 2")
    i = page.evaluate("E.tables.T_GFL2s.findIndex(r => r[1] === 121)")
    page.evaluate(f"selectionner('T_GFL2s', {i})")
    page.select_option(".bloc-modele-station select", "mmc")
    assert page.evaluate("E.tables.PARAMS_PAR_BUS['121'].modele") == "mmc"
    assert page.evaluate("defautsReglages('GFL', 121).W0") == 0.03
    assert "MMC" in page.text_content("#inspecteur")
    page.evaluate("calculerLoadflow()")
    page.wait_for_function("E.res.lf && E.res.lf.bus_dc.length === 2", timeout=60000)
    page.evaluate("selectionner('T_LinesDC', 0)")
    assert page.locator("#inspecteur select").last.input_value() == "monopole_sym"


@pytest.mark.parametrize(
    "reseau, attente, table, schema",
    [
        (
            "cigre_b4_dcs1",
            "E.tables && (E.tables.T_VSCs || []).length === 2",
            "T_GFL2s",
            "gfl_vdc_port",
        ),
        ("wscc9_gfm", None, "T_GFM1s", "gfm"),
        (
            "cigre_b4_mmc",
            "E.tables && (E.tables.T_VSCs || []).length === 11",
            "T_GFL2s",
            "mmc",
        ),
    ],
)
def test_fiche_au_double_clic(page, reseau, attente, table, schema, tmp_path):
    """Double-clic sur un convertisseur : la fiche s'ouvre, avec ses deux
    schemas (electrique, commande) et ses equations en MathML, rangees par
    sous-systeme ; Echap la ferme."""
    _charger(page, reseau, attente)
    if not page.evaluate(f"(E.tables.{table} || []).length"):
        pytest.skip(f"{reseau} n'a pas de {table}")
    # la boite du convertisseur (son groupe porte aussi sa source)
    sym = page.locator(f'#canevas [data-sel-table="{table}"] rect').first
    sym.dblclick(force=True)
    page.wait_for_function("E.fiche && !$('#fiche').hidden", timeout=60000)
    assert page.evaluate("E.fiche.schema") == schema
    # un schema-bloc par boucle : sommateurs, blocs de transfert, retours
    assert page.locator("#fiche .schema-bloc").count() >= 4
    assert page.locator("#fiche .schema-bloc circle.somme").count() >= 6
    titres = page.locator("#fiche [data-schema-bloc]").evaluate_all(
        "l => l.map(e => e.dataset.schemaBloc)"
    )
    assert any("Boucle de courant" in t for t in titres), titres
    # les equations restent, repliees
    assert page.locator("#fiche .fiche-equations math").count() >= 15
    page.set_viewport_size({"width": 1500, "height": 5200})
    page.locator("#fiche .fiche-page").screenshot(
        path=str(tmp_path / f"fiche_{table}.png")
    )
    page.keyboard.press("Escape")
    assert page.locator("#fiche").is_hidden()


@pytest.mark.parametrize(
    "outil,source", [("cote_pv", "PV"), ("cote_batterie", "batterie")]
)
def test_source_du_cote_continu(page, outil, source):
    """Batterie ou PV pose en cliquant un GFL a bus continu ideal : le GFL
    passe dans sa variante a port continu, recoit un bus continu, son
    condensateur et un DC/DC de la source, dessinee du cote continu."""
    _charger(page, "cigre_mv_dc", "E.tables && E.tables.Y_DER.length === 17")
    i = page.evaluate("E.tables.T_GFL1s.findIndex(r => r[1] === 14)")
    assert i >= 0
    page.click(f'button[data-outil="{outil}"]')
    page.locator('#canevas [data-sel-table="T_GFL1s"][data-bus="14"] rect').first.click(
        force=True
    )
    t = page.evaluate("E.tables")
    assert not any(r[1] == 14 for r in t["T_GFL1s"])
    assert any(r[1] == 14 for r in t["T_GFL2s"])
    vsc = next(r for r in t["T_VSCs"] if r[1] == 14)
    assert any(r[1] == vsc[2] and r[3] == source for r in t["T_DCDC"])
    assert any(r[1] == vsc[2] for r in t["T_ShuntsDC"])
    assert (
        page.locator(f'#canevas .symbole-source[data-source="{source}"]').count() >= 1
    )


def test_commande_du_convertisseur(page):
    """Le bloc Commande ecrit les tables : statisme Idc-Vdc (mode droopI) sur
    une station, statisme Q-V sur une autre, loi VSM sur un GFM ; le load flow
    converge."""
    _charger(page, "cigre_b4", B4)
    i = page.evaluate("E.tables.T_GFL2s.findIndex(r => r[1] === 111)")
    page.evaluate(f"selectionner('T_GFL2s', {i})")
    page.locator(".bloc-commande select").nth(0).select_option("droopI")
    assert page.evaluate("E.tables.T_VSCs.find(r => r[1] === 111)[3]") == "droopI"
    i = page.evaluate("E.tables.T_GFL2s.findIndex(r => r[1] === 122)")
    page.evaluate(f"selectionner('T_GFL2s', {i})")
    page.locator(".bloc-commande select").nth(1).select_option("QV")
    assert page.evaluate("E.tables.PARAMS_PAR_BUS[122].statisme_qv") == 0.05
    page.evaluate("calculerLoadflow()")
    page.wait_for_function("E.res && E.res.lf", timeout=60000)
    i = page.evaluate("E.tables.T_GFM3s.findIndex(r => r[1] === 121)")
    page.evaluate(f"selectionner('T_GFM3s', {i})")
    page.locator(".bloc-commande select").last.select_option("vsm")
    assert page.evaluate("E.tables.T_GFM4s.some(r => r[1] === 121)")


def test_editeur_de_blocs(page):
    """Onglet Blocs : un modele cree depuis un exemple, cable a la souris,
    compile en equations, puis pose sur un bus par la palette du reseau ; le
    reseau obtenu a un point de fonctionnement stable, ou l'on retrouve le
    pole du premier ordre du schema."""
    import numpy as np

    from g2egrid.robustesse import _spectre
    from g2egrid.tables import tables_depuis_json

    _charger(page, "smib")
    page.click('nav button[data-onglet="blocs"]')
    page.select_option("#blocs-exemple", "source_filtree")
    page.click("#btn-blocs-nouveau")
    assert page.locator("#blocs-svg .bloc").count() == 6
    # un gain insere entre le filtre de P et sa sortie, a la souris
    page.click('.blocs-palette button[data-bloc="gain"]')
    page.evaluate(
        "schemaCourant().liens = schemaCourant().liens.filter(l => l.vers !== 'op')"
    )
    page.evaluate("dessinerSchemaBlocs()")

    def relier(de, vers):
        s = page.locator(
            f'#blocs-svg .bloc-port-sortie[data-bloc="{de}"]'
        ).bounding_box()
        t = page.locator(
            f'#blocs-svg .bloc-port-entree[data-bloc="{vers}"]'
        ).bounding_box()
        page.mouse.move(s["x"] + 4, s["y"] + 4)
        page.mouse.down()
        page.mouse.move(t["x"] + 4, t["y"] + 4, steps=6)
        page.mouse.up()

    relier("fp", "g1")
    relier("g1", "op")
    liens = page.evaluate("schemaCourant().liens")
    assert {"de": "fp", "vers": "g1", "entree": 0} in liens
    assert {"de": "g1", "vers": "op", "entree": 0} in liens
    page.click("#btn-blocs-compiler")
    page.wait_for_function("B.compile !== null", timeout=30000)
    eq = page.locator("#blocs-equations").inner_text()
    assert "dx_fp/dt" in eq and "y_g1 = y_fp" in eq
    # pose sur le reseau
    page.click("#btn-blocs-poser")
    page.locator('#canevas circle[data-bus="2"]').first.click()
    t = page.evaluate("E.tables")
    assert len(t["T_Blocs"]) == 1
    modele = t["T_Blocs"][0][2]
    assert modele in t["MODELES_BLOCS"]
    b = int(t["T_Blocs"][0][1])
    t["Y_DER"] = [
        [*r[:5], 0.1 * t["Sb"], 0.0, *r[7:]] if int(r[0]) == b else r
        for r in t["Y_DER"]
    ]
    tau = t["MODELES_BLOCS"][modele]["parametres"]["tau"]
    lam = _spectre(tables_depuis_json(t)).lam
    assert float(lam.real.max()) < 1e-3
    assert np.min(np.abs(lam + 1 / tau)) < 1e-6 / tau


def test_gfl_complet_en_blocs(page):
    """L'exemple « GFL complet » (LCL, PLL, boucles de courant et de
    puissance, limitation de courant : 32 blocs) se cree, se compile et se pose ; le reseau obtenu a
    les memes valeurs propres qu'avec le GFL1 du moteur au meme endroit."""
    import numpy as np

    from g2egrid.robustesse import _spectre
    from g2egrid.tables import tables_depuis_json

    _charger(page, "smib")
    page.click('nav button[data-onglet="blocs"]')
    page.select_option("#blocs-exemple", "gfl_complet")
    page.click("#btn-blocs-nouveau")
    assert page.locator("#blocs-svg .bloc").count() == 33  # limiteurs compris
    page.click("#btn-blocs-compiler")
    page.wait_for_function("B.compile !== null", timeout=30000)
    assert page.evaluate("B.compile.etats.length") == 12
    page.click("#btn-blocs-poser")
    page.locator('#canevas circle[data-bus="2"]').first.click()
    t = page.evaluate("E.tables")
    lid, b, modele, Sn = t["T_Blocs"][0]
    t["Y_DER"] = [
        [*r[:5], 0.3 * Sn, 0.0, *r[7:]] if int(r[0]) == b else r for r in t["Y_DER"]
    ]
    lam_blocs = _spectre(tables_depuis_json(t)).lam
    t["T_Blocs"] = []
    t["T_GFL1s"] = [[1, b, Sn, 40e3]]
    lam_gfl1 = _spectre(tables_depuis_json(t)).lam

    def tri(lam):
        return np.sort_complex(lam[np.abs(lam) > 1e-6])

    a, c = tri(lam_gfl1), tri(lam_blocs)
    assert a.size == c.size
    assert np.max(np.abs(a - c) / np.maximum(1.0, np.abs(a))) < 1e-6


def test_puissance_de_court_circuit(page):
    """Onglet Point de fonctionnement : Scc et SCR de chaque noeud, dans la
    table et sur le schema (etiquettes et couleur par SCR)."""
    _charger(page, "cigre_mv")
    page.click('nav button[data-onglet="loadflow"]')
    page.click("#btn-scc")
    page.wait_for_function("E.res.scc", timeout=60000)
    n_bus = page.evaluate("E.tables.Y_DER.length")
    assert page.locator("#table-scc tbody tr").count() == n_bus
    assert "SCRIF min" in page.locator("#scc-etat").inner_text()
    textes = page.evaluate(
        "[...document.querySelectorAll('#canevas-scc text')].map(t => t.textContent)"
    )
    assert any(t.startswith("SCR ") for t in textes)
    assert any(t.endswith(" MVA") for t in textes)


def test_gfm_complet_en_blocs(page):
    """L'exemple « GFM complet » (VSM, statisme Q-V, boucles de tension et de
    courant) se cree, se compile et se pose en formeur (bus PV) ; le reseau a
    les memes valeurs propres qu'avec le GFM2 du moteur au meme endroit."""
    import numpy as np

    from g2egrid.robustesse import _spectre
    from g2egrid.tables import tables_depuis_json

    _charger(page, "smib")
    page.click('nav button[data-onglet="blocs"]')
    page.select_option("#blocs-exemple", "gfm_complet")
    page.click("#btn-blocs-nouveau")
    page.click("#btn-blocs-compiler")
    page.wait_for_function("B.compile !== null", timeout=30000)
    assert page.evaluate("B.compile.etats.length") == 14
    page.click("#btn-blocs-poser")
    page.locator('#canevas circle[data-bus="2"]').first.click()
    t = page.evaluate("E.tables")
    lid, b, modele, Sn = t["T_Blocs"][0]
    ligne = next(r for r in t["Y_DER"] if int(r[0]) == b)
    assert int(ligne[1]) in (1, 2)  # un formeur : slack ou PV, jamais PQ
    t["Y_DER"] = [
        [*r[:5], 0.3 * Sn, 0.0, *r[7:]] if int(r[0]) == b else r for r in t["Y_DER"]
    ]
    lam_blocs = _spectre(tables_depuis_json(t)).lam
    t["T_Blocs"] = []
    t["T_GFM2s"] = [[1, b, Sn, 40e3]]
    lam_gfm2 = _spectre(tables_depuis_json(t)).lam

    def tri(lam):
        return np.sort_complex(lam[np.abs(lam) > 1e-6])

    a, c = tri(lam_gfm2), tri(lam_blocs)
    assert a.size == c.size
    assert np.max(np.abs(a - c) / np.maximum(1.0, np.abs(a))) < 1e-6


def test_grandeurs_apres_un_evenement(page):
    """Onglet Evenements : un echelon de la tension amont fait passer le GFL
    au-dessus de 1,10 pu ; son relais (0,3 s) le fait declencher, et les
    tableaux des bus et des unites le disent."""
    _charger(page, "smib_gfl")
    page.click('nav button[data-onglet="evenements"]')
    page.wait_for_function("E.res.signaux", timeout=60000)
    page.select_option("#evt-genre", "echelon")
    page.fill("#evt-filtre", "E_ref")
    page.select_option("#evt-cibles", "IB1.E_ref")
    page.fill("#evt-instant", "0.1")
    page.fill("#evt-valeur", "0.12")
    page.click("#btn-evt-ajouter")
    assert page.locator("#evt-table tbody tr").count() == 1
    page.fill("#evt-retard", "0.3")
    page.fill("#evt-duree", "0.8")
    page.fill("#evt-pas", "0.1")
    page.fill("#evt-n", "161")
    page.click("#btn-evt-simuler")
    page.wait_for_function("EVT.res", timeout=300000)
    assert page.evaluate("EVT.res.declenchements.map(d => d.protection)") == [
        "surtension_ac:GFL2"
    ]
    page.wait_for_selector("#evt-rangee-v:not([hidden])")
    assert page.locator("#evt-table-bus tbody tr.hors-bornes").count() >= 1
    gfl = page.locator("#evt-table-unites tbody tr.sortie")
    assert gfl.count() == 1 and "GFL2" in gfl.text_content()
    assert page.locator("#evt-courbes-v polyline").count() >= 2
    assert "GFL2" in page.text_content("#evt-relais-liste")


def test_cascade_quasi_statique(page):
    """La cascade du cas d'etude, une a une : trois GFL ouverts au depart,
    quatre de plus declenchent (etudes/surtension_cascade)."""
    _charger(page, "wscc3sm_cigre3_surtension")
    page.click('nav button[data-onglet="evenements"]')
    page.wait_for_function("E.res.signaux", timeout=120000)
    page.fill("#evt-filtre", "GFL30")
    page.select_option("#evt-cibles", ["GFL303", "GFL304", "GFL305"])
    page.click("#btn-evt-ajouter")
    assert page.locator("#evt-table tbody tr").count() == 3
    page.select_option("#evt-ordre", "une_a_une")
    page.click("#btn-evt-cascade")
    page.wait_for_function("EVT.cascade", timeout=120000)
    assert page.evaluate("EVT.cascade.declenchees.length") == 4
    page.wait_for_selector("#evt-rangee-cascade:not([hidden])")
    assert "s'arrete" in page.text_content("#evt-cascade-table")


def test_onglet_documentation_au_premier_clic(page):
    """Le premier clic sur l'onglet Documentation : ``activerOnglet`` et
    ``ouvrirDoc`` s'appelaient l'un l'autre sans fin (pile epuisee, rien ne
    s'affichait). Le cadre recoit la documentation, ou la page qui dit
    comment la construire ; la fixture verifie qu'aucune erreur n'est levee."""
    page.click('nav button[data-onglet="doc"]')
    page.wait_for_function(
        "document.querySelector('#cadre-doc').getAttribute('src') || "
        "document.querySelector('#cadre-doc').srcdoc",
        timeout=15000,
    )
    assert page.locator("#onglet-doc").is_visible()
    # y revenir ne redemande rien et ne leve rien
    page.click('nav button[data-onglet="constructeur"]')
    page.click('nav button[data-onglet="doc"]')
    assert page.locator("#onglet-doc").is_visible()
