"""L'interface en deux langues (bouton FR / EN), dans Chromium.

* la langue de depart est celle du navigateur (francais s'il est en francais,
  anglais sinon), ``?lang=`` dans l'adresse l'emporte, et le choix fait au
  bouton est retenu ;
* en anglais, les textes de la page, les infobulles et les textes construits
  par les scripts (messages a trous compris) sont traduits ; revenir au
  francais rend les textes d'origine, a l'identique ;
* les donnees de l'utilisateur ne sont pas touchees, et un calcul fait en
  anglais donne les memes nombres qu'en francais.

Les fixtures (serveur, navigateur) sont celles de ``test_navigateur``.
"""

from __future__ import annotations

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("playwright.sync_api")

import test_navigateur as _nav

serveur = _nav.serveur
navigateur = _nav.navigateur
JETON = _nav.JETON


def _ouvrir(serveur, navigateur, locale: str, suffixe: str = ""):
    ctx = navigateur.new_context(
        viewport={"width": 1500, "height": 1100}, locale=locale
    )
    pg = ctx.new_page()
    erreurs: list[str] = []
    pg.on("pageerror", lambda e: erreurs.append(str(e)))
    pg.goto(f"{serveur}/?jeton={JETON}{suffixe}")
    pg.wait_for_function("typeof E !== 'undefined' && E.defauts", timeout=30000)
    return ctx, pg, erreurs


def _onglets(pg) -> list[str]:
    return [t.strip() for t in pg.locator("nav button[data-onglet]").all_inner_texts()]


def test_langue_du_navigateur_puis_adresse(serveur, navigateur):
    ctx, pg, err = _ouvrir(serveur, navigateur, "fr-FR")
    assert pg.evaluate("[I18N.langue, document.documentElement.lang]") == ["fr", "fr"]
    assert "Reseau" in _onglets(pg) and "Network" not in _onglets(pg)
    assert (
        pg.locator('[data-langue="fr"]').get_attribute("class").split().count("actif")
    )
    ctx.close()
    ctx, pg, err = _ouvrir(serveur, navigateur, "de-DE")
    assert pg.evaluate("I18N.langue") == "en", "un navigateur non francophone : anglais"
    assert "Network" in _onglets(pg)
    ctx.close()
    ctx, pg, err = _ouvrir(serveur, navigateur, "en-US", "&lang=fr")
    assert pg.evaluate("I18N.langue") == "fr", "?lang= l'emporte sur le navigateur"
    assert not err
    ctx.close()


def test_bascule_traduit_puis_rend_le_francais_a_l_identique(serveur, navigateur):
    ctx, pg, err = _ouvrir(serveur, navigateur, "fr-FR")
    _nav._charger(pg, "smib")
    pg.wait_for_timeout(1200)
    texte = "document.querySelector('header').innerText + document.querySelector('nav').innerText + document.querySelector('#onglet-constructeur').innerText"
    titres = "[...document.querySelectorAll('header [title], nav [title]')].map(e => e.title)"
    fr, fr_titres = pg.evaluate(texte), pg.evaluate(titres)
    pg.click('[data-langue="en"]')
    assert pg.evaluate("[I18N.langue, document.documentElement.lang]") == ["en", "en"]
    en = pg.evaluate(texte)
    onglets = _onglets(pg)
    for attendu in (
        "Network",
        "Operating point",
        "Modes",
        "Time domain",
        "Blocks",
        "Code",
    ):
        assert attendu in onglets, onglets
    for mot in (
        "Charger",
        "Enregistrer",
        "Outils",
        "Proprietes",
        "Rien de selectionne",
    ):
        assert mot not in en, mot
    for mot in ("Load", "Save", "TOOLS", "Nothing selected."):
        assert mot in en, mot
    # les infobulles suivent, les donnees de l'utilisateur non
    assert (
        pg.get_attribute("#btn-nouveau", "title") == "Start again from an empty network"
    )
    assert pg.input_value("#nom-reseau") == "SMIB"
    # un texte construit par un script, avec des nombres (gabarit a trous)
    assert (
        "3 buses, 1 SM, 0 GFM, 0 GFL"
        in pg.locator("#catalogue option:checked").inner_text()
    )
    # le choix est retenu d'une page a l'autre
    pg.reload()
    pg.wait_for_function("typeof E !== 'undefined' && E.defauts", timeout=30000)
    assert pg.evaluate("I18N.langue") == "en" and "Network" in _onglets(pg)
    _nav._charger(pg, "smib")
    pg.wait_for_timeout(1200)
    # retour au francais : le texte d'origine, a l'identique
    pg.click('[data-langue="fr"]')
    assert pg.evaluate(texte) == fr
    assert pg.evaluate(titres) == fr_titres
    assert not err, err
    ctx.close()


def test_un_calcul_en_anglais_donne_les_memes_nombres(serveur, navigateur):
    """La traduction ne touche que des textes : le tableau des modes porte
    les memes valeurs dans les deux langues, sous des en-tetes traduits."""
    valeurs = {}
    for langue in ("fr", "en"):
        ctx, pg, err = _ouvrir(serveur, navigateur, "fr-FR", f"&lang={langue}")
        _nav._charger(pg, "wscc9")
        pg.click('nav button[data-onglet="modal"]')
        pg.click("#btn-modal")
        pg.wait_for_function("E.res && E.res.modal", timeout=120000)
        pg.wait_for_selector("#table-modal tbody tr, #onglet-modal table tbody tr")
        lignes = pg.locator("#onglet-modal table tbody tr").all_inner_texts()[:12]
        valeurs[langue] = lignes
        etat = pg.inner_text("#mod-bornes") + pg.inner_text("#onglet-modal")
        assert ("states" in etat) == (langue == "en"), etat[:300]
        assert not err, err
        ctx.close()
    assert valeurs["fr"] == valeurs["en"] and valeurs["fr"]


def test_message_d_erreur_du_serveur_garde_son_prefixe_traduit(serveur, navigateur):
    """Un message compose « prefixe : message du serveur » : le prefixe est
    traduit, le message du serveur reste tel quel (il n'est pas au
    dictionnaire) -- rien n'est perdu."""
    ctx, pg, err = _ouvrir(serveur, navigateur, "en-GB")
    assert pg.evaluate("I18N.t('echec : division par zero dans le bloc b12')") == (
        "failed: division par zero dans le bloc b12"
    )
    assert pg.evaluate("I18N.t('bus 4 a 0.25 s')") == "bus 4 a 0.25 s", (
        "inconnu : inchange"
    )
    assert pg.evaluate("I18N.t('V max : bus 4 a 0.25 s')") == "V max: bus 4 at 0.25 s"
    assert "bus 4 a 0.25 s" in pg.evaluate("I18N.manquants()")
    assert not err, err
    ctx.close()
