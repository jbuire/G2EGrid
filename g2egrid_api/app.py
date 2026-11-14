"""L'interface web de g2egrid : une couche HTTP mince sur :mod:`analyse`.

Toutes les routes de calcul prennent **les tables du reseau** dans leur
corps, jamais un identifiant de reseau preenregistre. C'est ce qui fait que
le constructeur graphique n'est pas un cas particulier : un reseau dessine a
la main et un reseau du catalogue suivent exactement le meme chemin, et le
catalogue n'est qu'un moyen de remplir le formulaire de depart.

Lancement ::

    cd g2egrid
    uv run --extra api uvicorn g2egrid_api.app:app --port 8765 --reload
"""

from __future__ import annotations

import hmac
import logging
import os
import pathlib
import re
import uuid

from fastapi import FastAPI, HTTPException, Request
from fastapi.encoders import jsonable_encoder
from fastapi.exceptions import RequestValidationError
from fastapi.responses import (
    HTMLResponse,
    JSONResponse,
    PlainTextResponse,
    RedirectResponse,
    StreamingResponse,
)
from fastapi.staticfiles import StaticFiles

from g2egrid.noyau.progres import Interrompu

from . import analyse, taches
from .schemas import (
    CorpsApresEvenement,
    CorpsBalayage,
    CorpsBlocs,
    CorpsBlocsDegrouper,
    CorpsBlocsDepuisElement,
    CorpsBlocsRegrouper,
    CorpsBlocsRemplacer,
    CorpsBoucleMineure,
    CorpsCascade,
    CorpsContinuation,
    CorpsCourtCircuit,
    CorpsEquationsAjouter,
    CorpsEquationsElement,
    CorpsEquationsModele,
    CorpsEquationsModifier,
    CorpsEquationsParametre,
    CorpsEquationsSupprimer,
    CorpsFiche,
    CorpsFichier,
    CorpsGflAleatoires,
    CorpsImpact,
    CorpsImport,
    CorpsImportExterne,
    CorpsMachine,
    CorpsModal,
    CorpsMode,
    CorpsMulti,
    CorpsRapport,
    CorpsReglage,
    CorpsReparation,
    CorpsReseau,
    CorpsRobustesse,
    CorpsSimplification,
    CorpsStatismes,
    CorpsTemporel,
    CorpsUnSlack,
    CorpsZone,
)

WEB = pathlib.Path(__file__).parent / "web"

app = FastAPI(
    title="G2Egrid",
    description="Modelisation petit-signal dq : construction, load flow, analyse modale, temporel",
    version="0.1.0",
)

# ------------------------------------------------------------------ securite
# Un serveur local n'est pas un serveur prive : n'importe quelle page ouverte
# dans le navigateur peut envoyer des requetes a 127.0.0.1, et un nom de
# domaine qui se « rebinde » sur 127.0.0.1 peut meme en lire les reponses.
# Or l'atelier lit et ecrit des fichiers. D'ou deux gardes, actives des que
# G2EGRID_JETON est defini (``python -m g2egrid_api`` le tire au hasard) :
#
# * le JETON : l'URL ouverte au demarrage le porte (``/?jeton=...``) ; la page
#   d'accueil l'echange contre un cookie HttpOnly SameSite=Strict -- qu'une
#   autre origine ne peut ni lire ni envoyer -- et redirige vers ``/`` pour
#   l'oter de la barre d'adresse. Toute autre requete doit presenter le
#   cookie, ou l'en-tete ``X-G2egrid-Jeton`` (scripts) ;
# * l'HOTE : l'en-tete Host doit designer la machine locale (ou un nom de
#   G2EGRID_HOTES), ce qui ferme la porte au rebinding DNS.
#
# Sans G2EGRID_JETON (tests, ``uvicorn g2egrid_api.app:app`` a la main), rien ne
# change.
COOKIE_JETON = "g2egrid_jeton"
HOTES_LOCAUX = {"127.0.0.1", "localhost", "[::1]", "::1"}


def _meme_jeton(presente: str, jeton: str) -> bool:
    """Comparaison en temps constant, sur des octets.

    ``hmac.compare_digest`` refuse deux ``str`` non ASCII par un ``TypeError`` :
    un jeton ``?jeton=%C3%A9`` donnait une erreur 500 au lieu d'un refus.
    """
    return hmac.compare_digest(presente.encode(), jeton.encode())


def _hote_admis(hote: str) -> bool:
    nom = hote.rsplit(":", 1)[0] if not hote.endswith("]") else hote
    autres = {h.strip() for h in os.environ.get("G2EGRID_HOTES", "").split(",") if h}
    return nom in HOTES_LOCAUX or nom in autres


@app.middleware("http")
async def garde_locale(request: Request, suite):
    jeton = os.environ.get("G2EGRID_JETON")
    if not jeton:
        return await suite(request)
    if not _hote_admis(request.headers.get("host", "")):
        return PlainTextResponse("hote non admis", status_code=403)
    if request.url.path == "/" and "jeton" in request.query_params:
        if not _meme_jeton(request.query_params["jeton"], jeton):
            return PlainTextResponse("jeton invalide", status_code=403)
        # la langue demandee (?lang=fr|en) survit a la redirection ; rien
        # d'autre de l'adresse n'est repris
        langue = request.query_params.get("lang")
        cible = f"/?lang={langue}" if langue in ("fr", "en") else "/"
        r = RedirectResponse(cible, status_code=303)
        r.set_cookie(COOKIE_JETON, jeton, httponly=True, samesite="strict")
        return r
    presente = request.cookies.get(COOKIE_JETON) or request.headers.get(
        "x-g2egrid-jeton", ""
    )
    if presente and _meme_jeton(presente, jeton):
        return await suite(request)
    return HTMLResponse(
        "<!doctype html><meta charset=utf-8><title>g2egrid</title>"
        "<p>Cet atelier est protege par un jeton. Ouvrez l'adresse affichee "
        "dans la console au lancement (<code>python -m g2egrid_api</code>), "
        "elle se termine par <code>?jeton=...</code>.</p>"
        "<p>This workbench is protected by a token (<i>jeton</i>). Open the "
        "address printed in the console at start-up "
        "(<code>python -m g2egrid_api</code>): it ends with "
        "<code>?jeton=...</code>.</p>",
        status_code=403,
    )


def _tables(corps: CorpsReseau):
    """Les tables du corps de requete, ou une erreur 422 lisible."""
    try:
        return analyse.lire_tables(corps.tables)
    except Exception as err:
        raise HTTPException(422, f"tables illisibles : {err}") from err


def _calcul(fn, *a, **kw):
    """Execute un calcul metier ; une erreur de modele est une 422, pas une 500.

    Un reseau mal construit (pas de slack, charge a reactif negatif, boucle
    algebrique) n'est pas un bug du serveur : c'est une reponse a rendre a
    l'utilisateur, avec le message du modele.
    """
    try:
        return fn(*a, **kw)
    except (HTTPException, Interrompu):  # l'arret d'une tache n'est pas une erreur
        raise
    except KeyError as err:
        # str(KeyError(99)) vaut « 99 » : un message qui ne dit rien
        cle = err.args[0] if err.args else "?"
        raise HTTPException(
            422,
            f"reference introuvable : {cle!r}. Un bus, un equipement ou un "
            "parametre nomme dans les tables n'existe pas : verifiez les "
            "numeros de bus et les noms.",
        ) from err
    except FileNotFoundError as err:
        raise HTTPException(422, str(err) or "fichier ou reseau introuvable") from err
    except (ValueError, RuntimeError, ImportError, np_LinAlgError) as err:
        raise HTTPException(422, str(err) or type(err).__name__) from err
    except Exception as err:
        raise _erreur_interne(err) from err


def _erreur_interne(err: BaseException) -> HTTPException:
    """Une 500 SANS trace : un identifiant, et la trace dans le journal.

    La trace Python (chemins du disque, valeurs, lignes de code) n'a rien a
    faire dans la reponse : elle aide un attaquant et ne dit rien a
    l'utilisateur. Il garde l'identifiant, qui retrouve la trace dans la
    sortie du serveur."""
    ident = uuid.uuid4().hex[:8]
    journal.error("erreur interne %s", ident, exc_info=err)
    return HTTPException(
        500,
        f"erreur interne du serveur (reference {ident}) : la trace est dans le "
        "journal du serveur. Ce n'est pas une erreur de saisie.",
    )


journal = logging.getLogger("g2egrid_api")


@app.exception_handler(RequestValidationError)
async def corps_invalide(request: Request, err: RequestValidationError):
    """Un corps de requete refuse, en francais et en une phrase (l'interface
    affiche ``detail`` tel quel)."""
    return JSONResponse({"detail": _phrase_validation(err.errors())}, status_code=422)


@app.exception_handler(Exception)
async def erreur_non_prevue(request: Request, err: Exception):
    """Toute exception qui ne passe pas par :func:`_calcul` : 500 generique."""
    e = _erreur_interne(err)
    return JSONResponse({"detail": e.detail}, status_code=500)


_MOTS_VALIDATION = {
    "missing": "champ obligatoire absent",
    "greater_than": "doit etre strictement superieur a {gt}",
    "greater_than_equal": "doit etre superieur ou egal a {ge}",
    "less_than": "doit etre strictement inferieur a {lt}",
    "less_than_equal": "doit etre inferieur ou egal a {le}",
    "too_short": "doit contenir au moins {min_length} element(s)",
    "too_long": "doit contenir au plus {max_length} element(s)",
    "finite_number": "doit etre un nombre fini",
    "float_parsing": "doit etre un nombre",
    "float_type": "doit etre un nombre",
    "int_parsing": "doit etre un entier",
    "int_type": "doit etre un entier",
    "int_from_float": "doit etre un entier",
    "string_type": "doit etre du texte",
    "bool_parsing": "doit etre vrai ou faux",
    "dict_type": "doit etre un objet {{...}}",
    "list_type": "doit etre une liste",
    "json_invalid": "le corps n'est pas du JSON valide",
    "model_attributes_type": "doit etre un objet {{...}}",
}


def _phrase_validation(erreurs: list[dict]) -> str:
    morceaux = []
    for e in erreurs[:5]:
        champ = ".".join(str(x) for x in e.get("loc", ()) if x != "body") or "corps"
        modele = _MOTS_VALIDATION.get(e.get("type", ""))
        ctx = {
            k: (f"{v:g}" if isinstance(v, (int, float)) else v)
            for k, v in (e.get("ctx") or {}).items()
        }
        try:
            mot = modele.format(**ctx) if modele else e.get("msg", "")
        except (KeyError, IndexError):
            mot = e.get("msg", "valeur refusee")
        recu = e.get("input")
        morceaux.append(
            f"{champ} : {mot}"
            + (
                f" (recu : {recu!r})"
                if isinstance(recu, (int, float, str)) and len(repr(recu)) < 40
                else ""
            )
        )
    plus = f" (et {len(erreurs) - 5} autre(s))" if len(erreurs) > 5 else ""
    return "requete refusee - " + " ; ".join(morceaux) + plus


try:  # numpy.linalg.LinAlgError, sans importer numpy pour si peu
    from numpy.linalg import LinAlgError as np_LinAlgError
except Exception:  # pragma: no cover
    np_LinAlgError = ArithmeticError


# ------------------------------------------------------------------ catalogue
@app.get("/api/reseaux")
def reseaux():
    """Les reseaux de test disponibles, avec leur fiche."""
    return analyse.catalogue()


@app.get("/api/reseaux/{nom}")
def reseau(nom: str):
    """Les tables d'un reseau du catalogue, et un placement de ses bus."""
    try:
        return _calcul(analyse.depuis_catalogue, nom)
    except HTTPException as err:
        if err.status_code == 422 and "reseau inconnu" in str(err.detail):
            raise HTTPException(404, err.detail) from err
        raise


@app.get("/api/defauts")
def defauts():
    """Colonnes des tables, familles d'equipement, reglages par defaut."""
    return analyse.defauts()


# -------------------------------------------------------------------- calculs
@app.get("/api/equations")
def equations_liste():
    return analyse.equations_modele()


@app.post("/api/equations/reseau")
def equations_reseau(corps: CorpsReseau):
    return _calcul(analyse.equations_liste_reseau, _tables(corps))


@app.post("/api/equations/element")
def equations_element(corps: CorpsEquationsElement):
    return _calcul(
        analyse.equations_element,
        _tables(corps),
        corps.table,
        corps.index,
        corps.element,
    )


@app.post("/api/equations/modifier")
def equations_modifier(corps: CorpsEquationsModifier):
    return _calcul(
        analyse.equations_modifier,
        _tables(corps),
        corps.variable,
        corps.expression,
        corps.table,
        corps.index,
        corps.element,
        corps.genre,
    )


@app.post("/api/equations/ajouter")
def equations_ajouter(corps: CorpsEquationsAjouter):
    return _calcul(
        analyse.equations_ajouter,
        _tables(corps),
        corps.genre,
        corps.variable,
        corps.expression,
        corps.depart,
        corps.table,
        corps.index,
        corps.element,
    )


@app.post("/api/equations/supprimer")
def equations_supprimer(corps: CorpsEquationsSupprimer):
    return _calcul(
        analyse.equations_supprimer,
        _tables(corps),
        corps.variable,
        corps.table,
        corps.index,
        corps.element,
    )


@app.post("/api/equations/parametre")
def equations_parametre(corps: CorpsEquationsParametre):
    return _calcul(
        analyse.equations_parametre,
        _tables(corps),
        corps.nom,
        corps.valeur,
        corps.table,
        corps.index,
        corps.element,
    )


@app.post("/api/blocs/remplacer")
def blocs_remplacer(corps: CorpsBlocsRemplacer):
    return _calcul(
        analyse.blocs_remplacer,
        _tables(corps),
        corps.table,
        corps.index,
        corps.element,
        corps.detail,
    )


@app.post("/api/equations/revenir")
def equations_revenir(corps: CorpsEquationsElement):
    return _calcul(
        analyse.equations_revenir,
        _tables(corps),
        corps.table,
        corps.index,
        corps.element,
    )


@app.post("/api/equations/verifier")
def equations_verifier(corps: CorpsEquationsModele):
    return _calcul(analyse.equations_verifier, corps.modele)


@app.post("/api/equations/depuis_schema")
def equations_depuis_schema(corps: CorpsBlocs):
    return _calcul(analyse.equations_depuis_schema, corps.schema_bloc)


@app.post("/api/equations/vers_blocs")
def equations_vers_blocs(corps: CorpsEquationsModele):
    return _calcul(analyse.equations_vers_blocs, corps.modele)


@app.get("/api/equations/{nom}")
def equations_modele(nom: str):
    return _calcul(analyse.equations_modele, nom)


@app.get("/api/blocs/types")
def blocs_types():
    return analyse.blocs_types()


@app.post("/api/blocs/compiler")
def blocs_compiler(corps: CorpsBlocs):
    return _calcul(analyse.blocs_compiler, corps.schema_bloc)


@app.post("/api/blocs/depuis_element")
def blocs_depuis_element(corps: CorpsBlocsDepuisElement):
    t = None
    if corps.table is not None:
        t = _tables(CorpsReseau(tables=corps.tables or {}))
    return _calcul(
        analyse.blocs_depuis_element,
        corps.classe,
        corps.parametres,
        corps.detail,
        t,
        corps.table,
        corps.index,
        corps.f0,
    )


@app.post("/api/blocs/regrouper")
def blocs_regrouper(corps: CorpsBlocsRegrouper):
    return _calcul(
        analyse.blocs_regrouper, corps.schema_bloc, corps.ids, corps.nom, corps.chemin
    )


@app.post("/api/blocs/degrouper")
def blocs_degrouper(corps: CorpsBlocsDegrouper):
    return _calcul(analyse.blocs_degrouper, corps.schema_bloc, corps.id, corps.chemin)


@app.post("/api/blocs/aplatir")
def blocs_aplatir(corps: CorpsBlocs):
    return _calcul(analyse.blocs_aplatir, corps.schema_bloc)


def _modele_absent(fn, *a):
    """Un modele de la bibliotheque qui n'existe pas est une 404."""

    def appel():
        try:
            return fn(*a)
        except analyse.ModeleAbsent as err:  # par son TYPE, pas par son texte
            raise HTTPException(404, str(err)) from None

    return _calcul(appel)


@app.get("/api/blocs/bibliotheque")
def blocs_bibliotheque():
    return _calcul(analyse.bibliotheque_liste)


@app.get("/api/blocs/bibliotheque/{nom}")
def blocs_bibliotheque_lire(nom: str):
    return _modele_absent(analyse.bibliotheque_lire, nom)


@app.post("/api/blocs/bibliotheque")
def blocs_bibliotheque_enregistrer(corps: CorpsBlocs):
    return _calcul(analyse.bibliotheque_enregistrer, corps.schema_bloc)


@app.delete("/api/blocs/bibliotheque/{nom}")
def blocs_bibliotheque_supprimer(nom: str):
    return _modele_absent(analyse.bibliotheque_supprimer, nom)


@app.post("/api/fiche")
def fiche_element(c: CorpsFiche):
    return _calcul(analyse.fiche_element, _tables(c), c.table, c.index)


@app.post("/api/reglage/boucles")
def reglage_boucles(corps: CorpsReglage):
    return _calcul(analyse.reglage_boucles, corps.famille, corps.reglages, corps.f0)


@app.post("/api/machine/standard")
def machine_standard(corps: CorpsMachine):
    return _calcul(analyse.machine_standard, corps.valeurs, corps.f0, corps.sens)


@app.post("/api/disjoncteurs")
def disjoncteurs(corps: CorpsReseau):
    return _calcul(analyse.disjoncteurs, _tables(corps))


@app.post("/api/valider")
def valider(corps: CorpsReseau):
    problemes = _calcul(analyse.valider, _tables(corps), corps.positions)
    return {
        "ok": not any(p["gravite"] == "erreur" for p in problemes),
        "problemes": problemes,
    }


@app.post("/api/reparer")
def reparer(corps: CorpsReparation):
    """Applique une reparation proposee par /api/valider ; rend les tables."""
    t = _calcul(
        analyse.reparer, _tables(corps), corps.action, corps.args, corps.positions
    )
    return {"tables": analyse.tables_vers_json(t)}


@app.post("/api/topologie")
def topologie(corps: CorpsReseau):
    return _calcul(analyse.topologie, _tables(corps), corps.positions)


@app.post("/api/court_circuit")
def court_circuit(corps: CorpsCourtCircuit):
    return _calcul(
        analyse.court_circuit,
        _tables(corps),
        corps.gfm,
        corps.gfl,
        corps.charges,
        corps.dynamique,
        corps.f,
        corps.machine,
    )


@app.post("/api/loadflow")
def loadflow(corps: CorpsReseau):
    return _calcul(analyse.loadflow, _tables(corps))


@app.post("/api/modal")
def modal(corps: CorpsModal):
    return _calcul(analyse.modes, _tables(corps), corps.f_max, corps.n_etats)


@app.post("/api/mode")
def mode(corps: CorpsMode):
    return _calcul(
        analyse.detail_mode,
        _tables(corps),
        corps.indice,
        corps.n_etats,
        corps.n_sensibilites,
        corps.n_parametres,
    )


@app.post("/api/balayage")
def balayage(corps: CorpsBalayage):
    return _calcul(
        analyse.balayage,
        _tables(corps),
        corps.cible,
        corps.valeurs,
        corps.f_max,
        corps.n_modes,
    )


@app.post("/api/point/un_slack")
def point_un_slack(corps: CorpsUnSlack):
    return _calcul(analyse.un_slack_tout_pq, _tables(corps), corps.slack_bus)


@app.post("/api/t_eig")
def t_eig(corps: CorpsReseau):
    return _calcul(analyse.t_eig, _tables(corps))


@app.post("/api/robustesse")
def robustesse(corps: CorpsRobustesse):
    return _calcul(
        analyse.robustesse,
        _tables(corps),
        corps.f_min,
        corps.f_max,
        corps.zeta_max,
        corps.frac,
        corps.puissance,
        corps.deconnexion,
    )


# ------------------------------------------------------------------- taches
# Les calculs longs, en tache de fond (voir g2egrid_api.taches) : le corps est
# celui de la route synchrone du meme nom, et le resultat aussi.
@app.post("/api/impact_continu")
def impact_continu(corps: CorpsReseau):
    """Les modes que la dynamique du reseau continu deplace (continu rigide
    contre continu dynamique)."""
    return _calcul(analyse.impact_continu, _tables(corps))


@app.post("/api/elements")
def elements(corps: CorpsReseau):
    """Les elements assembles et leurs parametres."""
    return _calcul(analyse.elements_reseau, _tables(corps))


def _statismes(c: CorpsStatismes):
    return _calcul(
        analyse.balayage_statismes, _tables(c), c.debut, c.fin, c.pas, c.stations
    )


def _continuation(c: CorpsContinuation):
    return _calcul(
        analyse.continuation_parametre,
        _tables(c),
        c.element,
        c.parametre,
        c.debut,
        c.fin,
        c.pas,
        c.log,
        c.equilibre,
    )


def _boucle(c: CorpsBoucleMineure):
    return _calcul(analyse.boucle_mineure, _tables(c), c.element, c.f_min, c.f_max, c.n)


def _rapport(c: CorpsRapport):
    from .rapport import empreinte, rapport_etude

    t = _tables(c)
    page = _calcul(
        rapport_etude,
        t,
        analyse.construire,
        c.figures,
        c.n_modes,
        positions=c.positions,
    )
    nom = re.sub(r"[^0-9A-Za-z_-]", "_", t.nom or "reseau").strip("_") or "reseau"
    return {"html": page, "fichier": f"rapport_{nom}.html", "empreinte": empreinte(t)}


def _apres_evenement(c: CorpsApresEvenement):
    return _calcul(
        analyse.apres_evenement,
        _tables(c),
        c.evenements,
        c.duree,
        c.n,
        c.pas,
        c.limiteurs,
        c.relais,
        c.protections_dc,
        c.v_bas,
        c.v_haut,
        c.fenetre,
    )


def _cascade(c: CorpsCascade):
    return _calcul(
        analyse.cascade, _tables(c), c.initiales, c.v_max, c.familles, c.ordre
    )


@app.post("/api/apres_evenement")
def apres_evenement(corps: CorpsApresEvenement):
    """Tension de chaque bus, P et Q de chaque unite, frequence et
    declenchements des relais apres un evenement (simulation BDF2)."""
    return _apres_evenement(corps)


@app.post("/api/cascade")
def cascade(corps: CorpsCascade):
    """La cascade de surtension en regimes permanents successifs."""
    return _cascade(corps)


GENRES_TACHE: dict[str, tuple] = {
    "rapport": (CorpsRapport, _rapport),
    "robustesse": (CorpsRobustesse, robustesse),
    "statismes": (CorpsStatismes, _statismes),
    "continuation": (CorpsContinuation, _continuation),
    "impedance": (CorpsBoucleMineure, _boucle),
    "apres_evenement": (CorpsApresEvenement, _apres_evenement),
    "cascade": (CorpsCascade, _cascade),
}


@app.post("/api/taches/{genre}")
async def tache_lancer(genre: str, request: Request):
    if genre not in GENRES_TACHE:
        raise HTTPException(404, f"pas de tache {genre!r}")
    modele, fonction = GENRES_TACHE[genre]
    try:
        corps = modele.model_validate(await request.json())
    except Exception as err:
        raise HTTPException(422, f"corps illisible : {err}") from err
    try:
        return taches.lancer(genre, lambda: fonction(corps)).resume()
    except taches.TropDeTaches as err:
        raise HTTPException(429, str(err)) from err


@app.get("/api/taches")
def tache_liste():
    return taches.liste()


@app.get("/api/taches/{id_}")
def tache_etat(id_: str):
    t = taches.lire(id_)
    if t is None:
        raise HTTPException(404, "tache inconnue (trop ancienne ?)")
    return t.resume(avec_resultat=True)


@app.delete("/api/taches/{id_}")
def tache_annuler(id_: str):
    t = taches.annuler(id_)
    if t is None:
        raise HTTPException(404, "tache inconnue")
    return t.resume()


@app.get("/api/taches/{id_}/flux")
async def tache_flux(id_: str):
    """L'avancement en flux SSE : un evenement a chaque changement, le
    dernier avec le resultat ; puis le flux se ferme."""
    import asyncio
    import json

    t = taches.lire(id_)
    if t is None:
        raise HTTPException(404, "tache inconnue")

    async def flux():
        vue, silence = -1, 0
        while True:
            if t.revision != vue:
                vue, silence = t.revision, 0
                d = jsonable_encoder(t.resume(avec_resultat=True))
                yield f"data: {json.dumps(d)}\n\n"
                if t.terminee:
                    return
            else:
                silence += 1
                if silence % 75 == 0:  # toutes les 15 s : garder la connexion
                    yield ": attente\n\n"
            await asyncio.sleep(0.2)

    return StreamingResponse(
        flux(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


@app.post("/api/impact")
def impact(corps: CorpsImpact):
    return _calcul(analyse.impact, corps.reseau, corps.f_max, corps.par_poste)


@app.post("/api/code")
def code(corps: CorpsFichier):
    return _calcul(analyse.code_python, corps.reseau, corps.nom)


@app.post("/api/signaux")
def signaux(corps: CorpsReseau):
    return _calcul(analyse.signaux, _tables(corps))


@app.post("/api/temporel")
def temporel(corps: CorpsTemporel):
    return _calcul(
        analyse.temporel,
        _tables(corps),
        corps.entree,
        corps.amplitude,
        corps.duree,
        corps.t0,
        corps.n,
        corps.etats,
        corps.voie,
        corps.evenements,
        corps.mesures,
        corps.methode,
        corps.pas,
        corps.limiteurs,
        corps.protections_dc,
    )


@app.post("/api/exporter")
def exporter(corps: CorpsReseau):
    return {"python": _calcul(analyse.exporter, _tables(corps))}


# ------------------------------------------------------------------ fichiers
@app.get("/api/fichiers")
def fichiers():
    """Les reseaux enregistres au format natif, avec leur en-tete.

    Rend une LISTE, comme ``/api/reseaux`` : les deux repondent a la meme
    question -- quels reseaux puis-je ouvrir -- et l'interface les affiche
    cote a cote.
    """
    return _calcul(analyse.fichiers)


@app.post("/api/fichiers")
def enregistrer_fichier(corps: CorpsFichier):
    """Enregistre un reseau STRUCTURE : transport, zones et instances.

    A ne pas confondre avec ``/api/exporter`` : celui-la rend un texte a
    telecharger, celui-ci ecrit un fichier que le catalogue saura relire, donc
    que ``GET /api/reseaux/{nom}`` rouvrira tel quel.
    """
    return _calcul(analyse.enregistrer_fichier, corps.reseau, corps.nom)


@app.post("/api/fichiers/lire")
def lire_fichier(corps: CorpsFichier):
    """Un fichier ouvert par l'utilisateur, rendu sous la forme de lecture."""
    return _calcul(analyse.lire_fichier, corps.reseau)


@app.post("/api/fichiers/externe")
def importer_externe(corps: CorpsImportExterne):
    """Un fichier MATPOWER, PSS/E ou pandapower, lu par g2egrid.formats.interop."""
    return _calcul(analyse.importer_externe, corps.nom, corps.contenu, corps.charges)


@app.post("/api/fichiers/convertir")
def convertir_fichier(corps: CorpsFichier):
    """Le reseau du constructeur, au format ``.g2egrid.json``, a telecharger."""
    return _calcul(analyse.convertir_fichier, corps.reseau)


# --------------------------------------------------------------------- zones
@app.post("/api/zones/aplatir")
def zones_aplatir(corps: CorpsFichier):
    """Un reseau a zones, rendu a plat, sans passer par un fichier."""
    return _calcul(analyse.aplatir_reseau, corps.reseau)


@app.post("/api/reseau/importer")
def reseau_importer(corps: CorpsImport):
    """Importer un reseau (catalogue ou fichier) et le raccorder."""
    return _calcul(
        analyse.importer,
        corps.reseau,
        corps.source,
        corps.raccord,
        corps.garder_slack,
        corps.garder_ib,
    )


@app.post("/api/multi/ajouter_gfl")
def multi_ajouter_gfl(corps: CorpsGflAleatoires):
    """Des GFL ajoutes a des bus tires au hasard dans la distribution."""
    return _calcul(
        analyse.multi_ajouter_gfl, corps.reseau, corps.n, corps.graine, corps.part
    )


@app.post("/api/simplification/inventaire")
def simplification_inventaire(corps: CorpsReseau):
    """Les types d'elements du reseau, leurs etats, et les preregles."""
    return _calcul(analyse.simplification_inventaire, _tables(corps))


@app.post("/api/simplification/comparer")
def simplification_comparer(corps: CorpsSimplification):
    """Les modes sans puis avec les hypotheses de simplification."""
    return _calcul(
        analyse.simplification_comparer,
        _tables(corps),
        corps.simplifications,
        corps.preregles,
    )


@app.post("/api/multi/comparer")
def multi_comparer(corps: CorpsMulti):
    """Transport seul et distribution seule, confrontes au modele complet."""
    return _calcul(analyse.multi_comparer, corps.reseau, corps.instances)


@app.post("/api/multi/penetration")
def multi_penetration(corps: CorpsMulti):
    """Le balayage de penetration GFL."""
    return _calcul(
        analyse.multi_penetration,
        corps.reseau,
        corps.n_max,
        corps.pas,
        corps.graine,
        corps.part,
    )


@app.post("/api/zones/definir")
def zones_definir(corps: CorpsZone):
    """Faire d'un ensemble de bus une zone raccordee au transport."""
    return _calcul(analyse.definir_zone, corps.reseau, corps.bus, corps.nom)


# ------------------------------------------------------------------- exports
@app.post("/api/exporter/{quoi}")
def exporter_vers(quoi: str, corps: CorpsReseau):
    """``donnees``, ``notebook`` ou ``matlab``, avec ce que l'export y perd.

    Les avertissements sont dans la reponse, et l'interface doit les montrer :
    un ``.m`` dont les shunts sont partis en commentaire reste executable, et
    decrit pourtant un autre reseau que celui qu'on a dessine.
    """
    return _calcul(analyse.exporter_vers, _tables(corps), quoi)


@app.post("/api/cache/vider")
def vider():
    return {"vides": analyse.vider_cache()}


# ------------------------------------------------------------------- statique
@app.get("/")
def accueil():
    """La page, avec une empreinte de version sur ses fichiers statiques.

    LE CACHE DU NAVIGATEUR EST UN PIEGE DE DEVELOPPEMENT, pas un detail de
    confort : le fichier sur le disque etait a jour, la page servie aussi, et
    c'est le SCRIPT que le navigateur continuait a prendre dans son cache --
    si bien qu'on lit un comportement ancien en croyant tester le nouveau, et
    qu'on cherche un defaut la ou il n'y en a pas. Accrocher la date de
    derniere modification a l'URL suffit : elle change quand le fichier
    change, et ne change pas autrement.
    """
    html = (WEB / "index.html").read_text(encoding="utf-8")
    # TOUS les scripts et la feuille de style : equations.js en etait absent,
    # et le navigateur en gardait l'ancienne version
    for f in sorted([*WEB.glob("*.js"), *WEB.glob("*.css")]):
        v = int(f.stat().st_mtime)
        html = html.replace(f'/web/{f.name}"', f'/web/{f.name}?v={v}"')
    return HTMLResponse(html)


app.mount("/web", StaticFiles(directory=WEB), name="web")

# --------------------------------------------------------------- documentation
# La documentation Sphinx, servie telle que construite dans docs/_build/html.
# Montee meme si elle n'existe pas encore (check_dir=False) : l'onglet
# Documentation propose de la construire, et elle apparait sans redemarrer.
DOC = pathlib.Path(__file__).resolve().parent.parent / "docs"
DOC_HTML = DOC / "_build" / "html"
_construction_doc: dict = {"etat": "inconnu", "message": ""}


@app.get("/api/doc")
def doc_etat():
    """La documentation est-elle construite ? Et sa construction en cours ?"""
    return {
        "construite": (DOC_HTML / "index.html").exists(),
        "construction": _construction_doc["etat"],
        "message": _construction_doc["message"][-2000:],
    }


@app.post("/api/doc/construire")
def doc_construire():
    """Construit la documentation en tache de fond (Sphinx doit etre installe :
    ``uv sync --group docs``)."""
    import subprocess
    import sys
    import threading

    if _construction_doc["etat"] == "en cours":
        return doc_etat()
    try:
        import sphinx  # noqa: F401
    except ImportError as err:
        raise HTTPException(
            422,
            "Sphinx n'est pas installe : uv sync --group docs, puis relancer l'atelier",
        ) from err

    def construire():
        r = subprocess.run(
            [
                sys.executable,
                "-m",
                "sphinx",
                "-b",
                "html",
                "-q",
                str(DOC),
                str(DOC_HTML),
            ],
            capture_output=True,
            text=True,
        )
        _construction_doc["etat"] = "ok" if r.returncode == 0 else "echec"
        _construction_doc["message"] = r.stderr

    _construction_doc.update(etat="en cours", message="")
    threading.Thread(target=construire, daemon=True).start()
    return doc_etat()


_PAGE_DOC_ABSENTE = (
    "<!doctype html><meta charset=utf-8><title>g2egrid : documentation</title>"
    "<p>La documentation n'est pas encore construite. Dans l'atelier, onglet "
    "Documentation, lancez la construction (Sphinx doit etre installe : "
    "<code>uv sync --group docs</code>), puis rechargez cette page.</p>"
)


class _StatiqueDoc(StaticFiles):
    """Les fichiers de la doc ; une page claire (404) tant qu'elle n'est pas
    construite. ``StaticFiles(check_dir=False)`` seul levait une
    ``RuntimeError`` -- une 500 -- sur ``/doc/`` avant la construction."""

    async def check_config(self):
        """Pas de controle du dossier : il peut apparaitre apres le demarrage
        (construction depuis l'onglet Documentation)."""

    async def get_response(self, path, scope):
        if not (DOC_HTML / "index.html").exists():
            return HTMLResponse(_PAGE_DOC_ABSENTE, status_code=404)
        return await super().get_response(path, scope)


app.mount(
    "/doc", _StatiqueDoc(directory=DOC_HTML, html=True, check_dir=False), name="doc"
)


# --------------------------------------------------------------- taille limite
# Un corps de requete est lu en entier avant d'etre decode : sans plafond, une
# seule requete de quelques Go epuise la memoire du serveur. Le plus gros
# reseau du catalogue (ieee39_mixte, a plat) pese 83 ko en JSON (mesure :
# g2egrid/tests/test_securite_atelier.py) ; le plafond par defaut en laisse
# 8 Mio, soit une marge de l'ordre de 100 : un reseau MATPOWER de plusieurs
# milliers de bus tient encore. ``G2EGRID_TAILLE_MAX`` (en octets) le change.
TAILLE_MAX_DEFAUT = 8 * 1024 * 1024


def taille_max_requete() -> int:
    """Le plafond, en octets, du corps d'une requete (lu a chaque requete)."""
    try:
        v = int(os.environ.get("G2EGRID_TAILLE_MAX", ""))
    except ValueError:
        return TAILLE_MAX_DEFAUT
    return v if v > 0 else TAILLE_MAX_DEFAUT


class _TropGros(Exception):
    pass


class LimiteTaille:
    """Middleware ASGI : 413 si le corps depasse :func:`taille_max_requete`.

    On refuse sur l'en-tete ``Content-Length`` avant de lire quoi que ce soit,
    et on COMPTE les octets lus pour un corps sans longueur annoncee (transfert
    par morceaux), qui contournerait sinon l'en-tete."""

    def __init__(self, application):
        self.application = application

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            return await self.application(scope, receive, send)
        limite = taille_max_requete()
        annonce = dict(scope["headers"]).get(b"content-length", b"")
        if annonce.isdigit() and int(annonce) > limite:
            return await self._refuser(scope, receive, send, limite)
        recu = 0
        depasse = False
        commence = False

        async def recevoir():
            nonlocal recu, depasse
            m = await receive()
            if m["type"] == "http.request":
                recu += len(m.get("body", b""))
                if recu > limite:
                    # FastAPI avale cette exception en « 400 corps illisible » :
                    # on la leve pour arreter la lecture, et c'est ``envoyer``
                    # qui remplace la reponse de l'application par le 413
                    depasse = True
                    raise _TropGros
            return m

        async def envoyer(m):
            nonlocal commence
            if depasse:
                return
            commence = commence or m["type"] == "http.response.start"
            await send(m)

        try:
            await self.application(scope, recevoir, envoyer)
        except _TropGros:
            pass
        if depasse and not commence:
            await self._refuser(scope, receive, send, limite)

    @staticmethod
    async def _refuser(scope, receive, send, limite):
        r = JSONResponse(
            {
                "detail": f"requete trop volumineuse (plus de {limite} octets). "
                "Le plafond se regle par la variable G2EGRID_TAILLE_MAX."
            },
            status_code=413,
        )
        await r(scope, receive, send)


app.add_middleware(LimiteTaille)
