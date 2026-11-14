"""Les equations dynamiques modifiables (FR-052), sans HTTP.

Liste globale du reseau, modele a equations d'un element, edition (modifier,
ajouter, supprimer, revenir au natif), compilation a blanc, passerelles vers
l'editeur de blocs. Les fonctions levent
:class:`~g2egrid.elements.blocs.SchemaInvalide` (une ``ValueError``) avec un
message en francais : le transport en fait une 422.
"""

from __future__ import annotations

from g2egrid.elements import blocs, blocs_equations
from g2egrid.elements.blocs import SchemaInvalide
from g2egrid.reseau import equations_reseau as er
from g2egrid.reseau.tables import Tables, tables_vers_json


def _designation(table, index, element):
    if element:
        return element
    if table is None or index is None:
        raise SchemaInvalide("il faut l'element (element) ou sa ligne (table, index)")
    return (table, index)


def equations_liste_reseau(t: Tables) -> dict:
    """``{"lignes": [...]}`` : voir ``equations_reseau.equations_reseau``."""
    return {"lignes": er.equations_reseau(t)}


def equations_element(t: Tables, table=None, index=None, element=None) -> dict:
    """``{"modele": ...}`` : le modele a equations de l'element. Pour un
    element de ``T_Blocs``, en plus : ``natif`` (le nom de l'element natif
    qu'il remplace, ``None`` s'il n'en remplace aucun -- pas de retour au
    natif) et ``schema_bloc`` (son modele est un dessin, converti ici)."""
    el = _designation(table, index, element)
    rep = {"modele": er.modele_element(t, el)}
    _nom, tab, i = er._localiser(t, el)
    if tab == "T_Blocs":
        ligne = t.T_Blocs[i]
        brut = t.MODELES_BLOCS[ligne[2]]
        origine = (brut.get("origine") if isinstance(brut, dict) else None) or {}
        rep["natif"] = (
            origine.get("element") if origine.get("id_bloc") == int(ligne[0]) else None
        )
        rep["schema_bloc"] = not blocs_equations.est_modele_equations(brut)
    return rep


def _rendre(nouvelles: Tables, avant: Tables, element) -> dict:
    """``{tables, modele, table, index}`` : le modele de l'element dans les
    nouvelles tables et SON LIEU (un natif modifie est passe a la fin de
    ``T_Blocs``) -- l'interface suit l'element sans le deviner."""
    _nom, table, i = er._localiser(avant, element)
    if table != "T_Blocs" and len(nouvelles.T_Blocs) > len(avant.T_Blocs):
        cible = ("T_Blocs", len(nouvelles.T_Blocs) - 1)
    else:
        cible = (table, i)
    return {
        "tables": tables_vers_json(nouvelles),
        "modele": er.modele_element(nouvelles, cible),
        "table": cible[0],
        "index": cible[1],
    }


def equations_modifier(
    t: Tables, variable, expression, table=None, index=None, element=None, genre=None
) -> dict:
    el = _designation(table, index, element)
    return _rendre(er.modifier_equation(t, el, variable, expression, genre), t, el)


def equations_ajouter(
    t: Tables,
    genre,
    variable,
    expression,
    depart=None,
    table=None,
    index=None,
    element=None,
) -> dict:
    el = _designation(table, index, element)
    n = er.ajouter_equation(t, el, genre, variable, expression, depart)
    return _rendre(n, t, el)


def equations_supprimer(
    t: Tables, variable, table=None, index=None, element=None
) -> dict:
    el = _designation(table, index, element)
    return _rendre(er.supprimer_equation(t, el, variable), t, el)


def equations_parametre(
    t: Tables, nom, valeur, table=None, index=None, element=None
) -> dict:
    """``{tables, modele}`` : le parametre ``nom`` du MODELE de l'element vaut
    ``valeur``. Seulement pour un element deja decrit par un modele
    (``T_Blocs``) : les parametres d'un element natif se reglent par
    ``PARAMS_PAR_BUS`` (ses proprietes), sans le convertir."""
    el = _designation(table, index, element)
    _nom, tab, _i = er._localiser(t, el)
    if tab != "T_Blocs":
        raise SchemaInvalide(
            f"{_nom} est un element natif : ses parametres se reglent dans ses "
            "proprietes (PARAMS_PAR_BUS), il reste alors natif ; un parametre "
            "de modele ne se modifie ici que pour un element deja converti"
        )
    return _rendre(er.modifier_parametre(t, el, nom, valeur), t, el)


def blocs_remplacer(
    t: Tables, table=None, index=None, element=None, detail: str = "blocs"
) -> dict:
    """``{tables, modele, schema}`` : l'element natif est remplace, dans de
    nouvelles tables, par son schema-bloc (``modele`` : son nom dans
    ``MODELES_BLOCS``). Le serveur fait TOUT le remplacement (puissance d'une
    charge reportee au bus, ``Sn_VA``, origine gardee pour le retour au natif) :
    l'interface n'a pas a en reecrire les regles."""
    el = _designation(table, index, element)
    n, cle = er.remplacer_par_modele(t, el, "blocs", detail)
    return {
        "tables": tables_vers_json(n),
        "modele": cle,
        "schema": n.MODELES_BLOCS[cle],
        "table": "T_Blocs",
        "index": len(n.T_Blocs) - 1,
    }


def equations_revenir(t: Tables, table=None, index=None, element=None) -> dict:
    """``{tables, table, index}`` : l'element retrouve sa ligne native, dont
    le lieu est rendu."""
    n, (tab, i) = er.retour_au_natif(t, _designation(table, index, element))
    return {"tables": tables_vers_json(n), "table": tab, "index": i}


def equations_verifier(modele: dict) -> dict:
    """Compilation a blanc : ``{ok, message?, equations}``."""
    try:
        c = blocs.compiler(modele)
    except SchemaInvalide as err:
        return {"ok": False, "message": str(err), "equations": []}
    return {"ok": True, "equations": c.latex()}


def equations_depuis_schema(schema: dict) -> dict:
    return {"modele": blocs_equations.equations_depuis_schema(schema)}


def equations_vers_blocs(modele: dict) -> dict:
    return {"schema": blocs_equations.schema_depuis_equations(modele, "blocs")}
