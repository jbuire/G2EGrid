"""Export des tables vers un fichier reseau MATLAB.

POURQUOI UN EXPORT ALTERNATIF. La version MATLAB de G2ELIN lit un reseau par
``buildFromMatrices``, qui prend treize matrices et rien d'autre. Les fichiers
de ``matlab/G2Elin_V2/reseaux_test/reseaux_classiques/`` (``WSCC.m`` est le
modele) sont exactement cela : des donnees, puis un boilerplate. Ce module
refait ce fichier a partir d'un :class:`~g2egrid.reseau.tables.Tables`, ce qui permet
de rejouer cote MATLAB un reseau dessine ou calcule cote Python.

CE QUI NE PASSE PAS LA FRONTIERE, ET POURQUOI ON LE DIT. ``buildFromMatrices``
ne connait ni les shunts (``WSCC3SM_CIGRE3.m`` les ajoute apres coup, noeud par
noeud, par ``addShunt``), ni les charges a courant ou a puissance constante,
ni les reglages fins des machines. Plutot que de les perdre en silence, on les
ecrit EN COMMENTAIRE dans le fichier produit et on rend un avertissement par
omission : le fichier reste executable, et le lecteur sait ce qu'il lui manque.

LE CONTINU EST UN REFUS, PAS UN AVERTISSEMENT. ``T_BusDC``, ``T_LinesDC``,
``T_ShuntsDC`` et ``T_VSCs`` n'ont aucun equivalent dans ``buildFromMatrices``,
et un reseau a liaison continue ampute de sa liaison n'est pas une
approximation du reseau : c'est un autre reseau, qui converge et qui ment.

``T_Loads`` N'EST PAS LA MEME TABLE DE PART ET D'AUTRE. Cote MATLAB elle porte
les charges DYNAMIQUES ``[id, busID, P_W, Q_var]``, celles que Python fabrique
a partir des colonnes ``Pl_W``/``Ql_var`` de ``Y_DER`` ; la table ``T_Loads``
de Python, elle, ne contient que les charges SPECIALES (I et P). C'est la
seule correspondance non triviale de l'export, et ``WSCC.m`` la montre : ses
trois charges figurent a la fois dans ``Y_DER`` et dans ``T_Loads``.
"""

from __future__ import annotations

import json

import numpy as np

import re

from g2egrid.elements.bus import NOEUDS_DEFAUT, R_CAUSALITE
from g2egrid.reseau.tables import TABLES, Tables

__all__ = ["vers_matlab"]

#: Les tables passees a ``buildFromMatrices``, DANS SON ORDRE D'ARGUMENTS.
ORDRE_MATLAB = (
    "Y_DER",
    "T_Lines_SI",
    "T_Lines_pu",
    "T_Transfos",
    "T_IBs",
    "T_SMs",
    "T_GFM1s",
    "T_GFM2s",
    "T_GFM3s",
    "T_GFM4s",
    "T_GFL1s",
    "T_GFL2s",
    "T_Loads",
)

#: Les tables continues, qui n'ont pas d'equivalent MATLAB.
CONTINU = ("T_BusDC", "T_LinesDC", "T_ShuntsDC", "T_VSCs", "T_DCDC", "T_LoadsDC")

#: Le nombre de colonnes que ``buildFromMatrices`` lit. Au-dela, ce sont les
#: colonnes optionnelles Rt/Xt (machines) ou R2/L2 (convertisseurs), qui cote
#: MATLAB se reglent dans l'objet apres construction.
COLONNES_MATLAB = {nom: len(TABLES[nom]) for nom in ORDRE_MATLAB}
COLONNES_MATLAB["T_Loads"] = 4  # [id, busID, P_W, Q_var], sans la colonne Type


def _nom_fonction(nom: str) -> str:
    """Un identifiant MATLAB tire du nom du reseau."""
    ident = re.sub(r"[^0-9A-Za-z_]", "_", nom).strip("_") or "reseau"
    return ident if ident[0].isalpha() else f"N{ident}"


def _nombre(v) -> str:
    """Un nombre, ecrit sans perte et sans notation inutile."""
    if isinstance(v, str):
        return f"'{v}'"
    f = float(v)
    if f.is_integer() and abs(f) < 1e15:
        return str(int(f))
    return repr(f)


def _matrice(
    nom: str, entete: str, lignes: list[list[str]], indent: str = "    "
) -> list[str]:
    """Une affectation de matrice MATLAB, colonnes alignees comme ``WSCC.m``."""
    out = [f"{indent}% {nom} : {entete}"]
    if not lignes:
        out.append(f"{indent}{nom} = [];")
        return out
    larg = [max(len(r[c]) for r in lignes) for c in range(len(lignes[0]))]
    out.append(f"{indent}{nom} = [")
    for r in lignes:
        cells = "  ".join(v.rjust(larg[c]) for c, v in enumerate(r))
        out.append(f"{indent}    {cells} ;")
    out.append(f"{indent}];")
    return out


def _tronquer(nom: str, rows, avertissements: list[str]) -> list[list[str]]:
    """Les lignes d'une table, ramenees aux colonnes que MATLAB lit."""
    n = COLONNES_MATLAB[nom]
    trop = [r for r in rows if len(r) > n]
    if trop:
        avertissements.append(
            f"{nom} : {len(trop)} ligne(s) ont des colonnes optionnelles "
            f"(au-dela de {n}) que buildFromMatrices ne lit pas ; elles sont "
            "omises. Cote MATLAB, ces reglages se posent sur l'objet apres "
            "construction."
        )
    return [[_nombre(v) for v in r[:n]] for r in rows]


def vers_matlab(t: Tables) -> tuple[str, list[str]]:
    """Le source d'un fichier reseau MATLAB, et ce qu'il n'a pas pu emporter.

    Rend ``(source, avertissements)``. ``source`` est un ``.m`` complet, au
    format de ``matlab/G2Elin_V2/reseaux_test/reseaux_classiques/WSCC.m``, pret
    a etre appele ; ``avertissements`` est la liste, en clair, de ce que
    ``buildFromMatrices`` ne sait pas prendre et qui n'est donc present dans le
    fichier qu'en commentaire.

    Leve ``ValueError`` si le reseau a une partie continue : voir le module.
    """
    presentes = [nom for nom in CONTINU if getattr(t, nom)]
    if presentes:
        raise ValueError(
            f"{t.nom} : reseau a partie continue ({', '.join(presentes)}). "
            "buildFromMatrices ne prend que les treize tables alternatives ; "
            "exporter sans la liaison continue ne donnerait pas une "
            "approximation du reseau mais un autre reseau. Export refuse."
        )

    av: list[str] = []
    fn = _nom_fonction(t.nom)

    # Les charges dynamiques MATLAB : les colonnes Pl/Ql de Y_DER.
    consommateurs = [r for r in t.Y_DER if r[7] or r[8]]
    charges = [
        [_nombre(k), _nombre(int(r[0])), _nombre(r[7]), _nombre(r[8])]
        for k, r in enumerate(consommateurs, start=1)
    ]

    L = [
        f"function net = {fn}()",
        f"%{fn.upper()}  Reseau exporte de g2egrid.",
        "%",
        "%  Genere par g2egrid.formats.export.vers_matlab a partir des tables du reseau",
        f'%  "{" ".join(str(t.nom).split())}". Seule la section DONNEES change '
        "d'un reseau a l'autre ;",
        "%  le reste est le boilerplate de buildFromMatrices.",
        "%",
        f"%  Base systeme : {t.Sb / 1e6:g} MVA, {t.f0:g} Hz.",
        "%",
        "%  USAGE :",
        f"%    net = {fn}();",
        "%    net.displayLoadFlow();                   % tableau load flow",
        "%    [Sys, ev, LM] = net.buildGlobalModel();  % modele petit-signal",
        "",
        "    %% ===================== PARAMETRES =====================",
        f"    Parameters.fb = {_nombre(t.f0)};        % frequence reseau (Hz)",
        f"    Parameters.Sb = {_nombre(t.Sb)};     % puissance de base systeme (VA)",
        "",
        "    %% ===================== DONNEES =====================",
    ]

    for nom in ORDRE_MATLAB:
        L.append("")
        if nom == "T_Loads":
            L.extend(
                _matrice(
                    "T_Loads",
                    "[id, busID, P_W, Q_var]   (charges dynamiques : "
                    "colonnes Pl_W/Ql_var de Y_DER)",
                    charges,
                )
            )
            continue
        L.extend(
            _matrice(
                nom,
                "[" + ", ".join(TABLES[nom]) + "]",
                _tronquer(nom, getattr(t, nom), av),
            )
        )

    # --- ce que buildFromMatrices ne prend pas : en commentaire ---
    if t.T_Shunts:
        av.append(
            f"T_Shunts : {len(t.T_Shunts)} shunt(s) omis. buildFromMatrices "
            "n'a pas d'argument pour eux ; le bloc commente en fin de fichier "
            "les pose noeud par noeud, comme le fait WSCC3SM_CIGRE3.m."
        )
    if t.T_Loads:
        av.append(
            f"T_Loads : {len(t.T_Loads)} charge(s) speciale(s) (I ou P) omises. "
            "Cote MATLAB, T_Loads porte les charges a impedance constante ; "
            "les charges a courant ou a puissance constante n'ont pas "
            "d'equivalent dans buildFromMatrices."
        )
    for nom_reglage in ("SM_PARAMS", "GFM_PARAMS", "GFL_PARAMS"):
        d = getattr(t, nom_reglage)
        if d:
            av.append(
                f"{nom_reglage} : {len(d)} reglage(s) de famille omis "
                f"({', '.join(sorted(d))}). Le fichier produit utilise les "
                "valeurs par defaut des classes MATLAB."
            )
    if t.PARAMS_PAR_BUS:
        av.append(
            f"PARAMS_PAR_BUS : reglages de {len(t.PARAMS_PAR_BUS)} bus omis. "
            "buildFromMatrices ne sait pas regler un equipement seul."
        )

    if t.T_Shunts or t.T_Loads:
        L.append("")
        L.append(
            "    % --- NON REPRIS PAR buildFromMatrices ---------------------------"
        )
    if t.T_Shunts:
        L.append("    % Shunts [BusID, B_pu] : a poser apres construction, comme dans")
        L.append("    % WSCC3SM_CIGRE3.m :")
        L.append("    %   for s = 1:size(T_Shunts, 1)")
        L.append("    %       net.Nodes(T_Shunts(s, 1)).addShunt(T_Shunts(s, 2));")
        L.append("    %   end")
        L.append("    % T_Shunts = [")
        for r in t.T_Shunts:
            L.append(f"    %     {_nombre(int(r[1]))}  {_nombre(r[2])} ;")
        L.append("    % ];")
    if t.T_Loads:
        L.append("    % Charges speciales [ID, BusID, P_W, Q_var, Type] :")
        for r in t.T_Loads:
            vals = "  ".join(_nombre(v) for v in r)
            L.append(f"    %     {vals} ;")

    L += [
        "",
        "    %% ===================== CONSTRUCTION (boilerplate) ==================",
        "    net = Network(Parameters.fb, Parameters.Sb);",
        "    net.buildFromMatrices(Y_DER, T_Lines_SI, T_Lines_pu, T_Transfos, ...",
        "        T_IBs, T_SMs, T_GFM1s, T_GFM2s, T_GFM3s, T_GFM4s, T_GFL1s, "
        "T_GFL2s, T_Loads);",
        "",
        "    [R_Nodes, conv] = net.Load_flow(Y_DER, false);",
        "    if ~conv",
        f"        warning('{fn}:loadFlow', 'Load flow non converge.');",
        "    end",
        "end",
    ]
    return "\n".join(L) + "\n", av


def exporter(t: Tables) -> str:
    """Les tables, au format de la cellule de donnees des notebooks.

    Ce qui est dessine dans le constructeur redevient ainsi un reseau de
    test : meme format, meme ordre de colonnes, directement collable dans
    ``g2egrid/reseaux_test/``.
    """

    def val(v) -> str:
        if isinstance(v, str):
            # repr(), pas f"'{v}'" : une chaine avec une apostrophe ou un saut
            # de ligne ajoutait du code au script que l'on va executer
            return repr(v)
        if isinstance(v, (int, np.integer)) or (
            isinstance(v, float) and v.is_integer() and abs(v) < 1e6
        ):
            return str(int(v))
        return f"{v:.6g}"

    lignes = [
        f"NOM_RESEAU = {json.dumps(str(t.nom), ensure_ascii=False)}",
        "F_ZOOM = 5.0",
        "",
        f"f0 = {t.f0:g}",
        f"Sb = {t.Sb:g}",
        "",
        "# BusType : 1 = slack, 2 = PV, 3 = PQ",
    ]
    # le traitement des noeuds sans capacite, quand ce n'est pas le defaut :
    # les reseaux de parite MATLAB gardent la resistance de Node.m
    if t.R_NOEUDS != R_CAUSALITE:
        lignes[5:5] = [f"R_NOEUDS = {t.R_NOEUDS:g}"]
    if t.NOEUDS != NOEUDS_DEFAUT:
        lignes[5:5] = [f'NOEUDS = "{t.NOEUDS}"']
    for nom, colonnes in TABLES.items():
        rows = getattr(t, nom)
        lignes.append("")
        lignes.append(f"# {nom} : [{', '.join(colonnes)}]")
        if not rows:
            lignes.append(f"{nom} = []")
            continue
        lignes.append(f"{nom} = [")
        for r in rows:
            lignes.append("    (" + ", ".join(val(v) for v in r) + "),")
        lignes.append("]")
    for nom, d in (
        ("SM_PARAMS", t.SM_PARAMS),
        ("GFM_PARAMS", t.GFM_PARAMS),
        ("GFL_PARAMS", t.GFL_PARAMS),
    ):
        lignes.append("")
        lignes.append(f"{nom} = {dict(d)!r}")
    if t.PARAMS_PAR_BUS:
        lignes.append("")
        lignes.append(f"PARAMS_PAR_BUS = {dict(t.PARAMS_PAR_BUS)!r}")
    return "\n".join(lignes) + "\n"
