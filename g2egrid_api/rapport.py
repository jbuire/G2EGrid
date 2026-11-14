"""Le rapport d'etude : un fichier HTML autonome, qu'on peut archiver.

Un resultat de stabilite n'a de valeur que si l'on sait DE QUOI il est le
resultat : quel reseau exactement, calcule par quel code. Le rapport porte
donc, en tete :

* la VERSION du code -- le commit git et s'il y avait des modifications non
  enregistrees, a defaut la version du paquet ;
* l'EMPREINTE du reseau -- le SHA-256 de ses tables, ecrites sous forme
  canonique (cles triees) : deux rapports de meme empreinte portent sur le
  meme reseau, au bit pres ;
* les TABLES elles-memes, dans un ``<script type="application/json">`` : le
  bouton Ouvrir de l'atelier relit le rapport comme un reseau, et l'etude se
  refait.

Puis le reseau en chiffres, le load flow, le spectre, le rapport des modes
critiques (:func:`g2egrid.analyse.residus.rapport` : participations, forme,
commandabilite, observabilite, residus) et les figures que l'atelier a
envoyees (le schema, la carte des modes...), en SVG dans la page. Aucune
ressource exterieure : le fichier s'ouvre hors ligne, dans dix ans.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import html
import json
import pathlib
import re
import subprocess
from html.parser import HTMLParser

import numpy as np

from g2egrid.analyse import modal, residus
from g2egrid.reseau.tables import Tables, tables_vers_json

__all__ = [
    "assainir_svg",
    "empreinte",
    "rapport_etude",
    "version_code",
    "MARQUE_TABLES",
]

MARQUE_TABLES = "g2egrid-tables"
_WEB = pathlib.Path(__file__).parent / "web"


def empreinte(t: Tables) -> str:
    """SHA-256 des tables sous forme canonique."""
    brut = json.dumps(tables_vers_json(t), sort_keys=True, default=str)
    return hashlib.sha256(brut.encode()).hexdigest()


def version_code() -> dict:
    """Le commit git du depot (et s'il etait modifie), ou la version du paquet."""
    ici = pathlib.Path(__file__).resolve().parent
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ici,
            capture_output=True,
            text=True,
            timeout=5,
        ).stdout.strip()
        sale = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=ici,
            capture_output=True,
            text=True,
            timeout=5,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        commit, sale = "", ""
    try:
        from importlib.metadata import version

        paquet = version("g2egrid")
    except Exception:  # noqa: BLE001
        paquet = "inconnue"
    return {"commit": commit, "modifie": bool(sale), "paquet": paquet}


# ------------------------------------------------------------------ figures SVG
# Une figure vient du navigateur (la page serialise son propre SVG), donc d'un
# client qu'on ne croit pas : le rapport est ouvert plus tard, hors du jeton et
# de l'atelier, par quelqu'un qui lui fait confiance. Une regex qui « retire
# <script> et onclick="..." » laisse passer ``onload='..'``, ``onload=x``,
# ``<a href="javascript:..">``, ``<foreignObject>`` et un ``<svg/onload=..>`` :
# on ne FILTRE donc pas le texte, on le RECONSTRUIT. Le SVG est lu en jetons
# (``html.parser``), et seuls les elements et attributs d'une liste blanche
# sont reecrits, valeurs echappees ; tout le reste disparait.

#: elements de dessin admis, en minuscules -> graphie SVG
_SVG_ELEMENTS = {
    n.lower(): n
    for n in (
        "svg g defs title desc path text tspan line polyline polygon circle "
        "rect ellipse clipPath linearGradient radialGradient stop marker"
    ).split()
}
#: elements dont on garde le contenu mais pas la balise (un lien : son texte)
_SVG_DEROULES = {"a", "span", "div", "p", "b", "i", "em", "strong"}
#: attributs de PRESENTATION et de geometrie, en minuscules -> graphie SVG
_SVG_ATTRIBUTS = {
    n.lower(): n
    for n in (
        "x y x1 x2 y1 y2 cx cy r rx ry width height d points transform dx dy "
        "rotate textLength lengthAdjust fill fill-opacity fill-rule stroke "
        "stroke-width stroke-opacity stroke-dasharray stroke-dashoffset "
        "stroke-linecap stroke-linejoin stroke-miterlimit opacity "
        "font-family font-size font-weight font-style text-anchor "
        "dominant-baseline alignment-baseline letter-spacing text-decoration "
        "class id viewBox preserveAspectRatio xmlns version clip-path "
        "clipPathUnits visibility display offset stop-color stop-opacity "
        "gradientUnits gradientTransform markerWidth markerHeight refX refY "
        "orient markerUnits marker-start marker-end role aria-label "
        "pointer-events vector-effect shape-rendering"
    ).split()
}
#: proprietes CSS admises dans ``style`` (presentation seulement)
_CSS_PROPRIETES = {
    "cursor", "background", "background-color", "fill", "fill-opacity",
    "stroke", "stroke-width", "stroke-opacity", "stroke-dasharray", "opacity",
    "font-family", "font-size", "font-weight", "font-style", "text-anchor",
    "pointer-events", "display", "visibility",
}  # fmt: skip
_URL_INTERNE = re.compile(r"url\(\s*['\"]?#[A-Za-z_][\w.:-]*['\"]?\s*\)", re.I)
_VALEUR_SURE = re.compile(r"[^\x00-\x08\x0b\x0c\x0e-\x1f\x7f\\]*")
#: elements HTML sans balise fermante : ils ne changent pas la profondeur
_VIDES_HTML = {
    "img", "br", "hr", "input", "meta", "link", "area", "base", "col",
    "embed", "source", "track", "wbr", "param",
}  # fmt: skip


def _valeur_sure(v: str) -> bool:
    """Une valeur d'attribut sans schema actif ni ressource externe.

    On compare sur la valeur DEBARRASSEE des espaces et des caracteres de
    controle (``java\\tscript:``), et on n'admet ``url(...)`` que vers une
    ancre interne (``url(#id)``)."""
    if not _VALEUR_SURE.fullmatch(v):
        return False
    compact = re.sub(r"[\s\x00-\x20]+", "", v).lower()
    if any(m in compact for m in ("javascript:", "vbscript:", "data:", "expression(")):
        return False
    if "@import" in compact or "&#" in compact or "&{" in compact:
        return False
    sans_url = _URL_INTERNE.sub("", v)
    return "url(" not in re.sub(r"\s+", "", sans_url).lower()


def _style_sur(css: str) -> str:
    """Les seules declarations de presentation de ``css`` dont la valeur est sure."""
    gardees = []
    for decl in css.split(";"):
        prop, deux_points, val = decl.partition(":")
        prop, val = prop.strip().lower(), val.strip()
        if deux_points and prop in _CSS_PROPRIETES and val and _valeur_sure(val):
            gardees.append(f"{prop}:{val}")
    return ";".join(gardees)


class _AssainisseurSVG(HTMLParser):
    """Reecrit un SVG en ne gardant que la liste blanche (voir plus haut)."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.sortie: list[str] = []
        self.pile: list[str] = []  # elements SVG ouverts (graphie SVG)
        self.profondeur_svg = 0  # > 0 : on est dans un <svg>
        self.ignore = 0  # > 0 : dans un sous-arbre elimine

    def _attributs(self, attrs) -> str:
        sortie = []
        for nom, val in attrs:
            cle = _SVG_ATTRIBUTS.get((nom or "").lower())
            val = "" if val is None else val
            if cle is None and (nom or "").lower() != "style":
                continue  # on*, href, xlink:href, src, srcdoc, ...
            if (nom or "").lower() == "style":
                val = _style_sur(val)
                if not val:
                    continue
                cle = "style"
            elif not _valeur_sure(val) or (
                cle == "id"
                and not re.fullmatch(r"(?!g2egrid)[A-Za-z][\w-]{0,63}", val, re.I)
            ):
                continue  # un id « g2egrid... » masquerait les donnees du rapport
            sortie.append(f' {cle}="{html.escape(val, quote=True)}"')
        return "".join(sortie)

    def handle_starttag(self, tag, attrs, vide=False):
        tag = tag.lower()
        if self.ignore:
            if not vide and tag not in _VIDES_HTML:
                self.ignore += 1
            return
        if tag in _SVG_DEROULES:
            return
        nom = _SVG_ELEMENTS.get(tag)
        if nom is None or (nom != "svg" and not self.profondeur_svg):
            # element hors liste (script, foreignObject, style, use, image,
            # animate, iframe...) : lui et ses descendants disparaissent
            if not vide and tag not in _VIDES_HTML:
                self.ignore = 1
            return
        attrs_sures = self._attributs(attrs)
        if nom == "svg" and 'xmlns="' not in attrs_sures:
            attrs_sures += ' xmlns="http://www.w3.org/2000/svg"'
        self.sortie.append(f"<{nom}{attrs_sures}{'/' if vide else ''}>")
        if not vide:
            self.pile.append(nom)
            if nom == "svg":
                self.profondeur_svg += 1

    def handle_startendtag(self, tag, attrs):
        self.handle_starttag(tag, attrs, vide=True)

    def handle_endtag(self, tag):
        tag = tag.lower()
        if self.ignore:
            if tag not in _VIDES_HTML:
                self.ignore -= 1
            return
        nom = _SVG_ELEMENTS.get(tag)
        if nom is None or nom not in self.pile:
            return
        while self.pile:  # referme jusqu'a l'element vise
            dernier = self.pile.pop()
            self.sortie.append(f"</{dernier}>")
            if dernier == "svg":
                self.profondeur_svg -= 1
            if dernier == nom:
                break

    def handle_data(self, data):
        if self.profondeur_svg and not self.ignore:
            self.sortie.append(html.escape(data, quote=False))

    # commentaires, declarations, instructions de traitement, CDATA : rien
    def handle_comment(self, data):
        pass

    def handle_decl(self, decl):
        pass

    def handle_pi(self, data):
        pass

    def unknown_decl(self, data):
        pass

    def fin(self) -> str:
        self.close()
        while self.pile:
            self.sortie.append(f"</{self.pile.pop()}>")
        return "".join(self.sortie)


def assainir_svg(svg: str) -> str:
    """Le SVG ``svg`` reduit a des elements et attributs de dessin surs.

    Rien n'est execute ni charge : pas de script, de gestionnaire ``on*``, de
    ``javascript:``, de ``foreignObject``, de lien, de ``use`` ou d'image
    externes, de ``style`` avec ``url(`` hors ancre interne. Rend la chaine
    vide si le texte ne contient aucun ``<svg>``."""
    p = _AssainisseurSVG()
    p.feed(svg)
    return p.fin()


def _md_vers_html(texte: str) -> str:
    """Le petit Markdown de :meth:`RapportMode.texte` : titres, paragraphes,
    tableaux."""
    out, table = [], []

    def vider():
        if not table:
            return
        lignes = [r for r in table if not re.fullmatch(r"\|[-| ]+\|", r)]
        cellules = [[c.strip() for c in r.strip("|").split("|")] for r in lignes]
        tete, *corps = cellules
        out.append(
            "<table><tr>"
            + "".join(f"<th>{html.escape(c)}</th>" for c in tete)
            + "</tr>"
            + "".join(
                "<tr>" + "".join(f"<td>{html.escape(c)}</td>" for c in r) + "</tr>"
                for r in corps
            )
            + "</table>"
        )
        table.clear()

    for ligne in texte.splitlines():
        if ligne.startswith("|"):
            table.append(ligne)
            continue
        vider()
        if ligne.startswith("### "):
            out.append(f"<h3>{html.escape(ligne[4:])}</h3>")
        elif ligne.strip():
            out.append(f"<p>{html.escape(ligne)}</p>")
    vider()
    return "\n".join(out)


def _tableau(entetes, lignes) -> str:
    return (
        "<table><tr>"
        + "".join(f"<th>{html.escape(str(e))}</th>" for e in entetes)
        + "</tr>"
        + "".join(
            "<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in r) + "</tr>"
            for r in lignes
        )
        + "</table>"
    )


def _f(x, n=4) -> str:
    return f"{x:.{n}g}" if isinstance(x, (int, float, np.floating)) else str(x)


def rapport_etude(
    t: Tables,
    construire,
    figures: dict[str, str] | None = None,
    n_modes: int = 3,
    titre: str | None = None,
    positions: dict | None = None,
) -> str:
    """Le rapport HTML du reseau ``t``. ``construire`` : la fonction qui
    assemble et resout (celle de l'atelier, avec son cache). ``figures`` :
    ``{legende: svg}`` a inclure."""
    from g2egrid.noyau.progres import signaler

    v = version_code()
    emp = empreinte(t)
    nom = t.nom or "reseau"
    quand = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    sections = []

    # --- le reseau ----------------------------------------------------------
    signaler(0.05, "reseau")
    comptes = [
        (k, len(getattr(t, k)))
        for k in (
            "Y_DER", "T_Lines_SI", "T_Lines_pu", "T_Transfos", "T_IBs", "T_SMs",
            "T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s", "T_GFL1s", "T_GFL2s",
            "T_Loads", "T_Shunts", "T_BusDC", "T_VSCs",
        )
        if getattr(t, k, None)
    ]  # fmt: skip
    sections.append(
        "<h2>Le reseau</h2>"
        + _tableau(["table", "lignes"], comptes)
        + f"<p>f0 = {t.f0:g} Hz, Sb = {t.Sb / 1e6:g} MVA.</p>"
    )

    # --- le load flow ---------------------------------------------------------
    signaler(0.15, "load flow")
    a = construire(t)
    lf, net = a.lf, a.net
    idx = net._bus_index()
    types = {int(r[0]): int(r[1]) for r in t.Y_DER}
    lignes_lf = [
        (
            b,
            {1: "slack", 2: "PV", 3: "PQ"}.get(types.get(int(b), 3), "?"),
            _f(abs(lf.V[i]), 5),
            _f(float(np.degrees(np.angle(lf.V[i]))), 4),
            _f(lf.S[i].real, 4),
            _f(lf.S[i].imag, 4),
        )
        for b, i in sorted(idx.items())
    ]
    sections.append(
        "<h2>Point de fonctionnement</h2>"
        + f"<p>{lf.iterations} iterations, ecart {lf.mismatch:.1e} ; pertes "
        f"{lf.S.real.sum():.4g} pu.</p>"
        + _tableau(
            ["bus", "type", "|V| [pu]", "angle [deg]", "P [pu]", "Q [pu]"], lignes_lf
        )
    )

    # --- le spectre -----------------------------------------------------------
    signaler(0.35, "spectre")
    ss = a.etat_stationnaire()
    lam = ss.eig()
    jauge = np.abs(lam) < 1e-6 * max(1.0, float(np.abs(lam).max()))
    instables = int(np.sum((lam.real > 1e-6) & ~jauge))
    liste = modal.modes(ss, f_max=None, n_etats=3, oscillants=True)[:12]
    sections.append(
        "<h2>Spectre</h2>" + f"<p>{len(ss.states)} etats ; plus grande partie reelle "
        f"{float(lam.real[~jauge].max()):.4g} 1/s ; "
        + (
            f"<strong class='instable'>{instables} mode(s) instable(s)</strong>"
            if instables
            else "<strong class='stable'>stable</strong>"
        )
        + ".</p>"
        + _tableau(
            ["f [Hz]", "zeta [%]", "Re [1/s]", "etats dominants"],
            [
                (
                    _f(m.f, 4),
                    _f(m.zeta, 3),
                    _f(m.valeur.real, 4),
                    ", ".join(f"{n} {p:.0f} %" for n, p in m.etats),
                )
                for m in liste
            ],
        )
    )

    # --- les modes critiques ----------------------------------------------------
    signaler(0.6, "modes critiques")
    rap = residus.rapport(ss, n_modes=n_modes, f_min=0.01)
    sections.append(
        "<h2>Modes critiques</h2><p>Les moins amortis (mode de jauge exclu) : "
        "elements et etats qui y participent, forme modale relative a l'etat "
        "dominant, entrees qui les excitent, sorties qui les voient, residus "
        "les plus forts et coefficients du modele qui les deplacent.</p>"
        + "".join(
            f"<section class='mode'>{_md_vers_html(m.texte())}</section>" for m in rap
        )
    )

    # --- les figures ----------------------------------------------------------
    signaler(0.9, "figures")
    for legende, svg in (figures or {}).items():
        if not svg or "<svg" not in svg[:200]:
            continue
        svg = assainir_svg(svg)
        if not svg:
            continue
        sections.append(f"<h2>{html.escape(legende)}</h2><figure>{svg}</figure>")

    style_atelier = (
        (_WEB / "style.css").read_text(encoding="utf-8") if _WEB.exists() else ""
    )
    donnees = json.dumps(
        {
            "nom": nom,
            "tables": tables_vers_json(t),
            "positions": positions or {},
            "empreinte": emp,
        },
        default=str,
    ).replace("</", "<\\/")
    commit = (
        v["commit"][:12] + (" (modifie)" if v["modifie"] else "")
        if v["commit"]
        else "hors depot git"
    )
    return f"""<!doctype html>
<html lang="fr"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(titre or "Rapport g2egrid : " + nom)}</title>
<style>
{style_atelier}
body {{ max-width: 1100px; margin: 0 auto; padding: 24px; background: #fff; color: #1c2330;
  font-family: system-ui, sans-serif; }}
.cartouche {{ border: 1px solid #d0d5dd; border-radius: 8px; padding: 12px 16px; background: #f8f9fb; }}
.cartouche code {{ font-size: 12px; word-break: break-all; }}
table {{ border-collapse: collapse; margin: 8px 0 16px; font-size: 13px; }}
th, td {{ border: 1px solid #d0d5dd; padding: 3px 8px; text-align: left; }}
th {{ background: #eef1f5; }}
section.mode {{ border-left: 3px solid #1f4e79; padding-left: 12px; margin: 12px 0; }}
.instable {{ color: #bf2617; }} .stable {{ color: #2a7a4b; }}
figure svg {{ width: 100%; height: auto; border: 1px solid #eef1f5; }}
</style></head><body>
<h1>Rapport d'etude : {html.escape(nom)}</h1>
<div class="cartouche">
<p>Etabli le {quand} par g2egrid {html.escape(v["paquet"])}, code <code>{html.escape(commit)}</code>.</p>
<p>Empreinte du reseau (SHA-256 des tables) : <code>{emp}</code></p>
<p>Ce fichier contient les tables du reseau : le bouton <em>Ouvrir</em> de
l'atelier le relit comme un reseau.</p>
</div>
{"".join(sections)}
<script type="application/json" id="{MARQUE_TABLES}">{donnees}</script>
</body></html>
"""
