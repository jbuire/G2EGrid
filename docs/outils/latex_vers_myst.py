"""Convertit les chapitres LaTeX de ``documentation/`` en pages MyST pour Sphinx.

LA SOURCE RESTE LE LATEX. Les pages de ``docs/manuel/`` en sont TIREES par ce
script ; on les regenere apres chaque modification d'un chapitre, et on ne les
edite pas a la main (l'edition serait perdue a la conversion suivante) ::

    uvx --with pypandoc-binary --with pymupdf python docs/outils/latex_vers_myst.py

(a lancer depuis ``g2egrid/``). pandoc fait l'essentiel -- equations, tableaux,
listings, macros a arguments comme ``\\cartemodes`` --, et ce script ne fait que
le passage de son Markdown a celui de MyST :

* titres ``{#sec:x}`` -> cibles ``(sec-x)=`` ; renvois ``\\ref`` -> ``{ref}``,
  ``{numref}`` ou ``{eq}`` selon ce qu'ils designent ;
* equations numerotees -> blocs ``$$ ... $$ (eq-x)``, ``align`` -> ``aligned``
  (MathJax refuse un ``align`` dans un bloc d'affichage) ;
* encadres ``remarque`` / ``attention`` -> admonitions ;
* figures PDF -> PNG (PyMuPDF), rangees dans ``docs/_static/figures/`` ;
* la seule figure TikZ indispensable, la chaine de calcul, devient un schema
  Mermaid ; les autres sont retirees, la version PDF les garde.

Les deux-points des etiquettes LaTeX deviennent des tirets : ``sec:element``
se cite ``{ref}`sec-element``` depuis n'importe quelle page Markdown.
"""

from __future__ import annotations

import pathlib
import re

import pypandoc

ICI = pathlib.Path(__file__).resolve().parent
DOCS = ICI.parent
RACINE = DOCS.parent.parent
LATEX = RACINE / "documentation"
SORTIE = DOCS / "manuel"
FIGURES = DOCS / "_static" / "figures"
CHEMINS_FIGURES = (LATEX / "figures", RACINE / "tool_documentation" / "Figures")

#: (fichier, etiquette du chapitre dans main.tex, titre)
CHAPITRES = [
    ("01_introduction", "ch:intro", "Introduction"),
    ("02_principes", "ch:principes", "Principes de la modélisation"),
    ("03_modeles", "ch:modeles", "Bibliothèque de modèles"),
    ("04_v2", "ch:v2", "Architecture de G2Elin V2 (MATLAB)"),
    ("05_python", "ch:python", "Architecture du portage Python g2egrid"),
    ("06_usage", "ch:usage", "Mode d'emploi"),
    ("07_exemples", "ch:exemples", "Réseaux de test"),
    ("08_sources", "ch:sources", "Source de tension, source de courant"),
    ("09_transport_distribution", "ch:td", "Séparation transport / distribution"),
    ("10_transition", "ch:transition", "La transition, mode par mode"),
    ("11_validation", "ch:validation", "Validation contre G2Elin V1"),
    ("12_atelier", "ch:atelier", "L'atelier : deux voies, AC/DC, zones, interface"),
    ("A_per_unit", "ann:pu", "Annexe A — Système per-unit"),
    ("B_reglage", "ann:reglage", "Annexe B — Réglage automatique des correcteurs"),
]

#: macros du preambule de main.tex ; les encadres deviennent des
#: environnements que pandoc laisse passer en ``::: {.x}``
PREAMBULE = r"""
\newcommand{\code}[1]{\texttt{#1}}
\newcommand{\wb}{\omega_b}
\newcommand{\dt}[1]{\frac{\mathrm{d}#1}{\mathrm{d}t}}
\newcommand{\V}{V1}
\newcommand{\VD}{V2}
\newcommand{\up}[1]{\textsuperscript{#1}}
\newenvironment{remarque}[1][]{\begin{remarqueX}\textbf{#1}\ }{\end{remarqueX}}
\newenvironment{attention}[1][]{\begin{attentionX}\textbf{#1}\ }{\end{attentionX}}
"""

CHAINE_MERMAID = """```{mermaid}
flowchart LR
    tab["Tables<br/>Y_DER, T_*"] --> lf["Load flow<br/>V, S"]
    lf --> eq["Équilibre par élément<br/>x0, z0, u0"]
    eq --> lin["Jacobiennes<br/>par élément"]
    lin --> lft["Câblage LFT et<br/>élimination des algébriques"]
    lft --> ss["A, B, C, D"]
    ss --> modal["Analyse modale"]
```
"""


def etiquette(lbl: str) -> str:
    return lbl.strip().replace(":", "-").replace("_", "-")


# ------------------------------------------------------------------ figures
def figure_png(nom: str) -> str | None:
    """Convertit ``nom`` (PDF) en PNG dans ``_static/figures`` ; rend le chemin
    relatif a ``manuel/``."""
    import fitz  # PyMuPDF

    base = pathlib.Path(nom).name
    if not base.endswith(".pdf"):
        base += ".pdf"
    src = next((d / base for d in CHEMINS_FIGURES if (d / base).exists()), None)
    if src is None:
        return None
    FIGURES.mkdir(parents=True, exist_ok=True)
    dst = FIGURES / (pathlib.Path(base).stem + ".png")
    with fitz.open(src) as doc:
        doc[0].get_pixmap(dpi=130).save(dst)
    return f"../_static/figures/{dst.name}"


RE_FIGURE = re.compile(r"<figure(?P<attrs>[^>]*)>(?P<corps>.*?)</figure>", re.S)


def figures(md: str) -> str:
    def rempl(m):
        attrs, corps = m.group("attrs"), m.group("corps")
        ident = re.search(r'id="([^"]+)"', attrs)
        cap = re.search(r"<figcaption>(.*?)</figcaption>", corps, re.S)
        img = re.search(r'<(?:img|embed) src="([^"]+)"', corps)
        legende = cap.group(1).strip() if cap else ""
        if ident and ident.group(1) == "fig:chaine":
            return CHAINE_MERMAID + f"\n*{legende}*\n"
        if not img:
            return ""  # TikZ : gardee dans la version PDF seulement
        png = figure_png(img.group(1))
        if png is None:
            return ""
        nom = f":name: {etiquette(ident.group(1))}\n" if ident else ""
        return f"\n```{{figure}} {png}\n{nom}:width: 100%\n\n{legende}\n```\n"

    return RE_FIGURE.sub(rempl, md)


# ------------------------------------------------------------------- maths
def _equation(corps: str) -> str:
    lbl = re.search(r"\\label\{([^}]+)\}", corps)
    corps = re.sub(r"\\label\{[^}]+\}", "", corps)
    corps = re.sub(r"\\nonumber|\\notag", "", corps)
    corps = re.sub(r"\\begin\{equation\*?\}|\\end\{equation\*?\}", "", corps)
    for env, rempl in (("align", "aligned"), ("gather", "gathered")):
        corps = re.sub(rf"\\begin\{{{env}\*?\}}", rf"\\begin{{{rempl}}}", corps)
        corps = re.sub(rf"\\end\{{{env}\*?\}}", rf"\\end{{{rempl}}}", corps)
    corps = corps.strip()
    suffixe = f" ({etiquette(lbl.group(1))})" if lbl else ""
    return f"\n\n$$\n{corps}\n$${suffixe}\n\n"


def maths(md: str) -> str:
    return re.sub(r"\$\$(.+?)\$\$", lambda m: _equation(m.group(1)), md, flags=re.S)


# ---------------------------------------------------------------- renvois
RE_RENVOI = re.compile(
    r"\[[^\]]*?(?:\\\[[^\]]*\\\][^\]]*)?\]\(#([^)]+)\)"
    r"\{reference-type=\"(\w+)\" reference=\"[^\"]+\"\}"
)


def renvois(md: str) -> str:
    def rempl(m):
        cible, genre = m.group(1), m.group(2)
        e = etiquette(cible)
        if genre == "eqref" or cible.startswith("eq:"):
            return f"{{eq}}`{e}`"
        if cible == "fig:chaine":  # devenue un schema Mermaid, sans numero
            return "ci-dessous"
        if cible.startswith("fig:"):
            return f"{{numref}}`{e}`"
        if cible.startswith("tab:"):
            # un tableau HTML (cellules sur plusieurs lignes) n'est pas numerote :
            # un renvoi a titre explicite vaut pour les deux genres
            return f"@@TAB:{e}@@"  # resolu par main(), une fois tous les titres connus
        return f"{{ref}}`{e}`"

    return RE_RENVOI.sub(rempl, md)


def titres(md: str) -> str:
    return re.sub(
        r"^(#+ .*?) \{#([^}]+)\}\s*$",
        lambda m: f"({etiquette(m.group(2))})=\n{m.group(1)}",
        md,
        flags=re.M,
    )


# ---------------------------------------------------------------- tableaux
#: legende de chaque tableau, par etiquette : les renvois d'un chapitre a
#: l'autre citent le tableau par son titre
LEGENDES: dict[str, str] = {}

RE_TABLEAU = re.compile(r"^::: \{#(tab:[^}]+)\}\n(.*?)^:::\s*$", re.S | re.M)


def tableaux(md: str) -> str:
    """``::: {#tab:x}`` (tableau puis legende) -> directive ``table`` nommee."""

    def rempl(m):
        corps = m.group(2).strip()
        morceaux = corps.rsplit("\n\n", 1)
        if len(morceaux) == 2 and not morceaux[1].lstrip().startswith(("|", "<")):
            tableau, legende = morceaux
        else:
            tableau, legende = corps, ""
        legende = " ".join(legende.split())
        LEGENDES[etiquette(m.group(1))] = legende
        if tableau.lstrip().startswith("<"):  # tableau HTML : pas de directive
            tableau = re.sub(r"<caption>.*?</caption>\s*", "", tableau, flags=re.S)
            return f"({etiquette(m.group(1))})=\n\n*{legende}*\n\n{tableau}\n"
        return f":::{{table}} {legende}\n:name: {etiquette(m.group(1))}\n\n{tableau}\n:::\n"

    return RE_TABLEAU.sub(rempl, md)


# ------------------------------------------------------------- encadres
RE_ENCADRE = re.compile(r"^::: \{\.(remarque|attention)X\}\n(.*?)^:::\s*$", re.S | re.M)


def encadres(md: str) -> str:
    def rempl(m):
        genre, corps = m.group(1), m.group(2).strip()
        titre = re.match(r"\*\*(.+?)\*\*\s*", corps)
        classe = "note" if genre == "remarque" else "warning"
        if titre:
            corps = corps[titre.end() :]
            return f":::{{admonition}} {titre.group(1)}\n:class: {classe}\n\n{corps}\n:::\n"
        return f":::{{{classe}}}\n{corps}\n:::\n"

    return RE_ENCADRE.sub(rempl, md)


def _barres_en_math(ligne: str) -> str:
    """Dans une ligne de tableau, les ``|`` d'une formule casseraient les colonnes."""
    if not ligne.startswith("|"):
        return ligne
    return re.sub(r"\$[^$]*\$", lambda m: m.group(0).replace("|", "\\vert "), ligne)


def nettoyer(md: str) -> str:
    md = md.replace("’", "'")
    md = "\n".join(_barres_en_math(li) for li in md.split("\n"))
    # listings a attributs (``lstlisting`` avec legende) : on garde le langage
    md = re.sub(
        r"^``` \{.*$",  # l'espace distingue ``` {attributs} de ```{directive}
        lambda m: "```matlab" if "matlab" in m.group(0).lower() else "```python",
        md,
        flags=re.M,
    )
    # les ::: {.center} et autres conteneurs de mise en page
    md = re.sub(r"^::: \{?\.?(center|small|footnotesize)\}?\s*$\n", "", md, flags=re.M)
    md = re.sub(r"\n{3,}", "\n\n", md)
    return md.strip() + "\n"


# ------------------------------------------------------------ pre-traitement
def _arg(tex: str, i: int) -> int:
    """Indice qui suit l'argument ``{...}`` commencant a ``tex[i] == '{'``."""
    prof = 0
    for j in range(i, len(tex)):
        if tex[j] == "{":
            prof += 1
        elif tex[j] == "}":
            prof -= 1
            if prof == 0:
                return j + 1
    raise ValueError("accolade non fermee")


def sans_resizebox(tex: str) -> str:
    """``\\resizebox{l}{h}{X}`` -> ``X`` : pandoc perd le tableau qu'il enveloppe."""
    while (i := tex.find("\\resizebox{")) >= 0:
        a = _arg(tex, i + len("\\resizebox"))
        b = _arg(tex, a)
        c = _arg(tex, b)
        tex = tex[:i] + tex[b + 1 : c - 1] + tex[c:]
    return tex


# ----------------------------------------------------------------- pilote
def convertir(fichier: str, lbl: str, titre: str) -> str:
    tex = sans_resizebox(
        (LATEX / "chapitres" / f"{fichier}.tex").read_text(encoding="utf-8")
    )
    source = PREAMBULE + f"\n\\chapter{{{titre}}}\\label{{{lbl}}}\n" + tex
    md = pypandoc.convert_text(
        source,
        "commonmark_x+tex_math_dollars-raw_attribute",
        format="latex",
        extra_args=["--wrap=none", f"--resource-path={LATEX}"],
        cworkdir=str(LATEX),
    )
    for etape in (figures, maths, renvois, titres, tableaux, encadres, nettoyer):
        md = etape(md)
    entete = (
        "<!-- Page GENEREE depuis documentation/chapitres/"
        f"{fichier}.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->\n\n"
    )
    return entete + md


def _titre_court(legende: str) -> str:
    court = re.split(r"[.:]\s", legende + " ", maxsplit=1)[0].strip().rstrip(".")
    return court.replace("`", "")


def main() -> None:
    SORTIE.mkdir(exist_ok=True)
    pages = {f: convertir(f, lbl, titre) for f, lbl, titre in CHAPITRES}
    for fichier, md in pages.items():
        md = re.sub(
            r"@@TAB:([\w-]+)@@",
            lambda m: (
                f"« {{ref}}`{_titre_court(LEGENDES.get(m.group(1), 'ci-dessous'))} <{m.group(1)}>` »"
            ),
            md,
        )
        (SORTIE / f"{fichier}.md").write_text(md, encoding="utf-8")
        print(f"{fichier}.md")


if __name__ == "__main__":
    main()
