"""Configuration Sphinx de la documentation de g2egrid.

Construire la documentation (depuis ``g2egrid/``) ::

    uv run --group docs sphinx-build -b html docs docs/_build/html

puis ouvrir ``docs/_build/html/index.html``.

Trois sortes de pages cohabitent :

* les pages du GUIDE (``prise_en_main/``, ``guide/``), ecrites pour ce site ;
* le MANUEL (``manuel/``), tire des chapitres LaTeX de ``documentation/`` par
  ``docs/outils/latex_vers_myst.py`` -- la source reste le LaTeX ;
* la REFERENCE (``api/``), produite par autodoc a partir des docstrings.
"""

from __future__ import annotations

import os
import pathlib
import shutil
import sys

ICI = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ICI.parent))  # g2egrid et g2egrid_api, sans installation

import g2egrid

# -- Notebooks ------------------------------------------------------------------
# Les notebooks vivent hors de la documentation (reseaux_test/, etudes/) : on
# les copie dans _notebooks/ a chaque construction, et myst-nb les rend AVEC
# leurs sorties enregistrees, sans les executer.
NOTEBOOKS = ICI / "_notebooks"
shutil.rmtree(NOTEBOOKS, ignore_errors=True)
for dossier in ("reseaux_test", "etudes"):
    cible = NOTEBOOKS / dossier
    cible.mkdir(parents=True)
    for nb in sorted((ICI.parent / dossier).glob("*.ipynb")):
        # un lien ``../docs/x.md`` est juste depuis reseaux_test/ ou etudes/
        # (GitHub, Jupyter) ; la copie vit dans docs/_notebooks/<dossier>/
        texte = nb.read_text(encoding="utf-8").replace("](../docs/", "](../../")
        (cible / nb.name).write_text(texte, encoding="utf-8")
    if not any(cible.iterdir()):
        # un dossier sans notebook (la version publiee ne distribue pas les
        # etudes) : une page le dit, plutot qu'une section vide et un
        # avertissement de Sphinx
        (cible / "absent.md").write_text(
            f"# {dossier} : non distribue\n\nLes notebooks de `{dossier}/` ne font "
            "pas partie de cette distribution de l'outil ; ils sont dans son "
            "depot de developpement.\n",
            encoding="utf-8",
        )

# -- Projet -------------------------------------------------------------------
project = "G2Egrid"
author = "G2Elab"
copyright = "2026, G2Elab"
release = g2egrid.__version__
version = release
language = "fr"

# -- Extensions ---------------------------------------------------------------
extensions = [
    "myst_nb",  # l'analyseur de myst_parser, plus le rendu des notebooks
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxcontrib.mermaid",
]

source_suffix = {".rst": "restructuredtext", ".md": "myst-nb", ".ipynb": "myst-nb"}
nb_execution_mode = "off"
nb_merge_streams = True
root_doc = "index"
exclude_patterns = ["_build", "outils", "Thumbs.db", ".DS_Store"]
templates_path = ["_templates"]

myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
    "amsmath",
    "deflist",
    "attrs_inline",
]
myst_heading_anchors = 3
myst_dmath_double_inline = True
numfig = True
numfig_format = {
    "figure": "Figure %s",
    "table": "Tableau %s",
    "code-block": "Listing %s",
}

# -- autodoc ------------------------------------------------------------------
autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autodoc_member_order = "bysource"
# les docstrings du paquet sont en prose, pas au format NumPy/Google : napoleon
# ne sert qu'aux quelques sections "Parameters" qu'on y croise
napoleon_numpy_docstring = True
napoleon_google_docstring = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "sympy": ("https://docs.sympy.org/latest", None),
}

# -- Mathematiques et schemas ----------------------------------------------------
# MathJax et Mermaid viennent d'un CDN, sauf s'ils ont ete recuperes localement
# par docs/outils/vendoriser.py : la documentation se lit alors hors ligne.
VENDOR = ICI / "_static" / "vendor"
if (VENDOR / "tex-svg-full.js").exists():
    mathjax_path = "vendor/tex-svg-full.js"
if (VENDOR / "mermaid" / "mermaid.esm.min.mjs").exists():
    mermaid_use_local = "vendor/mermaid/mermaid.esm.min.mjs"
    mermaid_include_elk = False

# les macros du preambule LaTeX de documentation/main.tex
mathjax3_config = {
    "tex": {
        "macros": {
            "bm": ["\\boldsymbol{#1}", 1],
            "wb": "\\omega_b",
            "dt": ["\\frac{\\mathrm{d}#1}{\\mathrm{d}t}", 1],
        }
    }
}

# -- Sortie HTML --------------------------------------------------------------
# Read the Docs par defaut, comme VeraGrid ; G2EGRID_DOCS_THEME=furo pour l'autre
html_theme = os.environ.get("G2EGRID_DOCS_THEME", "sphinx_rtd_theme")
html_title = f"G2Egrid {release}"
html_short_title = "G2Egrid"
html_static_path = ["_static"]
html_css_files = ["g2egrid.css"]
html_logo = "_static/logo.svg"
html_favicon = "_static/logo.svg"
if html_theme == "sphinx_rtd_theme":
    html_theme_options = {
        "collapse_navigation": False,
        "navigation_depth": 3,
        "sticky_navigation": True,
        "includehidden": True,
        "titles_only": False,
        "style_external_links": True,
    }
    html_context = {
        # "Edit on GitHub" sur chaque page
        "display_github": True,
        "github_user": "jbuire",
        "github_repo": "g2elin_v2",
        "github_version": "master",
        "conf_py_path": "/g2egrid/docs/",
    }
html_show_sourcelink = False

# le manuel converti contient des renvois vers des tableaux sans etiquette et
# des ancres de la version PDF : on ne veut pas qu'ils masquent le reste
suppress_warnings = ["myst.header", "mystnb.unknown_mime_type"]


# -- Reexportations -------------------------------------------------------------
# ``g2egrid/__init__.py`` reexporte ``Network``, ``Element``, ``InitContext``...
# dans son ``__all__`` : la page du paquet les documentait une seconde fois,
# et toute annotation ``net: Network`` visait alors deux cibles
# (``g2egrid.Network`` et ``g2egrid.reseau.network.Network``). Chaque objet n'est
# documente que dans son module d'origine.
def _sauter_reexportations(app, what, name, obj, skip, options):
    module = getattr(obj, "__module__", None)
    if what == "module" and name and module and module.startswith("g2egrid"):
        parent = getattr(app.env, "temp_data", {}).get("autodoc:module")
        # tout PAQUET (g2egrid, g2egrid.elements, g2egrid_api.analyse...) : un nom
        # reexporte d'un sous-module n'est documente que la-bas. Sans cela,
        # g2egrid.elements.turbo_alternateur (sous-module ET fonction
        # reexportee) etait decrit deux fois.
        est_paquet = hasattr(sys.modules.get(parent or ""), "__path__")
        if parent and est_paquet and module != parent:
            return True
    return None


def setup(app):
    app.connect("autodoc-skip-member", _sauter_reexportations)
