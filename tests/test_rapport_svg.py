"""Les figures SVG du rapport sont RECONSTRUITES depuis une liste blanche.

L'ancien assainissement par regex (retirer ``<script>`` et ``onxxx="..."``)
laissait passer ``onload='..'`` (guillemets simples), ``onload=x`` (sans
guillemets), ``<a href="javascript:..">``, ``<foreignObject>`` et un
``<svg/onload=..>``. Ce fichier verifie, par attaque, que rien d'actif ne
survit, et qu'une figure reelle de l'atelier (clipPath, groupes, textes
accentues, ``style="cursor:pointer"``) ressort intacte.
"""

from __future__ import annotations

import re

import pytest

from g2egrid_api.rapport import assainir_svg

ATTAQUES = [
    "<svg xmlns='http://www.w3.org/2000/svg' onload='alert(1)'><circle r='1'/></svg>",
    '<svg onload=alert(1)><circle r="1"/></svg>',
    "<svg/onload=alert(1)><circle r='1'/></svg>",
    '<svg><circle r="1" onclick = "alert(1)" ONMOUSEOVER=\'x\' onfocus=y /></svg>',
    '<svg><a href="javascript:alert(1)"><text>x</text></a></svg>',
    '<svg><a xlink:href="JaVa&#x09;Script:alert(1)"><text>x</text></a></svg>',
    "<svg><script>alert(1)</script><circle r='1'/></svg>",
    "<svg><SCRIPT >alert(1)</SCRIPT ><circle r='1'/></svg>",
    "<svg><foreignObject><iframe srcdoc='<script>alert(3)</script>'></iframe>"
    "</foreignObject><circle r='1'/></svg>",
    "<svg><foreignObject><img src=x onerror=alert(1)></foreignObject><circle r='1'/></svg>",
    '<svg><use href="http://exemple.org/x.svg#a"/><use xlink:href="data:image/svg+xml,..."/></svg>',
    '<svg><image href="http://exemple.org/pixel.png"/></svg>',
    "<svg><style>@import url(http://exemple.org/x.css); circle{fill:red}</style></svg>",
    '<svg><circle r="1" style="fill:url(http://exemple.org/x)"/></svg>',
    '<svg><circle r="1" style="background:expression(alert(1))"/></svg>',
    '<svg><circle r="1" style="background:url(javascript:alert(1))"/></svg>',
    '<svg><circle r="1" fill="url(http://exemple.org/x#a)"/></svg>',
    '<svg><circle r="1" fill="javascript:alert(1)"/></svg>',
    '<svg><animate attributeName="href" values="javascript:alert(1)"/></svg>',
    '<svg><set attributeName="onmouseover" to="alert(1)"/></svg>',
    '<svg><circle r="1"><animate onbegin="alert(1)" attributeName="r"/></circle></svg>',
    "<svg><!-- <script>alert(1)</script> --><![CDATA[<script>alert(1)</script>]]></svg>",
    "<svg><text>&lt;script&gt;alert(1)&lt;/script&gt;</text></svg>",
    '<svg><circle r="1" id="g2egrid-tables"/></svg>',
    '<math><svg><circle r="1" onload="alert(1)"/></svg></math>',
    "<svg><circle r='1' ><script",
]
#: ce qui ne doit jamais subsister dans la sortie
INTERDITS = [
    r"<script",
    r"<foreignobject",
    r"<iframe",
    r"<img",
    r"<style",
    r"<use",
    r"<image",
    r"<animate",
    r"<set\b",
    r"<a[\s>]",
    r"\bon\w+\s*=",
    r"href",
    r"javascript:",
    r"expression\(",
    r"@import",
    r"srcdoc",
    r'id="g2egrid',
]


@pytest.mark.parametrize("attaque", ATTAQUES)
def test_rien_d_actif_ne_survit(attaque):
    sortie = assainir_svg(attaque)
    for motif in INTERDITS:
        assert not re.search(motif, sortie, re.I), (motif, sortie)
    # le texte d'une balise echappee reste du texte
    assert "<script" not in sortie.lower()
    # toute URL restante est une ancre interne
    for u in re.findall(r"url\([^)]*\)", sortie, re.I):
        assert re.fullmatch(r"url\(#[\w-]+\)", u), u


def test_texte_echappe_reste_du_texte():
    sortie = assainir_svg("<svg><text>&lt;b&gt; &amp; x</text></svg>")
    assert "<text>&lt;b&gt; &amp; x</text>" in sortie


def test_hors_svg_et_non_svg_vides():
    assert assainir_svg("<p>bonjour</p>") == ""
    assert assainir_svg("") == ""
    assert "bonjour" not in assainir_svg("<p>bonjour</p><svg><circle r='1'/></svg>")


FIGURE_REELLE = (
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 380" '
    'preserveAspectRatio="xMidYMid meet" style="background:#fff">'
    '<clipPath id="clip3"><rect x="60" y="20" width="900" height="300"/></clipPath>'
    '<g clip-path="url(#clip3)" transform="translate(1.5,2)">'
    '<path d="M0 0 L10 10 C 20 20, 30 30, 40 40" fill="none" stroke="#1f4e79" '
    'stroke-width="1.4" stroke-dasharray="4 2"/>'
    '<polyline points="0,0 5,5 10,0" fill="none" stroke="#c00"/>'
    '<polygon points="0,0 5,5 10,0" fill="#d98a00" fill-opacity="0.25"/>'
    '<line x1="0" y1="0" x2="9" y2="9" stroke="#000" stroke-width="2"/>'
    '<ellipse cx="1" cy="1" rx="2" ry="3"/>'
    '<circle cx="5" cy="5" r="3" fill="#f00" style="cursor:pointer"><title>mode 1 : '
    "f = 2,07 Hz</title></circle>"
    "</g>"
    '<text x="10" y="20" font-size="12" text-anchor="middle" '
    'font-family="system-ui, sans-serif"><tspan dy="1">&#951; = 9,8 % &amp; Re(&lambda;)</tspan>'
    "</text>"
    "</svg>"
)


def test_figure_reelle_intacte():
    """Memes elements, memes attributs, memes valeurs : seul l'ordre des
    attributs d'origine est conserve, aucune perte."""
    sortie = assainir_svg(FIGURE_REELLE)

    def jetons(s: str) -> list[str]:
        return sorted(re.findall(r'<[A-Za-z]+|[\w:-]+="[^"]*"', s))

    attendu = FIGURE_REELLE.replace("&#951;", "η").replace("&lambda;", "λ")
    assert jetons(sortie) == jetons(attendu)
    assert "η = 9,8 % &amp; Re(λ)" in sortie
    assert "<clipPath" in sortie and "viewBox=" in sortie
    assert 'style="cursor:pointer"' in sortie and 'style="background:#fff"' in sortie


def test_rapport_complet_sans_attaque(tmp_path):
    """Par ``rapport_etude`` : la page produite ne contient aucune des
    attaques, et la figure saine y est."""
    from g2egrid.catalogue import charger
    from g2egrid_api import analyse
    from g2egrid_api.rapport import rapport_etude

    figures = {f"attaque {i}": a for i, a in enumerate(ATTAQUES)}
    figures["saine"] = FIGURE_REELLE
    page = rapport_etude(charger("smib"), analyse.construire, figures, n_modes=1)
    corps = page.split("<h2>attaque 0</h2>", 1)[1].split(
        '<script type="application/json"'
    )[0]
    for motif in INTERDITS:
        assert not re.search(motif, corps, re.I), motif
    assert '<clipPath id="clip3">' in page
