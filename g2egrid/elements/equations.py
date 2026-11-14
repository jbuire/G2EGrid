"""Les equations des modeles, ecrites par SymPy a partir du code qui calcule.

POURQUOI -- une page d'equations tapee a la main derive des qu'on touche un
modele : un signe corrige, un limiteur ajoute, et la documentation decrit un
modele qui n'existe plus. Ici, les equations sont celles que
:func:`g2egrid.noyau.symbolic.compile_spec` derive et compile : le MEME appel a
``equations(s)`` de la classe, sur les memes symboles. Ce qui s'affiche est
donc, par construction, ce qui se calcule.

- :func:`equations` : un modele, en LaTeX et en texte, ligne par ligne ;
- :func:`documenter` : tous les modeles de :mod:`g2egrid.elements`, en
  Markdown (``$$ ... $$``) ;
- ``uv run python -m g2egrid.elements.equations > docs/equations_modeles.md``
  regenere la page de reference (un test verifie qu'elle est a jour).
"""

from __future__ import annotations

import inspect
from types import SimpleNamespace

import sympy as sp

from g2egrid.noyau.element import Element

__all__ = ["classes", "documenter", "equations"]


def classes() -> dict[str, type[Element]]:
    """Les classes d'elements exportees par :mod:`g2egrid.elements`, par nom.

    Les alias (``SM``, ``GFM``...) ne sont comptes qu'une fois, sous le nom
    de la classe.
    """
    from g2egrid import elements

    out: dict[str, type[Element]] = {}
    for nom in elements.__all__:
        c = getattr(elements, nom)
        if (
            inspect.isclass(c)
            and issubclass(c, Element)
            and c.__name__ == nom
            and (c.states or c.algebs or c.output_names())
        ):
            out[nom] = c
    return out


def _classe(nom_ou_classe) -> type[Element]:
    if inspect.isclass(nom_ou_classe):
        return nom_ou_classe
    c = classes().get(str(nom_ou_classe))
    if c is None:
        raise ValueError(f"modele {nom_ou_classe!r} inconnu : {sorted(classes())}")
    return c


def _symboles(spec):
    """Les symboles et l'espace de noms, exactement comme ``compile_spec``."""
    ns = SimpleNamespace()
    for noms in (spec.states, spec.algebs, spec.inputs, spec.params):
        for nm in noms:
            setattr(ns, nm, sp.Symbol(nm, real=True))
    return ns


def _latex(e) -> str:
    return sp.latex(e, mul_symbol="dot")


def equations(nom_ou_classe) -> dict:
    """Les equations d'un modele : ``dx/dt = f``, ``0 = g``, ``y = h``.

    Rend un dict JSON-compatible : les noms (parametres, etats, variables
    algebriques, entrees, sorties) et, pour chaque equation, la variable,
    le LaTeX et le texte SymPy de son membre droit.
    """
    c = _classe(nom_ou_classe)
    spec = c.spec()
    ns = _symboles(spec)
    f, g, h = spec.equations(ns)

    def lignes(variables, exprs, gauche):
        return [
            {
                "variable": v,
                "latex": f"{gauche(sp.Symbol(v))} = {_latex(sp.sympify(e))}",
                "texte": sp.sstr(sp.sympify(e)),
            }
            for v, e in zip(variables, exprs)
        ]

    doc = inspect.getdoc(c) or ""
    return {
        "classe": c.__name__,
        "modele": spec.name,
        "resume": doc.split("\n\n", 1)[0].replace("\n", " "),
        "parametres": list(spec.params),
        "etats": list(spec.states),
        "algebriques": list(spec.algebs),
        "entrees": list(spec.inputs),
        "sorties": list(spec.outputs),
        "f": lignes(spec.states, f, lambda s: rf"\frac{{d\,{_latex(s)}}}{{dt}}"),
        "g": lignes(spec.algebs, g, lambda s: "0"),
        "h": lignes(spec.outputs, h, _latex),
    }


def documenter(noms: list[str] | None = None) -> str:
    """Les equations de ``noms`` (tous les modeles par defaut), en Markdown."""
    noms = list(classes()) if noms is None else noms
    out = [
        "# Equations des modeles",
        "",
        "Page generee par `g2egrid.elements.equations` a partir des modeles eux-memes :",
        "ce sont les expressions que SymPy derive et compile, pas une copie.",
        "Pour la regenerer : `uv run python -m g2egrid.elements.equations >",
        "docs/equations_modeles.md`.",
        "",
        "Conventions : `x` les etats (`dx/dt = f`), `z` les variables",
        "algebriques (`0 = g`), `y` les sorties (`y = h`) ; les courants de",
        "port sont injectes dans le bus, en pu base systeme.",
        "",
    ]
    vus: set[str] = set()
    for nom in noms:
        d = equations(nom)
        out += [f"## {d['classe']}", ""]
        if d["modele"] in vus:
            out += [f"Memes equations que `{d['modele']}`.", ""]
            continue
        vus.add(d["modele"])
        if d["resume"]:
            out += [d["resume"], ""]
        for cle, titre in (
            ("parametres", "Parametres"),
            ("etats", "Etats"),
            ("algebriques", "Variables algebriques"),
            ("entrees", "Entrees"),
            ("sorties", "Sorties"),
        ):
            if d[cle]:
                out.append(f"- **{titre}** : " + ", ".join(f"`{n}`" for n in d[cle]))
        out.append("")
        for cle in ("f", "g", "h"):
            for ligne in d[cle]:
                out += ["$$", ligne["latex"], "$$", ""]
    return "\n".join(out)


if __name__ == "__main__":
    import sys

    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    print(documenter(), end="")
