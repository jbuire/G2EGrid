"""Les expressions des schemas-blocs sont lues par un analyseur SUR.

``sympy.parse_expr`` evalue le texte en Python : ``__import__('os')`` dans un
bloc ``fonction`` (ou dans un champ ``K``, ``tau``...) s'executait. Ce fichier
verifie que :

* l'arithmetique legitime (celle des exemples de l'editeur) passe et donne la
  meme expression SymPy que ``sympify`` du meme texte ;
* tout le reste (appels hors liste, attributs, indices, lambdas,
  comprehensions, chaines, exposants gigantesques, expressions geantes) est
  refuse avec un message, sans rien executer.
"""

from __future__ import annotations

import pytest
import sympy as sp

from g2egrid.blocs import SchemaInvalide, compiler, exemples
from tests.test_blocs import source_filtree


def _schema_expr(expr: str) -> dict:
    s = source_filtree()
    s["blocs"].append({"id": "g", "type": "fonction", "expr": expr, "n": 1})
    s["liens"].append({"de": "fp", "vers": "g"})
    return s


@pytest.mark.parametrize(
    "expr, attendu",
    [
        ("tau*u^2", "tau*y_fp**2"),
        ("u1 + 2.5*tau - u/3", "y_fp + 2.5*tau - y_fp/3"),
        ("-u**2 + (u - 1)*(u + 1)", "-y_fp**2 + (y_fp - 1)*(y_fp + 1)"),
        (
            "sin(u)*cos(2*u) + exp(-u) + sqrt(tau)",
            "sin(y_fp)*cos(2*y_fp) + exp(-y_fp) + sqrt(tau)",
        ),
        (
            "tanh(u) + atan2(u, tau) + log(tau) + tan(u)",
            "tanh(y_fp) + atan2(y_fp, tau) + log(tau) + tan(y_fp)",
        ),
        (
            "Abs(u) + abs(u) + max(u, tau) + min(u, tau) + pi",
            "2*Abs(y_fp) + Max(y_fp, tau) + Min(y_fp, tau) + pi",
        ),
    ],
)
def test_expressions_legitimes_inchangees(expr, attendu):
    c = compiler(_schema_expr(expr))
    j = c.algebriques.index("y_g")
    syms = {n: sp.Symbol(n) for n in ("y_fp", "tau")}
    voulu = sp.sympify(attendu, locals=syms)
    assert sp.simplify(c.g[j] - (sp.Symbol("y_g") - voulu)) == 0


@pytest.mark.parametrize(
    "expr",
    [
        "__import__('os').system('echo pwn')",
        "().__class__.__bases__[0].__subclasses__()",
        "u.__class__",
        "(lambda: 1)()",
        "[x for x in (1, 2)]",
        "{x: 1 for x in (1, 2)}",
        "(x for x in (1, 2))",
        "'a' + 'b'",
        "u if u else tau",
        "u and tau",
        "u < tau",
        "u[0]",
        "open('a')",
        "eval('1')",
        "exec('1')",
        "sin(u, x=1)",
        "sin(*[u])",
        "Symbol('z')",
        "Integer(3)",
        "u; import os",
        "u @ tau",
        "u // tau",
        "u % tau",
        "True",
        "None",
        "1j",
        "(u := 3)",
        "u**1000000",
        "9**9**9",
        "10**(10**6)",
        "sin(" * 60 + "u" + ")" * 60,
        "u+" * 600 + "u",
    ],
)
def test_expressions_malveillantes_refusees(expr, tmp_path, monkeypatch):
    marque = tmp_path / "pwn"
    monkeypatch.setattr("os.system", lambda *a, **k: marque.write_text("x"))
    with pytest.raises(SchemaInvalide):
        compiler(_schema_expr(expr))
    assert not marque.exists()


def test_expression_refusee_dit_pourquoi():
    with pytest.raises(SchemaInvalide, match="appel refuse"):
        compiler(_schema_expr("__import__('os').system('id')"))
    with pytest.raises(SchemaInvalide, match="Attribute refuse"):
        compiler(_schema_expr("u.real"))
    with pytest.raises(SchemaInvalide, match="Lambda refuse"):
        compiler(_schema_expr("(lambda: 1)"))
    with pytest.raises(SchemaInvalide, match="exposant constant"):
        compiler(_schema_expr("u**1000000"))


def test_champ_numerique_passe_par_le_meme_analyseur():
    """``K``, ``tau``... sont des expressions : un champ de bloc non plus ne
    doit pas executer de Python."""
    s = source_filtree()
    s["blocs"][2]["tau"] = "__import__('os').getcwd()"
    with pytest.raises(SchemaInvalide, match="parametres du schema"):
        compiler(s)


def test_exemples_de_l_editeur_compilent():
    """Les cinq schemas de l'editeur (dont le GFL et le GFM complets, avec
    leurs expressions de Park) compilent toujours."""
    for cle, ex in exemples().items():
        c = compiler(ex["schema"])
        assert c.etats or c.algebriques, cle
