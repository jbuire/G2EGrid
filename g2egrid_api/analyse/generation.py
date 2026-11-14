"""Le code Python d'un reseau, a copier dans un script.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import g2egrid.reseau.catalogue
from g2egrid.formats.export import exporter
from g2egrid.analyse.comparaison import comparer as comparer_versions
from g2egrid.reseau.tables import TABLES, Tables
from g2egrid.reseau.zones import aplatir
from .fichiers_reseaux import _reseau_depuis


# --------------------------------------------------------------------- code
def _tables_de(d: dict) -> Tables:
    """Les tables A PLAT d'une version : un nom du catalogue
    (``{"catalogue": "wscc9"}``) ou n'importe quelle forme de reseau."""
    if d.get("catalogue"):
        return g2egrid.reseau.catalogue.charger(str(d["catalogue"]))
    return aplatir(_reseau_depuis(d))


def comparer(a: dict, b: dict, f_max: float | None = None) -> dict:
    """Deux versions d'un reseau, cote a cote. Voir :mod:`g2egrid.analyse.comparaison`.

    ``a`` est la reference, ``b`` la version retouchee ; chacune sous l'une
    des formes de :func:`lire_fichier`, ou ``{"catalogue": nom}``. ``f_max``
    ne filtre que l'AFFICHAGE des paires : l'appariement porte toujours sur
    les spectres entiers.
    """
    ta, tb = _tables_de(a), _tables_de(b)
    c = comparer_versions(ta, tb)

    def garder(f):
        return f_max is None or f <= f_max

    paire = lambda m: {  # noqa: E731
        "f_a": m.f_a,
        "zeta_a": m.zeta_a,
        "f_b": m.f_b,
        "zeta_b": m.zeta_b,
        "re_a": m.lam_a.real,
        "im_a": m.lam_a.imag,
        "re_b": m.lam_b.real,
        "im_b": m.lam_b.imag,
        "ecart_f": m.ecart_f,
        "ecart_zeta": m.ecart_zeta,
        "dominant_a": m.dominant_a,
        "dominant_b": m.dominant_b,
        "degrade": m.degrade,
    }
    seul = lambda m: {  # noqa: E731
        "f": m.f,
        "zeta": m.zeta,
        "re": m.lam.real,
        "im": m.lam.imag,
        "dominant": m.dominant,
    }
    brut = lambda v: v if isinstance(v, (int, float, str, type(None))) else str(v)  # noqa: E731
    return {
        "noms": [ta.nom, tb.nom],
        "n_etats": list(c.n_etats),
        "max_re": list(c.max_re),
        "stable": [c.max_re[0] < 0, c.max_re[1] < 0],
        "ecart_tension": c.ecart_tension,
        "paires": [paire(m) for m in c.paires if garder(max(m.f_a, m.f_b))],
        "degrades": [paire(m) for m in c.degrades],
        "disparus": [seul(m) for m in c.disparus if garder(m.f)],
        "apparus": [seul(m) for m in c.apparus if garder(m.f)],
        "tables": [
            {
                "table": d.table,
                "colonnes": TABLES.get(d.table, []),
                "ajoutees": [[brut(v) for v in r] for r in d.ajoutees],
                "retirees": [[brut(v) for v in r] for r in d.retirees],
                "modifiees": [
                    {
                        "cle": brut(k),
                        "colonne": col,
                        "avant": brut(va),
                        "apres": brut(vb),
                    }
                    for k, col, va, vb in d.modifiees
                ],
            }
            for d in c.tables
        ],
        "point": [
            {"bus": bus, "V_a": v[0], "V_b": v[1], "theta_a": v[2], "theta_b": v[3]}
            for bus, v in sorted(c.point.items())
        ],
    }


def code_python(reseau: dict, nom_catalogue: str | None = None) -> dict:
    """Le script Python qui reprend le reseau du constructeur, pour developper.

    Le pont entre l'interface et le code : ce qu'on a dessine et regle ici se
    retrouve, dans un script ou un notebook, sous la forme d'un ``Atelier``
    -- l'objet qui offre en Python ce que proposent les onglets. Si le reseau
    est au catalogue, le script l'y relit par son nom ; sinon il embarque sa
    cellule de donnees.
    """
    r = _reseau_depuis(reseau)
    t = aplatir(r)
    if nom_catalogue:
        # le nom vient du client : un nom du catalogue (lettres, chiffres, _ et
        # -) et rien d'autre, ecrit en litteral Python par repr(). Un nom avec
        # un saut de ligne ajoutait sinon une ligne arbitraire au script que
        # l'utilisateur va executer.
        nom = g2egrid.reseau.catalogue._verifier_nom(nom_catalogue)
        entete = [
            "from g2egrid_api.atelier import Atelier",
            "",
            f"a = Atelier.catalogue({nom!r})",
        ]
    else:
        entete = [
            "from g2egrid_api.atelier import Atelier",
            "from g2egrid.reseau.tables import Tables",
            "import numpy as np, math",
            "",
            "# --- les tables du reseau, telles que l'interface les tient ---",
            exporter(t),
            "",
            "a = Atelier.depuis_globales(globals())",
        ]
    corps = [
        "",
        "# --- point de fonctionnement ---",
        "lf = a.loadflow()                # DataFrame : V, angle, P, Q par bus",
        "tr = a.transits()                # DataFrame : P, Q, pertes par branche",
        "",
        "# --- modes ---",
        "T_eig = a.T_eig()                # le tableau des modes des notebooks",
        "ss = a.etat()                    # le modele d'etat : ss.A, ss.B, ss.C",
        "",
        "# --- variations de parametre ---",
        '# a.balayage({"genre": "famille", "famille": "GFM_PARAMS", "cle": "mp"},',
        "#            np.linspace(0.005, 0.05, 20))",
        "",
        "# --- simulation ---",
        '# d = a.temporel("SM1.p_ref", 0.01, duree=2.0, mesures=["V_bus2"])',
    ]
    if r.instances:
        corps += [
            "",
            "# --- impact de la distribution sur le transport ---",
            "# im = a.impact()",
        ]
    return {"python": "\n".join(entete + corps) + "\n"}
