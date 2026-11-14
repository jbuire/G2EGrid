"""Les anciens chemins des modules de g2egrid, avant la reorganisation
en sous-paquets (2026-10-28).

``import g2egrid.modal``, ``from g2egrid import network`` ou
``from g2egrid.tables import assembler`` continuent de marcher : le chercheur
d'import ci-dessous rend, sous l'ancien nom, le MEME objet module que le
nouveau chemin (``g2egrid.analyse.modal``...). Un ``monkeypatch`` sur l'ancien
nom agit donc sur le module reel, et inversement. Les etudes, les notebooks
et les tests ecrits avant la reorganisation n'ont pas a changer.

``G2EGRID_AVERTIR_ANCIENS_CHEMINS=1`` fait emettre un ``DeprecationWarning`` a
chaque ancien chemin importe, pour migrer un script. Le jour ou plus rien ne
s'en sert, ce fichier et son appel dans ``g2egrid/__init__.py`` s'en vont.
"""

from __future__ import annotations

import importlib
import importlib.abc
import importlib.util
import os
import sys
import warnings

#: ancien nom -> nouveau nom
ANCIENS: dict[str, str] = {
    "g2egrid.symbolic": "g2egrid.noyau.symbolic",
    "g2egrid.element": "g2egrid.noyau.element",
    "g2egrid.loadflow": "g2egrid.noyau.loadflow",
    "g2egrid.limites": "g2egrid.noyau.limites",
    "g2egrid.integrateur": "g2egrid.noyau.integrateur",
    "g2egrid.blocs": "g2egrid.elements.blocs",
    "g2egrid.equations": "g2egrid.elements.equations",
    "g2egrid.fiches": "g2egrid.elements.fiches",
    "g2egrid.network": "g2egrid.reseau.network",
    "g2egrid.tables": "g2egrid.reseau.tables",
    "g2egrid.catalogue": "g2egrid.reseau.catalogue",
    "g2egrid.zones": "g2egrid.reseau.zones",
    "g2egrid.construction": "g2egrid.reseau.construction",
    "g2egrid.disjoncteurs": "g2egrid.reseau.disjoncteurs",
    "g2egrid.reglage": "g2egrid.reseau.reglage",
    "g2egrid.flux": "g2egrid.reseau.flux",
    "g2egrid.modal": "g2egrid.analyse.modal",
    "g2egrid.residus": "g2egrid.analyse.residus",
    "g2egrid.sensibilites": "g2egrid.analyse.sensibilites",
    "g2egrid.impedance": "g2egrid.analyse.impedance",
    "g2egrid.frequence": "g2egrid.analyse.frequence",
    "g2egrid.descripteur": "g2egrid.analyse.descripteur",
    "g2egrid.classement": "g2egrid.analyse.classement",
    "g2egrid.identification": "g2egrid.analyse.identification",
    "g2egrid.inertie": "g2egrid.analyse.inertie",
    "g2egrid.comparaison": "g2egrid.analyse.comparaison",
    "g2egrid.court_circuit": "g2egrid.analyse.court_circuit",
    "g2egrid.robustesse": "g2egrid.scenarios.robustesse",
    "g2egrid.continuation": "g2egrid.scenarios.continuation",
    "g2egrid.progres": "g2egrid.noyau.progres",
    "g2egrid.simplification": "g2egrid.scenarios.simplification",
    "g2egrid.domaines": "g2egrid.scenarios.domaines",
    "g2egrid.cascade": "g2egrid.scenarios.cascade",
    "g2egrid.multi_echelle": "g2egrid.multi_niveau.multi_echelle",
    "g2egrid.impact": "g2egrid.multi_niveau.impact",
    "g2egrid.separation": "g2egrid.multi_niveau.separation",
    "g2egrid.equivalents": "g2egrid.multi_niveau.equivalents",
    "g2egrid.interop": "g2egrid.formats.interop",
    "g2egrid.export": "g2egrid.formats.export",
    "g2egrid.notebook": "g2egrid.formats.notebook",
    "g2egrid.andes_croise": "g2egrid.validation.andes_croise",
}


class _Alias(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """Rend sous l'ancien nom le module importe sous le nouveau."""

    def find_spec(self, nom, chemin=None, cible=None):
        if nom in ANCIENS:
            return importlib.util.spec_from_loader(nom, self)
        return None

    def create_module(self, spec):
        if os.environ.get("G2EGRID_AVERTIR_ANCIENS_CHEMINS"):
            warnings.warn(
                f"{spec.name} s'appelle maintenant {ANCIENS[spec.name]}",
                DeprecationWarning,
                stacklevel=3,
            )
        module = importlib.import_module(ANCIENS[spec.name])
        self._spec_reelle = module.__spec__
        return module

    def exec_module(self, module):
        # l'import vient de poser l'ALIAS comme __spec__ du module reel :
        # le module garde sa vraie identite (pickle, inspect, Sphinx)
        module.__spec__ = self._spec_reelle


def installer() -> None:
    if not any(isinstance(f, _Alias) for f in sys.meta_path):
        sys.meta_path.insert(0, _Alias())


def module_ancien(nom: str):
    """``g2egrid.<nom>`` pour un ancien nom court (``g2egrid.modal``...)."""
    complet = f"g2egrid.{nom}"
    if complet in ANCIENS:
        return importlib.import_module(complet)
    raise AttributeError(nom)
