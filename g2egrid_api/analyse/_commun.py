"""Le socle commun de l'atelier : imports, constantes, tables, cache des assemblages.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Any
from g2egrid.reseau.network import Network
from g2egrid.reseau.tables import Tables, tables_vers_json
from g2egrid.temporel import ModeleNonLineaire


#: En dessous de cette pulsation [rad/s], un mode est tenu pour le mode de
#: JAUGE et non pour un mode physique. Un reseau sans bus infini n'a pas de
#: reference d'angle absolue : son spectre porte un zero exact, que l'algebre
#: rend a 1e-12 pres. Son amortissement vaut alors -Re/|lambda| = du bruit
#: d'arrondi divise par du bruit d'arrondi, soit +100 ou -100 % au hasard du
#: point de fonctionnement -- ce qui lui vaut d'arriver en tete du classement
#: par amortissement, d'occuper un lieu a lui seul, et de faire afficher
#: "max Re = +0.000" sur un reseau parfaitement stable. Le seuil est choisi
#: dans le fosse : sur wscc9_gfm, la jauge est a 4e-12 et le mode physique le
#: plus lent a 0,19 rad/s, dix ordres de grandeur plus haut.
JAUGE = 1e-6


#: Au-dela de cet amortissement [%], une paire complexe ne fait plus une
#: oscillation : elle decroit avant d'avoir boucle une periode. C'est la
#: valeur que nomme deja la docstring de ``modal.modes``.
ZETA_OSCILLANT = 95.0


#: Ce qu'on peut poser sur un bus, et dans quelle table ca atterrit.
EQUIPEMENTS: dict[str, dict[str, Any]] = {
    "SM": {"table": "T_SMs", "titre": "Machine synchrone", "couleur": "#1f4e79"},
    "GFM1": {
        "table": "T_GFM1s",
        "titre": "Grid-forming, statisme 1er ordre, bus continu ideal",
        "couleur": "#00743f",
    },
    "GFM2": {
        "table": "T_GFM2s",
        "titre": "Grid-forming, statisme 2e ordre, bus continu ideal",
        "couleur": "#00743f",
    },
    "GFM3": {
        "table": "T_GFM3s",
        "titre": "Grid-forming, statisme 1er ordre, bus continu dynamique",
        "couleur": "#00743f",
    },
    "GFM4": {
        "table": "T_GFM4s",
        "titre": "Grid-forming, statisme 2e ordre, bus continu dynamique",
        "couleur": "#00743f",
    },
    "GFL1": {
        "table": "T_GFL1s",
        "titre": "Grid-following, bus continu ideal",
        "couleur": "#bf7c17",
    },
    "GFL2": {
        "table": "T_GFL2s",
        "titre": "Grid-following, bus continu dynamique",
        "couleur": "#bf7c17",
    },
    "IB": {"table": "T_IBs", "titre": "Bus infini", "couleur": "#6b4c9a"},
    "Charge": {
        "table": "T_Loads",
        "titre": "Charge speciale (courant ou puissance constante)",
        "couleur": "#bf2617",
    },
    "Shunt": {
        "table": "T_Shunts",
        "titre": "Shunt : condensateur (B > 0) ou reactance (B < 0)",
        "couleur": "#2a8fb5",
    },
    "Blocs": {
        "table": "T_Blocs",
        "titre": "Element decrit par un schema-bloc (onglet Blocs)",
        "couleur": "#7a3e9d",
    },
}


#: Les dynamiques qu'on peut activer ou non, par famille d'equipement, pour
#: les cases a cocher de l'atelier. ``sans`` dit ce que devient le groupe
#: DECOCHE : ``"statique"`` (dx/dt = 0, la dynamique suit instantanement) ou
#: ``"fige"`` (x = x0, la dynamique est gelee) -- voir
#: :meth:`g2egrid.noyau.element.Element.reduire`. Les deux premiers groupes d'un GFM
#: et le premier d'un GFL ne passent pas par la reduction : ils choisissent la
#: TABLE (GFM1 a GFM4, GFL1 ou GFL2), qui porte deja la variante.
GROUPES_DYNAMIQUES: dict[str, list[dict]] = {
    "GFM": [
        {
            "cle": "statisme2",
            "titre": "Statisme du 2e ordre (inertie virtuelle)",
            "table": True,
        },
        {
            "cle": "dc",
            "titre": "Bus continu dynamique (batterie, reseau DC)",
            "table": True,
        },
        {
            "cle": "mesure",
            "titre": "Filtres de mesure de P et Q",
            "etats": ["pm", "qm"],
            "sans": "statique",
        },
        {
            "cle": "boucle_v",
            "titre": "Integrateurs de la boucle de tension",
            "etats": ["M_VLd", "M_VLq"],
            "sans": "fige",
        },
        {
            "cle": "boucle_i",
            "titre": "Integrateurs de la boucle de courant",
            "etats": ["M_CLd", "M_CLq"],
            "sans": "fige",
        },
    ],
    "GFL": [
        {
            "cle": "dc",
            "titre": "Regulation de la tension continue (Vdc/Q, bus continu dynamique)",
            "table": True,
        },
        {
            "cle": "externes",
            "titre": "Integrateurs des boucles externes",
            "etats": ["M_d", "M_q"],
            "sans": "fige",
        },
        {
            "cle": "boucle_i",
            "titre": "Integrateurs de la boucle de courant",
            "etats": ["M_CLd", "M_CLq"],
            "sans": "fige",
        },
        {
            "cle": "pll",
            "titre": "Integrateur de la PLL",
            "etats": ["M_pll"],
            "sans": "fige",
        },
    ],
    "SM": [
        {
            "cle": "stator",
            "titre": "Transitoires du stator",
            "etats": ["phi_d", "phi_q"],
            "sans": "statique",
        },
        {
            "cle": "amortisseurs",
            "titre": "Enroulements amortisseurs",
            "etats": ["phi_1d", "phi_1q", "phi_2q"],
            "sans": "statique",
        },
        {
            "cle": "avr",
            "titre": "Regulateur de tension (AVR, excitatrice)",
            "etats": ["e1", "e2", "efd", "e3"],
            "sans": "fige",
        },
        {
            "cle": "pss",
            "titre": "Stabilisateur (PSS)",
            "etats": ["Dw1", "v1", "v2", "vpss"],
            "sans": "fige",
        },
        {
            "cle": "vitesse",
            "titre": "Regulateur de vitesse (turbine)",
            "etats": ["Pm"],
            "sans": "fige",
        },
    ],
}


_CACHE: dict[str, "Assemble"] = {}


_CACHE_MAX = 8


@dataclass
class Assemble:
    """Un reseau construit, et ce qu'on en a deja calcule."""

    tables: Tables
    net: Network
    lf: Any = None
    ss: Any = None
    mnl: ModeleNonLineaire | None = None
    erreur: str = ""
    #: l'equilibre est-il deja resolu ? (partage par :func:`valider` et les
    #: calculs, pour ne pas le refaire)
    resolu: bool = False

    def equilibre(self) -> Network:
        """Le reseau, son point d'equilibre resolu (une fois)."""
        if not self.resolu:
            self.net.solve_equilibrium()
            self.resolu = True
        return self.net

    def etat_stationnaire(self):
        if self.ss is None:
            self.equilibre()
            self.ss = self.net.build_state_space()
        return self.ss

    def non_lineaire(self) -> ModeleNonLineaire:
        if self.mnl is None:
            self.etat_stationnaire()
            self.mnl = ModeleNonLineaire(self.net)
        return self.mnl


def _empreinte(t: Tables) -> str:
    return hashlib.sha256(
        json.dumps(tables_vers_json(t), sort_keys=True, default=str).encode()
    ).hexdigest()
