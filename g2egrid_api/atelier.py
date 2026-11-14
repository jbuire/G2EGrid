"""L'atelier en Python : ce que font les onglets de l'interface, en code.

On construit un reseau a la souris, on l'enregistre "Au catalogue", et on le
reprend ici pour developper -- ou l'inverse : on l'ecrit en code, et on
l'ouvre dans l'interface. Les deux passent par les MEMES fonctions
(``g2egrid_api.analyse``, ``g2egrid.reseau.construction``) : un resultat obtenu dans un
notebook est celui que montre l'onglet correspondant, au bit pres.

::

    from g2egrid_api.atelier import Atelier

    a = Atelier.catalogue("wscc9")          # ou Atelier.vide("demo")
    a.un_slack_tout_pq(slack_bus=1)         # onglet Point de fonctionnement
    a.loadflow()                            # DataFrame par bus
    a.transits()                            # DataFrame par branche
    T = a.T_eig()                           # onglet Modes, format des notebooks
    a.balayage({"genre": "famille", "famille": "GFM_PARAMS", "cle": "mp"},
               [0.01, 0.02, 0.05])          # onglet Variations
    a.temporel("SM1.p_ref", 0.01, duree=2)  # onglet Temporel
    a.impact()                              # onglet Transport / distribution
    a.comparer(Atelier.catalogue("wscc9"))  # onglet Comparer
    a.enregistrer("mon_reseau")             # visible dans l'interface

Construire en code (memes conventions que la palette) ::

    a = Atelier.vide("demo")
    b1 = a.bus(230e3)
    b2 = a.bus(230e3, Pl_W=80e6, Ql_var=20e6)
    a.ligne(b1, b2, R_pu=0.01, X_pu=0.1)
    a.equipement("IB", b1, slack=True)
    a.equipement("GFM1", b2, Sn_VA=100e6, P_W=50e6)

Un reseau a zones se tient sous sa forme structuree (``a.reseau``) ; les
calculs portent sur ses tables APLATIES (``a.tables``), comme dans
l'interface. La construction en code, elle, s'applique a un reseau a plat.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

import g2egrid.reseau.catalogue as cat
from g2egrid.reseau import construction
from g2egrid.reseau.tables import TABLES, Tables
from g2egrid.reseau.zones import Reseau, aplatir, reseau_vers_json

from . import analyse

__all__ = ["Atelier"]


class Atelier:
    def __init__(self, tables: Tables | None = None, reseau: Reseau | None = None):
        if (tables is None) == (reseau is None):
            raise ValueError("donner soit des tables, soit un reseau a zones")
        self.reseau = reseau
        self._tables = tables

    # ------------------------------------------------------------ ouvrir
    @classmethod
    def catalogue(cls, nom: str) -> "Atelier":
        """Un reseau du catalogue -- exemple, notebook, ou sauvegarde."""
        r = cat.charger_reseau(nom)
        return cls(reseau=r) if r.instances else cls(tables=aplatir(r))

    @classmethod
    def vide(
        cls, nom: str = "reseau", f0: float = 50.0, Sb: float = 100e6
    ) -> "Atelier":
        return cls(tables=construction.reseau_vide(nom, f0, Sb))

    @classmethod
    def depuis_json(cls, d: dict) -> "Atelier":
        """N'importe quelle forme : fichier .g2egrid.json, reponse de l'API..."""
        lu = analyse.lire_fichier(d)
        if lu["instances"]:
            return cls(reseau=analyse._reseau_depuis(lu))
        return cls(tables=analyse.tables_depuis_json(lu["tables"]))

    @classmethod
    def depuis_globales(cls, g: dict[str, Any]) -> "Atelier":
        """Les tables d'une cellule de donnees executee (``Y_DER = [...]`` ...)."""
        t = Tables(
            **{k: list(g.get(k, [])) for k in TABLES},
            **{
                k: dict(g.get(k, {}))
                for k in ("SM_PARAMS", "GFM_PARAMS", "GFL_PARAMS", "PARAMS_PAR_BUS")
            },
            f0=float(g.get("f0", 50.0)),
            Sb=float(g.get("Sb", 100e6)),
            nom=str(g.get("NOM_RESEAU", "reseau")),
            **({"NOEUDS": str(g["NOEUDS"])} if "NOEUDS" in g else {}),
            **({"R_NOEUDS": float(g["R_NOEUDS"])} if "R_NOEUDS" in g else {}),
        )
        return cls(tables=t)

    # ------------------------------------------------------------ tables
    @property
    def tables(self) -> Tables:
        """Les tables A PLAT : ce que resolvent tous les calculs."""
        return aplatir(self.reseau) if self.reseau is not None else self._tables

    def _a_plat(self) -> Tables:
        if self.reseau is not None:
            raise ValueError(
                "reseau a zones : la construction en code s'applique a un reseau "
                "a plat (a.reseau.transport, ou une zone de a.reseau.zones)"
            )
        return self._tables

    def table(self, nom: str) -> pd.DataFrame:
        """Une table de construction, en DataFrame -- l'onglet Tables."""
        return pd.DataFrame(getattr(self.tables, nom), columns=TABLES[nom])

    # ------------------------------------------------------ construction
    def bus(self, Un_V: float, **kw) -> int:
        return construction.poser_bus(self._a_plat(), Un_V, **kw)

    def ligne(self, de: int, vers: int, **kw) -> int:
        return construction.poser_ligne(self._a_plat(), de, vers, **kw)

    def transfo(self, de: int, vers: int, Sn_VA: float, **kw) -> int:
        return construction.poser_transfo(self._a_plat(), de, vers, Sn_VA, **kw)

    def equipement(self, famille: str, bus: int, **kw) -> int:
        return construction.poser_equipement(self._a_plat(), famille, bus, **kw)

    def charge(
        self, bus: int, P_W: float, Q_var: float = 0.0, genre: str = "Z"
    ) -> None:
        construction.poser_charge(self._a_plat(), bus, P_W, Q_var, genre)

    def shunt(self, bus: int, B_pu: float) -> int:
        return construction.poser_shunt(self._a_plat(), bus, B_pu)

    def equipements(self) -> pd.DataFrame:
        return pd.DataFrame(construction.equipements(self.tables))

    # ----------------------------------------- point de fonctionnement
    def point(self, bus: int, **kw) -> None:
        """``type`` (1 slack, 2 PV, 3 PQ), ``V_pu``, ``P_W``, ``Q_var`` d'un bus."""
        construction.fixer_point(self._a_plat(), bus, **kw)

    def un_slack_tout_pq(self, slack_bus: int) -> None:
        d = analyse.un_slack_tout_pq(self.tables, slack_bus)
        self._tables = analyse.tables_depuis_json(d["tables"])
        self.reseau = None

    def loadflow(self) -> pd.DataFrame:
        d = analyse.loadflow(self.tables)
        df = pd.DataFrame(d["bus"]).set_index("id")
        df.attrs.update(iterations=d["iterations"], mismatch=d["mismatch"])
        return df

    def transits(self) -> pd.DataFrame:
        return pd.DataFrame(analyse.loadflow(self.tables)["branches"]).set_index("nom")

    # ------------------------------------------------------------ modes
    def etat(self):
        """Le modele d'etat linearise (``StateSpace`` : A, B, C, D, noms)."""
        return analyse.construire(self.tables).etat_stationnaire()

    def T_eig(self) -> pd.DataFrame:
        """Le tableau des modes, au format des notebooks."""
        d = analyse.t_eig(self.tables)
        rows = []
        for m in d["modes"]:
            r = {
                "Re": m["re"],
                "Im": m["im"],
                "f_und [Hz]": m["f_und"],
                "f_dam [Hz]": m["f"],
                "zeta [%]": m["zeta"],
            }
            for j, e in enumerate(m["etats"], 1):
                r[f"Etat {j}"], r[f"%{j}"] = e["nom"], e["part"]
            rows.append(r)
        df = pd.DataFrame(rows)
        df.index.name = "Mode"
        return df

    def mode(self, indice: int) -> dict:
        """Forme modale et sensibilites d'un mode (indice du spectre)."""
        return analyse.detail_mode(self.tables, indice)

    # ------------------------------------------------------- variations
    def balayage(
        self, cible: dict, valeurs, f_max: float | None = None, n_modes: int = 12
    ) -> pd.DataFrame:
        d = analyse.balayage(
            self.tables, cible, list(np.asarray(valeurs, float)), f_max, n_modes
        )
        rows = []
        for p in d["pas"]:
            for j, m in enumerate(p.get("modes", [])):
                rows.append({"valeur": p["valeur"], "lieu": j, **m})
        return pd.DataFrame(rows)

    # ------------------------------------------------------ simulation
    def temporel(self, entree: str, amplitude: float = 0.01, **kw) -> dict:
        return analyse.temporel(self.tables, entree, amplitude, **kw)

    def validation(self, **kw) -> dict:
        return analyse.validation(self.tables, **kw)

    # ------------------------------------------- transport / distribution
    def impact(self, f_max: float = 5.0, par_poste: bool = True) -> dict:
        if self.reseau is None:
            raise ValueError("pas de zone : l'impact se mesure sur un reseau a zones")
        return analyse.impact(reseau_vers_json(self.reseau), f_max, par_poste)

    # ---------------------------------------------------------- comparer
    def comparer(
        self, reference: "Atelier", f_max: float | None = None
    ) -> pd.DataFrame:
        """Les modes de ``reference`` (A) face a ceux de ce reseau (B).

        Une ligne par mode apparie ; les modes presents d'un seul cote ont
        ``f_a`` ou ``f_b`` vide. Le detail -- tables modifiees, point de
        fonctionnement -- est dans ``df.attrs``, comme dans l'onglet Comparer.
        """
        d = analyse.comparer(reference.json(), self.json(), f_max)
        rows = [
            {
                k: m[k]
                for k in (
                    "f_a",
                    "zeta_a",
                    "f_b",
                    "zeta_b",
                    "ecart_f",
                    "ecart_zeta",
                    "dominant_b",
                )
            }
            for m in d["paires"]
        ]
        rows += [
            {"f_a": m["f"], "zeta_a": m["zeta"], "dominant_b": m["dominant"]}
            for m in d["disparus"]
        ]
        rows += [
            {"f_b": m["f"], "zeta_b": m["zeta"], "dominant_b": m["dominant"]}
            for m in d["apparus"]
        ]
        df = pd.DataFrame(rows)
        df.attrs.update(
            {k: d[k] for k in ("tables", "point", "max_re", "ecart_tension")}
        )
        return df

    # ------------------------------------------------------ enregistrer
    def json(self) -> dict:
        if self.reseau is not None:
            return reseau_vers_json(self.reseau)
        return analyse.convertir_fichier(
            {"tables": analyse.tables_vers_json(self._tables)}
        )

    def enregistrer(self, nom: str) -> dict:
        """Dans les sauvegardes du catalogue : l'interface l'ouvre aussitot."""
        return analyse.enregistrer_fichier(self.json(), nom)

    def __repr__(self) -> str:
        t = self.tables
        z = f", {len(self.reseau.instances)} instances de zone" if self.reseau else ""
        return f"<Atelier {t.nom!r} : {len(t.Y_DER)} bus{z}>"
