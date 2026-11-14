"""Cascades de declenchement en regime permanent (quasi statique).

Un producteur qui declenche change les tensions de ses voisins ; si l'un
d'eux passe son seuil de protection, il declenche a son tour, et ainsi de
suite. :func:`cascade_surtension` deroule cette suite de REGIMES PERMANENTS :

1. les unites ``initiales`` sont ouvertes (``Tables.OUVERTS``) ;
2. un load flow du reseau en service (:func:`g2egrid.reseau.disjoncteurs.en_service`)
   donne les nouvelles tensions ;
3. toutes les unites protegees encore en service dont le bus depasse
   ``v_max`` declenchent ENSEMBLE (``ordre="vague"``) -- ou seulement la plus
   haute d'entre elles (``ordre="une_a_une"``) ; retour en 2, jusqu'a ce
   qu'aucune ne depasse plus (cascade arretee), qu'il n'en reste plus, ou que
   le load flow ne converge plus.

C'est l'hypothese des relais temporises : entre deux vagues, le reseau a le
temps de rejoindre son nouveau regime permanent (les temporisations de
surtension, quelques centaines de millisecondes, depassent les transitoires
electriques). Les deux ordres bornent ce que font des relais a temps constant
de MEME reglage : s'ils demarrent ensemble et tiennent, ils declenchent
ensemble (la vague) ; mais un relais a retombee retombe sur le creux que
provoque le declenchement d'un voisin et repart de zero -- les declenchements
s'egrenent alors un par un, et chacun change la tension des suivants (une
unite qui aurait declenche dans la vague peut repasser sous le seuil). La
cascade depend de l'ordre. Les regulateurs restent ceux du load flow : les machines
tiennent leur tension (bus PV), le slack reprend l'ecart de puissance, les
prises des transformateurs ne bougent pas -- un regleur en charge agit en
dizaines de secondes, apres la cascade. La simulation temporelle, avec
:func:`g2egrid.temporel.protections_ac`, verifie ce deroule pas a pas.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from g2egrid.reseau.disjoncteurs import en_service, nom_element
from g2egrid.noyau.progres import signaler
from g2egrid.reseau.tables import Tables, assembler

__all__ = ["Cascade", "Etape", "cascade_surtension", "tensions_bus", "unites_protegees"]

#: tables des unites qu'un relais de surtension peut faire declencher
FAMILLES = {
    "GFL": ("T_GFL1s", "T_GFL2s"),
    "GFM": ("T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s"),
}


@dataclass
class Etape:
    """Un regime permanent de la cascade."""

    rang: int
    #: unites ouvertes a CE regime (initiales au rang 0, puis la vague d'avant)
    ouvertes: list[str]
    #: toutes les unites hors service a ce regime
    hors_service: list[str]
    #: tension de chaque bus en service (pu)
    V: dict[int, float]
    #: unites encore en service au-dessus du seuil : la vague suivante
    au_dessus: list[str]
    converge: bool = True

    @property
    def v_max(self) -> float:
        return max(self.V.values()) if self.V else float("nan")

    @property
    def bus_v_max(self) -> int | None:
        return max(self.V, key=self.V.get) if self.V else None


@dataclass
class Cascade:
    seuil: float
    initiales: list[str]
    #: le regime de depart, avant toute ouverture
    V0: dict[int, float]
    etapes: list[Etape] = field(default_factory=list)
    #: "arretee" (plus rien au-dessus du seuil), "totale" (plus aucune unite
    #: protegee en service), "divergence" (load flow non converge),
    #: "longueur" (n_max vagues)
    issue: str = ""

    @property
    def declenchees(self) -> list[str]:
        """Les unites declenchees PAR la cascade (initiales exclues)."""
        if not self.etapes:
            return []
        return [n for n in self.etapes[-1].hors_service if n not in self.initiales]


def unites_protegees(t: Tables, familles=("GFL",)) -> list[tuple[str, int, int, str]]:
    """``(table, ID, bus, nom)`` des unites que le relais surveille."""
    out = []
    for f in familles:
        if f not in FAMILLES:
            raise ValueError(f"famille {f!r} inconnue ; {sorted(FAMILLES)}")
        for table in FAMILLES[f]:
            for r in getattr(t, table):
                out.append(
                    (table, int(r[0]), int(r[1]), nom_element(t, table, int(r[0])))
                )
    return out


def tensions_bus(t: Tables) -> tuple[dict[int, float], bool]:
    """Le module de tension de chaque bus EN SERVICE apres load flow, et la
    convergence. Les disjoncteurs ouverts de ``t`` sont appliques."""
    t_es, _rap = en_service(t)
    lf = assembler(t_es).run_loadflow()
    bus = sorted(int(r[0]) for r in t_es.Y_DER)
    V = np.abs(np.asarray(lf.V)[: len(bus)])
    return {b: float(v) for b, v in zip(bus, V, strict=True)}, bool(lf.converged)


def cascade_surtension(
    t: Tables,
    initiales,
    v_max: float = 1.10,
    familles=("GFL",),
    n_max: int = 50,
    ordre: str = "vague",
) -> Cascade:
    """Deroule la cascade de surtension declenchee par l'ouverture des unites
    ``initiales`` (noms assembles, ``"GFL303"``, ou couples ``(table, ID)``).

    Une unite protegee declenche quand la tension de SON bus depasse
    ``v_max`` (strictement). Les unites deja ouvertes dans ``t.OUVERTS``
    restent ouvertes. ``ordre`` : ``"vague"`` (toutes celles au-dessus du
    seuil declenchent ensemble) ou ``"une_a_une"`` (la plus haute seule ;
    ``n_max`` compte alors des declenchements).
    """
    if ordre not in ("vague", "une_a_une"):
        raise ValueError(f"ordre {ordre!r} : 'vague' ou 'une_a_une'")
    unites = unites_protegees(t, familles)
    par_nom = {n: (tab, i) for tab, i, _b, n in unites}
    # toute unite manoeuvrable peut ouvrir la cascade, pas seulement les protegees
    from g2egrid.reseau.disjoncteurs import ouvrables

    tous = {o["nom"]: (o["table"], o["id"]) for o in ouvrables(t)}
    init: list[tuple[str, int]] = []
    noms_init: list[str] = []
    for x in initiales:
        if isinstance(x, str):
            if x not in tous:
                raise ValueError(f"unite {x!r} inconnue ; ouvrables : {sorted(tous)}")
            init.append(tous[x])
            noms_init.append(x)
        else:
            tab, i = x
            init.append((str(tab), int(i)))
            noms_init.append(nom_element(t, str(tab), int(i)))
    deja = [(a, int(b)) for a, b in t.OUVERTS]
    V0, ok0 = tensions_bus(t)
    if not ok0:
        raise RuntimeError("le load flow de depart ne converge pas")
    c = Cascade(seuil=float(v_max), initiales=noms_init, V0=V0)

    ouverts = list(dict.fromkeys(deja + init))
    vague = noms_init
    for rang in range(n_max):
        signaler(len(ouverts) / max(1, len(par_nom) + len(deja)), f"regime {rang}")
        t2 = t.copie()
        t2.OUVERTS = list(ouverts)
        hors = [n for n, (tab, i) in par_nom.items() if (tab, i) in set(ouverts)]
        hors += [n for n in noms_init if n not in hors]
        try:
            V, ok = tensions_bus(t2)
        except (RuntimeError, ValueError, np.linalg.LinAlgError):
            V, ok = {}, False
        if not ok:
            c.etapes.append(Etape(rang, vague, hors, V, [], converge=False))
            c.issue = "divergence"
            return c
        bus_de = {n: b for _tab, _i, b, n in unites}
        dessus = [
            n
            for n, (tab, i) in par_nom.items()
            if (tab, i) not in set(ouverts) and V.get(bus_de[n], 0.0) > v_max
        ]
        c.etapes.append(Etape(rang, vague, hors, V, dessus))
        if not dessus:
            restantes = [n for n, k in par_nom.items() if k not in set(ouverts)]
            c.issue = "totale" if not restantes else "arretee"
            return c
        if ordre == "une_a_une":
            dessus = [max(dessus, key=lambda n: V[bus_de[n]])]
        ouverts += [par_nom[n] for n in dessus]
        vague = dessus
    c.issue = "longueur"
    return c
