"""Le point de fonctionnement est-il dans le domaine des limiteurs ?

POURQUOI -- un convertisseur ne fournit pas plus que ``Imax`` (pu de son
``Sn``), et son onduleur ne fabrique pas plus que ``Mmax.v_dc``. Ces limites
sont des MODES DISCRETS (:mod:`g2egrid.elements.limiteurs`) : la simulation
temporelle les surveille, mais le modele linearise est TOUJOURS calcule en
mode libre, ou le limiteur n'existe pas. Tant que le point de fonctionnement
est dans les limites, c'est exact : le limiteur n'y agit pas. Au-dela -- on
demande 200 MW a un GFM de 130 MVA --, le convertisseur reel serait en
limitation, et son modele petit signal n'est plus celui que l'on calcule :
la consigne de courant n'y est plus tenue par la boucle de tension ou de
puissance, mais par le limiteur. Le spectre libre reste calculable, et il
peut sembler parfaitement stable ; il ne decrit simplement plus le systeme.
Les etudes qui poussent la puissance -- capacite d'accueil, robustesse,
continuation -- y tombent les premieres.

QUOI -- pour chaque element qui porte un limiteur (toute classe qui expose
``grandeurs_limitees``, c'est-a-dire tout :class:`~g2egrid.elements.limiteurs.LimiteurCourant` :
GFM, GFL, MMC et leurs variantes a port continu), on compare au point
d'equilibre la grandeur que lit la GARDE du limiteur a sa limite :

* courant : le module de la consigne de courant libre, qui vaut a l'equilibre
  le courant COTE CONVERTISSEUR ``|i_s|`` (et non le courant cote reseau
  ``|i_g|``, qui en differe du courant du condensateur du filtre), face a
  ``Imax``, les deux en pu de l'equipement ;
* modulation (si la saturation de tension est reglee, ``Mmax > 0``) :
  ``|e|/v_dc`` face a ``Mmax``.

Le critere est celui de la simulation temporelle quand elle demarre
(``ModeleNonLineaire.simuler`` : une garde deja franchie de plus de
:data:`TOLERANCE` bascule le limiteur avant d'integrer) : un depassement est
signale ICI si et seulement si la simulation commencerait en limitation.
``tests/test_limite_equilibre.py`` le verifie.

LIMITATION PAR IMPEDANCE VIRTUELLE -- la VI (``kRv > 0``) n'est PAS un mode :
elle est dans les equations du mode libre, donc dans le modele linearise. Un
point d'equilibre ou elle est declenchee (``I_vi < |i_s|``) reste dans le
domaine du modele ; il n'est hors limite que si ``|i_s| > Imax``, comme sans
VI. La grandeur comparee est la meme (``tests/test_limitation_impedance_virtuelle.py``).

LA FORME DU SIGNAL -- deux voies, pour deux lecteurs :

* l'attribut ``Network.depassements``, une liste de :class:`Depassement`
  (vide si tout va bien), mis a jour par chaque ``solve_equilibrium`` : ce
  que lisent les scripts d'etude, la continuation, la robustesse et
  l'atelier, sans analyser de texte ;
* un avertissement ``warnings.warn`` de categorie :class:`DepassementLimite`,
  pour l'utilisateur qui ne regarde pas l'attribut. Pas d'exception : une
  etude doit pouvoir explorer au-dela de la limite en le sachant. Pour la
  faire taire, ou au contraire la rendre bloquante ::

      warnings.simplefilter("ignore", DepassementLimite)
      warnings.simplefilter("error", DepassementLimite)
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

__all__ = [
    "TOLERANCE",
    "Depassement",
    "DepassementLimite",
    "depassements",
    "message",
    "signaler",
]

#: Ecart minimal (en pu de l'equipement, ou en indice de modulation) au-dela
#: duquel une limite est dite depassee. Ce n'est pas une marge nouvelle :
#: c'est le seuil auquel ``ModeleNonLineaire.simuler`` juge une garde deja
#: franchie au depart (``temporel.py``, ``g[4] * g[2] > 1e-9``), si bien que
#: les deux diagnostics s'accordent. Il est de l'ordre de la tolerance du
#: Newton d'equilibre (``solve_equilibrium(tol=1e-9)``) : il n'ecarte que
#: l'arrondi d'un point POSE exactement sur la limite, jamais un depassement
#: physique (1e-9 pu de 130 MVA, c'est 0,13 VA).
TOLERANCE = 1e-9

#: drapeau du limiteur -> (grandeur comparee, nom de la limite, unite)
_GRANDEURS = {
    "mode_lim": ("courant", "Imax", "pu"),
    "mode_vsat": ("modulation", "Mmax", ""),
}


class DepassementLimite(UserWarning):
    """Le point d'equilibre d'un convertisseur est au-dela de l'une de ses
    limites : le modele linearise (mode libre) ne le represente plus."""


@dataclass(frozen=True)
class Depassement:
    """Un limiteur dont la grandeur libre depasse la limite a l'equilibre."""

    #: nom de l'element (``GFM3``)
    element: str
    #: bus alternatif de l'element
    bus: int
    #: drapeau du limiteur : ``"mode_lim"`` (courant) ou ``"mode_vsat"``
    #: (modulation)
    drapeau: str
    #: grandeur libre a l'equilibre : ``|i_s|`` en pu de l'equipement, ou
    #: ``|e|/v_dc``
    valeur: float
    #: ``Imax`` (pu de l'equipement) ou ``Mmax``
    limite: float
    #: puissance assignee de l'equipement [VA], la base de ``valeur``
    Sn: float | None = None

    @property
    def grandeur(self) -> str:
        """``"courant"`` ou ``"modulation"``."""
        return _GRANDEURS.get(self.drapeau, (self.drapeau,))[0]

    @property
    def nom_limite(self) -> str:
        """``"Imax"`` ou ``"Mmax"``."""
        return _GRANDEURS.get(self.drapeau, ("", self.drapeau))[1]

    @property
    def rapport(self) -> float:
        """``valeur / limite`` : 1,25 pour un courant 25 % au-dessus."""
        return self.valeur / self.limite

    def __str__(self) -> str:
        unite = _GRANDEURS.get(self.drapeau, ("", "", ""))[2]
        u = f" {unite}" if unite else ""
        base = f", base Sn = {self.Sn / 1e6:.4g} MVA" if self.Sn and unite else ""
        return (
            f"{self.element} (bus {self.bus}) : {self.grandeur} {self.valeur:.4f}{u}"
            f" > {self.nom_limite} = {self.limite:.4g}{u}"
            f" (+{100 * (self.rapport - 1):.1f} %{base})"
        )

    def en_dict(self) -> dict:
        """Pour le JSON de l'atelier."""
        return {
            "element": self.element,
            "bus": self.bus,
            "grandeur": self.grandeur,
            "drapeau": self.drapeau,
            "limite_nom": self.nom_limite,
            "valeur": self.valeur,
            "limite": self.limite,
            "rapport": self.rapport,
            "Sn": self.Sn,
        }


def depassements(net) -> list[Depassement]:
    """Les limiteurs dont la grandeur libre depasse la limite au point
    d'equilibre de ``net`` (``solve_equilibrium`` deja fait).

    Ne compte que les limiteurs LIBRES (drapeau nul), ceux du modele
    linearise : un drapeau deja leve dirait que le modele est celui du mode
    limite, et le point y est coherent. Un convertisseur bloque n'a pas de
    limiteur surveille (voir ``grandeurs_limitees``). Un element sans point
    d'equilibre est ignore.
    """
    out: list[Depassement] = []
    for el in net.all_elements:
        lire = getattr(el, "grandeurs_limitees", None)
        if lire is None or el.x0 is None:
            continue
        for d, valeur, limite in lire(el.compiled(), el.x0, el.z0, el.u0, el.p):
            if el.p[el.params.index(d)] != 0.0:
                continue
            if valeur - limite > TOLERANCE:
                out.append(
                    Depassement(
                        element=el.name,
                        bus=int(el.buses[0]),
                        drapeau=d,
                        valeur=float(valeur),
                        limite=float(limite),
                        Sn=getattr(el, "Sn", None),
                    )
                )
    return out


def message(liste: list[Depassement], quoi: str = "") -> str:
    """Le texte de l'avertissement : une ligne par depassement."""
    tete = (f"{quoi} : " if quoi else "") + (
        "point d'equilibre au-dela d'une limite de convertisseur ; le modele "
        "linearise (mode libre) ne represente plus l'equipement, qui serait "
        "en limitation :"
    )
    return "\n".join([tete, *(f"  {d}" for d in liste)])


def signaler(liste: list[Depassement], quoi: str = "", stacklevel: int = 3) -> None:
    """Un avertissement :class:`DepassementLimite` si ``liste`` n'est pas
    vide. ``stacklevel`` = 3 designe l'appelant de la fonction qui signale
    (celui de ``solve_equilibrium``, de ``continuer``...)."""
    if liste:
        warnings.warn(message(liste, quoi), DepassementLimite, stacklevel=stacklevel)
