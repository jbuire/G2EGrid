"""Le repere de reference, en element a part.

D'ORDINAIRE, LE REPERE EST PORTE PAR LE SLACK -- ``theta_g``, l'angle que
toutes les unites retranchent du leur, est la sortie ``theta_out`` de
l'element pose au bus slack (``Network.cablage``), et la pulsation que voient
les elements passifs est la vitesse de ce meme element. C'est la convention de
G2Elin V1, et celle qui donne les matrices comparables terme a terme.

Elle a un defaut : le slack ne peut plus disparaitre. Declencher la machine de
reference, c'est emporter le repere avec elle. Ce bloc prend alors le relais :
un seul etat, l'angle du repere, qui tourne a la vitesse d'une machine qu'il
SUIT ::

    dtheta/dt = wb . (w_in - w_ref)        sorties : (theta_out = theta, w = w_in)

ou ``w_in`` est cablee sur la vitesse de la machine suivie (``Network
.repere_suit``), ou vaut 1 s'il ne reste aucune machine a suivre. Il reprend
EXACTEMENT les deux sorties que l'element de reference fournissait -- un angle
et une vitesse --, si bien que rien d'autre dans le reseau ne voit la
difference. Initialise a la valeur de ``theta_g`` a l'instant du changement,
il la prolonge sans saut : seule sa derivee change, et un repere n'est qu'un
choix de coordonnees -- la physique ne depend pas de la vitesse a laquelle on
le fait tourner, pourvu que les elements passifs et les unites voient la meme.

POURQUOI SUIVRE UNE MACHINE et non tourner a la vitesse nominale : apres une
perte de production, le reseau se stabilise a une frequence un peu differente
(statisme). Dans un repere fixe a 50 Hz, les grandeurs dq tourneraient alors
indefiniment, sans equilibre a lineariser. F. Kelada le mesure dans son
portage (l'erreur du modele linearise passe de 5 % a 18 % sur 1,5 s) ; son
``frame.py`` est l'origine de l'idee -- rien n'en est repris que l'idee.

Le repere n'a pas de position absolue : un reseau sans bus infini garde un
mode de jauge, comme avant.
"""

from __future__ import annotations

from g2egrid.noyau.element import Element, InitContext

__all__ = ["Repere"]


class Repere(Element):
    """Le repere dq global, porte par un bloc a un etat."""

    model_name = "Repere"
    ports = ()
    params = ("wb",)
    states = ("theta",)
    svc_inputs = ("w_in",)
    svc_outputs = ("theta_out", "w")

    def __init__(self, name: str = "Repere", wb: float = 100.0 * 3.141592653589793):
        super().__init__(name, [], wb=wb)
        #: valeur de l'angle au moment ou le repere prend la main
        self.theta0 = 0.0

    @staticmethod
    def equations(s):
        return [s.wb * (s.w_in - s.w_ref)], [], [s.theta, s.w_in]

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        return {"theta": float(self.theta0), "w_in": float(ctx.w)}
