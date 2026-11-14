"""Reseau type d'un centre de donnees alimente en MT, avec groupe de secours.

Le site assemble les blocs de :mod:`g2egrid.elements.datacenter` et le groupe de
:mod:`g2egrid.elements.groupe_electrogene` ::

    GRID --- L_AMONT --- [2] ---+--- T_DC --- [4] RED === bus continu [10] === OND --- [5] IT
     (20 kV, IB)      (jeu de    |                          |  (capa 2 pu)             (charge a p_ref)
                      barres MT) |                         BAT (DC/DC + batterie)
                                 +--- T_GE --- [3] GE (diesel + alternateur)

* ``GRID`` : le reseau amont, un bus infini (le slack) ;
* ``L_AMONT`` : la liaison dont l'ouverture ISOLE le site (perte du reseau) ;
* ``T_DC`` : transformateur MT/BT du centre de donnees ; bus ``[4]`` : le
  condensateur alternatif du redresseur (capacite du bus) ;
* ``RED`` : redresseur a diodes ; ``[10]`` : bus continu, dont la capacite est
  celle du condensateur du bus ; ``BAT`` : DC/DC et batterie ;
* ``OND`` : onduleur V/f constants, slack de l'ilot alternatif ``[5]`` ;
  ``IT`` : la charge informatique, pilotee par ``IT.p_ref`` ;
* ``T_GE`` / ``GE`` : le groupe electrogene de secours, en attente a chaud sur
  le jeu de barres MT (a vide, synchronise sur le reseau).

Base systeme ``Sb = 5 MVA`` ; la charge nominale est de 1,5 MW (``P_IT = 0,30``
pu), le groupe de 2,5 MVA (2,2 MW de moteur : ``Fmax = 0,88``) : il peut porter toute la charge,
avec une reserve de 30 %.

C'EST UN REGLAGE TYPE, pas un site : les valeurs sont des ordres de grandeur
coherents entre eux (tenue des tensions, partage des courants, stabilite
verifiee dans ``tests/test_datacenter.py``), a remplacer par les donnees d'un
site reel.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from g2egrid.elements import (
    IB,
    BatterieDCDC,
    BusDC,
    ChargeACPref,
    GroupeElectrogene,
    Line,
    OnduleurVF,
    RedresseurDiode,
    Transformer,
)
from g2egrid.elements.datacenter import reactif_commutation
from g2egrid.noyau.loadflow import DC_SLACK
from g2egrid.reseau.network import Network

__all__ = [
    "NOMS",
    "courbe_charge_ia",
    "reseau_datacenter",
    "scenario_perte_reseau",
]

#: Les noms des elements du site, pour s'y reporter sans les retaper.
NOMS = {
    "reseau": "GRID",
    "liaison": "L_AMONT",
    "groupe": "GE",
    "transfo_groupe": "T_GE",
    "transfo_site": "T_DC",
    "redresseur": "RED",
    "batterie": "BAT",
    "onduleur": "OND",
    "charge": "IT",
}


def reseau_datacenter(
    P_it: float = 0.30,
    cos_phi: float = 0.97,
    Sb: float = 5e6,
    avec_groupe: bool = True,
    avec_batterie: bool = True,
    reseau_amont: bool = True,
    lissage: float = 0.0,
    Kv_batterie: float = 20.0,
    Eb: float = 180.0,
    C_dc: float = 3.0,
    C_ac: float = 0.01,
    C_mt: float = 0.015,
    P_groupe_attente: float = 0.01,
    resoudre: bool = True,
    **params_groupe,
) -> Network:
    """Construit (et par defaut resout) le reseau type d'un centre de donnees.

    ``P_it`` : charge informatique au point de depart [pu de ``Sb``] ;
    ``cos_phi`` : son facteur de puissance. ``avec_groupe`` / ``avec_batterie``
    : retirent le groupe ou la batterie (comparaisons de configuration).
    ``reseau_amont=False`` : le site est ISOLE, le groupe est le slack et porte
    toute la charge ; c'est l'etat de regime apres la perte du reseau, autour
    duquel le modele linearise decrit le site alimente par le groupe.
    ``lissage`` : gain ``Kh`` du lissage rapide de la batterie (0 = aucun).
    ``Kv_batterie`` : statisme de tension de la batterie [pu courant/pu tension].
    ``Eb`` : energie utile de la batterie [pu.s] (180 : dix minutes a la charge
    nominale). ``C_ac`` : condensateur alternatif du redresseur [pu] (0,01 : 50 kvar ; un
    condensateur trop gros fait monter la tension d'un groupe a vide). ``C_mt`` : capacite des jeux de barres MT [pu] (0,015 : la charge d'un cable de
    5 km ; une capacite trop faible rend le circuit liaison-source resonant a plus de
    20 kHz, ce que l'integrateur paie tres cher a chaque manoeuvre). ``C_dc`` : condensateur du bus continu [pu] (3 : environ 20 ms de
    charge nominale). ``P_groupe_attente`` : puissance du groupe en attente
    [pu de ``Sb``]. ``params_groupe`` : parametres de
    :class:`~g2egrid.elements.GroupeElectrogene` (``f_tors``, ``Rd``...).

    ``resoudre`` lance le load flow mixte et l'equilibre ; il rend un reseau
    prets pour :func:`g2egrid.temporel.ModeleNonLineaire` et
    ``build_state_space``.
    """
    net = Network(f0=50.0, Sb=Sb)
    wb = net.wb
    tphi = float(np.tan(np.arccos(cos_phi)))
    if not reseau_amont and not avec_groupe:
        raise ValueError("un site isole a besoin du groupe electrogene")
    if reseau_amont:
        net.add_bus(1, 20e3, btype=1, V=1.0, C=C_mt)
    net.add_bus(2, 20e3, btype=3, C=C_mt)
    net.add_bus(4, 400.0, btype=3, C=C_ac)  # condensateur alternatif du redresseur
    net.add_bus(5, 400.0, btype=1, V=1.0, C=5e-2)  # ilot alternatif de l'onduleur
    net.add(BusDC(10, 700.0, C=C_dc, V=1.0))
    net.bus_specs_dc[10].btype = DC_SLACK
    if reseau_amont:
        net.add(IB(NOMS["reseau"], 1, wb=wb))
        net.add(Line(NOMS["liaison"], 1, 2, R=0.005, X=0.02, wb=wb))
    net.add(Transformer(NOMS["transfo_site"], 2, 4, R=0.004, X=0.14, wb=wb))
    net.add(RedresseurDiode(NOMS["redresseur"], 4, 10, Sn=2.0e6, Sb=Sb, wb=wb))
    net.add(OnduleurVF(NOMS["onduleur"], 5, 10, Sn=2.0e6, Sb=Sb, wb=wb))
    net.add(ChargeACPref(NOMS["charge"], 5, P=P_it, Q=P_it * tphi))
    net.vsc_modes[4] = "VdcQ"  # le pont tient la tension continue et son reactif
    net.vsc_modes[5] = "Vtheta"  # l'onduleur forme son ilot
    if avec_batterie:
        net.add(
            BatterieDCDC(
                NOMS["batterie"], 10, Kv=Kv_batterie, Kh=lissage, Eb=Eb, vref=0.92
            )
        )
    if avec_groupe:
        if reseau_amont:
            net.add_bus(3, 400.0, btype=2, V=1.0, Pg=P_groupe_attente, C=C_mt)
        else:
            net.add_bus(3, 400.0, btype=1, V=1.0, C=C_mt)  # le groupe est le slack
        net.add(Transformer(NOMS["transfo_groupe"], 3, 2, R=0.004, X=0.12, wb=wb))
        net.add(
            GroupeElectrogene(
                NOMS["groupe"], 3, Sn=2.5e6, Sb=Sb, wb=wb, **params_groupe
            )
        )
    # le repere dq global est porte par la source qui tient l'ilot MT : le
    # reseau amont, ou le groupe quand le site est isole. Sans cela, "le premier
    # bus slack" serait celui de l'onduleur (l'autre ilot), et le groupe
    # tournerait dans un repere qui n'est pas le sien.
    net.reference = NOMS["reseau"] if reseau_amont else NOMS["groupe"]
    if resoudre:
        net.run_loadflow()
        net.solve_equilibrium()
        # Le pont consomme le reactif de sa commutation (empietement) : on
        # l'impose au load flow, puis on refait l'equilibre. Deux passes
        # suffisent, le reactif ne change presque pas la puissance active.
        red = next(e for e in net.elements if e.name == NOMS["redresseur"])
        for _ in range(2):
            net.bus_specs[4].Ql = reactif_commutation(red)
            net.run_loadflow()
            net.solve_equilibrium()
    return net


# ----------------------------------------------------------------------
# courbe de charge
# ----------------------------------------------------------------------
def courbe_charge_ia(
    t,
    P_moy: float = 0.30,
    amplitude: float = 0.25,
    periode: float = 4.0,
    t_debut: float = 20.0,
    rampe: float = 15.0,
    intervalle_ckpt: float = 120.0,
    duree_ckpt: float = 6.0,
    profondeur_ckpt: float = 0.6,
    bruit: float = 0.01,
    graine: int = 0,
) -> np.ndarray:
    """Une courbe de puissance d'une grappe d'entrainement d'IA [pu de ``Sb``].

    Ce n'est pas une mesure : c'est un profil SYNTHETIQUE qui reproduit les
    trois traits que les exploitants decrivent pour ces charges :

    * a l'arret, la charge est stable a ``P_moy`` ;
    * l'entrainement demarre a ``t_debut`` et monte en ``rampe`` secondes vers un
      regime ALTERNE : a chaque iteration (``periode``) la grappe passe d'une
      phase de calcul (puissance haute) a une phase de communication/attente
      (puissance basse), en creneaux aux flancs raides de ``P_moy.(1 +- amplitude)`` ;
    * toutes les ``intervalle_ckpt`` secondes, une sauvegarde (checkpoint) fait
      tomber la charge a ``(1 - profondeur_ckpt)`` de sa valeur pendant
      ``duree_ckpt`` secondes ;
    * un bruit gaussien filtre de ``bruit`` (relatif) s'y superpose.

    La courbe vaut exactement ``P_moy`` en ``t[0]`` pour partir de l'equilibre.
    Pour utiliser une mesure, passer directement ses valeurs a
    :class:`g2egrid.temporel.Courbe`.
    """
    t = np.asarray(t, dtype=float)
    # creneau lisse de periode `periode`, apres le demarrage
    tau = np.maximum(t - t_debut, 0.0)
    env = np.clip(tau / max(rampe, 1e-9), 0.0, 1.0)
    creneau = np.tanh(3.0 * np.sin(2.0 * np.pi * tau / periode)) / np.tanh(3.0)
    forme = 1.0 + env * amplitude * creneau
    # sauvegardes periodiques : un creneau de duree `duree_ckpt`, flancs de 0,3 s,
    # qui commence une seconde apres le debut de chaque intervalle
    if intervalle_ckpt > 0:
        phase = np.mod(t - t_debut - rampe, intervalle_ckpt)
        fenetre = 0.5 * (
            np.tanh((phase - 1.0) / 0.3) - np.tanh((phase - 1.0 - duree_ckpt) / 0.3)
        )
        forme = forme * (1.0 - profondeur_ckpt * fenetre * (t >= t_debut + rampe))
    # bruit gaussien filtre (moyenne glissante sur ~0,5 s), d'ecart-type `bruit`
    rng = np.random.default_rng(graine)
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 1.0
    k = max(1, round(0.5 / dt))
    bruit_f = np.convolve(rng.standard_normal(len(t)), np.ones(k) / k, mode="same")
    bruit_f = bruit_f / max(bruit_f.std(), 1e-12) * bruit
    return P_moy * forme * (1.0 + env * bruit_f)


# ----------------------------------------------------------------------
# scenario
# ----------------------------------------------------------------------
@dataclass
class ScenarioPerteReseau:
    """Le scenario rendu par :func:`scenario_perte_reseau` : la fonction
    d'entrees et ses dates, pour les annoter sur une figure."""

    uext: object
    t_perte: float
    t_reprise: float
    duree_reprise: float
    v_batterie: np.ndarray


def scenario_perte_reseau(
    m,
    t,
    P_it,
    t_perte: float = 60.0,
    t_reprise: float = 75.0,
    duree_reprise: float = 60.0,
    v_batterie: float = 1.06,
    v_veille: float = 0.92,
    anticipation: float = 1.0,
) -> ScenarioPerteReseau:
    """Le scenario de reference : perte du reseau amont et reprise par le groupe.

    * ``t < t_perte - anticipation`` : le reseau alimente le site ; la batterie
      est en veille (``BAT.v_ref = v_veille``, sous la tension du bus continu) ;
    * ``t_perte - anticipation`` : la batterie releve sa consigne a
      ``v_batterie`` (en 20 ms), au-dessus de la tension de sortie du pont : les
      diodes se bloquent et la batterie porte seule la charge ;
    * ``t_perte`` : la liaison ``L_AMONT`` s'ouvre, le site est isole. Le groupe,
      deja en route et couple en attente, n'a plus que le reseau a tenir ;
    * ``t_reprise`` : la consigne de la batterie redescend lineairement a
      ``v_veille`` en ``duree_reprise`` : le redresseur, alimente par le groupe,
      reprend la charge PROGRESSIVEMENT (c'est la reprise de charge echelonnee d'un
      groupe de secours) ;
    * pendant tout le scenario la charge suit la courbe ``P_it``.

    POURQUOI LA BATTERIE PRECEDE L'OUVERTURE -- dans un site reel, la batterie
    prend la charge en quelques millisecondes quand le reseau disparait, et le
    groupe n'est couple que dix secondes apres, a vide. Le modele n'a pas de
    manoeuvre de fermeture : le groupe est couple des le depart, en attente. Si
    la liaison s'ouvrait AVANT que la batterie ait pris la charge, le groupe
    recevrait d'un coup les 0,6 pu de la charge, sa tension plongerait, et le
    pont de diodes se bloquerait puis reconduirait toutes les quelques
    millisecondes -- une oscillation qu'aucun site ne connait, puisque le groupe
    n'est pas la a cet instant. ``anticipation`` (1 s) transpose donc l'ordre
    reel : batterie d'abord, groupe ensuite. Mettre 0 retrouve l'ouverture seche
    (a ne faire que pour etudier ce choc : la simulation y est tres couteuse).

    ``m`` : un :class:`~g2egrid.temporel.ModeleNonLineaire`. Rend le scenario
    (:class:`ScenarioPerteReseau`), dont ``uext`` se passe a ``simuler``.
    """
    from g2egrid.temporel import Courbe, Manoeuvre

    t = np.asarray(t, dtype=float)
    t_batt = t_perte - anticipation
    if t_batt <= t[0]:
        raise ValueError("t_perte - anticipation doit suivre le debut de la simulation")
    if t_reprise <= t_perte:
        raise ValueError("t_reprise doit suivre t_perte")
    pts_t = [t[0], t_batt, t_batt + 0.02, t_reprise, t_reprise + duree_reprise]
    pts_v = [v_veille, v_veille, v_batterie, v_batterie, v_veille]
    v_ref = np.interp(t, pts_t, pts_v)
    # seuls les points de cassure sont declares a l'integrateur : une courbe
    # lineaire sur toute la grille le ferait redemarrer a chaque echantillon
    uext = m.scenario(
        Courbe(f"{NOMS['charge']}.p_ref", t, P_it),
        Courbe(f"{NOMS['batterie']}.v_ref", pts_t, pts_v, interpolation="lineaire"),
        Manoeuvre(t_perte, "ouverture", NOMS["liaison"]),
    )
    return ScenarioPerteReseau(uext, t_perte, t_reprise, duree_reprise, v_ref)
