"""Groupe electrogene de secours : moteur diesel + arbre + alternateur synchrone.

L'ELECTRIQUE EST CELUI DE ``SynchronousMachine`` (stator transitoire, rotor a
trois enroulements d'axe d et deux d'axe q, AVR a quatre etats) : ses 17 etats
gardent leur nom et leur ordre, et rien de l'electrique n'est reecrit. Ce qui
change est la MECANIQUE, que la machine standard ramenait a une masse et a un
regulateur de turbine de reseau. Un groupe diesel a une chaine d'entrainement
a part entiere, et c'est elle qui decide de ce que subit l'arbre ::

    consigne de vitesse        regulateur PI + statisme
          |                          |  u
          v                          v
    w_moteur  <--- masse moteur    actuateur (cremaillere)  ---> fuel   [butees Fmin/Fmax]
                       |                                          |
                       |                              retard d'injection (Pade 1er ordre)
                       |                                          |
                 couple d'arbre  <-- torsion                combustion + turbo --> Pm
                       |
                  masse alternateur  ---> couple electrique Ce

DEUX MASSES -- la masse du moteur (volant compris, inertie ``He``) et celle du
rotor de l'alternateur (``H``, la valeur de ``SynchronousMachine``) sont
reliees par un accouplement elastique de raideur ``Ksh`` et d'amortissement
``Dsh`` ::

    2 He . dDwe/dt  = Cmot - Carbre - KDe.Dwe
    2 H  . dDwr/dt  = Carbre - Ce - KD.Dwr
    d(dsh)/dt       = wb.(Dwe - Dwr)
    Carbre          = Ksh.dsh + Dsh.(Dwe - Dwr)

``Carbre`` est le couple TRANSMIS par l'accouplement : c'est la grandeur qui
fatigue l'arbre, et le constructeur la fournit comme sortie de service
(``C_arbre``, pu de la machine) pour :mod:`g2egrid.fatigue`. La raideur se donne
plus naturellement par la frequence du mode de torsion, ``f_tors`` (Hz), et son
amortissement relatif ``zeta_tors`` : le constructeur en tire ``Ksh`` et
``Dsh`` ::

    wn^2 = (wb.Ksh/2).(1/He + 1/H) ,   Dsh = 2.zeta.wn / (1/(2He) + 1/(2H))

REGULATEUR -- vitesse mesuree au volant, PI avec statisme permanent ``Rd``
autour de la consigne de puissance ``p_ref`` ::

    e = w_set - w_moteur - Rd.(fuel - p_ref) ,   u = Kp_g.e + xg ,   dxg/dt = Ki_g.e

``Rd = 0`` est le mode isochrone d'un groupe seul sur son reseau ; en parallele
avec un reseau rigide il faut un statisme (3 a 5 %), sinon la consigne de
vitesse et la frequence du reseau se battent.

BUTEES DE CARBURANT -- la cremaillere ne depasse pas ``Fmax`` (la puissance
maximale du moteur, pas celle de l'alternateur) ni ``Fmin``. C'est un MODE
discret, ``mode_fuel`` (cf. :mod:`g2egrid.elements.limiteurs`) : en butee,
l'actuateur est fige et l'integrateur du regulateur suit la butee par
retro-calcul (``Taw_g``), sans quoi il s'emballerait. A drapeau nul les
equations sont celles du regulateur libre ; le modele linearise est celui du
point de fonctionnement NON sature. C'est le phenomene qui compte pour une
reprise de charge brutale : le moteur ne peut pas donner plus que ``Fmax``, la
vitesse chute, et l'arbre voit un couple moteur plafonne.

LIMITES ASSUMEES -- le modele est MOYEN : ni pulsation de couple a la
frequence d'allumage des cylindres, ni harmonique 6 du redresseur dans le
couple electrique (le repere dq ne voit que le fondamental). Ces excitations
se rajoutent par l'entree ``dT_exc`` (couple moteur additionnel, pu), a
l'amplitude et a la frequence que l'on veut tester.

``Pm`` n'est plus la sortie d'un regulateur de turbine mais la puissance
mecanique du moteur apres combustion ; la constante ``mp``/``Tg`` de la
machine standard ne sert plus.
"""

from __future__ import annotations

import numpy as np

from g2egrid.elements.sync_machine import SM_DEFAULT, SynchronousMachine
from g2egrid.noyau.element import InitContext

__all__ = ["GROUPE_DEFAULT", "GroupeElectrogene"]

#: Jeu de parametres d'un groupe de ~2,5 MVA (alternateur 400 V, 1500 tr/min,
#: moteur diesel suralimente). Ce sont des ORDRES DE GRANDEUR de materiel de ce
#: type, pas les donnees d'un constructeur : a remplacer par la fiche du groupe
#: etudie. Les reactances de l'alternateur : Xd = 1,8 ; X'd ~ 0,25 ; X''d ~ 0,17.
GROUPE_DEFAULT = dict(
    SM_DEFAULT,
    # --- alternateur ---
    Ra=0.004,
    Ll=0.12,
    Lad=1.68,
    Laq=1.0,
    Lfd=0.141,
    L1d=0.081,
    L1q=0.5,
    L2q=0.25,
    Rfd=0.002,
    R1d=0.02,
    R1q=0.02,
    R2q=0.03,
    H=0.5,
    KD=0.5,
    # --- AVR (excitatrice a diodes tournantes, plus lente qu'une thyristor) ---
    Ka=200.0,
    Ta=0.01,
    K_PSS=0.0,
    # --- moteur diesel et chaine d'entrainement ---
    He=0.8,
    KDe=0.3,
    Rd=0.03,
    Kp_g=20.0,
    Ki_g=20.0,
    Ta_g=0.05,
    Td_c=0.01,
    Tb_c=0.03,
    Fmax=0.88,
    Fmin=0.0,
    Taw_g=0.1,
    mode_fuel=0.0,
)

_MOTEUR = (
    "He",
    "KDe",
    "Ksh",
    "Dsh",
    "Rd",
    "Kp_g",
    "Ki_g",
    "Ta_g",
    "Td_c",
    "Tb_c",
    "Fmax",
    "Fmin",
    "Taw_g",
    "mode_fuel",
)


class GroupeElectrogene(SynchronousMachine):
    """Alternateur synchrone entraine par un moteur diesel via un arbre elastique.

    22 etats : les 17 de ``SynchronousMachine`` plus ``Dwe`` (vitesse du
    moteur), ``dsh`` (torsion de l'arbre, rad electriques), ``xg``
    (integrateur du regulateur), ``fuel`` (cremaillere), ``xpd`` (retard
    d'injection). Entrees de service : ``p_ref`` (consigne de puissance du
    statisme), ``v_ref`` (consigne de tension de l'AVR), ``w_set`` (consigne de
    vitesse) et ``dT_exc`` (couple d'excitation additionnel, pu).
    """

    params = SynchronousMachine.params + _MOTEUR
    states = SynchronousMachine.states + ("Dwe", "dsh", "xg", "fuel", "xpd")
    svc_inputs = SynchronousMachine.svc_inputs + ("w_set", "dT_exc")
    svc_outputs = SynchronousMachine.svc_outputs + ("C_arbre", "w_moteur", "C_moteur")
    modes = SynchronousMachine.modes + ("mode_fuel",)

    def __init__(
        self,
        name: str,
        bus: int,
        Sn: float,
        Sb: float,
        wb: float,
        f_tors: float = 25.0,
        zeta_tors: float = 0.03,
        **params,
    ):
        """``f_tors`` (Hz) et ``zeta_tors`` fixent la raideur et l'amortissement
        de l'accouplement, sauf si ``Ksh`` / ``Dsh`` sont donnes. Les autres
        parametres reprennent :data:`GROUPE_DEFAULT`."""
        p = dict(GROUPE_DEFAULT)
        p.update(params)
        if p["He"] <= 0.0 or p["H"] <= 0.0:
            raise ValueError(f"{name}: inerties H et He doivent etre > 0")
        for k in ("Ta_g", "Td_c", "Tb_c", "Taw_g"):
            if p[k] <= 0.0:
                raise ValueError(f"{name}: {k} = {p[k]} ; constante de temps <= 0")
        if p["Fmax"] <= p["Fmin"]:
            raise ValueError(f"{name}: Fmax doit depasser Fmin")
        wn = 2.0 * np.pi * float(f_tors)
        if "Ksh" not in p:
            p["Ksh"] = 2.0 * wn**2 / (wb * (1.0 / p["He"] + 1.0 / p["H"]))
        if "Dsh" not in p:
            a = 1.0 / (2.0 * p["He"]) + 1.0 / (2.0 * p["H"])
            p["Dsh"] = 2.0 * float(zeta_tors) * wn / a
        super().__init__(name, bus, Sn=Sn, Sb=Sb, wb=wb, **p)
        self.f_tors, self.zeta_tors = float(f_tors), float(zeta_tors)

    # ------------------------------------------------------------------
    def gardes(self, cm, x, z, u, p) -> list[tuple[str, float, float, int]]:
        """Butees de la cremaillere, en plus de celles de l'AVR.

        Libre : on bute quand ``fuel`` atteint ``Fmax`` en montant ou ``Fmin``
        en descendant. En butee : on repart quand la derivee LIBRE de ``fuel``
        change de signe, c'est-a-dire quand le regulateur redemande une valeur
        a l'interieur de la plage.
        """
        out = super().gardes(cm, x, z, u, p)
        params = self.params
        j = params.index("mode_fuel")
        k = self.states.index("fuel")
        mode = p[j]
        if mode == 0.0:
            out.append(("mode_fuel", x[k] - p[params.index("Fmax")], 1.0, +1))
            out.append(("mode_fuel", x[k] - p[params.index("Fmin")], -1.0, -1))
            return out
        libre = np.array(p, dtype=float)
        libre[j] = 0.0
        dfuel = cm.f(x, z, u, libre)[k]
        out.append(("mode_fuel", dfuel, 0.0, -1 if mode > 0 else +1))
        return out

    # ------------------------------------------------------------------
    @staticmethod
    def equations(s):
        f, g, h = SynchronousMachine.equations(s)

        wre = s.Dwe + s.w_ref  # vitesse du moteur (volant)
        Ce = s.phi_d * s.i_q - s.phi_q * s.i_d
        dw = s.Dwe - s.Dwr
        C_arbre = s.Ksh * s.dsh + s.Dsh * dw
        C_mot = s.Pm / wre + s.dT_exc

        # regulateur de vitesse, statisme permanent autour de p_ref
        e = s.w_set - wre - s.Rd * (s.fuel - s.p_ref)
        u = s.Kp_g * e + s.xg
        # butee de cremaillere (mode_fuel = +1 haut, -1 bas, 0 libre)
        m = s.mode_fuel
        F_lim = s.Fmax * (m + m**2) / 2 + s.Fmin * (m**2 - m) / 2
        dxg = s.Ki_g * e + m**2 * (F_lim - u) / s.Taw_g
        dfuel = (1 - m**2) * (u - s.fuel) / s.Ta_g
        dxpd = 2.0 * (s.fuel - s.xpd) / s.Td_c  # Pade 1er ordre du retard
        y_inj = 2.0 * s.xpd - s.fuel

        # mecanique : masse alternateur et combustion (indices de SM : Dwr, Pm)
        f[6] = (C_arbre - Ce - s.KD * s.Dwr) / (2 * s.H)
        f[7] = (y_inj - s.Pm) / s.Tb_c
        f = f + [
            (C_mot - C_arbre - s.KDe * s.Dwe) / (2 * s.He),  # Dwe
            s.wb * dw,  # dsh
            dxg,
            dfuel,
            dxpd,
        ]
        h = list(h) + [C_arbre, wre, C_mot]
        return f, g, h

    # ------------------------------------------------------------------
    def initialize(self, ctx: InitContext) -> dict[str, float]:
        init = super().initialize(ctx)
        p = self.pdict
        Ce = init["phi_d"] * init["i_q"] - init["phi_q"] * init["i_d"]
        Pm = init["Pm"]
        if Pm > p["Fmax"] or Pm < p["Fmin"]:
            raise ValueError(
                f"{self.name}: la puissance du point de fonctionnement "
                f"({Pm:.4g} pu machine) sort de la plage du moteur "
                f"[{p['Fmin']:g}, {p['Fmax']:g}]"
            )
        init.update(
            Dwe=0.0,
            dsh=Ce / p["Ksh"],
            xg=Pm,
            fuel=Pm,
            xpd=Pm,
            w_set=ctx.w,
            dT_exc=0.0,
        )
        return init
