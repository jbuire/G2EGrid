"""Machine synchrone : portage FIDELE de ``matlab/G2Elin_V2/core/SM.m``.

Meme structure, memes equations, memes valeurs par defaut que le modele
MATLAB. Toute modification ici doit etre repercutee la-bas (et inversement) :
la table de correspondance est dans ``docs/equivalence_matlab.md``.

ETATS (17) -- ordre et noms identiques a SM.m
    phi_d, phi_q                 flux stator
    phi_fd, phi_1d               rotor axe d (champ + amortisseur)
    phi_1q, phi_2q               rotor axe q (deux amortisseurs)
    Dwr                          ecart de vitesse (pu)
    Pm                           puissance mecanique (regulateur de vitesse)
    Dw1, v1, v2, vpss            PSS (filtre + washout + 2 avances de phase)
    e1, e2, efd, e3              AVR (mesure, ampli, excitatrice, stabilisation)
    theta                        angle rotor / repere dq global

ALGEBRIQUES (6) : i_d, i_q, i_fd, i_1d, i_1q, i_2q
    relations flux/courant (matrice d'inductances). MATLAB les elimine
    symboliquement en inversant Ld_mat/Lq_mat ; les garder algebriques est
    strictement equivalent (complement de Schur) et numeriquement plus sain.

RESISTANCE DE CAUSALITE ``Rg``
    SM.m raccorde la machine au noeud a travers une resistance shunt
    ``Rg = 100`` pu : le courant INJECTE vaut ``ig = i - v/Rg`` et c'est ``ig``
    qui porte la chute ``Ra`` dans les equations statoriques. Le terme est
    reproduit ici a l'identique. Le mettre a une tres grande valeur
    (``Rg=1e9``) supprime ce shunt non physique.

REPERE
    L'axe d rotor est en avance de ``theta`` sur l'axe D global ::

        ved + j.veq = (vd + j.vq) . exp(-j.theta)
        iD  + j.iQ  = (igd + j.igq) . exp(+j.theta)

    MATLAB ecrit la meme rotation avec ``rot = theta_g - theta`` parce que son
    interface presente la tension de bus sous la forme (V, 0) + angle separe.

BASES
    Parametres en pu base MACHINE (Sn) ; le courant injecte est converti en pu
    base SYSTEME par ``kS = Sn/Sb`` (bases de tension identiques).
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from g2egrid.noyau.element import Element, InitContext

__all__ = ["SynchronousMachine", "SynchronousMachineTR", "SM_DEFAULT"]

# Valeurs par defaut : reprise LITTERALE du constructeur de matlab/G2Elin_V2/core/SM.m.
SM_DEFAULT = dict(
    # --- Electrical (jeu de Kundur, ex. 3.1 / 4.1) ---
    Ra=0.003,
    Ll=0.15,
    Lad=1.66,
    Laq=1.61,
    Lfd=0.165,
    L1d=0.1713,
    L1q=0.7252,
    L2q=0.125,
    Rfd=0.0006,
    R1d=0.0284,
    R1q=0.00619,
    R2q=0.02368,
    Rg=100.0,  # resistance de causalite aux bornes
    # --- Mechanical ---
    H=5.0,
    KD=1.0,
    # --- Governor (mp en POURCENT, comme SM.m) ---
    mp=0.5,
    Tg=0.2,
    # --- AVR ---
    Tr=0.02,
    Ka=300.0,
    Ta=0.001,
    Ke=1.0,
    Te=0.0001,
    Kfd=0.001,
    Tfd=0.1,
    # --- PSS ---
    K_PSS=20.0,
    T_LP=0.03,
    T_HP=2.0,
    T1n=0.05,
    T1d=0.02,
    T2n=3.0,
    T2d=5.4,
    # butees de l'amplificateur de l'AVR (e2, pu) et leur drapeau -- voir
    # g2egrid.elements.limiteurs ; a mode_e2 = 0 les equations sont celles de
    # SM.m, au bit pres
    VRmax=5.0,
    VRmin=-5.0,
    mode_e2=0.0,
)

# Ancien jeu de parametres du portage (conserve pour reference historique).
KUNDUR = dict(SM_DEFAULT)


class SynchronousMachine(Element):
    """Machine synchrone 6 enroulements + AVR 4 etats + PSS + regulateur."""

    ports = ("t",)
    params = (
        "Ra",
        "Ll",
        "Lad",
        "Laq",
        "Lfd",
        "L1d",
        "L1q",
        "L2q",
        "Rfd",
        "R1d",
        "R1q",
        "R2q",
        "Rg",
        "H",
        "KD",
        "mp",
        "Tg",
        "Tr",
        "Ka",
        "Ta",
        "Ke",
        "Te",
        "Kfd",
        "Tfd",
        "K_PSS",
        "T_LP",
        "T_HP",
        "T1n",
        "T1d",
        "T2n",
        "T2d",
        "kS",
        "wb",
        "theta_g0",
        "VRmax",
        "VRmin",
        "mode_e2",
    )
    states = (
        "phi_d",
        "phi_q",
        "phi_fd",
        "phi_1d",
        "phi_1q",
        "phi_2q",
        "Dwr",
        "Pm",
        "Dw1",
        "v1",
        "v2",
        "vpss",
        "e1",
        "e2",
        "efd",
        "e3",
        "theta",
    )
    algebs = ("i_d", "i_q", "i_fd", "i_1d", "i_1q", "i_2q")
    svc_inputs = ("p_ref", "v_ref")
    svc_outputs = (
        "P",
        "Q",
        "V",
        "Omega",
        "ved",
        "veq",
        "igd",
        "igq",
        "Ce",
        "theta_out",
    )

    def __init__(self, name: str, bus: int, Sn: float, Sb: float, wb: float, **params):
        p = dict(SM_DEFAULT)
        p.update(params)
        p["kS"] = Sn / Sb
        p["wb"] = wb
        p.setdefault("theta_g0", 0.0)  # fixe par initialize
        super().__init__(name, bus, **p)
        self.Sn, self.Sb = float(Sn), float(Sb)

    modes = ("mode_e2",)

    def gardes(self, cm, x, z, u, p) -> list[tuple[str, float, float, int]]:
        """Butees de l'amplificateur de l'AVR (cf. g2egrid.elements.limiteurs).

        Libre : on bute quand ``e2`` atteint ``VRmax`` en montant, ou ``VRmin``
        en descendant. En butee : on repart quand la derivee LIBRE de ``e2``
        change de signe -- l'amplificateur voudrait revenir vers l'interieur.
        """
        params = self.params
        j = params.index("mode_e2")
        k = self.states.index("e2")
        mode = p[j]
        if mode == 0.0:
            return [
                ("mode_e2", x[k] - p[params.index("VRmax")], 1.0, +1),
                ("mode_e2", x[k] - p[params.index("VRmin")], -1.0, -1),
            ]
        libre = np.array(p, dtype=float)
        libre[j] = 0.0
        de2 = cm.f(x, z, u, libre)[k]
        return [("mode_e2", de2, 0.0, -1 if mode > 0 else +1)]

    # Pas de shunt_admittance : Rg est interne au modele, dont le courant
    # injecte vaut deja ig = i - v/Rg. Le declarer aussi au Y-bus le comptait
    # deux fois -- la loi des noeuds etait violee de V/Rg au bus de chaque
    # machine et le spectre dependait du choix de la machine slack.

    # ------------------------------------------------------------------
    @staticmethod
    def equations(s):
        Ld, Lq = s.Lad + s.Ll, s.Laq + s.Ll
        Lffd, L11d = s.Lad + s.Lfd, s.Lad + s.L1d
        L11q, L22q = s.Laq + s.L1q, s.Laq + s.L2q
        mp = s.mp / 100  # SM.m : Governor.mp est en %

        # --- repere RESEAU -> repere ROTOR ---
        # Repere du rotor par rapport au repere global, comme tous les autres
        # equipements : theta - theta_g, decale de theta_g0 (sa valeur
        # d'equilibre) pour que l'equilibre reste celui du load flow. Pour la
        # machine de reference, theta_g EST theta : la rotation est constante et
        # le repere global suit son rotor, comme w_slack dans G2Elin V2.
        rot = s.theta - s.theta_g + s.theta_g0
        ct, st = sp.cos(rot), sp.sin(rot)
        ved = s.vd_t * ct + s.vq_t * st
        veq = -s.vd_t * st + s.vq_t * ct
        Et = sp.sqrt(s.vd_t**2 + s.vq_t**2)

        # --- courant INJECTE : la resistance de causalite Rg est aux bornes ---
        igd = s.i_d - ved / s.Rg
        igq = s.i_q - veq / s.Rg

        wr = s.Dwr + s.w_ref
        Ce = s.phi_d * s.i_q - s.phi_q * s.i_d
        Cm = s.Pm / wr
        DP = (wr - s.w_ref) / mp  # statisme du regulateur de vitesse

        # --- PSS : filtre passe-bas, washout, deux avances de phase ---
        dDw1 = (s.Dwr - s.Dw1) / s.T_LP
        dv1 = (s.T_HP * s.K_PSS * dDw1 - s.v1) / s.T_HP
        dv2 = (s.T1n * dv1 + s.v1 - s.v2) / s.T1d
        dvpss = (s.T2n * dv2 + s.v2 - s.vpss) / s.T2d

        # --- AVR : mesure, amplificateur, excitatrice, stabilisation ---
        de1 = (Et - s.e1) / s.Tr
        # en butee (mode_e2 = +1 ou -1), l'amplificateur est fige a sa limite
        de2 = (1 - s.mode_e2**2) * (
            (s.Ka * (s.v_ref - s.e1 - s.e3 + s.vpss) - s.e2) / s.Ta
        )
        defd = (s.Ke * s.e2 - s.efd) / s.Te
        de3 = (s.Kfd * defd - s.e3) / s.Tfd

        f = [
            # --- stator (transitoires inclus) ---
            s.wb * (ved + s.Ra * igd + wr * s.phi_q),
            s.wb * (veq + s.Ra * igq - wr * s.phi_d),
            # --- rotor ---
            (s.wb * s.Rfd / s.Lad) * s.efd - s.wb * s.Rfd * s.i_fd,
            s.wb * (-s.R1d * s.i_1d),
            s.wb * (-s.R1q * s.i_1q),
            s.wb * (-s.R2q * s.i_2q),
            # --- mecanique + regulateur de vitesse ---
            (Cm - Ce - s.KD * s.Dwr) / (2 * s.H),
            (s.p_ref - DP - s.Pm) / s.Tg,
            # --- PSS ---
            dDw1,
            dv1,
            dv2,
            dvpss,
            # --- AVR ---
            de1,
            de2,
            defd,
            de3,
            # --- angle rotor ---
            s.wb * s.Dwr,
        ]

        # --- relations flux / courant (convention generateur) ---
        g = [
            -s.phi_d + (-Ld * s.i_d + s.Lad * s.i_fd + s.Lad * s.i_1d),
            -s.phi_q + (-Lq * s.i_q + s.Laq * s.i_1q + s.Laq * s.i_2q),
            -s.phi_fd + (-s.Lad * s.i_d + Lffd * s.i_fd + s.Lad * s.i_1d),
            -s.phi_1d + (-s.Lad * s.i_d + s.Lad * s.i_fd + L11d * s.i_1d),
            -s.phi_1q + (-s.Laq * s.i_q + L11q * s.i_1q + s.Laq * s.i_2q),
            -s.phi_2q + (-s.Laq * s.i_q + s.Laq * s.i_1q + L22q * s.i_2q),
        ]

        # --- courant injecte : repere ROTOR -> repere GLOBAL, base systeme ---
        iD = s.kS * (igd * ct - igq * st)
        iQ = s.kS * (igd * st + igq * ct)
        h = [
            iD,
            iQ,
            ved * igd + veq * igq,  # P injecte (pu base machine)
            -ved * igq + veq * igd,  # Q injecte
            Et,
            wr,
            ved,
            veq,
            igd,
            igq,
            Ce,
            s.theta,
        ]
        return f, g, h

    # ------------------------------------------------------------------
    def initialize(self, ctx: InitContext) -> dict[str, float]:
        self.set_param("theta_g0", ctx.theta_g)
        p = self.pdict
        Ld, Lq = p["Lad"] + p["Ll"], p["Laq"] + p["Ll"]
        Lffd = p["Lad"] + p["Lfd"]
        Rg, Ra = p["Rg"], p["Ra"]

        V = ctx.v["t"]  # tension bus (pu systeme)
        S = ctx.s.get("t", 0j)  # puissance INJECTEE (systeme)
        Ig = (np.conj(S / V) if abs(V) > 1e-12 else 0j) / p["kS"]  # pu machine

        # Axe q aligne sur E = V + Ra.Ig + j.Lq.I, avec I = Ig + V/Rg le courant
        # MACHINE (celui des equations de flux) ; l'axe d est a -90 deg.
        E = V * (1 + 1j * Lq / Rg) + complex(Ra, Lq) * Ig
        theta = float(np.angle(E)) - np.pi / 2.0

        rot = np.exp(-1j * theta)
        e = V * rot  # (ved, veq)
        ig = Ig * rot  # (igd, igq)
        i = ig + e / Rg  # courant machine
        ved, veq = e.real, e.imag
        igd, igq = ig.real, ig.imag
        i_d, i_q = i.real, i.imag

        # stator a l'equilibre (wr = 1) : dphi/dt = 0
        phi_q = -(ved + Ra * igd)
        phi_d = veq + Ra * igq
        # amortisseurs a courant nul en regime permanent
        i_1d = i_1q = i_2q = 0.0
        i_fd = (phi_d + Ld * i_d) / p["Lad"]
        phi_fd = -p["Lad"] * i_d + Lffd * i_fd
        phi_1d = -p["Lad"] * i_d + p["Lad"] * i_fd
        phi_1q = -p["Laq"] * i_q
        phi_2q = -p["Laq"] * i_q

        efd = p["Lad"] * i_fd
        Ce = phi_d * i_q - phi_q * i_d
        Pm = Ce * ctx.w
        Et = abs(V)

        # AVR : de1 = de2 = defd = de3 = 0  ->  e1 = Et, e2 = efd/Ke, e3 = 0
        e2 = efd / p["Ke"]
        # PSS : Dwr = 0 -> Dw1 = v1 = v2 = vpss = 0
        return {
            "phi_d": phi_d,
            "phi_q": phi_q,
            "phi_fd": phi_fd,
            "phi_1d": phi_1d,
            "phi_1q": phi_1q,
            "phi_2q": phi_2q,
            "Dwr": 0.0,
            "Pm": Pm,
            "Dw1": 0.0,
            "v1": 0.0,
            "v2": 0.0,
            "vpss": 0.0,
            "e1": Et,
            "e2": e2,
            "efd": efd,
            "e3": 0.0,
            "theta": theta,
            "i_d": i_d,
            "i_q": i_q,
            "i_fd": i_fd,
            "i_1d": i_1d,
            "i_1q": i_1q,
            "i_2q": i_2q,
            "p_ref": Pm,
            "v_ref": Et + e2 / p["Ka"],
        }


class SynchronousMachineTR(SynchronousMachine):
    """Machine synchrone AVEC son transfo de raccordement, fidele a ``symSG``.

    Difference avec :class:`SynchronousMachine` : le transformateur de
    raccordement ``Rt``/``Lt`` est ICI, dans la machine, et son courant
    ``i_gd``/``i_gq`` est un ETAT exprime dans le **repere rotorique**.

    C'est la structure de ``symSG`` (G2Elin V1, la reference) :

        digd = (wb/Lt).(ved - vgd - Rt.igd + wr.Lt.igq)
        digq = (wb/Lt).(veq - vgq - Rt.igq - wr.Lt.igd)

    Le port de l'element est donc le noeud **RESEAU**, et non la borne de la
    machine : la tension aux bornes devient interne, portee par le shunt de
    causalite ``Rg`` (``alg2``/``alg3`` de symSG) ::

        ved = (i_d - i_gd).Rg        veq = (i_q - i_gq).Rg

    Avec :class:`SynchronousMachine`, la meme physique se modelise en trois
    objets -- machine, noeud terminal resistif, puis ``Transformer`` -- dont le
    courant est exprime dans le repere RESEAU. Les deux descriptions sont
    equivalentes, mais leurs matrices d'etat different alors d'une rotation de
    l'angle rotor, ce qui empeche toute comparaison terme a terme avec V1.

    L'ordre des etats reprend celui de ``symSG`` (``diffeqVec``), pour que les
    matrices A des deux outils se lisent ligne a ligne.
    """

    params = SynchronousMachine.params + ("Rt", "Lt")
    # ``ved``/``veq`` sont ALGEBRIQUES, comme dans symSG (algVec contient ved et
    # veq). La difference n'est pas cosmetique : la reference substitue dans son
    # jacobien les valeurs d'equilibre ed_0/eq_0 issues de ses propres formules,
    # et non ce que rendrait (i_d - i_gd).Rg au meme point. Les deux different
    # ici de 2,8 %, et cela se lit directement sur la ligne du regulateur de
    # tension. Les garder algebriques reproduit ce comportement.
    #: (``v_ed``/``v_eq`` sont les ``ved``/``veq`` de symSG ; le nom porte un
    #: souligne parce que la classe expose deja ``ved``/``veq`` en SORTIE.)
    algebs = SynchronousMachine.algebs + ("v_ed", "v_eq")
    states = (
        "i_gd",
        "i_gq",
        "phi_d",
        "phi_q",
        "phi_fd",
        "phi_1d",
        "phi_1q",
        "phi_2q",
        "Dwr",
        "theta",
        "Pm",
        "Dw1",
        "v1",
        "v2",
        "vpss",
        "e1",
        "e2",
        "efd",
        "e3",
    )

    #: conventions d'initialisation acceptees par ``init``
    INITS = ("physique", "symSG")

    def __init__(
        self,
        name: str,
        bus: int,
        Sn: float,
        Sb: float,
        wb: float,
        init: str = "physique",
        **params,
    ):
        """``init`` choisit la convention de calcul du point d'equilibre.

        ``'physique'`` (defaut) remonte du noeud reseau vers la borne machine a
        travers le transfo, puis ajoute le shunt Rg : le courant du transfo vaut
        alors exactement l'injection du load flow.

        ``'symSG'`` reproduit litteralement ``symSG`` de G2Elin V1, qui prend la
        tension ET la puissance de SON bus pour celles de la borne machine, puis
        RETIRE le courant du shunt au courant du transfo. Le transfo ne porte
        donc que ``I - V/Rg``, soit 0,1 pu de moins que le courant que le load
        flow avait fait passer par cette branche. C'est une incoherence de la
        reference -- l'equilibre obtenu ne verifie pas les equations du modele --
        mais il faut la rejouer pour comparer les deux matrices A terme a terme.
        Dans ce mode, ``raffiner_equilibre`` passe a False.

        Le bus de la machine, au sens de V1, est celui que ``network_form.m``
        cree DERRIERE le transfo. Poser :attr:`point_init` sur la tension et la
        generation de ce bus est donc necessaire a la fidelite ; a defaut,
        l'element retombe sur le noeud reseau, ce qui deplace l'angle rotor.
        """
        if init not in self.INITS:
            raise ValueError(
                f"init={init!r} inconnu (attendu : {', '.join(self.INITS)})"
            )
        self.init_mode = init
        # L'equilibre de symSG ne verifie pas ses propres equations : le
        # raffiner par Newton s'en ecarterait (cf. raffiner_equilibre).
        self.raffiner_equilibre = init != "symSG"
        params.setdefault("Rt", 0.0002)
        params.setdefault("Lt", 0.024)
        super().__init__(name, bus, Sn=Sn, Sb=Sb, wb=wb, **params)

    @staticmethod
    def equations(s):
        Ld, Lq = s.Lad + s.Ll, s.Laq + s.Ll
        Lffd, L11d = s.Lad + s.Lfd, s.Lad + s.L1d
        L11q, L22q = s.Laq + s.L1q, s.Laq + s.L2q
        mp = s.mp / 100

        # --- tension du noeud RESEAU, ramenee au repere rotor (alg10/alg11) ---
        #   symSG : vgd = vgd_g.cos(theta_g - theta) - vgq_g.sin(theta_g - theta)
        #           vgq = vgd_g.sin(theta_g - theta) + vgq_g.cos(theta_g - theta)
        #   theta_g est l'angle du slack, entree globale (nul a l'equilibre,
        #   mais sa DERIVEE ne l'est pas : c'est elle qui manquait au portage).
        rot = s.theta_g - s.theta
        ct, st = sp.cos(rot), sp.sin(rot)
        vgd = s.vd_t * ct - s.vq_t * st
        vgq = s.vd_t * st + s.vq_t * ct

        # --- tension AUX BORNES : portee par le shunt Rg (alg2/alg3) ---
        ved, veq = s.v_ed, s.v_eq
        Et = sp.sqrt(ved**2 + veq**2)

        wr = s.Dwr + s.w_ref
        Ce = s.phi_d * s.i_q - s.phi_q * s.i_d
        Cm = s.Pm / wr
        DP = (wr - s.w_ref) / mp

        dDw1 = (s.Dwr - s.Dw1) / s.T_LP
        dv1 = (s.T_HP * s.K_PSS * dDw1 - s.v1) / s.T_HP
        dv2 = (s.T1n * dv1 + s.v1 - s.v2) / s.T1d
        dvpss = (s.T2n * dv2 + s.v2 - s.vpss) / s.T2d

        de1 = (Et - s.e1) / s.Tr
        # en butee (mode_e2 = +1 ou -1), l'amplificateur est fige a sa limite
        de2 = (1 - s.mode_e2**2) * (
            (s.Ka * (s.v_ref - s.e1 - s.e3 + s.vpss) - s.e2) / s.Ta
        )
        defd = (s.Ke * s.e2 - s.efd) / s.Te
        de3 = (s.Kfd * defd - s.e3) / s.Tfd

        f = [
            # --- transfo de raccordement, repere ROTOR (symSG) ---
            (s.wb / s.Lt) * (ved - vgd - s.Rt * s.i_gd + wr * s.Lt * s.i_gq),
            (s.wb / s.Lt) * (veq - vgq - s.Rt * s.i_gq - wr * s.Lt * s.i_gd),
            # --- stator ---
            s.wb * (ved + s.Ra * s.i_d + wr * s.phi_q),
            s.wb * (veq + s.Ra * s.i_q - wr * s.phi_d),
            # --- rotor ---
            (s.wb * s.Rfd / s.Lad) * s.efd - s.wb * s.Rfd * s.i_fd,
            s.wb * (-s.R1d * s.i_1d),
            s.wb * (-s.R1q * s.i_1q),
            s.wb * (-s.R2q * s.i_2q),
            # --- mecanique, angle, regulateur de vitesse ---
            (Cm - Ce - s.KD * s.Dwr) / (2 * s.H),
            s.wb * s.Dwr,
            (s.p_ref - DP - s.Pm) / s.Tg,
            # --- PSS ---
            dDw1,
            dv1,
            dv2,
            dvpss,
            # --- AVR ---
            de1,
            de2,
            defd,
            de3,
        ]

        g = [
            -s.phi_d + (-Ld * s.i_d + s.Lad * s.i_fd + s.Lad * s.i_1d),
            -s.phi_q + (-Lq * s.i_q + s.Laq * s.i_1q + s.Laq * s.i_2q),
            -s.phi_fd + (-s.Lad * s.i_d + Lffd * s.i_fd + s.Lad * s.i_1d),
            -s.phi_1d + (-s.Lad * s.i_d + s.Lad * s.i_fd + L11d * s.i_1d),
            -s.phi_1q + (-s.Laq * s.i_q + L11q * s.i_1q + s.Laq * s.i_2q),
            -s.phi_2q + (-s.Laq * s.i_q + s.Laq * s.i_1q + L22q * s.i_2q),
            # alg2/alg3 de symSG : la borne machine est portee par le shunt Rg
            s.v_ed - (s.i_d - s.i_gd) * s.Rg,
            s.v_eq - (s.i_q - s.i_gq) * s.Rg,
        ]

        # --- courant injecte au noeud : repere ROTOR -> GLOBAL, base systeme ---
        #   rotation inverse, exp(+j.(theta - theta_g)), avec cos et sin de
        #   ``rot`` deja calcules : cos(theta - theta_g) = ct, sin = -st.
        iD = s.kS * (s.i_gd * ct + s.i_gq * st)
        iQ = s.kS * (-s.i_gd * st + s.i_gq * ct)
        h = [
            iD,
            iQ,
            ved * s.i_gd + veq * s.i_gq,
            -ved * s.i_gq + veq * s.i_gd,
            Et,
            wr,
            ved,
            veq,
            s.i_gd,
            s.i_gq,
            Ce,
            s.theta,
        ]
        return f, g, h

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        p = self.pdict
        Ld, Lq = p["Lad"] + p["Ll"], p["Laq"] + p["Ll"]
        Lffd = p["Lad"] + p["Lfd"]
        Rg, Ra, Rt, Lt = p["Rg"], p["Ra"], p["Rt"], p["Lt"]

        V = ctx.v["t"]  # tension du noeud RESEAU
        S = ctx.s.get("t", 0j)
        Inet = (np.conj(S / V) if abs(V) > 1e-12 else 0j) / p["kS"]

        if self.init_mode == "symSG":
            # V1 : la borne machine est le noeud que network_form.m a cree
            # DERRIERE le transfo, et le courant du load flow y est pris pour
            # celui de l'ARMATURE -- le transfo n'emporte alors que ce qui reste
            # une fois le shunt servi (cf. l'entete de la classe).
            Ve, Im = V, Inet
            Ig = Inet - Ve / Rg
        else:
            # remonter du noeud reseau vers la borne machine, a travers Rt + j.Lt
            Ig = Inet
            Ve = V + complex(Rt, Lt) * Ig
            Im = Ig + Ve / Rg  # courant machine
        E = Ve + complex(Ra, Lq) * Im
        theta = float(np.angle(E)) - np.pi / 2.0

        rot = np.exp(-1j * theta)
        e, ig, i = Ve * rot, Ig * rot, Im * rot
        if self.init_mode == "symSG":
            # Quirk assume de la reference : symSG ramene la TENSION au repere
            # rotor par exp(j.(alpha - theta)) alors que v0 porte deja exp(j.alpha)
            # -- les courants, eux, sont tournes correctement. La tension de borne
            # se retrouve donc en avance d'un alpha de trop, 1,4 degre ici.
            e = Ve * np.exp(1j * (float(np.angle(Ve)) - theta))
        ved, veq = e.real, e.imag
        i_d, i_q = i.real, i.imag

        phi_q = -(ved + Ra * i_d)
        phi_d = veq + Ra * i_q
        i_1d = i_1q = i_2q = 0.0
        i_fd = (phi_d + Ld * i_d) / p["Lad"]
        phi_fd = -p["Lad"] * i_d + Lffd * i_fd
        phi_1d = -p["Lad"] * i_d + p["Lad"] * i_fd
        phi_1q = -p["Laq"] * i_q
        phi_2q = -p["Laq"] * i_q

        efd = p["Lad"] * i_fd
        Ce = phi_d * i_q - phi_q * i_d
        Pm = Ce * ctx.w
        Et = abs(Ve)
        e2 = efd / p["Ke"]
        return {
            "i_gd": ig.real,
            "i_gq": ig.imag,
            "phi_d": phi_d,
            "phi_q": phi_q,
            "phi_fd": phi_fd,
            "phi_1d": phi_1d,
            "phi_1q": phi_1q,
            "phi_2q": phi_2q,
            "Dwr": 0.0,
            "theta": theta,
            "Pm": Pm,
            "Dw1": 0.0,
            "v1": 0.0,
            "v2": 0.0,
            "vpss": 0.0,
            "e1": Et,
            "e2": e2,
            "efd": efd,
            "e3": 0.0,
            "i_d": i_d,
            "i_q": i_q,
            "i_fd": i_fd,
            "i_1d": i_1d,
            "i_1q": i_1q,
            "i_2q": i_2q,
            "v_ed": ved,
            "v_eq": veq,
            "p_ref": Pm,
            "v_ref": Et + e2 / p["Ka"],
        }
