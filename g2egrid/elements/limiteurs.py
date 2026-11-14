"""Limiteurs et butees : les modes DISCRETS des equipements.

POURQUOI -- un convertisseur ne sait pas fournir plus que son courant
nominal, quelques dizaines de pourcent pres : sur un court-circuit, sa
consigne de courant part a plusieurs pu, et c'est sa limitation qui decide de
ce qu'il fait pendant le defaut -- donc de ce qui en sort. Sans elle, une
simulation de defaut d'un reseau riche en convertisseurs ne dit rien de
fiable. De meme, une excitatrice a une tension maximale : l'AVR y bute quand
la tension s'effondre.

COMMENT -- a la maniere de VeraGrid (sa "logique procedurale", evaluee aux
pas acceptes et jamais dans le residu de Newton), une limite n'est pas une
fonction non derivable dans les equations. C'est un MODE : un parametre
drapeau (``mode_lim``, ``mode_e2``) qui vaut 0 en fonctionnement libre et
change la forme des equations quand il vaut autre chose. La simulation
surveille des GARDES -- des fonctions continues de l'etat -- et bascule le
drapeau a l'instant precis ou une garde franchit zero
(:meth:`g2egrid.temporel.ModeleNonLineaire.simuler`). Entre deux basculements,
les equations restent lisses et le Newton ne voit jamais de coin.

A DRAPEAU NUL, LES EQUATIONS SONT EXACTEMENT CELLES D'AVANT : le terme ajoute
y est multiplie par zero, a la derivee pres. Le modele linearise, calcule en
mode libre, ne bouge donc pas d'un bit -- ni la parite avec le jumeau MATLAB.
Un point de fonctionnement AU-DELA d'une limite (un GFM surcharge) est
detecte au depart de la simulation, qui commence alors en mode limite.

Limitation de courant (GFM, GFL) : la consigne ``(isd_ref, isq_ref)`` est
ramenee DANS le disque ``|i| <= Imax`` (en pu de l'equipement), direction
conservee : ``k = min(1, Imax/|i_ref|)``. Les integrateurs de la boucle
EXTERNE suivent la consigne limitee par RETRO-CALCUL ::

    dM/dt = Ki.e + (i_lim - i_ref) / Taw

Un premier essai figeait simplement ces integrateurs (integration
conditionnelle) : sur un defaut franc aux bornes d'un GFM, la consigne
"libre" -- proportionnel faible plus integrateur fige -- redescendait sous la
limite pendant le defaut, le limiteur relachait, et la boucle de tension,
redevenue libre, demandait 25 pu au convertisseur. Avec le retro-calcul,
l'integrateur suit la limite : la consigne libre reste la vraie demande,
au-dessus de ``Imax`` tant que la tension manque, en dessous des qu'elle est
revenue. Et comme ``k`` et le terme de retro-calcul sont CONTINUS a la
frontiere, basculer n'y change rien a la dynamique : le mode ne sert qu'a
garder le Newton hors du coin du ``min``, la ou il n'a rien a faire.

La consigne libre n'est pas recalculee ici : elle se LIT sur le modele, par
l'equation de l'integrateur de la boucle de courant, ``dM_CL/dt =
KiCL.(i_ref - i_s)``, evaluee en mode libre. Rien n'est donc ecrit deux fois.

A L'EQUILIBRE, cette lecture dit si le point de fonctionnement est DANS le
domaine du mode libre. ``dM_CL/dt = 0`` y donne ``i_ref = i_s`` : la
grandeur que sature le limiteur est le courant COTE CONVERTISSEUR (``isd``,
``isq``, a travers ``L1``), en pu de l'equipement (base ``Sn``), et non le
courant cote reseau ``ig``, qui en differe du courant du condensateur du
filtre. :meth:`LimiteurCourant.grandeurs_limitees` rend cette grandeur et sa
limite, par le MEME calcul que la garde ; :mod:`g2egrid.noyau.limites` s'en sert
pour signaler un equilibre au-dela d'une limite, ou le modele linearise --
toujours en mode libre -- ne represente plus le convertisseur.
"""

from __future__ import annotations

import math

import numpy as np
import sympy as sp

__all__ = [
    "DUREE_MIN_MODE",
    "HYSTERESIS",
    "IMAX_BLOQUE",
    "LimiteurCourant",
    "LimiteurCourantTension",
    "MMAX",
    "REGLAGES_BUS_CONTINU",
    "depassement_vi",
    "impedance_virtuelle",
    "kRv_dimensionne",
    "limiter_courant",
    "limiter_tension",
]

#: Ecart relatif entre l'entree en limitation et la sortie. La dynamique etant
#: continue a la frontiere, il ne protege que du cout des basculements : sur
#: un defaut franc, la consigne, que la resonance du filtre LCL fait osciller,
#: franchissait la limite 70 fois en 30 ms sans lui ni duree minimale.
HYSTERESIS = 0.02
#: Duree minimale dans un mode avant d'en ressortir [s] -- la temporisation de
#: retour de tout limiteur reel, et la seconde garde contre le broutement.
DUREE_MIN_MODE = 5e-3
#: Le courant maximal d'un convertisseur BLOQUE (pu de l'equipement) : ses
#: interrupteurs ne conduisent plus, sa consigne de courant est ramenee a
#: (presque) zero par la limitation elle-meme -- memes equations, meme
#: retro-calcul des integrateurs, et la boucle de courant eteint le courant
#: en quelques millisecondes au lieu de le couper en zero seconde.
IMAX_BLOQUE = 1e-3
#: Regularisation du module de la consigne (une consigne nulle reste finie).
_EPS2 = 1e-12


def limiter_courant(s, isd_ref, isq_ref):
    """La consigne limitee ``(isd, isq)`` et les deux termes de retro-calcul
    ``(aw_d, aw_q)`` a ajouter aux derivees des integrateurs de la boucle
    externe. A ``mode_lim = 0`` : consigne inchangee, termes nuls."""
    module = sp.sqrt(isd_ref**2 + isq_ref**2 + _EPS2)
    k = 1 - s.mode_lim + s.mode_lim * sp.Min(1, s.Imax / module)
    isd_lim, isq_lim = isd_ref * k, isq_ref * k
    aw_d = s.mode_lim * (isd_lim - isd_ref) / s.Taw
    aw_q = s.mode_lim * (isq_lim - isq_ref) / s.Taw
    return isd_lim, isq_lim, aw_d, aw_q


# ---------------------------------------------------------------------------
# Limitation par IMPEDANCE VIRTUELLE (GFM a boucle de tension)
# ---------------------------------------------------------------------------
def depassement_vi(i, I_vi, dI_vi):
    """Le depassement LISSE ``ov(i)`` du seuil ``I_vi``, de classe C1 ::

        ov = 0                       i <= I_vi
        ov = (i - I_vi)^2 / (2.dI)   I_vi < i < I_vi + dI
        ov = i - I_vi - dI/2         i >= I_vi + dI

    Nul EXACTEMENT sous le seuil (et sa derivee aussi) : sous le seuil, les
    equations et le jacobien sont ceux du convertisseur sans impedance
    virtuelle, au bit pres. Une fonction ``max(0, i - I_vi)`` aurait un coin
    en ``I_vi`` ; une fonction partout positive (softplus) deplacerait le
    spectre sous le seuil. La parabole de raccord, de largeur ``dI``, donne
    un jacobien continu : le Newton et le petit signal n'ont pas de coin a
    traverser. Fonctionne sur des nombres comme sur des expressions SymPy.
    """
    d = i - I_vi
    return sp.Piecewise(
        (0, d <= 0),
        (d**2 / (2 * dI_vi), d < dI_vi),
        (d - dI_vi / 2, True),
    )


def impedance_virtuelle(s, isd, isq):
    """La chute de tension ``(dv_d, dv_q) = Zv.i`` de l'impedance virtuelle,
    a RETRANCHER de la consigne de la boucle de tension.

    ``Zv = Rv + j.Xv`` avec ``Rv = kRv.ov(|i|)`` et ``Xv = (X/R)_vi.Rv`` : une
    impedance virtuelle qui GRANDIT avec le depassement du seuil ``I_vi``
    (:func:`depassement_vi`) et qui est nulle en dessous. ``i`` est le
    courant cote convertisseur ``(isd, isq)``, celui que lit le limiteur
    circulaire (donc le meme ``Imax``). Un rapport ``X/R`` eleve garde la
    chute en quadrature avec le courant : elle abaisse la tension sans
    retarder l'angle, ce qui preserve la synchronisation.

    A ``kRv = 0`` (defaut) : chute nulle, equations inchangees au bit pres.
    """
    module = sp.sqrt(isd**2 + isq**2 + _EPS2)
    Rv = s.kRv * depassement_vi(module, s.I_vi, s.dI_vi)
    Xv = s.XsurR_vi * Rv
    return Rv * isd - Xv * isq, Rv * isq + Xv * isd


def kRv_dimensionne(Imax, I_vi=1.0, XsurR_vi=0.1, V=1.0, Zdefaut=0j, dI_vi=0.05):
    """Le gain ``kRv`` pour que le courant PERMANENT d'un defaut franc soit
    ``Imax`` (pu de l'equipement), a ``V`` pu de consigne de tension.

    En regime permanent du defaut, la boucle de tension tient ``ve + Zv.i =
    V`` et le reseau impose ``ve = Zdefaut.i`` : le courant ``I`` verifie
    ``|Zv(I) + Zdefaut| = V/I``, ``Zv(I) = kRv.ov(I).(1 + j.X/R)``.
    ``Zdefaut`` est l'impedance COMPLEXE (pu de l'equipement) entre la
    tension controlee et le defaut : filtre ``L2``, transformateur, defaut.
    Les deux impedances s'ajoutent en phaseurs et non en modules : une VI
    purement resistive devant un reseau inductif limite MOINS qu'une somme de
    modules ne le ferait (``sqrt(R^2 + X^2) < R + X``), et c'est pourquoi la
    VI doit etre dimensionnee avec l'impedance de defaut reelle. On resout
    pour ``I = Imax`` le module ``a = |Zv|`` par l'equation du second degre
    ``a^2 + 2 a Re(e^{j phi} conj(Zdefaut)) + |Zdefaut|^2 = (V/Imax)^2``.

    Ce calcul est STATIQUE : la dynamique de la boucle de tension
    (``tau_VL``) n'y entre pas. Le pic transitoire peut depasser ``Imax``
    avant que la VI ne s'etablisse -- c'est le role du limiteur circulaire,
    gardien de dernier recours.
    """
    ov = (
        Imax - I_vi - dI_vi / 2
        if Imax >= I_vi + dI_vi
        else (Imax - I_vi) ** 2 / (2 * dI_vi)
    )
    if Imax <= I_vi or ov <= 0:
        raise ValueError("Imax doit depasser I_vi pour qu'une VI soit utile")
    e = complex(1.0, XsurR_vi)
    e /= abs(e)
    zd = complex(Zdefaut)
    b = (e * zd.conjugate()).real
    c = abs(zd) ** 2 - (V / Imax) ** 2
    if c >= 0:
        return 0.0  # l'impedance de defaut seule limite deja le courant
    a = -b + math.sqrt(b * b - c)
    return a / (math.hypot(1.0, XsurR_vi) * ov)


def _valeur(el, nom: str, x, z, p) -> float:
    """La valeur de la variable ``nom`` de ``el``, ou qu'elle vive : etat,
    variable algebrique (dynamique rendue STATIQUE par
    :meth:`~g2egrid.noyau.element.Element.reduire`) ou parametre (dynamique FIGEE)."""
    if nom in el.states:
        return float(x[el.states.index(nom)])
    if nom in el.algebs:
        return float(z[el.algebs.index(nom)])
    return float(p[el.params.index(nom)])


class LimiteurCourant:
    """Les gardes de la limitation de courant (a melanger a un element).

    L'element doit avoir les variables ``isd``, ``isq`` et les parametres
    ``KiCL``, ``Imax``, ``mode_lim`` ; s'il a les etats ``M_CLd``, ``M_CLq``
    (integrateurs de la boucle de courant), c'est sur eux que se lit la
    consigne libre. C'est le cas des GFM, des GFL et des MMC de ce paquet.
    """

    modes = ("mode_lim",)
    #: un convertisseur BLOQUE (manoeuvre ``"blocage"``) reste en mode limite,
    #: a ``Imax`` quasi nul : sa garde n'est plus surveillee
    bloque: bool = False

    def _courant_libre(self, cm, x, z, u, p) -> float:
        """Le module de la consigne de courant LIBRE, ``|i_ref|`` avant le
        limiteur, en pu de l'equipement.

        Elle se lit sur l'integrateur de la boucle de courant :
        ``i_ref = (dM_CL/dt)/KiCL + i_s``, ``dM_CL/dt`` evalue en mode libre.
        Dans une variante ou cet integrateur n'est plus un etat
        (:meth:`~g2egrid.noyau.element.Element.reduire`, groupe ``boucle_i`` de
        l'atelier), on prend ``i_s`` : c'est exact quand l'integrateur est
        STATIQUE (sa contrainte ``0 = KiCL.(i_ref - i_s)`` tient a chaque
        instant) et a l'equilibre quand il est FIGE a sa valeur d'equilibre ;
        ailleurs, ``i_s`` suit la consigne a la vitesse de la boucle de
        courant. Avant, ces variantes faisaient echouer la garde.
        """
        params, etats = self.params, self.states
        i_d = _valeur(self, "isd", x, z, p)
        i_q = _valeur(self, "isq", x, z, p)
        if "M_CLd" in etats or "M_CLq" in etats:
            libre = np.array(p, dtype=float)
            libre[params.index("mode_lim")] = 0.0
            # la saturation de tension ajoute SON retro-calcul a ce meme
            # integrateur ; le laisser fausserait la consigne lue de awv/KiCL
            # (une fraction notable de pu des que Taw_v est court) : on
            # l'annule aussi
            if "mode_vsat" in params:
                libre[params.index("mode_vsat")] = 0.0
            f = cm.f(x, z, u, libre)
            KiCL = p[params.index("KiCL")]
            if "M_CLd" in etats:
                i_d = f[etats.index("M_CLd")] / KiCL + i_d
            if "M_CLq" in etats:
                i_q = f[etats.index("M_CLq")] / KiCL + i_q
        return math.hypot(i_d, i_q)

    def grandeurs_limitees(self, cm, x, z, u, p) -> list[tuple[str, float, float]]:
        """``(drapeau, grandeur libre, limite)`` de chaque limiteur REGLE.

        La grandeur est celle que lit la garde, par le meme calcul : le module
        de la consigne de courant libre face a ``Imax``, tous deux en pu de
        l'equipement. A l'equilibre, c'est ``|i_s|``, le courant cote
        convertisseur (voir l'en-tete du module). Un convertisseur BLOQUE n'a
        pas de limiteur surveille : rien.
        """
        if self.bloque:
            return []
        Imax = float(p[self.params.index("Imax")])
        return [("mode_lim", self._courant_libre(cm, x, z, u, p), Imax)]

    def gardes(self, cm, x, z, u, p) -> list[tuple[str, float, float, int]]:
        if self.bloque:
            return []
        params = self.params
        j = params.index("mode_lim")
        courant = self._courant_libre(cm, x, z, u, p)
        Imax = p[params.index("Imax")]
        if p[j] == 0.0:
            return [("mode_lim", courant - Imax, 1.0, +1)]
        return [("mode_lim", courant - Imax * (1 - HYSTERESIS), 0.0, -1)]


# ---------------------------------------------------------------------------
# Saturation de tension (modulation) des GFM et des GFL
# ---------------------------------------------------------------------------
#: Indice de modulation maximal d'une modulation vectorielle : la tension que
#: l'onduleur sait fabriquer, abs(e) <= MMAX.v_dc. Sur bus continu ideal, les
#: GFM et les GFL ont Mmax = 0 par defaut (saturation de tension inactive) ;
#: MMAX est la valeur a regler.
MMAX = 1.15
#: Les reglages par defaut d'un GFM dont la tension continue est une VARIABLE
#: (GFM3/GFM4, matching, port continu) : la saturation de tension y est
#: ACTIVE. Ce n'est pas un raffinement. Un GFM est une source de tension : il
#: fournit a son reseau alternatif la puissance que celui-ci appelle, et la
#: tire du bus continu quelle qu'en soit la tension. Sans borne sur la
#: modulation ``e/v_dc``, c'est une charge a puissance constante sur un
#: condensateur, qui le vide en un temps FINI : la simulation n'a plus de
#: solution. Cas type : un defaut alternatif aux bornes de la station qui
#: tient la tension d'une liaison continue, pendant que la station GFM de
#: l'autre bout soutient son reseau en puisant dans le continu
#: (``cigre_mv_dc``, defaut au bus 15 : le Newton lachait a 83 ms, v_dc a
#: 0,33 pu). Bornee, la modulation ne sait plus fabriquer la tension
#: alternative quand v_dc tombe vers abs(e)/Mmax : le transit s'arrete de
#: lui-meme et v_dc se tient (0,39 pu).
#:
#: Pas les GFL2 : un GFL tire du continu ce que ses consignes de courant
#: demandent, et la station qui tient v_dc le reduit d'elle-meme quand v_dc
#: baisse. Saturee, sa modulation lache au contraire la commande du courant
#: des que v_dc passe sous abs(e)/Mmax -- ce qui arrive sur tout defaut continu
#: ou blocage de la station d'en face (``cigre_b4_dcs1``) : le Newton lache
#: sur un defaut franc pole a pole, la liaison oscille apres un blocage.
REGLAGES_BUS_CONTINU = dict(Mmax=MMAX)


def limiter_tension(s, md, mq):
    """La modulation saturee ``(md, mq)``, ramenee dans le disque
    ``|m| <= Mmax`` direction conservee, et les deux termes de retro-calcul
    ``(aw_d, aw_q)`` (en indice de modulation par seconde) a ajouter, une
    fois multiplies par ``v_dc``, aux integrateurs de la boucle de COURANT.
    A ``mode_vsat = 0`` : modulation inchangee, termes nuls."""
    module = sp.sqrt(md**2 + mq**2 + _EPS2)
    k = 1 - s.mode_vsat + s.mode_vsat * sp.Min(1, s.Mmax / module)
    md_lim, mq_lim = md * k, mq * k
    # Un retro-calcul RAPIDE, a l'echelle de la boucle de courant (Taw_v =
    # tau_CL/2 par defaut, cf. GFM_DEFAULT et GFL_DEFAULT) : il garde la
    # demande libre de la boucle de courant collee a ce que l'onduleur sait
    # faire, et la saturation relache des que le courant est tenu. Sur le GFM
    # (tau_CL de 100 us), lent (dizaines de ms), il laisse l'integrateur de
    # courant s'emballer et la saturation ne relache plus ; entre les deux
    # (ms), la modulation libre colle au coin du Min et le Newton ne converge
    # plus.
    return (
        md_lim,
        mq_lim,
        s.mode_vsat * (md_lim - md) / s.Taw_v,
        s.mode_vsat * (mq_lim - mq) / s.Taw_v,
    )


class LimiteurCourantTension(LimiteurCourant):
    """Limitation de courant ET saturation de tension (GFM, GFL).

    La saturation de tension porte sur la modulation : l'onduleur ne sait pas
    fabriquer plus de ``Mmax.v_dc``. Sa garde lit la modulation LIBRE sur le
    modele, par l'equation de l'inductance cote onduleur ::

        e_d = (L1/wb).d isd/dt + ved + R1.isd - w.L1.isq

    evaluee a ``mode_vsat = 0`` (``w`` se lit sur ``d theta/dt``), comme la
    garde de courant lit la consigne libre sur la boucle de courant.

    Un GFL porte sa modulation en variables ALGEBRIQUES ``md``, ``mq`` -- la
    modulation libre, la saturation s'appliquant apres : la garde les lit
    directement.
    """

    modes = ("mode_lim", "mode_vsat")

    def _vdc_garde(self, x, z, u, p) -> float:
        if "vdc" in (*self.states, *self.algebs, *self.params):
            return _valeur(self, "vdc", x, z, p)
        if self.ports_dc:
            v = float(u[self.port_dc_input_idx(0)])
            b = p[self.params.index("bloque")] if "bloque" in self.params else 0.0
            return v * (1 - b) + b
        return 1.0

    def _vsat_reglee(self, p) -> bool:
        """La saturation de tension est-elle reglee (``Mmax > 0``) ?"""
        return "mode_vsat" in self.params and p[self.params.index("Mmax")] > 0.0

    def _modulation_libre(self, cm, x, z, u, p) -> float:
        """Le module de la modulation LIBRE, ``|e|/v_dc`` avant la saturation.

        Sur un GFL, ce sont les algebriques ``md``, ``mq``. Sur un GFM, elle se
        lit sur l'equation de l'inductance cote onduleur, ``d isd/dt`` et
        ``d theta/dt`` evalues a ``mode_vsat = 0`` ; un ``isd`` ou un ``theta``
        qui n'est plus un etat (variante reduite) y compte pour une derivee
        nulle -- exact a l'equilibre, et toujours pour une dynamique statique.
        """
        if "md" in self.algebs:  # GFL : la modulation libre est algebrique
            return math.hypot(z[self.algebs.index("md")], z[self.algebs.index("mq")])
        params, etats = self.params, self.states
        libre = np.array(p, dtype=float)
        libre[params.index("mode_vsat")] = 0.0
        f = cm.f(x, z, u, libre)
        P = lambda n: p[params.index(n)]  # noqa: E731
        X = lambda n: _valeur(self, n, x, z, p)  # noqa: E731
        D = lambda n: f[etats.index(n)] if n in etats else 0.0
        L1, R1, wb = P("L1_pu"), P("R1_pu"), P("wb")
        w = D("theta") / wb + u[self.wref_idx()]
        e_d = D("isd") * L1 / wb + X("ved") + R1 * X("isd") - w * L1 * X("isq")
        e_q = D("isq") * L1 / wb + X("veq") + R1 * X("isq") + w * L1 * X("isd")
        return math.hypot(e_d, e_q) / max(self._vdc_garde(x, z, u, p), 1e-6)

    def grandeurs_limitees(self, cm, x, z, u, p) -> list[tuple[str, float, float]]:
        """Celles de :class:`LimiteurCourant`, et la modulation libre face a
        ``Mmax`` quand la saturation de tension est reglee."""
        out = super().grandeurs_limitees(cm, x, z, u, p)
        if self.bloque or not self._vsat_reglee(p):
            return out
        Mmax = float(p[self.params.index("Mmax")])
        return out + [("mode_vsat", self._modulation_libre(cm, x, z, u, p), Mmax)]

    def gardes(self, cm, x, z, u, p) -> list[tuple[str, float, float, int]]:
        out = super().gardes(cm, x, z, u, p)
        if self.bloque or not self._vsat_reglee(p):
            return out  # saturation de tension non reglee
        j = self.params.index("mode_vsat")
        m = self._modulation_libre(cm, x, z, u, p)
        return out + [self._garde_vsat(m, p, j)]

    def _garde_vsat(self, m: float, p, j: int) -> tuple[str, float, float, int]:
        Mmax = p[self.params.index("Mmax")]
        if p[j] == 0.0:
            return ("mode_vsat", m - Mmax, 1.0, +1)
        return ("mode_vsat", m - Mmax * (1 - HYSTERESIS), 0.0, -1)
