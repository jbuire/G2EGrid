"""Reseau : load flow, equilibre, assemblage du modele d'etat global.

CHAINE DE CALCUL
    1. add_bus / add          declaration du reseau
    2. run_loadflow()         Newton-Raphson -> tensions complexes des bus
    3. solve_equilibrium()    par element : initialisation analytique puis
                              Newton sur [f; g] = 0 (entrees figees)
    4. build_state_space()    assemblage global -> (A, B, C, D) explicites

ASSEMBLAGE
Chaque element expose un modele DAE local (f, g, h) et des ports. Le cablage
est un simple routage de signaux ::

    entree tension de port  <- etat de tension du Bus
    entree isum du Bus      <- somme des courants de port (loi des noeuds)

On resout la boucle par LFT ::

    U = Ruext.uext + Ry.Y  et  Y = Hx.x + Hz.z + Hu.U
    => (I - Ry.Hu) U = Ruext.uext + Ry.Hx.x + Ry.Hz.z

Comme aucun element n'a de transmission directe tension -> courant, la
matrice (I - Ry.Hu) est bien conditionnee (pas de boucle algebrique).
Les variables algebriques restantes (relations flux/courant de la machine)
sont ensuite eliminees par complement de Schur.
"""

from __future__ import annotations

from collections import Counter

from dataclasses import dataclass, field

import numpy as np
import scipy.linalg
from scipy import sparse
from scipy.linalg import lu_factor, lu_solve
from scipy.linalg.lapack import dgecon
from scipy.optimize import root

from g2egrid.noyau import limites
from g2egrid.noyau.element import W_REF, Element, InitContext
from g2egrid.elements.branch import LineDC
from g2egrid.elements.bus import (
    NOEUDS_DEFAUT,
    R_CAUSALITE,
    Bus,
    BusDC,
    IdealBus,
    ResistiveBus,
)
from g2egrid.elements.source import SourceDC
from g2egrid.noyau.loadflow import (
    DC_P,
    DC_SLACK,
    PQ,
    SLACK,
    LoadFlowACDCResult,
    LoadFlowResult,
    ProblemeACDC,
    VSCSpec,
    run_loadflow,
)
from g2egrid.noyau.loadflow import build_ybus as ybus_creux

__all__ = ["Network", "StateSpace", "DAELineaire", "BusSpec", "BusSpecDC"]


@dataclass
class BusSpecDC:
    """Declaration d'un noeud CONTINU.

    Deliberement a part de :class:`BusSpec`, et range dans un autre
    dictionnaire : ``_bus_index``, ``build_ybus`` et ``run_loadflow`` ne
    voient litteralement rien de nouveau. Ce n'est pas de la prudence, c'est
    une garantie structurelle -- le load flow alternatif ne peut pas tomber
    sur un noeud continu, puisqu'il ne les lit pas.

    ``btype`` decide de ce que ``V`` veut dire. ``DC_P`` (2, le defaut) : la
    tension est une INCONNUE du load flow mixte, et ``V`` n'en est que
    l'amorce -- c'est le comportement du lot 3, et c'est pourquoi c'est le
    defaut : un reseau bati a la main, sans type declare, retrouve son point de
    fonctionnement d'avant. ``DC_SLACK`` (1) : la tension est IMPOSEE, et le
    noeud doit etre tenu, soit par un convertisseur en mode ``Vdc...`` -- dont
    la puissance alternative devient alors l'inconnue -- soit par une source
    continue, dont la f.e.m. se deduit du load flow comme celle d'un bus infini.

    ``Pinj`` est l'injection NETTE du noeud (pu base systeme), nommee ainsi
    pour ne pas la confondre avec le couple generation/charge du cote
    alternatif : sur un noeud continu, rien ne se compense, il n'y a qu'un
    solde.
    """

    id: int
    Un: float
    C: float = 0.0  # capacite pu (base systeme)
    V: float = 1.0
    btype: int = DC_P
    Pinj: float = 0.0  # injection nette pu (base systeme)


@dataclass
class BusSpec:
    id: int
    Un: float
    btype: int = PQ
    V: float = 1.0
    theta: float = 0.0
    Pg: float = 0.0  # pu base systeme
    Qg: float = 0.0
    Pl: float = 0.0
    Ql: float = 0.0
    C: float = 0.0  # susceptance shunt pu (capacite de bus / charge ligne)
    R: float = R_CAUSALITE  # resistance de causalite pu, utilisee si C == 0


@dataclass
class StateSpace:
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    states: list[str]
    inputs: list[str]
    outputs: list[str]
    #: etats ELIMINES par la reduction d'indice (courants en serie a un noeud
    #: ideal) : ils ne sont plus des etats, mais se reconstruisent par ``T``
    elimines: list[str] = field(default_factory=list)
    #: ``x_complet = T . x`` ; ``None`` quand rien n'a ete elimine
    T: np.ndarray | None = None
    #: les etats AVANT la reduction, dans l'ordre des lignes de ``T`` (vide
    #: sans reduction) : la base complete de :func:`g2egrid.modal.spectre_complet`
    etats_complets: list[str] = field(default_factory=list)
    #: inverse a gauche de ``T`` (``T_gauche . T = I``) dont le noyau porte
    #: les directions rapides des noeuds ideaux : ``w = w_r . T_gauche``
    T_gauche: np.ndarray | None = None

    def eig(self) -> np.ndarray:
        return np.linalg.eigvals(self.A)

    def marge_stabilite(self, tol: float = 1e-6, facteur: float = 10.0) -> float:
        """Pire partie reelle, apres retrait de ce que ``eig`` ne sait pas dire.

        Rend ``max_i (Re(lambda_i) - seuil_i)`` : negatif si tous les modes sont
        stables, positif sinon. Le seuil est PROPRE A CHAQUE VALEUR PROPRE ::

            seuil_i = max(tol, facteur . eps . ||A|| . kappa_i)

        ``kappa_i = ||x_i|| ||y_i|| / |y_i^H x_i|`` est le conditionnement de
        ``lambda_i`` (vecteurs propres a droite ``x`` et a gauche ``y``) : un
        algorithme QR remontant l'erreur inverse de ``A``, il deplace
        ``lambda_i`` de ``eps . ||A|| . kappa_i`` au plus. Un seuil unique du
        type ``rel . max|A|`` etait faux dans les deux sens : avec un bus de
        C = 1e-9 pu, ``max|A|`` atteint 3e11 et le seuil montait a plus de
        100 s^-1, masquant un mode instable de +61 s^-1 (ILoad) alors que son
        conditionnement (3,5e3) ne justifie que 0,4 s^-1.

        Le mode de jauge d'un reseau sans bus infini (faire tourner en bloc
        tous les angles ne change aucune equation : lambda = 0) a un
        conditionnement voisin de 1 et ne sort qu'a quelques ``eps . ||A||`` :
        c'est le seul que ce seuil absorbe. Une valeur propre (quasi) defective
        (kappa >= 1e8) n'est pas absorbee : on ne s'en excuse pas.

        Le calcul des vecteurs propres (3 fois le cout de ``eig``) n'est fait
        que s'il existe une valeur propre a partie reelle >= ``tol``.
        """
        w = self.eig()
        if not np.any(w.real >= tol):
            return float(np.max(w.real) - tol)
        A = self.A
        w, vl, vr = scipy.linalg.eig(A, left=True, right=True)
        nl = np.linalg.norm(vl, axis=0)
        nr = np.linalg.norm(vr, axis=0)
        croise = np.abs(np.einsum("ij,ij->j", vl.conj(), vr))
        with np.errstate(divide="ignore", invalid="ignore"):
            kappa = nl * nr / croise
        # borne superieure de ||A||_2 : sqrt(||A||_1 ||A||_inf)
        norme = float(np.sqrt(np.linalg.norm(A, 1) * np.linalg.norm(A, np.inf)))
        # au-dela de 1e8 la valeur propre est (quasi) defective : l'erreur ne
        # suit plus eps.kappa (elle va en racine de eps), on n'excuse rien
        erreur = np.where(
            kappa < 1e8, facteur * np.finfo(float).eps * norme * kappa, 0.0
        )
        return float(np.max(w.real - np.maximum(tol, erreur)))

    def is_stable(self, tol: float = 1e-6, facteur: float = 10.0) -> bool:
        """Toutes les valeurs propres a partie reelle negative, a leur erreur pres.

        Voir :meth:`marge_stabilite` : le seuil est celui de chaque valeur
        propre (``max(tol, facteur . eps . ||A|| . kappa_i)``), pas un seuil
        global sur ``max|A|``. Un reseau sans bus infini a un mode de jauge
        exactement nul dont l'arrondi donne quelques ``eps . ||A||`` : c'est le
        seul que le seuil absorbe.
        """
        return self.marge_stabilite(tol, facteur) < 0.0

    def iu(self, name: str) -> int:
        return self.inputs.index(name)

    def iy(self, name: str) -> int:
        return self.outputs.index(name)

    def to_control(self):
        """Objet ``control.StateSpace`` nomme (si python-control est installe)."""
        import control  # import optionnel

        sys = control.ss(self.A, self.B, self.C, self.D)
        sys.set_inputs(self.inputs)
        sys.set_outputs(self.outputs)
        sys.set_states(self.states)
        return sys

    def __repr__(self) -> str:
        mx = np.max(self.eig().real)
        return (
            f"<StateSpace {len(self.states)} etats, {len(self.inputs)} entrees, "
            f"{len(self.outputs)} sorties | max Re(lambda)={mx:+.3e}>"
        )


@dataclass
class Cablage:
    """Le cablage du reseau : qui envoie quoi a qui, et a quel indice.

    C'est la partie de l'assemblage qui ne depend QUE de la structure du
    reseau, pas de son point de fonctionnement : les indices globaux de
    chaque element, et les deux matrices de la fermeture LFT ::

        U = Ruext . uext + Ry . Y

    ``Ruext`` route les entrees externes (consignes, ``w_ref``), ``Ry``
    referme la boucle : tensions de port <- etats de tension des bus,
    courants de port -> loi des noeuds, angle et vitesse du repere <-
    element de reference.

    Les deux voies de resolution partent de la :

    * :meth:`Network.build_state_space` evalue les jacobiens au point de
      fonctionnement, ferme la boucle et elimine les variables algebriques
      par complement de Schur : c'est le modele d'etat lineaire ;
    * :mod:`g2egrid.temporel` garde les equations non lineaires et resout
      ``g(x, z, u) = 0`` et la fermeture LFT **ensemble**, par Newton, a
      chaque evaluation demandee par l'integrateur.
    """

    els: list
    cms: list
    ox: list[int]
    oz: list[int]
    ou: list[int]
    oy: list[int]
    nx: int
    nz: int
    nu: int
    ny: int
    Ruext: np.ndarray
    Ry: np.ndarray
    uext_names: list[str]
    state_names: list[str]
    #: indices, dans Y, des sorties exposees (hors courants de port)
    rows_sorties: np.ndarray
    output_names: list[str]

    def uext0(self) -> np.ndarray:
        """Les entrees externes au point de fonctionnement : chaque consigne
        vaut ce que l'element a retenu dans son ``u0``, et ``w_ref`` vaut 1."""
        v = np.zeros(len(self.uext_names))
        idx = {nm: i for i, nm in enumerate(self.uext_names)}
        v[idx[W_REF]] = 1.0
        for el in self.els:
            for j, s in enumerate(el.svc_inputs):
                nm = f"{el.name}.{s}"
                if nm in idx:
                    v[idx[nm]] = el.u0[j]
        return v

    def x0(self) -> np.ndarray:
        parts = [el.x0 for el in self.els if len(el.x0)]
        return np.concatenate(parts) if parts else np.zeros(0)

    def z0(self) -> np.ndarray:
        parts = [el.z0 for el in self.els if len(el.z0)]
        return np.concatenate(parts) if parts else np.zeros(0)

    def u0(self) -> np.ndarray:
        return np.concatenate([el.u0 for el in self.els])

    def iu(self, name: str) -> int:
        """L'indice de l'entree externe ``name`` dans ``uext``."""
        return self.uext_names.index(name)


class Network:
    """Conteneur du reseau et moteur d'assemblage."""

    def __init__(
        self,
        f0: float = 50.0,
        Sb: float = 1e6,
        noeuds_sans_capacite: str = NOEUDS_DEFAUT,
        r_causalite: float = R_CAUSALITE,
    ):
        """``noeuds_sans_capacite`` : ce que devient un bus declare avec ``C = 0``.

        ``"ideal"`` (defaut, FR-005) : un :class:`IdealBus`, loi des noeuds
        exacte sans element fictif ; les courants en serie y sont reduits
        algebriquement (voir :meth:`DAELineaire.etat`). En temporel, BDF2
        seulement (indice 2), que :meth:`ModeleNonLineaire.simuler` choisit
        de lui-meme.
        ``"resistance"`` (comme ``Node.m``) : un :class:`ResistiveBus`,
        tension ``R.isum`` a travers la resistance de causalite ``r_causalite``
        (pu systeme, les 100 pu de ``Node.m`` par defaut). Un bus qui donne
        son propre ``R`` le garde. C'est le mode des reseaux compares a MATLAB.
        """
        if noeuds_sans_capacite not in ("resistance", "ideal"):
            raise ValueError(
                f"noeuds_sans_capacite = {noeuds_sans_capacite!r} : "
                "'resistance' ou 'ideal'"
            )
        self.noeuds_sans_capacite = noeuds_sans_capacite
        if not r_causalite > 0:
            raise ValueError(f"r_causalite = {r_causalite!r} : une resistance > 0")
        self.r_causalite = float(r_causalite)
        self.f0 = float(f0)
        self.Sb = float(Sb)
        self.wb = 2.0 * np.pi * self.f0
        self.bus_specs: dict[int, BusSpec] = {}
        self.elements: list[Element] = []
        self.buses: dict[int, Bus] = {}
        #: noeuds CONTINUS, tenus a part des alternatifs (cf. BusSpecDC)
        self.bus_specs_dc: dict[int, BusSpecDC] = {}
        self.buses_dc: dict[int, BusDC] = {}
        #: tensions continues d'equilibre, remplies par ``_run_loadflow_acdc``
        self.v_dc: dict[int, float] = {}
        #: mode de conduite de chaque convertisseur, par numero de BUS
        #: ALTERNATIF -- la cle de ``T_VSCs``. Absent = le mode par defaut de
        #: la classe, celui qui ne tient pas la tension continue.
        self.vsc_modes: dict[int, str] = {}
        #: statisme tension continue / puissance des convertisseurs en mode
        #: ``droop``, par bus ALTERNATIF : pu de tension continue par pu de
        #: puissance, sur la base Sn du convertisseur
        self.vsc_statisme: dict[int, float] = {}
        #: statisme Q-V des convertisseurs (bus alternatif -> ``nq``, pu de
        #: tension par pu de reactif sur leur Sn) ; sur un bus PQ seulement
        self.vsc_qv: dict[int, float] = {}
        self.lf: LoadFlowResult | None = None
        self.w_sys = 1.0
        #: pulsation du repere dq global vue par les elements passifs (noeuds,
        #: branches, charges) : "nominale" (1 pu) ou "reference" (vitesse de
        #: l'element de reference, comme w_slack dans G2Elin V2)
        self.repere = "reference"
        #: NOM de l'element qui porte l'angle du repere, quand ce n'est pas
        #: celui du bus slack -- apres le declenchement de la machine de
        #: reference, c'est un :class:`~g2egrid.elements.Repere`. ``None`` = la
        #: regle de V1 (premier element du bus slack a exposer ``theta_out``).
        self.reference: str | None = None
        #: NOM de la machine dont un ``Repere`` suit la vitesse, ou dict
        #: ``{nom: poids}`` d'une moyenne ponderee (centre d'inertie, voir
        #: :meth:`reference_flottante`). ``None`` : il tourne a la vitesse
        #: nominale, s'il ne reste aucune machine.
        self.repere_suit: str | dict[str, float] | None = None
        #: les limiteurs dont la grandeur libre depasse la limite au point
        #: d'equilibre (:class:`g2egrid.noyau.limites.Depassement`), remis a jour par
        #: chaque :meth:`solve_equilibrium` ; vide quand tout est dans les
        #: limites. Le modele linearise, en mode libre, ne represente pas un
        #: convertisseur de cette liste.
        self.depassements: list[limites.Depassement] = []

    # ------------------------------------------------------------------
    # Declaration
    # ------------------------------------------------------------------
    def add(self, element: Element) -> Element:
        """Ajoute un element au reseau.

        Les bus passent par ici comme tout le reste -- ``net.add(Bus(...))`` --
        de sorte que la construction d'un reseau se lise comme les tables de
        ``buildFromMatrices`` cote MATLAB : d'abord les bus (``Y_DER``), puis
        les lignes, les transfos et les sources.
        """
        if isinstance(element, BusDC):
            return self._declarer_bus_dc(element)
        if isinstance(element, Bus):
            return self._declarer_bus(element)
        self._verifier_genre_ports(element)
        self.elements.append(element)
        return element

    def _verifier_genre_ports(self, el: Element) -> None:
        """Chaque port est-il branche sur un bus de SON genre ?

        Un port alternatif pose deux signaux sur son noeud, un port continu un
        seul. Les croiser ne produit pas une erreur la ou on les croise mais
        trois appels plus loin, dans la loi des noeuds, sous la forme d'un
        coefficient ecrit sur l'entree du VOISIN dans ``Ry`` -- une matrice
        fausse, pas une exception. D'ou cette verification, faite a la
        declaration puis refaite au cablage, qui nomme l'element et le bus.
        """
        n_ac = len(el.ports)
        for k, pt in enumerate(el.ports_tous()):
            b = el.buses[k]
            continu = k >= n_ac
            if b in (self.bus_specs_dc if continu else self.bus_specs):
                continue
            if b in (self.bus_specs if continu else self.bus_specs_dc):
                genre, autre = (
                    ("continu", "alternatif") if continu else ("alternatif", "continu")
                )
                raise ValueError(
                    f"{el.name}: le port {genre} '{pt}' pointe sur le bus {b}, "
                    f"qui est {autre}"
                )
            raise ValueError(f"{el.name}: bus {b} non declare")

    def _declarer_bus_dc(self, decl: BusDC) -> BusDC:
        """Enregistre un noeud continu et acheve sa construction.

        Les identifiants de bus sont partages entre les deux reseaux : le
        cablage repere un noeud par son seul numero (``bus_pos``), donc un
        numero reutilise d'un genre a l'autre ferait disparaitre un des deux
        noeuds sans rien dire.
        """
        if decl.bus in self.bus_specs:
            raise ValueError(
                f"{decl.name}: le bus {decl.bus} est deja declare comme alternatif"
            )
        si = decl.spec_si
        if si["C"] <= 0.0:
            raise ValueError(
                f"{decl.name}: C = {si['C']} ; un noeud continu sans capacite "
                "donnerait une equation algebrique sans variable pour la satisfaire"
            )
        spec = BusSpecDC(id=decl.bus, Un=decl.Un, C=si["C"], V=si["V"])
        self.bus_specs_dc[spec.id] = spec
        decl.finaliser(self.wb)
        self.buses_dc[spec.id] = decl
        return decl

    def _declarer_bus(self, decl: Bus) -> Bus:
        """Enregistre un bus declare et acheve sa construction.

        C'est ici que se fait la conversion SI -> pu : ``Bus`` porte des watts
        et des vars, comme ``Y_DER``, et seul le reseau connait la base ``Sb``.
        Un bus sans capacite shunt devient un noeud algebrique a resistance de
        causalite (``ResistiveBus``), comme ``setResistor(100)`` dans
        ``Network.m`` -- de :attr:`r_causalite` pu, sauf ``R`` propre au bus.
        """
        if decl.bus in self.bus_specs_dc:
            raise ValueError(
                f"{decl.name}: le bus {decl.bus} est deja declare comme continu"
            )
        si = dict(decl.spec_si)
        if si.get("R") is None:
            si["R"] = self.r_causalite
        spec = BusSpec(
            id=decl.bus,
            Un=decl.Un,
            btype=si["btype"],
            V=si["V"],
            theta=si["theta"],
            Pg=si["Pg"] / self.Sb,
            Qg=si["Qg"] / self.Sb,
            Pl=si["Pl"] / self.Sb,
            Ql=si["Ql"] / self.Sb,
            C=si["C"],
            R=si["R"],
        )
        self.bus_specs[spec.id] = spec
        if spec.C > 0:
            noeud = decl
        elif self.noeuds_sans_capacite == "ideal":
            noeud = IdealBus(spec.id, decl.Un, name=decl.name)
            noeud.spec_si = dict(si)
        else:
            noeud = ResistiveBus(spec.id, decl.Un, name=decl.name)
            noeud.spec_si = dict(si)
        noeud.finaliser(self.wb)
        self.buses[spec.id] = noeud
        return noeud

    def add_bus(self, bus_id: int, Un: float, **kw) -> BusSpec:
        """Declaration d'un bus en PU -- voie historique.

        Preferer ``net.add(Bus(...))``, qui suit les tables MATLAB et prend
        des watts. Ici ``Pg``/``Qg``/``Pl``/``Ql`` sont en **pu base systeme**,
        d'ou la reconversion en SI avant de rejoindre le chemin commun.
        """
        si = dict(kw)
        for cle in ("Pg", "Qg", "Pl", "Ql"):
            if cle in si:
                si[cle] = float(si[cle]) * self.Sb
        self.add(Bus(bus_id, Un, **si))
        return self.bus_specs[int(bus_id)]

    @property
    def all_elements(self) -> list[Element]:
        """Bus d'abord, puis les autres (ordre stable pour les indices).

        Les noeuds continus viennent apres les alternatifs : sur un reseau
        sans continu, la liste est mot pour mot celle d'avant le lot 3, donc
        les indices globaux -- et les matrices assemblees -- aussi.
        """
        return (
            [self.buses[b] for b in sorted(self.buses)]
            + [self.buses_dc[b] for b in sorted(self.buses_dc)]
            + self.elements
        )

    def _bus_index(self) -> dict[int, int]:
        return {b: i for i, b in enumerate(sorted(self.bus_specs))}

    # ------------------------------------------------------------------
    # 2) Load flow
    # ------------------------------------------------------------------
    def build_ybus(self, creux: bool = False):
        """Y-bus du reseau alternatif : dense par defaut, CSC si ``creux``.

        UN SEUL constructeur, :func:`g2egrid.noyau.loadflow.build_ybus` (creux) ; la
        forme dense n'en est que la conversion, pour les appelants qui
        indexent ou modifient la matrice (``Y[i, i] += y``). Le Newton AC
        (:meth:`run_loadflow`) prend le creux.
        """
        idx = self._bus_index()
        n = len(idx)
        branches = []
        shunts = []
        for el in self.elements:
            if len(el.ports) == 2:  # branche
                branches.append(
                    (
                        idx[el.buses[0]],
                        idx[el.buses[1]],
                        el.pdict["R"],
                        el.pdict["X"],
                        el.pdict.get("n", 1.0),
                    )
                )
        # capacites shunt des bus : elles sont physiques, donc dans le Y-bus
        for b, i in idx.items():
            shunts.append((i, 1j * self.bus_specs[b].C))
        # shunts permanents declares par les elements (resistance de
        # causalite d'un noeud sans capacite)
        for el in self.all_elements:
            y = el.shunt_admittance()
            if y is not None:
                shunts.append((idx[el.buses[0]], y))
        # elements qui SONT des shunts (reactances) : leur admittance est au
        # Y-bus, leur courant dans leurs propres etats (Element.admittance_ybus)
        for el in self.elements:
            y = el.admittance_ybus()
            if y is not None:
                shunts.append((idx[el.buses[0]], y))
        Y = ybus_creux(n, branches, shunts)
        return Y if creux else Y.toarray()

    def run_loadflow(self, **kw) -> LoadFlowResult:
        """Le load flow du reseau : alternatif seul, ou mixte s'il y a du continu.

        LA GARDE EST EN TETE, ET LE CORPS EST RESTE TEL QUEL -- pas une
        expression deplacee, pas une variable renommee. C'est delibere : un
        reseau sans noeud continu ne traverse pas une seule instruction
        nouvelle, et la non regression du load flow alternatif n'est alors pas
        une affaire de relecture mais de structure.

        ``tests/test_loadflow_acdc.py`` le verifie ensuite sur chaque reseau du
        catalogue, par deux filets. Le premier compare V a une reference
        produite depuis un export du code d'avant le lot, a 1e-11 pu et au
        meme nombre d'iterations. Le second recalcule les equations du load
        flow depuis le Y-bus, sans reference. Le filet comparait autrefois
        ``lf.V.tobytes()`` octet pour octet, mais l'octet ne survit pas d'une
        machine a l'autre (FR-001). C'est la garantie qui reste depuis que
        les garde-fous de texte du lot 3 sont tombes, et c'est la bonne : ce
        qu'on protege est le resultat, pas la mise en forme.
        """
        if self.bus_specs_dc:
            return self._run_loadflow_acdc(**kw)
        idx = self._bus_index()
        n = len(idx)
        Y = self.build_ybus(creux=True)
        btype = np.zeros(n, dtype=int)
        Vsp = np.ones(n)
        th = np.zeros(n)
        Psp = np.zeros(n)
        Qsp = np.zeros(n)
        for b, i in idx.items():
            s = self.bus_specs[b]
            btype[i], Vsp[i], th[i] = s.btype, s.V, s.theta
            Psp[i] = s.Pg - s.Pl
            Qsp[i] = s.Qg - s.Ql
        self.lf = run_loadflow(Y, btype, Vsp, th, Psp, Qsp, **kw)
        if not self.lf.converged:
            raise RuntimeError(
                f"load flow non converge (mismatch {self.lf.mismatch:.2e})"
            )
        return self.lf

    # ------------------------------------------------------------------
    # 3) Equilibre
    # ------------------------------------------------------------------
    def _element_injections(self) -> dict[int, complex]:
        """Puissance injectee attribuee a chaque element non-Bus (par id)."""
        idx = self._bus_index()
        # part connue (charges) par bus
        known: dict[int, complex] = {b: 0j for b in idx}
        unknown: dict[int, list[Element]] = {b: [] for b in idx}
        sched: dict[int, complex] = {}
        shunts: dict[int, complex] = {}
        for k, el in enumerate(self.elements):
            if len(el.ports) != 1:
                continue
            y = el.admittance_ybus()
            if y is not None:
                # Deja au Y-bus : lf.S la met a la charge des autres elements
                # du bus. On lui rend sa puissance, sans la compter dans
                # ``known`` -- sinon le generateur voisin la fournirait deux fois.
                V = self.lf.V[idx[el.buses[0]]]
                shunts[k] = -np.conj(y) * abs(V) ** 2
                continue
            s = getattr(el, "scheduled_injection", lambda: None)()
            if s is None:
                unknown[el.buses[0]].append(el)
            else:
                sched[k] = s
                known[el.buses[0]] += s
        out: dict[int, complex] = dict(shunts)
        for k, el in enumerate(self.elements):
            if k in sched:
                out[k] = sched[k]
        for b, els in unknown.items():
            if not els:
                continue
            # Injection nette du bus donnee par le load flow. Le Y-bus INCLUT
            # deja la capacite shunt : lf.S est donc la puissance fournie par
            # les elements, capacite comprise. Ne rien y rajouter (sinon la
            # capacite est comptee deux fois et la loi des noeuds sort a x2).
            rest = self.lf.S[idx[b]] - known[b]
            if len(els) > 1:
                raise RuntimeError(
                    f"bus {b}: {len(els)} elements sans injection programmee ; "
                    "impossible d'attribuer le load flow"
                )
            out[self.elements.index(els[0])] = rest
        return out

    # ------------------------------------------------------------------
    # Repere de reference flottant
    # ------------------------------------------------------------------
    def _poids_suivis(self) -> dict[str, float]:
        """``repere_suit`` sous forme ``{nom: poids}``, poids de somme 1."""
        suit = self.repere_suit
        if not suit:
            return {}
        if isinstance(suit, str):
            return {suit: 1.0}
        total = float(sum(suit.values()))
        if total <= 0.0:
            raise ValueError(f"poids du repere de somme {total:g} : rien a suivre")
        return {nom: float(w) / total for nom, w in suit.items()}

    def formeurs(self) -> list[Element]:
        """Les elements qui FONT tourner le reseau : un angle et une vitesse.

        Machines, grid-forming, bus infinis. Un grid-following n'en est pas :
        sa PLL suit la vitesse du reseau, elle ne la fait pas. Un
        :class:`~g2egrid.elements.Repere` non plus : il n'est qu'un choix de
        coordonnees.
        """
        from g2egrid.elements.repere import Repere

        return [
            e
            for e in self.elements
            if "theta_out" in e.svc_outputs
            and self.sortie_vitesse(e) is not None
            and not isinstance(e, Repere)
        ]

    def reference_flottante(
        self,
        suivre: str | dict[str, float] | None = "coi",
        nom: str = "Repere",
    ):
        """Un angle de reference FLOTTANT, detache du slack.

        D'ordinaire ``theta_g``, l'angle que toutes les unites retranchent du
        leur, est l'angle interne de l'element pose au bus slack (regle de V1).
        Le repere est alors LIE a cet element : le declencher, c'est emporter
        le repere avec lui, et le modele d'avant et celui d'apres n'ont plus le
        meme systeme de coordonnees.

        Ici le repere devient un element a part, un :class:`~g2egrid.elements.Repere`
        a un etat, ``dtheta_g/dt = wb.(w_in - 1)``, dont la vitesse ``w_in``
        est :

        * ``"coi"`` (defaut) -- celle du CENTRE D'INERTIE des formeurs
          (:meth:`formeurs`), ponderee par leur inertie ``2H.Sn`` quand ils en
          ont une, et par leur puissance ``Sn`` sinon (grid-forming, bus infini) ;
        * un nom -- celle d'une seule machine ;
        * un dict ``{nom: poids}`` -- une moyenne ponderee quelconque ;
        * ``None`` -- la vitesse nominale : un repere synchrone fixe.

        Le repere n'est qu'un choix de coordonnees : a ``w = 1`` pres, la
        physique n'en depend pas, et le spectre est celui du reseau lie au
        slack, a un mode de jauge pres (le repere a son propre angle absolu,
        qui ne rappelle rien). Mais aucun element n'y est plus indispensable :
        l'element du slack peut etre retire du reseau sans que le repere ne
        change -- le load flow garde un bus slack, le modele dynamique n'en a
        plus besoin.

        Rend le :class:`~g2egrid.elements.Repere` ajoute (ou celui qui existait
        deja sous ce nom, reconfigure).
        """
        from g2egrid.elements.repere import Repere

        rep = next((e for e in self.elements if e.name == nom), None)
        if rep is None:
            rep = self.add(Repere(nom, wb=self.wb))
        elif not isinstance(rep, Repere):
            raise ValueError(f"'{nom}' existe deja et n'est pas un Repere")
        if suivre == "coi":
            poids = {}
            for e in self.formeurs():
                kS = e.pdict.get("kS", 1.0)
                H = e.pdict.get("H")
                poids[e.name] = 2.0 * H * kS if H else kS
            suivre = poids or None
        self.reference = rep.name
        self.repere_suit = suivre
        for e in self.all_elements:  # le repere change : equilibre a refaire
            e.x0 = None
        return rep

    @staticmethod
    def sortie_vitesse(e: Element) -> int | None:
        """Indice, dans ``h``, de la sortie VITESSE d'un element (ou ``None``).

        Une machine l'appelle ``Omega``, un GFM ``w``, un bus infini ``omega``.
        Un grid-following n'en a pas au sens voulu ici : la vitesse de sa PLL
        suit le reseau, elle ne le fait pas tourner.
        """
        nom = next((n for n in ("Omega", "omega", "w") if n in e.svc_outputs), None)
        return None if nom is None else e.n_port_signaux() + e.svc_outputs.index(nom)

    def _reference(self, els) -> tuple[int | None, int | None, int | None]:
        """``(k_ref, ty_ref, tw_ref)`` : l'element qui porte le repere, l'indice
        de sa sortie d'angle ``theta_out`` et celui de sa sortie de vitesse.

        :attr:`reference` l'emporte quand il est pose ; sinon c'est la regle de
        V1 : le premier element non-Bus du bus slack qui expose ``theta_out``.
        """
        if self.reference is not None:
            for m, e in enumerate(els):
                if e.name == self.reference:
                    ty = e.n_port_signaux() + e.svc_outputs.index("theta_out")
                    return m, ty, self.sortie_vitesse(e)
            raise ValueError(f"element de reference '{self.reference}' introuvable")
        bus_slack = next(
            (sp.id for sp in self.bus_specs.values() if sp.btype == SLACK), None
        )
        if bus_slack is None:
            return None, None, None
        for m, e in enumerate(els):
            if isinstance(e, Bus) or bus_slack not in e.buses:
                continue
            if "theta_out" in e.svc_outputs:
                ty = e.n_port_signaux() + e.svc_outputs.index("theta_out")
                return m, ty, self.sortie_vitesse(e)
        return None, None, None

    def _element_reference(self, els=None):
        """Element qui porte l'angle du repere global, et l'indice de sa sortie.

        Rend ``(indice dans els, indice de la sortie)``, ou ``(None, None)``
        s'il n'y en a pas. Voir :meth:`_reference`.
        """
        if els is None:
            els = self.all_elements
        k, ty, _tw = self._reference(els)
        return k, ty

    def _angle_reference(self, idx, inj) -> float:
        """Valeur d'equilibre de ``theta_g``, l'angle du repere global.

        C'est l'angle INTERNE de l'element de reference, pas celui de son bus :
        pour un bus infini, l'angle de sa f.e.m. ; pour une machine, son angle
        rotor. C'est la regle de V1 (``theta_g_0 = alpha_s``, ou
        ``delta_s0 - pi/2`` quand le slack est une machine), et c'est aussi la
        seule valeur coherente avec le cablage : ``theta_g`` est branche sur la
        sortie ``theta_out`` de cet element, donc son point de fonctionnement
        doit valoir ce que cette sortie vaut a l'equilibre.

        L'element de reference est initialise une fois a ``theta_g = 0`` pour
        lire cette sortie. Aucun des modeles portes n'en depend pour son propre
        equilibre, donc une passe suffit.
        """
        k_ref, ty_ref = self._element_reference()
        if k_ref is None:
            return 0.0
        el = self.all_elements[k_ref]
        cm = el.compiled()
        ctx = self._contexte_init(el, idx, inj, 0.0)
        init = el.initialize(ctx)
        x0 = np.array([init[n] for n in cm.spec.states], dtype=float)
        z0 = np.array([init[n] for n in cm.spec.algebs], dtype=float)
        u0 = el.build_u0(init, ctx)
        y = (
            cm.hx(x0, z0, u0, el.p) @ x0
            + cm.hz(x0, z0, u0, el.p) @ z0
            + cm.hu(x0, z0, u0, el.p) @ u0
        )
        return float(y[ty_ref])

    def _contexte_init(self, el, idx, inj, theta_g: float, brut: bool = False):
        """Contexte d'initialisation d'un element, load flow et repere compris.

        ``brut=True`` ignore :attr:`Element.point_init`. Le point impose sert a
        rejouer les FORMULES d'equilibre de V1, qui lisent le bus que
        ``network_form.m`` a cree derriere le transfo ; le point de
        linearisation, lui, doit rester celui du vrai noeud reseau, sans quoi
        les tensions de port du modele ne seraient pas celles du reseau.
        """
        if isinstance(el, BusDC):
            return InitContext(
                v={"_self": self._tension_dc(el.bus)},
                w=self.w_sys,
                f0=self.f0,
                theta_g=theta_g,
            )
        if isinstance(el, Bus):
            return InitContext(
                v={"_self": self.lf.V[idx[el.bus]]},
                w=self.w_sys,
                f0=self.f0,
                theta_g=theta_g,
            )
        k = self.elements.index(el)
        v = {pt: self.lf.V[idx[b]] for pt, b in zip(el.ports, el.buses)}
        s = {}
        if len(el.ports) == 1 and k in inj:
            s = {el.ports[0]: inj[k]}
        if el.point_init is not None and len(el.ports) == 1 and not brut:
            v_imp, s_imp = el.point_init
            v = {el.ports[0]: complex(v_imp)}
            s = {el.ports[0]: complex(s_imp)}
        # ports continus : une tension REELLE, celle que porte ``v_dc``, et la
        # puissance que l'element connait de lui-meme (aucun load flow continu
        # ne la lui attribue avant le lot 4). Ajoutes APRES ``point_init``, qui
        # ne parle que du cote alternatif et reconstruit ``v`` de zero.
        p_dc = el.scheduled_injection_dc() if el.ports_dc else None
        for pk, pt in enumerate(el.ports_dc):
            v[pt] = self._tension_dc(el.buses[len(el.ports) + pk])
            if p_dc is not None:
                s[pt] = p_dc
        return InitContext(v=v, s=s, w=self.w_sys, f0=self.f0, theta_g=theta_g)

    # ------------------------------------------------------------------
    # 3 bis) Equilibre du reseau CONTINU
    # ------------------------------------------------------------------
    def _tension_dc(self, bus: int) -> float:
        """Tension d'essai du noeud continu ``bus``, avant ou apres Newton."""
        if bus in self.v_dc:
            return self.v_dc[bus]
        return self.bus_specs_dc[bus].V

    def _sortie_dc(self, el, ctx, k: int) -> float:
        """Courant que ``el`` injecte a son k-ieme noeud continu, au point ``ctx``.

        Le point vient des SEULES formules d'``initialize``, sans le Newton
        interne de :meth:`_point_element`. Ce n'est pas une approximation :
        l'equilibre analytique des convertisseurs a port continu est EXACT.
        Mesure : le Newton interne ne deplace pas d'un bit le ``x0`` d'un
        ``GFM3DC``, et de 2,2e-16 celui d'un ``GFL2DC``. Ce qu'on y gagne est
        le cout -- le Newton unifie evalue cette fonction quatre a six fois par
        convertisseur et par iteration.

        La sortie est LUE dans ``h``, jamais reconstruite : le courant continu
        d'un convertisseur vaut ``-(Gdc.v + md.isd + mq.isq)``, l'aboutissement
        de ses treize equations, et le reecrire ici serait un second modele a
        tenir en accord avec le premier.
        """
        cm = el.compiled()
        init = el.initialize(ctx)
        x0 = np.array([init[n] for n in cm.spec.states], dtype=float)
        z0 = np.array([init[n] for n in cm.spec.algebs], dtype=float)
        u0 = el.build_u0(init, ctx)
        return float(cm.h(x0, z0, u0, el.p)[el.port_dc_output_idx(k)])

    def _couplage_vsc(self, el):
        """La fonction ``p_dc(vm, p_ac, q_ac, v_dc)`` d'un convertisseur.

        C'EST LE COEUR DU LOT. La puissance qu'un VSC rend a son noeud continu
        n'est pas ``-P_ac / rendement`` : c'est ``-V_dc . i_port``, avec le
        ``i_port`` du MODELE DYNAMIQUE. Les pertes qui en sortent valent
        ``R1.|is|^2 + Gdc.V_dc^2`` -- celles du filtre et celles du
        condensateur, aux valeurs que le modele porte -- et non un coefficient
        libre que le load flow et la simulation interpreteraient differemment.

        C'est ce qui fait du point de fonctionnement mixte un POINT FIXE du
        modele non lineaire : ``residu_equilibre`` le mesure a 1e-8.

        ``theta_g = 0`` et une tension alternative posee sur l'axe reel : le
        courant continu est invariant par rotation du repere (verifie,
        ``i_dc(theta_g + 0.7) - i_dc`` vaut exactement zero), donc la phase ne
        lui manque pas.
        """

        def p_dc(vm, p_ac, q_ac, v_dc, _el=el):
            ctx = InitContext(
                v={_el.ports[0]: complex(vm, 0.0), _el.ports_dc[0]: float(v_dc)},
                s={_el.ports[0]: complex(p_ac, q_ac)},
                w=self.w_sys,
                f0=self.f0,
                theta_g=0.0,
            )
            return float(v_dc) * self._sortie_dc(_el, ctx, 0)

        return p_dc

    def _injection_dc_propre(self, el):
        """La fonction ``p(v)`` d'un element continu a un seul port.

        Source, charge resistive, charge a puissance constante : tous se
        resument, vus du load flow, a une puissance injectee fonction de la
        seule tension de leur noeud. On la demande au modele plutot que de la
        reecrire, pour la meme raison que :meth:`_couplage_vsc`.
        """

        def p(v, _el=el):
            ctx = InitContext(
                v={_el.ports_dc[0]: float(v)},
                w=self.w_sys,
                f0=self.f0,
                theta_g=0.0,
            )
            return float(v) * self._sortie_dc(_el, ctx, 0)

        return p

    def _elements_vsc(self) -> list:
        """Les convertisseurs : un element qui a des ports des DEUX genres."""
        return [el for el in self.elements if el.ports and el.ports_dc]

    def _reseau_dc(self, pos: dict[int, int]):
        """``(G, dctype, Vdc_sp, Pdc_sp, p_bus)`` du reseau continu.

        ``G`` est le laplacien des conductances des cables : le pendant exact
        du Y-bus, mais reel -- un reseau continu n'a ni angle ni susceptance.
        La capacite des cables n'y figure pas ; elle est deja versee par
        moities aux deux noeuds, ou elle est un ETAT et non une admittance.

        ``p_bus[j]`` rassemble ce que les elements a un seul port continu
        injectent au noeud ``j``, sous forme d'une fonction de sa tension.
        """
        n = len(pos)
        G = np.zeros((n, n))
        fonctions: list[list] = [[] for _ in range(n)]
        for el in self.elements:
            if not el.ports_dc or el.ports:
                continue
            if isinstance(el, LineDC):
                f, t = pos[el.buses[0]], pos[el.buses[1]]
                g = 1.0 / el.pdict["R"]
                G[f, f] += g
                G[t, t] += g
                G[f, t] -= g
                G[t, f] -= g
                continue
            if len(el.ports_dc) != 1:
                raise ValueError(
                    f"{el.name}: {len(el.ports_dc)} ports continus et aucun port "
                    "alternatif ; le load flow continu ne sait traiter que les "
                    "cables (deux ports) et les elements a un seul port"
                )
            fonctions[pos[el.buses[0]]].append(self._injection_dc_propre(el))

        dctype = np.zeros(n, dtype=int)
        Vdc_sp = np.ones(n)
        Pdc_sp = np.zeros(n)
        for b, j in pos.items():
            s = self.bus_specs_dc[b]
            dctype[j] = s.btype
            Vdc_sp[j] = s.V
            Pdc_sp[j] = s.Pinj

        def somme(fs):
            return (lambda v, _fs=fs: sum(f(v) for f in _fs)) if fs else None

        return G, dctype, Vdc_sp, Pdc_sp, [somme(fs) for fs in fonctions]

    def _specs_vsc(self, idx, pos, btype, Psp, Qsp) -> list[VSCSpec]:
        """Les :class:`~g2egrid.noyau.loadflow.VSCSpec`, un par convertisseur.

        LE MODE PAR DEFAUT est le premier de ``modes_dc`` qui ne tient PAS la
        tension continue -- ``PV`` pour un grid-forming, ``PQ`` pour un
        grid-following. Ce n'est pas un choix arbitraire : c'est le mode ou le
        convertisseur se comporte exactement comme au lot 3, sa puissance
        alternative venant du load flow et la tension continue sortant du bilan
        du noeud. Un reseau bati a la main, sans ``T_VSCs``, retrouve donc son
        point de fonctionnement d'avant.
        """
        out: list[VSCSpec] = []
        programmees = self._injections_programmees()
        for el in self._elements_vsc():
            if len(el.ports) != 1 or len(el.ports_dc) != 1:
                raise ValueError(
                    f"{el.name}: le load flow mixte ne traite que les "
                    "convertisseurs a un port alternatif et un port continu"
                )
            b_ac, b_dc = el.buses[0], el.buses[1]
            i = idx[b_ac]
            modes = getattr(el, "modes_dc", ())
            if not modes:
                raise ValueError(
                    f"{el.name}: la classe {type(el).__name__} ne declare aucun "
                    "mode continu (attribut de classe 'modes_dc')"
                )
            # sur le slack, un convertisseur qui sait former son ile la forme
            if btype[i] == SLACK and "Vtheta" in modes:
                defaut = "Vtheta"
            else:
                defaut = next(m for m in modes if not m.startswith("Vdc"))
            mode = self.vsc_modes.get(b_ac, defaut)
            if btype[i] == SLACK and mode != "Vtheta":
                raise ValueError(
                    f"{el.name}: son bus alternatif {b_ac} est le slack. Ses "
                    "deux puissances y seraient libres, alors qu'un VSC en tient "
                    "exactement deux : le declarer sur un bus PV ou PQ, ou, s'il "
                    "forme seul son ile alternative, en mode Vtheta."
                )
            if mode == "Vtheta" and btype[i] != SLACK:
                raise ValueError(
                    f"{el.name}: mode Vtheta (il forme son ile alternative), "
                    f"mais son bus {b_ac} n'est pas le slack de l'ile : le "
                    "declarer de type 1."
                )
            if mode not in modes:
                raise ValueError(
                    f"{el.name}: mode '{mode}' impossible pour "
                    f"{type(el).__name__} ; modes admis : {', '.join(modes)}"
                )
            if mode.startswith("Vdc") and self.bus_specs_dc[b_dc].btype != DC_SLACK:
                raise ValueError(
                    f"{el.name}: mode '{mode}', mais le noeud continu {b_dc} est "
                    f"de type {self.bus_specs_dc[b_dc].btype}. Un convertisseur "
                    "qui tient la tension de son noeud en fait le slack de son "
                    f"ile : declarer le bus {b_dc} de type {DC_SLACK}."
                )
            pente, vdc_ref = 0.0, 1.0
            if mode in ("droop", "droopI"):
                k = self.vsc_statisme.get(b_ac)
                if not k or k <= 0:
                    raise ValueError(
                        f"{el.name}: mode {mode} sans statisme : renseigner "
                        "'kdroop' (pu de tension continue par pu de puissance, "
                        "base Sn) dans PARAMS_PAR_BUS"
                    )
                if self.bus_specs_dc[b_dc].btype == DC_SLACK:
                    raise ValueError(
                        f"{el.name}: mode droop, mais le noeud continu {b_dc} est "
                        "de type 1 : un statisme laisse la tension libre, le "
                        "declarer de type 2 (sa V_pu est la consigne)."
                    )
                pente = (float(el.Sn) / self.Sb) / float(k)
                vdc_ref = float(self.bus_specs_dc[b_dc].V)
            kqv, v_ref_ac = 0.0, 1.0
            nq = self.vsc_qv.get(b_ac)
            if nq:
                if btype[i] != PQ:
                    raise ValueError(
                        f"{el.name}: statisme Q-V sur le bus {b_ac}, qui n'est pas "
                        "PQ : un bus PV tient deja sa tension, le slack forme son "
                        "ile ; le declarer de type 3"
                    )
                kqv = (float(el.Sn) / self.Sb) / float(nq)
                v_ref_ac = float(self.bus_specs[b_ac].V)
            p_aut, q_aut = programmees[b_ac]
            out.append(
                VSCSpec(
                    iac=i,
                    idc=pos[b_dc],
                    mode=mode,
                    p_dc=self._couplage_vsc(el),
                    p_ac=float(Psp[i] - p_aut),
                    q_ac=float(Qsp[i] - q_aut),
                    q_autres=float(q_aut),
                    nom=el.name,
                    p_autres=float(p_aut),
                    pente=pente,
                    vdc_ref=vdc_ref,
                    kqv=kqv,
                    v_ref_ac=v_ref_ac,
                )
            )
        return out

    def _injections_programmees(self) -> dict[int, tuple[float, float]]:
        """Par bus alternatif, la puissance que les elements PROGRAMMES y posent.

        La meme part ``known`` que :meth:`_element_injections` retranche avant
        d'attribuer le reliquat : c'est elle qui, otee de l'injection nette du
        bus, donne la puissance du convertisseur lui-meme.
        """
        out = {b: (0.0, 0.0) for b in self.bus_specs}
        for el in self.elements:
            if len(el.ports) != 1 or el.ports_dc or el.admittance_ybus() is not None:
                continue  # une reactance est au Y-bus, pas dans ``known``
            s = getattr(el, "scheduled_injection", lambda: None)()
            if s is None:
                continue
            p, q = out[el.buses[0]]
            out[el.buses[0]] = (p + complex(s).real, q + complex(s).imag)
        return out

    def _iles_dc(self, pos: dict[int, int], G: np.ndarray) -> list[list[int]]:
        """Les composantes connexes du graphe continu, en numeros de bus."""
        ids = sorted(pos)
        vus: set[int] = set()
        iles = []
        for b in ids:
            if b in vus:
                continue
            pile, ile = [b], []
            vus.add(b)
            while pile:
                c = pile.pop()
                ile.append(c)
                for d in ids:
                    if d not in vus and G[pos[c], pos[d]] != 0.0:
                        vus.add(d)
                        pile.append(d)
            iles.append(sorted(ile))
        return iles

    def _verifier_iles_dc(self, pos, G, vscs) -> None:
        """Une seule reference de tension par ile continue -- ou rien a dire.

        REGLE REPRISE DE VeraGrid (l'idee, pas le code) : un convertisseur est,
        vu du continu, un terminal a PUISSANCE constante. Une ile qui n'en
        contient que de ceux-la n'a aucune equation qui fixe son NIVEAU de
        tension, et son bilan de puissance ne boucle que par miracle -- les
        pertes des cables ne sont a la charge de personne. Il lui faut donc
        exactement une reference : un VSC en mode ``Vdc...``, ou une source
        continue.

        Une ile SANS convertisseur, elle, n'a rien a prouver : ses elements
        (resistance, source derriere R) ont une injection qui depend de leur
        tension, et fixent le niveau d'eux-memes. C'est le cas des reseaux
        purement continus du lot 3, qui traversent cette garde sans la voir.

        Une ile qui n'est pas reliee aux autres est simplement une ile de plus,
        et elle est jugee sur le meme critere : c'est ainsi qu'un reseau
        continu coupe en deux se signale ici, en nommant la moitie orpheline.
        """
        vsc_par_bus: dict[int, list] = {b: [] for b in pos}
        inv = {j: b for b, j in pos.items()}
        for c in vscs:
            vsc_par_bus[inv[c.idc]].append(c)
        tenue = {b: [] for b in pos}
        # une batterie a statisme PARTAGE le reglage de tension : elle suffit
        # a une ile qui n'a pas d'autre reference, sans en exclure une
        statisme = {b: [] for b in pos}
        for el in self.elements:
            if isinstance(el, SourceDC):
                tenue[el.buses[0]].append(el.name)
            elif getattr(el, "reference_tension_dc", False):
                statisme[el.buses[0]].append(el.name)

        iles = self._iles_dc(pos, G)
        for ile in iles:
            conv = [c for b in ile for c in vsc_par_bus[b]]
            if not conv:
                continue
            refs = [c.nom for c in conv if c.tient_vdc]
            refs += [nom for b in ile for nom in tenue[b]]
            # un statisme fixe lui aussi le niveau de tension (par sa consigne)
            refs_statisme = [c.nom for c in conv if c.statisme]
            ou = (
                f"ile continue {{{', '.join(str(b) for b in ile)}}}"
                if len(iles) > 1
                else f"reseau continu {{{', '.join(str(b) for b in ile)}}}"
            )
            if not refs and (refs_statisme or any(statisme[b] for b in ile)):
                continue
            if not refs:
                raise ValueError(
                    f"{ou} : aucune reference de tension. Ses "
                    f"{len(conv)} convertisseur(s) y sont tous a puissance "
                    "imposee, donc rien n'y fixe le niveau de tension. Passer "
                    "l'un d'eux en mode VdcQ ou VdcV, ou y raccorder une "
                    "source continue ou une batterie (convertisseur DC/DC)."
                )
            # PLUSIEURS references sont admises, a des noeuds distincts : comme
            # plusieurs bus PV du cote alternatif, chacune tient SA tension et
            # les cables portent la difference (CIGRE B4 : Cb-A1 a 1,01 pu,
            # Cm-B2 a 1,00 pu, dans le meme reseau continu). Deux references
            # au MEME noeud, elles, se contrediraient.
            par_noeud = Counter(inv[c.idc] for c in conv if c.tient_vdc)
            for b in ile:
                par_noeud[b] += len(tenue[b])
            doubles = [b for b, n in par_noeud.items() if n > 1]
            if doubles:
                raise ValueError(
                    f"{ou} : {par_noeud[doubles[0]]} references de tension au "
                    f"noeud continu {doubles[0]} ; un noeud n'en tient qu'une "
                    "(des references a des noeuds DISTINCTS, elles, sont admises)."
                )

        # Un noeud de type 1 sans tenant (ni VSC en mode Vdc..., ni source, ni
        # batterie) serait un slack IDEAL du load flow -- sa tension est
        # imposee, son bilan jamais ecrit -- alors que le modele dynamique n'y
        # met aucune source : le point calcule n'y serait pas un equilibre.
        # On refuse plutot que d'avertir : le resultat serait faux, pas
        # seulement discutable.
        for b in sorted(pos):
            tient = any(c.tient_vdc for c in vsc_par_bus[b]) or tenue[b] or statisme[b]
            if self.bus_specs_dc[b].btype == DC_SLACK and not tient:
                raise ValueError(
                    f"noeud continu {b} : declare de type 1 (slack) mais rien "
                    "ne le tient -- ni convertisseur en mode Vdc..., ni source "
                    "continue, ni batterie. Le load flow y poserait un slack "
                    "ideal que le modele dynamique n'a pas. Le declarer de "
                    "type 2 (passif) ou y raccorder un tenant de tension."
                )

    def _probleme_acdc(self):
        """Le systeme non lineaire du load flow mixte, monte mais PAS resolu.

        Le cote alternatif est monte EXACTEMENT comme dans
        :meth:`run_loadflow` -- meme Y-bus, memes types, memes consignes. Ce
        qui s'y ajoute ne le modifie pas, il s'y accroche : les noeuds continus
        apportent leurs inconnues et leurs equations en QUEUE du systeme, et
        les convertisseurs une colonne ``-1`` sur la ligne de leur bus. Un
        reseau sans noeud continu ne passe jamais ici.

        Separee de :meth:`_run_loadflow_acdc` pour que le jacobien soit
        eprouvable seul, sur le vrai reseau.
        """
        idx = self._bus_index()
        n = len(idx)
        Y = self.build_ybus()
        btype = np.zeros(n, dtype=int)
        Vsp = np.ones(n)
        th = np.zeros(n)
        Psp = np.zeros(n)
        Qsp = np.zeros(n)
        for b, i in idx.items():
            s = self.bus_specs[b]
            btype[i], Vsp[i], th[i] = s.btype, s.V, s.theta
            Psp[i] = s.Pg - s.Pl
            Qsp[i] = s.Qg - s.Ql

        pos = {b: j for j, b in enumerate(sorted(self.bus_specs_dc))}
        G, dctype, Vdc_sp, Pdc_sp, p_bus = self._reseau_dc(pos)
        vscs = self._specs_vsc(idx, pos, btype, Psp, Qsp)
        self._verifier_iles_dc(pos, G, vscs)
        return ProblemeACDC(
            Y, btype, Vsp, th, Psp, Qsp, G, dctype, Vdc_sp, Pdc_sp, vscs, p_bus_dc=p_bus
        )

    def _run_loadflow_acdc(self, **kw) -> LoadFlowACDCResult:
        """Resout le systeme de :meth:`_probleme_acdc` et range son resultat."""
        pos = {b: j for j, b in enumerate(sorted(self.bus_specs_dc))}
        self.lf = self._probleme_acdc().resoudre(**kw)
        self.v_dc = {b: float(self.lf.Vdc[j]) for b, j in pos.items()}
        if not self.lf.converged:
            raise RuntimeError(self._diagnostic_dc(pos))
        return self.lf

    def _diagnostic_dc(self, pos) -> str:
        """Le message d'echec du Newton mixte, rendu exploitable.

        La cause la plus frequente n'est pas numerique mais PHYSIQUE : un
        convertisseur se comporte, vu du continu, comme une charge a puissance
        constante (``i = -P/v``), et une source derriere une resistance ``R``
        ne peut en alimenter qu'au plus ``E^2/(4R)``. Au-dela, il n'existe
        aucune tension d'equilibre, et Newton n'a rien a trouver -- d'ou le
        rappel des tensions ou il s'est arrete, qui montrent l'effondrement,
        et la piste a verifier en premier.
        """
        entete = (
            f"load flow AC/DC non converge (mismatch {self.lf.mismatch:.2e} "
            f"apres {self.lf.iterations} iterations)"
        )
        lignes = [entete]
        if getattr(self.lf, "message", ""):
            lignes.append(f"  {self.lf.message}")
        for b in sorted(pos):
            lignes.append(f"  bus continu {b}: v = {self.v_dc[b]:+.4f} pu")
        lignes.append(
            "  piste : un port continu de convertisseur se comporte en charge "
            "a puissance constante ; verifier que la source continue peut la "
            "fournir (P_max = E^2 / 4R pour une f.e.m. derriere R)"
        )
        return "\n".join(lignes)

    def _element_injections_dc(self) -> dict[int, float]:
        """Puissance continue attribuee a chaque element, apres le load flow.

        Le pendant exact de :meth:`_element_injections`, cote continu : les
        elements qui connaissent leur puissance la declarent, et le RELIQUAT du
        noeud revient a celui qui ne la connait pas. Comme du cote alternatif,
        deux pretendants sur un meme noeud sont une erreur, pas un partage
        arbitraire.

        Une difference assumee : le reliquat n'a de sens qu'au SLACK continu.
        Sur un noeud de type 2, la tension est l'inconnue et le bilan est deja
        referme par le Newton -- il n'y a rien a attribuer, et chaque element y
        tient sa propre puissance de sa propre tension.
        """
        connu = {b: 0.0 for b in self.bus_specs_dc}
        inconnu: dict[int, list[int]] = {b: [] for b in self.bus_specs_dc}
        out: dict[int, float] = {}
        for k, el in enumerate(self.elements):
            if el.ports or len(el.ports_dc) != 1:
                continue
            b = el.buses[0]
            p = el.scheduled_injection_dc()
            if p is None and not (isinstance(el, SourceDC) or hasattr(el, "caler_dc")):
                # ni programmee ni ajustable (charge resistive) : sa puissance
                # est celle que lui impose la tension de son noeud
                if b in self.v_dc:
                    p = self._injection_dc_propre(el)(self.v_dc[b])
            if p is None:
                inconnu[b].append(k)
            else:
                out[k] = float(p)
                connu[b] += float(p)
        if self.lf is None or not isinstance(self.lf, LoadFlowACDCResult):
            return out

        pos = {b: j for j, b in enumerate(sorted(self.bus_specs_dc))}
        G, _t, _v, Pdc_sp, _p = self._reseau_dc(pos)
        v = np.array([self.v_dc[b] for b in sorted(pos)])
        sortant = v * (G @ v)
        for b, ks in inconnu.items():
            if not ks or self.bus_specs_dc[b].btype != DC_SLACK:
                continue
            if len(ks) > 1:
                raise RuntimeError(
                    f"noeud continu {b}: {len(ks)} elements sans injection "
                    "programmee ; impossible d'attribuer le load flow"
                )
            j = pos[b]
            out[ks[0]] = float(sortant[j] - Pdc_sp[j] - self.lf.Pdc[j] - connu[b])
        return out

    def _caler_sources_dc(self) -> None:
        """La f.e.m. des sources continues du SLACK, deduite du load flow.

        Le pendant de ce que fait deja ``IB.initialize`` du cote alternatif :
        au bus slack, la tension est une donnee et la f.e.m. s'en deduit, pas
        l'inverse. Un noeud continu de type 1 declare la meme chose -- sa
        tension est imposee -- donc la source qui l'alimente doit rendre
        exactement le reliquat, et son ``E`` suit.

        Cela n'arrive QUE sur un noeud de type 1, qui se declare explicitement.
        Un noeud de type 2 laisse la source avec le ``E`` qu'on lui a donne, et
        c'est sa tension qui s'ajuste : c'est le comportement du lot 3, et il
        ne bouge pas d'un chiffre.
        """
        attrib = self._element_injections_dc()
        for k, el in enumerate(self.elements):
            if k not in attrib:
                continue
            if hasattr(el, "caler_dc"):  # convertisseur DC/DC d'une batterie
                if self.bus_specs_dc[el.buses[0]].btype == DC_SLACK:
                    el.caler_dc(self.v_dc[el.buses[0]], attrib[k])
                continue
            if not isinstance(el, SourceDC):
                continue
            if self.bus_specs_dc[el.buses[0]].btype != DC_SLACK:
                continue
            v = self.v_dc[el.buses[0]]
            if abs(v) < 1e-12:
                raise RuntimeError(f"{el.name}: tension continue nulle au slack")
            i = attrib[k] / v / el.pdict["kS"]
            el.E = float(v + el.pdict["R"] * i)

    def _point_element(self, el, idx, inj, theta_g: float, tol: float):
        """``(x0, z0, u0)`` d'un element : initialisation analytique puis Newton.

        Extrait de :meth:`solve_equilibrium` parce que le Newton continu en a
        besoin lui aussi : pour connaitre le courant qu'un convertisseur rend
        a son noeud continu, il faut d'abord resoudre l'equilibre INTERNE de
        ce convertisseur a la tension d'essai. Le corps n'a pas bouge.
        """
        cm = el.compiled()
        ctx = self._contexte_init(el, idx, inj, theta_g)

        init = el.initialize(ctx)
        x0 = np.array([init[n] for n in cm.spec.states], dtype=float)
        z0 = np.array([init[n] for n in cm.spec.algebs], dtype=float)
        u0 = el.build_u0(init, self._contexte_init(el, idx, inj, theta_g, brut=True))

        if (cm.n_x + cm.n_z) and el.raffiner_equilibre:

            def res(w, _cm=cm, _u=u0, _p=el.p, _nx=cm.n_x):
                return _cm.residual(w[:_nx], w[_nx:], _u, _p)

            def jac(w, _cm=cm, _u=u0, _p=el.p, _nx=cm.n_x):
                x, z = w[:_nx], w[_nx:]
                top = np.hstack((_cm.fx(x, z, _u, _p), _cm.fz(x, z, _u, _p)))
                bot = np.hstack((_cm.gx(x, z, _u, _p), _cm.gz(x, z, _u, _p)))
                return np.vstack((top, bot))

            w0 = np.concatenate((x0, z0))
            sol = root(res, w0, jac=jac, method="hybr", tol=tol)
            if not sol.success:
                # Powell (hybr) s'arrete quand l'initialisation analytique est
                # loin -- un GFL dont la PLL part a quelques dizaines de degres
                # de son point, par exemple. Levenberg-Marquardt, plus lent
                # mais plus patient, part du meme point ; on ne l'essaie
                # qu'en second, pour que ce qui convergeait ne bouge pas.
                sol = root(res, w0, jac=jac, method="lm", tol=tol)
            if not sol.success:
                raise RuntimeError(f"{el.name}: equilibre non converge — {sol.message}")
            x0, z0 = sol.x[: cm.n_x], sol.x[cm.n_x :]
            r = (
                np.max(np.abs(cm.residual(x0, z0, u0, el.p)))
                if cm.n_x + cm.n_z
                else 0.0
            )
            if r > 1e-6:
                raise RuntimeError(f"{el.name}: residu d'equilibre {r:.2e}")
        return x0, z0, u0

    def solve_equilibrium(
        self, tol: float = 1e-9, pseudo_transitoire: bool = False
    ) -> None:
        """Le point de fonctionnement de chaque element, autour du load flow.

        IL N'Y A PLUS DE NEWTON CONTINU A PART. Le lot 3 en avait un, qui
        tournait ici sur les seules tensions continues ; le lot 4 l'a absorbe
        dans :meth:`_run_loadflow_acdc`, ou les tensions continues sont des
        inconnues comme les autres. La seule trace qu'il en reste est la garde
        ci-dessous : si le load flow a ete pose a la main (``net.lf = ...``)
        sans passer par le Newton unifie, les tensions continues manquent
        encore, et il faut les calculer avant de figer quoi que ce soit --
        sinon chaque convertisseur s'equilibrerait a une tension de port fausse
        et le modele serait linearise a cote de son point de fonctionnement.

        LES LIMITES DES CONVERTISSEURS sont verifiees au point trouve
        (:mod:`g2egrid.noyau.limites`) : chaque limiteur dont la grandeur libre
        depasse sa limite (``|i_s| > Imax``, ``|m| > Mmax``) est range dans
        :attr:`depassements` et signale par un avertissement
        :class:`~g2egrid.noyau.limites.DepassementLimite`. Le calcul continue : le
        modele linearise reste celui du mode libre, que l'on peut explorer en
        le sachant.
        """
        if self.lf is None:
            self.run_loadflow()
        if self.bus_specs_dc and not self.v_dc:
            self._run_loadflow_acdc()
        if self.bus_specs_dc:
            self._caler_sources_dc()
        idx = self._bus_index()
        inj = self._element_injections()

        theta_g = self._angle_reference(idx, inj)

        for el in self.all_elements:
            el.x0, el.z0, el.u0 = self._point_element(el, idx, inj, theta_g, tol)

        if pseudo_transitoire:
            self._equilibre_global()

        self.depassements = limites.depassements(self)
        limites.signaler(self.depassements)

    def _equilibre_global(self, tol: float = 1e-9) -> float:
        """Le point de chaque element, recale sur un equilibre du RESEAU ENTIER.

        L'initialisation par element est explicite : chacun calcule son point
        a partir du load flow. Elle donne un point fixe du modele complet quand
        tout est coherent -- c'est le cas du catalogue, verifie a 1e-8 --, mais
        pas quand un modele rejoue une convention qui ne verifie pas ses propres
        equations (les modes ``init='symSG'`` de V1), ni pour un reseau bati a
        la main dont une hypothese ne tient pas. La continuation
        pseudo-transitoire (:func:`g2egrid.noyau.integrateur.equilibre_pseudo_transitoire`,
        l'idee de VeraGrid) part de ce point et laisse le modele s'equilibrer
        lui-meme, consignes fixes. Rend le residu atteint.
        """
        from g2egrid.noyau.integrateur import equilibre_pseudo_transitoire
        from g2egrid.temporel import ModeleNonLineaire

        m = ModeleNonLineaire(self)
        x, z, u, r = equilibre_pseudo_transitoire(m, tol=tol)
        for el, xs, zs, us in zip(m.c.els, m._x(x), m._z(z), m._u(u)):
            el.x0, el.z0, el.u0 = np.array(xs), np.array(zs), np.array(us)
        return r

    # ------------------------------------------------------------------
    # 4) Assemblage du modele d'etat global
    # ------------------------------------------------------------------
    def cablage(self, injections: tuple[int, ...] = ()) -> Cablage:
        """Le cablage du reseau, independamment de son point de fonctionnement.

        ``injections`` : des bus ou l'on ajoute une entree externe de courant
        INJECTE (``inj<bus>.id``, ``inj<bus>.iq``), versee dans la loi des
        noeuds du bus -- de quoi lire l'impedance du reseau vue de ce bus
        (:mod:`g2egrid.analyse.impedance`). Voir :class:`Cablage`. Cette methode indexe les elements et construit
        ``Ruext`` et ``Ry`` ; c'est :meth:`build_state_space` qui evalue
        ensuite les jacobiens au point de fonctionnement, et
        :mod:`g2egrid.temporel` qui garde les equations non lineaires.
        """
        els = self.all_elements
        cms = [el.compiled() for el in els]
        for el in self.elements:
            self._verifier_genre_ports(el)

        # --- index globaux ---
        ox, oz, ou, oy = [], [], [], []
        nx = nz = nu = ny = 0
        state_names: list[str] = []
        for el, cm in zip(els, cms):
            ox.append(nx)
            nx += cm.n_x
            oz.append(nz)
            nz += cm.n_z
            ou.append(nu)
            nu += cm.n_u
            oy.append(ny)
            ny += cm.n_y
            state_names += [f"{el.name}.{n}" for n in cm.spec.states]

        # --- routage des signaux ---
        uext_names: list[str] = [W_REF]
        uext_idx: dict[str, int] = {W_REF: 0}
        suivi = self.repere_suit or None
        for el in els:
            for s in el.svc_inputs:
                nm = f"{el.name}.{s}"
                if el is not None and isinstance(el, Bus) and s.startswith("isum"):
                    continue  # pilote par la loi des noeuds
                if suivi is not None and el.name == self.reference and s == "w_in":
                    continue  # pilote par la vitesse de la machine suivie
                uext_idx[nm] = len(uext_names)
                uext_names.append(nm)

        Ruext = np.zeros((nu, len(uext_names)))
        Ry = np.zeros((nu, ny))
        bus_pos = {el.bus: k for k, el in enumerate(els) if isinstance(el, Bus)}

        # --- element de reference d'angle -------------------------------
        # theta_g est l'angle de l'element attache au bus slack (symSG et
        # consorts le recoivent en entree globale et tournent de
        # theta_g - theta). L'element de reference se le recoit LUI AUSSI :
        # son theta_g - theta vaut alors zero identiquement, si bien que sa
        # f.e.m. reste fixe dans le repere qu'il definit -- exactement ce que
        # fait symIB dans V1. Sans ce bouclage, l'angle de reference tournerait
        # deux fois : une fois sur la source, une fois sur tous les autres.
        # Le bouclage n'est pas algebrique : theta_out est un ETAT, donc le
        # chemin theta_g -> theta_out n'a pas de transmission directe.
        k_ref, ty_ref, tw_ref = self._reference(els)
        # un Repere qui suit une machine, ou la moyenne ponderee de plusieurs
        # (centre d'inertie) : son entree w_in <- somme des poids x vitesses
        suivis: list[tuple[int, int, float]] = []
        for nom, poids in self._poids_suivis().items():
            k_s = next((m for m, e in enumerate(els) if e.name == nom), None)
            if k_s is None:
                raise ValueError(f"machine suivie '{nom}' introuvable")
            tw_s = self.sortie_vitesse(els[k_s])
            if tw_s is None:
                raise ValueError(f"'{nom}' n'a pas de sortie de vitesse")
            suivis.append((k_s, tw_s, poids))

        for k, (el, cm) in enumerate(zip(els, cms)):
            su = ou[k]
            # entrees de service -> externes (sauf isum des bus)
            for j, s in enumerate(el.svc_inputs):
                nm = f"{el.name}.{s}"
                if nm in uext_idx:
                    Ruext[su + j, uext_idx[nm]] = 1.0
            # w_ref global
            if self.repere == "reference" and tw_ref is not None and el.passif:
                Ry[su + el.wref_idx(), oy[k_ref] + tw_ref] = 1.0
            else:
                Ruext[su + el.wref_idx(), uext_idx[W_REF]] = 1.0
            # theta_g <- angle de l'element de reference
            if k_ref is not None:
                Ry[su + el.thetaref_idx(), oy[k_ref] + ty_ref] = 1.0
            if suivis and k == k_ref:
                j_w = el.svc_inputs.index("w_in")
                for k_s, tw_s, poids in suivis:
                    Ry[su + j_w, oy[k_s] + tw_s] += poids
            # tensions de port <- etats de tension du bus, PAR GENRE : un port
            # alternatif prend les deux sorties vd_o/vq_o du noeud, un port
            # continu son unique v_o.
            n_ac = len(el.ports)
            for pk, b in enumerate(el.buses if not isinstance(el, Bus) else []):
                bi = bus_pos[b]
                if pk < n_ac:
                    ivd, ivq = el.port_input_idx(pk)
                    Ry[su + ivd, oy[bi] + 0] = 1.0  # vd_o du bus
                    Ry[su + ivq, oy[bi] + 1] = 1.0  # vq_o du bus
                else:
                    Ry[su + el.port_dc_input_idx(pk - n_ac), oy[bi] + 0] = 1.0  # v_o
            # isum du bus <- somme des courants de port a ce bus
            if isinstance(el, Bus):
                for m, (e2, cm2) in enumerate(zip(els, cms)):
                    if isinstance(e2, Bus):
                        continue
                    n_ac2 = len(e2.ports)
                    for pk, b2 in enumerate(e2.buses):
                        if b2 != el.bus:
                            continue
                        if pk < n_ac2:
                            iid, iiq = e2.port_output_idx(pk)
                            Ry[su + 0, oy[m] + iid] += 1.0
                            Ry[su + 1, oy[m] + iiq] += 1.0
                        else:
                            ii = e2.port_dc_output_idx(pk - n_ac2)
                            Ry[su + 0, oy[m] + ii] += 1.0

        # --- injections de courant d'essai aux bus demandes ---
        for b in injections:
            if b not in bus_pos:
                raise ValueError(f"injection : bus {b} inconnu")
            su = ou[bus_pos[b]]
            for j, axe in enumerate(("id", "iq")):
                uext_names.append(f"inj{b}.{axe}")
                col = np.zeros((nu, 1))
                col[su + j, 0] = 1.0  # isum_d / isum_q du bus
                Ruext = np.hstack((Ruext, col))

        # --- sorties exposees : tout sauf les courants de port ---
        out_names, rows = [], []
        for k, (el, cm) in enumerate(zip(els, cms)):
            for j, o in enumerate(cm.spec.outputs):
                if j < el.n_port_signaux():
                    continue  # courants de port (AC puis DC) : internes
                out_names.append(f"{el.name}.{o}")
                rows.append(oy[k] + j)

        return Cablage(
            els=els,
            cms=cms,
            ox=ox,
            oz=oz,
            ou=ou,
            oy=oy,
            nx=nx,
            nz=nz,
            nu=nu,
            ny=ny,
            Ruext=Ruext,
            Ry=Ry,
            uext_names=uext_names,
            state_names=state_names,
            rows_sorties=np.array(rows, dtype=int),
            output_names=out_names,
        )

    def jacobiens(self, c: Cablage | None = None) -> dict[str, np.ndarray]:
        """Les jacobiens de TOUS les elements, empiles bloc-diagonalement.

        ``Fx .. Hu`` : chaque element y verse ses neuf blocs locaux au point de
        fonctionnement, aux indices globaux du :class:`Cablage`. C'est
        l'ensemble des equations du reseau, element par element, AVANT toute
        fermeture : le cablage (``Ry``, ``Ruext``) les relie ensuite.
        """
        if any(el.x0 is None for el in self.all_elements):
            self.solve_equilibrium()
        if c is None:
            c = self.cablage()
        nx, nz, nu, ny = c.nx, c.nz, c.nu, c.ny
        J = {
            "Fx": np.zeros((nx, nx)),
            "Fz": np.zeros((nx, nz)),
            "Fu": np.zeros((nx, nu)),
            "Gx": np.zeros((nz, nx)),
            "Gz": np.zeros((nz, nz)),
            "Gu": np.zeros((nz, nu)),
            "Hx": np.zeros((ny, nx)),
            "Hz": np.zeros((ny, nz)),
            "Hu": np.zeros((ny, nu)),
        }
        for k, (el, cm) in enumerate(zip(c.els, c.cms)):
            x, z, u, p = el.x0, el.z0, el.u0, el.p
            lignes = {
                "F": (c.ox[k], cm.n_x),
                "G": (c.oz[k], cm.n_z),
                "H": (c.oy[k], cm.n_y),
            }
            cols = {
                "x": (c.ox[k], cm.n_x),
                "z": (c.oz[k], cm.n_z),
                "u": (c.ou[k], cm.n_u),
            }
            for bloc, (r0, nr) in lignes.items():
                for var, (c0, nc) in cols.items():
                    fn = getattr(cm, bloc.lower() + var)
                    J[bloc + var][r0 : r0 + nr, c0 : c0 + nc] = fn(x, z, u, p)
        return J

    def dae_lineaire(
        self, outputs: list[str] | None = None, injections: tuple[int, ...] = ()
    ) -> "DAELineaire":
        """Le DAE lineaire GLOBAL du reseau, avant elimination de quoi que ce soit.

        Chaque classe d'equipement donne l'ensemble de ses equations -- ``f``,
        ``g``, ``h`` --, le cablage les referme par la LFT, et il reste UN
        systeme pour tout le reseau ::

            dx/dt = AX.x + AZ.z + AU.u         (toutes les dynamiques)
                0 = BX.x + BZ.z + BU.u         (toutes les contraintes)
                y = CX.x + CZ.z + CU.u

        ``z`` rassemble les variables algebriques de tous les elements
        (relations flux/courant des machines, tensions des noeuds ideaux,
        etats rendus statiques par :meth:`Element.reduire`). La resolution est
        globale : :meth:`DAELineaire.etat` les elimine toutes ensemble, et
        reduit l'indice quand une contrainte ne porte que sur des etats.
        """
        c = self.cablage(injections)
        J = self.jacobiens(c)
        Fx, Fz, Fu = J["Fx"], J["Fz"], J["Fu"]
        Gx, Gz, Gu = J["Gx"], J["Gz"], J["Gu"]
        Hx, Hz, Hu = J["Hx"], J["Hz"], J["Hu"]
        Ruext, Ry = c.Ruext, c.Ry
        nu = c.nu

        # --- fermeture LFT : U = Ruext.uext + Ry.Y ---
        # Une factorisation LU, et non cond() puis inv() : sur un reseau de
        # quelques milliers d'etats, la SVD de cond() coutait a elle seule 40 s.
        # Le conditionnement est ESTIME sur la LU (LAPACK gecon, norme 1),
        # et les trois produits par M^-1 sont des resolutions triangulaires.
        # Ry et Ruext sont des matrices de CABLAGE : une poignee de 1 par
        # ligne (0,02 % de termes non nuls sur un reseau de 3000 etats). Leurs
        # produits se font en creux ; le reste reste dense.
        Rs = sparse.csr_matrix(Ry)
        M = -(Rs @ Hu)
        M[np.diag_indices(nu)] += 1.0
        lu, piv = lu_factor(M, check_finite=False, overwrite_a=False)
        rcond, _info = dgecon(lu, np.abs(M).sum(axis=0).max(), norm="1")
        if not rcond > 1e-12:
            raise RuntimeError(
                "boucle algebrique tension/courant : (I - Ry.Hu) singuliere"
            )
        resoudre = lambda B: lu_solve((lu, piv), B, check_finite=False)  # noqa: E731
        Wx, Wz, Wu = resoudre(Rs @ Hx), resoudre(Rs @ Hz), resoudre(Ruext)

        # --- sorties exposees (cf. Cablage), filtrees si on en demande moins ---
        garde = [
            i for i, nm in enumerate(c.output_names) if outputs is None or nm in outputs
        ]
        rows = c.rows_sorties[garde]
        noms_z = [
            f"{el.name}.{n}" for el, cm in zip(c.els, c.cms) for n in cm.spec.algebs
        ]
        return DAELineaire(
            AX=Fx + Fu @ Wx,
            AZ=Fz + Fu @ Wz,
            AU=Fu @ Wu,
            BX=Gx + Gu @ Wx,
            BZ=Gz + Gu @ Wz,
            BU=Gu @ Wu,
            CX=Hx[rows] + Hu[rows] @ Wx,
            CZ=Hz[rows] + Hu[rows] @ Wz,
            CU=Hu[rows] @ Wu,
            etats=list(c.state_names),
            algebriques=noms_z,
            entrees=list(c.uext_names),
            sorties=[c.output_names[i] for i in garde],
        )

    def build_state_space(
        self, outputs: list[str] | None = None, injections: tuple[int, ...] = ()
    ) -> StateSpace:
        """Le modele d'etat lineaire : jacobiens au point de fonctionnement,
        fermeture LFT, puis elimination GLOBALE des variables algebriques
        (:meth:`dae_lineaire`, :meth:`DAELineaire.etat`). ``injections`` :
        voir :meth:`cablage`."""
        return self.dae_lineaire(outputs, injections).etat()


def _rang_plein(M: np.ndarray, tol: float) -> bool:
    """``M`` carree est-elle de rang plein, a l'echelle pres ? Lignes et
    colonnes sont ramenees a un maximum unite (quelques passes de Ruiz) avant
    la SVD : le rang ne depend plus des unites des equations."""
    n = M.shape[0]
    if n == 0:
        return True
    A = np.abs(M)
    r = np.ones(n)
    c = np.ones(n)
    for _ in range(8):
        lr = (A * r[:, None] * c[None, :]).max(axis=1)
        if (lr <= 0).any():
            return False
        r /= np.sqrt(lr)
        lc = (A * r[:, None] * c[None, :]).max(axis=0)
        if (lc <= 0).any():
            return False
        c /= np.sqrt(lc)
    s = np.linalg.svd(M * r[:, None] * c[None, :], compute_uv=False)
    return bool(s.min() > tol * max(float(s.max()), 1.0))


def _pivots_stables(K: np.ndarray, seuil: float = 0.1) -> list[int]:
    """Les etats a eliminer, un par ligne de ``K`` (contraintes ``K.x = 0``).

    Un QR a pivot de colonne, a une difference pres : parmi les colonnes dont
    le residu atteint ``seuil`` fois le plus grand, on prend la PREMIERE dans
    l'ordre des etats, pas la plus grande. Les courants d'une contrainte de
    noeud ideal ont souvent des normes egales a l'arrondi pres, ou qui se
    croisent quand tourne l'angle du repere propre d'un convertisseur : le
    pivot strict basculait alors d'un point de fonctionnement a l'autre.
    Memes valeurs propres, mais pas les memes etats nommes -- ce qui casse le
    suivi des modes d'une continuation. Le seuil garde ``K[:, elim]`` bien
    conditionne.
    """
    R = np.array(K, dtype=float)
    libres = np.ones(R.shape[1], dtype=bool)
    choisis = []
    for _ in range(R.shape[0]):
        normes = np.where(libres, np.linalg.norm(R, axis=0), 0.0)
        if normes.max() <= 1e-12:
            raise RuntimeError("contraintes d'etat dependantes : reduction d'indice")
        j = int(np.flatnonzero(normes >= seuil * normes.max())[0])
        choisis.append(j)
        libres[j] = False
        q = R[:, j] / normes[j]
        R -= np.outer(q, q @ R)
    return choisis


@dataclass
class DAELineaire:
    """Le DAE lineaire global d'un reseau (voir :meth:`Network.dae_lineaire`)."""

    AX: np.ndarray
    AZ: np.ndarray
    AU: np.ndarray
    BX: np.ndarray
    BZ: np.ndarray
    BU: np.ndarray
    CX: np.ndarray
    CZ: np.ndarray
    CU: np.ndarray
    etats: list[str]
    algebriques: list[str]
    entrees: list[str]
    sorties: list[str]

    def etat(self, tol: float = 1e-10, indice_max: int = 3) -> StateSpace:
        """Elimine les variables algebriques, en reduisant l'indice s'il le faut.

        INDICE 1 -- ``BZ`` inversible : complement de Schur, ``z = -BZ^-1
        (BX.x + BU.u)``. C'est le cas de tous les reseaux dont chaque noeud a
        sa capacite ou sa resistance de causalite.

        INDICE 2 -- ``BZ`` singulier : certaines combinaisons ``N^T.g`` des
        contraintes ne dependent d'aucune variable algebrique. Elles
        contraignent les ETATS, ``K.x = 0`` -- deux courants de branche en
        serie a un noeud ideal : ``i1 - i2 = 0``. On les remplace par leur
        derivee, ``K.(AX.x + AZ.z + AU.u) = 0``, qui, elle, determine ``z`` (la
        tension du noeud est celle qui garde les deux courants egaux). La
        dynamique obtenue laisse invariante la variete ``K.x = 0`` ; on y
        restreint le modele en eliminant un etat par contrainte :
        ``x = T.x_r``. Les etats elimines sont nommes dans
        :attr:`StateSpace.elimines`, et ``T`` permet de les reconstruire.

        Une contrainte qui ferait intervenir une ENTREE demanderait sa derivee :
        elle est refusee, comme une contrainte qui ne porte sur rien (variable
        algebrique indeterminee).
        """
        AX, AZ, AU = self.AX, self.AZ, self.AU
        BX, BZ, BU = self.BX, self.BZ, self.BU
        CX, CZ, CU = self.CX, self.CZ, self.CU
        nx, nz = AX.shape[0], BZ.shape[0]
        contraintes: list[np.ndarray] = []
        rapides: list[np.ndarray] = []
        for _ in range(indice_max + 1):
            if nz == 0:
                break
            # Le rang se juge sur BZ EQUILIBRE (lignes et colonnes ramenees a
            # l'unite) : une bretelle de 0,001 pu a cote d'une charge donne
            # des coefficients de 1e7 et de 1, et la SVD brute voyait
            # singulier (1e-11 en relatif) un bloc qui ne l'est pas.
            if _rang_plein(BZ, tol):
                break
            U, s, _vh = np.linalg.svd(BZ)
            echelle = max(float(s.max()), 1.0)
            r = int(np.sum(s > tol * echelle))
            if r == nz:
                break
            Q, N = U[:, :r], U[:, r:]
            Kx, Ku = N.T @ BX, N.T @ BU
            ex = max(float(np.abs(BX).max(initial=0.0)), 1.0)
            if np.abs(Ku).max(initial=0.0) > tol * max(ex, np.abs(BU).max(initial=0.0)):
                raise RuntimeError(
                    "contrainte algebrique sur les ENTREES : sa derivee ferait "
                    "intervenir du/dt. Verifier les noeuds ideaux (IdealBus) : "
                    "un courant injecte a un noeud ideal ou ne se rencontrent "
                    "que des branches inductives y voit une impedance impropre "
                    "(s.L), sans modele d'etat. Donner une capacite au noeud, ou "
                    "passer en NOEUDS = 'resistance' ; l'impedance seule se lit "
                    "sur le DAE (DAELineaire.reponse_frequentielle)."
                )
            vides = np.abs(Kx).max(axis=1) <= tol * ex
            if vides.any():
                raise RuntimeError(
                    "variables algebriques indeterminees (contrainte vide) ; "
                    "suspectes : " + ", ".join(self._suspectes(N[:, vides]))
                )
            contraintes.append(Kx)
            # les directions ou pousse l'algebrique que g ne fixe pas (la
            # tension d'un noeud ideal) : celles du mode parasite de la
            # resistance de causalite quand R tend vers l'infini
            rapides.append(AZ @ _vh[r:].conj().T)
            BX = np.vstack((Q.T @ BX, Kx @ AX))
            BZ = np.vstack((Q.T @ BZ, Kx @ AZ))
            BU = np.vstack((Q.T @ BU, Kx @ AU))
        else:
            raise RuntimeError(
                f"indice du DAE superieur a {indice_max} : variables suspectes "
                + ", ".join(self._suspectes(np.linalg.svd(BZ)[0][:, r:]))
            )

        if nz:
            SX = np.linalg.solve(BZ, BX)
            SU = np.linalg.solve(BZ, BU)
            A, B = AX - AZ @ SX, AU - AZ @ SU
            Cm, Dm = CX - CZ @ SX, CU - CZ @ SU
        else:
            A, B, Cm, Dm = AX, AU, CX, CU

        if not contraintes:
            return StateSpace(
                A=A,
                B=B,
                C=Cm,
                D=Dm,
                states=list(self.etats),
                inputs=list(self.entrees),
                outputs=list(self.sorties),
            )

        # --- restriction a la variete K.x = 0 : un etat de moins par ligne ---
        K = np.vstack(contraintes)
        K = K / np.abs(K).max(axis=1, keepdims=True)
        elim = np.sort(_pivots_stables(K))
        reste = np.array([i for i in range(nx) if i not in set(elim)], dtype=int)
        T = np.zeros((nx, len(reste)))
        T[reste, np.arange(len(reste))] = 1.0
        T[elim] = -np.linalg.solve(K[:, elim], K[:, reste])
        # L'inverse a gauche de T dont le noyau est fait de ces directions
        # rapides : un mode lent du modele a resistance de causalite est
        # orthogonal (a gauche) au mode parasite en -R, et c'est la limite
        # R -> infini de ses participations que rend spectre_complet.
        try:
            T_gauche = np.linalg.inv(np.hstack([T, np.real(np.hstack(rapides))]))
            T_gauche = T_gauche[: len(reste)]
        except np.linalg.LinAlgError:
            T_gauche = None
        return StateSpace(
            A=A[reste] @ T,
            B=B[reste],
            C=Cm @ T,
            D=Dm,
            states=[self.etats[i] for i in reste],
            inputs=list(self.entrees),
            outputs=list(self.sorties),
            elimines=[self.etats[i] for i in elim],
            T=T,
            etats_complets=list(self.etats),
            T_gauche=T_gauche,
        )

    def reponse_frequentielle(
        self, f, entrees: list[str], sorties: list[str]
    ) -> np.ndarray:
        """``H(j.2.pi.f)`` des ``entrees`` vers les ``sorties``, forme
        ``(len(f), len(sorties), len(entrees))``, resolue sur le DAE lui-meme ::

            [s.I - AX, -AZ ; -BX, -BZ] . [x ; z] = [AU ; BU] . u

        sans eliminer ``z``. A l'indice 1, c'est exactement le modele d'etat
        de :meth:`etat`. Mais un courant injecte a un noeud ideal ou ne se
        rencontrent que des courants d'etat y voit ``Z(s) ~ s.L`` en haute
        frequence : une impedance IMPROPRE, sans modele d'etat (:meth:`etat`
        refuse la contrainte sur l'entree), et pourtant bien definie a toute
        frequence finie -- c'est elle que l'on lit ici.
        """
        iu = [self.entrees.index(n) for n in entrees]
        iy = [self.sorties.index(n) for n in sorties]
        nx, nz = self.AX.shape[0], self.BZ.shape[0]
        M = -np.block([[self.AX, self.AZ], [self.BX, self.BZ]])
        diag = np.r_[np.ones(nx), np.zeros(nz)]
        R = np.vstack([self.AU[:, iu], self.BU[:, iu]])
        C = np.hstack([self.CX[iy], self.CZ[iy]])
        D = self.CU[np.ix_(iy, iu)]
        f = np.atleast_1d(np.asarray(f, dtype=float))
        H = np.empty((f.size, len(iy), len(iu)), dtype=complex)
        for k, fk in enumerate(f):
            H[k] = C @ np.linalg.solve(M + np.diag(2j * np.pi * fk * diag), R) + D
        return H

    def _suspectes(self, N: np.ndarray, n: int = 6) -> list[str]:
        """Les equations algebriques qui portent le plus sur les directions ``N``."""
        poids = np.abs(N).sum(axis=1)
        return [self.algebriques[i] for i in np.argsort(-poids)[:n] if poids[i] > 0]
