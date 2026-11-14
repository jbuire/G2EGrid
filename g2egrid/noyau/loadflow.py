"""Load flow Newton-Raphson en coordonnees polaires.

Convention des types de bus (identique au framework MATLAB) ::

    1 = slack (V et theta imposes)
    2 = PV    (P et |V| imposes)
    3 = PQ    (P et Q imposes)

Le Y-bus est construit a partir des branches serie R-L avec rapport de
prise ``n`` cote ``from`` ::

    Y[f,f] += y/n^2   Y[f,t] -= y/n   Y[t,f] -= y/n   Y[t,t] += y

RESEAU MIXTE -- :class:`ProblemeACDC` et :func:`run_loadflow_acdc` ajoutent le
cote continu au MEME Newton. :func:`run_loadflow` n'en sait rien : un reseau
sans noeud continu n'emprunte que lui. Depuis FR-027 il est CREUX (Y-bus CSC,
jacobien de ``dSbus_dV``, SuperLU) ; le Newton mixte, lui, reste dense.

CE QUI EST GARANTI N'EST PLUS LE TEXTE MAIS LE RESULTAT. Le lot 3 pinait le
contenu de ce fichier ; le lot 4 doit y ecrire, donc la garantie a change de
nature. ``tests/test_loadflow_acdc.py`` la tient par deux filets
independants, sur chaque reseau du catalogue :

* une REFERENCE stockee. ``reference_loadflow.json`` (reseaux alternatifs)
  vient d'un export ``git archive`` du code d'avant le load flow mixte ;
  ``reference_loadflow_acdc.json`` (reseaux continus) a ete produite par
  ``developpement/regenerer_reference_loadflow.py``. Le filet exige
  ``max|V - V_ref| <= 1e-11`` pu (et de meme pour ``Vdc``), et le meme
  nombre d'iterations pour descendre sous 1e-10 ;
* les EQUATIONS, recalculees par le test depuis le Y-bus et la solution :
  P et Q a chaque bus, bilan de chaque noeud continu et de chaque station,
  sous 1e-10 pu, sans aucune reference.

CE QUI N'EST PLUS GARANTI : l'identite a l'octet. Elle a ete exigee jusqu'a
FR-001, et ne tenait que sur la machine qui avait produit la reference. Sur
une autre, l'arrondi de BLAS/LAPACK et des sin/cos vectorises de NumPy
deplace V de 1e-17 a 1,8e-13 pu, au meme nombre d'iterations. Une expression
deplacee qui ne change que l'arrondi (1 ulp) ne se voit donc plus. Une
consigne deplacee de 1e-6 se voit toujours, et un jacobien faux de 1 % aussi,
par les iterations.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

__all__ = [
    "DC_P",
    "DC_SLACK",
    "MODES_VSC",
    "LoadFlowACDCResult",
    "LoadFlowResult",
    "ProblemeACDC",
    "VSCSpec",
    "build_ybus",
    "run_loadflow",
    "run_loadflow_acdc",
]

SLACK, PV, PQ = 1, 2, 3

#: Types de noeud CONTINU. Volontairement numerotes a part de SLACK/PV/PQ : un
#: noeud continu n'a ni angle ni reactif, donc aucun des trois types
#: alternatifs ne s'y transpose. 1 = tension imposee, 2 = puissance imposee.
DC_SLACK, DC_P = 1, 2

#: Les modes d'un convertisseur : la grandeur continue tenue, puis
#: l'alternative. Regle reprise de VeraGrid (idee seulement, aucun code
#: copie) : un VSC tient EXACTEMENT deux grandeurs, ni une ni trois.
#:
#: Deux modes de plus pour les reseaux continus maillés (CIGRE B4) :
#:
#: * ``Vtheta`` -- le convertisseur FORME son ile alternative : il en tient le
#:   module et l'angle, il en est le slack. Ses deux puissances sont libres et
#:   sortent du bilan de l'ile ; vu du continu, c'est un terminal a puissance
#:   imposee par l'alternatif. C'est la station offshore d'un parc eolien ;
#: * ``droop`` -- statisme tension continue / puissance :
#:   ``p_ac = p_ref + s.(v_dc - v_dc,ref)``. Le convertisseur ne tient ni sa
#:   puissance ni sa tension continue, il les lie ; ``s`` (pu systeme par pu
#:   de tension) est la pente.
#: * ``droopI`` -- statisme COURANT continu / tension continue :
#:   ``i_dc = i_ref + s.(v_dc - v_dc,ref)``, d'ou
#:   ``p_ac = v_dc.(i_ref + s.(v_dc - v_dc,ref))`` avec ``i_ref = p_ref/v_dc,ref``.
#:
#: Cote alternatif, le convertisseur tient sa tension (bus PV), son reactif
#: (bus PQ) ou forme son ile (``Vtheta``, slack) -- c'est le type du bus qui le
#: dit -- ou suit un statisme Q-V (``kqv`` > 0 sur un bus PQ) :
#: ``q_ac = q_ref + kqv.(v_ref - |V|)``.
MODES_VSC = ("PQ", "PV", "VdcQ", "VdcV", "Vtheta", "droop", "droopI")


@dataclass
class LoadFlowResult:
    V: np.ndarray  # tensions complexes (pu, repere global)
    S: np.ndarray  # puissances nettes injectees par bus (pu)
    converged: bool
    iterations: int
    mismatch: float

    def vmag(self) -> np.ndarray:
        return np.abs(self.V)

    def theta(self) -> np.ndarray:
        return np.angle(self.V)


def build_ybus(n_bus: int, branches, shunts=()) -> sp.csc_array:
    """Y-bus CREUX (CSC), l'unique constructeur du Y-bus alternatif.

    ``branches`` : liste de ``(f, t, R, X, n)`` (indices 0-based, ``n`` le
    rapport de transformation cote ``f``) ; ``shunts`` : liste de ``(i, y)``,
    admittances a la terre du bus ``i`` (capacite de bus, resistance de
    causalite, reactance shunt). Les doublons (lignes en parallele, plusieurs
    shunts) se somment a la conversion.

    ``Network.build_ybus`` s'en sert, et en rend la forme dense par defaut ;
    le Newton de :func:`run_loadflow` travaille, lui, sur le creux.
    """
    br = np.asarray(list(branches), dtype=float).reshape(-1, 5)
    f = br[:, 0].astype(int)
    t = br[:, 1].astype(int)
    z = br[:, 2] + 1j * br[:, 3]
    if np.any(z == 0):
        raise ValueError("branche d'impedance serie nulle dans le Y-bus")
    n = br[:, 4]
    y = 1.0 / z
    lignes = [f, f, t, t]
    colonnes = [f, t, f, t]
    valeurs = [y / (n * n), -y / n, -y / n, y]
    sh = list(shunts)
    if sh:
        lignes.append(np.array([i for i, _ in sh], dtype=int))
        colonnes.append(lignes[-1])
        valeurs.append(np.array([v for _, v in sh], dtype=complex))
    return sp.coo_array(
        (
            np.concatenate(valeurs).astype(complex),
            (np.concatenate(lignes), np.concatenate(colonnes)),
        ),
        shape=(n_bus, n_bus),
    ).tocsc()


def jacobien_polaire(
    Y, Vc: np.ndarray, idx_th: np.ndarray, idx_v: np.ndarray
) -> sp.csc_array:
    """Jacobien polaire creux ``J = [[dP/dth, dP/dV], [dQ/dth, dQ/dV]]``.

    Les formules de ``dSbus_dV`` de MATPOWER, pour ``S = V.conj(Y.V)`` ::

        dS/dV_m[i,j] = V_i.conj(Y_ij.V_j/|V_j|) + [i=j].conj(I_i).V_i/|V_i|
        dS/dth [i,j] = j.V_i.conj([i=j].I_i - Y_ij.V_j)        (I = Y.V)

    ecrites directement sur les triplets (ligne, colonne, valeur) de ``Y`` :
    aucun produit creux-creux, donc un cout en O(nnz) et pas de surcout fixe
    notable sur les petits reseaux (les produits de matrices diagonales de
    SciPy coutaient 8 ms par load flow sur ieee14).

    Lignes : ``P`` des bus ``idx_th`` puis ``Q`` des bus ``idx_v`` ; colonnes :
    angles de ``idx_th`` puis modules de ``idx_v``. ``J.dx = mismatch`` avec
    ``mismatch = consigne - injection`` (meme sens que la boucle d'origine,
    gardee comme oracle dans ``tests/test_loadflow_acdc.py``).
    """
    n = Vc.size
    Yc = sp.coo_array(Y)
    Ib = Yc @ Vc
    Vn = Vc / np.abs(Vc)
    ligne, c0, y = Yc.row, Yc.col, Yc.data
    d = np.arange(n)
    lig = np.concatenate((ligne, d))
    col = np.concatenate((c0, d))
    # Les termes en Y_ij, puis ceux de la diagonale (sommes au meme triplet).
    dVm = np.concatenate((Vc[ligne] * np.conj(y * Vn[c0]), np.conj(Ib) * Vn))
    dVa = np.concatenate((1j * Vc[ligne] * np.conj(-y * Vc[c0]), 1j * Vc * np.conj(Ib)))
    nt, nv = len(idx_th), len(idx_v)
    pos_th = np.full(n, -1)
    pos_th[idx_th] = np.arange(nt)
    pos_v = np.full(n, -1)
    pos_v[idx_v] = np.arange(nv)
    r_th, r_v, c_th, c_v = pos_th[lig], pos_v[lig], pos_th[col], pos_v[col]
    blocs = (
        (r_th, c_th, dVa.real, 0, 0),  # dP/dth
        (r_th, c_v, dVm.real, 0, nt),  # dP/dV
        (r_v, c_th, dVa.imag, nt, 0),  # dQ/dth
        (r_v, c_v, dVm.imag, nt, nt),  # dQ/dV
    )
    L, C, D = [], [], []
    for pr, pc, val, ro, co in blocs:
        m = (pr >= 0) & (pc >= 0)
        L.append(pr[m] + ro)
        C.append(pc[m] + co)
        D.append(val[m])
    return sp.coo_array(
        (np.concatenate(D), (np.concatenate(L), np.concatenate(C))),
        shape=(nt + nv, nt + nv),
    ).tocsc()


def jacobien_polaire_dense(
    Y: np.ndarray,
    V: np.ndarray,
    th: np.ndarray,
    P: np.ndarray,
    Q: np.ndarray,
    idx_th: np.ndarray,
    idx_v: np.ndarray,
) -> np.ndarray:
    """Le meme jacobien, dense, vectorise (FR-027 etape 1) : voie des petits reseaux.

    Les quatre blocs sont formes sur tout le reseau d'un coup, puis on n'en
    garde que les lignes et colonnes des inconnues. Hors diagonale, les termes
    sont ceux de la double boucle d'origine, dans le meme ordre d'operations ;
    la diagonale est ensuite ecrasee par sa forme en P, Q. Cette voie est
    CELLE D'AVANT FR-027 : sous ``SEUIL_CREUX`` bus elle rend les memes
    tensions qu'avant (a 1e-15 pu pres, l'ecart venant de l'ordre de somme du
    Y-bus), et ``ProblemeACDC.jacobien``, qui en est la copie terme a terme,
    reste comparable a elle.
    """
    G, B = Y.real, Y.imag
    cos = np.cos(th[:, None] - th[None, :])
    sin = np.sin(th[:, None] - th[None, :])
    a = G * sin - B * cos
    b = G * cos + B * sin
    VV = V[:, None] * V[None, :]
    H = VV * a
    N = V[:, None] * b
    M = -VV * b
    L = V[:, None] * a
    d = np.arange(len(V))
    Gd, Bd = G[d, d], B[d, d]
    H[d, d] = -Q - Bd * V**2
    N[d, d] = P / V + Gd * V
    M[d, d] = P - Gd * V**2
    L[d, d] = Q / V - Bd * V
    return np.block(
        [
            [H[np.ix_(idx_th, idx_th)], N[np.ix_(idx_th, idx_v)]],
            [M[np.ix_(idx_v, idx_th)], L[np.ix_(idx_v, idx_v)]],
        ]
    )


#: Nombre de bus a partir duquel ``run_loadflow`` passe au Newton creux.
#: Mesure (grilles m x m de ``developpement/bench_loadflow.py``, mediane de 30) : le dense gagne
#: de 0,5 a 1 ms sur 3 a 40 bus (le creux paie ~1 ms de construction de
#: matrices SciPy), les deux se croisent vers 100 bus, et le creux est
#: ensuite 3 fois (146 bus) a 15 fois (423 bus) plus rapide, puis 200 fois
#: sur 4 900 bus (le dense y coute 74 s). Sous le seuil, la voie dense
#: garde l'arithmetique du Newton mixte (``ProblemeACDC``, dense), ce qu'un
#: test du continu (``test_residu_equilibre_liaison_mixte``, residu amplifie
#: 3e5 fois, egalite a 1e-12 relatif entre load flow alternatif et mixte)
#: exige : avec le creux il echoue (dx 3,2206e-8 contre 3,2200e-8).
SEUIL_CREUX = 100


def run_loadflow(
    Y,
    btype: np.ndarray,
    Vsp: np.ndarray,
    theta_sp: np.ndarray,
    Psp: np.ndarray,
    Qsp: np.ndarray,
    tol: float = 1e-10,
    max_iter: int = 40,
    amorti: bool = False,
    polissage: bool = False,
    methode: str = "auto",
) -> LoadFlowResult:
    """Newton-Raphson polaire. ``Psp``/``Qsp`` sont les injections NETTES (pu).

    ``methode`` : ``"creux"`` (Y-bus CSC, jacobien de
    :func:`jacobien_polaire`, factorisation SuperLU), ``"dense"`` (jacobien de
    :func:`jacobien_polaire_dense`, ``np.linalg.solve``) ou ``"auto"`` (creux
    a partir de ``SEUIL_CREUX`` bus). ``Y`` peut etre dense ou creux dans les
    deux cas. Les deux voies resolvent les memes equations : tensions a 1e-12
    pu et memes iterations sur le catalogue (``tests/test_loadflow_creux.py``).

    ``iterations`` compte les evaluations du mismatch jusqu'a ``mismatch <
    tol`` (convention d'origine). ``converged`` est juge sur le mismatch du
    point RENDU, y compris quand le dernier pas de la boucle tombe sur la
    solution (avant FR-027 il restait faux : mismatch perime).

    ``polissage`` : apres convergence, un pas de Newton de plus ramene le
    mismatch du point rendu de la tolerance d'arret (jusqu'a 3,2e-11 sur le
    catalogue) vers le plancher d'arrondi, sans compter une iteration ;
    ``mismatch`` est celui du point poli. DESACTIVE par defaut : le pas
    deplace V de 1,7e-11 pu (wscc9_cigre), au-dela de la tolerance du filet de
    non-regression (1e-11 pu). Le plancher, lui, vaut 1e-14 a 1e-12 selon le
    reseau (1,2e-12 sur ieee39_mixte) : le polissage ne donne pas 1e-15.

    ``amorti`` : si le pas plein n'a pas reduit la norme 2 du mismatch, il est
    divise par deux (jusqu'a 8 fois) et l'on garde le meilleur essai. Sans
    effet tant que le pas plein reduit le mismatch, mais desactive par defaut :
    un Newton non amorti peut remonter avant de converger (wscc3sm_*), et
    l'amortir changerait alors le compte des iterations.
    """
    if methode not in ("auto", "creux", "dense"):
        raise ValueError(f"methode inconnue : {methode!r}")
    n = len(btype)
    creux = methode == "creux" or (methode == "auto" and n >= SEUIL_CREUX)
    if creux:
        Yop = sp.csc_array(Y)
    else:
        Yop = Y.toarray() if sp.issparse(Y) else np.asarray(Y)

    V = np.array(Vsp, dtype=float).copy()
    th = np.array(theta_sp, dtype=float).copy()

    pv = np.where(btype == PV)[0]
    pq = np.where(btype == PQ)[0]
    idx_th = np.concatenate((pv, pq))  # inconnues d'angle
    idx_v = pq  # inconnues de module
    nt = len(idx_th)

    def evaluer():
        Vc = V * np.exp(1j * th)
        S = Vc * np.conj(Yop @ Vc)
        mis = np.concatenate((Psp[idx_th] - S.real[idx_th], Qsp[idx_v] - S.imag[idx_v]))
        return Vc, S, mis

    def pas(Vc, S, mis):
        if creux:
            J = jacobien_polaire(Yop, Vc, idx_th, idx_v)
            try:
                return splu(J).solve(mis)
            except RuntimeError as e:  # SuperLU : "Factor is exactly singular"
                raise np.linalg.LinAlgError(f"jacobien singulier : {e}") from e
        J = jacobien_polaire_dense(Yop, V, th, S.real, S.imag, idx_th, idx_v)
        return np.linalg.solve(J, mis)

    def norme(mis):
        return float(np.max(np.abs(mis))) if mis.size else 0.0

    def avancer(alpha, dx, V0, th0):
        th[idx_th] = th0[idx_th] + alpha * dx[:nt]
        V[idx_v] = V0[idx_v] + alpha * dx[nt:]
        return evaluer()

    Vc, S, mis = evaluer()
    mismatch = norme(mis)
    pas_faits = 0
    while mismatch >= tol and pas_faits < max_iter:
        dx = pas(Vc, S, mis)
        V0, th0 = V.copy(), th.copy()
        Vc, S, mis_new = avancer(1.0, dx, V0, th0)
        if amorti:
            n0 = np.linalg.norm(mis)
            essais = [(np.linalg.norm(mis_new), 1.0)]
            alpha = 1.0
            while essais[-1][0] >= n0 and len(essais) <= 8:
                alpha *= 0.5
                Vc, S, mis_new = avancer(alpha, dx, V0, th0)
                essais.append((np.linalg.norm(mis_new), alpha))
            if essais[-1][0] >= n0:  # aucun essai ne descend : le moins mauvais
                alpha = min(essais)[1]
                Vc, S, mis_new = avancer(alpha, dx, V0, th0)
        mis = mis_new
        mismatch = norme(mis)
        pas_faits += 1
    converged = mismatch < tol
    it = min(pas_faits + 1, max_iter)

    if converged and polissage and mismatch > 0.0:
        dx = pas(Vc, S, mis)
        th[idx_th] += dx[:nt]
        V[idx_v] += dx[nt:]
        Vc, S, mis = evaluer()
        mismatch = norme(mis)

    return LoadFlowResult(
        V=Vc, S=S, converged=converged, iterations=it, mismatch=mismatch
    )


# ==========================================================================
# Load flow AC/DC unifie
# ==========================================================================
@dataclass
class VSCSpec:
    """Un convertisseur, vu du load flow : deux noeuds et deux grandeurs tenues.

    C'est la SEULE chose que ce module sait d'un convertisseur. Il n'importe
    aucun element, ne connait ni ``md`` ni ``isd`` : il recoit ``p_dc``, une
    fonction, et s'en sert comme d'une boite noire. Le couplage physique est
    monte par :meth:`g2egrid.reseau.network.Network._couplage_vsc`, qui, lui, a acces
    au modele dynamique.

    ``p_dc(vm, p_ac, q_ac, v_dc)`` rend la puissance INJECTEE au noeud continu
    (pu base systeme), au point alternatif ``(vm, p_ac, q_ac)`` et a la tension
    continue ``v_dc``. Elle est negative pour un convertisseur qui tire du
    continu pour alimenter l'alternatif.

    LES QUATRE ARGUMENTS SONT LES SEULS DONT ELLE DEPEND, et c'est ce qui rend
    le bloc de couplage du jacobien si maigre : le courant de port se deduit du
    MODULE de la tension alternative et de la puissance complexe, pas de leur
    phase -- tourner le repere ne change pas un courant continu. Mesure :
    ``i_dc(theta_g + 0.7) - i_dc`` vaut exactement zero.

    ``p_ac``/``q_ac`` sont les puissances PROGRAMMEES du convertisseur lui-meme,
    deja nettes de ce que les autres elements de son bus consomment. En mode
    ``Vdc...``, ``p_ac`` n'est que le point de depart : le Newton lui ajoute
    l'inconnue additive. ``q_autres`` sert au cas inverse -- sur un bus PV le
    reactif du convertisseur n'est pas programme, il vaut ``Q_bus - q_autres``.
    """

    iac: int
    idc: int
    mode: str
    p_dc: Callable[[float, float, float, float], float]
    p_ac: float = 0.0
    q_ac: float = 0.0
    q_autres: float = 0.0
    nom: str = "VSC"
    #: ce que les autres elements du bus posent en actif : sur le slack
    #: (``Vtheta``) le P du convertisseur est ``P_bus - p_autres``
    p_autres: float = 0.0
    #: statisme (``droop``) : pente ``s`` en pu systeme par pu de tension
    #: continue, et tension de consigne ; ``p_ac`` est alors ``p_ref``
    pente: float = 0.0
    vdc_ref: float = 1.0
    #: statisme Q-V (bus PQ) : pente ``kqv`` en pu systeme de reactif par pu de
    #: tension alternative, et consigne ``v_ref`` ; ``q_ac`` est ``q_ref``
    kqv: float = 0.0
    v_ref_ac: float = 1.0

    def __post_init__(self) -> None:
        if self.mode not in MODES_VSC:
            raise ValueError(
                f"{self.nom}: mode '{self.mode}' inconnu "
                f"(attendu : {', '.join(MODES_VSC)})"
            )

    @property
    def tient_vdc(self) -> bool:
        """Le convertisseur tient-il la tension continue de son noeud ?"""
        return self.mode.startswith("Vdc")

    @property
    def forme(self) -> bool:
        """Le convertisseur forme-t-il son ile alternative (slack) ?"""
        return self.mode == "Vtheta"

    @property
    def statisme(self) -> bool:
        """Le convertisseur suit-il un statisme sur la tension continue
        (puissance ``droop`` ou courant ``droopI``) ?"""
        return self.mode in ("droop", "droopI")

    def ecart_statisme(self, vdc: float) -> tuple[float, float]:
        """``(p_ac - p_ref, dp_ac/dv_dc)`` du statisme a la tension ``vdc``."""
        dv = vdc - self.vdc_ref
        if self.mode == "droopI":
            i_ref = self.p_ac / self.vdc_ref
            return vdc * (i_ref + self.pente * dv) - self.p_ac, i_ref + self.pente * (
                vdc + dv
            )
        return self.pente * dv, self.pente


@dataclass
class LoadFlowACDCResult(LoadFlowResult):
    """Le resultat du Newton unifie : le cote alternatif, plus le continu.

    Herite de :class:`LoadFlowResult` sans rien lui retirer, de sorte que tout
    ce qui lit deja ``lf.V``, ``lf.S`` ou ``lf.converged`` -- l'equilibre par
    element, les notebooks, l'interface -- continue de le lire sans savoir
    qu'il y a du continu derriere.
    """

    #: tensions continues (pu), dans l'ordre des noeuds continus tries
    Vdc: np.ndarray = field(default_factory=lambda: np.zeros(0))
    #: puissance que les convertisseurs injectent a chaque noeud continu (pu)
    Pdc: np.ndarray = field(default_factory=lambda: np.zeros(0))
    #: correction du P programme des VSC en mode ``Vdc...`` (pu), dans l'ordre
    #: de la liste des VSC en mode ``Vdc...``
    P_vsc: np.ndarray = field(default_factory=lambda: np.zeros(0))
    #: pourquoi le Newton s'est arrete avant d'avoir converge ("" sinon)
    message: str = ""


class ProblemeACDC:
    """Le systeme non lineaire du load flow mixte : inconnues, residu, jacobien.

    POURQUOI UNE CLASSE, ET PAS UNE FONCTION -- le jacobien de ce systeme est
    la piece qu'il faut pouvoir eprouver SEULE. Un jacobien faux ne fait pas
    diverger un Newton : il le fait converger plus lentement. Sur les reseaux
    simples on ne voit rien, et l'erreur ressort sur le cas qui compte, six
    heures plus tard. D'ou :meth:`residu` et :meth:`jacobien` exposees, et un
    test qui les confronte aux differences finies avant meme que le Newton
    existe.

    LES INCONNUES, dans cet ordre ::

        [ theta (bus PV et PQ) | V (bus PQ) | V_dc (noeuds de type 2)
          | p (VSC en mode Vdc...) ]

    ``p`` est la correction ADDITIVE du P programme du convertisseur. Le bus
    alternatif d'un VSC garde le type que lui donne ``Y_DER`` -- PQ pour un
    grid-following, PV pour un grid-forming -- et c'est ce qui garantit la non
    regression : le bloc alternatif/alternatif du jacobien est celui d'avant,
    aux memes lignes et aux memes colonnes, plus une colonne ``-1``.

    LES EQUATIONS, dans le meme ordre ::

        P programme - P(V, theta)        aux bus PV et PQ
        Q programme - Q(V, theta)        aux bus PQ
        R_j                              aux noeuds continus de type 2 et a
                                         ceux que tient un VSC en mode Vdc...

    avec ``R_j = Pinj_j + somme_c P_dc(c) - V_j . somme_k G_jk (V_j - V_k)``.
    Un noeud continu de type 1 alimente par une source n'apporte NI inconnue NI
    equation : c'est le slack du reseau continu, et la source y encaisse le
    reliquat, exactement comme le bus infini du cote alternatif. Un noeud de
    type 1 tenu par un VSC, lui, garde son equation : sa tension est imposee,
    mais le bilan du noeud reste a satisfaire, et c'est lui qui fixe le P du
    convertisseur.

    LE COUPLAGE est tout entier dans ``P_dc(c)``, et il vient du MODELE
    DYNAMIQUE : ``P_dc = -V_dc . (Gdc.V_dc + md.isd + mq.isq)``, ou ``md``,
    ``mq``, ``isd`` et ``isq`` sortent des formules d'initialisation du
    convertisseur au point alternatif courant. Les pertes valent donc
    ``R1.|is|^2 + Gdc.V_dc^2`` -- ce que le modele fait vraiment, jamais un
    rendement pose a la main. C'est ce qui fait du point de fonctionnement
    mixte un POINT FIXE du modele non lineaire, et non un point voisin ou le
    reseau linearise n'est deja plus a l'equilibre.

    CONVENTION DU JACOBIEN -- :meth:`jacobien` rend ``J = -dF/dx``, soit la
    derivee des grandeurs CALCULEES et non celle du residu. C'est la convention
    de :func:`run_loadflow` (``J.dx = mismatch``), et c'est elle qui permet d'y
    recopier ``J1..J4`` sans changer un signe.
    """

    def __init__(
        self,
        Y: np.ndarray,
        btype: np.ndarray,
        Vsp: np.ndarray,
        theta_sp: np.ndarray,
        Psp: np.ndarray,
        Qsp: np.ndarray,
        Gdc: np.ndarray,
        dctype: np.ndarray,
        Vdc_sp: np.ndarray,
        Pdc_sp: np.ndarray,
        vscs: list[VSCSpec],
        p_bus_dc: list | None = None,
        pas: float = 1e-6,
    ):
        self.Y = np.asarray(Y, dtype=complex)
        self.btype = np.asarray(btype, dtype=int)
        self.Vsp = np.asarray(Vsp, dtype=float)
        self.theta_sp = np.asarray(theta_sp, dtype=float)
        self.Psp = np.asarray(Psp, dtype=float)
        self.Qsp = np.asarray(Qsp, dtype=float)
        self.Gdc = np.asarray(Gdc, dtype=float)
        self.dctype = np.asarray(dctype, dtype=int)
        self.Vdc_sp = np.asarray(Vdc_sp, dtype=float)
        self.Pdc_sp = np.asarray(Pdc_sp, dtype=float)
        self.vscs = list(vscs)
        n_dc = len(self.dctype)
        self.p_bus_dc = list(p_bus_dc) if p_bus_dc is not None else [None] * n_dc
        self.pas = float(pas)

        self.idx_th = np.concatenate(
            (np.where(self.btype == PV)[0], np.where(self.btype == PQ)[0])
        ).astype(int)
        self.idx_v = np.where(self.btype == PQ)[0].astype(int)
        self.idx_dcv = np.where(self.dctype == DC_P)[0].astype(int)
        #: les VSC qui portent une inconnue de puissance
        self.idx_p = [k for k, c in enumerate(self.vscs) if c.tient_vdc]
        tenus = {self.vscs[k].idc for k in self.idx_p}
        self.idx_dceq = np.array(sorted(set(self.idx_dcv.tolist()) | tenus), dtype=int)

        # position d'un bus dans chaque liste d'inconnues, pour ne pas l'y
        # rechercher a chaque ligne du jacobien
        self._pos_th = {int(b): a for a, b in enumerate(self.idx_th)}
        self._pos_v = {int(b): a for a, b in enumerate(self.idx_v)}
        self._pos_dcv = {int(b): a for a, b in enumerate(self.idx_dcv)}
        self._pos_p = {k: a for a, k in enumerate(self.idx_p)}

        nt, nv = len(self.idx_th), len(self.idx_v)
        self.n_x = nt + nv + len(self.idx_dcv) + len(self.idx_p)
        self.n_f = nt + nv + len(self.idx_dceq)
        if self.n_x != self.n_f:
            raise RuntimeError(
                f"load flow AC/DC mal pose : {self.n_x} inconnue(s) pour "
                f"{self.n_f} equation(s). Un noeud continu de type 1 doit etre "
                "tenu par un VSC en mode Vdc... ou par une source continue, et "
                "un seul VSC peut tenir un noeud donne."
            )

    # ------------------------------------------------------------------
    def depart(self) -> np.ndarray:
        """Le point de depart : les consignes du cote alternatif, les amorces du
        cote continu, et zero pour les corrections de puissance."""
        return np.concatenate(
            (
                self.theta_sp[self.idx_th],
                self.Vsp[self.idx_v],
                self.Vdc_sp[self.idx_dcv],
                np.zeros(len(self.idx_p)),
            )
        )

    def decoder(self, x: np.ndarray):
        """``x`` -> ``(V, theta, V_dc, p)``, chaque vecteur COMPLET.

        Les composantes qui ne sont pas des inconnues gardent leur consigne :
        c'est ce qui permet d'ecrire les formules du jacobien pour un bus
        quelconque, qu'il soit PV, PQ ou slack.
        """
        nt, nv, ndv = len(self.idx_th), len(self.idx_v), len(self.idx_dcv)
        th = self.theta_sp.copy()
        V = self.Vsp.copy()
        vdc = self.Vdc_sp.copy()
        th[self.idx_th] = x[:nt]
        V[self.idx_v] = x[nt : nt + nv]
        vdc[self.idx_dcv] = x[nt + nv : nt + nv + ndv]
        return V, th, vdc, x[nt + nv + ndv :]

    def injections(self, V: np.ndarray, th: np.ndarray):
        """``(P, Q)`` injectes a chaque bus alternatif -- comme run_loadflow."""
        Vc = V * np.exp(1j * th)
        S = Vc * np.conj(self.Y @ Vc)
        return S.real, S.imag

    def point_vsc(
        self, k: int, V, vdc, p, Q, P=None
    ) -> tuple[float, float, float, float]:
        """Les quatre scalaires dont depend le ``p_dc`` du VSC ``k``.

        ``q_ac`` est le seul qui demande un choix. Sur un bus PQ il est
        PROGRAMME, et le prendre tel quel evite de faire dependre le bloc
        continu des angles ; sur un bus PV il est libre, et il faut le lire
        dans le load flow -- d'ou le terme de chaine ``dP_dc/dq . dQ/dx`` du
        jacobien, sans lequel les derivees de ce bloc seraient fausses des que
        le convertisseur est un grid-forming.
        """
        c = self.vscs[k]
        p_ac = c.p_ac + (p[self._pos_p[k]] if c.tient_vdc else 0.0)
        if c.forme:
            # le slack de son ile : sa puissance est ce que l'ile appelle
            p_ac = P[c.iac] - c.p_autres
        elif c.statisme:
            p_ac = c.p_ac + c.ecart_statisme(float(vdc[c.idc]))[0]
        if self.btype[c.iac] == PQ:
            q_ac = c.q_ac + c.kqv * (c.v_ref_ac - V[c.iac])
        else:
            q_ac = Q[c.iac] - c.q_autres
        return float(V[c.iac]), float(p_ac), float(q_ac), float(vdc[c.idc])

    def injections_dc(self, V, vdc, p, Q, P=None) -> np.ndarray:
        """Puissance injectee a chaque noeud continu, convertisseurs compris."""
        inj = self.Pdc_sp.astype(float).copy()
        for j, f in enumerate(self.p_bus_dc):
            if f is not None:
                inj[j] += float(f(float(vdc[j])))
        for k, c in enumerate(self.vscs):
            inj[c.idc] += float(c.p_dc(*self.point_vsc(k, V, vdc, p, Q, P)))
        return inj

    # ------------------------------------------------------------------
    def residu(self, x: np.ndarray) -> np.ndarray:
        """``F(x)``, dans l'ordre des equations annonce par la docstring."""
        V, th, vdc, p = self.decoder(x)
        P, Q = self.injections(V, th)
        Psp = self.Psp.copy()
        for k in self.idx_p:
            Psp[self.vscs[k].iac] += p[self._pos_p[k]]
        Qsp = self.Qsp.copy()
        for c in self.vscs:
            if c.statisme:
                Psp[c.iac] += c.ecart_statisme(float(vdc[c.idc]))[0]
            if c.kqv and self.btype[c.iac] == PQ:
                Qsp[c.iac] += c.kqv * (c.v_ref_ac - V[c.iac])
        inj = self.injections_dc(V, vdc, p, Q, P)
        flux = vdc * (self.Gdc @ vdc)
        return np.concatenate(
            (
                Psp[self.idx_th] - P[self.idx_th],
                Qsp[self.idx_v] - Q[self.idx_v],
                (inj - flux)[self.idx_dceq],
            )
        )

    # ------------------------------------------------------------------
    def _derivees_vsc(self, k: int, pt) -> dict[str, float]:
        """Les une a trois derivees du ``p_dc`` d'un VSC, par difference centree.

        NUMERIQUES ET ASSUMEES : ``p_dc`` est l'aboutissement des formules
        d'initialisation d'un convertisseur a treize etats, et en ecrire la
        derivee a la main serait un SECOND modele, a tenir en accord avec le
        premier a chaque retouche. Le bloc alternatif/alternatif, lui, reste
        analytique : c'est lui qui porte la convergence, et c'est lui qui ne
        doit pas bouger.

        Une seule derivee par grandeur REELLEMENT inconnue : sur un bus PV le
        module de tension est impose (pas de ``dvm``), sur un bus PQ le reactif
        est programme (pas de ``dq``), un VSC qui ne tient pas la tension
        continue n'a pas d'inconnue de puissance (pas de ``dp``), et un noeud
        continu de type 1 a sa tension imposee (pas de ``dvdc``). D'ou une a
        trois derivees, jamais quatre.
        """
        c = self.vscs[k]
        besoins = [0] if self.btype[c.iac] == PQ else [2]
        if self.btype[c.iac] == PQ and c.kqv:
            besoins.append(2)  # le reactif suit la tension (statisme Q-V)
        if c.tient_vdc or c.forme or c.statisme:
            besoins.append(1)
        if self.dctype[c.idc] == DC_P:
            besoins.append(3)
        noms = ("dvm", "dp", "dq", "dvdc")
        base = list(pt)
        out: dict[str, float] = {}
        for i in besoins:
            h = self.pas * max(1.0, abs(base[i]))
            plus, moins = list(base), list(base)
            plus[i] += h
            moins[i] -= h
            out[noms[i]] = (c.p_dc(*plus) - c.p_dc(*moins)) / (2.0 * h)
        return out

    def jacobien(self, x: np.ndarray) -> np.ndarray:
        """``J = -dF/dx``, au sens de :func:`run_loadflow` : ``J.dx = mismatch``."""
        V, th, vdc, p = self.decoder(x)
        P, Q = self.injections(V, th)
        nt, nv = len(self.idx_th), len(self.idx_v)
        ndv = len(self.idx_dcv)
        J = np.zeros((self.n_f, self.n_x))

        # --- blocs alternatif/alternatif -------------------------------
        # LES FORMULES DE run_loadflow, terme pour terme, mais ecrites pour
        # TOUS les couples (i, j) et non pour les seules lignes retenues. La
        # generalisation ne change pas une expression ; elle donne en plus les
        # lignes dQ_i/dx des bus PV, dont le bloc de couplage a besoin pour sa
        # regle de chaine. Copie plutot que factorisation avec run_loadflow :
        # run_loadflow reste intouchee, et cela vaut mieux qu'une ligne
        # economisee.
        G, B = self.Y.real, self.Y.imag
        cos = np.cos(th[:, None] - th[None, :])
        sin = np.sin(th[:, None] - th[None, :])
        VV = V[:, None] * V[None, :]
        A = G * sin - B * cos
        M = G * cos + B * sin
        Vinv = np.where(V != 0.0, 1.0 / np.where(V != 0.0, V, 1.0), 0.0)
        dPdth, dPdV = VV * A, V[:, None] * M
        dQdth, dQdV = -VV * M, V[:, None] * A
        np.fill_diagonal(dPdth, -Q - np.diag(B) * V**2)
        np.fill_diagonal(dPdV, P * Vinv + np.diag(G) * V)
        np.fill_diagonal(dQdth, P - np.diag(G) * V**2)
        np.fill_diagonal(dQdV, Q * Vinv - np.diag(B) * V)

        if nt:
            J[:nt, :nt] = dPdth[np.ix_(self.idx_th, self.idx_th)]
            J[nt : nt + nv, :nt] = dQdth[np.ix_(self.idx_v, self.idx_th)]
        if nv:
            J[:nt, nt : nt + nv] = dPdV[np.ix_(self.idx_th, self.idx_v)]
            J[nt : nt + nv, nt : nt + nv] = dQdV[np.ix_(self.idx_v, self.idx_v)]

        # --- la colonne -E : F_P = Psp + p - P, donc -dF/dp = -1 -------
        for k in self.idx_p:
            iac = self.vscs[k].iac
            J[self._pos_th[iac], nt + nv + ndv + self._pos_p[k]] = -1.0

        # --- lignes continues : F_j = inj_j - V_j.(Gdc.V)_j ------------
        flux = self.Gdc @ vdc
        dvsc = [
            self._derivees_vsc(k, self.point_vsc(k, V, vdc, p, Q, P))
            for k in range(len(self.vscs))
        ]
        # statisme : la ligne P du bus alternatif depend de v_dc,
        # F_P = Psp + s.(v_dc - v_ref) - P, donc -dF/dv_dc = -s
        pentes = {
            k: c.ecart_statisme(float(vdc[c.idc]))[1]
            for k, c in enumerate(self.vscs)
            if c.statisme
        }
        for k, c in enumerate(self.vscs):
            if c.statisme and c.idc in self._pos_dcv and c.iac in self._pos_th:
                J[self._pos_th[c.iac], nt + nv + self._pos_dcv[c.idc]] -= pentes[k]
            # statisme Q-V : F_Q = Qsp + kqv.(v_ref - V) - Q, -dF/dV = kqv + dQ/dV
            if c.kqv and c.iac in self._pos_v:
                a = nt + self._pos_v[c.iac]
                J[a, a] += c.kqv
        dbus = []
        for j, f in enumerate(self.p_bus_dc):
            if f is None:
                dbus.append(0.0)
                continue
            v = float(vdc[j])
            h = self.pas * max(1.0, abs(v))
            dbus.append((f(v + h) - f(v - h)) / (2.0 * h))

        for r, j in enumerate(self.idx_dceq):
            ligne = nt + nv + r
            # cables, puis injection propre au noeud (source, charge)
            for kk, k in enumerate(self.idx_dcv):
                d = (flux[j] if k == j else 0.0) + vdc[j] * self.Gdc[j, k]
                if k == j:
                    d -= dbus[j]
                J[ligne, nt + nv + kk] += d
            # convertisseurs raccordes a ce noeud
            for k, c in enumerate(self.vscs):
                if c.idc != j:
                    continue
                d = dvsc[k]
                if "dvdc" in d:
                    J[ligne, nt + nv + self._pos_dcv[c.idc]] -= d["dvdc"]
                if "dp" in d and c.tient_vdc:
                    J[ligne, nt + nv + ndv + self._pos_p[k]] -= d["dp"]
                elif "dp" in d and c.statisme and c.idc in self._pos_dcv:
                    J[ligne, nt + nv + self._pos_dcv[c.idc]] -= d["dp"] * pentes[k]
                elif "dp" in d and c.forme:
                    # regle de chaine : le P du formeur est celui de son bus
                    if nt:
                        J[ligne, :nt] -= d["dp"] * dPdth[c.iac, self.idx_th]
                    if nv:
                        J[ligne, nt : nt + nv] -= d["dp"] * dPdV[c.iac, self.idx_v]
                if "dvm" in d:
                    J[ligne, nt + self._pos_v[c.iac]] -= d["dvm"]
                if "dq" in d and self.btype[c.iac] == PQ:
                    # statisme Q-V : q_ac = q_ref + kqv.(v_ref - V)
                    if c.iac in self._pos_v:
                        J[ligne, nt + self._pos_v[c.iac]] += d["dq"] * c.kqv
                elif "dq" in d and nt:
                    # regle de chaine : sur un bus PV le reactif du
                    # convertisseur n'est pas programme, il suit le load flow
                    J[ligne, :nt] -= d["dq"] * dQdth[c.iac, self.idx_th]
                    if nv:
                        J[ligne, nt : nt + nv] -= d["dq"] * dQdV[c.iac, self.idx_v]
        return J

    # ------------------------------------------------------------------
    def resoudre(self, tol: float = 1e-12, max_iter: int = 40) -> LoadFlowACDCResult:
        """Le Newton lui-meme : ``J.dx = F`` jusqu'a ce que ``F`` soit nul.

        Porte par la classe et non par :func:`run_loadflow_acdc` pour que le
        probleme puisse etre construit sans etre resolu -- c'est ainsi que le
        test du jacobien l'attrape, sur le VRAI reseau et non sur un cas
        fabrique pour lui.
        """
        x = self.depart()
        mismatch = np.inf
        it = 0
        message = ""
        for it in range(1, max_iter + 1):
            F = self.residu(x)
            mismatch = float(np.max(np.abs(F))) if F.size else 0.0
            if mismatch < tol:
                break
            try:
                x = x + np.linalg.solve(self.jacobien(x), F)
            except np.linalg.LinAlgError:
                message = (
                    "jacobien singulier a l'iteration "
                    f"{it} : le systeme mixte n'a pas de solution isolee "
                    "(noeud continu sans tenant, ile sans reference de "
                    "tension, ou convertisseurs a puissance imposee tous "
                    "ensemble)"
                )
                break

        V, th, vdc, p = self.decoder(x)
        Vc = V * np.exp(1j * th)
        S = Vc * np.conj(self.Y @ Vc)
        Pdc = np.zeros(len(self.dctype))
        for k, c in enumerate(self.vscs):
            Pdc[c.idc] += float(c.p_dc(*self.point_vsc(k, V, vdc, p, S.imag, S.real)))
        return LoadFlowACDCResult(
            V=Vc,
            S=S,
            converged=mismatch < tol,
            iterations=it,
            mismatch=mismatch,
            Vdc=vdc,
            Pdc=Pdc,
            P_vsc=p,
            message=message,
        )


def run_loadflow_acdc(
    Y: np.ndarray,
    btype: np.ndarray,
    Vsp: np.ndarray,
    theta_sp: np.ndarray,
    Psp: np.ndarray,
    Qsp: np.ndarray,
    Gdc: np.ndarray,
    dctype: np.ndarray,
    Vdc_sp: np.ndarray,
    Pdc_sp: np.ndarray,
    vscs: list[VSCSpec],
    p_bus_dc=None,
    tol: float = 1e-12,
    max_iter: int = 40,
    pas: float = 1e-6,
) -> LoadFlowACDCResult:
    """Newton unifie sur les inconnues alternatives ET continues.

    Un seul systeme, un seul jacobien, une seule convergence. Enchainer deux
    Newton ferait converger chacun autour d'un point que l'autre a deja quitte,
    et les pertes du convertisseur -- qui dependent des deux cotes a la fois --
    ne boucleraient jamais : c'est la lecon des 0,08 pu de la version Kelada.

    POURQUOI ``tol`` EST PLUS SERRE QU'AU LOAD FLOW ALTERNATIF (1e-12 contre
    1e-10) -- un residu de puissance n'est pas un residu de courant. La loi des
    noeuds du modele dynamique s'ecrit ``dv/dt = (wb/C).isum`` : sur un noeud de
    capacite 0,16 pu, un reliquat de 1e-10 ressort a 2e-7 dans ``dx/dt``, et
    ``residu_equilibre`` -- qui mesure precisement si le point de
    linearisation est un point fixe -- le verrait. Une iteration de plus coute
    quelques millisecondes et fait tomber le residu a 5e-15 sur ``hvdc2t`` ;
    le mismatch de :func:`run_loadflow`, lui, ne change pas, ce chemin ne le
    traverse pas.

    LIMITE CONNUE (FR-001) -- 1e-12 est au niveau du plancher d'arrondi du
    residu des que le reseau a des branches raides. Une tension ne se
    represente qu'a un ulp pres (2,2e-16 pres de 1 pu), et un ulp sur ``V_j``
    deplace le bilan du noeud de ``G_jk.ulp``. Mesure en poursuivant le Newton
    apres convergence : 3e-13 a 1,3e-12 sur ``cigre_b4`` et ``cigre_b4_mmc``
    (abs(G_dc) = 4,3e3 pu), 2e-13 a 6e-13 sur ``cigre_mv_dc`` (abs(Y) = 1,9e3 pu).
    Deux consequences :

    * sur ces reseaux, le nombre d'iterations depend de l'arrondi, donc de la
      machine (6 a 8 sous un bruit de 4 ulp), mais pas la solution (4e-14) ;
    * au-dela, le Newton atteint le point fixe et n'en sort plus. ``hvdc2t``
      avec un cable de 30 m (abs(G) = 2,7e5 pu) reste a 2,3e-11 des la 5e
      iteration et finit en « non converge ». Voir
      ``tests/test_loadflow_acdc.py::test_cable_continu_court_converge``
      (``xfail`` strict).
    """
    return ProblemeACDC(
        Y,
        btype,
        Vsp,
        theta_sp,
        Psp,
        Qsp,
        Gdc,
        dctype,
        Vdc_sp,
        Pdc_sp,
        vscs,
        p_bus_dc=p_bus_dc,
        pas=pas,
    ).resoudre(tol=tol, max_iter=max_iter)
