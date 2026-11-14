"""Sensibilites parametriques : quel parametre deplace un mode, et de combien.

La sensibilite aux coefficients de ``A`` (:func:`g2egrid.analyse.modal.sensibilite`)
designe une EQUATION -- le couple d'etats ``(i, j)`` dont le coefficient
compte --, mais pas le bouton sur lequel agir : un coefficient de ``A``
melange plusieurs parametres, et le complement de Schur qui elimine les
variables algebriques les etale sur toute la matrice. Ce module remonte
jusqu'au PARAMETRE.

PRINCIPE -- on ne derive pas ``A``, on derive les equations. Le reseau est
un DAE lineaire global (:meth:`g2egrid.reseau.network.Network.dae_lineaire`), et sa
forme la plus simple garde les entrees des elements ``U`` comme inconnues ::

        d/dt x = Fx.x + Fz.z + Fu.U
            0  = Gx.x + Gz.z + Gu.U
            0  = Ry.(Hx.x + Hz.z + Hu.U) - U

Un faisceau ``J - lambda.E``, avec ``E = diag(I, 0, 0)``, dont les valeurs
propres finies sont celles de ``A``. Un parametre d'un element ne touche que
les neuf blocs de CET element, et c'est la toute l'astuce : pour une valeur
propre simple, de vecteurs propres a droite ``v`` et a gauche ``w`` du
faisceau ::

        dlambda/dp = w^T . (dJ/dp) . v  /  (w^T . E . v)

``v`` et ``w`` se deduisent de ceux de ``A`` en resolvant les equations
algebriques, une fois pour tout le reseau ; ``dJ/dp`` est la derivee
SYMBOLIQUE du jacobien de l'element (:meth:`CompiledModel.derivees_parametres`),
compilee une fois par classe. La sensibilite a TOUS les parametres de TOUS
les elements coute donc une factorisation, puis un produit scalaire creux par
element.

CE QUI EST TENU FIXE -- le point de fonctionnement. ``dlambda/dp`` mesure
l'effet du parametre sur la DYNAMIQUE autour du meme equilibre : c'est la
bonne question pour un gain de regulation, qui ne deplace pas l'equilibre.
Pour un parametre qui le deplace (une reactance de ligne, une consigne),
``methode="differences"`` avec ``equilibre=True`` refait load flow et
equilibre a chaque perturbation et donne la sensibilite TOTALE. Les
drapeaux des modes discrets (limiteurs, ``Element.modes``) sont aussi tenus
a leur valeur : c'est ce qui garde les derivees symboliques d'un GFM courtes
(voir :meth:`CompiledModel.derivees_parametres`).

La valeur rendue en premier est la sensibilite NORMALISEE ``p.dlambda/dp`` :
le deplacement du mode pour une variation RELATIVE du parametre (1 = 100 %).
Elle rend comparables un gain de 1e3 et une constante de temps de 1e-3.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from g2egrid.analyse import modal

__all__ = ["Sensibilite", "sensibilites", "tableau"]


@dataclass
class Sensibilite:
    """La sensibilite d'un mode ``lambda`` a un parametre ``p`` d'un element."""

    element: str
    parametre: str
    valeur: float  # p, au point de fonctionnement
    mode: complex  # lambda
    dlambda: complex  # dlambda/dp

    @property
    def normalisee(self) -> complex:
        """``p . dlambda/dp`` : deplacement du mode pour ``dp/p = 1``."""
        return self.valeur * self.dlambda

    @property
    def df(self) -> float:
        """``p . df/dp`` [Hz] : effet d'une variation relative sur la frequence."""
        return float(self.normalisee.imag / (2.0 * np.pi))

    @property
    def dzeta(self) -> float:
        """``p . dzeta/dp`` [%] : effet d'une variation relative sur l'amortissement.

        ``zeta = -sigma/|lambda|`` donne ``dzeta = omega.(sigma.domega -
        omega.dsigma)/|lambda|^3``.
        """
        lam, d = self.mode, self.normalisee
        s, w = lam.real, lam.imag
        m = abs(lam)
        if m == 0.0:
            return 0.0
        return float(100.0 * w * (s * d.imag - w * d.real) / m**3)

    def __repr__(self) -> str:
        n = self.normalisee
        return (
            f"<Sensibilite {self.element}.{self.parametre}={self.valeur:.4g} : "
            f"p.dlambda/dp = {n.real:+.4g}{n.imag:+.4g}j>"
        )


def _choisir_mode(ss, mode) -> int:
    """L'indice du mode dans ``eig(A)``, qu'on le designe par un indice, une
    valeur propre (la plus proche est retenue) ou un :class:`modal.Mode`."""
    lam = np.linalg.eigvals(ss.A)
    if isinstance(mode, modal.Mode):
        mode = mode.valeur
    if isinstance(mode, (int, np.integer)):
        if not -len(lam) <= mode < len(lam):
            raise IndexError(f"mode {mode} hors du spectre ({len(lam)} modes)")
        return int(mode) % len(lam)
    return int(np.argmin(np.abs(lam - complex(mode))))


def sensibilites(
    net,
    mode,
    elements: list[str] | None = None,
    parametres: list[str] | None = None,
    methode: str = "analytique",
    equilibre: bool = False,
    pas: float = 1e-6,
    plancher: float = 1e-3,
) -> list[Sensibilite]:
    """Les sensibilites du mode ``mode`` a chaque parametre, les plus fortes d'abord.

    ``mode`` : indice dans ``eig(A)`` (celui de :func:`g2egrid.analyse.modal.modes`),
    valeur propre approchee, ou :class:`g2egrid.analyse.modal.Mode`.
    ``elements`` / ``parametres`` : restreindre aux noms donnes.
    ``methode`` :

    * ``"analytique"`` (defaut) -- derivees symboliques, point de
      fonctionnement fixe. Exacte, et d'un seul coup pour tout le reseau.
    * ``"differences"`` -- differences finies centrees sur le modele
      reconstruit, au pas ``h = pas . max(|p0|, plancher)`` (relatif a la
      valeur du parametre, avec un plancher ABSOLU pour p0 ~ 0). Si ``p0 - h``
      change de signe (p0 < h : un parametre physique positif ne doit pas
      devenir negatif), la difference est unilaterale d'ordre 2. Avec ``equilibre=True``, load flow et equilibre sont
      refaits a chaque perturbation : la sensibilite TOTALE, point de
      fonctionnement compris. Un parametre que l'initialisation recalcule
      (``Kpdc`` auto-regle, impedance d'une charge a puissance constante) y
      sort nul, puisqu'il est ecrase.

    Tri par ``|p.dlambda/dp|`` decroissant.
    """
    ss = net.build_state_space()
    k = _choisir_mode(ss, mode)
    if methode == "analytique":
        if ss.elimines:
            raise NotImplementedError(
                "reseau reduit en indice (noeuds ideaux) : utiliser "
                "methode='differences'"
            )
        out = _analytique(net, ss, k, elements, parametres)
    elif methode == "differences":
        out = _differences(net, ss, k, elements, parametres, equilibre, pas, plancher)
    else:
        raise ValueError(f"methode {methode!r} : 'analytique' ou 'differences'")
    out.sort(key=lambda s: -abs(s.normalisee))
    return out


def _retenu(el, nom, elements, parametres) -> bool:
    return (elements is None or el.name in elements) and (
        parametres is None or nom in parametres
    )


def _analytique(net, ss, k, elements, parametres) -> list[Sensibilite]:
    c = net.cablage()
    J = net.jacobiens(c)
    nz, nu = c.nz, c.nu
    Ry = c.Ry

    lam, Vr, Wl = modal._decomposition(ss)
    lam_k = complex(lam[k])
    vx, wx = Vr[:, k], Wl[k]  # wx . vx = 1 (bi-orthonormalisation)

    # le faisceau etendu : d = x, a = (z, U)
    Jda = np.hstack((J["Fz"], J["Fu"]))
    Jad = np.vstack((J["Gx"], Ry @ J["Hx"]))
    Jaa = np.block(
        [
            [J["Gz"], J["Gu"]],
            [Ry @ J["Hz"], Ry @ J["Hu"] - np.eye(nu)],
        ]
    )
    if nz + nu:
        va = -np.linalg.solve(Jaa, Jad @ vx)
        wa = -np.linalg.solve(Jaa.T, Jda.T @ wx)
    else:
        va = wa = np.zeros(0, dtype=complex)
    vz, vU = va[:nz], va[nz:]
    wz, wU = wa[:nz], wa[nz:]
    # les lignes U du faisceau valent Ry.(dH).v : leur poids a gauche est Ry^T.wU
    r = Ry.T @ wU
    denom = complex(wx @ vx)

    out: list[Sensibilite] = []
    for m, (el, cm) in enumerate(zip(c.els, c.cms)):
        if elements is not None and el.name not in elements:
            continue
        # les drapeaux des modes discrets (limiteurs) sont figes a leur valeur
        params = cm.spec.params
        figes = {d: el.p[params.index(d)] for d in el.modes if d in params}
        ip, li, co, valeurs = cm.derivees_parametres(figes)
        if not len(ip):
            continue
        ox, oz, ou, oy = c.ox[m], c.oz[m], c.ou[m], c.oy[m]
        gauche = np.concatenate(
            (wx[ox : ox + cm.n_x], wz[oz : oz + cm.n_z], r[oy : oy + cm.n_y])
        )
        droite = np.concatenate(
            (vx[ox : ox + cm.n_x], vz[oz : oz + cm.n_z], vU[ou : ou + cm.n_u])
        )
        d = valeurs(el.x0, el.z0, el.u0, el.p)
        dl = np.zeros(cm.n_p, dtype=complex)
        np.add.at(dl, ip, gauche[li] * d * droite[co])
        dl /= denom
        for j, nom in enumerate(cm.spec.params):
            if not _retenu(el, nom, elements, parametres):
                continue
            out.append(Sensibilite(el.name, nom, float(el.p[j]), lam_k, complex(dl[j])))
    return out


def _differences(net, ss, k, elements, parametres, equilibre, pas, plancher=1e-3):
    lam_k = complex(np.linalg.eigvals(ss.A)[k])

    def spectre():
        if equilibre:
            net.run_loadflow()
            net.solve_equilibrium()
        return np.linalg.eigvals(net.build_state_space().A)

    out: list[Sensibilite] = []
    for el in net.all_elements:
        if el.x0 is None:
            continue
        for j, nom in enumerate(el.params):
            if not _retenu(el, nom, elements, parametres):
                continue
            p0 = float(el.p[j])
            h = pas * max(abs(p0), plancher)
            if (p0 > 0.0 and p0 - h <= 0.0) or p0 == 0.0:
                sens, lam_pas = +1.0, (0.0, 1.0, 2.0)  # unilateral vers le haut
            elif p0 < 0.0 and p0 + h >= 0.0:
                sens, lam_pas = -1.0, (0.0, 1.0, 2.0)  # unilateral vers le bas
            else:
                sens, lam_pas = None, (+1.0, -1.0)  # centre
            vals = []
            for m in lam_pas:
                el.set_param(nom, p0 + (m if sens is None else sens * m) * h)
                try:
                    lam = spectre()
                finally:
                    el.set_param(nom, p0)
                vals.append(lam[np.argmin(np.abs(lam - lam_k))])
            if sens is None:
                d = (vals[0] - vals[1]) / (2 * h)
            else:  # (-3 f0 + 4 f1 - f2) / 2h, d'ordre 2, signe selon le sens
                d = sens * (-3 * vals[0] + 4 * vals[1] - vals[2]) / (2 * h)
            out.append(Sensibilite(el.name, nom, p0, lam_k, complex(d)))
    if equilibre:  # rendre au reseau son point de fonctionnement d'origine
        net.run_loadflow()
        net.solve_equilibrium()
    return out


def tableau(sens: list[Sensibilite], n: int = 10) -> str:
    """Les ``n`` premieres sensibilites, en texte aligne."""
    lignes = [
        f"{'element.parametre':<28}{'valeur':>12}{'p.dRe':>12}{'p.df [Hz]':>12}"
        f"{'p.dzeta [%]':>13}"
    ]
    for s in sens[:n]:
        lignes.append(
            f"{s.element + '.' + s.parametre:<28}{s.valeur:>12.4g}"
            f"{s.normalisee.real:>12.4g}{s.df:>12.4g}{s.dzeta:>13.4g}"
        )
    return "\n".join(lignes)
