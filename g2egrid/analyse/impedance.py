"""Stabilite par les impedances : un equipement contre le reseau qui l'entoure.

Les etudes de reseaux riches en convertisseurs parlent d'impedances : chaque
convertisseur a son admittance dq ``Y_e(s)`` (:mod:`g2egrid.analyse.frequence`), le
reste du reseau son impedance ``Z_r(s)`` vue du bus de raccordement, et la
stabilite de leur association se juge sur le gain de boucle

    L(s) = Z_r(s) . Y_e(s)      source de COURANT (GFL) contre le reseau
    L(s) = Z_e(s) . Y_r(s)      source de TENSION (GFM, SM) contre le reseau

par le critere de Nyquist generalise. Les deux formes ont les memes zeros --
ceux de ``Y_r + Y_e``, les modes du reseau complet -- mais pas les memes
poles : on prend celle dont la boucle ouverte est stable. Une source de
courant est stable a vide de tension imposee (``Y_e`` sans pole instable), une
source de tension l'est en court-circuit... de courant impose (``Z_e``) ; le
reste du reseau, lui, est d'ordinaire stable sous tension imposee (``Y_r``) et
le reste branche sur une source de courant aussi. Avec ``P`` poles instables
en boucle ouverte, le nombre de modes instables du reseau complet est
``Z = P - N``, ``N`` le nombre de tours de ``det(I + L(jw))`` autour de
l'origine dans le sens trigonometrique quand ``w`` parcourt l'axe imaginaire.
Comme le modele d'etat donne ``Z`` directement, :func:`boucle_mineure` en
deduit ``P`` : le critere ne s'applique tel quel (stable si et seulement si
``N = 0``) que lorsque ``P = 0``, ce qu'il dit.

D'OU VIENNENT LES DEUX IMPEDANCES. Tout vient du modele d'etat GLOBAL, au meme
point de fonctionnement :

* ``Z_t(s)`` : l'impedance de tout le reseau, equipement compris, vue du bus --
  on injecte un courant d'essai dans la loi des noeuds du bus
  (:meth:`Network.build_state_space` avec ``injections``) et on lit sa
  tension ;
* ``Y_e(s)`` : l'admittance d'entree de l'equipement seul,
  :func:`g2egrid.analyse.frequence.admittance_dq` ;
* le reste du reseau est en parallele de l'equipement au bus, donc
  ``Z_r = (Z_t^-1 - Y_e)^-1``.

Ce partage est exact tant que l'equipement n'est couple au reste QUE par sa
tension de bornes. Il ne l'est pas pour l'element de reference d'angle (il
donne ``theta_g`` a tous) : :func:`boucle_mineure` le refuse.

Les poles de ``Z_t`` sont les modes du reseau complet : le critere de Nyquist
sur ``L`` doit donc dire instable exactement quand le modele d'etat l'est --
c'est ce que verifient les tests.

On y ajoute la PASSIVITE d'un equipement : ``Y_e`` est passive a une frequence
quand sa partie hermitienne ``(Y + Y^H)/2`` est definie positive ; une bande ou
elle ne l'est pas est une bande ou l'equipement peut fournir de l'energie a une
oscillation, la ou une resonance du reseau devient dangereuse.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.linalg as sl

from g2egrid.analyse.frequence import admittance_dq

__all__ = [
    "impedance_bus",
    "BoucleMineure",
    "boucle_mineure",
    "nyquist_generalise",
    "passivite",
    "bandes",
]


def _bus_element(net, bus: int):
    from g2egrid.elements.bus import Bus

    for el in net.all_elements:
        if isinstance(el, Bus) and el.bus == bus:
            return el
    raise ValueError(f"bus {bus} introuvable")


def impedance_bus(net, bus: int, f: np.ndarray) -> np.ndarray:
    """L'impedance dq de TOUT le reseau vue du bus ``bus`` : tension du bus
    sur courant injecte, forme ``(len(f), 2, 2)``, en pu base systeme.

    Lue sur le DAE lineaire (:meth:`DAELineaire.reponse_frequentielle`) : a
    un noeud ideal ou ne se rencontrent que des branches inductives,
    l'impedance croit comme ``s.L`` et n'a pas de modele d'etat."""
    if any(el.x0 is None for el in net.all_elements):
        net.solve_equilibrium()
    b = _bus_element(net, bus)
    sorties = [f"{b.name}.vd_o", f"{b.name}.vq_o"]
    dae = net.dae_lineaire(outputs=sorties, injections=(bus,))
    return dae.reponse_frequentielle(f, [f"inj{bus}.id", f"inj{bus}.iq"], sorties)


def nyquist_generalise(L: np.ndarray, f: np.ndarray) -> dict:
    """Le critere de Nyquist generalise sur le gain de boucle ``L(j.2.pi.f)``
    (forme ``(len(f), 2, 2)``, ``f`` croissantes et positives).

    Les systemes dq ont des coefficients reels : ``L(-jw)`` est le conjugue de
    ``L(jw)``. Le chemin parcourt donc les frequences negatives (conjugues,
    dans l'ordre inverse), puis les positives, et se referme en ligne droite
    -- au-dela de ``f[-1]`` le gain de boucle doit etre negligeable. Rend le
    nombre de tours de l'origine (``tours``, sens trigonometrique), la plus petite
    distance des lieux caracteristiques (valeurs propres de ``L``) au point
    ``-1``, et la frequence ou elle est atteinte."""
    d = np.linalg.det(np.eye(2)[None] + L)
    chemin = np.concatenate([np.conj(d[::-1]), d])
    arg = np.unwrap(np.angle(chemin))
    fermeture = np.angle(chemin[0] / chemin[-1])
    tours = (arg[-1] - arg[0] + fermeture) / (2 * np.pi)
    lieux = np.linalg.eigvals(L)
    dist = np.abs(lieux + 1).min(axis=1)
    k = int(np.argmin(dist))
    return {
        "tours": int(np.round(tours)),
        "tours_brut": float(tours),
        "marge": float(dist[k]),
        "f_marge": float(f[k]),
        "lieux": lieux,
        "det": d,
    }


def passivite(Y: np.ndarray) -> np.ndarray:
    """La plus petite valeur propre de la partie hermitienne de chaque
    matrice de ``Y`` : negative, l'equipement n'est pas passif a cette
    frequence."""
    H = 0.5 * (Y + np.conj(np.swapaxes(Y, -1, -2)))
    return np.linalg.eigvalsh(H)[:, 0]


def bandes(f: np.ndarray, masque: np.ndarray) -> list[tuple[float, float]]:
    """Les intervalles de frequence ou ``masque`` est vrai."""
    out, debut, prec = [], None, None
    for fk, m in zip(f, masque):
        if m and debut is None:
            debut = fk
        if not m and debut is not None:
            out.append((float(debut), float(prec)))
            debut = None
        prec = fk
    if debut is not None:
        out.append((float(debut), float(f[-1])))
    return out


@dataclass
class BoucleMineure:
    """Un equipement contre le reste du reseau, frequence par frequence."""

    element: str
    bus: int
    f: np.ndarray
    Y_e: np.ndarray  # admittance d'entree de l'equipement
    Z_r: np.ndarray  # impedance du reste du reseau vue du bus
    Z_t: np.ndarray  # impedance de tout le reseau vue du bus
    L: np.ndarray  # gain de boucle (forme ``forme``)
    forme: str  # "courant" : Z_r.Y_e ; "tension" : Z_e.Y_r
    instables: int  # modes instables du reseau complet (modele d'etat), Z
    poles_equipement: int  # poles instables de Y_e (equipement sous tension imposee)
    nyquist: dict = field(default_factory=dict)
    non_passif: list[tuple[float, float]] = field(default_factory=list)

    @property
    def poles_ouverts(self) -> int:
        """``P = Z + N`` : les poles instables de la boucle ouverte que le
        nombre de tours et le modele d'etat impliquent."""
        return self.instables + self.nyquist["tours"]

    @property
    def coherent(self) -> bool:
        """``P >= 0`` : un nombre de poles ne peut pas etre negatif. Faux, le
        comptage des tours est a reprendre (grille, ou gain de boucle non
        negligeable au-dela de la derniere frequence)."""
        return self.poles_ouverts >= 0

    @property
    def applicable(self) -> bool:
        """Boucle ouverte stable : le critere se lit tel quel."""
        return self.poles_ouverts == 0

    @property
    def stable(self) -> bool:
        """Verdict du critere (s'il est applicable) : aucun tour."""
        return self.nyquist["tours"] == 0

    def resume(self) -> str:
        n = self.nyquist
        return (
            f"{self.element} au bus {self.bus} ({self.forme}) : {n['tours']} tour(s), "
            f"{self.instables} mode(s) instable(s) au complet"
            + (
                ""
                if self.applicable
                else f", boucle ouverte instable (P={self.poles_ouverts})"
            )
            + f" ; marge {n['marge']:.3g} a {n['f_marge']:.1f} Hz ; non passif sur "
            + (
                ", ".join(f"{a:.1f}-{b:.1f} Hz" for a, b in self.non_passif)
                or "aucune bande"
            )
        )


def _frequences_critiques(net, el, bus: int, port: int) -> np.ndarray:
    """Les valeurs propres ou det(I + L) a ses poles et ses zeros, connues
    par le modele d'etat : modes du reseau complet, modes de l'equipement,
    zeros de transmission de ``Z_t`` et modes du RESTE du reseau.

    Ces derniers sont les zeros de det(Y_t - Y_e), c'est-a-dire de
    det(I - Z_t.Y_e) : les modes du bouclage de ``Z_t`` par ``-Y_e``, dont la
    matrice d'etat s'ecrit directement. Un mode du reste peut etre amorti a
    quelques 1/s a 16 kHz (smib_gfl) : une largeur de 0,3 Hz, qu'aucune grille
    reguliere ne voit, et qui suffit a fausser le nombre de tours."""
    sortie = [
        f"{_bus_element(net, bus).name}.vd_o",
        f"{_bus_element(net, bus).name}.vq_o",
    ]
    ss = net.build_state_space(outputs=sortie, injections=(bus,))
    iu = [ss.inputs.index(f"inj{bus}.id"), ss.inputs.index(f"inj{bus}.iq")]
    iy = [ss.outputs.index(n) for n in sortie]
    A, B, C, D = ss.A, ss.B[:, iu], ss.C[iy], ss.D[np.ix_(iy, iu)]
    Ae, Be, Ce, De = el.compiled().linearize(el.x0, el.z0, el.u0, el.p)[:4]
    vu, io = list(el.port_input_idx(port)), list(el.port_output_idx(port))
    # Y_e = -(courant de port / tension de port), cf. admittance_dq
    Bv, Ci, Dvi = Be[:, vu], -Ce[io], -De[np.ix_(io, vu)]
    out = [np.linalg.eigvals(A)]
    if Ae.size:
        out.append(np.linalg.eigvals(Ae))
        A_reste = np.block([[A + B @ Dvi @ C, B @ Ci], [Bv @ C, Ae]])
    else:
        A_reste = A + B @ Dvi @ C
    out.append(np.linalg.eigvals(A_reste))
    n = A.shape[0]
    M = np.block([[A, B], [C, D]])
    N = np.zeros_like(M)
    N[:n, :n] = np.eye(n)
    z = sl.eigvals(M, N)
    out.append(z[np.isfinite(z)])
    return np.concatenate(out)


def _autour(lam: np.ndarray, f_min: float, f_max: float) -> np.ndarray:
    """Des frequences resserrees autour de chaque ``Im lam / 2.pi`` de
    ``[f_min, f_max]``, a l'echelle de sa largeur ``|Re lam| / 2.pi``."""
    ecarts = np.array([0, 0.1, 0.3, 0.6, 1, 2, 4, 8, 16, 40, 100])
    ecarts = np.concatenate([-ecarts[:0:-1], ecarts])
    out = []
    for x in lam:
        f0 = abs(x.imag) / (2 * np.pi)
        if not f_min < f0 < f_max:
            continue
        w = max(abs(x.real) / (2 * np.pi), 1e-7 * f0)
        out.append(f0 + w * ecarts)
    if not out:
        return np.zeros(0)
    f = np.concatenate(out)
    return f[(f > f_min) & (f < f_max)]


def boucle_mineure(
    net,
    element: str,
    f: np.ndarray | None = None,
    port: int = 0,
    forme: str | None = None,
) -> BoucleMineure:
    """Le critere de la boucle mineure pour l'equipement ``element`` (nom)
    raccorde a son bus : gain de boucle, Nyquist generalise et bandes de
    non-passivite de l'equipement. ``f`` : frequences [Hz], croissantes
    (defaut : 0,01 Hz a 20 kHz, logarithmique) ; la grille est RAFFINEE la ou
    la phase de det(I + L) saute, et PROLONGEE vers le haut tant que le gain
    de boucle n'est pas negligeable -- ``BoucleMineure.f`` est la grille
    effectivement employee. ``forme`` : ``"courant"``
    (``Z_r.Y_e``) ou ``"tension"`` (``Z_e.Y_r``) ; par defaut ``"courant"``
    pour un GFL, ``"tension"`` pour les autres."""
    if f is None:
        f = np.geomspace(0.01, 20000.0, 2400)
    f = np.asarray(f, dtype=float)
    if any(el.x0 is None for el in net.all_elements):
        net.solve_equilibrium()
    el = next((e for e in net.all_elements if e.name == element), None)
    if el is None:
        raise ValueError(f"equipement {element!r} introuvable")
    if el.name == getattr(net, "reference", None):
        raise ValueError(
            f"{element} est la reference d'angle : il est couple a tout le reseau "
            "par theta_g, pas seulement par sa tension de bornes"
        )
    bus = el.buses[port]
    if forme is None:
        forme = "courant" if el.name.startswith("GFL") else "tension"
    if forme not in ("courant", "tension"):
        raise ValueError(f"forme {forme!r} : 'courant' ou 'tension'")

    def gain(fs):
        Y_e = admittance_dq(el, fs, port)
        Z_t = impedance_bus(net, bus, fs)
        Y_r = np.linalg.inv(Z_t) - Y_e
        Z_r = np.linalg.inv(Y_r)
        L = Z_r @ Y_e if forme == "courant" else np.linalg.inv(Y_e) @ Y_r
        return Y_e, Z_t, Z_r, L

    # GRILLE RAFFINEE LA OU LA PHASE SAUTE. Pres d'un mode peu amorti du
    # reseau complet, det(I + L) frole l'origine et sa phase tourne de pi en
    # une fraction de hertz : entre deux echantillons trop espaces, le
    # deroulage de phase ne sait plus dans quel sens elle a tourne. On insere
    # donc des frequences (moyenne geometrique) partout ou la phase varie de
    # plus de 0,3 rad d'un point au suivant.
    #
    # ET PROLONGEE VERS LE HAUT TANT QUE LE GAIN DE BOUCLE COMPTE. Le contour
    # se referme en ligne droite au-dela de la derniere frequence, ce qui
    # suppose L negligeable. Ce n'est pas toujours le cas a 20 kHz : sur
    # smib_gfl, une resonance du filtre fait encore deux tours entre 20 et
    # 100 kHz, et s'arreter a 20 kHz comptait -2 tours (P = -2, impossible)
    # pour un reseau stable. On ajoute donc des decades tant que la plus grande
    # valeur singuliere de L depasse 0,05 (jusqu'a 100 MHz).
    def raffiner(f, Y_e, Z_t, Z_r, L):
        for _ in range(30):
            d = np.linalg.det(np.eye(2)[None] + L)
            saut = np.abs(np.angle(d[1:] / d[:-1])) > 0.3
            if not saut.any() or f.size > 50000:
                break
            neuves = np.sqrt(f[:-1][saut] * f[1:][saut])
            ajout = gain(neuves)
            ordre = np.argsort(np.concatenate([f, neuves]))
            f = np.concatenate([f, neuves])[ordre]
            Y_e, Z_t, Z_r, L = (
                np.concatenate([x, y])[ordre] for x, y in zip((Y_e, Z_t, Z_r, L), ajout)
            )
        return f, Y_e, Z_t, Z_r, L

    # les poles et zeros connus de det(I + L) d'abord, puis le raffinement
    crit = _autour(_frequences_critiques(net, el, bus, port), f[0], f[-1])
    f = np.unique(np.concatenate([f, crit]))
    f, Y_e, Z_t, Z_r, L = raffiner(f, *gain(f))
    while np.linalg.svd(L[-1], compute_uv=False).max() > 0.05 and f[-1] < 1e8:
        haut = np.geomspace(f[-1], 10 * f[-1], 61)[1:]
        haut = np.unique(
            np.concatenate(
                [
                    haut,
                    _autour(_frequences_critiques(net, el, bus, port), f[-1], haut[-1]),
                ]
            )
        )
        ajout = gain(haut)
        f = np.concatenate([f, haut])
        Y_e, Z_t, Z_r, L = (
            np.concatenate([x, y]) for x, y in zip((Y_e, Z_t, Z_r, L), ajout)
        )
        f, Y_e, Z_t, Z_r, L = raffiner(f, Y_e, Z_t, Z_r, L)
    cm = el.compiled()
    A_e = cm.linearize(el.x0, el.z0, el.u0, el.p)[0]
    lam_e = np.linalg.eigvals(A_e) if A_e.size else np.zeros(0)
    lam = np.linalg.eigvals(net.build_state_space().A)
    instables = int(np.sum(lam.real > 1e-6 * np.maximum(1.0, np.abs(lam))))
    return BoucleMineure(
        element=element,
        bus=bus,
        f=f,
        Y_e=Y_e,
        Z_r=Z_r,
        Z_t=Z_t,
        L=L,
        forme=forme,
        instables=instables,
        poles_equipement=int(
            np.sum(lam_e.real > 1e-6 * np.maximum(1.0, np.abs(lam_e)))
        ),
        nyquist=nyquist_generalise(L, f),
        non_passif=bandes(f, passivite(Y_e) < 0),
    )
