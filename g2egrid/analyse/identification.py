"""Identification modale sur une trajectoire simulee, et appariement au spectre.

Le modele d'etat de :meth:`Network.build_state_space` et le modele non
lineaire de :mod:`g2egrid.temporel` decrivent le MEME reseau, avec le meme
cablage et les memes equations, par deux chemins qui ne partagent pas une
ligne de code de resolution. Le contrat de validation en trois epreuves du
plan dit comment s'en servir :

1. le point de fonctionnement est un point fixe du modele non lineaire
   (``ModeleNonLineaire.residu_equilibre`` ~ 0) ;
2. l'ecart entre reponse lineaire et reponse non lineaire croit comme le
   CARRE de l'amplitude de la perturbation ;
3. les modes lus DIRECTEMENT sur la trajectoire non lineaire sont ceux de la
   matrice ``A``.

Ce module est la troisieme epreuve, et c'est la seule qui ait un sens
physique : les deux premieres disent que le modele lineaire est coherent avec
lui-meme au voisinage de son point de linearisation ; celle-ci dit que ``A``
decrit vraiment ce que fait le systeme non lineaire quand on le perturbe.
C'est aussi, mot pour mot, ce qu'on fait sur un enregistrement de reseau
reel -- d'ou le fait que :func:`identifier` ne demande rien d'autre qu'un
tableau de temps et un tableau de mesures.

DEUX ESTIMATEURS, VOLONTAIREMENT DISTINCTS
    :func:`identifier` porte *matrix pencil* (Hua & Sarkar) et *ERA* (Juang &
    Pappa). Ils partagent la preparation du signal (fenetre, offset,
    verification du pas) et la post-analyse (residus, energie, residu de
    reconstruction), mais pas l'algebre : le pencil exploite l'invariance par
    decalage dans les COLONNES de la matrice de Hankel, via ses vecteurs
    singuliers a droite ; ERA l'exploite dans les LIGNES, via une seconde
    matrice de Hankel decalee. Ils n'ont donc pas les memes defauts -- le
    pencil est plus sobre et plus precis a rang exact, ERA supporte mieux un
    rang surestime -- et pouvoir croiser les deux est la moitie de l'interet
    de l'exercice. Un mode que les deux methodes rendent a la meme place
    existe ; un mode que l'une seule voit est un artefact de troncature.

CE QU'IL FAUT SAVOIR AVANT DE S'EN SERVIR
    * **le pas doit etre uniforme** -- les deux methodes estiment un
      operateur de decalage d'un pas, qui n'existe que si le pas est le meme
      partout. Radau rend des pas qui varient de plusieurs decades. Voir
      :func:`identifier` et :func:`reechantillonner` ;
    * **l'offset doit partir** -- apres un echelon le signal ne revient pas a
      zero, il tend vers un AUTRE equilibre ;
    * **la fenetre ne doit pas commencer a l'evenement** -- le modele porte
      des modes a 1e6-1e7 rad/s qui sont eteints en quelques microsecondes
      mais qui, vus a 5 ms de pas, se replient n'importe ou.

    Ces trois points sont traites la ou ils se posent, avec le detail du
    pourquoi.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import linear_sum_assignment

from g2egrid.analyse.modal import Mode, amortissement, frequence

__all__ = [
    "Appariement",
    "Identification",
    "ModeIdentifie",
    "Paire",
    "RETARD_DEFAUT",
    "apparier",
    "identifier",
    "reconstruire",
    "reechantillonner",
    "sur_trajectoire",
]

#: Retard par defaut entre le debut du signal fourni et le debut de la fenetre
#: d'identification [s]. Le signal est cense commencer a l'evenement ; 5 ms,
#: c'est 300 fois la constante de temps du mode a 1e5 rad/s le plus lent des
#: modes "rapides" d'un reseau dq (resistances de causalite, troncons de
#: cable courts), et moins de 1 % de la periode d'un mode electromecanique a
#: 2 Hz. La fenetre perd donc le transitoire rapide sans rien perdre de ce
#: qu'on cherche.
RETARD_DEFAUT = 5e-3

#: Tolerance relative sur la regularite du pas (cf. :func:`_pas_uniforme`).
_TOL_PAS = 1e-6

#: Deux valeurs propres plus proches que ceci EN RELATIF sont le meme mode.
#: Ne sert qu'a reconnaitre une paire conjuguee, qui est exacte au bit pres :
#: ``numpy.linalg.eig`` d'une matrice REELLE rend ``(a + bj)`` et ``(a - bj)``
#: avec les memes ``a`` et ``b``, et les Hankel manipules ici sont reels.
_TOL_CONJUGUE = 1e-8


# ----------------------------------------------------------------------
# le resultat
# ----------------------------------------------------------------------
@dataclass
class ModeIdentifie:
    """Un mode lu sur le signal, avec de quoi le comparer a un ``Mode``.

    ``f`` et ``zeta`` sont calcules par :func:`g2egrid.analyse.modal.frequence` et
    :func:`g2egrid.analyse.modal.amortissement`, les MEMES fonctions que pour le
    spectre de ``A``. Ce n'est pas de la coquetterie : il existe au moins
    trois definitions courantes de l'amortissement (``-Re/|lam|``,
    ``-Re/Im``, le decrement logarithmique) et deux de la frequence (amortie
    et propre). Deux d'entre elles qui divergeraient de part et d'autre de la
    comparaison donneraient un ecart de quelques pour cent parfaitement
    stable, parfaitement reproductible, et parfaitement faux.

    ``amplitude`` et ``phase`` sont le module et l'argument du residu
    complexe, voie par voie, pour la reconstruction ::

        y_c(t) = somme_i  amplitude[c] . exp(j.phase[c]) . exp(lam_i . t)

    ``t`` y compte a partir du DEBUT DE LA FENETRE d'identification, pas de
    l'origine des temps de la simulation (cf. :func:`reconstruire`). Un mode
    reel et son conjugue apparaissent tous deux dans la liste, chacun avec la
    moitie de l'amplitude physique : l'oscillation visible d'une paire
    conjuguee vaut ``2.amplitude.exp(Re(lam).t).cos(Im(lam).t + phase)``.

    ``energie`` est la part de ce mode dans l'energie des contributions
    modales, en %, sur la fenetre et toutes voies confondues. C'est le
    critere de tri, et le seul qui separe un mode physique d'un mode de
    troncature : ce dernier ajuste du bruit, donc porte une energie minuscule
    (typiquement moins de 0,1 %), meme quand sa frequence est plausible.
    """

    valeur: complex
    f: float  # frequence amortie [Hz]
    zeta: float  # amortissement [%]
    amplitude: np.ndarray  # (n_voies,)
    phase: np.ndarray  # (n_voies,) [deg]
    energie: float  # [%] de l'energie modale totale
    voies: list[str] = field(default_factory=list)

    @property
    def stable(self) -> bool:
        return self.valeur.real < 0

    @property
    def oscillant(self) -> bool:
        return abs(self.valeur.imag) > 1e-6

    @property
    def residu(self) -> np.ndarray:
        """Le residu complexe par voie, ``amplitude . exp(j.phase)``."""
        return self.amplitude * np.exp(1j * np.radians(self.phase))

    @property
    def voie_dominante(self) -> str:
        if not len(self.amplitude) or not self.voies:
            return ""
        return self.voies[int(np.argmax(self.amplitude))]

    def __repr__(self) -> str:
        return (
            f"<ModeIdentifie {self.valeur:+.4g} : {self.f:.4g} Hz, "
            f"zeta = {self.zeta:.3g} %, energie = {self.energie:.3g} %>"
        )


@dataclass
class Identification:
    """Ce que :func:`identifier` a lu sur le signal, et de quoi en juger.

    ``residu`` est le garde-fou honnete : la norme de l'ecart entre le signal
    identifie et sa resynthese, rapportee a la norme du signal. Un rang trop
    bas rend des modes parfaitement plausibles qui ne reconstruisent rien, et
    rien dans les valeurs propres seules ne le dit. Au-dela de quelques
    pour-cent, les modes rendus ne sont pas ceux du signal.

    CE QU'IL NE DIT PAS -- une fenetre trop COURTE. Sur wscc9 a rang impose
    16, 0,6 s de fenetre donnent un residu de 0,0002 %, dix fois meilleur
    que les 0,0033 % de la fenetre de 4 s, et un amortissement faux de 56 %.
    C'est logique et il faut le savoir : sur une demi-periode, une
    exponentielle amortie et une exponentielle plus raide de residu plus
    grand sont indiscernables, et le modele identifie explique parfaitement
    ce qu'on lui a montre. Le residu mesure l'adequation au signal FOURNI,
    pas la pertinence de ce qu'on a choisi de fournir. Ce qui garde de cette
    erreur-la est le nombre de periodes du mode le plus lent dans la
    fenetre, et rien d'autre.

    ``valeurs_singulieres`` sert a relire le choix du rang apres coup : c'est
    la ou se voit le decrochement entre ce que porte le signal et ce qui
    n'est que du bruit.
    """

    modes: list[ModeIdentifie]
    methode: str
    rang: int
    valeurs_singulieres: np.ndarray
    residu: float
    t: np.ndarray  # grille (absolue) effectivement utilisee
    y: np.ndarray  # (n_t, n_voies), offset deja retranche
    offset: np.ndarray  # (n_voies,) ce qui a ete retranche
    voies: list[str] = field(default_factory=list)

    @property
    def t0(self) -> float:
        """L'origine des residus : le premier instant de la fenetre."""
        return float(self.t[0])

    def reconstruction(self) -> np.ndarray:
        """Le signal resynthetise sur la grille d'identification."""
        return reconstruire(self.modes, self.t - self.t0)

    def oscillants(self, energie_min: float = 0.0) -> list[ModeIdentifie]:
        """Les modes oscillants de partie imaginaire POSITIVE, par energie.

        L'identification rend les deux membres de chaque paire conjuguee ;
        ils portent la meme information. Pour lire un tableau, on n'en garde
        qu'un, par convention celui de partie imaginaire positive -- la meme
        que ``modal.modes(..., oscillants=True)``.
        """
        return [
            m for m in self.modes if m.valeur.imag > 1e-6 and m.energie >= energie_min
        ]

    def __repr__(self) -> str:
        return (
            f"<Identification {self.methode} : rang {self.rang}, "
            f"{len(self.modes)} modes, residu {100 * self.residu:.3g} %, "
            f"{len(self.t)} points sur [{self.t[0]:g}, {self.t[-1]:g}] s>"
        )


# ----------------------------------------------------------------------
# preparation du signal
# ----------------------------------------------------------------------
def _pas_uniforme(t: np.ndarray) -> float:
    """Le pas de ``t``, ou une erreur qui dit quoi faire.

    Matrix pencil comme ERA estiment un operateur de decalage d'UN pas :
    ``y[k+1] = somme z_i . r_i . z_i^k``. Cet operateur n'existe que si le
    pas est le meme partout. Or ``Radau`` choisit ses pas, et sur un reseau
    dq ils varient de plusieurs decades le long d'une meme trajectoire --
    resserres dans le transitoire, relaches ensuite. Identifier sur cette
    grille-la ne donne pas un resultat degrade : il donne un resultat faux,
    et silencieusement, parce que chaque ligne de la matrice de Hankel y
    porte un decalage temporel different.

    On IMPOSE donc l'uniformite plutot que de la rattraper : la trajectoire
    se demande avec ``simuler(..., t_eval=np.linspace(...))``, qui rend
    exactement les instants voulus sans rien couter a l'integration (ce sont
    des interpolants d'ordre eleve de Radau, pas des pas supplementaires).
    C'est plus propre que de reinterpoler apres coup, ou l'on n'a plus que la
    trajectoire echantillonnee. :func:`reechantillonner` reste la pour les
    enregistrements dont on ne maitrise pas la grille.
    """
    t = np.asarray(t, dtype=float).ravel()
    if t.size < 4:
        raise ValueError(f"{t.size} points : il en faut au moins 4 pour identifier")
    dt = np.diff(t)
    if dt.min() <= 0.0:
        raise ValueError("les instants ne sont pas strictement croissants")
    moyen = float(dt.mean())
    ecart = float(np.abs(dt - moyen).max() / moyen)
    if ecart > _TOL_PAS:
        raise ValueError(
            f"pas non uniforme : il varie de {100 * ecart:.3g} % autour de "
            f"{moyen:g} s (de {dt.min():g} a {dt.max():g} s). Matrix pencil "
            "et ERA supposent un echantillonnage uniforme : simuler avec "
            "t_eval=np.linspace(...), ou passer par "
            "identification.reechantillonner()."
        )
    return moyen


def reechantillonner(
    t: np.ndarray, y: np.ndarray, n: int | None = None, dt: float | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Remet ``(t, y)`` sur une grille uniforme, par spline cubique.

    Utile quand on n'a pas la main sur la grille : un enregistrement, ou une
    trajectoire deja simulee sans ``t_eval``. Quand on l'a, mieux vaut
    demander directement les bons instants (cf. :func:`_pas_uniforme`).

    POURQUOI PAS UNE INTERPOLATION LINEAIRE -- elle serait tentante, et elle
    est mauvaise ici pour une raison precise : sur un pas ``h``, elle se
    comporte comme un filtre de gain ``sinc`` et d'erreur ``O(h^2.y'')``, ce
    qui retranche de l'amplitude a chaque sommet d'oscillation. Sur une
    grille a pas VARIABLE, elle en retranche donc davantage la ou les pas
    sont longs -- c'est-a-dire dans la queue de la trajectoire, la ou
    l'oscillation est deja petite. Le signal reconstruit decroit ainsi plus
    vite que le vrai, et l'amortissement identifie sort trop grand, d'un
    biais systematique que rien dans le resultat ne signale. La spline
    cubique, en ``O(h^4)``, rend ce biais negligeable devant les autres
    erreurs pour le meme prix.

    ``n`` fixe le nombre de points (defaut : autant qu'il y en avait),
    ``dt`` le pas ; l'un ou l'autre, pas les deux.
    """
    t = np.asarray(t, dtype=float).ravel()
    y = np.asarray(y, dtype=float)
    y = y.reshape(len(t), -1)
    if n is not None and dt is not None:
        raise ValueError("reechantillonner : donner n OU dt, pas les deux")
    if dt is not None:
        n = int(np.floor((t[-1] - t[0]) / float(dt))) + 1
        tu = t[0] + float(dt) * np.arange(n)
    else:
        tu = np.linspace(t[0], t[-1], int(n) if n else len(t))
    return tu, CubicSpline(t, y, axis=0)(tu)


def _fenetre(
    t: np.ndarray, y: np.ndarray, debut: float | None, fin: float | None
) -> tuple[np.ndarray, np.ndarray]:
    """Coupe la fenetre d'identification.

    POURQUOI NE PAS COMMENCER A L'EVENEMENT -- un reseau dq porte des modes a
    1e6-1e7 rad/s : resistances de causalite, troncons de cable courts,
    boucles de courant d'electronique. Un echelon de consigne les excite peu,
    mais il les excite, et ils sont eteints en quelques microsecondes. Vus a
    un pas de quelques millisecondes, ils ne sont pas "petits" : ils sont
    REPLIES, et le premier point de la fenetre porte un residu qui ne
    correspond a aucune frequence representable. La matrice de Hankel en
    herite quelques lignes aberrantes qui gonflent les premieres valeurs
    singulieres, deplacent la troncature, et font sortir des poles de rebut
    au voisinage du cercle unite -- la ou l'on cherche justement les modes
    electromecaniques. Attendre ``RETARD_DEFAUT`` apres l'evenement supprime
    le probleme a la source : ce qui est mort n'a plus besoin d'etre
    identifie.

    CE QUE LA MESURE DIT REELLEMENT -- sur wscc9, ce retard ne change
    RIEN : de 0 a 300 ms de decalage, le residu de reconstruction va de
    0,0035 % a 0,0017 % et les deux modes electromecaniques ne bougent pas
    de plus de 0,15 %. La raison est que l'echelon y passe par la consigne
    d'un regulateur de vitesse, qui n'excite les modes a 1e6 rad/s qu'a un
    niveau negligeable, et que ``simuler`` decoupe l'integration a
    l'instant de l'echelon -- Radau y redemarre proprement au lieu de
    franchir la discontinuite avec un pas trop grand. Le retard par defaut
    est donc ici une assurance, pas un correctif : il reste justifie pour un
    court-circuit, une manoeuvre d'organe ou un enregistrement reel, ou le
    contenu rapide est d'un tout autre ordre.
    """
    t = np.asarray(t, dtype=float).ravel()
    y = np.asarray(y, dtype=float).reshape(len(t), -1)
    a = t[0] + RETARD_DEFAUT if debut is None else float(debut)
    b = t[-1] if fin is None else float(fin)
    garde = (t >= a - 1e-12) & (t <= b + 1e-12)
    if garde.sum() < 4:
        raise ValueError(
            f"fenetre [{a:g}, {b:g}] s : {int(garde.sum())} points sur les "
            f"{len(t)} fournis ([{t[0]:g}, {t[-1]:g}] s)"
        )
    return t[garde], y[garde]


def _retrancher_offset(y: np.ndarray, offset: bool) -> tuple[np.ndarray, np.ndarray]:
    """Ramene le signal a zero en fin de fenetre.

    Apres un echelon de consigne, le systeme ne revient pas d'ou il vient :
    il s'arrete sur un AUTRE equilibre. Le signal contient donc une
    composante continue, et un modele en somme d'exponentielles ne peut la
    representer qu'en s'inventant un pole en zero -- c'est la seconde
    solution, legitime, mais elle coute un rang, et surtout elle laisse ce
    pole capter aussi la derive residuelle des modes tres lents, ce qui
    pollue son voisinage. Retrancher la valeur finale est plus robuste : on
    enleve exactement ce qu'on sait etre la, et le modele identifie n'a plus
    a representer que le transitoire.

    On retranche la MOYENNE des derniers points plutot que le dernier point :
    avec du bruit de mesure, le dernier point porte tout le bruit, et cette
    erreur-la se retrouve integralement dans l'offset. Une moyenne sur 2 % de
    la fenetre le divise par plusieurs, et reste courte devant la periode du
    mode le plus lent qu'on cherche.
    """
    if not offset:
        return y, np.zeros(y.shape[1])
    queue = max(1, len(y) // 50)
    val = y[-queue:].mean(axis=0)
    return y - val, val


# ----------------------------------------------------------------------
# algebre commune : Hankel, rang, residus
# ----------------------------------------------------------------------
def _hankel(y: np.ndarray, n_lignes: int, n_cols: int, k0: int = 0) -> np.ndarray:
    """Hankel par BLOCS de lignes : ``H[i.C + c, j] = y[k0 + i + j, c]``.

    Le multi-voies n'est pas un luxe. D'abord parce qu'il conditionne mieux
    le probleme : a nombre d'echantillons egal, trois voies donnent trois
    fois plus de lignes a la matrice de Hankel pour le meme nombre
    d'inconnues. Ensuite et surtout parce qu'un mode n'est pas observable de
    la meme facon partout : un mode inter-machines ou deux machines oscillent
    en opposition est presque invisible sur la somme de leurs vitesses, et
    parfaitement visible sur leur difference. Identifier sur une seule voie,
    c'est accepter de rater les modes dont cette voie se trouve etre un noeud.

    L'empilement par blocs de lignes convient aux deux methodes : il ne
    change pas l'espace ligne (donc les vecteurs singuliers a droite du
    pencil), et il donne a ERA la structure de Hankel par blocs de ses
    parametres de Markov.
    """
    n_t, n_voies = y.shape
    if k0 + n_lignes + n_cols - 1 > n_t:
        raise ValueError(
            f"Hankel {n_lignes}x{n_cols} decale de {k0} : il faudrait "
            f"{k0 + n_lignes + n_cols - 1} echantillons, il y en a {n_t}"
        )
    i = np.arange(n_lignes)
    j = np.arange(n_cols)
    bloc = y[k0 + i[:, None] + j[None, :], :]  # (n_lignes, n_cols, n_voies)
    return bloc.transpose(0, 2, 1).reshape(n_lignes * n_voies, n_cols)


def _forme_hankel(n_t: int, n_voies: int, methode: str) -> int:
    """Le parametre de forme du Hankel : ``L`` pour le pencil, ``alpha`` pour ERA.

    Les deux methodes decoupent les memes ``N`` echantillons en un Hankel par
    blocs, et le seul choix libre est l'endroit de la coupe : le pencil a
    ``C.(N - L)`` lignes pour ``L + 1`` colonnes, ERA a ``C.alpha`` lignes
    pour ``N - alpha`` colonnes. Dans les deux cas, la regle qui sort de la
    mesure est la meme : **le Hankel doit etre carre**. D'ou ::

        pencil : C.(N - L) = L + 1    ->  L     ~ C.N / (C + 1)
        ERA    : C.alpha   = N - alpha ->  alpha ~ N / (C + 1)

    Sur les trois vitesses rotor de wscc9 (``C = 3``, ``N = 781``), a rang
    impose 16, le residu de reconstruction en fonction de ``L`` vaut 1,30 %
    a ``N/4``, 0,54 % a ``N/3``, 0,020 % a ``N/2``, 0,004 % a ``2N/3``,
    0,003 % a ``3N/4`` -- la valeur carree -- puis remonte a 0,032 % a
    ``0,9 N``. ERA suit exactement la courbe symetrique : 0,003 % a
    ``alpha = N/4``, la valeur carree, et 2,3 % a ``alpha = 0,9 N``. Les
    deux estimateurs se degradent donc de la meme facon -- quand le Hankel
    devient haut et mince -- et pour la meme raison : la troncature au rang
    ``M`` n'a plus assez de colonnes independantes pour separer ``M`` modes.

    Le nombre de VOIES entre dans la regle, et c'est ce qui la distingue de
    la recommandation ``N/3 <= L <= N/2`` de Hua & Sarkar : celle-ci est
    ecrite pour une seule voie, ou ``C = 1`` la ramene precisement a ``L ~
    N/2``, c'est-a-dire au meme Hankel carre. A trois voies, s'y tenir
    triplerait le nombre de lignes sans toucher aux colonnes.
    """
    if methode == "pencil":
        L = int(round((n_voies * n_t - 1) / (n_voies + 1)))
        return int(np.clip(L, 2, n_t - 2))
    return int(np.clip(round(n_t / (n_voies + 1)), 2, n_t - 2))


def _rang(s: np.ndarray, rang: int | None, seuil: float, rang_max: int) -> int:
    """Le rang de troncature, impose ou choisi sur les valeurs singulieres.

    Le rang fixe le NOMBRE de modes identifies : c'est le seul parametre qui
    compte vraiment, et le seul qu'on ne sache pas deduire du signal sans
    convention. Celle retenue est la plus lisible : on garde les valeurs
    singulieres qui valent plus de ``seuil`` fois la plus grande. Le rapport
    est sans dimension, il ne depend donc ni de l'unite du signal ni de
    l'amplitude de la perturbation -- deux choses qui changent a chaque
    etude.

    LE DEFAUT EST DELIBEREMENT BAS. Sur les trois vitesses rotor de wscc9
    apres un echelon de 0,01 pu, le spectre singulier descend en marches --
    1 ; 0,385 ; 0,320 ; 0,101 ; 0,081 ; ... ; 1,8e-4 ; 6,1e-5 ; 2,2e-5 ;
    1,0e-5 ; 6,7e-6 -- puis s'aplatit sur un plancher a 4e-6, qui est le
    bruit de l'integration a ``rtol = 1e-6``. Un seuil a 1e-2 coupe a 9, et
    rate le critere de validation sur les deux modes electromecaniques
    (-4,6 % en frequence sur l'un, +28 % d'amortissement sur l'autre) ; a
    1e-3 il coupe a 12, ce qui les trouve a 4 % d'amortissement pres ; a
    1e-4 il coupe a 14, juste au-dessus du plancher, et l'erreur tombe a
    0,5 % au plus.

    SUR-ESTIMER LE RANG EST LE MOINDRE MAL. Descendre le seuil a 1e-6, donc
    plonger dans le plancher de bruit et retenir 68 modes au lieu de 14,
    deplace les deux modes electromecaniques de moins de 0,2 % en frequence
    et 0,3 % en amortissement : les modes surnumeraires se partagent le
    bruit sans toucher a ce qui porte de l'energie. Le manque, lui, ne
    pardonne pas -- a rang 8, le mode a 2,2 Hz ne sort pas du tout. Sur un
    enregistrement reel, ou le plancher est plus haut, ce seuil se remonte.
    """
    if rang is not None:
        r = int(rang)
        if not 1 <= r <= rang_max:
            raise ValueError(f"rang {r} hors de [1, {rang_max}]")
        return r
    if s[0] <= 0.0:
        raise ValueError("signal identiquement nul : rien a identifier")
    r = int(np.count_nonzero(s > float(seuil) * s[0]))
    return int(np.clip(r, 1, rang_max))


def _puissances(z: np.ndarray, k: np.ndarray) -> np.ndarray:
    """``z_i^k``, calcule par ``exp(k.log z)`` avec garde-fou d'exposant.

    Rien n'empeche une identification bruitee de rendre un pole legerement a
    l'exterieur du cercle unite ; eleve a la puissance de quelques milliers
    d'echantillons, il deborde. On plafonne donc l'exposant reel : la colonne
    correspondante sature au lieu de rendre des ``inf`` qui contamineraient
    tout le moindres carres.
    """
    with np.errstate(divide="ignore"):
        log_z = np.log(z.astype(complex))
    expo = np.outer(k, log_z)
    return np.exp(np.clip(expo.real, -700.0, 700.0) + 1j * expo.imag)


def _residus(y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, float]:
    """Les residus par voie, au sens des moindres carres, et le residu relatif.

    Les poles une fois connus, le probleme redevient LINEAIRE en les
    amplitudes : ``y[k, c] = somme_i R[i, c] . z_i^k``. C'est le seul endroit
    ou l'on utilise tous les echantillons de toutes les voies ensemble, et
    c'est ce qui donne un sens au residu de reconstruction : il compare la
    MEILLEURE resynthese possible avec ces poles-la au signal qu'on lui a
    donne.
    """
    n_t = len(y)
    Z = _puissances(z, np.arange(n_t))
    R, *_ = np.linalg.lstsq(Z, y.astype(complex), rcond=None)
    ecart = y - (Z @ R).real
    norme = float(np.linalg.norm(y))
    return R, float(np.linalg.norm(ecart) / norme) if norme > 0 else 0.0


def _energies(R: np.ndarray, z: np.ndarray, n_t: int) -> np.ndarray:
    """Part de chaque mode dans l'energie des contributions modales [%].

    L'energie de la contribution d'un mode, sommee sur la fenetre et sur les
    voies, ne depend ni de l'unite ni du nombre d'echantillons une fois
    rapportee au total. Elle trie ce qui explique le signal avant ce qui
    l'ajuste : un mode de troncature a une amplitude minuscule OU une
    decroissance si rapide qu'il ne vit que sur deux echantillons, et les
    deux se voient ici, alors que ni la frequence ni l'amortissement ne les
    trahissent.
    """
    Z = _puissances(z, np.arange(n_t))  # (n_t, n_modes)
    e = (np.abs(Z) ** 2).sum(axis=0) * (np.abs(R) ** 2).sum(axis=1)
    total = float(e.sum())
    return 100.0 * e / total if total > 0 else e


# ----------------------------------------------------------------------
# les deux estimateurs
# ----------------------------------------------------------------------
def _poles_pencil(
    y: np.ndarray, L: int, rang: int | None, seuil: float
) -> tuple[np.ndarray, np.ndarray, int]:
    """Matrix pencil (Hua & Sarkar 1990), multi-voies.

    Une matrice de Hankel ``Y`` de ``L + 1`` colonnes se factorise
    ``Y = Z1 . diag(r) . Z2`` avec ``Z2[i, j] = z_i^j`` : deux sous-matrices
    obtenues en otant a ``Y`` sa derniere puis sa premiere colonne sont donc
    liees par ``Y2 = Z1 . diag(r) . diag(z) . Z2'``, et les ``z_i`` sont les
    valeurs propres generalisees du faisceau ``(Y2, Y1)``.

    On ne forme pas ce faisceau : la SVD de ``Y`` tronquee au rang ``M`` rend
    ``V``, dont les ``M`` premieres colonnes engendrent l'espace ligne, et le
    decalage de colonne devient un decalage de LIGNE sur ``V``. Les ``z_i``
    sont alors les valeurs propres de ``pinv(V1) . V2``, ou ``V1`` et ``V2``
    sont ``V[:, :M]`` prive respectivement de sa derniere et de sa premiere
    ligne. C'est la forme "SVD" de Hua & Sarkar : elle ne manipule qu'une
    matrice ``M x M``, et la troncature y fait tout le debruitage.

    ``L``, le nombre de colonnes moins une, est le parametre du faisceau ;
    :func:`_forme_hankel` dit comment il se choisit.
    """
    n_t = y.shape[0]
    Y = _hankel(y, n_t - L, L + 1)
    U, s, Vh = np.linalg.svd(Y, full_matrices=False)
    M = _rang(s, rang, seuil, min(L, Y.shape[0] - 1))
    V = Vh[:M].conj().T  # (L + 1, M)
    z = np.linalg.eigvals(np.linalg.pinv(V[:-1]) @ V[1:])
    return z, s, M


def _poles_era(
    y: np.ndarray, alpha: int, beta: int, rang: int | None, seuil: float
) -> tuple[np.ndarray, np.ndarray, int]:
    """ERA (Juang & Pappa 1985), multi-voies.

    Meme point de depart -- une matrice de Hankel de reponse libre -- mais le
    decalage est pris dans les LIGNES, sur une SECONDE matrice de Hankel
    ``H1`` decalee d'un echantillon. La SVD tronquee de ``H0`` fournit une
    base equilibree de l'espace d'etat, et la realisation minimale s'y ecrit
    ::

        A_r = S^-1/2 . U^T . H1 . V . S^-1/2

    dont les valeurs propres sont les ``z_i``. La racine carree de ``S`` de
    part et d'autre n'est pas une commodite : c'est elle qui rend la
    realisation equilibree, donc aussi bien conditionnee en observabilite
    qu'en commandabilite. Une realisation non equilibree donne les MEMES
    valeurs propres en arithmetique exacte, et des valeurs propres
    sensiblement differentes en flottant des que le rang monte.

    Difference pratique avec le pencil : ERA lit le decalage sur ``H1``, qui
    est une matrice de donnees, la ou le pencil le lit sur ``V``, qui est
    deja debruitee par la troncature. A forme de Hankel egale, les deux
    rendent pourtant la meme chose a la precision ou on les lit (sur wscc9,
    les memes 0,013 % de residu et les memes ecarts a 0,01 % pres) : c'est
    rassurant, et c'est ce qui donne son sens au croisement des deux. Ce qui
    les separe est la FORME optimale du Hankel, qui est inverse de l'une a
    l'autre (cf. :func:`_forme_hankel`) : un ``L`` bon pour le pencil est
    mauvais pour ERA. Les faire tourner toutes deux avec le meme ``L``
    impose, au lieu de laisser chacune prendre le sien, est donc le meilleur
    moyen de croire a un desaccord qui n'est qu'un mauvais reglage.
    """
    H0 = _hankel(y, alpha, beta, 0)
    H1 = _hankel(y, alpha, beta, 1)
    U, s, Vh = np.linalg.svd(H0, full_matrices=False)
    M = _rang(s, rang, seuil, min(H0.shape) - 1)
    racine = 1.0 / np.sqrt(s[:M])
    A_r = (racine[:, None] * (U[:, :M].T @ H1 @ Vh[:M].T)) * racine[None, :]
    return np.linalg.eigvals(A_r), s, M


# ----------------------------------------------------------------------
# l'entree principale
# ----------------------------------------------------------------------
def identifier(
    t: np.ndarray,
    y: np.ndarray,
    methode: str = "pencil",
    debut: float | None = None,
    fin: float | None = None,
    rang: int | None = None,
    seuil_rang: float = 1e-4,
    offset: bool = True,
    L: int | None = None,
    voies: list[str] | None = None,
) -> Identification:
    """Les modes portes par un signal de reponse libre.

    ``y`` est ``(n_t,)`` ou ``(n_t, n_voies)`` : plusieurs voies valent mieux
    qu'une (cf. :func:`_hankel`). ``t`` doit etre a pas uniforme -- demander
    la trajectoire avec ``simuler(..., t_eval=np.linspace(...))`` -- et
    commencer A l'evenement : la fenetre demarre ``RETARD_DEFAUT`` plus tard,
    sauf si ``debut`` en decide autrement.

    ``methode`` vaut ``"pencil"`` ou ``"era"``. Les faire tourner toutes deux
    sur le meme signal est le meilleur controle disponible : ce qu'elles
    rendent en commun est dans le signal, le reste est dans la troncature.

    ``rang`` force le nombre de valeurs singulieres retenues, donc le nombre
    de modes ; a defaut il se choisit sur ``seuil_rang`` (cf. :func:`_rang`).
    ``offset=False`` desactive le retrait de la valeur finale, ce qui n'a de
    sens que sur un signal dont on sait qu'il retourne a zero.

    Le passage au continu est ``lam = log(z) / dt``. ``log`` prend sa
    determination principale : une frequence au-dela de Nyquist ne peut pas
    etre distinguee de son repliee, et c'est la repliee qui sort. C'est une
    raison de plus de ne pas commencer la fenetre a l'evenement.
    """
    if methode not in ("pencil", "era"):
        raise ValueError(f"methode {methode!r} inconnue ('pencil' ou 'era')")
    t = np.asarray(t, dtype=float).ravel()
    y = np.asarray(y, dtype=float).reshape(len(t), -1)
    noms = list(voies) if voies is not None else [f"y{i}" for i in range(y.shape[1])]
    if len(noms) != y.shape[1]:
        raise ValueError(f"{len(noms)} noms de voies pour {y.shape[1]} colonnes")

    tw, yw = _fenetre(t, y, debut, fin)
    dt = _pas_uniforme(tw)
    yw, val = _retrancher_offset(yw, offset)

    n_t = len(yw)
    Lp = int(L) if L is not None else _forme_hankel(n_t, yw.shape[1], methode)
    if not 2 <= Lp <= n_t - 2:
        raise ValueError(f"L = {Lp} hors de [2, {n_t - 2}]")
    if methode == "pencil":
        z, s, M = _poles_pencil(yw, Lp, rang, seuil_rang)
    else:
        z, s, M = _poles_era(yw, Lp, n_t - Lp, rang, seuil_rang)

    # Un pole nul (mode eteint en un pas, ou rang superieur au rang reel du
    # signal) n'a pas d'image dans le plan continu : log(0) = -inf. On
    # l'ecarte plutot que de propager des infinis dans les residus.
    z = z[np.abs(z) > 1e-300]
    lam = np.log(z) / dt

    R, residu = _residus(yw, z)
    energie = _energies(R, z, n_t)
    f = frequence(lam)
    zeta = amortissement(lam)

    modes = [
        ModeIdentifie(
            valeur=complex(lam[i]),
            f=float(f[i]),
            zeta=float(zeta[i]),
            amplitude=np.abs(R[i]),
            phase=np.degrees(np.angle(R[i])),
            energie=float(energie[i]),
            voies=noms,
        )
        for i in range(len(lam))
    ]
    modes.sort(key=lambda m: -m.energie)
    return Identification(
        modes=modes,
        methode=methode,
        rang=M,
        valeurs_singulieres=s,
        residu=residu,
        t=tw,
        y=yw,
        offset=val,
        voies=noms,
    )


def reconstruire(modes: list[ModeIdentifie], t: np.ndarray) -> np.ndarray:
    """Resynthetise le signal a partir des modes identifies.

    ``t`` compte a partir du DEBUT DE LA FENETRE d'identification, pas de
    l'origine des temps de la simulation : les residus y sont references
    (``Identification.t0`` donne le decalage). Ce n'est pas un detail de
    convention -- ramener les residus a l'origine des temps reviendrait a les
    multiplier par ``exp(-lam.t0)``, c'est-a-dire a faire exploser ceux des
    modes amortis, et a perdre en chiffres significatifs exactement ce que le
    mode a d'amortissement.

    Rend ``(len(t), n_voies)``, reel : les parties imaginaires se compensent
    entre paires conjuguees.
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    if not modes:
        raise ValueError("aucun mode a reconstruire")
    lam = np.array([m.valeur for m in modes], dtype=complex)
    R = np.array([m.residu for m in modes], dtype=complex)  # (n_modes, n_voies)
    expo = np.outer(t, lam)
    E = np.exp(np.clip(expo.real, -700.0, 700.0) + 1j * expo.imag)
    return (E @ R).real


def sur_trajectoire(
    tr,
    voies: Sequence[str],
    evenement: float | None = None,
    methode: str = "pencil",
    **kw,
) -> Identification:
    """Identifier directement sur une :class:`g2egrid.temporel.Trajectoire`.

    C'est la porte d'entree de la validation du petit signal par l'EMT :
    simuler, identifier, apparier au spectre de ``A``. :func:`identifier`
    travaille sur des tableaux ``(t, y)``, ce qui est le bon niveau pour un
    enregistrement venu d'ailleurs ; mais l'appelant normal a une
    trajectoire et des noms d'etats, et lui laisser refaire a chaque fois la
    decoupe et l'empilement est le meilleur moyen que deux appelants ne la
    fassent pas pareil.

    ``evenement`` est la date de la perturbation : tout ce qui PRECEDE est
    jete. Avant l'evenement le signal est plat, et une fenetre qui
    l'engloberait demanderait au modele en somme d'exponentielles de
    representer aussi le saut -- ce qu'il ne peut pas, et il paierait en
    poles de rebut. A defaut, la trajectoire est prise telle quelle, ce qui
    suppose qu'elle commence a l'evenement.

    Le pas doit etre uniforme, donc la trajectoire demandee avec un
    ``t_eval`` regulier ; :func:`identifier` le verifie et le dit.
    """
    voies = list(voies)
    inconnues = [n for n in voies if n not in tr.etats]
    if inconnues:
        raise ValueError(
            f"etat(s) absent(s) de la trajectoire : {', '.join(inconnues)}"
        )
    garde = (
        np.ones(len(tr.t), dtype=bool)
        if evenement is None
        else tr.t >= float(evenement) - 1e-12
    )
    if garde.sum() < 4:
        raise ValueError(
            f"evenement a {evenement} s : il ne reste que {int(garde.sum())} "
            f"points sur [{tr.t[0]:g}, {tr.t[-1]:g}] s"
        )
    y = np.column_stack([tr.etat(n)[garde] for n in voies])
    return identifier(tr.t[garde], y, methode=methode, voies=voies, **kw)


# ----------------------------------------------------------------------
# appariement au spectre du modele lineaire
# ----------------------------------------------------------------------
@dataclass
class Paire:
    """Un mode identifie en face d'un mode du modele lineaire.

    ``ecart_f`` est en % de la frequence PREDITE, ``ecart_zeta`` en POINTS
    d'amortissement (zeta est deja un pourcentage : un mode predit a 8,7 % et
    identifie a 9,4 % fait ``+0,7`` point), et ``ecart_zeta_rel`` le meme en
    % de la valeur predite (ici ``+8 %``). Les deux sont donnes parce qu'ils
    ne servent pas a la meme chose : le point d'amortissement est ce qui se
    compare d'un mode a l'autre et ce qui s'ecrit dans un cahier des charges
    ("au moins 5 %") ; le relatif est ce qui se compare d'un modele a
    l'autre, et c'est donc lui qu'on tolere a quelques dizaines de %.
    """

    identifie: ModeIdentifie
    lineaire: Mode
    distance: float  # sans dimension, cf. apparier
    ecart_f: float  # [%] de la frequence predite
    ecart_zeta: float  # [points] d'amortissement
    ecart_zeta_rel: float  # [%] de l'amortissement predit


@dataclass
class Appariement:
    """Le resultat de :func:`apparier` : ce qui se correspond, et le reste.

    Les deux listes de laisses-pour-compte ne sont pas un detail
    d'implementation : un mode predit que l'identification ne trouve pas est
    soit inobservable sur les voies choisies, soit inexcite par la
    perturbation, soit absent du systeme reel -- trois diagnostics tres
    differents, et aucun ne se lit sur les paires.
    """

    paires: list[Paire]
    identifies_seuls: list[ModeIdentifie]
    lineaires_seuls: list[Mode]

    def tableau(self) -> str:
        """Le tableau "predit / identifie / ecart", en texte."""
        lignes = [
            f"{'predit f [Hz]':>13} {'zeta [%]':>9} | "
            f"{'identifie f':>12} {'zeta [%]':>9} | "
            f"{'d f [%]':>8} {'d zeta [pt]':>11} {'energie [%]':>11}"
        ]
        for p in self.paires:
            lignes.append(
                f"{p.lineaire.f:13.4f} {p.lineaire.zeta:9.3f} | "
                f"{p.identifie.f:12.4f} {p.identifie.zeta:9.3f} | "
                f"{p.ecart_f:8.2f} {p.ecart_zeta:11.3f} "
                f"{p.identifie.energie:11.3f}"
            )
        for m in self.lineaires_seuls:
            lignes.append(
                f"{m.f:13.4f} {m.zeta:9.3f} | {'-':>12} {'-':>9} | "
                f"{'non identifie':>32}"
            )
        for m in self.identifies_seuls:
            lignes.append(
                f"{'-':>13} {'-':>9} | {m.f:12.4f} {m.zeta:9.3f} | "
                f"{'non predit':>21} {m.energie:11.3f}"
            )
        return "\n".join(lignes)


def _valeur(m) -> complex:
    """La valeur propre d'un ``Mode``, d'un ``ModeIdentifie`` ou d'un complexe."""
    return complex(getattr(m, "valeur", m))


def _replier(modes: list) -> list:
    """Un seul representant par paire conjuguee, celui de partie imaginaire >= 0.

    CONVENTION -- ``modal.modes(..., oscillants=True)`` ne garde que la
    moitie superieure du spectre ; l'identification, elle, rend les deux
    membres de chaque paire. Comparer les deux listes telles quelles
    laisserait la moitie des modes identifies orphelins, et un appariement
    global n'est pas fait pour digerer autant de figurants. On replie donc
    les DEUX listes de la meme facon, ce qui rend :func:`apparier`
    indifferent a la convention de son appelant : que les conjugues soient
    presents ou non de chaque cote, le resultat est le meme.

    Le repliement ne fusionne que ce qui est conjugue au bit pres (cf.
    ``_TOL_CONJUGUE``) : deux modes physiquement voisins -- les 2,01 et
    2,20 Hz de wscc9 sont a 9 % l'un de l'autre -- restent bien distincts.
    """
    gardes: list = []
    replies: list[complex] = []
    for m in modes:
        v = _valeur(m)
        w = v.conjugate() if v.imag < 0 else v
        for k, deja in enumerate(replies):
            if abs(w - deja) <= _TOL_CONJUGUE * max(abs(deja), 1.0):
                if v.imag > _valeur(gardes[k]).imag:
                    gardes[k] = m  # on prefere le representant a Im >= 0
                break
        else:
            replies.append(w)
            gardes.append(m)
    return gardes


def apparier(
    identifies: list[ModeIdentifie],
    lineaires: list[Mode],
    seuil: float = 0.25,
    energie_min: float = 0.0,
) -> Appariement:
    """Apparie globalement les modes identifies aux modes predits.

    APPARIEMENT GLOBAL, PAS GLOUTON -- prendre a chaque tour la paire la plus
    proche encore disponible est tentant et faux : sur deux modes voisins, le
    glouton verrouille d'abord la paire la mieux ajustee, quitte a forcer
    l'autre mode a aller chercher un partenaire lointain. Le cout total y
    perd, et surtout le diagnostic est inverse -- on conclut qu'un mode
    predit "n'est pas la" alors qu'il est juste passe apres. L'algorithme
    hongrois (``scipy.optimize.linear_sum_assignment``) minimise la somme des
    distances, ce qui est exactement la question posee.

    DISTANCE SANS DIMENSION -- ``|lam_i - lam_l| / |lam_l|``. Un ecart absolu
    de 0,1 rad/s est enorme sur un mode a 0,05 Hz et invisible sur un mode a
    300 Hz : comparer des distances absolues revient a n'apparier
    correctement que les modes rapides. Rapporter au module de la valeur
    propre PREDITE (et non a une moyenne des deux) garde une reference qui ne
    depend pas de la qualite de l'identification, et mesure d'un seul nombre
    ce qui se separerait sinon en un ecart de frequence et un ecart
    d'amortissement -- les deux erreurs deplacent ``lam`` dans le plan
    complexe, et c'est leur somme geometrique qui dit si c'est le meme mode.

    ``seuil`` est la distance au-dela de laquelle on REFUSE d'apparier : le
    defaut 0,25 laisse passer un mode identifie a 20 % de sa frequence
    predite, ce qui est deja beaucoup, et rejette les rapprochements de
    hasard. ``energie_min`` [%] ecarte d'abord les modes de troncature, qui
    sont nombreux et qui, a seuil egal, monopoliseraient des partenaires.
    """
    ident = _replier([m for m in identifies if m.energie >= energie_min])
    lin = _replier(list(lineaires))
    if not ident or not lin:
        return Appariement([], list(ident), list(lin))

    vi = np.array([_valeur(m) for m in ident])
    vl = np.array([_valeur(m) for m in lin])
    d = np.abs(vi[:, None] - vl[None, :]) / np.maximum(np.abs(vl)[None, :], 1e-12)

    li, lj = linear_sum_assignment(d)
    paires, pris_i, pris_j = [], set(), set()
    for i, j in zip(li, lj):
        if d[i, j] > seuil:
            continue
        pris_i.add(int(i))
        pris_j.add(int(j))
        mi, ml = ident[int(i)], lin[int(j)]
        f_l, z_l = float(ml.f), float(ml.zeta)
        paires.append(
            Paire(
                identifie=mi,
                lineaire=ml,
                distance=float(d[i, j]),
                ecart_f=100.0 * (mi.f - f_l) / f_l if f_l > 1e-12 else np.inf,
                ecart_zeta=mi.zeta - z_l,
                ecart_zeta_rel=(
                    100.0 * (mi.zeta - z_l) / abs(z_l) if abs(z_l) > 1e-12 else np.inf
                ),
            )
        )
    paires.sort(key=lambda p: p.lineaire.f)
    return Appariement(
        paires=paires,
        identifies_seuls=[m for k, m in enumerate(ident) if k not in pris_i],
        lineaires_seuls=[m for k, m in enumerate(lin) if k not in pris_j],
    )
