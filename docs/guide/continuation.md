# Continuation : jusqu'où le réseau reste stable

Une étude modale répond pour **un** point de fonctionnement ; la question
pratique est souvent « jusqu'où ? » : jusqu'à quel gain de PLL, quelle
puissance d'un parc, et quel mode lâche en premier.
{mod}`g2egrid.scenarios.continuation` fait varier un paramètre $\mu$, suit chaque mode d'un
pas à l'autre et raffine par dichotomie les points où un mode traverse l'axe
imaginaire.

## Un statisme qui déstabilise : bifurcation de Hopf

```python
import numpy as np
from g2egrid import catalogue
from g2egrid import continuation as cont

net = catalogue.construire("smib_gfm")
mp0 = next(e for e in net.all_elements if e.name == "GFM2").pdict["mp"]   # 0,005

c = cont.continuer(cont.par_parametre(net, "GFM2", "mp"), mp0 * np.geomspace(1, 100, 10))
print(c.resume())
b = c.bifurcations[0]
print(b.genre, round(b.mu, 4), round(b.f, 1), b.sens, b.elements[:2])
```

```text
Hopf a mu = 0.0758397 (24.9 Hz, destabilise) : GFM2 88 %, IB1 12 %, L1_1_2 0 %
hopf 0.0758 24.9 1 [('GFM2', 88.01426362535187), ('IB1', 11.862106290476234)]
```

Au-delà d'un statisme P-f de 7,6 %, une oscillation à 24,9 Hz naît entre le
GFM et le bus infini.

`par_parametre(net, element, nom)` rend la fabrique $\mu \mapsto$ réseau où
`element.nom = mu` ; `equilibre=True` refait load flow et équilibre à chaque pas
(pour un paramètre qui déplace le point de fonctionnement). N'importe quelle
fonction $\mu \mapsto$ `Network` résolu convient aussi.

Le résultat, {class}`~g2egrid.scenarios.continuation.Continuation` :

| Champ | Contenu |
|---|---|
| `mu`, `lam` | les valeurs calculées, et le spectre suivi `(n_mu, n_modes)` : la colonne `j` est le même mode d'un bout à l'autre |
| `bifurcations` | `genre` (`"hopf"` ou `"noeud-col"`), `mu`, `valeur`, `f`, `elements`, `sens` (+1 : se déstabilise quand $\mu$ croît) |
| `marge` | le premier $\mu$ qui déstabilise |
| `arret_a` | le $\mu$ sans point de fonctionnement qui a arrêté la suite (`None` si elle est allée au bout) |
| `tensions` | $\lvert V\rvert$ de chaque bus à chaque pas |
| `depassements`, `hors_limite`, `mu_hors_limite` | les convertisseurs au-delà de leur limite à chaque pas, le masque des pas concernés, et le premier $\mu$ hors limite (voir plus bas) |

Sur un tirage de la campagne où le GFL est stable, le gain intégral de sa PLL
le déstabilise vers 54 Hz : la valeur trouvée, Kipll = 255,18, encadre le
changement de stabilité au millième (test `tests/test_continuation.py`).

## Le suivi des modes

D'un pas au suivant, les modes sont appariés par l'algorithme hongrois sur la
distance des valeurs propres et la ressemblance des participations, comme pour
la {doc}`robustesse <robustesse>` : un mode garde son indice même quand deux modes
se croisent. Un mode dont la partie réelle change de signe entre deux pas est
une bifurcation, raffinée par dichotomie jusqu'à `tol` (relatif) : de **Hopf**
pour une paire complexe (une oscillation naît à $\operatorname{Im}\lambda/2\pi$),
**nœud-col** pour une valeur propre réelle. Le mode de jauge est écarté.

## Quand le point d'équilibre disparaît

Si, à un pas, le load flow ne converge plus, la suite s'arrête au dernier point
convergé : `mu` s'arrête là, `arret_a` donne la valeur refusée et `resume()` le
dit (« suite arrêtée : pas de point d'équilibre à mu = … »). La suite ne cherche
pas la frontière entre les deux pas.

## Les limites des convertisseurs

Pousser la puissance finit par sortir un convertisseur de son courant maximal.
Le modèle linéarisé, calculé en mode libre, ne le voit pas : au-delà de
$I_{max}$, le convertisseur réel serait en limitation, et le spectre calculé ne
le décrit plus, même s'il paraît parfaitement stable. La continuation
**marque** donc chaque pas où un convertisseur sort de ses limites (voir
{mod}`g2egrid.noyau.limites` et la section *Hors des limites des convertisseurs* de
{doc}`analyse_modale`). Ici, on monte la consigne du GFM de 130 MVA du bus 3 de
`wscc9_gfm` :

```python
from g2egrid.reseau.tables import assembler

def gfm3_a(p_mw):
    """wscc9_gfm, le GFM du bus 3 (130 MVA) à p_mw MW."""
    t = catalogue.charger("wscc9_gfm")
    t.Y_DER = [r if r[0] != 3 else (*r[:5], p_mw * 1e6, *r[6:]) for r in t.Y_DER]
    net = assembler(t)
    net.solve_equilibrium()
    return net

c = cont.continuer(gfm3_a, [85, 120, 150, 180, 210])
print(c.resume())
print(c.hors_limite, round(c.mu_hors_limite, 2))
print(c.depassements[-1][0])
```

```text
aucune bifurcation
hors limite des convertisseurs a partir de mu = 159.346 (2 point(s) sur 5) : GFM3 (courant)
[False False False  True  True] 159.35
GFM3 (bus 3) : courant 1.5745 pu > Imax = 1.2 pu (+31.2 %, base Sn = 130 MVA)
```

Au-delà de 159,35 MW, le courant côté convertisseur du GFM3 dépasse
$I_{max} = 1{,}2$ pu ; à 210 MW, il en demande 1,57 pu. Le spectre du dernier
pas est stable ($\max \operatorname{Re}\lambda \approx -10^{-11}$, le mode de
jauge) : c'est celui d'un GFM qui saurait fournir ce courant.

| Champ | Contenu |
|---|---|
| `depassements[k]` | les {class}`~g2egrid.noyau.limites.Depassement` du pas `k` : élément, bus, grandeur (`"courant"` ou `"modulation"`), valeur, limite |
| `hors_limite` | masque des pas qui en ont au moins un |
| `mu_hors_limite` | le premier $\mu$ hors limite, raffiné par dichotomie entre le dernier pas dans les limites et le premier dehors (`mu[0]` si le départ est déjà dehors) |
| `Bifurcation.hors_limite` | l'une des bornes de la dichotomie qui encadre la traversée est hors limite |

**Marquer, pas s'arrêter.** Sans point d'équilibre, la suite s'arrête.
Au-delà de $I_{max}$, il en reste un, celui du modèle libre : il dit ce que ferait un convertisseur plus
gros, et où la commande perdrait sa stabilité s'il l'était. C'est une
information de dimensionnement, que l'arrêt cacherait. `marge` reste le
premier $\mu$ qui déstabilise, sans égard aux limites : à comparer à
`mu_hors_limite`. Un seul avertissement
{class}`~g2egrid.noyau.limites.DepassementLimite` est émis, à la fin de la suite,
au lieu d'un par pas.

## Suivre l'avancement

La boucle signale son avancement par {mod}`g2egrid.noyau.progres` ; voir
{doc}`automatiser`.

## Dans l'atelier

Onglet **Limites et impédances** : *Continuation d'un paramètre* (élément et paramètre choisis dans des
listes, échelle logarithmique, max Re du spectre et lieu des modes que le
paramètre déplace). Elle tourne en tâche de fond, avec une barre
d'avancement et un bouton *Arrêter*.
