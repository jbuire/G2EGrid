# Stabilité par les impédances : la boucle mineure

Les études de réseaux riches en convertisseurs raisonnent souvent en
**impédances** : chaque convertisseur a une admittance dq $Y_e(s)$, le reste du
réseau une impédance $Z_r(s)$ vue du bus de raccordement, et la stabilité de
l'association se juge sur le gain de boucle par le **critère de Nyquist
généralisé**. {mod}`g2egrid.analyse.impedance` le calcule à partir du même modèle d'état
que l'analyse modale, ce qui permet de **confronter** les deux lectures.

## Utilisation

```python
from g2egrid import catalogue
from g2egrid.analyse.impedance import boucle_mineure

net = catalogue.construire("smib_gfl")
b = boucle_mineure(net, "GFL2")          # le convertisseur, contre tout le reste
print(b.resume())
print(b.forme, b.applicable, b.stable, b.coherent)
```

```text
GFL2 au bus 2 (courant) : 0 tour(s), 0 mode(s) instable(s) au complet ; marge 3.68e-05 a 22882.5 Hz ; non passif sur 0.0-815.8 Hz
courant True True True
```

Aucun tour : le critère s'applique et conclut à la stabilité, comme le modèle
d'état. La marge, minuscule, désigne un mode du réseau complet presque non
amorti à 22,9 kHz, dans le filtre du convertisseur : stable, mais à surveiller
dès qu'on touche aux résistances du filtre.

| Attribut | Contenu |
|---|---|
| `f`, `L` | la grille de fréquences effectivement employée et le gain de boucle $L(j2\pi f)$, forme `(n, 2, 2)` |
| `Y_e`, `Z_r`, `Z_t` | admittance du convertisseur, impédance du reste, impédance de tout le réseau vue du bus |
| `nyquist` | `tours` (sens trigonométrique), `marge` $= \min\lvert 1 + \lambda_i(L)\rvert$, `f_marge`, `lieux` (valeurs propres de $L$), `det` |
| `instables` | modes instables du réseau complet, lus sur le modèle d'état ($Z$) |
| `poles_ouverts` | $P = Z + N$, les pôles instables de la boucle ouverte |
| `applicable`, `stable`, `coherent` | $P = 0$ ; aucun tour ; $P \ge 0$ |
| `non_passif` | les bandes de fréquence où le convertisseur n'est pas passif |

## Les deux formes du gain de boucle

$$
L(s) = Z_r(s)\,Y_e(s) \quad\text{(source de courant : GFL)}, \qquad
L(s) = Z_e(s)\,Y_r(s) \quad\text{(source de tension : GFM, SM)}.
$$

Les deux ont les mêmes zéros, ceux de $Y_r + Y_e$, c'est-à-dire les modes du
réseau complet, mais pas les mêmes pôles : on prend celle dont la **boucle ouverte
est stable**. Par défaut, `forme="courant"` pour un GFL et `"tension"` pour les
autres ; `forme=` l'impose. Avec $P$ pôles instables en boucle ouverte, le
réseau complet a $Z = P - N$ modes instables, $N$ étant le nombre de tours de
$\det(I + L(j\omega))$ autour de l'origine. Comme le modèle d'état donne $Z$
directement, `boucle_mineure` en déduit $P$ : le critère ne se lit tel quel
(stable si et seulement si $N = 0$) que lorsque `applicable` est vrai.

Les impédances viennent toutes du modèle d'état global : $Z_t$ par un courant
d'essai injecté dans la loi des nœuds du bus
({meth}`~g2egrid.reseau.network.Network.build_state_space` avec `injections=`), $Y_e$
par {func}`g2egrid.analyse.frequence.admittance_dq`, et $Z_r = (Z_t^{-1} - Y_e)^{-1}$.
L'élément de référence d'angle est refusé : il est couplé à tout le réseau par
$\theta_g$, pas seulement par sa tension de bornes.

## Un convertisseur qui déstabilise le réseau

Sur un tirage de la campagne (`etudes/campagne`), le GFL du bus 2 de WSCC9
crée une oscillation à 78 Hz. Le critère la voit :

```python
import sys
sys.path.insert(0, "etudes/campagne")
import campagne

t = campagne.construire(campagne.tirer(28, "wscc9"))
net28 = campagne.modele_complet(t)[0]
b28 = boucle_mineure(net28, "GFL2")
print(b28.resume())
```

```text
GFL2 au bus 2 (courant) : -2 tour(s), 2 mode(s) instable(s) au complet ; marge 0.128 a 954.7 Hz ; non passif sur 0.0-644.4 Hz
```

Deux tours dans le sens horaire, deux modes instables (une paire) et une boucle
ouverte stable : le convertisseur et le reste du réseau, stables chacun de leur
côté, sont instables ensemble.

## La grille de fréquences

Le nombre de tours est un comptage, et un comptage sur une grille peut se
tromper : près d'un mode peu amorti, $\det(I + L)$ tourne de $\pi$ en une
fraction de hertz. Trois précautions le rendent indépendant de la grille :

1. la grille est **resserrée** autour des pôles et zéros de $\det(I + L)$,
   connus par le modèle d'état : modes du réseau complet, de l'équipement, du
   *reste* du réseau (le bouclage de $Z_t$ par $-Y_e$) et zéros de transmission
   de $Z_t$. Sur smib_gfl, un mode du reste à 16 kHz n'est amorti que sur
   0,3 Hz ;
2. elle est **raffinée** là où la phase saute de plus de 0,3 rad d'un point au
   suivant ;
3. elle est **prolongée** vers le haut tant que $\sigma_{\max}(L)$ dépasse
   0,05 : le contour se referme en ligne droite au-delà de la dernière
   fréquence, ce qui suppose le gain de boucle négligeable.

`b.coherent` signale un $P$ négatif, signe d'un comptage à reprendre.

## Passivité

```python
from g2egrid.analyse.impedance import passivite, bandes

p = passivite(b.Y_e)                 # plus petite valeur propre de (Y + Y^H)/2
print(bandes(b.f, p < 0)[:3])        # les bandes non passives
```

Une bande où le convertisseur n'est pas passif est une bande où il peut fournir
de l'énergie à une oscillation : c'est là qu'une résonance du réseau devient
dangereuse.

## Dans l'atelier

Onglet **Limites et impédances**, carte *Boucle mineure* : lieux
caractéristiques sur tout l'axe (échelle radiale compressée autour de $-1$),
$\lvert Z\rvert$ du reste et du convertisseur, bandes non passives, verdict.
