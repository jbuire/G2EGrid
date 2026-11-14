# Sensibilités paramétriques

La table des modes dit **qui** oscille ; la sensibilité dit **sur quoi agir**.
g2egrid en offre deux niveaux :

| Fonction | Rend | Désigne |
|---|---|---|
| {func}`g2egrid.analyse.modal.sensibilite` | $\partial\lambda / \partial A_{ij}$ | une **équation** : l'effet de l'état $j$ sur la dérivée de l'état $i$ |
| {func}`g2egrid.analyse.sensibilites.sensibilites` | $\partial\lambda / \partial p$ | un **paramètre** d'équipement : un gain, une inertie, une réactance |

Un coefficient de $A$ mélange plusieurs paramètres, et le complément de Schur
qui élimine les variables algébriques les étale sur toute la matrice. Le
second niveau remonte jusqu'au bouton à tourner.

## Utilisation

```python
from g2egrid import catalogue, modal
from g2egrid.analyse.sensibilites import sensibilites, tableau

net = catalogue.construire("wscc9")
ss = net.build_state_space()
m = modal.modes(ss, f_max=5, oscillants=True)[0]      # le mode le moins amorti

sens = sensibilites(net, m)       # tous les paramètres de tous les équipements
print(tableau(sens, n=6))
```

```text
element.parametre                 valeur       p.dRe   p.df [Hz]  p.dzeta [%]
SM3.H                                  5      0.3936     -0.8081       0.3546
SM3.wb                               377      0.1664      0.4478       -2.947
SM3.mp                               0.5      0.6536      -0.394       -3.129
SM3.Tg                               0.2       1.426     -0.3226       -8.933
T3_3_9.n                               1      0.3706     -0.3075       -1.444
SM2.H                                  5     -0.3014     -0.2994         3.33
```

Le mode peut être désigné par son indice dans `eig(A)`, par une valeur propre
approchée (la plus proche est retenue) ou par un {class}`~g2egrid.analyse.modal.Mode`.
`elements=[...]` et `parametres=[...]` restreignent le calcul.

Chaque {class}`~g2egrid.analyse.sensibilites.Sensibilite` donne :

- `dlambda` : $\partial\lambda/\partial p$ ;
- `normalisee` : $p\,\partial\lambda/\partial p$, le déplacement du mode pour une
  variation **relative** du paramètre ($\Delta p / p = 1$). C'est elle qui
  classe les résultats : elle rend comparables un gain de $10^3$ et une
  constante de temps de $10^{-3}$ ;
- `df` : $p\,\partial f/\partial p$ en Hz, et `dzeta` : $p\,\partial\zeta/\partial p$
  en points d'amortissement.

Dans l'atelier, le détail d'un mode (onglet *Modes*, clic sur un mode) affiche
ce tableau à côté de la sensibilité aux coefficients de $A$.

## Principe

On ne dérive pas $A$ : on dérive les **équations**. Le réseau est un DAE
linéaire global ; en gardant les entrées $\bm U$ des équipements comme
inconnues, il s'écrit sous forme d'un faisceau $J - \lambda E$ :

$$
\begin{aligned}
\dot{\bm x} &= F_x\,\bm x + F_z\,\bm z + F_u\,\bm U \\
\bm 0 &= G_x\,\bm x + G_z\,\bm z + G_u\,\bm U \\
\bm 0 &= R_y\,(H_x\,\bm x + H_z\,\bm z + H_u\,\bm U) - \bm U
\end{aligned}
\qquad E = \mathrm{diag}(I, 0, 0)
$$

dont les valeurs propres finies sont celles de $A$. Un paramètre d'un
équipement ne touche que les neuf blocs **de cet équipement**. Pour une valeur
propre simple de vecteurs propres à droite $\bm v$ et à gauche $\bm w$ du
faisceau :

$$
\frac{\partial\lambda}{\partial p} = \frac{\bm w^\top \, \dfrac{\partial J}{\partial p} \, \bm v}{\bm w^\top E\, \bm v}.
$$

$\bm v$ et $\bm w$ se déduisent de ceux de $A$ en résolvant les équations
algébriques, une fois pour tout le réseau. $\partial J / \partial p$ est la
dérivée **symbolique** du jacobien de l'élément
({meth}`~g2egrid.noyau.symbolic.CompiledModel.derivees_parametres`), compilée une fois
par classe. La sensibilité à tous les paramètres de tous les équipements coûte
donc une factorisation, puis un produit scalaire creux par équipement.

## Point de fonctionnement fixe ou non

La méthode analytique tient le **point de fonctionnement fixe** : elle mesure
l'effet du paramètre sur la dynamique autour du même équilibre. C'est la bonne
question pour un gain de régulation, qui ne déplace pas l'équilibre.

Pour un paramètre qui le déplace (une réactance de ligne, une consigne),
`methode="differences", equilibre=True` refait load flow et équilibre à chaque
perturbation et donne la sensibilité **totale** :

```python
sensibilites(net, m, methode="differences", equilibre=True, elements=["L1_4_5"])
```

Un paramètre que l'initialisation recalcule (`Kpdc` auto-réglé, impédance d'une
charge à puissance constante) y sort nul, puisqu'il est écrasé.

Les **drapeaux des modes discrets** (`mode_lim`, `mode_vsat` des limiteurs de
courant et de tension, `mode_e2` de la butée d'excitatrice) sont eux aussi
tenus à leur valeur, 0 en mode libre. Les dérivées symboliques sont prises
drapeaux remplacés par leur valeur : les `Min` des limiteurs disparaissent, et
celles d'un GFM, dont les deux limiteurs s'empilent autour du seuil de
l'impédance virtuelle, restent assez courtes pour se calculer en quelques
secondes au lieu de ne plus rendre. La sensibilité à un drapeau reste exacte :
elle se dérive les autres figés, puis à sa valeur. Sous la limite, celle de
`mode_lim` est nulle, puisque la consigne n'est pas réduite.

La méthode analytique est vérifiée contre les différences finies, à point de
fonctionnement fixe, à $10^{-6}$ près en relatif (`tests/test_dynamiques.py`).
