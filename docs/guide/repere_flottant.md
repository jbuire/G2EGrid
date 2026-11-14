# Repère de référence flottant

## Le repère lié au slack

Toutes les tensions et tous les courants du réseau sont exprimés dans un
repère $dq$ global d'angle $\theta_g$. Chaque équipement qui a un angle propre
(rotor, PLL, statisme) tourne ses grandeurs de $\theta - \theta_g$, et les
éléments passifs écrivent leurs couplages $\omega L i$ avec la pulsation du
repère ({doc}`../reperes_et_angles`).

Par défaut — c'est la règle de G2Elin V1 —, $\theta_g$ est l'angle **interne
de l'équipement du bus slack**. Le repère est alors lié à cet équipement :
le déclencher, c'est emporter le repère avec lui.

## Un repère qui ne dépend de personne

```python
net.reference_flottante("coi")        # ou Tables.REFERENCE = "coi"
```

ajoute au réseau un élément {class}`~g2egrid.elements.Repere` à un seul état :

$$
\dot\theta_g = \omega_b\,(\omega_{\mathrm{in}} - 1),
$$

qui fournit exactement les deux sorties qu'on lisait sur l'équipement de
référence — un angle et une vitesse. $\omega_{\mathrm{in}}$ est au choix :

| `suivre` | $\omega_{\mathrm{in}}$ |
|---|---|
| `"coi"` (défaut) | vitesse du **centre d'inertie** des formeurs : $\sum_k w_k \omega_k / \sum_k w_k$, avec $w_k = 2H_k S_{n,k}$, ou $S_{n,k}$ pour un formeur sans inertie (GFM, bus infini) |
| `None` | vitesse nominale : repère synchrone fixe (`Tables.REFERENCE = "nominale"`) |
| `"SM2"` | vitesse d'une machine |
| `{"SM1": 2, "GFM3": 1}` | moyenne pondérée quelconque |

Les *formeurs* sont les équipements qui ont un angle et une vitesse propres —
machines, grid-forming, bus infinis ({meth}`~g2egrid.reseau.network.Network.formeurs`).
Un grid-following n'en est pas : sa PLL suit le réseau, elle ne le fait pas
tourner.

## Ce qui change, ce qui ne change pas

Un repère n'est qu'un choix de coordonnées. Le spectre avec repère flottant
est celui du repère lié au slack, à $10^{-11}$ près, **plus un mode de
jauge** : le repère a son propre angle absolu, que rien ne rappelle.

Ce qui change : aucun équipement n'est plus indispensable au modèle
dynamique. La machine du slack peut être déclenchée — par un disjoncteur des
tables ou une manœuvre en simulation — sans que le repère ne change :

```python
t2 = t.copie(REFERENCE="coi", OUVERTS=[("T_SMs", 1)])   # la machine du slack ouverte
net = assembler(t2)
net.run_loadflow()          # le load flow garde un slack : le plus gros formeur restant
ss = net.build_state_space()
```

En simulation, quand une machine suivie par le centre d'inertie est
déclenchée, elle en sort et les poids des autres sont renormalisés.

```{note}
Pourquoi suivre une vitesse plutôt que tourner à 50 Hz ? Après une perte de
production, le réseau se stabilise à une fréquence un peu différente
(statisme). Dans un repère fixe, les grandeurs $dq$ tourneraient alors
indéfiniment, sans équilibre autour duquel linéariser.
```
