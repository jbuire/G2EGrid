# Choisir les dynamiques

Un modèle d'équipement porte toutes ses dynamiques : filtre, boucles de
régulation, bus continu, statisme du second ordre... Pour une étude donnée,
certaines sont inutiles, ou gênantes : trop rapides pour la question posée,
ou tenues pour parfaitement régulées. g2egrid permet de les **supprimer** sans
écrire un nouveau modèle.

## Deux façons de dire « cette dynamique vaut zéro »

Pour un état $x_i$ d'équation $\dot x_i = f_i(x, z, u)$ :

`statique` — $\dot x_i = 0$
: La dynamique est jugée **infiniment rapide** : l'état suit instantanément son
  équilibre. Il devient une variable algébrique et son équation une
  contrainte $0 = f_i(x, z, u)$. C'est la perturbation singulière classique.

`fige` — $x_i = x_{i,0}$
: La dynamique est jugée **infiniment lente**, ou parfaitement régulée : l'état
  reste à sa valeur d'équilibre (ou à une valeur imposée). Il devient un
  **paramètre** du modèle, et son équation disparaît.

Les deux ne donnent pas le même modèle. Sur le statisme du second ordre d'un
grid-forming, $\dot{\Delta\omega} = \omega_c\,(m_p (p^\star - p_m) - \Delta\omega)$ :

- `Dw` **statique** donne $\Delta\omega = m_p (p^\star - p_m)$, le statisme du
  premier ordre ;
- `Dw` **figé** donne $\Delta\omega = 0$ : plus de statisme du tout.

## En code

```python
from g2egrid.elements import GFM4

gfm = GFM4("GFM1", 1, Sn=100e6, Sb=100e6, wb=net.wb)
gfm.reduire(statique=["Dw"], fige=["vdc", "idc"])       # rend l'élément lui-même
net.add(gfm)

gfm.dynamiques()      # {'isd': 'dynamique', ..., 'vdc': 'fige', 'Dw': 'statique', ...}
```

- `fige` accepte aussi un dict `{état: valeur}` : `fige={"vdc": 1.0}` impose la
  valeur au lieu de celle de l'équilibre ;
- les réductions **se cumulent** d'un appel à l'autre ;
- un nom d'état inconnu lève une erreur, sauf avec `strict=False` ;
- l'élément change de classe — une sous-classe engendrée et mise en cache, au
  modèle compilé propre —, `isinstance(gfm, GFM4)` reste vrai, et tout le
  reste (équilibre, linéarisation, simulation, sensibilités) le traite comme
  n'importe quel élément. Deux éléments réduits de la même façon partagent le
  même code compilé.

Un état figé devient un paramètre comme un autre : sa sensibilité se lit avec
les autres ({doc}`sensibilites`).

## Les variantes deviennent des sélections

GFM1 à GFM3 ne sont que des sélections de dynamiques de GFM4. Le test
`tests/test_dynamiques.py` le vérifie sur le spectre complet, à $10^{-8}$ près :

| Variante | Sélection sur GFM4 |
|---|---|
| GFM1 (statisme 1{sup}`er` ordre, bus continu idéal) | `statique=["Dw"], fige=["vdc", "idc"]` |
| GFM2 (statisme 2{sup}`e` ordre, bus continu idéal) | `fige=["vdc", "idc"]` |
| GFM3 (statisme 1{sup}`er` ordre, bus continu dynamique) | `statique=["Dw"]` |

Les classes GFM1…GFM3 restent disponibles : leur ordre d'états est celui des
classes MATLAB, ce que demandent les comparaisons terme à terme.

D'autres sélections utiles sur la machine synchrone :

| But | Sélection |
|---|---|
| puissance mécanique constante (sans régulateur de vitesse) | `fige=["Pm"]` |
| tension d'excitation constante (sans AVR) | `fige=["e1", "e2", "efd", "e3"]` |
| sans PSS | `fige=["Dw1", "v1", "v2", "vpss"]` avec `K_PSS = 0` |

```{warning}
Figer n'est pas débrancher : figer un état interne d'une chaîne (le filtre
d'entrée du PSS, `Dw1`, seul) laisse le reste de la chaîne en place, qui voit
alors une entrée constante là où elle attendait un signal filtré — et peut
déstabiliser le réseau. Figer la chaîne entière, ou annuler son gain.
```

## Dans les tables et l'atelier

La clé `dynamiques` d'un réglage fait la même sélection, pour une **famille**
(`SM_PARAMS`, `GFM_PARAMS`, `GFL_PARAMS`) ou pour l'équipement d'un **bus**
(`PARAMS_PAR_BUS`). Le bus complète la famille, clé par clé :

```python
t2 = t.copie(
    GFM_PARAMS={"dynamiques": {"Dw": "statique", "vdc": "fige", "idc": "fige"}},
    PARAMS_PAR_BUS={3: {"dynamiques": {"Dw": "dynamique"}}},   # le GFM du bus 3 garde Dw
)
```

| Valeur | Effet |
|---|---|
| `"statique"` | $\dot x = 0$ |
| `"fige"` | $x = x_0$ |
| un nombre | $x$ figé à cette valeur |
| `"dynamique"` | rien — rend à un bus une dynamique que sa famille supprime |

Une sélection de famille vaut pour des variantes qui n'ont pas toutes les
mêmes états (`Dw` n'existe pas dans GFM1) : chaque équipement ignore les noms
qu'il n'a pas. Mais un nom qu'**aucun** équipement du réseau ne porte est une
faute de frappe, et l'assemblage la signale. Dans l'atelier, la sélection se
saisit dans les « Réglages fins » d'un équipement.

```{admonition} Limite
:class: note

Rendre statique un **intégrateur** de correcteur ($\dot M = K_i\,e$) donne la
contrainte $0 = e$, qui ne porte pas sur $M$ : le système est alors d'indice
supérieur. La résolution globale sait réduire l'indice quand la contrainte
porte sur des états ({doc}`resolution_globale`), mais ce n'est en général pas
ce qu'on veut : figer l'intégrateur est plus sûr.
```
