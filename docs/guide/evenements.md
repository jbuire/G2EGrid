# Après un événement : tensions, puissances, relais et cascades

Que devient le réseau après un déclenchement, un défaut éliminé, un échelon
de consigne ? L'onglet **Événements** de l'atelier répond par une simulation
non linéaire de tout le réseau, puis, pour **chaque bus** et **chaque
unité**, l'état avant l'événement, les extrêmes du transitoire et le régime
atteint. Les relais de surtension des convertisseurs peuvent y être armés,
et la cascade de surtension peut aussi se dérouler en régimes permanents
successifs, sans dynamique.

## Les événements

| Genre | Effet | Valeurs |
|---|---|---|
| `ouverture` | le disjoncteur de l'élément s'ouvre (ligne, transformateur, unité, charge) | — |
| `blocage` | le convertisseur se bloque | — |
| `defaut` / `elimination` | court-circuit R + jX au bus, puis son retrait | R, X [pu] |
| `charge` | la charge varie de dP, dQ (relatif) | dP, dQ |
| `saut_de_phase` | la tension du bus tourne de l'angle donné | angle [deg] |
| `echelon` / `rampe` | une consigne (entrée du modèle) saute, ou monte en rampe | amplitude [pu], durée [s] |

Choisir plusieurs éléments à la fois (Ctrl-clic dans la liste, avec un filtre
sur le nom) donne un événement par élément, au même instant : c'est ainsi
qu'on fait déclencher simultanément plusieurs GFL.

## Ce que rend la simulation

{func}`g2egrid_api.analyse.apres_evenement` simule en BDF2. Le pas vaut 1 ms
par défaut : on regarde des secondes de régulation. Un défaut proche
demande un pas plus fin. La fonction rend :

- **par bus** : `V0`, `Vmin` et `Vmax` avec leurs instants, `Vfin` (moyenne
  sur la dernière fenêtre, 0,2 s par défaut), `dV`, le temps passé hors des
  bornes `[v_bas, v_haut]`, et `V_lf`, la tension du **load flow** après les
  mêmes ouvertures. C'est le régime vers lequel la simulation doit tendre,
  aux régulations de fréquence près. `V_lf` n'existe que si le scénario se
  ramène à des ouvertures et des blocages : une consigne ou un échelon de
  charge n'ont pas d'équivalent dans les tables ;
- **par unité** (machine, convertisseur, charge) : P et Q avant et après, et
  l'instant où elle est sortie ;
- **la fréquence** : l'enveloppe (min, max, moyenne) des fréquences filtrées
  des bus des machines et des GFM. À défaut, elle porte sur tous les bus,
  sauf ceux d'une unité sortie ou d'un défaut : à l'instant de la manœuvre,
  l'angle de ces bus saute, et sa dérivée n'est pas une fréquence ;
- **les relais** : démarrages, retombées, déclenchements.

```python
from g2egrid.reseau.catalogue import charger
from g2egrid_api import analyse

# un GFL derriere un bus infini : la tension amont monte de 0,12 pu a 0,1 s
t = charger("smib_gfl")
r = analyse.apres_evenement(
    t,
    [{"genre": "echelon", "instant": 0.1, "entree": "IB1.E_ref", "amplitude": 0.12}],
    duree=1.0,
    relais={"v_max": 1.10, "retard": 0.3, "retombee": 0.98, "familles": ["GFL"]},
)
print([(d["protection"], round(d["t"], 2)) for d in r["declenchements"]])
# [('surtension_ac:GFL2', 0.4)]
bus2 = {b["bus"]: b for b in r["bus"]}[2]
print(round(bus2["Vmax"], 3), round(bus2["t_hors"], 1))  # 1.122 0.9
```

Sur le cas d'étude (plus bas), la même fonction suit une cascade : trois GFL
ouverts à 0,5 s, relais à 1,10 pu, 3 s et retombée 0,999. GFL314, 313 et 311
déclenchent à 3,50, 6,51 et 9,51 s.

## Relais de surtension des convertisseurs

{func}`g2egrid.temporel.protections_ac` arme un relais `surtension_ac` sur
chaque convertisseur sans port continu des familles demandées (GFL, GFM). La
mesure est `|ve|`, la tension aux bornes du convertisseur. Le seuil et la
temporisation se règlent. L'action est l'ouverture ou le blocage.

Le **rapport de retombée** (`retombee`) fait d'un relais un vrai relais à
temps constant. Le franchissement du seuil ne fait que démarrer la
temporisation. Si la mesure repasse sous `retombee · seuil` avant le terme,
le relais retombe et repart de zéro. Il ne déclenche que si la surtension a
tenu toute la temporisation. Sans retombée (`None`), le relais est verrouillé
dès le premier dépassement, si bref soit-il.

Le réglage compte quand la tension de service est près du seuil. Avec une
retombée de 0,98 et un seuil de 1,10 pu, le relais ne retombe que sous
1,078 pu. Si le réseau vit à 1,095 pu, il ne retombe jamais, et le moindre
pic le fait déclencher 3 s plus tard. Il faut alors une retombée proche de 1.

## La cascade quasi statique

{func}`g2egrid.scenarios.cascade.cascade_surtension` (bouton **Cascade quasi
statique**) part des unités ouvertes par le scénario. Elle fait déclencher
celles dont le bus dépasse le seuil, refait un load flow, et recommence
jusqu'à ce que plus rien ne dépasse. Deux ordres :

- `"vague"` : toutes les unités au-dessus du seuil déclenchent ensemble ;
- `"une_a_une"` : la plus haute d'abord, puis nouveau load flow.

Ils bornent ce que font des relais de même réglage. Des relais **à retombée**
suivent le un-à-un : le déclenchement d'une unité fait plonger un instant la
tension des voisines, leurs relais retombent, et seul le plus haut tient sa
temporisation. La cascade peut dépendre fortement de l'ordre : sur le cas
d'étude, 21 GFL déclenchent en vagues, et 4 un à un, comme en simulation.

```python
from g2egrid.scenarios.cascade import cascade_surtension

t = charger("wscc3sm_cigre3_surtension")
c = cascade_surtension(t, ["GFL303", "GFL304", "GFL305"], v_max=1.10, ordre="une_a_une")
c.issue, c.declenchees          # ('arretee', ['GFL307', 'GFL311', 'GFL313', 'GFL314'])
[e.ouvertes for e in c.etapes]  # un régime permanent par étape
```

## Le cas d'étude

`etudes/surtension_cascade` (réseau `wscc3sm_cigre3_surtension` au
catalogue) : WSCC à trois machines, trois distributions CIGRE MV, 33 GFL à
tan φ = −0,3, régleurs à 1,02 pu, câbles 110 kV. Il montre trois choses :

- la déconnexion de GFL qui absorbent du réactif fait monter la tension
  partout où le chemin vers la source a un X/R supérieur à 1/|tan φ| ;
- la cascade bascule sur quelques dixièmes de pour cent de marge ;
- un transport à X/R élevé concentre la hausse dans la zone qui perd son
  réactif.
