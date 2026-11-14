# Groupes nucléaires : fiches de machines, lignes d'arbre, réseaux des thèses

Cette page décrit ce que le paquet tire de deux thèses sur le nucléaire et le
réseau — celle de Boudot sur l'intégration de parcs de petits réacteurs
modulaires (SMR), celle d'A.-L. Mazauric sur la flexibilité du nucléaire — et
du rapport EPRI 1013460 sur les lignes d'arbre :

- un **catalogue de machines** ({mod}`g2egrid.reseau.machines`) : 19 fiches
  d'alternateurs, converties à la demande dans les paramètres du modèle ;
- deux **lignes d'arbre de groupes nucléaires**, à ajouter à celle de l'IEEE
  First Benchmark Model, et la clé `arbre` qui les met sur une machine d'un
  réseau décrit par ses tables ;
- deux **réseaux du catalogue** : `ieee39_boudot` et `kundur_mazauric`.

L'étude `etudes/fatigue_nucleaire/` s'en sert pour chiffrer la fatigue
mécanique des lignes d'arbre.

## Le catalogue de machines

Une fiche porte les données d'une fiche constructeur (réactances en pu de la
machine, constantes de temps à circuit ouvert), l'inertie, la filière et la
source :

```python
from g2egrid.reseau.machines import fiche_machine, lister_machines

f = fiche_machine("smr_245")
f.Sn_VA, f.H, f.standard["Xd"]        # 245,5 MVA ; 4,624 s ; 1,71 pu
f.reglages(f0=50.0)                   # paramètres de Park + H, prêts pour SM_PARAMS
[m.nom for m in lister_machines("nucleaire")]
```

| Fiches | Machines | Source |
|---|---|---|
| `smr_25`, `smr_35`, `smr_51`, `smr_75`, `smr_100`, `smr_245` | turbo-alternateurs de SMR, 25 à 245,5 MVA | thèse Boudot, tableau IV.1 (d'après Anderson et Fouad) |
| `ne39_g1` à `ne39_g10` | les dix machines de l'IEEE 39 nœuds : cinq nucléaires (2, 3, 6, 8, 9), trois thermiques (4, 5, 7), une hydraulique (10), l'interconnexion (1) | thèse Boudot, tableaux A.1 à A.5 |
| `kundur_900_a`, `kundur_900_b` | machines de 900 MVA du réseau à quatre machines | thèse Mazauric, tableau III-4 |
| `nucleaire_epri_1500` | groupe nucléaire réel de 1500 MVA, 1800 tr/min | EPRI 1013460, annexe A |

Remplacer une machine d'un réseau par une autre tient en une ligne. La
« modularisation » de la thèse de Boudot — une centrale remplacée par le
nombre minimal de SMR — s'écrit par exemple :

```python
from g2egrid.catalogue import charger

t = charger("ieee39_boudot")
smr = fiche_machine("smr_245")
t.T_SMs = [(i, b, 3 * smr.Sn_VA if b == 37 else s) for i, b, s in t.T_SMs]
t.PARAMS_PAR_BUS[37] |= smr.reglages(f0=60.0)     # trois SMR à la place de la centrale 8
```

### Ce qui est traduit, et comment

Les thèses emploient des régulations normalisées que le modèle de machine du
paquet n'a pas telles quelles. Le champ `regulation` d'une fiche en donne
l'équivalent :

| Régulation de la source | Dans le modèle | Exact ? |
|---|---|---|
| AVR ST1, avance-retard $(1+sT_c)/(1+sT_b)$ | boucle de stabilisation : $T_{fd} = T_c$, $K_{fd} = (T_b - T_c)/K_a$ | oui |
| AVR EXST1, boucle $sK_f/(1+sT_f)$ | $K_{fd} = K_f$, $T_{fd} = T_f$ | oui |
| stabilisateur conventionnel | gain, passe-haut, deux avances de phase, terme à terme | oui |
| régulateur de vitesse IEEEG1 (vapeur) | premier ordre : même statisme, $T_g$ = retard moyen (1,26 s) | non : étages non distingués |
| régulateur de vitesse IEEEG3 (hydraulique) | premier ordre : même statisme, $T_g$ = 52,4 s | non : pas de réponse inverse |

La première ligne n'est pas un détail : sans la réduction du gain transitoire
de l'AVR, l'IEEE 39 nœuds paramétré par la thèse est instable. Avec elle, tous
ses modes électromécaniques sont amortis à plus de 11 %.

La machine hydraulique (10) est à pôles saillants : la thèse ne donne ni
$x'_q$ ni $T'_{q0}$. Le modèle a deux enroulements d'axe q ; le premier est
rendu inopérant ($x'_q = 0{,}98\,x_q$).

## Les lignes d'arbre

Aucune des deux thèses ne décrit la ligne d'arbre de ses machines : elles en
donnent l'inertie totale $H$, pas sa répartition ni les raideurs. Les lignes
viennent donc de la littérature ({mod}`g2egrid.elements.turbo_alternateur`) :

| Ligne | Machine | Masses | Modes de torsion [Hz] |
|---|---|---|---|
| `ARBRE_NUCLEAIRE_EPRI` (`"nucleaire_epri"`) | groupe nucléaire réel, 1500 MVA, 1800 tr/min, 60 Hz ; $H$ = 3,75 s | 9 : corps HP et trois corps BP en deux demi-corps, alternateur | 8,15 / 15,83 / 22,05 / 23,71 ; 50,6 / 52,5 / 54,5 ; 133,9 |
| `ARBRE_AM600` (`"am600"`) | projet AM600, 600 MWe, 1500 tr/min, 50 Hz ; $H$ = 4,70 s | 3 : HP, BP, alternateur | 18,87 / 28,80 |
| `ARBRE_FBM` (`"fbm"`) | IEEE First Benchmark Model, 892 MVA, 3600 tr/min (thermique) | 6 | 15,71 / 20,21 / 25,55 / 32,28 / 47,46 |

La ligne EPRI est publiée en unités physiques (kg·m², N·m/rad) ;
{func}`~g2egrid.elements.turbo_alternateur.ligne_depuis_physique` la convertit en pu. La même
conversion, appliquée aux données physiques du FBM que donne le même rapport,
redonne ses valeurs en pu de 1977 : c'est ce qui la valide. Ses quatre
premiers modes sont ceux de la ligne entière ; les trois suivants opposent les
deux flux d'un même corps basse pression.

```{admonition} Deux réserves sur l'AM600
:class: warning
L'article ne donne ni la puissance apparente ni le partage de puissance entre
corps : 720 MVA et 30 % / 70 % sont des hypothèses. Et ses propres données
donnent 18,87 et 28,80 Hz, alors qu'il annonce 18,2 et 23,9 Hz : le second
mode n'est pas retrouvé. C'est un ordre de grandeur d'une ligne courte et
raide, pas une référence.
```

### Prêter une ligne à une machine dont on ne connaît que $H$

{func}`~g2egrid.elements.turbo_alternateur.ligne_a_l_echelle` ramène une ligne à l'inertie d'une
autre machine sans changer ses fréquences de torsion ni ses déformées : toutes
les inerties, raideurs et amortissements sont multipliés par le même facteur
(et les raideurs corrigées du changement de fréquence du réseau, pour garder
les fréquences en hertz). **C'est une hypothèse de travail** : elle suppose
que la machine cible a la même répartition d'inertie et de raideur que la
machine d'origine.

Dans un réseau décrit par ses tables, la clé `arbre` des réglages d'une
machine fait ce travail :

```python
t = charger("ieee39_boudot")
t.PARAMS_PAR_BUS[38]["arbre"] = "nucleaire_epri"   # la centrale 9 reçoit la ligne à 9 masses
net = assembler(t)
```

Le `H` de la machine garde son sens de fiche constructeur — l'inertie de toute
la ligne — et la ligne nommée le répartit entre ses masses. Les couples des
tronçons deviennent des sorties de service (`SM38.C_LP3B_GEN`...) : ce sont
les grandeurs de la fatigue ({mod}`g2egrid.fatigue`). `ligne.arbres` nomme les
tronçons qui sont des accouplements entre rotors, ceux qu'il faut juger.

```{admonition} Stabilisateur et torsion
:class: warning
Un stabilisateur à entrée vitesse sans filtre de torsion peut rendre instable
le premier mode de torsion. C'est le cas sur l'IEEE 39 nœuds avec les
stabilisateurs de la thèse et la ligne EPRI (mode à 8,2 Hz, trois groupes sur
cinq). Allonger le filtre passe-bas `T_LP` du stabilisateur (0,1 s) le
stabilise, au prix d'un amortissement électromécanique moindre : voir l'étude.
```

## Les deux réseaux du catalogue

| Réseau | Contenu |
|---|---|
| `ieee39_boudot` | l'IEEE 39 nœuds (60 Hz) paramétré par la thèse de Boudot : puissances nominales, données de chaque alternateur, AVR, stabilisateurs et régulateurs de vitesse des tableaux A.1 à A.5. Topologie, charges et plan de production d'`ieee39` |
| `kundur_mazauric` | le réseau à quatre machines passé à 50 Hz par la thèse de Mazauric : GEN_A1 thermique (référence), GEN_A2 et GEN_B2 nucléaires (statisme 4 %, constante de 10 s qui tient la rampe de 5 %Pn/min), GEN_B1 hydraulique, 500 MW de renouvelable au nœud de la charge A. Plan de production du tableau III-9 |

Ils sont écrits par `developpement/reseaux_theses.py`, qui dit aussi ce qui
n'est pas repris. Pour `kundur_mazauric` : le stabilisateur STAB1 de la thèse
est réglé pour une excitatrice ESAC8B dont elle ne donne pas les paramètres ;
posé sur le régulateur de tension du paquet, il rend instable un mode
d'excitation à 11 Hz. GEN_B1 porte donc le stabilisateur par défaut du paquet,
et lui seul en a un, comme dans la thèse.

## Références

- Boudot, thèse sur l'intégration de parcs de SMR dans les réseaux
  électriques : chapitres III et IV, annexe A.
- A.-L. Mazauric, thèse sur la flexibilité du nucléaire et la stabilité du
  réseau : chapitre III, annexe A.
- EPRI, *Torsional Interaction Between Electrical Network Phenomena and
  Turbine-Generator Shafts: Plant Vulnerability*, rapport 1013460, 2006.
- T. Mudau, R. M. Field, « Rotordynamic Analysis of the AM600
  Turbine-Generator Shaftline », Trans. Korean Nuclear Society, automne 2018.
