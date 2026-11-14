# Défauts et protections du réseau continu

Un défaut continu ne ressemble à aucun défaut alternatif : les condensateurs
des stations et des câbles se déchargent en quelques centaines de
microsecondes, la tension continue s'effondre, et rien ne passe par zéro
pour aider à couper. Ce qui décide de la suite, ce sont les **protections** :
blocage des convertisseurs en moins d'une milliseconde, ouverture des
disjoncteurs à courant continu en quelques millisecondes.

## Les manœuvres

| Manœuvre | Effet |
|---|---|
| `defaut_dc` | court-circuit pôle à pôle au nœud continu `cible` (`"bus221"`), branche R-L ({class}`~g2egrid.elements.DefautDC`), `R = 1e-4`, `X = 1e-3` pu par défaut |
| `elimination_dc` | retire le défaut |
| `ouverture` d'un câble (`CDC1_201_221`) | le disjoncteur à courant continu l'ouvre |
| `blocage` | le convertisseur `cible` se bloque : sa limitation de courant passe à `Imax ≈ 0` et y reste, la boucle de courant éteint le courant en quelques millisecondes, et il ne prend plus rien à son nœud continu |

Le blocage ne retire pas l'élément : couper en zéro seconde le courant
inductif du transformateur ferait une surtension qui effondrerait le poste.
Le modèle est celui d'un convertisseur capable de bloquer un défaut continu
(pont complet, ou isolé par ses disjoncteurs), pas les diodes de roue libre
d'un demi-pont.

## Les protections

Une {class}`~g2egrid.temporel.Protection` est une mesure, un seuil, un retard
et les manœuvres qu'elle déclenche. Pendant la simulation, son seuil est une
**garde** de l'intégrateur, comme celles des limiteurs : l'intégration
s'arrête à l'instant exact du franchissement, les manœuvres sont posées
`retard` plus tard, et l'intégration reprend.

| Protection | Mesure | Action par défaut |
|---|---|---|
| `surintensite_dc` | `\|i\|` d'un câble continu | ouverture du câble |
| `sous_tension_dc` | tension d'un nœud continu | blocage des convertisseurs du nœud |
| `surintensite_convertisseur` | `\|is\|` d'un convertisseur | blocage |
| `surtension_ac` | `\|ve\|` aux bornes d'un convertisseur | blocage (parc dont l'île a perdu son formeur) |

{func}`~g2egrid.temporel.protections_dc` construit un plan type :

- un disjoncteur par câble, à 3 fois son courant d'équilibre (au moins le
  courant moyen des câbles), 3 ms ;
- chaque station se bloque sous 0,8 pu continu ou au-delà de 2 pu de
  courant, 0,2 ms ;
- chaque source sans port continu (éolien, photovoltaïque) décroche au-delà
  de 1,3 pu alternatif, 1 ms.

```python
from g2egrid import catalogue
from g2egrid.temporel import ModeleNonLineaire, Manoeuvre, protections_dc

net = catalogue.construire("cigre_b4_dcs1")
m = ModeleNonLineaire(net)
plan = protections_dc(net)
for p in plan:
    print(f"{p.nom:40s} seuil {p.seuil:.3g}  retard {1e3 * p.retard:.1f} ms")
sc = m.scenario(Manoeuvre(0.02, "defaut_dc", "bus221"))
tr = m.simuler((0, 0.06), uext_fn=sc, methode="BDF2", pas=2.5e-5, protections=plan)
for t, nom, mans in tr.declenchements:
    print(f"{1e3 * t:7.3f} ms  {nom:28s} -> {mans}")
```

```text
surintensite_dc:CDC1_201_221             seuil 11.8  retard 3.0 ms
sous_tension_dc:bus201                   seuil 0.8  retard 0.2 ms
sous_tension_dc:bus221                   seuil 0.8  retard 0.2 ms
surintensite_convertisseur:GFL101        seuil 2  retard 0.2 ms
surintensite_convertisseur:GFL121        seuil 2  retard 0.2 ms
 20.176 ms  sous_tension_dc:bus221       -> [('blocage', 'GFL121')]
 22.763 ms  sous_tension_dc:bus201       -> [('blocage', 'GFL101')]
```

Défaut franc au nœud C1 de la liaison DCS1 : Cm-C1 se bloque 0,18 ms après le
défaut, Cm-A1 quand l'onde a parcouru les 200 km de câble. Le courant de
décharge reste sous le seuil du disjoncteur, qui ne s'ouvre pas.

Sur le CIGRE B4 maillé, un défaut au nœud B4 fait ouvrir les disjoncteurs des
deux câbles qui le portent, bloque les stations voisines, puis, de proche en
proche, celles dont le nœud continu s'effondre : le plan type est un point de
départ, pas un plan sélectif. Le blocage seul de Cb-B2 (1,5 GW), lui, ne
déclenche rien d'autre que la station elle-même.

## Les charges sous 0,7 pu

Une puissance ou un courant constants n'ont pas de sens à tension nulle : le
courant de consigne y perd sa direction. Sous `Vbt = 0,7` pu, une charge
`PLoad` ou `ILoad` devient l'impédance qui appelle la même puissance à
`Vbt`, et revient à son modèle au-dessus de `Vbt` (+2 %) : un mode discret
(`mode_bt`), surveillé comme les limiteurs. À drapeau nul, les équations
sont celles d'avant.

## Un défaut alternatif aux bornes d'une station

Un défaut **alternatif** peut vider un réseau continu. Sur `cigre_mv_dc`, un
défaut au bus 15 prive de réseau la station GFL qui tient la tension de la
liaison MVDC. Pendant ce temps, la station GFM du bus 16 soutient son propre
réseau, lui aussi déprimé, en puisant dans la liaison, que plus rien ne
recharge. Pour le continu, c'est une charge à puissance constante sur un
condensateur : sans borne sur la modulation $e/v_{dc}$, $v_{dc}$ s'effondre
en un temps fini et la simulation n'a plus de solution (le Newton lâchait à
83 ms, $v_{dc}$ à 0,33 pu).

C'est la **saturation de tension** qui l'arrête : $\lvert m\rvert \le M_{max}$,
l'onduleur ne sait plus fabriquer sa tension alternative quand $v_{dc}$ tombe
vers $\lvert e\rvert/M_{max}$, et le transit s'arrête de lui-même. Elle est
donc **active par défaut** ($M_{max} = 1{,}15$) sur tout GFM dont la tension
continue est une variable : `T_GFM3s` et `T_GFM4s` (y compris la commande par
appariement), monolithiques ou à port continu. Sur le défaut au bus 15,
$v_{dc}$ descend à 0,39 pu et revient à 1 pu en 0,3 s environ. Un réglage de
famille ou de bus l'emporte, comme pour tout paramètre.

Pas sur les GFL : un GFL ne tire du continu que ce que ses consignes de
courant demandent, et la station qui tient $v_{dc}$ le réduit d'elle-même
quand $v_{dc}$ baisse. Saturée, sa modulation lâche au contraire la commande
du courant dès que $v_{dc}$ passe sous $\lvert e\rvert/M_{max}$, ce qui arrive
sur tout défaut continu ou blocage de la station d'en face : sur
`cigre_b4_dcs1`, la liaison oscille après un blocage, et le calcul n'aboutit
plus sur un défaut franc pôle à pôle. Comme sur bus continu idéal, elle reste
à régler (`Mmax = 0` par défaut).

La boucle de tension continue des GFL2 est proportionnelle par défaut
(`Kidc = 0`, comme `GFL2.m`). Son état `M_d` est alors le courant du point de
fonctionnement : la limitation de courant ne le rétro-calcule pas, sans quoi
$v_{dc}$ garderait l'écart accumulé pendant la limitation (4 % sur ce défaut,
quand la station réglée avec la saturation de tension passe en limitation à
l'élimination).

## Dans l'atelier

Onglet **Temporel** : les manœuvres `defaut dc`, `elimination dc` et
`blocage` dans la frise (nœuds continus et convertisseurs dans les listes,
câbles continus parmi les éléments à ouvrir) ; la case **protections DC**
arme le plan type. Les déclenchements s'affichent sous les courbes, avec
leur instant et ce qu'ils ont manœuvré.

## Ouvrir un câble ou perdre une station en régime permanent

Un câble continu s'ouvre comme une ligne : `OUVERTS=[("T_LinesDC", 2)]`
(nom assemblé `CDC2_202_211`). {func}`g2egrid.reseau.disjoncteurs.en_service` rend
alors, dans `Rapport.iles_dc`, les îles continues (nœuds, convertisseurs,
références de tension) et **refuse** une île qui garde des convertisseurs sans
rien pour tenir sa tension. Quand une île alternative meurt faute de formeur,
les nœuds continus devenus orphelins sont retirés avec leurs shunts et leurs
câbles pendants ; un nœud de type 1 qui perd son seul tenant passe en type 2.

Perdre une station n'est pas anodin pour le spectre : sur `hvdc2t` et
`mtdc3t`, le convertisseur qui reste en régulation de Vdc reprend la puissance
et le spectre devient instable (de +0,7 à +10 s⁻¹). Ce n'est pas un artefact
d'assemblage : le même réseau intact, avec la puissance de la station forcée à
zéro (`hvdc2t`), l'est aussi. `Rapport.messages` nomme la station perdue, ce
qu'elle tenait et ce qui reste.
