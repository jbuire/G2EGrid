# Repères, angles, et ce que « modèle local » veut dire

Ce document explique la seule différence de fond qui reste entre les deux
implémentations : **le choix du repère dq dans lequel les équipements
s'échangent leurs signaux**. C'est aussi ce qui décide si les modèles
d'équipement sont *locaux* ou non.

## 1. Le problème

Chaque équipement a un repère naturel qui lui est propre :

| Équipement | Repère naturel | Angle |
|---|---|---|
| Machine synchrone | celui du rotor | `theta` (angle rotor) |
| GFM | celui imposé par le statisme | `theta` |
| GFL | celui de la PLL | `theta_pll` |
| Ligne, transfo, charge, nœud | aucun (passifs) | — |

Écrire les équations d'un équipement dans son repère propre est ce qui rend le
modèle simple : `veq_ref = 0` pour un GFM, `veq = 0` à l'équilibre pour une
PLL, l'axe q aligné sur la f.e.m. pour une machine.

Mais deux équipements raccordés au même nœud doivent s'échanger une tension et
un courant **exprimés dans le même repère**. Il faut donc un repère commun, et
chaque équipement doit tourner ses grandeurs entre son repère propre et ce
repère commun. Toute la question est : **qui définit le repère commun ?**

## 2. Convention MATLAB actuelle : le repère est celui du slack

Un élément est désigné `isSlack`. Son angle devient la référence de tout le
réseau :

```
      slack (SM ou GFM)                        élément non-slack
  ┌──────────────────────┐                ┌──────────────────────────┐
  │ rot = theta_rotor    │   theta_slack  │ rot = theta_g - theta    │
  │ (CONSTANTE)          │ ─────────────► │ (theta_g = angle slack)  │
  │ sorties : i, theta,w │   w_slack      │ sorties : i              │
  └──────────────────────┘ ─────────────► └──────────────────────────┘
```

Concrètement, dans `SM.m` :

```matlab
if obj.isSlack
    ugVec = [vgd_g  vgq_g];
    rot = obj.theta_rotor;      % CONSTANTE : le slack définit le repère
else
    ugVec = [theta_g  vgd_g  vgq_g];
    rot = theta_g - theta;      % l'angle du slack est une ENTRÉE
end
```

et dans `Network.buildGlobalModel`, le routage par nom envoie la sortie
`theta_slack` du slack vers l'entrée `theta_slack` de tous les autres.

**Conséquences.**

1. **Les modèles ne sont pas locaux.** Un élément non-slack a une entrée
   supplémentaire (`theta_g`) qui vient d'un *autre* élément. Ses équations ne
   se referment pas sur ses propres variables.
2. **Deux variantes par modèle.** `SM.m` et `GFM1..4.m` contiennent chacun un
   `if obj.isSlack ... else ... end` qui change le nombre d'entrées, le nombre
   de sorties, les équations de rotation *et* le jeu de contraintes
   d'équilibre. C'est deux modèles à maintenir par équipement.
3. **Le modèle dépend de la topologie.** Changer de slack, ou retirer le slack
   d'une étude de sous-réseau, oblige à reconstruire les modèles de tous les
   équipements — pas seulement le câblage.
4. **Un GFL ne peut jamais être slack** (il n'a pas d'angle propre au sens
   « référence »), donc tout réseau 100 % GFL est hors du domaine du code.
5. Le repère commun tourne à la vitesse **du rotor du slack**, donc à une
   vitesse variable pendant un transitoire. C'est licite (repère machine de
   Kundur) mais il faut garder en tête que « l'angle » affiché pour un élément
   est mesuré par rapport à cette référence mobile.

## 3. Convention du portage Python : repère synchrone commun

Le repère commun est **fixe et défini a priori** : il tourne à `w_ref` (1 pu,
soit `2.pi.f0`). Il n'appartient à personne.

```
     élément A                repère commun (w_ref)                élément B
  ┌────────────────┐        vd, vq  ▲   ▲  vd, vq        ┌────────────────┐
  │ rot = -theta_A │ ───────────────┘   └─────────────── │ rot = -theta_B │
  │ theta_A : ÉTAT │        id, iq                id, iq │ theta_B : ÉTAT │
  └────────────────┘                                     └────────────────┘
```

Chaque équipement porte **son propre angle absolu** comme état et fait sa
rotation lui-même :

```python
ct, st = sp.cos(s.theta), sp.sin(s.theta)
ved =  s.vd_t * ct + s.vq_t * st      # repère commun -> repère propre
veq = -s.vd_t * st + s.vq_t * ct
...
iD = s.kS * (igd * ct - igq * st)     # repère propre -> repère commun
iQ = s.kS * (igd * st + igq * ct)
```

Et l'angle s'intègre sur l'**écart** de vitesse au repère :

```
dtheta/dt = wb . (w_element - w_ref)
```

**C'est la définition d'un modèle local :** les équations d'un équipement ne
font intervenir que ses propres états, ses propres paramètres, les tensions de
ses ports et le scalaire `w_ref` qui définit le repère. Aucun signal ne vient
d'un autre équipement.

Conséquences pratiques :

1. **Un seul modèle par équipement**, pas de variante `isSlack`.
2. **Aucun élément n'est privilégié.** Un réseau 100 % GFM, 100 % GFL, avec ou
   sans bus infini, se traite avec exactement les mêmes modèles.
3. **Composable** : on peut extraire un sous-réseau, remplacer une source,
   assembler deux réseaux — les modèles d'équipement ne changent pas.
4. Les angles affichés sont **absolus** (par rapport au repère synchrone), donc
   directement comparables entre eux et au load flow.

## 4. Les deux conventions donnent la même physique

Le passage de l'une à l'autre est une rotation globale d'angle
`theta_slack(t)`. Formellement, si `x_MATLAB` et `x_Python` désignent les mêmes
grandeurs dans les deux repères :

```
x_Python = R(+theta_slack) . x_MATLAB          pour les grandeurs dq
theta_i,Python = theta_i,MATLAB + theta_slack  pour les angles
```

Les deux systèmes sont donc reliés par un changement de variables régulier :
**les valeurs propres sont les mêmes**. Ce qui change à la lecture :

| | MATLAB | Python |
|---|---|---|
| Angle d'un élément | relatif au slack | absolu (repère synchrone) |
| Angle du slack | constant par construction | état comme les autres |
| Entrées d'un élément non-slack | `theta_g`, `vgd`, `vgq` | `vd_t`, `vq_t` |
| Mode de jauge `λ = 0` | présent | présent |

### Sur le mode de jauge

Il est **présent des deux côtés**, et c'est normal : aucun de ces réseaux ne
possède de référence d'angle absolue. Décaler tous les angles d'une constante
laisse la physique inchangée.

- Côté Python, le vecteur propre associé est « tous les `theta` décalés de la
  même quantité », ce qui se lit directement dans les facteurs de
  participation.
- Côté MATLAB, la même symétrie existe : le slack n'utilise pas son propre état
  `theta` dans ses équations (il utilise la constante `theta_rotor`), et les
  autres n'utilisent que la différence `theta_g - theta`. Décaler tous les
  `theta` ne change donc rien non plus. Le mode `λ = 0` est simplement porté
  par d'autres états.

Sur le SMIB, ce mode est porté à 100 % par `IB.theta` côté Python. Ce n'est pas
une instabilité : c'est l'information « ce réseau n'a pas de référence d'angle
absolue ».

## 5. Ce qu'il faudrait pour rendre les modèles MATLAB locaux

Le changement est mécanique, mais il touche tous les modèles à angle et
l'assembleur :

1. **`SM.m`, `GFM1.m` … `GFM4.m`** : supprimer les branches `if obj.isSlack`.
   Toujours utiliser `rot = -theta` (avec `theta` état), et remplacer
   `dtheta = wb*w` (angle absolu, exclu de l'équilibre) par
   `dtheta = wb*(w - w_ref)`, qui s'annule à l'équilibre — la ligne rentre
   alors dans le système résolu par `fsolve` et il n'y a plus de contrainte
   ad hoc à écrire.
2. **Interfaces** : supprimer l'entrée `theta_g` et les sorties
   `theta_slack` / `w_slack`. `w_ref` reste une entrée externe commune
   (constante à 1 pu), déjà présente partout.
3. **`Network.buildGlobalModel`** : supprimer le routage `THETA`/`W` de
   `parsePort` et le repérage `slackThetaOut` / `slackWout`. Le câblage se
   réduit à « tension du nœud → entrées de tension » et « courants → loi des
   nœuds », c'est-à-dire ce que fait déjà le reste de la fonction.
4. **Contraintes d'équilibre** : pour chaque élément, remplacer les
   contraintes « angle = référence » par la condition naturelle
   `p = P_lf, q = Q_lf`, l'angle étant alors déterminé par l'équilibre.
5. **Effet attendu sur les résultats** : valeurs propres inchangées, un mode
   `λ = 0` qui change simplement d'états porteurs, et les angles affichés qui
   deviennent absolus.

Je ne l'ai pas fait : **aucun interpréteur MATLAB n'est disponible dans cet
environnement**, et cette modification touche `Network.buildGlobalModel` et six
classes d'équipement sans que je puisse exécuter un seul test. Les trois
corrections de bugs, elles, sont locales à quelques lignes et vérifiables à la
lecture ; ce refactor ne l'est pas. Dis-moi si tu veux que je le fasse quand
même — je peux le livrer sur une branche séparée avec, en regard, les valeurs
propres attendues calculées côté Python pour chaque réseau de test, ce qui te
donne le jeu de comparaison pour valider une fois chez toi.
