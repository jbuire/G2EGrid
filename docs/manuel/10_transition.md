<!-- Page GENEREE depuis documentation/chapitres/10_transition.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ch-transition)=
# La transition, mode par mode
Un système électrique ne passe pas des machines synchrones aux onduleurs d'un coup : les convertisseurs arrivent d'abord par la distribution, en petites unités, puis gagnent le transport. Ce chapitre suit cette progression sur un même réseau et regarde ce qu'elle fait au spectre : quels modes disparaissent, lesquels naissent, et à quelle condition l'un d'eux devient dangereux.

## Le réseau et l'échelle des scénarios

Le point de départ est le réseau `wscc3sm_cigre3` (section {ref}`sec-wscc3sm-cigre3`) : la WSCC 9 bus à *trois machines synchrones*, dont chacun des trois bus de charge alimente un réseau de distribution CIGRE MT complet. 429 états, 54 bus, aucun convertisseur. On lui fait subir sept étages, en ajoutant d'abord des grid-following dans un départ, puis dans les trois, puis en passant les plus gros producteurs en grid-forming, et enfin en remplaçant les machines du transport.

À chaque étage, la production distribuée **efface** la production synchrone : G2 et G3 réduisent leur consigne au prorata, de sorte que la charge et le transit du slack ne bougent pas. La transition se lit donc sur le taux d'onduleurs, et non sur le point de fonctionnement.

:::{table} Les sept scénarios. Modes compte les modes oscillants sous 40 Hz, Électrom. ceux dont l'état dominant est une vitesse de rotor.
:name: tab-transition-resume

|  | Scenario | Etats | SM | GFM | GFL | IBR \[%\] | Modes | Electrom. | $\zeta_{\min}$ \[%\] |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| S0 | 3 SM, aucun convertisseur | 429 | 3 | 0 | 0 | 0.0 | 6 | 2 | 8.1 |
| S1 | 3 GFL dans D1 | 465 | 3 | 0 | 3 | 2.4 | 9 | 2 | 8.2 |
| S2 | 11 GFL dans D1 | 561 | 3 | 0 | 11 | 5.6 | 17 | 2 | 8.4 |
| S3 | 11 GFL par depart | 825 | 3 | 0 | 33 | 12.9 | 39 | 2 | 8.9 |
| S4 | 3 GFM par depart | 834 | 3 | 9 | 24 | 12.9 | 54 | 3 | 9.0 |
| S5 | G2 en GFM | 830 | 2 | 10 | 24 | 43.1 | 54 | 2 | 9.7 |
| S6 | G2 et G3 en GFM | 826 | 1 | 11 | 24 | 62.6 | 54 | 1 | 15.1 |
:::

## Où vont les modes

```{figure} ../_static/figures/transition_modes.png
:name: fig-transition-modes
:width: 100%

Migration des modes le long de la transition, de S0 (violet) à S6 (jaune).
```

**Les modes électromécaniques se raréfient, et ceux qui restent s'amortissent.** Trois machines portent deux oscillations inter-machines, deux machines une seule, une machine plus qu'une : celle où elle oscille contre le reste du système. Et leur amortissement *augmente* avec la pénétration : 8,1 à 8,9 % pour le mode de G3, 11,3 à 11,9 % pour celui de G1 contre G2, entre S0 et S3. Ce n'est pas un effet de la commande des onduleurs, c'est un effet de dispatch : la production synchrone qui s'efface rapproche les angles internes du réseau, et l'oscillation entre rotors porte moins d'énergie.

:::{admonition} Les GFM deviennent l'autre extrémité de l'oscillation
:class: note

À S4, un *troisième* mode électromécanique apparaît, à 1,94 Hz et 20,7 % d'amortissement, alors que le transport porte toujours ses trois machines. Les GFM de la distribution viennent de prendre un comportement de source de tension : la machine G1 oscille désormais contre eux. Un grid-forming ne supprime pas les oscillations de rotor, il s'y substitue comme partenaire.
:::

**De nouvelles familles naissent.** Le nombre de modes oscillants sous 40 Hz passe de 6 à 54 : boucles de statisme et filtres de puissance des GFM, PLL et boucles de courant des GFL. Le spectre se densifie entre 5 et 40 Hz, une bande où un système purement synchrone n'avait presque rien.

```{figure} ../_static/figures/transition_familles.png
:name: fig-transition-familles
:width: 100%

À gauche, les modes par famille d'équipement (celle de l'état qui participe le plus). À droite, l'amortissement de tous les modes de la bande.
```

## Le mode dangereux

Avec les réglages par défaut de la bibliothèque, *aucun* des sept scénarios n'est instable. Mais ces réglages ne sont pas anodins : `GFM_DEFAULT` donne aux GFM un statisme $m_p = 0{,}001$, soit **0,1 %**, qui en fait une source de tension quasi rigide. Un statisme de groupe réel vaut 2 à 5 %.

:::{table} Amortissement du mode le moins amorti (sous 40 Hz), selon le statisme des GFM de *distribution*. Les trois colonnes sont identiques.
:name: tab-transition-statisme

| Statisme \[%\] | $f$ \[Hz\] |     S4 |     S5 |     S6 |
|---------------:|-----------:|-------:|-------:|-------:|
|           0.10 |        2.0 |  +9.01 |  +9.74 | +15.10 |
|           0.20 |        2.0 |  +9.02 |  +9.87 | +15.49 |
|           0.50 |       11.7 |  +8.99 |  +9.94 | +14.47 |
|           0.75 |       14.7 |  +7.09 |  +7.09 |  +7.09 |
|           1.00 |       17.2 |  +1.19 |  +1.19 |  +1.19 |
|           1.25 |       19.3 |  -3.89 |  -3.89 |  -3.89 |
|           1.50 |       21.0 |  -8.33 |  -8.33 |  -8.33 |
|           2.00 |       23.8 | -15.58 | -15.58 | -15.58 |
|           3.00 |       27.5 | -25.49 | -25.49 | -25.49 |
:::

Trois colonnes identiques à la deuxième décimale : que le transport porte trois machines synchrones (S4), deux (S5) ou une seule (S6), le mode le moins amorti du système est le même, et il change de signe au même endroit, à **1,06 %** de statisme. C'est la signature d'un mode qui ne voit pas le transport.

```{figure} ../_static/figures/transition_statisme.png
:name: fig-transition-statisme
:width: 100%

À gauche, le lieu du mode le moins amorti quand le statisme des GFM de distribution croît de 0,1 à 3 % ; en dessous de 0,75 %, le mode le moins amorti est encore l'électromécanique. À droite, son amortissement : le seuil ne dépend pas du transport.
```

(distribution-ou-transport)=
### Distribution ou transport ?
Le même réglage, 2 % de statisme, appliqué aux seuls GFM de la distribution puis aux seuls GFM du transport, tranche la question.

:::{table} Le même statisme de 2 %, selon où on l'applique, sur le scénario S6.
:name: tab-transition-ou

| Reglage (S6) | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant | $\max \mathrm{Re}\,\lambda$ |
|:---|---:|---:|:---|:---|
| defaut (0,1 % partout) | 2.02 | +15.10 | SM1.Dwr | $-9{,}8\cdot 10^{-11}$ |
| 2 % en distribution | 23.81 | -15.58 | GFM209.pm | $2{,}4\cdot 10^{1}$ |
| 2 % en transport | 2.02 | +18.06 | SM1.Dwr | $-9{,}6\cdot 10^{-11}$ |
| 2 % partout | 23.81 | -15.58 | GFM209.pm | $2{,}4\cdot 10^{1}$ |
:::

Sur les GFM du **transport**, 2 % de statisme *améliore* l'amortissement, de 15,1 à 18,1 %. Sur ceux de la **distribution**, le système devient instable : trois modes à 23,8 Hz, un par départ, à $-15{,}6\,\%$ d'amortissement. Les états qui y participent le plus le disent : `GFM209.pm` et `GFM209.theta`, avec `GFM207.pm` --- deux GFM du même départ, qui oscillent l'un contre l'autre par leur boucle de statisme.

**La cause est le rapport $R/X$.** Une ligne du CIGRE MT a $R/X = 0{,}82$ en moyenne ; une ligne de la WSCC, $0{,}16$. Le statisme $P$--$f$ suppose un réseau inductif, où la puissance active répond à l'angle. Dans un réseau résistif, l'angle commande surtout le réactif, et la boucle $P$--$f$ se referme à l'envers : plus elle est souple, plus elle amplifie. C'est un résultat connu des réseaux de distribution, et le modèle le retrouve sans qu'on le lui ait dit.

## Ce qu'un modèle du transport seul ne voit pas

Le chapitre {ref}`ch-td` a montré qu'un modèle réduit au transport, $\bm{A}_{TT}$, ignore la dynamique de la distribution. Le mode instable en est l'illustration la plus brutale, puisqu'il n'est fait que d'états de distribution :

- **$\bm{A}_{TT}$ seul est stable** alors que le système complet ne l'est pas. Il ne contient aucun état de GFM de distribution : il ne peut pas voir un mode qui n'est fait que de ceux-là.

- **Le quasi statique est instable, mais pour de mauvaises raisons** : il annonce un mode à $+410\,\mathrm{s}^{-1}$, sans rapport avec les $+23{,}6\,\mathrm{s}^{-1}$ du vrai. La substitution $\bm{A}_{qs} = \bm{A}_{TT} - \bm{A}_{TD}\bm{A}_{DD}^{-1}\bm{A}_{DT}$ remplace les états de la distribution par leur équilibre : quand $\bm{A}_{DD}$ porte elle-même des modes instables, cet équilibre n'existe pas et le résultat n'a de sens ni comme prédiction, ni comme alerte.

C'est la justification la plus directe du chapitre précédent : le quasi statique suffit tant que la distribution est rapide *et amortie*.

## Ce qu'il faut retenir

1.  **La transition, vue des modes électromécaniques, est une bonne nouvelle** : ils se raréfient et s'amortissent, tant que les onduleurs sont correctement réglés.

2.  **Le spectre se densifie** entre 5 et 40 Hz, une bande que la surveillance d'un système synchrone ignorait.

3.  **Le danger ne vient pas du taux d'onduleurs, mais du couple emplacement / réglage** : 2 % de statisme est sain en transport et fatal en distribution, et le seuil ne dépend pas de ce que fait le transport.

4.  **Un modèle du transport seul ne voit rien venir.**

## Outils

Côté Python, tout est dans le notebook `transition_ibr.ipynb` du dossier `python/etudes/`, qui construit les scénarios avec `g2egrid.reseau.tables` : `Tables` porte les tables du réseau, `Tables.copie` en tire une variante, `assembler` construit le `Network`, et `Tables.PARAMS_PAR_BUS` règle autrement les convertisseurs de la distribution et ceux du transport. Côté MATLAB, `etudes/transition_ibr.m` refait la même étude à partir de `WSCC3SM_CIGRE3.m`, dont l'argument optionnel décrit le scénario ; les deux outils donnent les mêmes chiffres à $5\cdot 10^{-4}$ près sur les modes lents.

`etudes/transition_ibr_emt.m` construit en outre le modèle Simulink (Simscape Electrical) des cas extrêmes par `SimulinkGenerator`, avec son excitation et sa sonde déjà branchées. Il l'enregistre *sans* le simuler : sur ce réseau de 54 bus et près de deux cents blocs, la seule génération demande une dizaine de minutes et la simulation bien davantage.

:::{admonition} Un défaut du générateur Simulink, révélé par la taille du réseau
:class: note

`calculatePositions` indexait la position des blocs par le *rang* du nœud dans la table, et non par son identifiant. Les deux coïncident tant que les bus sont numérotés de 1 à $N$ : c'est le cas de tous les réseaux sur lesquels le générateur avait tourné jusqu'ici. Dès qu'un réseau de distribution porte des numéros supérieurs à 100, la position d'une ligne devient introuvable et la génération s'arrête. C'est la même confusion rang / identifiant que celle corrigée dans `buildFromMatrices` (section {ref}`sec-v2-multimachine`).
:::

```python
from g2egrid.reseau.tables import Tables, assembler

t = T0.copie(nom="S6")                       # les tables du reseau de depart
for prefixe, _amont, K in DEPARTS:
    raccorder(t, prefixe, K, 11, part_gfm=3) # 11 producteurs, dont 3 GFM
effacer(t)                                   # G2 et G3 s'effacent
en_gfm(t, 2, 200e6); en_gfm(t, 3, 130e6)     # le transport bascule
t.PARAMS_PAR_BUS = {b: dict(mp=0.02) for b in bus_dist}   # statisme realiste

net = assembler(t)
net.run_loadflow(); net.solve_equilibrium()
ss = net.build_state_space()
```
