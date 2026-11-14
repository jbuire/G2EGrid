<!-- Page GENEREE depuis documentation/chapitres/07_exemples.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ch-exemples)=
# Réseaux de test
Quatorze réseaux de test accompagnent l'outil, un notebook chacun dans le dossier `g2egrid/reseaux_test/`. Neuf reprennent les réseaux classiques de V2, trois les réseaux élémentaires de V1, et les deux derniers couplent un réseau de transport à des réseaux de distribution. Le tableau « {ref}`Les réseaux de test <tab-resume-reseaux>` » les résume. Toutes les figures et tous les tableaux de ce chapitre sont produits par le script `generer.py` du dossier `documentation/figures/`, qui rejoue les notebooks : les chiffres sont ceux que l'on obtient en les exécutant.

:::{table} Les réseaux de test. Stable tient compte du mode de jauge (section {ref}`sec-jauge`), qui explique les valeurs propres de l'ordre de $10^{-12}$ des réseaux sans bus infini.
:name: tab-resume-reseaux

| Reseau | Bus | Lignes | Transfos | SM | GFM | GFL | Charges | Etats | $\max \mathrm{Re}\,\lambda$ | Stable |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--:|
| SMIB | 3 | 1 | 1 | 1 | 0 | 0 | 0 | 28 | $+0.0e+00$ | oui |
| SMIB elementaire SM | 2 | 1 | 0 | 1 | 0 | 0 | 1 | 32 | $+0.0e+00$ | oui |
| SMIB elementaire GFM | 2 | 1 | 0 | 0 | 1 | 0 | 1 | 28 | $+0.0e+00$ | oui |
| SMIB elementaire GFL | 2 | 1 | 0 | 0 | 0 | 1 | 1 | 27 | $+0.0e+00$ | oui |
| Trois machines + bus infini | 5 | 2 | 3 | 3 | 0 | 0 | 1 | 72 | $+0.0e+00$ | oui |
| WSCC 9 bus | 9 | 6 | 3 | 3 | 0 | 0 | 3 | 93 | $-4.7e-12$ | oui |
| WSCC 9 bus GFM | 9 | 6 | 3 | 1 | 2 | 0 | 3 | 85 | $-4.1e-12$ | oui |
| Kundur deux zones | 11 | 8 | 4 | 4 | 0 | 0 | 2 | 114 | $+1.6e-12$ | oui |
| IEEE 14 bus | 14 | 20 | 0 | 5 | 0 | 0 | 11 | 181 | $-1.8e-14$ | oui |
| IEEE 33 bus | 33 | 32 | 0 | 0 | 0 | 0 | 32 | 195 | $+0.0e+00$ | oui |
| IEEE 39 bus | 39 | 46 | 0 | 10 | 0 | 0 | 21 | 402 | $-8.7e-13$ | oui |
| CIGRE MT | 15 | 12 | 2 | 0 | 0 | 3 | 13 | 147 | $+0.0e+00$ | oui |
| WSCC 9 bus + CIGRE | 24 | 18 | 6 | 1 | 1 | 4 | 16 | 230 | $+1.2e-11$ | oui |
| WSCC 3 SM + 3 CIGRE | 54 | 42 | 12 | 3 | 0 | 0 | 42 | 429 | $+2.6e-11$ | oui |
:::

Chaque section donne la topologie du réseau, sa carte des valeurs propres et ses modes oscillants de fréquence inférieure à la bande de zoom du notebook, classés par amortissement croissant. Sur les cartes, le panneau de gauche montre tous les modes, l'axe réel en échelle logarithmique symétrique ; celui de droite zoome sur les modes lents et les étiquette par l'élément qui y participe le plus. Sauf mention contraire, les machines gardent leurs réglages par défaut (tableau « {ref}`Paramètres par défaut de la machine synchrone (SM_DEFAULT, identiques à SM.m) <tab-sm-defaut>` »), qui ne sont pas ceux des publications d'origine : les modes obtenus ne se comparent donc pas à ceux de la littérature.

(smib-machine-synchrone-sur-bus-infini)=
## SMIB : machine synchrone sur bus infini
Une machine de 1000 MVA, derrière son transformateur 20/230 kV et une ligne haute tension, alimente un bus infini de 10 GW. C'est le cas d'école du mode électromécanique local. Le transformateur est déclaré dans `T_Transfos` et la table `T_SMs` n'a que trois colonnes. Le mode électromécanique sort à 1,80 Hz et 13 % d'amortissement ; un mode lié au circuit d'excitation apparaît vers 0,92 Hz, amorti à 49 %.

```{figure} ../_static/figures/modes_smib.png
:name: fig-modes-smib
:width: 100%

SMIB : carte des valeurs propres.
```

:::{table} SMIB : modes oscillants lents, par amortissement croissant.
:name: tab-modes-smib

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant       |
|:------------------------|-----------:|--------------:|:--------------------|
| $-1.491 \pm j\,11.279$  |      1.795 |          13.1 | `SM1.Dwr` (35 %)    |
| $-3.272 \pm j\,5.809$   |      0.925 |          49.1 | `SM1.phi_fd` (26 %) |
:::

## Réseaux élémentaires de V1

Les trois réseaux `smib_sm`, `smib_gfm` et `smib_gfl` reprennent les préréglages de `ref_SMIB.m` : un bus infini à 20 kV, une ligne de 1 km, un équipement de production de 0,5 MW et une charge de 1 kW. L'équipement porte son transformateur de raccordement ($R_t = 0{,}0002$, $X_t = 0{,}024$ pu), déclaré par les colonnes `Rt_pu`, `Xt_pu` de sa table. Ces notebooks montrent l'usage courant ; la reproduction exacte de la matrice de V1 demande en plus ses conventions de load flow et d'initialisation (chapitre {ref}`ch-validation`).

**Machine à transformateur intégré.** La `SynchronousMachineTR` alimente aussi une charge de 100 kW à ses bornes ($R_g = 10$ pu), avec un PSS de gain 2 et un amortissement mécanique nul. Son mode électromécanique est à 2,1 Hz et 15 % d'amortissement.

```{figure} ../_static/figures/modes_smib_sm.png
:name: fig-modes-smib-sm
:width: 100%

SMIB élémentaire, machine synchrone : carte des valeurs propres.
```

:::{table} SMIB élémentaire, machine synchrone : modes oscillants lents, par amortissement croissant.
:name: tab-modes-smib-sm

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant       |
|:------------------------|-----------:|--------------:|:--------------------|
| $-2.077 \pm j\,13.224$  |      2.105 |          15.5 | `SM2.Dwr` (46 %)    |
| $-2.085 \pm j\,2.354$   |      0.375 |          66.3 | `SM2.phi_fd` (32 %) |
| $-0.563 \pm j\,0.035$   |      0.006 |          99.8 | `SM2.v1` (40 %)     |
:::

**Grid-forming.** Un `GFM3` à statisme du premier ordre et bus continu dynamique remplace la machine. Son seul mode oscillant lent, dominé par le filtre de puissance active, est à 4,5 Hz et amorti à 50 %.

```{figure} ../_static/figures/modes_smib_gfm.png
:name: fig-modes-smib-gfm
:width: 100%

SMIB élémentaire, grid-forming : carte des valeurs propres.
```

:::{table} SMIB élémentaire, grid-forming : modes oscillants lents, par amortissement croissant.
:name: tab-modes-smib-gfm

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant    |
|:------------------------|-----------:|--------------:|:-----------------|
| $-16.036 \pm j\,28.147$ |      4.480 |          49.5 | `GFM2.pm` (50 %) |
:::

**Grid-following.** Un `GFL2` à PLL et régulation de tension continue remplace la machine. Le mode lent est porté par l'intégrale de la boucle de tension continue, à 4,7 Hz et 72 % d'amortissement.

```{figure} ../_static/figures/modes_smib_gfl.png
:name: fig-modes-smib-gfl
:width: 100%

SMIB élémentaire, grid-following : carte des valeurs propres.
```

:::{table} SMIB élémentaire, grid-following : modes oscillants lents, par amortissement croissant.
:name: tab-modes-smib-gfl

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant     |
|:------------------------|-----------:|--------------:|:------------------|
| $-30.466 \pm j\,29.824$ |      4.747 |          71.5 | `GFL2.M_d` (49 %) |
:::

## Trois machines sur bus infini

Trois machines de 150 MVA, chacune derrière son transformateur 20/230 kV, alimentent un bus de charge de 300 MW relié par deux lignes parallèles à un bus infini. Les modes électromécaniques se regroupent vers 2,2 Hz avec 9 % d'amortissement : ce sont les oscillations des machines entre elles et de l'ensemble contre le réseau. Les modes d'excitation, vers 0,6--0,7 Hz, sont bien amortis.

```{figure} ../_static/figures/modes_three_machine_ib.png
:name: fig-modes-three-machine-ib
:width: 100%

Trois machines sur bus infini : carte des valeurs propres.
```

:::{table} Trois machines sur bus infini : modes oscillants lents, par amortissement croissant.
:name: tab-modes-three-machine-ib

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant       |
|:------------------------|-----------:|--------------:|:--------------------|
| $-1.315 \pm j\,14.181$  |      2.257 |           9.2 | `SM3.Dwr` (22 %)    |
| $-1.279 \pm j\,13.613$  |      2.167 |           9.4 | `SM3.Dwr` (12 %)    |
| $-1.352 \pm j\,14.194$  |      2.259 |           9.5 | `SM1.Dwr` (23 %)    |
| $-3.243 \pm j\,4.640$   |      0.738 |          57.3 | `SM1.phi_fd` (9 %)  |
| $-2.988 \pm j\,3.786$   |      0.603 |          62.0 | `SM1.phi_fd` (13 %) |
| $-2.987 \pm j\,3.718$   |      0.592 |          62.6 | `SM3.phi_fd` (15 %) |
:::

## WSCC 9 bus

Le réseau WSCC d'Anderson et Fouad compte trois machines de 247,5, 192 et 128 MVA, chacune derrière son transformateur élévateur, un anneau de six lignes à 230 kV et trois charges. Ses deux modes électromécaniques inter-machines sortent à 2,20 et 2,01 Hz, amortis à 9 et 11 %. Le mode quasi apériodique à $-0{,}601$ s$^{-1}$ appartient au stabilisateur de puissance.

```{figure} ../_static/figures/modes_wscc9.png
:name: fig-modes-wscc9
:width: 100%

WSCC 9 bus : carte des valeurs propres.
```

:::{table} WSCC 9 bus : modes oscillants lents, par amortissement croissant.
:name: tab-modes-wscc9

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant       |
|:------------------------|-----------:|--------------:|:--------------------|
| $-1.207 \pm j\,13.823$  |      2.200 |           8.7 | `SM3.Dwr` (25 %)    |
| $-1.431 \pm j\,12.638$  |      2.011 |          11.3 | `SM1.Dwr` (17 %)    |
| $-2.348 \pm j\,11.567$  |      1.841 |          19.9 | `SM1.Dwr` (9 %)     |
| $-4.408 \pm j\,8.123$   |      1.293 |          47.7 | `SM1.phi_fd` (12 %) |
| $-3.464 \pm j\,5.405$   |      0.860 |          54.0 | `SM2.phi_fd` (12 %) |
| $-3.164 \pm j\,4.120$   |      0.656 |          60.9 | `SM3.phi_fd` (15 %) |
| $-0.601 \pm j\,0.031$   |      0.005 |          99.9 | `SM3.v1` (28 %)     |
:::

## WSCC 9 bus à onduleurs grid-forming

Les machines G2 et G3 sont remplacées par des `GFM1`, à statisme du premier ordre et bus continu idéal, qui gardent leur filtre par défaut ; leurs transformateurs élévateurs restent dans `T_Transfos`. Il ne reste qu'un mode électromécanique, celui de G1 contre le reste du réseau, à 1,98 Hz et 15 % d'amortissement. La machine de référence étant une machine synchrone, la plus grande valeur propre, à $-4\cdot 10^{-12}$, est le mode de jauge.

```{figure} ../_static/figures/modes_wscc9_gfm.png
:name: fig-modes-wscc9-gfm
:width: 100%

WSCC 9 bus à onduleurs grid-forming : carte des valeurs propres.
```

:::{table} WSCC 9 bus à onduleurs grid-forming : modes oscillants lents, par amortissement croissant.
:name: tab-modes-wscc9-gfm

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant       |
|:------------------------|-----------:|--------------:|:--------------------|
| $-1.860 \pm j\,12.459$  |      1.983 |          14.8 | `SM1.Dwr` (34 %)    |
| $-3.989 \pm j\,6.676$   |      1.063 |          51.3 | `SM1.phi_fd` (27 %) |
| $-29.329 \pm j\,0.099$  |      0.016 |         100.0 | `GFM3.pm` (25 %)    |
:::

## Kundur deux zones

Deux zones de deux machines de 900 MVA chacune sont reliées par deux lignes parallèles de 220 km. Les bus 7 et 9 portent chacun une charge et un banc de condensateurs de 200 et 350 Mvar, déclaré dans `T_Shunts`. Avec les réglages par défaut des machines (PSS de gain 20, régulateur de tension rapide), les modes électromécaniques sortent entre 1,8 et 2,1 Hz, amortis de 10 à 17 % : le mode inter-zones faiblement amorti du livre de Kundur n'apparaît pas tel quel.

```{figure} ../_static/figures/modes_kundur2area.png
:name: fig-modes-kundur2area
:width: 100%

Kundur deux zones : carte des valeurs propres.
```

:::{table} Kundur deux zones : modes oscillants lents, par amortissement croissant.
:name: tab-modes-kundur2area

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant       |
|:------------------------|-----------:|--------------:|:--------------------|
| $-1.269 \pm j\,13.303$  |      2.117 |           9.5 | `SM4.Dwr` (19 %)    |
| $-1.288 \pm j\,13.259$  |      2.110 |           9.7 | `SM2.Dwr` (18 %)    |
| $-1.667 \pm j\,11.313$  |      1.800 |          14.6 | `SM1.Dwr` (9 %)     |
| $-1.944 \pm j\,11.516$  |      1.833 |          16.6 | `SM3.Dwr` (9 %)     |
| $-4.195 \pm j\,7.213$   |      1.148 |          50.3 | `SM3.phi_fd` (8 %)  |
| $-4.482 \pm j\,7.601$   |      1.210 |          50.8 | `SM1.phi_fd` (7 %)  |
| $-3.359 \pm j\,5.016$   |      0.798 |          55.6 | `SM2.phi_fd` (11 %) |
| $-3.338 \pm j\,4.941$   |      0.786 |          56.0 | `SM4.phi_fd` (12 %) |
:::

## IEEE 14 bus

Deux générateurs et trois compensateurs synchrones alimentent onze charges. Les données sont sur une base de tension unique de 100 kV, qui conserve les valeurs en pu de MATPOWER ; les transformateurs de la source sont des lignes en pu sans leur prise. Le condensateur de 19 Mvar du bus 9 est déclaré dans `T_Shunts`, et la demande réactive négative du bus 4 est représentée par une charge de 1 % de $P$ et un shunt. Les cinq modes électromécaniques s'étalent de 1,8 à 2,2 Hz, avec 14 à 24 % d'amortissement.

```{figure} ../_static/figures/modes_ieee14.png
:name: fig-modes-ieee14
:width: 100%

IEEE 14 bus : carte des valeurs propres.
```

:::{table} IEEE 14 bus : modes oscillants lents, par amortissement croissant.
:name: tab-modes-ieee14

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant       |
|:------------------------|-----------:|--------------:|:--------------------|
| $-1.786 \pm j\,12.765$  |      2.032 |          13.9 | `SM8.Dwr` (9 %)     |
| $-2.075 \pm j\,14.120$  |      2.247 |          14.5 | `SM2.Dwr` (26 %)    |
| $-2.198 \pm j\,12.844$  |      2.044 |          16.9 | `SM8.Dwr` (14 %)    |
| $-2.472 \pm j\,13.035$  |      2.075 |          18.6 | `SM3.Dwr` (25 %)    |
| $-2.726 \pm j\,11.104$  |      1.767 |          23.8 | `SM1.Dwr` (6 %)     |
| $-4.409 \pm j\,8.533$   |      1.358 |          45.9 | `SM1.Pm` (5 %)      |
| $-3.404 \pm j\,6.253$   |      0.995 |          47.8 | `SM8.phi_fd` (24 %) |
| $-2.918 \pm j\,5.188$   |      0.826 |          49.0 | `SM6.phi_fd` (27 %) |
:::

## IEEE 33 bus

Le réseau de distribution radial de Baran et Wu, à 12,66 kV, porte 3,7 MW de charge sur 32 bus, alimentés par un bus infini au poste source. Ses lignes n'ont pas de capacité : tous les nœuds sont algébriques et la dynamique ne vient que des branches et des charges. Il n'a donc aucun mode oscillant lent, et toutes ses valeurs propres sont réelles ou très rapides.

```{figure} ../_static/figures/modes_ieee33.png
:name: fig-modes-ieee33
:width: 100%

IEEE 33 bus : carte des valeurs propres.
```

(tab-modes-ieee33)=

**

<table>
<thead>
<tr>
<th style="text-align: left;">Valeur propre <span class="math inline"><em>λ</em></span></th>
<th style="text-align: right;"><span class="math inline"><em>f</em></span> [Hz]</th>
<th style="text-align: right;"><span class="math inline"><em>ζ</em></span> [%]</th>
<th style="text-align: left;">Etat dominant</th>
</tr>
</thead>
<tbody>
<tr>
<td colspan="4" style="text-align: center;">aucun mode oscillant sous 5 Hz</td>
</tr>
</tbody>
</table>

## IEEE 39 bus

Le réseau New England compte dix machines, dont une grosse machine au bus 39 qui représente le reste du réseau interconnecté, 46 branches et 21 charges. Avec 402 états, c'est le plus gros des réseaux de test. Ses modes électromécaniques forment un amas entre 1,9 et 2,3 Hz, amortis à 8--10 %.

```{figure} ../_static/figures/modes_ieee39.png
:name: fig-modes-ieee39
:width: 100%

IEEE 39 bus : carte des valeurs propres.
```

:::{table} IEEE 39 bus : modes oscillants lents, par amortissement croissant.
:name: tab-modes-ieee39

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant     |
|:------------------------|-----------:|--------------:|:------------------|
| $-1.150 \pm j\,13.568$  |      2.159 |           8.4 | `SM37.Dwr` (20 %) |
| $-1.188 \pm j\,13.985$  |      2.226 |           8.5 | `SM33.Dwr` (19 %) |
| $-1.165 \pm j\,13.406$  |      2.134 |           8.7 | `SM32.Dwr` (22 %) |
| $-1.263 \pm j\,14.126$  |      2.248 |           8.9 | `SM35.Dwr` (20 %) |
| $-1.203 \pm j\,13.061$  |      2.079 |           9.2 | `SM39.Dwr` (16 %) |
| $-1.215 \pm j\,12.874$  |      2.049 |           9.4 | `SM36.Dwr` (10 %) |
| $-1.222 \pm j\,12.176$  |      1.938 |          10.0 | `SM38.Dwr` (24 %) |
| $-1.259 \pm j\,12.342$  |      1.964 |          10.1 | `SM31.Dwr` (12 %) |
:::

## CIGRE moyenne tension

Le réseau européen de distribution CIGRE TF C6.04.02, en exploitation radiale, est alimenté par un réseau amont à 110 kV à travers deux transformateurs de 25 MVA. Trois GFL y sont raccordés : un aérogénérateur de 1,5 MVA et deux installations photovoltaïques de 1 MVA. Ses données suivent pandapower. Ses treize charges, déclarées dans `Y_DER`, deviennent des charges R-L et portent le modèle à 147 états. Le réseau n'a aucun mode oscillant sous 5 Hz : les modes des GFL, de leurs PLL et de leurs boucles de courant sont plus rapides.

```{figure} ../_static/figures/modes_cigre_mv.png
:name: fig-modes-cigre-mv
:width: 100%

CIGRE moyenne tension : carte des valeurs propres.
```

(tab-modes-cigre-mv)=

**

<table>
<thead>
<tr>
<th style="text-align: left;">Valeur propre <span class="math inline"><em>λ</em></span></th>
<th style="text-align: right;"><span class="math inline"><em>f</em></span> [Hz]</th>
<th style="text-align: right;"><span class="math inline"><em>ζ</em></span> [%]</th>
<th style="text-align: left;">Etat dominant</th>
</tr>
</thead>
<tbody>
<tr>
<td colspan="4" style="text-align: center;">aucun mode oscillant sous 5 Hz</td>
</tr>
</tbody>
</table>

(sec-wscc9-cigre)=
## WSCC 9 bus et distribution CIGRE
Réseau de transport WSCC dont les trois groupes sont de trois technologies : G1 reste une machine synchrone et sert de référence d'angle, G2 devient un `GFM1` sur un bus PV à 1,00 pu et G3 un `GFL1` sur un bus PQ. Au bus de charge 5, un poste 230/110 kV alimente un réseau de distribution CIGRE MT complet, avec ses deux départs, ses treize charges et ses trois GFL.

**Échelle de puissance.** Un réseau CIGRE porte 25,8 MW de charge quand le bus 5 en consomme 125. Le poste alimente donc $K = 4$ réseaux CIGRE identiques, agrégés en un seul : toutes les puissances sont multipliées par $K$ et toutes les impédances divisées par $K$, si bien que les tensions, les pertes relatives et les dynamiques en pu restent celles d'un seul réseau CIGRE. La distribution porte alors 103 MW de charge et 11 MW de production ; les 21,6 MW restants demeurent une charge agrégée au bus 5, au facteur de puissance d'origine.

**Plan de tension.** La prise du poste 230/110 kV ($U_2 = 1{,}025 \times
110$ kV) remet le 110 kV à 1,02 pu, la tension de la source amont du CIGRE d'origine ; les transformateurs 110/20 kV restent sans prise. Le transport reste entre 1,01 et 1,04 pu, la moyenne tension entre 0,98 et 1,02 pu.

**Numérotation.** Tout ce qui appartient à la distribution porte un numéro supérieur ou égal à 100, poste 230/110 kV compris ; le bus 5 est le point de raccordement. C'est la convention qu'utilise le chapitre {ref}`ch-td` pour séparer les deux réseaux. Avec 230 états, ce réseau est le deuxième plus gros des réseaux de test. Son mode électromécanique, celui de G1 contre le reste, sort à 1,95 Hz et 15 % d'amortissement.

```{figure} ../_static/figures/modes_wscc9_cigre.png
:name: fig-modes-wscc9-cigre
:width: 100%

WSCC 9 bus et distribution CIGRE : carte des valeurs propres.
```

:::{table} WSCC 9 bus et distribution CIGRE : modes oscillants lents, par amortissement croissant.
:name: tab-modes-wscc9-cigre

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant       |
|:------------------------|-----------:|--------------:|:--------------------|
| $-1.912 \pm j\,12.252$  |      1.950 |          15.4 | `SM1.Dwr` (33 %)    |
| $-4.162 \pm j\,7.074$   |      1.126 |          50.7 | `SM1.phi_fd` (27 %) |
| $-34.462 \pm j\,1.561$  |      0.248 |          99.9 | `SM1.Dw1` (26 %)    |
| $-30.315 \pm j\,0.051$  |      0.008 |         100.0 | `GFM2.pm` (37 %)    |
:::

(sec-wscc3sm-cigre3)=
## WSCC 9 bus et trois distributions CIGRE
Le réseau WSCC d'Anderson et Fouad tel quel -- trois machines synchrones, leur dispatch d'origine, aucun convertisseur --, dont *chacun* des trois bus de charge alimente un réseau CIGRE MT complet par un poste 230/110 kV. C'est le point de départ de la transition du chapitre {ref}`ch-transition` : un système entièrement synchrone, auquel on ajoute ensuite des convertisseurs.

**Échelle de puissance.** Un réseau CIGRE porte 25,8 MW de charge ; les bus 5, 6 et 8 en consomment 125, 90 et 100. Les postes alimentent donc $K = 4$, $3$ et $3$ réseaux CIGRE identiques, agrégés en un seul chacun. La charge résiduelle du bus WSCC est ce que la distribution ne porte pas, **en $P$ comme en $Q$** : chaque bus de charge garde exactement sa consommation WSCC, le point de fonctionnement du transport est celui du réseau d'origine, et les modes électromécaniques sortent bien à 2,21 et 2,01 Hz comme dans le tableau « {ref}`WSCC 9 bus <tab-modes-wscc9>` ».

**Plan de tension.** Chaque poste source porte un gradin de condensateurs (35, 25 et 25 Mvar, déclarés dans `T_Shunts`) qui compense la consommation réactive de ses transformateurs ; sans lui, la moyenne tension descendrait à 0,90 pu. Aucune prise hors nominale n'est nécessaire. Le transport reste entre 1,018 et 1,040 pu, la moyenne tension entre 0,954 et 1,038 pu.

**Numérotation.** La distribution du bus 5 porte les numéros 100 à 114, celle du bus 6 les numéros 200 à 214, celle du bus 8 les numéros 300 à 314 ; le poste 230/110 kV appartient à la distribution. Un numéro supérieur ou égal à 100 signale donc la distribution, comme dans `wscc9_cigre`.

```{figure} ../_static/figures/modes_wscc3sm_cigre3.png
:name: fig-modes-wscc3sm-cigre3
:width: 100%

WSCC 9 bus et trois distributions CIGRE : carte des valeurs propres.
```

:::{table} WSCC 9 bus et trois distributions CIGRE : modes oscillants lents, par amortissement croissant.
:name: tab-modes-wscc3sm-cigre3

| Valeur propre $\lambda$ | $f$ \[Hz\] | $\zeta$ \[%\] | Etat dominant       |
|:------------------------|-----------:|--------------:|:--------------------|
| $-1.125 \pm j\,13.870$  |      2.208 |           8.1 | `SM3.Dwr` (24 %)    |
| $-1.428 \pm j\,12.607$  |      2.006 |          11.3 | `SM1.Dwr` (18 %)    |
| $-2.294 \pm j\,11.645$  |      1.853 |          19.3 | `SM1.phi_fd` (8 %)  |
| $-4.323 \pm j\,8.135$   |      1.295 |          46.9 | `SM1.phi_fd` (12 %) |
| $-3.409 \pm j\,5.255$   |      0.836 |          54.4 | `SM2.phi_fd` (12 %) |
| $-3.103 \pm j\,3.881$   |      0.618 |          62.4 | `SM3.phi_fd` (15 %) |
| $-0.585 \pm j\,0.059$   |      0.009 |          99.5 | `SM3.v1` (23 %)     |
:::
