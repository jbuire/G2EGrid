<!-- Page GENEREE depuis documentation/chapitres/08_sources.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ch-sources)=
# Source de tension, source de courant
Une machine synchrone et un onduleur *grid-forming* tiennent la tension du réseau ; un onduleur *grid-following* y injecte un courant. Cette distinction, centrale pour l'exploitation d'un réseau à forte part d'électronique de puissance, n'est pas une convention de vocabulaire : elle se lit directement sur l'admittance des modèles. Ce chapitre la mesure sur les huit classes d'équipement.

## L'admittance d'un équipement

Chaque équipement a deux entrées de port ($v_d$, $v_q$) et deux sorties de port ($i_d$, $i_q$, courant *injecté* au nœud, section {ref}`sec-element`). En figeant ses autres entrées -- ses consignes, $\omega_{ref}$ et $\theta_g$ --, sa linéarisation autour d'un point de fonctionnement donne une matrice $2\times 2$ : 

$$
\begin{bmatrix} \Delta i_d \\ \Delta i_q \end{bmatrix}
 = -\,\bm{Y}(s) \begin{bmatrix} \Delta v_d \\ \Delta v_q \end{bmatrix},
\qquad
\bm{Y}(s) = -\left[\bm{C}(s\bm{I}-\bm{A})^{-1}\bm{B} + \bm{D}\right],
$$

 où $\bm{A}$, $\bm{B}$, $\bm{C}$, $\bm{D}$ sont les matrices de l'équipement seul, après élimination de ses variables algébriques (section {ref}`sec-assemblage`). Le signe fait de $\bm{Y}$ l'admittance *d'entrée*, courant compté vers l'équipement, ce qui rend aux mots leur sens usuel :

|  | $\Vert\bm{Y}\Vert$ | réaction à une variation de tension du réseau |
|:---|:---|:---|
| Source de tension | grand | débite un fort courant pour tenir sa tension |
| Source de courant | petit | garde son courant, laisse la tension varier |
:::

### Deux directions, pas un scalaire

Une matrice dq n'est pas un nombre, et sa décomposition en vecteur complexe ne suffit pas : une charge à puissance constante, par exemple, n'a qu'une composante inverse. On trace donc le gain $\Vert \bm{Y} e\Vert$ dans deux directions du plan dq, prises au point de fonctionnement :

- la direction **module**, portée par la tension de bornes : une perturbation dans cette direction est une variation du module de la tension ;

- la direction **angle**, orthogonale : une variation de l'angle.

Il faut bien prendre la norme $\Vert \bm{Y}e\Vert$ et non la projection $e^{\!\top}\bm{Y}e$ : sous une impédance inductive, une variation de module appelle un courant en quadrature, que la projection annulerait.

:::{admonition} Pourquoi la direction angle est molle pour tout le monde
:class: note

Un équipement regardé *seul* n'a aucune référence d'angle absolu : son angle interne est un intégrateur libre. Sa direction angle est donc molle en basse fréquence, machine synchrone comprise. C'est le réseau, à travers les autres équipements, qui fixe cet angle -- c'est aussi l'origine du mode de jauge (section {ref}`sec-jauge`). La nature de la source se lit dans la direction module.
:::

## Banc d'essai

Les huit équipements sont pris à la même puissance nominale que la base système (1 MVA, 20 kV, 50 Hz), seuls à leur nœud face à un réseau infini, au même point de fonctionnement : $V = 1\angle 0$ pu, $P = 0{,}5$ pu, $Q = 0$. Les réglages sont ceux par défaut des classes (chapitre {ref}`ch-modeles`) ; la machine à transformateur intégré porte le transformateur de `ref_SMIB.m` ($R_t = 0{,}0002$, $X_t = 0{,}024$ pu). Seul le point de fonctionnement compte : la linéarisation d'un équipement ne dépend pas du reste du réseau.

## Résultats

```{figure} ../_static/figures/sources_admittance.png
:name: fig-sources
:width: 100%

Admittance des équipements vue de leurs bornes, dans la direction du module de tension (gauche) et de l'angle (droite).
```

(tab-sources)=

**

<table>
<thead>
<tr>
<th style="text-align: left;"></th>
<th style="text-align: left;"></th>
<th colspan="5" style="text-align: center;"><span class="math inline">‖<strong>Y</strong><em>e</em><sub><em>m</em><em>o</em><em>d</em></sub>‖</span> [pu] a</th>
<th style="text-align: right;"><span class="math inline">1/‖<strong>Y</strong><em>e</em><sub><em>m</em><em>o</em><em>d</em></sub>‖</span></th>
<th style="text-align: left;"></th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><span>3-7</span>(lr)<span>8-8</span> Classe</td>
<td style="text-align: left;">Modele</td>
<td style="text-align: right;">0.01 Hz</td>
<td style="text-align: right;">0.1 Hz</td>
<td style="text-align: right;">1 Hz</td>
<td style="text-align: right;">10 Hz</td>
<td style="text-align: right;">100 Hz</td>
<td style="text-align: right;">a 0,1 Hz [pu]</td>
<td style="text-align: left;">Comportement</td>
</tr>
<tr>
<td style="text-align: left;"><code>SM</code></td>
<td style="text-align: left;">machine synchrone</td>
<td style="text-align: right;">224.20</td>
<td style="text-align: right;">37.74</td>
<td style="text-align: right;">3.90</td>
<td style="text-align: right;">4.00</td>
<td style="text-align: right;">3.10</td>
<td style="text-align: right;">0.0265</td>
<td style="text-align: left;">tension</td>
</tr>
<tr>
<td style="text-align: left;"><code>SMTR</code></td>
<td style="text-align: left;">machine + transfo</td>
<td style="text-align: right;">35.53</td>
<td style="text-align: right;">25.23</td>
<td style="text-align: right;">3.81</td>
<td style="text-align: right;">3.67</td>
<td style="text-align: right;">2.82</td>
<td style="text-align: right;">0.0396</td>
<td style="text-align: left;">tension</td>
</tr>
<tr>
<td style="text-align: left;"><code>GFM1</code></td>
<td style="text-align: left;">GFM statisme 1er ordre, DC ideal</td>
<td style="text-align: right;">90.90</td>
<td style="text-align: right;">90.90</td>
<td style="text-align: right;">91.58</td>
<td style="text-align: right;">102.44</td>
<td style="text-align: right;">67.35</td>
<td style="text-align: right;">0.0110</td>
<td style="text-align: left;">tension</td>
</tr>
<tr>
<td style="text-align: left;"><code>GFM2</code></td>
<td style="text-align: left;">GFM statisme 2e ordre, DC ideal</td>
<td style="text-align: right;">90.90</td>
<td style="text-align: right;">90.90</td>
<td style="text-align: right;">91.92</td>
<td style="text-align: right;">98.60</td>
<td style="text-align: right;">67.38</td>
<td style="text-align: right;">0.0110</td>
<td style="text-align: left;">tension</td>
</tr>
<tr>
<td style="text-align: left;"><code>GFM3</code></td>
<td style="text-align: left;">GFM statisme 1er ordre, DC dynamique</td>
<td style="text-align: right;">90.90</td>
<td style="text-align: right;">90.90</td>
<td style="text-align: right;">91.58</td>
<td style="text-align: right;">102.44</td>
<td style="text-align: right;">67.35</td>
<td style="text-align: right;">0.0110</td>
<td style="text-align: left;">tension</td>
</tr>
<tr>
<td style="text-align: left;"><code>GFM4</code></td>
<td style="text-align: left;">GFM statisme 2e ordre, DC dynamique</td>
<td style="text-align: right;">90.90</td>
<td style="text-align: right;">90.90</td>
<td style="text-align: right;">91.92</td>
<td style="text-align: right;">98.60</td>
<td style="text-align: right;">67.38</td>
<td style="text-align: right;">0.0110</td>
<td style="text-align: left;">tension</td>
</tr>
<tr>
<td style="text-align: left;"><code>GFL1</code></td>
<td style="text-align: left;">GFL boucles P/Q, DC ideal</td>
<td style="text-align: right;">0.50</td>
<td style="text-align: right;">0.50</td>
<td style="text-align: right;">0.49</td>
<td style="text-align: right;">0.23</td>
<td style="text-align: right;">0.25</td>
<td style="text-align: right;">2.0063</td>
<td style="text-align: left;">courant</td>
</tr>
<tr>
<td style="text-align: left;"><code>GFL2</code></td>
<td style="text-align: left;">GFL boucle Vdc/Q, DC dynamique</td>
<td style="text-align: right;">8.91</td>
<td style="text-align: right;">8.91</td>
<td style="text-align: right;">8.82</td>
<td style="text-align: right;">5.13</td>
<td style="text-align: right;">1.41</td>
<td style="text-align: right;">0.1122</td>
<td style="text-align: left;">intermediaire</td>
</tr>
</tbody>
</table>

**Les grid-forming** affichent $\Vert\bm{Y}\Vert \approx 91$ pu jusqu'à 1 Hz, soit une impédance de 0,011 pu. C'est, à 5 % près, l'inverse de l'impédance de leur filtre de raccordement, $|R_2 + jX_2| = 0{,}0104$ pu : la boucle de tension tient le condensateur, et le réseau ne voit plus que cette branche. Les quatre variantes GFM1 à GFM4 se superposent : ni l'ordre du statisme ni la dynamique du bus continu ne changent la raideur.

**La machine synchrone** atteint 224 pu à 0,01 Hz -- c'est le régulateur de tension qui tient la tension de bornes --, puis s'effondre dès 0,5 Hz, quand l'AVR ne suit plus, jusqu'aux réactances transitoires vers 2,5 pu. La `SynchronousMachineTR` plafonne bien plus bas, à 36 pu : son transformateur de raccordement est en série et limite à lui seul l'admittance à $1/X_t \approx 42$ pu.

**Le grid-following à boucles P/Q** reste à 0,5 pu, soit une impédance de 2 pu. Ce n'est pas zéro : ses boucles tiennent $P$ et $Q$, donc son courant suit $\overline{S/V}$ et $\Vert\bm{Y}\Vert = S/V^2 = 0{,}5$. Son minimum, 0,11 pu vers 34 Hz, est la source de courant pure : les boucles externes ne suivent plus, la boucle de courant tient encore.

**Le grid-following à régulation de tension continue** est nettement plus raide, 8,9 pu : en tenant $v_{dc}$, il ajuste sa puissance à la tension alternative. Il reste dix fois moins raide qu'un grid-forming.

:::{admonition} Les deux résonances
:class: note

Le pic vers 50 Hz des machines et des GFM n'est pas une résonance physique : à 50 Hz dans le repère dq, l'une des deux composantes de la perturbation tombe à 0 Hz en grandeurs de phase, où l'impédance d'une inductance s'annule. Le pic à 1,5 kHz des GFL est, lui, la résonance du filtre LCL, que la boucle de tension d'un grid-forming amortit. Au-delà de 2 kHz, tous les équipements convergent vers l'admittance de leur inductance de raccordement : aucune commande ne suit plus, il ne reste que la physique du filtre.
:::

(sec-sources-transport)=
## À réglages comparables
Les valeurs par défaut des classes ne sont *pas* comparables d'une famille à l'autre, et une partie des écarts ci-dessus en vient :

- la boucle de courant est réglée à 0,1 ms pour un GFM contre 10 ms pour un GFL, soit un facteur 100 sur la bande passante ;

- le statisme d'un GFM vaut 0,1 % contre 0,5 % pour le gouverneur de la machine, et son statisme de tension 0,1 % contre 4 % dans les codes de réseau ;

- surtout, un GFM par défaut n'a pour tout raccordement que son filtre de 0,01 pu, alors qu'une machine de transport est derrière un transformateur élévateur de 0,15 pu.

On reprend donc le calcul avec un **profil commun de gros groupe de transport**, appliqué à la machine comme aux convertisseurs :

| Grandeur | Valeur | Appliqué à |
|:---|:---|:---|
| Transformateur élévateur | $R = 0{,}005$, $X = 0{,}15$ pu | machine (`Rt`, `Lt`), convertisseurs (`R2_pu`, `L2_pu`) |
| Statisme de fréquence | 4 % | machine (gouverneur) et GFM (`mp`) |
| Inertie | $H = 5$ s | machine ; GFM par $\omega_f = 1/(2 m_p H)$ |
| Régulation de tension | AVR $K_a = 300$ | machine |
| Statisme de tension | 4 % ou $1/K_a = 0{,}33$ % | GFM (`nq`) |
| Boucle de courant | $\tau = 1$ ms | tous les convertisseurs |
:::

L'équivalence $H = 1/(2 m_p \omega_f)$ est celle d'un statisme du premier ordre filtré, lu comme une équation du mouvement : c'est la façon usuelle de donner à un GFM l'inertie d'une machine. Les boucles externes du GFL (P/Q à 100 ms, PLL à 50 ms) restent à leurs valeurs par défaut, faute d'équivalent côté machine.

```{figure} ../_static/figures/sources_transport.png
:name: fig-sources-transport
:width: 100%

Même profil de gros groupe de transport. Le repère pointillé est <span class="math inline">1/<em>X</em><sub><em>T</em></sub> = 6, 7</span> pu.
```

(tab-sources-transport)=

**

<table>
<thead>
<tr>
<th style="text-align: left;"></th>
<th style="text-align: left;"></th>
<th colspan="5" style="text-align: center;"><span class="math inline">‖<strong>Y</strong><em>e</em><sub><em>m</em><em>o</em><em>d</em></sub>‖</span> [pu] a</th>
<th style="text-align: right;"><span class="math inline">1/‖<strong>Y</strong><em>e</em><sub><em>m</em><em>o</em><em>d</em></sub>‖</span></th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><span>3-7</span>(lr)<span>8-8</span> Equipement</td>
<td style="text-align: left;">Reglage</td>
<td style="text-align: right;">0.01 Hz</td>
<td style="text-align: right;">0.1 Hz</td>
<td style="text-align: right;">1 Hz</td>
<td style="text-align: right;">10 Hz</td>
<td style="text-align: right;">100 Hz</td>
<td style="text-align: right;">a 0,1 Hz [pu]</td>
</tr>
<tr>
<td style="text-align: left;">SMTR</td>
<td style="text-align: left;">machine + transfo, AVR</td>
<td style="text-align: right;">6.47</td>
<td style="text-align: right;">6.40</td>
<td style="text-align: right;">3.30</td>
<td style="text-align: right;">2.57</td>
<td style="text-align: right;">1.91</td>
<td style="text-align: right;">0.156</td>
</tr>
<tr>
<td style="text-align: left;">GFM1 (nq=4%)</td>
<td style="text-align: left;">GFM, statisme de tension 4 %</td>
<td style="text-align: right;">5.24</td>
<td style="text-align: right;">5.29</td>
<td style="text-align: right;">6.39</td>
<td style="text-align: right;">7.09</td>
<td style="text-align: right;">3.11</td>
<td style="text-align: right;">0.189</td>
</tr>
<tr>
<td style="text-align: left;">GFM1 (nq=1/Ka)</td>
<td style="text-align: left;">GFM, tension aussi raide que l'AVR</td>
<td style="text-align: right;">6.50</td>
<td style="text-align: right;">6.51</td>
<td style="text-align: right;">6.63</td>
<td style="text-align: right;">7.09</td>
<td style="text-align: right;">3.11</td>
<td style="text-align: right;">0.154</td>
</tr>
<tr>
<td style="text-align: left;">GFM3 (nq=4%)</td>
<td style="text-align: left;">GFM a bus continu dynamique</td>
<td style="text-align: right;">5.24</td>
<td style="text-align: right;">5.29</td>
<td style="text-align: right;">6.39</td>
<td style="text-align: right;">7.09</td>
<td style="text-align: right;">3.11</td>
<td style="text-align: right;">0.189</td>
</tr>
<tr>
<td style="text-align: left;">GFL1</td>
<td style="text-align: left;">GFL boucles P/Q</td>
<td style="text-align: right;">0.50</td>
<td style="text-align: right;">0.50</td>
<td style="text-align: right;">0.49</td>
<td style="text-align: right;">0.24</td>
<td style="text-align: right;">0.27</td>
<td style="text-align: right;">1.988</td>
</tr>
<tr>
<td style="text-align: left;">GFL2</td>
<td style="text-align: left;">GFL boucle Vdc/Q</td>
<td style="text-align: right;">10.38</td>
<td style="text-align: right;">10.38</td>
<td style="text-align: right;">10.17</td>
<td style="text-align: right;">4.44</td>
<td style="text-align: right;">0.35</td>
<td style="text-align: right;">0.096</td>
</tr>
</tbody>
</table>

**La machine et le grid-forming ne se distinguent plus en basse fréquence.** Tous deux plafonnent à l'inverse de leur impédance de raccordement : 6,40 pu pour la machine et 6,51 pu pour le GFM à tension raide, soit $|Z| = 0{,}155$ pu, le transformateur élévateur lui-même. Le GFM à statisme de tension de 4 % descend à 5,29 pu, soit 0,19 pu : son statisme s'ajoute simplement à la réactance de couplage, $|Z| \approx X_T + n_q$. C'est la seule différence en régime lent, et c'est un choix de réglage, pas une propriété de la technologie.

**Ils se séparent en revanche entre 0,5 et 10 Hz.** La machine y creuse jusqu'à 1,7 pu vers 2,6 Hz, soit $|Z| = 0{,}58$ pu : son AVR ne suit plus et il ne reste que ses réactances transitoires. Le GFM, dont la boucle de tension est réglée à 15 ms, ne descend jamais sous 5,6 pu. La vraie différence entre les deux n'est donc pas la raideur statique, mais la *bande dans laquelle elle est tenue*.

**Le grid-following reste à 0,5 pu**, treize fois plus mou, et cela ne dépend d'aucun réglage. Avec la même boucle de courant à 1 ms que le GFM, l'écart à 34 Hz reste de deux décades : le contraste en moyenne fréquence est bien structurel, pas un artefact de bande passante.

:::{admonition} Deux chiffres à ne pas surinterpréter
:class: warning

**Les 91 pu de la figure {numref}`fig-sources`** tiennent au filtre de 0,01 pu des réglages par défaut. Rapporté à un raccordement réaliste de 0,15 pu, le GFM tombe à 6 pu, du même ordre que la machine : l'écart de deux décades entre GFM et machine de la figure {numref}`fig-sources` était surtout un écart d'impédance de raccordement.

**Les 10,4 pu du `GFL2`** ne sont pas une raideur de source de tension. Ils viennent de sa boucle de tension continue et dépendent fortement de son réglage : en passant $\tau_{VL}$ de 100 à 20 puis 5 ms, le gain tombe à 2,9 puis 1,1 pu. C'est une conséquence du couplage entre puissance et tension continue, à examiner au cas par cas.
:::

## Ce que cela implique

Une source de tension s'oppose aux variations de tension du réseau : c'est elle qui tient le réseau et qui permet à un îlot d'exister. Une source de courant subit la tension : elle a besoin d'une référence extérieure, et sa PLL doit pouvoir s'y accrocher. Un réseau entièrement composé de grid-following n'a personne pour tenir sa tension -- c'est l'argument central en faveur d'une part minimale de grid-forming, et c'est aussi pourquoi le réseau `wscc9_gfm.ipynb` garde une machine synchrone au nœud slack.

Ces courbes sont également la matière première de l'**analyse de stabilité par impédances** : le système bouclé équipement/réseau est stable si le rapport $\bm{Z}_{r\acute{e}seau}\,\bm{Y}_{\acute{e}quipement}$ satisfait le critère de Nyquist généralisé. Le mode de jauge et l'analyse modale du chapitre {ref}`ch-exemples` donnent la même information par les valeurs propres du modèle global ; l'approche fréquentielle a l'avantage de se mesurer sur un équipement réel, sans modèle de son intérieur.

## Outils

Côté Python, les fonctions sont dans `g2egrid.analyse.frequence` : `admittance_dq`, `impedance_dq`, `axes_tension`, `gain`, `vecteur_complexe`. Les figures de ce chapitre sont produites par le notebook `bode_equipements.ipynb` du dossier `python/etudes/`. Côté MATLAB, le script `etudes/bode_equipements.m` de V2 refait le même calcul avec les classes de V2 à partir de `getLinearModel`, trace les mêmes courbes et exporte un CSV que le notebook superpose aux siennes.

```python
from g2egrid.analyse.frequence import admittance_dq, axes_tension, gain

f = np.logspace(-2, 4, 400)          # Hz
Y = admittance_dq(el, f)             # (400, 2, 2), pu base systeme
e_mod, e_ang = axes_tension(el)      # directions module et angle
raideur = gain(Y, e_mod)             # ||Y e|| dans la direction du module
```

**Accord entre les deux versions.** Sept équipements sur huit coïncident à $10^{-6}$ près, terme à terme sur toute la bande. La `SMTR` fait exception en très basse fréquence, jusqu'à 36 % d'écart sous 0,1 Hz : `SMTR.m` rejoue l'équilibre de V1, qui ne vérifie pas ses propres équations (section {ref}`sec-v1-equilibres`), alors que Python raffine le sien par Newton. Les deux outils linéarisent donc autour de deux points voisins mais différents, et le gain de l'AVR amplifie l'écart. En posant `init='symSG'`, Python rejoue la convention de V1 et l'écart retombe à $10^{-6}$ : la différence est bien dans le point d'équilibre, pas dans les équations.
