<!-- Page GENEREE depuis documentation/chapitres/11_validation.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ch-validation)=
# Validation contre G2Elin V1
V1 étant la référence, les deux nouvelles versions ont été confrontées à elle sur les trois réseaux élémentaires de `ref_SMIB.m`. Ce chapitre décrit la méthode, les résultats, puis les conventions de V1 qu'il a fallu comprendre pour y parvenir et les défauts de V2 que la comparaison a corrigés. Les outils se trouvent dans `python/developpement/comparaison_matlab/`.

## Méthode

Les trois outils construisent le même réseau à partir des mêmes paramètres : `ref_SMIB.m` pour V1, `SMIB_V1.m` pour V2, `smib_v1.py` pour Python. Les exports MATLAB sont produits par `export_ss.m` et `export_v2.m`, qui écrivent la matrice $\bm{A}$, les noms d'états et, pour V1, la solution du load flow et les structures d'initialisation. Le script `comparer_trois.py` compare ensuite les matrices à deux niveaux.

**Le spectre** est invariant par changement de base d'état et se compare toujours. Les modes sont appariés par affectation optimale (algorithme de Kuhn-Munkres) sur la distance relative dans le plan complexe. Un appariement par plus proche voisin fabrique de faux écarts quand deux modes réels d'un outil font face à une paire complexe de l'autre.

**La matrice $\bm{A}$ terme à terme** n'a de sens que si les trois outils ont le même jeu d'états dans le même repère. Leurs lignes et colonnes sont remises dans l'ordre canonique de V1 : bus infini, équipement, nœuds, ligne, charge. L'écart maximal est normalisé par $\max|A_{ij}|$ de la référence, et l'on compte les coefficients qui en diffèrent de plus de $10^{-6}$ en relatif.

## Résultats

:::{table} Écarts à V1 sur les réseaux élémentaires. Aucun coefficient de $\bm{A}$ ne diffère de plus de $10^{-6}$ en relatif.
:name: tab-validation

| Préréglage | Comparaison | États | Spectre méd. | Spectre max | Modes $<10^{-3}$ | Écart $\bm{A}$ |
|:---|:---|---:|---:|---:|---:|---:|
| SM | V1 $\to$ V2 | 30 | $1{,}1\cdot10^{-12}$ | $2{,}5\cdot10^{-6}$ | 30/30 | $1{,}4\cdot10^{-7}$ |
| SM | V1 $\to$ Python | 30 | $2{,}8\cdot10^{-15}$ | $2{,}5\cdot10^{-6}$ | 30/30 | $4{,}9\cdot10^{-11}$ |
| GFM | V1 $\to$ V2 | 26 | $1{,}6\cdot10^{-13}$ | $6{,}2\cdot10^{-9}$ | 26/26 | $6{,}7\cdot10^{-9}$ |
| GFM | V1 $\to$ Python | 26 | $1{,}3\cdot10^{-15}$ | $3{,}3\cdot10^{-12}$ | 26/26 | $5{,}5\cdot10^{-15}$ |
| GFL | V1 $\to$ V2 | 25 | $8{,}3\cdot10^{-16}$ | $1{,}9\cdot10^{-14}$ | 25/25 | $2{,}9\cdot10^{-14}$ |
| GFL | V1 $\to$ Python | 25 | $5{,}2\cdot10^{-15}$ | $1{,}5\cdot10^{-14}$ | 25/25 | $3{,}1\cdot10^{-14}$ |
:::

Les trois versions produisent la même matrice $\bm{A}$ : à la précision machine pour quatre comparaisons sur six, à sept chiffres significatifs pour la machine côté V2. Le seul écart encore lisible, $2{,}5\cdot 10^{-6}$ sur le spectre de la machine, porte sur deux pôles réels de la référence, $-0{,}5154$ et $-0{,}7787$, sur le point de fusionner : dans cette configuration l'écart relatif amplifie une perturbation infime de $\bm{A}$.

## Ce qu'il a fallu comprendre de V1

Aucun des écarts initiaux ne venait des paramètres. Tous tenaient à des conventions de V1, dont plusieurs sont des incohérences qu'il faut rejouer pour retrouver sa matrice. La démarche qui a fonctionné est toujours la même : exporter les valeurs intermédiaires de V1 et les comparer nombre à nombre.

### Le load flow a quatre nœuds

`network_form.m` ne place pas l'équipement sur le nœud du réseau : il lui crée un nœud à lui, derrière le transformateur de raccordement, et c'est là que vont la consigne de tension et le type de bus de `Y_DER`. Le réseau réellement résolu par V1 compte donc quatre nœuds :

:::

Poser directement $V = 1$ sur les nœuds 1 et 2 donne un tout autre point de fonctionnement : sur une ligne de $2{,}5 \cdot 10^{-4}$ pu, imposer deux modules égaux force un transit réactif de 0,4 pu et fait tourner l'angle rotor de plus de trente degrés. Les scripts de comparaison résolvent ce load flow à quatre nœuds, puis imposent sa solution au réseau réduit qui porte le modèle petit-signal.

### Deux conventions de capacité de ligne

Le load flow de V1 place $B/2$ à chaque extrémité de ligne, alors que son modèle de nœud applique la capacité entière à chaque nœud. La moitié qui les sépare vaut $2 \cdot 10^{-4}$ pu de réactif, exactement l'écart qui apparaissait sur le courant réactif de l'équipement.

### L'angle du repère global

$\theta_{g,0}$ est l'angle interne de l'élément de référence, et non celui du bus slack (section {ref}`sec-reperes`). Utiliser l'angle du bus décalait la rotation de tous les équipements de 0,7 degré.

(sec-v1-equilibres)=
### Des équilibres qui ne vérifient pas leurs équations
**La machine.** `ref_SMIB.m` prend le courant du load flow pour celui de l'armature, si bien que le transformateur ne porte que $I - V/R_g$, soit 0,1 pu de moins que le courant qui traverse réellement cette branche. Il ramène aussi la tension au repère rotorique par $e^{j(\alpha - \theta)}$ alors que $v_0$ porte déjà $e^{j\alpha}$ : la tension de borne se retrouve en avance d'un angle $\alpha$ de trop, 1,4 degré ici.

**Le grid-following.** Les intégrales des correcteurs sont divisées par leur propre gain : $M_{d,0} = I_{s,d,0}/K_{id}$ avec $K_{id} = -77\,783$, donc $M_{d,0} \approx 0$ au lieu de $I_{s,d,0}$. C'est l'initialisation d'un schéma-bloc où le gain suit l'intégrateur, transposée à un modèle d'état où il le précède.

**Conséquence.** Raffiner ces équilibres par Newton s'en écarte. Les modes `init='symSG'`, `init='symGFL'`, `PointInit` et `InitV1` rejouent ces formules sans raffinement.

(sec-v1-algebriques)=
### Les variables algébriques ne sont pas une commodité d'écriture
C'est le point qui a fermé les derniers écarts. `symSG` déclare algébriques les courants de la machine et sa tension de borne ; `symGFL` déclare algébriques les rapports cycliques et la pulsation de PLL. V1 substitue dans ses jacobiennes ses propres valeurs d'équilibre pour ces grandeurs, et non celles que rendraient leurs équations au même point. Ses équilibres étant incohérents, les deux diffèrent.

:::{table} Valeurs substituées par V1 et valeurs rendues par leurs équations.
:name: tab-algebriques

| Grandeur | Substituée par V1 | Rendue par son équation |
|:---|---:|---:|
| courant $i_q$ de la machine | 0,37585 | 0,30284 |
| tension de borne $v_{ed}$ de la machine | 0,64225 | 0,66050 |
| pulsation $\omega_{pll}$ du GFL | 1,0 | 0,95638 |
:::

Écrire ces grandeurs comme des expressions revient à corriger la référence, et fait apparaître 3 % d'écart sur la boucle mécanique, sur celle du régulateur de tension et sur les coefficients en $v_{dc}$. Les garder algébriques les reproduit exactement.

## Défauts de V2 corrigés par la comparaison

**Équilibre des convertisseurs posé au mauvais endroit.** Les six classes de convertisseurs imposaient que $p$ et $q$, mesurés aux bornes du condensateur, valent les consignes du load flow, qui sont celles du nœud. Entre les deux, le transformateur de raccordement consomme 0,006 pu de réactif : le convertisseur injectait 0,0187 pu de courant réactif au lieu de 0,0127. Elles contraignent désormais le courant injecté au nœud. Au passage, les modes nuls de WSCC et de Kundur sont retombés de $10^{-4}$ à $10^{-10}$.

**Boucle de tension continue du GFM recalée d'office.** `GFM3.m` réglait sa boucle continue sur la conductance effective $G_{\mathrm{eff}} = G_{dc} - P_0$, là où V1 la règle sur $G_{dc}$. Un `tau_DC` vide désactive ce recalage.

**Boucle de tension continue du GFL deux fois trop molle.** `GFL2.m` la réglait en $3\tau/t_r$ au lieu de $2\zeta\omega_n\tau - 1$, soit un gain proportionnel deux fois plus petit, et forçait l'intégrale à zéro. L'écart maximal sur le spectre valait alors 1,00, soit un mode entièrement déplacé.

**Machine à transformateur intégré absente.** `SM.m` n'a que 17 états et laisse son transformateur dehors. La classe `SMTR` ajoute la structure de `symSG`, sans laquelle la comparaison terme à terme était impossible.

**Résistance de causalité parasite.** Déclarer la machine sur un nœud terminal ajoutait la résistance de causalité de 100 pu en parallèle du shunt de la machine, et décalait de 12 % les pôles rapides du transformateur.

(sec-v2-multimachine)=
## Réseaux multi-machines : réconcilier V2 et Python
Les trois réseaux ci-dessus ont un *bus infini* pour référence d'angle. Sur un réseau où la référence est une *machine synchrone*, les deux versions donnaient au contraire des modes lents différents : 1,5 % d'écart sur le WSCC à trois machines, et jusqu'à 22 % dès qu'un GFM entrait en jeu (mode électromécanique à 14,13 rad/s côté Python contre 12,47 côté V2). Les équipements pris séparément coïncidaient pourtant à $10^{-5}$ près (chapitre {ref}`ch-sources`) : l'écart venait de l'assemblage et des points de fonctionnement, et quatre défauts distincts s'y cachaient.

### Les quatre défauts

**Le shunt $R_g$ compté deux fois (les deux versions).** Le courant injecté par une machine vaut déjà $i - v/R_g$ : le shunt est interne au modèle. Les deux outils l'ajoutaient pourtant au $Y$-bus du load flow. La loi des nœuds du modèle dynamique était alors violée de $V/R_g$ -- 0,01 pu, soit 1 MW sur une base de 100 MVA -- au bus de chaque machine, et le point de linéarisation n'était pas un équilibre. Le shunt a été retiré des deux $Y$-bus ; le test `test_loi_des_noeuds_avec_machine_de_reference` vérifie désormais la loi des nœuds au bus de chaque machine.

**Un repère global incohérent (Python).** Les convertisseurs se repéraient par $\theta - \theta_g$, c'est-à-dire par rapport au rotor de la machine de référence, mais les nœuds, lignes et charges écrivaient leurs couplages $dq$ à la pulsation nominale, et la machine tournait ses grandeurs de son seul angle $\theta$. Trois repères pour un même réseau : sans conséquence avec un bus infini, dont l'angle ne bouge pas, mais faux avec une machine. Le repère suit désormais la référence partout (section {ref}`sec-reperes`), et le test `test_spectre_independant_de_la_machine_de_reference` le contrôle.

**La vitesse du slack branchée comme un courant (V2).** Dans `Node.m`, les colonnes de $\bm{B}$ étaient rangées dans l'ordre $[\omega, i_{gd}, i_{gq}, i_{sh,d}, i_{sh,q}]$ alors que les étiquettes annonçaient $[i_{gd}, i_{gq}, \omega_{slack}, i_{sh,d}, i_{sh,q}]$, et l'assemblage câble par étiquette. La vitesse du slack arrivait donc sur la colonne du courant d'axe $q$, avec un coefficient $\omega_b/C$ au lieu du couplage de repère $-\omega_b v_{d,0}$. Là encore, invisible avec un bus infini.

**L'angle rotor du slack calculé sans $R_g$ (V2).** `SM.m` figeait l'angle du rotor de la machine slack à $\pi/2 - \arg(E_q)$ avec $E_q = V + (R_a + jL_q)I_{LF}$, sans le courant $V/R_g$ du shunt. L'angle était décalé de $L_q|V|/R_g$, soit 0,6 degré, et `fsolve` compensait sur la tension d'interface. La formule retenue est celle du portage : $E_q = V(1 + jL_q/R_g) + (R_a + jL_q)I_{LF}$.

### Ce qui reste, et pourquoi

Une dernière différence subsiste, mais c'est une convention assumée : V2 calcule une **fréquence système** à partir du statisme collectif des groupes (`Network.SystemOmega_pu`), en prenant pour consigne le $P_g$ du slack déclaré dans `Y_DER`. Ce $P_g$ valant zéro dans les réseaux de test, tout ce que produit le slack est lu comme un déséquilibre : V2 linéarise alors à 0,9993 pu là où Python linéarise à 1,00. Il suffit de déclarer la production réelle du slack pour que les deux coïncident.

:::{table} Écart relatif entre V2 et le portage Python sur les modes lents (sous 5 Hz), avant et après les corrections.
:name: tab-v2-multimachine

| Réseau | Avant | Après | Après, à $\omega_{sys} = 1$ |
|:---|---:|---:|---:|
| WSCC 9 bus (3 machines) | $1{,}5\cdot 10^{-2}$ | $4{,}6\cdot 10^{-4}$ | $8{,}2\cdot 10^{-8}$ |
| WSCC 9 bus à GFM | $2{,}2\cdot 10^{-1}$ | $9{,}7\cdot 10^{-5}$ | -- |
| WSCC 9 bus + CIGRE | $1{,}8\cdot 10^{-1}$ | $1{,}1\cdot 10^{-4}$ | -- |
:::

Les corrections ne touchent pas les trois réseaux élémentaires de V1 : le tableau « {ref}`Écarts à V1 sur les réseaux élémentaires <tab-validation>` » reste inchangé, avec zéro coefficient différent. Le $Y$-bus de V2 garde d'ailleurs une option, `ResistanceCausaliteLF`, que `SMIB_V1.m` met à faux pour rejouer le load flow de V1, qui ignore les résistances de causalité des nœuds.

Enfin, la séparation transport / distribution du chapitre {ref}`ch-td` donne désormais les mêmes chiffres dans les deux outils, jusqu'à la troisième décimale des tableaux « {ref}`Modes lents du transport <tab-td-modes>` » et « {ref}`Écart maximal sur les modes lents du transport, selon l'ordre et la méthode de réduction du modèle de distribution (146 états à l'origine) <tab-td-reduction>` » : c'est la vérification la plus complète dont on dispose sur un réseau de 230 états.
