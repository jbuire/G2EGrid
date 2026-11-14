<!-- Page GENEREE depuis documentation/chapitres/12_atelier.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ch-atelier)=
# L'atelier : deux voies, AC/DC, zones, interface
Les chapitres précédents décrivent un modèle : des équations, une matrice $\bm{A}$, et des scripts qui la construisent. Celui-ci décrit ce qui a été bâti autour, et qui tient en une idée : *un même réseau, décrit une seule fois, doit pouvoir être résolu de deux manières, comporter du continu, se composer par zones, et se manipuler aussi bien en code qu'à la souris*.

Chacune de ces quatre extensions aurait pu se faire en dupliquant quelque chose --- une seconde description du réseau pour le non linéaire, une seconde table de colonnes pour l'interface, une seconde copie de la distribution pour chaque départ. C'est précisément ce qu'il a fallu éviter, et la plupart des défauts racontés ici viennent d'un endroit où la duplication avait déjà eu lieu sans qu'on s'en aperçoive.

## Deux voies sur le même câblage

Un élément de `g2egrid` déclare un système différentiel-algébrique : $f$ pour ses états, $g$ pour ses variables algébriques, $h$ pour ses sorties. Le *câblage* referme l'ensemble en une boucle fermée, sous forme LFT : 

$$
\bm{U} = \bm{R}_{\text{uext}}\,\bm{u}_{\text{ext}} + \bm{R}_y\,\bm{Y},
$$

 les tensions de port venant des nœuds et les courants de port y retournant par la loi des nœuds.

Jusqu'ici, ce câblage vivait à l'intérieur de `build_state_space` et n'existait que linéarisé. Il en a été extrait --- classe `Cablage`, méthode `Network.cablage()` --- pour que **les deux voies partent exactement de la même description** :

- la voie **linéaire** élimine les variables algébriques par complément de Schur et rend $\bm{A}$ ;

- la voie **non linéaire** (module `temporel`) garde $f$, $g$, $h$ et le câblage tels quels, et résout à chaque pas le système couplé 

$$
\begin{bmatrix} g(x,z,u) \\ u - \bm{R}_{\text{uext}}\,u_{\text{ext}} - \bm{R}_y\,h(x,z,u)\end{bmatrix} = \bm{0}
$$

 en $(z,u)$ par un Newton à jacobien analytique, avant d'intégrer par Radau.

Les jacobiens sont *analytiques* : chaque élément déclare $f$, $g$, $h$ en SymPy, qui les dérive une fois pour toutes ; `lambdify` (avec élimination des sous-expressions communes) en fait des fonctions NumPy, compilées une fois par classe d'élément. Le système couplé, lui, est très creux -- 144 inconnues et un jacobien plein à 1,1 % sur `wscc9` : chaque élément ne voit que ses variables, et $\bm{R}_y$ n'est qu'une matrice d'incidence. Il est donc assemblé creux, sur une structure calculée une fois, et résolu par un Newton *modifié* : la factorisation LU creuse (`splu`) est gardée d'une résolution à la suivante, et refaite seulement quand la convergence ralentit. Le dogleg de MINPACK, en matrice pleine, ne sert plus que de repli -- et aucune réponse n'est acceptée sans que le résidu soit vérifié nul, car le « succès » de MINPACK ne porte que sur la taille du pas. Sur `wscc9`, une simulation EMT va ainsi 4 à 5 fois plus vite qu'avec MINPACK seul, pour un écart aux trajectoires inférieur à l'erreur d'intégration.

L'intérêt n'est pas d'avoir une simulation : c'est d'avoir *deux chemins qui ne partagent aucune ligne de code de résolution*. Leur accord ne peut donc pas venir d'une erreur commune, et c'est ce qui rend la validation du §{ref}`sec-atelier-validation` probante.

:::{admonition} La matrice $\bm{A}$ instantanée
:class: note

Le jacobien du système couplé donne, par le théorème des fonctions implicites, la matrice $\bm{A}$ *à l'instant courant*. Fournie à Radau, elle fait passer une seconde simulée de 34 s à 0,6 s de calcul sur un modèle à 61 états. Exposée à l'utilisateur sous le nom `linearisation(m, tr, t)`, elle répond à une question que le modèle linéarisé à l'équilibre ne pose pas : *ce mode est-il déjà dangereux pendant le transitoire*, quand l'angle rotor et la tension se sont éloignés de leur valeur nominale.
:::

(sec-atelier-validation)=
## Valider le petit signal par l'EMT
Trois épreuves, de coût croissant et de portée différente.

(point-fixe.)=
##### Point fixe.
Au point de fonctionnement, $\dot{x}$ doit être nul et les variables algébriques doivent déjà vérifier leurs équations. Sinon les points de fonctionnement calculés élément par élément ne forment pas un point fixe du *réseau* : le défaut est invisible en modal --- on linéarise quand même, autour d'un point qui n'existe pas --- et fatal en temporel, où la trajectoire dérive sans qu'on ait rien perturbé.

(le-carré-de-lamplitude.)=
##### Le carré de l'amplitude.
L'écart entre réponse linéaire et réponse non linéaire doit croître comme $a^2$ : c'est la définition d'une linéarisation au premier ordre. Doubler l'amplitude doit donc quadrupler l'écart. Un rapport proche de $2$ dirait que l'écart est du premier ordre, c'est-à-dire que $\bm{A}$ n'est pas la bonne matrice.

(identification-modale.)=
##### Identification modale.
On identifie les modes *directement sur la trajectoire non linéaire* --- matrix pencil et [era]{.smallcaps}, module `identification` --- et on les apparie aux valeurs propres de $\bm{A}$ par affectation optimale.

:::{table} Les trois épreuves, échelon de 0,01 pu sur une consigne de puissance. Le rapport attendu est exactement $4$.
:name: tab-atelier-validation

| Réseau  |     $\max\vert \dot{x}\vert $ |     Écart algébrique | Rapport $a^2$ | Résidu ident. |
|:--------|--------------------:|---------------------:|--------------:|--------------:|
| `smib`  | $1{,}5\cdot10^{-7}$ | $1{,}2\cdot10^{-11}$ |     $4{,}012$ |      0,0127 % |
| `wscc9` | $6{,}9\cdot10^{-9}$ | $1{,}1\cdot10^{-12}$ |     $3{,}995$ |      0,0125 % |
:::

(tab-atelier-identification)=

**

<table>
<thead>
<tr>
<th colspan="2" style="text-align: center;">Prédit</th>
<th colspan="2" style="text-align: center;">Identifié</th>
<th colspan="2" style="text-align: center;">Écart</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: right;"><span>1-2</span>(lr)<span>3-4</span>(l)<span>5-6</span> <span class="math inline"><em>f</em></span> [Hz]</td>
<td style="text-align: right;"><span class="math inline"><em>ζ</em></span> [%]</td>
<td style="text-align: right;"><span class="math inline"><em>f</em></span> [Hz]</td>
<td style="text-align: right;"><span class="math inline"><em>ζ</em></span> [%]</td>
<td style="text-align: right;"><span class="math inline"><em>f</em></span></td>
<td style="text-align: right;"><span class="math inline"><em>ζ</em></span></td>
</tr>
<tr>
<td style="text-align: right;"><span class="math inline">2, 2000</span></td>
<td style="text-align: right;"><span class="math inline">8, 696</span></td>
<td style="text-align: right;"><span class="math inline">2, 2000</span></td>
<td style="text-align: right;"><span class="math inline">8, 679</span></td>
<td style="text-align: right;">+0<span>,</span>00 %</td>
<td style="text-align: right;">-0<span>,</span>20 %</td>
</tr>
<tr>
<td style="text-align: right;"><span class="math inline">2, 0115</span></td>
<td style="text-align: right;"><span class="math inline">11, 251</span></td>
<td style="text-align: right;"><span class="math inline">2, 0130</span></td>
<td style="text-align: right;"><span class="math inline">11, 304</span></td>
<td style="text-align: right;">+0<span>,</span>08 %</td>
<td style="text-align: right;">+0<span>,</span>47 %</td>
</tr>
</tbody>
</table>

:::{admonition} Le résidu de reconstruction ne prévient pas
:class: warning

L'identification rend un résidu, qui mesure à quel point les modes trouvés reconstruisent le signal fourni. Il est tentant d'y voir un indicateur de confiance. **Il ne l'est pas**, et c'est la mesure la plus utile de ce travail.

Sur une fenêtre de 0,6 s au lieu de 4 s, le résidu vaut 0,0002 % --- dix fois *mieux* que sur la fenêtre correcte --- pendant que l'amortissement identifié est faux de 56 %. Sur une demi-période, une exponentielle amortie et une exponentielle plus raide de résidu plus grand sont indiscernables. Le même piège se referme sur une seule voie : `SM1.Dwr` seule reconstruit à 0,0091 % et rend le mode 2,2 Hz amorti à 17,4 % au lieu de 8,7 %, parce qu'elle en est presque un nœud.

Le résidu mesure l'adéquation au signal *fourni*, pas la pertinence de ce qu'on a choisi de fournir. Les remèdes sont d'allonger la fenêtre et d'identifier sur plusieurs voies ; l'interface affiche les deux avertissements.
:::

## Le continu

Un port alternatif porte deux signaux ($v_d, v_q$ en entrée, $i_d, i_q$ en sortie), un port continu un seul. Toute l'extension tient dans cette asymétrie, et dans la règle qui l'a encadrée : *le chemin alternatif existant ne devait pas bouger d'un bit*.

### Le bus continu est le condensateur du convertisseur

C'est le point de conception central. Le modèle `GFM3` portait ses propres états $v_{dc}$ et $i_{dc}$ ; sa variante `GFM3DC` les perd (15 états $\rightarrow$ 13), lit $v_{dc}$ comme une *entrée de port* et rend $i_{dc}$ en sortie. La loi des nœuds du bus redonne alors, à l'identique, l'équation qui était dans le modèle.

L'équivalence se vérifie exactement. En écrivant les deux équations continues de `GFM3`, 

$$
\frac{\mathrm{d}v_{dc}}{\mathrm{d}t} = \frac{\omega_b}{C_{dc}}\bigl(i_{dc} - G_{dc}v_{dc} - m_d i_{sd} - m_q i_{sq}\bigr),
  \qquad
  \frac{\mathrm{d}i_{dc}}{\mathrm{d}t} = \frac{i_{dc}^{\text{ref}} - i_{dc}}{T_{dc}},
$$

 avec $i_{dc}^{\text{ref}} = p^{\text{ref}}/v_{dc}^{\text{ref}} + K_{p,dc}(v_{dc}^{\text{ref}} - v_{dc})$, la première est reproduite telle quelle par un bus continu de capacité $C_{dc}$. La seconde, en factorisant $K_{p,dc}$, devient 

$$
\frac{\mathrm{d}i_{dc}}{\mathrm{d}t} = \frac{K_{p,dc}}{T_{dc}}
  \Bigl[\underbrace{v_{dc}^{\text{ref}} + \frac{p^{\text{ref}}}{K_{p,dc}\,v_{dc}^{\text{ref}}}}_{E}
  - v_{dc} - \underbrace{\tfrac{1}{K_{p,dc}}}_{R}\,i_{dc}\Bigr],
$$

 soit exactement une f.é.m. derrière R-L avec 

$$
R = \frac{1}{K_{p,dc}}, \qquad X = \frac{\omega_b\,T_{dc}}{K_{p,dc}}, \qquad
  E = v_{dc}^{\text{ref}} + \frac{p^{\text{ref}}}{K_{p,dc}\,v_{dc}^{\text{ref}}}.
$$

Mesuré : `GFM3DC` + `BusDC`($C = C_{dc}$) + `SourceDC` a le même spectre que `GFM3` monolithique --- 24 états de part et d'autre, écart apparié $2{,}9\cdot10^{-10}$ sur un spectre allant à $3{,}2\cdot10^{5}$, soit $5{,}4\cdot10^{-12}$ en relatif. Le point de linéarisation coïncide à quinze chiffres.

:::{admonition} Ce qui n'est pas transposable
:class: note

`GFL2DC` n'est délibérément pas comparé à `GFL2`. Le côté continu de `GFL2` est une source de *courant* filtrée, sans réaction à $v_{dc}$, là où `SourceDC` est une f.é.m. derrière R-L. Seul le grid-forming avait une boucle exactement transposable.
:::

### Le load flow AC/DC

Un seul Newton sur $[\theta_{\text{AC}},\; V_{\text{AC}},\; V_{\text{DC}},\;
P_{\text{ac}}$ des [vsc]{.smallcaps} en mode $V_{dc}]$. Le bus alternatif d'un convertisseur garde le type que lui donne `Y_DER` ; en mode $V_{dc}$, sa puissance programmée devient une inconnue *additive*, si bien que le bloc alternatif/alternatif du jacobien est celui d'avant, copié sans modification, plus une colonne $-1$. C'est ce qui garantit la non-régression.

Deux modes s'ajoutent aux quatre de départ pour les réseaux continus maillés : `Vtheta`, où la station forme seule son île alternative et en est le slack --- sa puissance $P_{\text{bus}} - p_{\text{autres}}$ entre dans le bilan continu par la règle de chaîne $\partial P_{dc}/\partial p \cdot \partial P/\partial x$ ---, et `droop`, où $p_{ac} = p_{ref} + s\,(v_{dc} - v_{dc,ref})$ apporte une colonne $-s$ à la ligne $P$ du bus et un terme $s\,\partial P_{dc}/\partial p$ au bilan continu. Une île continue peut avoir plusieurs références, à des nœuds distincts ({doc}`../guide/cigre_b4`).

Le couplage n'invente aucun coefficient : $P_{dc} = -V_{dc}\,i_{\text{port}}$, avec le $i_{\text{port}}$ du *modèle dynamique*, dont $m_d, m_q, i_{sd}, i_{sq}$ sortent des formules d'initialisation au point alternatif courant. Les pertes du convertisseur sont donc celles de ses équations, et le point de fonctionnement AC/DC est un point fixe du modèle non linéaire --- résidu $3{,}1\cdot10^{-11}$ sur `hvdc2t`, $2{,}1\cdot10^{-11}$ sur `mtdc3t`.

:::{admonition} Deux termes qu'il ne faut pas oublier
:class: warning

Deux erreurs de conception ont été trouvées à l'écriture du jacobien, et toutes deux donnaient un résultat *presque* juste.

La décomposition par blocs naïve met $0$ sous $\theta_{\text{AC}}$. C'est faux dès qu'un grid-forming est en jeu : sur un bus [pv]{.smallcaps}, le réactif du convertisseur n'est pas programmé, il suit le load flow, et il manque donc la règle de chaîne $\partial P_{dc}/\partial q \cdot \partial Q/\partial x$.

De même, écrire les pertes $R_1|i_s|^2 + G_{dc}V_{dc}^2$ oublie le $R_2|i_g|^2$ de la branche située entre le condensateur, où le modèle mesure $p$, et le nœud réseau. Sans ce terme, le bilan de puissance boucle à $10^{-3}$ ; avec, à $4\cdot10^{-15}$.
:::

(non-régression.)=
##### Non-régression.
Les quatorze réseaux alternatifs du catalogue rendent un load flow *identique bit à bit* --- `lf.V.tobytes()` et nombre d'itérations --- vérifié contre une extraction propre de `git archive`. La garantie a été délibérément déplacée du texte vers le résultat : geler le source de `run_loadflow` aurait interdit d'y écrire la garde `if not self.bus_specs_dc`, qui est justement la bonne façon de ne rien changer.

## Les zones

Un réseau de transport auquel se raccordent trois distributions identiques s'écrivait à plat : 54 bus numérotés à la main, la même distribution recopiée trois fois à trois numérotations différentes. Une *zone* se décrit une fois, en numérotation locale et *sans slack*, et s'*instancie* --- chaque instance pouvant représenter $K$ départs identiques, les puissances multipliées par $K$ et les impédances divisées par $K$.

Le transformateur de poste appartient à l'*instance* et non à la zone : sa puissance dépend de $K$ et son primaire du bus amont. Les transformateurs de départ, eux, appartiennent à la zone.

:::{table} `wscc3sm_cigre3` réécrit en zones, comparé au notebook à plat. L'écart n'est pas « sous $10^{-12}$ » : il est exactement nul.
:name: tab-atelier-zones

| Comparaison                                          |                Écart |
|:-----------------------------------------------------|---------------------:|
| Tables, ligne à ligne (18 tables, 54 bus)            |   0 ligne différente |
| Éléments assemblés (99, plus 54 nœuds)               |           même ordre |
| États (429 de part et d'autre)                       |           mêmes noms |
| Matrice $\bm{A}$, pour $\max\vert A\vert  = 6{,}62\cdot10^{6}$ | $0{,}000\cdot10^{0}$ |
| $\mathrm{Re}(\lambda)$                               | $0{,}000\cdot10^{0}$ |
:::

Le zéro exact n'est pas un heureux hasard, c'est une contrainte de conception : le notebook écrit $K\cdot P$ et $R/K$, la renumérotation fait les mêmes opérations *dans le même ordre*, et la somme des charges d'une zone parcourt les lignes de `Y_DER` dans l'ordre --- ajouter les $0{,}0$ des bus sans charge est exact en [ieee]{.smallcaps} 754.

## L'atelier

L'interface web (`g2egrid_api/`) tient en une page, un fichier JavaScript et une feuille de style : aucune chaîne de construction, aucun cadriciel, les schémas et les graphes en [svg]{.smallcaps} construit à la main.

Sa règle d'architecture est que **ce qu'on dessine [est]{.smallcaps} les tables**. L'interface ne manipule aucune structure intermédiaire : elle édite directement `Y_DER`, `T_Lines_pu`, `T_SMs`..., c'est-à-dire exactement ce que déclarent les notebooks de `reseaux_test/` et les fichiers `.m` de `G2Elin V2`. Le corollaire est que tout ce que fait l'interface se fait en code, par les mêmes fonctions : `g2egrid_api/analyse.py` ne sait pas ce qu'est [http]{.smallcaps}, et `app.py` n'est qu'une couche de transport.

:::{table} Les onglets de l'atelier.
:name: tab-atelier-onglets

| Onglet | Ce qu'on y fait |
|:---|:---|
| Réseau | poser bus, machines, convertisseurs, nœuds continus, zones ; relier, déplacer, éditer. La validation tourne en continu sous le schéma. |
| Tables | la description de référence, table par table, en miroir du schéma. |
| Point de fonctionnement | type de chaque bus (slack, PV, PQ) et consignes ; « un slack, tout le reste en PQ » ; load flow et transits sur le schéma. |
| Modes | le tableau `T_eig` et la carte des valeurs propres, axes au choix (Re/Im, $f$/$\zeta$) ; cliquer un mode ouvre sa *forme modale* et ses *sensibilités*. |
| Variations | faire varier un réglage et suivre les modes : lieu des racines, marge de stabilité, valeur critique. |
| Temporel | échelon ou scénario daté, résolu par le modèle linéaire, non linéaire, ou les deux ; mesures physiques. |
| Validation | les trois épreuves du §{ref}`sec-atelier-validation`. |
| Transport / distribution | ce que la dynamique de la distribution change au transport, face à un équivalent statique au même point de fonctionnement. |
| Comparer | deux versions du réseau : spectres appariés en entier (algorithme hongrois), modes déplacés, disparus, apparus ; cellules modifiées ; déplacement du load flow. |
| Code | le script Python qui reprend le réseau (un `Atelier`), exports notebook et MATLAB. |
:::

### Lire une forme modale

La table des modes dit *qui* oscille ; elle ne dit pas comment les états se placent les uns par rapport aux autres. C'est pourtant là que se lit la nature d'un mode. Le diagramme polaire place chaque état participant par son module relatif au plus grand et sa phase relative à l'état le plus participant. Deux vitesses rotor en phase : tout le groupe bouge ensemble, mode local. À 180 degres : les machines s'opposent, mode inter-machines. Sur `wscc9`, le mode à 2,2 Hz montre `SM3.Dwr` à 0 degres et `SM2.Dwr` à 172 degres.

### Trois nombres qui avaient l'air justes

Les défauts les plus coûteux de ce travail n'ont pas été des plantages, mais des nombres faux qui se lisaient comme des nombres justes. Ils méritent d'être consignés, parce que chacun se reproduirait volontiers ailleurs.

(le-mode-de-jauge.)=
##### Le mode de jauge.
Un réseau dont la référence d'angle n'est tenue par aucun bus infini porte un zéro *exact*, que l'algèbre rend à $10^{-12}$ près. Son amortissement vaut $-\mathrm{Re}/|\lambda|$, soit du bruit d'arrondi divisé par du bruit d'arrondi : il arrivait en tête du classement par amortissement et faisait afficher $\max\mathrm{Re} = +0{,}00\cdot10^{0}$. Les liaisons continues ont rendu le défaut criant, car elles séparent les îles et chacune apporte sa jauge.

:::{table} Marge de stabilité avant et après exclusion des modes de jauge.
:name: tab-atelier-jauge

| Réseau   | États | Modes de jauge | $\max\mathrm{Re}$ physique |
|:---------|------:|---------------:|---------------------------:|
| `smib`   |    28 |              1 |                 $-0{,}187$ |
| `wscc9`  |    93 |              1 |                 $-0{,}180$ |
| `hvdc2t` |    50 |              3 |                 $-14{,}20$ |
| `mtdc3t` |    77 |              5 |                 $-0{,}880$ |
:::

Leur nombre est désormais *affiché* plutôt que tu : il dit combien d'îles flottent sans référence absolue.

(la-consigne-perturbée-par-défaut.)=
##### La consigne perturbée par défaut.
L'entrée `w_ref` vient en tête de la liste des entrées parce que c'est l'entrée *globale*, la pulsation du repère. Un échelon de 0,01 pu dessus déplace le repère de tout le réseau d'un demi-hertz : sur `smib`, les deux voies s'écartent alors de 386 %, contre 0,15 % sur une vraie consigne d'équipement. Le réglage par défaut faisait donc conclure que le modèle linéaire était faux.

(une-cellule-dassemblage-qui-avait-divergé.)=
##### Une cellule d'assemblage qui avait divergé.
Le générateur de notebooks portait sa propre copie de l'assemblage. Elle ignorait `T_Loads`, `PARAMS_PAR_BUS` et les quatre tables continues, et surtout accumulait les $B_{pu}$ de `T_Lines_pu` *sans le changement de base ligne $\rightarrow$ système*. Sur `wscc9` l'écart était nul --- la base de ligne y est la base système et les longueurs valent 1 --- mais sur `kundur2area` le facteur atteignait 110 sur la susceptance d'une ligne. Un notebook bâti dessus décrivait un autre réseau que celui de ses propres tables, sans que rien ne le signale.

:::{admonition} Le remède est toujours le même
:class: note

Ces trois défauts, et deux autres du même ordre --- le contrat de colonnes dupliqué entre le paquet et l'interface, la fonction d'export que le paquet allait chercher *dans* l'interface --- se soignent identiquement : une seule source, et tout ce qui la lit suit. C'est ce qui a fait apparaître les quatre tables continues dans l'export, dans l'onglet Tables et dans le constructeur sans une ligne de code supplémentaire.
:::
