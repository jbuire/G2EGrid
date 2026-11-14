<!-- Page GENEREE depuis documentation/chapitres/02_principes.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ch-principes)=
# Principes de la modélisation
Ce chapitre décrit la méthode commune à V2 et à Python, indépendamment des équipements. La figure ci-dessous en donne la vue d'ensemble : les tables de données décrivent le réseau, le load flow fixe le point de fonctionnement, chaque équipement calcule son équilibre et sa linéarisation, et l'assemblage referme les connexions pour produire le modèle global.

```{mermaid}
flowchart LR
    tab["Tables<br/>Y_DER, T_*"] --> lf["Load flow<br/>V, S"]
    lf --> eq["Équilibre par élément<br/>x0, z0, u0"]
    eq --> lin["Jacobiennes<br/>par élément"]
    lin --> lft["Câblage LFT et<br/>élimination des algébriques"]
    lft --> ss["A, B, C, D"]
    ss --> modal["Analyse modale"]
```

*Chaîne de calcul commune à V2 et à Python.*

(sec-element)=
## L'élément, système algébro-différentiel
Chaque équipement est un système d'équations algébro-différentielles 

$$
\dot{\bm{x}} = \bm{f}(\bm{x}, \bm{z}, \bm{u}), \qquad
\bm{0} = \bm{g}(\bm{x}, \bm{z}, \bm{u}), \qquad
\bm{y} = \bm{h}(\bm{x}, \bm{z}, \bm{u}),
$$ (eq-dae)

 où $\bm{x}$ regroupe les états, $\bm{z}$ les variables algébriques (courants d'une machine tirés de ses flux, rapports cycliques d'un convertisseur), $\bm{u}$ les entrées et $\bm{y}$ les sorties. Les entrées et les sorties se répartissent en deux familles, qui jouent un rôle distinct dans l'assemblage.

**Entrées et sorties de service.** Ce sont les consignes et les mesures propres à l'équipement : puissance de référence $p^\star$, tension de référence $v^\star$ d'un régulateur, puissance active mesurée. Elles restent des entrées et des sorties du modèle global, ce qui permet d'y appliquer un échelon ou d'y lire une réponse.

**Entrées et sorties de réseau.** Un équipement raccordé à un nœud reçoit la tension de ce nœud $(v_d, v_q)$ et y renvoie le courant qu'il injecte $(i_d, i_q)$. Une branche a deux ports, un par extrémité. Tous les équipements reçoivent aussi deux entrées globales : la fréquence du système $\omega_{\mathrm{ref}}$ et l'angle du repère de référence $\theta_g$ (section {ref}`sec-reperes`). Ces signaux sont internes au réseau : l'assemblage les relie entre eux et ils disparaissent du modèle final.

(sec-reperes)=
## Repères de référence
### Repère global et repères locaux

Toutes les tensions de nœud et tous les courants de branche sont exprimés dans un unique repère $dq$ *global*, qui suit le rotor de l'élément de référence (section suivante) : son angle est $\theta_g$ et sa pulsation $\omega_g$, celle de cet élément. Pour un bus infini, $\omega_g$ est la pulsation nominale ; pour une machine synchrone slack, c'est la vitesse de son rotor. Un phaseur de tension de nœud s'y écrit $v_d + jv_q = Ve^{j\theta}$, et les éléments sans angle propre -- nœuds, lignes, transformateurs, charges -- écrivent leurs termes de couplage $dq$ ($\omega L i$, $\omega C v$) avec $\omega_g$. Chaque équipement qui possède un angle propre travaille en revanche dans son repère *local* : le repère rotorique d'une machine, le repère de la PLL d'un GFL, le repère de la tension interne d'un GFM. On passe de l'un à l'autre par une rotation de l'angle $\theta_g - \theta$ : 

$$
\begin{bmatrix} v_d^{\mathrm{loc}} \\ v_q^{\mathrm{loc}} \end{bmatrix} =
\begin{bmatrix} \cos(\theta_g-\theta) & -\sin(\theta_g-\theta) \\
                \sin(\theta_g-\theta) &  \cos(\theta_g-\theta) \end{bmatrix}
\begin{bmatrix} v_d \\ v_q \end{bmatrix},
\qquad
\begin{bmatrix} i_d \\ i_q \end{bmatrix} =
\begin{bmatrix} \cos(\theta-\theta_g) & -\sin(\theta-\theta_g) \\
                \sin(\theta-\theta_g) &  \cos(\theta-\theta_g) \end{bmatrix}
\begin{bmatrix} i_d^{\mathrm{loc}} \\ i_q^{\mathrm{loc}} \end{bmatrix}.
$$ (eq-rotation)

 La tension du nœud entre dans l'équipement par la première relation ; le courant qu'il injecte en ressort par la seconde, qui en est l'inverse exact.

### L'angle de référence $\theta_g$

$\theta_g$ est l'angle *interne* de l'élément de référence, c'est-à-dire du premier équipement raccordé au bus slack qui expose une sortie d'angle : l'angle de la force électromotrice pour un bus infini, l'angle rotor pour une machine. Ce n'est pas l'angle du bus slack lui-même. L'entrée $\theta_g$ de tous les équipements est câblée sur la sortie d'angle de cet élément, y compris l'élément de référence lui-même : pour lui, $\theta_g - \theta$ vaut identiquement zéro et sa force électromotrice reste fixe dans le repère qu'il définit. Sa valeur d'équilibre est celle que prend cette sortie au point de fonctionnement, ce qui reproduit la règle de V1 ($\theta_{g,0} = \alpha_s$ pour un bus infini, $\delta_{s,0} - \pi/2$ pour une machine).

:::{admonition} Un repère cohérent ne dépend pas du choix de la référence
:class: note

Les deux règles vont ensemble : les équipements à angle propre se repèrent par $\theta - \theta_g$, et les éléments passifs utilisent $\omega_g$. Si l'on applique l'une sans l'autre -- rotations relatives à $\theta_g$ mais couplages à la pulsation nominale --, le repère n'est cohérent qu'avec un bus infini, dont l'angle ne bouge pas. Avec une machine pour référence, le spectre dépend alors de la machine choisie comme slack, et un GFM se repère sur un rotor qui oscille alors que le réseau ne le suit pas (20 % d'écart sur le mode électromécanique du WSCC à GFM, section {ref}`sec-v2-multimachine`). Le test `test_spectre_independant_de_la_machine_de_reference` vérifie qu'en changeant de machine slack au même point de fonctionnement, le spectre ne bouge pas à $10^{-8}$ près.
:::

(sec-jauge)=
### Mode de jauge
Faire tourner en bloc tous les angles et toutes les grandeurs $dq$ d'un réseau laisse ses équations inchangées. Un réseau qui n'est pas ancré à un bus infini possède donc une valeur propre exactement nulle, associée à cette rotation d'ensemble. Avec un repère cohérent, elle sort numériquement vers $10^{-12}\,\mathrm{s}^{-1}$ ; une valeur plus grande, de l'ordre de $10^{-4}$, signale au contraire un repère incohérent (remarque ci-dessus). Ce n'est pas une instabilité. La méthode `is_stable` de Python en tient compte en comparant la plus grande partie réelle à $\max(10^{-6}, 10^{-9}\max|A_{ij}|)$.

(sec-noeuds)=
## Les nœuds
Un nœud reçoit les courants injectés par tous les équipements qui s'y raccordent et leur renvoie sa tension. Deux cas se présentent selon qu'une capacité shunt est disponible au nœud.

**Nœud capacitif.** Quand les lignes qui aboutissent au nœud lui apportent une charge capacitive, ou qu'un banc de condensateurs y est déclaré, la tension devient un état : 

$$
\frac{\mathrm{d}v_d}{\mathrm{d}t} = \frac{\omega_b}{C}\left(i_{\Sigma,d} + \omega\, C\, v_q\right), \qquad
\frac{\mathrm{d}v_q}{\mathrm{d}t} = \frac{\omega_b}{C}\left(i_{\Sigma,q} - \omega\, C\, v_d\right),
$$ (eq-noeud)

 où $i_\Sigma$ est la somme des courants injectés : c'est la loi des nœuds. La capacité d'un nœud vaut la moitié de la susceptance totale de chaque ligne qui y aboutit (schéma en $\pi$), plus d'éventuels shunts.

**Nœud algébrique.** Sans capacité, la tension ne peut pas être un état. Par défaut, le nœud est **idéal** : sa tension est une variable algébrique et son équation est la loi des nœuds elle-même, $i_\Sigma = 0$, exacte et sans élément fictif. Quand les courants qui s'y rejoignent sont tous des états (deux branches en série), la contrainte est d'indice 2 : l'analyse linéaire la réduit (section {ref}`sec-noeuds`), et en temporel seul BDF2 l'intègre. MATLAB pose au contraire une résistance de causalité $R = 100$ pu, et la tension suit $v = R\,i_\Sigma$. Ce shunt n'a pas de sens physique ; versé au $Y$-bus du load flow, il consomme $|V|^2/R$, soit $0{,}01\,S_b$ par nœud, jusqu'à 77 % de la charge d'un départ HTA. Python le garde en option (`NOEUDS = "resistance"`) pour les réseaux comparés à MATLAB.

(sec-loadflow)=
## Load flow
Le point de fonctionnement est calculé par un load flow de Newton-Raphson en coordonnées polaires. Chaque bus est déclaré slack (type 1, module et angle imposés), PV (type 2, puissance active et module imposés) ou PQ (type 3, puissances imposées). L'injection nette d'un bus vaut $P_g - P_l$ et $Q_g - Q_l$. La matrice d'admittance contient les branches, les capacités de nœud et les résistances de causalité des nœuds sans capacité, que le modèle dynamique du nœud contient aussi. Le shunt $R_g$ d'une machine synchrone n'y figure pas : il est interne au modèle, dont le courant injecté vaut déjà $i - v/R_g$. L'y ajouter le compterait deux fois, et la loi des nœuds du modèle dynamique serait violée de $V/R_g$ au bus de chaque machine.

:::{admonition} Deux conventions de capacité de ligne dans V1
:class: warning

Le load flow de V1 place $B/2$ à chaque extrémité d'une ligne, alors que son modèle de nœud applique la capacité *entière* à chacun des deux nœuds. Pour reproduire V1, il faut donc résoudre le load flow avec la moitié de la capacité que l'on donne ensuite au modèle dynamique. Les deux nouvelles versions appliquent la même capacité au load flow et au modèle ; seuls les scripts de comparaison rejouent la convention de V1 (chapitre {ref}`ch-validation`).
:::

À la fin du load flow, la puissance injectée par chaque bus est répartie entre les équipements qui y sont raccordés. Les charges ont une injection programmée $-(P + jQ)$ ; la part restante revient au seul équipement de production du bus. Pour que chaque puissance du load flow ait un équipement qui l'absorbe, l'assemblage des notebooks crée une charge R-L pour toute consommation déclarée dans `Y_DER`, et ajoute au load flow la consommation des charges spéciales de `T_Loads` (section {ref}`sec-tables-charges`). Un équipement déclaré sans sa puissance dans le load flow, ou l'inverse, laisserait le modèle dynamique hors équilibre au point calculé.

(sec-equilibre)=
## Équilibre des équipements
Chaque équipement reçoit du load flow la tension de son nœud et la puissance qu'il y injecte, et en déduit son état d'équilibre : angle interne, courants, flux, intégrales des correcteurs, consignes. Python calcule d'abord cet équilibre par des formules analytiques, puis le raffine par une méthode de Newton sur $[\bm{f}; \bm{g}] = \bm{0}$ à entrées fixées. V2 résout le même système par `fsolve`, sous la contrainte que la puissance injectée *au nœud* vaille celle du load flow.

:::{admonition} Conventions de V1 et raffinement
:class: note

L'équilibre calculé par V1 ne vérifie pas toujours ses propres équations (chapitre {ref}`ch-validation`). Pour reproduire exactement sa matrice $\bm{A}$, les équipements concernés disposent d'un mode d'initialisation qui rejoue ses formules et désactive le raffinement : `init='symSG'` et `init='symGFL'` en Python, `PointInit` et `InitV1` en MATLAB. Le fonctionnement par défaut reste l'équilibre exact.
:::

(sec-linearisation)=
## Linéarisation d'un équipement
Autour de l'équilibre, le système {eq}`eq-dae` se linéarise en 

$$
\Delta\dot{\bm{x}} = F_x\Delta\bm{x} + F_z\Delta\bm{z} + F_u\Delta\bm{u}, \quad
\bm{0} = G_x\Delta\bm{x} + G_z\Delta\bm{z} + G_u\Delta\bm{u}, \quad
\Delta\bm{y} = H_x\Delta\bm{x} + H_z\Delta\bm{z} + H_u\Delta\bm{u},
$$

 où $F_x = \partial\bm{f}/\partial\bm{x}$ et ainsi de suite, évalués à l'équilibre. Les jacobiennes sont obtenues symboliquement, une fois pour toutes par modèle, puis évaluées numériquement. Les variables algébriques s'éliminent par complément de Schur : 

$$
\bm{A} = F_x - F_z G_z^{-1} G_x, \quad
\bm{B} = F_u - F_z G_z^{-1} G_u, \quad
\bm{C} = H_x - H_z G_z^{-1} G_x, \quad
\bm{D} = H_u - H_z G_z^{-1} G_u.
$$ (eq-schur)

 Lorsque les variables algébriques sont données par des expressions explicites $\bm{z} = \bm{e}(\bm{x}, \bm{z}, \bm{u})$, on a $G_z = I - E_z$ et $G_x = -E_x$, d'où $\bm{A} = F_x + F_z (I - E_z)^{-1} E_x$ : c'est la forme utilisée par les classes MATLAB qui éliminent leurs algébriques à la main.

:::{admonition} Pourquoi garder des variables algébriques
:class: note

Remplacer une variable algébrique par son expression donne le même modèle *si* l'équilibre est cohérent. Quand il ne l'est pas, les jacobiennes diffèrent : $F_z$ est évalué avec la valeur d'équilibre de $\bm{z}$, et non avec celle que rendrait son expression. V1 déclare algébriques les courants de la machine, sa tension de borne, et les rapports cycliques et la pulsation de PLL des GFL. Les nouvelles versions font de même, ce qui est indispensable pour retrouver sa matrice $\bm{A}$ (section {ref}`sec-v1-algebriques`).
:::

(sec-assemblage)=
## Assemblage du modèle global
Les modèles linéarisés des éléments sont empilés en blocs diagonaux : $F_x = \mathrm{diag}(F_x^{(1)}, F_x^{(2)}, \dots)$, et de même pour les autres jacobiennes. Il reste à relier les entrées de chaque élément aux sorties des autres, ce qu'exprime une transformation fractionnaire linéaire (LFT).

Toutes les entrées des éléments forment le vecteur $\bm{U}$, toutes les sorties le vecteur $\bm{Y}$. Le routage s'écrit 

$$
\bm{U} = R_{\mathrm{ext}}\,\bm{u}_{\mathrm{ext}} + R_y\,\bm{Y},
$$

 où $R_y$ encode trois familles de connexions : les tensions de port de chaque équipement lues sur les états (ou la sortie algébrique) de son nœud, la somme des courants de port de chaque nœud, et l'angle $\theta_g$ lu sur l'élément de référence. $R_{\mathrm{ext}}$ sélectionne les entrées de service et $\omega_{\mathrm{ref}}$. En reportant $\bm{Y} = H_x\bm{x} + H_z\bm{z} + H_u\bm{U}$, on obtient 

$$
\bm{U} = (I - R_y H_u)^{-1}\left(R_{\mathrm{ext}}\,\bm{u}_{\mathrm{ext}} +
R_y H_x\,\bm{x} + R_y H_z\,\bm{z}\right),
$$

 qui existe tant qu'aucune boucle algébrique ne relie une tension à un courant. C'est le rôle des capacités de nœud et des résistances de causalité. Il suffit alors de substituer $\bm{U}$ dans $F$, $G$ et $H$, puis d'éliminer les variables algébriques globales par {eq}`eq-schur`.

(sec-modale)=
## Analyse modale
Pour chaque valeur propre $\lambda_k$ de $\bm{A}$, de vecteurs propres à droite $\bm{v}_k$ et à gauche $\bm{w}_k$ ($\bm{W} = \bm{V}^{-1}$), les notebooks et `displayModalAnalysis` calculent : 

$$
f_{\mathrm{amort},k} = \frac{|\mathrm{Im}\,\lambda_k|}{2\pi}, \qquad
f_{\mathrm{propre},k} = \frac{|\lambda_k|}{2\pi}, \qquad
\zeta_k = -\frac{\mathrm{Re}\,\lambda_k}{|\lambda_k|}, \qquad
p_{ki} = \frac{|w_{ki}|\,|v_{ik}|}{\sum_j |w_{kj}|\,|v_{jk}|}.
$$

 Le facteur de participation $p_{ki}$ mesure la part de l'état $i$ dans le mode $k$ ; les trois états les plus participants suffisent en général à attribuer un mode à un équipement et à une boucle. Les modes sont triés par partie réelle décroissante, les plus lents d'abord.
