<!-- Page GENEREE depuis documentation/chapitres/09_transport_distribution.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ch-td)=
# Séparation transport / distribution
Un modèle petit-signal de réseau de transport s'arrête d'ordinaire au poste source : la distribution y est une charge. Ce chapitre montre ce que cette simplification coûte, et comment la remplacer par un équivalent dynamique juste. Il met en œuvre la note *Impact du réseau de distribution sur le modèle d'état du réseau de transport* sur le réseau `wscc9_cigre` (section {ref}`sec-wscc9-cigre`), dont le bus 5 alimente quatre réseaux CIGRE agrégés.

## Partition du modèle

On sépare les états du modèle global en transport $x_T$ et distribution $x_D$ : 

$$
\begin{bmatrix}\dot x_T\\ \dot x_D\end{bmatrix} =
\begin{bmatrix}\bm{A}_{TT} & \bm{A}_{TD}\\ \bm{A}_{DT} & \bm{A}_{DD}\end{bmatrix}
\begin{bmatrix}x_T\\ x_D\end{bmatrix} +
\begin{bmatrix}\bm{B}_T\\ \bm{B}_D\end{bmatrix} u .
$$

 Ici 84 états de transport et 146 de distribution. $\bm{A}_{TD}$ est l'action de la distribution sur le transport, $\bm{A}_{DT}$ celle du transport sur la distribution. Trois modèles du transport en découlent :

| Modèle | Expression | Hypothèse |
|:---|:---|:---|
| transport seul | $\bm{A}_{TT}$ | la distribution n'a aucun effet dynamique |
| élimination exacte | $\bm{A}_{TT} + \bm{A}_{TD}(s\bm{I}-\bm{A}_{DD})^{-1}\bm{A}_{DT}$ | aucune |
| quasi statique | $\bm{A}_{TT} - \bm{A}_{TD}\bm{A}_{DD}^{-1}\bm{A}_{DT}$ | distribution rapide |
:::

(sec-td-interface)=
### Une interface étroite
Les blocs de couplage sont creux, et leur structure se lit directement : la distribution n'excite que les *deux* états de tension du bus 5, et ne voit du transport que quatre états. Son action sur le transport passe donc entièrement par le courant qu'elle injecte au poste, et le transport ne lui parle que par la tension de ce poste : 

$$
\dot x_D = \bm{A}_D x_D + \bm{B}_D \begin{bmatrix} v_{PCC} \\ \theta_g \\ \omega_g \end{bmatrix},
\qquad
i_{PCC} = \bm{C}_D x_D + \bm{D}_D \begin{bmatrix} v_{PCC} \\ \theta_g \\ \omega_g \end{bmatrix}.
$$

 Aux deux composantes $dq$ de la tension du poste s'ajoutent l'angle et la pulsation de la machine de référence : ils définissent le repère $dq$ commun (section {ref}`sec-reperes`), dont dépendent les rotations des PLL de la distribution et les termes de couplage de ses lignes et de ses charges. Un équivalent qui les ignorerait ne verrait pas le transport tourner. V2 et le portage Python trouvent exactement la même interface.

Le nœud du poste intègre $\dot v = (\omega_b/C)\sum i$ : le courant injecté se retrouve donc en divisant la ligne correspondante de $\bm{A}_{TD}$ par $\omega_b/C$. C'est ce qui donne l'**admittance dynamique** de la distribution vue du poste, $\bm{Y}_D(s)$, tracée plus bas.

## Ce que coûte le transport seul

:::{table} Modes lents du transport : modèle complet, $\bm{A}_{TT}$ seul et quasi statique. L'écart est la distance dans le plan complexe au mode du modèle complet.
:name: tab-td-modes

| Complet | $A_{TT}$ seul | Ecart | Quasi statique | Ecart |
|:---|:---|:---|:---|:---|
| $-0{,}0000$ | $-0{,}0000$ | $1{,}9\cdot 10^{-12}$ | $-0{,}0000$ | $2{,}9\cdot 10^{-11}$ |
| $-0{,}1846$ | $-0{,}1853$ | $6{,}5\cdot 10^{-4}$ | $-0{,}1846$ | $1{,}5\cdot 10^{-7}$ |
| $-0{,}5152$ | $-0{,}5214 + j\,0{,}0411$ | $4{,}2\cdot 10^{-2}$ | $-0{,}5152$ | $5{,}3\cdot 10^{-6}$ |
| $-0{,}6053$ | $-0{,}5214 - j\,0{,}0411$ | $9{,}3\cdot 10^{-2}$ | $-0{,}6053$ | $2{,}9\cdot 10^{-5}$ |
| $-1{,}9115 - j\,12{,}2517$ | $-1{,}9446 - j\,12{,}2971$ | $5{,}6\cdot 10^{-2}$ | $-1{,}9120 - j\,12{,}2532$ | $1{,}6\cdot 10^{-3}$ |
| $-1{,}9115 + j\,12{,}2517$ | $-1{,}9446 + j\,12{,}2971$ | $5{,}6\cdot 10^{-2}$ | $-1{,}9120 + j\,12{,}2532$ | $1{,}6\cdot 10^{-3}$ |
| $-2{,}4735$ | $-2{,}4789$ | $5{,}3\cdot 10^{-3}$ | $-2{,}4735$ | $1{,}4\cdot 10^{-6}$ |
| $-4{,}1623 - j\,7{,}0736$ | $-4{,}1596 - j\,7{,}0965$ | $2{,}3\cdot 10^{-2}$ | $-4{,}1595 - j\,7{,}0756$ | $3{,}4\cdot 10^{-3}$ |
| $-4{,}1623 + j\,7{,}0736$ | $-4{,}1596 + j\,7{,}0965$ | $2{,}3\cdot 10^{-2}$ | $-4{,}1595 + j\,7{,}0756$ | $3{,}4\cdot 10^{-3}$ |
:::

**$\bm{A}_{TT}$ seul est faux.** Figer $x_D$ ne supprime pas la distribution : cela la remplace par une injection de courant constante, qui ne réagit plus à la tension. Les modes lents se déplacent jusqu'à $9\cdot 10^{-2}\,\mathrm{s}^{-1}$, et deux modes réels ($-0{,}52$ et $-0{,}61$) fusionnent en une oscillation parasite. En régime permanent, l'erreur est plus parlante encore : après un échelon de consigne de la machine, la tension du poste se stabilise 4,5 fois trop haut (figure {numref}`fig-td-echelons`).

**Le quasi statique suffit ici.** Le mode le plus lent de la distribution seule est à $-29\,\mathrm{s}^{-1}$, quinze fois plus rapide que les modes électromécaniques : l'hypothèse de la note est vérifiée, et l'erreur tombe à $3\cdot 10^{-3}$. C'est la version dynamiquement correcte de la charge agrégée.

:::{admonition} L'élimination exacte n'est pas une matrice d'état
:class: note

$\bm{A}_{\mathrm{eff}}(s)$ dépend de $s$ : la distribution y apparaît comme une fonction de transfert, pas comme une matrice constante. Les valeurs propres du modèle complet sont exactement les solutions de $\det(\lambda\bm{I} - \bm{A}_{\mathrm{eff}}(\lambda)) = 0$ ; on le vérifie sur les neuf modes lents en calculant la plus petite valeur singulière de $\lambda\bm{I} - \bm{A}_{\mathrm{eff}}(\lambda)$, qui vaut $10^{-20}$ en relatif.
:::

## Réduire la distribution

Le modèle d'interface a 146 états. On l'équilibre (valeurs singulières de Hankel), puis on élimine les états les moins commandables et observables de deux façons : la **troncature** les supprime, la **résidualisation** les suppose à l'équilibre. La première est exacte en haute fréquence, la seconde conserve le gain statique -- et à l'ordre 0, elle redonne exactement le modèle quasi statique.

:::{table} Écart maximal sur les modes lents du transport, selon l'ordre et la méthode de réduction du modèle de distribution (146 états à l'origine).
:name: tab-td-reduction

| Ordre | Troncature                      | Residualisation      |
|------:|:--------------------------------|:---------------------|
|     0 | --                              | $3{,}4\cdot 10^{-3}$ |
|     4 | $9{,}3\cdot 10^{-2}$            | $3{,}4\cdot 10^{-3}$ |
|     8 | $1{,}8\cdot 10^{-2}$            | $2{,}9\cdot 10^{-3}$ |
|    12 | $1{,}6\cdot 10^{-2}$            | $3{,}0\cdot 10^{-3}$ |
|    16 | $7{,}7\cdot 10^{-3}$ (instable) | $2{,}9\cdot 10^{-3}$ |
|    24 | $8{,}4\cdot 10^{-3}$ (instable) | $7{,}5\cdot 10^{-4}$ |
|    32 | $9{,}5\cdot 10^{-3}$            | $3{,}0\cdot 10^{-4}$ |
|    48 | $8{,}1\cdot 10^{-4}$ (instable) | $2{,}1\cdot 10^{-5}$ |
:::

Les modes qui nous intéressent sont lents devant la distribution : c'est le gain statique qui compte, et la troncature seule est donc inadaptée. Elle fait même apparaître des instabilités parasites, alors que la résidualisation conserve la stabilité à tous les ordres et converge. À l'ordre 32, le modèle de distribution est 4,5 fois plus petit et l'erreur onze fois plus faible qu'avec le quasi statique.

## Validation

```{figure} ../_static/figures/td_modes.png
:name: fig-td-modes
:width: 100%

Modes lents du transport selon le modèle de distribution.
```

```{figure} ../_static/figures/td_admittance.png
:name: fig-td-admittance
:width: 100%

Admittance de la distribution vue du poste, complète, réduite à l'ordre 32 et quasi statique. Le quasi statique est une constante : il ne peut pas reproduire la résonance du filtre des GFL ni celle des câbles.
```

```{figure} ../_static/figures/td_echelons.png
:name: fig-td-echelons
:width: 100%

Échelons de consigne, côté transport (en haut) et côté distribution (en bas). <span class="math inline"><strong>A</strong><sub><em>T</em><em>T</em></sub></span> seul se trompe de régime permanent.
```

:::{table} Écart maximal après 0,1 s, rapporté à l'amplitude de la réponse complète. $\bm{A}_{TT}$ seul n'a pas de courant de poste : ses colonnes $P$ et $Q$ n'ont pas de sens.
:name: tab-td-essais

| Echelon | Modele | vitesse SM1 \[pu\] | \|V_PCC\| \[pu\] | P soutire \[MW\] | Q soutire \[Mvar\] |
|:---|:---|---:|---:|---:|---:|
| machine | reduit (ordre 32) | $3{,}0\cdot 10^{-5}$ | $2{,}6\cdot 10^{-4}$ | $5{,}7\cdot 10^{-4}$ | $6{,}0\cdot 10^{-3}$ |
| machine | quasi statique | $2{,}4\cdot 10^{-4}$ | $2{,}1\cdot 10^{-3}$ | $6{,}0\cdot 10^{-3}$ | $3{,}6\cdot 10^{-2}$ |
| machine | A_TT seul | $5{,}1\cdot 10^{-2}$ | $2{,}2\cdot 10^{0}$ | -- | -- |
| eolien GFL107 | reduit (ordre 32) | $6{,}5\cdot 10^{-3}$ | $7{,}3\cdot 10^{-4}$ | $6{,}0\cdot 10^{-3}$ | $1{,}9\cdot 10^{-3}$ |
| eolien GFL107 | quasi statique | $2{,}4\cdot 10^{-1}$ | $2{,}2\cdot 10^{-2}$ | $5{,}1\cdot 10^{-2}$ | $4{,}7\cdot 10^{-2}$ |
:::

L'essai le plus instructif est le second. Tant que la perturbation vient du transport, le quasi statique suffit ; mais dès qu'elle naît *dans* la distribution -- un échelon de consigne de l'éolien --, elle traverse justement la dynamique qu'il a supprimée, et l'erreur sur la réponse de vitesse de la machine atteint 24 %. Le modèle réduit d'ordre 32 reste sous $7\cdot 10^{-3}$.

Pendant les premières millisecondes, les deux modèles réduits répondent instantanément là où le complet met quelques millisecondes : c'est le prix de la résidualisation, qui remplace les états rapides par leur équilibre. D'où la mesure des écarts après 0,1 s.

## Ce qu'il faut retenir

1.  **Le transport seul est faux** : il fige les courants de la distribution, déplace les modes lents et se trompe de régime permanent.

2.  **L'interface est étroite** : tension du poste en entrée, courant du poste en sortie, plus l'angle et la pulsation de la machine de référence.

3.  **Le quasi statique suffit pour les modes et les perturbations du transport** quand la distribution est nettement plus rapide. Il ne suffit plus pour une perturbation née dans la distribution, ni avec une distribution à dynamique lente propre : régulateurs de tension lents, charges moteurs, GFM de distribution.

4.  **Pour réduire, résidualiser plutôt que tronquer.**

## Outils

Côté Python, tout est dans `g2egrid/separation.py` : `separer` partitionne et construit le modèle d'interface, `Separation.a_effectif` et `a_quasi_statique` donnent les deux modèles du transport, `troncature_equilibree` réduit, `Separation.modele_global` reconnecte et `reponse_indicielle` simule. Le notebook `python/etudes/transport_distribution.ipynb` produit les figures de ce chapitre. Côté MATLAB, `matlab/G2Elin_V2/etudes/transport_distribution.m` refait la même étude avec `separerTD`, `modeleGlobalTD`, `reductionEquilibree` et `reponseIndicielle` ; les deux outils donnent les mêmes chiffres à trois décimales.

```python
from g2egrid.multi_niveau.separation import separer, troncature_equilibree

# tout element dont un bus depasse 100 appartient a la distribution
elements_D = {el.name for el in net.all_elements if any(b >= 100 for b in el.buses)}
sep = separer(ss, lambda etat: etat.split(".")[0] in elements_D)

red, hsv = troncature_equilibree(sep.distribution, 32)   # residualisation
A_reduit, B_reduit, etats = sep.modele_global(red)
```
