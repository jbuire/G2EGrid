<!-- Page GENEREE depuis documentation/chapitres/01_introduction.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ch-intro)=
# Introduction
## Objet de l'outil

G2Egrid construit le modèle linéarisé, dit *petit-signal*, d'un réseau électrique quelconque autour d'un point de fonctionnement. Le réseau peut mêler des machines synchrones, des convertisseurs *grid-forming* (GFM) et *grid-following* (GFL), des bus infinis, des lignes, des transformateurs et des charges. Le résultat est une représentation d'état 

$$
\dot{\bm{x}} = \bm{A}\,\bm{x} + \bm{B}\,\bm{u}, \qquad
  \bm{y} = \bm{C}\,\bm{x} + \bm{D}\,\bm{u},
$$

 dont on tire l'analyse modale : valeurs propres, fréquences, amortissements et facteurs de participation. Cette analyse permet d'identifier les modes mal amortis, d'attribuer chacun d'eux aux équipements qui y participent et de régler les correcteurs en conséquence.

Tous les modèles sont écrits dans des repères $dq$ tournants et conservent la dynamique des inductances et des capacités du réseau. Ils vont donc au-delà du modèle électromécanique classique : les modes électriques rapides des lignes, des filtres LCL et des boucles de courant apparaissent dans le spectre au même titre que les modes électromécaniques.

## Trois versions de l'outil

L'outil existe en trois générations, qui doivent produire la même matrice $\bm{A}$ pour un même réseau.

:::{table} Les trois versions de G2Egrid.
:name: tab-versions

| Version | Nature | Rôle |
|:---|:---|:---|
| V1 | MATLAB procédural : fonctions `symSG`, `symGFM`, `symGFL`, `symIB`, `symLine`, `symNode`, `symLoad`, assemblage par `assoc_matrices`. Dossiers `Functions/` et `Symbolic/`. | la référence |
| V2 | MATLAB orienté objet : une classe par équipement, une classe `Network` qui construit le réseau depuis des tables, calcule le load flow et assemble le modèle. Génération automatique du modèle Simulink. | doit suivre V1 |
| Python | paquet `g2egrid` : éléments déclaratifs dont les équations sont écrites une seule fois en calcul symbolique, compilées en fonctions numériques, puis assemblées par la même méthode que V2. | doit suivre V1 |
:::

La version V1 fait foi. Ses fichiers `Functions/` et `Symbolic/` ne sont jamais modifiés ; V2 et Python sont alignés sur elle. Le chapitre {ref}`ch-validation` montre que les trois versions produisent aujourd'hui la même matrice $\bm{A}$, coefficient par coefficient, sur les réseaux élémentaires de référence.

Les deux nouvelles versions partagent la même organisation. Un réseau se décrit par des tables de données identiques d'une version à l'autre : une table de bus `Y_DER`, puis une table par famille d'équipements. Chaque équipement est un objet qui connaît ses équations, calcule son propre point d'équilibre et fournit son modèle linéarisé. La classe réseau se charge du load flow, câble les équipements entre eux et élimine les variables algébriques.

## Organisation du dépôt

:::{table} Arborescence utile du dépôt `g2elin_v2`.
:name: tab-depot

| Dossier | Contenu |
|:---|:---|
| `matlab/G2Elin_V1/` | la référence, à ne pas modifier |
| `matlab/G2Elin_V2/core/` | classes de V2 : `Network`, `Node`, `TransmissionLine`, `Transformer`, `IB`, `SM`, `SMTR`, `GFM1` à `GFM4`, `GFL1`, `GFL2`, `GenericLoad`, `AdvancedLoad`, `SimulinkGenerator` |
| `matlab/G2Elin_V2/reseaux_test/` | réseaux classiques et tests élémentaires |
| `g2egrid/g2egrid/` | paquet Python : noyau (`element.py`, `symbolic.py`, `loadflow.py`, `network.py`) et bibliothèque d'équipements (`elements/`) |
| `g2egrid/reseaux_test/` | un notebook par réseau de test, tous au même format |
| `python/etudes/` | notebooks d'analyse : nature des sources (chapitre {ref}`ch-sources`), séparation transport / distribution (chapitre {ref}`ch-td`), transition vers les onduleurs (chapitre {ref}`ch-transition`) |
| `g2egrid/tests/` | suite `pytest` du noyau |
| `python/developpement/` | outils de travail : comparaison V1/V2/Python, générateur de notebooks, exemples, anciennes versions |
| `documentation/` | ce document, ses figures et le script qui les produit |
:::

## Prise en main rapide

Le plus court chemin pour obtenir un modèle est de copier l'un des notebooks de `g2egrid/reseaux_test/` et d'en modifier les tables. Chaque notebook enchaîne sept cellules de code : imports, données, assemblage, load flow, modèle d'état, tableau modal et carte des valeurs propres. Seule la cellule de données change d'un réseau à l'autre. Côté MATLAB, les fichiers du dossier

`matlab/G2Elin_V2/reseaux_test/reseaux_classiques/`
:::

jouent le même rôle, avec les mêmes tables.

```python
net = Network(f0=f0, Sb=Sb)
# ... une boucle net.add(...) par table (cellule d'assemblage des notebooks)
lf = net.run_loadflow()          # load flow de Newton
net.solve_equilibrium()          # equilibre element par element
ss = net.build_state_space()     # modele d'etat global
lam = ss.eig()                   # valeurs propres
```

```matlab
net = Network(fb, Sb);
net.buildFromMatrices(Y_DER, T_Lines_SI, T_Lines_pu, T_Transfos, T_IBs, ...
    T_SMs, T_GFM1s, T_GFM2s, T_GFM3s, T_GFM4s, T_GFL1s, T_GFL2s, T_Loads);
net.Load_flow(Y_DER, false);
[Sys, EigenValues, LinearModel] = net.buildGlobalModel();
```

## Notations

Les grandeurs sont en per-unit (annexe {ref}`ann-pu`) sauf mention contraire. $\omega_b= 2\pi f_0$ désigne la pulsation de base, $\omega$ une pulsation en pu. Un phaseur de tension de bus s'écrit $v_d + j v_q = V e^{j\theta}$ dans le repère global. L'indice $g$ désigne le côté réseau d'un équipement, l'indice $e$ la borne interne (condensateur de filtre ou borne machine), l'indice $s$ le côté convertisseur. Les courants de port sont comptés positifs lorsqu'ils sont *injectés* dans le nœud. Les noms d'états sont ceux du code, écrits en `police machine`.
