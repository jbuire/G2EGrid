<!-- Page GENEREE depuis documentation/chapitres/06_usage.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ch-usage)=
# Mode d'emploi
Ce chapitre explique comment décrire un réseau, lancer le calcul et lire les résultats. Les tables sont les mêmes en MATLAB et en Python : en MATLAB ce sont des matrices, une ligne par élément ; dans les notebooks ce sont des listes de tuples, alignées colonne par colonne.

## Les tables d'entrée

(y-der-les-bus)=
### `Y_DER` : les bus
Une ligne par bus : `[BusID, BusType, V_base, V_pu, Theta_rad, Pg_W, Qg_var, Pl_W, Ql_var]`.

:::{table} Colonnes de `Y_DER`.
:name: tab-yder

| Colonne | Signification |
|:---|:---|
| `BusID` | numéro du bus, entier, repris par toutes les autres tables |
| `BusType` | slack, 2 PV, 3 PQ |
| `V_base` | tension nominale entre phases, en volts |
| `V_pu`, `Theta_rad` | module et angle de départ ; imposés pour un slack, module imposé pour un PV |
| `Pg_W`, `Qg_var` | production du bus, en W et var |
| `Pl_W`, `Ql_var` | charge classique du bus, en W et var : devient une charge R-L (`RLLoad`) |
:::

Le bus slack porte l'équipement de référence d'angle, en général le bus infini ou une machine synchrone. Un bus où se raccorde un convertisseur est habituellement de type PQ ; un GFM peut aussi être PV. La consommation déclarée dans `Y_DER` sert à la fois au load flow et au modèle dynamique : l'assemblage crée une charge R-L à chaque bus où `Pl_W` ou `Ql_var` est non nul. Une charge R-L exige $Q_l > 0$ ; l'assemblage s'arrête avec un message explicite sinon.

### Lignes, shunts et transformateurs

:::{table} Tables de branches et de shunts.
:name: tab-branches

| Table | Colonnes et unités |
|:---|:---|
| `T_Lines_SI` | `[ID, From, To, Length_km, R_ohm/km, X_ohm/km, B_S/km]` : valeurs linéiques, multipliées par la longueur |
| `T_Lines_pu` | `[ID, From, To, Length_km, R_pu, X_pu, B_pu, Un_V, Sb_VA]` : valeurs *totales* en pu sur $(U_n, S_{b,\mathrm{ligne}})$, longueur informative |
| `T_Shunts` | `[ID, BusID, B_pu]` : susceptance shunt en pu sur $S_b$, positive pour un condensateur ; Python seulement |
| `T_Transfos` | `[ID, From, To, Sn_VA, U1_V, U2_V, R_pu, X_pu]` : impédance en pu sur $(S_n, U_1)$, ramenée à la base système à l'assemblage. Une **prise hors nominale** se déclare par des tensions $U_1$, $U_2$ différentes de celles des bus : le rapport vaut $n = (U_1/U_2)/(U_{n,i}/U_{n,j})$, et vaut 1 quand elles coïncident |
:::

### Sources, équipements et charges

:::{table} Tables d'équipements.
:name: tab-equipements

| Table | Colonnes et unités |
|:---|:---|
| `T_IBs` | `[ID, BusID, Sn_VA, E_pu, theta_deg, R_pu, X_pu]` : impédance en pu sur $S_n$ ; `E` et `theta` sont recalculés par le load flow |
| `T_SMs` | `[ID, BusID, Sn_VA]` ou `[ID, BusID, Sn_VA, Rt_pu, Xt_pu]` |
| `T_GFM1s` à `T_GFM4s` | `[ID, BusID, Sn_VA, Vdc_V]` ou `[..., Rt_pu, Xt_pu]` |
| `T_GFL1s`, `T_GFL2s` | `[ID, BusID, Sn_VA, Vdc_V]` ou `[..., Rt_pu, Xt_pu]` |
| `T_Loads` | `[ID, BusID, P_W, Q_var, Type]` : charge spéciale, `'I'` courant constant (`ILoad`) ou `'P'` puissance constante (`PLoad`) |
:::

Les colonnes optionnelles `Rt_pu`, `Xt_pu` intègrent le transformateur de raccordement à l'équipement (section {ref}`sec-transfo-der`). La colonne `Vdc_V` fixe la base continue des convertisseurs en MATLAB ; les modèles Python travaillent directement en pu et ne l'utilisent pas.

(sec-tables-charges)=
### Charges classiques et charges spéciales
Les charges se déclarent à deux endroits, selon leur loi de charge :

- **charge classique** : les colonnes `Pl_W`, `Ql_var` de `Y_DER`. Elle devient une charge R-L à impédance constante, nommée `Load<bus>`. C'est le cas courant, et il ne demande aucune autre table ;

- **charge spéciale** : une ligne de `T_Loads`, de type `'I'` (courant constant) ou `'P'` (puissance constante), nommée `LoadP<id>_<bus>` ou `LoadI<id>_<bus>`. Sa consommation *s'ajoute* à celle de `Y_DER` dans le load flow : on ne la recopie pas dans `Y_DER`.

Un même bus peut porter les deux. Dans l'exemple suivant, le bus 7 consomme au total 1,5 MW et 0,3 Mvar : 1 MW à impédance constante et 0,5 MW à puissance constante.

```python
Y_DER = [ # [BusID, BusType, V_base, V_pu, Theta_rad, Pg_W, Qg_var,   Pl_W, Ql_var]
        (     7,       3,   20e3,  1.00,        0,    0,      0,    1.0e6,  0.2e6),
]
T_Loads = [ # [ID, BusID, P_W, Q_var, Type]
    (  1,  7, 0.5e6, 0.1e6, 'P'),
]
```

Les notebooks de test n'ont que des charges classiques : `T_Loads` y est vide, prêt à recevoir des charges spéciales. En MATLAB, `buildFromMatrices` garde l'ancienne convention : la consommation de `Y_DER` ne sert qu'au load flow, et `T_Loads = [ID, BusID, P_W, Q_var]` crée les `GenericLoad` qui lui correspondent.

### Réglages des équipements

Dans les notebooks, les réglages communs à tous les équipements d'une famille se passent par un dictionnaire placé juste après sa table : `SM_PARAMS`, `GFM_PARAMS`, `GFL_PARAMS`. Un dictionnaire vide garde les valeurs par défaut du chapitre {ref}`ch-modeles`. En MATLAB, on modifie les objets après `buildFromMatrices`, en parcourant `net.Elements`.

```python
T_SMs = [
    (  1,  1, 247.5e6),
]
SM_PARAMS = dict(H=6.5, K_PSS=0.0)      # toutes les machines du reseau

T_GFM3s = [
    (  1,  2, 1e6, 29394, 0.0002, 0.024),   # transfo integre
]
GFM_PARAMS = dict(mp=0.005, tau_VL=15e-3)
```

## Déroulé d'un notebook

Tous les notebooks de `g2egrid/reseaux_test/` suivent la même structure, celle de `cigre_mv.ipynb`. Pour un nouveau réseau, on copie un notebook proche et on ne réécrit que la présentation et la cellule de données.

1.  **Imports.** Le paquet est chargé depuis le dossier parent.

2.  **Données.** `NOM_RESEAU`, la bande de zoom `F_ZOOM` de la carte des modes, `f0`, `Sb`, puis toutes les tables, vides ou non.

3.  **Assemblage.** Calcul des capacités de bus, puis une boucle `net.add(...)` par table. Cette cellule ne se modifie pas.

4.  **Load flow.** Tableau des tensions et des injections de chaque bus.

5.  **Modèle d'état.** Équilibre, assemblage et liste des états par élément.

6.  **Tableau modal.** DataFrame pandas `T_eig`, une ligne par valeur propre triée par partie réelle décroissante : `Re`, `Im`, `f_und [Hz]`, `f_dam [Hz]`, `zeta [%]` et les trois états les plus participants avec leur participation en %. L'index `Mode` est l'indice $k$ de `lam[k]`. Les modes instables sont surlignés en rouge, les modes oscillants amortis à moins de 5 % en orange. Un second tableau, `T_lents`, isole les modes oscillants sous `F_ZOOM` par amortissement croissant.

7.  **Carte des valeurs propres.** Tous les modes, puis un zoom sur les modes lents étiquetés par l'élément dominant.

## Déroulé d'un script MATLAB

Les réseaux de `reseaux_test/reseaux_classiques/` sont des fonctions qui rendent un `Network` prêt, load flow compris. On les utilise ainsi :

```matlab
addpath(genpath('matlab'));
net = WSCC();
net.displayLoadFlow();
[Sys, EigenValues, LinearModel] = net.buildGlobalModel();
T = net.displayModalAnalysis(Sys.A, EigenValues, Sys.StateName);
```

## Lire le tableau modal

Chaque ligne du tableau est un mode. Une valeur propre réelle négative est un mode apériodique, de constante de temps $-1/\mathrm{Re}\,\lambda$. Une paire complexe apparaît sur deux lignes consécutives, de même partie réelle. Les modes électromécaniques d'une machine se reconnaissent à la participation de `Dwr` et `theta`, entre 0,1 et 3 Hz ; ceux du régulateur de tension à `phi_fd`, `efd` et `e2` ; ceux d'une PLL à `M_pll` et `theta_pll` ; les modes électriques du réseau, très rapides, aux courants de branche et aux tensions de nœud. Un amortissement inférieur à 5 % signale en général un réglage à revoir.

La première ligne d'un réseau sans bus infini est le mode de jauge, à $\lambda \approx 0$ : il est attendu et ne traduit aucune instabilité (section {ref}`sec-jauge`).

## Pièges fréquents

:::{admonition} Avant de conclure à une instabilité
:class: warning

**Charge à $Q$ négatif ou nul dans `Y_DER`.** Une charge R-L ne peut pas consommer de réactif négatif sans devenir instable, et un réactif nul lui donnerait une inductance nulle. Déclarer les condensateurs dans `T_Shunts`, ou la charge dans `T_Loads` en type `'P'` ou `'I'`, qui accepte tout $Q$.

**Charge comptée deux fois.** Une charge de `T_Loads` s'ajoute à `Y_DER` : la recopier aussi dans `Pl_W`, `Ql_var` doublerait la consommation.

**Valeurs en pu d'une ligne.** `T_Lines_pu` attend des valeurs totales, et `T_Lines_SI` des valeurs linéiques. Déclarer des valeurs totales avec une longueur de 25 km dans `T_Lines_SI` rend la ligne 25 fois trop longue.

**Base d'un transformateur.** Son impédance est donnée sur sa propre puissance $S_n$. Un transformateur de 900 MVA à 15 % vaut 1,67 % sur une base de 100 MVA.

**Transformateur compté deux fois.** Un convertisseur déclaré à six colonnes porte déjà son transformateur ; ne pas en ajouter un autre dans `T_Transfos`.

**Pas d'élément de référence.** Le bus slack doit porter un équipement doté d'un angle : bus infini, machine synchrone ou GFM.
:::
