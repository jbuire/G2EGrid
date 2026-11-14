# Décrire un équipement par un schéma-bloc

Un nouvel équipement n'a pas toujours besoin d'une classe Python
({doc}`ajouter_un_equipement`). S'il se décrit par des blocs de commande
(gains, sommes, premiers ordres, PI, intégrateurs, fonctions...) reliés entre
eux, comme sous Simulink, le module {mod}`g2egrid.elements.blocs` le compile en équations
et en fait un élément du réseau.

## Le schéma

Un schéma est un dictionnaire JSON :

```python
schema = {
    "nom": "source_pi",
    "parametres": {"Kp": 0.5, "Ki": 20.0, "tm": 0.02},
    "blocs": [
        {"id": "pr", "type": "consigne", "nom": "P_ref"},
        {"id": "qr", "type": "consigne", "nom": "Q_ref"},
        {"id": "e", "type": "somme", "signes": "+-"},
        {"id": "c", "type": "pi", "Kp": "Kp", "Ki": "Ki"},
        {"id": "m", "type": "premier_ordre", "K": 1, "tau": "tm"},
        {"id": "op", "type": "sortie", "signal": "P"},
        {"id": "oq", "type": "sortie", "signal": "Q"},
    ],
    "liens": [
        {"de": "pr", "vers": "e", "entree": 0}, {"de": "m", "vers": "e", "entree": 1},
        {"de": "e", "vers": "c"}, {"de": "c", "vers": "m"},
        {"de": "c", "vers": "op"}, {"de": "qr", "vers": "oq"},
    ],
}
```

| type | entrées | champs | ce qu'il calcule |
|---|---|---|---|
| `entree` | 0 | `signal` : `vd`, `vq`, `V`, `w` | tension ou fréquence du bus |
| `consigne` | 0 | `nom` | entrée de service, calée au load flow |
| `constante` | 0 | `valeur` | une constante |
| `gain` | 1 | `K` | $K u$ |
| `somme` | une par signe | `signes` (`"+-"`...) | $\pm u_1 \pm u_2 \dots$ |
| `produit`, `division` | 2 | | $u_1 u_2$, $u_1/u_2$ |
| `fonction` | `n` | `expr`, `n` | expression en `u` (ou `u1`, `u2`...) |
| `integrateur` | 1 | `K` | $K/s$ |
| `premier_ordre` | 1 | `K`, `tau` | $K/(1+\tau s)$ |
| `pi` | 1 | `Kp`, `Ki` | $K_p + K_i/s$ |
| `derivee` | 1 | `K`, `tau` | $K s/(1+\tau s)$ |
| `park` | 3 : $x_d$, $x_q$, $\theta$ | | sorties `d`, `q` : passage au repère d'angle $\theta$ |
| `park_inverse` | 3 : $d$, $q$, $\theta$ | | sorties `d`, `q` : retour au repère global |
| `rl_dq` | 5 : $v_{1d}$, $v_{1q}$, $v_{2d}$, $v_{2q}$, $\omega$ | `R`, `L` | courant de la branche R-L dans un repère tournant à $\omega$ |
| `c_dq` | 5 : $i_{1d}$, $i_{1q}$, $i_{2d}$, $i_{2q}$, $\omega$ | `C` | tension d'un condensateur chargé par $i_1 - i_2$ |
| `limiteur_courant` | 2 : $i_d^*$, $i_q^*$ | `Imax`, `Taw` | sorties `d`, `q` : la consigne ramenée dans $\lvert i\rvert \le I_{max}$ ; `aw_d`, `aw_q` : le rétro-calcul $(i_{lim} - i^*)/T_{aw}$ |
| `saturation_tension` | 2 : $m_d^*$, $m_q^*$ | `Mmax`, `Taw` | sorties `d`, `q` : la modulation ramenée dans $\lvert m\rvert \le M_{max}$ ; `aw_d`, `aw_q` : le rétro-calcul |
| `pi_aw` | 2 : $e$, $aw$ | `Kp`, `Ki` | PI dont l'intégrale reçoit le rétro-calcul : $\dot x = K_i e + aw$, $y = K_p e + x$ |
| `pid` | 1 | `Kp`, `Ki`, `Kd`, `Tf` ($T_f > 0$) | $K_p + K_i/s + K_d s/(1+T_f s)$ ; états `x_<id>_1` (intégrale), `x_<id>_2` (filtre de la dérivée) |
| `avance_retard` | 1 | `K`, `T1`, `T2` ($T_2 > 0$) | $K(1+T_1 s)/(1+T_2 s)$ ; un état `x_<id>` |
| `second_ordre` | 1 | `K`, `wn`, `zeta` | $K\omega_n^2/(s^2 + 2\zeta\omega_n s + \omega_n^2)$ ; états `x_<id>_1` ($y$), `x_<id>_2` ($\dot y$) |
| `fonction_transfert` | 1 | `num`, `den` | $\mathrm{num}(s)/\mathrm{den}(s)$, listes de coefficients (puissances décroissantes ; nombres ou expressions des paramètres), propre ; forme compagne, états `x_<id>_1`… `x_<id>_n` |
| `saturation` | 1 | `min`, `max` | $\min \le y \le \max$ ; **libre** au point de fonctionnement ($y = u$), drapeau `mode_<id>` basculé par les gardes de la simulation temporelle, comme un limiteur |
| `implicite` | `n` | `expr`, `n` | $y$ tel que `expr`$(u_1..u_n, y) = 0$ : une contrainte algébrique, pas une affectation |
| `sous_systeme` | ses `port_entree` | `nom`, `schema` | un schéma imbriqué (voir plus bas) |
| `port_entree`, `port_sortie` | 0, 1 | `index`, `nom` | l'interface d'un sous-système (admis seulement dedans) |
| `sortie` | 1 | `signal` : `P`, `Q` ou `id`, `iq` ; `theta`, `w` | ce que l'élément injecte au réseau ; son angle et sa vitesse |

Chaque type porte une `categorie` (`sources`, `operations`, `dynamique`,
`repere_dq`, `limiteurs`, `structure`, `sorties`) : la palette de l'atelier
est construite depuis `GET /api/blocs/types`, sans liste en dur.

Une `consigne` a deux champs facultatifs, `calee` (vrai par défaut) et
`valeur` (0 par défaut). Une consigne non calée est une entrée de service
fixée à `valeur` ; l'initialisation exige toujours **exactement deux**
consignes calées.

Un champ numérique est un nombre ou une expression des `parametres`. Les
quatre derniers blocs ont DEUX sorties, `d` et `q` : un lien qui part de la
seconde porte `"sortie": 1`. Les entrées `entree` lisent `vd`, `vq`, `V`, `w`
(fréquence du repère du réseau) et `theta_g` (angle de l'élément de
référence).

Deux paramètres sont **réservés** et fournis par le réseau : `wb`, la
pulsation de base ($2\pi f_0$), et `kS = Sn/Sb`, tiré de la colonne `Sn_VA`
de `T_Blocs`. Un schéma s'en sert sans les déclarer.

### Ce que la compilation refuse

`compiler(schema)` lève `SchemaInvalide` (une `ValueError`), avec la raison et
le bloc en cause, plutôt que de rendre un modèle qui échouera plus loin :

- un schéma **mal formé** : `blocs` ou `liens` qui ne sont pas des listes
  d'objets, `parametres` qui n'est pas un objet, sous-système sans schéma ;
- une **valeur** qui n'est pas un nombre fini : paramètre `NaN` ou infini,
  booléen, constante hors de portée (`1e999`) ;
- un **diviseur nul** : `tau = 0` d'un `premier_ordre` ou d'une `derivee`,
  `L = 0` d'un `rl_dq`, `C = 0` d'un `c_dq`, et tout paramètre nul au
  dénominateur d'une équation (`wb*(...)/L1_pu` avec `L1_pu = 0`) — une
  valeur négative reste admise, c'est un modèle instable mais bien posé ;
- une `saturation` dont `min` dépasse `max`, un `pid` ou un `avance_retard`
  dont la constante de temps du dénominateur n'est pas strictement positive ;
- une boucle algébrique sans état, une entrée libre, un nom inconnu dans une
  expression, autre chose que deux consignes calées.

## Un GFL complet, bloc à bloc

{func}`g2egrid.elements.blocs.schema_gfl` dessine un grid-following entier (33 blocs) :

- le **filtre LCL** dans le repère de la PLL : `rl_dq` (R1, L1), `c_dq` (C),
  `rl_dq` (R2, L2), tournant à $\omega_{pll}$ ;
- la **PLL** : `park` de la tension du bus, PI sur $v_{eq}$, plus
  $\omega_{ff}$, puis l'intégrateur $\dot\theta = \omega_b(\omega_{pll} - \omega)$ ;
- les **boucles de puissance** : PI sur $P$ et $Q$ mesurées
  ($p = v_e \cdot i_g$, `fonction`) qui donnent les consignes de courant ;
- la **limitation de courant** : la consigne $(i_d^*, i_q^*)$ des PI de
  puissance passe par `limiteur_courant`, et ces PI (`pi_aw`) suivent la
  consigne limitée par rétro-calcul ;
- les **boucles de courant** : PI, anticipation de $v_e$ et découplage
  $\omega_{ff} L_1 i$ ;
- la **saturation de tension** (`vsat`) entre la boucle de courant et le
  filtre : la modulation est ramenée dans $\lvert m\rvert \le M_{max} = 1{,}15$,
  avec rétro-calcul sur les PI de courant ;
- le **retour au réseau** : `park_inverse`, puis le gain `kS`.

Ce sont exactement les équations du `GFL1` du moteur (bus continu idéal ;
limitation et saturation sont au repos au point de fonctionnement) : les tests vérifient que les valeurs propres d'un
réseau sont les mêmes, à $10^{-7}$ près, avec l'un ou l'autre. C'est le point
de départ pour essayer une autre PLL, une autre boucle de courant ou un
filtre L seul. Comme les GFL du moteur, `park` et `park_inverse` retranchent
d'eux-mêmes l'angle $\theta_g$ de l'élément de référence, et l'initialisation
part de l'angle de la tension du bus pour placer la PLL sur la bonne branche.

## Ce que le compilateur en fait

Chaque bloc qui n'est pas une source reçoit une variable algébrique `y_<id>`,
chaque bloc dynamique un état `x_<id>`. Les boucles algébriques (une somme
rebouclée par un gain) se traitent donc comme le reste : l'élément est un
système algébro-différentiel. Des sorties `P`, `Q` donnent le courant injecté
$i = \overline{S/V_{mes}}$ sur la tension **mesurée**, filtrée au premier ordre
(`tau_mes`, 2 ms par défaut) : calculé sur la tension instantanée, une
injection à puissance constante a une conductance incrémentale négative qui
déstabiliserait la moindre résonance du réseau.

Les deux **consignes** sont les réglages que l'initialisation cale pour que
l'élément produise la puissance de son bus au load flow (`Pg`, `Qg` de
`Y_DER`, comme un générateur).

```python
from g2egrid.elements.blocs import compiler

c = compiler(schema)
c.etats          # ['x_c', 'x_m', 'vd_mes', 'vq_mes']
for e in c.latex():
    print(e["texte"])   # dx_c/dt = Ki*y_e ...
```

## Un GFM complet, bloc à bloc

{func}`g2egrid.elements.blocs.schema_gfm` dessine un grid-forming à **machine synchrone
virtuelle** (43 blocs) :

- la **VSM** : $2H\,\dot{\Delta\omega} = P_{ref} - p_m - D\,\Delta\omega$, sur la
  puissance mesurée filtrée ($\omega_{f,p}$), puis
  $\dot\theta = \omega_b \Delta\omega$ ; la fréquence du convertisseur est
  $\omega + \Delta\omega$ ;
- le **statisme Q-V** : $v_{ref} = E_{ref} - n_q\, q_m$ ;
- l'**impédance virtuelle** (`zvd`, `zvq`, nulle par défaut) : la chute
  $R_v i_g + L_v(\tfrac{1}{\omega_b}\tfrac{d}{dt} + j\omega) i_g$ retranchée
  des deux consignes de tension, dynamique ou quasi statique
  ({doc}`lois_gfm`) ;
- la **boucle de tension** : PI sur $v_{ref} - v_{ed}$ et $0 - v_{eq}$, plus
  l'anticipation du courant de sortie ($K_{ffi}\, i_g$) et le découplage du
  condensateur ($\mp \omega_{ff} C v_e$) ; elle donne les consignes de
  courant ;
- la **limitation de courant** entre la boucle de tension et la boucle de
  courant, avec rétro-calcul sur les PI de tension ;
- la **saturation de tension** (`vsat`) entre la boucle de courant et le
  filtre : la modulation est ramenée dans $\lvert m\rvert \le M_{max} = 1{,}15$
  (modulation vectorielle), avec rétro-calcul sur les PI de courant ;
- la **boucle de courant**, le filtre LCL dans le repère de la VSM et le
  retour au réseau, comme pour le GFL ;
- deux **sorties de service**, `theta` et `w` (l'angle et la vitesse de la
  VSM).

Ce sont exactement les équations du `GFM2` du moteur (statisme du 2ᵉ ordre),
avec $H = 1/(2 m_p \omega_c)$ et $D = 1/m_p$ : les tests vérifient les mêmes
valeurs propres, à $10^{-7}$ près, sur réseau fort, sur réseau faible, en PV
comme en PQ. Le GFM du moteur a trois consignes ($p_{ref}$, $q_{ref}$,
$v_{e,ref}$) ; le load flow n'en peut caler que deux, et la dynamique ne voit
que $E_{ref} = v_{e,ref} + n_q q_{ref}$ : le schéma n'a donc que $P_{ref}$ et
$E_{ref}$.

Un élément qui écrit `theta` et `w` est un **formeur** pour le moteur : posé
au bus slack, il porte le repère du réseau ($\theta_g$ est son angle) et le
centre d'inertie le compte (avec le poids $2H\,S_n$ s'il a un paramètre `H`).
Le GFM à schéma-bloc tient donc seul un îlot, comme le `GFM2`. Dans l'atelier,
« Poser sur le réseau » le place sur un bus PV (ou slack s'il est le premier
formeur de son île), et un modèle sans `theta` sur un bus PQ.

## La limitation de courant

Comme dans le moteur, la limitation n'est pas un `min` dans les équations :
c'est un **mode**. Chaque `limiteur_courant` ajoute un paramètre drapeau
`mode_<id>`, nul en fonctionnement libre : la consigne passe telle quelle et
le rétro-calcul est nul, si bien que le modèle linéarisé, les valeurs propres
et l'initialisation sont ceux du régime libre. En simulation temporelle, le
limiteur a ses **gardes** : il bascule en mode limité quand le module de la
consigne libre dépasse $I_{max}$, et en ressort sous $I_{max}(1 - 2\,\%)$.
En mode limité, la consigne est ramenée dans $\lvert i\rvert \le I_{max}$,
direction conservée, et `aw` ramène l'intégrale des PI amont sur la consigne
limitée.

Contrairement aux GFM, GFL et MMC du moteur, un élément à schéma-bloc n'est
pas encore vérifié **à l'équilibre** : `solve_equilibrium` ne signale pas un
point de fonctionnement où sa consigne libre dépasse déjà $I_{max}$ (voir la
section *Hors des limites des convertisseurs* de {doc}`analyse_modale`). En
simulation temporelle, ses gardes sont surveillées dès $t = 0$, comme celles
du moteur.

La **saturation de tension** (`saturation_tension`) est le même mécanisme
appliqué à la modulation, comme `limiter_tension` du moteur : drapeau
`mode_vsat` dans le GFM et le GFL complets, $M_{max} = 0$ la désactive (c'est
le défaut des GFL du moteur, et des GFM sur bus continu idéal ; les exemples
la règlent à 1,15, comme le moteur pour un GFM sur bus continu dynamique :
voir {doc}`protections_dc`). Sur
un défaut franc aux bornes d'un GFM, elle s'enclenche avant la limitation de
courant, quelques millisecondes ; il faut alors un pas d'intégration fin
(20 µs en BDF2). Sur un GFL, elle s'enclenche à l'élimination du défaut : la
tension revient, la boucle de courant demande plus que l'onduleur ne sait
faire ; elle relâche quelques millisecondes plus tard (pas de 100 µs).

Son rétro-calcul doit être **rapide**, de l'ordre de la boucle de courant :
$T_{aw,v} = \tau_{CL}/2$, le défaut du moteur comme des exemples (50 µs pour
le GFM, 5 ms pour le GFL, dont la boucle de courant est plus lente). Avec
un rétro-calcul lent (20 ms, l'ancien défaut du moteur), l'intégrateur de
courant ($K_{iCL} \approx 4\cdot10^5$ s⁻¹) s'emballe pendant le défaut, la
modulation reste collée à $M_{max}$ dans une direction qui ne commande plus le
courant, la saturation ne relâche plus et le courant monte à 4 pu (un test le
vérifie) ; sur `wscc9_gfm`, le Newton lâchait même dès 51 ms. Des valeurs
intermédiaires (2 ms, 0,2 ms) le font lâcher aussi : la modulation libre colle
au coin du `min`.

Dans les deux exemples, le limiteur s'appelle `lim` : son drapeau est donc
`mode_lim`, comme celui des GFL et des GFM du moteur. La manœuvre `blocage`
le trouve (elle ramène `Imax` à presque zéro et fige le mode limité), et les
tests vérifient que, sur un défaut, le limiteur du schéma-bloc bascule aux
**mêmes instants** que celui du `GFL1` ou du `GFM2`, avec la même trajectoire
de courant ($10^{-9}$ à $10^{-6}$ pu près). Saturation de tension comprise
(GFM2 et GFL1 réglés de même), un défaut franc aux bornes du GFM donne les
mêmes 33 basculements, à 1 µs près, et le même courant à $3\cdot10^{-4}$ pu
près ; un défaut proche du GFL, les mêmes 6 basculements aux mêmes instants
et le même courant à $10^{-7}$ pu près.

## Dans un réseau

Les modèles vivent dans `Tables.MODELES_BLOCS`, les éléments dans la table
`T_Blocs` (`ID`, `BusID`, `Modele`, `Sn_VA`). Les paramètres du schéma valent par
défaut ; `PARAMS_PAR_BUS` les surcharge bus par bus.

```python
from g2egrid.reseau.tables import Tables, assembler

t = Tables(
    Y_DER=[
        (1, 1, 20e3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (2, 3, 20e3, 1.0, 0.0, 20e6, 5e6, 0.0, 0.0),
    ],
    T_Lines_pu=[(1, 1, 2, 1.0, 0.01, 0.05, 0.0, 20e3, 100e6)],
    T_IBs=[(1, 1, 1e9, 1.0, 0.0, 0.001, 0.01)],
    T_Shunts=[(1, 2, 0.01)],
)
t.MODELES_BLOCS["source_pi"] = schema
t.T_Blocs = [(1, 2, "source_pi", 10e6)]   # ID, bus, modele, Sn_VA
t.PARAMS_PAR_BUS[2] = {"Ki": 40.0}
net = assembler(t)
```

Pour mettre un schéma-bloc **à la place d'un élément natif** du réseau,
`remplacer_par_modele` fait tout le remplacement (voir plus bas) : il n'y a
pas à déplacer la ligne à la main.

## Sous-systèmes

Un bloc `sous_systeme` contient un schéma de même forme que le premier niveau
(`blocs`, `liens`, positions comprises). Ses **ports** sont des blocs
`port_entree` et `port_sortie` (`index` à partir de 0, `nom`) : le
sous-système a autant d'entrées que de `port_entree` et autant de sorties que
de `port_sortie`, dans l'ordre des `index`, nommées par leurs `nom`.
L'imbrication va jusqu'à 8 niveaux. Les paramètres du modèle sont visibles
partout (pas de masque) et les identifiants de blocs ne contiennent pas `__`.

Un sous-système est de l'**organisation**, pas de la physique : `compiler`
dissout d'abord la structure avec `aplatir` (chaque bloc intérieur devient
`<id_sous_système>__<id_intérieur>`, récursivement ; les ports disparaissent
et les liens se raccordent), si bien que `compiler(s)` et
`compiler(aplatir(s))` sont la même compilation. Les noms de variables
portent ces identifiants (`x_regulateur__c`).

`regrouper(schema, ids, nom)` range des blocs d'un même niveau dans un nouveau
sous-système : un port d'entrée par signal qui entre (un seul pour une
sortie source qui alimente plusieurs blocs), un port de sortie par signal
interne utilisé dehors. `degrouper(schema, id)` est son inverse. Les deux
acceptent `chemin=[...]` (les sous-systèmes qui mènent au niveau à éditer).
Elles ne modifient jamais le schéma reçu et n'exigent pas qu'il compile.

```python
from g2egrid.elements.blocs import aplatir, compiler, degrouper, regrouper

groupe = regrouper(schema, ["e", "c"], "regulateur", id="reg")
sous = next(b for b in groupe["blocs"] if b["type"] == "sous_systeme")
ports = [p["nom"] for p in sous["schema"]["blocs"] if p["type"] == "port_entree"]
print(sous["id"], ports)    # reg ['pr', 'm']
assert "x_reg__c" in compiler(groupe).etats
assert compiler(aplatir(groupe)).etats == compiler(groupe).etats
rendu = degrouper(groupe, "reg")
assert {b["id"] for b in rendu["blocs"]} == {b["id"] for b in schema["blocs"]}
```

## Éditer un élément de la bibliothèque : équations → blocs

`depuis_element(classe, parametres=None, detail="blocs", f0=50.0)` écrit le
schéma-bloc d'une classe **native** (`"GFL1"`, `"GFM2"`, `"GFMdVOC"`,
`"IB"`, `"SynchronousMachine"`, `"RLLoad"`...) à partir de ses équations
$\dot x = f$, $0 = g$, $i = h$. Seul un élément à un port alternatif `t` est
pris en charge : un élément à port continu (GFM3DC, MMC...) ou à plusieurs
ports (branche) est refusé avec un message clair.

- un `integrateur` par état, nommé d'après l'état, dont le gain $K$ est le
  facteur paramétrique de l'équation (`wb/L1_pu`, `KiCL`...) ;
- un **sous-système par équation d'état** (`d<état>/dt`) qui calcule la
  dérivée en blocs élémentaires (`somme` avec ses signes, `gain`, `produit`,
  `division`, `fonction` pour sin, cos, sqrt, min, max). `detail="equations"`
  met à la place **un seul bloc `fonction` par équation** ;
- chaque variable algébrique devient un sous-système qui la calcule quand son
  équation se résout en elle (linéaire en elle : toutes les classes de la
  bibliothèque), sinon un bloc `implicite` ;
- les entrées du réseau en blocs `entree`, les courants injectés en `sortie`
  `id`, `iq`, l'angle et la vitesse des formeurs en `theta`, `w`. Les
  **autres sorties de service** (P, Q, Vt, vdc...) **ne sont pas reprises** :
  ce sont des mesures qu'aucune équation du modèle n'emploie ;
- les entrées de service natives en `consigne`, deux calées, les autres non
  calées à une valeur fixe (le choix par classe est `blocs_inverse.CONSIGNES` :
  `p_ref` et `ve_ref` calées et `q_ref` à 0 pour un GFM, puisque la
  dynamique ne voit que `ve_ref + nq.q_ref`) ;
- les paramètres natifs en `parametres` (`kS` et `wb` restent réservés).

**Le modèle généré décrit le fonctionnement LIBRE.** Les drapeaux de mode
(`mode_lim`, `mode_vsat`, `mode_e2`) deviennent des paramètres à 0 et les
gardes de limitation du natif **ne sont pas reprises** : la limitation de
courant et la saturation de tension du GFL ou du GFM natif n'existent pas dans
le schéma généré (on peut les ajouter avec les blocs `limiteur_courant`,
`saturation_tension` ou `saturation`). Il en va de même de l'impédance
virtuelle de limitation (`kRv`) quand `kRv = 0` à la génération.
Conséquence à connaître : ces paramètres (`kRv`, les `mode_*`) restent listés
dans `parametres`, mais plus aucune équation ne les emploie. Mettre ensuite
`kRv = 0.3` dans le schéma **ne change rien** ; pour une impédance virtuelle
de limitation, générer le schéma à partir d'un élément où `kRv > 0`. Le bouton
« Remplacer par son modèle à blocs » de l'atelier rend donc un élément qui ne
se limite plus en simulation temporelle.

Seules les classes à un port alternatif `t` et à au moins deux entrées de
service (les deux consignes calées au load flow) s'écrivent en blocs ; les
autres (`ZLoad`, sans consigne ; `ReactanceShunt`, à argument obligatoire du
constructeur ; tout élément à port continu) sont refusées avec la raison. Un
paramètre inconnu de la classe, ou une valeur qui n'est pas un nombre fini, est
refusé aussi.

**Amorce de l'initialisation.** Le schéma porte
`"amorce": {"classe": "GFL1", "correspondance": {...}, "recales": {...}}` :
`classe_element` appelle d'abord l'`initialize` natif (mêmes paramètres, même
contexte) et s'en sert comme point de départ du Newton ; s'il échoue (schéma
trop modifié), le départ d'avant reprend la main. `recales` liste les
paramètres que l'initialisation native recalcule au point de fonctionnement
(`Kpdc` d'un GFL2, `theta_g0` d'une machine, `R` et `X` d'une charge, la
fréquence `w_set` d'un bus infini) : tant qu'ils valent leur valeur de départ,
ils prennent celle du natif.

Posé à la place du natif (`T_Blocs`), l'élément a le même load flow et les
mêmes valeurs propres (écart $\le 10^{-8}$ relatif, mode de jauge exclu) ; un
gain édité dans le schéma donne le spectre du natif avec ce paramètre changé :

```python
import numpy as np
from g2egrid.elements.blocs import depuis_element
from g2egrid.reseau.catalogue import charger

t = charger("wscc9_cigre6")
gfl = depuis_element("GFL1", f0=t.f0)       # 60 Hz : les gains sont regles pour lui
gfl["parametres"]["Kipll"] *= 2             # on edite un gain du schema
ligne = t.T_GFL1s.pop(0)
t.MODELES_BLOCS["gfl_edite"] = gfl
t.T_Blocs = [(900, ligne[1], "gfl_edite", ligne[2])]   # meme bus, meme Sn_VA
net = assembler(t)
net.run_loadflow()
net.solve_equilibrium()
lam = np.linalg.eigvals(net.build_state_space().A)
print(len(gfl["blocs"]), "blocs au premier niveau,", len(lam), "valeurs propres")
```

Dans l'atelier, `POST /api/blocs/depuis_element` accepte soit
`{classe, parametres?, detail?}`, soit `{tables, table, index, detail?}` (l'élément
d'un réseau, avec **ses** paramètres).

### Remplacer un élément du réseau par son schéma-bloc

L'exemple ci-dessus déplace la ligne à la main, ce qui suffit pour un
convertisseur. `remplacer_par_modele` le fait pour tout élément remplaçable,
et garde de quoi revenir en arrière :

```python
from g2egrid.reseau.equations_reseau import remplacer_par_modele, revenir_au_natif

t = charger("smib_gfl")
t2, nom = remplacer_par_modele(t, ("T_GFL2s", 0), "blocs")
assert t2.T_GFL2s == [] and t2.T_Blocs[0][2] == nom      # "GFL2_GFL2"
assert "blocs" in t2.MODELES_BLOCS[nom]                  # un dessin, editable
assert revenir_au_natif(t2, ("T_Blocs", 0)).T_GFL2s == t.T_GFL2s
```

- le modèle est écrit avec les paramètres de l'élément et rangé sous
  `<classe>_<élément>` ; il garde son `origine` (la ligne native), donc le
  retour au natif marche aussi après avoir enregistré puis rouvert le réseau ;
- la puissance d'une **charge** (`T_Loads`, ou les colonnes `Pl`, `Ql` de
  `Y_DER`) est reportée à son bus en génération négative : le load flow ne
  bouge pas. Il lui faut un bus à elle (un autre élément dont l'injection
  n'est pas programmée au même bus est refusé) ;
- `revenir_au_natif` rend la ligne d'origine ; si le bus a été retouché
  entre-temps (consigne de tension, bus ajoutés), seules les puissances
  déplacées reviennent, les retouches restent.

Dans l'atelier, c'est le bouton « Remplacer par son modèle à blocs » :
`POST /api/blocs/remplacer` `{tables, table, index, detail?}` rend
`{tables, modele, schema, table, index}`.

## La bibliothèque de modèles

Les modèles — schémas-blocs et modèles à équations — se rangent dans une
bibliothèque persistante, partagée entre les réseaux : un fichier JSON par
modèle dans `modeles/`, sous le dossier des sauvegardes du catalogue. Le nom
du fichier est le `nom` du modèle, validé comme un nom de réseau (lettres,
chiffres, `_`, `-`, 64 caractères au plus, aucun chemin). Un modèle mal formé
ou de plus de 2 Mo est refusé.

| route | corps | réponse |
|---|---|---|
| `GET /api/blocs/bibliotheque` | | `[{nom, titre, modifie}]` |
| `GET /api/blocs/bibliotheque/{nom}` | | `{schema}` |
| `POST /api/blocs/bibliotheque` | `{schema_bloc}` | `{nom}` |
| `DELETE /api/blocs/bibliotheque/{nom}` | | `{nom}` |
| `POST /api/blocs/regrouper` | `{schema_bloc, ids, nom, chemin?}` | `{schema}` |
| `POST /api/blocs/degrouper` | `{schema_bloc, id, chemin?}` | `{schema}` |
| `POST /api/blocs/aplatir` | `{schema_bloc}` | `{schema}` |
| `POST /api/blocs/remplacer` | `{tables, table, index, detail?}` | `{tables, modele, schema, table, index}` |

Une erreur de schéma est une 422 avec un message en français ; un modèle qui
n'existe pas dans la bibliothèque est une 404. En code,
`g2egrid_api.analyse.bibliotheque_enregistrer(schema)`, `bibliotheque_lire(nom)`,
`bibliotheque_liste()` et `bibliotheque_supprimer(nom)` font la même chose sans
HTTP.

## Dans l'atelier

L'onglet **Blocs** est l'éditeur : une palette de blocs, un canevas où l'on
glisse les blocs et où l'on tire un lien depuis la sortie d'un bloc (point
noir) jusqu'à l'entrée d'un autre, les champs du bloc choisi et les
paramètres du modèle. **Compiler en équations** affiche les variables et les
équations que le moteur assemblera. **Poser sur le réseau** choisit l'outil
« Élément à schéma-bloc » de la palette du réseau : un clic sur un bus y
raccorde un élément du modèle, sur son propre bus derrière un transformateur.
Les modèles sont enregistrés avec le réseau.
