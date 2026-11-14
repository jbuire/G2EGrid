# Reseaux de test

Tous les notebooks suivent le format de `cigre_mv.ipynb` : une presentation du
reseau, un bloc de donnees fait uniquement de tables (`Y_DER`, `T_Lines_SI`,
`T_Lines_pu`, `T_Shunts`, `T_Transfos`, `T_IBs`, `T_SMs`, `T_GFM1s` a
`T_GFM4s`, `T_GFL1s`, `T_GFL2s`, `T_Loads`), puis les memes cellules
d'assemblage, de load flow, de modele d'etat et d'analyse modale. Pour un
nouveau reseau, il suffit de copier un notebook et de changer les tables.

| Notebook | Reseau | Etats | Equivalent MATLAB V2 |
|---|---|---|---|
| `smib.ipynb` | machine synchrone sur bus infini | 28 | `SMIB.m` |
| `smib_sm.ipynb` | SMIB elementaire de V1, machine a transfo integre | 32 | `SMIB_V1.m` |
| `smib_gfm.ipynb` | SMIB elementaire de V1, grid-forming | 28 | `SMIB_V1.m` |
| `smib_gfl.ipynb` | SMIB elementaire de V1, grid-following | 27 | `SMIB_V1.m` |
| `three_machine_ib.ipynb` | trois machines sur bus infini | 72 | `ThreeMachineInfiniteBus.m` |
| `wscc9.ipynb` | WSCC 9 bus, trois machines | 93 | `WSCC.m` |
| `wscc9_gfm.ipynb` | WSCC 9 bus, deux onduleurs grid-forming | 85 | `WSCC_GFM.m` |
| `kundur2area.ipynb` | Kundur deux zones, quatre machines | 114 | `Kundur2Area.m` |
| `ieee14.ipynb` | IEEE 14 bus | 181 | `IEEE14.m` |
| `ieee33.ipynb` | IEEE 33 bus, distribution radiale | 195 | `IEEE33.m` |
| `ieee39.ipynb` | IEEE 39 bus, New England | 402 | `IEEE39.m` |
| `cigre_mv.ipynb` | CIGRE MT europeen, trois GFL, 13 charges | 147 | aucun (donnees pandapower) |
| `wscc9_cigre.ipynb` | WSCC 9 bus (SM + GFM + GFL) et distribution CIGRE au bus 5 | 230 | `WSCC_CIGRE.m` |
| `wscc3sm_cigre3.ipynb` | WSCC 9 bus, 3 SM, une distribution CIGRE a chacun des 3 bus de charge | 429 | `WSCC3SM_CIGRE3.m` |
| `datacenter.ipynb` | centre de donnees (redresseur a diodes, onduleur V/f, batterie) et groupe de secours : courbe de charge de 10 minutes, perte du reseau, fatigue de l'arbre | 61 (avec reseau amont) / 54 (site isole) | aucun (blocs propres a la version Python) |

`datacenter.ipynb` fait exception : il n'est pas bati sur des tables mais par le code de
`g2egrid.datacenter` (`reseau_datacenter()`), parce que ses blocs n'existent pas cote MATLAB ;
il dure plusieurs minutes (simulations non lineaires de 2 a 10 minutes de temps simule).

Les nombres d'etats sont ceux de V2, sauf deux cas :

- les trois SMIB elementaires utilisent la charge `RLLoad` (4 etats), la ou
  `SMIB_V1.m` reproduit la charge a 2 etats de V1 ;
- l'IEEE 14 declare en `T_Shunts` le condensateur que `IEEE14.m` fusionne dans
  une charge a Q negatif, ce qui donne une capacite, donc deux etats, au bus 9.

## Transformateur de raccordement

| Equipement | Transfo inclus ? |
|---|---|
| Bus infini (`T_IBs`) | oui : colonnes `R_pu`, `X_pu` |
| Machine synchrone, table a 3 colonnes | non : le declarer dans `T_Transfos` |
| Machine synchrone, table a 5 colonnes `[..., Rt_pu, Xt_pu]` | oui : `SynchronousMachineTR`, courant du transfo en repere rotorique |
| GFM et GFL, table a 4 colonnes | non : la seconde branche du filtre LCL garde ses valeurs par defaut, et le transfo se declare dans `T_Transfos` |
| GFM et GFL, table a 6 colonnes `[..., Rt_pu, Xt_pu]` | oui : `Rt`, `Xt` deviennent la seconde branche du filtre LCL, comme dans V1 |

`T_Transfos` reste disponible pour tout autre transformateur : postes sources,
transfos elevateurs declares a part.

## Charges

| Declaration | Modele cree | Loi de charge |
|---|---|---|
| colonnes `Pl_W`, `Ql_var` de `Y_DER` | `RLLoad` nommee `Load<bus>` | impedance constante (charge classique) |
| `T_Loads` : `(ID, BusID, P_W, Q_var, 'I')` | `ILoad` | courant constant |
| `T_Loads` : `(ID, BusID, P_W, Q_var, 'P')` | `PLoad` | puissance constante |

Les charges classiques ne demandent aucune autre table : l'assemblage cree une
charge R-L a chaque bus dont `Pl_W` ou `Ql_var` est non nul. Une charge
speciale de `T_Loads` s'ajoute a `Y_DER`, y compris dans le load flow : ne pas
la recopier dans `Pl_W`, `Ql_var`. Un meme bus peut porter les deux.

Dans les notebooks fournis, `T_Loads` est vide : toutes les charges sont
classiques.

## Reseaux couples transport / distribution

Deux reseaux couplent un transport a de la distribution. Le second, plus gros,
sert de point de depart a l'etude de la transition.

### wscc3sm_cigre3

La WSCC 9 bus **telle quelle** -- trois machines synchrones, dispatch d'origine,
aucun convertisseur --, dont chacun des trois bus de charge alimente un CIGRE MT
complet (`K = 4`, `3`, `3` reseaux agreges). La charge residuelle du bus WSCC est
ce que la distribution ne porte pas, **en P comme en Q** : chaque bus de charge
garde exactement sa consommation WSCC, donc le point de fonctionnement du
transport est celui du reseau d'origine et les modes electromecaniques sortent a
2,21 et 2,01 Hz, comme dans `wscc9`.

Chaque poste source porte un gradin de condensateurs (`T_Shunts` : 35, 25 et
25 Mvar) qui compense la consommation reactive de ses transformateurs ; sans lui
la moyenne tension descendrait a 0,90 pu. Numerotation : 100-114 pour la
distribution du bus 5, 200-214 pour celle du bus 6, 300-314 pour celle du bus 8.

C'est le reseau de `etudes/transition_ibr.ipynb`, qui lui ajoute des GFL puis
des GFM.

### wscc9_cigre

`wscc9_cigre.ipynb` branche un reseau CIGRE MT complet sur le bus de charge 5 du
WSCC, a travers un poste 230/110 kV a prise (`U2 = 1,025 x 110 kV`). Ses trois
groupes sont de trois technologies : G1 machine synchrone et reference d'angle,
G2 `GFM1` sur bus PV, G3 `GFL1` sur bus PQ.

- **Echelle.** Le bus 5 consomme 125 MW, un reseau CIGRE en porte 25,8 : le
  poste alimente `K = 4` reseaux CIGRE identiques, agreges en un seul
  (puissances x K, impedances / K, donc memes tensions et memes dynamiques en
  pu qu'un seul CIGRE). Les 21,6 MW restants restent une charge au bus 5.
- **Numerotation.** Tout ce qui appartient a la distribution porte un numero
  >= 100, poste compris ; c'est cette convention qui sert a separer transport et
  distribution dans `etudes/transport_distribution.ipynb`.
- **Plan de tension.** Transport entre 1,01 et 1,04 pu, moyenne tension entre
  0,98 et 1,02 pu.

### Charges a Q negatif

Une charge `RLLoad` est une branche R-L : a Q negatif, son inductance serait
negative et ses poles instables. L'assemblage refuse donc `Ql_var <= 0` dans
`Y_DER`. Les bancs de condensateurs se declarent dans `T_Shunts` (les donnees de
Kundur, de l'IEEE 14 et de l'IEEE 39 le font), ou la charge se declare dans
`T_Loads` en type `'P'` ou `'I'`, qui n'ont pas d'inductance et acceptent tout Q.
