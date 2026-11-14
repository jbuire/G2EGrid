# Puissance de court-circuit et SCR à chaque nœud

{mod}`g2egrid.analyse.court_circuit` estime la puissance de court-circuit $S_{cc}$ de
chaque nœud alternatif, et le rapport de court-circuit (SCR) des nœuds qui
portent des GFL.

## Le calcul statique

C'est la méthode de la CEI 60909-0 pour le courant de court-circuit initial
symétrique $I_k''$ : un défaut triphasé franc, une source de tension $c \cdot 1$ pu
au point de défaut, le réseau réduit à son Y-bus de branches (lignes et
transformateurs ; capacités et charges négligées), les sources de tension
derrière leur impédance. Ce n'est **pas** un calcul normatif complet : ni $c_{min}$,
ni courant de crête, ni courant coupé, ni composantes symétriques, et les
rapports de tensions nominales sont supposés égaux à 1.

| équipement | hypothèse |
|---|---|
| bus infini | derrière son $R + jX$ |
| machine synchrone, `machine="soustransitoire"` (défaut) | derrière sa réactance **sous-transitoire** $x''_d = L_l + 1/(1/L_{ad} + 1/L_{fd} + 1/L_{1d})$, avec les facteurs de correction ci-dessous |
| machine synchrone, `machine="transitoire"` | derrière sa réactance transitoire $x'_d = L_l + L_{ad}L_{fd}/(L_{ad}+L_{fd})$, sans facteur de correction (l'ancien calcul) |
| GFM, `gfm="tension"` (défaut) | source de tension derrière son filtre ($R_1+R_2$, $L_1+L_2$) |
| GFM, `gfm="courant"` | source de courant limitée à $I_{max} S_n$ |
| GFL, `gfl="aucun"` (défaut) | rien : c'est la convention du SCR |
| GFL, `gfl="courant"` | source de courant $I_{max} S_n$ (CEI 60909-0:2016) |

La réactance transitoire ignore l'amortisseur de l'axe d : elle **sous-estime**
$S_{cc}$ près d'une machine (sur le catalogue, $S_{cc}$ médiane de 1,01 à
1,21 fois plus grande avec x''d et les facteurs, jusqu'à 1,40 aux bornes d'une machine).
La réactance sous-transitoire est celle de la CEI 60909 pour $I_k''$.

### Facteurs de correction d'impédance

Avec $c_{max} = c$ (l'argument `c`, 1 par défaut ; 1,1 pour le calcul normatif
de $I_k''$ maximal), et $\sin\varphi_{rG}$ tiré de `cos_phi` (0,9 par défaut :
valeur courante, le catalogue ne porte pas le facteur de puissance nominal des
machines). `correction_cei=False` les désactive.

| où | facteur |
|---|---|
| générateur raccordé directement au réseau | $K_G = c / (1 + x''_d \sin\varphi_{rG})$, sur son impédance |
| bloc centrale (machine seule sur son bus, derrière son transformateur élévateur, intégré à la machine ou dans `T_Transfos`), régleur en charge | $K_S = c / (1 + \lvert x''_d - x_T\rvert \sin\varphi_{rG})$, sur l'impédance du bloc (machine **et** transformateur) |
| transformateur de réseau à deux enroulements | $K_T = 0{,}95\,c / (1 + 0{,}6\,x_T)$ |

Un défaut **aux bornes** de la machine d'un bloc centrale se calcule avec $K_G$
sur la machine et le transformateur sans correction (les autres blocs gardent
$K_S$). Les facteurs appliqués sont rendus dans `hypotheses["facteurs_cei"]`.

Les formules sont celles de l'implémentation de pandapower
(`pandapower.shortcircuit`), avec laquelle le module s'accorde à $10^{-5}$ près
sur WSCC 9 bus et $10^{-9}$ près sur un réseau à générateur direct et
transformateur de réseau (tests `test_*_contre_pandapower`) ; elles n'ont **pas**
été relues dans la norme elle-même : *à vérifier*. Le bloc centrale **sans**
régleur en charge ($K_{SO}$) n'est pas traité.

Avec $Z = Y^{-1}$, le courant de défaut au nœud $i$ est

$$
I_k(i) = \frac{c + \sum_j |Z_{ij}|\, I_j}{|Z_{ii}|},
\qquad S_{cc}(i) = I_k(i)\, S_b ,
$$

le second terme étant celui des sources de courant $I_j$. Le rapport $X/R$ est
celui de $Z_{ii}$.

## Le SCR

Pour un nœud qui porte des GFL de puissance nominale $S_n$ :
$\mathrm{SCR} = S_{cc}/S_n$. Le réseau est dit **faible** sous 3 et **très
faible** sous 2 (NERC, *Integrating inverter-based resources into low short
circuit strength systems*, 2017). Des GFL **voisins** s'affaiblissent l'un
l'autre : le SCRIF (*SCR with interaction factors*, même guide) en tient compte
nœud par nœud,

$$
\mathrm{SCRIF}_i = \frac{S_{cc,i}}{P_i + \sum_{j \ne i} \mathrm{IF}_{ji}\, P_j},
\qquad \mathrm{IF}_{ji} = \left|\frac{Z_{ji}}{Z_{ii}}\right| ,
$$

$\mathrm{IF}_{ji}$ étant la variation de tension au nœud $j$ pour une variation
au nœud $i$. Le WSCR de l'ERCOT ($\sum S_{cc,i} P_i / (\sum P_i)^2$) n'est pas
rendu pour tout le réseau : il est fait pour une grappe de producteurs
proches, et sur un grand réseau à GFL dispersés il divise le SCR par leur
nombre.

Le SCR critique d'un GFL dépend de son réglage : avec le réglage par défaut du
moteur (PLL de 50 ms), il est voisin de 3.

## La vue par l'impédance

`dynamique=True` ajoute $S_b/\sigma_{max}(Z_{dq}(j 2\pi f))$ : l'impédance dq de
**tout** le réseau linéarisé au point de fonctionnement, vue de chaque nœud, à
$f = 1$ Hz dans le repère dq (autour de 50 Hz). Elle compte ce que le calcul
statique ignore : les boucles des convertisseurs, la PLL des GFL et les
charges. C'est la force que « voit » un nouvel équipement.

```python
from g2egrid.reseau.catalogue import charger
from g2egrid.analyse.court_circuit import puissance_court_circuit

r = puissance_court_circuit(charger("cigre_mv"), dynamique=True)
for n in r.noeuds:
    print(n.bus, n.Scc / 1e6, n.XR, n.SCR, n.Scc_dyn)
print("SCR min", r.SCR_min, "SCRIF min", r.SCRIF_min)
```

## Dans l'atelier

Onglet **Point de fonctionnement**, carte « Puissance de court-circuit à
chaque nœud » : les hypothèses (GFM, GFL, charges, vue par l'impédance), puis
**Estimer Scc**. Le schéma colore chaque nœud à GFL par son SCR (rouge foncé
sous 2, rouge de 2 à 3, orange de 3 à 5, vert au-delà) et étiquette les
autres de leur $S_{cc}$ ; la table donne $U_n$, $S_{cc}$, $I_k$, $X/R$, la
puissance des GFL, le SCR et le SCRIF.

## Des réseaux faibles au catalogue

Trois réseaux, stables au point livré mais proches de leur limite, ont toutes
leurs centrales GFL raccordées sous SCR = 3. Ils sont reconstruits et vérifiés
par `etudes/reseaux_faibles/construire.py`, dont le README détaille chaque cas.

| réseau | genre | centrales GFL | SCR | SCRIF |
|---|---|---|---|---|
| `wscc9_faible` | transport, 9 bus | 2 (machines des bus 2 et 3 remplacées) | 1,47 et 2,19 | 1,02 et 1,10 |
| `ieee39_faible` | transport, 39 bus | 5 (machines des bus 32 à 36 remplacées) | 2,07 à **3,06** | 0,76 à 1,38 |
| `cigre_mv_faible` | distribution | 3 PV en bout de départs, amont de 100 MVA | 1,88 à 2,65 | 0,87 à 1,20 |

`cigre_mv_faible` a 16 % de pertes de branche (4,1 MW pour 25,8 MW de charge) :
c'est voulu, les trois PV (28,8 MW) produisent plus que la charge et absorbent
du réactif ; voir le README de `etudes/reseaux_faibles`.

Ces valeurs sont celles du calcul par défaut (x''d et facteurs de la CEI 60909).
Les réseaux ont été construits avec x'd (`machine="transitoire"`) : 1,30 et 1,95
pour `wscc9_faible`, 2,00 à 2,89 pour `ieee39_faible`, 1,85 à 2,61 pour
`cigre_mv_faible`. Au bus 32 de `ieee39_faible`, le SCR de la centrale passe de
2,89 à 3,06 : elle n'est plus « faible » au sens de la convention par défaut.
