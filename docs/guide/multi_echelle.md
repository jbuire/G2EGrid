# Transport et distribution : modèles réduits

Faut-il le transport pour étudier la distribution, et la distribution pour
étudier le transport ? {mod}`g2egrid.multi_niveau.multi_echelle` construit de quoi y
répondre, et l'onglet **Transport / distribution** de l'atelier en donne les
mêmes outils.

## Les réseaux

| Réseau | Construction |
|---|---|
| `wscc9_cigre6` | WSCC9 ; chaque charge (bus 5, 6, 8) remplacée par **deux** réseaux CIGRE MV à GFL, chacun agrégeant $K = 3$ départs |
| `ieee14_ieee33` | IEEE14 ; les charges des bus 9, 10, 13, 14 remplacées par des départs IEEE33, trois GFL (PV) en bout d'antenne |

Chaque distribution est une **zone à elle** : on peut y ajouter des GFL
indépendamment. Les charges sont tirées bus par bus (±20 %, total conservé)
pour que deux zones au même poste ne soient pas identiques, puis mises à
l'échelle pour reprendre 98 % de la charge du poste. Les deux réseaux sont au
catalogue.

```python
from g2egrid.multi_niveau.multi_echelle import wscc9_cigre6, ajouter_gfl, evaluer

r = wscc9_cigre6()
r, ajouts = ajouter_gfl(r, n=10, graine=1)   # 10 GFL à des bus tirés au hasard
pal = evaluer(r)                              # les modèles réduits contre le complet
pal.transport.ecart_max(), pal.transport_hybride.ecart_max()
[c.ecart_max() for c in pal.distribution]
```

## Les modèles réduits

**Transport seul** — chaque distribution devient la puissance qu'elle soutire,
lue sur le load flow du modèle complet (même point de fonctionnement) :

- en charge **R-L** (`comparer_transport`) ;
- en **hybride** (`comparer_transport(..., hybride=True)`) : la production GFL
  en puissance constante, le reste en R-L (`equivalent_hybride`).

**Distribution seule** (`distribution_seule`, `comparer_distribution`) — une
instance avec son transformateur de poste ; au bus amont, un bus infini
derrière l'impédance de Thévenin du reste du réseau
(`impedance_thevenin` : Y-bus, charges en admittance, machines derrière $X'_d$,
GFM derrière leur filtre, GFL ignorés). `zth=False` : bus infini raide.

## La partition de l'état

Le modèle complet est un seul modèle d'état, $\dot x = A x + B u$. Ses états
se partagent entre ceux du transport, $x_T$, et ceux de la distribution,
$x_D$ :

$$
\frac{d}{dt}\begin{bmatrix} x_T \\ x_D \end{bmatrix} =
\begin{bmatrix} A_{TT} & A_{TD} \\ A_{DT} & A_{DD} \end{bmatrix}
\begin{bmatrix} x_T \\ x_D \end{bmatrix} + B u .
$$

Étudier le transport seul, c'est garder $\dot x_T = A_{TT} x_T + A_{TD} x_D +
B_T u$ en faisant passer $x_D$ **du côté des entrées** : les modes du
transport sont alors les valeurs propres de $A_{TT}$, et le couplage
$A_{TD}$, $A_{DT}$ est exactement ce que la séparation néglige. Deux façons
de traiter l'autre côté :

- **figé** — $x_D$ est une entrée constante : modes de $A_{TT}$ ;
- **quasi statique** — $\dot x_D = 0$, donc $x_D = -A_{DD}^{-1}(A_{DT} x_T +
  B_D u)$ : modes du complément de Schur $A_{TT} - A_{TD} A_{DD}^{-1}
  A_{DT}$. Les équivalents de charge (R-L, hybride) en sont des
  approximations.

Symétriquement pour la distribution. C'est la perturbation singulière : un
côté RAPIDE se rend quasi statique, un côté LENT se fige. Sur WSCC9 + 6 CIGRE,
le transport avec la distribution quasi statique est à 0,4 % du complet
(8,8 % figée) ; la distribution avec le transport figé à 1,3 % (32 % quasi
statique).

```python
from g2egrid.multi_niveau.multi_echelle import partition, comparer_partition

P = partition(r)                      # (Partition, spectre complet)
c = comparer_partition(r, "transport", quasi_statique=True, part=P)
print(100 * c.ecart_max())
```

## L'appariement des modes

Les deux modèles comparés ont la même numérotation (on découpe le réseau
aplati), donc les mêmes noms d'états. Les modes sont appariés par
l'algorithme hongrois sur la ressemblance de leurs participations et la
distance des valeurs propres (voir {doc}`robustesse`). On retient les modes
que la partie étudiée porte à plus de 50 %, oscillants, de fréquence propre
$|\lambda|/2\pi \le 100$ Hz et amortis à moins de 95 %. Les modes que le
modèle réduit n'a pas (`absents()`) se comptent à part : côté distribution,
ce sont les modes **collectifs** partagés entre départs voisins d'un même
poste.

## Le balayage de pénétration

`penetration(r, n_max, pas, graine)` ajoute des GFL par paquets et, à chaque
palier, confronte les modèles réduits au complet ; `modes_interaction` liste
les modes portés à la fois par le transport et la distribution.

L'étude complète du papier se rejoue par :

```bash
cd g2egrid
uv run python etudes/transport_distribution/etude.py      # resultats.json
uv run --extra plot python etudes/transport_distribution/figures.py
uv run python etudes/transport_distribution/tables_latex.py
cd etudes/transport_distribution/papier && pdflatex papier.tex
```
