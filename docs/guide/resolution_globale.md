# Résolution globale et nœuds idéaux

## Toutes les équations, puis une seule résolution

Chaque classe d'équipement donne l'ensemble de ses équations
$(\bm f, \bm g, \bm h)$. Le câblage les referme par une LFT
({doc}`../manuel/02_principes`), et il reste **un** système pour tout le
réseau :

$$
\begin{aligned}
\dot{\bm x} &= A_X\,\bm x + A_Z\,\bm z + A_U\,\bm u \\
\bm 0 &= B_X\,\bm x + B_Z\,\bm z + B_U\,\bm u \\
\bm y &= C_X\,\bm x + C_Z\,\bm z + C_U\,\bm u
\end{aligned}
$$

où $\bm z$ rassemble les variables algébriques de **tous** les équipements :
relations flux/courant des machines, tensions des nœuds idéaux, états rendus
statiques ({doc}`dynamiques`).

```python
dae = net.dae_lineaire()          # DAELineaire : AX ... CU, noms des variables
dae.algebriques[:5]               # ['SM1.i_d', 'SM1.i_q', ...]
ss = dae.etat()                   # = net.build_state_space()
```

{meth}`~g2egrid.reseau.network.DAELineaire.etat` élimine $\bm z$ **globalement**.
Quand $B_Z$ est inversible (indice 1), c'est un complément de Schur :
$A = A_X - A_Z B_Z^{-1} B_X$.

## Le nœud idéal

Un bus sans capacité shunt ne peut pas avoir sa tension pour état. Deux
traitements, au choix (`Network(noeuds_sans_capacite=...)`, `Tables.NOEUDS`) :

| Traitement | Load flow | Modes | Temporel | Usage |
|---|---|---|---|---|
| `"ideal"` (défaut) | exact | exacts, sans mode parasite | BDF2 seulement, choisi de lui-même | tous les réseaux |
| `"resistance"` | shunt fictif $1/R$ | un mode parasite rapide par nœud, vers $-R$ | Radau ou BDF2 | parité MATLAB (`Node.m`) |

La **résistance de causalité** (`ResistiveBus`) pose $v = R\,i_\Sigma$ : un
shunt sans sens physique, qui n'existe que pour rendre le bloc causal. Versé
au load flow, il y consomme $|V|^2/R$ à chaque nœud. MATLAB pose
$R = 100$ pu (`setResistor(100)`, `R_NODE_M`) : jusqu'à 77 % de la charge en
pertes fictives sur `ieee33`, et le mode le moins amorti d'`ieee14_ieee33`
déplacé de 14,4 à 5,5 % (constat FR-005, `tests/test_noeuds_sans_capacite.py`).
Une résistance assez grande pour ne plus fausser le point de fonctionnement
($10^5$ pu) rend le modèle mille fois plus raide et bruite les valeurs
propres proches de zéro. D'où le nœud idéal par défaut ; les réseaux comparés
à MATLAB (`smib`, `wscc9`, `wscc9_gfm`, `wscc9_cigre`, `kundur2area`) gardent
`NOEUDS = "resistance"` (fichier natif : bloc `reglages` ; notebook : cellule
de données). La valeur de $R$ se règle par `Network(r_causalite=...)`,
`Tables.R_NOEUDS` ou bus par bus (`Bus(..., R=...)`).

Le **nœud idéal** (`IdealBus`) s'en passe : sa tension est une variable
algébrique, et son équation est la loi des nœuds elle-même, exacte :

$$
0 = i_{\Sigma,d}, \qquad 0 = i_{\Sigma,q}.
$$

```python
net = Network(f0=50.0, Sb=1e6, noeuds_sans_capacite="ideal")
# ou, par les tables :
t2 = t.copie(NOEUDS="ideal")
```

## Réduction algébrique des courants en série

Quand tous les équipements d'un nœud idéal ont leur courant pour **état** —
une ligne et un transformateur en série —, la loi des nœuds ne porte que sur
des états : $i_1 - i_2 = 0$. Elle ne dit rien de la tension du nœud, et $B_Z$
est singulier : le système est d'**indice 2**.

La résolution globale le traite ainsi :

1. une décomposition en valeurs singulières de $B_Z$ isole les combinaisons
   $N^\top \bm g$ des contraintes qui ne dépendent d'aucune variable
   algébrique : ce sont des contraintes sur les états, $K\bm x = 0$ ;
2. on les remplace par leur **dérivée**, $K(A_X\bm x + A_Z\bm z + A_U\bm u) = 0$,
   qui, elle, détermine $\bm z$ — la tension du nœud est celle qui garde les
   deux courants égaux ;
3. la dynamique obtenue laisse invariante la variété $K\bm x = 0$ : on y
   restreint le modèle en **éliminant un état par contrainte**,
   $\bm x = T\,\bm x_r$. Le choix des états éliminés est un pivot de QR
   **stable** : parmi les colonnes de $K$ dont le résidu atteint 10 % du plus
   grand, la première dans l'ordre des états. Un pivot strict tranchait les
   égalités à l'arrondi près, et le choix basculait d'un point de
   fonctionnement à l'autre : mêmes valeurs propres, mais d'autres états
   nommés, et une continuation (balayage des statismes de CIGRE B4) qui ne
   retrouvait plus ses modes. Les facteurs de participation se lisent sur les
   états gardés : un mode partagé entre un état éliminé et ses voisins se
   reporte sur ces derniers.

```python
ss = net.build_state_space()
ss.elimines        # ['L1.i_d', 'L1.i_q'] : ils ne sont plus des états
ss.T               # x_complet = T . x_reduit
```

Deux branches $R_1 + jX_1$ et $R_2 + jX_2$ en série donnent **exactement** le
spectre d'une branche $(R_1 + R_2) + j(X_1 + X_2)$ ($10^{-11}$ près), et un
nœud idéal est la limite d'une capacité qui s'annule, sans le mode très
rapide qu'elle introduit.

```{admonition} Limites
:class: warning

- **En temporel, BDF2 seulement.** BDF2 résout états et contraintes
  ensemble ; sa matrice de Newton reste inversible pour un pas fini, indice 2
  compris. Radau résout $\bm g = 0$ à $\bm x$ figé, ce qui n'a pas de
  solution quand la contrainte ne porte que sur des états : `simuler`
  (`methode="auto"`, le défaut) et l'atelier prennent BDF2 d'eux-mêmes, et
  `methode="Radau"` est refusé. Les variables algébriques de l'intégration sont gardées
  (`Trajectoire.algebriques`) : les mesures se relisent sans re-résolution.
  Sur `ieee14_ieee33` (140 nœuds idéaux), 0,3 s s'intègrent en 12 s, comme
  avec la résistance de causalité, et les deux réponses s'accordent à
  $6\cdot10^{-4}$ près en relatif.
- **Chaque reprise est réinitialisée** (`g2egrid.noyau.integrateur.reinitialiser`) :
  après une manœuvre ou un basculement de limiteur, deux Euler implicites
  de pas $+h$ puis $-h$ ($h = 10^{-8}$ s) rendent un point cohérent sans
  faire avancer le temps. Retirer une branche parcourue par un courant
  (élimination d'un défaut, ouverture) viole la loi d'un nœud où ne se
  rencontrent que des branches inductives : le premier pas applique
  l'impulsion de tension de la limite idéale, dont le saut de courant
  conserve le flux. Sans cette réinitialisation, l'intervalle repartait de
  la tension d'avant l'événement, et un GFM en schéma-bloc s'écartait de
  celui du moteur de 1,25 pu ; avec elle, les deux s'accordent à
  $3\cdot10^{-8}$ pu, sur les 34 basculements d'un défaut franc. Sur
  `wscc9_faible`, les basculements des limiteurs sont ceux d'une résistance
  de causalité de $10^5$ pu à $10^{-5}$ s près ; les 100 pu de `Node.m`
  chargeaient assez les nœuds pour masquer une saturation de modulation de
  5 ms à l'apparition du défaut.
- La linéarisation **le long d'une trajectoire**
  (`g2egrid.temporel.linearisation`) est refusée : linéariser le réseau
  (`Network.build_state_space`), qui réduit les courants en série.
- Une contrainte qui ferait intervenir une **entrée** demanderait sa dérivée :
  elle est refusée par le modèle d'état. C'est le cas d'un courant d'essai
  injecté à un nœud idéal où ne se rencontrent que des branches inductives :
  l'impédance y croît comme $sL$, elle est impropre. L'impédance vue d'un bus
  (`impedance_bus`) et la $S_{cc}$ dynamique se lisent donc sur le DAE
  linéaire lui-même (`DAELineaire.reponse_frequentielle`), exact à toute
  fréquence finie ; la boucle mineure (`boucle_mineure`), qui demande les
  modèles d'état de $Z_t$, reste refusée à un tel nœud.
- Les sensibilités analytiques ({doc}`sensibilites`) ne traitent pas encore un
  modèle réduit en indice : utiliser `methode="differences"`.
```
