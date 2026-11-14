# Robustesse des modes : puissance et N-1

Quels modes bougent le plus quand le point de fonctionnement change ou qu'un
élément déclenche ? {mod}`g2egrid.scenarios.robustesse` refait l'analyse modale sur des
variantes du réseau et suit chaque mode du cas de base.

```python
from g2egrid.reseau.catalogue import charger
from g2egrid.scenarios.robustesse import etudier, resume

e = etudier(charger("ieee14"), f_max=5.0, zeta_max=30, frac=0.1)
print(resume(e))
m = e.modes[0]                         # le plus sensible
m.score("puissance"), m.score("deconnexion"), m.zeta_min, m.pire.perturbation.libelle
```

## Les variantes

| Genre | Variante | Qui reprend la puissance |
|---|---|---|
| `puissance` | une unité non slack produit `frac` de plus | le slack |
| `puissance` | une charge (ou toute la charge) monte de `frac` | les unités réglantes (SM, GFM), au prorata de $S_n$ |
| `deconnexion` | un élément manœuvrable est ouvert (ligne, transformateur, unité, charge, shunt) | les autres unités réglantes, au prorata de $S_n$ |
| `ilotage` | une ouverture qui laisse plusieurs îles vivantes | rendue, mais hors score |

Par convention, un GFM est à batterie et reprend sa part ; un GFL est un
PV, qui reste à sa puissance MPPT et ne reprend rien.

Les ouvertures passent par {func}`g2egrid.reseau.disjoncteurs.en_service` : une île
sans formeur est mise hors tension, une nouvelle référence est désignée si le
slack déclenche.

## L'appariement des modes

Pour chaque variante, les modes physiques ($\operatorname{Im}\lambda \ge 0$)
sont appariés à ceux du cas de base par l'algorithme hongrois, sur le coût

$$
c_{ij} = \Bigl(1 - \sum_{k \in \text{communs}} \sqrt{p_{ik}\, p_{jk}}\Bigr)
       + \min\Bigl(\frac{|\lambda_j - \lambda_i|}{|\lambda_i|},\ 10\Bigr),
$$

où $p_{ik}$ est la participation de l'état $k$ au mode $i$, comparée sur les
états que les deux modèles partagent (par leur nom). La ressemblance des
participations suit un mode quand deux modes voisins se croisent, là où la
seule distance des valeurs propres se trompe. Un mode dont moins de 50 % de
la participation reste sur les états communs **disparaît** : il était porté
par l'élément ouvert.

## Ce qu'on lit

- `score(genre)` : le plus fort $|\Delta\lambda| / |\lambda|$ ;
- `zeta_min` : l'amortissement le plus bas atteint ;
- `pire` : la perturbation qui désamortit le plus ;
- `disparitions` ; au niveau de l'étude, `instables` (variantes instables,
  avec les éléments qui portent le mode instable) et `echecs` (sans point de
  fonctionnement).

## Les limites des convertisseurs

Une variante qui pousse la puissance peut sortir un convertisseur de son
courant maximal $I_{max}$ (ou de sa modulation maximale $M_{max}$). Son spectre
est alors celui du modèle libre, qui ne représente plus ce convertisseur :
réel, il serait en limitation (voir {mod}`g2egrid.noyau.limites` et la section
*Hors des limites des convertisseurs* de {doc}`analyse_modale`). L'étude la
**marque** :

| Champ | Contenu |
|---|---|
| `Etude.depassements` | les convertisseurs hors limite au point de **base** |
| `Etude.hors_limite` | les variantes hors limite, chacune avec ses {class}`~g2egrid.noyau.limites.Depassement` |
| `Effet.hors_limite` | l'effet vient d'une variante hors limite |

Sur `wscc9_gfm`, le GFM du bus 3 (130 MVA) poussé à 150 MW
($\lvert i_s\rvert = 1{,}13$ pu), seule la variante « GFM3 P +10 % » sort de la
limite :

```python
t = charger("wscc9_gfm")
t.Y_DER = [r if r[0] != 3 else (*r[:5], 150e6, *r[6:]) for r in t.Y_DER]
e = etudier(t, deconnexion=False)
print(resume(e))
```

```text
WSCC9_GFM : 2 modes, 6 variantes (puissance +10%, N-1)
  HORS LIMITE apres GFM3 P +10% : GFM3 (bus 3) : courant 1.2415 pu > Imax = 1.2 pu (+3.5 %, base Sn = 130 MVA)
- 1.146 Hz, zeta 47.8 % [SM1 100%, GFM3 0%, GFM2 0%] : score P 0.013, N-1 0.000, zeta min 47.0 %
    charge totale +10% : df -0.017 Hz, dzeta +0.78 pts, score 0.013
    GFM2 P +10% : df +0.011 Hz, dzeta -0.85 pts, score 0.011
    GFM3 P +10% : df +0.011 Hz, dzeta -0.78 pts, score 0.010 [hors limite]
...
```

La variante n'est **ni retirée ni sortie des scores**. Qu'une perturbation de
10 % ou un N-1 sorte un convertisseur de ses limites est en soi un résultat
de robustesse, que le retrait cacherait ; et retirer ces effets des scores
changerait en silence les classements des études déjà faites. Pour les
écarter, filtrer sur `Effet.hors_limite`. Un seul avertissement
{class}`~g2egrid.noyau.limites.DepassementLimite` résume l'étude.

L'étude se lance aussi depuis l'onglet **Robustesse** de l'atelier. Les
constats sur le catalogue sont résumés dans la page de
{doc}`l'atelier <../prise_en_main/atelier>`, section *Robustesse des modes*.
