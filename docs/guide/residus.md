# Résidus, grammiens et rapport par mode critique

Les participations disent quels **états** font un mode. Elles ne disent ni par
où l'exciter, ni où le voir, ni ce qui en déplace la valeur propre. C'est
pourtant ce qui sert à placer un stabilisateur, un amortisseur de
convertisseur ou une mesure. {mod}`g2egrid.analyse.residus` y répond, depuis le même
modèle d'état $(A, B, C, D)$.

## Le rapport d'un mode critique

```python
from g2egrid import catalogue, residus

net = catalogue.construire("wscc9")
ss = net.build_state_space()
r = residus.rapport(ss, n_modes=2, f_min=0.5, f_max=3.0)    # les deux moins amortis
print(r[0].texte())
```

```text
### Mode 59 : -1.207 + 13.82j  (2.2 Hz, zeta 8.70 %)

Elements : SM3 72 %, SM2 26 %, SM1 2 %

Etats : SM3.Dwr 25 %, SM3.Pm 12 %, SM3.theta 10 %, SM2.Dwr 9 %, SM2.Pm 4 %

Forme modale : SM3.Dwr 1.00/-0 deg, SM3.Pm 69.76/+105 deg, SM3.theta 27.17/-95 deg, SM2.Dwr 0.45/+172 deg, SM2.Pm 31.09/-82 deg

Commandable par : SM3.v_ref 100 %, SM2.v_ref 90 %, SM3.p_ref 21 %, w_ref 19 %

Observable sur : T3_3_9.p_f 100 %, SM3.Ce 79 %, SM3.P 78 %, SM3.igd 75 %

Indice de dominance : 14.6

| sortie | entree | abs(R) | phase |
|---|---|---|---|
| T3_3_9.p_f | SM3.v_ref | 5.51 | +145 |
| T3_3_9.p_f | SM2.v_ref | 4.94 | -42 |
| SM3.Ce | SM3.v_ref | 4.33 | -36 |
...

Coefficients de A les plus influents : dSM3.Dwr/dSM3.Pm (26), dSM2.Dwr/dSM2.Pm (9.6), ...
```

Le mode inter-machines à 2,2 Hz oppose SM3 et SM2 (phases 0 et 172°). Il
s'excite le mieux par les consignes de tension des deux machines, se voit le
mieux sur le transit du transformateur T3, et c'est la paire (`T3_3_9.p_f`,
`SM3.v_ref`) qui y porte le plus fort résidu.

Chaque {class}`~g2egrid.analyse.residus.RapportMode` porte, en plus du texte :
`elements` (participations cumulées par élément), `etats`, `forme` (module et
phase relatifs à l'état dominant), `entrees` (commandabilité), `sorties`
(observabilité), `residus` (les plus forts), `sensibilites` (coefficients de
$A$) et `dominance`. `critere="dominance"` classe les modes par leur présence aux
entrées et sorties retenues (`entrees=`, `sorties=`) plutôt que par leur
amortissement.

## Résidus

$$
H(s) = D + \sum_k \frac{R_k}{s - \lambda_k}, \qquad R_k = (C v_k)(w_k B).
$$

Le résidu $R_k[i, j]$ est le poids du mode $k$ dans le transfert de l'entrée
$j$ vers la sortie $i$. Grand, le mode est à la fois commandable par $j$ et
observable sur $i$, et c'est là que se boucle un amortisseur : au premier ordre,
un gain $K$ de $i$ vers $j$ déplace la valeur propre de $R_k[i, j]\,K$.

```python
R = residus.residus(ss, entrees=["SM3.v_ref"], sorties=["T3_3_9.p_f"])   # (modes, 1, 1)
d = residus.dominance(ss)                 # sigma_max(R_k) / |Re lambda_k|
print(R.shape, round(float(d[r[0].indice]), 2))
```

La commandabilité $\lvert w_k B_j\rvert$ et l'observabilité $\lvert C_i v_k\rvert$
dépendent de la normalisation des vecteurs propres, mais pas leur produit : on
les rend en % du plus grand, entrée contre entrée dans un même mode
({func}`~g2egrid.analyse.residus.commandabilite`, {func}`~g2egrid.analyse.residus.observabilite`).
Les résidus portent les unités des entrées et des sorties : on ne compare que
des paires homogènes.

## Grammiens et valeurs de Hankel

```python
g = residus.grammiens(ss)
print(len(g.lam), "modes stables ;", len(g.exclus), "écarté(s) (jauge) ;",
      "ordre à 1e-3 :", g.ordre(1e-3))
Wc, Wo = g.physique()                     # dans la base des états
```

```text
92 modes stables ; 1 écarté(s) (jauge) ; ordre à 1e-3 : 67
```

Les grammiens n'existent que pour un $A$ stable. Or un réseau sans bus infini
porte un mode de jauge nul, et un tirage de campagne peut être instable. On les
calcule donc dans la base modale, sur la seule partie stable, où ils ont une
forme close :

$$
W_c[i, j] = -\frac{(B_m B_m^H)_{ij}}{\lambda_i + \bar\lambda_j}, \qquad
W_o[i, j] = -\frac{(C_m^H C_m)_{ij}}{\bar\lambda_i + \lambda_j},
$$

avec $B_m = W B$ et $C_m = C V$. Les modes écartés sont rendus dans `exclus`.
Les valeurs singulières de Hankel (`g.hankel`) ne dépendent pas de la base ;
leur décroissance dit combien d'états un modèle réduit équilibré devrait garder
(`g.ordre(tol)`). Les tests les comparent à la résolution de Lyapunov de SciPy.

## Dans l'atelier

Le **rapport d'étude** (bouton *Rapport*) contient le rapport des trois modes
critiques ; voir {doc}`automatiser`.
