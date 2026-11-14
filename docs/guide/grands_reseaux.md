# Grands réseaux : les modes critiques sans le spectre complet

Le modèle d'état $A$ s'obtient en éliminant les variables algébriques par
complément de Schur ({doc}`resolution_globale`). Cette élimination **remplit** la
matrice, et `eig(A)` coûte $O(n^3)$ pour tout le spectre, alors qu'on ne veut
que la poignée de modes proches de l'axe imaginaire. {mod}`g2egrid.analyse.descripteur`
garde le DAE linéaire tel quel, creux :

$$
E\,\dot x = J\,x, \qquad E = \begin{bmatrix} I & 0 \\ 0 & 0 \end{bmatrix}, \qquad
J = \begin{bmatrix} A_X & A_Z \\ B_X & B_Z \end{bmatrix}.
$$

Les valeurs propres **finies** du faisceau sont exactement celles de $A$ ; les
autres sont infinies, ce sont les contraintes.

## Les modes autour d'un point

```python
from g2egrid import catalogue
from g2egrid import descripteur as dd

net = catalogue.construire("ieee39")
d = dd.descripteur(net)
print(d.n, "inconnues dont", d.nx, "états ; densité de J :", f"{100 * d.densite:.2f} %")

m = dd.modes_proches(d, sigma=-0.5 + 2j * 3.14159 * 1.0, k=6)    # autour de 1 Hz
for lam, p in zip(m.lam, m.p):
    print(f"{lam.real:8.4f} {lam.imag / 6.2832:+7.3f} Hz   {d.etats[p.argmax()]}")
```

```text
462 inconnues dont 402 états ; densité de J : 0.94 %
 -3.4448  +0.838 Hz   SM31.phi_fd
 -3.5185  +0.876 Hz   SM34.phi_fd
 -3.6315  +0.946 Hz   SM38.phi_fd
 -3.2271  +0.742 Hz   SM37.phi_fd
 -3.1044  +0.701 Hz   SM32.phi_fd
 -3.7793  +0.972 Hz   SM39.phi_fd
```

Arnoldi sur $(J - \sigma E)^{-1} E$ rend les valeurs propres les plus proches de
$\sigma$ ; une factorisation LU creuse par point. Les vecteurs à gauche, pour les
participations, viennent d'une itération inverse sur la même LU : les
participations sont celles de {func}`g2egrid.analyse.modal.spectre`, sans avoir formé
$A$. Chaque paire trouvée est vérifiée par son résidu ; $\sigma$ est écarté
d'office du mode de jauge, sur lequel il tomberait exactement en $0$.

## Couvrir une zone du plan complexe

```python
import numpy as np

m = dd.balayage(d, f_max=3.0, profondeur=0.5)      # |Re| <= 0,5 s^-1, f <= 3 Hz
physiques = m.lam[np.abs(m.lam) > 1e-6]            # sans le mode de jauge (0)
print(len(m), "modes ; le plus critique :", np.round(physiques[0], 4))
```

```text
12 modes ; le plus critique : (-0.177-0j)
```

Les modes sont triés par partie réelle décroissante ; le premier est le mode de
jauge, exactement nul, d'un réseau sans bus infini, qu'on écarte ici.

{func}`~g2egrid.analyse.descripteur.balayage` place des points le long de
$\operatorname{Re} = \sigma_r$ et double `k` en chaque point tant que les valeurs
trouvées ne couvrent pas le disque nécessaire : le rectangle demandé est couvert
**entièrement**, multiplicités comprises, quel que soit le nombre de modes qu'il
contient.

| Réseau | États | `balayage` (0–3 Hz) | `eig(A)` dense |
|---|---|---|---|
| ieee39 | 402 | 0,3 s | < 0,1 s |
| ieee39_mixte | 3410 | 1,5 s | 11 à 20 s |

Sur un petit réseau, le calcul dense reste plus rapide ; le faisceau creux
s'impose à partir de quelques milliers d'états.

## Le load flow n'est plus le goulot

Le Newton alternatif de {func}`g2egrid.noyau.loadflow.run_loadflow` est creux depuis
FR-027 : Y-bus en CSC, jacobien polaire de `dSbus_dV` (MATPOWER) formé en
$O(\text{nnz})$, factorisation SuperLU. En dessous de 100 bus (`SEUIL_CREUX`),
le dense vectorisé reste plus rapide (le creux paie environ 1 ms de
construction) et garde l'arithmétique du Newton mixte ; `methode="creux"` ou
`"dense"` force l'une ou l'autre. Mesure (`developpement/bench_loadflow.py`,
`OPENBLAS_NUM_THREADS=1`, médiane) :

| Cas | Bus | Dense (avant) | Auto (après) |
|---|---|---|---|
| ieee14_ieee33 | 146 | 10 ms | 3,7 ms |
| ieee39_mixte | 423 | 125 ms | 8 ms |
| grille synthétique 45×45 | 2 025 | 5,0 s | 0,12 s |
| grille synthétique 70×70 | 4 900 | 73,7 s | 0,38 s |

Les tensions sont celles de l'ancien code à $10^{-13}$ pu (ieee39_mixte), avec
les mêmes itérations.
