# Analyse modale et simulation

Toutes les fonctions de ce chapitre prennent un
{class}`~g2egrid.reseau.network.StateSpace` — celui de `net.build_state_space()`.

## Spectre et modes

```python
from g2egrid import catalogue, modal

net = catalogue.construire("kundur2area")
ss = net.build_state_space()

ss.is_stable()                              # tient compte du mode de jauge
lam, Vr, P = modal.spectre(ss)              # valeurs propres, vecteurs, participations
f, zeta = modal.frequence(lam), modal.amortissement(lam)
for m in modal.modes(ss, f_max=2, oscillants=True):
    print(m.indice, f"{m.f:.3f} Hz", f"{m.zeta:.1f} %", m.etats)
```

{func}`~g2egrid.analyse.modal.modes` classe les modes du moins amorti au plus amorti et
filtre en fréquence (`f_max`), en amortissement (`zeta_max`) et sur les seules
paires oscillantes. `m.indice` est le rang du mode dans `eig(A)` : c'est lui
qu'attendent les fonctions suivantes.

```{note}
Un réseau sans bus infini a une valeur propre nulle : faire tourner tous les
angles ensemble ne change rien à ses équations. Ce **mode de jauge** n'est pas
une instabilité (chapitre {doc}`../manuel/02_principes`).
```

## Hors des limites des convertisseurs

Les limiteurs des convertisseurs (courant $\lvert i\rvert \le I_{max}$,
modulation $\lvert m\rvert \le M_{max}$) sont des **modes discrets**,
surveillés par la simulation temporelle ({mod}`g2egrid.elements.limiteurs`).
Le modèle linéarisé, lui, est toujours calculé en mode **libre**. Tant que le
point de fonctionnement est dans les limites, c'est exact ; au-delà, le
convertisseur réel serait en limitation, et le spectre calculé ne le décrit
plus.

`solve_equilibrium` le vérifie ({mod}`g2egrid.noyau.limites`). Chaque limiteur dont
la grandeur dépasse sa limite au point d'équilibre est rangé dans
`net.depassements`, et signalé par un avertissement
{class}`~g2egrid.noyau.limites.DepassementLimite`. Sur `wscc9_gfm`, en demandant
200 MW au GFM du bus 3 (130 MVA) :

```text
DepassementLimite: point d'equilibre au-dela d'une limite de convertisseur ; le modele
linearise (mode libre) ne represente plus l'equipement, qui serait en limitation :
  GFM3 (bus 3) : courant 1.4998 pu > Imax = 1.2 pu (+25.0 %, base Sn = 130 MVA)
```

Le calcul continue : on peut explorer le modèle libre au-delà de la limite,
en le sachant. `warnings.simplefilter("error", DepassementLimite)` rend
l'avertissement bloquant, `"ignore"` le fait taire ; `net.depassements` reste
rempli dans les deux cas.

- **La grandeur comparée** est celle que lit la garde du limiteur, par le même
  calcul (`grandeurs_limitees`) : le module de la consigne de courant libre,
  qui vaut à l'équilibre le courant **côté convertisseur** $\lvert i_s\rvert$
  (à travers $L_1$), en pu de $S_n$. Ce n'est pas le courant côté réseau
  $\lvert i_g\rvert$, qui en diffère du courant du condensateur du filtre.
  La modulation libre $\lvert e\rvert / v_{dc}$ face à $M_{max}$ n'est
  comparée que si la saturation de tension est réglée ($M_{max} > 0$).
- **Le seuil**, $10^{-9}$, est celui de la simulation pour une garde déjà
  franchie au départ : un dépassement est signalé si et seulement si la
  simulation, limiteurs surveillés, commencerait en limitation.
- **Les familles couvertes** sont celles qui héritent de
  `LimiteurCourant` : GFM, GFL, MMC, et leurs variantes à port continu. Les
  limiteurs des schémas-blocs ({doc}`schemas_blocs`) ne sont pas encore
  vérifiés à l'équilibre.

La {doc}`continuation <continuation>` et la {doc}`robustesse <robustesse>`
marquent leurs points et variantes hors limite ; la validation de l'atelier
en fait un avis.

## Forme modale et sensibilités

```python
k = m.indice
modal.forme_modale(ss, k, n=6)       # (état, module relatif, phase [deg])
modal.sensibilite_top(ss, k, n=8)    # coefficients A[i, j] qui déplacent le mode
```

La forme modale dit **qui** oscille et en quelle phase : deux vitesses rotor à
180° signent un mode inter-machines. La sensibilité aux coefficients de $A$
désigne une **équation** ; pour remonter jusqu'au **paramètre**, voir
{doc}`sensibilites`.

## Grands réseaux

Au-delà de quelques milliers d'états, calculer tout le spectre coûte cher.
{func}`~g2egrid.analyse.modal.valeurs_propres_partielles` rend les `k` valeurs propres
les plus proches d'un point `sigma` du plan complexe (Arnoldi en
shift-and-invert). {mod}`g2egrid.multi_niveau.separation` sépare un modèle en parties lente
et rapide, ou le réduit par troncature équilibrée.

## Réponses linéaires

```python
import numpy as np
from g2egrid.multi_niveau.separation import reponse_indicielle

t = np.linspace(0, 5, 1001)
u = np.zeros(len(ss.inputs)); u[ss.iu("SM1.p_ref")] = 0.01
x, parts = modal.reponse_indicielle_modale(ss, u, t, par_mode=True)   # mode par mode
```

## Simulation non linéaire

Le même réseau, les mêmes équations, **sans linéarisation** :

```python
from g2egrid.temporel import ModeleNonLineaire, Evenement, Manoeuvre, mesures

mnl = ModeleNonLineaire(net)
tr = mnl.simuler((0.0, 5.0), uext_fn=mnl.echelon("SM1.p_ref", 0.05, t0=0.1),
                 t_eval=np.linspace(0, 5, 501))
tr.etat("SM1.Dwr")                                   # une trajectoire d'état
mes = mesures(mnl, tr, uext_fn=mnl.echelon("SM1.p_ref", 0.05, t0=0.1))  # V, f, P, Q
```

Un scénario enchaîne plusieurs {class}`~g2egrid.temporel.Evenement` ; une
{class}`~g2egrid.temporel.Manoeuvre` ouvre une ligne ou déclenche une machine
en cours de simulation. La comparaison des deux voies (linéaire et non
linéaire) est décrite au chapitre {doc}`../manuel/12_atelier`.
