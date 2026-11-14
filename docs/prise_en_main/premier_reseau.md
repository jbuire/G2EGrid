# Premier réseau

Ce tutoriel construit une machine synchrone raccordée à un réseau amont
(*SMIB*), calcule son point de fonctionnement, son modèle d'état et ses modes,
puis passe au catalogue des réseaux de test.

## Déclarer le réseau

Un réseau se déclare comme dans les tables MATLAB : d'abord les bus, puis les
branches et les équipements.

```python
from g2egrid import Network
from g2egrid.elements import SynchronousMachine, Line, VoltageSource

net = Network(f0=50.0, Sb=1e6)
net.add_bus(1, 20e3, btype=2, V=1.0, Pg=0.8, C=1e-3)   # PV : la machine
net.add_bus(2, 20e3, btype=1, V=1.0, C=1e-3)           # slack : réseau amont
net.add(SynchronousMachine("SM1", 1, Sn=1e6, Sb=1e6, wb=net.wb))
net.add(Line("L1", 1, 2, R=0.01, X=0.20, wb=net.wb))
net.add(VoltageSource("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
```

- `btype` : 1 = slack, 2 = PV, 3 = PQ ;
- `add_bus` prend les puissances **en pu** base système ; `net.add(Bus(...))`
  les prend en watts, comme la table `Y_DER` ;
- `C` est la capacité shunt du bus (pu). Elle fait de la tension du bus un
  **état**, et la loi des nœuds est exacte. Un bus sans capacité devient un
  nœud à résistance de causalité, ou un nœud idéal
  ({doc}`../guide/resolution_globale`).

## Point de fonctionnement et modèle d'état

```python
net.run_loadflow()            # Newton-Raphson polaire
net.solve_equilibrium()       # équilibre de chaque équipement
ss = net.build_state_space()  # A, B, C, D
print(ss)                     # <StateSpace 26 etats, 6 entrees, 23 sorties | ...>
```

`solve_equilibrium` est appelé de lui-même par `build_state_space` s'il ne
l'a pas été. Le modèle d'état est **nommé** :

```python
ss.states[:4]                  # ['bus1.vd', 'bus1.vq', 'bus2.vd', 'bus2.vq']
ss.iu("SM1.p_ref")             # indice de l'entrée
ss.iy("SM1.P")                 # indice de la sortie
ss.eig()                       # valeurs propres
sys = ss.to_control()          # objet python-control (option "control")
```

## Les modes

```python
from g2egrid import modal

for m in modal.modes(ss, f_max=5, oscillants=True):
    print(f"{m.f:6.3f} Hz  zeta = {m.zeta:5.1f} %  {m.dominant}")
```

```text
 2.065 Hz  zeta =   9.8 %  SM1.Dwr
 0.621 Hz  zeta =  59.2 %  SM1.phi_fd
```

Le premier est le mode électromécanique de la machine contre le réseau ; le
second est lié à son circuit d'excitation. Chaque {class}`~g2egrid.analyse.modal.Mode` donne sa valeur propre, sa fréquence,
son amortissement et les états qui y participent le plus. Le
{doc}`../guide/analyse_modale` va plus loin : forme modale, sensibilités,
réponses.

## Les réseaux de test

Quatorze réseaux sont livrés, décrits par des tables (un notebook chacun dans
`reseaux_test/`) :

```python
from g2egrid import catalogue

[f.nom for f in catalogue.lister()]
# ['smib', 'smib_sm', 'smib_gfm', 'smib_gfl', 'three_machine_ib', 'wscc9',
#  'wscc9_gfm', 'kundur2area', 'ieee14', 'ieee33', 'ieee39', 'cigre_mv',
#  'wscc9_cigre', 'wscc3sm_cigre3', 'hvdc2t', 'mtdc3t', 'wscc3sm_zones']

net = catalogue.construire("wscc9")     # assemblé, load flow et équilibre faits
ss = net.build_state_space()
for m in modal.modes(ss, f_max=5, oscillants=True)[:3]:
    print(f"{m.f:6.3f} Hz  zeta = {m.zeta:5.1f} %  {m.dominant}")
```

```text
 2.200 Hz  zeta =   8.7 %  SM3.Dwr
 2.011 Hz  zeta =  11.3 %  SM1.Dwr
 1.841 Hz  zeta =  19.9 %  SM1.Dwr
```

Pour modifier un réseau du catalogue, on part de ses **tables**
({doc}`../guide/construire`) :

```python
t = catalogue.charger("wscc9")          # les tables
t2 = t.copie(GFM_PARAMS={"mp": 0.02})   # une variante
```

Les résultats de ces réseaux sont commentés au chapitre
{doc}`../manuel/07_exemples`.
