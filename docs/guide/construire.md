# Construire un réseau

Un réseau se décrit de deux façons, qui aboutissent au même objet
{class}`~g2egrid.reseau.network.Network` :

- **en code**, élément par élément (`net.add(...)`) ;
- **par ses tables** (`Y_DER`, `T_Lines_SI`, `T_SMs`...), exactement celles des
  notebooks de `reseaux_test/` et des scripts MATLAB, assemblées par
  {func}`g2egrid.reseau.tables.assembler`. C'est la voie du catalogue et de l'atelier.

## En code

```python
from g2egrid import Network
from g2egrid.elements import Bus, Line, Transformer, SynchronousMachine, RLLoad, IB

net = Network(f0=50.0, Sb=100e6)
net.add(Bus(1, 20e3, btype=2, V=1.0, Pg=400e6))          # watts et vars, comme Y_DER
net.add(Bus(2, 230e3, btype=3, Pl=100e6, Ql=20e6, C=0.05))  # C : capacité shunt, pu
net.add(Bus(3, 230e3, btype=1, V=1.0, C=0.05))
net.add(Transformer("T1", 1, 2, R=0.0, X=0.02, wb=net.wb))
net.add(Line("L1", 2, 3, R=0.01, X=0.1, wb=net.wb))
net.add(SynchronousMachine("SM1", 1, Sn=500e6, Sb=net.Sb, wb=net.wb))
net.add(RLLoad("Load2", 2, P=1.0, Q=0.2, wb=net.wb))    # pu base système
net.add(IB("IB3", 3, R=1e-4, X=1e-3, wb=net.wb))
net.run_loadflow()
ss = net.build_state_space()
```

Les éléments disponibles sont listés dans {mod}`g2egrid.elements` ; leurs
équations au chapitre {doc}`../manuel/03_modeles` et dans la page générée
{doc}`../equations_modeles`. Les paramètres des machines et des convertisseurs
ont des valeurs par défaut (`SM_DEFAULT`, `GFM_DEFAULT`, `GFL_DEFAULT`) ; ceux
qu'on donne au constructeur les remplacent.

Les options du réseau lui-même :

| Option | Valeurs | Effet |
|---|---|---|
| `Network(noeuds_sans_capacite=...)`, `Tables.NOEUDS` | `"ideal"` (défaut), `"resistance"` | ce que devient un bus à `C = 0` ({doc}`resolution_globale`) |
| `Network(r_causalite=...)`, `Tables.R_NOEUDS` | `100` pu (défaut, `Node.m`) | la résistance de causalité du mode `"resistance"` ({doc}`resolution_globale`) |
| `net.reference_flottante(...)` | `"coi"`, `None`, un nom, un dict | repère détaché du slack ({doc}`repere_flottant`) |
| `net.repere` | `"reference"` (défaut), `"nominale"` | pulsation vue par les éléments passifs |

## Par les tables

```python
from g2egrid import catalogue
from g2egrid.reseau.tables import assembler

t = catalogue.charger("wscc9")       # un objet Tables
net = assembler(t)
net.run_loadflow()
```

{class}`~g2egrid.reseau.tables.Tables` porte une liste de tuples par table, dans l'ordre
de colonnes de {data}`g2egrid.reseau.tables.TABLES`, plus des réglages :

| Champ | Contenu |
|---|---|
| `Y_DER` | bus : `BusID, Type, Un_V, V_pu, Theta_rad, Pg_W, Qg_var, Pl_W, Ql_var` |
| `T_Lines_SI`, `T_Lines_pu`, `T_Transfos`, `T_Shunts` | branches et shunts |
| `T_SMs`, `T_GFM1s`…`T_GFM4s`, `T_GFL1s`, `T_GFL2s`, `T_IBs`, `T_Loads` | équipements |
| `T_BusDC`, `T_LinesDC`, `T_ShuntsDC`, `T_VSCs` | le réseau continu |
| `SM_PARAMS`, `GFM_PARAMS`, `GFL_PARAMS` | réglages d'une famille entière |
| `PARAMS_PAR_BUS` | réglages propres à l'équipement d'un bus ; l'emportent sur la famille |
| `OUVERTS` | disjoncteurs ouverts, `(table, ID)` — voir {mod}`g2egrid.reseau.disjoncteurs` |
| `REFERENCE` | `"slack"` (défaut), `"coi"`, `"nominale"` — {doc}`repere_flottant` |
| `NOEUDS` | `"resistance"` (défaut), `"ideal"` — {doc}`resolution_globale` |

Une variante se fait par copie, sans toucher à l'original :

```python
t2 = t.copie(GFM_PARAMS={"mp": 0.02, "tau_VL": 20e-3},
             PARAMS_PAR_BUS={3: {"H": 6.0}},
             REFERENCE="coi")
```

La clé `dynamiques` d'un réglage choisit les dynamiques des équipements
concernés ({doc}`dynamiques`). Le détail des tables et le déroulé d'un notebook
sont au chapitre {doc}`../manuel/06_usage`.

## Point de fonctionnement

```python
net.run_loadflow()          # Newton-Raphson ; AC/DC unifié s'il y a du continu
net.solve_equilibrium()     # chaque élément : formules analytiques puis Newton
net.solve_equilibrium(pseudo_transitoire=True)   # recalage sur l'équilibre global
```

Le load flow attribue à chaque équipement la puissance qu'il injecte : les
charges connaissent la leur, le reliquat du bus revient à l'unique générateur
qui s'y trouve. Le détail est au chapitre {doc}`../manuel/02_principes`.
