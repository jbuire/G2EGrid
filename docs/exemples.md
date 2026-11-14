# Exemples par fonctionnalité

Un exemple court, **exécutable tel quel**, pour chaque fonctionnalité de
g2egrid, de la construction d'un réseau au rapport d'étude. Chaque section
renvoie à la page du guide qui l'explique en détail. Les blocs se suivent : on
les exécute dans l'ordre, depuis le dossier `g2egrid/` du dépôt
(`uv run python`). La suite de tests les exécute tous
(`tests/test_docs_exemples.py`), si bien qu'un exemple qui ne tournerait plus
ferait échouer la suite.

```python
import numpy as np
from g2egrid import catalogue, modal
```

## Construire et charger

### Un réseau en code

```python
from g2egrid import Network
from g2egrid.elements import SynchronousMachine, Line, VoltageSource

net = Network(f0=50.0, Sb=1e6)
net.add_bus(1, 20e3, btype=2, V=1.0, Pg=0.8, C=1e-3)   # PV : la machine
net.add_bus(2, 20e3, btype=1, V=1.0, C=1e-3)           # slack : le réseau amont
net.add(SynchronousMachine("SM1", 1, Sn=1e6, Sb=1e6, wb=net.wb))
net.add(Line("L1", 1, 2, R=0.01, X=0.20, wb=net.wb))
net.add(VoltageSource("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
net.run_loadflow()
ss = net.build_state_space()
print(ss)
```

Voir {doc}`prise_en_main/premier_reseau` et {doc}`guide/construire`.

### Un réseau du catalogue, ses tables, une variante

```python
noms = [f.nom for f in catalogue.lister()]
t = catalogue.charger("wscc9")                    # les tables (Y_DER, T_SMs...)
t2 = t.copie(SM_PARAMS={"H": 6.0})                # une variante : toutes les inerties
net = catalogue.construire("wscc9")               # assemblé, load flow et équilibre faits
print(len(noms), "réseaux ;", len(t.Y_DER), "bus dans wscc9")
```

### Le point de fonctionnement

```python
lf = net.run_loadflow()
print(f"{lf.iterations} itérations, |V| de {np.abs(lf.V).min():.3f} à {np.abs(lf.V).max():.3f} pu")
```

## Analyse modale

### Modes, participations, forme modale

```python
ss = net.build_state_space()
for m in modal.modes(ss, f_max=5, oscillants=True)[:3]:
    print(f"{m.f:6.3f} Hz  zeta = {m.zeta:5.1f} %  {m.dominant}")
k = modal.modes(ss, f_max=5, oscillants=True)[0].indice
print(modal.forme_modale(ss, k, n=3))     # (état, module relatif, phase [deg])
```

Voir {doc}`guide/analyse_modale`.

### Le paramètre qui déplace un mode

```python
from g2egrid.analyse.sensibilites import sensibilites, tableau

m = modal.modes(ss, f_max=5, oscillants=True)[0]
print(tableau(sensibilites(net, m), n=3))
```

Voir {doc}`guide/sensibilites`.

### Résidus, commandabilité, observabilité : le rapport d'un mode critique

```python
from g2egrid import residus

for r in residus.rapport(ss, n_modes=1, f_min=0.5, f_max=3.0):
    print(r.elements[:2], r.entrees[:2], r.sorties[:2])
```

Voir {doc}`guide/residus`.

## Choisir les dynamiques

### Rendre des dynamiques algébriques (dx/dt = 0)

```python
from g2egrid.scenarios.simplification import comparer, preregle

tk = catalogue.charger("kundur2area")
c = comparer(tk, preregle(tk, "reseau_quasi_statique"))
print(c.resume().splitlines()[0])
```

Voir {doc}`guide/dynamiques` et {doc}`guide/simplification`.

### Nœuds idéaux et repère flottant

```python
from g2egrid.reseau.tables import assembler

net_ideal = assembler(t.copie(NOEUDS="ideal", REFERENCE="coi"))
net_ideal.run_loadflow()
ss_ideal = net_ideal.build_state_space()
print(len(ss_ideal.states), "états,", len(ss_ideal.elimines), "éliminés")
```

Voir {doc}`guide/resolution_globale` et {doc}`guide/repere_flottant`.

## Ce qui fait bouger les modes

### Puissance et N-1

```python
from g2egrid.scenarios.robustesse import etudier

e = etudier(catalogue.charger("wscc9"), f_max=3.0, zeta_max=30, deconnexion=False)
print(e.n_variantes, "variantes ; mode le plus sensible :", round(e.modes[0].f, 3), "Hz")
```

Voir {doc}`guide/robustesse`.

### Continuation d'un paramètre, bifurcation de Hopf

```python
from g2egrid import continuation as cont

net_gfl = catalogue.construire("smib_gfl")
suivi = cont.continuer(cont.par_parametre(net_gfl, "GFL2", "Kipll"), np.geomspace(10, 1000, 6))
print(suivi.resume(), "; max Re :", round(float(suivi.lam.real.max()), 3))
```

Voir {doc}`guide/continuation`.

## Convertisseurs et interactions

### Lois de synchronisation des GFM : statisme, dVOC, appariement

```python
tg = catalogue.charger("wscc9_gfm")
lignes = [r for k in ("T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s") for r in getattr(tg, k)]
tg = tg.copie(T_GFM1s=lignes, T_GFM2s=[], T_GFM3s=[], T_GFM4s=[],
              GFM_PARAMS=dict(tg.GFM_PARAMS, loi="dvoc"))
net_voc = assembler(tg)
net_voc.solve_equilibrium()
print([type(e).__name__ for e in net_voc.all_elements if e.name.startswith("GFM")],
      net_voc.build_state_space().is_stable())
```

Voir {doc}`guide/lois_gfm`.

### Critère d'impédance d'un convertisseur (boucle mineure)

```python
from g2egrid.analyse.impedance import boucle_mineure

b = boucle_mineure(catalogue.construire("smib_gfl"), "GFL2")
print(b.resume())
```

Voir {doc}`guide/impedance`.

## Grands réseaux

### Les modes critiques sans le spectre complet

```python
from g2egrid import descripteur as dd

d = dd.descripteur(catalogue.construire("ieee39"))
m = dd.balayage(d, f_max=3.0, profondeur=0.5)
print(f"J : {d.n} inconnues, densité {100 * d.densite:.2f} % ; {len(m)} modes dans la zone")
```

Voir {doc}`guide/grands_reseaux`.

### Transport et distribution

```python
from g2egrid.multi_niveau.multi_echelle import wscc9_cigre6, partition, comparer_partition

r = wscc9_cigre6()
c = comparer_partition(r, "transport", quasi_statique=True, part=partition(r))
print(f"écart des modes du transport réduit : {100 * c.ecart_max():.2f} %")
```

Voir {doc}`guide/multi_echelle`.

## Simuler dans le temps

```python
from g2egrid.temporel import ModeleNonLineaire

mnl = ModeleNonLineaire(catalogue.construire("smib"))
tr = mnl.simuler((0.0, 2.0), uext_fn=mnl.echelon("SM1.p_ref", 0.05, t0=0.1),
                 t_eval=np.linspace(0, 2, 201))
print(f"vitesse rotor en fin de simulation : {tr.etat('SM1.Dwr')[-1]:.2e} pu")
```

Voir {doc}`guide/analyse_modale`, section *Simulation non linéaire*.

## Échanger avec d'autres outils

```python
from g2egrid import interop

texte = interop.vers_matpower(catalogue.charger("wscc9"))   # un cas MATPOWER
t_mp = interop.depuis_matpower(texte, nom="wscc9_mp")       # et retour
t_raw = interop.depuis_psse("tests/donnees/interop/case9_v33.raw")
t_pp = interop.depuis_pandapower("tests/donnees/interop/pp_simple.json")
print(len(t_mp.Y_DER), len(t_raw.Y_DER), len(t_pp.Y_DER), t_pp.avertissements[:1])
```

Voir {doc}`guide/interoperabilite`.

## Exporter

```python
from g2egrid.formats.export import exporter, vers_matlab

cellule = exporter(t)                    # la cellule de données d'un notebook
source_m, perdu = vers_matlab(t)         # un script MATLAB G2Elin V2, et ce qu'il perd
print(len(cellule.splitlines()), "lignes Python ;", len(source_m.splitlines()), "lignes .m")
```

## Automatiser : l'atelier sans l'interface

Tout ce que fait l'atelier se fait par les mêmes fonctions, en Python :

```python
from g2egrid_api import analyse

print(analyse.valider(t)[:1])                                     # rien : réseau valide
print(analyse.loadflow(t)["iterations"], "itérations")
co = analyse.continuation_parametre(t, "SM1", "H", 2.0, 6.0, pas=5)   # l'onglet Limites
print(len(co["mu"]), "pas ;", co["bifurcations"])
```

Le suivi d'un calcul long, et son arrêt, passent par {mod}`g2egrid.noyau.progres` :

```python
from g2egrid import progres

avancement = []
with progres.suivi(lambda f, msg: avancement.append(round(f, 2))):
    etudier(catalogue.charger("smib"), deconnexion=False)
print(avancement[:4])
```

Le rapport d'étude autonome (version du code, empreinte des tables) :

```python
from g2egrid_api.rapport import empreinte, rapport_etude

page = rapport_etude(t, analyse.construire)
print(empreinte(t)[:16], len(page) > 10000)
```

Voir {doc}`guide/automatiser`.
