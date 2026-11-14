# Interopérabilité : MATPOWER, PSS/E, pandapower

{mod}`g2egrid.formats.interop` lit les formats des autres outils **sans les installer** :
le `.m` de MATPOWER, le RAW de PSS/E et le JSON de pandapower sont lus tels
quels.

| Format | Lecture | Écriture |
|---|---|---|
| MATPOWER (`.m`, ou dict `mpc` de PYPOWER) | {func}`~g2egrid.formats.interop.depuis_matpower` | {func}`~g2egrid.formats.interop.vers_matpower` |
| PSS/E RAW, versions 33 et 34 | {func}`~g2egrid.formats.interop.depuis_psse` | |
| pandapower (JSON de `pp.to_json`, ou l'objet) | {func}`~g2egrid.formats.interop.depuis_pandapower` | |

{func}`~g2egrid.formats.interop.importer` choisit le lecteur d'après l'extension.

## Lire

```python
from g2egrid import interop

t = interop.depuis_matpower("tests/donnees/interop/case9.m")
print(len(t.Y_DER), "bus,", len(t.T_SMs), "machines,", len(t.T_Lines_pu), "lignes")

t_raw = interop.importer("tests/donnees/interop/case9_v34.raw")
t_pp = interop.depuis_pandapower("tests/donnees/interop/pp_mv_oberrhein.json")
print(len(t_pp.Y_DER), "bus,", len(t_pp.T_GFL1s), "générateurs statiques")
for a in t_pp.avertissements:
    print(" -", a)
```

```text
9 bus, 3 machines, 9 lignes
179 bus, 153 générateurs statiques
 - transformateur 114 : couplage 150 deg ignore (angles decales)
 - transformateur 142 : couplage 150 deg ignore (angles decales)
```

Le résultat est un objet {class}`~g2egrid.reseau.tables.Tables` ordinaire : il s'assemble,
se calcule, s'enregistre comme n'importe quel réseau du catalogue.

**Ce qui est repris** : bus, lignes, transformateurs à deux enroulements (prises
comprises), charges, shunts, unités et consignes. **La dynamique ne l'est pas** :
aucun de ces formats ne la porte (PSS/E la met dans un `.dyr`, dont les modèles
ne sont pas ceux de g2egrid). Les unités prennent les réglages par défaut du
catalogue ; `unites=` choisit leur table (machines synchrones par défaut), et
pour pandapower `sgen=` celle des générateurs statiques (GFL à puissance imposée
par défaut). Ce qui n'a pas d'équivalent (déphaseur, transformateur à trois
enroulements, magnétisation, couplage horaire) est ignoré **et signalé** dans
`t.avertissements`.

Particularités de pandapower : les interrupteurs bus-bus fermés fusionnent les
bus ; un interrupteur de ligne ouvert ne coupe qu'une extrémité, et la ligne à
vide devient l'admittance d'entrée de son schéma en π au bus qui l'alimente ;
l'impédance des transformateurs est rapportée au côté BT, comme pandapower le
fait ; `ext_grid` devient un bus infini.

## Les options qui comptent

| Option | Défaut | Effet |
|---|---|---|
| `charges` | `"Z"` | `"Z"` : une charge à $P \ge 0$, $Q > 0$ devient une R-L à impédance constante (convention du catalogue) ; `"P"` : toutes à puissance constante |
| `x_sur_r_max` | 100 | plancher $R \ge X/100$ : une branche sans résistance n'amortit rien ; `None` le retire |
| `x_min` | 0 | réactance minimale des branches (bretelles), pu base 100 MVA |
| `un_defaut` | 138 kV | tension de base des bus MATPOWER sans `BASE_KV` |

## Fidélité du load flow

Les tests (`tests/test_interop.py`) comparent le load flow importé à celui de
l'outil d'origine, calculé une fois et gardé dans `tests/donnees/interop/` :

| Cas | Écart sur $\lvert V\rvert$ |
|---|---|
| MATPOWER case9, case57 | $10^{-15}$ |
| PSS/E (case9 réécrit en v33 et v34) | $10^{-15}$ |
| pandapower `example_simple`, `mv_oberrhein` | $5\cdot 10^{-6}$ (la magnétisation des transformateurs, ignorée) |

Pour cette comparaison, `x_sur_r_max=None` et `t.NOEUDS = "ideal"` : la
résistance de causalité des bus sans capacité ({doc}`resolution_globale`) et le
plancher $X/R$ sont des choix de modélisation dynamique, pas des données.

```python
import numpy as np
from g2egrid.reseau.tables import assembler

t = interop.depuis_matpower("tests/donnees/interop/case57.m", x_sur_r_max=None)
t.NOEUDS = "ideal"
lf = assembler(t).run_loadflow()
print(f"{np.abs(lf.V).min():.4f} {np.abs(lf.V).max():.4f}")
```

```text
0.9359 1.0598
```

## Écrire

```python
from g2egrid import catalogue

texte = interop.vers_matpower(catalogue.charger("ieee39"), "ieee39_mpc.m")
retour = interop.depuis_matpower("ieee39_mpc.m", x_sur_r_max=None)
print(texte.splitlines()[0], "->", len(retour.Y_DER), "bus")
```

```text
function mpc = IEEE39 -> 39 bus
```

L'aller-retour redonne le même load flow (à $10^{-8}$). Le réseau continu et la
dynamique ne passent pas : l'export de l'atelier le dit dans ses avertissements.

## Dans l'atelier

Le bouton **Ouvrir** accepte `.m`, `.raw` et le JSON de pandapower (reconnu à sa
signature) ; les avertissements s'affichent dans la barre d'état. L'onglet
**Code** exporte un *Cas MATPOWER*.

Non couverts : CGMES (CIM) et l'export FMU.
