# Topologie continue et MMC moyen

Deux ajouts au format, tous deux optionnels : un fichier qui ne les emploie
pas garde exactement le sens qu'il avait.

## La topologie d'un câble continu

Un réseau continu de g2egrid est décrit par un conducteur ÉQUIVALENT, à la
tension pôle à pôle de ses nœuds. Les fiches de câble, elles, donnent des
paramètres PAR PÔLE. Deux colonnes optionnelles de `T_LinesDC`, en fin de
ligne, font le passage :

| Colonne | Valeurs | Défaut |
|---|---|---|
| `Config` | `equivalent`, `monopole_sym`, `bipole`, `monopole_metallique`, `monopole_terre` | `equivalent` |
| `Circuits` | nombre de circuits identiques en parallèle | 1 |

| Configuration | R et L | C |
|---|---|---|
| `equivalent` | telles quelles | telle quelle |
| `monopole_sym`, `bipole` | 2r, 2l (la boucle traverse les deux pôles) | c/2 (deux capacités pôle-terre en série) |
| `monopole_metallique` | 2r, 2l (pôle et retour) | c (seul le pôle est sous tension) |
| `monopole_terre` | r, l | c |

Un bipôle est supposé équilibré, son retour métallique sans courant, et ses
deux demi-stations agrégées : son équivalent est celui du monopôle symétrique.
`n` circuits divisent R et L par `n` et multiplient C par `n`.

```python
from g2egrid import catalogue
from g2egrid.reseau.tables import cable_dc

t = catalogue.charger("cigre_b4")
print(t.T_LinesDC[0])
print(cable_dc(t.T_LinesDC[0]))
```

```text
(1, 201, 221, 200.0, 0.011, 2.615, 0.2185, 'monopole_sym', 1)
(1, 201, 221, 200.0, 0.022, 5.23, 0.10925)
```

Le CIGRE B4 du catalogue décrit maintenant ses lignes par pôle (tableau 6
de la brochure), avec leur configuration : c'est l'assembleur qui en tire
l'équivalent. Dans l'atelier, l'inspecteur d'un câble continu propose les
deux champs.

## Le MMC moyen

Une station de `T_GFL2s` raccordée par `T_VSCs` est par défaut un VSC à
deux niveaux, avec son filtre LCL. Le réglage `modele` (dans `GFL_PARAMS`
pour toutes les stations, dans `PARAMS_PAR_BUS` pour une seule) en fait un
convertisseur modulaire multiniveaux, représenté par son modèle MOYEN :

| `modele` | Classe | États | Ce qui fixe la tension continue des bras |
|---|---|---|---|
| `vsc` (défaut) | `GFL2DC` | 12 | (pas de bras : un condensateur continu) |
| `mmc` | `MMC` | 11 | l'énergie des sous-modules, laissée libre |
| `mmc_energie` | `MMCEnergie` | 13 | une boucle de courant continu, sous une boucle d'énergie |

Le MMC n'a pas de bus continu : il a l'ÉNERGIE `W` stockée dans ses
sous-modules, que les deux côtés alimentent (`dW/dt = v_z.idc - p_ac`).
Côté alternatif, un courant à travers `L1 = Lc + Larm/2` ; côté continu, un
courant à travers `Ldc = 2.Larm/3`. Ses modes de station sont ceux du VSC
(`VdcQ`, `PQ`, `droop`), et sa boucle externe aussi. Les équations sont
dans {doc}`../equations_modeles`.

```python
from g2egrid.reseau.tables import assembler

d = catalogue.charger("cigre_b4_dcs1")
for modele in ("vsc", "mmc", "mmc_energie"):
    net = assembler(d.copie(GFL_PARAMS=dict(d.GFL_PARAMS, modele=modele)))
    print(modele, [(e.name, type(e).__name__, len(e.states))
                   for e in net.elements if e.ports and e.ports_dc])
```

```text
vsc [('GFL101', 'GFL2DC', 12), ('GFL121', 'GFL2DC', 12)]
mmc [('GFL101', 'MMC', 11), ('GFL121', 'MMC', 11)]
mmc_energie [('GFL101', 'MMCEnergie', 13), ('GFL121', 'MMCEnergie', 13)]
```

L'initialisation est analytique et exacte : le load flow mixte prend les
pertes de la station au modèle (`R1.|i|² + Rdc.idc²`), et le point de
fonctionnement est un point fixe du modèle non linéaire à 1e-7 près.

### Énergie libre ou régulée

**Libre** (`mmc`, le modèle moyen de référence des brochures B4.57/B4.58) :
la tension continue des bras est la somme des tensions des condensateurs,
`v_z = sqrt(W/W0)`. Vu du réseau continu, la station est une inductance en
série avec la capacité équivalente de ses sous-modules (`2.W0` secondes) :
c'est elle qui tient la tension du nœud.

**Régulée** (`mmc_energie`) : une boucle d'énergie fixe la consigne de
courant continu, une boucle de courant continu fabrique `v_z`. L'énergie
est découplée de la tension continue — la station ne montre plus sa
capacité au réseau, dont les nœuds ne gardent que celle des câbles.

```python
from g2egrid import domaines

for modele in ("mmc", "mmc_energie"):
    net = assembler(d.copie(GFL_PARAMS=dict(d.GFL_PARAMS, modele=modele)))
    net.run_loadflow()
    ms = domaines.classer_modes(net, oscillants=True)
    print(modele)
    for m in sorted(ms, key=lambda m: m.zeta)[:2]:
        print(m.texte())
```

```text
mmc
   -9.470   189.37 Hz  zeta   0.8 %  DC    (DC 96 %)  DC 201 96 %, AC 2 2 %, AC 21 2 %
  -16.126   200.31 Hz  zeta   1.3 %  DC    (DC 96 %)  DC 201 96 %, AC 21 3 %, AC 2 1 %
mmc_energie
   -35.129    73.88 Hz  zeta   7.5 %  DC    (DC 94 %)  DC 201 94 %, AC 2 5 %, AC 21 0 %
 -464.649   174.62 Hz  zeta  39.0 %  AC    (DC 4 %)  AC 21 96 %, DC 201 4 %, AC 2 0 %
```

Sur DCS1, l'énergie libre laisse une résonance du câble (190 à 200 Hz) très
peu amortie : la capacité des sous-modules, derrière `Ldc`, en fait un
circuit LC que seules les résistances amortissent. Régulée, l'énergie
supprime ce circuit. Mais sur le système complet, c'est l'inverse : les
statismes raides (Cb-B1, Cb-B2, Cm-B3, `K = 0,2 pu`) n'agissent plus que sur
la capacité des lignes, et leur boucle part instable (`tests/test_mmc.py`).
C'est le résultat connu des réseaux MTDC à MMC régulés en énergie ; la
variante du catalogue `cigre_b4_mmc` laisse donc l'énergie libre, et elle est
stable.

### Deux choix de réglage, et pourquoi

* **La mesure de tension continue est filtrée** (`T_vdc`, 5 ms). Sans filtre,
  le proportionnel de la boucle de tension agit encore à la fréquence de
  résonance du nœud ; à travers l'intégrateur d'énergie, il y donne à la
  station une conductance négative, et DCS1 part instable vers 200 Hz.
* **Pas d'anticipation de tension alternative** (`Kffv = 0`). À 1, le
  système complet est instable vers 195 Hz.

### Défauts et blocage

Le blocage (manœuvre `blocage`, ou plan de {doc}`protections_dc`) coupe la
tension qui entraîne le courant continu : les bras bloqués l'éteignent en
1 ms environ, le courant restant continu au nœud, et l'énergie des
sous-modules ne bouge plus. Les diodes de roue libre, qui empêcheraient la
tension continue de s'inverser, ne sont pas représentées.

```python
from g2egrid.temporel import Manoeuvre, ModeleNonLineaire, protections_dc

net = assembler(d.copie(GFL_PARAMS=dict(d.GFL_PARAMS, modele="mmc")))
net.run_loadflow()
m = ModeleNonLineaire(net)
sc = m.scenario(Manoeuvre(0.02, "defaut_dc", "bus221"))
tr = m.simuler((0, 0.05), uext_fn=sc, methode="BDF2", pas=2.5e-5,
               protections=protections_dc(net))
for t, nom, actions in tr.declenchements:
    print(f"{1e3 * t:6.2f} ms  {nom:32s} {actions}")
```

```text
 20.16 ms  sous_tension_dc:bus221           [('blocage', 'GFL121')]
 31.38 ms  surintensite_dc:CDC1_201_221     [('ouverture', 'CDC1_201_221')]
 32.11 ms  sous_tension_dc:bus201           [('blocage', 'GFL101')]
```

Le MMC en énergie libre tient la tension de Bm-A1 bien plus longtemps que
le VSC (bloqué à 22,8 ms) : sa capacité équivalente est cent fois celle du
condensateur continu d'un VSC. C'est donc la surintensité du câble qui
déclenche d'abord, puis la station.
