# Le CIGRE B4 DC Grid Test System

Le réseau de test continu de référence (CIGRE WG B4.57/B4.58, Vrana *et al.*,
Electra 270, 2013) est au catalogue sous trois formes :

| Nom | Contenu |
|---|---|
| `cigre_b4` | le système complet : réseaux alternatifs A et B, îles offshore C, D, E, F, onze stations AC/DC, DCS1, DCS2 et DCS3 |
| `cigre_b4_mmc` | le même, stations terrestres en MMC moyen à énergie libre (voir {doc}`mmc`) |
| `cigre_b4_dcs1` | le cas d'étude « DCS1 » : la liaison point à point A1–C1 entre deux réseaux extérieurs |

Les trois sont construits par `etudes/cigre_b4/construire.py`, à partir des
tableaux publiés avec le modèle de load flow DIgSILENT PowerFactory du réseau.

## Le load flow, contre la référence

```python
from g2egrid import catalogue
from g2egrid.reseau.tables import assembler

net = assembler(catalogue.charger("cigre_b4"))
lf = net.run_loadflow()
v = dict(zip(sorted(net.bus_specs_dc), lf.Vdc))
print({b: round(float(v[b]), 4) for b in (202, 211, 213, 241, 251)})
```

```text
{202: 1.01, 211: 1.0005, 213: 0.9868, 241: 1.0067, 251: 1.0069}
```

Les treize tensions continues tombent à moins de $10^{-3}$ pu de celles de
PowerFactory (1,01 ; 1,0002 ; 0,9861 ; 1,0065 ; 1,0068 pour les mêmes nœuds :
Bb-A1, Bb-B1, Bb-B2, E1, Bm-F1). Sur DCS1, Bm-C1 est à 1,0108 pu des deux
côtés. Les écarts restants viennent des pertes des stations, qui sont celles
des modèles de convertisseur de g2egrid et non celles du modèle de référence
(`tests/test_cigre_b4.py`).

## Ce que le réseau a fait entrer dans le moteur

Deux modes de convertisseur (`T_VSCs.Mode`) :

| Mode | La station tient | Pour |
|---|---|---|
| `Vtheta` | le module et l'angle de son bus alternatif : elle **forme** son île, elle en est le slack ; sa puissance est celle que l'île appelle | les stations offshore Cm-C1, Cb-D1, Cm-E1, Cm-F1 (GFM) |
| `droop` | une relation : $p_{ac} = p_{ref} + \dfrac{S_n}{S_b\,k_{droop}}\,(v_{dc} - v_{dc,ref})$ | les statismes Cb-B1, Cb-B2, Cm-B3 (GFL) |

En `droop`, $p_{ref}$ est le `Pg_W` du bus de la station, $v_{dc,ref}$ le
`V_pu` de son nœud continu (de type 2) et `kdroop` (pu de tension continue
par pu de puissance, base $S_n$) se règle dans `PARAMS_PAR_BUS`. Le gain de
la boucle continue du GFL suit le même statisme, $K_{p,dc} = 1/(k_{droop}\,v_{e,d})$ :
load flow et modèle dynamique ont la même pente.

Le partage de puissance dépend des statismes. {func}`g2egrid.scenarios.continuation.par_statisme`
les multiplie tous (ou ceux de quelques stations) par $\mu$ :

```python
from g2egrid import continuation as cont

fab = cont.par_statisme(catalogue.charger("cigre_b4"))
c = cont.continuer(fab, [0.5, 1.0, 2.0])
for mu in (0.5, 1.0, 2.0):
    p = fab.mesures[mu]["p_mw"]
    print(mu, {b: round(p[b]) for b in (111, 113, 114, 102, 112)})
print(c.resume())
```

```text
0.5 {111: 1520, 113: 1421, 114: 717, 102: -1511, 112: -251}
1.0 {111: 1505, 113: 1541, 114: 754, 102: -1632, 112: -281}
2.0 {111: 1501, 113: 1615, 114: 775, 102: -1711, 112: -299}
aucune bifurcation
```

Plus raides ($\mu < 1$), les trois stations en statisme (bus 111, 113, 114)
s'écartent davantage de leur consigne et soulagent les deux références de
tension (Cb-A1 au bus 102, Cm-B2 au bus 112). Dans l'atelier : bloc
*Statisme Vdc-P* des propriétés d'une station, et carte *Statismes Vdc-P* de
l'onglet *Limites et impédances*.

Une île continue peut aussi avoir **plusieurs références de tension**, à des
nœuds distincts — Cb-A1 tient 1,01 pu et Cm-B2 1,00 pu dans le même réseau,
comme deux bus PV du côté alternatif. Deux références au même nœud restent
refusées.

## Ce qui est représenté autrement

- **Bipôles et monopôles symétriques** : un conducteur équivalent à la tension
  pôle à pôle (boucle $2r$, $2l$, capacité entre pôles $c/2$) ; les deux pôles
  d'un bipôle forment une station de puissance double. Les lignes gardent
  leurs données par pôle et déclarent leur configuration (`Config`,
  `Circuits` de `T_LinesDC`, voir {doc}`mmc`).
- **DC/DC idéaux** : Cd-B1 (rapport 1) confond Bb-B1x et Bb-B1 ; Cd-E1
  (rapport 0,5, exactement celui des bases) vaut 1 en pu. DCS2 est ramené à la
  base 800 kV et fusionné avec DCS3 au nœud E1.
- **Stations** : VSC à filtre LCL de g2egrid (le modèle moyen de MMC est dans
  la variante `cigre_b4_mmc`), réglés au plus près d'un MMC : réactance de phase 0,15 pu, condensateur
  de filtre réduit, boucle de courant à 2 ms, sans anticipation de tension.
  Avec le filtre LCL et l'anticipation par défaut, DCS1 est instable dès le
  départ : la résonance du filtre avec le transformateur et le réseau C part à
  230 Hz.
- **Productions** : machines synchrones à terre, GFL à puissance imposée en
  mer, bus infinis de 30 GVA pour les réseaux extérieurs.

## Petit signal

```python
from g2egrid import modal

net = catalogue.construire("cigre_b4")
ss = net.build_state_space()
print(len(ss.states), "états ; stable :", ss.is_stable())
for m in modal.modes(ss, f_max=1000, oscillants=True)[:2]:
    print(f"{m.f:7.1f} Hz  zeta = {m.zeta:5.1f} %  {m.dominant}")
```

```text
348 états ; stable : True
  780.0 Hz  zeta =   0.9 %  bus10.vq
  680.4 Hz  zeta =   1.2 %  bus10.vq
```

Les modes les moins amortis sont des résonances des lignes 380 kV du réseau B
avec l'impédance du réseau extérieur. Une île offshore n'a pas de référence
d'angle absolue : chacune apporte un mode de jauge, écarté du verdict.
