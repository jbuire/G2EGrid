# Modes alternatifs, continus et couplés

Dans un réseau à liaisons continues, la première question qu'on pose à un
mode n'est pas sa classe (angle, tension, convertisseur...) mais son **côté** :
vit-il dans un réseau alternatif, dans le réseau continu, ou traverse-t-il les
convertisseurs ? {mod}`g2egrid.scenarios.domaines` répond à partir des participations.

## Le domaine d'un état, puis d'un mode

| États | Domaine |
|---|---|
| nœuds et câbles continus, DC/DC, charges continues | DC |
| boucle de tension continue d'un GFL à port continu (`M_d`) ; condensateur et courant continus internes (GFM3, GFM4, GFL2) | DC |
| filtre, PLL, boucles de courant et de puissance des convertisseurs ; machines, réseau alternatif, charges | AC |

Un mode est **AC** si sa part continue est sous 10 %, **DC** si sa part
alternative l'est, **AC/DC** (couplé) sinon. Chaque état appartient aussi à
une **île** — composante connexe du réseau alternatif ou continu, nommée par
son plus petit numéro de bus — et un mode que deux îles portent chacune à
10 % au moins est **inter-îles** : il traverse une liaison continue.

```python
from g2egrid import catalogue, domaines

net = catalogue.construire("cigre_b4")
ms = domaines.classer_modes(net, oscillants=True)
print(domaines.bilan(ms))
for m in [m for m in ms if m.couple][:4]:
    print(m.texte())
```

```text
{'AC': 179, 'DC': 9, 'AC/DC': 12, 'inter-iles': 18}
   -1.410     1.94 Hz  zeta  11.5 %  AC    (DC 0 %)  AC 1 88 %, AC 10 12 %, DC 201 0 %
   -1.440     1.92 Hz  zeta  11.9 %  AC    (DC 0 %)  AC 10 89 %, AC 1 11 %, DC 202 0 %
  -25.128     0.35 Hz  zeta  99.6 %  AC/DC (DC 45 %)  AC 10 54 %, DC 202 45 %, AC 1 1 %
  -30.705    82.17 Hz  zeta   5.9 %  AC/DC (DC 76 %)  DC 202 76 %, AC 1 12 %, AC 10 11 %
```

Sur le CIGRE B4, le mode électromécanique du réseau A (1,94 Hz) est porté à
12 % par le réseau B, qui ne lui est relié que par le continu : les stations
en statisme transmettent l'oscillation. Les résonances des câbles (40 à
214 Hz) sont des modes DC ; les modes AC/DC mêlent la boucle de tension
continue des stations et les nœuds continus.

## Ce que le continu fait aux modes alternatifs

{func}`~g2egrid.scenarios.domaines.impact_continu` compare le réseau à son double où le
continu est **rigide** — tensions des nœuds figées, câbles instantanés
({data}`~g2egrid.scenarios.domaines.CONTINU_RIGIDE`), avec les outils de
{doc}`simplification` :

```python
c = domaines.impact_continu(catalogue.charger("cigre_b4"))
print("\n".join(c.resume(n=0).splitlines()[:2]))
```

```text
348 → 322 etats ; max Re -0.1865 → -0.1863
retrouves 181, retires 12, perdus 7, nouveaux 3 ; ecart max 14.67 %, dzeta max 9.78 pts
```

Un mode **retrouvé mais déplacé** est un mode alternatif que la dynamique du
continu influence ; un mode **retiré** est un mode du continu lui-même.

## Dans l'atelier

Onglet **Modes**, pour un réseau à continu : le bilan par domaine s'affiche
avec la marge ; un filtre *Domaine* (AC, DC, AC/DC, couplés ou inter-îles) ;
sur la carte, un rond pour un mode alternatif, un carré pour un mode continu,
un losange pour un mode couplé ; dans `T_eig`, le domaine et la part continue
(`*` : inter-îles). La carte **Impact du continu sur les modes** liste les
modes que le continu déplace le plus, ceux du continu, et ceux qu'on perd en le
rendant rigide.
