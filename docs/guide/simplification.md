# Hypothèses de simplification du calcul des modes

Les modèles de g2egrid sont complets : courants des lignes, flux statoriques,
filtres et boucles des convertisseurs... Une étude de stabilité se fait
pourtant souvent sur un modèle **réduit**, où les dynamiques rapides sont
tenues pour instantanées. La question est alors : *qu'a-t-on perdu ?*
{mod}`g2egrid.scenarios.simplification` y répond en confrontant les modes du modèle
complet à ceux du modèle simplifié.

## Déclarer les hypothèses

Le champ `SIMPLIFICATIONS` des tables dit, **par type d'équipement** ou
**par équipement**, le sort de chaque état :

```python
from g2egrid.reseau.catalogue import charger

t = charger("kundur2area")
t = t.copie(SIMPLIFICATIONS={
    # tous les courants de ligne : dx/dt = 0 (régime sinusoïdal établi)
    "Line": {"i_d": "statique", "i_q": "statique"},
    # les flux statoriques de toutes les machines...
    "SynchronousMachine": {"phi_d": "statique", "phi_q": "statique"},
    # ... sauf ceux de SM2, qui restent dynamiques
    "SM2": {"phi_d": "dynamique", "phi_q": "dynamique"},
})
```

- `"statique"` : $\dot x = 0$, l'état devient une variable algébrique
  (voir {doc}`dynamiques`) ;
- `"fige"` : $x = x_0$, l'état garde sa valeur d'équilibre ;
- `"dynamique"` : rien — c'est ce qui permet à un équipement de rendre une
  dynamique que son type supprime.

L'ordre de priorité, état par état : le **type**, puis la famille et le bus
(`..._PARAMS["dynamiques"]`, `PARAMS_PAR_BUS`), puis l'**équipement**. Un nom
qui n'est ni un équipement ni un type du réseau est une erreur. Le champ
voyage avec le fichier `.g2egrid.json`.

## Les préréglages

{data}`g2egrid.scenarios.simplification.PREREGLES` rassemble les hypothèses classiques ;
{func}`~g2egrid.scenarios.simplification.preregle` les restreint aux types présents :

| Préréglage | Hypothèse |
|---|---|
| `reseau_quasi_statique` | réseau en phaseurs : courants des branches et des charges, tensions des bus algébriques |
| `flux_stator` | flux statoriques des machines algébriques (modèle « RMS ») |
| `amortisseurs` | plus de subtransitoire |
| `excitation` | mesure et excitatrice instantanées, retour de stabilisation supprimé |
| `pss` | stabilisateurs instantanés |
| `turbine` | turbine et régulateur de vitesse instantanés |
| `filtres_convertisseurs` | filtres LC des convertisseurs algébriques |
| `boucles_courant` | boucles de courant instantanées (intégrateur **et** inductance côté convertisseur) |
| `boucles_tension_gfm` | boucles de tension et de courant des GFM instantanées |
| `pll` | correcteur PI des PLL instantané |
| `mesures_puissance` | filtres de mesure de puissance instantanés |
| `bus_dc` | bus continu parfaitement régulé |

Toutes les combinaisons ne sont pas des modèles : rendre algébrique
l'intégrateur d'une boucle de courant *sans* l'inductance qu'elle commande
impose le courant, état dynamique, par une contrainte — le système devient
d'indice 2 et l'assemblage le refuse. Les préréglages associent ce qui doit
l'être.

## Ce qu'on perd

```python
from g2egrid.scenarios.simplification import comparer, preregle

c = comparer(t, preregle(t, "reseau_quasi_statique", "flux_stator"))
print(c.resume())
```

{func}`~g2egrid.scenarios.simplification.comparer` apparie les modes des deux modèles
(par leurs participations sur les états communs, puis leurs valeurs propres)
et range chaque mode du modèle complet :

retrouvé
: le modèle simplifié l'a ; on lit son déplacement $|\Delta\lambda|/|\lambda|$
  et $\Delta\zeta$ ;

retiré
: porté à plus de 50 % par les états supprimés. C'est le prix **attendu** de
  l'hypothèse ; il ne coûte rien tant que ce mode est bien amorti ;

perdu
: porté par les états gardés, mais introuvable : l'hypothèse a changé plus
  que ce qu'elle supprime ;

nouveau
: un mode du modèle simplifié sans équivalent dans le complet.

Une simplification est **sûre** pour la stabilité
({meth}`~g2egrid.scenarios.simplification.ComparaisonSimplification.sure`) quand elle ne
retire ni ne perd aucun mode instable ou mal amorti ($\zeta < 5\,\%$ sous
100 Hz), n'en invente aucun, et ne change le signe d'aucune partie réelle.
Les hypothèses classiques ne le sont pas toujours : sur WSCC9, rendre seul le
retour tachymétrique de l'AVR ($e_3$) algébrique rend le réseau instable ; sur
le SMIB à GFL (`smib_gfl`), des boucles de courant instantanées créent une
instabilité à haute fréquence que le modèle complet n'a pas.

## Dans l'atelier

L'onglet **Simplifications** liste les préréglages applicables au réseau,
puis chaque type d'équipement avec un choix *dynamique / dx/dt = 0 / figé*
par état ; déplier un type pour régler ses équipements un par un.
**Comparer les modes** trace les deux spectres dans le plan complexe et
dresse le tableau des modes, du moins amorti au plus amorti, avec les
alertes. La case **appliquer à tous les calculs** fait prendre le modèle
simplifié par les autres onglets (Modes, Variations, Robustesse,
Temporel...) ; un réseau enregistré la garde.
