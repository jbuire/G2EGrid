# Équations dynamiques modifiables

Les équations `dx/dt = f(x, u)` de chaque équipement sont lisibles et
modifiables, élément par élément ou dans une **liste globale** du réseau
(onglet Équations de l'atelier ; tout ce qu'il fait se fait aussi en code).

## Le modèle à équations

Un élément s'écrit en texte, nom par nom : états (`dx/dt = ...`), variables
algébriques (`z = ...`, ou `0 = résidu` si elles ne s'explicitent pas),
consignes, paramètres, sorties (`id`, `iq`). Les expressions passent par
l'analyseur sûr des schémas-blocs (arithmétique et fonctions usuelles, rien
d'autre n'est jamais évalué). Les entrées du réseau s'appellent `vd`, `vq`,
`V`, `w`, `theta_g` ; `wb` et `kS` sont fournis par le réseau.

```python
from g2egrid.elements.blocs import compiler, equations_depuis_element

m = equations_depuis_element("GFL1")
assert m["genre"] == "equations"
print([e["nom"] for e in m["etats"]][:4])
print(m["etats"][-1]["nom"], "=", m["etats"][-1]["derivee"])
assert len(compiler(m).etats) == 12
```

Les variables algébriques natives restent des signaux nommés (`w_pll`, `md`,
`mq` pour un GFL1) : on ne les développe pas dans les dérivées.

## La liste globale

```python
import g2egrid.catalogue as cat
from g2egrid.elements.blocs import equations_reseau

t = cat.charger("wscc9_gfm")
lignes = equations_reseau(t)
gfm = [x for x in lignes if x["element"] == "GFM2" and x["genre"] == "etat"]
print(len(lignes), "equations ;", len(gfm), "etats pour GFM2")
ligne = next(x for x in lignes if x["classe"] == "Line")
print(ligne["modifiable"], "-", ligne["pourquoi"])
```

Une ligne par équation de chaque élément assemblé, branches et nœuds compris.
Ces derniers ne se modifient pas ; `pourquoi` donne la raison.

## Modifier une équation

```python
import copy
from g2egrid.elements.blocs import modifier_equation, revenir_au_natif

t = copy.deepcopy(cat.charger("smib_gfl"))
t.T_GFL1s, t.T_GFL2s = t.T_GFL2s, []
# on double Kipll dans dM_pll/dt de l'element GFL2 (GFL1 natif, bus 2)
t2 = modifier_equation(t, "GFL2", "M_pll", "2*Kipll*veq")
assert t2.T_GFL1s == [] and len(t2.T_Blocs) == 1   # remplace par son modele
# les tables recues n'ont pas bouge ; une erreur les laisse intactes
t3 = revenir_au_natif(t2, ("T_Blocs", 0))
assert t3.T_GFL1s == t.T_GFL1s and t3.T_Blocs == []
```

L'élément natif est d'abord remplacé par son modèle à équations (avec ses
paramètres), rangé dans `MODELES_BLOCS` sous `<classe>_<élément>` et posé dans
`T_Blocs` ; l'équation y est ensuite remplacée. Le spectre obtenu est celui du
natif avec `Kipll` doublé (écart de l'ordre de 1e-14). Autres opérations
(`g2egrid.reseau.equations_reseau`) :

- `ajouter_equation` (nouvel état ou nouvelle variable, avec sa valeur de
  départ), `supprimer_equation` (si plus personne ne l'emploie) ;
- `modifier_parametre` : la valeur d'un paramètre du modèle. Un schéma-bloc
  **reste** un schéma-bloc, seul son paramètre change ;
- `remplacer_par_modele(t, element, "equations" | "blocs")` : le remplacement
  seul, sans rien modifier ;
- `revenir_au_natif` (ou `retour_au_natif`, qui rend aussi l'emplacement de la
  ligne rendue).

Le modèle garde son **origine** (la ligne native, et pour une charge la
puissance reportée à son bus) : le retour au natif rend les tables de départ,
à l'identique, y compris après avoir enregistré puis rouvert le réseau. Si le
bus a été retouché entre-temps, seules les puissances déplacées reviennent ;
si l'identifiant de la ligne a été repris, elle en reçoit un neuf.

Modifier une équation d'un élément décrit par un **schéma-bloc** n'efface pas
le dessin : l'élément reçoit un modèle à équations à part (`<nom>_eq`), le
schéma reste dans `MODELES_BLOCS`.

Une modification que la compilation refuse (nom inconnu, valeur non finie,
paramètre nul au dénominateur d'une équation...) laisse les tables intactes
et dit pourquoi.

Limites : un modèle à équations décrit le fonctionnement **libre** (les
drapeaux de limiteur deviennent des paramètres à 0) ; les éléments à port
continu, les branches et les nœuds ne se remplacent pas ; un élément cité par
`SIMPLIFICATIONS` ou ouvert par un disjoncteur non plus.

## Vers les blocs, et retour

`schema_depuis_equations(modele)` donne le schéma-bloc (onglet Blocs) ;
`equations_depuis_schema(schema)` fait l'inverse en substituant les sorties
intermédiaires de blocs sans état tant que l'expression reste lisible.

## L'API

| route | corps | réponse |
|---|---|---|
| `POST /api/equations/reseau` | `{tables}` | `{lignes}` |
| `POST /api/equations/element` | `{tables, table, index}` | `{modele}` ; pour une ligne de `T_Blocs`, en plus `natif` (l'élément natif remplacé, ou `null`) et `schema_bloc` |
| `POST /api/equations/modifier` | `{tables, table, index, variable, expression}` | `{tables, modele, table, index}` |
| `POST /api/equations/ajouter` | `{..., genre, variable, expression, depart?}` | `{tables, modele, table, index}` |
| `POST /api/equations/supprimer` | `{..., variable}` | `{tables, modele, table, index}` |
| `POST /api/equations/parametre` | `{..., nom, valeur}` | `{tables, modele, table, index}` ; un élément natif est refusé (ses paramètres se règlent par `PARAMS_PAR_BUS`) |
| `POST /api/equations/revenir` | `{tables, table, index}` | `{tables, table, index}` |
| `POST /api/equations/verifier` | `{modele}` | `{ok, message?, equations}` |
| `POST /api/equations/depuis_schema` | `{schema_bloc}` | `{modele}` |
| `POST /api/equations/vers_blocs` | `{modele}` | `{schema}` |

`table` et `index` de la réponse donnent l'emplacement de l'élément dans les
nouvelles tables (un natif modifié est passé dans `T_Blocs`, un élément rendu
à son modèle natif a retrouvé sa table). Dans l'atelier, la section
**Équations** d'un équipement liste aussi les **paramètres** du modèle :
modifiables sur place pour un élément converti, en lecture seule pour un natif.
