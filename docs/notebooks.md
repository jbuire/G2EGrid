# Notebooks

Les notebooks du dépôt, rendus **avec leurs sorties enregistrées** (ils ne sont
pas réexécutés à la construction de la documentation). Chaque réseau de test a
le sien, au format de `cigre_mv.ipynb` : tables de données, assemblage, load
flow, modèle d'état, table des modes et carte des valeurs propres. Les études
reprennent les chapitres {doc}`manuel/08_sources`,
{doc}`manuel/09_transport_distribution` et {doc}`manuel/10_transition`.

Pour les exécuter : `uv run --extra plot jupyter lab` depuis `g2egrid/`.

```{toctree}
:maxdepth: 1
:caption: Réseaux de test
:glob:

_notebooks/reseaux_test/*
```

```{toctree}
:maxdepth: 1
:caption: Études
:glob:

_notebooks/etudes/*
```
