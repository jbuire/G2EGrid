# G2Egrid

[![Licence : GPL v3](https://img.shields.io/badge/Licence-GPLv3-blue.svg)](LICENSE)

*[English version](README.md)*

**G2Egrid** est un outil Python d'analyse de stabiltié de réseaux électriques. Chaque équipement
déclare un modèle DAE symbolique ; le réseau assemble le tout en un modèle
d'état global (A, B, C, D), dont on tire analyse modale, sensibilités,
impédances, continuation et simulation temporelle.

Il est développé au G2Elab (Univ. Grenoble Alpes, CNRS, Grenoble INP) et
succède à l'outil MATLAB G2ELIN.

## Ce qu'il fait

- **Modèles** : machines synchrones (avec AVR, PSS, régulateur de vitesse,
  lignes d'arbre multi-masses), convertisseurs grid-forming (statisme,
  oscillateur virtuel, appariement), convertisseurs grid-following, stations
  MMC et VSC, réseaux continus, lignes, transformateurs, charges ; équipements
  de l'utilisateur en **schémas-blocs** ou en **équations modifiables**
  `dx/dt = f(x, u)`.
- **Point de fonctionnement** : load flow AC et AC/DC, limites des
  convertisseurs.
- **Analyses** : modes, facteurs de participation et formes modales,
  sensibilités, résidus, impédances dq et critère de Nyquist généralisé,
  puissance de court-circuit et SCR, continuation d'un paramètre et
  bifurcations de Hopf, robustesse aux variations de puissance et au N-1,
  simplification de modèle, interaction transport / distribution.
- **Temporel** : simulation linéaire et non linéaire d'événements et de
  scénarios.
- **Formats** : lecture MATPOWER, PSS/E RAW (v33, v34) et pandapower ;
  export MATPOWER, MATLAB et notebook.
- **Un catalogue de 34 réseaux** : SMIB, WSCC 9 nœuds, Kundur deux zones,
  IEEE 14, 33, 39 et 118, CIGRE MT, réseau continu CIGRE B4, Nordic 44, et
  leurs variantes.
- **Un atelier web** qui tourne en local : construire un réseau sur un schéma
  unifilaire ou dans des tables, lancer toutes les analyses, éditer les
  schémas-blocs et les équations. Son interface est disponible **en français
  et en anglais** (bouton FR / EN).

## Démarrage

Il faut [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/jbuire/g2elin_v2.git
cd g2elin_v2/g2egrid
uv sync --all-extras --all-groups
uv run pytest -q --ignore=tests/e2e        # la suite de tests (un quart d'heure)
```

G2Egrid s'utilise depuis un clone du dépôt : le catalogue de réseaux
(`reseaux/`, `reseaux_test/`) est à côté du paquet, pas dedans. Pour une
installation non éditable, la variable `G2EGRID_DONNEES` désigne le dossier
d'un clone.

L'atelier (un serveur local et une page de navigateur, sur votre machine
seulement) :

```bash
./lancer_atelier.sh            # Linux, macOS ; lancer_atelier.bat sous Windows
# ou : uv run --extra api python -m g2egrid_api
```

Un premier modèle en quelques lignes :

```python
from g2egrid.reseau.catalogue import construire

net = construire("wscc9")              # load flow et equilibre compris
ss = net.build_state_space()           # (A, B, C, D) global
print(ss.eig()[:5], ss.is_stable())
```

## Organisation

| Chemin | Contenu |
|---|---|
| `g2egrid/` | le paquet : `noyau/` (DAE, compilation symbolique, load flow), `elements/` (bibliothèque d'équipements, schémas-blocs), `reseau/` (`Network`, tables, catalogue, zones), `temporel.py`, `analyse/`, `scenarios/`, `multi_niveau/`, `formats/`, `validation/` |
| `g2egrid_api/` | l'atelier : serveur FastAPI (`app.py`), calculs sans HTTP (`analyse/`), interface web (`web/`) |
| `tests/` | la suite pytest ; `tests/e2e/` pilote l'atelier dans Chromium |
| `reseaux/`, `reseaux_test/` | le catalogue de réseaux (`*.g2egrid.json`) et les notebooks de réseaux |
| `docs/` | la documentation Sphinx (guide, théorie, équations des modèles, API) |

Le code, ses commentaires et la documentation sont écrits en français ;
l'interface de l'atelier est bilingue. `g2egrid_api/README.md` décrit
l'atelier et son API HTTP ; `docs/guide/index.md` est l'entrée du guide.

Construire la documentation en local :

```bash
uv run --group docs --extra api sphinx-build -b html docs docs/_build/html
```

## Citer G2Egrid

G2Egrid est un outil de recherche : si vous l'utilisez dans une publication,
merci de le citer. Les métadonnées de citation sont dans
[`CITATION.cff`](CITATION.cff).

Pour la méthode de modélisation petit-signal sous-jacente, citer aussi :
F. Kelada, J. Buire, N. HadjSaid, « Revisiting small-signal modeling for
analyzing fast dynamic interactions in converter-dominated power systems »,
*IEEE Transactions on Power Systems*, 40(5):3819-3832, 2025.

## Licence

G2Egrid est un logiciel libre, distribué sous
[licence publique générale GNU, version 3 ou ultérieure](LICENSE)
(`GPL-3.0-or-later`) : vous pouvez le redistribuer et le modifier selon les
termes de cette licence. Il est fourni **sans aucune garantie**.
