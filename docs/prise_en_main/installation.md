# Installation

g2egrid demande Python 3.10 ou plus récent. Le dépôt utilise
[uv](https://docs.astral.sh/uv/) pour gérer l'environnement : c'est la voie
recommandée, et la seule qu'emploient les scripts du dépôt.

## Avec uv (recommandé)

```bash
git clone https://github.com/jbuire/g2elin_v2.git
cd g2elin_v2/g2egrid
uv sync                     # dépendances de base : numpy, scipy, sympy, pandas...
uv run pytest               # la suite de tests
```

Les groupes et options facultatifs :

| Commande | Ajoute | Pour |
|---|---|---|
| `uv sync --extra plot` | matplotlib | figures des notebooks et des études |
| `uv sync --extra control` | python-control | `StateSpace.to_control()` |
| `uv sync --extra api` | FastAPI, uvicorn | l'atelier web ({doc}`atelier`) |
| `uv sync --all-extras` | tout ce qui précède | le plus simple pour travailler |
| `uv sync --group validation` | ANDES | validation croisée ({doc}`../guide/validation`) |
| `uv sync --group docs` | Sphinx et extensions | cette documentation |
| `uv sync --all-extras --group e2e` | Playwright | les tests de l'atelier dans un navigateur |

Le groupe `dev` (installé par défaut) contient `httpx`, qu'emploient les tests des
routes de l'atelier.

## L'atelier comme commande (pipx)

```bash
pipx install --editable "./python[api]"
g2egrid-atelier                  # --port, --sans-navigateur, --sans-jeton
```

L'installation **éditable** garde les exemples du catalogue (`reseaux/`,
`reseaux_test/`), qui ne font pas partie du paquet. Pour une installation non
éditable, la variable `G2EGRID_DONNEES` désigne le dossier `g2egrid/` d'un dépôt
cloné.

## En conteneur (Docker)

```bash
cd g2elin_v2/g2egrid
docker build -t g2egrid-atelier .
docker run --rm -p 127.0.0.1:8765:8765 g2egrid-atelier
```

L'adresse à ouvrir, jeton compris, est dans les journaux du conteneur. Le serveur
écoute sur toutes les interfaces *du conteneur* ; c'est la publication du port
sur `127.0.0.1` qui le garde sur la machine locale.

## Avec pip

```bash
cd g2elin_v2/g2egrid
pip install -e ".[plot,control]"
```

## Construire cette documentation

```bash
cd g2elin_v2/g2egrid
uv run --group docs sphinx-build -b html docs docs/_build/html
```

Les exemples des pages du guide consacrées aux nouvelles fonctionnalités, et ceux
de {doc}`../exemples`, sont exécutés par la suite de tests
(`tests/test_docs_exemples.py`) : une page dont un exemple ne tourne plus fait
échouer `uv run pytest`.

La page d'accueil est alors `docs/_build/html/index.html`. Une partie des pages
(le *manuel*) est tirée des chapitres LaTeX de `documentation/` ; après une
modification d'un chapitre, on la régénère par :

```bash
uvx --with pypandoc-binary --with pymupdf python docs/outils/latex_vers_myst.py
```

## Vérifier l'installation

```python
import g2egrid
from g2egrid import catalogue

print(g2egrid.__version__)
net = catalogue.construire("smib")
print(net.build_state_space())
```
