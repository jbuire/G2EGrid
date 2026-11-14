# L'atelier web

L'atelier est une interface graphique locale : on y construit un réseau au
glisser-déposer et on le fait passer par le load flow, l'analyse modale, les
balayages de paramètres et la simulation temporelle. Pour le lancer, double-cliquer
`lancer_atelier.bat` (Windows) ou `lancer_atelier.sh` (macOS, Linux) à la racine
du dépôt, ou :

```bash
cd g2egrid
uv run --extra api python -m g2egrid_api       # --port, --sans-navigateur
```

puis ouvrir <http://localhost:8765>. Son fonctionnement interne (deux voies,
continu, zones) est décrit au chapitre {doc}`../manuel/12_atelier`.

```{include} ../../g2egrid_api/README.md
:start-line: 1
:relative-docs: ../docs/
:heading-offset: 1
```
