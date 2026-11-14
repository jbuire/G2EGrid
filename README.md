# G2Egrid

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

*[Version française](README.fr.md)*

**G2Egrid** is a Python tool for stability analysis of power system. Each piece of equipment
declares a symbolic DAE model; the network assembles them into one global
state-space model (A, B, C, D), from which modal analysis, sensitivities,
impedances, continuation and time-domain simulation are derived.

It is developed at G2Elab (Univ. Grenoble Alpes, CNRS, Grenoble INP) and is the
Python successor of the MATLAB tool G2ELIN.

## What it does

- **Models**: synchronous machines (with AVR, PSS, governor, multi-mass
  shafts), grid-forming converters (droop, virtual oscillator, matching),
  grid-following converters, MMC and VSC stations, DC grids, lines,
  transformers, loads; user-defined equipment as **block diagrams** or as
  **editable equations** `dx/dt = f(x, u)`.
- **Operating point**: AC and AC/DC load flow, converter limits.
- **Analyses**: modes, participation factors and mode shapes, sensitivities,
  residues, dq impedances and the generalised Nyquist criterion, short-circuit
  power and SCR, parameter continuation and Hopf bifurcations, robustness to
  power variations and N-1, model simplification, transmission/distribution
  interaction.
- **Time domain**: linear and non-linear simulation of events and scenarios.
- **Formats**: MATPOWER, PSS/E RAW (v33, v34) and pandapower import; MATPOWER,
  MATLAB and notebook export.
- **A catalogue of 34 networks**: SMIB, WSCC 9-bus, Kundur two-area, IEEE 14,
  33, 39 and 118, CIGRE MV, CIGRE B4 DC grid, Nordic 44, and variants.
- **A web workbench** that runs locally: build a network on a single-line
  diagram or in tables, run every analysis, edit block diagrams and equations.
  Its interface is available **in French and in English** (FR / EN button).

## Quick start

[uv](https://docs.astral.sh/uv/) is required.

```bash
git clone https://github.com/jbuire/g2elin_v2.git
cd g2elin_v2/g2egrid
uv sync --all-extras --all-groups
uv run pytest -q --ignore=tests/e2e        # the test suite (about 15 minutes)
```

G2Egrid is used from a clone of the repository: the network catalogue
(`reseaux/`, `reseaux_test/`) sits next to the package, not inside it. For a
non-editable install, point the `G2EGRID_DONNEES` environment variable to
the folder of a clone.

The workbench (a local server and a browser page, on your machine only):

```bash
./lancer_atelier.sh            # Linux, macOS ; lancer_atelier.bat on Windows
# or: uv run --extra api python -m g2egrid_api
```

A first model in a few lines:

```python
from g2egrid.reseau.catalogue import construire

net = construire("wscc9")              # load flow and equilibrium included
ss = net.build_state_space()           # global (A, B, C, D)
print(ss.eig()[:5], ss.is_stable())
```

## Layout

| Path | Content |
|---|---|
| `g2egrid/` | the package: `noyau/` (DAE, symbolic compilation, load flow), `elements/` (equipment library, block diagrams), `reseau/` (`Network`, tables, catalogue, zones), `temporel.py`, `analyse/`, `scenarios/`, `multi_niveau/`, `formats/`, `validation/` |
| `g2egrid_api/` | the workbench: FastAPI server (`app.py`), computations without HTTP (`analyse/`), web interface (`web/`) |
| `tests/` | the pytest suite; `tests/e2e/` drives the workbench in Chromium |
| `reseaux/`, `reseaux_test/` | the network catalogue (`*.g2egrid.json`) and the network notebooks |
| `docs/` | the Sphinx documentation (guide, theory, model equations, API) |

The code, its comments and the documentation are written in French; the
workbench interface is bilingual. `g2egrid_api/README.md` describes the
workbench and its HTTP API; `docs/guide/index.md` is the entry point of the
user guide.

Build the documentation locally:

```bash
uv run --group docs --extra api sphinx-build -b html docs docs/_build/html
```

## Citing G2Egrid

G2Egrid is a research tool: if you use it in a publication, please cite it.
Citation metadata are in [`CITATION.cff`](CITATION.cff).

```bibtex
@software{g2egrid,
  author  = {Buire, Jerome},
  title   = {{G2Egrid}: small-signal modelling of power systems
             with a high share of power electronics},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/jbuire/g2elin_v2},
  license = {GPL-3.0-or-later}
}
```

For the underlying small-signal modelling method, please also cite:
F. Kelada, J. Buire, N. HadjSaid, "Revisiting small-signal modeling for
analyzing fast dynamic interactions in converter-dominated power systems",
*IEEE Transactions on Power Systems*, 40(5):3819-3832, 2025.

## License

G2Egrid is free software, distributed under the
[GNU General Public License, version 3 or later](LICENSE)
(`GPL-3.0-or-later`): you may redistribute and modify it under the terms of
that license. It comes with **no warranty**.
