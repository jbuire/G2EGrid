"""Le notebook d'un reseau, fabrique a partir de ses tables.

POURQUOI CE MODULE EXISTE. Les cellules communes des notebooks de
``reseaux_test/`` (imports, assemblage, load flow, modele d'etat, tableau
modal, carte des valeurs propres) vivaient dans
``developpement/generateur_notebooks/build_reseaux.py``, c'est-a-dire dans un
outil, hors du paquet. Consequence : rien de ce qui est INSTALLE ne savait
fabriquer un notebook, et l'interface web ne pouvait rendre qu'une cellule de
donnees a coller a la main. Les cellules descendent donc ici, et
``build_reseaux.py`` n'est plus qu'un appelant.

LA CELLULE D'ASSEMBLAGE EST ``tables.assembler``, ET NON SA COPIE. C'etait la
troisieme version du meme code -- une dans ``tables.py``, une dans
``build_reseaux.py``, une dans ``common_cells.py`` -- et elles avaient deja
diverge : ``common_cells.BUILD_TABLES`` ne connait ni ``T_Loads``, ni
``PARAMS_PAR_BUS``, ni les tables continues, et ramene les ``B_pu`` de
``T_Lines_pu`` sans changement de base. Un notebook bati dessus decrit un autre
reseau que le meme fichier lu par ``catalogue.charger``. Une seule definition,
et elle est dans ``g2egrid/tables.py``.

Les notebooks DEJA ECRITS de ``reseaux_test/`` gardent leur cellule
d'assemblage developpee : ils sont mis en forme a la main, ils se lisent, et
``test_assembleur_identique_au_notebook`` verifie justement que cette cellule
et ``tables.assembler`` construisent le meme reseau. Ce module sert a ce qu'on
ecrit MAINTENANT.
"""

from __future__ import annotations

import pathlib
from typing import Any

from g2egrid.formats.export import exporter
from g2egrid.reseau.tables import Tables

__all__ = [
    "ASSEMBLAGE",
    "CARTE",
    "CARTE_PLOTLY",
    "IMPORTS",
    "LOADFLOW",
    "MODAL",
    "MODELE_ETAT",
    "cellule_donnees",
    "cellules",
    "ecrire",
    "notebook",
]


IMPORTS = """\
import sys
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import display
import plotly.graph_objects as go


# Le paquet g2egrid est a la racine du projet (dossier parent de reseaux_test/).
ROOT = pathlib.Path.cwd()
ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))

from g2egrid.tables import TABLES, Tables, assembler  # noqa: E402

np.set_printoptions(precision=4, suppress=True)
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 20)"""

TITRE_CONSTRUCTION = """\
## 1. Construction du reseau

### 1.1 Donnees et assemblage"""

ASSEMBLAGE = """\
# --- assemblage : les tables ci-dessus, puis g2egrid.tables.assembler -------
#
# L'assemblage n'est pas recopie ici. Il l'a ete longtemps, et il y avait donc
# autant de versions que de notebooks : celle qui lisait T_Loads, celle qui
# ignorait PARAMS_PAR_BUS, celle qui ne changeait pas la base des B_pu de
# T_Lines_pu. Lire l'assemblage, c'est lire g2egrid/tables.py.
#
# globals().get : une cellule de donnees ne declare que les tables dont son
# reseau a besoin ; les autres prennent une table vide, exactement comme le
# fait g2egrid.catalogue.charger.
t = Tables(
    **{nom: list(globals().get(nom, [])) for nom in TABLES},
    f0=f0,
    Sb=Sb,
    SM_PARAMS=dict(globals().get("SM_PARAMS", {})),
    GFM_PARAMS=dict(globals().get("GFM_PARAMS", {})),
    GFL_PARAMS=dict(globals().get("GFL_PARAMS", {})),
    PARAMS_PAR_BUS=dict(globals().get("PARAMS_PAR_BUS", {})),
    nom=NOM_RESEAU,
)
net = assembler(t)

n_gfm = sum(len(getattr(t, f"T_GFM{i}s")) for i in (1, 2, 3, 4))
n_rl = sum(1 for r in t.Y_DER if r[7] or r[8])
print(
    f"{NOM_RESEAU} : {len(t.Y_DER)} bus, "
    f"{len(t.T_Lines_SI) + len(t.T_Lines_pu)} ligne(s), "
    f"{len(t.T_Transfos)} transfo(s), {len(t.T_IBs)} bus infini, "
    f"{len(t.T_SMs)} machine(s), {n_gfm} GFM, "
    f"{len(t.T_GFL1s) + len(t.T_GFL2s)} GFL, {n_rl} charge(s) R-L, "
    f"{len(t.T_Loads)} charge(s) speciale(s)"
)
if t.T_BusDC:
    print(
        f"  continu : {len(t.T_BusDC)} noeud(s), {len(t.T_LinesDC)} cable(s), "
        f"{len(t.T_VSCs)} convertisseur(s) a port continu"
    )"""

TITRE_LOADFLOW = "### 1.2 Load flow"

LOADFLOW = """\
lf = net.run_loadflow()
print(
    f"Load flow : converge en {lf.iterations} iterations (mismatch {lf.mismatch:.2e})\\n"
)

print(
    f"{'bus':>4} {'|V| [pu]':>10} {'angle [deg]':>12} "
    f"{'P inj [pu]':>12} {'Q inj [pu]':>12}"
)
print("-" * 54)
for b, i in net._bus_index().items():
    print(
        f"{b:>4} {abs(lf.V[i]):>10.4f} {np.degrees(np.angle(lf.V[i])):>12.4f} "
        f"{lf.S[i].real:>12.4f} {lf.S[i].imag:>12.4f}"
    )
print(
    f"\\nSomme des injections : P = {lf.S.real.sum():+.4f} pu, "
    f"Q = {lf.S.imag.sum():+.4f} pu  (= pertes)"
)"""

TITRE_MODELE = """\
## 2. Modele d'etat

Equilibre element par element (Newton sur `[f; g] = 0`), puis assemblage
global : cablage des ports par LFT, elimination des variables algebriques
par complement de Schur."""

MODELE_ETAT = """\
net.solve_equilibrium()
ss = net.build_state_space()

print(ss)
print(f"\\n{len(ss.states)} etats, {len(ss.inputs)} entrees, {len(ss.outputs)} sorties")
print(f"stable : {ss.is_stable()}   max Re(lambda) = {ss.eig().real.max():+.4e}")

print("\\nEtats par element :")
seen = []
for s in ss.states:
    owner = s.split(".")[0]
    if owner not in seen:
        seen.append(owner)
for owner in seen:
    noms = [s.split(".", 1)[1] for s in ss.states if s.startswith(owner + ".")]
    print(f"  {owner:<12} {len(noms):>3} : {', '.join(noms)}")"""

TITRE_MODAL = """\
### 2.1 Tableau modal

Une ligne par valeur propre, triee par partie reelle decroissante (modes les
plus lents d'abord), comme `displayModalAnalysis` en MATLAB : une paire
conjuguee occupe deux lignes. L'index `Mode` est l'indice `k` de `lam[k]`.

| Colonne | Definition |
|---|---|
| `Re`, `Im` | parties reelle [1/s] et imaginaire [rad/s] de la valeur propre |
| `f_und [Hz]` | frequence propre non amortie, `abs(lambda) / 2pi` |
| `f_dam [Hz]` | frequence d'oscillation, `abs(Im) / 2pi` |
| `zeta [%]` | amortissement, `-Re / abs(lambda)` |
| `Etat j`, `%j` | j-ieme etat le plus participant et sa participation `abs(v_ki . w_ik)`, normalisee a 100 % par mode |"""

MODAL = """\
lam, Vr = np.linalg.eig(ss.A)
Wl = np.linalg.inv(Vr)
Pf = np.abs(Vr.T) * np.abs(Wl)  # Pf[k, i] : participation de l'etat i au mode k
Pf = Pf / np.maximum(Pf.sum(axis=1, keepdims=True), 1e-30)
top3 = np.argsort(-Pf, axis=1)[:, :3]  # les trois etats les plus participants
k = np.arange(len(lam))

mag = np.abs(lam)
zeta = np.where(mag > 1e-10, -lam.real / np.where(mag > 1e-10, mag, 1.0), 0.0)

T_eig = pd.DataFrame(
    {
        "Re": lam.real,
        "Im": lam.imag,
        "f_und [Hz]": mag / (2 * np.pi),
        "f_dam [Hz]": np.abs(lam.imag) / (2 * np.pi),
        "zeta [%]": 100 * zeta,
    }
)
for j in range(3):
    T_eig[f"Etat {j + 1}"] = [ss.states[i] for i in top3[:, j]]
    T_eig[f"%{j + 1}"] = 100 * Pf[k, top3[:, j]]
T_eig = T_eig.iloc[np.argsort(-lam.real, kind="stable")]
T_eig.index.name = "Mode"



T_eig"""

TITRE_CARTE = "### 2.2 Carte des valeurs propres"

CARTE = """\
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))

ax[0].scatter(lam.real, lam.imag / (2 * np.pi), s=22, marker="x")
ax[0].axvline(0.0, color="k", lw=0.8)
ax[0].set_xlabel("Re($\\\\lambda$) [1/s]")
ax[0].set_ylabel("Im($\\\\lambda$)/2$\\\\pi$ [Hz]")
ax[0].set_title(f"{NOM_RESEAU} - tous les modes ({len(lam)} etats)")
ax[0].grid(True, alpha=0.3)

zoom = np.abs(lam.imag) / (2 * np.pi) <= F_ZOOM
ax[1].scatter(lam[zoom].real, lam[zoom].imag / (2 * np.pi), s=32, marker="x")
for k in np.flatnonzero(zoom & (lam.imag > 1e-6)):
    ax[1].annotate(
        ss.states[top3[k, 0]].split(".")[0],
        (lam[k].real, lam[k].imag / (2 * np.pi)),
        fontsize=7,
        xytext=(4, 3),
        textcoords="offset points",
    )
ax[1].axvline(0.0, color="k", lw=0.8)
ax[1].set_xlabel("Re($\\\\lambda$) [1/s]")
ax[1].set_ylabel("Im($\\\\lambda$)/2$\\\\pi$ [Hz]")
ax[1].set_title(f"modes lents (|f| < {F_ZOOM} Hz), etiquetes par element dominant")
ax[1].grid(True, alpha=0.3)

fig.tight_layout()
plt.show()"""

CARTE_PLOTLY = """\
fig = go.Figure()

fig.add_trace(go.Scatter(x=T_eig["Re"], y=T_eig["Im"], mode="markers"))

fig.update_xaxes(range=[-1000, 10])
fig.update_yaxes(range=[-10, 1000])
fig.show()"""


def cellule_donnees(t: Tables) -> str:
    """La cellule de donnees d'un reseau, au format des notebooks."""
    return exporter(t)


def cellules(intro_md: str, donnees_src: str) -> list[tuple[str, str]]:
    """Les cellules du notebook, comme ``(type, source)``.

    Rendre des couples plutot que des objets ``nbformat`` permet de composer
    un notebook, d'en comparer les cellules ou d'en extraire une seule sans
    dependre de ``nbformat``, qui n'est installe que pour le developpement.
    """
    return [
        ("markdown", intro_md),
        ("code", IMPORTS),
        ("markdown", TITRE_CONSTRUCTION),
        ("code", donnees_src),
        ("code", ASSEMBLAGE),
        ("markdown", TITRE_LOADFLOW),
        ("code", LOADFLOW),
        ("markdown", TITRE_MODELE),
        ("code", MODELE_ETAT),
        ("markdown", TITRE_MODAL),
        ("code", MODAL),
        ("markdown", TITRE_CARTE),
        ("code", CARTE),
        ("code", CARTE_PLOTLY),
    ]


def notebook(intro_md: str, donnees_src: str) -> Any:
    """Un ``nbformat.NotebookNode`` complet, non execute."""
    import nbformat  # noqa: PLC0415

    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(src)
        if genre == "markdown"
        else nbformat.v4.new_code_cell(src)
        for genre, src in cellules(intro_md, donnees_src)
    ]
    nb.metadata.setdefault(
        "kernelspec",
        {"display_name": "Python 3", "language": "python", "name": "python3"},
    )
    return nb


def ecrire(
    chemin: str | pathlib.Path,
    intro_md: str,
    donnees_src: str,
    executer: bool = False,
    timeout: int = 1800,
) -> pathlib.Path:
    """Ecrit le notebook, execute si on le demande, et rend son chemin."""
    import nbformat  # noqa: PLC0415

    chemin = pathlib.Path(chemin)
    nb = notebook(intro_md, donnees_src)
    if executer:
        from nbconvert.preprocessors import ExecutePreprocessor  # noqa: PLC0415

        ExecutePreprocessor(timeout=timeout, kernel_name="python3").preprocess(
            nb, {"metadata": {"path": str(chemin.parent)}}
        )
    chemin.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(nb, chemin)
    return chemin
