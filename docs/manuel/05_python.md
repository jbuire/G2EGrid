<!-- Page GENEREE depuis documentation/chapitres/05_python.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ch-python)=
# Architecture du portage Python g2egrid
Le paquet `g2egrid` reprend les modèles de V2 avec une architecture différente. Un équipement n'est plus une classe qui assemble ses propres matrices : c'est une *déclaration* de ses variables et de ses équations symboliques. Le noyau en dérive les jacobiennes, les compile en fonctions numériques une seule fois par classe, puis s'occupe de l'équilibre et de l'assemblage pour tous les équipements de la même façon.

## Organisation du paquet

:::{table} Modules du paquet `g2egrid`.
:name: tab-modules

| Module | Rôle |
|:---|:---|
| `element.py` | classe de base `Element`, contexte d'initialisation `InitContext`, conventions de ports et d'entrées globales |
| `symbolic.py` | déclaration `ModelSpec`, compilation `compile_spec` en modèle `CompiledModel`, cache par classe |
| `loadflow.py` | load flow de Newton-Raphson polaire `run_loadflow`, types de bus `SLACK`, `PV`, `PQ` |
| `network.py` | classe `Network` (déclaration, $Y$-bus, load flow, équilibre, assemblage) et résultat `StateSpace` |
| `frequence.py` | réponse fréquentielle d'un équipement : `reponse`, `admittance_dq`, `impedance_dq`, `axes_tension`, `gain` (chapitre {ref}`ch-sources`) |
| `separation.py` | séparation transport / distribution : `separer`, `Separation.a_effectif`, `a_quasi_statique`, `modele_global`, `troncature_equilibree`, `reponse_indicielle` (chapitre {ref}`ch-td`) |
| `elements/` | bibliothèque : `bus.py`, `branch.py`, `source.py`, `load.py`, `sync_machine.py`, `gfm.py`, `gfl.py` |
:::

(déclarer-un-équipement-la-classe-element)=
## Déclarer un équipement : la classe `Element`
Une sous-classe d'`Element` déclare, au niveau de la classe, six tuples de noms et deux méthodes :

:::{table} Ce que déclare une classe d'équipement.
:name: tab-declaration

| Attribut ou méthode | Contenu |
|:---|:---|
| `ports` | noms des ports réseau : `("t",)` pour un élément shunt, `("f", "t")` pour une branche |
| `params` | paramètres numériques, passés au constructeur |
| `states`, `algebs` | états et variables algébriques |
| `svc_inputs`, `svc_outputs` | entrées et sorties de service |
| `equations(s)` | méthode statique qui reçoit un espace de noms de symboles et rend les listes $(\bm{f}, \bm{g}, \bm{h})$ |
| `initialize(ctx)` | calcule l'équilibre à partir du load flow et rend un dictionnaire de toutes les variables et des entrées de service |
:::

Chaque port `p` ajoute automatiquement deux entrées `vd_p`, `vq_p` (la tension du nœud) et deux sorties `id_p`, `iq_p` (le courant injecté). Tous les éléments reçoivent en plus les entrées globales `w_ref` et `theta_g`. Le vecteur d'entrée complet est donc $[\text{service}, \omega_{\mathrm{ref}}, \theta_g, \text{tensions de port}]$, et le vecteur de sortie $[\text{courants de port}, \text{service}]$ : `h` doit rendre les courants de port en premier.

Trois hooks optionnels complètent la déclaration. `shunt_admittance` rend un shunt permanent à ajouter au $Y$-bus, quand le nœud absorbe un courant que l'élément ne compte pas dans ses sorties : c'est le cas du nœud résistif, pas du $R_g$ de la machine, déjà retranché de son courant injecté. L'attribut de classe `passif` marque les éléments sans angle propre (nœuds, branches, charges), dont l'entrée `w_ref` est câblée sur la pulsation de l'élément de référence (`Network.repere = "reference"`, défaut). `scheduled_injection` rend la puissance injectée connue d'avance, $-(P + jQ)$ pour une charge ; un élément qui rend `None` reçoit le reliquat du load flow à son nœud. `set_param` permet à `initialize` de fixer un paramètre qui dépend du point de fonctionnement, comme l'impédance d'une charge.

```python
class ZLoad(Element):
    ports = ("t",)
    params = ("R", "X", "wb")
    states = ("i_d", "i_q")
    svc_outputs = ("p", "q")

    def __init__(self, name, bus, P, Q, wb):
        super().__init__(name, bus, R=1.0, X=1.0, wb=wb)  # R, X fixes a l'equilibre
        self.P, self.Q = float(P), float(Q)

    def scheduled_injection(self):
        return -complex(self.P, self.Q)                  # une charge injecte du negatif

    @staticmethod
    def equations(s):
        f = [(s.wb / s.X) * (s.vd_t - s.R * s.i_d + s.w_ref * s.X * s.i_q),
             (s.wb / s.X) * (s.vq_t - s.R * s.i_q - s.w_ref * s.X * s.i_d)]
        h = [-s.i_d, -s.i_q,                              # courants de port d'abord
             s.vd_t * s.i_d + s.vq_t * s.i_q,
             s.vq_t * s.i_d - s.vd_t * s.i_q]
        return f, [], h

    def initialize(self, ctx):
        V = ctx.v["t"]
        Z = abs(V) ** 2 / np.conj(complex(self.P, self.Q))
        self.set_param("R", Z.real)
        self.set_param("X", Z.imag)
        i = V / Z
        return {"i_d": i.real, "i_q": i.imag}
```

## Compilation symbolique

`Element.compiled()` construit un `ModelSpec` à partir de la déclaration et le passe à `compile_spec`. Celle-ci crée un symbole SymPy par nom, appelle `equations`, vérifie que les listes ont la bonne longueur, dérive les neuf jacobiennes $F_x, F_z, F_u, G_x, \dots, H_u$ et les transforme en fonctions NumPy par `lambdify`, avec élimination des sous-expressions communes. Les paramètres restent symboliques : la compilation est donc faite une fois par *classe*, et tous les équipements du même type la partagent. Le résultat est mis en cache sous le nom du modèle ; deux classes aux équations identiques, comme `Line` et `Transformer`, partagent ce nom par l'attribut `model_name`. `clear_cache()` vide le cache.

Le `CompiledModel` expose les fonctions `f`, `g`, `h` et leurs jacobiennes, toutes de signature `(x, z, u, p)`, le résidu `residual` et une méthode `linearize` qui applique le complément de Schur {eq}`eq-schur` à un élément isolé.

## La classe `Network`

### Déclaration

`Network(f0, Sb)` crée un réseau vide. `add(element)` y ajoute un équipement ou un bus. Un `Bus` est déclaré en unités SI, comme une ligne de `Y_DER` : le réseau convertit ses puissances en pu et le remplace par un `ResistiveBus` si sa capacité est nulle. Les autres équipements doivent être raccordés à des bus déjà déclarés. La méthode historique `add_bus` accepte des puissances en pu.

### Load flow

`build_ybus()` forme la matrice d'admittance à partir des branches, des capacités de bus et des shunts permanents des équipements. `run_loadflow()` résout le load flow et rend un `LoadFlowResult` dont les champs `V` et `S` donnent la tension complexe et l'injection nette de chaque bus, dans l'ordre croissant des numéros de bus (`_bus_index()`), avec `converged`, `iterations` et `mismatch`.

### Équilibre

`solve_equilibrium()` commence par répartir l'injection de chaque bus entre ses équipements (`_element_injections`) et par déterminer l'angle d'équilibre $\theta_{g,0}$ de l'élément de référence. Pour chaque élément, il construit ensuite un `InitContext` portant la tension et la puissance du port, la fréquence et $\theta_{g,0}$. Il appelle `initialize`, assemble le vecteur d'entrée d'équilibre et raffine l'état par `scipy.optimize.root` sur $[\bm{f}; \bm{g}] = \bm{0}$. Un résidu supérieur à $10^{-6}$ arrête le calcul.

Deux attributs d'élément modifient ce comportement. `raffiner_equilibre = False` garde tel quel l'équilibre rendu par `initialize`. `point_init = (V, S)` remplace la tension et la puissance lues sur le load flow *pour le calcul de l'équilibre seulement* : le point de linéarisation garde la tension du vrai nœud.

### Assemblage

`build_state_space(outputs=None)` évalue les jacobiennes de chaque élément à son équilibre, les empile, construit le routage décrit à la section {ref}`sec-assemblage` et élimine les variables algébriques. Il rend un `StateSpace` : matrices `A`, `B`, `C`, `D`, noms `states`, `inputs` et `outputs` préfixés par le nom de l'élément (`SM1.phi_d`, `bus2.vd`), et les méthodes suivantes.

:::{table} Méthodes du résultat `StateSpace`.
:name: tab-statespace

| Méthode | Rôle |
|:---|:---|
| `eig()` | valeurs propres de $\bm{A}$ |
| `is_stable(tol, rel)` | stabilité, mode de jauge toléré (section {ref}`sec-jauge`) |
| `iu(nom)`, `iy(nom)` | indice d'une entrée ou d'une sortie par son nom |
| `to_control()` | objet `control.StateSpace` nommé, si `python-control` est installé |
:::

Les entrées externes du modèle global sont `w_ref` et toutes les entrées de service ; les sorties sont toutes les sorties de service, ou la sélection passée dans `outputs`.

## Modes d'initialisation des équipements

:::{table} Options d'initialisation propres à certains équipements.
:name: tab-options-py

| Classe | Option | Effet |
|:---|:---|:---|
| `SynchronousMachineTR` | `init="physique"` (défaut) | remonte du nœud à la borne à travers le transformateur puis ajoute le courant du shunt |
|  | `init="symSG"` | formules de V1, sans raffinement ; utiliser avec `point_init` |
| `GFL2` | `init="physique"` (défaut) | intégrales à la valeur qui annule leur dérivée |
|  | `init="symGFL"` | intégrales divisées par leur gain, sans raffinement |
| `GFM3`, `GFL2` | `Kpdc` fourni | désactive le recalage du gain continu au point de fonctionnement |
:::

## Exemple minimal hors notebook

```python
import numpy as np
from g2egrid import Network, SLACK, PV
from g2egrid.elements import Bus, IB, Line, SynchronousMachine

net = Network(f0=50.0, Sb=100e6)
net.add(Bus(1, 20e3, btype=PV, V=1.0, Pg=80e6, C=0.01))
net.add(Bus(2, 20e3, btype=SLACK, V=1.0, C=0.01))
net.add(SynchronousMachine("SM1", 1, Sn=100e6, Sb=100e6, wb=net.wb))
net.add(Line("L12", 1, 2, R=0.01, X=0.2, wb=net.wb))
net.add(IB("IB2", 2, R=1e-3, X=1e-2, wb=net.wb))

net.run_loadflow()
net.solve_equilibrium()
ss = net.build_state_space()
print(len(ss.states), ss.is_stable(), ss.eig().real.max())
```

## Tests

La suite `g2egrid/tests/test_g2egrid.py` vérifie le noyau et les modèles : nombre d'états et de variables algébriques de chaque classe, cohérence de l'équilibre avec le load flow, loi des nœuds, gains statiques, comportements de référence des convertisseurs. Elle se lance par `uv run pytest` depuis `g2egrid/`.
