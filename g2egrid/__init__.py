"""g2egrid — modelisation petit-signal dq de reseaux electriques.

Portage Python de G2ELIN, avec deux changements d'architecture :

1. Couche symbolique-numerique hybride (approche ANDES) : chaque modele est
   declare une fois en SymPy, les jacobiens sont derives analytiquement puis
   ``lambdify``-es en numpy, et le resultat est mis en cache **par classe**.
   La version MATLAB re-executait le calcul symbolique a chaque appel.

2. Loi des noeuds exacte, sans resistance de causalite. Les bus portent leur
   capacite shunt physique : la tension devient un etat, la loi des noeuds
   est exacte, et le modele global est une EDO propre (pas de boucle
   algebrique, pas de shunt fictif de 0,01 pu).

Exemple minimal (SMIB) :

    from g2egrid import Network
    from g2egrid.elements import SynchronousMachine, Line, VoltageSource

    net = Network(f0=50, Sb=1e6)
    net.add_bus(1, 20e3, btype=2, V=1.0, Pg=0.8)      # PV : machine
    net.add_bus(2, 20e3, btype=1, V=1.0)              # slack : reseau amont
    net.add(SynchronousMachine("SM1", 1, Sn=1e6, Sb=1e6, wb=net.wb))
    net.add(Line("L1", 1, 2, R=0.01, X=0.2, wb=net.wb))
    net.add(VoltageSource("GRID", 2, R=1e-4, X=1e-3, wb=net.wb))
    net.run_loadflow()
    ss = net.build_state_space()
    print(ss, ss.is_stable())
"""

# les anciens chemins (g2egrid.modal, g2egrid.network...) avant toute chose :
# voir g2egrid/_compat.py
from g2egrid import _compat

_compat.installer()

from g2egrid.noyau.element import Element, InitContext
from g2egrid.noyau.loadflow import PQ, PV, SLACK, LoadFlowResult, run_loadflow
from g2egrid.reseau.network import BusSpec, Network, StateSpace
from g2egrid.noyau.symbolic import CompiledModel, ModelSpec, clear_cache, compile_spec

__version__ = "0.1.0"

__all__ = [
    "Network",
    "StateSpace",
    "BusSpec",
    "Element",
    "InitContext",
    "ModelSpec",
    "CompiledModel",
    "compile_spec",
    "clear_cache",
    "run_loadflow",
    "LoadFlowResult",
    "SLACK",
    "PV",
    "PQ",
    "__version__",
]


def __getattr__(nom: str):
    """``g2egrid.modal`` et les autres anciens noms de modules (``_compat``)."""
    return _compat.module_ancien(nom)
