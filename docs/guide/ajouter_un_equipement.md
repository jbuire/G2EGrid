# Ajouter un équipement

Une classe d'équipement **déclare** ses variables et ses équations ; les
jacobiens et le code numérique sont dérivés automatiquement, une fois par
classe (chapitre {doc}`../manuel/05_python`).

```python
from g2egrid import Element

class MonElement(Element):
    ports = ("t",)                       # un port alternatif, raccordé à un bus
    params = ("R", "X", "wb")
    states = ("i_d", "i_q")
    algebs = ()                          # variables algébriques éventuelles
    svc_inputs = ("consigne",)           # consignes, exposées dans le modèle global
    svc_outputs = ("p",)                 # mesures, exposées aussi

    @staticmethod
    def equations(s):                    # symbolique (SymPy)
        f = [(s.wb / s.X) * (s.vd_t - s.R * s.i_d + s.w_ref * s.X * s.i_q),
             (s.wb / s.X) * (s.vq_t - s.R * s.i_q - s.w_ref * s.X * s.i_d)]
        h = [-s.i_d, -s.i_q, s.vd_t * s.i_d + s.vq_t * s.i_q]
        return f, [], h                  # (f, g, h)

    def initialize(self, ctx):           # numérique, par instance
        I = ctx.v["t"] / complex(self.pdict["R"], self.pdict["X"])
        return {"i_d": I.real, "i_q": I.imag, "consigne": 0.0}
```

Les règles :

- chaque port `p` ajoute les entrées `vd_p`, `vq_p` (tension du bus) et les
  sorties `id_p`, `iq_p` (courant **injecté** dans le bus) ; tous les éléments
  reçoivent aussi `w_ref` et `theta_g` ;
- `h` rend d'abord les courants de chaque port, puis les sorties de service ;
- `initialize` rend la valeur d'équilibre de **toutes** les variables et des
  entrées de service, à partir de la tension et de la puissance que le load
  flow attribue à l'élément (`ctx.v`, `ctx.s`) ; un Newton la raffine ensuite ;
- les conventions de signe et de repère sont celles de {mod}`g2egrid.noyau.element`.

La nouvelle classe bénéficie aussitôt de tout le reste : sélection des
dynamiques ({doc}`dynamiques`), sensibilités ({doc}`sensibilites`), simulation.
Ses équations s'ajoutent à la page {doc}`../equations_modeles` en la
régénérant :

```bash
uv run python -m g2egrid.elements.equations > docs/equations_modeles.md
```
