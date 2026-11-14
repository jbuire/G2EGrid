# Lois de synchronisation des grid-forming

Les quatre GFM du catalogue (`GFM1` à `GFM4`) se synchronisent par **statisme**
P-f, du premier ou du second ordre. Deux autres lois, très étudiées dans la
littérature, partagent le même filtre LCL, les mêmes boucles internes de tension
et de courant et la même limitation de courant :

| Loi | Classe | Bus continu | Synchronisation |
|---|---|---|---|
| `"statisme"` (défaut) | `GFM1`…`GFM4` | idéal ou dynamique | $\omega = \omega_{ref} + m_p (p_{ref} - p)$ |
| `"dvoc"` | {class}`~g2egrid.elements.gfm.GFMdVOC` | idéal (`T_GFM1s`, `T_GFM2s`) | oscillateur virtuel dispatchable |
| `"matching"` | {class}`~g2egrid.elements.gfm.GFMMatching` | dynamique (`T_GFM3s`, `T_GFM4s`) | $\omega = \omega_{ref}\, v_{dc}/v_{dc,ref}$ |

La loi se choisit par le réglage `loi`, pour une famille (`GFM_PARAMS`) ou un bus
(`PARAMS_PAR_BUS`).

## Oscillateur virtuel dispatchable (dVOC)

```python
from g2egrid import catalogue
from g2egrid.reseau.tables import assembler

t = catalogue.charger("smib_gfm")
gfm = [r for k in ("T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s") for r in getattr(t, k)]
t_voc = t.copie(T_GFM1s=gfm, T_GFM2s=[], T_GFM3s=[], T_GFM4s=[],
                GFM_PARAMS=dict(t.GFM_PARAMS, loi="dvoc", eta=0.02, alpha=5.0))
net = assembler(t_voc)
net.solve_equilibrium()
print([type(e).__name__ for e in net.all_elements if e.name.startswith("GFM")],
      net.build_state_space().is_stable())
```

```text
['GFMdVOC'] True
```

Forme polaire pour un réseau inductif (Colombino, Groß, Dörfler) :

$$
\omega = \omega_{ref} + \eta\Bigl(\frac{p_{ref}}{v_{ref}^2} - \frac{p}{V_r^2}\Bigr), \qquad
\dot V_r = \omega_b\,\eta\Bigl[\frac{q_{ref}}{v_{ref}^2} - \frac{q}{V_r^2}
         + \alpha\,\frac{v_{ref}^2 - V_r^2}{v_{ref}^2}\Bigr] V_r .
$$

$V_r$, l'amplitude de l'oscillateur, est la consigne de la boucle de tension.
Près du nominal, $\eta$ est la pente du statisme P-f et $1/(2\alpha)$ celle du
statisme Q-V équivalent ; la dynamique de tension est intrinsèque à
l'oscillateur.

## Commande par appariement (matching control)

```python
t_mc = t.copie(T_GFM1s=[], T_GFM2s=[], T_GFM3s=gfm, T_GFM4s=[],
               GFM_PARAMS=dict(t.GFM_PARAMS, loi="matching"))
net = assembler(t_mc)
net.solve_equilibrium()
ss = net.build_state_space()
print(ss.is_stable(), [n for n in ss.states if n.startswith("GFM")][:8])
```

```text
True ['GFM2.isd', 'GFM2.isq', 'GFM2.igd', 'GFM2.igq', 'GFM2.ved', 'GFM2.veq', 'GFM2.vdc', 'GFM2.idc']
```

La fréquence du convertisseur **est** la tension continue (Arghir, Jouini,
Dörfler) : le condensateur continu tient le rôle de l'inertie, et le régulateur
de la source continue celui du régulateur de vitesse,
$i_{dc,ref} = p_{ref}/v_{dc,ref} + K_{pdc}(v_{dc,ref} - v_{dc})$. Par défaut,
$K_{pdc} = 1/m_{dc}$, avec un statisme $m_{dc}$ de 2 %. Une source continue qui ne
suit plus se lit immédiatement en fréquence.

## Limitation de courant par impédance virtuelle

Tous les GFM à boucle de tension (`GFM1`…`GFM4`, `GFM3DC`/`GFM4DC`, dVOC,
appariement) peuvent limiter leur courant par une **impédance virtuelle** (VI),
en complément — ou à la place — du limiteur circulaire de
{mod}`g2egrid.elements.limiteurs`. Au-delà du seuil $I_{vi}$, une impédance
$Z_v = R_v + jX_v$ s'ôte de la consigne de tension :

$$
v_{ref}' = v_{ref} - Z_v\,i_s,\qquad R_v = k_{Rv}\,\mathrm{ov}(|i_s|),\qquad
X_v = (X/R)_{vi}\,R_v ,
$$

où $i_s$ est le courant côté convertisseur (celui que lit le limiteur circulaire)
et $\mathrm{ov}$ un dépassement **lissé de classe C1** : nul sous $I_{vi}$,
parabolique sur $[I_{vi}, I_{vi}+\delta]$, linéaire au-delà. Le jacobien est
continu : pas de coin pour le Newton ni pour le petit signal.

| Paramètre | Défaut | Rôle |
|---|---|---|
| `kRv` | 0 (VI **inactive**) | gain de $R_v$ ; 0 laisse les équations de l'ancien modèle au bit près |
| `I_vi` | 1,0 pu | seuil de déclenchement |
| `dI_vi` | 0,05 pu | largeur du raccord parabolique |
| `XsurR_vi` | 0,1 | rapport $X/R$ virtuel |

Ils se règlent comme les autres paramètres de convertisseur (`GFM_PARAMS` ou
`PARAMS_PAR_BUS`) ; l'atelier n'a pas de champ dédié (il ne montre pas non plus
`Imax`). Sous le seuil, $A$, $B$, $C$, $D$ et le statisme sont identiques au bit
près à ceux sans VI ; au-dessus, le modèle linéarisé **inclut** la VI (elle est
dans les équations, pas dans un mode) et l'équilibre reste exact : le repère
local suit le phaseur $v_e + Z_v i_s$ au lieu de $v_e$.

`kRv_dimensionne` donne le $k_{Rv}$ qui fixe le courant **permanent** de défaut
franc à $I_{max}$, à partir de l'impédance de défaut complexe vue par le
convertisseur (filtre, transformateur, défaut) :

```python
from g2egrid.elements.limiteurs import kRv_dimensionne

kRv = kRv_dimensionne(Imax=1.2, I_vi=1.0, XsurR_vi=0.1, Zdefaut=0.0005 + 0.046j)
print(round(kRv, 2))
```

```text
4.7
```

**Ce que la VI ne fait pas, et un constat de conception.** (1) La boucle de
tension fait 15 ms (réglage MATLAB) : sur un défaut franc la VI seule tient le
courant permanent (1,19 pu pour $I_{max}$ = 1,2) mais pas le pic transitoire
(1,67 pu) ; le limiteur circulaire reste le gardien (pic 1,21 pu) et la VI lui
épargne la moitié de ses basculements. (2) Le rapport $X/R$ « élevé » de la
littérature (5 à 10) **n'est pas utilisable ici** : une réactance virtuelle
algébrique dans le repère dq est un gyrateur qui retire de l'amortissement au
mode LC du filtre (instable dès $X/R \ge 1$ pour $R_v$ = 0,3 pu, avec ce réglage
des boucles ; en simulation de défaut, $X/R$ = 0,25 perd déjà le Newton à
l'élimination). Il faudrait une boucle de tension nettement plus rapide, ou une
vraie inductance virtuelle $L_v\,di/dt$. Les chiffres et leurs scripts sont
dans `tests/test_limitation_impedance_virtuelle.py`.

Références : le principe de l'impédance virtuelle pour limiter le courant d'un
GFM est ancien (Paquette et Divan, IEEE Trans. Industry Applications, 2015 ;
Qoria, Gruson, Colas, Guillaud, L2EP : *à vérifier* pour les détails de
chaque publication) ; la forme à dépassement lissé C1 est un choix de ce
paquet.

## Une loi incompatible

```python
try:
    assembler(t.copie(GFM_PARAMS=dict(t.GFM_PARAMS, loi="matching"), T_GFM1s=gfm,
                      T_GFM2s=[], T_GFM3s=[], T_GFM4s=[]))
except ValueError as err:
    print(err)
```

```text
GFM du bus 2 : loi 'matching' impossible avec GFM1 (bus continu ideal) ; possibles : statisme, dvoc
```

Les équations des deux modèles figurent dans {doc}`../equations_modeles`.

## Impédance virtuelle de stabilisation (réseau fort)

Face à un réseau **fort** (peu d'impédance entre le GFM et le reste du réseau), la
boucle de puissance d'un grid-forming a un gain élevé, en $1/X$ : un GFM à statisme
du second ordre (inertie virtuelle) oscille et diverge, à quelques hertz. L'impédance
virtuelle ajoute de l'impédance **dans la commande** : la consigne de la boucle de
tension perd la chute

$$
\underline v_{e,ref} = \underline e - R_v\,\underline i_g - L_v\Bigl(\frac{1}{\omega_b}\frac{d}{dt}
+ j\omega\Bigr)\underline i_g ,
$$

$\underline e$ étant la f.é.m. virtuelle, sur l'axe $d$ du GFM (statisme Q-V ou
oscillateur dVOC). Vu du réseau, le GFM est une source derrière $R_v + jL_v$ de plus.
Trois réglages, communs à tous les GFM (`GFM_PARAMS` ou `PARAMS_PAR_BUS`), en pu de
l'équipement :

| Réglage | Défaut | Rôle |
|---|---|---|
| `Lv_pu` | 0 | inductance virtuelle |
| `Rv_pu` | 0 | résistance virtuelle ; $R_v \approx L_v/5$ amortit aussi le mode du flux de l'inductance totale, vers 50 Hz dans le repère dq |
| `zv_dyn` | 1 | 1 : forme **dynamique** ; 0 : forme **quasi statique** $(R_v + j\omega L_v)\,\underline i_g$ |

La forme dynamique n'a besoin ni de dérivée ni d'état de plus : la chute sur
l'inductance de sortie du filtre la donne exactement,
$L_v(\tfrac{1}{\omega_b}\tfrac{d}{dt} + j\omega)\,\underline i_g = \tfrac{L_v}{L_2}
(\underline v_e - \underline v_g - R_2\,\underline i_g)$. La forme quasi statique, la
plus répandue, ignore le $d/dt$ : avec la boucle de tension par défaut, elle la
déstabilise. Les deux ont le même équilibre ; le repère du GFM s'y aligne sur la
f.é.m. virtuelle (`ve_ref` vaut $|\underline e|$, et $v_{eq} \neq 0$).

```python
import numpy as np


def pire_mode(x_ib, **reglages):
    """Le mode le moins amorti : GFM a statisme du 2e ordre (inertie virtuelle)
    face a un bus infini de reactance x_ib (pu, R = X/10)."""
    ib = list(t.T_IBs[0])
    ib[5], ib[6] = x_ib / 10, x_ib
    t_zv = t.copie(T_GFM1s=[], T_GFM2s=gfm, T_GFM3s=[], T_GFM4s=[], T_IBs=[tuple(ib)],
                   GFM_PARAMS=dict(t.GFM_PARAMS, mp=0.02, **reglages))
    net = assembler(t_zv)
    net.solve_equilibrium()
    lam = np.linalg.eigvals(net.build_state_space().A)
    lam = lam[np.abs(lam) > 1e-6]
    pire = lam[np.argmax(lam.real)]
    return f"{pire.real:+6.1f} 1/s ({abs(pire.imag) / (2 * np.pi):5.1f} Hz)"


for nom, reg in [("sans", {}),
                 ("Lv = 0,1, Rv = 0,02 pu", dict(Lv_pu=0.1, Rv_pu=0.02)),
                 ("idem, quasi statique", dict(Lv_pu=0.1, Rv_pu=0.02, zv_dyn=0.0))]:
    print(f"{nom:24s} fort : {pire_mode(0.01, **reg)}   faible : {pire_mode(0.5, **reg)}")
```

```text
sans                     fort :   +9.2 1/s (  7.5 Hz)   faible :   -7.8 1/s (  2.2 Hz)
Lv = 0,1, Rv = 0,02 pu   fort :   -1.5 1/s (  4.4 Hz)   faible :   -8.4 1/s (  1.9 Hz)
idem, quasi statique     fort : +281.8 1/s (105.6 Hz)   faible :   +8.3 1/s ( 58.2 Hz)
```

- **Réseau fort** ($x = 0{,}01$ pu) : sans impédance virtuelle, le mode de la boucle
  de puissance diverge à 7,5 Hz ; 0,1 pu d'inductance virtuelle dynamique le ramène
  dans le demi-plan gauche.
- **Réseau faible** ($x = 0{,}5$ pu) : elle ne gêne pas.
- **Forme quasi statique**, à la même valeur : un mode de la boucle de tension
  diverge, vers 100 Hz, sur les deux réseaux.

Le schéma-bloc `gfm_complet` ({doc}`schemas_blocs`) porte la même impédance virtuelle
(blocs `zvd`, `zvq`, réglages `Rv`, `Lv`, `zv_dyn`).
