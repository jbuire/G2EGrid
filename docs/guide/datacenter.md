# Un centre de données et son groupe de secours : blocs, courbe de charge, fatigue

Un centre de données alimenté par le réseau MT se décrit par une chaîne de
**blocs élémentaires**, chacun un {class}`~g2egrid.Element` ordinaire. Le même
modèle sert à la simulation **non linéaire** ({class}`~g2egrid.temporel.ModeleNonLineaire`)
et au modèle d'état **linéarisé** (`Network.build_state_space`) :

```
GRID --- L_AMONT --- [2] ---+--- T_DC --- [4] RED === bus continu [10] === OND --- [5] IT
 (20 kV)             MT     |                          |
                            |                         BAT
                            +--- T_GE --- [3] GE (diesel + alternateur)
```

| Bloc | Classe | Ce qu'il modélise |
|---|---|---|
| transformateur | `Transformer` | MT/BT, branche R-L (déjà dans le noyau) |
| condensateur alternatif | capacité du bus (`add_bus(..., C=...)`) | filtre d'entrée du redresseur |
| redresseur à diodes | {class}`~g2egrid.elements.RedresseurDiode` | pont triphasé, modèle moyen dq, 1 état |
| condensateur continu | `BusDC` (`C=`) | le condensateur du bus continu |
| onduleur V/f constants | {class}`~g2egrid.elements.OnduleurVF` | tension et fréquence de consigne, sans statisme, 13 états |
| charge IT | {class}`~g2egrid.elements.ChargeACPref` | puissance constante, entrée `p_ref` |
| DC/DC + batterie | {class}`~g2egrid.elements.BatterieDCDC` | soutien du bus continu, lissage, état de charge |
| groupe de secours | {class}`~g2egrid.elements.GroupeElectrogene` | alternateur + moteur diesel + arbre élastique, 22 états |

Leurs équations sont dans {doc}`../equations_modeles`. Les mêmes blocs
existent en MATLAB (`matlab/G2Elin_V2/core/datacenter/`, charge composite
`DataCenter`), vérifiés contre ce portage ({doc}`../equivalence_matlab`, §I),
avec une bibliothèque Simulink de blocs élémentaires et du bloc composite.

## Le redresseur à diodes

Modèle **moyen** du pont triphasé : la tension continue moyenne est
proportionnelle à l'amplitude de la tension alternative, diminuée de la chute de
commutation (non dissipative), et le fondamental du courant alternatif suit la
tension, au déplacement de commutation près ::

$$
V_r = k_v\,|V| ,\qquad
\frac{L_t}{\omega_b}\frac{di_{dc}}{dt} = g(V_r - v_{dc}) - (R_{dc}+R_{com})\,i_{dc} ,\qquad
P_{ac} = (V_r - R_{com}\,i_{dc})\,i_{dc}
$$

En ohms, la chute de commutation est `3·Xc/π`. En pu, avec pour base continue la
tension redressée à vide `3√2/π·U` et la même puissance de base, l'impédance de
base continue vaut `18/π²` fois l'alternative : `Rcom = π·Xcom/6` (`Xcom` en pu
alternatif), et les deux inductances de commutation en série côté continu
ajoutent `π²·Xcom/9` à l'inductance de lissage, `Lt = Ld + π²·Xcom/9`. Exemple :
pont de 2 MVA sous 400 V, `Xcom = 0,03` pu ; `Vd0 = 540 V`, `Id_base = 3,70 kA`,
`Xc = 2,4 mΩ`, chute `3·Xc·Id/π = 8,5 V = 0,0157 pu = π·0,03/6`.

Le pont consomme du réactif : le fondamental du courant est déphasé par
l'empiètement `μ`, que donne `Vd/Vd0 = (V_r − R_com·i_dc)/V_r = (1 + cos μ)/2`.
Pour un pont à diodes, `tan φ = (2μ − sin 2μ)/(1 − cos 2μ)` (forme exacte,
vérifiée dans les tests par la série de Fourier du courant de phase ;
l'approximation `cos φ ≈ Vd/Vd0` surestime `tan φ` de 6 %). `Xcom` est exprimé
dans la base alternative du pont, `kv` jouant le rôle d'un transformateur idéal
entre le bus et le pont. Le modèle garde `Q = tphi·P` avec `tphi` calé au point
de fonctionnement ;
{func}`~g2egrid.elements.datacenter.reactif_commutation` en donne la valeur
physique, que {func}`~g2egrid.datacenter.reseau_datacenter` impose au load flow
(`tphi = 0,1455` à pleine charge pour `Xcom = 0,03`).

La conduction est une caractéristique **lisse** (`g` = maximum lisse de zéro sur
`eps_v`) : un mode discret blocage/conduction a été essayé d'abord, et ses
basculements répétés quand `V_r` et `v_dc` sont proches coûtaient des dizaines de
milliers de résolutions en 30 ms. Pont conducteur, le modèle est exactement
celui d'une inductance en série avec `Rd`, et le linéarisé est celui du point
conducteur.

Le courant alternatif est une fonction **algébrique** de `idc` et de la
tension (les diodes commutent sur la tension instantanée, sans retard). Un filtre
du premier ordre sur ce courant rend instable, vers 500 Hz, la résonance entre la
capacité du bus alternatif et l'inductance du transformateur. Le bus alternatif
du redresseur doit donc avoir une capacité : c'est justement le condensateur
alternatif du bloc.

`kv` et `tphi` ne sont pas des données mais le **calage** : `initialize` les
règle pour que le pont délivre au point du load flow la tension continue du bus et
le réactif qu'il y consomme. Pour le load flow mixte, le pont est un VSC en mode
`VdcQ` (`net.vsc_modes[bus] = "VdcQ"`).

## L'onduleur V/f constants

Un {class}`~g2egrid.elements.GFM3DC` sans statisme (`mp = nq = 0`) : mêmes filtre
LCL, boucles de tension et de courant et limitations, mais la fréquence est une
consigne `w_set` et le module de la tension une consigne `ve_ref`. L'îlot qu'il
alimente est indépendant de la fréquence du réseau amont ; au load flow il en est
le slack (mode `Vtheta`). Pour le bus continu, c'est une charge à puissance
constante.

## La batterie

Une source de courant du premier ordre (`Tc`) dont la consigne est ::

$$
i_{libre} = \frac{p_{ref}}{v} + K_v\,g(v_{ref} - v) + K_h\,(v_f - v) ,\qquad
i_{ref} = \text{sat}(i_{libre}) ,\qquad
\dot{soc} = -P_{bat}/E_b
$$

* le **statisme est unilatéral** : la batterie ne soutient le bus que sous
  `v_ref`. Un statisme bilatéral partagerait durablement la charge avec le
  redresseur, et la batterie se viderait ;
* **relever `v_ref`** au-dessus de la tension de sortie du pont fait prendre la
  charge à la batterie ; **l'abaisser progressivement** la rend au redresseur :
  c'est la reprise échelonnée du groupe ;
* `Kh` lisse : la batterie n'absorbe que la composante rapide de la tension
  continue (`vf` est la tension filtrée), jamais la puissance moyenne ;
* l'état de charge `soc` est un intégrateur pur : le linéarisé a un mode nul, et
  la batterie ne peut **démarrer** qu'en veille (le point d'équilibre est la veille).

## Le groupe électrogène

L'électrique est celui de {class}`~g2egrid.elements.SynchronousMachine` (17 états,
inchangés). La mécanique est remplacée par une chaîne diesel complète ::

```
régulateur PI + statisme -> actuateur (Ta_g) -> retard d'injection (Padé, Td_c)
   -> combustion + turbo (Tb_c) -> Pm -> masse moteur (He) -- accouplement (Ksh, Dsh) -- masse alternateur (H)
```

Le couple transmis par l'accouplement, `C_arbre = Ksh·dsh + Dsh·(Dwe − Dwr)`, est
une sortie de service : c'est la grandeur de la fatigue. La raideur se donne par
la fréquence du mode de torsion, `f_tors`, et son amortissement relatif :

```python
import numpy as np
from g2egrid.elements import GroupeElectrogene

ge = GroupeElectrogene("GE", 3, Sn=2.5e6, Sb=5e6, wb=2 * np.pi * 50,
                       f_tors=25.0, zeta_tors=0.03)
print(f"raideur Ksh = {ge.pdict['Ksh']:.1f} pu/rad, amortissement Dsh = {ge.pdict['Dsh']:.2f} pu")
```

La butée de carburant (`Fmax`, `Fmin`) est un **mode discret** (`mode_fuel`, comme
les butées de l'AVR) avec rétro-calcul de l'intégrateur. L'entrée `dT_exc` ajoute un
couple moteur : c'est par elle que se superposent la pulsation d'allumage des
cylindres ou l'harmonique 6 du redresseur, que le modèle moyen ne contient pas.

## Construire le site

```python
from g2egrid.datacenter import reseau_datacenter, courbe_charge_ia, scenario_perte_reseau

net = reseau_datacenter()                        # load flow mixte + équilibre
ilot = reseau_datacenter(reseau_amont=False)     # le groupe seul porte le site
ss = net.build_state_space()                     # linéarisé
print(ss, ss.is_stable())
```

`reseau_datacenter` ne construit qu'un réglage **type** (5 MVA de base, 1,5 MW de
charge, groupe de 2,5 MVA) : à remplacer par les données d'un site. Les options
(`avec_batterie`, `lissage`, `C_dc`, `C_ac`, `C_mt`, paramètres du groupe)
permettent les comparaisons de configuration.

## Une courbe de charge de plusieurs minutes

`ChargeACPref` a une entrée absolue `IT.p_ref`. {class}`~g2egrid.temporel.Courbe`
la pilote par une série `(t, P)` — ici le profil synthétique
`courbe_charge_ia`, ou une mesure :

```python
import numpy as np
from g2egrid.temporel import ModeleNonLineaire, Courbe, simuler_lineaire

t = np.arange(0.0, 60.0 + 1e-9, 0.05)              # 60 s ici ; plusieurs minutes de même
P = courbe_charge_ia(t, P_moy=0.30, amplitude=0.25, periode=4.0, t_debut=10.0, rampe=5.0)
m = ModeleNonLineaire(ilot)
sc = m.scenario(Courbe("IT.p_ref", t, P))          # valeur absolue de la consigne
sortie = ["GE.C_arbre"]
tr_nl = m.simuler((0, t[-1]), uext_fn=sc, t_eval=t, rtol=1e-4, atol=1e-6, sorties=sortie)
tr_li = simuler_lineaire(m, t, sc, sorties=sortie)  # même entrée, modèle linéarisé
ecart = np.abs(tr_nl.sorties["GE.C_arbre"] - tr_li.sorties["GE.C_arbre"]).max()
print(f"couple d'arbre : écart linéaire / non linéaire {ecart:.3f} pu sur {np.ptp(tr_nl.sorties['GE.C_arbre']):.3f} pu")
```

Trois interpolations : `"pchip"` (défaut, dérivable, **sans** découpage de
l'intégration), `"lineaire"` et `"palier"` (chaque point est alors une date de
redémarrage de l'intégrateur : à réserver aux courbes courtes ou aux consignes
déclarées par leurs seuls points de cassure).
`simuler_lineaire` discrétise exactement (bloqueur d'ordre un) : le pas ne sert
qu'à échantillonner, une simulation de 10 minutes coûte une fraction de seconde,
contre plusieurs minutes pour le non linéaire.

## Perte du réseau et reprise de charge

`scenario_perte_reseau` enchaîne : la batterie prend la charge, la liaison amont
s'ouvre, la consigne de la batterie redescend et le redresseur, alimenté par le
groupe, reprend la charge progressivement. Le modèle n'a pas de manœuvre de
fermeture, donc le groupe est couplé en attente dès le départ ; la batterie
précède l'ouverture (`anticipation`) pour transposer l'ordre réel — batterie
d'abord, groupe ensuite. Ouvrir la liaison avant donnerait au groupe d'un coup
toute la charge, ce qu'aucun site ne connaît.

**Coût** : une ouverture de ligne coûte environ 40 s de calcul (une résonance LC
excitée par l'interruption du courant, que l'intégrateur résout fidèlement), puis
la simulation est rapide (~0,3 s de calcul par seconde simulée avec une charge
lisse). La capacité des jeux de barres MT (`C_mt`, 0,015 pu, la charge d'un câble
de 5 km) la borne ; avec 1e-3 pu, la même ouverture coûte plusieurs minutes.

## Fatigue mécanique de l'arbre

{mod}`g2egrid.fatigue` ::

```python
from g2egrid.fatigue import evaluer_fatigue

r = evaluer_fatigue(t, tr_nl.sorties["GE.C_arbre"], Sn_VA=2.5e6)
print(r, "marge d'amplitude pour 25 ans :", round(r.marge_amplitude(), 2))
```

Chaîne : couple d'arbre → contrainte de cisaillement `16·C/(π·d³)` → comptage
**rainflow** (ASTM E1049, avec hystérésis `seuil_pu`) → correction de **Goodman** →
courbe S-N de **Basquin** à coude de Haibach → règle de **Miner**.

* L'indicateur à retenir est `marge_amplitude()` : le facteur par lequel on peut
  multiplier toutes les amplitudes de contrainte avant que la vie, la sollicitation
  répétée, tombe sous 25 ans. Le dommage absolu varie comme (contrainte)^7 avec
  les données matériau.
* Les **données matériau par défaut** (`COURBE_ACIER_ARBRE`, `ARBRE_DEFAUT`) sont des
  ordres de grandeur, pas celles d'un arbre précis : seules les comparaisons sont
  robustes.
* Pas de pulsation d'allumage ni d'harmonique 6 dans le couple : voir `dT_exc`.

Le notebook `reseaux_test/datacenter.ipynb` déroule l'ensemble sur dix minutes de
courbe de charge, la perte du réseau, et les balayages (lissage par la batterie,
amortissement de torsion, amplitude et période de la charge) au linéarisé.

## Tests

`tests/test_datacenter.py` (bilan de puissance du pont, jacobiens contre
différences finies, linéarisé contre non linéaire, mode de torsion à `f_tors`,
butée de carburant, scénario d'îlotage court) et `tests/test_fatigue.py` (exemple
de la norme ASTM, Miner sur une sinusoïde, pentes de Basquin, marge d'amplitude).
