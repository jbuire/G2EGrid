# Équivalence des modèles Python / MATLAB

Ce document est la référence de correspondance entre `g2egrid/g2egrid/elements/`
et `matlab/G2Elin_V2/core/`. Il liste, modèle par modèle, les états, les valeurs par
défaut et **tous les écarts qui subsistent**.

Statut : les équipements sont des portages fidèles (mêmes équations, mêmes
états, mêmes noms d'états, mêmes valeurs par défaut). Les trois incohérences
relevées dans le code MATLAB (§B, §D, §H) ont été **corrigées côté MATLAB** ;
le portage Python utilisait déjà la forme physique. Ce qui reste sont des
différences d'**assemblage** (§A, §C, §E) et deux réglages à revoir (§F, §G).

Le choix du repère — et donc le caractère *local* ou non des modèles
d'équipement — est traité à part dans
[`reperes_et_angles.md`](reperes_et_angles.md).

## Table des modèles

| MATLAB | Python | États | Algébriques | Statut |
|---|---|---|---|---|
| `SM.m` | `SynchronousMachine` | 17 | 6 | équivalent |
| `IB.m` | `IB` (alias `VoltageSource`) | 3 | 0 | équivalent, voir §A |
| `TransmissionLine.m` | `Line` | 2 | 0 | équivalent |
| `Transformer.m` | `Transformer` | 2 | 0 | équivalent |
| `GenericLoad.m` | `RLLoad` | 4 | 0 | équivalent, voir §D |
| `Node.m` (C > 0) | `Bus` | 2 | 0 | équivalent |
| `Node.m` (C = 0) | `ResistiveBus` | 0 | 0 | équivalent, voir §C |
| `GFM1.m` | `GFM1` | 13 | 0 | équivalent |
| `GFM2.m` | `GFM2` | 14 | 0 | équivalent |
| `GFM3.m` | `GFM3` | 15 | 0 | équivalent |
| `GFM4.m` | `GFM4` | 16 | 0 | équivalent |
| `GFL1.m` | `GFL1` | 12 | 0 | équivalent |
| `GFL2.m` | `GFL2` | 14 | 0 | équivalent, voir §F |
| `datacenter/DataCenter.m` | site de `g2egrid.datacenter` (sans groupe) : `Transformer`, `Bus`, `RedresseurDiode`, `BusDC`, `BatterieDCDC`, `OnduleurVF`, `ChargeACPref` | 23 | 0 | équivalent, voir §I |

Les variables algébriques `md`, `mq`, `w` des convertisseurs sont substituées
en ligne des deux côtés (ODE pure). Pour la machine, MATLAB inverse les
matrices d'inductance symboliquement ; le portage garde les six courants
comme variables algébriques et les élimine par complément de Schur — c'est le
même modèle.

## Machine synchrone — `SM.m` ↔ `SynchronousMachine`

États, dans l'ordre :

```
phi_d phi_q phi_fd phi_1d phi_1q phi_2q   flux (6 enroulements)
Dwr Pm                                    mécanique + régulateur de vitesse
Dw1 v1 v2 vpss                            PSS
e1 e2 efd e3                              AVR (mesure, ampli, excitatrice, stab.)
theta                                     angle rotor
```

Valeurs par défaut (`SM_DEFAULT`), reprises littéralement du constructeur :

| Groupe | Paramètres |
|---|---|
| Électrique | `Ra=0.003, Ll=0.15, Lad=1.66, Laq=1.61, Lfd=0.165, L1d=0.1713, L1q=0.7252, L2q=0.125, Rfd=0.0006, R1d=0.0284, R1q=0.00619, R2q=0.02368, Rg=100` |
| Mécanique | `H=5, KD=1` |
| Régulateur de vitesse | `mp=0.5` (en **%**, comme MATLAB), `Tg=0.2` |
| AVR | `Tr=0.02, Ka=300, Ta=0.001, Ke=1, Te=1e-4, Kfd=0.001, Tfd=0.1` |
| PSS | `K_PSS=20, T_LP=0.03, T_HP=2, T1n=0.05, T1d=0.02, T2n=3, T2d=5.4` |

La résistance de causalité `Rg` est portée à l'identique : le courant injecté
vaut `ig = i - v/Rg` et c'est `ig` qui porte la chute `Ra`. Elle figure aussi
dans le Y-bus du load flow, comme dans `Ybus_creation`.

Les butées `VSmin`/`VSmax` du PSS et `Efmin`/`Efmax` de l'AVR de MATLAB
n'interviennent pas dans le modèle linéarisé, et ne sont pas portées sous cette
forme. Le modèle Python a en revanche la butée de l'amplificateur de l'AVR
(`VRmax`/`VRmin` = ±5 pu sur `e2`, drapeau `mode_e2`) : c'est un mode discret
de la simulation temporelle (voir {mod}`g2egrid.elements.limiteurs`), nul à
l'équilibre, si bien que le linéarisé reste celui de MATLAB au bit près.

## Convertisseurs

Gains PI auto-réglés depuis les temps de réponse, formules identiques à
`computeControllerGains()` :

```
KpVL = 6/(K0_VL·tau_VL)                  K0_VL   = wb/C_pu
KiVL = (3/(zeta_VL·tau_VL))²/K0_VL
KpCL = (6·tau0_CL - tau_CL)/(K0_CL·tau_CL)   K0_CL   = vdc0/R1_pu
KiCL = (3/(zeta_CL·tau_CL))²·tau0_CL/K0_CL   tau0_CL = L1_pu/(wb·R1_pu)
Kppll = 2·zeta_PLL·wn/wb, Kipll = wn²/wb     wn      = 3/(zeta_PLL·tau_PLL)
```

Défauts GFM (`GFM_DEFAULT`) : `R1_pu=0.01, L1_pu=0.1, C_pu=0.1, R2_pu=0.003,
L2_pu=0.01, mp=0.001, nq=0.001, wf_p=wf_q=0.1·2π·50, wff=Kffi=Kffv=1,
tau_VL=15 ms, tau_CL=0.1 ms, zeta=0.707` ; GFM3/4 ajoutent
`Cdc_pu=0.1, Gdc_pu=17.28e-6, Tdc=1e-4, tau_DC=5 ms`.

Défauts GFL (`GFL_DEFAULT`) : `R1_pu=0.015, L1_pu=0.1, C_pu=0.11, R2_pu=0.003,
L2_pu=0.01, tau_CL=10 ms, tau_PL=tau_QL=0.1 s, tau_PLL=50 ms, zeta=0.7` ;
GFL2 ajoute `Cdc_pu=0.2, Gdc_pu=17e-6, Tdc=1e-3, tau_VL=0.1 s`.

Comme dans MATLAB, `Kpdc` de GFM3/GFM4/GFL2 est recalculé **au point de
fonctionnement** (correction de la charge à puissance constante).

## Écarts d'assemblage

### §A — Référence d'angle (seule différence de fond restante)

MATLAB désigne un élément *slack* dont l'angle devient la référence du réseau :
tous les autres éléments reçoivent `theta_slack` en entrée. Le portage utilise
un **repère synchrone commun** tournant à `w_ref`, que personne ne possède, et
chaque élément porte son propre angle absolu comme état.

Les deux formulations sont reliées par une rotation globale d'angle
`theta_slack(t)` : **les valeurs propres sont les mêmes**, seuls les angles
affichés diffèrent d'une constante et le mode de jauge (λ = 0, présent des deux
côtés) change d'états porteurs.

C'est ce choix qui décide si les modèles d'équipement sont *locaux*. Explication
détaillée, conséquences et plan de migration côté MATLAB :
[`reperes_et_angles.md`](reperes_et_angles.md).

### §B — Couplage de fréquence `w_ref` — ✅ corrigé côté MATLAB

Côté MATLAB, `w_slack` est rebouclé du slack vers tous les éléments, mais les
colonnes correspondantes de `B` n'étaient pas homogènes :

| Modèle MATLAB | ∂f/∂w_pu écrit | valeur physique |
|---|---|---|
| `SM.m` | `wb·(...)` (cohérent) | ✔ |
| `TransmissionLine.m` | `i_lq` ❌ | `wb·i_lq` |
| `Transformer.m` | `i_2q` ❌ | `wb·i_2q` |
| `Node.m` | `v_gq0` ❌ | `wb·v_gq0` |
| `GenericLoad.m` | `wb·i_cq0` (cohérent) | ✔ |
| `IB.m` | `i_q0·wb` (cohérent) | ✔ |

Il manquait donc un facteur `wb` sur les lignes, transformateurs et nœuds.
Comme `w_slack` est un signal **interne** (rebouclé), l'écart se retrouvait
dans la matrice `A` du système bouclé, pas seulement dans `B`.

**Corrigé** dans `TransmissionLine.m`, `Transformer.m` et `Node.m` : les trois
colonnes portent maintenant le facteur `wb`, comme le faisaient déjà `SM.m`,
`GenericLoad.m` et `IB.m`. Le portage Python utilisait déjà la forme physique.

### §C — Nœuds et shunts dans le load flow

Règle identique des deux côtés : un nœud reçoit la moitié de la charge
capacitive de chaque ligne qui y aboutit ; s'il n'en reçoit aucune, il devient
un nœud algébrique à résistance de causalité `R = 100` pu.

**Différence** : MATLAB met dans le Y-bus du load flow la charge des lignes et
le `Rg` des machines, mais **pas** la résistance de causalité des nœuds. Le
point de fonctionnement ne satisfait alors pas exactement la loi des nœuds du
modèle dynamique (écart de `V/100`, soit 1 % du courant nominal), et la
linéarisation se fait hors équilibre. Le portage injecte aussi ce shunt dans
le Y-bus (`Element.shunt_admittance`), ce qui referme l'équilibre à 1e-9 près
(test `test_loi_des_noeuds`).

### §D — Charge : lignes 3-4 de `B` dans `GenericLoad.m` — ✅ corrigé

```matlab
B = [0,     0,    wb*i_cq0,  1/L_pu,      0;
     0,     0,   -wb*i_cd0,       0, 1/L_pu;
     1/Tau, 0,           0,   1/Tau,      0;     % <- 1/Tau sur vgd
     0, 1/Tau,           0,       0,  1/Tau];    % <- 1/Tau sur vgq
```

Les états `DP_filt`/`DQ_filt` sont des premiers ordres sur les entrées
`DP_load`/`DQ_load` : ils ne dépendent pas de la tension du bus. Les termes
`B(3,4)` et `B(4,5)` étaient un report des lignes 1-2 et injectaient la tension
du bus dans le filtre de puissance.

**Corrigé** dans `GenericLoad.m` (les deux termes sont mis à zéro). Le portage
Python n'a jamais eu ces termes.

### §E — Point d'application de la puissance du load flow (convertisseurs)

À l'équilibre, MATLAB impose `p = P_lf` et `q = Q_lf` sur la puissance
**aux bornes du condensateur du filtre** (`p = ved·igd + veq·igq`), alors que
le load flow est une grandeur **au bus**. L'écart vaut la perte dans `R2`/`L2`
(≈ 0,3 % de la puissance nominale avec les défauts). Le portage cale
l'équilibre sur l'injection au bus, ce qui est cohérent avec le load flow.

### §F — GFL2 : réglage DC instable au-delà de ~0,3 pu

Avec les défauts portés à l'identique (`Cdc_pu = 0.2`, `tau_VL = 0.1 s`,
`tau_CL = 10 ms`), le GFL2 seul est **instable à forte puissance** :

| P injecté | max Re(λ) du modèle élément, `Cdc=0.2` | `Cdc=1.0` |
|---|---|---|
| 0.1 pu | 0 | 0 |
| 0.2 pu | 0 | 0 |
| 0.4 pu | +10.8 | 0 |
| 0.8 pu | +259 | 0 |

Explication : la charge à puissance constante donne au bus DC un pôle instable
en `P·wb/Cdc ≈ 1250 rad/s` à `P = 0,8`, alors que la boucle de courant qui doit
le rattraper est réglée à `1/tau_CL = 100 rad/s`. Le recalage CPL de `Kpdc`
place le pôle voulu à `-3/tau_VL` **en supposant la boucle de courant
instantanée** — hypothèse fausse ici. `Cdc_pu ≥ 1` (ou `tau_CL` plus court)
règle le problème. Le commentaire de `GFL2.m` (« `Cdc_pu >= 0.2` requis »)
signale déjà le sujet ; la valeur retenue n'est pas suffisante.

### §G — PSS : un réglage inadapté, pas une erreur de signe (FR-007)

**Verdict.** La chaîne Δω → PSS → AVR ne contient **aucune erreur de signe ni
de convention**, ni dans le portage, ni dans `SM.m`. Le PSS désamortit parce
que ses constantes ne compensent pas le retard de la chaîne
excitation-machine : à la fréquence des modes locaux du catalogue (≈ 2 Hz), le
transfert GEP(s) retarde de 95 à 115°, et le PSS n'apporte que ≈ 0°. Son
couple s'oppose donc en partie à l'écart de vitesse. Trois causes se
combinent :

1. un PSS **sans avance de phase utile** : avance 0,05/0,02 s, mais passe-bas
   0,03 s et **retard** 3/5,4 s, soit −3° à +3° entre 1 et 2,5 Hz ;
2. un AVR dont le **retour de stabilisation** (`Kfd = 0,001`, `Tfd = 0,1`,
   soit `Ka·Kfd/Tfd = 3`) divise le gain transitoire par 4 et ajoute de 35 à
   70° de retard entre 1 et 2,5 Hz. Sans lui, le même PSS amortit ;
3. un régulateur de vitesse très raide (statisme 0,5 %, `Tg = 0,2 s`), qui
   pousse les modes locaux vers 2 Hz. Le mode rotorique du SMIB du catalogue
   est à 0,9 Hz avec un statisme de 5 % et `Tg = 5 s`, 1,7 Hz avec le réglage
   du catalogue. Or c'est au-delà de ≈ 1,5 Hz que GEP retarde de plus de 90°.

Les constantes ressemblent à celles du PSS de l'exemple 12.6 de Kundur
(*Power System Stability and Control*, 1994 : `K = 20`, `T_W = 10 s`, avance
0,05/0,02, retard 3/5,4), associées à une excitatrice statique sans retour de
stabilisation (`KA = 200`, `TR = 0,01`). Ces valeurs sont **à vérifier dans
l'ouvrage**. *Supposé, non vérifié* : l'AVR et le PSS de `SM.m` reprennent
les valeurs par défaut des blocs Simulink « Excitation System » et « Generic
Power System Stabilizer ». La liste des paramètres de V1 (`Tb`, `Tc`,
`Efmin = −11,5`, `Efmax = 11,5`, `Kp`) a la forme du masque du premier. Une
recherche web sur la documentation MathWorks de 2020b donne `Kf, Tf = 0,001 ;
0,1` pour le premier et `30 ms, 20, 2 s, [0,05 0,02], [3 5,4]` pour le second,
mais les pages n'ont pas pu être relues. Les scripts V1
(`matlab/G2Elin_V1/WSCC/script_WSCC.m`) utilisaient les mêmes constantes avec
`K_PSS = 2`.

#### La chaîne, telle qu'elle est écrite

Identique, au caractère près, dans `sync_machine.py` (l. 238-256),
`SM.m` (l. 312-341) et `symSG.m` de V1 (l. 89-101) :

```
Dwr  = wr - w_ref                      # w_ref : entrée externe, 1 pu  =>  Δω = ω − 1
dDw1  = (Dwr - Dw1)/T_LP                         # transducteur 1/(1+s·T_LP)
dv1   = (T_HP·K_PSS·dDw1 - v1)/T_HP              # K·s·T_HP/(1+s·T_HP)
dv2   = (T1n·dv1 + v1 - v2)/T1d                  # (1+s·T1n)/(1+s·T1d)
dvpss = (T2n·dv2 + v2 - vpss)/T2d                # (1+s·T2n)/(1+s·T2d)
de2   = (Ka·(v_ref - e1 - e3 + vpss) - e2)/Ta    # vpss ajouté à la consigne
defd  = (Ke·e2 - efd)/Te
de3   = (Kfd·defd - e3)/Tfd                      # Kfd·s/(1+s·Tfd)·efd, retranché
```

C'est le schéma IEEE PSS1A (norme IEEE 421.5) sans filtre de torsion ni
butées, avec un gain positif sur Δω = ω − 1 et une sortie ajoutée à la
consigne de l'AVR. C'est aussi le schéma de Kundur (gain, washout,
compensation de phase, sortie au point de sommation de l'excitatrice), et
celui de l'`IEEEST` d'ANDES en mode 1 (`omega - 1`, sortie ajoutée à `vi`).
L'AVR est un `EXST1` d'ANDES (`TB = TC = 0`, `KC = 0`), au pôle de
l'excitatrice `Te = 1e-4 s` près. Une inversion pourrait se cacher à six
endroits. Aucun n'est en cause :

| Endroit | Ce qui est écrit | Vérifié par |
|---|---|---|
| Δω | `wr − w_ref`, soit ω − 1 | `test_entree_du_pss_est_w_moins_1` |
| gain, washout, avances | tous positifs | `test_pss_est_le_schema_pss1a_a_gain_positif` (1e-10) |
| point de sommation de l'AVR | `+vpss`, comme `v_ref` ; `−e1`, `−e3` | `test_sortie_du_pss_ajoutee_a_la_consigne_de_l_avr` |
| `efd` → champ, couple, mouvement | GEP de gain statique positif (+2,1), phase → 0 en basse fréquence | Heffron-Phillips indépendant, ci-dessous |
| angle rotor (couple synchronisant) | une inversion rendrait instable la machine sans régulateur | croisement ANDES sans régulateur (`test_andes.py`) |
| convention de signe d'ensemble | ANDES (IEEE) voit le même désamortissement, et le cas inverse instable | croisement ANDES avec régulateurs, ci-dessous |

#### Deuxième voie : GEP(s)

GEP(s) est le transfert de la consigne de l'AVR vers le couple électrique,
vitesse et angle de **toutes** les machines figés (modèle linéaire complet,
lignes et colonnes `Dwr`, `theta` retirées). Phase en degrés :

| f (Hz) | 0,3 | 0,5 | 0,8 | 1,0 | 1,3 | 1,6 | 2,0 | 2,5 |
|---|---|---|---|---|---|---|---|---|
| `smib` SM1 | −6,6 | −14,8 | −38,7 | −60,9 | −87,2 | −100,5 | −108,2 | −112,1 |
| `wscc9` SM1…SM3 | −4…−11 | −12…−25 | −33…−59 | −52…−78 | −74…−92 | −86…−98 | −98…−105 | −107…−109 |
| `kundur2area` SM1…SM4 | −6 | −14…−15 | −36…−39 | −52…−57 | −69…−75 | −83…−87 | −97…−100 | −105…−107 |
| Kundur 12.6, AVR par défaut | −6 | −14 | −35…−37 | −51…−56 | −68…−75 | −84…−89 | −98…−101 | −107…−109 |
| les quatre réseaux, `Kfd = 0` | −3…−7 | −7…−13 | −12…−21 | −15…−27 | −20…−34 | −26…−41 | −33…−49 | −41…−58 |
| **PSS par défaut** | +10,4 | +6,3 | +3,8 | +2,8 | +1,6 | +0,5 | −1,0 | −3,2 |

Sur le SMIB du catalogue, PSS coupé, le mode local est à 1,6816 Hz
(ζ = 19,054 %). On y trouve GEP = 1,388 ∠ −102,6°, le PSS par défaut
∠ +0,2°, et `Re(GEP·PSS) = −3,49` pu de couple par pu de vitesse : un
amortissement **négatif**. Trois recoupements :

- **Heffron-Phillips**, écrit à la main sans rien de g2egrid (modèle à un
  axe, bus infini, P = 0,9 pu sous 1,03 pu). On trouve K1 = 0,81, K2 = 1,07,
  K3 = 0,39, K6 = 0,49, et GEP(1,68 Hz) = 1,354 ∠ −91,4°. Le modèle complet
  donne −102,6° : l'écart de 11° vient des amortisseurs et du stator.
  **Sans** le retour `Kfd`, les deux modèles donnent −26,8° et −32,7° ;
- **sensibilité de la valeur propre** au gain, en `K_PSS = 0`. On trouve
  dλ/dK = +0,0352 + 0,0377j par unité de gain, aux différences finies comme
  par le résidu exact (accord à 1e-3). L'approximation de l'équation du
  mouvement, −GEP·PSS/(4H), donne +0,0255 + 0,0372j. Partie réelle
  **positive** : un gain positif pousse le mode vers la droite ;
- **sur chaque machine** de `smib`, `wscc9` et `kundur2area`, la phase de
  GEP·PSS au mode local est comprise entre −97,6 et −106,7°.

#### Croisement ANDES avec AVR et PSS actifs

Le réseau est celui de Kundur 12.6 (`andes_croise.kundur_deux_zones`), avec
la fiche machine du livre, et l'AVR et le PSS de `SM_DEFAULT`. Le régulateur
de vitesse est neutralisé des deux côtés : il n'a pas d'homologue dans la
traduction. Côté ANDES, on a `EXST1` et `IEEEST`, voir
`g2egrid/andes_croise.py`. La machine du catalogue ne peut pas servir ici :
X''d = 0,230 ≠ X''q = 0,250, ce que le GENROU refuse, et `vers_andes` lève
désormais une erreur. Le tableau donne ζ en %, ANDES / g2egrid :

| Cas | inter-zones | local 1 | local 2 | max \|Δf\| |
|---|---|---|---|---|
| sans AVR ni PSS | 3,62 / 3,97 (0,54 Hz) | 8,88 / 8,87 | 8,67 / 8,66 | 0,01 % |
| AVR seul | −0,22 / −1,70 (0,62 Hz) | 4,44 / 4,48 | 4,48 / 4,58 | 1,06 % |
| AVR + PSS par défaut | 24,03 / 22,36 (0,59 Hz) | **1,30 / 1,57** | **1,04 / 1,32** | 0,72 % |
| AVR + PSS inversé (K = −20) | −20,23 / −21,35 | −7,62 / −8,09 | −8,01 / −8,47 | 0,78 % |

- **Modes locaux** : l'écart reste sous 0,5 point dans les quatre cas, sous
  le critère de 1 point. Les deux outils voient le PSS par défaut faire
  passer les modes locaux de 4,5 % à 1,0-1,6 %, amortir fortement le mode
  inter-zones (0,6 Hz, où GEP ne retarde que d'environ 20°) et, inversé,
  tout déstabiliser. L'AVR et le PSS sont ainsi validés hors du jumeau
  MATLAB. Leur comportement est un fait du réglage, pas d'une implémentation
  qui serait propre à G2ELIN.
- **Mode inter-zones** : l'écart va jusqu'à 1,7 point, au-delà du critère.
  Il vient de la **f.é.m. de vitesse** du stator. Le GENROU d'ANDES l'écrit
  sans la vitesse (`Flux0`), ce paquet avec (`wr·phi`), et l'AVR l'amplifie,
  puisque la tension qu'il régule porte ce terme. L'écart vaut 0,2, 1,05 et
  1,5 point pour `Ka` = 30, 100 et 300. Une fois ce terme retiré du stator
  (machine de test `MachineSansFemVitesse`), tous les modes reviennent à
  moins de 0,5 point (0,13 / 0,14 / 0,20 / 0,31). C'est aussi l'essentiel
  des 0,35 point publiés sans régulateur, attribués jusqu'ici à « la
  dynamique du réseau et du stator » : c'est le terme de vitesse du stator
  (0,13 point une fois retiré).

#### Des constantes qui amortiraient

Le tableau donne le mode électromécanique le moins amorti entre 0,1 et 3 Hz,
en ζ (%), sur huit configurations : quatre réseaux, chacun avec le régulateur
de vitesse du catalogue (« cat ») puis avec un statisme de 5 % et `Tg = 5 s`
(« 5 % »).

| Réglage | smib cat / 5 % | wscc9 cat / 5 % | kundur2area cat / 5 % | SMIB de test cat / 5 % |
|---|---|---|---|---|
| PSS coupé | 19,1 / **−12,0** | 14,4 / 6,1 | 14,6 / **−3,3** | 15,4 / 3,8 |
| défaut | 13,1 / **−2,2** | 8,7 / 0,2 | 9,5 / 0,6 | 9,8 / **−2,0** |
| (1) avance : K = 20, `T1n = T2n = 0,15`, `T1d = T2d = 0,02` | 39,5 / 29,6 | 29,3 / 39,3 | 39,7 / 25,7 | 48,4 / 54,1 |
| (2) PSS par défaut, AVR sans retour (`Kfd = 0`) | 21,5 / 17,6 | 14,0 / 22,2 | 19,2 / 21,4 | 12,7 / 18,1 |
| PSS des études (K = 10, `0,15/0,05` ×2) | 17,7 / 27,1 | 12,1 / 13,3 | 13,6 / 16,5 | 12,9 / 13,3 |

Les valeurs négatives sont des instabilités. Les trois réglages proposés
sont stables sur les huit configurations, spectre complet. Le réglage (1)
apporte +39 à +78° entre 0,3 et 2,5 Hz, ce qui ramène GEP·PSS entre −36° et
+35° pour toutes les machines des quatre réseaux : du couple d'amortissement
à toutes les fréquences de la bande. ANDES le confirme sur Kundur 12.6, à
24,7 / 23,7 % pour le mode inter-zones et de 48 à 59 % pour les autres modes
de la bande. Son gain haute fréquence est élevé
(`K·(T1n/T1d)·(T2n/T2d)` ≈ 1 100), ce qui n'est vérifié ni pour le bruit ni
pour la torsion. Le réglage (2) garde le PSS de Kundur et retire à l'AVR
le retour qui en faussait la compensation. ANDES le confirme aussi
(35,0 / 33,2 % locaux, 20,5 / 19,2 % inter-zones), mais il reste en deçà du
PSS coupé sur le SMIB de test raide (12,7 contre 15,4 %). Le PSS « corrigé »
de `wscc_saturation` et `inertie_virtuelle` n'apporte que +32 à +47° : il
stabilise les cas à 5 %, mais, avec le régulateur du catalogue, il reste
sous le PSS coupé.

Le choix revient à l'utilisateur (FR-015, profils de réglage). Toute
correction déplace l'invariant du SMIB de test (2,07 Hz, 9,8 %). Elle doit
être faite des deux côtés, dans `SM.m` comme ici.

#### Reproduire

```bash
cd g2egrid
uv sync --all-extras --group validation
uv run pytest -q tests/test_pss_avr.py   # 32 réussis, 3 xfail (comportement attendu, strict)
```

Les tests chiffrent tout ce qui précède : conventions, GEP, Heffron-Phillips,
sensibilité, tableau FR-007, croisement ANDES et f.é.m. de vitesse.
L'ancienne mesure sur le SMIB de test reste valable (`test_g2egrid.py`,
`test_pss_agit_sur_le_mode_electromecanique`) :

| K_PSS | mode local | ζ |
|---|---|---|
| 0 | 1,95 Hz | 15,3 % |
| +20 (défaut) | 2,06 Hz | **9,8 %** |
| −20 | 1,85 Hz | 22,8 % |

### §H — `IB.m` : `wb = obj.W_ref` — ✅ corrigé

`IB.m` posait `wb = obj.W_ref` (soit 1 pu) là où les commentaires du même code
(« couplage … en rad/s », « Δω_ref en p.u. × wb = rad/s ») indiquent `2πf`.
Avec `wb = 1`, la dynamique électrique de la source est 314 fois trop lente
(pôle à `-R/X = -0,1 s⁻¹` au lieu de `-31,4 s⁻¹`).

**Corrigé** : `wb = obj.wb` (pulsation de base) et `w0 = obj.W_ref` (vitesse
d'équilibre en pu), les termes de couplage devenant `w0·wb`. Effet mesuré sur
le SMIB : négligeable sur le mode électromécanique (1,796 Hz dans les deux
cas), l'impédance de la source étant très faible — mais les modes propres de
la source, eux, passent de 0,1 à 31,4 s⁻¹.

### §I — Centre de données : `core/datacenter/` ↔ `g2egrid.elements.datacenter`

Les blocs élémentaires MATLAB (`RedresseurDiode`, `OnduleurVF`, `ChargeIT`,
`BatterieDCDC`, `TransfoRL`, `CondensateurAC/DC`) portent les équations
non linéaires du portage Python et se linéarisent par différences finies. La
charge `DataCenter` les assemble. Deux différences de découpage, sans effet
sur le spectre :

- **Îlot de l'onduleur** : MATLAB l'écrit dans le repère propre de
  l'onduleur, Python dans le repère global. Python a donc un état de plus,
  l'angle de l'onduleur, qui donne un mode nul.
- **Filtres de mesure de P et Q de l'onduleur** : ils sont sans effet à
  statisme nul et ne sont pas repris en MATLAB. Python a donc deux modes de
  plus, à −wf.

`matlab/G2Elin_V2/reseaux_test/datacenter/test_datacenter.m` le vérifie
contre `reference_python.json`. L'équilibre est identique à 1e-14 près, et
les 32 valeurs propres du réseau bus infini – ligne – centre, assemblé par
`Network.m`, à 1e-9 près.

## Comptage d'états sur les réseaux de test

| Réseau | Python | `reseaux_classiques/README.md` |
|---|---|---|
| SMIB | **28** | 28 ✔ |
| WSCC 9 bus | 93 | 87 |
| Kundur two-area | 114 | 110 |

Les écarts WSCC et Kundur valent exactement **2 états par charge**
(3 charges → 6, 2 charges → 4) : les chiffres du README datent d'avant
l'ajout des états `DP_filt`/`DQ_filt` à `GenericLoad`. Un `buildGlobalModel`
relancé aujourd'hui côté MATLAB devrait donner 93 et 114.
