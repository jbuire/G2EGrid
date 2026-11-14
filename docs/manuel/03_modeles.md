<!-- Page GENEREE depuis documentation/chapitres/03_modeles.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ch-modeles)=
# Bibliothèque de modèles
Ce chapitre donne, pour chaque équipement, ses équations, ses états, ses paramètres par défaut et la façon dont son équilibre est calculé. Les équations sont celles du code Python ; les classes MATLAB de V2 en sont l'équivalent exact, et le tableau « {ref}`Correspondance des modèles entre les trois versions <tab-correspondance>` » fait le lien entre les trois versions.

:::{table} Correspondance des modèles entre les trois versions.
:name: tab-correspondance

| Équipement | Python (`g2egrid.elements`) | MATLAB V2 | V1 | États |
|:---|:---|:---|---:|---:|
| Nœud capacitif | `Bus` | `Node` ($C>0$) | `symNode` | 2 |
| Nœud algébrique | `ResistiveBus` | `Node` ($C=0$) | -- | 0 |
| Ligne | `Line` | `TransmissionLine` | `symLine` | 2 |
| Transformateur | `Transformer` | `Transformer` | -- | 2 |
| Bus infini | `IB` | `IB` | `symIB` | 3 |
| Charge R-L + échelon | `RLLoad` | `GenericLoad` | -- | 4 |
| Charge à impédance | `ZLoad` | `AdvancedLoad` (mode Z) | `symLoad` | 2 |
| Charge à courant constant | `ILoad` | `AdvancedLoad` (mode I) | -- | 2 |
| Charge à puissance constante | `PLoad` | `AdvancedLoad` (mode P) | -- | 2 |
| Machine synchrone | `SynchronousMachine` | `SM` | -- | 17 |
| Machine à transfo intégré | `SynchronousMachineTR` | `SMTR` | `symSG` | 19 |
| GFM statisme 1^er^ ordre, DC idéal | `GFM1` | `GFM1` | -- | 13 |
| GFM statisme 2^e^ ordre, DC idéal | `GFM2` | `GFM2` | -- | 14 |
| GFM statisme 1^er^ ordre, DC dynamique | `GFM3` | `GFM3` | `symGFM` | 15 |
| GFM statisme 2^e^ ordre, DC dynamique | `GFM4` | `GFM4` | -- | 16 |
| GFL boucles P/Q, DC idéal | `GFL1` | `GFL1` | -- | 12 |
| GFL boucle $v_{dc}$/Q, DC dynamique | `GFL2` | `GFL2` | `symGFL` | 14 |
:::

## Nœuds

Le modèle des nœuds est donné à la section {ref}`sec-noeuds`. En Python, un bus se déclare par `Bus(bus_id, Un, btype=..., V=..., theta=..., Pg=..., Qg=..., Pl=..., Ql=..., C=...)`, avec les puissances en watts et en vars et la capacité en pu. `Network.add` convertit les puissances sur la base $S_b$ et, si $C = 0$, remplace le bus par un nœud idéal (`IdealBus`, loi des nœuds exacte, le défaut) ou par un `ResistiveBus` de résistance `Network.r_causalite` (100 pu, mode `"resistance"`). En MATLAB, `Node` reçoit sa capacité de `addShuntCapacitance`, appelée par `buildFromMatrices` pour chaque ligne, et sa résistance de causalité de `setResistor(100)` si aucune ligne ne lui en apporte ; les réseaux comparés à MATLAB gardent ce traitement (`NOEUDS = "resistance"`).

## Lignes et transformateurs

Une branche est une impédance série $R + jX$, de courant $i$ orienté du port `f` vers le port `t`, avec un rapport de transformation $n$ côté `f` (égal à 1 pour une ligne) : 

$$
\begin{aligned}
\frac{\mathrm{d}i_d}{\mathrm{d}t} &= \frac{\omega_b}{X}\left(\frac{v_{d,f}}{n} - v_{d,t} - R\,i_d + \omega X i_q\right), \\
\frac{\mathrm{d}i_q}{\mathrm{d}t} &= \frac{\omega_b}{X}\left(\frac{v_{q,f}}{n} - v_{q,t} - R\,i_q - \omega X i_d\right).
\end{aligned}
$$

 Les courants injectés valent $-i/n$ au port `f` et $+i$ au port `t`. La susceptance d'une ligne n'est pas portée par la branche : elle est reportée pour moitié sur chacun de ses deux nœuds.

**Bases.** Dans les tables, une ligne `T_Lines_SI` est donnée en $\Omega$/km et S/km, et une ligne `T_Lines_pu` en valeurs *totales* en pu sur sa propre base $(U_n, S_{b,\mathrm{ligne}})$ ; la longueur y est informative. Un transformateur `T_Transfos` est en pu sur $(S_n, U_1)$ ; l'assemblage le ramène à la base système par $Z_{\mathrm{sys}} = Z_{\mathrm{pu}}\,(S_b/S_n)\,(U_1/U_{\text{nœud}})^2$.

## Bus infini

Le bus infini est une force électromotrice $E\angle\delta$ derrière une impédance $R + jX$, avec un état d'angle qui intègre l'écart entre la fréquence qu'il impose et celle du système : 

$$
\begin{aligned}
\frac{\mathrm{d}i_d}{\mathrm{d}t} &= \frac{\omega_b}{X}\left(e_d - v_d - R\,i_d\right) + \omega_b\,\omega\, i_q, &
e_d &= E\cos(\delta + \theta - \theta_g), \\
\frac{\mathrm{d}i_q}{\mathrm{d}t} &= \frac{\omega_b}{X}\left(e_q - v_q - R\,i_q\right) - \omega_b\,\omega\, i_d, &
e_q &= E\sin(\delta + \theta - \theta_g), \\
\frac{\mathrm{d}\theta}{\mathrm{d}t} &= \omega_b\left(\omega_{\mathrm{set}} - \omega_{\mathrm{ref}}\right).
\end{aligned}
$$

 États : `i_d`, `i_q`, `theta`. Entrées de service : `E_ref`, `delta_ref`, `w_set`. L'impédance, en pu sur $S_n$, est ramenée à la base système par $k_S = S_n/S_b$. À l'équilibre, $E = V + (R+jX)\,I$ avec $I = \overline{S/V}$. Quand le bus infini est l'élément de référence, $\theta_g = \theta$ et sa force électromotrice ne tourne pas.

## Charges

### Charge R-L avec échelon de puissance (`RLLoad`, `GenericLoad`)

Une branche R-L série, dont l'impédance est figée à la tension du load flow, complétée par une injection de courant commandée en puissance et filtrée au premier ordre. Cette seconde voie fournit les entrées de service $\Delta p$ et $\Delta q$, qui appliquent un échelon de puissance propre sans modifier l'impédance : 

$$
\begin{aligned}
\frac{\mathrm{d}i_d}{\mathrm{d}t} &= \frac{\omega_b}{X}\left(v_d - R\,i_d + \omega X i_q\right), &
\frac{\mathrm{d}p_f}{\mathrm{d}t} &= \frac{\Delta p - p_f}{T_f}, \\
\frac{\mathrm{d}i_q}{\mathrm{d}t} &= \frac{\omega_b}{X}\left(v_q - R\,i_q - \omega X i_d\right), &
\frac{\mathrm{d}q_f}{\mathrm{d}t} &= \frac{\Delta q - q_f}{T_f},
\end{aligned}
$$

 avec un courant appelé total $i + (p_f v_d + q_f v_q,\; p_f v_q - q_f v_d)/|v|^2$. L'impédance vaut $Z = V_0^2/\overline{S}$, soit $R = P/(P^2+Q^2)$ et $X = Q/(P^2+Q^2)$ à 1 pu. $T_f = 1$ ms par défaut.

:::{admonition} Charges à puissance réactive négative
:class: warning

Pour $Q < 0$, la réactance $X$ serait négative et les pôles de la branche instables, de l'ordre de $+10^4\,\mathrm{s}^{-1}$. Un banc de condensateurs se déclare donc dans `T_Shunts`, et non comme une charge à $Q$ négatif. La classe MATLAB `GenericLoad` calcule $X = Z\sqrt{1-\cos^2\varphi}$, toujours positive : elle reste stable mais modélise alors une charge inductive dont le load flow, lui, voit un $Q$ négatif.
:::

### Charge à impédance constante (`ZLoad`, `AdvancedLoad` mode Z)

La même branche R-L, sans la voie d'échelon : deux états `i_d`, `i_q`. C'est le modèle `symLoad` de V1, à utiliser pour toute comparaison avec la référence.

(sec-charges-speciales)=
### Charges à courant constant et à puissance constante (`ILoad`, `PLoad`)
Ces deux charges n'ont pas d'inductance : elles appellent un courant de consigne qui ne dépend que de la tension à leurs bornes, et le courant réel le suit au premier ordre, comme `AdvancedLoad` en modes I et P : 

$$
T_f\,\frac{\mathrm{d}i_d}{\mathrm{d}t} = k(v)\,\frac{P\,v_d + Q\,v_q}{|v|^2} - i_d, \qquad
T_f\,\frac{\mathrm{d}i_q}{\mathrm{d}t} = k(v)\,\frac{P\,v_q - Q\,v_d}{|v|^2} - i_q,
$$

 avec $P = P_0 + \Delta p$ et $Q = Q_0 + \Delta q$. Le terme en $1/|v|^2$ est le courant $\overline{S/V}$ qui consomme exactement $P + jQ$ ; le facteur $k(v)$ fixe la loi de charge : 

$$
\text{\texttt{PLoad}} : \; k = 1, \quad S(V) = S_0 ; \qquad
\text{\texttt{ILoad}} : \; k = \frac{|v|}{V_0}, \quad S(V) = S_0\,\frac{|V|}{V_0},
$$

 où $V_0$ est la tension du load flow, recalée à l'initialisation pour que la charge consomme $S_0$ au point de fonctionnement. En régime permanent, la puissance active varie donc comme $|V|^n$ avec $n = 2$ pour `ZLoad`, $n = 1$ pour `ILoad` et $n = 0$ pour `PLoad` ; un test du paquet le vérifie. $T_f = 10$ ms par défaut, la valeur de `AdvancedLoad`.

Faute d'inductance, un réactif négatif ne pose aucun problème de stabilité. En revanche, une charge à puissance constante présente une résistance incrémentale négative : sur un réseau faible, elle peut à elle seule déstabiliser le point de fonctionnement, ce qui est précisément le phénomène qu'elle sert à étudier.

(sec-sm)=
## Machine synchrone
```{figure} ../_static/figures/SM_1.png
:name: fig-sm
:width: 100%

Machine synchrone et ses régulateurs (figure reprise de la documentation de V1).
```

Le modèle est celui d'une machine à six enroulements (stator $d$ et $q$, champ, un amortisseur d'axe $d$, deux d'axe $q$), avec un régulateur de vitesse, un régulateur de tension à quatre états et un stabilisateur de puissance. Les paramètres électriques par défaut sont ceux de l'exemple de Kundur.

### Équations électriques

Les flux sont reliés aux courants par les matrices d'inductances, en convention générateur, avec $L_d = L_{ad} + L_l$ et $L_q = L_{aq} + L_l$ : 

$$
\begin{aligned}
\psi_d &= -L_d\, i_d + L_{ad}\, i_{fd} + L_{ad}\, i_{1d}, &
\psi_q &= -L_q\, i_q + L_{aq}\, i_{1q} + L_{aq}\, i_{2q}, \\
\psi_{fd} &= -L_{ad}\, i_d + (L_{ad}+L_{fd})\, i_{fd} + L_{ad}\, i_{1d}, &
\psi_{1q} &= -L_{aq}\, i_q + (L_{aq}+L_{1q})\, i_{1q} + L_{aq}\, i_{2q}, \\
\psi_{1d} &= -L_{ad}\, i_d + L_{ad}\, i_{fd} + (L_{ad}+L_{1d})\, i_{1d}, &
\psi_{2q} &= -L_{aq}\, i_q + L_{aq}\, i_{1q} + (L_{aq}+L_{2q})\, i_{2q}. 
\end{aligned}
$$

 Les six courants sont les variables algébriques du modèle. Les flux évoluent selon, avec $\omega_r = \Delta\omega_r + \omega_{\mathrm{ref}}$ : 

$$
\begin{aligned}
\frac{\mathrm{d}\psi_d}{\mathrm{d}t} &= \omega_b\left(v_{ed} + R_a\, i_{gd} + \omega_r\,\psi_q\right), &
\frac{\mathrm{d}\psi_{fd}}{\mathrm{d}t} &= \frac{\omega_bR_{fd}}{L_{ad}}\, e_{fd} - \omega_bR_{fd}\, i_{fd}, &
\frac{\mathrm{d}\psi_{1d}}{\mathrm{d}t} &= -\omega_bR_{1d}\, i_{1d}, \\
\frac{\mathrm{d}\psi_q}{\mathrm{d}t} &= \omega_b\left(v_{eq} + R_a\, i_{gq} - \omega_r\,\psi_d\right), &
\frac{\mathrm{d}\psi_{1q}}{\mathrm{d}t} &= -\omega_bR_{1q}\, i_{1q}, &
\frac{\mathrm{d}\psi_{2q}}{\mathrm{d}t} &= -\omega_bR_{2q}\, i_{2q}.
\end{aligned}
$$

 La tension de borne $v_e$ est celle du nœud ramenée au repère rotorique par {eq}`eq-rotation`. Une résistance shunt $R_g$ est placée aux bornes : le courant injecté vaut $i_g = i - v_e/R_g$, et c'est lui qui porte la chute $R_a$. Avec la valeur par défaut $R_g = 100$ pu ce shunt est négligeable ; il sert de résistance de causalité. Il est interne au modèle -- le courant injecté en tient déjà compte -- et ne figure donc pas au $Y$-bus du load flow (section {ref}`sec-loadflow`). La rotation {eq}`eq-rotation` porte l'angle $\theta - \theta_g + \theta_{g,0}$, où $\theta_{g,0}$ est la valeur d'équilibre de $\theta_g$, figée à l'initialisation : pour la machine de référence, la rotation est alors constante et le repère global suit son rotor.

### Mécanique et régulateur de vitesse

$$
\begin{aligned}
\frac{\mathrm{d}\Delta\omega_r}{\mathrm{d}t} &= \frac{1}{2H}\left(\frac{P_m}{\omega_r} - C_e - K_D\,\Delta\omega_r\right),
\qquad C_e = \psi_d\, i_q - \psi_q\, i_d, \\
\frac{\mathrm{d}\theta}{\mathrm{d}t} &= \omega_b\,\Delta\omega_r, \qquad
\frac{\mathrm{d}P_m}{\mathrm{d}t} = \frac{1}{T_g}\left(p^\star - \frac{\omega_r - \omega_{\mathrm{ref}}}{m_p} - P_m\right),
\end{aligned}
$$

 où le statisme $m_p$ est donné *en pourcent* dans les paramètres.

### Stabilisateur de puissance et régulateur de tension

Le PSS est un filtre passe-bas, un filtre passe-haut et deux avances de phase : 

$$
\begin{aligned}
\frac{\mathrm{d}\Delta\omega_1}{\mathrm{d}t} &= \frac{\Delta\omega_r - \Delta\omega_1}{T_{LP}}, &
\frac{\mathrm{d}v_1}{\mathrm{d}t} &= \frac{T_{HP}K_{PSS}\,\dot{\Delta\omega}_1 - v_1}{T_{HP}}, \\
\frac{\mathrm{d}v_2}{\mathrm{d}t} &= \frac{T_{1n}\dot{v}_1 + v_1 - v_2}{T_{1d}}, &
\frac{\mathrm{d}v_{pss}}{\mathrm{d}t} &= \frac{T_{2n}\dot{v}_2 + v_2 - v_{pss}}{T_{2d}}.
\end{aligned}
$$

 Le régulateur de tension enchaîne mesure, amplificateur, excitatrice et boucle de stabilisation : 

$$
\begin{aligned}
\frac{\mathrm{d}e_1}{\mathrm{d}t} &= \frac{E_t - e_1}{T_r}, &
\frac{\mathrm{d}e_2}{\mathrm{d}t} &= \frac{K_a(v^\star - e_1 - e_3 + v_{pss}) - e_2}{T_a}, \\
\frac{\mathrm{d}e_{fd}}{\mathrm{d}t} &= \frac{K_e\, e_2 - e_{fd}}{T_e}, &
\frac{\mathrm{d}e_3}{\mathrm{d}t} &= \frac{K_{fd}\,\dot{e}_{fd} - e_3}{T_{fd}},
\end{aligned}
$$

 avec $E_t = |v_e|$ la tension mesurée aux bornes.

### Machine à transformateur intégré (`SynchronousMachineTR`, `SMTR`)

Cette variante, fidèle à `symSG` de V1, place le transformateur de raccordement $R_t + jX_t$ *dans* la machine. Son courant $i_g$ devient un état exprimé dans le repère rotorique, et le port de l'élément est le nœud réseau : 

$$
\begin{aligned}
\frac{\mathrm{d}i_{gd}}{\mathrm{d}t} &= \frac{\omega_b}{X_t}\left(v_{ed} - v_{gd} - R_t\, i_{gd} + \omega_r X_t\, i_{gq}\right), \\
\frac{\mathrm{d}i_{gq}}{\mathrm{d}t} &= \frac{\omega_b}{X_t}\left(v_{eq} - v_{gq} - R_t\, i_{gq} - \omega_r X_t\, i_{gd}\right),
\end{aligned}
$$

 où $v_g$ est la tension du nœud ramenée au repère rotorique. La tension de borne $v_e$ devient interne, portée par le shunt $R_g$ : $v_{ed} = (i_d - i_{gd})\,R_g$ et $v_{eq} = (i_q - i_{gq})\,R_g$, et la chute $R_a$ porte alors le courant d'armature $i$. Les états sont ordonnés comme dans `symSG` : `i_gd`, `i_gq`, les six flux, `Dwr`, `theta`, `Pm`, les quatre états du PSS et les quatre du régulateur de tension, soit 19 états. Le shunt $R_g$ étant derrière le transformateur, il ne figure pas dans le $Y$-bus du load flow.

Décrire la même physique par une machine simple, un nœud terminal et un `Transformer` séparé donne un modèle équivalent, mais le courant du transformateur y est exprimé dans le repère réseau. Les deux matrices $\bm{A}$ diffèrent alors d'une rotation de l'angle rotor : même spectre, coefficients différents.

### Paramètres par défaut et équilibre

:::{table} Paramètres par défaut de la machine synchrone (`SM_DEFAULT`, identiques à `SM.m`).
:name: tab-sm-defaut

| Électrique | Mécanique et vitesse | Régulateur de tension | PSS |
|:---|:---|:---|:---|
| $R_a = 0{,}003$ | $H = 5$ s | $T_r = 20$ ms | $K_{PSS} = 20$ |
| $L_l = 0{,}15$ | $K_D = 1$ | $K_a = 300$ | $T_{LP} = 30$ ms |
| $L_{ad} = 1{,}66$, $L_{aq} = 1{,}61$ | $m_p = 0{,}5$ % | $T_a = 1$ ms | $T_{HP} = 2$ s |
| $L_{fd} = 0{,}165$, $R_{fd} = 0{,}0006$ | $T_g = 0{,}2$ s | $K_e = 1$, $T_e = 0{,}1$ ms | $T_{1n} = 50$ ms, $T_{1d} = 20$ ms |
| $L_{1d} = 0{,}1713$, $R_{1d} = 0{,}0284$ |  | $K_{fd} = 0{,}001$ | $T_{2n} = 3$ s, $T_{2d} = 5{,}4$ s |
| $L_{1q} = 0{,}7252$, $R_{1q} = 0{,}00619$ |  | $T_{fd} = 0{,}1$ s |  |
| $L_{2q} = 0{,}125$, $R_{2q} = 0{,}02368$ |  |  |  |
| $R_g = 100$ |  |  |  |
:::

L'équilibre suit la construction du phaseur classique. L'axe $q$ est aligné sur $E = V + R_a I_g + jL_q I$, où $I = I_g + V/R_g$ est le courant d'armature, et $\theta_0 = \arg E - \pi/2$. Les tensions et courants sont ramenés au repère rotorique, les courants d'amortisseurs sont nuls, les flux s'en déduisent par les équations statoriques à $\dot\psi = 0$, puis $e_{fd} = L_{ad}\, i_{fd}$, $P_m = C_e$, $e_1 = E_t$, $e_2 = e_{fd}/K_e$, $e_3 = 0$ et $v^\star = E_t + e_2/K_a$.

(sec-gfm)=
## Convertisseurs grid-forming
```{figure} ../_static/figures/GFM_1.png
:name: fig-gfm
:width: 100%

Convertisseur grid-forming : filtre LCL, boucles internes de tension et de courant, statisme (figure reprise de la documentation de V1).
```

Les quatre classes `GFM1` à `GFM4` partagent un filtre LCL, deux boucles internes cascadées et un statisme puissance-fréquence et puissance réactive-tension. Elles diffèrent par l'ordre du statisme et par le modèle du bus continu (tableau « {ref}`Correspondance des modèles entre les trois versions <tab-correspondance>` »).

### Filtre LCL

Le convertisseur est raccordé par une inductance $R_1 + jL_1$, un condensateur $C$ aux bornes duquel on mesure $v_e$, puis une seconde branche $R_2 + jL_2$ vers le nœud. Dans le repère local tournant à $\omega$ : 

$$
\begin{aligned}
\frac{\mathrm{d}i_{sd}}{\mathrm{d}t} &= \frac{\omega_b}{L_1}\left(m_d v_{dc} - v_{ed} - R_1 i_{sd} + \omega L_1 i_{sq}\right), &
\frac{\mathrm{d}i_{gd}}{\mathrm{d}t} &= \frac{\omega_b}{L_2}\left(v_{ed} - v_{gd} - R_2 i_{gd} + \omega L_2 i_{gq}\right), \\
\frac{\mathrm{d}i_{sq}}{\mathrm{d}t} &= \frac{\omega_b}{L_1}\left(m_q v_{dc} - v_{eq} - R_1 i_{sq} - \omega L_1 i_{sd}\right), &
\frac{\mathrm{d}i_{gq}}{\mathrm{d}t} &= \frac{\omega_b}{L_2}\left(v_{eq} - v_{gq} - R_2 i_{gq} - \omega L_2 i_{gd}\right), \\
\frac{\mathrm{d}v_{ed}}{\mathrm{d}t} &= \frac{\omega_b}{C}\left(i_{sd} - i_{gd} + \omega C v_{eq}\right), &
\frac{\mathrm{d}v_{eq}}{\mathrm{d}t} &= \frac{\omega_b}{C}\left(i_{sq} - i_{gq} - \omega C v_{ed}\right). 
\end{aligned}
$$

 La seconde branche $R_2 + jL_2$ représente soit l'inductance côté réseau du filtre, soit le transformateur de raccordement (section {ref}`sec-transfo-der`).

### Boucles de commande

Les puissances mesurées sont filtrées au premier ordre : $\dot p_m = \omega_{fp}(p - p_m)$ et $\dot q_m = \omega_{fq}(q - q_m)$, avec $p = v_{ed} i_{gd} + v_{eq} i_{gq}$ et $q = -v_{ed} i_{gq} + v_{eq} i_{gd}$. Le statisme fixe la pulsation et la consigne de tension : 

$$
\begin{aligned}
\text{1\textsuperscript{er} ordre (GFM1, GFM3)} &: \quad \omega = \omega_{\mathrm{ref}} + m_p\,(p^\star - p_m), \\
\text{2\textsuperscript{e} ordre (GFM2, GFM4)} &: \quad \omega = \omega_{\mathrm{ref}} + \Delta\omega, \quad
\frac{\mathrm{d}\Delta\omega}{\mathrm{d}t} = m_p\,\omega_c\,(p^\star - p_m) - \omega_c\,\Delta\omega, \\
&\quad \frac{\mathrm{d}\theta}{\mathrm{d}t} = \omega_b\,(\omega - \omega_{\mathrm{ref}}), \qquad
v_{ed}^\star = v_e^\star + n_q\,(q^\star - q_m), \quad v_{eq}^\star = 0.
\end{aligned}
$$

 La boucle de tension produit la consigne de courant, avec découplage et compensation du courant réseau : 

$$
\begin{aligned}
i_{sd}^\star &= K_{p,VL}(v_{ed}^\star - v_{ed}) + M_{VLd} + K_{ffi}\, i_{gd} - \omega_{ff} C v_{eq}, &
\frac{\mathrm{d}M_{VLd}}{\mathrm{d}t} &= K_{i,VL}(v_{ed}^\star - v_{ed}), \\
i_{sq}^\star &= K_{p,VL}(v_{eq}^\star - v_{eq}) + M_{VLq} + K_{ffi}\, i_{gq} + \omega_{ff} C v_{ed}, &
\frac{\mathrm{d}M_{VLq}}{\mathrm{d}t} &= K_{i,VL}(v_{eq}^\star - v_{eq}),
\end{aligned}
$$

 et la boucle de courant produit les rapports cycliques : 

$$
\begin{aligned}
m_d &= \frac{1}{v_{dc}}\left(K_{p,CL}(i_{sd}^\star - i_{sd}) + K_{ffv} v_{ed} - \omega_{ff} L_1 i_{sq} + M_{CLd}\right), &
\frac{\mathrm{d}M_{CLd}}{\mathrm{d}t} &= K_{i,CL}(i_{sd}^\star - i_{sd}), \\
m_q &= \frac{1}{v_{dc}}\left(K_{p,CL}(i_{sq}^\star - i_{sq}) + K_{ffv} v_{eq} + \omega_{ff} L_1 i_{sd} + M_{CLq}\right), &
\frac{\mathrm{d}M_{CLq}}{\mathrm{d}t} &= K_{i,CL}(i_{sq}^\star - i_{sq}).
\end{aligned}
$$

### Bus continu (GFM3, GFM4)

$$
\frac{\mathrm{d}v_{dc}}{\mathrm{d}t} = \frac{\omega_b}{C_{dc}}\left(i_{dc} - G_{dc} v_{dc} - m_d i_{sd} - m_q i_{sq}\right), \qquad
\frac{\mathrm{d}i_{dc}}{\mathrm{d}t} = \frac{i_{dc}^\star - i_{dc}}{T_{dc}}, \qquad
i_{dc}^\star = \frac{p^\star}{v_{dc}^\star} + K_{p,dc}\,(v_{dc}^\star - v_{dc}).
$$

 Pour GFM1 et GFM2, $v_{dc} = 1$ et ces deux états disparaissent.

### Paramètres et équilibre

:::{table} Paramètres par défaut des GFM (`GFM_DEFAULT` et `GFM_DC_DEFAULT`).
:name: tab-gfm-defaut

| Filtre | Commande | Bus continu (GFM3, GFM4) |
|:---|:---|:---|
| $R_1 = 0{,}01$, $L_1 = 0{,}1$ | $m_p = 0{,}001$, $n_q = 0{,}001$ | $C_{dc} = 0{,}1$ |
| $C = 0{,}1$ | $\omega_{fp} = \omega_{fq} = 0{,}1\,\omega_b$ | $G_{dc} = 17{,}28 \cdot 10^{-6}$ |
| $R_2 = 0{,}003$, $L_2 = 0{,}01$ | $\omega_{ff} = K_{ffi} = K_{ffv} = 1$ | $T_{dc} = 0{,}1$ ms |
|  | $\tau_{VL} = 15$ ms, $\tau_{CL} = 0{,}1$ ms, $\zeta = 0{,}707$ | $\tau_{DC} = 5$ ms |
|  | $\omega_c = 0{,}1\,\omega_b$ (GFM2, GFM4) |  |
:::

Les gains des correcteurs sont calculés à partir des temps de réponse (annexe {ref}`ann-reglage`). Un gain passé explicitement au constructeur n'est pas recalculé. L'équilibre aligne l'axe $d$ local sur la tension du condensateur : $I_g = \overline{S/V}$, $V_e = V + (R_2 + jL_2) I_g$, $\theta_0 = \arg V_e$, $v_{eq} = 0$. On en déduit $i_s = i_g + j\omega C v_e$, les rapports cycliques, les intégrales des correcteurs et les consignes $p^\star = p_0$, $q^\star = q_0$, $v_e^\star = v_{ed,0}$.

(sec-gfl)=
## Convertisseurs grid-following
```{figure} ../_static/figures/GFL_1.png
:name: fig-gfl
:width: 100%

Convertisseur grid-following : filtre LCL, PLL, boucles externes et boucle de courant (figure reprise de la documentation de V1).
```

Un GFL se synchronise sur le réseau par une boucle à verrouillage de phase et commande le courant qu'il injecte. Il ne forme pas la tension : c'est toute la différence avec le GFM.

### PLL et filtre

Le filtre LCL obéit aux mêmes équations que pour le GFM, en remplaçant $\omega$ par la pulsation de la PLL $\omega_{pll}$ et l'angle $\theta$ par $\theta_{pll}$. La PLL se verrouille sur la tension du condensateur, dont elle annule la composante $q$ : 

$$
\omega_{pll} = M_{pll} + K_{p,pll}\, v_{eq} + \omega_{ff}, \qquad
\frac{\mathrm{d}M_{pll}}{\mathrm{d}t} = K_{i,pll}\, v_{eq}, \qquad
\frac{\mathrm{d}\theta_{pll}}{\mathrm{d}t} = \omega_b\,(\omega_{pll} - \omega_{\mathrm{ref}}).
$$

 Les rapports cycliques $m_d$, $m_q$ et la pulsation $\omega_{pll}$ sont des variables algébriques du modèle, comme dans `symGFL`.

### Boucles externes

**GFL1** régule les puissances, avec un bus continu idéal ($v_{dc} = 1$) : 

$$
i_{sd}^\star = K_{pp}(p^\star - p) + M_d, \quad \frac{\mathrm{d}M_d}{\mathrm{d}t} = K_{ip}(p^\star - p), \qquad
i_{sq}^\star = K_{pq}(q - q^\star) + M_q, \quad \frac{\mathrm{d}M_q}{\mathrm{d}t} = K_{iq}(q - q^\star).
$$

 **GFL2** régule la tension continue sur l'axe $d$ et la puissance réactive sur l'axe $q$, avec un bus continu dynamique : 

$$
i_{sd}^\star = K_{p,dc}(v_{dc} - v_{dc}^\star) + M_d, \quad \frac{\mathrm{d}M_d}{\mathrm{d}t} = K_{i,dc}(v_{dc} - v_{dc}^\star),
$$

 la boucle $q$ étant celle de GFL1, et le bus continu celui du GFM3 avec $\dot i_{dc} = (i_{dc}^\star - i_{dc})/T_{dc}$. La boucle de courant est celle du GFM, sans compensation du courant réseau.

### Paramètres et équilibre

:::{table} Paramètres par défaut des GFL (`GFL_DEFAULT` et `GFL2_DEFAULT`).
:name: tab-gfl-defaut

| Filtre | Commande | GFL2 |
|:---|:---|:---|
| $R_1 = 0{,}015$, $L_1 = 0{,}1$ | $\tau_{CL} = 10$ ms, $\zeta_{CL} = 0{,}7$ | $C_{dc} = 0{,}2$, $G_{dc} = 17 \cdot 10^{-6}$ |
| $C = 0{,}11$ | $\tau_{PL} = \tau_{QL} = 100$ ms | $T_{dc} = 1$ ms |
| $R_2 = 0{,}003$, $L_2 = 0{,}01$ | $\tau_{PLL} = 50$ ms, $\zeta_{PLL} = 0{,}7$ | $\tau_{VL} = 100$ ms, $\zeta_{VL} = 0{,}7$ |
|  | $K_{ffv} = \omega_{ff} = 1$ | $K_{iq} = 100$, $K_{pq} = 0$ |
:::

L'équilibre aligne la PLL sur la tension du condensateur, comme pour le GFM, et pose les intégrales à la valeur qui annule leur dérivée : $M_d = i_{sd}$, $M_q = i_{sq}$, $M_{pll} = \omega - \omega_{ff}$. Le mode `init='symGFL'` (`InitV1` en MATLAB) les divise au contraire par leur gain, comme V1 (section {ref}`sec-v1-equilibres`).

(sec-transfo-der)=
## Transformateur de raccordement des équipements
Un équipement de production est presque toujours raccordé par un transformateur. Selon l'équipement et la table utilisée, ce transformateur est porté par l'équipement ou déclaré à part.

:::{table} Où déclarer le transformateur de raccordement.
:name: tab-transfo-der

| Déclaration | Transformateur |
|:---|:---|
| `T_IBs` | inclus : l'impédance $R + jX$ du bus infini |
| `T_SMs` à 3 colonnes | non inclus : le déclarer dans `T_Transfos`, avec un nœud terminal pour la machine |
| `T_SMs` à 5 colonnes `[..., Rt_pu, Xt_pu]` | inclus : machine `SynchronousMachineTR` ou `SMTR`, courant du transformateur dans le repère rotorique |
| `T_GFMxs`, `T_GFLxs` à 4 colonnes | non inclus : la seconde branche du filtre garde ses valeurs par défaut ($R_2 = 0{,}003$, $L_2 = 0{,}01$), le transformateur se déclare dans `T_Transfos` |
| `T_GFMxs`, `T_GFLxs` à 6 colonnes `[..., Rt_pu, Xt_pu]` | inclus : $R_t$ et $X_t$ remplacent $R_2$ et $L_2$, comme dans `symGFM` et `symGFL` |
:::

`T_Transfos` reste disponible pour tout autre transformateur : postes sources, transformateurs élevateurs déclarés à part. Intégrer le transformateur à l'équipement économise un nœud et, pour la machine, reproduit exactement la structure de V1.
