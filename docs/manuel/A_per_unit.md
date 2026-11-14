<!-- Page GENEREE depuis documentation/chapitres/A_per_unit.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ann-pu)=
# Annexe A --- Système per-unit
Toutes les grandeurs des modèles sont en per-unit. Les bases alternatives et continues sont celles de la fonction `PU_calc` de V1, reprises par V2 et par Python.

## Bases alternatives

Pour une puissance de base $S_b$ et une tension nominale entre phases $U_l$, les bases sont définies en valeurs crêtes de phase, ce qui rend la transformation de Park conservative en amplitude : 

$$
U_b = \sqrt{\tfrac{2}{3}}\,U_l, \qquad
I_b = \sqrt{2}\,\frac{S_b}{\sqrt{3}\,U_l}, \qquad
Z_b = \frac{U_b}{I_b} = \frac{U_l^2}{S_b}, \qquad
L_b = \frac{Z_b}{\omega_b}, \qquad
C_b = \frac{1}{Z_b\,\omega_b}.
$$

 Une réactance en pu est donc aussi une inductance en pu : $X_{\mathrm{pu}} =
L_{\mathrm{pu}}$ à la fréquence nominale. C'est pourquoi les équations écrivent $\omega_b/X$ ou $\omega_b/L$ indifféremment.

## Bases continues

Pour le bus continu d'un convertisseur, les bases se déduisent des bases alternatives de façon à conserver la puissance : 

$$
U_{b,dc} = 2\,U_b, \qquad
I_{b,dc} = \tfrac{3}{4}\,I_b, \qquad
Z_{b,dc} = \tfrac{8}{3}\,Z_b, \qquad
C_{b,dc} = \tfrac{3}{8}\,C_b.
$$

 Pour un convertisseur à 18 kV et 1 MVA, $Z_b = 324\,\Omega$ et une capacité de 50 mF vaut $C_{dc} = 50 \cdot 10^{-3} / C_{b,dc} \approx 13\,570$ pu : la constante de temps du bus continu est alors très grande devant celle des boucles, ce qui impose de régler sa boucle de tension en conséquence.

## Bases des équipements et base système

Chaque équipement est paramétré sur sa propre puissance nominale $S_n$. Son courant injecté est ramené à la base système par $k_S = S_n/S_b$, les bases de tension étant celles du nœud. Les tables de lignes en pu portent leur propre base $(U_n, S_{b,\mathrm{ligne}})$, et les transformateurs la leur $(S_n, U_1)$ ; l'assemblage les ramène tous à la base système : 

$$
Z_{\mathrm{sys}} = Z_{\mathrm{pu}}\,\frac{Z_{b,\mathrm{propre}}}{Z_{b,\mathrm{sys}}},
\qquad
B_{\mathrm{sys}} = B_{\mathrm{pu}}\,\frac{Z_{b,\mathrm{sys}}}{Z_{b,\mathrm{propre}}}.
$$
