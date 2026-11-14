<!-- Page GENEREE depuis documentation/chapitres/B_reglage.tex par docs/outils/latex_vers_myst.py : ne pas editer. -->

(ann-reglage)=
# Annexe B --- Réglage automatique des correcteurs
Les gains des correcteurs PI se calculent à partir d'un temps de réponse à 5 % $\tau$ et d'un amortissement $\zeta$, par identification à un système du second ordre de pulsation propre $\omega_n = 3/(\zeta\tau)$. Les formules ci-dessous sont celles de `computeControllerGains` en MATLAB et des constructeurs Python ; un gain fourni explicitement n'est jamais recalculé.

## Boucle de courant

L'inductance du filtre vue par le correcteur est un premier ordre de gain $K_0 = v_{dc,0}/R_1$ et de constante de temps $\tau_0 = L_1/(\omega_bR_1)$, avec $v_{dc,0} = 1$ : 

$$
K_{p,CL} = \frac{6\,\tau_0 - \tau_{CL}}{K_0\,\tau_{CL}}, \qquad
K_{i,CL} = \left(\frac{3}{\zeta_{CL}\,\tau_{CL}}\right)^2 \frac{\tau_0}{K_0}.
$$

## Boucle de tension du GFM

Le condensateur du filtre est un intégrateur de gain $K_0 = \omega_b/C$ : 

$$
K_{p,VL} = \frac{6}{K_0\,\tau_{VL}}, \qquad
K_{i,VL} = \left(\frac{3}{\zeta_{VL}\,\tau_{VL}}\right)^2 \frac{1}{K_0}.
$$

## PLL du GFL

Linéarisée à tension nominale, la PLL est un second ordre de gain $\omega_bV_0$ avec $V_0 = 1$ : 

$$
\omega_n = \frac{3}{\zeta_{PLL}\,\tau_{PLL}}, \qquad
K_{i,pll} = \frac{\omega_n^2}{\omega_b}, \qquad
K_{p,pll} = \frac{2\,\zeta_{PLL}\,\omega_n}{\omega_b}.
$$

## Boucles externes du GFL1

Les boucles de puissance sont des intégrateurs purs réglés au premier ordre : $K_{ip} = 3/\tau_{PL}$, $K_{iq} = 3/\tau_{QL}$ et $K_{pp} = K_{pq} = 0$.

## Boucle de tension continue

Le bus continu est un premier ordre de gain $K_0 = 1/G_{dc}$ et de constante de temps $\tau_0 = C_{dc}/(G_{dc}\,\omega_b)$. Mais au point de fonctionnement, la puissance $P_0$ extraite par le convertisseur se comporte comme une conductance négative : la conductance effective vaut $G_{\mathrm{eff}} = G_{dc} - P_0/v_{dc,0}^2$. Les classes recalent donc le gain proportionnel au point de fonctionnement pour placer le pôle du bus continu à $-3/\tau$ : 

$$
\begin{aligned}
\text{GFM3, GFM4} &: \quad K_{p,dc} = \frac{3\,C_{dc}}{\tau_{DC}\,\omega_b} - G_{\mathrm{eff}}, \\
\text{GFL2} &: \quad K_{p,dc} = \frac{1}{m_{d,0}}\left(\frac{3\,C_{dc}}{\tau_{VL}\,\omega_b} - G_{\mathrm{eff}}\right),
\qquad K_{i,dc} = 0.
\end{aligned}
$$

 V1 règle ces mêmes boucles sur $G_{dc}$ au lieu de $G_{\mathrm{eff}}$, et, pour le GFL, en $2\zeta\omega_n\tau_0 - 1$ avec un terme intégral non nul : 

$$
K_{p,dc}^{V1{}} = \frac{2\,\zeta\,\omega_n\,C_{dc}}{\omega_b} - G_{dc}, \qquad
K_{i,dc}^{V1{}} = \frac{\omega_n^2\,C_{dc}}{\omega_b}.
$$

 Pour reproduire V1, on passe ces gains explicitement et on désactive le recalage (`tau_DC = []` ou `tau_VL = []` en MATLAB, `Kpdc` fourni en Python).
