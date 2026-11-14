# Tests et validation

## La suite de tests

```bash
cd g2egrid
uv run pytest
```

Elle vérifie, entre autres :

- le nombre d'états de chaque modèle, conforme à son homologue MATLAB ;
- les jacobiens symboliques contre les différences finies ($10^{-9}$ à $10^{-11}$) ;
- la loi des nœuds à l'équilibre ($10^{-15}$) ;
- les modes de référence des réseaux de test ;
- les équivalences exactes de ce guide : variantes GFM comme sélections de
  GFM4, branches en série à un nœud idéal, repère flottant, sensibilités
  analytiques contre différences finies (`tests/test_dynamiques.py`).

La validation terme à terme contre G2Elin V1 et V2 (MATLAB) est décrite au
chapitre {doc}`../manuel/11_validation`, et la correspondance des classes dans
{doc}`../equivalence_matlab`.

## Validation croisée contre ANDES

Le jumeau MATLAB vérifie ce paquet ligne à ligne, mais c'est un jumeau : même
architecture, mêmes choix de modélisation. {mod}`g2egrid.validation.andes_croise` confronte
le paquet à un outil **indépendant**, [ANDES](https://github.com/CURENT/andes),
sur le réseau à deux zones de Kundur (exemple 12.6).

| | ANDES (GENROU) | g2egrid (Park) | écart |
|---|---|---|---|
| Load flow | | | 6e-7 pu, 2e-6 rad |
| Mode inter-zones | 0,5427 Hz, ζ 3,62 % | 0,5427 Hz, ζ 3,97 % | 0,35 point |
| Premier mode local | 1,0835 Hz, ζ 8,88 % | 1,0835 Hz, ζ 8,87 % | 0,01 point |
| Second mode local | 1,1190 Hz, ζ 8,67 % | 1,1190 Hz, ζ 8,66 % | 0,01 point |

Les régulateurs sont neutralisés des deux côtés. Sans régulateur de vitesse,
ANDES garde le **couple** mécanique constant, alors que ce paquet garde la
**puissance** constante ($C_m = P_m/\omega$) ; linéarisée, cette dernière
apporte un amortissement $D = P_m$, qu'on donne donc au GENROU. L'écart
résiduel sur le mode inter-zones vient de la dynamique du stator et du réseau,
que ce paquet représente et qu'ANDES tient pour algébrique.

ANDES (GPL-3.0) n'est qu'un outil de vérification ; il s'installe à part :

```bash
uv sync --group validation
uv run pytest tests/test_andes.py
```
