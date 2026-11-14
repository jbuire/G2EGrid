# Automatiser : calculs longs, atelier par script, rapport d'étude

## Suivre et interrompre un calcul long

Une étude de robustesse refait load flow et spectre pour chaque variante ; une
continuation, pour chaque pas. Leurs boucles signalent leur avancement par
{mod}`g2egrid.noyau.progres`, que le moteur appelle sans rien savoir de qui écoute :

```python
from g2egrid import catalogue, progres
from g2egrid.scenarios.robustesse import etudier

def afficher(fraction, message):
    print(f"{100 * fraction:5.1f} %  {message}")

with progres.suivi(afficher):
    e = etudier(catalogue.charger("smib"), deconnexion=False)
```

```text
  0.0 %  SM1 P +10%
 50.0 %  charge totale +10%
```

Pour interrompre, le rappel lève {class}`~g2egrid.noyau.progres.Interrompu` : le calcul
s'arrête au prochain point de contrôle de sa boucle.

```python
n = 0
def arreter_au_deuxieme(fraction, message):
    global n
    n += 1
    if n == 2:
        raise progres.Interrompu("assez")

try:
    with progres.suivi(arreter_au_deuxieme):
        etudier(catalogue.charger("smib"), deconnexion=False)
except progres.Interrompu as err:
    print("interrompu :", err)
```

```text
interrompu : assez
```

Sans rappel installé, `progres.signaler` ne coûte qu'une lecture de variable.

## L'atelier par script

Tout ce que fait l'atelier passe par {mod}`g2egrid_api.analyse`, qui ne sait pas ce
qu'est HTTP : les mêmes fonctions s'appellent en Python.

```python
from g2egrid_api import analyse

t = catalogue.charger("wscc9")
pbs = analyse.valider(t.copie(Y_DER=[(r[0], 2 if r[1] == 1 else r[1], *r[2:]) for r in t.Y_DER]))
print(pbs[0]["message"], "->", pbs[0]["reparation"]["libelle"])
r = pbs[0]["reparation"]
t_repare = analyse.reparer(t, r["action"], r["args"])

print(analyse.elements_reseau(t)[18]["nom"])              # listes de choix
co = analyse.continuation_parametre(t, "SM1", "H", 2.0, 6.0, pas=5)   # onglet Limites
print(len(co["mu"]), "pas")
```

```text
aucun bus slack (type 1) dans le reseau -> Prendre SM du bus 1 comme slack
SM1
5 pas
```

Par HTTP, le serveur lancé avec `python -m g2egrid_api` exige son jeton : un
script le présente dans l'en-tête `X-G2egrid-Jeton` (le jeton est rangé dans
`~/.cache/g2egrid/atelier-<port>.jeton`). Les calculs longs passent par les
**tâches** :

```bash
JETON=$(cat ~/.cache/g2egrid/atelier-8765.jeton)
H="X-G2egrid-Jeton: $JETON"
# lancer : le corps est celui de la route synchrone (tables du réseau + options)
ID=$(curl -s -H "$H" -H 'Content-Type: application/json' \
     -d @corps_pv.json http://127.0.0.1:8765/api/taches/pv | jq -r .id)
curl -N -H "$H" http://127.0.0.1:8765/api/taches/$ID/flux       # avancement (SSE)
curl -s -X DELETE -H "$H" http://127.0.0.1:8765/api/taches/$ID  # arrêter
```

| Genre de tâche | Corps | Calcul |
|---|---|---|
| `robustesse` | tables, `frac`, `puissance`, `deconnexion`... | {doc}`robustesse` |
| `continuation` | tables, `element`, `parametre`, `debut`, `fin`, `pas`, `log`, `equilibre` | {doc}`continuation` |
| `impedance` | tables, `element`, `f_min`, `f_max`, `n` | {doc}`impedance` |
| `rapport` | tables, `positions`, `figures` (SVG), `n_modes` | ci-dessous |

Le flux SSE envoie un événement à chaque changement (`etat`, `fraction`,
`message`), le dernier avec le `resultat` ; `GET /api/taches/{id}` relit une
tâche terminée (les vingt dernières sont gardées).

## Le rapport d'étude

Un résultat de stabilité ne vaut que si l'on sait **de quoi** il est le
résultat. {func}`g2egrid_api.rapport.rapport_etude` produit un fichier HTML
autonome (aucune ressource extérieure) qui porte :

- la **version du code** : le commit git, et s'il y avait des modifications non
  enregistrées ;
- l'**empreinte** du réseau : le SHA-256 de ses tables sous forme canonique.
  Deux rapports de même empreinte portent sur le même réseau, au bit près ;
- le réseau, le load flow, le spectre, le rapport des modes critiques
  ({doc}`residus`) et les figures que l'atelier a envoyées ;
- les **tables** elles-mêmes : le bouton *Ouvrir* de l'atelier relit le rapport
  comme un réseau.

```python
from g2egrid_api.rapport import empreinte, rapport_etude

page = rapport_etude(t, analyse.construire, n_modes=3)
with open("rapport_wscc9.html", "w", encoding="utf-8") as f:
    f.write(page)
print(empreinte(t)[:16], "...")
```

Dans l'atelier : bouton **Rapport** de l'en-tête.
