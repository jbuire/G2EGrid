# Interface web

Un atelier de reseaux : on construit un reseau au glisser-deposer, et on le
fait passer par le load flow, l'analyse modale et la simulation temporelle --
les deux voies, lineaire et non lineaire.

Pour le lancer, double-cliquer `lancer_atelier.bat` a la racine du depot
(`lancer_atelier.sh` sous macOS ou Linux). Il demarre le serveur local et ouvre
le navigateur des que la page repond. S'il tourne deja, il ouvre seulement la
page. Il faut [uv](https://docs.astral.sh/uv/) ; le premier lancement
installe l'environnement. En ligne de commande :

```bash
cd g2egrid
uv run --extra api python -m g2egrid_api            # --port, --sans-navigateur
```
puis l'adresse affichee dans la console : elle porte le **jeton d'acces**
(`http://127.0.0.1:8765/?jeton=...`, voir *Securite* ci-dessous).

## Le principe

**Ce qu'on dessine EST les tables.** L'interface ne manipule aucune structure
intermediaire : elle edite directement `Y_DER`, `T_Lines_pu`, `T_Transfos`,
`T_SMs`... c'est-a-dire exactement ce que declarent les notebooks de
`reseaux_test/` et les fichiers `.m` de `matlab/G2Elin_V2/`. Dans l'onglet
Reseau, le bouton **Creation par tables** montre cette description en direct
(et **Creation graphique** ramene au schema) ; l'onglet **Code** en rend la
cellule de donnees, collable telle quelle dans un notebook.

Corollaire : **tout ce que fait l'interface se fait en code**, par les memes
fonctions. `g2egrid_api/analyse/` ne sait pas ce qu'est HTTP ; `app.py` n'est
qu'une couche de transport.

```python
from g2egrid_api import analyse
from g2egrid.catalogue import charger

t = charger("wscc9")  # les tables du notebook
analyse.loadflow(t)["bus"][:3]
analyse.modes(t, f_max=5.0)["modes"][0]  # le mode le moins amorti
analyse.detail_mode(t, 59)["forme"]  # qui bouge, et en quelle phase
analyse.temporel(t, "SM1.p_ref", 0.01, duree=1.0, voie="les_deux")["ecart"]
analyse.noms_mesures(t)  # ce qu'on peut relever sur ce reseau
analyse.balayage(
    t, {"genre": "famille", "famille": "GFM_PARAMS", "cle": "mp"}, [0.01, 0.02, 0.05]
)["pas"][0]["max_re"]
analyse.validation(t)["carre"]["rapport"]  # 4, si A est bien la bonne
analyse.robustesse(t)["modes"][0]  # le mode le plus sensible (N-1, +10 % P)
analyse.importer(
    analyse.depuis_catalogue("wscc9"),  # importer ieee14 dans wscc9
    {"catalogue": "ieee14"},
    {"mode": "transfo", "bus_importe": 1, "bus_existant": 5},
)["import"]
print(analyse.exporter(t))  # la cellule de donnees
```

## Les onglets

L'ordre des onglets est celui d'une etude : on construit, on fixe le point de
fonctionnement, on lit les modes, on fait varier, on simule.

| Onglet | Ce qu'on y fait |
|---|---|
| **Reseau** | construire le reseau, caracteristiques nominales : bus, lignes, transformateurs, machines, convertisseurs, continu, zones. Schema unifilaire (**Creation graphique**) ou tables (**Creation par tables**, ci-dessous). |
| *Creation par tables* | les **matrices de construction** (`Y_DER`, `T_Lines_pu`, `T_SMs`...), editables, **en miroir du schema** : modifier une cellule redessine le schema, poser un equipement ajoute ses lignes ; cliquer une ligne la montre sur le schema, et ce qui est selectionne dans l'onglet Reseau est surligne ici avec les lignes qui touchent le meme bus. |
| **Point de fonctionnement** | le type de chaque bus (slack, PV, PQ) et ses consignes ; **« tout le reste en PQ »** autour d'un slack choisi (les PV passent en PQ avec le reactif qu'ils fournissaient : le point ne bouge pas) ; load flow, plan de tension et **transits** sur le schema. |
| **Modes** | le **T_eig** des notebooks (Re, Im, f propre, f amortie, zeta, trois etats participants) et la carte des valeurs propres, **axes au choix** (Re / Im, Re / f, f / zeta...). Cliquer un mode ouvre sa forme modale et ses sensibilites. |
| **Variations** | faire varier un parametre et suivre les modes : on choisit un **element** (SM2, L3_5_7, la charge du bus 5...), puis l'un de ses **parametres** (colonnes de sa table, consignes P/Q/V, reglages) ; De et A reprennent leur valeur par defaut -- de 0,5 a 1,5 fois la valeur actuelle -- a chaque nouveau parametre. Lieu des racines dans le plan (Re, Im), marge, valeur critique. |
| **Robustesse** | quels modes bougent le plus quand la puissance varie ou qu'un element declenche (N-1) : plan complexe du cas de base et de ses variantes, classement des modes, detail par perturbation. |
| **Limites et impedances** | **continuation d'un parametre** (element et parametre au choix : max Re du spectre et lieu des modes que le parametre deplace, bifurcations raffinees) et **boucle mineure** d'un convertisseur (lieux caracteristiques autour de -1, echelle compressee, impedances du reste et du convertisseur, bandes non passives, verdict). |
| **Simplifications** | les hypotheses de simplification du calcul des modes : rendre des dynamiques algebriques (dx/dt = 0), par type d'equipement ou par equipement ; les modes sans puis avec. |
| **Temporel** | echelon ou scenario (avec sa **frise** : glisser un evenement change son instant), voie lineaire / non lineaire, mesures physiques, **rejeu** sur le schema colore par la tension simulee. |
| **Evenements** | les grandeurs electriques **apres** un evenement transitoire ou de consigne : tension de chaque bus, P et Q des unites, frequence ; relais de surtension et cascades en regimes permanents successifs. |
| **Transport / distribution** | ce que la **dynamique de la distribution** change au transport : modes deplaces, part de la distribution dans chaque mode, modes crees ou degrades, et poste par poste. **Ajouter des GFL au hasard** dans la distribution ; **modele complet ou modeles reduits** (transport seul avec ses distributions en charges R-L ou hybrides, chaque distribution seule derriere un bus infini et son impedance de Thevenin) ; **balayage de penetration**. Reseaux a zones : `wscc9_cigre6`, `ieee14_ieee33`, `wscc3sm_zones`. Voir le guide *Transport et distribution*. |
| **Equations** | la **liste globale** des equations dynamiques `dx/dt = f(x, u)` de tous les equipements, **modifiables sur place** (filtres par equipement et par genre, recherche, LaTeX). Modifier l'equation d'un equipement natif le remplace par son modele a equations ; « Revenir au modele natif » le rend. La meme liste, restreinte a un equipement et avec ses parametres, est dans sa fiche et ses proprietes. Voir [`docs/guide/equations_modifiables.md`](../docs/guide/equations_modifiables.md). |
| **Blocs** | l'editeur de **schemas-blocs** : palette de blocs elementaires (somme, gain, integrateur, premier ordre, PI, PID...), liens, **sous-systemes** (regrouper, entrer, remonter), compilation en equations, bibliotheque de modeles partagee entre les reseaux. « Editer en blocs » ouvre un equipement du reseau (GFL, GFM, machine...) en blocs ; « Remplacer par son modele a blocs » le met a sa place. Voir [`docs/guide/schemas_blocs.md`](../docs/guide/schemas_blocs.md). |
| **Code** | le script Python qui reprend le reseau (un `Atelier`), et les exports notebook / MATLAB / **MATPOWER**. |
| **Documentation** | cette documentation (Sphinx), servie par l'atelier ; le bouton *Voir dans la documentation* des equations d'un modele y mene. Si elle n'est pas construite, un bouton la construit. |

## Les cartes des modes

Pour un reseau a continu, chaque mode a un **domaine** (`g2egrid.domaines`) :
AC, DC, ou AC/DC couple (part continue entre 10 et 90 %), et ses iles ; un
mode porte par deux iles que le continu relie est *inter-iles*. Le bilan
s'affiche avec la marge ; le filtre *Domaine* ne garde qu'un cote ; sur la
carte, rond = AC, carre = DC, losange = AC/DC ; la table `T_eig` gagne les
colonnes *Dom.* (`*` : inter-iles) et *DC %*. La carte **Impact du continu
sur les modes** compare au reseau dont le continu est rigide (tensions
figees, cables instantanes) : modes alternatifs deplaces, modes du continu,
modes perdus (`POST /api/impact_continu`).

Les cartes des onglets Modes, Variations et Robustesse partagent les memes
regles : **echelle lineaire**, abscisse et ordonnee en **traits epais**, axes
Re = 0 et Im = 0 traces, et **tous les modes** affiches -- c'est l'utilisateur
qui zoome. Au-dessus de chaque carte, **X min, X max, Y min, Y max** fixent
l'echelle (une case vide reste automatique) ; tirer un rectangle sur la carte
les remplit, double-cliquer ou *Tout voir* revient a l'automatique. Variations
et Robustesse tracent le plan complexe (Re en abscisse, Im en ordonnee), avec
les rayons d'iso-amortissement a 5 et 10 %.

## Construire : palette, parametres, copier-coller

La palette de l'onglet Reseau est **une seule liste**, alternatif et continu
ensemble, rangee par role :

- **Reseau** : bus, bus continu, ligne, transformateur, cable continu, reseau
  amont ;
- **Production** : machine synchrone ; **GFM batterie**, **GFM
  photovoltaique**, **GFM sur bus DC** ; **GFL photovoltaique**, **GFL
  batterie**, **GFL sur bus DC** -- chacun pose avec sa chaine continue et
  porte sur son symbole le pictogramme de sa source (voir *Dessiner du
  continu*) ; **DC/DC** (deux cas : *batterie*, qui soutient la tension par
  un statisme, ou *PV*, qui tient sa puissance) ;
- **Consommation et compensation** : charge, charge continue (puissance ou
  resistance), inductance, condensateur, condensateur continu ;
- **Raccordement** : *raccorder un GFM/GFL a un bus continu*.

Apres chaque pose, l'outil revient a **Selection**. En Selection, **maintenir
le clic dans le vide** trace un rectangle : les bus qu'il prend (alternatifs et
continus) forment un groupe, cercle sur le schema, et glisser l'un d'eux
deplace tout le groupe (qui reste dans le cadre). Echap ou un clic dans le vide
defait le groupe.

**Convention : un GFM est a batterie, un GFL est un PV.** Un GFM tient sa
tension : son bus est pose en **PV** (ou slack s'il est le premier formeur
de son ile), comme une machine synchrone. Un GFL produit a sa puissance MPPT
sans tenir la tension : son bus est pose en **PQ**, et il ne devient jamais
le slack. La verification du reseau le signale (avis) quand un GFM est en PQ
ou un GFL en PV. Un convertisseur raccorde au continu suit, lui, le mode de
`T_VSCs`.

Dans les **proprietes** d'une unite :

- **Dynamiques** : une case par groupe de dynamiques (statisme du 2e ordre,
  bus continu, filtre de mesure, boucles de tension et de courant, PLL,
  stator, amortisseurs, AVR, PSS...). Decocher un groupe le supprime
  (`dx/dt = 0` ou `x = x0`, voir le guide *Supprimer une dynamique*) ; pour
  les GFM et GFL, les cas « statisme » et « bus continu » font passer l'unite
  d'une variante a l'autre (GFM1 a GFM4, GFL1/GFL2).
- **Parametres** : toutes les valeurs de la famille, pre-remplies de leurs
  valeurs par defaut ; une valeur modifiee s'affiche en gras et ne vaut que
  pour ce bus. Le reglage des boucles par intention et les donnees
  constructeur y sont regroupes.

**Copier-coller** : selectionner une machine, *Copier* (Ctrl+C), puis *Coller*
(Ctrl+V) et cliquer un bus. La machine arrive avec son transformateur, ses
parametres, ses dynamiques et sa consigne (P, V).

**Importer un reseau** : le bouton *Importer un reseau* (palette, Edition)
ouvre un panneau : on choisit un reseau du catalogue ou un fichier, et la
facon de le raccorder -- par un **transformateur**, une **ligne**, la
**fusion** de deux bus de meme tension, ou *aucun* pour tracer soi-meme les
branches. Le bus du reseau importe est pre-choisi (son bus infini, sinon son
slack) ; le bus existant se choisit dans la liste ou d'un clic sur le schema.
Le reseau importe se pose a droite du schema, avec (`g2egrid.zones.importer_reseau`) :

- une **renumerotation** par un prefixe, la centaine au-dessus des numeros
  deja pris (zones instanciees comprises) : ieee14 importe dans wscc9 devient
  les bus 101 a 114 ;
- **une seule reference**, celle du reseau existant : ses bus infinis sont
  retires (le reseau existant devient son amont) et ses slacks passent en PV,
  avec la production que leur donne le load flow du reseau importe seul ;
  deux cases gardent l'un ou l'autre ;
- ses **reglages de famille** reportes bus par bus, pour ne pas se meler a
  ceux de l'existant ;
- les shunts convertis si `Sb` differe ; une autre frequence prend celle du
  reseau existant (reactances et pu conservees, message a l'appui).

Une fusion est refusee si les tensions nominales different, ou si le bus
porterait deux unites. Sans raccordement, le reseau importe reste une ile
sans reference, signalee par la verification jusqu'a ce qu'on le raccorde.

**Definir une zone** : l'outil *Definir une zone* trace un rectangle ; les bus
qu'il contient deviennent une zone de distribution, instanciee a la meme
place (le reseau aplati ne change pas). La selection doit se raccorder au
reste par une seule branche et ne contenir ni le slack ni de convertisseur
raccorde au continu.

## Robustesse des modes

Un mode lu en un seul point de fonctionnement ne dit pas s'il est robuste. L'onglet
**Robustesse** (et `g2egrid.robustesse.etudier`) refait l'analyse modale sur des variantes du
reseau :

- **puissance** : chaque unite non slack pousse sa production de 10 %,
  chaque charge monte de 10 %, puis toute la charge ensemble ; une charge qui
  monte est reprise par les unites reglantes -- machines synchrones et GFM a
  batterie, au prorata de leur puissance nominale ; les GFL (PV) restent a
  leur puissance MPPT ;
- **deconnexions** : chaque element manoeuvrable est ouvert a son tour (N-1) ;
  la puissance perdue est reprise de meme. Une ouverture qui coupe le reseau
  en plusieurs iles vivantes est un *ilotage* : rendue, mais hors score.

Les modes de chaque variante sont apparies a ceux du cas de base par
affectation optimale sur la ressemblance de leurs **participations** (sur les
etats communs) et la distance des valeurs propres. Un mode porte a plus de
50 % par l'element ouvert **disparait** avec lui. Le score d'un mode est son
plus fort deplacement relatif `|d lambda| / |lambda|`, par genre de
perturbation ; *zeta min* est l'amortissement le plus bas atteint.

**Lire l'onglet.** La carte est le plan complexe (Re en abscisse, Im en
ordonnee, echelle lineaire, tous les modes) : un rond noir par mode du cas de
base, et autour le meme mode dans chaque variante -- rond bleu pour une
variation de puissance, carre orange pour une deconnexion. Cliquer un mode
relie ses variantes a lui : une variante qui part vers la droite le
desamortit, et passer l'axe Re = 0, c'est devenir instable ; les rayons
pointilles sont les iso-amortissements a 5 et 10 %. Un nuage serre dit un
mode robuste ; une variante isolee designe l'element dont le mode depend.
Pour les modes electromecaniques, zoomer pres de l'origine (rectangle, ou
bornes X/Y). Le tableau classe les modes du plus sensible au moins sensible,
et le detail donne chaque perturbation (survol : les memes chiffres sur la
carte).

Ce qui en ressort sur le catalogue (bande 0,05-5 Hz, zeta <= 30 %) :

- **+10 % de puissance deplace peu les modes electromecaniques** : score
  inferieur a 0,02, amortissement a moins de 1,5 point pres, sur `wscc9`,
  `kundur2area`, `ieee14`, `ieee39`.
- **Les deconnexions dominent** : la perte d'une machine recompose les modes
  (jusqu'a 0,7 Hz et 40 points d'amortissement), et l'ouverture d'une ligne
  peut desamortir fortement : sur `ieee14`, ouvrir `L1_1_2` fait passer le
  mode a 2,03 Hz de 13,9 % a 7,1 %.
- Des variantes n'ont **pas de point de fonctionnement** : sur
  `kundur2area`, perdre SM3 ou SM4 (zone 2, importatrice) fait diverger le
  load flow. Sur `wscc9_gfm`, quand les GFM etaient declares en PQ a Q = 0,
  aucune ligne ne pouvait s'ouvrir sans effondrement de tension ; ils sont
  desormais en PV (GFM a batterie) et tout le N-1 passe.
- Sur `wscc9_cigre`, ouvrir le transformateur d'un poste de distribution
  rend un **GFL** instable (GFL3, ou GFL113 dans la zone) : ce sont les
  modes des convertisseurs, pas les electromecaniques, qui sont fragiles au
  N-1.

## Du constructeur au code : l'Atelier

Tout ce que font les onglets se fait en Python, par les **memes fonctions** --
un T_eig calcule dans un notebook est celui de l'onglet Modes, au bit pres ::

    from g2egrid_api.atelier import Atelier

    a = Atelier.catalogue("mon_reseau")    # enregistre depuis l'interface ("Au catalogue")
    a.un_slack_tout_pq(slack_bus=1)
    a.loadflow(); a.transits()             # DataFrames
    T = a.T_eig()
    ss = a.etat()                          # A, B, C, D
    a.balayage({"genre": "famille", "famille": "GFM_PARAMS", "cle": "mp"}, [0.01, 0.02])
    a.impact()                             # reseau a zones

Et dans l'autre sens : construire en code, avec les conventions de la palette,
puis ouvrir le reseau dans l'interface ::

    a = Atelier.vide("demo")
    b1 = a.bus(230e3); b2 = a.bus(230e3, Pl_W=80e6, Ql_var=20e6)
    a.ligne(b1, b2, R_pu=0.01, X_pu=0.1)
    a.equipement("IB", b1, slack=True)
    a.equipement("GFM1", b2, Sn_VA=100e6, P_W=50e6)
    a.enregistrer("demo")                  # visible aussitot dans le catalogue

L'onglet **Code** produit le script correspondant au reseau ouvert.

## Etats bruts et mesures

Ce que rend le modele est un vecteur d'**etats** : des flux de machine, des
courants de branche, des tensions dq de noeud. Ce qu'on compare a un
enregistrement, c'est un module de tension, une frequence, un transit. Les
**mesures** sont cette seconde couche, reconstruite le long de la
trajectoire : `V_busN`, `angle_busN`, `f_busN` et `fm_busN` (filtree sur une
periode), `P_<element>` et `Q_<element>`, et les quatre transits d'une
branche.

Elles ne sont rendues que sur la voie **non lineaire** : elles se
reconstruisent a partir des `(x, z, u, y)` re-resolus a chaque instant, ce
qu'une exponentielle de matrice sur les etats ne fournit pas.

De meme, **la voie lineaire ne suit qu'un echelon unique** : sa reponse est
calculee en forme close. Des que le scenario compte plusieurs evenements, ou
une rampe, elle est omise -- et l'interface le dit -- plutot que rendue
fausse.

Attention a la consigne perturbee : `w_ref` vient en tete de la liste des
entrees parce que c'est l'entree **globale**, la pulsation du repere. Un
echelon de 0,01 pu dessus deplace le repere de tout le reseau d'un
demi-hertz : sur `smib` les deux voies s'ecartent alors de 386 %, ce qui ne
dit rien du modele et tout de la perturbation choisie. L'interface
preselectionne donc une vraie consigne d'equipement (0,15 % sur le meme
reseau).

## Valider le petit signal par l'EMT

(En Python : `analyse.validation` ; l'onglet correspondant a ete retire de
l'atelier.)

Trois epreuves de cout croissant, qui ne disent pas la meme chose.

1. **Point fixe.** Au point de fonctionnement, `dx/dt` doit etre nul et les
   variables algebriques doivent deja verifier leurs equations. Sinon les
   points de fonctionnement par element ne forment pas un point fixe du
   *reseau* : invisible en modal -- on linearise quand meme, autour d'un
   point qui n'existe pas -- et fatal en temporel.
2. **Le carre de l'amplitude.** L'ecart entre les deux voies doit croitre
   comme `a^2` : c'est la definition d'une linearisation au premier ordre.
   Doubler l'amplitude doit quadrupler l'ecart. Mesure : **3,995** sur
   `wscc9`, **4,012** sur `smib`. Un rapport proche de 2 dirait que `A`
   n'est pas la bonne.
3. **Identification.** On identifie les modes directement sur la
   trajectoire non lineaire et on les apparie aux valeurs propres de `A`.
   C'est la seule epreuve qui ait un sens physique : les deux chemins ne
   partagent aucune ligne de code de resolution, si bien que leur accord ne
   peut pas venir d'une erreur commune. Sur `wscc9`, les deux modes
   electromecaniques (2,2000 Hz / 8,696 % et 2,0115 Hz / 11,251 %) sont
   retrouves a mieux que 0,1 % en frequence et 0,5 % en amortissement.

Deux avertissements que l'interface affiche, parce que les taire serait
pire que de ne rien dire :

- **Le residu de reconstruction ne previent pas** des deux facons de se
  tromper qui comptent. Sur une fenetre trop courte il vaut 0,0002 %, dix
  fois *mieux* que sur une fenetre correcte, pendant que l'amortissement est
  faux de 56 % : sur une demi-periode, une exponentielle amortie et une
  exponentielle plus raide de residu plus grand sont indiscernables. Meme
  piege sur une seule voie placee au noeud d'un mode. Le remede est
  d'allonger la duree et de prendre plusieurs voies.
- **Les seuils de 2 % et 15 % visent les modes electromecaniques**, les
  moins amortis. Un mode tres amorti les depasse sans que le modele y soit
  pour rien : il decroit avant d'avoir boucle une periode et ne laisse
  presque rien dans la fenetre. Les modes amortis a plus de `ZETA_OSCILLANT`
  (95 %) sont d'ailleurs ecartes d'office -- sur `wscc9`, `modal.modes`
  rendait autrement un mode a 0,005 Hz amorti a 99,9 %, qu'aucune
  identification ne peut lire.

## Le balayage, et le mode de jauge

Le balayage repond a la question que ni une carte des modes ni une
simulation ne posent : **a partir de quelle valeur** ce reglage rend le
reseau instable. Sur `smib_gfm`, un balayage du statisme `mp` de 0,5 a 10 %
en 30 points donne une valeur critique de **7,58 %**, le mode qui traverse
etant a 25 Hz.

Chaque pas assemble le reseau et resout son equilibre : de l'ordre de
0,3 s/pas sur `wscc9_gfm`, 1,8 s/pas sur `wscc3sm_cigre3` (54 bus). Le
nombre de valeurs est donc borne a 60 ; au-dela, c'est un script qu'il faut.

Les modes sont **suivis d'un pas au suivant** par un appariement global
(`linear_sum_assignment`) : sans lui, un lieu des racines n'est qu'un nuage,
`eig` ne promettant aucun ordre stable. Le controle est qu'une paire
conjuguee traverse l'axe a la **meme** valeur, ce que le test verifie.

Un pas qui ne converge pas est garde comme resultat, avec son message :
c'est souvent la fin du domaine de validite, et l'effacer ferait croire que
le balayage s'arrete sans raison.

**Le mode de jauge est ecarte.** Un reseau sans bus infini n'a pas de
reference d'angle absolue : son spectre porte un zero exact, rendu a 1e-12
pres. Son amortissement vaut alors `-Re/|lambda|`, du bruit d'arrondi divise
par du bruit d'arrondi -- soit +100 ou -100 % au hasard du point de
fonctionnement. Il arrivait donc en tete du classement par amortissement,
occupait un lieu a lui seul, et faisait afficher `max Re = +0.000` sur un
reseau parfaitement stable (`wscc9_gfm`, ou le vrai maximum est -2,3). Le
seuil `JAUGE` est pris dans le fosse : la jauge y est a 4e-12, le mode
physique le plus lent a 0,19 rad/s.

## Dessiner du continu

La palette est **une seule liste**, alternatif et continu ensemble, rangee
par role : *Reseau* (bus, **bus continu** -- un carre, pour qu'on ne le
confonde pas d'un coup d'oeil avec un noeud alternatif --, ligne,
transformateur, **cable continu**, reseau amont), *Production* (machine
synchrone, GFM et GFL avec leur source, **DC/DC**), *Consommation et
compensation* (charge, **charge continue**, inductance, condensateur,
**condensateur continu**) et *Raccordement* (**raccorder un convertisseur** a
un bus continu).

### GFM et GFL avec leur source

Un GFM ou un GFL se choisit **avec sa source d'energie** :

| Outil | Pose | Commande accordee |
|---|---|---|
| GFM batterie | GFM3, bus continu + condensateur, DC/DC batterie | `PV` (statisme P-f), la batterie tient la tension continue |
| GFM photovoltaique | GFM3, bus continu + condensateur, DC/DC PV | `VdcV` |
| GFM sur bus DC | GFM3, bus continu + condensateur | `VdcV` tant qu'il est seul |
| GFL photovoltaique | GFL2, bus continu + condensateur, DC/DC PV | `VdcQ` |
| GFL batterie | GFL2, bus continu + condensateur, DC/DC batterie | `PQ`, la batterie tient la tension continue |
| GFL sur bus DC | GFL2, bus continu + condensateur | `VdcQ` tant qu'il est seul |

Cliquer un bus y raccorde le convertisseur par un **transformateur** ;
cliquer dans le vide le pose seul, son bus alternatif a relier par une ligne
ou un transfo (le convertisseur « a bornes ouvertes », ci-dessous). Le
DC/DC fait la moitie de `Sb` et injecte 10 % de `Sb` ; derriere une
batterie, le convertisseur regle 10 % de `Sb`.

Le symbole porte le **pictogramme de sa source** dans une pastille en haut
a droite -- une pile verte, un panneau sous un soleil, ou le carre violet du
bus continu --, et l'infobulle la nomme. La source se lit dans l'ile
continue du convertisseur : une batterie s'il y en a une, sinon un PV, sinon
le bus continu ; passer un DC/DC de PV a batterie dans ses proprietes change
donc le pictogramme. Un convertisseur a bus continu **ideal** (GFM1, GFM2,
GFL1, ceux du catalogue) suit la convention de l'editeur : un GFM est a
batterie, un GFL photovoltaique.

Une exception : un **GFM batterie qui serait le premier formeur de son ile**
doit en tenir la reference d'angle, ce qu'un convertisseur raccorde au
continu ne peut pas faire (regle ci-dessous). Il se pose alors a bus continu
ideal (GFM1, la batterie par convention), en slack. Les autres variantes,
posees dans une ile sans reference, le signalent dans la barre d'etat.

### Raccorder a la main

Ce dernier est le seul qui demande une explication. Une ligne de `T_VSCs`
n'est pas une branche : c'est une **declaration**. Elle dit que le GFM3 (ou
GFM4, ou GFL2) pose sur tel bus alternatif devient sa variante a port
continu, dont le condensateur n'est plus un etat interne mais le bus continu
lui-meme. On clique donc le bus alternatif du convertisseur, **puis** son bus
continu -- l'ordre de lecture de la table.

Le mode se choisit tout seul, **par ile continue** : une ile continue a
besoin d'exactement une reference de tension. Un convertisseur la tient
(`VdcV` pour un grid-forming, `VdcQ` pour un grid-following) et les autres
tiennent leur puissance -- sauf quand l'ile ne compte qu'un convertisseur et
une **batterie** : la batterie tient alors la tension continue et le
convertisseur regle sa puissance (`PV` pour un GFM, statisme P-f ; `PQ` pour
un GFL). Un DC/DC **PV**, a puissance constante, ne tient jamais la tension.
Le bus continu de la reference passe en type 1, les autres en type 2.

### GFM / GFL a bornes ouvertes

Les outils **GFM sur bus DC** et **GFL sur bus DC**, d'un clic dans le
vide, posent un convertisseur **avec ses deux cotes ouverts** : son bus
alternatif (PV pour un GFM, PQ pour un GFL -- jamais le slack) et son bus
continu, deja declares dans `T_VSCs`, rien de branche ni d'un cote ni de
l'autre. On branche ensuite a la main :

- cote alternatif, une **ligne** ou un **transformateur** vers le reseau ;
- cote continu, un **DC/DC batterie**, un **DC/DC PV** ou un **cable
  continu** vers un autre bus continu (un autre convertisseur, un reseau
  continu existant).

La commande suit ce qu'on branche (regle ci-dessus) ; passer un DC/DC de PV
a batterie dans ses proprietes la refait, et la barre d'etat dit ce qui a
change. Le bus continu recoit d'office un **condensateur de 0,2 pu** : sans
lui, la boucle de courant du DC/DC et la regulation du GFM se renvoient la
tension du bus, et ce mode part instable (Re = +57 1/s pour un GFM a
batterie qui exporte 10 MW). Les combinaisons choisies par l'editeur sont
verifiees stables par `tests/test_dc_ouvert.py`, et les six outils de la
palette, posees sur un reseau, par `tests/e2e/test_navigateur.py`.

Deux regles s'appliquent, que l'interface fait respecter pendant qu'on
dessine plutot que de les laisser decouvrir au premier load flow :

- **un convertisseur ne peut pas etre sur le slack** -- sauf s'il FORME seul
  son ile (mode `Vtheta`, la station d'un parc offshore). Il tient
  exactement deux grandeurs, et sur la reference d'angle les deux puissances
  sont libres. C'est la regle qui surprend le plus, parce que le premier
  equipement pose devient justement la reference ;
- **la reference d'angle se compte par ILE.** Une liaison continue transmet
  de la puissance mais aucun angle : chaque cote garde la sienne. Poser un
  equipement dans une ile qui n'en a pas encore l'y installe.

Les bus alternatifs et continus **partagent l'espace de numeros** -- c'est le
contrat des tables, et le cablage s'appuie dessus pour retrouver un noeud
sans avoir a connaitre son genre.

### Statisme Vdc-P

Le mode d'une station se choisit dans ses proprietes (ligne de `T_VSCs`),
parmi ceux que sa classe admet : `VdcQ`, `PQ`, `droop` pour un GFL a port
continu ; `PV`, `VdcV`, `Vtheta` pour un GFM. Le choix met les tables en
coherence aussitot, et la barre d'etat dit ce qui a change :

| Mode | Ce que l'atelier ajuste |
|---|---|
| `droop` | le noeud continu passe en type 2 (sa `V_pu` devient la consigne) ; un statisme de 5 % s'il n'y en a pas |
| `VdcQ`, `VdcV` | le noeud continu passe en type 1 (tenu par la station) |
| `Vtheta` | le bus alternatif devient le slack de son ile ; en quittant ce mode, il redevient PV (GFM) ou PQ (GFL) |

Une station en `droop` montre, sous ses champs, le bloc **Statisme Vdc-P** :
`k` (en % : pu de tension continue par pu de puissance, sur la puissance Sn de
la station ; rangé dans `PARAMS_PAR_BUS[bus].kdroop`), `P_ref` (le `Pg` de son
bus) et `Vdc_ref` (la `V_pu` de son noeud continu), et sa caracteristique
`P = P_ref + (Sn/k).(Vdc - Vdc_ref)`, avec le point du dernier load flow. Le
bloc apparait aussi quand on selectionne le convertisseur lui-meme.

Une ile continue peut avoir plusieurs references de tension a des noeuds
distincts, et plusieurs statismes : l'accord automatique des commandes ne
retire plus une reference posee (il n'en installe une que s'il n'y en a
aucune, ni statisme, ni batterie).

**Balayer les statismes** : onglet *Limites et impedances*, carte
*Statismes Vdc-P*. Le `k` de toutes les stations en statisme (ou d'une
seule) est multiplie par un facteur, en echelle logarithmique ; a chaque pas,
load flow, equilibre et spectre. Trois traces : la puissance de chaque
convertisseur a port continu (le partage que fixent les statismes ; traits
pleins pour les stations en statisme), les tensions continues, et la plus
grande partie reelle du spectre, avec les bifurcations eventuelles. Le calcul
passe par les taches de fond (`POST /api/taches/statismes`).

## Lire une forme modale

Le diagramme polaire place chaque etat participant par son **module relatif**
au plus grand et sa **phase relative a l'etat le plus participant**. Deux
vitesses rotor en phase : les machines bougent ensemble, mode local. A 180
degres : elles s'opposent, mode inter-machines. Sur `wscc9`, le mode a 2,2 Hz
montre `SM3.Dwr` a 0 degre et `SM2.Dwr` a 172 degres -- SM2 et SM3 oscillent
l'une contre l'autre.

Le selecteur « Etats » restreint a un suffixe (`Dwr`, `theta`, ...) : une
vitesse en pu et une puissance en pu n'ont pas le meme ordre de grandeur, et
tout tracer ensemble ecrase le diagramme.

La table de sensibilite se lit comme le modele : `A[i, j]` est l'effet de
l'etat **colonne** sur la derivee de l'etat **ligne**. C'est un coefficient
d'equation -- un statisme, un gain de boucle -- donc un endroit ou agir, pas
une paire d'etats symetrique.

## Les zones

Un reseau de transport auquel se raccordent plusieurs distributions
identiques se decrit **par zones** : la distribution une fois, en numerotation
locale et sans slack, puis **instanciee** a chaque bus du transport ou elle se
raccorde -- chaque instance pouvant valoir K departs identiques.

Ouvrir un tel reseau (`wscc3sm_zones`) montre le **transport**, ou chaque
instance est un bloc en pointille qui dit sa zone, son prefixe de
numerotation et son K. **Double-cliquer un bloc ouvre la zone** : le schema
montre alors ses bus locaux, son point de raccordement (PCC) cerne de
pointilles, et le fil d'Ariane rappelle a combien d'instances la zone sert.
**Une modification dans une zone vaut pour toutes ses instances**, chacune a
son echelle K : doubler une charge de la zone la double a chacun des trois
postes, fois 4, 3 et 3.

Les calculs portent toujours sur le reseau **aplati**, jamais sur la vue
ouverte : lances depuis une zone de 15 bus, le load flow resout les 54 bus du
reseau et l'analyse modale ses 429 etats. Le serveur re-aplatit apres chaque
edition, et un calcul attend que ce soit fait.

Un clic simple sur un bloc selectionne l'**instance** : on y regle son bus
amont, son prefixe et son K. La zone elle-meme ne s'edite qu'en l'ouvrant,
puisqu'elle vaut pour toutes.

## Disjoncteurs, ilotage et manoeuvres

**En regime permanent**, chaque ligne, transfo, unite, charge et shunt a son
disjoncteur : la case *en service* de l'inspecteur l'ouvre. L'etat est une
liste `OUVERTS = [(table, ID), ...]` rangee avec les reglages ; les tables
elles-memes ne bougent pas, refermer rend le reseau d'origine. L'assembleur
(`g2egrid.disjoncteurs.en_service`) retire ce qui est ouvert, met HORS TENSION
une ile qui n'a plus de formeur de reseau (machine, grid-forming, bus infini
-- un grid-following ne fait que suivre), et donne un slack a chaque ile
vivante (son slack s'il est encore la, sinon son plus gros formeur). Le
schema montre les disjoncteurs ouverts (carres rouges evides, pointilles) et
estompe ce qui est hors tension ; la validation dit ce qui s'est passe.

**En simulation**, une *manoeuvre* du scenario de l'onglet Temporel change le
reseau a son instant : ouverture d'un element, variation d'une charge
(impedance changee a la tension d'equilibre), saut de phase d'un bus. Rien
n'est recalcule par un load flow : les elements gardent leurs consignes,
leurs etats sont repris par leur nom (`Manoeuvre`,
`ModeleNonLineaire.appliquer`). Declencher la machine de reference passe le
repere a un bloc a part (`Repere`) qui suit la plus grosse machine restante
-- l'idee du `frame.py` de F. Kelada. Les grandeurs d'un element retire
valent `null` apres son retrait.

```python
from g2egrid.temporel import Manoeuvre, ModeleNonLineaire

m = ModeleNonLineaire(net)
tr = m.simuler((0, 1), uext_fn=m.scenario(Manoeuvre(0.1, "ouverture", "SM2")))
```

## Defauts, limiteurs, integrateurs

**Defaut** : une manoeuvre `defaut` pose un court-circuit R-L au bus choisi
(`R = 1e-3`, `X = 1e-2` pu par defaut), `elimination` le retire.

**Defaut continu** : `defaut_dc` pose un court-circuit pole a pole au noeud
continu choisi (`R = 1e-4`, `X = 1e-3` pu par defaut), `elimination_dc` le
retire ; l'`ouverture` d'un cable continu est l'action d'un disjoncteur a
courant continu ; `blocage` bloque un convertisseur (limitation de courant a
zero, decouple de son noeud continu). La case **protections DC** arme le
plan type de `g2egrid.temporel.protections_dc` -- disjoncteur sur chaque
cable, blocage des stations sur sous-tension continue ou surintensite,
decrochage des parcs sur surtension alternative -- surveille pendant la
simulation ; les declenchements s'affichent sous les courbes. Voir
`docs/guide/protections_dc.md`.

**Charges sous 0,7 pu** : une charge a puissance ou courant constant devient
une impedance (mode discret `mode_bt`), sans quoi elle perdrait toute
solution a tension nulle.

**Limiteurs** (`g2egrid.elements.limiteurs`) : chaque GFM et GFL limite sa
consigne de courant a `Imax` (1,2 pu de l'equipement), les integrateurs de sa
boucle externe la suivent par retro-calcul (`Taw`) ; l'amplificateur de l'AVR
des machines bute a `VRmax` / `VRmin`. Ce sont des MODES : un drapeau qui
bascule a l'instant ou une garde franchit zero, de sorte que les equations
restent lisses entre deux basculements -- et qu'a drapeau nul elles sont
exactement celles d'avant (petit signal et parite MATLAB inchanges). Les
basculements sont rendus avec la simulation et affiches sous les courbes. La
case *limiteurs* les desactive, pour comparer.

**Integrateurs** : `Radau` (pas adaptatif, ordre 5) pour la validation fine du
petit signal ; `BDF2` (`g2egrid.integrateur`) pour les defauts et manoeuvres :
pas fixe, etats et contraintes resolus ENSEMBLE a chaque pas (jacobien creux a
structure fixe, LU gardee), les resonances plus rapides que le pas amorties au
lieu d'etre suivies. *auto* choisit BDF2 des qu'il y a une manoeuvre. Sur un
defaut franc, Radau suit une resonance a ~20 kHz et n'avance plus ; BDF2 a
100 us simule 0,6 s de wscc9_gfm en ~20 s.

**Equilibre** : `net.solve_equilibrium(pseudo_transitoire=True)` recale le
point de chaque element sur un equilibre du reseau entier (continuation
pseudo-transitoire, l'idee de VeraGrid), quand l'initialisation explicite n'en
donne pas un.

## Regler par intention

Dans l'inspecteur d'un **GFM** ou d'un **GFL**, *Reglage des boucles* donne
chaque boucle (tension, courant, PLL, puissance) par son temps de reponse et
son amortissement ; les gains qui en sortent sont LUS sur le modele
(placement de poles, `wn = 3/(zeta.tau)`), et `g2egrid.reglage.
intention_depuis_gains` fait le chemin inverse pour un jeu arrive avec des
gains. Les `tau_*` et `zeta_*` se balaient comme tout reglage de famille.

Pour une **machine synchrone**, *Donnees constructeur* saisit la fiche
(Xd, X'd, X''d, T'd0, T''d0 et l'axe q) et la convertit en parametres du
modele de Park (`g2egrid.reglage.machine_depuis_standard`, definitions
classiques de Kundur ; le jeu par defaut redonne son exemple 4.1).

## Langue de l'interface (FR / EN)

Le bouton **FR / EN** de l'en-tete bascule toute l'interface entre le francais
et l'anglais, sans recharger la page ni perdre le reseau en cours. La langue de
depart est celle du navigateur (francais s'il est en francais, anglais sinon) ;
le choix fait au bouton est retenu, et `?lang=fr` ou `?lang=en` dans l'adresse
l'emporte.

L'atelier est ecrit en francais : `web/i18n.js` traduit ce que la page affiche
(textes, infobulles, boites de confirmation) a partir du dictionnaire
`web/i18n_en.js`, qui va du texte francais a l'anglais. Une entree a trous
(`"bus {0} a {1} s"`) vaut pour tous les nombres et tous les noms. Un texte
absent du dictionnaire reste en francais : `I18N.manquants()`, dans la console
du navigateur, liste ceux que la page a montres -- c'est la liste a completer
quand on ajoute un texte a l'interface.

Restent en francais : la documentation (onglet Documentation), le rapport
d'etude, le code exporte, et une partie des messages d'erreur du serveur (leur
prefixe est traduit, le detail peut ne pas l'etre). Les nombres, les noms
d'equipements et les donnees saisies ne sont jamais modifies.

## Securite

Le serveur n'ecoute que la machine locale, et `python -m g2egrid_api` tire a
chaque lancement un **jeton** aleatoire. L'adresse ouverte le porte ; la page
d'accueil l'echange contre un cookie `HttpOnly`, `SameSite=Strict`, puis
redirige vers `/` pour l'oter de la barre d'adresse. Aucune autre page du
navigateur, qui ne le connait pas, ne peut appeler l'atelier ; l'en-tete
`Host` doit designer la machine locale (rebinding DNS). Un script presente le
jeton dans l'en-tete `X-G2egrid-Jeton` ; il est range, lisible du seul
utilisateur, dans `~/.cache/g2egrid/atelier-<port>.jeton`, ce qui permet a un
second double-clic de rouvrir la page du serveur deja lance. `--sans-jeton`
retire la garde (conteneur isole seulement). Sans la variable `G2EGRID_JETON`
-- `uvicorn g2egrid_api.app:app` lance a la main, tests --, rien n'est exige.

## Annuler, reparer

**Annuler / retablir** (Ctrl+Z, Ctrl+Y ou Ctrl+Maj+Z, et les fleches de
l'en-tete) : l'historique garde des instantanes complets du reseau, pris a
chaque edition ; les editions rapprochees (une saisie) se fondent, et un
reseau ouvert est une entree a lui seul.

**Reparer.** Un probleme de validation qui a une correction sure porte un
bouton : prendre le plus gros formeur de l'ile comme slack, n'en garder qu'un,
passer un bus en PV ou en PQ, raccorder un bus isole au bus de meme tension le
plus proche sur le schema, mettre une charge a Ql <= 0 a puissance constante,
ajouter un condensateur a un noeud continu qui n'en a pas. La correction
s'annule par Ctrl+Z. En code : `analyse.valider(t)[k]["reparation"]`, puis
`analyse.reparer(t, action, args)`.

**Convertisseurs hors limite.** Sur un reseau sans erreur qui porte des
convertisseurs, la validation calcule aussi le point de fonctionnement (load
flow et equilibre, gardes en cache pour les calculs qui suivent) et rend un
**avis** par convertisseur dont le courant cote convertisseur depasse `Imax`
(ou la modulation `Mmax`, si elle est reglee) : le modele petit signal, en
mode libre, ne le represente pas. Le detail chiffre est sous
`analyse.valider(t)[k]["depassement"]`. Voir `g2egrid.limites`.

## Le mode sur le schema

Ouvrir un mode (onglet Modes) affiche, sous la carte, le **schema unifilaire
teinte par sa participation** : chaque bus prend la participation cumulee des
elements qui y sont raccordes, chaque branche la sienne. Une **fleche** porte la
phase de l'etat dominant du bus, relative a celle de tout le mode ; animees,
toutes les fleches tournent ensemble (un tour en 4 s) : deux fleches opposees,
ce sont deux groupes qui oscillent l'un contre l'autre.

## Calculs longs et rapport

Les calculs longs -- robustesse, continuation, boucle mineure,
rapport -- partent en **tache de fond** : une barre d'avancement s'affiche dans
l'en-tete, alimentee par le moteur lui-meme (`g2egrid.progres`), et le bouton
**Arreter** les interrompt au prochain point de controle.

Le bouton **Rapport** telecharge un **rapport d'etude HTML autonome** : version
du code (commit git, modifie ou non), empreinte SHA-256 des tables, reseau,
load flow, spectre, rapport des modes critiques et figures deja calculees. Il
contient les tables : **Ouvrir** le relit comme un reseau.

## Confort

**Brouillon automatique.** Chaque edition range le reseau, zones comprises,
dans le stockage local du navigateur. Rien ne part au serveur. La page le
rouvre au demarrage, avec un lien *Repartir du catalogue*. Ce n'est pas une
sauvegarde : seuls le fichier `.g2egrid.json` et le catalogue durent. Le
stockage d'un navigateur peut etre vide, plein ou interdit, et l'atelier
marche sans lui.

**Lieu des racines parcouru.** Sous le lieu d'un balayage, un curseur (ou
&#9654;) parcourt les valeurs une a une. Le tableau donne, a chaque pas, les
**participations de chaque mode regroupees par element** et l'etat dominant.
Le trait montre ou passe un mode ; le tableau montre qui oscille. On y voit
par exemple un mode de machine devenir un mode de convertisseur a l'approche
de l'instabilite, alors que le trait reste continu. En code :
`analyse.balayage(...)["pas"][k]["modes"][j]["part"]`.

**Equations du modele.** L'inspecteur de chaque element ouvre les equations
de sa classe (`dx/dt = f`, `0 = g`, `y = h`), avec un bouton *Copier en
LaTeX*. Elles sont ecrites par SymPy a partir des memes expressions que
celles qu'il compile : ce qui s'affiche est ce qui se calcule. La page
[`docs/equations_modeles.md`](../docs/equations_modeles.md) les reunit toutes.
Un test verifie qu'elle est a jour ; pour la regenerer :

```bash
uv run python -m g2egrid.equations > docs/equations_modeles.md
```

## Ouvrir, enregistrer

- **Enregistrer** telecharge un vrai `.g2egrid.json`, celui que le catalogue
  relit : on peut le deposer dans `g2egrid/reseaux/`.
- **Au catalogue** l'ecrit cote serveur, dans `g2egrid/reseaux/sauvegardes/`,
  et il apparait aussitot dans la liste. **Pas dans `g2egrid/reseaux/`** : ce
  dossier porte les exemples livres avec l'outil, suivis par git et couverts
  par les tests -- une sauvegarde y cassait la non-regression.
- **Ouvrir** accepte un fichier natif, la reponse d'une lecture, ou l'ancien
  JSON a plat du constructeur. Le format est l'affaire du serveur ; un
  fichier qui ne contient aucun bus est signale, plutot qu'ouvert en silence.
- **Ouvrir** lit aussi les fichiers des autres outils : **MATPOWER** (`.m`),
  **PSS/E** (`.raw`, versions 33 et 34) et **pandapower** (son JSON, reconnu a
  sa signature), par `g2egrid.interop`. Ce qui n'a pas pu etre repris
  (dephaseur, transformateur a trois enroulements, couplage horaire...)
  s'affiche dans la barre d'etat. Un **rapport d'etude** se rouvre aussi.

## Les symboles du schema

Le schema est un **schema unifilaire** : il se lit par ses symboles, sans
legende.

| Symbole | Equipement |
|---|---|
| cercle marque d'une sinusoide | machine synchrone |
| carre marque **GFM**, variante en indice | grid-forming (GFM1 a GFM4) |
| carre marque **GFL**, variante en indice | grid-following (GFL1, GFL2) |
| carre hachure | reseau amont (bus infini) |
| deux cercles enlaces, au milieu de la branche | transformateur |
| fleche qui quitte le bus | charge -- celles de `Y_DER` comme celles de `T_Loads` |
| deux plaques et une terre | gradin de condensateurs |
| petit carre violet | bus continu |

Un bus est un simple noeud, colore par son type (slack, PV, PQ) -- ou par
son module de tension dans l'onglet Load flow. Les symboles partent du noeud
**a l'oppose de ses voisins**, la ou il y a du vide, et se repartissent en
eventail quand plusieurs equipements partagent un bus.

**Cliquer un symbole selectionne son equipement, le glisser deplace son
bus.** Les outils qui attendent un bus (ligne, transformateur, raccordement
continu) prennent le clic sur un symbole pour un clic sur son bus. Une branche
se vise par une bande de 12 px : un trait de 2 px ne se clique pas.

## Conventions du constructeur

- **Un bus nu est toujours PQ.** C'est l'equipement qui designe le slack : le
  premier qu'on pose devient la reference d'angle du reseau. Faire du premier
  bus dessine le slack donnerait un reseau dont la reference n'est portee par
  aucune source, et le modele d'etat n'aurait pas de repere.
- **Poser une machine ou un convertisseur sur un bus lui cree son bus terminal
  et son transformateur de raccordement**, comme le veulent les reseaux de
  test : une unite vit sur son propre bus, derriere un transfo. Les charges
  speciales et les gradins, eux, se posent directement sur le bus vise.
- Les reglages fins d'un equipement (gains, constantes de temps) s'editent
  dans le panneau de proprietes, sous « Reglages fins » : ils atterrissent
  dans `PARAMS_PAR_BUS`, qui l'emporte sur les reglages de famille.

## Routes

Toutes les routes de calcul prennent **les tables** en corps, jamais un
identifiant de reseau preenregistre : un reseau dessine a la main et un reseau
du catalogue suivent exactement le meme chemin.

| Route | Rend |
|---|---|
| `GET /api/reseaux` | les fiches du catalogue |
| `GET /api/reseaux/{nom}` | ses tables et un placement de ses bus |
| `GET /api/defauts` | colonnes des tables, familles d'equipement, reglages par defaut |
| `POST /api/valider` | les problemes du reseau, avec leur gravite et leur reparation quand il y en a une, et les convertisseurs hors limite au point de fonctionnement |
| `POST /api/reparer` | applique une reparation (`action`, `args`) ; rend les tables |
| `POST /api/elements` | les elements assembles et leurs parametres |
| `POST /api/fiche` | la fiche d'un element (`table`, `index`) : ses equations ordonnees par sous-systeme, la valeur de ses parametres et, si le load flow converge, son point d'equilibre |
| `POST /api/taches/{genre}` | lance une tache de fond : `robustesse`, `statismes`, `continuation`, `impedance`, `rapport`, `apres_evenement`, `cascade` ; rend son identifiant |
| `GET /api/taches` | les taches en cours et recentes |
| `GET /api/taches/{id}` | l'etat d'une tache, et son resultat si elle est finie |
| `GET /api/taches/{id}/flux` | l'avancement en flux SSE ; le dernier evenement porte le resultat |
| `DELETE /api/taches/{id}` | arrete la tache au prochain point de controle |
| `POST /api/topologie` | bus, branches, positions (placement automatique) |
| `POST /api/loadflow` | tensions, angles, injections, pertes |
| `POST /api/court_circuit` | la puissance de court-circuit et le SCR de chaque noeud alternatif (machines en x''d et facteurs CEI 60909 par defaut, ou en x'd) |
| `POST /api/modal` | valeurs propres, amortissements, etats participants |
| `POST /api/mode` | forme modale et sensibilites d'UN mode, designe par son indice dans le spectre |
| `POST /api/balayage` | lieu des modes quand un parametre varie (60 valeurs au plus), avec les participations par element a chaque pas |
| `POST /api/signaux` | entrees, etats, sorties **et mesures** du modele |
| `POST /api/temporel` | reponse a un echelon ou a un scenario, par l'une des deux voies ou les deux, avec les mesures demandees |
| `POST /api/apres_evenement` | tension de chaque bus, P et Q de chaque unite et frequence apres un evenement |
| `POST /api/cascade` | la cascade de surtension en regimes permanents successifs |
| `POST /api/robustesse` | les modes qui bougent quand la puissance varie ou qu'un element declenche (aussi en tache de fond) |
| `POST /api/exporter` | la cellule de donnees |
| `GET /api/fichiers` | les fichiers natifs : exemples et sauvegardes, avec leur origine |
| `POST /api/fichiers` | enregistre un reseau dans les sauvegardes |
| `POST /api/fichiers/lire` | n'importe quelle forme de reseau, ramenee a la forme de lecture |
| `POST /api/fichiers/convertir` | le reseau au format `.g2egrid.json`, a telecharger |
| `POST /api/fichiers/externe` | un fichier MATPOWER, PSS/E ou pandapower (`nom`, `contenu`), sous la forme de lecture, avec ses avertissements |
| `POST /api/point/un_slack` | un slack, tout le reste en PQ, au meme point de fonctionnement |
| `POST /api/t_eig` | le tableau des modes des notebooks |
| `POST /api/impact` | l'impact de la distribution sur le transport |
| `POST /api/multi/comparer`, `/api/multi/penetration`, `/api/multi/ajouter_gfl` | transport seul et distribution seule face au modele complet ; balayage de penetration GFL ; GFL ajoutes au hasard dans la distribution |
| `POST /api/simplification/inventaire`, `/api/simplification/comparer` | les types d'elements, leurs etats et les prereglages ; les modes sans puis avec les hypotheses |
| `POST /api/reglage/boucles` | les boucles d'un GFM/GFL : tau, zeta, et les gains qui en sortent |
| `POST /api/machine/standard` | fiche constructeur <-> parametres du modele de Park |
| `POST /api/disjoncteurs` | les elements manoeuvrables, et ce que l'ouverture des disjoncteurs a fait (iles, bus hors tension, nouveau slack) |
| `POST /api/code` | le script Python de l'Atelier |
| `POST /api/zones/aplatir`, `/api/zones/definir` | un reseau a zones, rendu a plat ; un ensemble de bus devient une zone raccordee au transport |
| `POST /api/reseau/importer` | importe un reseau (catalogue ou fichier) et le raccorde |
| `POST /api/exporter/{quoi}` | `donnees`, `notebook`, `matlab` ou `matpower`, avec ses avertissements |
| `GET /api/equations` | les modeles dont on peut lire les equations |
| `GET /api/equations/{nom}` | les equations d'un modele, en LaTeX et en texte SymPy |
| `POST /api/equations/reseau`, `/element` | la liste globale des equations du reseau ; le modele a equations d'un element |
| `POST /api/equations/modifier`, `/ajouter`, `/supprimer`, `/parametre`, `/revenir` | l'edition : rend les nouvelles tables et l'emplacement de l'element (`table`, `index`) |
| `POST /api/equations/verifier`, `/depuis_schema`, `/vers_blocs` | compilation a blanc d'un modele ; passages entre schema-bloc et modele a equations |
| `GET /api/blocs/types` | les types de blocs et leurs champs (la palette en est construite) |
| `POST /api/blocs/compiler` | les variables et les equations d'un schema-bloc |
| `POST /api/blocs/regrouper`, `/degrouper`, `/aplatir` | sous-systemes |
| `POST /api/blocs/depuis_element`, `/remplacer` | le schema-bloc d'une classe ou d'un element du reseau ; le meme, mis a la place de l'element |
| `GET`, `POST /api/blocs/bibliotheque`, `GET`, `DELETE /api/blocs/bibliotheque/{nom}` | la bibliotheque de modeles |
| `GET /api/doc`, `POST /api/doc/construire` | la documentation est-elle construite ; la construire en tache de fond |
| `POST /api/cache/vider` | vide le cache des reseaux assembles |

## Notes de mise en oeuvre

- **Aucune dependance de build, aucun framework** : une page, un fichier JS,
  un CSS. Les schemas et les graphes sont du SVG construit a la main.
- **Rien n'est persiste cote serveur** : le navigateur detient le reseau et le
  renvoie entier a chaque appel. Enregistrer telecharge un JSON, ouvrir en
  relit un.
- Les reseaux assembles sont **mis en cache** cote serveur, indexes par une
  empreinte de leurs tables : assembler un gros reseau et resoudre son
  equilibre prend une dizaine de secondes, et chaque onglet redemande le meme.
- Chaque edition incremente un compteur de **version** ; un resultat calcule
  pour une version perimee n'est jamais affiche.
- **Teste dans un vrai navigateur** : `tests/e2e/` lance un serveur avec jeton
  et rejoue dans Chromium (Playwright) les gestes principaux -- porte fermee
  sans jeton, annuler / retablir, reparer, mode sur le schema, continuation
  en tache, ouverture d'un `.m`, rapport telecharge puis rouvert, frise et rejeu.
  `uv sync --all-extras --group e2e`, puis `uv run pytest tests/e2e`.
