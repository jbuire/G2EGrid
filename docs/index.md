# Bienvenue dans la documentation de g2egrid

```{image} _static/banniere.svg
:alt: g2egrid
:width: 100%
```

**g2egrid** est un outil d'analyse petit-signal et de simulation temporelle des
réseaux électriques à machines synchrones et à convertisseurs. Il construit un
réseau de machines, de convertisseurs grid-forming et grid-following, de bus
infinis, de lignes, de transformateurs et de charges ; il en calcule le load
flow ; il linéarise le tout en **un seul modèle d'état en boucle fermée** pour
l'analyse modale (valeurs propres, participations, sensibilités, lieux des
racines) ; et il intègre ce même modèle, non linéaire, dans le temps.

C'est le portage Python de G2Elin, développé au G2Elab. Chaque équipement y
déclare **une seule fois** son modèle algébro-différentiel en SymPy : les
jacobiens sont dérivés analytiquement puis compilés, une fois par classe, et
le réseau assemble ces modèles par une loi des nœuds exacte. Les modèles sont
des portages fidèles des classes MATLAB de G2Elin V2, validés contre elles
terme à terme, et contre l'outil indépendant ANDES.

```{figure} _static/figures/modes_kundur2area.png
:alt: Carte des valeurs propres du réseau de Kundur
:width: 100%

Carte des valeurs propres du réseau à deux zones de Kundur : tout le spectre
(à gauche, axe réel en échelle logarithmique symétrique) et les modes lents,
étiquetés par l'élément qui y participe le plus (à droite).
```

Ce qu'on peut faire avec :

- **Construire ou charger un réseau** — vingt-quatre réseaux de test (SMIB, WSCC
  9 bus, Kundur deux zones, IEEE 14, 33, 39 et 57 bus, Nordic 44, CIGRE MV,
  couplages transport-distribution, liaisons à courant continu, CIGRE B4 DC
  Grid Test System), un réseau
  dessiné dans l'atelier, ou un fichier MATPOWER, PSS/E ou pandapower, avec des
  disjoncteurs sur chaque branche et chaque unité.
- **Calculer le point de fonctionnement** — load flow de Newton, alternatif ou
  alternatif/continu unifié, équilibre exact de chaque équipement.
- **Choisir les modèles** — machine synchrone à AVR, PSS et régulateur de
  vitesse ; quatre grid-forming à statisme, plus le dVOC et la commande par
  appariement ({doc}`guide/lois_gfm`), deux grid-following ; et, pour chacun, **quelles
  dynamiques garder** ({doc}`guide/dynamiques`).
- **Analyser les modes** — cartes de valeurs propres, participations, formes
  modales, réponses, et **quel paramètre physique déplace un mode**
  ({doc}`guide/sensibilites`).
- **Simuler dans le temps** — le modèle non linéaire, avec échelons de
  consigne, ouvertures de lignes et déclenchements, et la réponse linéarisée
  superposée.
- **Séparer transport et distribution** — ce que la dynamique de la
  distribution change aux modes du transport ({doc}`manuel/09_transport_distribution`).
- **Chercher les limites** — suivre un paramètre jusqu'à la bifurcation de
  Hopf, repérer les convertisseurs hors limite ({doc}`guide/continuation`).
- **Lire la stabilité en impédances** — boucle mineure d'un convertisseur,
  Nyquist généralisé, passivité ({doc}`guide/impedance`).
- **Situer un mode** — résidus, commandabilité, observabilité, grammiens, rapport
  par mode critique ({doc}`guide/residus`).
- **Passer à l'échelle** — modes critiques d'un réseau de milliers d'états sans le
  spectre complet ({doc}`guide/grands_reseaux`).
- **Échanger** — lire MATPOWER, PSS/E et pandapower, écrire MATPOWER
  ({doc}`guide/interoperabilite`).

Un exemple court et exécutable pour **chaque** fonctionnalité :
{doc}`exemples`.

g2egrid s'utilise comme bibliothèque Python, depuis les notebooks, ou dans
l'atelier web. Pour commencer : {doc}`prise_en_main/installation`, puis
{doc}`prise_en_main/premier_reseau`.

## Comment les pièces s'assemblent

```{mermaid}
flowchart TD
    subgraph Donnees["Description du réseau"]
        TAB["tables<br/>Y_DER, T_Lines, T_SMs..."]
        CAT["catalogue<br/>réseaux de test"]
        DISJ["disjoncteurs · zones"]
    end
    subgraph Noyau["Noyau g2egrid"]
        EL["elements<br/>SM, GFM, GFL, IB, lignes, charges"]
        SYM["symbolic<br/>SymPy → jacobiens compilés"]
        NET["network<br/>load flow, équilibre, câblage LFT"]
        DAE["DAE global<br/>résolution, réduction d'indice"]
    end
    subgraph Analyses["Analyses"]
        MOD["modal<br/>modes, participations"]
        SENS["sensibilites<br/>dλ/dp"]
        SEP["separation · impact<br/>transport / distribution"]
        TEMP["temporel<br/>simulation non linéaire"]
        LIM["continuation · impedance<br/>limites, Hopf, boucle mineure"]
        RES["residus · descripteur<br/>modes critiques, grands réseaux"]
    end
    IO["interop<br/>MATPOWER, PSS/E, pandapower"]
    API["g2egrid_api<br/>FastAPI + atelier web"]

    CAT --> TAB
    DISJ --> TAB
    IO --> TAB
    TAB --> NET
    EL --> SYM
    SYM --> NET
    NET --> DAE
    DAE --> MOD
    DAE --> SENS
    DAE --> SEP
    NET --> TEMP
    DAE --> LIM
    DAE --> RES
    MOD --> API
    LIM --> API
    RES --> API
    SENS --> API
    TEMP --> API
    SEP --> API

    click TAB "api/auto/g2egrid.reseau.tables.html" "g2egrid.reseau.tables"
    click CAT "api/auto/g2egrid.reseau.catalogue.html" "g2egrid.reseau.catalogue"
    click DISJ "api/auto/g2egrid.reseau.disjoncteurs.html" "g2egrid.reseau.disjoncteurs"
    click EL "api/auto/g2egrid.elements.html" "g2egrid.elements"
    click SYM "api/auto/g2egrid.noyau.symbolic.html" "g2egrid.noyau.symbolic"
    click NET "api/auto/g2egrid.reseau.network.html" "g2egrid.reseau.network"
    click DAE "guide/resolution_globale.html" "résolution globale"
    click MOD "api/auto/g2egrid.analyse.modal.html" "g2egrid.analyse.modal"
    click SENS "guide/sensibilites.html" "sensibilités"
    click SEP "api/auto/g2egrid.multi_niveau.separation.html" "g2egrid.multi_niveau.separation"
    click TEMP "api/auto/g2egrid.temporel.html" "g2egrid.temporel"
    click LIM "guide/continuation.html" "continuation et impédances"
    click RES "guide/residus.html" "résidus et grands réseaux"
    click IO "guide/interoperabilite.html" "interopérabilité"
    click API "prise_en_main/atelier.html" "atelier web"
```

Les **tables** décrivent le réseau comme les notebooks et les scripts MATLAB ;
le **catalogue** en fournit vingt-quatre, et les **disjoncteurs** en retirent ce
qui est ouvert. Chaque **élément** déclare ses équations, que la couche
**symbolique** dérive et compile. Le **réseau** calcule le load flow et
l'équilibre, puis câble les éléments en un **DAE global**, résolu d'un bloc
({doc}`guide/resolution_globale`). Le modèle d'état qui en sort alimente
l'**analyse modale**, les **sensibilités** et la séparation
**transport/distribution** ; le même câblage, gardé non linéaire, alimente la
**simulation temporelle**. L'**atelier web** expose le tout. Un clic sur un
bloc ouvre sa page.

## État d'avancement

| Fonction | État |
|---|---|
| Load flow AC et AC/DC unifié | fait |
| Modèles SM, GFM1–4, GFL1–2, IB, lignes, transformateurs, charges Z/I/P | fait, validés contre MATLAB V1/V2 |
| GFM à oscillateur virtuel (dVOC) et à appariement (matching) | fait |
| Analyse modale, formes modales, réponses | fait |
| Sélection des dynamiques ($\dot x = 0$, $x = x_0$) | fait |
| Résolution globale, nœuds idéaux, réduction d'indice | fait — analyse linéaire |
| Repère flottant (centre d'inertie) | fait |
| Sensibilités paramétriques analytiques | fait |
| Simulation non linéaire, manœuvres, limiteurs | fait |
| Transport / distribution, zones | fait |
| Continuation, Hopf, limites des convertisseurs | fait |
| Critère d'impédance (boucle mineure, Nyquist généralisé, passivité) | fait — vérifié contre le modèle d'état |
| Résidus, grammiens, rapport par mode critique | fait |
| Faisceau creux et shift-invert pour les grands réseaux | fait |
| Import MATPOWER, PSS/E RAW 33-34, pandapower ; export MATPOWER | fait — CGMES non couvert |
| Atelier web | fait — quatorze onglets, interface en français et en anglais, tâches de fond, jeton d'accès, rapport d'étude |

Ce tableau décrit ce qui est implémenté **et testé** dans ce dépôt ; les
limites de chaque fonction sont données sur sa page.

Certaines pages renvoient à `etudes/` (les études de recherche et leurs
papiers) et à `developpement/` (les outils du portage depuis MATLAB) : ces
dossiers vivent dans l'arbre de développement de l'outil et ne font pas partie
de sa version publiée.

## Contact

- Signaler un problème ou proposer une évolution : section
  [Issues](https://github.com/jbuire/g2elin_v2/issues) du dépôt.
- Le code source : [github.com/jbuire/g2elin_v2](https://github.com/jbuire/g2elin_v2).

## Sommaire

```{toctree}
:maxdepth: 2
:caption: Prise en main

prise_en_main/installation
prise_en_main/premier_reseau
exemples
manuel/06_usage
prise_en_main/atelier
```

```{toctree}
:maxdepth: 2
:caption: Guide

guide/index
guide/construire
guide/analyse_modale
guide/dynamiques
guide/resolution_globale
guide/repere_flottant
guide/sensibilites
guide/robustesse
guide/continuation
guide/residus
guide/impedance
guide/court_circuit
guide/lois_gfm
guide/grands_reseaux
guide/multi_echelle
guide/simplification
guide/domaines
guide/evenements
guide/protections_dc
guide/mmc
guide/cigre_b4
guide/datacenter
guide/nucleaire
guide/interoperabilite
guide/automatiser
guide/ajouter_un_equipement
guide/schemas_blocs
guide/equations_modifiables
```

```{toctree}
:maxdepth: 2
:caption: Théorie et modèles

manuel/01_introduction
manuel/02_principes
reperes_et_angles
manuel/03_modeles
equations_modeles
manuel/A_per_unit
manuel/B_reglage
```

```{toctree}
:maxdepth: 2
:caption: Études et notebooks

manuel/07_exemples
manuel/08_sources
manuel/09_transport_distribution
manuel/10_transition
notebooks
```

```{toctree}
:maxdepth: 2
:caption: Architecture et validation

manuel/05_python
manuel/12_atelier
manuel/04_v2
manuel/11_validation
equivalence_matlab
guide/validation
```

```{toctree}
:maxdepth: 1
:caption: Référence

api/index
```

## Avertissement

Les résultats des réseaux de test sont obtenus avec les réglages par défaut
des équipements, qui ne sont pas ceux des publications d'origine : ils ne se
comparent pas directement à la littérature. Les marques citées dans la
documentation ou le code appartiennent à leurs propriétaires respectifs.
