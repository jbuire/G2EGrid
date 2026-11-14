# Guide de l'utilisateur

Ce guide suit le déroulé d'une étude : on **construit** un réseau, on choisit
le niveau de détail de ses **dynamiques** et la façon de le **résoudre**, on
fixe son **repère**, puis on lit ses **modes** et ce qui les déplace.

| Page | Pour |
|---|---|
| {doc}`construire` | déclarer un réseau en code ou par ses tables, les réglages par famille et par bus |
| {doc}`analyse_modale` | modes, participations, forme modale, réponses, simulation non linéaire |
| {doc}`dynamiques` | supprimer une dynamique : $\dot x = 0$ ou $x = x_0$ ; les variantes comme sélections |
| {doc}`resolution_globale` | le DAE global, les nœuds idéaux et la réduction des courants en série |
| {doc}`repere_flottant` | un angle de référence détaché du slack |
| {doc}`sensibilites` | quel paramètre déplace un mode, et de combien |
| {doc}`robustesse` | quels modes bougent quand la puissance varie ou qu'un élément déclenche (N-1) |
| {doc}`continuation` | jusqu'où le réseau reste stable : suivi d'un paramètre, Hopf, limites des convertisseurs |
| {doc}`residus` | par où exciter un mode, où le voir : résidus, grammiens, rapport par mode critique |
| {doc}`court_circuit` | puissance de court-circuit, X/R et SCR de chaque nœud ; WSCR ; réseaux faibles |
| {doc}`impedance` | la stabilité vue des impédances : boucle mineure, Nyquist généralisé, passivité |
| {doc}`lois_gfm` | GFM à statisme, à oscillateur virtuel (dVOC) et à appariement (matching) ; impédance virtuelle |
| {doc}`grands_reseaux` | les modes critiques d'un grand réseau sans le spectre complet |
| {doc}`multi_echelle` | transport et distribution : modèles réduits, GFL au hasard, pénétration |
| {doc}`simplification` | les hypothèses dx/dt = 0 par type ou par équipement, et ce qu'elles font perdre |
| {doc}`evenements` | après un événement : tension de chaque bus, P et Q des unités, relais de surtension des convertisseurs, cascade quasi statique (onglet Événements) |
| {doc}`protections_dc` | défauts continus, blocage des convertisseurs, disjoncteurs DC et plans de protection, défaut alternatif aux bornes d'une station |
| {doc}`domaines` | modes alternatifs, continus et couplés AC/DC ; îles ; impact du continu |
| {doc}`mmc` | topologie des câbles continus (monopôle, bipôle, circuits) et modèle moyen du MMC, énergie libre ou régulée |
| {doc}`cigre_b4` | le CIGRE B4 DC Grid Test System : stations offshore formeuses, statismes Vdc-P, références multiples |
| {doc}`datacenter` | un centre de données (redresseur à diodes, onduleur V/f, batterie) et son groupe de secours : courbe de charge de plusieurs minutes, perte du réseau, fatigue de l'arbre |
| {doc}`nucleaire` | groupes nucléaires : catalogue de machines (SMR, IEEE 39 nœuds, quatre machines), lignes d'arbre multi-masses, réseaux des thèses de Boudot et de Mazauric |
| {doc}`interoperabilite` | lire MATPOWER, PSS/E et pandapower, écrire MATPOWER |
| {doc}`automatiser` | suivre et interrompre un calcul long, l'atelier par script, le rapport d'étude |
| {doc}`ajouter_un_equipement` | écrire une nouvelle classe d'équipement |
| {doc}`schemas_blocs` | décrire un équipement par un schéma-bloc, compilé en équations (onglet Blocs) |
| {doc}`equations_modifiables` | lire et modifier les équations dynamiques d'un équipement ou de tout le réseau |
| {doc}`validation` | les tests et la validation croisée |

Un exemple court pour chaque fonctionnalité : {doc}`../exemples`.

Les notions de fond — élément DAE, repères, load flow, linéarisation,
assemblage LFT — sont exposées au chapitre {doc}`../manuel/02_principes`.
