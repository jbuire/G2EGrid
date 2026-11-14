"""Les corps de requete de l'API.

Ils sont volontairement permissifs sur les tables (``dict`` brut) : la
validation du reseau est un resultat a rendre a l'utilisateur, pas un rejet
au parsing. Un reseau en cours de construction est presque toujours
incomplet, et l'interface doit pouvoir le poster pour savoir ce qui lui
manque.
"""

from __future__ import annotations

from pydantic import BaseModel, Field

#: un flottant fini : ni NaN ni infini (JSON de Python les laisse passer)
FINI = {"allow_inf_nan": False}


class CorpsReseau(BaseModel):
    """Le reseau tel que le constructeur le manipule."""

    tables: dict = Field(default_factory=dict)
    positions: dict | None = None


class CorpsReparation(CorpsReseau):
    """Une reparation proposee par la validation, a appliquer."""

    action: str
    args: dict = Field(default_factory=dict)


class CorpsFichier(BaseModel):
    """Un reseau STRUCTURE : exactement le contenu d'un ``.g2egrid.json``.

    Ce corps porte le reseau tel que le FICHIER le decrit -- transport, zones,
    instances -- et non ses tables a plat. C'est le seul endroit de l'API ou la
    distinction compte : l'aplatir avant de l'envoyer perdrait precisement ce
    qu'on veut enregistrer, a savoir qu'une meme distribution est decrite une
    fois et instanciee trois.

    ``reseau`` reste un ``dict`` brut pour la meme raison que ``tables`` :
    c'est ``g2egrid.reseau.zones.reseau_depuis_json`` qui en juge, et son refus est un
    message a rendre a l'utilisateur, pas un rejet au parsing.
    """

    reseau: dict = Field(default_factory=dict)
    #: le nom du fichier ; a defaut, celui que porte le reseau lui-meme
    nom: str | None = None


class CorpsModal(CorpsReseau):
    #: frequence maximale des modes rendus [Hz] ; a defaut, tous
    f_max: float | None = Field(None, gt=0, **FINI)
    n_etats: int = Field(3, ge=0, le=1000)


class CorpsTemporel(CorpsReseau):
    entree: str
    amplitude: float = Field(0.01, **FINI)
    #: duree simulee, en secondes, strictement positive
    duree: float = Field(1.0, gt=0, le=1e6, **FINI)
    t0: float = Field(0.0, ge=0, **FINI)
    #: nombre d'instants de calcul, au moins 2
    n: int = Field(401, ge=2, le=1_000_000)
    etats: list[str] | None = None
    voie: str = "les_deux"
    #: scenario a plusieurs perturbations ; a defaut, l'echelon ci-dessus
    evenements: list[dict] | None = None
    #: grandeurs physiques a relever en plus des etats (voie non lineaire)
    mesures: list[str] | None = None
    #: integrateur de la voie non lineaire : auto, Radau ou BDF2 (pas fixe)
    methode: str = "auto"
    pas: float | None = Field(None, gt=0, **FINI)
    limiteurs: bool = True
    #: le plan de protection type du reseau continu (g2egrid.temporel.protections_dc)
    protections_dc: bool = False


class CorpsApresEvenement(CorpsReseau):
    """Les grandeurs electriques apres un evenement (:func:`analyse.apres_evenement`)."""

    #: le scenario, au format de CorpsTemporel.evenements (au moins un)
    evenements: list[dict]
    duree: float = 5.0
    n: int = 501
    #: pas de BDF2 [s] ; null : 1/(400 f0)
    pas: float | None = 1e-3
    limiteurs: bool = True
    #: relais de surtension des convertisseurs : {"v_max", "retard",
    #: "retombee" (null : verrouille), "action", "familles"} ; null : aucun
    relais: dict | None = None
    protections_dc: bool = False
    v_bas: float = 0.9
    v_haut: float = 1.1
    #: la fenetre finale sur laquelle on moyenne le regime atteint [s]
    fenetre: float = 0.2


class CorpsCascade(CorpsReseau):
    """La cascade de surtension en regimes permanents (:func:`analyse.cascade`)."""

    initiales: list[str]
    v_max: float = 1.10
    familles: list[str] = Field(default_factory=lambda: ["GFL"])
    #: "vague" (ensemble) ou "une_a_une" (la plus haute d'abord)
    ordre: str = "vague"


class CorpsImport(BaseModel):
    """Importer un reseau dans celui qu'on construit.

    ``source`` est ``{"catalogue": nom}`` ou le contenu d'un fichier ;
    ``raccord`` vaut ``{"mode": "fusion" | "ligne" | "transfo", "bus_importe":
    b, "bus_existant": B}`` ou ``None``.
    """

    reseau: dict
    source: dict
    raccord: dict | None = None
    garder_slack: bool = False
    garder_ib: bool = False


class CorpsGflAleatoires(BaseModel):
    """Ajouter des GFL au hasard dans les distributions d'un reseau a zones."""

    reseau: dict
    n: int = 5
    graine: int | None = None
    part: float = 0.05


class CorpsMulti(BaseModel):
    """Confronter les modeles reduits (transport seul, distribution seule)."""

    reseau: dict
    instances: list[int] | None = None
    n_max: int = 30
    pas: int = 10
    graine: int = 1
    part: float = 0.05


class CorpsZone(BaseModel):
    """Definir une zone : le reseau (a plat ou a zones), ses bus, son nom."""

    reseau: dict
    bus: list[int]
    nom: str = "zone"


class CorpsMode(CorpsReseau):
    """Le detail d'un mode : son indice DANS LE SPECTRE, pas dans la table."""

    #: de 0 au nombre de modes moins un (la borne haute depend du spectre :
    #: elle est verifiee par ``analyse.detail_mode``)
    indice: int = Field(ge=0)
    n_etats: int = Field(8, ge=0, le=1000)
    n_sensibilites: int = Field(8, ge=0, le=1000)
    #: nombre de parametres d'equipement rendus, par sensibilite normalisee
    n_parametres: int = Field(10, ge=0, le=1000)


class CorpsBalayage(CorpsReseau):
    """Un parametre a faire varier, et les valeurs a essayer.

    ``cible`` designe ou vit le parametre : ``{"genre": "famille",
    "famille": "GFM_PARAMS", "cle": "mp"}``, ``{"genre": "bus", "bus": 5,
    "cle": "mp"}`` ou ``{"genre": "table", "table": "T_Lines_pu",
    "ligne": 0, "colonne": "X_pu"}``.
    """

    cible: dict
    valeurs: list[float]
    f_max: float | None = None
    n_modes: int = 12


class CorpsValidation(CorpsReseau):
    """Les trois epreuves du petit signal face a l'EMT."""

    entree: str | None = None
    amplitude: float = 0.01
    duree: float = 4.0
    dt: float = 0.005
    t_evenement: float = 0.1
    voies: list[str] | None = None
    f_max: float = 6.0


class CorpsUnSlack(CorpsReseau):
    slack_bus: int


class CorpsRobustesse(CorpsReseau):
    """Etude de robustesse des modes : variations de puissance et N-1."""

    f_min: float = 0.05
    f_max: float = 5.0
    zeta_max: float | None = 30.0
    frac: float = 0.1
    puissance: bool = True
    deconnexion: bool = True


class CorpsSimplification(CorpsReseau):
    """Hypotheses de simplification du calcul des modes : ``{type ou element:
    {etat: "statique" | "fige" | "dynamique"}}``, et des preregles."""

    simplifications: dict[str, dict[str, str]] = Field(default_factory=dict)
    preregles: list[str] = Field(default_factory=list)


class CorpsImpact(BaseModel):
    reseau: dict = Field(default_factory=dict)
    f_max: float = 5.0
    par_poste: bool = True


class CorpsComparer(BaseModel):
    """Deux versions : chacune un reseau (toute forme) ou ``{"catalogue": nom}``."""

    a: dict = Field(default_factory=dict)
    b: dict = Field(default_factory=dict)
    f_max: float | None = None


class CorpsReglage(BaseModel):
    famille: str
    reglages: dict = Field(default_factory=dict)
    f0: float = 50.0


class CorpsMachine(BaseModel):
    valeurs: dict = Field(default_factory=dict)
    f0: float = 50.0
    sens: str = "vers_standard"


class CorpsStatismes(CorpsReseau):
    """Balayage du statisme Vdc-P des stations en mode droop."""

    debut: float = 0.25
    fin: float = 4.0
    pas: int = 9
    stations: list[int] | None = None


class CorpsContinuation(CorpsReseau):
    """Continuation sur un parametre d'element."""

    element: str
    parametre: str
    debut: float
    fin: float
    pas: int = 15
    log: bool = False
    equilibre: bool = False


class CorpsBoucleMineure(CorpsReseau):
    """Critere de la boucle mineure d'un convertisseur."""

    element: str
    f_min: float = 0.1
    f_max: float = 5000.0
    n: int = 800


class CorpsImportExterne(BaseModel):
    """Un fichier d'un autre outil : son nom (l'extension decide) et son texte."""

    nom: str
    contenu: str
    charges: str = "Z"


class CorpsRapport(CorpsReseau):
    """Le rapport d'etude : figures SVG de l'atelier a inclure."""

    figures: dict[str, str] = Field(default_factory=dict)
    n_modes: int = 3


class CorpsFiche(CorpsReseau):
    """La fiche d'un element : sa ligne de table."""

    table: str
    #: rang de la ligne dans la table, a partir de 0
    index: int = Field(ge=0)


class CorpsBlocs(BaseModel):
    """Un schema-bloc a compiler (voir :mod:`g2egrid.elements.blocs`)."""

    schema_bloc: dict = Field(default_factory=dict)


class CorpsBlocsRegrouper(CorpsBlocs):
    """Regrouper des blocs d'un niveau en sous-systeme."""

    ids: list[str] = Field(default_factory=list)
    nom: str = ""
    #: identifiants de sous-systemes qui menent au niveau a editer (facultatif :
    #: l'interface peut aussi envoyer le seul niveau affiche)
    chemin: list[str] = Field(default_factory=list)


class CorpsBlocsDegrouper(CorpsBlocs):
    """Dissoudre un sous-systeme dans son niveau."""

    id: str = ""
    chemin: list[str] = Field(default_factory=list)


class CorpsBlocsDepuisElement(BaseModel):
    """L'element de la bibliotheque a ecrire en blocs : une classe native
    (``classe``, ``parametres``), ou un element d'un reseau (``tables``,
    ``table``, ``index``) avec ses parametres."""

    classe: str | None = None
    parametres: dict | None = None
    detail: str = "blocs"
    f0: float | None = Field(default=None, gt=0)
    tables: dict | None = None
    table: str | None = None
    index: int | None = Field(default=None, ge=0)


class CorpsEquationsElement(CorpsReseau):
    """Un element du reseau : sa ligne (``table``, ``index``) ou son nom
    (``element``, qui leve l'ambiguite d'une ligne de ``Y_DER``)."""

    table: str | None = None
    index: int | None = Field(default=None, ge=0)
    element: str | None = None


class CorpsEquationsModifier(CorpsEquationsElement):
    variable: str
    expression: str
    #: ``etat``, ``algebrique`` ou ``sortie`` (nom a la fois etat et sortie)
    genre: str | None = None


class CorpsEquationsAjouter(CorpsEquationsElement):
    genre: str
    variable: str
    expression: str
    depart: float | None = None


class CorpsEquationsSupprimer(CorpsEquationsElement):
    variable: str


class CorpsEquationsParametre(CorpsEquationsElement):
    """Un parametre du modele d'un element deja converti (``T_Blocs``)."""

    nom: str
    valeur: float = Field(**FINI)


class CorpsBlocsRemplacer(CorpsEquationsElement):
    """Remplacer un element natif par son schema-bloc."""

    detail: str = "blocs"


class CorpsEquationsModele(BaseModel):
    modele: dict = Field(default_factory=dict)


class CorpsCourtCircuit(CorpsReseau):
    """Puissance de court-circuit a chaque noeud (voir
    :mod:`g2egrid.analyse.court_circuit`)."""

    gfm: str = "tension"
    gfl: str = "aucun"
    charges: bool = False
    dynamique: bool = False
    f: float = 1.0
    machine: str = "soustransitoire"
