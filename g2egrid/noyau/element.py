"""Classe de base des equipements.

CONVENTIONS (valables pour tout le paquet)
------------------------------------------
::

    Repere      : un unique repere dq global tournant a w0 = 2*pi*f0 (synchrone).
                  Un phaseur de bus V/_theta s'ecrit  vd = V cos(theta),
                  vq = V sin(theta). Chaque machine porte son propre angle
                  interne relatif a ce repere -> pas de "slack qui definit le
                  repere", donc pas de convention d'angle rotor fragile.
    Courants    : le courant d'un PORT est celui INJECTE dans le bus
                  (positif = l'element fournit du courant au reseau).
                  Generateur : i > 0 ; charge : i < 0 ; branche : -i cote from,
                  +i cote to. La loi des noeuds est donc  somme(i_ports) = 0.
    Puissance   : S = V conj(I)  =>  P = vd.id + vq.iq ,  Q = vq.id - vd.iq
    Per-unit    : base systeme Sb ; les impedances sont en pu base systeme.

Chaque port ``p`` ajoute automatiquement au modele ::

    entrees  : vd_<p>, vq_<p>      (tension du bus, imposee par le reseau)
    sorties  : id_<p>, iq_<p>      (courant injecte, vu par la loi des noeuds)

Tous les elements recoivent aussi l'entree globale ``w_ref`` (frequence
systeme en pu).

PORTS CONTINUS -- un port DC ``p`` (declare dans ``ports_dc``) n'apporte
qu'UN signal de chaque cote : entree ``v_<p>``, sortie ``i_<p>``. Ces signaux
s'ajoutent EN QUEUE de ceux des ports alternatifs, jamais au milieu : c'est ce
qui garantit que ``port_input_idx``, ``port_output_idx``, ``wref_idx`` et
``thetaref_idx`` gardent exactement la meme valeur qu'avant, et qu'aucun
modele existant ne voit un seul de ses signaux se decaler. Cote sorties, les
courants DC viennent donc APRES les courants alternatifs mais AVANT les
sorties de service : les ``n_port_signaux()`` premieres sorties sont les
courants de port, internes au cablage, et le reste est expose.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from g2egrid.noyau.symbolic import CompiledModel, ModelSpec, compile_spec

__all__ = [
    "Element",
    "InitContext",
    "port_v_names",
    "port_i_names",
    "port_v_name_dc",
    "port_i_name_dc",
    "W_REF",
    "THETA_REF",
]

W_REF = "w_ref"  # entree globale : frequence systeme (pu)
# Entree globale : angle du repere de reference, c'est-a-dire l'angle de
# l'element slack. Les DG qui possedent leur propre angle tournent de
# ``theta_g - theta``, comme ``theta_g`` dans symSG / symGFM / symGFL.
# Les noeuds et les lignes ne l'utilisent pas (leurs ug_vec n'en ont pas non
# plus cote V1) : ils en heritent indirectement, par le courant injecte.
THETA_REF = "theta_g"


def port_v_names(port: str) -> tuple[str, str]:
    return (f"vd_{port}", f"vq_{port}")


def port_i_names(port: str) -> tuple[str, str]:
    return (f"id_{port}", f"iq_{port}")


def port_v_name_dc(port: str) -> str:
    return f"v_{port}"


def port_i_name_dc(port: str) -> str:
    return f"i_{port}"


@dataclass
class InitContext:
    """Point de fonctionnement fourni par le load flow, pour l'initialisation.

    v : dict port -> tension complexe du bus (pu, repere global, V*exp(j.theta))
    s : dict port -> puissance complexe INJECTEE par l'element (pu base systeme)

    Pour un port CONTINU, il n'y a pas de phaseur : ``v[port]`` est la tension
    continue du bus, un REEL, et ``s[port]`` la puissance active injectee,
    reelle elle aussi. Les elements a port DC lisent donc ces deux entrees
    comme des flottants ; ``float(np.real(...))`` les accepte des deux genres.
    """

    v: dict[str, complex]
    s: dict[str, complex] = field(default_factory=dict)
    w: float = 1.0
    f0: float = 50.0
    theta_g: float = 0.0  # angle du repere de reference (element slack)

    @property
    def wb(self) -> float:
        return 2.0 * np.pi * self.f0


class Element:
    """Base declarative d'un equipement (modele DAE + ports reseau).

    Une sous-classe declare (au niveau CLASSE) :
        ports, params, states, algebs, svc_inputs, svc_outputs
    et implemente :
        equations(s) -> (f, g, h)   [statique, symbolique]
        initialize(ctx) -> dict     [numerique, par instance]

    ``h`` doit renvoyer, DANS CET ORDRE : les courants de chaque port
    alternatif (id, iq), puis celui de chaque port continu (i), puis les
    sorties de service.
    """

    # --- declaration (a surcharger) --------------------------------------
    # Nom du modele compile. Deux classes qui partagent EXACTEMENT les memes
    # equations (ex. Line et Transformer) peuvent partager ce nom pour ne
    # compiler qu'une fois.
    model_name: str | None = None
    #: raffiner le point d'equilibre par Newton apres ``initialize`` ?
    #: Les elements qui rejouent une convention d'initialisation de G2Elin V1
    #: (``init='symSG'`` et consorts) la mettent a False : la reference calcule
    #: son equilibre par formules directes, et celui-ci NE VERIFIE PAS ses
    #: propres equations. Le raffiner s'en ecarterait, et les deux matrices A
    #: ne seraient plus comparables terme a terme.
    raffiner_equilibre: bool = True
    #: point de fonctionnement impose a ``initialize``, sous la forme (V, S)
    #: complexes en pu systeme. Quand il est pose, il remplace ce que le load
    #: flow attribue a l'element -- utile pour rejouer la convention de G2Egrid
    #: V1, qui initialise chaque DER sur la ligne de ``bus_sol`` de SON bus,
    #: celui que ``network_form.m`` cree derriere le transfo de raccordement,
    #: et y lit la GENERATION, sans en deduire la charge du bus.
    point_init: tuple[complex, complex] | None = None
    #: element passif (noeud, branche, charge) : ses termes de couplage dq
    #: suivent la pulsation du repere global (voir Network.repere)
    passif: bool = False
    #: le modele vit du cote CONTINU du reseau (marqueur declaratif ; le
    #: cablage, lui, raisonne sur le genre de chaque port)
    dc: bool = False
    ports: tuple[str, ...] = ("t",)
    #: ports CONTINUS, un seul signal de chaque cote (``v_<p>`` / ``i_<p>``).
    #: Les bus correspondants viennent APRES ceux de ``ports`` dans ``buses``.
    ports_dc: tuple[str, ...] = ()
    #: les modes de conduite admis par ce convertisseur, au sens de
    #: ``g2egrid.noyau.loadflow.MODES_VSC`` : les DEUX grandeurs qu'il tient, la
    #: continue puis l'alternative. Vide = l'element n'est pas un
    #: convertisseur, ou ne se declare pas dans ``T_VSCs``. C'est le seul
    #: endroit ou la question se decide : un mode que le modele ne peut pas
    #: tenir doit etre refuse par l'assembleur, pas decouvert au Newton.
    modes_dc: tuple[str, ...] = ()
    params: tuple[str, ...] = ()
    states: tuple[str, ...] = ()
    algebs: tuple[str, ...] = ()
    svc_inputs: tuple[str, ...] = ()
    svc_outputs: tuple[str, ...] = ()

    def __init__(self, name: str, buses: int | list[int], **params):
        if isinstance(buses, int):
            buses = [buses]
        attendus = len(self.ports) + len(self.ports_dc)
        if len(buses) != attendus:
            raise ValueError(
                f"{name}: {len(buses)} bus fourni(s) pour {attendus} port(s)"
            )
        self.name = name
        self.buses = list(buses)
        missing = set(self.params) - set(params)
        if missing:
            raise ValueError(f"{name}: parametre(s) manquant(s) {sorted(missing)}")
        self.p = np.array([float(params[k]) for k in self.params])
        self.pdict = {k: float(params[k]) for k in self.params}
        # rempli par Network.solve_equilibrium
        self.x0: np.ndarray | None = None
        self.z0: np.ndarray | None = None
        self.u0: np.ndarray | None = None

    # --- noms complets ----------------------------------------------------
    @classmethod
    def input_names(cls) -> tuple[str, ...]:
        pv = [n for pt in cls.ports for n in port_v_names(pt)]
        pv += [port_v_name_dc(pt) for pt in cls.ports_dc]
        return tuple(cls.svc_inputs) + (W_REF, THETA_REF) + tuple(pv)

    @classmethod
    def output_names(cls) -> tuple[str, ...]:
        pi = [n for pt in cls.ports for n in port_i_names(pt)]
        pi += [port_i_name_dc(pt) for pt in cls.ports_dc]
        return tuple(pi) + tuple(cls.svc_outputs)

    @classmethod
    def spec(cls) -> ModelSpec:
        return ModelSpec(
            name=cls.model_name or cls.__name__,
            params=cls.params,
            states=cls.states,
            algebs=cls.algebs,
            inputs=cls.input_names(),
            outputs=cls.output_names(),
            equations=cls.equations,
        )

    @classmethod
    def compiled(cls) -> CompiledModel:
        """Modele compile, mis en cache au niveau de la CLASSE."""
        return compile_spec(cls.spec())

    # --- indices pratiques ------------------------------------------------
    @classmethod
    def port_input_idx(cls, k: int) -> tuple[int, int]:
        """Indices (vd, vq) du port k dans le vecteur d'entree."""
        base = len(cls.svc_inputs) + 2  # +2 pour w_ref et theta_g
        return base + 2 * k, base + 2 * k + 1

    @classmethod
    def port_output_idx(cls, k: int) -> tuple[int, int]:
        """Indices (id, iq) du port k dans le vecteur de sortie."""
        return 2 * k, 2 * k + 1

    @classmethod
    def port_dc_input_idx(cls, k: int) -> int:
        """Indice de ``v_<p>`` du k-ieme port CONTINU, dans le vecteur d'entree."""
        base = len(cls.svc_inputs) + 2  # +2 pour w_ref et theta_g
        return base + 2 * len(cls.ports) + k

    @classmethod
    def port_dc_output_idx(cls, k: int) -> int:
        """Indice de ``i_<p>`` du k-ieme port CONTINU, dans le vecteur de sortie."""
        return 2 * len(cls.ports) + k

    @classmethod
    def ports_tous(cls) -> tuple[str, ...]:
        """Tous les ports, alternatifs D'ABORD puis continus -- l'ordre de ``buses``."""
        return tuple(cls.ports) + tuple(cls.ports_dc)

    @classmethod
    def n_port_signaux(cls) -> int:
        """Nombre de signaux de port : 2 par port alternatif, 1 par port continu.

        C'est le SEUL nombre qui dit ou s'arretent les courants de port dans
        ``h`` et ou commencent les sorties de service. Tout ``2 * len(ports)``
        ecrit en dur ailleurs laisserait fuir ``i_<p>`` dans les sorties
        exposees et decalerait la matrice C de tout ce qui suit.
        """
        return 2 * len(cls.ports) + len(cls.ports_dc)

    @classmethod
    def wref_idx(cls) -> int:
        return len(cls.svc_inputs)

    @classmethod
    def thetaref_idx(cls) -> int:
        """Indice de l'entree globale ``theta_g`` (angle du slack)."""
        return len(cls.svc_inputs) + 1

    # --- a surcharger -----------------------------------------------------
    @staticmethod
    def equations(s):
        """Equations symboliques. Renvoie (f, g, h)."""
        raise NotImplementedError

    def initialize(self, ctx: InitContext) -> dict[str, float]:
        """Valeurs initiales de TOUTES les variables et entrees de service.

        Renvoie un dict {nom: valeur}. Les entrees de service dont la valeur
        est imposee par le point de fonctionnement (p_ref, v_ref, ...) sont
        calculees analytiquement ici.
        """
        raise NotImplementedError

    def shunt_admittance(self) -> complex | None:
        """Admittance shunt pu (base systeme) a ajouter au Y-bus du load flow.

        A declarer quand le NOEUD absorbe un courant que le modele de l'element
        ne compte pas dans ses sorties : c'est le cas du noeud resistif, dont la
        tension vaut ``R.isum``. Sinon le point de fonctionnement du load flow
        ne satisfait pas la loi des noeuds du modele dynamique. A l'inverse, un
        shunt deja soustrait du courant injecte par l'element (``Rg`` de la
        machine synchrone) ne doit PAS y figurer : il serait compte deux fois.
        ``None`` = aucun shunt.
        """
        return None

    #: les parametres DRAPEAUX des modes discrets de l'element (limiteurs,
    #: butees) : 0 = libre. Voir :mod:`g2egrid.elements.limiteurs`.
    modes: tuple[str, ...] = ()

    def gardes(self, cm, x, z, u, p) -> list[tuple[str, float, float, int]]:
        """Les gardes des modes discrets, au point ``(x, z, u)`` de l'element.

        Chaque garde est ``(drapeau, g, valeur_apres, sens)`` : quand ``g``
        franchit zero dans le sens ``sens`` (+1 en montant, -1 en descendant),
        le parametre ``drapeau`` prend ``valeur_apres``. ``p`` est le vecteur
        de parametres COURANT (celui du mode en cours). Aucune garde par defaut.
        """
        return []

    def admittance_ybus(self) -> complex | None:
        """Admittance a inscrire au Y-bus, pour un element qui EST un shunt.

        L'inverse exact de :meth:`shunt_admittance` : ici le courant du shunt
        est bien une sortie de l'element (il a ses propres etats), mais le load
        flow doit le voir comme une admittance -- sa puissance suit ``|V|^2``,
        elle n'est pas programmee. ``Network`` l'inscrit donc au Y-bus, et
        attribue a l'element, a l'initialisation, la puissance
        ``-conj(y).|V|^2`` qu'il y consomme, SANS la retrancher du reliquat du
        generateur voisin : le Y-bus l'a deja mise a sa charge. ``None`` = pas
        un shunt.
        """
        return None

    def scheduled_injection(self) -> complex | None:
        """Puissance complexe INJECTEE connue d'avance (pu base systeme).

        ``None`` = inconnue : l'element recupere le reliquat du load flow a
        son bus (cas d'un generateur). Une charge, elle, connait sa
        consommation et renvoie ``-(P + jQ)``.
        """
        return None

    def scheduled_injection_dc(self) -> float | None:
        """Puissance active INJECTEE au port continu, quand l'element la connait.

        Le pendant DC de :meth:`scheduled_injection`, mais sans load flow
        derriere : tant que le lot 4 n'a pas unifie le Newton AC/DC, aucune
        puissance n'est ATTRIBUEE a un element du cote continu. Ne la declarent
        donc que les elements qui la tiennent d'eux-memes -- une charge a
        puissance constante. ``None`` = l'element se debrouille de son propre
        point de depart.
        """
        return None

    # --- helpers ----------------------------------------------------------
    def set_param(self, name: str, value: float) -> None:
        """Met a jour un parametre numerique (utilisable depuis ``initialize``).

        Utile pour les elements dont les parametres derivent du point de
        fonctionnement, p.ex. l'impedance equivalente d'une charge a
        puissance constante calculee a la tension du load flow.
        """
        if name not in self.params:
            raise KeyError(f"{self.name}: parametre inconnu '{name}'")
        self.pdict[name] = float(value)
        self.p[self.params.index(name)] = float(value)

    def build_u0(self, init: dict[str, float], ctx: InitContext) -> np.ndarray:
        """Vecteur d'entree d'equilibre : service + w_ref + tensions de port."""
        u = np.zeros(len(self.input_names()))
        for i, nm in enumerate(self.svc_inputs):
            u[i] = init.get(nm, 0.0)
        u[self.wref_idx()] = ctx.w
        u[self.thetaref_idx()] = ctx.theta_g
        n_ac = len(self.ports)
        for k, pt in enumerate(self.ports_tous()):
            if k < n_ac:
                iv, iq = self.port_input_idx(k)
                u[iv] = complex(ctx.v[pt]).real
                u[iq] = complex(ctx.v[pt]).imag
            else:
                u[self.port_dc_input_idx(k - n_ac)] = float(np.real(ctx.v[pt]))
        return u

    # --- selection des dynamiques ------------------------------------------
    def reduire(
        self,
        statique=(),
        fige=(),
        *,
        strict: bool = True,
    ) -> "Element":
        """Supprime des dynamiques de l'element : ``dx/dt = 0`` ou ``x = x0``.

        Deux facons de dire qu'une dynamique vaut zero, qui ne donnent pas le
        meme modele :

        ``statique``  -- ``dx/dt = 0`` : l'etat devient une variable ALGEBRIQUE
            et son equation differentielle une contrainte ``0 = f(x, z, u)``.
            C'est la perturbation singuliere classique : la dynamique est jugee
            infiniment rapide, la grandeur suit instantanement son equilibre.
            ``Dw`` statique dans ``GFM2`` donne ``Dw = mp.(p_ref - pm)``, le
            statisme du 1er ordre de ``GFM1``.
        ``fige``      -- ``x = x0`` : l'etat devient un PARAMETRE, qui vaut sa
            valeur d'equilibre (ou celle donnee, si ``fige`` est un dict
            ``{nom: valeur}``). La dynamique est jugee infiniment lente, ou
            parfaitement regulee : ``vdc`` fige a 1 dans ``GFM3`` donne le bus
            continu IDEAL de ``GFM1``.

        Les variantes d'une meme famille ne sont alors plus des classes a part
        mais des SELECTIONS de dynamiques du modele complet : ``GFM4`` avec
        ``Dw`` statique et ``vdc``, ``idc`` figes a exactement le spectre de
        ``GFM1`` (``tests/test_dynamiques.py``).

        L'element change de classe -- une sous-classe engendree et mise en
        cache, dont le modele compile porte un nom propre --, si bien que tout
        le reste du paquet (cablage, equilibre, linearisation, simulation) le
        traite comme n'importe quel element. ``isinstance`` reste vrai pour la
        classe d'origine. Les reductions se cumulent d'un appel a l'autre.

        ``strict=False`` ignore les noms que l'element n'a pas : c'est ce qui
        permet de donner une meme selection a toute une famille (``GFM_PARAMS``)
        dont les variantes n'ont pas toutes les memes etats.

        Rend l'element lui-meme, pour chainer : ``net.add(GFM4(...).reduire(...))``.
        """
        if isinstance(statique, str):
            statique = (statique,)
        if isinstance(fige, str):
            fige = (fige,)
        valeurs = dict(fige) if isinstance(fige, dict) else {n: None for n in fige}
        statique = tuple(statique)
        connus = set(self.states)
        inconnus = [n for n in (*statique, *valeurs) if n not in connus]
        if inconnus and strict:
            raise ValueError(
                f"{self.name}: pas d'etat {', '.join(map(repr, inconnus))} ; "
                f"etats de {type(self).__name__} : {', '.join(self.states)}"
            )
        statique = tuple(n for n in statique if n in connus)
        valeurs = {n: v for n, v in valeurs.items() if n in connus}
        doublons = set(statique) & set(valeurs)
        if doublons:
            raise ValueError(
                f"{self.name}: {sorted(doublons)} a la fois statique(s) et fige(s)"
            )
        if not statique and not valeurs:
            return self

        base = getattr(type(self), "_classe_complete", type(self))
        deja_s = tuple(getattr(type(self), "_statiques", ()))
        deja_f = tuple(getattr(type(self), "_figes", ()))
        cls = _variante(base, deja_s + statique, deja_f + tuple(valeurs))
        pdict = dict(self.pdict)
        valeurs_fixees = dict(getattr(self, "_valeurs_figees", {}))
        for n, v in valeurs.items():
            pdict[n] = np.nan if v is None else float(v)
            if v is not None:
                valeurs_fixees[n] = float(v)
        self.__class__ = cls
        self._valeurs_figees = valeurs_fixees
        self.pdict = {k: float(pdict[k]) for k in cls.params}
        self.p = np.array([self.pdict[k] for k in cls.params])
        # le point de fonctionnement n'a plus les bonnes dimensions : a refaire
        self.x0 = self.z0 = self.u0 = None
        return self

    @classmethod
    def dynamiques(cls) -> dict[str, str]:
        """Le sort de chaque etat du modele COMPLET : ``dynamique``,
        ``statique`` ou ``fige``."""
        base = getattr(cls, "_classe_complete", cls)
        st, fi = set(getattr(cls, "_statiques", ())), set(getattr(cls, "_figes", ()))
        return {
            n: "statique" if n in st else "fige" if n in fi else "dynamique"
            for n in base.states
        }

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.name} @ bus {self.buses}>"


_VARIANTES: dict[tuple, type] = {}


def _variante(base: type, statiques: tuple[str, ...], figes: tuple[str, ...]) -> type:
    """La sous-classe de ``base`` ou ``statiques`` et ``figes`` n'ont plus de
    dynamique. Mise en cache : deux elements reduits de la meme facon partagent
    la meme classe, donc le meme modele compile."""
    cle = (base, tuple(sorted(statiques)), tuple(sorted(figes)))
    if cle in _VARIANTES:
        return _VARIANTES[cle]
    retire = set(statiques) | set(figes)
    garde = [i for i, n in enumerate(base.states) if n not in retire]
    # l'ordre du modele complet, pas celui des appels : la variante ne depend
    # que de l'ENSEMBLE des dynamiques retirees
    i_stat = [i for i, n in enumerate(base.states) if n in statiques]
    n_fige = tuple(n for n in base.states if n in figes)
    equations_completes = base.equations

    def equations(s, _eq=equations_completes, _g=tuple(garde), _s=tuple(i_stat)):
        f, g, h = _eq(s)
        return [f[i] for i in _g], list(g) + [f[i] for i in _s], h

    morceaux = []
    if i_stat:
        morceaux.append("statique:" + ",".join(base.states[i] for i in i_stat))
    if n_fige:
        morceaux.append("fige:" + ",".join(n_fige))
    nom_modele = f"{base.model_name or base.__name__}[{'|'.join(morceaux)}]"

    def initialize(self, ctx, _base=base, _figes=n_fige):
        init = _base.initialize(self, ctx)
        imposees = getattr(self, "_valeurs_figees", {})
        for n in _figes:
            v = imposees.get(n, init[n])
            self.set_param(n, v)
            init[n] = v
        return init

    cls = type(
        base.__name__,
        (base,),
        {
            "__doc__": base.__doc__,
            "__module__": base.__module__,
            "model_name": nom_modele,
            "states": tuple(base.states[i] for i in garde),
            "algebs": tuple(base.algebs) + tuple(base.states[i] for i in i_stat),
            "params": tuple(base.params) + n_fige,
            "equations": staticmethod(equations),
            "initialize": initialize,
            "_classe_complete": base,
            "_statiques": tuple(base.states[i] for i in i_stat),
            "_figes": n_fige,
        },
    )
    _VARIANTES[cle] = cls
    return cls
