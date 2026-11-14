"""Classer une instabilite (un mode) par type de stabilite.

La classification suit celle de Hatziargyriou et al., « Definition and
Classification of Power System Stability -- Revisited & Extended », IEEE
TPWRS 36(4), 2021 (le Task Force IEEE/PES). Elle est PHYSIQUE : une valeur
propre ne porte pas son etiquette, on la construit par le diagnostic
« mode -> etats dominants -> mecanisme -> categorie ».

1. Chaque etat appartient a un GROUPE physique (rotor d'une machine,
   regulateur de vitesse, excitation, PLL d'un GFL, boucle de courant...).
2. ARBRE DE DECISION (:func:`decider`), sur la part de participation des
   groupes, la frequence du mode et la part commune de son energie
   cinetique :

   * les COMMANDES de convertisseurs portent une part notable du mode
     (>= ``SEUIL_CONVERTISSEUR``) -> stabilite pilotee par convertisseur,
     lente si ``f < 10 Hz``, rapide au-dela (§B5 : c'est la FREQUENCE du
     phenomene qui separe les deux, pas la boucle qui y participe -- une PLL
     peut porter une interaction rapide) ;
   * sinon, le mode est porte par les machines, les charges ou le seul
     reseau : angle (rotor, PSS, amortisseurs ; local au-dessus de 1 Hz,
     inter-zones en dessous), frequence (regulateurs de vitesse, ou mode
     commun : toutes les machines ensemble), tension (excitation, champ,
     charges) -- ou, quand le reseau et les stators dominent sans commande
     de convertisseur ni machine, resonance ELECTRIQUE (§B4.2).

   La boucle de puissance du GFM (statisme P-f, filtre de puissance) est
   comptee avec les regulateurs de vitesse des machines : elle gere la
   puissance comme eux. Elle n'est donc pas une « commande » au sens du
   seuil ci-dessus.

   La resonance TORSIONNELLE (§B4.1) demande les modeles d'arbre, absents
   ici ; l'instabilite de tension a long terme, les regleurs en charge, et
   l'harmonique, le decoupage : aucune de ces dynamiques n'est modelisee, et
   un mode ne peut reveler que les etats qu'on a ecrits.
3. Le type SECOND est garde (un mode peut etre hybride : un convertisseur
   qui degrade un mode electromecanique).

::

    from g2egrid.analyse.classement import classer
    c = classer(participations, etats, f)   # {"categorie": ..., "parts": ...}
"""

from __future__ import annotations

import re

import numpy as np

#: les categories, dans l'ordre d'affichage
CATEGORIES = (
    "angle",
    "frequence",
    "tension",
    "convertisseur_lent",
    "convertisseur_rapide",
    "resonance",
)
LIBELLES = {
    "angle": "angle (électromécanique, synchronisation)",
    "frequence": "fréquence (réglage de vitesse)",
    "tension": "tension (excitation, réglage de tension, charges)",
    "convertisseur_lent": "pilotée par convertisseur, interaction lente (< 10 Hz)",
    "convertisseur_rapide": "pilotée par convertisseur, interaction rapide (> 10 Hz)",
    "resonance": "résonance électrique (réseau, flux statoriques, sans commande)",
    "point_de_fonctionnement": "pas de point de fonctionnement (tension statique)",
}

#: etat d'un element -> groupe, par famille d'element
_SM = {
    "Dwr": "SM_rotor",
    "theta": "SM_rotor",
    "Pm": "SM_gouverneur",
    "Dw1": "SM_pss",
    "v1": "SM_pss",
    "v2": "SM_pss",
    "vpss": "SM_pss",
    "e1": "SM_excitation",
    "e2": "SM_excitation",
    "efd": "SM_excitation",
    "e3": "SM_excitation",
    "phi_fd": "SM_champ",
    "phi_1d": "SM_amortisseurs",
    "phi_1q": "SM_amortisseurs",
    "phi_2q": "SM_amortisseurs",
    "phi_d": "SM_stator",
    "phi_q": "SM_stator",
    "i_gd": "reseau",
    "i_gq": "reseau",
}
_GFM = {
    "theta": "GFM_synchronisation",
    "pm": "GFM_synchronisation",
    "Dw": "GFM_synchronisation",
    "qm": "GFM_tension",
    "M_VLd": "GFM_tension",
    "M_VLq": "GFM_tension",
    "ved": "GFM_tension",
    "veq": "GFM_tension",
    "M_CLd": "GFM_courant",
    "M_CLq": "GFM_courant",
    "isd": "GFM_courant",
    "isq": "GFM_courant",
    "igd": "GFM_courant",
    "igq": "GFM_courant",
    "vdc": "GFM_continu",
    "idc": "GFM_continu",
}
_GFL = {
    "M_pll": "GFL_pll",
    "theta_pll": "GFL_pll",
    "M_d": "GFL_puissance_active",
    "M_q": "GFL_puissance_reactive",
    "M_CLd": "GFL_courant",
    "M_CLq": "GFL_courant",
    "isd": "GFL_courant",
    "isq": "GFL_courant",
    "igd": "GFL_courant",
    "igq": "GFL_courant",
    "ved": "GFL_courant",
    "veq": "GFL_courant",
    "vdc": "GFL_continu",
    "idc": "GFL_continu",
}
_CHARGE = {"dp_f": "charge", "dq_f": "charge"}

CATEGORIE_DE = {
    "SM_rotor": "angle",
    "SM_pss": "angle",
    "SM_amortisseurs": "angle",
    # boucle de puissance du GFM : l'equivalent du regulateur de vitesse des
    # machines (elle gere la puissance), donc de la frequence
    "GFM_synchronisation": "frequence",
    "SM_gouverneur": "frequence",
    "SM_excitation": "tension",
    "SM_champ": "tension",
    "GFM_tension": "tension",
    "GFL_puissance_reactive": "tension",
    "charge": "tension",
    "GFL_pll": "convertisseur_lent",
    "GFL_puissance_active": "convertisseur_lent",
    "GFM_courant": "convertisseur_rapide",
    "GFL_courant": "convertisseur_rapide",
    "GFM_continu": "convertisseur_rapide",
    "GFL_continu": "convertisseur_rapide",
    "SM_stator": "resonance",
    "reseau": "resonance",
}


def groupe(etat: str) -> str:
    """Le groupe physique d'un etat ``"ELEMENT.etat"``."""
    nom, _, x = etat.partition(".")
    famille = re.match(r"[A-Za-z]+", nom)
    famille = famille.group(0) if famille else nom
    if famille in ("SM", "SynchronousMachine"):
        return _SM.get(x, "SM_rotor")
    if famille.startswith("GFM"):
        return _GFM.get(x, "GFM_synchronisation")
    if famille.startswith("GFL"):
        return _GFL.get(x, "GFL_courant")
    if famille.startswith("Load") and x in _CHARGE:
        return _CHARGE[x]
    return "reseau"


def parts_par_groupe(p: np.ndarray, etats: list[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    for v, e in zip(p, etats):
        g = groupe(e)
        out[g] = out.get(g, 0.0) + float(v)
    return out


def part_commune(
    v: np.ndarray, etats: list[str], inerties: dict[str, float]
) -> float | None:
    """La part COMMUNE de l'energie cinetique du mode de vecteur propre
    ``v`` : ``|sum M_k v_k|^2 / (sum M_k * sum M_k |v_k|^2)``, sur les angles
    des machines et des GFM (``inerties`` : nom d'element -> M = H.S). Elle
    vaut 1 quand toutes les unites oscillent ensemble (un mode de
    FREQUENCE), 0 quand elles oscillent les unes contre les autres (un mode
    d'ANGLE)."""
    num, m_tot, e_tot = 0j, 0.0, 0.0
    for x, e in zip(v, etats):
        nom, _, etat = e.partition(".")
        if etat != "theta" or nom not in inerties:
            continue
        m = inerties[nom]
        num += m * x
        m_tot += m
        e_tot += m * abs(x) ** 2
    if m_tot == 0 or e_tot == 0:
        return None
    return float(abs(num) ** 2 / (m_tot * e_tot))


#: au-dela, un mode electromecanique commun est un mode de frequence
SEUIL_COMMUN = 0.5
#: part des commandes de convertisseurs au-dela de laquelle le mode est
#: pilote par convertisseur
SEUIL_CONVERTISSEUR = 0.2
#: frontiere entre interactions lente et rapide (§B5 : « typiquement moins
#: de 10 Hz » contre « dizaines a centaines de Hz, voire kHz »)
F_LENTE = 10.0
#: en dessous, un mode d'angle oscillant est inter-zones (0,1-1 Hz), au-dessus
#: local (1-2 Hz)
F_INTER_ZONES = 1.0
#: les groupes d'etats des COMMANDES de convertisseurs. La boucle de
#: puissance du GFM n'en fait pas partie : elle gere la puissance comme le
#: regulateur de vitesse d'une machine, et compte avec lui (frequence)
GROUPES_COMMANDE = {
    "GFM_tension",
    "GFM_courant",
    "GFM_continu",
    "GFL_pll",
    "GFL_puissance_active",
    "GFL_puissance_reactive",
    "GFL_courant",
    "GFL_continu",
}
GROUPES_ANGLE = {"SM_rotor", "SM_pss", "SM_amortisseurs"}
GROUPES_FREQUENCE = {"SM_gouverneur", "GFM_synchronisation"}
GROUPES_TENSION = {"SM_excitation", "SM_champ", "charge"}
GROUPES_PASSIFS = {"reseau", "SM_stator"}


def decider(
    groupes: dict[str, float],
    f: float,
    commun: float | None = None,
    seuil_convertisseur: float = SEUIL_CONVERTISSEUR,
    f_lente: float = F_LENTE,
) -> dict:
    """L'arbre de decision de la classification 2021, a partir des parts de
    participation PAR GROUPE (normalisees ou non), de la frequence ``f`` du
    mode [Hz] et de sa part commune (:func:`part_commune`).

    Rend la categorie, la sous-categorie, la seconde categorie et les
    attributs qui ont decide (part des commandes, nature apériodique ou
    oscillante). Ne demande que ce que la campagne enregistre deja : on peut
    reclasser une campagne sans la recalculer. ``seuil_convertisseur`` et
    ``f_lente`` sont des frontieres OPERATIONNELLES, pas des constantes
    physiques : les faire varier mesure la robustesse du classement."""
    tot = sum(groupes.values()) or 1.0
    g = {k: v / tot for k, v in groupes.items()}
    conv = sum(v for k, v in g.items() if k in GROUPES_COMMANDE)
    ang = sum(v for k, v in g.items() if k in GROUPES_ANGLE)
    freq = sum(v for k, v in g.items() if k in GROUPES_FREQUENCE)
    tens = sum(v for k, v in g.items() if k in GROUPES_TENSION)
    pas = sum(v for k, v in g.items() if k in GROUPES_PASSIFS)
    if commun is not None and commun >= SEUIL_COMMUN:
        # la participation rotor d'un mode commun est celle du centre
        # d'inertie : de la frequence, pas de l'angle
        freq, ang = freq + ang, 0.0
    nature = "aperiodique" if f < 0.01 else "oscillant"
    lent = f < f_lente
    conv_cat = "convertisseur_lent" if lent else "convertisseur_rapide"
    # les candidats et leur poids : la commande compte pour la categorie
    # convertisseur de la bande du mode
    poids = {
        conv_cat: conv,
        "angle": ang,
        "frequence": freq,
        "tension": tens,
        "resonance": pas,
    }
    if conv >= seuil_convertisseur:
        cat = conv_cat
    else:
        # sans commande notable : le mecanisme electromecanique, de tension
        # ou de frequence qui domine ; le reseau seul fait une resonance
        machines = {k: poids[k] for k in ("angle", "frequence", "tension")}
        meca = max(machines, key=machines.get)
        cat = meca if machines[meca] >= pas or machines[meca] >= 0.2 else "resonance"
    sous = None
    if cat == "angle":
        sous = (
            "inter_zones"
            if 0.01 <= f < F_INTER_ZONES
            else ("aperiodique" if f < 0.01 else "local")
        )
    elif cat == "frequence":
        sous = "court_terme"
    elif cat == "tension":
        sous = "court_terme"
    elif cat == "resonance":
        sous = "electrique"
    reste = sorted(
        ((k, v) for k, v in poids.items() if k != cat), key=lambda kv: -kv[1]
    )
    seconde = reste[0][0] if reste and reste[0][1] >= 0.2 else None
    return {
        "categorie": cat,
        "sous_categorie": sous,
        "seconde": seconde,
        "nature": nature,
        "part_commande": round(conv, 4),
        "parts": {k: round(v, 4) for k, v in poids.items() if v >= 0.01},
    }


def classer(
    p: np.ndarray,
    etats: list[str],
    f: float | None = None,
    commun: float | None = None,
) -> dict:
    """La categorie du mode de participations ``p`` (sur ``etats``) et de
    frequence ``f`` [Hz], selon :func:`decider`, avec les parts par groupe.
    ``commun`` (voir :func:`part_commune`) fait d'un mode electromecanique
    ou toutes les unites oscillent ensemble un mode de frequence."""
    g = parts_par_groupe(p, etats)
    tot = sum(g.values()) or 1.0
    d = decider(g, 0.0 if f is None else float(f), commun)
    return {
        "commun": None if commun is None else round(commun, 4),
        **d,
        "groupes": {
            k: round(v / tot, 4)
            for k, v in sorted(g.items(), key=lambda kv: -kv[1])
            if v / tot >= 0.02
        },
    }


# ------------------------------------------------------- procedure complete
#: les composants LENTS (electromecanique, boucles externes, PLL) et RAPIDES
#: (reseau, stators, boucles internes) : la seconde lecture lente/rapide
GROUPES_LENTS_CONVERTISSEUR = {
    "GFL_pll",
    "GFL_puissance_active",
    "GFL_puissance_reactive",
    "GFL_continu",
    "GFM_continu",
    "GFM_synchronisation",
    "SM_rotor",
    "SM_amortisseurs",
    "SM_gouverneur",
    "SM_pss",
    "SM_excitation",
    "SM_champ",
    "charge",
}
GROUPES_RAPIDES_CONVERTISSEUR = {
    "GFL_courant",
    "GFM_courant",
    "GFM_tension",
    "reseau",
    "SM_stator",
}
#: les groupes d'etats par FAMILLE de la procedure (le reseau et les stators
#: forment le support passif)
FAMILLES_PROCEDURE = {
    "commandes": GROUPES_COMMANDE,
    "electromecanique": {
        "SM_rotor",
        "SM_pss",
        "SM_amortisseurs",
        "SM_gouverneur",
        "GFM_synchronisation",
    },
    "tension": GROUPES_TENSION,
    "passif": GROUPES_PASSIFS,
}
#: les commandes portent le mode quand elles en portent au moins cette part
TAU_P = 0.5
#: le support passif porte une resonance quand il porte au moins cette part
TAU_RES = 0.5
#: le mecanisme de chaque groupe de commande dominant
MECANISME_COMMANDE = {
    "GFL_pll": "PLL-réseau",
    "GFL_courant": "boucle de courant-réseau",
    "GFM_courant": "boucle de courant-réseau",
    "GFM_tension": "boucle de tension-réseau",
    "GFL_puissance_active": "boucles de puissance",
    "GFL_puissance_reactive": "boucles de puissance",
    "GFM_continu": "bus continu",
    "GFL_continu": "bus continu",
}
#: le nom court du support de chaque groupe d'etats
SUPPORT = {
    "SM_rotor": "rotor",
    "SM_pss": "PSS",
    "SM_amortisseurs": "amortisseurs",
    "SM_gouverneur": "turbine",
    "SM_excitation": "excitation",
    "SM_champ": "champ",
    "SM_stator": "stator",
    "GFM_synchronisation": "angle GFM",
    "GFM_tension": "tension GFM",
    "GFM_courant": "courant GFM",
    "GFM_continu": "continu",
    "GFL_pll": "PLL",
    "GFL_puissance_active": "P GFL",
    "GFL_puissance_reactive": "Q GFL",
    "GFL_courant": "courant GFL",
    "GFL_continu": "continu",
    "charge": "charges",
    "reseau": "réseau",
}


def _parts(groupes: dict[str, float]) -> dict[str, float]:
    tot = sum(groupes.values()) or 1.0
    return {
        f: sum(v for k, v in groupes.items() if k in g) / tot
        for f, g in FAMILLES_PROCEDURE.items()
    }


def procedure(
    groupes: dict[str, float],
    f: float,
    commun: float | None = None,
    sensibilites: dict[str, float] | None = None,
    consequence: dict[str, float] | None = None,
    regime: str = "linéaire",
    tau_p: float = TAU_P,
    tau_res: float = TAU_RES,
    f_lente: float = F_LENTE,
) -> dict:
    """La procedure de classement complete d'un mode instable : classe
    principale, sous-classe, mecanisme, support, consequence, sensibilite et
    regime. Elle donne la priorite au MECANISME pour les deux classes de 2021
    definies par lui (pilotee par convertisseur, resonance), puis a la
    CONSEQUENCE systeme pour les trois classes historiques :

    A. mode de perte : oscillant (paire complexe) ou aperiodique ;
    B. les commandes de convertisseurs reglent-elles la marge (famille
       d'equipements de plus grande sensibilite, le reseau etant le support
       de toutes les interactions), ou portent-elles la majorite du mode
       (``tau_p``) ? -> pilotee par convertisseur, lente ou rapide ;
    C. un mode oscillant porte par le reseau et les stators (``tau_res``) ->
       resonance electrique, avec sa paire d'elements en interaction ;
    D. la famille electromecanique domine -> frequence si les unites bougent
       ensemble (``commun``), angle sinon (synchronisme d'un GFM quand son
       angle porte le mode) ;
    E. la famille de tension domine -> tension ;

    sinon le mode est hybride.

    ``sensibilites`` : famille -> sensibilite (``commandes``, ``angle``,
    ``frequence``, ``tension``, ``reseau``, ``charge``) ; sans elles, l'etape
    B ne lit que la participation (seuil de 20 %, comme :func:`decider`)."""
    q = _parts(groupes)
    tot = sum(groupes.values()) or 1.0
    g = {k: v / tot for k, v in groupes.items()}
    aperiodique = f < 0.01
    perte = (
        "perte d'équilibre ou de synchronisme"
        if aperiodique
        else "amortissement négatif"
    )
    support = [
        SUPPORT.get(k, k) for k, _ in sorted(g.items(), key=lambda kv: -kv[1])[:2]
    ]
    s = None
    if sensibilites:
        fam = {
            k: sensibilites.get(k, 0.0)
            for k in ("commandes", "angle", "frequence", "tension", "reseau")
        }
        tot_s = sum(fam.values())
        s = {k: v / tot_s for k, v in fam.items()} if tot_s > 0 else None
    cons = max(consequence, key=consequence.get) if consequence else None
    out = {
        "nature": "apériodique" if aperiodique else "oscillant",
        "perte": perte,
        "support": support,
        "consequence": cons,
        "sensibilite": max(s, key=s.get) if s else None,
        "regime": regime,
        "parts": {k: round(v, 4) for k, v in q.items()},
    }
    # B. les commandes. Elles pilotent le mode quand elles en portent la
    # majorite (``tau_p``), quand elles reglent la marge plus que toute autre
    # famille, reseau compris, ou quand elles y prennent part (20 %) et
    # dominent la sensibilite parmi les seules familles d'EQUIPEMENTS : le
    # reseau est le support de toutes les interactions, et la PLL d'un
    # reseau faible reste une interaction de commande meme si la force du
    # reseau regle sa marge
    eq = (
        {
            k: sensibilites.get(k, 0.0)
            for k in ("commandes", "angle", "frequence", "tension")
        }
        if sensibilites
        else {}
    )
    ctrl_eq = sum(eq.values()) > 0 and max(eq, key=eq.get) == "commandes"
    ctrl_tout = s is not None and max(s, key=s.get) == "commandes"
    ctrl_s = ctrl_tout or (ctrl_eq and q["commandes"] >= SEUIL_CONVERTISSEUR)
    ctrl_p = q["commandes"] >= (tau_p if s is not None else SEUIL_CONVERTISSEUR)
    if ctrl_s or ctrl_p:
        lent_f = f < f_lente
        lents = sum(v for k, v in g.items() if k in GROUPES_LENTS_CONVERTISSEUR)
        rapides = sum(v for k, v in g.items() if k in GROUPES_RAPIDES_CONVERTISSEUR)
        lent_c = lents > rapides
        commande = max(
            (k for k in g if k in GROUPES_COMMANDE), key=lambda k: g[k], default=None
        )
        return out | {
            "classe": "convertisseur_lent" if lent_f else "convertisseur_rapide",
            "sous_classe": "hybride lent/rapide"
            if lent_f != lent_c
            else ("lente" if lent_f else "rapide"),
            "mecanisme": MECANISME_COMMANDE.get(commande, "commandes"),
            "critere": "sensibilité" if ctrl_s else "participation",
        }
    # C. la resonance : un echange oscillant porte par le reseau
    if not aperiodique and q["passif"] >= tau_res:
        autres = sorted(
            ((k, v) for k, v in g.items() if k not in GROUPES_PASSIFS),
            key=lambda kv: -kv[1],
        )
        partenaire = SUPPORT.get(autres[0][0], autres[0][0]) if autres else "réseau"
        return out | {
            "classe": "resonance",
            "sous_classe": "électrique",
            "mecanisme": f"réseau-{partenaire}",
            "critere": "participation",
        }
    # D. l'electromecanique : frequence si mouvement commun, angle sinon
    em, te = q["electromecanique"], q["tension"]
    if em >= te and em >= q["passif"]:
        if commun is not None and commun >= SEUIL_COMMUN:
            return out | {
                "classe": "frequence",
                "sous_classe": "court terme",
                "mecanisme": "mouvement commun des unités",
                "critere": "participation",
            }
        gfm = g.get("GFM_synchronisation", 0.0)
        sm = sum(g.get(k, 0.0) for k in ("SM_rotor", "SM_pss", "SM_amortisseurs"))
        if gfm > sm:
            # synchronisme d'un GFM : angle, ou angle/convertisseur quand sa
            # loi de puissance regle la marge plus que le reseau
            loi = s is not None and s["frequence"] > s["reseau"]
            return out | {
                "classe": "angle",
                "sous_classe": "synchronisme GFM",
                "mecanisme": "loi de puissance du GFM"
                if loi
                else "couple synchronisant",
                "seconde": "convertisseur_lent" if loi else None,
                "critere": "participation",
            }
        sous = (
            "apériodique"
            if aperiodique
            else ("inter-zones" if f < F_INTER_ZONES else "local")
        )
        return out | {
            "classe": "angle",
            "sous_classe": sous,
            "mecanisme": "couple synchronisant"
            if aperiodique
            else "couple d'amortissement",
            "critere": "participation",
        }
    # E. la tension
    if te >= q["passif"]:
        return out | {
            "classe": "tension",
            "sous_classe": "court terme",
            "mecanisme": "excitation et charges",
            "critere": "participation",
        }
    return out | {
        "classe": "hybride",
        "sous_classe": None,
        "mecanisme": "indéterminé",
        "critere": "aucun",
    }
