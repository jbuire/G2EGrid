"""Les controles numeriques, le point de fonctionnement (un slack) et le tableau T_eig.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import numpy as np
from g2egrid.analyse import identification, modal
from g2egrid.reseau import construction
from g2egrid.reseau.tables import Tables, tables_vers_json
from g2egrid.temporel import Evenement
from ._commun import JAUGE, ZETA_OSCILLANT
from .assemblage import construire
from .resultats import _domaines_modes, _etats_par_defaut, temporel


# ---------------------------------------------------------------- validation
def validation(
    t: Tables,
    entree: str | None = None,
    amplitude: float = 0.01,
    duree: float = 4.0,
    dt: float = 0.005,
    t_evenement: float = 0.1,
    voies: list[str] | None = None,
    f_max: float = 6.0,
) -> dict:
    """Le modele petit signal decrit-il ce que fait le modele non lineaire ?

    Trois epreuves de cout croissant, qui ne disent pas la meme chose.

    1. **POINT FIXE.** Au point de fonctionnement, ``dx/dt`` doit etre nul et
       les variables algebriques doivent deja verifier leurs equations. Si ce
       n'est pas le cas, les points de fonctionnement par element ne forment
       pas un point fixe du RESEAU : invisible en modal -- on linearise
       quand meme, autour d'un point qui n'existe pas -- et fatal en
       temporel, ou la trajectoire part en derivant sans qu'on ait rien
       perturbe.

    2. **LE CARRE DE L'AMPLITUDE.** L'ecart entre reponse lineaire et
       reponse non lineaire doit croitre comme ``a^2`` : c'est la definition
       meme d'une linearisation au premier ordre. On mesure donc l'ecart a
       deux amplitudes et on en prend le rapport, qui doit valoir 4 quand on
       double. Un rapport proche de 2 dirait que l'ecart est du premier
       ordre, c'est-a-dire que la matrice ``A`` n'est pas la bonne ; un
       rapport qui s'emballe, qu'on est sorti du domaine de validite.

    3. **IDENTIFICATION.** On identifie les modes directement sur la
       trajectoire non lineaire et on les apparie aux valeurs propres de
       ``A``. C'est la seule epreuve qui ait un sens PHYSIQUE : les deux
       chemins ne partagent aucune ligne de code de resolution, si bien que
       leur accord ne peut pas venir d'une erreur commune.

    Le rapport dit aussi si un convertisseur est au-dela de sa limite au
    point de fonctionnement (``limites``, :mod:`g2egrid.noyau.limites`) : la
    simulation, limiteurs surveilles, y demarre en limitation, et le modele
    lineaire -- en mode libre -- ne peut plus la suivre ; les epreuves 2 et 3
    echouent alors pour une raison qui n'est pas une erreur de linearisation.

    ``entree`` vaut par defaut la premiere consigne de puissance venue. Ce
    n'est pas un detail : ``w_ref`` est l'entree globale, la pulsation du
    repere, et un echelon dessus deplace le repere de tout le reseau -- les
    deux voies s'ecartent alors de centaines de % sans que le modele y soit
    pour rien.
    """
    a = construire(t)
    ss = a.etat_stationnaire()
    m = a.non_lineaire()

    if entree is None:
        entree = next(
            (n for n in ss.inputs if n.endswith(".p_ref")),
            next((n for n in ss.inputs if "." in n), ss.inputs[0]),
        )
    if entree not in ss.inputs:
        raise ValueError(f"entree inconnue : {entree}")

    # --- 1. point fixe -------------------------------------------------
    r_x, r_alg = m.residu_equilibre()

    # --- 2. l'ecart croit-il comme a^2 ? -------------------------------
    # Deux amplitudes dans un rapport 2, sur une duree courte : ce qu'on
    # mesure est une PENTE en log, elle n'a pas besoin du regime etabli.
    etats = _etats_par_defaut(ss.states)
    duree_a2 = min(duree, 1.0)
    ecarts = []
    for amp in (amplitude, 2.0 * amplitude):
        d = temporel(
            t,
            entree,
            amp,
            duree=duree_a2,
            n=int(duree_a2 / dt) + 1,
            etats=etats,
            voie="les_deux",
        )
        ecarts.append(d["ecart"]["absolu"])
    rapport = ecarts[1] / ecarts[0] if ecarts[0] > 0 else float("nan")

    # --- 3. identification sur la trajectoire --------------------------
    t_eval = np.arange(0.0, duree + 0.5 * dt, dt)
    tr = m.simuler(
        (0.0, duree),
        uext_fn=m.scenario(Evenement(t_evenement, entree, amplitude)),
        t_eval=t_eval,
    )
    if voies is None:
        voies = [n for n in ss.states if n.endswith((".Dwr", ".Dw", ".w"))][:4]
    if not voies:
        voies = etats[:2]
    idn = identification.sur_trajectoire(tr, voies, evenement=t_evenement)
    # Deux exclusions, pour ne pas reprocher a l'EMT de ne pas retrouver ce
    # qui n'est pas la :
    #  - le mode de JAUGE, qui n'est pas un mode (cf. :data:`JAUGE`) ;
    #  - les modes amortis a plus de ``ZETA_OSCILLANT``. Le filtre
    #    ``oscillants`` ne regarde que la partie imaginaire, si bien qu'il
    #    garde des paires amorties a 99,9 % -- sur wscc9, un mode a 0,005 Hz
    #    porte par l'excitation de SM3. Une telle paire est une exponentielle
    #    decroissante, pas une oscillation : elle ne laisse aucune periode
    #    dans une fenetre de quelques secondes, et aucune identification ne
    #    peut l'y lire. L'exiger ferait echouer la validation sur une
    #    propriete du signal, pas du modele.
    lineaires = [
        x
        for x in modal.modes(ss, f_max=f_max, zeta_max=ZETA_OSCILLANT, oscillants=True)
        if abs(x.valeur) > JAUGE
    ]
    ap = identification.apparier(idn.modes, lineaires)

    return {
        "entree": entree,
        "amplitude": amplitude,
        "point_fixe": {
            "residu_x": r_x,
            "residu_algebrique": r_alg,
            "verdict": bool(r_x < 1e-6 and r_alg < 1e-6),
        },
        "limites": {
            "depassements": [d.en_dict() for d in a.net.depassements],
            "verdict": not a.net.depassements,
        },
        "carre": {
            "ecart_a": ecarts[0],
            "ecart_2a": ecarts[1],
            "rapport": rapport,
            # 4 est la valeur exacte ; on accepte largement, parce que
            # l'ecart a la petite amplitude finit par toucher le bruit de
            # l'integrateur et que le rapport s'en ressent
            "verdict": bool(3.0 <= rapport <= 5.0),
        },
        "identification": {
            "voies": list(voies),
            "rang": int(idn.rang),
            "residu": float(idn.residu),
            "paires": [
                {
                    "f_predit": pr.lineaire.f,
                    "zeta_predit": pr.lineaire.zeta,
                    "f_identifie": pr.identifie.f,
                    "zeta_identifie": pr.identifie.zeta,
                    "ecart_f": pr.ecart_f,
                    "ecart_zeta": pr.ecart_zeta,
                    "ecart_zeta_rel": pr.ecart_zeta_rel,
                    "dominant": pr.lineaire.dominant,
                }
                for pr in sorted(ap.paires, key=lambda x: -x.lineaire.f)
            ],
            "predits_seuls": [
                {"f": x.f, "zeta": x.zeta, "dominant": x.dominant}
                for x in ap.lineaires_seuls
            ],
        },
    }


# ----------------------------------------------------- point de fonctionnement
def un_slack_tout_pq(t: Tables, slack_bus: int) -> dict:
    """Le reseau decrit avec UN slack et tout le reste en PQ.

    C'est la facon la plus courante de poser un point de fonctionnement --
    une reference, et des injections imposees partout ailleurs. Les bus PV
    passent en PQ avec le REACTIF qu'ils fournissaient au load flow courant :
    le point de fonctionnement ne bouge pas, seule change sa description.
    Rend les tables modifiees ; rien n'est ecrit tant que l'appelant ne les
    garde pas.
    """
    a = construire(t)
    idx = a.net._bus_index()
    lf_bus = {int(b): complex(a.lf.S[i]) for b, i in idx.items()}
    if int(slack_bus) not in lf_bus:
        raise ValueError(f"bus {slack_bus} inconnu")
    t2 = t.copie()
    construction.un_slack_tout_pq(t2, int(slack_bus), lf_bus)
    return {"tables": tables_vers_json(t2)}


# ------------------------------------------------------------------- T_eig
def t_eig(t: Tables, n_etats: int = 3) -> dict:
    """Le tableau des modes au format des notebooks (``T_eig``).

    Toutes les valeurs propres physiques -- pas seulement les paires
    oscillantes d'une bande : c'est le tableau de reference, celui qu'on
    exporte et qu'on compare d'une version a l'autre. Colonnes : Re, Im,
    frequence non amortie et amortie, amortissement, et les trois etats les
    plus participants avec leur part. Trie par partie reelle decroissante :
    le moins stable en tete. Les modes de jauge sont retires, et comptes.
    """
    a = construire(t)
    ss = a.etat_stationnaire()
    lam, _, p = modal.spectre(ss)
    p = p / np.maximum(p.sum(axis=1, keepdims=True), 1e-30)
    garde = np.abs(lam) > JAUGE
    noms = list(ss.states)
    domaine, bilan_dom = _domaines_modes(a.net, ss)
    lignes = []
    for k in np.where(garde)[0]:
        top = np.argsort(-p[k])[:n_etats]
        mag = abs(lam[k])
        lignes.append(
            {
                "indice": int(k),
                "re": float(lam[k].real),
                "im": float(lam[k].imag),
                "f_und": float(mag / (2 * np.pi)),
                "f": float(abs(lam[k].imag) / (2 * np.pi)),
                "zeta": float(100 * -lam[k].real / mag),
                "etats": [{"nom": noms[i], "part": float(100 * p[k, i])} for i in top],
                **domaine(k),
            }
        )
    lignes.sort(key=lambda x: -x["re"])
    return {
        "domaines": bilan_dom,
        "n_etats": len(noms),
        "jauge": int((~garde).sum()),
        "max_re": lignes[0]["re"] if lignes else 0.0,
        "stable": bool(lignes and lignes[0]["re"] < 0),
        "modes": lignes,
    }
