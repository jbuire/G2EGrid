"""Validation croisee contre ANDES : le meme reseau, deux outils independants.

POURQUOI -- le jumeau MATLAB V2 valide ce paquet ligne a ligne, mais c'est un
jumeau : meme architecture, memes choix de modelisation, meme auteur. VeraGrid
se valide contre ANDES (CURENT, H. Cui), outil symbolique-numerique de
reference pour le petit signal. On fait de meme : un reseau decrit UNE fois,
dans les tables de ce paquet, est traduit en elements ANDES
(:func:`vers_andes`), et l'on compare ce qui doit coincider.

- Le **load flow** doit coincider exactement : memes donnees, memes equations
  (lignes en pi, transformateurs, charges a puissance constante, shunts).
- Les **modes electromecaniques** doivent coincider a la difference de
  modele pres : ANDES utilise le GENROU, stator et reseau ALGEBRIQUES ; ce
  paquet le modele de Park avec transitoires de stator, et un reseau
  dynamique (lignes R-L, capacites de noeud). Les regulateurs sont
  neutralises des deux cotes (excitation et puissance mecanique constantes,
  pas de stabilisateur) : on compare les machines et le reseau, pas des
  regulateurs qui ne sont pas les memes.

COUPLE OU PUISSANCE -- sans regulateur de vitesse, le GENROU d'ANDES garde
le COUPLE mecanique constant ; la machine de ce paquet, la PUISSANCE
(``Cm = Pm/w``). Linearise, ``Pm/w`` apporte un amortissement ``D = Pm`` (pu
machine) que le couple constant n'a pas : sans correction, les modes locaux
differaient de 0,45 point d'amortissement, frequences identiques. On donne
donc au GENROU ce ``D`` (converti en base systeme, comme ANDES le fait de
ses propres parametres), et l'on compare deux modeles mecaniquement
equivalents. Mesure sur Kundur 12.6 : frequences identiques a 1e-4 pres,
modes locaux a 0,01 point d'amortissement, mode inter-zones a 0,35 point --
le mode qui fait transiter la puissance sur la longue liaison, le plus
sensible a la dynamique du reseau et du stator que seul ce paquet represente.

AVEC REGULATEURS (``regulateurs=True``) -- l'AVR et le PSS de la machine du
paquet ont chacun un homologue ANDES de MEME structure, ce qui permet de les
valider hors du jumeau MATLAB (FR-007) :

- l'AVR (mesure ``Tr``, amplificateur ``Ka/(1+sTa)``, excitatrice
  ``Ke/(1+sTe)``, retour de stabilisation ``Kfd.s/(1+sTfd)`` pris sur
  ``efd`` et retranche a l'entree) est l'``EXST1`` d'ANDES a avance-retard
  neutre (``TB = TC = 0``), sans terme de commutation (``KC = 0``), butees
  ecartees. Seul ecart de structure : ANDES n'a pas le pole de l'excitatrice
  ``1/(1+sTe)``, a ``Te = 1e-4`` s (1e4 rad/s, cent fois au-dela de tout
  mode electromecanique) ; ``KA = Ka.Ke`` ;
- le PSS (passe-bas ``T_LP``, washout ``T_HP``, deux avances-retards, gain,
  entree ``w - 1``, sortie ajoutee a la consigne de l'AVR) est l'``IEEEST``
  d'ANDES en mode 1 (ecart de vitesse ``omega - 1``), filtre ``A1 = T_LP``
  (les autres ``A`` nuls : ANDES replie l'etat a constante de temps nulle),
  ``T1..T4`` les avances-retards, ``KS = K_PSS``, ``T5 = T6 = T_HP``.

Le regulateur de vitesse, lui, n'est pas traduit : il doit etre neutralise
(``mp`` immense, voir :data:`AVR_PSS`), sinon :func:`vers_andes` refuse.

F.E.M. DE VITESSE -- le GENROU d'ANDES ecrit son stator sans la vitesse
(``Flux0`` : ``psid = vq + ra.Iq``), ce paquet avec (``wr.phi``). C'est
l'essentiel de l'ecart du mode inter-zones : remplacer ici ``wr.phi`` par
``phi`` le ramene de 0,35 a 0,13 point sans regulateur. L'AVR l'amplifie,
puisque la tension mesuree porte ce terme : 0,2 / 1,05 / 1,5 point pour
``Ka`` = 30 / 100 / 300, et de nouveau 0,1 point une fois le terme retire
(voir ``tests/test_pss_avr.py``). Les modes locaux, eux, restent a moins de
0,3 point avec ou sans regulateurs.

Le reseau de reference est celui de l'exemple 12.6 de Kundur -- deux zones de
deux machines de 900 MVA reliees par une longue liaison, le cas d'ecole des
oscillations inter-zones --, construit d'apres le livre par
:func:`kundur_deux_zones`.

ANDES (GPL-3.0) n'est qu'un outil de VERIFICATION : il n'est ni embarque ni
recopie, s'installe a part (groupe ``validation`` de ``pyproject.toml``), et
ce module ne l'importe qu'a l'appel.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import linear_sum_assignment

from g2egrid.reseau.reglage import machine_depuis_standard, standard_depuis_machine
from g2egrid.reseau.tables import Tables, assembler

__all__ = [
    "AVR_PSS",
    "KUNDUR_12_6",
    "SANS_REGULATEURS",
    "comparer",
    "kundur_deux_zones",
    "vers_andes",
]

#: Kundur, exemple 12.6 : la fiche des quatre machines (base 900 MVA, 20 kV).
KUNDUR_12_6 = dict(
    Xd=1.8, Xq=1.7, Xl=0.2, Xd_p=0.3, Xq_p=0.55, Xd_pp=0.25, Xq_pp=0.25,
    Ra=0.0025, Td0_p=8.0, Tq0_p=0.4, Td0_pp=0.03, Tq0_pp=0.05,
)  # fmt: skip

#: Ce qui neutralise les regulateurs de la machine du paquet : un AVR de gain
#: negligeable (excitation constante), un statisme immense (puissance
#: mecanique constante), pas de PSS ni d'amortissement ; la resistance de
#: causalite aux bornes rendue negligeable.
SANS_REGULATEURS = dict(Ka=1e-6, mp=1e6, K_PSS=0.0, KD=0.0, Rg=1e4)

#: Ce que la traduction ANDES reproduit AVEC regulateurs : l'AVR et le PSS de
#: :data:`~g2egrid.elements.sync_machine.SM_DEFAULT` restent actifs, le
#: regulateur de vitesse et l'amortissement mecanique sont neutralises comme
#: dans :data:`SANS_REGULATEURS`. ``K_PSS`` se surcharge a la demande.
AVR_PSS = {"mp": 1e6, "KD": 0.0, "Rg": 1e4}


def kundur_deux_zones(
    avec_regulateurs: bool = False, reglage: dict | None = None
) -> Tables:
    """Le reseau a deux zones de Kundur (exemple 12.6), d'apres le livre.

    Base 100 MVA ; lignes de 230 kV a r = 0,0001, x = 0,001, b = 0,00175 pu/km ;
    transfos 20/230 kV de x = 0,15 pu sur 900 MVA ; charges et condensateurs
    aux bus 7 et 9. Un petit condensateur (0,01 pu) a chaque bus de machine
    lui donne une capacite : sans elle, le noeud serait resistif, et sa
    resistance de causalite n'a pas d'equivalent dans ANDES.

    ``reglage`` : parametres communs aux quatre machines, poses par-dessus la
    fiche 12.6 (par exemple :data:`AVR_PSS`). Il l'emporte sur
    ``avec_regulateurs``.
    """
    Sb, f0 = 100e6, 60.0
    Y_DER = [
        (1, 1, 20e3, 1.03, math.radians(20.2), 700e6, 0, 0, 0),
        (2, 2, 20e3, 1.01, 0, 700e6, 0, 0, 0),
        (3, 2, 20e3, 1.03, 0, 719e6, 0, 0, 0),
        (4, 2, 20e3, 1.01, 0, 700e6, 0, 0, 0),
        *[(b, 3, 230e3, 1.0, 0, 0, 0, 0, 0) for b in (5, 6, 8, 10, 11)],
        (7, 3, 230e3, 1.0, 0, 0, 0, 967e6, 100e6),
        (9, 3, 230e3, 1.0, 0, 0, 0, 1767e6, 100e6),
    ]
    Y_DER.sort(key=lambda r: r[0])
    lignes = [(5, 6, 25), (6, 7, 10), (7, 8, 110), (7, 8, 110), (8, 9, 110), (8, 9, 110),
              (9, 10, 10), (10, 11, 25)]  # fmt: skip
    T_Lines_pu = [
        (k + 1, i, j, km, 1e-4 * km, 1e-3 * km, 1.75e-3 * km, 230e3, Sb)
        for k, (i, j, km) in enumerate(lignes)
    ]
    T_Transfos = [
        (k + 1, g, b, 900e6, 20e3, 230e3, 0.0, 0.15)
        for k, (g, b) in enumerate(((1, 5), (2, 6), (3, 11), (4, 10)))
    ]
    T_Shunts = [(1, 7, 2.0), (2, 9, 3.5), *[(2 + g, g, 0.01) for g in (1, 2, 3, 4)]]
    T_SMs = [(g, g, 900e6) for g in (1, 2, 3, 4)]
    park = machine_depuis_standard(KUNDUR_12_6, f0)
    if reglage is None:
        reglage = {} if avec_regulateurs else SANS_REGULATEURS
    par_bus = {}
    for g, H in ((1, 6.5), (2, 6.5), (3, 6.175), (4, 6.175)):
        par_bus[g] = {**park, "H": H, **reglage}
    return Tables(
        Y_DER=Y_DER,
        T_Lines_pu=T_Lines_pu,
        T_Transfos=T_Transfos,
        T_Shunts=T_Shunts,
        T_SMs=T_SMs,
        PARAMS_PAR_BUS=par_bus,
        f0=f0,
        Sb=Sb,
        nom="kundur_12_6",
    )


def vers_andes(t: Tables, D: dict[int, float] | None = None, regulateurs: bool = False):
    """Le meme reseau en elements ANDES : Bus, Line, PQ, Shunt, PV/Slack, GENROU.

    Couvre ce dont la validation a besoin -- ``Y_DER`` (bus, charges,
    consignes), ``T_Lines_pu``, ``T_Lines_SI``, ``T_Transfos``, ``T_Shunts``
    et ``T_SMs`` (machine convertie en fiche standard par
    :func:`g2egrid.reseau.reglage.standard_depuis_machine`). Les charges deviennent
    des impedances constantes en dynamique (``p2z = q2z = 1``), comme les
    charges R-L du paquet.

    ``D`` : amortissement mecanique de chaque machine, par bus, en pu de la
    machine (voir l'en-tete : l'equivalent du couple a puissance constante).

    ``regulateurs`` : ajoute a chaque machine son AVR (``EXST1``) et, si
    ``K_PSS`` n'est pas nul, son PSS (``IEEEST``) -- correspondance dans
    l'en-tete. Refuse une machine dont le regulateur de vitesse est actif
    ou dont l'excitatrice a une constante de temps non negligeable : la
    traduction ne les reproduit pas.

    EFFET DE BORD D'ANDES -- voir :func:`importer_andes`. La classe
    ``Piecewise`` d'ANDES n'est remise en place que le temps de construire
    son systeme (c'est la qu'il genere son code).
    """
    import sympy

    andes = importer_andes()
    from andes.thirdparty.sympymod import FixPiecewise

    origine = sympy.Piecewise
    sympy.Piecewise = FixPiecewise
    try:
        return _systeme_andes(andes, t, D, regulateurs)
    finally:
        sympy.Piecewise = origine


def importer_andes():
    """``import andes``, sans laisser SymPy modifie apres coup.

    A son import, ANDES remplace ``sympy.Piecewise`` par sa propre classe
    (``andes/core/symprocessor.py``), dont le code lambdifie appelle
    ``__trues``, ``__zeros``... des noms qu'ANDES est seul a fournir. Tout
    modele de ce paquet compile ensuite avec un ``Piecewise`` echouerait
    (``NameError: name '__trues' is not defined``), et resterait casse dans
    le cache de compilation. On rend donc a SymPy sa classe des la fin de
    l'import -- y compris quand ANDES avait deja ete importe ailleurs, d'ou
    ``OldPiecewise``, ou ANDES range l'original. Tout import d'ANDES de ce
    paquet et de ses tests passe par ici.
    """
    import sympy

    origine = getattr(sympy, "OldPiecewise", sympy.Piecewise)
    import andes

    sympy.Piecewise = getattr(sympy, "OldPiecewise", origine)
    return andes


def _systeme_andes(andes, t: Tables, D, regulateurs: bool):
    """Corps de :func:`vers_andes`, ``sympy.Piecewise`` a la mode d'ANDES."""
    andes.config_logger(stream_level=40)
    from g2egrid.elements.sync_machine import SM_DEFAULT

    ss = andes.System(default_config=True, no_output=True)
    ss.config.mva = t.Sb / 1e6
    ss.config.freq = t.f0
    UN = {int(r[0]): float(r[2]) for r in t.Y_DER}
    for bid, _ty, Un, V, th, *_ in t.Y_DER:
        ss.add(
            "Bus", dict(idx=int(bid), name=f"bus{int(bid)}", Vn=Un / 1e3, v0=V, a0=th)
        )
    for lid, i, j, _km, r, x, b, Un, Sb_l in t.T_Lines_pu:
        ss.add("Line", dict(idx=f"L{int(lid)}", bus1=int(i), bus2=int(j), r=r, x=x, b=b,
                            Sn=Sb_l / 1e6, Vn1=Un / 1e3, Vn2=Un / 1e3))  # fmt: skip
    for lid, i, j, km, r_km, x_km, b_km in t.T_Lines_SI:
        Zb = UN[int(i)] ** 2 / t.Sb
        ss.add("Line", dict(idx=f"LS{int(lid)}", bus1=int(i), bus2=int(j),
                            r=r_km * km / Zb, x=x_km * km / Zb, b=b_km * km * Zb,
                            Sn=t.Sb / 1e6, Vn1=UN[int(i)] / 1e3, Vn2=UN[int(j)] / 1e3))  # fmt: skip
    for tid, i, j, Sn, U1, U2, R, X in t.T_Transfos:
        ss.add("Line", dict(idx=f"T{int(tid)}", bus1=int(i), bus2=int(j), r=R, x=X, b=0.0,
                            Sn=Sn / 1e6, Vn1=U1 / 1e3, Vn2=U2 / 1e3, trans=1, tap=1.0))  # fmt: skip
    for sid, b, B in t.T_Shunts:
        ss.add("Shunt", dict(idx=f"S{int(sid)}", bus=int(b), b=B, g=0.0,
                             Sn=t.Sb / 1e6, Vn=UN[int(b)] / 1e3))  # fmt: skip
    for bid, _ty, _Un, _V, _th, _Pg, _Qg, Pl, Ql in t.Y_DER:
        if Pl or Ql:
            ss.add(
                "PQ",
                dict(idx=f"PQ{int(bid)}", bus=int(bid), p0=Pl / t.Sb, q0=Ql / t.Sb),
            )
    ss.PQ.config.p2z = ss.PQ.config.q2z = 1.0
    ss.PQ.config.p2p = ss.PQ.config.q2q = 0.0
    for bid, ty, _Un, V, th, Pg, *_ in t.Y_DER:
        machines = [r for r in t.T_SMs if int(r[1]) == int(bid)]
        if not machines:
            continue
        (_id, b, Sn) = machines[0][:3]
        gen = f"G{int(b)}"
        if int(ty) == 1:
            ss.add("Slack", dict(idx=gen, bus=int(b), v0=V, a0=th, p0=Pg / t.Sb,
                                 Sn=Sn / 1e6, Vn=UN[int(b)] / 1e3))  # fmt: skip
        else:
            ss.add("PV", dict(idx=gen, bus=int(b), v0=V, p0=Pg / t.Sb,
                              Sn=Sn / 1e6, Vn=UN[int(b)] / 1e3))  # fmt: skip
        p = {**SM_DEFAULT, **t.SM_PARAMS, **t.PARAMS_PAR_BUS.get(int(b), {})}
        s = standard_depuis_machine(p, t.f0)
        # le GENROU suppose X''d = X''q ; sinon ANDES ne s'initialise pas (il le
        # dit, mais poursuit, et ses modes sont pris hors equilibre). C'est le
        # cas de la machine par defaut : X''d = 0,230, X''q = 0,250.
        if abs(s["Xd_pp"] - s["Xq_pp"]) > 1e-6:
            raise ValueError(
                f"machine du bus {int(b)} : X''d = {s['Xd_pp']:.4f} != X''q = "
                f"{s['Xq_pp']:.4f}, le GENROU d'ANDES impose leur egalite"
            )
        ss.add(
            "GENROU",
            dict(
                idx=f"GR{int(b)}", bus=int(b), gen=gen, Sn=Sn / 1e6, Vn=UN[int(b)] / 1e3,
                M=2 * p["H"], D=(D or {}).get(int(b), 0.0), ra=s["Ra"], xl=s["Xl"], xd=s["Xd"], xq=s["Xq"],
                xd1=s["Xd_p"], xq1=s["Xq_p"], xd2=s["Xd_pp"], xq2=s["Xq_pp"],
                Td10=s["Td0_p"], Tq10=s["Tq0_p"], Td20=s["Td0_pp"], Tq20=s["Tq0_pp"],
                S10=0.0, S12=0.0,
            ),
        )  # fmt: skip
        if regulateurs:
            _ajouter_regulateurs(ss, int(b), p)
    ss.setup()
    return ss


#: butees ecartees : le petit signal se linearise hors butee des deux cotes
_LARGE = 99.0


def _ajouter_regulateurs(ss, b: int, p: dict) -> None:
    """``EXST1`` (+ ``IEEEST``) pour la machine du bus ``b`` (voir l'en-tete)."""
    if p["mp"] < 1e3:
        raise ValueError(
            f"machine du bus {b} : regulateur de vitesse actif (mp = {p['mp']:g} %), "
            "non traduit vers ANDES -- le neutraliser (mp = 1e6, voir AVR_PSS)"
        )
    if p["Te"] > 1e-3:
        raise ValueError(
            f"machine du bus {b} : Te = {p['Te']:g} s, l'EXST1 n'a pas de pole "
            "d'excitatrice"
        )
    ss.add(
        "EXST1",
        {
            "idx": f"AVR{b}", "syn": f"GR{b}", "TR": p["Tr"],
            "VIMAX": _LARGE, "VIMIN": -_LARGE, "TC": 0.0, "TB": 0.0,
            "KA": p["Ka"] * p["Ke"], "TA": p["Ta"], "VRMAX": _LARGE,
            "VRMIN": -_LARGE, "KC": 0.0, "KF": p["Kfd"], "TF": p["Tfd"],
        },
    )  # fmt: skip
    if p["K_PSS"] == 0.0:
        return
    ss.add(
        "IEEEST",
        {
            "idx": f"PSS{b}", "avr": f"AVR{b}", "MODE": 1,
            "A1": p["T_LP"], "A2": 0.0, "A3": 0.0, "A4": 0.0, "A5": 0.0, "A6": 0.0,
            "T1": p["T1n"], "T2": p["T1d"], "T3": p["T2n"], "T4": p["T2d"],
            "T5": p["T_HP"], "T6": p["T_HP"], "KS": p["K_PSS"],
            "LSMAX": _LARGE, "LSMIN": -_LARGE,
        },
    )  # fmt: skip


def _modes(lam, fmin, fmax, zmax):
    lam = np.asarray(lam)
    f = np.abs(lam.imag) / (2 * np.pi)
    z = -lam.real / np.maximum(np.abs(lam), 1e-30) * 100
    garde = (lam.imag > 1e-6) & (f >= fmin) & (f <= fmax) & (z <= zmax)
    return lam[garde]


def comparer(
    t: Tables,
    fmin: float = 0.2,
    fmax: float = 3.0,
    zmax: float = 30.0,
    regulateurs: bool = False,
) -> dict:
    """Load flow et modes electromecaniques, ici et dans ANDES.

    Rend l'ecart de tension (module et angle) bus par bus, et les modes de la
    bande ``[fmin, fmax]`` Hz amortis a moins de ``zmax`` %, apparies par
    l'algorithme hongrois : ``(f, zeta)`` des deux cotes. ``regulateurs`` :
    voir :func:`vers_andes`.
    """
    net = assembler(t)
    lf = net.run_loadflow()
    net.solve_equilibrium()
    # l'equivalent en couple de la puissance mecanique constante (en-tete),
    # plus l'amortissement mecanique KD de la machine
    D = {
        int(e.buses[0]): float(e.x0[e.states.index("Pm")])
        + float(e.pdict.get("KD", 0.0))
        for e in net.elements
        if "Pm" in e.states
    }
    ss = vers_andes(t, D=D, regulateurs=regulateurs)
    ss.PFlow.run()
    if not ss.PFlow.converged:
        raise RuntimeError("load flow ANDES non converge")
    idx = net._bus_index()
    ecarts = []
    for k, b in enumerate(ss.Bus.idx.v):
        V = lf.V[idx[int(b)]]
        ecarts.append(
            (int(b), abs(abs(V) - ss.Bus.v.v[k]), abs(np.angle(V) - ss.Bus.a.v[k]))
        )
    ss.EIG.run()
    lam_a = _modes(ss.EIG.mu, fmin, fmax, zmax)
    lam_g = _modes(net.build_state_space().eig(), fmin, fmax, zmax)
    cout = np.abs(lam_a[:, None] - lam_g[None, :]) / np.abs(lam_a[:, None])
    i, j = linear_sum_assignment(cout)

    def fz(lam):
        return float(abs(lam.imag) / (2 * np.pi)), float(-lam.real / abs(lam) * 100)

    paires = [
        {
            "andes": fz(lam_a[a]),
            "g2egrid": fz(lam_g[b]),
            "ecart_relatif": float(cout[a, b]),
        }
        for a, b in zip(i, j)
    ]
    paires.sort(key=lambda p: p["andes"][0])
    return {
        "ecart_V": max(e[1] for e in ecarts),
        "ecart_angle": max(e[2] for e in ecarts),
        "bus": ecarts,
        "modes": paires,
        "n_andes": len(lam_a),
        "n_g2egrid": len(lam_g),
    }
