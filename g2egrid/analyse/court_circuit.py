"""Puissance de court-circuit a chaque noeud, et force du reseau vue des GFL.

DEUX ESTIMATIONS, PARCE QU'ELLES NE REPONDENT PAS A LA MEME QUESTION.

**Statique** (:func:`puissance_court_circuit`) -- un defaut triphase franc, le
reseau ramene a son Y-bus de branches (lignes R + jX et transformateurs,
capacites et charges negligees), les sources de tension derriere leur
impedance. C'est la methode de la CEI 60909-0 (courant de court-circuit
initial symetrique I_k'' par la source de tension equivalente au point de
defaut), avec ses facteurs de correction d'impedance, mais ce N'EST PAS un
calcul normatif complet : pas de c_min, pas de courant de crete ni de courant
coupe, pas de composantes symetriques, et les rapports de tension nominale
sont supposes egaux a 1 (U_n/U_rG = 1, U_rTLV/U_rTHV = rapport du reseau).

* bus infini derriere son ``R + jX`` ;
* machine synchrone derriere sa reactance SOUS-TRANSITOIRE (defaut,
  ``machine="soustransitoire"``)
  ``x''d = Ll + 1/(1/Lad + 1/Lfd + 1/L1d)`` (plus son transfo integre s'il y en
  a un) : c'est la reactance de la CEI 60909 pour I_k''. La reactance
  TRANSITOIRE ``x'd = Ll + Lad.Lfd/(Lad + Lfd)`` (``machine="transitoire"``)
  ignore l'amortisseur et sous-estime Scc pres d'une machine (Scc de 1,01 a 1,21 fois
  plus grande en mediane sur le catalogue avec x''d, jusqu'a 1,40) ; elle
  rend l'ancien calcul, sans facteur de correction ;
* GFM, au choix : source de tension derriere son filtre (``"tension"``, ce
  qu'il est en petits signaux), source de courant limitee a ``Imax.Sn``
  (``"courant"``, ce qu'il est pendant un vrai defaut), ou rien ;
* GFL : rien (``"aucun"``, la convention du SCR) ou source de courant
  ``Imax.Sn`` (``"courant"``, la CEI 60909-0:2016 pour les convertisseurs).

Facteurs de correction d'impedance (CEI 60909-0, avec ``c_max = c``), appliques
avec la reactance sous-transitoire sauf ``correction_cei=False`` :

* generateur raccorde directement au reseau,
  ``K_G = c / (1 + x''d.sin(phi_rG))`` sur son impedance ;
* bloc centrale (machine + transformateur elevateur, seul sur son bus : transfo
  integre de la machine ou ``T_Transfos`` dont le bus bas ne porte que la
  machine), regleur en charge,
  ``K_S = c / (1 + |x''d - x_T|.sin(phi_rG))`` sur l'impedance du bloc, avec
  ``x_T`` la reactance du transfo sur sa propre base ;
* transformateur de reseau a deux enroulements (``T_Transfos`` hors bloc
  centrale), ``K_T = 0,95.c / (1 + 0,6.x_T)``.

``sin(phi_rG)`` vient de ``cos_phi`` (0,9 par defaut : valeur courante, le
catalogue ne porte pas le facteur de puissance nominal des machines). Les
formules ont ete confrontees a l'implementation de pandapower
(``pandapower.shortcircuit``), pas relues dans la norme : « a verifier ». Le
cas du bloc centrale SANS regleur en charge (``K_SO``) n'est pas traite : le
regleur est suppose present. Avec ``c = 1`` (defaut), ``c_max`` de la norme
(1,1 en general) n'est pas celui des facteurs : passer ``c = 1.1`` pour le
calcul normatif de I_k'' maximal.

Avec ``Z = Y^-1``, le courant de defaut au noeud ``i`` est, en pu,

    I_k(i) = ( c + sum_j abs(Z_ij) . I_j ) / abs(Z_ii)

-- le terme des sources de courant ``I_j`` est celui de la norme -- et
``S_cc(i) = I_k(i) . Sb`` (tension prefault 1 pu, ``c = 1`` par defaut).
Le rapport ``X/R`` est celui de ``Z_ii``.

**Vue par l'impedance** (``dynamique=True``) -- ``Sb / sigma_max(Z_dq(j2pi.f))``,
l'impedance dq de TOUT le reseau linearise au point de fonctionnement, vue du
noeud, a basse frequence ``f`` dans le repere dq (c'est-a-dire autour de
50 Hz). Elle compte ce que la statique ignore : les boucles des
convertisseurs, la PLL des GFL, les charges. C'est la force que « voit » un
nouvel equipement, au sens de l'etude ``etudes/accueil_gfl_dc``.

LE SCR. Pour chaque noeud qui porte des GFL (``T_GFL1s``, ``T_GFL2s`` et les
elements a schema-bloc de ``T_Blocs``, des sources de courant), ``SCR =
S_cc / Sn_GFL`` ; en dessous de 3 le reseau est dit faible, en dessous de 2
tres faible (NERC, *Integrating inverter-based resources into low short circuit
strength systems*, 2017). Des GFL VOISINS s'affaiblissent l'un l'autre : le
SCRIF (*SCR with interaction factors*, meme guide) en tient compte noeud par
noeud,

    SCRIF_i = S_cc,i / (P_i + sum_{j != i} IF_ji . P_j),   IF_ji = abs(Z_ji / Z_ii),

``IF_ji`` la variation de tension au noeud ``j`` pour une variation au noeud
``i``, ``P`` la puissance nominale des GFL. On ne rend PAS le WSCR de l'ERCOT
(``sum S_cc,i P_i / (sum P_i)^2``) pour tout le reseau : il est fait pour une
grappe de producteurs proches, et sur un grand reseau a GFL disperses il
divise le SCR par leur nombre.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

import numpy as np

from g2egrid.reseau.tables import Tables, assembler

__all__ = [
    "SccNoeud",
    "RapportScc",
    "MODES_GFM",
    "MODES_GFL",
    "MODES_MACHINE",
    "puissance_court_circuit",
]

#: ce que deviennent les GFM et les GFL dans le calcul statique
MODES_GFM = ("tension", "courant", "aucun")
MODES_GFL = ("aucun", "courant")
#: la reactance des machines synchrones : x''d (CEI 60909, defaut) ou x'd
MODES_MACHINE = ("soustransitoire", "transitoire")
#: seuils usuels du SCR
SCR_FAIBLE, SCR_TRES_FAIBLE = 3.0, 2.0


@dataclass
class SccNoeud:
    """La puissance de court-circuit d'un noeud alternatif."""

    bus: int
    Un: float  # tension nominale [V]
    Scc: float  # puissance de court-circuit [VA]
    Ik: float  # courant de court-circuit [A]
    XR: float  # rapport X/R de l'impedance de Thevenin
    R_th: float  # impedance de Thevenin [pu base Sb]
    X_th: float
    S_gfl: float  # puissance nominale des GFL raccordes au noeud [VA]
    S_gfm: float  # idem des GFM [VA]
    SCR: float | None  # Scc / S_gfl, None sans GFL
    SCRIF: float | None = None  # SCR avec facteurs d'interaction
    Scc_dyn: float | None = None  # vue par l'impedance dq [VA]

    @property
    def faible(self) -> bool:
        return self.SCR is not None and self.SCR < SCR_FAIBLE


@dataclass
class RapportScc:
    noeuds: list[SccNoeud]
    hypotheses: dict

    @property
    def SCR_min(self) -> float | None:
        v = [n.SCR for n in self.noeuds if n.SCR is not None]
        return min(v) if v else None

    @property
    def SCRIF_min(self) -> float | None:
        v = [n.SCRIF for n in self.noeuds if n.SCRIF is not None]
        return min(v) if v else None

    def noeud(self, bus: int) -> SccNoeud:
        return next(n for n in self.noeuds if n.bus == int(bus))

    def vers_json(self) -> dict:
        def propre(x):
            return None if x is None or not math.isfinite(x) else float(x)

        return {
            "noeuds": [
                {
                    k: (propre(v) if isinstance(v, float) else v)
                    for k, v in asdict(n).items()
                }
                for n in self.noeuds
            ],
            "SCR_min": propre(self.SCR_min),
            "SCRIF_min": propre(self.SCRIF_min),
            "hypotheses": self.hypotheses,
            "seuils": {"faible": SCR_FAIBLE, "tres_faible": SCR_TRES_FAIBLE},
        }


def _reglages(t: Tables, famille: dict, defaut: dict, b: int) -> dict:
    return defaut | famille | t.PARAMS_PAR_BUS.get(int(b), {})


def _x_machine(p: dict, machine: str) -> float:
    """La reactance de la machine [pu base machine] : ``x''d`` (CEI 60909, avec
    l'amortisseur de l'axe d) ou ``x'd`` (sans)."""
    if machine == "soustransitoire":
        return p["Ll"] + 1.0 / (1.0 / p["Lad"] + 1.0 / p["Lfd"] + 1.0 / p["L1d"])
    return p["Ll"] + p["Lad"] * p["Lfd"] / (p["Lad"] + p["Lfd"])


def _blocs_centrale(t: Tables) -> dict[int, tuple[int, int]]:
    """Les blocs centrale de ``T_Transfos`` : ``{bus de la machine: (id du
    transfo, bus haut)}``. Un bloc est une machine seule sur son bus (ni autre
    source, ni autre branche) derriere un transformateur."""
    branches: dict[int, list[tuple[str, int, int]]] = {}
    for nom, table in (("L", t.T_Lines_SI), ("L", t.T_Lines_pu), ("T", t.T_Transfos)):
        for r in table:
            i, j = int(r[1]), int(r[2])
            branches.setdefault(i, []).append((nom, int(r[0]), j))
            branches.setdefault(j, []).append((nom, int(r[0]), i))
    sources: dict[int, int] = {}
    for table in (
        t.T_IBs,
        t.T_SMs,
        t.T_GFL1s,
        t.T_GFL2s,
        t.T_GFM1s,
        t.T_GFM2s,
        t.T_GFM3s,
        t.T_GFM4s,
        t.T_Blocs,
    ):
        for r in table:
            sources[int(r[1])] = sources.get(int(r[1]), 0) + 1
    blocs = {}
    for r in t.T_SMs:
        b = int(r[1])
        if sources.get(b) == 1 and len(branches.get(b, [])) == 1:
            nom, tid, haut = branches[b][0]
            if nom == "T":
                blocs[b] = (tid, haut)
    return blocs


def _sources(t: Tables, gfm: str, gfl: str, machine: str, facteurs: dict | None):
    """Les sources du calcul statique : admittances a la terre (sources de
    tension) et courants injectes (sources de courant), en pu base Sb, par
    bus ; et la puissance nominale des GFL et des GFM de chaque bus.

    ``facteurs`` (ou ``None`` : aucune correction) porte ``c``, ``cos_phi`` et,
    en sortie, les facteurs ``K_G``/``K_S`` appliques a chaque machine
    (``facteurs["KG"]``, ``facteurs["KS"]``) ainsi que ``facteurs["kT"]`` des
    transformateurs, rempli par l'appelant."""
    from g2egrid.elements.gfl import GFL_DEFAULT
    from g2egrid.elements.gfm import GFM_DEFAULT
    from g2egrid.elements.sync_machine import SM_DEFAULT

    Sb = t.Sb
    y: dict[int, complex] = {}
    i: dict[int, float] = {}
    s_gfl: dict[int, float] = {}
    s_gfm: dict[int, float] = {}

    def ajoute(d, b, v):
        d[int(b)] = d.get(int(b), 0) + v

    blocs = _blocs_centrale(t) if facteurs is not None else {}
    for r in t.T_IBs:
        b, Sn = int(r[1]), float(r[2])
        ajoute(y, b, (Sn / Sb) / complex(float(r[5]), float(r[6])))
    for _id, b, Sn, *Tr in t.T_SMs:
        p = _reglages(t, t.SM_PARAMS, SM_DEFAULT, b)
        xd = _x_machine(p, machine)
        z = complex(p["Ra"], xd)
        if Tr:
            z += complex(float(Tr[0]), float(Tr[1]))
        if facteurs is not None:
            sin_phi = math.sqrt(1.0 - facteurs["cos_phi"] ** 2)
            c = facteurs["c"]
            if Tr:  # bloc centrale a transfo integre : x_T en base machine
                k = c / (1.0 + abs(xd - float(Tr[1])) * sin_phi)
                facteurs["KS"][int(b)] = k
            elif int(b) in blocs:  # transfo de T_Transfos : x_T sur sa base
                tid = blocs[int(b)][0]
                x_t = next(float(r[7]) for r in t.T_Transfos if int(r[0]) == tid)
                k = c / (1.0 + abs(xd - x_t) * sin_phi)
                facteurs["KS"][int(b)] = k
                facteurs["bloc"][tid] = k  # le transfo du bloc porte le meme K_S
                # defaut AUX BORNES de la machine : K_G sur elle, transfo nu
                facteurs["bornes"][int(b)] = (tid, k, c / (1.0 + xd * sin_phi))
            else:  # generateur direct
                k = c / (1.0 + xd * sin_phi)
                facteurs["KG"][int(b)] = k
            z *= k
        ajoute(y, b, (float(Sn) / Sb) / z)
    for table in ("T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s"):
        for _id, b, Sn, _vdc, *Tr in getattr(t, table):
            ajoute(s_gfm, b, float(Sn))
            p = _reglages(t, t.GFM_PARAMS, GFM_DEFAULT, b)
            if Tr:
                p = p | dict(R2_pu=float(Tr[0]), L2_pu=float(Tr[1]))
            if gfm == "tension":
                z = complex(p["R1_pu"] + p["R2_pu"], p["L1_pu"] + p["L2_pu"])
                ajoute(y, b, (float(Sn) / Sb) / z)
            elif gfm == "courant":
                ajoute(i, b, float(p.get("Imax", 1.2)) * float(Sn) / Sb)
    for table in ("T_GFL1s", "T_GFL2s"):
        for r in getattr(t, table):
            b, Sn = int(r[1]), float(r[2])
            ajoute(s_gfl, b, Sn)
            if gfl == "courant":
                p = _reglages(t, t.GFL_PARAMS, GFL_DEFAULT, b)
                ajoute(i, b, float(p.get("Imax", 1.2)) * Sn / Sb)
    for r in t.T_Blocs:
        b = int(r[1])
        Sn = float(r[3]) if len(r) > 3 and r[3] else Sb
        ajoute(s_gfl, b, Sn)
        if gfl == "courant":
            ajoute(i, b, 1.2 * Sn / Sb)
    return y, i, s_gfl, s_gfm


def _charges(t: Tables) -> dict[int, complex]:
    """Les charges en impedance constante, a 1 pu [pu base Sb]."""
    y: dict[int, complex] = {}
    for r in t.Y_DER:
        s = complex(float(r[7]), float(r[8])) / t.Sb
        if s:
            y[int(r[0])] = y.get(int(r[0]), 0j) + s.conjugate()
    for r in t.T_Loads:
        s = complex(float(r[2]), float(r[3])) / t.Sb
        y[int(r[1])] = y.get(int(r[1]), 0j) + s.conjugate()
    return y


def _scc_dynamique(net, bus: list[int], f: float) -> dict[int, float]:
    """``1 / sigma_max(Z_dq(j2pi.f))`` a chaque bus [pu base Sb] : l'impedance
    de tout le reseau linearise vue de chaque noeud, en un seul modele d'etat
    (une injection de courant d'essai par bus). Lue sur le DAE lineaire,
    pas sur le modele d'etat : a un noeud ideal ou ne se rencontrent que des
    branches inductives, l'impedance est impropre."""
    from g2egrid.analyse.impedance import _bus_element

    if any(el.x0 is None for el in net.all_elements):
        net.solve_equilibrium()
    elements = {b: _bus_element(net, b) for b in bus}
    sorties = [
        n for b in bus for n in (f"{elements[b].name}.vd_o", f"{elements[b].name}.vq_o")
    ]
    dae = net.dae_lineaire(outputs=sorties, injections=tuple(bus))
    entrees = [f"inj{b}.{c}" for b in bus for c in ("id", "iq")]
    Z = dae.reponse_frequentielle(f, entrees, sorties)[0]
    out = {}
    for k, b in enumerate(bus):
        z = Z[2 * k : 2 * k + 2, 2 * k : 2 * k + 2]
        out[b] = 1.0 / float(np.linalg.svd(z, compute_uv=False)[0])
    return out


def puissance_court_circuit(
    t: Tables,
    *,
    gfm: str = "tension",
    gfl: str = "aucun",
    charges: bool = False,
    c: float = 1.0,
    dynamique: bool = False,
    f: float = 1.0,
    machine: str = "soustransitoire",
    correction_cei: bool = True,
    cos_phi: float = 0.9,
) -> RapportScc:
    """La puissance de court-circuit de chaque noeud alternatif de ``t``.

    ``gfm`` (:data:`MODES_GFM`) et ``gfl`` (:data:`MODES_GFL`) disent ce que
    deviennent les convertisseurs ; ``charges`` ajoute les charges en
    impedance constante ; ``c`` est le facteur de tension de la CEI 60909.
    ``machine`` (:data:`MODES_MACHINE`) choisit la reactance des machines
    synchrones : ``"soustransitoire"`` (defaut, x''d, celle de la CEI 60909
    pour I_k'') ou ``"transitoire"`` (x'd, l'ancien calcul, sans facteur de
    correction). ``correction_cei`` applique avec x''d les facteurs K_G, K_S et
    K_T de la CEI 60909-0 (voir le module), ``cos_phi`` est le facteur de
    puissance nominal des machines qu'ils utilisent.
    ``dynamique`` ajoute la puissance vue par l'impedance dq du reseau
    linearise a ``f`` Hz (repere dq), qui demande un load flow et un
    equilibre."""
    if gfm not in MODES_GFM:
        raise ValueError(f"gfm = {gfm!r} : {', '.join(MODES_GFM)}")
    if gfl not in MODES_GFL:
        raise ValueError(f"gfl = {gfl!r} : {', '.join(MODES_GFL)}")
    if machine not in MODES_MACHINE:
        raise ValueError(f"machine = {machine!r} : {', '.join(MODES_MACHINE)}")
    if not 0.0 < cos_phi <= 1.0:
        raise ValueError(f"cos_phi = {cos_phi} : dans ]0, 1]")
    corriger = correction_cei and machine == "soustransitoire"
    facteurs = (
        {"c": c, "cos_phi": cos_phi, "KG": {}, "KS": {}, "bloc": {}, "bornes": {}}
        if corriger
        else None
    )
    net = assembler(t)
    idx = net._bus_index()
    n = len(idx)
    if n == 0:
        raise ValueError("aucun bus alternatif")
    # Y-bus des branches seules : lignes, transformateurs (rapport compris)
    Y = np.zeros((n, n), dtype=complex)
    y_src, i_src, s_gfl, s_gfm = _sources(t, gfm, gfl, machine, facteurs)
    # K_T des transformateurs de reseau ; le transfo d'un bloc centrale porte le
    # K_S de son bloc (calcule avec la machine)
    k_trafo: dict[str, float] = {}
    k_t: dict[str, float] = {}
    if corriger:
        for tid, i, j, _Sn, _U1, _U2, _R, X in t.T_Transfos:
            nom = f"T{int(tid)}_{i}_{j}"
            if int(tid) in facteurs["bloc"]:
                k_trafo[nom] = facteurs["bloc"][int(tid)]
            else:
                k_trafo[nom] = k_t[nom] = 0.95 * c / (1.0 + 0.6 * float(X))

    def estampe(M, a, b, yb, nn):
        M[a, a] += yb / (nn * nn)
        M[a, b] -= yb / nn
        M[b, a] -= yb / nn
        M[b, b] += yb

    stamps = {}  # transfo -> (a, b, 1/Z nu, n), pour les defauts aux bornes
    for el in net.elements:
        if len(el.ports) == 2 and {"R", "X"} <= set(el.pdict):
            a, b = idx[el.buses[0]], idx[el.buses[1]]
            z_nu = complex(el.pdict["R"], el.pdict["X"])
            yb = 1.0 / (z_nu * k_trafo.get(el.name, 1.0))
            nn = el.pdict.get("n", 1.0)
            estampe(Y, a, b, yb, nn)
            stamps[el.name] = (a, b, 1.0 / z_nu, nn)
    for b, v in y_src.items():
        if b in idx:
            Y[idx[b], idx[b]] += v
    if charges:
        for b, v in _charges(t).items():
            if b in idx:
                Y[idx[b], idx[b]] += v
    # Une ile sans source de tension rend Y singuliere : une admittance de
    # fuite infime la regularise sans rien changer ailleurs (Scc ~ 0 la-bas).
    Y += 1e-9 * np.eye(n)
    Z = np.linalg.inv(Y)
    # Defaut aux bornes d'une machine de bloc centrale : la machine porte K_G et
    # son transfo n'est pas corrige (CEI 60909-0, 6.7.1) ; les autres blocs gardent
    # K_S. Une inversion de plus par bloc, reseau de quelques dizaines de noeuds.
    Z_bornes: dict[int, np.ndarray] = {}
    for b, (tid, k_s, k_g) in (facteurs["bornes"] if corriger else {}).items():
        if b not in idx:
            continue
        nom = next(m for m in stamps if m.startswith(f"T{tid}_"))
        a_t, b_t, y_nu, nn = stamps[nom]
        Ya = Y.copy()
        estampe(Ya, a_t, b_t, y_nu - y_nu / k_s, nn)  # transfo nu
        Ya[idx[b], idx[b]] += y_src[b] * (k_s / k_g - 1.0)  # machine en K_G
        Z_bornes[idx[b]] = np.linalg.inv(Ya)
    i_inj = np.zeros(n)
    for b, v in i_src.items():
        if b in idx:
            i_inj[idx[b]] += v

    Un = {int(r[0]): float(r[2]) for r in t.Y_DER}
    dyn = {}
    if dynamique:
        if not net.run_loadflow().converged:
            raise RuntimeError("load flow non converge : pas de Scc dynamique")
        dyn = _scc_dynamique(net, sorted(idx), f)
    noeuds = []
    for b, k in sorted(idx.items()):
        Zk = Z_bornes.get(k, Z)
        zii = Zk[k, k]
        ik = (c + float(np.abs(Zk[k]) @ i_inj)) / abs(zii)
        Scc = ik * t.Sb
        U = Un.get(b, float("nan"))
        sg = s_gfl.get(b, 0.0)
        noeuds.append(
            SccNoeud(
                bus=b,
                Un=U,
                Scc=Scc,
                Ik=Scc / (math.sqrt(3) * U) if U > 0 else float("nan"),
                XR=zii.imag / zii.real if zii.real > 0 else float("inf"),
                R_th=float(zii.real),
                X_th=float(zii.imag),
                S_gfl=sg,
                S_gfm=s_gfm.get(b, 0.0),
                SCR=Scc / sg if sg > 0 else None,
                Scc_dyn=dyn[b] * t.Sb if b in dyn else None,
            )
        )
    # SCR avec facteurs d'interaction : les GFL voisins chargent le noeud
    voisins = [(idx[b], s) for b, s in s_gfl.items() if b in idx and s > 0]
    for x in noeuds:
        if x.S_gfl > 0:
            k = idx[x.bus]
            charge = sum(
                s * (1.0 if j == k else abs(Z[j, k] / Z[k, k])) for j, s in voisins
            )
            x.SCRIF = x.Scc / charge
    return RapportScc(
        noeuds=noeuds,
        hypotheses=dict(
            gfm=gfm,
            gfl=gfl,
            charges=charges,
            c=c,
            dynamique=dynamique,
            f=f,
            machine=machine,
            correction_cei=corriger,
            cos_phi=cos_phi,
            # les facteurs appliques, par bus (K_G, K_S) ou par transfo (K_T)
            facteurs_cei=(
                {
                    "KG": {str(b): float(v) for b, v in facteurs["KG"].items()},
                    "KS": {str(b): float(v) for b, v in facteurs["KS"].items()},
                    "kT": k_t,
                }
                if corriger
                else {}
            ),
        ),
    )
