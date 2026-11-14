"""FR-005 : ce que coutent les resistances de causalite au point de
fonctionnement, reseau par reseau sur tout le catalogue.

LE POURQUOI PHYSIQUE -- un bus sans capacite shunt (aucune ligne n'y apporte
de charge capacitive) n'a pas d'etat de tension. Pour donner au bloc une forme
causale, ``Node.m`` lui pose une resistance vers la terre
(:class:`~g2egrid.elements.bus.ResistiveBus`, ``NOEUDS = "resistance"``) : la
tension vaut ``R.isum``, R EN BASE SYSTEME. Ce shunt n'existe
pas dans le reseau reel, mais il est verse au Y-bus du load flow (sans quoi
l'equilibre ne serait pas un point fixe du modele) : il consomme ``|V|^2/R``
par noeud, quelle que soit la puissance qui transite au noeud.

``Node.m`` pose ``setResistor(100)`` : 1 % de ``Sb`` par noeud, 1 MW a
100 MVA. Sur un depart HTA dont les charges font 0,1 MW par noeud, il pesait
autant que la charge : 77 % de la charge sur ``ieee33``, ou le slack
fournissait 0,31 pu de trop. Une resistance de 1e5 pu ramenait tout le
catalogue sous 0,1 %, mais rendait le modele mille fois plus raide. Depuis
FR-005, le defaut est le NOEUD IDEAL (``NOEUDS = "ideal"``) : la loi des
noeuds exacte, sans shunt, donc sans pertes fictives ; en temporel, BDF2.

LE SEUIL -- 0,5 % de la charge totale. C'est un ordre de grandeur sous les
pertes reelles d'un reseau (5,5 % sur le depart de Baran et Wu d'origine :
202,7 kW pour 3,7 MW) et sous l'incertitude des donnees de charge : en
dessous, le shunt fictif ne change pas les conclusions d'une etude de point de
fonctionnement. Un reseau sans charge (``smib``, ou la machine debite dans le
bus infini) se rapporte a sa production declaree. Un reseau ajoute au
catalogue au-dela du seuil echoue franchement ; ``AU_DELA`` permettrait d'en
porter un, mesure a l'appui, en ``xfail(strict=True)``.

LES RESEAUX DE PARITE MATLAB sont exclus du seuil : leurs chiffres sont
compares a V2 (ou leur equivalence avec V2 est declaree), qui pose la
resistance de 100 pu de ``Node.m``. Ils la gardent par un reglage explicite
(``NOEUDS = "resistance"`` : reglages du fichier natif, cellule de donnees
du notebook), et le second test verifie qu'ils l'ont.
"""

import pytest

from g2egrid.catalogue import RESEAUX, charger
from g2egrid.elements.bus import NOEUDS_DEFAUT, R_NODE_M, IdealBus, ResistiveBus
from g2egrid.tables import assembler

#: Seuil : pertes fictives / charge totale.
SEUIL = 0.005


#: Reseaux dont un chiffre est compare a MATLAB (V1 ou V2) : ils gardent la
#: resistance de causalite de ``Node.m``, et donc ses pertes fictives.
PARITE_MATLAB = {
    # equivalence DECLAREE (notebook, et 28 etats des deux cotes selon le
    # README des reseaux classiques V2), sans comparaison chiffree publiee ;
    # 1 noeud resistif, 1,2 % de la production (Sb = 1000 MVA)
    "smib": "SMIB.m (V2), equivalent declare du notebook",
    # tab-v2-multimachine : ecart V2 / Python de 8,2e-8 sur les modes lents
    "wscc9": "WSCC.m, ecart V2 / Python publie (docs/manuel/11_validation.md)",
    # tab-v2-multimachine : ecart de 9,7e-5
    "wscc9_gfm": "WSCC_GFM.m, ecart V2 / Python publie",
    # tab-v2-multimachine (1,1e-4) et chapitre transport / distribution, dont
    # les tableaux sont identiques dans V2 et Python a la troisieme decimale
    "wscc9_cigre": "WSCC_CIGRE.m, ecart V2 / Python et chapitre T/D publies",
    # Kundur2Area.m (V2, 114 etats des deux cotes) ; sous le seuil de toute
    # facon (0,15 %)
    "kundur2area": "Kundur2Area.m (V2), et ecarts publies contre ANDES",
}

#: Reseaux au-dela du seuil aujourd'hui : pertes fictives (pu systeme) et part
#: de la charge, mesurees avec le traitement par defaut.
AU_DELA: dict[str, str] = {
    # vide depuis FR-005 (a) : a R = 100 pu, 15 reseaux y figuraient, de 1 %
    # (petits reseaux a 100 MVA) a 77 % (ieee33) de la charge
}


def _pertes_et_charge(nom: str) -> tuple[float, float, int]:
    """Pertes des shunts fictifs, charge de reference (pu systeme) et nombre
    de noeuds resistifs, au point de fonctionnement du traitement par defaut.

    La charge comprend celles des bus (``Y_DER`` et ``T_Loads``, deja fondues
    dans ``bus_specs``), les productions declarees negatives (un groupe qui
    absorbe est une charge) et les charges continues de ``T_LoadsDC``.
    """
    t = charger(nom)
    net = assembler(t)
    lf = net.run_loadflow()
    assert lf.converged, f"{nom} : load flow non converge"
    ix = net._bus_index()
    noeuds = [el for el in net.all_elements if isinstance(el, ResistiveBus)]
    pertes = sum(abs(lf.V[ix[el.bus]]) ** 2 / el.pdict["R"] for el in noeuds)
    specs = net.bus_specs.values()
    charge = (
        sum(s.Pl for s in specs)
        + sum(max(-s.Pg, 0.0) for s in specs)
        + sum(float(r[2]) for r in t.T_LoadsDC) / t.Sb
    )
    if charge <= 0.0:
        charge = sum(max(s.Pg, 0.0) for s in specs)
    return pertes, charge, len(noeuds)


def _cas():
    for nom in RESEAUX:
        if nom in PARITE_MATLAB:
            continue
        marques = []
        if nom in AU_DELA:
            # raises : seul le depassement du seuil est attendu ; un load flow
            # qui casserait ne doit pas passer pour le constat
            marques.append(
                pytest.mark.xfail(
                    strict=True,
                    raises=AssertionError,
                    reason=f"FR-005 : {AU_DELA[nom]}",
                )
            )
        yield pytest.param(nom, marks=marques, id=nom)


@pytest.mark.parametrize("nom", list(_cas()))
def test_pertes_fictives_sous_le_seuil(nom):
    """Les shunts de causalite consomment moins de 0,5 % de la charge."""
    pertes, charge, n = _pertes_et_charge(nom)
    assert charge > 0.0, f"{nom} : ni charge ni production declaree"
    assert pertes <= SEUIL * charge, (
        f"{nom} : {n} noeud(s) resistif(s), pertes fictives {pertes:.4g} pu "
        f"= {100 * pertes / charge:.2f} % de la charge ({charge:.4g} pu), "
        f"seuil {100 * SEUIL:g} %"
    )


@pytest.mark.parametrize("nom", sorted(PARITE_MATLAB))
def test_parite_matlab_garde_la_resistance_de_node_m(nom):
    """Exception explicite : un reseau de parite MATLAB garde, a chaque noeud
    sans capacite, la resistance de causalite de ``Node.m``. Un traitement des
    noeuds qui la changerait deplacerait des comparaisons publiees : ce test
    le fait voir au lieu de le laisser passer."""
    t = charger(nom)
    net = assembler(t)
    sans_capacite = {b for b, s in net.bus_specs.items() if s.C == 0.0}
    resistifs = {
        el.bus: el.pdict["R"] for el in net.all_elements if isinstance(el, ResistiveBus)
    }
    assert set(resistifs) == sans_capacite
    assert all(r == R_NODE_M for r in resistifs.values())


def test_les_listes_ne_citent_que_des_reseaux_du_catalogue():
    """Une liste qui nommerait un reseau disparu ou renomme ne protegerait
    plus rien ; un reseau ne peut pas etre a la fois exempte et en echec."""
    assert set(PARITE_MATLAB) <= set(RESEAUX)
    assert set(AU_DELA) <= set(RESEAUX)
    assert not set(PARITE_MATLAB) & set(AU_DELA)


@pytest.mark.parametrize("nom", ["ieee33", "wscc3sm_cigre3", "wscc9_cigre6"])
def test_hors_parite_le_noeud_ideal_par_defaut(nom):
    """Hors parite MATLAB, chaque noeud sans capacite est un noeud ideal :
    aucun shunt fictif."""
    t = charger(nom)
    assert t.NOEUDS == NOEUDS_DEFAUT == "ideal"
    net = assembler(t)
    assert not [el for el in net.all_elements if isinstance(el, ResistiveBus)]
    assert [el for el in net.all_elements if isinstance(el, IdealBus)]


def test_le_reglage_traverse_les_deux_formats():
    """``NOEUDS`` et ``R_NOEUDS`` voyagent avec les tables (atelier) et dans
    le fichier natif ; le defaut n'est pas ecrit."""
    from g2egrid.catalogue import charger_reseau
    from g2egrid.tables import tables_depuis_json, tables_vers_json
    from g2egrid.zones import reseau_depuis_json, reseau_vers_json

    t = charger("wscc9")
    assert (t.NOEUDS, t.R_NOEUDS) == ("resistance", R_NODE_M)
    d = tables_vers_json(t)
    assert d["NOEUDS"] == "resistance" and "R_NOEUDS" not in d
    t2 = tables_depuis_json(dict(d, R_NOEUDS=1e3))
    assert (t2.NOEUDS, t2.R_NOEUDS) == ("resistance", 1e3)
    assert tables_vers_json(t2)["R_NOEUDS"] == 1e3
    assert "NOEUDS" not in tables_vers_json(charger("ieee14"))

    r = charger_reseau("wscc9")
    j = reseau_vers_json(r)
    assert j["reglages"]["NOEUDS"] == "resistance"
    assert reseau_depuis_json(j).transport.NOEUDS == "resistance"
    r2 = charger_reseau("wscc9_faible")
    assert "NOEUDS" not in reseau_vers_json(r2)["reglages"]
    assert reseau_depuis_json(reseau_vers_json(r2)).transport.NOEUDS == "ideal"


def test_un_bus_garde_sa_propre_resistance():
    """``Bus(R=...)`` l'emporte sur le reglage du reseau ; ``r_causalite``
    doit etre positive."""
    from g2egrid.elements import Bus
    from g2egrid.network import Network

    net = Network(f0=50.0, Sb=1e6, noeuds_sans_capacite="resistance", r_causalite=1e3)
    net.add(Bus(1, 20e3, btype=1))
    net.add(Bus(2, 20e3, R=50.0))
    assert net.buses[1].pdict["R"] == 1e3 and net.buses[2].pdict["R"] == 50.0
    assert net.bus_specs[1].R == 1e3
    with pytest.raises(ValueError, match="r_causalite"):
        Network(r_causalite=0.0)
