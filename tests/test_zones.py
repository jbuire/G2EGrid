"""Les zones interconnectees, les fichiers de reseau et les exports.

LE TEST DECISIF EST ``test_zones_donnent_la_meme_matrice_a``. Le reste mesure
des morceaux ; celui-la mesure la seule chose qui compte : un reseau ecrit en
zones et le meme reseau ecrit a plat doivent donner LE MEME MODELE. Il est
ecrit en trois assertions de cout croissant -- les tables, puis les noms et
l'ordre des elements, puis la matrice ``A`` -- pour que l'echec dise OU ca
diverge et pas seulement QU'il diverge. Comparer ``A`` d'abord couterait le
meme temps et n'apprendrait rien : une permutation de deux lignes de ``Y_DER``
et une erreur d'echelle sur ``K`` s'y ressemblent.
"""

from __future__ import annotations

import json
import re

import numpy as np
import pytest

from g2egrid import catalogue
from g2egrid.catalogue import charger, charger_reseau
from g2egrid.export import vers_matlab
from g2egrid.tables import TABLES, Tables, assembler
from g2egrid.zones import (
    Instance,
    Reseau,
    TransfoPoste,
    Zone,
    aplatir,
    renumeroter,
    reseau_depuis_json,
    reseau_vers_json,
)

MATLAB = (
    catalogue.DOSSIER.parents[1]
    / "matlab"
    / "G2Elin_V2"
    / "reseaux_test"
    / "reseaux_classiques"
)


# ==========================================================================
# de quoi fabriquer de petites zones
# ==========================================================================
def transport_minimal(**modifs) -> Tables:
    """Deux bus 20 kV, un slack, de quoi raccorder une zone."""
    champs = dict(
        Y_DER=[
            (1, 1, 20e3, 1.0, 0, 0, 0, 0, 0),
            (2, 3, 20e3, 1.0, 0, 0, 0, 1e6, 0.5e6),
        ],
        T_Lines_pu=[(1, 1, 2, 1, 0.01, 0.1, 0.02, 20e3, 100e6)],
        f0=50.0,
        Sb=100e6,
        nom="transport",
    )
    return Tables(**(champs | modifs))


def zone_minimale(**modifs) -> Zone:
    """Un PCC et un bus de charge, sans slack."""
    champs = dict(
        Y_DER=[
            (0, 3, 20e3, 1.0, 0, 0, 0, 0, 0),
            (1, 3, 20e3, 1.0, 0, 0, 0, 2e6, 0.4e6),
        ],
        T_Lines_SI=[(1, 0, 1, 2.0, 0.5, 0.7, 4e-5)],
        f0=50.0,
        Sb=100e6,
        nom="zone",
    )
    return Zone(tables=Tables(**(champs | modifs)), bus_pcc=0, nom="zone")


def zone_complete() -> Tables:
    """Une zone qui touche TOUTES les tables, pour eprouver l'echelle K."""
    return Tables(
        Y_DER=[
            (0, 3, 20e3, 1.0, 0, 0, 0, 0, 0),
            (1, 3, 20e3, 1.0, 0, 1e6, 0.2e6, 2e6, 0.4e6),
        ],
        T_Lines_SI=[(1, 0, 1, 2.0, 0.5, 0.7, 4e-5)],
        T_Lines_pu=[(2, 0, 1, 1, 0.01, 0.1, 0.02, 20e3, 50e6)],
        T_Shunts=[(1, 1, 0.3)],
        T_Transfos=[(1, 0, 1, 25e6, 20e3, 20e3, 0.002, 0.1)],
        T_IBs=[(1, 1, 5e6, 1.0, 0, 0.001, 0.01)],
        T_SMs=[(1, 1, 10e6)],
        T_GFM1s=[(1, 1, 3e6, 1e3)],
        T_GFM3s=[(2, 1, 4e6, 40e3)],
        T_GFL2s=[(3, 1, 5e6, 40e3)],
        T_Loads=[(1, 1, 0.7e6, 0.1e6, "P")],
        T_BusDC=[(0, 1, 40e3, 1.0, 1e6)],
        T_LinesDC=[(1, 0, 1, 10.0, 0.05, 0.2, 0.3)],
        T_ShuntsDC=[(1, 0, 0.02)],
        T_VSCs=[(1, 1, 0, "Vdc")],
        f0=50.0,
        Sb=100e6,
        nom="complete",
    )


# ==========================================================================
# LE TEST DECISIF
# ==========================================================================
@pytest.fixture(scope="module")
def paire_wscc3sm():
    """Les tables de ``wscc3sm_cigre3`` et celles de sa reecriture en zones.

    L'une vient du notebook ecrit a plat, l'autre du fichier natif
    ``g2egrid/reseaux/wscc3sm_zones.g2egrid.json`` : la chaine complete est donc
    eprouvee, du fichier aux tables, et pas seulement :func:`aplatir`.
    """
    return charger("wscc3sm_cigre3"), charger("wscc3sm_zones")


def test_zones_tables_identiques_au_reseau_a_plat(paire_wscc3sm):
    """1re assertion : les tables, ligne a ligne.

    La moins chere des trois, et la seule qui designe la ligne fautive. Un
    ecart ici veut dire que la numerotation, l'ordre d'emission ou une regle
    d'echelle K ne reproduit pas ce que le notebook ecrivait a la main.
    """
    plat, zones = paire_wscc3sm
    for nom in TABLES:
        a = [tuple(r) for r in getattr(plat, nom)]
        b = [tuple(r) for r in getattr(zones, nom)]
        assert len(a) == len(b), f"{nom} : {len(a)} lignes a plat, {len(b)} en zones"
        for k, (la, lb) in enumerate(zip(a, b)):
            assert la == lb, f"{nom}, ligne {k}\n  a plat : {la}\n  zones  : {lb}"
    assert plat.PARAMS_PAR_BUS == zones.PARAMS_PAR_BUS
    assert (plat.f0, plat.Sb) == (zones.f0, zones.Sb)


def test_zones_memes_elements_dans_le_meme_ordre(paire_wscc3sm):
    """2e assertion : les noms et l'ordre des elements assembles.

    L'ordre d'insertion fixe l'ordre des etats : deux reseaux aux memes tables
    mais assembles dans un autre ordre ont la meme physique et une matrice
    ``A`` permutee. C'est ce que cette assertion separe de la suivante.
    """
    plat, zones = paire_wscc3sm
    a = [el.name for el in assembler(plat).all_elements]
    b = [el.name for el in assembler(zones).all_elements]
    assert len(a) == len(b), f"{len(a)} elements a plat, {len(b)} en zones"
    for k, (x, y) in enumerate(zip(a, b)):
        assert x == y, f"element de rang {k} : '{x}' a plat, '{y}' en zones"


def test_zones_donnent_la_meme_matrice_a(paire_wscc3sm):
    """3e assertion : la matrice ``A`` a 1e-12. La seule qui prouve.

    Les deux precedentes peuvent passer alors que le modele differe : une
    charge residuelle mal retranchee laisse les tables identiques en forme et
    change le point de fonctionnement, donc la linearisation. Le seuil est
    absolu et volontairement severe -- les deux chemins font les memes
    operations sur les memes flottants, l'ecart attendu est nul.
    """
    plat, zones = paire_wscc3sm
    modeles = []
    for t in (plat, zones):
        net = assembler(t)
        net.run_loadflow()
        net.solve_equilibrium()
        modeles.append(net.build_state_space())
    A, B = modeles[0].A, modeles[1].A
    assert A.shape == B.shape, f"{A.shape[0]} etats a plat, {B.shape[0]} en zones"
    assert modeles[0].states == modeles[1].states
    ecart = np.abs(A - B).max()
    assert ecart < 1e-12, f"ecart maximal sur A : {ecart:.3e}"


# ==========================================================================
# l'echelle K, table par table
# ==========================================================================
def test_echelle_k_table_par_table():
    """K departs identiques : puissances x K, impedances / K, table par table.

    Chaque regle est verifiee sur sa propre colonne, parce qu'elles ne sont
    justement pas toutes les memes -- ``T_Lines_pu`` ne voit bouger que sa
    base, et ``T_VSCs`` ne bouge pas du tout.
    """
    K = 3
    t = zone_complete()
    r = renumeroter(t, 100, K)

    # Y_DER : Pg, Qg, Pl, Ql
    assert r.Y_DER[1][5:9] == (1e6 * K, 0.2e6 * K, 2e6 * K, 0.4e6 * K)
    assert r.Y_DER[1][2] == 20e3  # Un ne bouge pas
    # T_Lines_SI : R/K, X/K, B*K
    assert r.T_Lines_SI[0][4:7] == (0.5 / K, 0.7 / K, 4e-5 * K)
    assert r.T_Lines_SI[0][3] == 2.0  # la longueur ne bouge pas
    # T_Lines_pu : Sb_VA seulement
    assert r.T_Lines_pu[0][4:8] == (0.01, 0.1, 0.02, 20e3)
    assert r.T_Lines_pu[0][8] == 50e6 * K
    # T_Shunts, T_Transfos, DER
    assert r.T_Shunts[0][2] == 0.3 * K
    assert r.T_Transfos[0][3] == 25e6 * K
    assert r.T_Transfos[0][6:8] == (0.002, 0.1)  # les pu du transfo ne bougent pas
    assert r.T_IBs[0][2] == 5e6 * K
    assert r.T_SMs[0][2] == 10e6 * K
    assert r.T_GFM1s[0][2] == 3e6 * K
    assert r.T_GFM3s[0][2] == 4e6 * K
    assert r.T_GFL2s[0][2] == 5e6 * K
    # T_Loads : P, Q ; le type reste un type
    assert r.T_Loads[0][2:4] == (0.7e6 * K, 0.1e6 * K)
    assert r.T_Loads[0][4] == "P"
    # le continu
    assert r.T_BusDC[0][4] == 1e6 * K
    assert r.T_BusDC[0][2] == 40e3
    assert r.T_LinesDC[0][4:7] == (0.05 / K, 0.2 / K, 0.3 * K)
    assert r.T_ShuntsDC[0][2] == 0.02 * K
    # T_VSCs : aucune grandeur physique, seuls les numeros bougent
    assert r.T_VSCs[0] == (101, 101, 100, "Vdc")


def test_renumerotation_decale_bus_et_identifiants():
    """Les identifiants prennent le prefixe, sinon deux instances se marchent
    dessus : le nom d'un element assemble est bati sur son ID."""
    r = renumeroter(zone_complete(), 200, 1)
    assert [int(x[0]) for x in r.Y_DER] == [200, 201]
    assert r.T_Lines_SI[0][:3] == (201, 200, 201)
    assert r.T_Transfos[0][:3] == (201, 200, 201)
    assert [int(x[0]) for x in r.T_BusDC] == [200]
    assert r.T_LinesDC[0][:3] == (201, 200, 201)


def test_echelle_k_egale_un_ne_change_rien():
    """K = 1 doit laisser les valeurs BIT a bit : c'est le cas courant."""
    t = zone_complete()
    r = renumeroter(t, 100, 1)
    for nom in TABLES:
        for avant, apres in zip(getattr(t, nom), getattr(r, nom)):
            assert tuple(avant)[3:] == tuple(apres)[3:]


def test_params_par_bus_decale_les_cles_sans_toucher_aux_valeurs():
    """Un reglage ne depend pas du nombre de departs agreges."""
    t = zone_complete()
    t.PARAMS_PAR_BUS = {1: {"H": 4.5}}
    r = renumeroter(t, 300, 7)
    assert r.PARAMS_PAR_BUS == {301: {"H": 4.5}}


def test_reglages_de_famille_replies_par_bus():
    """Deux zones peuvent regler leurs GFM autrement : un GFM_PARAMS global ne
    saurait pas les distinguer, PARAMS_PAR_BUS si."""
    zt = zone_complete()
    zt.GFM_PARAMS = {"Kp": 0.02}
    zt.GFL_PARAMS = {"Ti": 0.01}
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": Zone(tables=zt, bus_pcc=0, nom="z")},
        instances=[Instance("z", 100, 2, K=1, transfo=None)],
    )
    t = aplatir(r)
    assert t.GFM_PARAMS == {}  # rien ne remonte au global
    assert t.PARAMS_PAR_BUS[101] == {"Kp": 0.02, "Ti": 0.01}


# ==========================================================================
# le raccordement : transfo de poste, fusion, charge residuelle
# ==========================================================================
def test_transfo_de_poste_appartient_a_l_instance():
    """Sa puissance depend de K, son primaire est le bus amont : la zone ne
    peut donc pas le porter."""
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": zone_minimale()},
        instances=[
            Instance(
                "z",
                100,
                2,
                K=5,
                transfo=TransfoPoste(10e6, 20e3, 20e3, 0.004, 0.11),
            )
        ],
    )
    t = aplatir(r)
    assert t.T_Transfos[0] == (100, 2, 100, 10e6 * 5, 20e3, 20e3, 0.004, 0.11)
    assert len(t.Y_DER) == 4  # 2 de transport + 2 de zone, rien de fusionne


def test_sans_transfo_le_pcc_fusionne_avec_le_bus_amont():
    """Charges et shunts s'additionnent, et il ne reste qu'un bus."""
    z = zone_minimale()
    z.tables.Y_DER[0] = (0, 3, 20e3, 1.0, 0, 0, 0, 0.5e6, 0.1e6)
    z.tables.T_Shunts = [(1, 0, 0.05)]
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": z},
        instances=[Instance("z", 100, 2, K=2, transfo=None)],
    )
    t = aplatir(r)
    assert [int(x[0]) for x in t.Y_DER] == [1, 2, 101]  # le PCC a disparu
    bus2 = t.bus(2)
    assert bus2[7] == 1e6 + 0.5e6 * 2  # charge du transport + celle du PCC, x K
    assert bus2[8] == 0.5e6 + 0.1e6 * 2
    assert t.T_Shunts[0] == (101, 2, 0.05 * 2)  # le shunt du PCC va au bus amont
    assert t.T_Lines_SI[0][1] == 2  # la ligne de la zone part du bus amont


def test_charge_residuelle_conserve_la_consommation_du_bus_amont():
    """C'est ce qui permet de reproduire ``wscc3sm_cigre3`` : le bus WSCC garde
    exactement sa consommation, la distribution en reprend une part."""
    z = zone_minimale()
    P, Q = z.charge_totale()
    assert (P, Q) == (2e6, 0.4e6)
    K = 3
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": z},
        instances=[
            Instance(
                "z",
                100,
                2,
                K=K,
                transfo=TransfoPoste(10e6, 20e3, 20e3, 0.004, 0.11),
                charge_residuelle=True,
            )
        ],
    )
    t = aplatir(r)
    assert t.bus(2)[7] == 1e6 - K * P
    assert t.bus(2)[8] == 0.5e6 - K * Q
    assert sum(x[7] for x in t.Y_DER) == 1e6  # la charge totale ne bouge pas


def test_ordre_d_emission_transport_puis_instances():
    """Table par table, et les instances dans l'ordre : c'est ce contrat qui
    fait que ``A`` est comparable terme a terme."""
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": zone_minimale()},
        instances=[
            Instance("z", 100, 2, transfo=TransfoPoste(1e6, 20e3, 20e3, 0, 0.1)),
            Instance("z", 200, 2, transfo=TransfoPoste(1e6, 20e3, 20e3, 0, 0.1)),
        ],
    )
    t = aplatir(r)
    assert [int(x[0]) for x in t.Y_DER] == [1, 2, 100, 101, 200, 201]
    assert [int(x[0]) for x in t.T_Transfos] == [100, 200]


def test_poste_emis_avant_ce_qu_il_alimente():
    """Le transfo et le gradin du poste passent devant les transfos et les
    shunts de la zone, comme on les lit sur un schema."""
    z = zone_minimale()
    z.tables.T_Transfos = [(1, 0, 1, 5e6, 20e3, 20e3, 0.001, 0.05)]
    z.tables.T_Shunts = [(1, 1, 0.1)]
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": z},
        instances=[
            Instance(
                "z",
                100,
                2,
                transfo=TransfoPoste(1e6, 20e3, 20e3, 0, 0.1),
                shunt_pcc=0.2,
            )
        ],
    )
    t = aplatir(r)
    assert [int(x[0]) for x in t.T_Transfos] == [100, 101]
    assert t.T_Shunts == [(100, 100, 0.2), (101, 101, 0.1)]


# ==========================================================================
# les validations qui doivent echouer
# ==========================================================================
def test_une_zone_ne_peut_pas_avoir_de_slack():
    z = zone_minimale()
    z.tables.Y_DER[0] = (0, 1, 20e3, 1.0, 0, 0, 0, 0, 0)
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": z},
        instances=[Instance("z", 100, 2, transfo=None)],
    )
    with pytest.raises(ValueError, match="slack"):
        aplatir(r)


def test_prefixe_sous_cent_refuse():
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": zone_minimale()},
        instances=[Instance("z", 50, 2, transfo=None)],
    )
    with pytest.raises(ValueError, match="prefixe"):
        aplatir(r)


def test_prefixes_identiques_refuses():
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": zone_minimale()},
        instances=[
            Instance("z", 100, 2, transfo=None),
            Instance("z", 100, 2, transfo=None),
        ],
    )
    with pytest.raises(ValueError, match="deja pris"):
        aplatir(r)


def test_prefixes_qui_se_chevauchent_refuses():
    """Deux plages qui se recouvrent donneraient deux bus de meme numero, et
    le second ecraserait le premier sans un mot."""
    grande = zone_minimale()
    grande.tables.Y_DER.append((5, 3, 20e3, 1.0, 0, 0, 0, 1e6, 0.1e6))
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": grande},
        instances=[
            Instance("z", 100, 2, transfo=None),
            Instance("z", 104, 2, transfo=None),
        ],
    )
    with pytest.raises(ValueError, match="chevauchent"):
        aplatir(r)


def test_bus_amont_inconnu_refuse():
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": zone_minimale()},
        instances=[Instance("z", 100, 42, transfo=None)],
    )
    with pytest.raises(ValueError, match="bus amont 42"):
        aplatir(r)


def test_f0_different_refuse():
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": zone_minimale(f0=60.0)},
        instances=[Instance("z", 100, 2, transfo=None)],
    )
    with pytest.raises(ValueError, match="f0"):
        aplatir(r)


def test_sb_different_refuse():
    """Les colonnes en pu base systeme changeraient de sens a l'aplatissement."""
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": zone_minimale(Sb=25e6)},
        instances=[Instance("z", 100, 2, transfo=None)],
    )
    with pytest.raises(ValueError, match="Sb"):
        aplatir(r)


def test_fusion_exige_des_tensions_egales():
    z = zone_minimale()
    z.tables.Y_DER[0] = (0, 3, 110e3, 1.0, 0, 0, 0, 0, 0)
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": z},
        instances=[Instance("z", 100, 2, transfo=None)],
    )
    with pytest.raises(ValueError, match="Un doivent"):
        aplatir(r)


def test_reglage_de_famille_du_transport_qui_deborde_refuse():
    """Sans cette garde, un ``GFM_PARAMS`` de transport s'appliquerait aussi
    aux convertisseurs de la zone, sans que rien ne le dise."""
    zt = zone_complete()
    r = Reseau(
        transport=transport_minimal(GFM_PARAMS={"Kp": 0.05}),
        zones={"z": Zone(tables=zt, bus_pcc=0, nom="z")},
        instances=[Instance("z", 100, 2, transfo=None)],
    )
    with pytest.raises(ValueError, match="GFM_PARAMS"):
        aplatir(r)


def test_zone_inconnue_refusee():
    r = Reseau(
        transport=transport_minimal(),
        zones={"z": zone_minimale()},
        instances=[Instance("absente", 100, 2, transfo=None)],
    )
    with pytest.raises(ValueError, match="zone inconnue"):
        aplatir(r)


# ==========================================================================
# le fichier .g2egrid.json
# ==========================================================================
@pytest.mark.parametrize("nom", ["smib", "wscc9", "wscc3sm_zones"])
def test_aller_retour_json(nom):
    """Un reseau ecrit puis relu doit rendre les memes tables, a l'identique.

    Les trois reseaux ne mesurent pas la meme chose : ``smib`` et ``wscc9``
    sont des reseaux a plat convertis (le fichier n'a qu'une piece, le
    transport), ``wscc3sm_zones`` a une zone instanciee trois fois, avec
    transfos de poste, gradins et charges residuelles.
    """
    r = charger_reseau(nom)
    d = reseau_vers_json(r)
    # le format doit survivre a un passage par du JSON reel, pas seulement par
    # un dict python : c'est la que les cles entieres deviennent des chaines.
    bis = reseau_depuis_json(json.loads(json.dumps(d)))
    a, b = aplatir(r), aplatir(bis)
    for table in TABLES:
        assert [tuple(x) for x in getattr(a, table)] == [
            tuple(x) for x in getattr(b, table)
        ], table
    assert (a.f0, a.Sb, a.nom) == (b.f0, b.Sb, b.nom)
    assert a.PARAMS_PAR_BUS == b.PARAMS_PAR_BUS
    assert len(bis.zones) == len(r.zones)
    assert len(bis.instances) == len(r.instances)


def test_fichier_natif_decrit_bien_des_zones():
    """``wscc3sm_zones`` doit vraiment etre ecrit en zones, et non a plat."""
    r = charger_reseau("wscc3sm_zones")
    assert len(r.zones) == 1
    assert len(r.instances) == 3
    assert [i.K for i in r.instances] == [4, 3, 3]
    assert [i.bus_amont for i in r.instances] == [5, 6, 8]
    assert all(i.charge_residuelle for i in r.instances)
    assert len(r.transport.Y_DER) == 9  # le transport seul
    assert len(aplatir(r).Y_DER) == 54  # et 54 bus une fois a plat


def test_le_fichier_natif_passe_avant_le_notebook():
    """La resolution est un contrat : ``g2egrid/reseaux/`` d'abord."""
    assert catalogue.chemin_natif("wscc9") is not None
    assert catalogue.chemin_natif("wscc3sm_cigre3") is None
    assert "wscc3sm_zones" in catalogue.RESEAUX


def test_charger_rend_toujours_des_tables_aplaties():
    """Le contrat de ``charger`` ne bouge pas, quelle que soit la source."""
    for nom in ("wscc9", "wscc3sm_zones", "wscc3sm_cigre3"):
        t = charger(nom)
        assert isinstance(t, Tables)
        assert t.Y_DER


def test_fiche_marche_sur_les_deux_sources():
    f = catalogue.fiche("wscc3sm_zones")
    assert f.bus == 54
    assert f.sm == 3
    assert f.f0 == 60.0
    assert f.f_zoom == 5.0
    assert f.titre == "WSCC3SM_ZONES"


def test_charger_reseau_sur_un_notebook_rend_une_seule_piece():
    """Un notebook est deja a plat : il devient un reseau sans zone, pour que
    le code qui manipule des reseaux structures marche sur les deux sources."""
    r = charger_reseau("wscc3sm_cigre3")
    assert r.zones == {}
    assert r.instances == []
    assert len(r.transport.Y_DER) == 54


def test_format_inconnu_refuse():
    with pytest.raises(ValueError, match="format"):
        reseau_depuis_json({"format": "autre", "version": 1})


def test_version_future_refusee():
    with pytest.raises(ValueError, match="version"):
        reseau_depuis_json({"format": "g2egrid", "version": 99})


# ==========================================================================
# l'export MATLAB
# ==========================================================================
def _matrices_m(source: str) -> dict[str, list[list[float]]]:
    """Les matrices d'un fichier reseau MATLAB, blancs et commentaires enleves.

    Comparer les deux fichiers caractere par caractere n'aurait aucun sens :
    ``WSCC.m`` ecrit ``16.5e3`` la ou l'export ecrit ``16500``, et les deux
    disent le meme nombre. Ce qui doit coincider, ce sont les VALEURS.
    """
    lignes = []
    for ligne in source.splitlines():
        i = ligne.find("%")
        lignes.append(ligne if i < 0 else ligne[:i])
    out: dict[str, list[list[float]]] = {}
    for m in re.finditer(r"(\w+)\s*=\s*\[(.*?)\]\s*;", "\n".join(lignes), re.S):
        corps = [r.split() for r in m.group(2).split(";")]
        out[m.group(1)] = [[float(c) for c in r] for r in corps if r]
    return out


def test_export_matlab_de_wscc9_equivaut_a_wscc_m():
    """Le ``.m`` produit doit porter les memes tables que ``WSCC.m``.

    ``WSCC.m`` est le modele du format : c'est lui qui dit quelles matrices
    ``buildFromMatrices`` attend, et dans quel etat. Le seul point non trivial
    est ``T_Loads``, qui cote MATLAB porte les charges dynamiques, c'est-a-dire
    les colonnes ``Pl_W``/``Ql_var`` de ``Y_DER``.
    """
    reference = _matrices_m((MATLAB / "WSCC.m").read_text(encoding="utf-8"))
    source, avertissements = vers_matlab(charger("wscc9"))
    produit = _matrices_m(source)
    assert avertissements == []
    for nom, attendu in reference.items():
        assert nom in produit, f"{nom} absente de l'export"
        assert produit[nom] == attendu, nom
    assert produit["T_Loads"] == [
        [1, 5, 125e6, 50e6],
        [2, 6, 90e6, 30e6],
        [3, 8, 100e6, 35e6],
    ]


def test_export_matlab_refuse_le_continu():
    """Un reseau a liaison continue ampute de sa liaison n'est pas une
    approximation du reseau : c'est un autre reseau, qui converge et qui ment."""
    with pytest.raises(ValueError, match="continue"):
        vers_matlab(charger("hvdc2t"))


def test_export_matlab_avertit_de_ce_qu_il_laisse():
    """Shunts et charges speciales partent en commentaire, avec un mot."""
    t = charger("wscc9")
    t.T_Shunts = [(1, 5, 0.2)]
    t.T_Loads = [(1, 7, 0.5e6, 0.1e6, "P")]
    t.GFM_PARAMS = {"Kp": 0.02}
    source, av = vers_matlab(t)
    assert len(av) == 3
    assert any("T_Shunts" in a for a in av)
    assert any("T_Loads" in a for a in av)
    assert any("GFM_PARAMS" in a for a in av)
    # ce qui est omis reste LISIBLE dans le fichier, en commentaire
    assert "addShunt" in source
    assert _matrices_m(source)["T_Loads"] == [
        [1, 5, 125e6, 50e6],
        [2, 6, 90e6, 30e6],
        [3, 8, 100e6, 35e6],
    ]


def test_export_matlab_nomme_sa_fonction_comme_le_fichier():
    source, _ = vers_matlab(charger("wscc3sm_zones"))
    assert source.startswith("function net = WSCC3SM_ZONES()")
    assert source.rstrip().endswith("end")


# ==========================================================================
# le notebook complet
# ==========================================================================
def test_notebook_exporte_est_relu_par_le_catalogue(tmp_path, monkeypatch):
    """Un reseau exporte en notebook doit revenir par le catalogue.

    C'est la boucle complete de l'interface : des tables, un notebook, et le
    meme reseau de l'autre cote. Les valeurs passent par la cellule de
    donnees, qui les ecrit a six chiffres significatifs -- ``wscc9`` n'en
    demande pas davantage, et l'egalite est donc exacte.
    """
    from g2egrid import notebook as nbmod

    t = charger("wscc9")
    nbmod.ecrire(
        tmp_path / "essai.ipynb",
        "# Essai",
        nbmod.cellule_donnees(t),
        executer=False,
    )
    monkeypatch.setattr(catalogue, "DOSSIER", tmp_path)
    monkeypatch.setattr(catalogue, "DOSSIER_RESEAUX", tmp_path / "vide")
    relu = charger("essai")
    for table in TABLES:
        assert [tuple(x) for x in getattr(relu, table)] == [
            tuple(x) for x in getattr(t, table)
        ], table
    assert (relu.f0, relu.Sb) == (t.f0, t.Sb)


def test_cellule_d_assemblage_du_notebook_appelle_l_assembleur():
    """La cellule ASSEMBLAGE ne recopie plus ``tables.assembler``, elle
    l'appelle : c'etait la troisieme copie du meme code."""
    from g2egrid import notebook as nbmod

    assert "assembler(t)" in nbmod.ASSEMBLAGE
    assert "net.add(" not in nbmod.ASSEMBLAGE
    assert "from g2egrid.tables import" in nbmod.IMPORTS


def test_cellules_du_notebook_dans_l_ordre():
    from g2egrid import notebook as nbmod

    cs = nbmod.cellules("# intro", "Y_DER = []")
    assert [genre for genre, _ in cs] == [
        "markdown",
        "code",
        "markdown",
        "code",
        "code",
        "markdown",
        "code",
        "markdown",
        "code",
        "markdown",
        "code",
        "markdown",
        "code",
        "code",
    ]
    assert cs[0][1] == "# intro"
    assert cs[3][1] == "Y_DER = []"


# ----------------------------------------------------------------------
# definir une zone a partir d'une selection de bus
# ----------------------------------------------------------------------
def _spectre_lent(t):
    net = assembler(t)
    net.solve_equilibrium()
    lam = net.build_state_space().eig()
    return np.sort_complex(lam[np.abs(lam) > 1e-6])


def test_extraire_une_zone_ne_change_pas_le_spectre():
    from g2egrid.zones import extraire_zone

    r = charger_reseau("ieee14")
    r2 = extraire_zone(r, [8], "z")
    assert list(r2.zones) == ["z"]
    (inst,) = r2.instances
    assert inst.bus_amont == 7 and inst.prefixe == 100
    a, b = _spectre_lent(aplatir(r)), _spectre_lent(aplatir(r2))
    assert len(a) == len(b)
    assert np.max(np.abs(a - b) / np.abs(a)) < 1e-8


def test_extraire_une_zone_refuse_le_slack_et_plusieurs_branches():
    from g2egrid.zones import extraire_zone

    r = charger_reseau("ieee14")
    with pytest.raises(ValueError, match="slack"):
        extraire_zone(r, [1], "z")
    with pytest.raises(ValueError, match="exactement une"):
        extraire_zone(r, [4, 5], "z")


# ----------------------------------------------------------------------
# importer un reseau dans un autre
# ----------------------------------------------------------------------
def test_importer_renumerote_et_passe_le_slack_en_pv():
    from g2egrid.zones import _production_des_slacks, importer_reseau

    r, src = charger_reseau("wscc9"), charger_reseau("ieee14")
    raccord = {"mode": "transfo", "bus_importe": 1, "bus_existant": 5}
    n, info = importer_reseau(r, src, raccord=raccord)
    assert info["prefixe"] == 100
    assert info["bus"][1] == 101 and info["bus"][14] == 114
    t = aplatir(n)
    y = {int(x[0]): x for x in t.Y_DER}
    assert [b for b, x in y.items() if int(x[1]) == 1] == [1]  # une seule reference
    pg = _production_des_slacks(src.transport)[1]
    assert int(y[101][1]) == 2 and y[101][5] == pytest.approx(pg)
    assert (t.T_Transfos[-1][1], t.T_Transfos[-1][2]) == (5, 101)
    net = assembler(t)
    assert net.run_loadflow().converged


def test_importer_retire_le_bus_infini_et_convertit_les_shunts():
    from g2egrid.zones import importer_reseau

    r, src = charger_reseau("wscc9"), charger_reseau("ieee33")
    n, info = importer_reseau(
        r, src, raccord={"mode": "ligne", "bus_importe": 1, "bus_existant": 8}
    )
    t = aplatir(n)
    assert not t.T_IBs
    assert int(next(x for x in t.Y_DER if int(x[0]) == 101)[1]) == 3
    assert any("Sb" in m for m in info["messages"])
    k = src.transport.Sb / r.transport.Sb
    avant = {int(x[0]): x[2] for x in src.transport.T_Shunts}
    for x in t.T_Shunts:
        if int(x[0]) > 100:
            assert x[2] == pytest.approx(avant[int(x[0]) - 100] * k)


def test_importer_par_fusion_et_ses_refus():
    from g2egrid.zones import importer_reseau

    r = charger_reseau("wscc9")
    with pytest.raises(ValueError, match="deux unites"):
        importer_reseau(
            r,
            charger_reseau("wscc9"),
            raccord={"mode": "fusion", "bus_importe": 1, "bus_existant": 1},
        )
    with pytest.raises(ValueError, match="transformateur"):
        importer_reseau(
            r,
            charger_reseau("ieee14"),
            raccord={"mode": "fusion", "bus_importe": 1, "bus_existant": 5},
        )
    # meme tension : le bus 4 d'un wscc9 importe devient le bus 4 existant
    n, info = importer_reseau(
        r,
        charger_reseau("wscc9"),
        raccord={"mode": "fusion", "bus_importe": 4, "bus_existant": 4},
    )
    assert info["bus"][4] == 4 and info["bus"][5] == 105
    t = aplatir(n)
    assert [int(x[0]) for x in t.Y_DER].count(4) == 1
    assert assembler(t).run_loadflow().converged
