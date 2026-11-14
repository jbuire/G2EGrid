"""Le catalogue vu de l'atelier : liste des reseaux, lecture, positions du schema.

Partie de :mod:`g2egrid_api.analyse` (ancien ``analyse.py``, decoupe le 2026-10-28)."""

from __future__ import annotations

import g2egrid.reseau.catalogue
from g2egrid.reseau.catalogue import RESEAUX, fiche
from g2egrid.reseau.tables import Tables, tables_vers_json
from g2egrid.reseau.zones import Instance, Reseau, Zone, aplatir
from ..layout import placer


# ------------------------------------------------------------------ catalogue
def catalogue() -> list[dict]:
    """Les fiches des reseaux de test, sans les assembler."""
    out = []
    # les sauvegardes sont relues a CHAQUE appel : RESEAUX est fige a
    # l'import, et un reseau enregistre depuis l'interface n'apparaissait dans
    # la liste qu'apres un redemarrage du serveur
    for nom in dict.fromkeys((*RESEAUX, *g2egrid.reseau.catalogue.sauvegardes())):
        try:
            f = fiche(nom)
        except Exception:  # noqa: BLE001  -- notebook absent ou illisible
            continue
        out.append(
            {
                "nom": f.nom,
                "titre": f.titre,
                "bus": f.bus,
                "lignes": f.lignes,
                "transfos": f.transfos,
                "sm": f.sm,
                "gfm": f.gfm,
                "gfl": f.gfl,
                "charges": f.charges,
                "f0": f.f0,
                "Sb": f.Sb,
                "f_zoom": f.f_zoom,
            }
        )
    return out


def depuis_catalogue(nom: str) -> dict:
    """Un reseau du catalogue, pret a ouvrir dans le constructeur."""
    return _reponse_reseau(g2egrid.reseau.catalogue.charger_reseau(nom))


def _reponse_reseau(r: Reseau) -> dict:
    """Un reseau sous SES DEUX FORMES : a plat, et piece par piece.

    ``tables`` et ``positions`` sont le reseau APLATI, et ce sont eux que
    lisent l'onglet Tables, le schema et toutes les routes de calcul -- ce
    contrat ne bouge pas. ``transport``, ``zones`` et ``instances`` sont la
    MEME chose non aplatie, et c'est tout l'objet de cette fonction : une
    distribution instanciee trois fois n'a qu'une description, et l'editer une
    fois doit suffire. Rendre seulement le reseau a plat revenait a montrer
    trois copies d'une meme zone -- c'est-a-dire a perdre exactement ce que les
    zones apportent.

    Un reseau ecrit a plat (un notebook) n'est pas un cas particulier : il rend
    ``zones = {}``, ``instances = []`` et un ``transport`` egal a ses tables.
    L'interface lit donc toujours la meme forme, sans avoir a demander
    laquelle elle tient.
    """
    t = aplatir(r)
    pos_transport = _positions_connues(r.positions, r.transport)
    return {
        "nom": r.nom,
        "f_zoom": r.f_zoom,
        "tables": tables_vers_json(t),
        # Les positions du reseau a plat sont DERIVEES, jamais stockees : les
        # ranger dans le fichier en ferait une quatrieme copie, qui serait
        # fausse des qu'une zone bouge. Quand il n'y a aucune instance, le
        # reseau a plat EST le transport, et le meme placement sert aux deux.
        "positions": pos_transport if not r.instances else _positions(t),
        "transport": {
            "tables": tables_vers_json(r.transport),
            "positions": pos_transport,
        },
        "zones": {
            nom: {
                "nom": nom,
                "bus_pcc": z.bus_pcc,
                "tables": tables_vers_json(z.tables),
                "positions": _positions_connues(z.positions, z.tables),
            }
            for nom, z in r.zones.items()
        },
        "instances": [
            {
                "zone": i.zone,
                "prefixe": i.prefixe,
                "bus_amont": i.bus_amont,
                "K": i.K,
                "transfo": None if i.transfo is None else [float(v) for v in i.transfo],
                "charge_residuelle": bool(i.charge_residuelle),
                "shunt_pcc": i.shunt_pcc,
                "bus": _bus_d_instance(r.zones[i.zone], i),
            }
            for i in r.instances
        ],
    }


def _bus_d_instance(z: Zone, inst: Instance) -> dict:
    """Les numeros qu'une instance occupe DANS LE RESEAU A PLAT.

    C'est ce qui relie les deux formes : sans cela l'interface voit 54 bus a
    plat et ne sait pas lequel appartient a quelle instance, donc quelle zone
    ouvrir quand on en clique un. La regle est celle de ``zones.aplatir`` --
    le prefixe sur tout, sauf le PCC qui fusionne avec le bus amont quand il
    n'y a pas de transfo de poste. La redire ici est un risque de derive, et
    c'est pourquoi ``test_api`` ne la croit pas : il verifie que ces numeros
    sont exactement ceux que le ``Y_DER`` aplati porte en plus du transport.
    """
    ac = [
        int(ligne[0]) + inst.prefixe
        for ligne in z.tables.Y_DER
        if not (inst.transfo is None and int(ligne[0]) == z.bus_pcc)
    ]
    return {
        "ac": ac,
        "dc": [int(ligne[0]) + inst.prefixe for ligne in z.tables.T_BusDC],
        "pcc": inst.bus_amont if inst.transfo is None else z.bus_pcc + inst.prefixe,
    }


def _noeuds(t: Tables) -> list[str]:
    """Tous les noeuds d'un jeu de tables, continus compris, dans l'ordre."""
    return [str(int(r[0])) for r in t.Y_DER] + [str(int(r[0])) for r in t.T_BusDC]


def _positions_connues(stockees: dict | None, t: Tables) -> dict[str, list[float]]:
    """Les positions du fichier si elles suffisent, un placement neuf sinon.

    Un fichier peut porter un dessin fait a la main, et le recalculer le
    jetterait ; mais un dessin INCOMPLET -- un bus ajoute depuis la derniere
    sauvegarde -- est pire qu'aucun, parce qu'un noeud sans position atterrit
    au centre du canevas, par-dessus les autres. On ne panache donc pas : ou
    bien le dessin couvre tout et on le garde, ou bien on replace tout.

    Les positions devenues orphelines (un bus supprime) sont ecartees, pour
    qu'un fichier ancien ne fasse pas apparaitre des noeuds qui n'existent
    plus.
    """
    noeuds = _noeuds(t)
    p = {str(k): [float(v[0]), float(v[1])] for k, v in (stockees or {}).items()}
    if noeuds and set(noeuds) <= set(p):
        return {b: p[b] for b in noeuds}
    return _positions(t)


def _positions(t: Tables) -> dict[str, list[float]]:
    """Un placement automatique de TOUS les noeuds, continus compris.

    Les bus continus doivent entrer dans le meme placement que les
    alternatifs, et non etre ajoutes apres : un noeud sans position se
    retrouve au centre du canevas, et sur une liaison a deux terminaux les
    deux s'empilent au meme endroit. Les aretes du placement comprennent
    donc les cables continus ET les raccordements de T_VSCs -- ces derniers
    ne sont pas des branches du modele, mais ce sont eux qui disent quel
    noeud continu va pres de quel poste alternatif.
    """
    bus = [int(r[0]) for r in t.Y_DER] + [int(r[0]) for r in t.T_BusDC]
    aretes = [(int(r[1]), int(r[2])) for r in t.T_Lines_SI]
    aretes += [(int(r[1]), int(r[2])) for r in t.T_Lines_pu]
    aretes += [(int(r[1]), int(r[2])) for r in t.T_Transfos]
    aretes += [(int(r[1]), int(r[2])) for r in t.T_LinesDC]
    aretes += [(int(r[1]), int(r[2])) for r in t.T_VSCs]
    return {str(b): [x, y] for b, (x, y) in placer(bus, aretes).items()}


def _iles(ids: list[int], voisins: dict[int, set[int]]) -> list[set[int]]:
    """Les composantes connexes du graphe ALTERNATIF.

    Le parcours partait autrefois du premier bus et declarait isole tout ce
    qu'il n'atteignait pas. C'etait juste tant qu'un reseau n'avait qu'une
    ile ; des qu'une liaison continue en separe deux, la seconde entiere
    etait denoncee comme isolee.
    """
    vus: set[int] = set()
    out: list[set[int]] = []
    for depart in ids:
        if depart in vus:
            continue
        ile, pile = set(), [depart]
        while pile:
            b = pile.pop()
            if b in ile:
                continue
            ile.add(b)
            pile += [v for v in voisins.get(b, ()) if v not in ile]
        vus |= ile
        out.append(ile)
    return out
