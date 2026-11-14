"""Le rapport d'etude de l'atelier : autonome, date de son code, et
relisible comme un reseau."""

from __future__ import annotations

import json
import re

import pytest

pytest.importorskip("fastapi")

from g2egrid.catalogue import charger  # noqa: E402
from g2egrid_api import analyse  # noqa: E402
from g2egrid_api.rapport import MARQUE_TABLES, empreinte, rapport_etude  # noqa: E402


def test_rapport_autonome_et_relisible():
    t = charger("wscc9")
    svg = '<svg viewBox="0 0 10 10"><circle r="3" onclick="x()"/><script>y()</script></svg>'
    page = rapport_etude(
        t, analyse.construire, {"Schema": svg}, positions={"1": [0.1, 0.2]}
    )
    assert empreinte(t) in page and "stable" in page and "Modes critiques" in page
    # aucune ressource exterieure, aucun script actif venu des figures
    assert not re.search(r"(src|href)=\"https?://", page)
    assert "onclick" not in page and "y()" not in page
    m = re.search(
        rf'<script type="application/json" id="{MARQUE_TABLES}">(.*?)</script>',
        page,
        re.S,
    )
    d = json.loads(m.group(1))
    assert d["empreinte"] == empreinte(t) and d["positions"] == {"1": [0.1, 0.2]}
    t2 = analyse.tables_depuis_json(d["tables"])
    assert empreinte(t2) == empreinte(t)


def test_empreinte_change_avec_le_reseau():
    t = charger("wscc9")
    t2 = t.copie()
    t2.Y_DER = [(*r[:7], r[7] * 1.01, r[8]) for r in t.Y_DER]
    assert empreinte(t) != empreinte(t2)
