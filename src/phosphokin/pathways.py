"""Pathway over-representation for phosphoproteins (gene-set files in GMT format).

Reactome and MSigDB provide GMT downloads, e.g.
https://reactome.org/download/current/ReactomePathways.gmt.zip
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom

from .enrichment import bh_adjust


def read_gmt(path: str | Path) -> dict[str, set[str]]:
    """Read a GMT file into ``{set_name: {genes}}``."""
    sets: dict[str, set[str]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = {g for g in parts[2:] if g}
    return sets


def overrepresentation(query, background, gene_sets: dict[str, set[str]],
                       min_size: int = 5, max_size: int = 500) -> pd.DataFrame:
    """Hypergeometric test of ``query`` genes against each gene set.

    Only genes in ``background`` (e.g. all proteins with a quantified site) are
    counted, which avoids inflating enrichment for proteins that were never
    measurable.
    """
    bg = set(background)
    q = set(query) & bg
    rows = []
    for name, genes in gene_sets.items():
        g = genes & bg
        if not (min_size <= len(g) <= max_size):
            continue
        k = len(q & g)
        p = hypergeom.sf(k - 1, len(bg), len(g), len(q)) if k else 1.0
        rows.append({"gene_set": name, "set_size": len(g), "overlap": k,
                     "expected": len(q) * len(g) / len(bg), "p": p,
                     "genes": ";".join(sorted(q & g))})
    res = pd.DataFrame(rows, columns=["gene_set", "set_size", "overlap", "expected", "p", "genes"])
    res["padj"] = bh_adjust(res["p"]) if len(res) else []
    return res.sort_values("p").reset_index(drop=True)
