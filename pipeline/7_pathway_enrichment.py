"""Step 7 (optional) - Pathway over-representation of regulated phosphoproteins.

For each comparison, tests whether human genes carrying up- (or down-)
regulated sites are over-represented in gene sets from a GMT file (e.g. Reactome),
using all genes with a quantified, human-mapped site as the background.

Example
  python 7_pathway_enrichment.py --gmt ../data/ReactomePathways.gmt

Outputs
  results/pathway_enrichment.tsv
"""

from __future__ import annotations

import argparse

import pandas as pd

from _common import RESULTS, out
from phosphokin.pathways import overrepresentation, read_gmt


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gmt", required=True, help="gene-set file (GMT, human gene symbols)")
    ap.add_argument("--status", default=str(RESULTS / "site_status.tsv"))
    ap.add_argument("--sites", default=str(RESULTS / "sites_human.tsv"))
    args = ap.parse_args()

    sets = read_gmt(args.gmt)
    genes = pd.read_csv(args.sites, sep="\t").dropna(subset=["human_gene"])[["site_id", "human_gene"]]
    status = pd.read_csv(args.status, sep="\t").merge(genes, on="site_id")

    parts = []
    for comp, st in status.groupby("comparison"):
        background = st["human_gene"].unique()
        for direction in ("up", "down"):
            res = overrepresentation(st.loc[st["status"] == direction, "human_gene"].unique(), background, sets)
            res.insert(0, "direction", direction)
            res.insert(0, "comparison", comp)
            parts.append(res)
    pd.concat(parts).to_csv(out(RESULTS / "pathway_enrichment.tsv"), sep="\t", index=False)
    print(f"tested {len(sets):,} gene sets across {status['comparison'].nunique()} comparisons")


if __name__ == "__main__":
    main()
