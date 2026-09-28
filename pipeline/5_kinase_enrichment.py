"""Step 5 - Kinase motif enrichment for every comparison.

For each comparison, tests whether each kinase is favoured more often among
up- (or down-) regulated sites than among unregulated sites (one-sided Fisher's
exact test, Benjamini-Hochberg correction), following Johnson et al. 2023.

Outputs
  results/kinase_enrichment.tsv   one row per comparison x kinase
"""

from __future__ import annotations

import argparse

import pandas as pd

from _common import RESULTS, out
from phosphokin.enrichment import kinase_enrichment


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--status", default=str(RESULTS / "site_status.tsv"))
    ap.add_argument("--favored", default=str(RESULTS / "favored_kinases.tsv.gz"))
    ap.add_argument("--padj", type=float, default=0.1)
    args = ap.parse_args()

    status = pd.read_csv(args.status, sep="\t").dropna(subset=["status"])
    fav = pd.read_csv(args.favored, sep="\t")

    parts = []
    for comp, st in status.groupby("comparison"):
        res = kinase_enrichment(st, fav)
        res.insert(0, "comparison", comp)
        parts.append(res)
        top = res[res["padj"] <= args.padj].head(3)
        hits = ", ".join(f"{k} ({d})" for k, d in zip(top["kinase"], top["direction"])) or "-"
        print(f"{comp:24s} up={res['n_up'].iat[0]:5d} down={res['n_down'].iat[0]:5d}  top: {hits}")
    pd.concat(parts).to_csv(out(RESULTS / "kinase_enrichment.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
