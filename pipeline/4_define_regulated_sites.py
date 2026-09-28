"""Step 4 - Label each site as up / down / unregulated for every comparison.

Two modes:

* ``--limma FILE [FILE ...]``  differential results (columns: site_id, logFC, P.Value);
  one comparison per file (named after the file).
* ``--tissue-contrast``  (demo) for each tissue, compare a site's intensity with
  its median across the other tissues in which it was detected. Intensities are
  median-centred per tissue first. Sites need detection in the tissue and in at
  least ``--min-other`` other tissues.

Outputs
  results/site_status.tsv   site_id, comparison, log2fc, status
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import RESULTS, out
from phosphokin.enrichment import classify_sites


def tissue_contrasts(intensities: pd.DataFrame, min_other: int = 3) -> pd.DataFrame:
    x = intensities.set_index("site_id")
    x = (x - x.median()) * np.log2(10)             # log10 -> log2, median-centred per tissue
    rows = []
    for tissue in x.columns:
        others = x.drop(columns=tissue)
        n_other = others.notna().sum(axis=1)
        keep = x[tissue].notna() & (n_other >= min_other)
        fc = x.loc[keep, tissue] - others.loc[keep].median(axis=1)
        rows.append(pd.DataFrame({"site_id": fc.index, "comparison": tissue, "log2fc": fc.to_numpy()}))
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limma", nargs="*", help="limma-style tables (site_id, logFC, P.Value)")
    ap.add_argument("--tissue-contrast", action="store_true")
    ap.add_argument("--intensities", default=str(RESULTS / "intensities.tsv"))
    ap.add_argument("--min-other", type=int, default=3)
    ap.add_argument("--fc", type=float, default=2.0, help="|log2FC| threshold for up/down")
    ap.add_argument("--p", type=float, default=0.05, help="p-value threshold (limma mode)")
    args = ap.parse_args()

    if args.tissue_contrast:
        res = tissue_contrasts(pd.read_csv(args.intensities, sep="\t"), args.min_other)
        res["status"] = classify_sites(res["log2fc"], fc_threshold=args.fc)
    elif args.limma:
        parts = []
        for f in args.limma:
            t = pd.read_csv(f, sep="\t" if f.endswith((".tsv", ".txt")) else ",")
            parts.append(pd.DataFrame({"site_id": t["site_id"], "comparison": Path(f).stem,
                                       "log2fc": t["logFC"],
                                       "status": classify_sites(t["logFC"], t["P.Value"], args.fc, args.p)}))
        res = pd.concat(parts, ignore_index=True)
    else:
        ap.error("choose --tissue-contrast or --limma FILES")

    res.to_csv(out(RESULTS / "site_status.tsv"), sep="\t", index=False)
    summary = res.groupby(["comparison", "status"]).size().unstack(fill_value=0)
    print(summary.to_string())


if __name__ == "__main__":
    main()
