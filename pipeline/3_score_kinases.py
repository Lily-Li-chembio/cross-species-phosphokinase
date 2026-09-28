"""Step 3 - Score every human-mapped S/T site against 303 Ser/Thr kinases.

Uses the kinome-atlas matrices and reference-phosphoproteome percentiles of
Johnson et al. (Nature 2023). For each site the kinases are ranked by percentile;
the top 15 are the site's "biochemically favoured" kinases, as in the paper.

Outputs
  results/favored_kinases.tsv.gz   site x top-15 kinases (percentile, rank)
  results/reference_check.tsv      how well our percentiles reproduce the published ones
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from _common import KINOME_MATRICES, KINOME_REFERENCE, RESULTS, out
from phosphokin.kinase_atlas import KinomeAtlas


def reference_check(atlas: KinomeAtlas, ref: pd.DataFrame, frac: float = 0.1, seed: int = 0) -> pd.DataFrame:
    """Hold out reference sites, rebuild the percentile curves, compare to published values."""
    rng = np.random.default_rng(seed)
    test = rng.random(len(ref)) < frac
    atlas_cv = KinomeAtlas(kinases=atlas.kinases, log_matrix=atlas.log_matrix)
    atlas_cv.set_reference(ref[~test])
    held = ref[test]
    pred = atlas_cv.percentile(held["SITE_+/-7_AA"])
    pub = held[[f"{k}_percentile" for k in atlas.kinases]].to_numpy(float)
    err = np.abs(pred - pub)
    top = lambda a: np.argsort(-a, axis=1)[:, :15]
    overlap = np.array([len(set(x) & set(y)) for x, y in zip(top(pred), top(pub))])
    return pd.DataFrame({"metric": ["held-out sites", "mean |percentile error|", "99th pct |error|",
                                    "mean top-15 overlap (of 15)"],
                         "value": [int(test.sum()), err.mean(), np.percentile(err, 99), overlap.mean()]})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", default=str(RESULTS / "sites_human.tsv"))
    ap.add_argument("--matrices", default=str(KINOME_MATRICES))
    ap.add_argument("--reference", default=str(KINOME_REFERENCE))
    ap.add_argument("--top", type=int, default=15)
    args = ap.parse_args()

    ref = pd.read_csv(args.reference)
    atlas = KinomeAtlas.from_tables(
        pd.read_excel(args.matrices, sheet_name="ser_thr_all_norm_scaled_matrice", index_col=0), ref)
    print(f"loaded {len(atlas.kinases)} kinase matrices and {len(ref):,} reference sites")

    check = reference_check(atlas, ref)
    check.to_csv(out(RESULTS / "reference_check.tsv"), sep="\t", index=False)
    print(check.to_string(index=False))

    sites = pd.read_csv(args.sites, sep="\t").dropna(subset=["human_window"])
    fav = atlas.favored_kinases(sites, window_col="human_window", id_col="site_id", top_n=args.top)
    fav.to_csv(out(RESULTS / "favored_kinases.tsv.gz"), sep="\t", index=False)
    print(f"scored {fav['site_id'].nunique():,} S/T sites -> {len(fav):,} site-kinase pairs")


if __name__ == "__main__":
    main()
