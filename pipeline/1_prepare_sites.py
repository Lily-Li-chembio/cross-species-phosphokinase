"""Step 1 - Build a clean phosphosite table from a MaxQuant-style site table.

Works with MaxQuant's ``Phospho (STY)Sites.txt`` and with the mouse tissue atlas
(Giansanti et al., Nat. Methods 2022, Supplementary Table 2), which uses the
same columns. Keeps singly-phosphorylated, well-localized sites.

Outputs
  results/sites.tsv         one row per site: site_id, accession, gene, residue, position, window
  results/intensities.tsv   site_id x sample log10 intensities (if intensity columns exist)
"""

from __future__ import annotations

import argparse

import pandas as pd

from _common import MOUSE_ATLAS, RESULTS, out
from phosphokin.sequences import center_window


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", default=str(MOUSE_ATLAS), help="site table (.xlsx or tab-separated .txt)")
    ap.add_argument("--sheet", default="Supplementary Table 2 - Organs", help="Excel sheet name")
    ap.add_argument("--min-prob", type=float, default=0.75, help="minimum localization probability")
    ap.add_argument("--intensity-prefix", default="Intensity (Log10) ",
                    help="prefix of per-sample intensity columns")
    args = ap.parse_args()

    if args.input.endswith((".xlsx", ".xls")):
        raw = pd.read_excel(args.input, sheet_name=args.sheet)
    else:
        raw = pd.read_csv(args.input, sep="\t", low_memory=False)
    print(f"read {len(raw):,} site rows")

    df = raw[raw["Localization prob"] >= args.min_prob]
    if "Multiplicity" in df.columns:          # keep singly phosphorylated peptides only
        df = df[df["Multiplicity"].astype(str).str.endswith("_1")]
    df = df[df["Amino acid"].isin(["S", "T", "Y"])]
    df = df.dropna(subset=["Position", "Leading proteins", "Sequence window"])

    sites = pd.DataFrame({
        "accession": df["Leading proteins"].astype(str).str.split(";").str[0],
        "gene": df["Gene names"].astype("string").str.split(";").str[0],
        "residue": df["Amino acid"],
        "position": df["Position"].astype(int),
        "probability": df["Localization prob"],
        # MaxQuant windows are ±15 and may list several proteins separated by ';'
        "window": df["Sequence window"].astype(str).str.split(";").str[0].map(center_window),
    }, index=df.index)
    sites["site_id"] = sites["accession"] + "_" + sites["residue"] + sites["position"].astype(str)
    sites = sites.dropna(subset=["gene"]).drop_duplicates("site_id")

    icols = [c for c in df.columns if c.startswith(args.intensity_prefix)]
    if icols:
        inten = df.loc[sites.index, icols].copy()
        inten.columns = [c[len(args.intensity_prefix):] for c in icols]
        inten.insert(0, "site_id", sites["site_id"])
        inten.to_csv(out(RESULTS / "intensities.tsv"), sep="\t", index=False)

    cols = ["site_id", "accession", "gene", "residue", "position", "probability", "window"]
    sites[cols].to_csv(out(RESULTS / "sites.tsv"), sep="\t", index=False)
    print(f"kept {len(sites):,} singly-phosphorylated sites (prob >= {args.min_prob}) "
          f"on {sites['accession'].nunique():,} proteins; {len(icols)} samples")


if __name__ == "__main__":
    main()
