"""Step 2 - Transfer mouse phosphosites to their human orthologous residues.

Pairs each mouse gene with a human ortholog, aligns the two proteins
(global Needleman-Wunsch / Gotoh, BLOSUM62) and keeps sites whose aligned human
residue is phospho-compatible (S<->S/T, T<->S/T, Y<->Y).

Outputs
  results/sites_human.tsv
"""

from __future__ import annotations

import argparse
import time

import pandas as pd

from _common import HUMAN_UNIPROT, MOUSE_UNIPROT, RESULTS, out
from phosphokin.orthology import (map_sites_to_human, orthologs_by_symbol, orthologs_from_table,
                                  orthologs_gprofiler)
from phosphokin.sequences import load_uniprot


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", default=str(RESULTS / "sites.tsv"))
    ap.add_argument("--mouse-uniprot", default=str(MOUSE_UNIPROT))
    ap.add_argument("--human-uniprot", default=str(HUMAN_UNIPROT))
    ap.add_argument("--orthologs", default="symbol",
                    help="'symbol' (offline), 'gprofiler' (web), or a TSV with columns source_gene, human_gene")
    ap.add_argument("--method", choices=["alignment", "window"], default="alignment")
    args = ap.parse_args()

    sites = pd.read_csv(args.sites, sep="\t")
    mouse = load_uniprot(args.mouse_uniprot)
    human = load_uniprot(args.human_uniprot)

    if args.orthologs == "symbol":
        orth = orthologs_by_symbol(sites["gene"], human)
    elif args.orthologs == "gprofiler":
        orth = orthologs_gprofiler(sites["gene"])
    else:
        orth = orthologs_from_table(pd.read_csv(args.orthologs, sep="\t"), "source_gene", "human_gene")
    print(f"{len(orth):,} of {sites['gene'].nunique():,} genes have a human ortholog ({args.orthologs})")

    t0 = time.time()
    mapped = map_sites_to_human(sites, mouse, human, orth, method=args.method, progress=True)
    ok = mapped["human_site"].notna()
    print(f"mapped {ok.sum():,} / {len(mapped):,} sites to human in {time.time() - t0:.0f}s; "
          f"median ±7 window identity {mapped.loc[ok, 'window_identity'].median():.2f}")

    mapped["human_site_id"] = mapped["human_accession"] + "_" + mapped["human_site"]
    mapped.to_csv(out(RESULTS / "sites_human.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
