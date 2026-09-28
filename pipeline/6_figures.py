"""Step 6 - Figures: a volcano plot per chosen comparison and a summary heatmap.

Outputs
  docs/figures/volcano_<comparison>.png
  docs/figures/heatmap.png
"""

from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from _common import FIGURES, RESULTS, out  # noqa: E402
from phosphokin.plots import enrichment_heatmap, volcano  # noqa: E402

DEMO_TISSUES = ["Frontal.lobe", "Hippocampus", "Cerebellum", "Spinal.cord", "Testes", "Thymus",
                "Spleen", "Lung", "Urinary.bladder", "Heart", "Skeletal.muscle", "Liver", "Kidneys"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--enrichment", default=str(RESULTS / "kinase_enrichment.tsv"))
    ap.add_argument("--volcano", nargs="*", default=["Hippocampus", "Urinary.bladder"])
    ap.add_argument("--heatmap", nargs="*", default=DEMO_TISSUES)
    ap.add_argument("--padj", type=float, default=0.1)
    ap.add_argument("--top", type=int, default=4, help="kinases per comparison in the heatmap")
    args = ap.parse_args()

    res = pd.read_csv(args.enrichment, sep="\t")
    by_comp = {c: t for c, t in res.groupby("comparison")}

    for comp in args.volcano:
        fig, ax = plt.subplots(figsize=(5.2, 4.4), dpi=150)
        volcano(by_comp[comp], title=f"{comp.replace('.', ' ')} vs other tissues", padj_cutoff=args.padj, ax=ax)
        fig.tight_layout()
        fig.savefig(out(FIGURES / f"volcano_{comp}.png"))
        plt.close(fig)

    chosen = {c: by_comp[c] for c in args.heatmap if c in by_comp}
    fig, ax = plt.subplots(figsize=(7.2, 7.6), dpi=150)
    enrichment_heatmap(chosen, padj_cutoff=args.padj, top=args.top, ax=ax)
    ax.set_title("Kinase motifs enriched in tissue-specific phosphosites", fontsize=11, loc="left")
    fig.tight_layout()
    fig.savefig(out(FIGURES / "heatmap.png"))
    plt.close(fig)
    print(f"wrote figures to {FIGURES}")


if __name__ == "__main__":
    main()
