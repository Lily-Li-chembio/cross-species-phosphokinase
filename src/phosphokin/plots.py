"""Figures: kinase-enrichment volcano plot and a multi-condition heatmap."""

from __future__ import annotations

import numpy as np
import pandas as pd

UP = "#e34948"        # enriched among up-regulated sites
DOWN = "#2a78d6"      # enriched among down-regulated sites
NEUTRAL = "#b9b8b3"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e6e5e1"


def _style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_2)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK_2, labelsize=9)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def _label_points(ax, rows: pd.DataFrame) -> None:
    """Direct labels with a simple greedy nudge so labels do not overlap."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    placed = []
    for _, r in rows.iterrows():
        for dy in (3, -11, 13, -21, 23, -31, 33):
            txt = ax.annotate(r["kinase"], (r["log2_ff"], r["y"]), xytext=(5, dy),
                              textcoords="offset points", fontsize=8, color=INK)
            box = txt.get_window_extent(renderer).expanded(1.05, 1.15)
            if not any(box.overlaps(b) for b in placed):
                placed.append(box)
                break
            txt.remove()


def volcano(enrichment: pd.DataFrame, title: str = "", padj_cutoff: float = 0.1,
            n_labels: int = 8, ax=None):
    """Volcano plot of :func:`phosphokin.enrichment.kinase_enrichment` output.

    x: log2 frequency factor (negative = enriched among down-regulated sites),
    y: -log10 adjusted p-value.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(5.2, 4.4), dpi=150)
    df = enrichment.dropna(subset=["log2_ff", "padj"]).copy()
    df["y"] = -np.log10(df["padj"].clip(lower=1e-300))
    sig = (df["padj"] <= padj_cutoff) & (df["log2_ff"].abs() > 0)
    color = np.where(~sig, NEUTRAL, np.where(df["direction"] == "up", UP, DOWN))
    _style(ax)
    ax.scatter(df["log2_ff"], df["y"], s=22, c=color, edgecolors="white", linewidths=0.6, zorder=3)
    ax.axhline(-np.log10(padj_cutoff), color=INK_2, lw=0.8, ls=(0, (3, 3)), zorder=2)
    ax.axvline(0, color=INK_2, lw=0.8, zorder=2)
    _label_points(ax, df[sig].sort_values("padj").head(n_labels))
    ax.set_xlabel("log2 frequency factor  (← down-regulated · up-regulated →)", color=INK_2, fontsize=9)
    ax.set_ylabel("−log10 adjusted p", color=INK_2, fontsize=9)
    if title:
        ax.set_title(title, color=INK, fontsize=11, loc="left")
    return ax


def enrichment_heatmap(results: dict[str, pd.DataFrame], kinases: list[str] | None = None,
                       padj_cutoff: float = 0.1, top: int = 5, ax=None):
    """Heatmap of signed log2 frequency factors (kinases x conditions).

    ``results`` maps a condition name to its enrichment table. Cells with
    ``padj > padj_cutoff`` are drawn neutral. By default shows the ``top`` most
    significant kinases of each condition.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    if kinases is None:
        chosen: list[str] = []
        for res in results.values():
            for k in res[res["padj"] <= padj_cutoff].sort_values("padj")["kinase"].head(top):
                if k not in chosen:
                    chosen.append(k)
        kinases = chosen
    conds = list(results)
    mat = np.full((len(kinases), len(conds)), np.nan)
    for j, c in enumerate(conds):
        r = results[c].set_index("kinase")
        for i, k in enumerate(kinases):
            if k in r.index and r.at[k, "padj"] <= padj_cutoff:
                mat[i, j] = r.at[k, "log2_ff"]

    cmap = LinearSegmentedColormap.from_list("div", [DOWN, "#f0efec", UP])
    cmap.set_bad("#f7f6f3")
    lim = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1.0
    if ax is None:
        _, ax = plt.subplots(figsize=(0.55 * len(conds) + 2.2, 0.28 * len(kinases) + 1.6), dpi=150)
    im = ax.imshow(np.ma.masked_invalid(mat), cmap=cmap, norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
    ax.set_xticks(range(len(conds)), [c.replace(".", " ") for c in conds], rotation=45, ha="right",
                  fontsize=8, color=INK_2)
    ax.set_yticks(range(len(kinases)), kinases, fontsize=8, color=INK_2)
    ax.set_xticks(np.arange(-0.5, len(conds)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(kinases)), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("signed log2 frequency factor", fontsize=8, color=INK_2)
    cb.ax.tick_params(labelsize=7, colors=INK_2)
    cb.outline.set_visible(False)
    return ax
