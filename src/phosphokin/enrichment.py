"""Kinase motif enrichment (Johnson et al., Nature 2023, "Kinase enrichment analysis").

For each kinase, compare how often it is among the top-ranked (biochemically
favoured) kinases for *up*- or *down*-regulated sites versus *unregulated* sites:

    frequency factor = (fav_up / n_up) / (fav_unreg / n_unreg)

Significance is a one-sided Fisher's exact test per direction, with
Benjamini-Hochberg correction across kinases. When a contingency cell is zero
the frequency factor uses the Haldane correction (+0.5 to every cell).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def bh_adjust(pvalues) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaNs are kept as NaN)."""
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    q = p[ok]
    n = len(q)
    if n == 0:
        return out
    order = np.argsort(q)
    scaled = q[order] * n / np.arange(1, n + 1)
    scaled = np.minimum.accumulate(scaled[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(scaled, 1.0)
    out[ok] = adj
    return out


def classify_sites(log2fc, pvalue=None, fc_threshold: float = 1.0,
                   p_threshold: float = 0.05) -> np.ndarray:
    """Label sites 'up', 'down' or 'unregulated'.

    Up/down need ``|log2FC| >= fc_threshold`` and, if p-values are given,
    ``p < p_threshold``. Everything else is 'unregulated'.
    """
    fc = np.asarray(log2fc, dtype=float)
    sig = np.ones(fc.shape, bool) if pvalue is None else (np.asarray(pvalue, float) < p_threshold)
    labels = np.full(fc.shape, "unregulated", dtype=object)
    labels[(fc >= fc_threshold) & sig] = "up"
    labels[(fc <= -fc_threshold) & sig] = "down"
    labels[np.isnan(fc)] = None
    return labels


def _freq_factor(a, b, c, d):
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return (a / (a + b)) / (c / (c + d))


def kinase_enrichment(site_status: pd.DataFrame, favored: pd.DataFrame,
                      id_col: str = "site_id", status_col: str = "status") -> pd.DataFrame:
    """Kinase enrichment among up- and down-regulated sites.

    Parameters
    ----------
    site_status
        One row per site with ``id_col`` and ``status_col`` in {'up', 'down', 'unregulated'}.
        Only sites present in ``favored`` (i.e. scoreable S/T sites) are counted.
    favored
        Output of :meth:`KinomeAtlas.favored_kinases` (``id_col``, ``kinase``).

    Returns
    -------
    One row per kinase with counts, ``log2_ff_up/down``, one-sided p-values for
    enrichment, BH-adjusted p-values, and the most significant direction.
    """
    scored = set(favored[id_col])
    st = site_status[site_status[id_col].isin(scored)].drop_duplicates(id_col)
    groups = {s: set(st.loc[st[status_col] == s, id_col]) for s in ("up", "down", "unregulated")}
    n = {s: len(v) for s, v in groups.items()}
    by_kinase = favored.groupby("kinase")[id_col].agg(set)

    rows = []
    for kinase, sites in by_kinase.items():
        fav = {s: len(sites & groups[s]) for s in groups}
        row = {"kinase": kinase, "n_up": n["up"], "n_down": n["down"], "n_unregulated": n["unregulated"],
               "fav_up": fav["up"], "fav_down": fav["down"], "fav_unregulated": fav["unregulated"]}
        c, d = fav["unregulated"], n["unregulated"] - fav["unregulated"]
        for s in ("up", "down"):
            a, b = fav[s], n[s] - fav[s]
            if n[s] == 0 or n["unregulated"] == 0:
                row[f"log2_ff_{s}"], row[f"p_{s}"] = np.nan, np.nan
                continue
            row[f"log2_ff_{s}"] = np.log2(_freq_factor(a, b, c, d))
            row[f"p_{s}"] = fisher_exact([[a, b], [c, d]], alternative="greater")[1]
        rows.append(row)

    res = pd.DataFrame(rows)
    for s in ("up", "down"):
        res[f"padj_{s}"] = bh_adjust(res[f"p_{s}"])
    # Most significant side (as in the paper's volcano plots): down is shown as negative.
    use_down = res["padj_down"].fillna(1) < res["padj_up"].fillna(1)
    res["direction"] = np.where(use_down, "down", "up")
    res["log2_ff"] = np.where(use_down, -res["log2_ff_down"], res["log2_ff_up"])
    res["padj"] = np.where(use_down, res["padj_down"], res["padj_up"])
    return res.sort_values("padj").reset_index(drop=True)
