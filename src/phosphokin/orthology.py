"""Transfer phosphosites from a model organism (e.g. mouse) to human orthologs.

Workflow
--------
1. Pair each source gene with a human ortholog (g:Profiler, a user table such as
   MGI/Ensembl, or a same-symbol fallback).
2. Pick one representative protein per gene (Swiss-Prot first).
3. Globally align each source/human protein pair and move every site to the
   aligned human residue. A site is transferred only if the human residue is
   phospho-compatible (S<->S/T, T<->S/T, Y<->Y).

An alignment-free alternative (``method="window"``) looks for the human position
whose ±7 residue window best matches the source window. It is useful for very
long proteins and as a cross-check of the alignment.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd

from .align import aligned_identity, align_position_map
from .sequences import canonical_by_gene, site_window

COMPATIBLE = {"S": "ST", "T": "ST", "Y": "Y"}


# --------------------------------------------------------------------------- #
# Ortholog tables
# --------------------------------------------------------------------------- #
def orthologs_by_symbol(source_genes: Iterable[str], human_proteins: pd.DataFrame) -> pd.DataFrame:
    """Same-symbol orthologs (mouse ``Braf`` -> human ``BRAF``).

    A quick, offline approximation: most mouse/human one-to-one orthologs share
    a symbol. Prefer :func:`orthologs_gprofiler` or a curated table when you can.
    """
    human = set(human_proteins["gene"].dropna().astype(str))
    rows = [(g, g.upper()) for g in pd.unique(pd.Series(list(source_genes)).dropna().astype(str))
            if g.upper() in human]
    return pd.DataFrame(rows, columns=["source_gene", "human_gene"])


def orthologs_from_table(table: pd.DataFrame, source_col: str, human_col: str) -> pd.DataFrame:
    """Use a curated ortholog table (e.g. MGI ``HOM_MouseHumanSequence.rpt``)."""
    out = table[[source_col, human_col]].dropna().astype(str)
    out.columns = ["source_gene", "human_gene"]
    return out.drop_duplicates().reset_index(drop=True)


def orthologs_gprofiler(source_genes: Iterable[str], organism: str = "mmusculus") -> pd.DataFrame:
    """Orthologs from the g:Profiler web service (needs ``pip install gprofiler-official``)."""
    from gprofiler import GProfiler  # optional dependency, imported lazily

    genes = list(pd.unique(pd.Series(list(source_genes)).dropna().astype(str)))
    res = GProfiler(return_dataframe=True).orth(genes, organism=organism, target="hsapiens")
    res = res[res["name"].notna() & (res["name"] != "N/A")]
    out = res[["incoming", "name"]].rename(columns={"incoming": "source_gene", "name": "human_gene"})
    return out.drop_duplicates().reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Site transfer
# --------------------------------------------------------------------------- #
def _window_identity(w1: str, w2: str) -> float:
    pairs = [(a, b) for a, b in zip(w1, w2) if a != "_" and b != "_"]
    return sum(a == b for a, b in pairs) / len(pairs) if pairs else 0.0


def _best_window_match(window: str, residue: str, human_seq: str, flank: int) -> tuple[int, float]:
    """Human position (1-based) whose window best matches ``window``; 0 if ambiguous."""
    seq = np.frombuffer(human_seq.encode(), dtype=np.uint8)
    padded = np.concatenate([np.full(flank, ord("_"), np.uint8), seq, np.full(flank, ord("_"), np.uint8)])
    windows = np.lib.stride_tricks.sliding_window_view(padded, 2 * flank + 1)
    query = np.frombuffer(window.encode(), dtype=np.uint8)
    ok = np.isin(seq, np.frombuffer(COMPATIBLE[residue].encode(), dtype=np.uint8))
    if not ok.any():
        return 0, 0.0
    valid = (windows != ord("_")) & (query != ord("_"))
    ident = ((windows == query) & valid).sum(1) / np.maximum(valid.sum(1), 1)
    ident = np.where(ok, ident, -1.0)
    best = ident.max()
    hits = np.flatnonzero(ident == best)
    return (int(hits[0]) + 1, float(best)) if len(hits) == 1 else (0, float(best))


def _check_source_position(seq: str, pos: int, window: str, flank: int) -> int:
    """Confirm a site against its measured sequence window; relocate it if the
    reference sequence version is shifted. Returns the 1-based position or 0."""
    if 1 <= pos <= len(seq) and site_window(seq, pos, flank) == window:
        return pos
    core = window.strip("_")
    offset = len(window) - len(window.lstrip("_"))       # padding before the core
    hits = [h for h in _find_all(seq, core)]
    if len(hits) == 1:
        return hits[0] + (flank - offset) + 1
    return 0


def _find_all(seq: str, sub: str) -> list[int]:
    out, i = [], seq.find(sub)
    while i != -1:
        out.append(i)
        i = seq.find(sub, i + 1)
    return out


def map_sites_to_human(sites: pd.DataFrame, source_proteins: pd.DataFrame,
                       human_proteins: pd.DataFrame, orthologs: pd.DataFrame,
                       method: str = "alignment", flank: int = 7,
                       gap_open: int = 11, gap_extend: int = 1,
                       max_length: int = 8000, progress: bool = False) -> pd.DataFrame:
    """Transfer source-organism sites to human.

    Parameters
    ----------
    sites
        Needs ``accession`` (source UniProt), ``gene``, ``residue`` and ``position``.
        If a ``window`` column (the measured ±7 sequence) is present, each site
        is checked against the reference sequence and re-located when the
        sequence version differs (the corrected position is in ``source_position``).
    source_proteins, human_proteins
        Protein tables from :func:`phosphokin.sequences.load_uniprot`.
    orthologs
        Two columns ``source_gene, human_gene``.
    method
        ``"alignment"`` (global alignment; default) or ``"window"``. With
        ``"alignment"``, proteins longer than ``max_length`` fall back to ``"window"``.

    Returns
    -------
    ``sites`` plus ``human_gene, human_accession, human_position, human_site,
    human_window, window_identity, protein_identity, mapping``. Unmapped sites
    keep NaN in the human columns.
    """
    if method not in {"alignment", "window"}:
        raise ValueError("method must be 'alignment' or 'window'")
    src_seq = dict(zip(source_proteins["accession"], source_proteins["sequence"]))
    human = canonical_by_gene(human_proteins)
    h_acc = dict(zip(human["gene"], human["accession"]))
    h_seq = dict(zip(human["gene"], human["sequence"]))
    orth = orthologs.drop_duplicates("source_gene").set_index("source_gene")["human_gene"]

    out = sites.copy()
    cols = ["human_gene", "human_accession", "human_position", "human_site",
            "human_window", "window_identity", "protein_identity", "mapping"]
    found: dict = {}
    has_window = "window" in out.columns

    groups = out.groupby(["accession", "gene"], sort=False).groups
    n_done = 0
    for (acc, gene), idx in groups.items():
        hg = orth.get(gene)
        s_seq, hs = src_seq.get(acc), h_seq.get(hg) if hg is not None else None
        n_done += 1
        if progress and n_done % 500 == 0:
            print(f"  {n_done}/{len(groups)} proteins")
        if s_seq is None or hs is None:
            continue
        use_align = method == "alignment" and max(len(s_seq), len(hs)) <= max_length
        if use_align:
            pos_map, _ = align_position_map(s_seq, hs, gap_open, gap_extend)
            prot_id = aligned_identity(s_seq, hs, pos_map)
        for i in idx:
            res, pos = out.at[i, "residue"], int(out.at[i, "position"])
            if has_window:
                pos = _check_source_position(s_seq, pos, str(out.at[i, "window"]), flank)
            elif pos > len(s_seq) or s_seq[pos - 1] != res:
                pos = 0
            if not pos:
                continue  # site does not match this version of the source sequence
            s_win = site_window(s_seq, pos, flank)
            if use_align:
                hp = int(pos_map[pos])
                label = "alignment"
            else:
                hp, _ = _best_window_match(s_win, res, hs, flank)
                prot_id, label = np.nan, "window"
            if hp == 0 or hs[hp - 1] not in COMPATIBLE[res]:
                continue
            h_win = site_window(hs, hp, flank)
            found[i] = (pos, hg, h_acc[hg], hp, f"{hs[hp - 1]}{hp}", h_win,
                        round(_window_identity(s_win, h_win), 3), prot_id, label)
    mapped = pd.DataFrame.from_dict(found, orient="index", columns=["source_position"] + cols)
    out = out.drop(columns=[c for c in ["source_position"] + cols if c in out.columns]).join(mapped)
    out["human_position"] = out["human_position"].astype("Int64")
    out["source_position"] = out["source_position"].astype("Int64")
    return out
