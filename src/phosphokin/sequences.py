"""Protein sequence tables and phosphosite sequence windows."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# Column names used across the package for a protein table.
PROTEIN_COLUMNS = ["accession", "gene", "reviewed", "sequence"]

_UNIPROT_RENAME = {
    "Entry": "accession",
    "Gene Names (primary)": "gene",
    "Gene Names": "gene",
    "Reviewed": "reviewed",
    "Sequence": "sequence",
}


def load_uniprot(path: str | Path) -> pd.DataFrame:
    """Load a UniProt export (TSV, CSV or Excel) into a tidy protein table.

    The export needs the columns *Entry*, *Gene Names (primary)*, *Reviewed* and
    *Sequence* (the default names from https://www.uniprot.org downloads).
    Returns a table with columns ``accession, gene, reviewed, sequence`` where
    ``reviewed`` is a boolean (Swiss-Prot entries are True).
    """
    path = Path(path)
    name = path.name.lower()
    if name.endswith((".xlsx", ".xls")):
        df = pd.read_excel(path)
    elif name.endswith((".tsv", ".tsv.gz", ".tab")):
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)
    df = df.rename(columns=_UNIPROT_RENAME)
    missing = [c for c in PROTEIN_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing UniProt columns: {missing}")
    df = df[PROTEIN_COLUMNS].copy()
    df["reviewed"] = df["reviewed"].astype(str).str.lower().eq("reviewed")
    # "Gene Names" can hold several synonyms; the first one is the primary symbol.
    df["gene"] = df["gene"].astype("string").str.split(r"[ ;]").str[0]
    return df.dropna(subset=["sequence"]).reset_index(drop=True)


def canonical_by_gene(proteins: pd.DataFrame) -> pd.DataFrame:
    """One protein per gene symbol: reviewed (Swiss-Prot) first, then the longest."""
    p = proteins.dropna(subset=["gene"]).copy()
    p["_len"] = p["sequence"].str.len()
    p = p.sort_values(["gene", "reviewed", "_len"], ascending=[True, False, False])
    return p.drop_duplicates("gene").drop(columns="_len").reset_index(drop=True)


def site_window(sequence: str, position: int, flank: int = 7, pad: str = "_") -> str:
    """Return the residues ``position-flank .. position+flank`` (1-based, inclusive).

    Positions outside the protein are filled with ``pad`` so the window always has
    length ``2*flank + 1`` and the phosphosite sits in the middle.
    """
    if not 1 <= position <= len(sequence):
        raise ValueError(f"position {position} outside sequence of length {len(sequence)}")
    start, end = position - 1 - flank, position + flank
    left = pad * max(0, -start)
    right = pad * max(0, end - len(sequence))
    return left + sequence[max(0, start):min(len(sequence), end)] + right


def center_window(window: str, flank: int = 7) -> str:
    """Trim a longer centred window (e.g. MaxQuant's ±15 'Sequence window') to ±flank."""
    mid = len(window) // 2
    return window[mid - flank: mid + flank + 1]
