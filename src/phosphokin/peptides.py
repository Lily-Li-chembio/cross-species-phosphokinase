"""From modified peptides with localization probabilities to protein phosphosites.

Search engines such as MaxQuant report phosphopeptides with per-residue
localization probabilities written inline, e.g. ``AAS(0.98)PT(0.02)K``.
This module turns such strings into protein-level sites like ``S123`` with the
probability attached, and filters out poorly localized sites.
"""

from __future__ import annotations

import re

import pandas as pd

_MOD = re.compile(r"([A-Z])\(([0-9.eE+-]+)\)")


def parse_probability_sequence(annotated: str) -> tuple[str, list[tuple[str, int, float]]]:
    """Split ``AAS(0.98)PT(0.02)K`` into the plain peptide and its candidate sites.

    Returns
    -------
    peptide
        The unmodified peptide, e.g. ``AASPTK``.
    sites
        ``(residue, position_in_peptide, probability)`` for every annotated
        residue with probability > 0. Positions are 1-based.
    """
    plain: list[str] = []
    sites: list[tuple[str, int, float]] = []
    pos = 0
    last = 0
    for match in _MOD.finditer(annotated):
        chunk = annotated[last:match.start()]
        plain.append(chunk)
        pos += len(chunk)
        residue, prob = match.group(1), float(match.group(2))
        plain.append(residue)
        pos += 1
        if prob > 0:
            sites.append((residue, pos, prob))
        last = match.end()
    plain.append(annotated[last:])
    peptide = "".join(plain)
    if not peptide.isalpha():
        raise ValueError(f"could not parse peptide annotation: {annotated!r}")
    return peptide, sites


def locate_peptide(peptide: str, sequence: str) -> list[int]:
    """All 0-based start offsets of ``peptide`` in ``sequence`` (overlaps allowed)."""
    starts, i = [], sequence.find(peptide)
    while i != -1:
        starts.append(i)
        i = sequence.find(peptide, i + 1)
    return starts


def peptides_to_sites(peptides: pd.DataFrame, proteins: pd.DataFrame,
                      annotated_col: str = "annotated_peptide",
                      accession_col: str = "accession",
                      min_probability: float = 0.75,
                      sequence_overrides: dict[str, str] | None = None) -> pd.DataFrame:
    """Map annotated phosphopeptides onto their proteins.

    Parameters
    ----------
    peptides
        One row per (peptide, protein) pair with an annotated peptide string.
    proteins
        Protein table from :func:`phosphokin.sequences.load_uniprot`.
    min_probability
        Keep sites whose localization probability is at least this value
        (0.75 is the common "class I" cutoff).
    sequence_overrides
        Optional ``{accession: sequence}`` for isoforms that are not the UniProt
        canonical sequence (the peptide is looked up in this sequence instead).

    Returns
    -------
    DataFrame with ``accession, gene, site, residue, position, probability,
    peptide`` plus all original columns. Peptides that cannot be found in their
    protein are reported in ``DataFrame.attrs['unmapped']``.
    """
    seqs = dict(zip(proteins["accession"], proteins["sequence"]))
    genes = dict(zip(proteins["accession"], proteins["gene"]))
    if sequence_overrides:
        seqs.update(sequence_overrides)

    rows, unmapped = [], []
    for rec in peptides.to_dict("records"):
        acc = rec[accession_col]
        peptide, sites = parse_probability_sequence(rec[annotated_col])
        seq = seqs.get(acc)
        starts = locate_peptide(peptide, seq) if seq else []
        if not starts:
            unmapped.append((acc, peptide))
            continue
        for start in starts:
            for residue, p_pos, prob in sites:
                if prob < min_probability:
                    continue
                position = start + p_pos
                rows.append({**rec, "accession": acc, "gene": genes.get(acc),
                             "peptide": peptide, "residue": residue,
                             "position": position, "site": f"{residue}{position}",
                             "probability": prob})
    out = pd.DataFrame(rows)
    out.attrs["unmapped"] = unmapped
    return out
