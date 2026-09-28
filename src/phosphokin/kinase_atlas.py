"""Score phosphosites with the Ser/Thr kinome atlas (Johnson et al., Nature 2023).

The atlas provides, for each of 303 Ser/Thr kinases, a position-specific matrix
over residues -5..-1 and +1..+4 (20 amino acids plus phospho-S/T/Y written in
lower case). A site's raw score for a kinase is the product of the matrix values
of the residues around it. Because raw scores are not comparable between
kinases, each score is converted to a *percentile* of that kinase's score
distribution over the human reference phosphoproteome, and kinases are then
ranked by percentile for every site.

Percentiles here are obtained by interpolating against the published reference
table (Supplementary Table 3), separately for serine- and threonine-centred
sites. This reproduces the published percentiles (held-out mean error
~0.001 percentile points) without re-deriving the kinase-specific S/T
preference factors.

The atlas tables are not redistributed with this package; download them from
the paper's supplementary information (see README).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

POSITIONS = [-5, -4, -3, -2, -1, 1, 2, 3, 4]
RESIDUES = list("PGACSTVILMFYWHKRQNDEsty")
_UNKNOWN = len(RESIDUES)          # '_' padding, X, U ... contribute a neutral factor (log2 = 0)
_FLANK = 7                        # windows are ±7 around the site (15 residues)

_LUT = np.full(256, _UNKNOWN, dtype=np.int64)
for _i, _aa in enumerate(RESIDUES):
    _LUT[ord(_aa)] = _i


def encode_windows(windows) -> np.ndarray:
    """Encode ±7 windows (15 characters, site in the middle) as residue indices at POSITIONS."""
    w = pd.Series(windows).astype(str).to_numpy()
    if len(w) == 0:
        return np.zeros((0, len(POSITIONS)), dtype=np.int64)
    lengths = {len(x) for x in w}
    if lengths != {2 * _FLANK + 1}:
        raise ValueError(f"windows must be {2 * _FLANK + 1} residues long, got lengths {sorted(lengths)}")
    arr = np.frombuffer("".join(w).encode("ascii"), dtype=np.uint8).reshape(len(w), 2 * _FLANK + 1)
    return _LUT[arr[:, [_FLANK + p for p in POSITIONS]]]


@dataclass
class KinomeAtlas:
    """Kinase matrices plus the reference phosphoproteome used for percentiles."""

    kinases: list[str]
    log_matrix: np.ndarray                     # kinases x positions x (residues + unknown), log2
    reference: dict = field(default_factory=dict)  # residue -> (list of x arrays, list of y arrays)

    # ------------------------------------------------------------------ #
    @classmethod
    def from_tables(cls, matrices: pd.DataFrame, reference_sites: pd.DataFrame | None = None,
                    window_col: str = "SITE_+/-7_AA") -> "KinomeAtlas":
        """Build from the matrix table (kinases x '<pos><aa>' columns) and, optionally,
        the reference-site table with ``<KINASE>_percentile`` columns."""
        kinases = [str(k) for k in matrices.index]
        logm = np.zeros((len(kinases), len(POSITIONS), len(RESIDUES) + 1))
        for pi, p in enumerate(POSITIONS):
            for ai, aa in enumerate(RESIDUES):
                col = f"{p}{aa}"
                if col not in matrices.columns:
                    raise ValueError(f"matrix column {col!r} missing")
                logm[:, pi, ai] = np.log2(matrices[col].to_numpy(dtype=float))
        atlas = cls(kinases=kinases, log_matrix=logm)
        if reference_sites is not None:
            atlas.set_reference(reference_sites, window_col=window_col)
        return atlas

    @classmethod
    def from_supplementary(cls, matrix_xlsx: str | Path, reference_csv: str | Path | None = None,
                           sheet: str = "ser_thr_all_norm_scaled_matrice") -> "KinomeAtlas":
        """Load Supplementary Tables 2 (matrices) and 3 (reference sites) of Johnson et al. 2023."""
        matrices = pd.read_excel(matrix_xlsx, sheet_name=sheet, index_col=0)
        ref = pd.read_csv(reference_csv) if reference_csv is not None else None
        return cls.from_tables(matrices, ref)

    # ------------------------------------------------------------------ #
    def score(self, windows) -> np.ndarray:
        """log2 raw scores, shape (n_sites, n_kinases)."""
        idx = encode_windows(windows)
        out = np.zeros((len(idx), len(self.kinases)))
        for pi in range(len(POSITIONS)):
            out += self.log_matrix[:, pi, :].T[idx[:, pi]]
        return out

    def set_reference(self, reference_sites: pd.DataFrame, window_col: str = "SITE_+/-7_AA") -> None:
        """Store score -> published-percentile curves for S- and T-centred reference sites."""
        pct_cols = [f"{k}_percentile" for k in self.kinases]
        missing = [c for c in pct_cols if c not in reference_sites.columns]
        if missing:
            raise ValueError(f"reference table lacks {len(missing)} percentile columns, e.g. {missing[:3]}")
        windows = reference_sites[window_col].astype(str)
        center = windows.str[_FLANK].str.upper()
        scores = self.score(windows)
        pcts = reference_sites[pct_cols].to_numpy(dtype=float)
        self.reference = {}
        for res in ("S", "T"):
            m = (center == res).to_numpy()
            xs, ys = [], []
            for k in range(len(self.kinases)):
                x, y = scores[m, k], pcts[m, k]
                ok = ~np.isnan(y)
                ux, inv = np.unique(x[ok], return_inverse=True)
                xs.append(ux)
                ys.append(np.bincount(inv, y[ok]) / np.bincount(inv))
            self.reference[res] = (xs, ys)

    def percentile(self, windows) -> np.ndarray:
        """Percentile of each site for each kinase, shape (n_sites, n_kinases).

        Sites centred on Y (not covered by the Ser/Thr atlas) get NaN.
        """
        if not self.reference:
            raise RuntimeError("no reference loaded; call set_reference() first")
        windows = pd.Series(windows).astype(str).reset_index(drop=True)
        scores = self.score(windows)
        center = windows.str[_FLANK].str.upper().to_numpy()
        out = np.full(scores.shape, np.nan)
        for res, (xs, ys) in self.reference.items():
            m = center == res
            if not m.any():
                continue
            for k in range(len(self.kinases)):
                out[m, k] = np.interp(scores[m, k], xs[k], ys[k])
        return out

    def rank(self, percentiles: np.ndarray) -> np.ndarray:
        """Rank of each kinase for each site (1 = highest percentile)."""
        order = np.argsort(-np.nan_to_num(percentiles, nan=-1.0), axis=1, kind="stable")
        ranks = np.empty_like(order)
        rows = np.arange(order.shape[0])[:, None]
        ranks[rows, order] = np.arange(1, order.shape[1] + 1)
        return ranks

    def favored_kinases(self, sites: pd.DataFrame, window_col: str = "human_window",
                        id_col: str = "site_id", top_n: int = 15) -> pd.DataFrame:
        """Long table of the top-``top_n`` kinases for every S/T site.

        Columns: ``site_id, kinase, percentile, rank``.
        """
        st = sites[sites[window_col].astype(str).str[_FLANK].str.upper().isin(["S", "T"])]
        st = st.drop_duplicates(id_col)
        pct = self.percentile(st[window_col])
        rnk = self.rank(pct)
        site_idx, kin_idx = np.nonzero(rnk <= top_n)
        return pd.DataFrame({
            id_col: st[id_col].to_numpy()[site_idx],
            "kinase": np.asarray(self.kinases)[kin_idx],
            "percentile": pct[site_idx, kin_idx],
            "rank": rnk[site_idx, kin_idx],
        }).sort_values([id_col, "rank"]).reset_index(drop=True)
