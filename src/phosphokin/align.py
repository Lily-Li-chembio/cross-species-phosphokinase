"""Global protein alignment (Needleman-Wunsch with affine gaps, Gotoh).

A small, dependency-free implementation used to transfer phosphosite positions
between orthologous proteins. Each row of the dynamic-programming matrix is
computed with NumPy, so a typical pair of ~800-residue proteins aligns in tens
of milliseconds.

Gap cost for a gap of length L is ``gap_open + (L - 1) * gap_extend``.
"""

from __future__ import annotations

import numpy as np

# BLOSUM62 (NCBI). Order of rows/columns:
_B62_ORDER = "ARNDCQEGHILKMFPSTWYVBZX*"
_B62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""

BLOSUM62 = np.array([[int(v) for v in line.split()] for line in _B62_ROWS.strip().splitlines()],
                    dtype=np.int32)

# Map every byte to a BLOSUM62 index; unknown letters (U, O, J, lowercase ...) -> X.
_LUT = np.full(256, _B62_ORDER.index("X"), dtype=np.int64)
for _i, _aa in enumerate(_B62_ORDER):
    _LUT[ord(_aa)] = _i
    _LUT[ord(_aa.lower())] = _i

_NEG = np.int64(-(10 ** 9))


def _encode(seq: str) -> np.ndarray:
    return _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def align_position_map(seq_a: str, seq_b: str, gap_open: int = 11,
                       gap_extend: int = 1) -> tuple[np.ndarray, int]:
    """Globally align two protein sequences and return a residue position map.

    Parameters
    ----------
    seq_a, seq_b
        Protein sequences (one-letter code).
    gap_open, gap_extend
        Positive gap penalties; a gap of length L costs ``gap_open + (L-1)*gap_extend``.

    Returns
    -------
    pos_map
        Integer array of length ``len(seq_a) + 1``. ``pos_map[p]`` is the 1-based
        position in ``seq_b`` aligned to 1-based position ``p`` of ``seq_a``, or 0
        if that residue is aligned to a gap. ``pos_map[0]`` is unused (0).
    score
        The optimal alignment score.
    """
    if gap_extend > gap_open:
        raise ValueError("gap_extend must not exceed gap_open")
    a, b = _encode(seq_a), _encode(seq_b)
    n, m = len(a), len(b)
    pos_map = np.zeros(n + 1, dtype=np.int64)
    if n == 0 or m == 0:
        return pos_map, -(gap_open + (max(n, m) - 1) * gap_extend) if max(n, m) else 0

    go, ge = np.int64(gap_open), np.int64(gap_extend)
    j_idx = np.arange(m + 1, dtype=np.int64)

    # Row 0 (only b consumed -> horizontal gap)
    H_prev = np.empty(m + 1, dtype=np.int64)
    H_prev[0] = 0
    H_prev[1:] = -(go + (j_idx[1:] - 1) * ge)
    X_prev = np.full(m + 1, _NEG, dtype=np.int64)

    # Traceback tables
    tb_h = np.zeros((n + 1, m + 1), dtype=np.int8)     # state of H: 0=M, 1=X, 2=Y
    tb_h0 = np.zeros((n + 1, m + 1), dtype=np.int8)    # state of max(M, X): 0=M, 1=X
    tb_x = np.zeros((n + 1, m + 1), dtype=np.int8)     # X: 0=open from H above, 1=extend
    y_from = np.zeros((n + 1, m + 1), dtype=np.int32)  # Y: column where the gap opened

    sub = BLOSUM62[a][:, b].astype(np.int64)           # n x m substitution scores

    for i in range(1, n + 1):
        M = np.full(m + 1, _NEG, dtype=np.int64)
        M[1:] = H_prev[:-1] + sub[i - 1]

        x_open = H_prev - go
        x_ext = X_prev - ge
        X = np.maximum(x_open, x_ext)
        tb_x[i] = (x_ext > x_open)
        X[0] = -(go + (i - 1) * ge)
        tb_x[i, 0] = 1 if i > 1 else 0

        H0 = np.maximum(M, X)
        tb_h0[i] = (X > M)

        # Horizontal gaps: Y[j] = max_{k<j} H0[k] - go - ge*(j-1-k)
        T = H0 + ge * j_idx
        run = np.maximum.accumulate(T)
        arg = np.maximum.accumulate(np.where(T == run, j_idx, 0))
        Y = np.full(m + 1, _NEG, dtype=np.int64)
        Y[1:] = run[:-1] - go - ge * (j_idx[1:] - 1)
        y_from[i, 1:] = arg[:-1]

        H = np.maximum(H0, Y)
        tb_h[i] = np.where(Y > H0, 2, tb_h0[i])
        H_prev, X_prev = H, X

    score = int(H_prev[m])

    # Traceback
    i, j, state = n, m, int(tb_h[n, m])
    while i > 0 and j > 0:
        if state == 0:          # match/mismatch
            pos_map[i] = j
            i, j = i - 1, j - 1
            state = int(tb_h[i, j])
        elif state == 1:        # residue of a against a gap
            ext = tb_x[i, j]
            i -= 1
            state = 1 if ext else int(tb_h[i, j])
        else:                   # residues of b against gaps
            j = int(y_from[i, j])
            state = int(tb_h0[i, j])
    return pos_map, score


def aligned_identity(seq_a: str, seq_b: str, pos_map: np.ndarray) -> float:
    """Fraction of residues of ``seq_a`` aligned to an identical residue in ``seq_b``."""
    same = sum(1 for p in range(1, len(seq_a) + 1)
               if pos_map[p] and seq_a[p - 1] == seq_b[pos_map[p] - 1])
    return same / max(len(seq_a), 1)
