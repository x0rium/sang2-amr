"""6-signal resonance scan for genomic pair analysis.

Universal pair-scoring: R = geometric_mean(A, D, MI, P, X, G).
Works on aligned matrices (rows = sequences, columns = positions).

Validated in exploration/004:
  - Nucleotide level (|Σ|=4): weak P, R up to 0.25
  - Codon level (|Σ|=64): strong P, R up to 0.53 (ATG = magic bytes)
  - Recommendation: use codon level for CDS regions

Signals:
  A — frequency (log-scaled)
  D — directionality (asymmetry of forward vs reverse pair)
  MI — PPMI (positive pointwise mutual information)
  P — positional stability (column concentration)
  X — context stability (predictability of neighbors)
  G — gap consistency (what fills between b1 and b2)
"""

from __future__ import annotations

import numpy as np


def _norm01(x: float) -> float:
    if x != x:  # NaN
        return 0.0
    return max(0.0, min(1.0, x))


def _stability_from_hist(hist: np.ndarray) -> float:
    """1 - normalized entropy of a frequency histogram."""
    total = int(hist.sum())
    if total <= 1:
        return 0.0
    p = hist[hist > 0].astype(np.float64) / total
    h = float(-(p * np.log2(p)).sum())
    unique = int((hist > 0).sum())
    if unique < 2:
        return 1.0
    h_max = float(np.log2(unique))
    return 1.0 - (h / h_max)


def cds_to_codons(mat: np.ndarray) -> np.ndarray:
    """Convert nucleotide matrix to codon matrix (each 3 nt → index 0..63)."""
    n_rows, n_cols = mat.shape
    n_codons = n_cols // 3
    trimmed = mat[:, :n_codons * 3].reshape(n_rows, n_codons, 3).astype(np.int64)
    return (trimmed[:, :, 0] * 16 + trimmed[:, :, 1] * 4 + trimmed[:, :, 2]).astype(np.uint8)


def scan_and_score(
    mat: np.ndarray,
    alphabet_size: int = 4,
    k_range: tuple[int, int] = (1, 12),
    min_count: int = 10,
) -> list[tuple[int, int, int, dict[str, float], float, int]]:
    """Scan all (b1, b2, k) pairs and compute 6 signals + resonance R.

    Args:
        mat: aligned matrix (n_rows × n_cols), uint8
        alphabet_size: 4 for nucleotides, 64 for codons
        k_range: (min_k, max_k) pair offset range
        min_count: minimum pair count to include

    Returns list of (b1, b2, k, signals_dict, R, count) sorted by R desc.
    """
    n_rows, n_cols = mat.shape
    mat_i = mat.astype(np.int64)

    byte_hist = np.bincount(mat.ravel(), minlength=alphabet_size).astype(np.float64)
    byte_total = byte_hist.sum()
    p_byte = byte_hist / max(byte_total, 1)

    results = []

    for k in range(k_range[0], k_range[1] + 1):
        if k >= n_cols:
            continue
        n_pairs_per_row = n_cols - k
        total_positions = n_rows * n_pairs_per_row

        b1_flat = mat_i[:, :n_pairs_per_row].ravel()
        b2_flat = mat_i[:, k:k + n_pairs_per_row].ravel()
        pair_key = b1_flat * alphabet_size + b2_flat
        col_idx = np.tile(np.arange(n_pairs_per_row, dtype=np.int64), n_rows)
        row_idx = np.repeat(np.arange(n_rows, dtype=np.int64), n_pairs_per_row)

        global_counts = np.bincount(pair_key, minlength=alphabet_size**2)

        for pk in range(alphabet_size**2):
            count = int(global_counts[pk])
            if count < min_count:
                continue
            b1_val = pk // alphabet_size
            b2_val = pk % alphabet_size

            mask = pair_key == pk
            if mask.sum() == 0:
                continue

            sel_cols = col_idx[mask]
            sel_rows = row_idx[mask]

            # A — log-scaled frequency
            A = _norm01(np.log1p(count * 10) / np.log1p(total_positions))

            # D — directionality
            rev_pk = b2_val * alphabet_size + b1_val
            rev_count = int(global_counts[rev_pk])
            denom = count + rev_count
            D = _norm01(abs(count - rev_count) / denom if denom > 0 else 0.0)

            # MI — PPMI
            p_ab = count / max(total_positions, 1)
            p_a = float(p_byte[b1_val])
            p_b = float(p_byte[b2_val])
            if p_ab > 0 and p_a > 0 and p_b > 0:
                pmi = float(np.log2(p_ab / (p_a * p_b)))
            else:
                pmi = 0.0
            MI = _norm01(1.0 - float(np.exp(-max(pmi, 0.0) / 4.0)))

            # P — positional stability
            col_hist = np.bincount(sel_cols, minlength=n_pairs_per_row)
            P = _stability_from_hist(col_hist)

            # X — context stability
            left_mask = sel_cols > 0
            right_mask = sel_cols + k + 1 < n_cols
            left_hist = np.zeros(alphabet_size, dtype=np.int64)
            if left_mask.sum() > 0:
                left_hist = np.bincount(
                    mat_i[sel_rows[left_mask], sel_cols[left_mask] - 1],
                    minlength=alphabet_size,
                )
            right_hist = np.zeros(alphabet_size, dtype=np.int64)
            if right_mask.sum() > 0:
                rb = mat_i[sel_rows[right_mask], sel_cols[right_mask] + k]
                if rb.max() < alphabet_size:
                    right_hist = np.bincount(rb, minlength=alphabet_size)
            X = 0.5 * (_stability_from_hist(left_hist) + _stability_from_hist(right_hist))

            # G — gap consistency
            if k > 1:
                gap_hist = np.zeros(alphabet_size, dtype=np.int64)
                for g in range(1, min(k, 4)):
                    gap_pos = sel_cols + g
                    valid = gap_pos < n_cols
                    if valid.sum() > 0:
                        gap_hist += np.bincount(
                            mat_i[sel_rows[valid], gap_pos[valid]],
                            minlength=alphabet_size,
                        )
                G = _stability_from_hist(gap_hist)
            else:
                G = 1.0

            signals = {"A": A, "D": D, "MI": MI, "P": P, "X": X, "G": G}

            vals = np.array([A, D, MI, P, X, G])
            R = float(np.exp(np.log(vals).mean())) if (vals > 0).all() else 0.0

            results.append((b1_val, b2_val, k, signals, R, count))

        del b1_flat, b2_flat, pair_key, col_idx, row_idx, global_counts

    results.sort(key=lambda x: -x[4])
    return results
