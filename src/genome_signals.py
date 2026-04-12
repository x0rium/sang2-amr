"""Genomic-specific signals for SANG2-AMR.

Functions that embed DNA domain knowledge: strand symmetry, codon tables,
replication geometry. Not universal — specific to nucleotide sequences.

Validated in exploration/002b (GC-skew, Chargaff-2, ρ, RSCU, IR density).
"""

from __future__ import annotations

import numpy as np
from collections import Counter

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ALPHABET_SIZE = 4
NT_MAP = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3,
          ord("a"): 0, ord("c"): 1, ord("g"): 2, ord("t"): 3}
NT_CHARS = "ACGT"
COMPLEMENT = np.array([3, 2, 1, 0], dtype=np.uint8)  # A↔T, C↔G


def seq_to_array(seq: str) -> np.ndarray:
    """Convert DNA string to uint8 array (A=0, C=1, G=2, T=3)."""
    arr = np.array([NT_MAP.get(ord(c), 255) for c in seq], dtype=np.uint8)
    arr[arr >= 4] = 0  # replace ambiguous bases
    return arr


def reverse_complement(data: np.ndarray) -> np.ndarray:
    """Reverse complement of nucleotide array."""
    return COMPLEMENT[data[::-1]]


# ---------------------------------------------------------------------------
# GC-skew (replication geometry)
# ---------------------------------------------------------------------------
def gc_skew(data: np.ndarray, window: int = 1000) -> tuple[np.ndarray, np.ndarray]:
    """Sliding GC-skew: (G-C)/(G+C) per window.

    Validated in 002b: finds E. coli oriC with 0.1% error (3 kb on 4.6 Mb).
    """
    n = len(data)
    n_windows = n // window
    positions = np.arange(n_windows) * window + window // 2
    skew = np.empty(n_windows)
    for i in range(n_windows):
        seg = data[i * window:(i + 1) * window]
        g = int(np.sum(seg == 2))
        c = int(np.sum(seg == 1))
        denom = g + c
        skew[i] = (g - c) / denom if denom > 0 else 0.0
    return positions, skew


def gc_skew_cumulative(data: np.ndarray, window: int = 1000) -> tuple[np.ndarray, np.ndarray]:
    """Cumulative GC-skew. Minimum ≈ origin, maximum ≈ terminus."""
    positions, skew_vals = gc_skew(data, window)
    return positions, np.cumsum(skew_vals)


# ---------------------------------------------------------------------------
# Chargaff-2 deviation (windowed)
# ---------------------------------------------------------------------------
def chargaff2_deviation(data: np.ndarray, window: int = 5000,
                        step: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Per-window Chargaff's 2nd parity rule deviation.

    Returns (positions, deviations). High deviation = ssDNA origin, phage,
    or foreign functional block.
    Validated in 002b: conjugal region vs rest, p=0.001.
    """
    if step is None:
        step = window // 2
    n = len(data)
    n_windows = max(1, (n - window) // step + 1)
    positions = np.empty(n_windows)
    deviations = np.empty(n_windows)
    for i in range(n_windows):
        start = i * step
        seg = data[start:start + window]
        counts = np.bincount(seg, minlength=4).astype(np.float64)
        total = counts.sum()
        if total == 0:
            deviations[i] = 0.0
        else:
            freq = counts / total
            deviations[i] = abs(freq[0] - freq[3]) + abs(freq[1] - freq[2])
        positions[i] = start + window / 2
    return positions, deviations


# ---------------------------------------------------------------------------
# Inverted repeat density (IS-element detection)
# ---------------------------------------------------------------------------
def inverted_repeat_density(data: np.ndarray, min_len: int = 8,
                            max_gap: int = 500, window: int = 5000,
                            step: int = 1000) -> tuple[np.ndarray, np.ndarray]:
    """Count inverted repeat k-mer seeds per window.

    An IR: seq[i:i+k] == revcomp(seq[j:j+k]) within gap distance.
    Validated in 002b: higher near transposases (636 vs 552, trend correct).
    """
    n = len(data)
    k = min_len
    n_windows = max(1, (n - window) // step + 1)
    positions = np.empty(n_windows)
    ir_counts = np.empty(n_windows)

    # Precompute forward k-mer hashes
    fwd_keys = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        fwd_keys += data[i:n - k + 1 + i].astype(np.int64) * (4 ** (k - 1 - i))

    # Precompute reverse-complement k-mer hashes
    rc_hashes = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        rc_hashes += COMPLEMENT[data[k - 1 - i:n - i]].astype(np.int64) * (4 ** (k - 1 - i))

    for wi in range(n_windows):
        w_start = wi * step
        w_end = min(w_start + window, n - k + 1)
        fwd_set = set(fwd_keys[w_start:w_end].tolist())
        search_end = min(w_end + max_gap, n - k + 1)
        rc_set = set(rc_hashes[w_start:search_end].tolist())
        ir_counts[wi] = len(fwd_set & rc_set)
        positions[wi] = w_start + window / 2

    return positions, ir_counts


# ---------------------------------------------------------------------------
# RSCU (Relative Synonymous Codon Usage)
# ---------------------------------------------------------------------------
_SYNONYMOUS_GROUPS: dict[str, list[int]] | None = None


def _get_synonymous_groups() -> dict[str, list[int]]:
    global _SYNONYMOUS_GROUPS
    if _SYNONYMOUS_GROUPS is None:
        from Bio.Data.CodonTable import standard_dna_table
        groups: dict[str, list[int]] = {}
        for codon, aa in standard_dna_table.forward_table.items():
            idx = NT_MAP[ord(codon[0])] * 16 + NT_MAP[ord(codon[1])] * 4 + NT_MAP[ord(codon[2])]
            groups.setdefault(aa, []).append(idx)
        for codon in standard_dna_table.stop_codons:
            idx = NT_MAP[ord(codon[0])] * 16 + NT_MAP[ord(codon[1])] * 4 + NT_MAP[ord(codon[2])]
            groups.setdefault("*", []).append(idx)
        _SYNONYMOUS_GROUPS = groups
    return _SYNONYMOUS_GROUPS


def rscu(data: np.ndarray) -> np.ndarray:
    """Relative Synonymous Codon Usage — 64-element vector.

    RSCU(codon) = observed / expected_if_uniform_within_AA_group.
    """
    n = len(data)
    n_codons = n // 3
    if n_codons < 10:
        return np.ones(64)
    codons = data[:n_codons * 3].reshape(n_codons, 3).astype(np.int64)
    codon_idx = codons[:, 0] * 16 + codons[:, 1] * 4 + codons[:, 2]
    counts = np.bincount(codon_idx, minlength=64).astype(np.float64)

    rscu_vec = np.ones(64)
    for aa, group in _get_synonymous_groups().items():
        group_count = sum(counts[c] for c in group)
        n_syn = len(group)
        if group_count > 0:
            expected = group_count / n_syn
            for c in group:
                rscu_vec[c] = counts[c] / expected
    return rscu_vec


def codon_bias_distance(data: np.ndarray, ref_rscu: np.ndarray) -> float:
    """Euclidean distance between window RSCU and reference RSCU."""
    return float(np.sqrt(np.sum((rscu(data) - ref_rscu) ** 2)))


# ---------------------------------------------------------------------------
# Codon utilities
# ---------------------------------------------------------------------------
def to_codons(data: np.ndarray) -> np.ndarray:
    """Convert nucleotide array to codon index array (0..63)."""
    n_codons = len(data) // 3
    if n_codons == 0:
        return np.array([], dtype=np.int64)
    trimmed = data[:n_codons * 3].reshape(n_codons, 3).astype(np.int64)
    return trimmed[:, 0] * 16 + trimmed[:, 1] * 4 + trimmed[:, 2]


def codon_name(idx: int) -> str:
    """Convert codon index (0-63) to 3-letter string."""
    return NT_CHARS[idx // 16] + NT_CHARS[(idx // 4) % 4] + NT_CHARS[idx % 4]
