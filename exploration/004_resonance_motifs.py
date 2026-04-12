#!/usr/bin/env python3
"""Exploration 004: Resonance scan on E. coli K-12 CDS.

Resonance = geometric mean of 6 independent pair signals (A, D, MI, P, X, G).
Works on aligned matrix: rows = different sequences, columns = positions.

Two levels tested:
  A. Nucleotide (alphabet=4): rows = CDS aligned by start codon
  B. Codon (alphabet=64): rows = CDS converted to codons, aligned

Expected findings:
  - Top-R pairs at nucleotide level: positionally conserved dinucleotides
  - Eigendistances at k=3 (codon structure)
  - Top-R codon pairs: start/stop codons, conserved AA pairs

Success criteria: ≥3 of top-20 R-pairs are biologically interpretable.
"""

from __future__ import annotations

import numpy as np
from collections import Counter
from pathlib import Path
import sys

from Bio import SeqIO

DATA_DIR = Path(__file__).parent / "data"
NT_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
NT_CHARS = "ACGT"

# Codon table for interpretation
CODON_TO_AA = {}
def _init_codon_table():
    from Bio.Data.CodonTable import standard_dna_table
    for codon, aa in standard_dna_table.forward_table.items():
        idx = NT_MAP[codon[0]] * 16 + NT_MAP[codon[1]] * 4 + NT_MAP[codon[2]]
        CODON_TO_AA[idx] = aa
    for codon in standard_dna_table.stop_codons:
        idx = NT_MAP[codon[0]] * 16 + NT_MAP[codon[1]] * 4 + NT_MAP[codon[2]]
        CODON_TO_AA[idx] = '*'
_init_codon_table()

def codon_name(idx: int) -> str:
    c2 = idx % 4
    c1 = (idx // 4) % 4
    c0 = idx // 16
    return NT_CHARS[c0] + NT_CHARS[c1] + NT_CHARS[c2]

def codon_aa(idx: int) -> str:
    return CODON_TO_AA.get(idx, '?')


# =========================================================================
# Load data
# =========================================================================
def load_cds_matrix(gb_path: str, fasta_path: str, max_genes: int = 500,
                    target_len: int = 300) -> tuple[np.ndarray, list[str]]:
    """Build aligned CDS matrix: rows = genes, cols = nucleotide positions.

    Selects + strand genes of length ≥ target_len, truncates to target_len.
    Returns (matrix uint8 shape (N, target_len), gene_names).
    """
    genome_rec = SeqIO.read(fasta_path, "fasta")
    genome = str(genome_rec.seq).upper()

    gb_rec = SeqIO.read(gb_path, "genbank")
    rows = []
    names = []
    for feat in gb_rec.features:
        if feat.type != 'CDS':
            continue
        loc = feat.location
        if loc is None or loc.strand != 1:
            continue
        start, end = int(loc.start), int(loc.end)
        if end - start < target_len:
            continue
        seq = genome[start:start + target_len]
        row = np.array([NT_MAP.get(c, 0) for c in seq], dtype=np.uint8)
        rows.append(row)
        name = feat.qualifiers.get('gene', feat.qualifiers.get('product', ['?']))[0]
        names.append(name)
        if len(rows) >= max_genes:
            break

    mat = np.stack(rows)
    return mat, names


def cds_to_codons(mat: np.ndarray) -> np.ndarray:
    """Convert nucleotide matrix to codon matrix (each 3 nt → one codon index 0..63)."""
    n_rows, n_cols = mat.shape
    n_codons = n_cols // 3
    trimmed = mat[:, :n_codons * 3].reshape(n_rows, n_codons, 3).astype(np.int64)
    return (trimmed[:, :, 0] * 16 + trimmed[:, :, 1] * 4 + trimmed[:, :, 2]).astype(np.uint8)


# =========================================================================
# Resonance engine (adapted from auto-reverser/src/resonance.py)
# =========================================================================
def _norm01(x: float) -> float:
    if x != x:
        return 0.0
    return max(0.0, min(1.0, x))


def _stability_from_hist(hist: np.ndarray) -> float:
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


def scan_and_score(mat: np.ndarray, alphabet_size: int,
                   k_range: tuple[int, int] = (1, 12),
                   min_count: int = 10):
    """Scan all (b1, b2, k) pairs and compute 6 signals + resonance R.

    Returns list of (b1, b2, k, signals_dict, R, count) sorted by R desc.
    """
    n_rows, n_cols = mat.shape
    mat_i = mat.astype(np.int64)

    # Global byte histogram
    byte_hist = np.bincount(mat.ravel(), minlength=alphabet_size).astype(np.float64)
    byte_total = byte_hist.sum()
    p_byte = byte_hist / max(byte_total, 1)

    results = []

    for k in range(k_range[0], k_range[1] + 1):
        if k >= n_cols:
            continue
        n_pairs_per_row = n_cols - k

        # Flatten pairs
        b1 = mat_i[:, :n_pairs_per_row].ravel()
        b2 = mat_i[:, k:k + n_pairs_per_row].ravel()
        pair_key = b1 * alphabet_size + b2
        col_idx = np.tile(np.arange(n_pairs_per_row, dtype=np.int64), n_rows)

        # Global pair counts
        global_counts = np.bincount(pair_key, minlength=alphabet_size**2)

        # For each pair that passes min_count
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

            # A — log-scaled frequency
            total_positions = n_rows * n_pairs_per_row
            max_possible = total_positions  # theoretical max
            A = _norm01(np.log1p(count * 10) / np.log1p(max_possible))

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
            ppmi = max(pmi, 0.0)
            MI = _norm01(1.0 - float(np.exp(-ppmi / 4.0)))

            # P — positional stability (column concentration)
            col_hist = np.bincount(sel_cols, minlength=n_pairs_per_row)
            P = _stability_from_hist(col_hist)

            # X — context stability (left + right byte)
            left_mask = sel_cols > 0
            right_mask = sel_cols + k + 1 < n_cols

            # Flatten positions to get row indices
            sel_rows = np.repeat(np.arange(n_rows, dtype=np.int64), n_pairs_per_row)[mask]

            left_hist = np.zeros(alphabet_size, dtype=np.int64)
            if left_mask.sum() > 0:
                left_bytes = mat_i[sel_rows[left_mask], sel_cols[left_mask] - 1]
                left_hist = np.bincount(left_bytes, minlength=alphabet_size)

            right_hist = np.zeros(alphabet_size, dtype=np.int64)
            if right_mask.sum() > 0:
                right_bytes = mat_i[sel_rows[right_mask], sel_cols[right_mask] + k]
                if right_bytes.max() < alphabet_size:
                    right_hist = np.bincount(right_bytes, minlength=alphabet_size)

            X = 0.5 * (_stability_from_hist(left_hist) + _stability_from_hist(right_hist))

            # G — gap consistency (for k > 1)
            if k > 1:
                gap_hist = np.zeros(alphabet_size, dtype=np.int64)
                for g in range(1, min(k, 4)):  # limit gap scan for speed
                    gap_pos = sel_cols + g
                    valid = gap_pos < n_cols
                    if valid.sum() > 0:
                        gap_bytes = mat_i[sel_rows[valid], gap_pos[valid]]
                        gap_hist += np.bincount(gap_bytes, minlength=alphabet_size)
                G = _stability_from_hist(gap_hist)
            else:
                G = 1.0

            signals = {'A': A, 'D': D, 'MI': MI, 'P': P, 'X': X, 'G': G}

            # Geometric mean
            vals = np.array([A, D, MI, P, X, G])
            if (vals <= 0).any():
                R = 0.0
            else:
                R = float(np.exp(np.log(vals).mean()))

            results.append((b1_val, b2_val, k, signals, R, count))

    results.sort(key=lambda x: -x[4])
    return results


# =========================================================================
# Main
# =========================================================================
def main():
    print("=" * 70)
    print("EXPLORATION 004: Resonance Scan on E. coli K-12 CDS")
    print("=" * 70)

    # Load aligned CDS matrix
    print("\n[1] Building aligned CDS matrix...")
    mat_nt, gene_names = load_cds_matrix(
        str(DATA_DIR / "ecoli_k12.gb"),
        str(DATA_DIR / "ecoli_k12.fasta"),
        max_genes=500, target_len=300,
    )
    print(f"  Matrix: {mat_nt.shape[0]} genes × {mat_nt.shape[1]} nt positions")

    mat_codon = cds_to_codons(mat_nt)
    print(f"  Codon matrix: {mat_codon.shape[0]} genes × {mat_codon.shape[1]} codon positions")

    # ===================================================================
    # A. Nucleotide-level resonance
    # ===================================================================
    print("\n" + "=" * 70)
    print("[A] Nucleotide-Level Resonance (alphabet=4)")
    print("=" * 70)
    print("  k_range=(1, 15), min_count=50\n")

    nt_results = scan_and_score(mat_nt, alphabet_size=4, k_range=(1, 15), min_count=50)

    print(f"  Total pairs with R > 0: {sum(1 for r in nt_results if r[4] > 0)}")
    print(f"  Total pairs scanned: {len(nt_results)}")

    print(f"\n  Top 30 resonant nucleotide pairs:")
    print(f"  {'Pair':>6s}  {'k':>3s}  {'R':>7s}  {'A':>5s} {'D':>5s} {'MI':>5s} {'P':>5s} {'X':>5s} {'G':>5s}  {'Count':>6s}  Note")
    print("  " + "-" * 80)

    for b1, b2, k, sigs, R, count in nt_results[:30]:
        pair_str = f"{NT_CHARS[b1]}{NT_CHARS[b2]}"
        note = ""
        if k == 3:
            note = "codon-spacing"
        elif k == 6:
            note = "2-codon"
        elif k == 9:
            note = "3-codon"
        elif k == 1:
            note = "adjacent"
        elif k == 2:
            note = "intra-codon"
        print(f"  {pair_str:>6s}  {k:>3d}  {R:7.4f}  "
              f"{sigs['A']:5.3f} {sigs['D']:5.3f} {sigs['MI']:5.3f} "
              f"{sigs['P']:5.3f} {sigs['X']:5.3f} {sigs['G']:5.3f}  "
              f"{count:>6d}  {note}")

    # Eigendistance analysis
    print(f"\n  Eigendistance analysis (best k per pair):")
    eigen: dict[str, tuple[int, float]] = {}
    for b1, b2, k, sigs, R, count in nt_results:
        pair = f"{NT_CHARS[b1]}{NT_CHARS[b2]}"
        if pair not in eigen or R > eigen[pair][1]:
            eigen[pair] = (k, R)

    k_counts = Counter(k for k, r in eigen.values())
    print(f"  Distribution of eigendistances:")
    for k_val in sorted(k_counts.keys()):
        n = k_counts[k_val]
        bar = "#" * n
        note = ""
        if k_val == 3:
            note = " ← codon"
        elif k_val == 1:
            note = " ← adjacent"
        elif k_val % 3 == 0:
            note = f" ← {k_val//3}-codon"
        print(f"    k={k_val:2d}: {n:2d} pairs {bar}{note}")

    # ===================================================================
    # B. Codon-level resonance
    # ===================================================================
    print("\n" + "=" * 70)
    print("[B] Codon-Level Resonance (alphabet=64)")
    print("=" * 70)
    print("  k_range=(1, 8), min_count=20\n")

    codon_results = scan_and_score(mat_codon, alphabet_size=64, k_range=(1, 8), min_count=20)

    print(f"  Total pairs with R > 0: {sum(1 for r in codon_results if r[4] > 0)}")
    print(f"  Total pairs scanned: {len(codon_results)}")

    print(f"\n  Top 30 resonant codon pairs:")
    print(f"  {'Pair':>10s}  {'AA':>4s}  {'k':>3s}  {'R':>7s}  {'A':>5s} {'D':>5s} {'MI':>5s} {'P':>5s} {'X':>5s} {'G':>5s}  {'Count':>5s}")
    print("  " + "-" * 80)

    for b1, b2, k, sigs, R, count in codon_results[:30]:
        c1_name = codon_name(b1)
        c2_name = codon_name(b2)
        aa_pair = f"{codon_aa(b1)}{codon_aa(b2)}"
        print(f"  {c1_name}-{c2_name}  {aa_pair:>4s}  {k:>3d}  {R:7.4f}  "
              f"{sigs['A']:5.3f} {sigs['D']:5.3f} {sigs['MI']:5.3f} "
              f"{sigs['P']:5.3f} {sigs['X']:5.3f} {sigs['G']:5.3f}  "
              f"{count:>5d}")

    # Codon eigendistances
    print(f"\n  Codon eigendistance distribution:")
    codon_eigen: dict[str, tuple[int, float]] = {}
    for b1, b2, k, sigs, R, count in codon_results:
        if R <= 0:
            continue
        pair = f"{codon_name(b1)}-{codon_name(b2)}"
        if pair not in codon_eigen or R > codon_eigen[pair][1]:
            codon_eigen[pair] = (k, R)

    codon_k_counts = Counter(k for k, r in codon_eigen.values())
    for k_val in sorted(codon_k_counts.keys()):
        n = codon_k_counts[k_val]
        bar = "#" * min(n, 50)
        print(f"    k={k_val}: {n:4d} pairs {bar}")

    # Biological interpretation of top codon pairs
    print(f"\n  Biological patterns in top-50 codon pairs:")

    # Check for start codon involvement
    atg_idx = NT_MAP['A'] * 16 + NT_MAP['T'] * 4 + NT_MAP['G']
    atg_pairs = [(b1, b2, k, R) for b1, b2, k, _, R, _ in codon_results[:100]
                 if b1 == atg_idx or b2 == atg_idx]
    if atg_pairs:
        print(f"  ATG (start) in top-100:")
        for b1, b2, k, R in atg_pairs[:5]:
            print(f"    {codon_name(b1)}-{codon_name(b2)} k={k} R={R:.4f} ({codon_aa(b1)}{codon_aa(b2)})")

    # Check for stop codon involvement
    stop_idxs = [NT_MAP[c[0]]*16 + NT_MAP[c[1]]*4 + NT_MAP[c[2]] for c in ['TAA', 'TAG', 'TGA']]
    stop_pairs = [(b1, b2, k, R) for b1, b2, k, _, R, _ in codon_results[:100]
                  if b1 in stop_idxs or b2 in stop_idxs]
    if stop_pairs:
        print(f"  Stop codons in top-100:")
        for b1, b2, k, R in stop_pairs[:5]:
            print(f"    {codon_name(b1)}-{codon_name(b2)} k={k} R={R:.4f} ({codon_aa(b1)}{codon_aa(b2)})")

    # Same-AA codon pairs (synonymous codon coupling)
    same_aa = [(b1, b2, k, R) for b1, b2, k, _, R, _ in codon_results[:100]
               if codon_aa(b1) == codon_aa(b2) and b1 != b2]
    if same_aa:
        print(f"  Synonymous codon pairs (same AA, different codon) in top-100:")
        for b1, b2, k, R in same_aa[:5]:
            print(f"    {codon_name(b1)}-{codon_name(b2)} k={k} R={R:.4f} ({codon_aa(b1)}-{codon_aa(b2)})")

    # ===================================================================
    # C. Signal independence check
    # ===================================================================
    print("\n" + "=" * 70)
    print("[C] Signal Independence")
    print("=" * 70)

    # Correlation matrix of 6 signals at nucleotide level
    if len(nt_results) > 10:
        sig_names = ['A', 'D', 'MI', 'P', 'X', 'G']
        sig_data = np.array([[r[3][s] for s in sig_names] for r in nt_results])
        corr = np.corrcoef(sig_data.T)
        print(f"\n  Nucleotide-level signal correlation:")
        print(f"  {'':>4s}", end="")
        for s in sig_names:
            print(f" {s:>6s}", end="")
        print()
        for i, s1 in enumerate(sig_names):
            print(f"  {s1:>4s}", end="")
            for j, s2 in enumerate(sig_names):
                val = corr[i, j]
                marker = "*" if abs(val) > 0.3 and i != j else " "
                print(f" {val:5.2f}{marker}", end="")
            print()

    if len(codon_results) > 10:
        sig_data_c = np.array([[r[3][s] for s in sig_names] for r in codon_results])
        corr_c = np.corrcoef(sig_data_c.T)
        print(f"\n  Codon-level signal correlation:")
        print(f"  {'':>4s}", end="")
        for s in sig_names:
            print(f" {s:>6s}", end="")
        print()
        for i, s1 in enumerate(sig_names):
            print(f"  {s1:>4s}", end="")
            for j, s2 in enumerate(sig_names):
                val = corr_c[i, j]
                marker = "*" if abs(val) > 0.3 and i != j else " "
                print(f" {val:5.2f}{marker}", end="")
            print()

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    checks = []

    # Check 1: Eigendistance k=3 prominent at nucleotide level
    k3_fraction = k_counts.get(3, 0) / max(sum(k_counts.values()), 1)
    ok = k3_fraction > 0.1 or k_counts.get(3, 0) >= 2
    checks.append(("Eigendistance k=3 prominent", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] k=3 eigendistance: {k_counts.get(3,0)} pairs ({k3_fraction*100:.0f}%)")

    # Check 2: ≥3 of top-20 nucleotide pairs biologically interpretable
    # (k=3 or k=multiples of 3 → codon-related)
    top20_bio = sum(1 for _, _, k, _, R, _ in nt_results[:20] if k % 3 == 0 and R > 0)
    ok = top20_bio >= 3
    checks.append(("≥3 top-20 bio-interpretable", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Top-20 nt pairs with codon-related k: {top20_bio}")

    # Check 3: Codon resonance produces diverse results
    ok = len(codon_results) > 50
    checks.append(("Codon resonance produces >50 pairs", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Codon pairs with R>0: {len([r for r in codon_results if r[4] > 0])}")

    # Check 4: ATG appears in top codon pairs
    ok = len(atg_pairs) > 0
    checks.append(("ATG in top-100 codon pairs", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] ATG pairs in top-100: {len(atg_pairs)}")

    # Check 5: Signals are reasonably independent
    if len(nt_results) > 10:
        max_offdiag = max(abs(corr[i, j]) for i in range(6) for j in range(6) if i != j)
        ok = max_offdiag < 0.7
        checks.append(("Signal independence (max |ρ| < 0.7)", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] Max off-diagonal |ρ|: {max_offdiag:.2f}")

    passed = sum(1 for _, ok in checks if ok)
    print(f"\n  Result: {passed}/{len(checks)} checks passed")
    print(f"  {'SUCCESS — proceed to 005' if passed >= 3 else 'NEEDS INVESTIGATION'}")


if __name__ == "__main__":
    main()
