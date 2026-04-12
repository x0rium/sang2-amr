#!/usr/bin/env python3
"""Exploration 001: SANG2 signals on E. coli K-12 (U00096.3).

Goal: Verify that SANG2 signal algorithms, adapted from byte (256-symbol)
to nucleotide (4-symbol) alphabet, produce biologically meaningful results.

Key hypothesis: MI half-life ≈ 3 for coding regions (codon structure).

Success criteria (from exploration/README.md):
  - MI half-life for CDS in range 2.5–4.0
  - Autocorrelation peak at lag=3 in coding regions
  - H_ratio in CODING range (0.60–0.85) for CDS
  - Hurst ≈ 0.6–0.7 (persistent)
  - ΔH negative at gene starts (convergent = structural)
"""

from __future__ import annotations

import sys
import os
import numpy as np
from collections import Counter
from pathlib import Path
from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
from Bio import SeqIO

DATA_DIR = Path(__file__).parent / "data"
FASTA_PATH = DATA_DIR / "ecoli_k12.fasta"
GB_PATH = DATA_DIR / "ecoli_k12.gb"

NT_MAP = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3,
          ord('a'): 0, ord('c'): 1, ord('g'): 2, ord('t'): 3}

ALPHABET_SIZE = 4
H_MAX = np.log2(ALPHABET_SIZE)  # 2.0 bits


def load_genome() -> np.ndarray:
    """Load E. coli K-12 as uint8 array (A=0, C=1, G=2, T=3)."""
    rec = SeqIO.read(str(FASTA_PATH), "fasta")
    seq = str(rec.seq).upper()
    arr = np.array([NT_MAP.get(ord(c), 255) for c in seq], dtype=np.uint8)
    valid = arr < 4
    if not valid.all():
        n_bad = int((~valid).sum())
        print(f"  WARNING: {n_bad} non-ACGT bases, replacing with random nt")
        arr[~valid] = np.random.randint(0, 4, size=n_bad).astype(np.uint8)
    return arr


@dataclass
class GeneRegion:
    start: int
    end: int
    strand: int  # 1 or -1
    kind: str    # 'CDS', 'rRNA', 'tRNA', etc.
    name: str


def load_annotations() -> list[GeneRegion]:
    """Load gene annotations from GenBank."""
    rec = SeqIO.read(str(GB_PATH), "genbank")
    regions = []
    for feat in rec.features:
        if feat.type in ('CDS', 'rRNA', 'tRNA'):
            loc = feat.location
            if loc is None:
                continue
            name = feat.qualifiers.get('gene', feat.qualifiers.get('product', ['?']))[0]
            regions.append(GeneRegion(
                start=int(loc.start),
                end=int(loc.end),
                strand=int(loc.strand) if loc.strand else 1,
                kind=feat.type,
                name=name,
            ))
    return regions


# ---------------------------------------------------------------------------
# Signal functions adapted for 4-letter alphabet
# ---------------------------------------------------------------------------

def entropy(counts) -> float:
    """Shannon entropy H = -sum P(x) log2 P(x) in bits."""
    if isinstance(counts, (Counter, dict)):
        values = np.fromiter(counts.values(), dtype=np.float64)
    elif isinstance(counts, np.ndarray):
        values = counts.astype(np.float64, copy=False)
    else:
        values = np.fromiter(counts, dtype=np.float64)
    total = values.sum()
    if total <= 1:
        return 0.0
    probs = values / total
    probs = probs[probs > 0]
    return float(-np.sum(probs * np.log2(probs)))


def sdr_kmer(sequence: np.ndarray, k: int = 4):
    """SDR on k-mer frequencies (alphabet too small for single-nt SDR).

    Returns (structural_kmers, theta, sdr_ratio).
    sdr_ratio = |structural| / |observed kmers|.
    """
    n = len(sequence)
    if n < k:
        return set(), 0.0, 0.0

    # Build k-mer keys as base-4 integers
    keys = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        keys += sequence[i:n - k + 1 + i].astype(np.int64) * (4 ** (k - 1 - i))

    freq = Counter(int(x) for x in keys)
    total = sum(freq.values())
    if total == 0:
        return set(), 0.0, 0.0

    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = (lo + hi) / 2
        mass = sum(c for c in freq.values() if c / total > mid)
        if mass / total > 0.50:
            lo = mid
        else:
            hi = mid
    theta = (lo + hi) / 2
    structural = {x for x, c in freq.items() if c / total > theta}
    sdr_ratio = len(structural) / max(len(freq), 1)
    return structural, theta, sdr_ratio


def otsu_threshold(values) -> tuple[float, float]:
    """Otsu's method: optimal binary split."""
    vals = sorted(values)
    n = len(vals)
    if n < 2:
        return 0.0, 0.0
    total_sum = float(sum(vals))
    best_t = vals[0]
    best_var = 0.0
    left_sum = 0.0
    for i in range(1, n):
        left_sum += vals[i - 1]
        right_sum = total_sum - left_sum
        w0 = i / n
        w1 = (n - i) / n
        m0 = left_sum / i
        m1 = right_sum / (n - i)
        between = w0 * w1 * (m0 - m1) ** 2
        if between > best_var:
            best_var = between
            best_t = (vals[i - 1] + vals[i]) / 2
    return float(best_t), float(best_var)


def mi_at_lag(data: np.ndarray, lag: int) -> float:
    """Mutual information I(X_i; X_{i+lag}) for 4-letter alphabet."""
    n = len(data)
    if n <= lag + 1:
        return 0.0
    marg_hist = np.bincount(data, minlength=ALPHABET_SIZE)
    h_marg = entropy(marg_hist)

    a = data[:n - lag].astype(np.int64)
    b = data[lag:].astype(np.int64)
    joint_key = a * ALPHABET_SIZE + b  # 0..15
    joint_hist = np.bincount(joint_key, minlength=ALPHABET_SIZE**2)
    h_joint = entropy(joint_hist)
    return max(2 * h_marg - h_joint, 0.0)


def mi_curve(data: np.ndarray, max_lag: int = 64) -> tuple[np.ndarray, np.ndarray]:
    """MI(lag) for lag = 1..max_lag."""
    n = len(data)
    max_lag = max(1, min(max_lag, n // 4))
    lags = np.arange(1, max_lag + 1)
    mi = np.array([mi_at_lag(data, int(k)) for k in lags], dtype=np.float64)
    return lags, mi


def mi_halflife(data: np.ndarray, max_lag: int = 64) -> float:
    """MI half-life: lag at which MI decays to 1/2 of start.

    Log-linear regression: log(MI) ~ lag; half-life = -ln(2)/slope.
    """
    lags, mi = mi_curve(data, max_lag=max_lag)
    if len(mi) < 3 or mi[0] <= 0:
        return 0.0
    mi_safe = np.clip(mi, 1e-10, None)
    log_mi = np.log(mi_safe)
    A = np.vstack([lags.astype(np.float64), np.ones(len(lags))]).T
    slope, _ = np.linalg.lstsq(A, log_mi, rcond=None)[0]
    if slope >= -1e-6:
        return float("inf")
    return float(-np.log(2) / slope)


def delta_h(data: np.ndarray, k: int = 3) -> np.ndarray:
    """Bidirectional entropy ΔH for 4-letter alphabet.

    ΔH[i] = H(next | past_k) - H(prev | future_k)
    Negative = convergent (structural), Positive = divergent (content).
    """
    n = len(data)
    if n < 2 * k + 2:
        return np.zeros(n, dtype=np.float64)

    arr = data.astype(np.int64)
    dh = np.zeros(n, dtype=np.float64)

    # Hash contexts as base-4 integers
    def hash_ctx(start, length):
        out = np.zeros(length, dtype=np.int64)
        for s in range(k):
            out = out * ALPHABET_SIZE + arr[start + s:start + s + length]
        return out

    def context_h(contexts, targets):
        n_pos = contexts.size
        if n_pos == 0:
            return np.zeros(0, dtype=np.float64)
        uniq, inverse = np.unique(contexts, return_inverse=True)
        n_ctx = uniq.size
        flat = inverse.astype(np.int64) * ALPHABET_SIZE + targets
        joint = np.bincount(flat, minlength=n_ctx * ALPHABET_SIZE).reshape(n_ctx, ALPHABET_SIZE).astype(np.float64)
        row_sums = joint.sum(axis=1, keepdims=True)
        row_sums_safe = np.maximum(row_sums, 1.0)
        p = joint / row_sums_safe
        with np.errstate(divide="ignore", invalid="ignore"):
            log_p = np.where(p > 0, np.log2(p), 0.0)
        h_per_ctx = -(p * log_p).sum(axis=1)
        return h_per_ctx[inverse]

    # Forward: H(next | past_k)
    fwd_ctx = hash_ctx(0, n - k)
    fwd_targets = arr[k:]
    h_fwd = context_h(fwd_ctx, fwd_targets)

    # Backward: H(prev | future_k)
    bwd_ctx = hash_ctx(1, n - k)
    bwd_targets = arr[:n - k]
    h_bwd = context_h(bwd_ctx, bwd_targets)

    lo, hi = k, n - k
    dh[lo:hi] = h_fwd[:hi - lo] - h_bwd[lo:hi]
    return dh


def autocorrelation(data: np.ndarray, max_lag: int = 128):
    """Autocorrelation + peak detection."""
    arr = data.astype(np.float64)
    n = len(arr)
    max_lag = max(2, min(max_lag, n // 4))
    centered = arr - arr.mean()
    norm = float(centered @ centered)
    lags = np.arange(max_lag)

    if norm == 0.0:
        return lags, np.zeros(max_lag), []

    ac = np.empty(max_lag, dtype=np.float64)
    for lag in range(max_lag):
        ac[lag] = (centered[:n - lag] @ centered[lag:]) / norm

    peaks = []
    for i in range(2, max_lag - 1):
        if ac[i] > ac[i - 1] and ac[i] > ac[i + 1] and ac[i] > 0.005:
            peaks.append((int(i), float(ac[i])))
    return lags, ac, peaks


def h_ratio(data: np.ndarray) -> float:
    """H / H_max for 4-letter alphabet."""
    freq = np.bincount(data, minlength=ALPHABET_SIZE)
    h = entropy(freq)
    return float(h / H_MAX)


def regime_classify(data: np.ndarray) -> str:
    """Classify regime for genomic data (adapted thresholds for H_max=2.0)."""
    r = h_ratio(data)
    if r < 0.25:
        return "SPARSE"        # homopolymer tracts
    if r < 0.60:
        return "BIASED"        # extreme GC bias
    if r < 0.85:
        return "CODING"        # coding regions with codon bias
    return "BALANCED"          # intergenic, high variability


def hurst(data: np.ndarray) -> float:
    """Hurst exponent via R/S analysis."""
    arr = data.astype(np.float64)
    n = len(arr)
    if n < 64:
        return 0.5

    windows = []
    rs_values = []
    w = 8
    while w < n // 4:
        rs_list = []
        for start in range(0, n - w, w):
            seg = arr[start:start + w]
            mean = seg.mean()
            cumdev = np.cumsum(seg - mean)
            R = float(cumdev.max() - cumdev.min())
            S = float(seg.std(ddof=1))
            if S > 0:
                rs_list.append(R / S)
        if rs_list:
            windows.append(w)
            rs_values.append(float(np.mean(rs_list)))
        w = int(w * 1.5)

    if len(windows) < 3:
        return 0.5
    log_w = np.log(windows)
    log_rs = np.log(rs_values)
    A_mat = np.vstack([log_w, np.ones(len(log_w))]).T
    H, _ = np.linalg.lstsq(A_mat, log_rs, rcond=None)[0]
    return float(np.clip(H, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("EXPLORATION 001: SANG2 Signals on E. coli K-12 (U00096.3)")
    print("=" * 70)

    # Load data
    print("\n[1] Loading genome...")
    genome = load_genome()
    print(f"  Genome length: {len(genome):,} nt")
    print(f"  GC content: {np.mean((genome == 1) | (genome == 2)):.4f}")

    print("\n[2] Loading annotations...")
    annotations = load_annotations()
    cds_regions = [r for r in annotations if r.kind == 'CDS']
    print(f"  Total annotated features: {len(annotations)}")
    print(f"  CDS: {len(cds_regions)}")

    # -----------------------------------------------------------------------
    # A1. Global Shannon Entropy
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("[A1] Shannon Entropy")
    print("=" * 70)
    h_global = entropy(np.bincount(genome, minlength=4))
    print(f"  H(genome)     = {h_global:.4f} bits")
    print(f"  H_max         = {H_MAX:.4f} bits")
    print(f"  H/H_max       = {h_global/H_MAX:.4f}")

    # Sliding window entropy
    window = 300
    step = 300
    n_windows = (len(genome) - window) // step
    h_values = np.empty(n_windows)
    for i in range(n_windows):
        start = i * step
        seg = genome[start:start + window]
        h_values[i] = entropy(np.bincount(seg, minlength=4))

    print(f"  Sliding window (w={window}, step={step}):")
    print(f"    mean H = {h_values.mean():.4f}, std = {h_values.std():.4f}")
    print(f"    min H  = {h_values.min():.4f}, max H = {h_values.max():.4f}")

    # -----------------------------------------------------------------------
    # A2. SDR on k-mers
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("[A2] SDR (Structural Detection Rule) on 4-mers")
    print("=" * 70)
    structural, theta, sdr_ratio = sdr_kmer(genome, k=4)
    print(f"  Theta          = {theta:.6f}")
    print(f"  Structural 4-mers: {len(structural)} / {4**4} possible")
    print(f"  SDR ratio      = {sdr_ratio:.4f}")

    # Decode top structural k-mers
    nt_chars = "ACGT"
    def decode_kmer(key, k=4):
        s = ""
        for _ in range(k):
            s = nt_chars[key % 4] + s
            key //= 4
        return s

    print(f"  Top structural 4-mers (carrying ~50% mass):")
    # Get frequencies for structural k-mers
    n = len(genome)
    keys = np.zeros(n - 3, dtype=np.int64)
    for i in range(4):
        keys += genome[i:n - 3 + i].astype(np.int64) * (4 ** (3 - i))
    freq = Counter(int(x) for x in keys)
    struct_sorted = sorted(structural, key=lambda x: freq.get(x, 0), reverse=True)
    for km in struct_sorted[:10]:
        print(f"    {decode_kmer(km)} : {freq[km]:>7,} ({freq[km]/(n-3)*100:.2f}%)")

    # -----------------------------------------------------------------------
    # A3. Otsu on GC-content per window
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("[A3] Otsu Threshold on GC-content")
    print("=" * 70)
    gc_values = []
    for i in range(n_windows):
        start = i * step
        seg = genome[start:start + window]
        gc = float(np.mean((seg == 1) | (seg == 2)))
        gc_values.append(gc)

    gc_thresh, gc_var = otsu_threshold(gc_values)
    below = sum(1 for g in gc_values if g < gc_thresh)
    above = len(gc_values) - below
    print(f"  Otsu threshold = {gc_thresh:.4f}")
    print(f"  Inter-class var = {gc_var:.6f}")
    print(f"  Windows below: {below}, above: {above}")
    print(f"  GC stats: mean={np.mean(gc_values):.4f}, std={np.std(gc_values):.4f}")

    # -----------------------------------------------------------------------
    # A5. MI half-life (KEY EXPERIMENT)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("[A5] MI Half-life (KEY HYPOTHESIS)")
    print("=" * 70)
    print("  Prediction: MI half-life ≈ 3 for coding DNA (codon triplets)")
    print()

    # Sample CDS regions
    np.random.seed(42)
    sample_size = min(200, len(cds_regions))
    sampled_cds = [cds_regions[i] for i in np.random.choice(len(cds_regions), sample_size, replace=False)]

    cds_halflifes = []
    for reg in sampled_cds:
        seg = genome[reg.start:reg.end]
        if len(seg) < 100:
            continue
        hl = mi_halflife(seg, max_lag=32)
        if 0 < hl < 100:  # filter degenerate
            cds_halflifes.append(hl)

    print(f"  CDS regions sampled: {sample_size}, valid: {len(cds_halflifes)}")
    if cds_halflifes:
        arr_hl = np.array(cds_halflifes)
        print(f"  MI half-life (CDS):")
        print(f"    mean   = {arr_hl.mean():.2f}")
        print(f"    median = {np.median(arr_hl):.2f}")
        print(f"    std    = {arr_hl.std():.2f}")
        print(f"    range  = [{arr_hl.min():.2f}, {arr_hl.max():.2f}]")
        in_range = np.sum((arr_hl >= 2.5) & (arr_hl <= 4.0))
        print(f"    in [2.5, 4.0]: {in_range}/{len(arr_hl)} ({in_range/len(arr_hl)*100:.1f}%)")

    # Intergenic comparison
    intergenic_halflifes = []
    sorted_cds = sorted(cds_regions, key=lambda r: r.start)
    for i in range(len(sorted_cds) - 1):
        gap_start = sorted_cds[i].end
        gap_end = sorted_cds[i + 1].start
        if gap_end - gap_start >= 100:
            seg = genome[gap_start:gap_end]
            hl = mi_halflife(seg, max_lag=32)
            if 0 < hl < 100:
                intergenic_halflifes.append(hl)
        if len(intergenic_halflifes) >= 200:
            break

    if intergenic_halflifes:
        arr_ig = np.array(intergenic_halflifes)
        print(f"\n  MI half-life (intergenic, n={len(intergenic_halflifes)}):")
        print(f"    mean   = {arr_ig.mean():.2f}")
        print(f"    median = {np.median(arr_ig):.2f}")
        print(f"    std    = {arr_ig.std():.2f}")

    # MI curve for one representative CDS
    print(f"\n  MI decay curve for a representative CDS (dnaA, first gene):")
    dnaa = [r for r in cds_regions if 'dnaA' in r.name.lower() or 'dnaa' in r.name.lower()]
    if dnaa:
        seg = genome[dnaa[0].start:dnaa[0].end]
        lags, mi_vals = mi_curve(seg, max_lag=20)
        for lag, mi_val in zip(lags, mi_vals):
            bar = "#" * int(mi_val * 200)
            print(f"    lag {lag:2d}: MI={mi_val:.4f} {bar}")

    # -----------------------------------------------------------------------
    # A6. Autocorrelation
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("[A6] Autocorrelation + Peak Detection")
    print("=" * 70)

    # CDS sample
    big_cds = [r for r in cds_regions if r.end - r.start >= 500]
    if big_cds:
        sample_ac = big_cds[:5]
        for reg in sample_ac:
            seg = genome[reg.start:reg.end]
            _, ac_vals, peaks = autocorrelation(seg, max_lag=30)
            peak_str = ", ".join(f"lag={p[0]}({p[1]:.4f})" for p in peaks[:5])
            print(f"  {reg.name[:15]:15s} [{reg.end-reg.start:5d} nt]: peaks = {peak_str}")

    # Check lag-3 prevalence across CDS
    lag3_count = 0
    total_checked = 0
    for reg in cds_regions:
        if reg.end - reg.start < 200:
            continue
        seg = genome[reg.start:reg.end]
        _, _, peaks = autocorrelation(seg, max_lag=20)
        total_checked += 1
        if any(p[0] == 3 for p in peaks):
            lag3_count += 1
        if total_checked >= 500:
            break

    print(f"\n  Lag-3 peak in CDS: {lag3_count}/{total_checked} ({lag3_count/max(total_checked,1)*100:.1f}%)")
    print(f"  (Expected: high % due to codon triplet structure)")

    # -----------------------------------------------------------------------
    # A4. Bidirectional entropy (ΔH) around gene starts
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("[A4] Bidirectional Entropy (ΔH) around gene starts")
    print("=" * 70)
    print("  Prediction: ΔH < 0 (convergent) at promoters/RBS, ΔH > 0 after gene end")

    # Average ΔH profile around gene starts
    flank = 100  # nt before/after start codon
    dh_profiles = []
    for reg in cds_regions[:300]:
        if reg.strand != 1:
            continue
        s = reg.start - flank
        e = reg.start + flank
        if s < 0 or e >= len(genome):
            continue
        seg = genome[s:e]
        dh = delta_h(seg, k=5)
        dh_profiles.append(dh)

    if dh_profiles:
        avg_dh = np.mean(dh_profiles, axis=0)
        print(f"  Averaged ΔH around {len(dh_profiles)} gene starts (+ strand):")
        print(f"  Position relative to start codon:")
        positions = [-50, -30, -20, -10, -5, 0, 5, 10, 20, 30, 50]
        for pos in positions:
            idx = pos + flank
            if 0 <= idx < len(avg_dh):
                val = avg_dh[idx]
                marker = "◀ convergent" if val < -0.01 else ("▶ divergent" if val > 0.01 else "  neutral")
                print(f"    {pos:+4d}: ΔH = {val:+.4f} {marker}")

    # -----------------------------------------------------------------------
    # A7. H/H_max Regime Classification
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("[A7] H/H_max Regime Classification")
    print("=" * 70)

    # Classify windows across genome
    regime_counts = Counter()
    for i in range(n_windows):
        start = i * step
        seg = genome[start:start + window]
        r = regime_classify(seg)
        regime_counts[r] += 1

    print(f"  Window classification (w={window}):")
    for r in ["SPARSE", "BIASED", "CODING", "BALANCED"]:
        c = regime_counts.get(r, 0)
        pct = c / n_windows * 100
        bar = "#" * int(pct)
        print(f"    {r:10s}: {c:5d} ({pct:5.1f}%) {bar}")

    # H_ratio for CDS vs intergenic
    cds_ratios = []
    for reg in sampled_cds:
        seg = genome[reg.start:reg.end]
        if len(seg) >= 50:
            cds_ratios.append(h_ratio(seg))

    ig_ratios = []
    for i in range(len(sorted_cds) - 1):
        gap_start = sorted_cds[i].end
        gap_end = sorted_cds[i + 1].start
        if gap_end - gap_start >= 50:
            seg = genome[gap_start:gap_end]
            ig_ratios.append(h_ratio(seg))
        if len(ig_ratios) >= 200:
            break

    if cds_ratios:
        print(f"\n  H_ratio (CDS, n={len(cds_ratios)}): mean={np.mean(cds_ratios):.4f}, std={np.std(cds_ratios):.4f}")
    if ig_ratios:
        print(f"  H_ratio (intergenic, n={len(ig_ratios)}): mean={np.mean(ig_ratios):.4f}, std={np.std(ig_ratios):.4f}")

    # -----------------------------------------------------------------------
    # A8. Hurst Exponent
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("[A8] Hurst Exponent")
    print("=" * 70)
    print("  Prediction: H ≈ 0.58–0.67 (persistent, Peng et al. 1992)")

    # Global Hurst on first 100k
    h_global_hurst = hurst(genome[:100_000])
    print(f"  Hurst (first 100kb): {h_global_hurst:.4f}")

    # CDS vs intergenic Hurst
    cds_hursts = []
    for reg in cds_regions:
        if reg.end - reg.start >= 500:
            h_val = hurst(genome[reg.start:reg.end])
            cds_hursts.append(h_val)
        if len(cds_hursts) >= 100:
            break

    ig_hursts = []
    for i in range(len(sorted_cds) - 1):
        gap_start = sorted_cds[i].end
        gap_end = sorted_cds[i + 1].start
        if gap_end - gap_start >= 500:
            ig_hursts.append(hurst(genome[gap_start:gap_end]))
        if len(ig_hursts) >= 50:
            break

    if cds_hursts:
        print(f"  Hurst (CDS, n={len(cds_hursts)}): mean={np.mean(cds_hursts):.4f}, std={np.std(cds_hursts):.4f}")
    if ig_hursts:
        print(f"  Hurst (intergenic, n={len(ig_hursts)}): mean={np.mean(ig_hursts):.4f}, std={np.std(ig_hursts):.4f}")

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY — Hypothesis Validation")
    print("=" * 70)

    checks = []

    # Check 1: MI half-life
    if cds_halflifes:
        med_hl = np.median(cds_halflifes)
        ok = 2.5 <= med_hl <= 4.0
        checks.append(("MI half-life ∈ [2.5, 4.0]", med_hl, ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] MI half-life (CDS median) = {med_hl:.2f}  (target: [2.5, 4.0])")

    # Check 2: Lag-3 autocorrelation
    lag3_pct = lag3_count / max(total_checked, 1)
    ok = lag3_pct > 0.5
    checks.append(("Lag-3 autocorr > 50% CDS", lag3_pct, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Lag-3 autocorr in CDS = {lag3_pct*100:.1f}%  (target: >50%)")

    # Check 3: H_ratio in CODING range
    if cds_ratios:
        mean_hr = np.mean(cds_ratios)
        ok = 0.60 <= mean_hr <= 0.85
        checks.append(("H_ratio(CDS) ∈ [0.60, 0.85]", mean_hr, ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] H_ratio (CDS mean) = {mean_hr:.4f}  (target: [0.60, 0.85])")

    # Check 4: Hurst persistent
    if cds_hursts:
        mean_hurst = np.mean(cds_hursts)
        ok = 0.55 <= mean_hurst <= 0.75
        checks.append(("Hurst(CDS) ∈ [0.55, 0.75]", mean_hurst, ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] Hurst (CDS mean) = {mean_hurst:.4f}  (target: [0.55, 0.75])")

    # Check 5: ΔH convergent before gene start
    if dh_profiles:
        dh_before = avg_dh[flank - 10:flank]  # 10 nt before start
        mean_dh_before = float(dh_before.mean())
        ok = mean_dh_before < 0
        checks.append(("ΔH < 0 before gene start", mean_dh_before, ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] ΔH (10nt before start) = {mean_dh_before:+.4f}  (target: < 0)")

    passed = sum(1 for _, _, ok in checks if ok)
    print(f"\n  Result: {passed}/{len(checks)} checks passed")
    print(f"  {'SUCCESS' if passed >= 3 else 'NEEDS INVESTIGATION'}: "
          f"{'proceed to 002' if passed >= 3 else 'review failed checks'}")


if __name__ == "__main__":
    main()
