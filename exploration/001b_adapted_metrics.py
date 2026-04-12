#!/usr/bin/env python3
"""Exploration 001b: Adapted metrics after initial findings.

Findings from 001:
  1. MI half-life doesn't work for 4-letter alphabet — MI oscillates with
     period 3 (codon structure) instead of decaying exponentially.
  2. Autocorrelation lag-3 is present but often weaker than harmonics (lag 9).
  3. H/H_max ≈ 0.99 for everything — 4-symbol alphabet can't produce low H.

Adapted metrics:
  1. MI Periodicity Score: MI(3)/mean(MI(1),MI(2),MI(4)) — detects codon signal
  2. Autocorrelation: check for "strongest peak at multiple of 3"
  3. k-mer entropy (k=3): effective alphabet of 64 codons → more separation
  4. Codon-level MI half-life: treat codons (not nucleotides) as symbols
"""

from __future__ import annotations

import numpy as np
from collections import Counter
from pathlib import Path

from Bio import SeqIO

DATA_DIR = Path(__file__).parent / "data"

NT_MAP = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3,
          ord('a'): 0, ord('c'): 1, ord('g'): 2, ord('t'): 3}

ALPHABET_SIZE = 4


def load_genome() -> np.ndarray:
    rec = SeqIO.read(str(DATA_DIR / "ecoli_k12.fasta"), "fasta")
    seq = str(rec.seq).upper()
    arr = np.array([NT_MAP.get(ord(c), 255) for c in seq], dtype=np.uint8)
    arr[arr >= 4] = np.random.randint(0, 4, size=int((arr >= 4).sum())).astype(np.uint8)
    return arr


def load_cds():
    from dataclasses import dataclass

    @dataclass
    class Region:
        start: int
        end: int
        strand: int
        name: str

    rec = SeqIO.read(str(DATA_DIR / "ecoli_k12.gb"), "genbank")
    regions = []
    for feat in rec.features:
        if feat.type == 'CDS':
            loc = feat.location
            if loc is None:
                continue
            name = feat.qualifiers.get('gene', feat.qualifiers.get('product', ['?']))[0]
            regions.append(Region(
                start=int(loc.start), end=int(loc.end),
                strand=int(loc.strand) if loc.strand else 1, name=name,
            ))
    return regions


def entropy(counts) -> float:
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


# ---------------------------------------------------------------------------
# ADAPTED METRIC 1: MI Periodicity Score
# ---------------------------------------------------------------------------
def mi_at_lag(data: np.ndarray, lag: int) -> float:
    n = len(data)
    if n <= lag + 1:
        return 0.0
    marg_hist = np.bincount(data, minlength=ALPHABET_SIZE)
    h_marg = entropy(marg_hist)
    a = data[:n - lag].astype(np.int64)
    b = data[lag:].astype(np.int64)
    joint_key = a * ALPHABET_SIZE + b
    joint_hist = np.bincount(joint_key, minlength=ALPHABET_SIZE**2)
    h_joint = entropy(joint_hist)
    return max(2 * h_marg - h_joint, 0.0)


def mi_periodicity_score(data: np.ndarray, period: int = 3) -> float:
    """Ratio of MI at target period vs surrounding lags.

    Score > 1.0 means the target period is over-represented in MI structure.
    For coding DNA with period=3, we expect score > 1.0 from codon structure.
    """
    max_check = period * 4
    mis = [mi_at_lag(data, lag) for lag in range(1, max_check + 1)]
    if not mis or max(mis) == 0:
        return 0.0

    # MI at multiples of period
    period_mis = [mis[i] for i in range(period - 1, len(mis), period)]
    # MI at non-multiples
    other_mis = [mis[i] for i in range(len(mis)) if (i + 1) % period != 0]

    if not other_mis or np.mean(other_mis) == 0:
        return float('inf') if period_mis else 0.0

    return float(np.mean(period_mis) / np.mean(other_mis))


def mi_lag3_ratio(data: np.ndarray) -> float:
    """Simple ratio: MI(lag=3) / mean(MI(lag=1), MI(lag=2), MI(lag=4), MI(lag=5)).

    > 1.0 indicates codon periodicity.
    """
    mi = [mi_at_lag(data, lag) for lag in range(1, 7)]
    if len(mi) < 6:
        return 0.0
    baseline = np.mean([mi[0], mi[1], mi[3], mi[4]])  # lags 1,2,4,5
    if baseline == 0:
        return float('inf') if mi[2] > 0 else 0.0
    return float(mi[2] / baseline)


# ---------------------------------------------------------------------------
# ADAPTED METRIC 2: Autocorrelation with "multiple of 3" check
# ---------------------------------------------------------------------------
def autocorrelation(data: np.ndarray, max_lag: int = 30):
    arr = data.astype(np.float64)
    n = len(arr)
    max_lag = max(2, min(max_lag, n // 4))
    centered = arr - arr.mean()
    norm = float(centered @ centered)
    if norm == 0:
        return np.zeros(max_lag), []
    ac = np.empty(max_lag, dtype=np.float64)
    for lag in range(max_lag):
        ac[lag] = (centered[:n - lag] @ centered[lag:]) / norm
    peaks = []
    for i in range(2, max_lag - 1):
        if ac[i] > ac[i - 1] and ac[i] > ac[i + 1] and ac[i] > 0.003:
            peaks.append((int(i), float(ac[i])))
    return ac, peaks


def strongest_peak_is_multiple_of_3(data: np.ndarray) -> bool:
    """Check if the strongest autocorrelation peak is at a multiple of 3."""
    _, peaks = autocorrelation(data, max_lag=30)
    if not peaks:
        return False
    best = max(peaks, key=lambda p: p[1])
    return best[0] % 3 == 0


def any_peak_at_multiple_of_3(data: np.ndarray) -> bool:
    """Check if any significant peak is at a multiple of 3."""
    _, peaks = autocorrelation(data, max_lag=30)
    return any(p[0] % 3 == 0 for p in peaks)


# ---------------------------------------------------------------------------
# ADAPTED METRIC 3: k-mer Entropy
# ---------------------------------------------------------------------------
def kmer_entropy(data: np.ndarray, k: int = 3) -> float:
    """Entropy of k-mer distribution. Effective alphabet = 4^k."""
    n = len(data)
    if n < k:
        return 0.0
    keys = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        keys += data[i:n - k + 1 + i].astype(np.int64) * (4 ** (k - 1 - i))
    hist = np.bincount(keys, minlength=4**k)
    h = entropy(hist)
    h_max = np.log2(4**k)  # k * 2.0
    return float(h / h_max)


# ---------------------------------------------------------------------------
# ADAPTED METRIC 4: Codon-level MI half-life
# ---------------------------------------------------------------------------
def to_codons(data: np.ndarray) -> np.ndarray:
    """Convert nucleotide array to codon array (base-4 triplets → 0..63)."""
    n = len(data)
    n_codons = n // 3
    if n_codons == 0:
        return np.array([], dtype=np.int64)
    trimmed = data[:n_codons * 3].reshape(n_codons, 3).astype(np.int64)
    return trimmed[:, 0] * 16 + trimmed[:, 1] * 4 + trimmed[:, 2]


def codon_mi_at_lag(codons: np.ndarray, lag: int) -> float:
    """MI between codons at given lag (codon-level, alphabet=64)."""
    n = len(codons)
    if n <= lag + 1:
        return 0.0
    marg = np.bincount(codons, minlength=64)
    h_marg = entropy(marg)
    a = codons[:n - lag].astype(np.int64)
    b = codons[lag:].astype(np.int64)
    joint = a * 64 + b
    h_joint = entropy(np.bincount(joint, minlength=64 * 64))
    return max(2 * h_marg - h_joint, 0.0)


def codon_mi_halflife(data: np.ndarray, max_lag: int = 32) -> float:
    """MI half-life on codon level (alphabet=64, one symbol = one codon)."""
    codons = to_codons(data)
    n = len(codons)
    if n < 10:
        return 0.0
    max_lag = min(max_lag, n // 4)
    lags = np.arange(1, max_lag + 1)
    mi = np.array([codon_mi_at_lag(codons, int(k)) for k in lags])
    if len(mi) < 3 or mi[0] <= 0:
        return 0.0
    mi_safe = np.clip(mi, 1e-10, None)
    log_mi = np.log(mi_safe)
    A = np.vstack([lags.astype(np.float64), np.ones(len(lags))]).T
    slope, _ = np.linalg.lstsq(A, log_mi, rcond=None)[0]
    if slope >= -1e-6:
        return float("inf")
    return float(-np.log(2) / slope)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("EXPLORATION 001b: Adapted Metrics")
    print("=" * 70)

    genome = load_genome()
    cds_regions = load_cds()

    np.random.seed(42)
    sample_idx = np.random.choice(len(cds_regions), min(300, len(cds_regions)), replace=False)
    sampled_cds = [cds_regions[i] for i in sample_idx]

    # Sort for intergenic gaps
    sorted_cds = sorted(cds_regions, key=lambda r: r.start)

    # ===================================================================
    # METRIC 1: MI Periodicity Score
    # ===================================================================
    print("\n" + "=" * 70)
    print("[M1] MI Periodicity Score (replaces MI half-life)")
    print("=" * 70)
    print("  Score = MI(multiples of 3) / MI(non-multiples of 3)")
    print("  Coding DNA: score > 1.0 (codon structure)")
    print()

    cds_scores = []
    for reg in sampled_cds:
        seg = genome[reg.start:reg.end]
        if len(seg) < 100:
            continue
        score = mi_periodicity_score(seg, period=3)
        if 0 < score < 100:
            cds_scores.append(score)

    ig_scores = []
    for i in range(len(sorted_cds) - 1):
        s, e = sorted_cds[i].end, sorted_cds[i + 1].start
        if e - s >= 100:
            score = mi_periodicity_score(genome[s:e], period=3)
            if 0 < score < 100:
                ig_scores.append(score)
        if len(ig_scores) >= 200:
            break

    if cds_scores:
        arr = np.array(cds_scores)
        print(f"  MI Periodicity (CDS, n={len(cds_scores)}):")
        print(f"    mean   = {arr.mean():.3f}")
        print(f"    median = {np.median(arr):.3f}")
        print(f"    > 1.0  : {np.sum(arr > 1.0)}/{len(arr)} ({np.sum(arr > 1.0)/len(arr)*100:.1f}%)")

    if ig_scores:
        arr = np.array(ig_scores)
        print(f"  MI Periodicity (intergenic, n={len(ig_scores)}):")
        print(f"    mean   = {arr.mean():.3f}")
        print(f"    median = {np.median(arr):.3f}")
        print(f"    > 1.0  : {np.sum(arr > 1.0)}/{len(arr)} ({np.sum(arr > 1.0)/len(arr)*100:.1f}%)")

    # Simple lag-3 ratio
    print("\n  MI lag-3 ratio (MI(3) / mean(MI(1,2,4,5))):")
    cds_r3 = []
    for reg in sampled_cds:
        seg = genome[reg.start:reg.end]
        if len(seg) < 100:
            continue
        r = mi_lag3_ratio(seg)
        if 0 < r < 100:
            cds_r3.append(r)

    ig_r3 = []
    for i in range(len(sorted_cds) - 1):
        s, e = sorted_cds[i].end, sorted_cds[i + 1].start
        if e - s >= 100:
            r = mi_lag3_ratio(genome[s:e])
            if 0 < r < 100:
                ig_r3.append(r)
        if len(ig_r3) >= 200:
            break

    if cds_r3:
        arr = np.array(cds_r3)
        print(f"    CDS (n={len(cds_r3)}): mean={arr.mean():.3f}, median={np.median(arr):.3f}")
    if ig_r3:
        arr = np.array(ig_r3)
        print(f"    Intergenic (n={len(ig_r3)}): mean={arr.mean():.3f}, median={np.median(arr):.3f}")

    # ===================================================================
    # METRIC 2: Autocorrelation — multiples of 3
    # ===================================================================
    print("\n" + "=" * 70)
    print("[M2] Autocorrelation — multiples of 3")
    print("=" * 70)

    strongest_m3 = 0
    any_m3 = 0
    total = 0
    for reg in cds_regions:
        if reg.end - reg.start < 200:
            continue
        seg = genome[reg.start:reg.end]
        if strongest_peak_is_multiple_of_3(seg):
            strongest_m3 += 1
        if any_peak_at_multiple_of_3(seg):
            any_m3 += 1
        total += 1
        if total >= 500:
            break

    print(f"  CDS (n={total}):")
    print(f"    Strongest peak at 3k: {strongest_m3} ({strongest_m3/max(total,1)*100:.1f}%)")
    print(f"    Any peak at 3k:       {any_m3} ({any_m3/max(total,1)*100:.1f}%)")

    # ===================================================================
    # METRIC 3: k-mer Entropy (k=3)
    # ===================================================================
    print("\n" + "=" * 70)
    print("[M3] Codon-level Entropy (k=3 mer entropy / H_max)")
    print("=" * 70)
    print("  H_max(3-mer) = log2(64) = 6.0 bits")

    cds_ke = []
    for reg in sampled_cds:
        seg = genome[reg.start:reg.end]
        if len(seg) >= 50:
            cds_ke.append(kmer_entropy(seg, k=3))

    ig_ke = []
    for i in range(len(sorted_cds) - 1):
        s, e = sorted_cds[i].end, sorted_cds[i + 1].start
        if e - s >= 50:
            ig_ke.append(kmer_entropy(genome[s:e], k=3))
        if len(ig_ke) >= 200:
            break

    if cds_ke:
        arr = np.array(cds_ke)
        print(f"  3-mer H_ratio (CDS, n={len(cds_ke)}): mean={arr.mean():.4f}, std={arr.std():.4f}")
    if ig_ke:
        arr = np.array(ig_ke)
        print(f"  3-mer H_ratio (intergenic, n={len(ig_ke)}): mean={arr.mean():.4f}, std={arr.std():.4f}")

    # Otsu on 3-mer entropy for regime separation
    if cds_ke and ig_ke:
        all_ke = cds_ke + ig_ke
        vals = sorted(all_ke)
        n_v = len(vals)
        total_sum = sum(vals)
        best_t, best_var = vals[0], 0.0
        left_sum = 0.0
        for i_v in range(1, n_v):
            left_sum += vals[i_v - 1]
            right_sum = total_sum - left_sum
            w0, w1 = i_v / n_v, (n_v - i_v) / n_v
            m0, m1 = left_sum / i_v, right_sum / (n_v - i_v)
            bv = w0 * w1 * (m0 - m1) ** 2
            if bv > best_var:
                best_var = bv
                best_t = (vals[i_v - 1] + vals[i_v]) / 2
        thresh, var = best_t, best_var
        print(f"  Otsu threshold (3-mer H_ratio): {thresh:.4f}, var={var:.6f}")

    # ===================================================================
    # METRIC 4: Codon-level MI half-life
    # ===================================================================
    print("\n" + "=" * 70)
    print("[M4] Codon-level MI Half-life")
    print("=" * 70)
    print("  Working on codon alphabet (64 symbols) instead of nt (4 symbols)")
    print("  Prediction: shorter half-life for CDS (structured codon usage)")

    cds_chl = []
    for reg in sampled_cds:
        seg = genome[reg.start:reg.end]
        if len(seg) < 150:
            continue
        hl = codon_mi_halflife(seg, max_lag=20)
        if 0 < hl < 100:
            cds_chl.append(hl)

    ig_chl = []
    for i in range(len(sorted_cds) - 1):
        s, e = sorted_cds[i].end, sorted_cds[i + 1].start
        if e - s >= 150:
            hl = codon_mi_halflife(genome[s:e], max_lag=20)
            if 0 < hl < 100:
                ig_chl.append(hl)
        if len(ig_chl) >= 100:
            break

    if cds_chl:
        arr = np.array(cds_chl)
        print(f"  Codon MI half-life (CDS, n={len(cds_chl)}):")
        print(f"    mean   = {arr.mean():.2f}")
        print(f"    median = {np.median(arr):.2f}")
        print(f"    std    = {arr.std():.2f}")
    if ig_chl:
        arr = np.array(ig_chl)
        print(f"  Codon MI half-life (intergenic, n={len(ig_chl)}):")
        print(f"    mean   = {arr.mean():.2f}")
        print(f"    median = {np.median(arr):.2f}")
        print(f"    std    = {arr.std():.2f}")

    # Show MI decay for dnaA at codon level
    dnaa = [r for r in cds_regions if 'dnaa' in r.name.lower()]
    if dnaa:
        seg = genome[dnaa[0].start:dnaa[0].end]
        codons = to_codons(seg)
        print(f"\n  Codon-level MI decay for dnaA ({len(codons)} codons):")
        for lag in range(1, min(16, len(codons) // 4)):
            mi = codon_mi_at_lag(codons, lag)
            bar = "#" * int(mi * 40)
            print(f"    lag {lag:2d}: MI={mi:.4f} {bar}")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print("\n" + "=" * 70)
    print("SUMMARY — Adapted Metrics")
    print("=" * 70)

    checks = []

    # M1: MI periodicity
    if cds_scores and ig_scores:
        cds_med = np.median(cds_scores)
        ig_med = np.median(ig_scores)
        ok = cds_med > 1.0 and cds_med > ig_med
        checks.append(("MI periodicity: CDS > 1.0 and > intergenic", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] MI period-3 score: CDS={cds_med:.3f}, IG={ig_med:.3f}")

    # M2: Autocorrelation multiples of 3
    ok = any_m3 / max(total, 1) > 0.5
    checks.append(("Any peak at lag=3k in >50% CDS", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Autocorr peak at 3k: {any_m3/max(total,1)*100:.1f}%")

    # M3: 3-mer entropy separates CDS from intergenic
    if cds_ke and ig_ke:
        cds_mean = np.mean(cds_ke)
        ig_mean = np.mean(ig_ke)
        ok = cds_mean < ig_mean  # CDS should have lower entropy (codon bias)
        checks.append(("3-mer H_ratio: CDS < intergenic", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] 3-mer H_ratio: CDS={cds_mean:.4f}, IG={ig_mean:.4f}")

    # M4: Codon MI half-life differs
    if cds_chl and ig_chl:
        cds_med_chl = np.median(cds_chl)
        ig_med_chl = np.median(ig_chl)
        ok = abs(cds_med_chl - ig_med_chl) > 0.5
        checks.append(("Codon MI half-life: CDS ≠ intergenic", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] Codon MI HL: CDS={cds_med_chl:.2f}, IG={ig_med_chl:.2f}")

    passed = sum(1 for _, ok in checks if ok)
    print(f"\n  Result: {passed}/{len(checks)} adapted checks passed")


if __name__ == "__main__":
    main()
