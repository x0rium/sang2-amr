#!/usr/bin/env python3
"""Exploration 005: Ab initio AMR detection on pKpQIL (KPC-2 plasmid).

THE test: can SANG2 signals find blaKPC WITHOUT knowing what to look for?

Method:
  1. Compute all validated signals on sliding windows along the plasmid
  2. Compute z-scores relative to E. coli chromosome baseline
  3. Combine into composite anomaly score
  4. Rank regions, check if blaKPC is in top-10

Signals used (all validated in 001-004):
  - Hurst exponent (deviation from host)
  - 3-mer entropy (deviation from host)
  - Autocorrelation period-3 strength
  - MI(lag=1) (dinucleotide bias)
  - Signal gradient (GC + 3-mer H)
  - Chargaff-2 deviation
  - ρ dinucleotide distance from host

Ground truth: blaKPC-2 at 6620-7502 (882 bp).

Success criterion: blaKPC region in top-10 anomalous windows.
"""

from __future__ import annotations

import numpy as np
from collections import Counter
from pathlib import Path
from dataclasses import dataclass

from Bio import SeqIO

DATA_DIR = Path(__file__).parent / "data"
NT_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
ALPHABET_SIZE = 4
NT_CHARS = "ACGT"
COMPLEMENT = np.array([3, 2, 1, 0], dtype=np.uint8)


def load_seq(path: str) -> np.ndarray:
    rec = SeqIO.read(path, "fasta")
    seq = str(rec.seq).upper()
    arr = np.array([NT_MAP.get(c, 0) for c in seq], dtype=np.uint8)
    arr[arr >= 4] = 0
    return arr


@dataclass
class AnnotatedRegion:
    start: int
    end: int
    name: str
    kind: str  # 'amr', 'transposase', 'other'


def load_ground_truth(gb_path: str) -> list[AnnotatedRegion]:
    rec = SeqIO.read(gb_path, "genbank")
    regions = []
    for f in rec.features:
        if f.type != 'CDS':
            continue
        loc = f.location
        if loc is None:
            continue
        name = f.qualifiers.get('gene', f.qualifiers.get('product', ['?']))[0]
        name_lower = name.lower()
        if any(k in name_lower for k in ['kpc', 'carbapenem', 'beta-lactam', 'bla']):
            kind = 'amr'
        elif any(k in name_lower for k in ['transpos', 'resolv', 'integr', 'is']):
            kind = 'transposase'
        else:
            kind = 'other'
        regions.append(AnnotatedRegion(int(loc.start), int(loc.end), name, kind))
    return regions


# =========================================================================
# Signal functions (all validated in 001-004)
# =========================================================================

def entropy(counts) -> float:
    values = np.array(list(counts.values()) if isinstance(counts, (dict, Counter)) else counts, dtype=np.float64)
    total = values.sum()
    if total <= 1:
        return 0.0
    probs = values[values > 0] / total
    return float(-np.sum(probs * np.log2(probs)))


def gc_content(data: np.ndarray) -> float:
    return float(np.mean((data == 1) | (data == 2)))


def kmer_entropy_3(data: np.ndarray) -> float:
    n = len(data)
    if n < 3:
        return 0.0
    keys = data[:-2].astype(np.int64) * 16 + data[1:-1].astype(np.int64) * 4 + data[2:].astype(np.int64)
    hist = np.bincount(keys, minlength=64)
    h = entropy(hist)
    return float(h / np.log2(64))


def hurst(data: np.ndarray) -> float:
    arr = data.astype(np.float64)
    n = len(arr)
    if n < 64:
        return 0.5
    windows, rs_values = [], []
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


def autocorr_period3(data: np.ndarray, max_lag: int = 30) -> float:
    arr = data.astype(np.float64)
    n = len(arr)
    max_lag = min(max_lag, n // 4)
    if max_lag < 6:
        return 0.0
    centered = arr - arr.mean()
    norm = float(centered @ centered)
    if norm == 0:
        return 0.0
    ac = np.empty(max_lag, dtype=np.float64)
    for lag in range(max_lag):
        ac[lag] = (centered[:n - lag] @ centered[lag:]) / norm
    ac_abs = np.abs(ac[1:])
    total = ac_abs.sum()
    if total == 0:
        return 0.0
    p3 = sum(ac_abs[i] for i in range(2, len(ac_abs), 3))
    return float(p3 / total)


def mi_lag1(data: np.ndarray) -> float:
    n = len(data)
    if n < 10:
        return 0.0
    h_m = entropy(np.bincount(data, minlength=4))
    a, b = data[:-1].astype(np.int64), data[1:].astype(np.int64)
    h_j = entropy(np.bincount(a * 4 + b, minlength=16))
    return max(2 * h_m - h_j, 0.0)


def chargaff2_dev(data: np.ndarray) -> float:
    counts = np.bincount(data, minlength=4).astype(np.float64)
    total = counts.sum()
    if total == 0:
        return 0.0
    freq = counts / total
    return abs(freq[0] - freq[3]) + abs(freq[1] - freq[2])


def rho_vector(data: np.ndarray) -> np.ndarray:
    n = len(data)
    if n < 10:
        return np.ones(16)
    mono = np.bincount(data, minlength=4).astype(np.float64) / n
    di_keys = data[:-1].astype(np.int64) * 4 + data[1:].astype(np.int64)
    di_freq = np.bincount(di_keys, minlength=16).astype(np.float64) / (n - 1)
    rho = np.ones(16)
    for a in range(4):
        for b in range(4):
            exp = mono[a] * mono[b]
            if exp > 0:
                rho[a * 4 + b] = di_freq[a * 4 + b] / exp
    return rho


def rho_distance(data: np.ndarray, ref_rho: np.ndarray) -> float:
    return float(np.sqrt(np.sum((rho_vector(data) - ref_rho) ** 2)))


# =========================================================================
# Composite anomaly scorer
# =========================================================================

@dataclass
class WindowScore:
    position: int  # center of window
    start: int
    end: int
    scores: dict  # signal_name → z-score
    composite: float  # combined anomaly score


def compute_host_baseline(host_seq: np.ndarray, window: int = 2000, step: int = 1000):
    """Compute mean ± std for each signal on host genome windows."""
    n = len(host_seq)
    n_windows = (n - window) // step + 1
    # Sample up to 500 windows for speed
    indices = np.linspace(0, n_windows - 1, min(500, n_windows), dtype=int)

    signals = {
        'gc': [], 'h3': [], 'hurst': [], 'ac3': [],
        'mi1': [], 'c2': [],
    }

    for idx in indices:
        start = idx * step
        seg = host_seq[start:start + window]
        signals['gc'].append(gc_content(seg))
        signals['h3'].append(kmer_entropy_3(seg))
        signals['hurst'].append(hurst(seg))
        signals['ac3'].append(autocorr_period3(seg))
        signals['mi1'].append(mi_lag1(seg))
        signals['c2'].append(chargaff2_dev(seg))

    baseline = {}
    for name, vals in signals.items():
        arr = np.array(vals)
        baseline[name] = (arr.mean(), max(arr.std(), 1e-6))

    # Host ρ-vector (global)
    baseline['rho_ref'] = rho_vector(host_seq)

    return baseline


def scan_plasmid(plasmid_seq: np.ndarray, baseline: dict,
                 window: int = 2000, step: int = 500) -> list[WindowScore]:
    """Scan plasmid with sliding window, compute z-scores from host baseline."""
    n = len(plasmid_seq)
    n_windows = max(1, (n - window) // step + 1)
    rho_ref = baseline['rho_ref']

    results = []
    for i in range(n_windows):
        start = i * step
        end = start + window
        seg = plasmid_seq[start:end]
        center = start + window // 2

        # Compute signals
        raw = {
            'gc': gc_content(seg),
            'h3': kmer_entropy_3(seg),
            'hurst': hurst(seg),
            'ac3': autocorr_period3(seg),
            'mi1': mi_lag1(seg),
            'c2': chargaff2_dev(seg),
        }

        # Z-scores relative to host
        z_scores = {}
        for name, val in raw.items():
            mean, std = baseline[name]
            z_scores[name] = abs(val - mean) / std

        # ρ-distance (already a distance, just normalize)
        rho_dist = rho_distance(seg, rho_ref)
        z_scores['rho'] = rho_dist  # raw distance, no z-score needed

        # Signal gradient (requires neighbors)
        z_scores['gradient'] = 0.0  # computed after

        # Composite: geometric mean of z-scores (zero-safe)
        z_vals = np.array(list(z_scores.values()))
        z_vals = np.clip(z_vals, 0.01, None)  # floor to avoid log(0)
        composite = float(np.exp(np.log(z_vals).mean()))

        results.append(WindowScore(
            position=center, start=start, end=end,
            scores=z_scores, composite=composite,
        ))

    # Add gradient signal (diff between adjacent windows)
    for i in range(1, len(results)):
        gc_diff = abs(results[i].scores.get('gc', 0) - results[i - 1].scores.get('gc', 0))
        h3_diff = abs(results[i].scores.get('h3', 0) - results[i - 1].scores.get('h3', 0))
        results[i].scores['gradient'] = gc_diff + h3_diff

    # Recompute composite with gradient
    for r in results:
        z_vals = np.array(list(r.scores.values()))
        z_vals = np.clip(z_vals, 0.01, None)
        r.composite = float(np.exp(np.log(z_vals).mean()))

    results.sort(key=lambda w: -w.composite)
    return results


# =========================================================================
# Main
# =========================================================================
def main():
    print("=" * 70)
    print("EXPLORATION 005: Ab Initio AMR Detection on pKpQIL")
    print("=" * 70)

    # Load data
    print("\n[1] Loading data...")
    chr_seq = load_seq(str(DATA_DIR / "ecoli_k12.fasta"))
    plas_seq = load_seq(str(DATA_DIR / "pkpqil.fasta"))
    ground_truth = load_ground_truth(str(DATA_DIR / "pkpqil.gb"))

    print(f"  Host chromosome: {len(chr_seq):,} nt (E. coli K-12)")
    print(f"  Target plasmid: {len(plas_seq):,} nt (pKpQIL)")
    print(f"  Ground truth features: {len(ground_truth)}")

    # Show ground truth
    print(f"\n  Ground truth annotation:")
    kpc_region = None
    for r in ground_truth:
        marker = ""
        if r.kind == 'amr':
            marker = " <<<< TARGET"
            kpc_region = r
        elif r.kind == 'transposase':
            marker = " [transposase]"
        print(f"    {r.start:>6}-{r.end:>6} ({r.end-r.start:4d} bp) {r.name:35s}{marker}")

    if kpc_region is None:
        print("  ERROR: No KPC gene found in annotation!")
        return

    # ===================================================================
    # Compute host baseline
    # ===================================================================
    print("\n[2] Computing host baseline (E. coli K-12)...")
    baseline = compute_host_baseline(chr_seq, window=2000, step=2000)
    for name, val in baseline.items():
        if name != 'rho_ref':
            mean, std = val
            print(f"    {name:10s}: mean={mean:.5f}, std={std:.5f}")

    # ===================================================================
    # Scan plasmid
    # ===================================================================
    print("\n[3] Scanning plasmid (window=2000, step=500)...")
    results = scan_plasmid(plas_seq, baseline, window=2000, step=500)
    print(f"  Windows scanned: {len(results)}")

    # ===================================================================
    # Results
    # ===================================================================
    print("\n" + "=" * 70)
    print("[4] Top 20 Anomalous Regions")
    print("=" * 70)

    print(f"\n  {'Rank':>4s} {'Position':>10s} {'Composite':>10s} "
          f"{'GC':>5s} {'H3':>5s} {'Hurst':>5s} {'AC3':>5s} {'MI1':>5s} {'C2':>5s} {'ρ':>5s}  "
          f"Overlaps")
    print("  " + "-" * 95)

    kpc_rank = None
    proximity = 1000  # bp tolerance for overlap

    for rank, w in enumerate(results[:20], 1):
        # Check overlap with ground truth
        overlaps = []
        for r in ground_truth:
            if w.start <= r.end + proximity and w.end >= r.start - proximity:
                overlaps.append(f"{r.name[:20]}({r.kind})")

        overlap_str = ", ".join(overlaps) if overlaps else ""

        # Check if this window overlaps KPC
        if (w.start <= kpc_region.end + proximity and
                w.end >= kpc_region.start - proximity):
            if kpc_rank is None:
                kpc_rank = rank
            overlap_str = "** " + overlap_str

        print(f"  {rank:4d} {w.start:>5d}-{w.end:>5d} {w.composite:10.3f} "
              f"{w.scores.get('gc',0):5.1f} {w.scores.get('h3',0):5.1f} "
              f"{w.scores.get('hurst',0):5.1f} {w.scores.get('ac3',0):5.1f} "
              f"{w.scores.get('mi1',0):5.1f} {w.scores.get('c2',0):5.1f} "
              f"{w.scores.get('rho',0):5.2f}  {overlap_str}")

    # ===================================================================
    # KPC-specific analysis
    # ===================================================================
    print("\n" + "=" * 70)
    print("[5] KPC Region Analysis")
    print("=" * 70)

    print(f"  KPC annotation: {kpc_region.start}-{kpc_region.end} ({kpc_region.end-kpc_region.start} bp)")
    print(f"  KPC rank in anomaly list: {kpc_rank if kpc_rank else '>20'}")

    # Find ALL windows overlapping KPC
    kpc_windows = [w for w in results
                   if w.start <= kpc_region.end + 500 and w.end >= kpc_region.start - 500]
    if kpc_windows:
        best_kpc = kpc_windows[0]
        kpc_ranks = [i + 1 for i, w in enumerate(results) if w in kpc_windows]
        print(f"  Windows overlapping KPC: {len(kpc_windows)}")
        print(f"  Best KPC window rank: {kpc_ranks[0]}")
        print(f"  Best KPC composite: {best_kpc.composite:.3f}")
        print(f"  KPC window scores:")
        for name, val in best_kpc.scores.items():
            print(f"    {name:10s}: {val:.3f}")

    # ===================================================================
    # Per-feature recall
    # ===================================================================
    print("\n" + "=" * 70)
    print("[6] Per-Feature Recall (top-N anomalous windows)")
    print("=" * 70)

    for top_n in [5, 10, 20]:
        top_windows = results[:top_n]
        detected = set()
        for w in top_windows:
            for r in ground_truth:
                if w.start <= r.end + proximity and w.end >= r.start - proximity:
                    detected.add(r.name)

        amr_detected = sum(1 for r in ground_truth if r.kind == 'amr' and r.name in detected)
        tnp_detected = sum(1 for r in ground_truth if r.kind == 'transposase' and r.name in detected)
        total_amr = sum(1 for r in ground_truth if r.kind == 'amr')
        total_tnp = sum(1 for r in ground_truth if r.kind == 'transposase')

        print(f"\n  Top-{top_n}:")
        print(f"    AMR recall: {amr_detected}/{total_amr}")
        print(f"    Transposase recall: {tnp_detected}/{total_tnp}")
        print(f"    Features detected: {sorted(detected)}")

    # ===================================================================
    # Individual signal contribution
    # ===================================================================
    print("\n" + "=" * 70)
    print("[7] Which signals drive KPC detection?")
    print("=" * 70)

    # For each signal, rank by that signal alone and find KPC rank
    signal_names = ['gc', 'h3', 'hurst', 'ac3', 'mi1', 'c2', 'rho', 'gradient']
    for sname in signal_names:
        sorted_by_sig = sorted(results, key=lambda w: -w.scores.get(sname, 0))
        kpc_rank_sig = None
        for i, w in enumerate(sorted_by_sig, 1):
            if (w.start <= kpc_region.end + proximity and
                    w.end >= kpc_region.start - proximity):
                kpc_rank_sig = i
                break
        total = len(sorted_by_sig)
        print(f"  {sname:10s}: KPC rank = {kpc_rank_sig if kpc_rank_sig else '>all'}/{total}")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    checks = []

    # Check 1: KPC in top-10
    ok = kpc_rank is not None and kpc_rank <= 10
    checks.append(("KPC in top-10 anomalous", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] KPC rank: {kpc_rank if kpc_rank else '>20'} (target: ≤10)")

    # Check 2: KPC in top-5
    ok5 = kpc_rank is not None and kpc_rank <= 5
    checks.append(("KPC in top-5", ok5))
    print(f"  [{'PASS' if ok5 else 'FAIL'}] KPC in top-5: {kpc_rank if kpc_rank else '>20'}")

    # Check 3: At least 1 transposase in top-10
    top10_tnp = sum(1 for w in results[:10]
                    for r in ground_truth
                    if r.kind == 'transposase' and
                    w.start <= r.end + proximity and w.end >= r.start - proximity)
    ok = top10_tnp > 0
    checks.append(("Transposase in top-10", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Transposases in top-10 windows: {top10_tnp}")

    # Check 4: Composite anomaly > median for KPC region
    if kpc_windows:
        median_composite = np.median([w.composite for w in results])
        best_kpc_comp = kpc_windows[0].composite
        ok = best_kpc_comp > median_composite
        checks.append(("KPC composite > median", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] KPC composite ({best_kpc_comp:.3f}) > median ({median_composite:.3f})")

    passed = sum(1 for _, ok in checks if ok)
    print(f"\n  Result: {passed}/{len(checks)} checks passed")
    print(f"  {'SUCCESS' if passed >= 3 else 'NEEDS INVESTIGATION'}")


if __name__ == "__main__":
    main()
