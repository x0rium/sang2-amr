#!/usr/bin/env python3
"""Exploration 002: Signal contrast — E. coli K-12 chromosome vs K. pneumoniae plasmid.

Goal: Verify that SANG2 signals can distinguish chromosomal from plasmid DNA,
and detect structural features (transposases, conjugal regions, IS-elements)
as signal anomalies.

Learnings from 001 applied:
  - Use 3-mer entropy instead of nucleotide H/H_max for regime classification
  - Use autocorrelation periodicity instead of MI half-life
  - Keep Hurst, ΔH, SDR as-is (they worked)

Success criteria (from exploration/README.md):
  ≥ 3 of 8 signals statistically differ between chromosome and plasmid (p < 0.01).

Datasets:
  - Chromosome: E. coli K-12 MG1655 (U00096.3, 4.6 Mb)
  - Plasmid: pCAV1392-131 (CP011577.1, 131 kb, K. pneumoniae)
"""

from __future__ import annotations

import numpy as np
from collections import Counter
from pathlib import Path
from scipy import stats

from Bio import SeqIO

DATA_DIR = Path(__file__).parent / "data"
ALPHABET_SIZE = 4
H_MAX = np.log2(ALPHABET_SIZE)

NT_MAP = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3,
          ord('a'): 0, ord('c'): 1, ord('g'): 2, ord('t'): 3}


def load_seq(fasta_path: str) -> np.ndarray:
    rec = SeqIO.read(fasta_path, "fasta")
    seq = str(rec.seq).upper()
    arr = np.array([NT_MAP.get(ord(c), 255) for c in seq], dtype=np.uint8)
    arr[arr >= 4] = np.random.randint(0, 4, size=int((arr >= 4).sum())).astype(np.uint8)
    return arr


def load_features(gb_path: str):
    from dataclasses import dataclass

    @dataclass
    class Feature:
        start: int
        end: int
        kind: str
        name: str

    rec = SeqIO.read(gb_path, "genbank")
    feats = []
    for f in rec.features:
        if f.type in ('CDS', 'gene', 'mobile_element', 'repeat_region'):
            loc = f.location
            if loc is None:
                continue
            name = f.qualifiers.get('gene',
                   f.qualifiers.get('product',
                   f.qualifiers.get('mobile_element_type', ['?'])))[0]
            feats.append(Feature(
                start=int(loc.start), end=int(loc.end),
                kind=f.type, name=name,
            ))
    return feats


# ---------------------------------------------------------------------------
# Signal functions (from 001, proven to work)
# ---------------------------------------------------------------------------

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


def gc_content(data: np.ndarray) -> float:
    return float(np.mean((data == 1) | (data == 2)))


def kmer_entropy(data: np.ndarray, k: int = 3) -> float:
    n = len(data)
    if n < k:
        return 0.0
    keys = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        keys += data[i:n - k + 1 + i].astype(np.int64) * (4 ** (k - 1 - i))
    hist = np.bincount(keys, minlength=4**k)
    h = entropy(hist)
    h_max = np.log2(4**k)
    return float(h / h_max)


def sdr_kmer(data: np.ndarray, k: int = 4):
    n = len(data)
    if n < k:
        return 0.0
    keys = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        keys += data[i:n - k + 1 + i].astype(np.int64) * (4 ** (k - 1 - i))
    freq = Counter(int(x) for x in keys)
    total = sum(freq.values())
    if total == 0:
        return 0.0
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
    return len(structural) / max(len(freq), 1)


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


def autocorr_period3_strength(data: np.ndarray, max_lag: int = 30) -> float:
    """Fraction of autocorrelation power at multiples of 3."""
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
    # Power at 3k vs total
    ac_abs = np.abs(ac[1:])  # skip lag 0
    total = ac_abs.sum()
    if total == 0:
        return 0.0
    period3 = sum(ac_abs[i] for i in range(2, len(ac_abs), 3))  # lags 3,6,9,...
    return float(period3 / total)


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


def delta_h_mean(data: np.ndarray, k: int = 5) -> float:
    """Mean absolute ΔH — measures directional structure."""
    n = len(data)
    if n < 2 * k + 2:
        return 0.0
    arr = data.astype(np.int64)

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
        row_sums = np.maximum(joint.sum(axis=1, keepdims=True), 1.0)
        p = joint / row_sums
        with np.errstate(divide="ignore", invalid="ignore"):
            log_p = np.where(p > 0, np.log2(p), 0.0)
        h_per_ctx = -(p * log_p).sum(axis=1)
        return h_per_ctx[inverse]

    fwd_ctx = hash_ctx(0, n - k)
    fwd_targets = arr[k:]
    h_fwd = context_h(fwd_ctx, fwd_targets)

    bwd_ctx = hash_ctx(1, n - k)
    bwd_targets = arr[:n - k]
    h_bwd = context_h(bwd_ctx, bwd_targets)

    lo, hi = k, n - k
    dh = h_fwd[:hi - lo] - h_bwd[lo:hi]
    return float(np.mean(np.abs(dh)))


def otsu_threshold(values) -> tuple[float, float]:
    vals = sorted(values)
    n = len(vals)
    if n < 2:
        return 0.0, 0.0
    total_sum = sum(vals)
    best_t, best_var = vals[0], 0.0
    left_sum = 0.0
    for i in range(1, n):
        left_sum += vals[i - 1]
        right_sum = total_sum - left_sum
        w0, w1 = i / n, (n - i) / n
        m0, m1 = left_sum / i, right_sum / (n - i)
        bv = w0 * w1 * (m0 - m1) ** 2
        if bv > best_var:
            best_var = bv
            best_t = (vals[i - 1] + vals[i]) / 2
    return float(best_t), float(best_var)


# ---------------------------------------------------------------------------
# Windowed signal computation
# ---------------------------------------------------------------------------
def compute_window_signals(data: np.ndarray, window: int = 2000, step: int = 500):
    """Compute all signals on sliding windows. Returns dict of arrays."""
    n = len(data)
    n_windows = max(1, (n - window) // step + 1)

    signals = {
        'gc': np.empty(n_windows),
        'h_nt': np.empty(n_windows),
        'h_3mer': np.empty(n_windows),
        'sdr_ratio': np.empty(n_windows),
        'hurst': np.empty(n_windows),
        'ac_period3': np.empty(n_windows),
        'mi_lag1': np.empty(n_windows),
        'dh_abs': np.empty(n_windows),
        'position': np.empty(n_windows),
    }

    for i in range(n_windows):
        start = i * step
        seg = data[start:start + window]
        signals['gc'][i] = gc_content(seg)
        signals['h_nt'][i] = entropy(np.bincount(seg, minlength=4)) / H_MAX
        signals['h_3mer'][i] = kmer_entropy(seg, k=3)
        signals['sdr_ratio'][i] = sdr_kmer(seg, k=4)
        signals['hurst'][i] = hurst(seg)
        signals['ac_period3'][i] = autocorr_period3_strength(seg)
        signals['mi_lag1'][i] = mi_at_lag(seg, 1)
        signals['dh_abs'][i] = delta_h_mean(seg, k=5)
        signals['position'][i] = start + window / 2

    return signals


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("EXPLORATION 002: Chromosome vs Plasmid Signal Contrast")
    print("=" * 70)

    # Load data
    print("\n[1] Loading sequences...")
    chr_seq = load_seq(str(DATA_DIR / "ecoli_k12.fasta"))
    plas_seq = load_seq(str(DATA_DIR / "pkpc_cav1321.fasta"))
    print(f"  Chromosome: {len(chr_seq):,} nt (E. coli K-12)")
    print(f"  Plasmid:    {len(plas_seq):,} nt (pCAV1392-131, K. pneumoniae)")

    plas_feats = load_features(str(DATA_DIR / "pkpc_cav1321.gb"))
    n_tnp = sum(1 for f in plas_feats if 'transposase' in f.name.lower())
    n_tra = sum(1 for f in plas_feats if f.name.lower().startswith('tra') or 'conjugal' in f.name.lower())
    print(f"  Plasmid features: {len(plas_feats)} total, {n_tnp} transposases, {n_tra} conjugal transfer")

    # ===================================================================
    # PART A: Global signal comparison
    # ===================================================================
    print("\n" + "=" * 70)
    print("[PART A] Global Signal Comparison")
    print("=" * 70)

    # Compute signals on equivalent-sized chunks
    # Sample chromosome in chunks matching plasmid size
    chunk_size = len(plas_seq)
    np.random.seed(42)

    # Chromosome: random non-overlapping chunks
    n_chr_chunks = len(chr_seq) // chunk_size
    chr_chunks = [chr_seq[i * chunk_size:(i + 1) * chunk_size] for i in range(n_chr_chunks)]

    # Plasmid: one big sequence, we can also use overlapping windows
    plas_window = 5000
    plas_step = 2000
    n_plas_windows = (len(plas_seq) - plas_window) // plas_step + 1

    print(f"\n  Chromosome: {n_chr_chunks} chunks of {chunk_size:,} bp")
    print(f"  Plasmid: {n_plas_windows} windows of {plas_window} bp (step {plas_step})")

    # Global stats on whole sequences
    print(f"\n  --- Global stats ---")

    for name, seq in [("Chromosome", chr_seq), ("Plasmid", plas_seq)]:
        gc = gc_content(seq)
        h = entropy(np.bincount(seq, minlength=4))
        h3 = kmer_entropy(seq, k=3)
        sdr_r = sdr_kmer(seq, k=4)
        hr = hurst(seq[:100_000] if len(seq) > 100_000 else seq)
        ac3 = autocorr_period3_strength(seq[:50_000] if len(seq) > 50_000 else seq)
        mi1 = mi_at_lag(seq[:50_000] if len(seq) > 50_000 else seq, 1)
        print(f"\n  {name}:")
        print(f"    GC content   = {gc:.4f}")
        print(f"    H (nt)       = {h:.4f} bits")
        print(f"    H/H_max (3m) = {h3:.4f}")
        print(f"    SDR ratio    = {sdr_r:.4f}")
        print(f"    Hurst        = {hr:.4f}")
        print(f"    AC period-3  = {ac3:.4f}")
        print(f"    MI(lag=1)    = {mi1:.5f}")

    # ===================================================================
    # PART B: Statistical comparison on matched windows
    # ===================================================================
    print("\n" + "=" * 70)
    print("[PART B] Statistical Comparison (windowed)")
    print("=" * 70)

    window = 5000
    step = 2000

    print(f"  Window: {window} bp, step: {step} bp")
    print("  Computing signals on chromosome windows...")
    chr_sigs = compute_window_signals(chr_seq, window=window, step=step)
    print(f"    {len(chr_sigs['gc'])} windows")

    print("  Computing signals on plasmid windows...")
    plas_sigs = compute_window_signals(plas_seq, window=window, step=step)
    print(f"    {len(plas_sigs['gc'])} windows")

    # Statistical tests
    sig_names = ['gc', 'h_3mer', 'sdr_ratio', 'hurst', 'ac_period3', 'mi_lag1', 'dh_abs']
    sig_labels = ['GC content', '3-mer H_ratio', 'SDR ratio', 'Hurst', 'AC period-3', 'MI(lag=1)', '|ΔH| mean']

    print(f"\n  {'Signal':20s} {'Chr mean':>10s} {'Plas mean':>10s} {'Diff':>8s} {'p-value':>12s} {'Sig?':>5s}")
    print("  " + "-" * 70)

    sig_results = []
    for sname, label in zip(sig_names, sig_labels):
        c_vals = chr_sigs[sname]
        p_vals = plas_sigs[sname]
        c_mean = c_vals.mean()
        p_mean = p_vals.mean()
        # Mann-Whitney U (non-parametric, different sample sizes)
        u_stat, p_val = stats.mannwhitneyu(c_vals, p_vals, alternative='two-sided')
        sig = p_val < 0.01
        sig_results.append((label, c_mean, p_mean, p_val, sig))
        diff = p_mean - c_mean
        print(f"  {label:20s} {c_mean:10.5f} {p_mean:10.5f} {diff:+8.5f} {p_val:12.2e} {'***' if sig else ''}")

    n_sig = sum(1 for _, _, _, _, s in sig_results if s)
    print(f"\n  Significant (p < 0.01): {n_sig}/{len(sig_results)}")

    # ===================================================================
    # PART C: Sliding signal profile along plasmid
    # ===================================================================
    print("\n" + "=" * 70)
    print("[PART C] Sliding Signal Profile Along Plasmid")
    print("=" * 70)
    print("  Looking for regime shifts that correspond to transposases/IS elements")

    # Compute detailed profile
    detail_window = 2000
    detail_step = 500
    plas_detail = compute_window_signals(plas_seq, window=detail_window, step=detail_step)

    # Compute chromosome baseline (mean ± 2σ for each signal)
    chr_detail = compute_window_signals(chr_seq[:500_000], window=detail_window, step=detail_step)

    print(f"\n  Anomaly detection: plasmid windows outside chr mean ± 2σ")
    print(f"  Chromosome baseline: {len(chr_detail['gc'])} windows of {detail_window} bp")
    print(f"  Plasmid profile: {len(plas_detail['gc'])} windows of {detail_window} bp")

    for sname, label in zip(sig_names, sig_labels):
        c_mean = chr_detail[sname].mean()
        c_std = chr_detail[sname].std()
        p_vals = plas_detail[sname]
        anomalous = np.sum((p_vals < c_mean - 2 * c_std) | (p_vals > c_mean + 2 * c_std))
        pct = anomalous / len(p_vals) * 100
        print(f"  {label:20s}: {anomalous:3d}/{len(p_vals)} anomalous ({pct:.1f}%)")

    # ===================================================================
    # PART D: Detect transposase regions via signal discontinuities
    # ===================================================================
    print("\n" + "=" * 70)
    print("[PART D] Transposase Detection via Signal Discontinuities")
    print("=" * 70)

    # Transposase locations from annotation
    tnp_regions = [(f.start, f.end, f.name) for f in plas_feats if 'transposase' in f.name.lower()]
    print(f"  Annotated transposases: {len(tnp_regions)}")
    for s, e, name in tnp_regions[:5]:
        print(f"    {s:>6}-{e:>6} ({e-s:4d} bp) {name}")
    if len(tnp_regions) > 5:
        print(f"    ... and {len(tnp_regions)-5} more")

    # Find signal discontinuities (gradient magnitude)
    print(f"\n  Signal gradient peaks (top discontinuities in GC + 3-mer H):")
    gc_profile = plas_detail['gc']
    h3_profile = plas_detail['h_3mer']
    positions = plas_detail['position']

    # Combined gradient: |Δ(gc)| + |Δ(h3)| normalized
    gc_grad = np.abs(np.diff(gc_profile))
    h3_grad = np.abs(np.diff(h3_profile))
    gc_grad_norm = gc_grad / (gc_grad.max() + 1e-10)
    h3_grad_norm = h3_grad / (h3_grad.max() + 1e-10)
    combined_grad = gc_grad_norm + h3_grad_norm

    # Top 15 discontinuity positions
    top_idx = np.argsort(combined_grad)[::-1][:15]
    top_positions = positions[top_idx].astype(int)

    # Check if they're near annotated transposases
    hits = 0
    proximity = 3000  # bp tolerance
    print(f"\n  Top 15 signal discontinuities (proximity tolerance = {proximity} bp):")
    for idx in top_idx:
        pos = int(positions[idx])
        grad_val = combined_grad[idx]
        # Check if near a transposase
        near_tnp = ""
        for ts, te, tn in tnp_regions:
            if abs(pos - ts) < proximity or abs(pos - te) < proximity:
                near_tnp = f"  <-- near {tn} ({ts}-{te})"
                hits += 1
                break
        print(f"    pos {pos:>7}: gradient = {grad_val:.4f}{near_tnp}")

    print(f"\n  Hits near transposases: {hits}/15")

    # Reverse: for each transposase, is there a nearby gradient peak?
    gradient_threshold = np.percentile(combined_grad, 90)  # top 10% gradients
    tnp_detected = 0
    print(f"\n  Transposase recall (gradient > 90th percentile within {proximity} bp):")
    for ts, te, tn in tnp_regions:
        detected = False
        for i, grad in enumerate(combined_grad):
            pos = int(positions[i])
            if (abs(pos - ts) < proximity or abs(pos - te) < proximity) and grad > gradient_threshold:
                detected = True
                break
        if detected:
            tnp_detected += 1
        status = "DETECTED" if detected else "MISSED"
        print(f"    {ts:>6}-{te:>6} {tn[:30]:30s} [{status}]")

    print(f"\n  Transposase recall: {tnp_detected}/{len(tnp_regions)} ({tnp_detected/max(len(tnp_regions),1)*100:.0f}%)")

    # ===================================================================
    # PART E: Conjugal transfer region as coherent block
    # ===================================================================
    print("\n" + "=" * 70)
    print("[PART E] Conjugal Transfer Region — Coherent Signal Block")
    print("=" * 70)

    tra_regions = [(f.start, f.end) for f in plas_feats if 'conjugal' in f.name.lower()]
    if tra_regions:
        tra_start = min(s for s, e in tra_regions)
        tra_end = max(e for s, e in tra_regions)
        print(f"  Conjugal region: {tra_start:,}-{tra_end:,} ({tra_end-tra_start:,} bp, {len(tra_regions)} genes)")

        # Signal inside vs outside conjugal region
        tra_mask = (plas_detail['position'] >= tra_start) & (plas_detail['position'] <= tra_end)
        non_tra_mask = ~tra_mask

        if tra_mask.sum() > 2 and non_tra_mask.sum() > 2:
            print(f"\n  {'Signal':20s} {'tra region':>10s} {'rest':>10s} {'p-value':>12s}")
            print("  " + "-" * 55)
            for sname, label in zip(sig_names, sig_labels):
                tra_vals = plas_detail[sname][tra_mask]
                rest_vals = plas_detail[sname][non_tra_mask]
                _, p_val = stats.mannwhitneyu(tra_vals, rest_vals, alternative='two-sided')
                print(f"  {label:20s} {tra_vals.mean():10.5f} {rest_vals.mean():10.5f} {p_val:12.2e}")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    checks = []

    # Check 1: ≥ 3 signals statistically different
    ok = n_sig >= 3
    checks.append(("≥3 signals differ (p<0.01)", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Significant signals: {n_sig}/7  (target: ≥3)")

    # Check 2: GC content differs
    gc_diff = abs(gc_content(chr_seq) - gc_content(plas_seq))
    ok = gc_diff > 0.01
    checks.append(("GC content differs", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] GC difference: {gc_diff:.4f}  (target: >0.01)")

    # Check 3: Hurst differs
    chr_hurst = chr_detail['hurst'].mean()
    plas_hurst = plas_detail['hurst'].mean()
    hurst_diff = abs(chr_hurst - plas_hurst)
    ok = hurst_diff > 0.01
    checks.append(("Hurst differs", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Hurst: chr={chr_hurst:.4f}, plas={plas_hurst:.4f}, diff={hurst_diff:.4f}")

    # Check 4: Transposase recall > 30%
    recall = tnp_detected / max(len(tnp_regions), 1)
    ok = recall > 0.3
    checks.append(("Transposase recall > 30%", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Transposase recall: {recall*100:.0f}%  (target: >30%)")

    # Check 5: Conjugal region signal differs from rest
    if tra_regions and tra_mask.sum() > 2:
        gc_tra = plas_detail['gc'][tra_mask].mean()
        gc_rest = plas_detail['gc'][non_tra_mask].mean()
        ok = abs(gc_tra - gc_rest) > 0.01
        checks.append(("Conjugal region GC differs", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] Conjugal GC: tra={gc_tra:.4f}, rest={gc_rest:.4f}")

    passed = sum(1 for _, ok in checks if ok)
    print(f"\n  Result: {passed}/{len(checks)} checks passed")
    print(f"  {'SUCCESS — proceed to 003' if passed >= 3 else 'NEEDS INVESTIGATION'}")


if __name__ == "__main__":
    main()
