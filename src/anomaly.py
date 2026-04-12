"""Composite anomaly scorer for ab initio AMR detection.

THE CORE of SANG2-AMR: combines all validated signals into a single
anomaly score per sliding window. Regions with high composite score
are candidates for novel AMR genes, mobile elements, or foreign DNA.

Validated in exploration/005: blaKPC-2 found at rank 8/37 on pKpQIL.

Method:
  1. Compute host baseline (mean ± std for each signal on reference genome)
  2. Scan target sequence with sliding windows
  3. Per window: compute z-scores relative to host baseline
  4. Composite = geometric mean of |z-scores|
  5. Rank windows by composite → top = most anomalous

Signals (10 total):
  Window-level (from signals.py):
    gc, h3, hurst, ac3, mi1, c2, rho, gradient
  Resonance-derived (codon-level pair analysis):
    cpd  — codon pair divergence (JSD from host profile)
    rscu — codon usage deviation from host
"""

from __future__ import annotations

import numpy as np

from ir import WindowScore
from signals import (
    gc_content, kmer_entropy, hurst,
    mi_at_lag, autocorr_period3, chargaff2_dev,
    rho_vector, rho_distance,
    codon_pair_profile, codon_pair_distance, codon_mi_profile,
)
from genome_signals import rscu, codon_bias_distance


def compute_host_baseline(
    host_seq: np.ndarray,
    window: int = 2000,
    step: int = 1000,
    max_windows: int = 500,
) -> dict:
    """Compute signal statistics (mean, std) on host genome windows.

    Returns dict with:
      - signal_name → (mean, std) for scalar signals
      - 'rho_ref' → global ρ-vector of host
      - 'codon_pair_ref' → global codon-pair frequency profile
      - 'rscu_ref' → global RSCU vector of host CDS
    """
    n = len(host_seq)
    n_windows = (n - window) // step + 1
    indices = np.linspace(0, n_windows - 1, min(max_windows, n_windows), dtype=int)

    collectors: dict[str, list[float]] = {
        "gc": [], "h3": [], "hurst": [], "ac3": [],
        "mi1": [], "c2": [], "cpd": [], "rscu_d": [], "codon_mi": [],
        "rho_raw": [],
    }

    # Global host profiles for codon-level comparison
    host_cp_ref = codon_pair_profile(host_seq)
    host_rscu_ref = rscu(host_seq)

    for idx in indices:
        start = int(idx) * step
        seg = host_seq[start:start + window]
        collectors["gc"].append(gc_content(seg))
        collectors["h3"].append(kmer_entropy(seg, k=3))
        collectors["hurst"].append(hurst(seg))
        collectors["ac3"].append(autocorr_period3(seg))
        collectors["mi1"].append(mi_at_lag(seg, 1))
        collectors["c2"].append(chargaff2_dev(seg))
        collectors["cpd"].append(codon_pair_distance(seg, host_cp_ref))
        collectors["rscu_d"].append(codon_bias_distance(seg, host_rscu_ref))
        collectors["codon_mi"].append(codon_mi_profile(seg))
        collectors["rho_raw"].append(rho_distance(seg, rho_vector(host_seq)))

    baseline: dict = {}
    for name, vals in collectors.items():
        arr = np.array(vals)
        baseline[name] = (float(arr.mean()), max(float(arr.std()), 1e-6))

    # rho baseline stats for z-scoring
    rho_arr = np.array(collectors["rho_raw"])
    baseline["rho_stats"] = (float(rho_arr.mean()), max(float(rho_arr.std()), 1e-6))

    baseline["rho_ref"] = rho_vector(host_seq)
    baseline["codon_pair_ref"] = host_cp_ref
    baseline["rscu_ref"] = host_rscu_ref
    return baseline


def scan_sequence(
    target_seq: np.ndarray,
    baseline: dict,
    window: int = 2000,
    step: int = 500,
) -> list[WindowScore]:
    """Scan target sequence and rank windows by composite anomaly score.

    10 signals: 7 window-level + 3 resonance-derived codon-level.
    Composite = geometric mean of all z-scores.

    Returns list sorted by composite descending (most anomalous first).
    """
    n = len(target_seq)
    n_windows = max(1, (n - window) // step + 1)
    rho_ref = baseline["rho_ref"]
    cp_ref = baseline["codon_pair_ref"]
    rscu_ref = baseline["rscu_ref"]

    results: list[WindowScore] = []

    for i in range(n_windows):
        start = i * step
        end = min(start + window, n)
        seg = target_seq[start:end]
        center = start + window // 2

        # Window-level signals
        raw = {
            "gc": gc_content(seg),
            "h3": kmer_entropy(seg, k=3),
            "hurst": hurst(seg),
            "ac3": autocorr_period3(seg),
            "mi1": mi_at_lag(seg, 1),
            "c2": chargaff2_dev(seg),
        }

        # Z-scores relative to host
        z_scores: dict[str, float] = {}
        for name, val in raw.items():
            mean, std = baseline[name]
            z_scores[name] = abs(val - mean) / std

        # ρ-distance (dinucleotide level) — z-scored like everything else
        rho_val = rho_distance(seg, rho_ref)
        rho_mean, rho_std = baseline.get("rho_stats", (0.4, 0.1))
        z_scores["rho"] = abs(rho_val - rho_mean) / rho_std

        # Resonance-derived codon-level signals (z-scored)
        cpd_val = codon_pair_distance(seg, cp_ref)
        cpd_mean, cpd_std = baseline["cpd"]
        z_scores["cpd"] = abs(cpd_val - cpd_mean) / cpd_std

        rscu_d_val = codon_bias_distance(seg, rscu_ref)
        rscu_mean, rscu_std = baseline["rscu_d"]
        z_scores["rscu_d"] = abs(rscu_d_val - rscu_mean) / rscu_std

        cmi_val = codon_mi_profile(seg)
        cmi_mean, cmi_std = baseline["codon_mi"]
        z_scores["codon_mi"] = abs(cmi_val - cmi_mean) / cmi_std

        # Gradient placeholder (post-pass)
        z_scores["gradient"] = 0.0

        results.append(WindowScore(
            position=center, start=start, end=end,
            scores=z_scores, composite=0.0,
        ))

    # Post-pass: signal gradient between adjacent windows
    for i in range(1, len(results)):
        gc_diff = abs(results[i].scores.get("gc", 0) - results[i - 1].scores.get("gc", 0))
        h3_diff = abs(results[i].scores.get("h3", 0) - results[i - 1].scores.get("h3", 0))
        results[i].scores["gradient"] = gc_diff + h3_diff

    # TWO-TRACK SCORING: geometric resonance + rank fusion.
    # Each catches different anomaly patterns:
    #   Geometric: "everything slightly off" → recent HGT from distant organism
    #   Rank fusion: "one signal screaming" → close-relative gene with one extreme feature
    # Final score = best of both tracks.
    all_signal_names = [
        "gc", "h3", "hurst", "ac3", "mi1", "c2", "rho", "gradient",
        "cpd", "rscu_d", "codon_mi",
    ]

    n_total = len(results)
    if n_total == 0:
        return results

    # --- Track A: Geometric mean (SANG2 resonance philosophy) ---
    # All signals must contribute. Low values drag score down.
    # Catches: broadly anomalous regions.
    for r in results:
        z_vals = np.array([r.scores.get(s, 0.01) for s in all_signal_names])
        z_vals = np.clip(z_vals, 0.01, None)
        r.scores["_geomean"] = float(np.exp(np.log(z_vals).mean()))

    # --- Track B: Rank fusion (top-3 mean rank) ---
    # Each signal votes independently. One strong signal enough.
    # Catches: one-signal-dominant anomalies.
    n_signals = len(all_signal_names)
    rank_matrix = np.zeros((n_total, n_signals), dtype=np.float64)
    for j, sig_name in enumerate(all_signal_names):
        vals = np.array([r.scores.get(sig_name, 0.0) for r in results])
        order = np.argsort(-vals)
        ranks = np.empty_like(order, dtype=np.float64)
        ranks[order] = np.arange(1, n_total + 1, dtype=np.float64)
        rank_matrix[:, j] = ranks

    top_k = min(3, n_signals)
    for i, r in enumerate(results):
        sorted_ranks = np.sort(rank_matrix[i])
        mean_top3 = float(sorted_ranks[:top_k].mean())
        r.scores["_rankfusion"] = n_total / mean_top3

    # --- Combine: best of both tracks ---
    # Normalize both to same scale (percentile rank within each track)
    geomeans = np.array([r.scores["_geomean"] for r in results])
    rankfusions = np.array([r.scores["_rankfusion"] for r in results])

    geo_order = np.argsort(-geomeans)
    geo_ranks = np.empty_like(geo_order, dtype=np.float64)
    geo_ranks[geo_order] = np.arange(1, n_total + 1, dtype=np.float64)

    rf_order = np.argsort(-rankfusions)
    rf_ranks = np.empty_like(rf_order, dtype=np.float64)
    rf_ranks[rf_order] = np.arange(1, n_total + 1, dtype=np.float64)

    for i, r in enumerate(results):
        best_rank = min(geo_ranks[i], rf_ranks[i])
        r.composite = n_total / best_rank  # higher = better
        r.scores["_track"] = "geometric" if geo_ranks[i] <= rf_ranks[i] else "rankfusion"

    # Deduplicate overlapping windows
    results.sort(key=lambda w: -w.composite)
    deduped: list[WindowScore] = []
    used_positions: set[int] = set()
    merge_distance = step * 2
    for w in results:
        center = (w.start + w.end) // 2
        is_dup = any(abs(center - p) < merge_distance for p in used_positions)
        if not is_dup:
            deduped.append(w)
            used_positions.add(center)

    return deduped
