"""Pipeline orchestration: wires L0-L5 into end-to-end workflows.

Three modes:
  A. Single isolate — find novel AMR on one genome
  B. Metagenome — cluster contigs, find AMR per cluster
  C. Cohort — N isolates with MIC data (stub)
"""

from __future__ import annotations

import numpy as np

from ir import Contig, AMRCandidate, GenomicRegion, SignalVector, WindowScore
from io_module import contig_to_array, load_fasta
from signals import (
    entropy, gc_content, kmer_entropy, hurst, mi_at_lag,
    autocorr_period3, autocorr_peaks, chargaff2_dev, sdr_kmer,
    otsu_threshold,
)
from anomaly import compute_host_baseline, scan_sequence
from context import apply_context_boost
from cluster import build_feature_matrix, cluster_hierarchical


def _compute_signal_vector(data: np.ndarray) -> SignalVector:
    """Compute all signals for a region."""
    return SignalVector(
        entropy=entropy(np.bincount(data, minlength=4)),
        sdr_ratio=sdr_kmer(data, k=4),
        delta_h=0.0,  # scalar summary not meaningful; use per-position delta_h separately
        autocorr_period3=autocorr_period3(data),
        autocorr_peaks=autocorr_peaks(data),
        kmer_entropy_3=kmer_entropy(data, k=3),
        regime=_classify_regime(data),
        hurst=hurst(data) if len(data) >= 64 else 0.5,
        gc_content=gc_content(data),
        chargaff2_dev=chargaff2_dev(data),
        mi_lag1=mi_at_lag(data, 1),
    )


def _classify_regime(data: np.ndarray) -> str:
    """Classify regime using 3-mer entropy + Otsu (data-driven)."""
    h3 = kmer_entropy(data, k=3)
    if h3 < 0.5:
        return "SPARSE"
    if h3 < 0.85:
        return "BIASED"
    if h3 < 0.95:
        return "CODING"
    return "BALANCED"


def _windows_to_candidates(
    windows: list[WindowScore],
    contig_id: str,
    top_n: int = 20,
) -> list[AMRCandidate]:
    """Convert top anomalous windows to AMRCandidate objects."""
    candidates = []
    for w in windows[:top_n]:
        evidence = []
        for sig_name, z_val in sorted(w.scores.items(), key=lambda x: -x[1]):
            if z_val > 1.5:
                evidence.append(f"{sig_name}: z={z_val:.1f}")

        region = GenomicRegion(
            contig_id=contig_id,
            start=w.start,
            end=w.end,
            strand="+",
            kind="unknown",
        )
        candidates.append(AMRCandidate(
            region=region,
            evidence=evidence,
            confidence=min(w.composite / 3.0, 1.0),  # normalize to [0, 1]
        ))
    return candidates


# =========================================================================
# Mode A: Single Isolate
# =========================================================================
def run_single_isolate(
    target_contigs: list[Contig],
    host_contigs: list[Contig] | None = None,
    window: int = 2000,
    step: int = 500,
    top_n: int = 20,
) -> list[AMRCandidate]:
    """Mode A: Find anomalous regions on a single genome/plasmid.

    If host_contigs provided, uses them as baseline.
    Otherwise, uses the target itself as baseline (self-comparison).
    """
    # Build host baseline
    if host_contigs:
        host_seq = np.concatenate([contig_to_array(c) for c in host_contigs])
    else:
        host_seq = np.concatenate([contig_to_array(c) for c in target_contigs])

    baseline = compute_host_baseline(host_seq, window=window, step=window)

    # Scan each target contig
    all_candidates: list[AMRCandidate] = []
    for contig in target_contigs:
        seq = contig_to_array(contig)
        if len(seq) < window:
            continue
        windows = scan_sequence(seq, baseline, window=window, step=step)
        windows = apply_context_boost(windows)
        candidates = _windows_to_candidates(windows, contig.id, top_n=top_n)
        all_candidates.extend(candidates)

    # Sort by confidence
    all_candidates.sort(key=lambda c: -c.confidence)
    return all_candidates[:top_n]


# =========================================================================
# Mode B: Metagenome
# =========================================================================
def run_metagenome(
    contigs: list[Contig],
    n_clusters: int | None = None,
    min_contig_len: int = 2000,
) -> dict:
    """Mode B: Cluster contigs by signal profile, find AMR per cluster.

    Returns dict with:
      'clusters': list of cluster labels
      'feature_matrix': numpy array (n_contigs, 22)
      'candidates_per_cluster': dict[int, list[AMRCandidate]]
    """
    # Filter short contigs
    valid = [(c, contig_to_array(c)) for c in contigs
             if len(c.sequence) >= min_contig_len]
    if not valid:
        return {"clusters": [], "feature_matrix": np.array([]), "candidates_per_cluster": {}}

    contig_list, seq_list = zip(*valid)
    X = build_feature_matrix(list(seq_list))

    # Auto-detect k if not provided (Otsu on GC for initial estimate)
    if n_clusters is None:
        gc_vals = X[:, 0].tolist()
        if len(gc_vals) >= 4:
            _, var = otsu_threshold(gc_vals)
            n_clusters = max(2, min(10, int(np.sqrt(len(gc_vals) / 10))))
        else:
            n_clusters = 2

    labels = cluster_hierarchical(X, n_clusters)

    # Per-cluster AMR scan
    candidates_per_cluster: dict[int, list[AMRCandidate]] = {}
    for k in range(n_clusters):
        mask = labels == k
        cluster_contigs = [contig_list[i] for i in range(len(contig_list)) if mask[i]]
        if cluster_contigs:
            candidates = run_single_isolate(cluster_contigs, top_n=10)
            candidates_per_cluster[k] = candidates

    return {
        "clusters": labels.tolist(),
        "feature_matrix": X,
        "candidates_per_cluster": candidates_per_cluster,
    }


# =========================================================================
# Mode C: Cohort (stub)
# =========================================================================
def run_cohort(
    isolates: list[list[Contig]],
    mic_values: list[float] | None = None,
) -> dict:
    """Mode C: Cohort analysis — correlate signals with phenotype.

    STUB: Not yet validated. Requires MIC data.
    """
    raise NotImplementedError(
        "Mode C (cohort) is not yet implemented. "
        "Requires validation on PATRIC/BV-BRC MIC data (see docs/data-sources.md)."
    )
