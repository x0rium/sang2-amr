#!/usr/bin/env python3
"""Exploration 006: Reference-free metagenomic clustering.

In silico mock metagenome: fragment 5 genomes, mix, cluster back by
SANG2 signal profiles alone (no reference DB, no alignment).

Species (wide GC range):
  - E. coli K-12:     4.6 Mb, GC=50.8%
  - B. subtilis 168:  4.2 Mb, GC=43.5%
  - S. aureus 8325:   2.8 Mb, GC=32.9%
  - P. aeruginosa:    6.3 Mb, GC=66.6%
  - K. pneumoniae (plasmid): 131 kb, GC=51.7%

Method:
  1. Fragment each genome into contigs (2-10 kb, random sizes)
  2. Mix all contigs, discard labels
  3. Compute signal vector per contig (GC, 3-mer H, Hurst, ρ-vector, MI, AC)
  4. Cluster with k-means (k=5) and hierarchical
  5. Evaluate: ARI, NMI vs ground truth

Success criteria: ARI > 0.5, NMI > 0.6.
"""

from __future__ import annotations

import numpy as np
from collections import Counter
from pathlib import Path
from dataclasses import dataclass

from Bio import SeqIO
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

DATA_DIR = Path(__file__).parent / "data"
NT_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
ALPHABET_SIZE = 4


def load_seq(path: str) -> np.ndarray:
    rec = SeqIO.read(path, "fasta")
    seq = str(rec.seq).upper()
    arr = np.array([NT_MAP.get(c, 0) for c in seq], dtype=np.uint8)
    arr[arr >= 4] = 0
    return arr


# =========================================================================
# Signal functions (compact, all validated)
# =========================================================================

def entropy(hist: np.ndarray) -> float:
    total = hist.sum()
    if total <= 1:
        return 0.0
    p = hist[hist > 0].astype(np.float64) / total
    return float(-np.sum(p * np.log2(p)))


def gc_content(data: np.ndarray) -> float:
    return float(np.mean((data == 1) | (data == 2)))


def kmer_entropy_3(data: np.ndarray) -> float:
    n = len(data)
    if n < 3:
        return 0.0
    keys = data[:-2].astype(np.int64) * 16 + data[1:-1].astype(np.int64) * 4 + data[2:].astype(np.int64)
    h = entropy(np.bincount(keys, minlength=64))
    return float(h / np.log2(64))


def kmer_entropy_4(data: np.ndarray) -> float:
    n = len(data)
    if n < 4:
        return 0.0
    keys = (data[:-3].astype(np.int64) * 64 + data[1:-2].astype(np.int64) * 16
            + data[2:-1].astype(np.int64) * 4 + data[3:].astype(np.int64))
    h = entropy(np.bincount(keys, minlength=256))
    return float(h / np.log2(256))


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
    A_m = np.vstack([log_w, np.ones(len(log_w))]).T
    H, _ = np.linalg.lstsq(A_m, log_rs, rcond=None)[0]
    return float(np.clip(H, 0.0, 1.0))


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


def autocorr_period3(data: np.ndarray) -> float:
    arr = data.astype(np.float64)
    n = len(arr)
    max_lag = min(30, n // 4)
    if max_lag < 6:
        return 0.0
    centered = arr - arr.mean()
    norm = float(centered @ centered)
    if norm == 0:
        return 0.0
    ac_abs_sum = 0.0
    p3_sum = 0.0
    for lag in range(1, max_lag):
        val = abs(float((centered[:n - lag] @ centered[lag:]) / norm))
        ac_abs_sum += val
        if lag % 3 == 0:
            p3_sum += val
    return float(p3_sum / ac_abs_sum) if ac_abs_sum > 0 else 0.0


def mi_lag1(data: np.ndarray) -> float:
    n = len(data)
    if n < 10:
        return 0.0
    h_m = entropy(np.bincount(data, minlength=4))
    a, b = data[:-1].astype(np.int64), data[1:].astype(np.int64)
    h_j = entropy(np.bincount(a * 4 + b, minlength=16))
    return max(2 * h_m - h_j, 0.0)


# =========================================================================
# Contig simulation
# =========================================================================

@dataclass
class Contig:
    sequence: np.ndarray
    source: str
    source_id: int
    length: int


def fragment_genome(genome: np.ndarray, source_name: str, source_id: int,
                    min_len: int = 2000, max_len: int = 10000,
                    max_contigs: int = 100, rng: np.random.Generator = None) -> list[Contig]:
    """Fragment genome into random-length contigs."""
    if rng is None:
        rng = np.random.default_rng(42)

    contigs = []
    pos = 0
    n = len(genome)
    while pos < n - min_len and len(contigs) < max_contigs:
        length = rng.integers(min_len, max_len + 1)
        end = min(pos + length, n)
        if end - pos >= min_len:
            contigs.append(Contig(
                sequence=genome[pos:end],
                source=source_name,
                source_id=source_id,
                length=end - pos,
            ))
        pos = end
    return contigs


# =========================================================================
# Feature extraction
# =========================================================================

def extract_features(contig: Contig) -> np.ndarray:
    """Extract signal vector for one contig.

    Features (22 total):
      - GC content (1)
      - 3-mer entropy (1)
      - 4-mer entropy (1)
      - Hurst exponent (1)
      - MI(lag=1) (1)
      - Autocorrelation period-3 strength (1)
      - ρ-vector (16 dinucleotide odds ratios)
    """
    data = contig.sequence
    features = [
        gc_content(data),
        kmer_entropy_3(data),
        kmer_entropy_4(data),
        hurst(data),
        mi_lag1(data),
        autocorr_period3(data),
    ]
    features.extend(rho_vector(data).tolist())
    return np.array(features, dtype=np.float64)


# =========================================================================
# Main
# =========================================================================
def main():
    print("=" * 70)
    print("EXPLORATION 006: Reference-Free Metagenomic Clustering")
    print("=" * 70)

    rng = np.random.default_rng(42)

    # Load genomes
    print("\n[1] Loading genomes...")
    genomes = {
        'E_coli': load_seq(str(DATA_DIR / "ecoli_k12.fasta")),
        'B_subtilis': load_seq(str(DATA_DIR / "mock_B_subtilis.fasta")),
        'S_aureus': load_seq(str(DATA_DIR / "mock_S_aureus.fasta")),
        'P_aeruginosa': load_seq(str(DATA_DIR / "mock_P_aeruginosa.fasta")),
    }

    for name, seq in genomes.items():
        gc = gc_content(seq)
        print(f"  {name:15s}: {len(seq):>10,} bp, GC={gc:.1%}")

    # Fragment into contigs
    print("\n[2] Fragmenting genomes into contigs (2-10 kb)...")
    all_contigs: list[Contig] = []
    for i, (name, seq) in enumerate(genomes.items()):
        contigs = fragment_genome(seq, name, i, min_len=2000, max_len=10000,
                                  max_contigs=80, rng=rng)
        all_contigs.extend(contigs)
        print(f"  {name:15s}: {len(contigs)} contigs")

    # Shuffle
    rng.shuffle(all_contigs)
    print(f"\n  Total contigs: {len(all_contigs)}")
    print(f"  Labels hidden. Now clustering blind.")

    # Extract features
    print("\n[3] Extracting signal vectors (22 features per contig)...")
    X = np.array([extract_features(c) for c in all_contigs])
    y_true = np.array([c.source_id for c in all_contigs])
    source_names = {c.source_id: c.source for c in all_contigs}

    print(f"  Feature matrix: {X.shape}")
    print(f"  Feature names: GC, H3, H4, Hurst, MI1, AC3, ρ(AA..TT) × 16")

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # ===================================================================
    # Clustering approaches
    # ===================================================================
    n_species = len(genomes)

    # A. K-means
    print(f"\n" + "=" * 70)
    print(f"[A] K-means (k={n_species})")
    print("=" * 70)

    kmeans = KMeans(n_clusters=n_species, random_state=42, n_init=20)
    y_km = kmeans.fit_predict(X_scaled)

    ari_km = adjusted_rand_score(y_true, y_km)
    nmi_km = normalized_mutual_info_score(y_true, y_km)
    print(f"  ARI = {ari_km:.4f}")
    print(f"  NMI = {nmi_km:.4f}")

    # Confusion-like: for each cluster, what's the majority species?
    print(f"\n  Cluster composition:")
    for k in range(n_species):
        mask = y_km == k
        cluster_sources = Counter(all_contigs[i].source for i in range(len(all_contigs)) if mask[i])
        total = mask.sum()
        majority = cluster_sources.most_common(1)[0]
        purity = majority[1] / total
        comp_str = ", ".join(f"{s}:{c}" for s, c in cluster_sources.most_common())
        print(f"    Cluster {k}: n={total:3d}, majority={majority[0]} ({purity:.0%}), [{comp_str}]")

    # B. Hierarchical (Ward's linkage)
    print(f"\n" + "=" * 70)
    print(f"[B] Hierarchical Clustering (Ward)")
    print("=" * 70)

    dist_matrix = pdist(X_scaled, metric='euclidean')
    Z = linkage(dist_matrix, method='ward')
    y_hc = fcluster(Z, t=n_species, criterion='maxclust')
    y_hc -= 1  # 0-indexed

    ari_hc = adjusted_rand_score(y_true, y_hc)
    nmi_hc = normalized_mutual_info_score(y_true, y_hc)
    print(f"  ARI = {ari_hc:.4f}")
    print(f"  NMI = {nmi_hc:.4f}")

    print(f"\n  Cluster composition:")
    for k in range(n_species):
        mask = y_hc == k
        cluster_sources = Counter(all_contigs[i].source for i in range(len(all_contigs)) if mask[i])
        total = mask.sum()
        if total == 0:
            continue
        majority = cluster_sources.most_common(1)[0]
        purity = majority[1] / total
        comp_str = ", ".join(f"{s}:{c}" for s, c in cluster_sources.most_common())
        print(f"    Cluster {k}: n={total:3d}, majority={majority[0]} ({purity:.0%}), [{comp_str}]")

    # C. Just GC + ρ-vector (minimal features)
    print(f"\n" + "=" * 70)
    print(f"[C] K-means on GC + ρ-vector only (17 features)")
    print("=" * 70)

    # GC = feature 0, ρ = features 6:22
    X_rho = np.hstack([X[:, 0:1], X[:, 6:]])
    X_rho_scaled = StandardScaler().fit_transform(X_rho)
    kmeans_rho = KMeans(n_clusters=n_species, random_state=42, n_init=20)
    y_rho = kmeans_rho.fit_predict(X_rho_scaled)

    ari_rho = adjusted_rand_score(y_true, y_rho)
    nmi_rho = normalized_mutual_info_score(y_true, y_rho)
    print(f"  ARI = {ari_rho:.4f}")
    print(f"  NMI = {nmi_rho:.4f}")

    # D. Just GC alone
    print(f"\n" + "=" * 70)
    print(f"[D] K-means on GC content only (1 feature, baseline)")
    print("=" * 70)

    X_gc = X[:, 0:1]
    kmeans_gc = KMeans(n_clusters=n_species, random_state=42, n_init=20)
    y_gc = kmeans_gc.fit_predict(X_gc)

    ari_gc = adjusted_rand_score(y_true, y_gc)
    nmi_gc = normalized_mutual_info_score(y_true, y_gc)
    print(f"  ARI = {ari_gc:.4f}")
    print(f"  NMI = {nmi_gc:.4f}")

    # ===================================================================
    # Feature importance
    # ===================================================================
    print(f"\n" + "=" * 70)
    print("[E] Feature Importance (leave-one-out ARI drop)")
    print("=" * 70)

    feature_names = ['GC', 'H3', 'H4', 'Hurst', 'MI1', 'AC3'] + \
                    [f'ρ_{chr(65+a//4)}{chr(65+a%4)}' for a in range(16)]

    base_ari = ari_km  # Full feature ARI
    print(f"  Full ARI: {base_ari:.4f}\n")
    drops = []
    for i in range(X_scaled.shape[1]):
        X_loo = np.delete(X_scaled, i, axis=1)
        km_loo = KMeans(n_clusters=n_species, random_state=42, n_init=10)
        y_loo = km_loo.fit_predict(X_loo)
        ari_loo = adjusted_rand_score(y_true, y_loo)
        drop = base_ari - ari_loo
        drops.append((feature_names[i], drop, ari_loo))

    drops.sort(key=lambda x: -x[1])
    print(f"  {'Feature':>10s} {'ARI drop':>10s} {'ARI w/o':>10s}  Impact")
    for name, drop, ari_loo in drops[:10]:
        impact = "HIGH" if drop > 0.05 else ("medium" if drop > 0.01 else "low")
        print(f"  {name:>10s} {drop:+10.4f} {ari_loo:10.4f}  {impact}")

    # ===================================================================
    # Pairwise species confusion
    # ===================================================================
    print(f"\n" + "=" * 70)
    print("[F] Pairwise Species Separability")
    print("=" * 70)

    species_list = sorted(source_names.values())
    species_ids = {name: sid for sid, name in source_names.items()}

    print(f"  {'Species A':>15s} vs {'Species B':>15s}: ARI (k=2, full features)")
    for i in range(len(species_list)):
        for j in range(i + 1, len(species_list)):
            s1, s2 = species_list[i], species_list[j]
            id1, id2 = species_ids[s1], species_ids[s2]
            mask = (y_true == id1) | (y_true == id2)
            if mask.sum() < 10:
                continue
            X_pair = X_scaled[mask]
            y_pair = y_true[mask]
            km_pair = KMeans(n_clusters=2, random_state=42, n_init=10)
            y_pred = km_pair.fit_predict(X_pair)
            ari_pair = adjusted_rand_score(y_pair, y_pred)
            hard = " ← HARD" if ari_pair < 0.5 else ""
            print(f"  {s1:>15s} vs {s2:>15s}: ARI={ari_pair:.4f}{hard}")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print(f"\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    best_ari = max(ari_km, ari_hc)
    best_nmi = max(nmi_km, nmi_hc)
    best_name = "K-means" if ari_km >= ari_hc else "Hierarchical"

    print(f"\n  {'Method':>20s} {'ARI':>8s} {'NMI':>8s}")
    print("  " + "-" * 40)
    print(f"  {'K-means (all 22)':>20s} {ari_km:8.4f} {nmi_km:8.4f}")
    print(f"  {'Hierarchical':>20s} {ari_hc:8.4f} {nmi_hc:8.4f}")
    print(f"  {'GC + ρ only (17)':>20s} {ari_rho:8.4f} {nmi_rho:8.4f}")
    print(f"  {'GC only (1)':>20s} {ari_gc:8.4f} {nmi_gc:8.4f}")

    checks = []

    ok = best_ari > 0.5
    checks.append(("ARI > 0.5", ok))
    print(f"\n  [{'PASS' if ok else 'FAIL'}] Best ARI = {best_ari:.4f} ({best_name}) (target: >0.5)")

    ok = best_nmi > 0.6
    checks.append(("NMI > 0.6", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Best NMI = {best_nmi:.4f} (target: >0.6)")

    ok = best_ari > ari_gc + 0.05
    checks.append(("Full > GC-only by >0.05 ARI", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Full ARI ({best_ari:.4f}) > GC-only ({ari_gc:.4f}) + 0.05")

    ok = ari_rho > ari_gc
    checks.append(("ρ adds value over GC alone", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] GC+ρ ARI ({ari_rho:.4f}) > GC-only ({ari_gc:.4f})")

    passed = sum(1 for _, ok in checks if ok)
    print(f"\n  Result: {passed}/{len(checks)} checks passed")
    print(f"  {'SUCCESS' if passed >= 3 else 'NEEDS INVESTIGATION'}")


if __name__ == "__main__":
    main()
