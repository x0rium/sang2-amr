"""Feature extraction and clustering for metagenomic binning (L4).

Extracts 22-dimensional signal vector per contig, clusters with
k-means or hierarchical clustering.

Validated in exploration/006: ARI=0.875 on 4-species mock metagenome.
"""

from __future__ import annotations

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

from signals import (
    gc_content, kmer_entropy, hurst, mi_at_lag,
    autocorr_period3, rho_vector,
)


def extract_features(data: np.ndarray) -> np.ndarray:
    """Extract 22-element signal vector for one contig.

    Features:
      [0]    GC content
      [1]    3-mer entropy
      [2]    4-mer entropy
      [3]    Hurst exponent
      [4]    MI(lag=1)
      [5]    Autocorrelation period-3 strength
      [6:22] ρ-vector (16 dinucleotide odds ratios)
    """
    features = [
        gc_content(data),
        kmer_entropy(data, k=3),
        kmer_entropy(data, k=4),
        hurst(data),
        mi_at_lag(data, 1),
        autocorr_period3(data),
    ]
    features.extend(rho_vector(data).tolist())
    return np.array(features, dtype=np.float64)


FEATURE_NAMES = (
    ["GC", "H3", "H4", "Hurst", "MI1", "AC3"]
    + [f"rho_{i:02d}" for i in range(16)]
)


def build_feature_matrix(contigs: list[np.ndarray]) -> np.ndarray:
    """Extract features for multiple contigs. Returns (n_contigs, 22) matrix."""
    return np.array([extract_features(c) for c in contigs])


def cluster_kmeans(X: np.ndarray, k: int, n_init: int = 20) -> np.ndarray:
    """K-means clustering on feature matrix. Returns labels."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    km = KMeans(n_clusters=k, random_state=42, n_init=n_init)
    return km.fit_predict(X_scaled)


def cluster_hierarchical(X: np.ndarray, k: int) -> np.ndarray:
    """Hierarchical clustering (Ward's linkage). Returns labels (0-indexed)."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    dist = pdist(X_scaled, metric="euclidean")
    Z = linkage(dist, method="ward")
    return fcluster(Z, t=k, criterion="maxclust") - 1
