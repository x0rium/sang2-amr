"""Universal signal functions for SANG2-AMR.

Pure, stateless functions operating on numpy uint8 arrays.
No genomic domain knowledge — works on any integer alphabet.
All thresholds derived from data (SDR / Otsu), zero magic numbers.

Validated on E. coli K-12 in exploration/001-005.
"""

from __future__ import annotations

import numpy as np
from collections import Counter


# ---------------------------------------------------------------------------
# A1. Shannon Entropy
# ---------------------------------------------------------------------------
def entropy(counts) -> float:
    """Shannon entropy H = -sum P(x) log2 P(x) in bits.

    Accepts Counter, dict, ndarray, or any iterable of counts.
    """
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
# GC Content
# ---------------------------------------------------------------------------
def gc_content(data: np.ndarray) -> float:
    """GC fraction for nucleotide array (C=1, G=2)."""
    return float(np.mean((data == 1) | (data == 2)))


# ---------------------------------------------------------------------------
# k-mer Entropy
# ---------------------------------------------------------------------------
def kmer_entropy(data: np.ndarray, k: int = 3, alphabet_size: int = 4) -> float:
    """Normalized entropy of k-mer distribution: H(k-mers) / H_max.

    Validated in 001b: 3-mer entropy separates CDS (0.966) from IG (0.914).
    """
    n = len(data)
    if n < k:
        return 0.0
    keys = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        keys += data[i:n - k + 1 + i].astype(np.int64) * (alphabet_size ** (k - 1 - i))
    hist = np.bincount(keys, minlength=alphabet_size**k)
    h = entropy(hist)
    h_max = np.log2(alphabet_size**k)
    return float(h / h_max) if h_max > 0 else 0.0


# ---------------------------------------------------------------------------
# A8. Hurst Exponent (R/S rescaled-range)
# ---------------------------------------------------------------------------
def hurst(data: np.ndarray) -> float:
    """Hurst exponent via R/S analysis.

    H > 0.5: persistent (structured), H = 0.5: random, H < 0.5: anti-persistent.
    Validated in 001: CDS=0.594, IG=0.624. Requires ≥64 data points.
    """
    arr = data.astype(np.float64)
    n = len(arr)
    if n < 64:
        return 0.5
    windows: list[int] = []
    rs_values: list[float] = []
    w = 8
    while w < n // 4:
        rs_list: list[float] = []
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
    A = np.vstack([log_w, np.ones(len(log_w))]).T
    H, _ = np.linalg.lstsq(A, log_rs, rcond=None)[0]
    return float(np.clip(H, 0.0, 1.0))


# ---------------------------------------------------------------------------
# A5. Mutual Information at lag
# ---------------------------------------------------------------------------
def mi_at_lag(data: np.ndarray, lag: int, alphabet_size: int = 4) -> float:
    """I(X_i; X_{i+lag}) in bits.

    Validated in 002: MI(lag=1) is the strongest chr/plas discriminator (p=5.5e-21).
    """
    n = len(data)
    if n <= lag + 1:
        return 0.0
    marg_hist = np.bincount(data, minlength=alphabet_size)
    h_marg = entropy(marg_hist)
    a = data[:n - lag].astype(np.int64)
    b = data[lag:].astype(np.int64)
    joint_key = a * alphabet_size + b
    joint_hist = np.bincount(joint_key, minlength=alphabet_size**2)
    h_joint = entropy(joint_hist)
    return max(2 * h_marg - h_joint, 0.0)


# ---------------------------------------------------------------------------
# A6. Autocorrelation period-3 strength
# ---------------------------------------------------------------------------
def autocorr_period3(data: np.ndarray, max_lag: int = 30) -> float:
    """Fraction of autocorrelation power at lag multiples of 3.

    Replaces MI half-life (G15). Validated in 001: 97.2% CDS have lag-3k peaks.
    """
    arr = data.astype(np.float64)
    n = len(arr)
    max_lag = min(max_lag, n // 4)
    if max_lag < 6:
        return 0.0
    centered = arr - arr.mean()
    norm = float(centered @ centered)
    if norm == 0:
        return 0.0
    total = 0.0
    period3 = 0.0
    for lag in range(1, max_lag):
        val = abs(float((centered[:n - lag] @ centered[lag:]) / norm))
        total += val
        if lag % 3 == 0:
            period3 += val
    return float(period3 / total) if total > 0 else 0.0


def autocorr_peaks(data: np.ndarray, max_lag: int = 30,
                   threshold: float = 0.005) -> list[int]:
    """Find autocorrelation peak lags above threshold."""
    arr = data.astype(np.float64)
    n = len(arr)
    max_lag = min(max_lag, n // 4)
    if max_lag < 4:
        return []
    centered = arr - arr.mean()
    norm = float(centered @ centered)
    if norm == 0:
        return []
    ac = np.empty(max_lag, dtype=np.float64)
    for lag in range(max_lag):
        ac[lag] = (centered[:n - lag] @ centered[lag:]) / norm
    peaks = []
    for i in range(2, max_lag - 1):
        if ac[i] > ac[i - 1] and ac[i] > ac[i + 1] and ac[i] > threshold:
            peaks.append(i)
    return peaks


# ---------------------------------------------------------------------------
# Dinucleotide odds ratio (ρ-vector)
# ---------------------------------------------------------------------------
def rho_vector(data: np.ndarray, alphabet_size: int = 4) -> np.ndarray:
    """Dinucleotide odds ratio: ρ(XY) = f(XY) / (f(X) * f(Y)).

    Returns alphabet_size^2 element vector. ρ=1.0 = no bias.
    Validated in 002b: 9/16 dinucleotides differ between E. coli chr and K. pneumoniae plasmid.
    """
    n = len(data)
    if n < 10:
        return np.ones(alphabet_size**2)
    mono = np.bincount(data, minlength=alphabet_size).astype(np.float64) / n
    di_keys = data[:-1].astype(np.int64) * alphabet_size + data[1:].astype(np.int64)
    di_freq = np.bincount(di_keys, minlength=alphabet_size**2).astype(np.float64) / (n - 1)
    rho = np.ones(alphabet_size**2)
    for a in range(alphabet_size):
        for b in range(alphabet_size):
            exp = mono[a] * mono[b]
            if exp > 0:
                rho[a * alphabet_size + b] = di_freq[a * alphabet_size + b] / exp
    return rho


def rho_distance(data: np.ndarray, ref_rho: np.ndarray,
                 alphabet_size: int = 4) -> float:
    """Euclidean distance between window ρ-vector and reference ρ-vector."""
    return float(np.sqrt(np.sum((rho_vector(data, alphabet_size) - ref_rho) ** 2)))


# ---------------------------------------------------------------------------
# Chargaff-2 deviation (scalar)
# ---------------------------------------------------------------------------
def chargaff2_dev(data: np.ndarray) -> float:
    """Chargaff's second parity rule deviation: |f(A)-f(T)| + |f(C)-f(G)|.

    Validated in 005: #1 single-signal AMR detector (KPC rank 1/37).
    """
    counts = np.bincount(data, minlength=4).astype(np.float64)
    total = counts.sum()
    if total == 0:
        return 0.0
    freq = counts / total
    return float(abs(freq[0] - freq[3]) + abs(freq[1] - freq[2]))


# ---------------------------------------------------------------------------
# A4. Bidirectional entropy (ΔH)
# ---------------------------------------------------------------------------
def delta_h(data: np.ndarray, k: int = 5, alphabet_size: int = 4) -> np.ndarray:
    """Bidirectional entropy ΔH[i] = H(next|past_k) - H(prev|future_k).

    Negative = convergent (structural). Positive = divergent (content).
    Validated in 001: convergent at promoters (-30..-10 nt before gene start).
    """
    n = len(data)
    if n < 2 * k + 2:
        return np.zeros(n, dtype=np.float64)
    arr = data.astype(np.int64)
    dh = np.zeros(n, dtype=np.float64)

    def hash_ctx(start, length):
        out = np.zeros(length, dtype=np.int64)
        for s in range(k):
            out = out * alphabet_size + arr[start + s:start + s + length]
        return out

    def context_h(contexts, targets):
        if contexts.size == 0:
            return np.zeros(0, dtype=np.float64)
        uniq, inverse = np.unique(contexts, return_inverse=True)
        n_ctx = uniq.size
        flat = inverse.astype(np.int64) * alphabet_size + targets
        joint = np.bincount(flat, minlength=n_ctx * alphabet_size).reshape(
            n_ctx, alphabet_size).astype(np.float64)
        row_sums = np.maximum(joint.sum(axis=1, keepdims=True), 1.0)
        p = joint / row_sums
        with np.errstate(divide="ignore", invalid="ignore"):
            log_p = np.where(p > 0, np.log2(p), 0.0)
        h_per_ctx = -(p * log_p).sum(axis=1)
        return h_per_ctx[inverse]

    fwd_ctx = hash_ctx(0, n - k)
    h_fwd = context_h(fwd_ctx, arr[k:])
    bwd_ctx = hash_ctx(1, n - k)
    h_bwd = context_h(bwd_ctx, arr[:n - k])

    lo, hi = k, n - k
    dh[lo:hi] = h_fwd[:hi - lo] - h_bwd[lo:hi]
    return dh


# ---------------------------------------------------------------------------
# A2. SDR (Structural Detection Rule) on k-mers
# ---------------------------------------------------------------------------
def sdr_kmer(data: np.ndarray, k: int = 4, alphabet_size: int = 4) -> float:
    """SDR ratio: fraction of k-mers that are 'structural' (carry ~50% mass).

    Validated in 001: 91/256 4-mers structural, SDR ratio=35.5%.
    """
    n = len(data)
    if n < k:
        return 0.0
    keys = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        keys += data[i:n - k + 1 + i].astype(np.int64) * (alphabet_size ** (k - 1 - i))
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
    structural = sum(1 for c in freq.values() if c / total > theta)
    return structural / max(len(freq), 1)


# ---------------------------------------------------------------------------
# A3. Otsu Threshold
# ---------------------------------------------------------------------------
def otsu_threshold(values) -> tuple[float, float]:
    """Otsu's method: optimal binary split. Returns (threshold, inter-class variance).

    Zero magic numbers: threshold derived purely from data distribution.
    """
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


# ---------------------------------------------------------------------------
# Codon-pair statistics (resonance-derived, for anomaly scoring)
# ---------------------------------------------------------------------------
def codon_pair_profile(data: np.ndarray) -> np.ndarray:
    """Frequency vector of adjacent codon pairs (64×64 = 4096 elements).

    Computes on reading frame starting at position 0.
    Returns normalized frequency vector (sums to 1).
    """
    n = len(data)
    n_codons = n // 3
    if n_codons < 2:
        return np.zeros(4096, dtype=np.float64)
    codons = data[:n_codons * 3].reshape(n_codons, 3).astype(np.int64)
    codon_idx = codons[:, 0] * 16 + codons[:, 1] * 4 + codons[:, 2]
    pairs = codon_idx[:-1].astype(np.int64) * 64 + codon_idx[1:].astype(np.int64)
    hist = np.bincount(pairs, minlength=4096).astype(np.float64)
    total = hist.sum()
    return hist / total if total > 0 else hist


def codon_pair_distance(data: np.ndarray, ref_profile: np.ndarray) -> float:
    """Jensen-Shannon divergence between window's codon-pair profile and reference.

    JSD is symmetric, bounded [0, 1], and handles zero counts gracefully.
    Higher = more different from host = more anomalous.
    """
    p = codon_pair_profile(data)
    q = ref_profile
    # Smooth to avoid log(0)
    eps = 1e-10
    p = p + eps
    q = q + eps
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)
    with np.errstate(divide="ignore", invalid="ignore"):
        kl_pm = np.sum(np.where(p > 0, p * np.log2(p / m), 0.0))
        kl_qm = np.sum(np.where(q > 0, q * np.log2(q / m), 0.0))
    return float(0.5 * kl_pm + 0.5 * kl_qm)


def codon_mi_profile(data: np.ndarray) -> float:
    """Average PPMI of adjacent codon pairs — measures codon-codon association.

    High MI in CDS = structured codon usage. Foreign genes may have different MI.
    """
    n = len(data)
    n_codons = n // 3
    if n_codons < 10:
        return 0.0
    codons = data[:n_codons * 3].reshape(n_codons, 3).astype(np.int64)
    codon_idx = codons[:, 0] * 16 + codons[:, 1] * 4 + codons[:, 2]

    marg = np.bincount(codon_idx, minlength=64).astype(np.float64)
    h_marg = entropy(marg)

    pairs = codon_idx[:-1].astype(np.int64) * 64 + codon_idx[1:].astype(np.int64)
    h_joint = entropy(np.bincount(pairs, minlength=4096))
    return max(2 * h_marg - h_joint, 0.0)
