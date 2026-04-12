#!/usr/bin/env python3
"""Exploration 002b: Genomic-specific signals for L0/L1/L3.

Resonance stays universal (6 pair signals). Genomic domain knowledge goes here.

NEW L0 signals (window-level):
  S1. GC-skew          — (G-C)/(G+C), cumulative → replication origin
  S2. Chargaff-2 dev   — strand symmetry violation → ssDNA origin, phage
  S3. Dinucleotide ρ   — obs/exp per dinucleotide → CpG suppression, HGT
  S4. Codon Usage Bias — deviation from host → HGT marker

NEW L1 signals (structure-level):
  S5. Inverted repeat density — IS-element flanks, terminators

NEW L3 signals (anomaly-level):
  S6. CUB anomaly score — per-window deviation from host codon usage

Test on: E. coli K-12 chromosome + pCAV1392-131 plasmid.
"""

from __future__ import annotations

import numpy as np
from collections import Counter
from pathlib import Path
from scipy import stats

from Bio import SeqIO

DATA_DIR = Path(__file__).parent / "data"
ALPHABET_SIZE = 4

NT_MAP = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3,
          ord('a'): 0, ord('c'): 1, ord('g'): 2, ord('t'): 3}
NT_CHARS = "ACGT"
COMPLEMENT = np.array([3, 2, 1, 0], dtype=np.uint8)  # A↔T, C↔G


def load_seq(fasta_path: str) -> np.ndarray:
    rec = SeqIO.read(fasta_path, "fasta")
    seq = str(rec.seq).upper()
    arr = np.array([NT_MAP.get(ord(c), 255) for c in seq], dtype=np.uint8)
    arr[arr >= 4] = np.random.randint(0, 4, size=int((arr >= 4).sum())).astype(np.uint8)
    return arr


def load_cds(gb_path: str):
    from dataclasses import dataclass

    @dataclass
    class Region:
        start: int
        end: int
        strand: int
        name: str
        kind: str

    rec = SeqIO.read(gb_path, "genbank")
    regions = []
    for feat in rec.features:
        if feat.type in ('CDS', 'gene'):
            loc = feat.location
            if loc is None:
                continue
            product = feat.qualifiers.get('product', [''])[0]
            name = feat.qualifiers.get('gene', [product])[0]
            regions.append(Region(
                start=int(loc.start), end=int(loc.end),
                strand=int(loc.strand) if loc.strand else 1,
                name=name, kind=feat.type,
            ))
    return regions


# =========================================================================
# S1. GC-skew
# =========================================================================
def gc_skew(data: np.ndarray, window: int = 1000) -> tuple[np.ndarray, np.ndarray]:
    """Sliding GC-skew: (G-C)/(G+C) per window.

    Returns (positions, skew_values).
    Sign change = replication origin/terminus.
    """
    n = len(data)
    n_windows = n // window
    positions = np.arange(n_windows) * window + window // 2
    skew = np.empty(n_windows)
    for i in range(n_windows):
        seg = data[i * window:(i + 1) * window]
        g = np.sum(seg == 2)
        c = np.sum(seg == 1)
        denom = g + c
        skew[i] = (g - c) / denom if denom > 0 else 0.0
    return positions, skew


def gc_skew_cumulative(data: np.ndarray, window: int = 1000) -> tuple[np.ndarray, np.ndarray]:
    """Cumulative GC-skew. Minimum = origin, maximum = terminus."""
    positions, skew = gc_skew(data, window)
    return positions, np.cumsum(skew)


# =========================================================================
# S2. Chargaff-2 deviation
# =========================================================================
def chargaff2_deviation(data: np.ndarray, window: int = 5000) -> tuple[np.ndarray, np.ndarray]:
    """Chargaff's second parity rule deviation per window.

    For each window: sum of |freq(X) - freq(complement(X))| for X in {A,C}.
    Perfect Chargaff-2 → 0. ssDNA regions, phage → high deviation.
    """
    n = len(data)
    step = window // 2
    n_windows = max(1, (n - window) // step + 1)
    positions = np.empty(n_windows)
    deviations = np.empty(n_windows)

    for i in range(n_windows):
        start = i * step
        seg = data[start:start + window]
        counts = np.bincount(seg, minlength=4).astype(np.float64)
        total = counts.sum()
        if total == 0:
            deviations[i] = 0.0
        else:
            freq = counts / total
            # |freq(A) - freq(T)| + |freq(C) - freq(G)|
            deviations[i] = abs(freq[0] - freq[3]) + abs(freq[1] - freq[2])
        positions[i] = start + window / 2

    return positions, deviations


# =========================================================================
# S3. Dinucleotide odds ratio (ρ)
# =========================================================================
def dinucleotide_rho(data: np.ndarray) -> dict[str, float]:
    """Dinucleotide odds ratio: ρ(XY) = f(XY) / (f(X) * f(Y)).

    ρ = 1.0: no bias.
    ρ < 1.0: suppressed (CpG in vertebrates, CG in many bacteria).
    ρ > 1.0: over-represented.

    Classic HGT marker: genome-wide ρ profile differs between species.
    """
    n = len(data)
    if n < 2:
        return {}

    mono = np.bincount(data, minlength=4).astype(np.float64) / n
    di_keys = data[:-1].astype(np.int64) * 4 + data[1:].astype(np.int64)
    di_counts = np.bincount(di_keys, minlength=16).astype(np.float64)
    di_freq = di_counts / (n - 1)

    rho = {}
    for a in range(4):
        for b in range(4):
            label = NT_CHARS[a] + NT_CHARS[b]
            expected = mono[a] * mono[b]
            if expected > 0:
                rho[label] = float(di_freq[a * 4 + b] / expected)
            else:
                rho[label] = 0.0
    return rho


def rho_vector(data: np.ndarray) -> np.ndarray:
    """16-element ρ vector for a sequence. For distance computation."""
    rho = dinucleotide_rho(data)
    return np.array([rho.get(NT_CHARS[a] + NT_CHARS[b], 1.0)
                     for a in range(4) for b in range(4)])


def rho_distance(data: np.ndarray, reference_rho: np.ndarray) -> float:
    """Euclidean distance between window ρ-vector and reference ρ-vector."""
    window_rho = rho_vector(data)
    return float(np.sqrt(np.sum((window_rho - reference_rho) ** 2)))


# =========================================================================
# S4. Codon Usage Bias (RSCU-based)
# =========================================================================
CODON_TABLE = {
    # AA: list of codon indices (base-4 triplet → 0..63)
    # Simplified: group by amino acid
}

def _build_synonymous_groups():
    """Build codon → amino acid mapping, return groups of synonymous codons."""
    from Bio.Data.CodonTable import standard_dna_table
    groups = {}
    for codon, aa in standard_dna_table.forward_table.items():
        idx = sum(NT_MAP[ord(c)] * (4 ** (2 - i)) for i, c in enumerate(codon))
        groups.setdefault(aa, []).append(idx)
    # Add stop codons
    for codon in standard_dna_table.stop_codons:
        idx = sum(NT_MAP[ord(c)] * (4 ** (2 - i)) for i, c in enumerate(codon))
        groups.setdefault('*', []).append(idx)
    return groups

SYNONYMOUS_GROUPS = _build_synonymous_groups()


def rscu(data: np.ndarray) -> np.ndarray:
    """Relative Synonymous Codon Usage — 64-element vector.

    RSCU(codon) = observed(codon) / expected_if_uniform_within_AA_group.
    RSCU = 1.0 for all → no bias. Deviation = codon usage bias.
    """
    n = len(data)
    n_codons = n // 3
    if n_codons < 10:
        return np.ones(64)

    codons = data[:n_codons * 3].reshape(n_codons, 3).astype(np.int64)
    codon_idx = codons[:, 0] * 16 + codons[:, 1] * 4 + codons[:, 2]
    counts = np.bincount(codon_idx, minlength=64).astype(np.float64)

    rscu_vec = np.ones(64)
    for aa, group in SYNONYMOUS_GROUPS.items():
        group_count = sum(counts[c] for c in group)
        n_syn = len(group)
        if group_count > 0:
            expected = group_count / n_syn
            for c in group:
                rscu_vec[c] = counts[c] / expected if expected > 0 else 1.0
    return rscu_vec


def codon_bias_distance(data: np.ndarray, ref_rscu: np.ndarray) -> float:
    """Euclidean distance between window RSCU and reference RSCU."""
    w_rscu = rscu(data)
    return float(np.sqrt(np.sum((w_rscu - ref_rscu) ** 2)))


# =========================================================================
# S5. Inverted Repeat Density
# =========================================================================
def reverse_complement(data: np.ndarray) -> np.ndarray:
    """Reverse complement of a nucleotide array."""
    return COMPLEMENT[data[::-1]]


def inverted_repeat_density(data: np.ndarray, min_len: int = 8, max_len: int = 30,
                            max_gap: int = 500, window: int = 5000, step: int = 1000) -> tuple[np.ndarray, np.ndarray]:
    """Count inverted repeats (IR) per window.

    An IR is: seq[i:i+L] == reverse_complement(seq[j:j+L]) where j > i+L
    and j - (i+L) <= max_gap.

    Uses k-mer hashing for efficiency. Returns (positions, ir_counts).
    """
    n = len(data)
    n_windows = max(1, (n - window) // step + 1)
    positions = np.empty(n_windows)
    ir_counts = np.empty(n_windows)

    k = min_len  # search for exact matches of length k as seeds

    # Precompute forward k-mer hashes
    fwd_keys = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        fwd_keys += data[i:n - k + 1 + i].astype(np.int64) * (4 ** (k - 1 - i))

    # Precompute reverse complement k-mer hashes
    rc_data = COMPLEMENT[data]
    rc_keys = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        rc_keys += rc_data[i:n - k + 1 + i].astype(np.int64) * (4 ** (k - 1 - i))
    # Reverse the rc_keys positions: rc at position j means the RC of seq[j:j+k]
    # We want: fwd_keys[i] == rc_keys_reversed[j] where rc of seq[j:j+k] read backwards
    # Actually: IR means seq[i:i+k] == revcomp(seq[j:j+k])
    # revcomp(seq[j:j+k]) as a forward read = complement(seq[j+k-1]) ... complement(seq[j])
    # So we need to hash the reverse complement starting from each position
    # Let's recompute properly:
    rc_seq = reverse_complement(data)  # full RC of genome
    # IR: seq[i:i+k] matches rc_seq[n-j-k : n-j] = rc_seq at position (n-j-k)
    # Easier: for each position i, hash seq[i:i+k]. For each position j, hash revcomp(seq[j:j+k]).
    # revcomp(seq[j:j+k]) = complement(seq[j+k-1:j-1:-1])
    # Let's just hash the reversed complement sequences

    # Build RC hashes: for position j, hash of revcomp(seq[j:j+k])
    rc_hashes = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        rc_hashes += COMPLEMENT[data[k - 1 - i:n - i]].astype(np.int64) * (4 ** (k - 1 - i))

    for wi in range(n_windows):
        w_start = wi * step
        w_end = min(w_start + window, n - k + 1)
        count = 0

        # Get forward hashes in this window
        fwd_slice = fwd_keys[w_start:w_end]
        # Build set for O(1) lookup
        fwd_set = set(fwd_slice.tolist())

        # Check RC hashes in nearby region (window + max_gap)
        search_end = min(w_end + max_gap, n - k + 1)
        rc_slice = rc_hashes[w_start:search_end]

        # Count matches (this is approximate — counts k-mer IR seeds)
        rc_set = set(rc_slice.tolist())
        count = len(fwd_set & rc_set)

        positions[wi] = w_start + window / 2
        ir_counts[wi] = count

    return positions, ir_counts


# =========================================================================
# S6. CUB Anomaly Score (L3)
# =========================================================================
def cub_anomaly_profile(data: np.ndarray, host_rscu: np.ndarray,
                        window: int = 3000, step: int = 500) -> tuple[np.ndarray, np.ndarray]:
    """Per-window CUB anomaly: Euclidean distance from host RSCU.

    High score = foreign codon usage → HGT candidate.
    """
    n = len(data)
    n_windows = max(1, (n - window) // step + 1)
    positions = np.empty(n_windows)
    scores = np.empty(n_windows)

    for i in range(n_windows):
        start = i * step
        seg = data[start:start + window]
        scores[i] = codon_bias_distance(seg, host_rscu)
        positions[i] = start + window / 2

    return positions, scores


# =========================================================================
# Main experiment
# =========================================================================
def main():
    print("=" * 70)
    print("EXPLORATION 002b: Genomic-Specific Signals (L0/L1/L3)")
    print("=" * 70)

    np.random.seed(42)

    chr_seq = load_seq(str(DATA_DIR / "ecoli_k12.fasta"))
    plas_seq = load_seq(str(DATA_DIR / "pkpc_cav1321.fasta"))
    plas_feats = load_cds(str(DATA_DIR / "pkpc_cav1321.gb"))

    print(f"  Chromosome: {len(chr_seq):,} nt")
    print(f"  Plasmid:    {len(plas_seq):,} nt")

    tnp_regions = [(f.start, f.end, f.name) for f in plas_feats
                   if 'transposase' in f.name.lower()]
    tra_regions = [(f.start, f.end) for f in plas_feats
                   if 'conjugal' in f.name.lower()]
    tra_start = min(s for s, e in tra_regions) if tra_regions else 0
    tra_end = max(e for s, e in tra_regions) if tra_regions else 0

    # ===================================================================
    # S1. GC-skew
    # ===================================================================
    print("\n" + "=" * 70)
    print("[S1] GC-skew")
    print("=" * 70)

    pos_chr, cumskew_chr = gc_skew_cumulative(chr_seq, window=5000)
    min_idx = np.argmin(cumskew_chr)
    max_idx = np.argmax(cumskew_chr)
    print(f"  Chromosome GC-skew (cumulative):")
    print(f"    Minimum (≈ oriC) at position: {int(pos_chr[min_idx]):,}")
    print(f"    Maximum (≈ terC) at position: {int(pos_chr[max_idx]):,}")
    print(f"    Known oriC ≈ 3,925,744; known terC ≈ 1,588,774")
    orc_err = abs(int(pos_chr[min_idx]) - 3_925_744)
    ter_err = abs(int(pos_chr[max_idx]) - 1_588_774)
    print(f"    oriC error: {orc_err:,} bp ({orc_err/len(chr_seq)*100:.1f}%)")
    print(f"    terC error: {ter_err:,} bp ({ter_err/len(chr_seq)*100:.1f}%)")

    pos_plas, cumskew_plas = gc_skew_cumulative(plas_seq, window=2000)
    print(f"\n  Plasmid GC-skew (cumulative):")
    print(f"    Range: [{cumskew_plas.min():.2f}, {cumskew_plas.max():.2f}]")
    print(f"    Std: {cumskew_plas.std():.3f} (chromosome: {cumskew_chr.std():.3f})")

    # ===================================================================
    # S2. Chargaff-2 deviation
    # ===================================================================
    print("\n" + "=" * 70)
    print("[S2] Chargaff-2 Deviation")
    print("=" * 70)

    pos_c2_chr, dev_chr = chargaff2_deviation(chr_seq, window=10000)
    pos_c2_plas, dev_plas = chargaff2_deviation(plas_seq, window=5000)

    print(f"  Chromosome (10kb windows): mean={dev_chr.mean():.5f}, std={dev_chr.std():.5f}")
    print(f"  Plasmid (5kb windows):     mean={dev_plas.mean():.5f}, std={dev_plas.std():.5f}")

    # Chargaff-2 in conjugal vs rest of plasmid
    if tra_regions:
        tra_mask = (pos_c2_plas >= tra_start) & (pos_c2_plas <= tra_end)
        if tra_mask.sum() > 0 and (~tra_mask).sum() > 0:
            _, p = stats.mannwhitneyu(dev_plas[tra_mask], dev_plas[~tra_mask])
            print(f"  Conjugal region: {dev_plas[tra_mask].mean():.5f} vs rest: {dev_plas[~tra_mask].mean():.5f} (p={p:.2e})")

    # ===================================================================
    # S3. Dinucleotide odds ratio (ρ)
    # ===================================================================
    print("\n" + "=" * 70)
    print("[S3] Dinucleotide Odds Ratio (ρ)")
    print("=" * 70)

    rho_chr = dinucleotide_rho(chr_seq)
    rho_plas = dinucleotide_rho(plas_seq)

    print(f"  {'Dinuc':>5s} {'Chr ρ':>8s} {'Plas ρ':>8s} {'Δ':>8s}  Note")
    print("  " + "-" * 55)
    for dn in ['AA', 'TT', 'AT', 'TA', 'CG', 'GC', 'CC', 'GG', 'AC', 'CA', 'AG', 'GA', 'TC', 'CT', 'TG', 'GT']:
        c, p = rho_chr[dn], rho_plas[dn]
        diff = p - c
        note = ""
        if dn == 'CG' and c < 0.85:
            note = "← CpG suppressed"
        elif dn == 'TA' and c < 0.85:
            note = "← TA suppressed"
        elif abs(diff) > 0.05:
            note = "← large Δ"
        print(f"  {dn:>5s} {c:8.4f} {p:8.4f} {diff:+8.4f}  {note}")

    # ρ-distance for sliding windows along plasmid
    chr_rho_vec = rho_vector(chr_seq)
    rho_window = 5000
    rho_step = 1000
    n_rho_w = (len(plas_seq) - rho_window) // rho_step + 1
    rho_distances = np.empty(n_rho_w)
    rho_positions = np.empty(n_rho_w)
    for i in range(n_rho_w):
        start = i * rho_step
        seg = plas_seq[start:start + rho_window]
        rho_distances[i] = rho_distance(seg, chr_rho_vec)
        rho_positions[i] = start + rho_window / 2

    print(f"\n  ρ-distance from chromosome (sliding 5kb along plasmid):")
    print(f"    mean = {rho_distances.mean():.4f}, std = {rho_distances.std():.4f}")
    print(f"    max  = {rho_distances.max():.4f} at pos {int(rho_positions[np.argmax(rho_distances)]):,}")

    # Check ρ near transposases
    tnp_rho_dists = []
    rest_rho_dists = []
    for i, pos in enumerate(rho_positions):
        near_tnp = any(abs(pos - ts) < 3000 or abs(pos - te) < 3000 for ts, te, _ in tnp_regions)
        if near_tnp:
            tnp_rho_dists.append(rho_distances[i])
        else:
            rest_rho_dists.append(rho_distances[i])

    if tnp_rho_dists and rest_rho_dists:
        _, p = stats.mannwhitneyu(tnp_rho_dists, rest_rho_dists, alternative='greater')
        print(f"    Near transposases: {np.mean(tnp_rho_dists):.4f} vs rest: {np.mean(rest_rho_dists):.4f} (p={p:.2e})")

    # ===================================================================
    # S4. Codon Usage Bias (RSCU)
    # ===================================================================
    print("\n" + "=" * 70)
    print("[S4] Codon Usage Bias (RSCU)")
    print("=" * 70)

    # Host baseline: E. coli CDS regions
    chr_feats = load_cds(str(DATA_DIR / "ecoli_k12.gb"))
    chr_cds = [f for f in chr_feats if f.kind == 'CDS']

    # Concatenate first 500 CDS for host RSCU
    host_codons = []
    for f in chr_cds[:500]:
        seg = chr_seq[f.start:f.end]
        if len(seg) >= 30:
            host_codons.append(seg)
    host_concat = np.concatenate(host_codons)
    host_rscu = rscu(host_concat)

    print(f"  Host RSCU computed from {len(host_codons)} E. coli CDS")

    # Plasmid global RSCU
    plas_rscu = rscu(plas_seq)
    rscu_dist_global = float(np.sqrt(np.sum((plas_rscu - host_rscu) ** 2)))
    print(f"  Plasmid global RSCU distance from host: {rscu_dist_global:.3f}")

    # Per-CDS RSCU distance on plasmid
    plas_cds = [f for f in plas_feats if f.kind == 'CDS']
    cds_dists = []
    tnp_dists = []
    tra_dists = []
    other_dists = []
    for f in plas_cds:
        seg = plas_seq[f.start:f.end]
        if len(seg) < 60:
            continue
        d = codon_bias_distance(seg, host_rscu)
        cds_dists.append((f.name, f.start, d))
        if 'transposase' in f.name.lower():
            tnp_dists.append(d)
        elif 'conjugal' in f.name.lower() or f.name.lower().startswith('tra'):
            tra_dists.append(d)
        else:
            other_dists.append(d)

    print(f"\n  Per-CDS RSCU distance from E. coli host:")
    print(f"    All CDS (n={len(cds_dists)}):       mean={np.mean([d for _,_,d in cds_dists]):.3f}")
    if tnp_dists:
        print(f"    Transposases (n={len(tnp_dists)}):   mean={np.mean(tnp_dists):.3f}")
    if tra_dists:
        print(f"    Conjugal (n={len(tra_dists)}):       mean={np.mean(tra_dists):.3f}")
    if other_dists:
        print(f"    Other (n={len(other_dists)}):         mean={np.mean(other_dists):.3f}")

    # Top CUB deviants
    cds_dists.sort(key=lambda x: -x[2])
    print(f"\n  Top 10 CUB deviants (furthest from E. coli codon usage):")
    for name, start, d in cds_dists[:10]:
        is_tnp = "TNP" if 'transposase' in name.lower() else "   "
        print(f"    {is_tnp} {start:>6}: d={d:.3f}  {name[:40]}")

    # ===================================================================
    # S5. Inverted Repeat Density
    # ===================================================================
    print("\n" + "=" * 70)
    print("[S5] Inverted Repeat Density")
    print("=" * 70)

    ir_pos, ir_counts = inverted_repeat_density(
        plas_seq, min_len=8, max_gap=500, window=5000, step=2000)

    print(f"  Plasmid IR density (8-mer seeds, gap ≤ 500):")
    print(f"    mean = {ir_counts.mean():.1f}, std = {ir_counts.std():.1f}")
    print(f"    max  = {ir_counts.max():.0f} at pos {int(ir_pos[np.argmax(ir_counts)]):,}")

    # Compare near transposases vs elsewhere
    tnp_ir = []
    rest_ir = []
    for i, pos in enumerate(ir_pos):
        near = any(abs(pos - ts) < 3000 or abs(pos - te) < 3000 for ts, te, _ in tnp_regions)
        if near:
            tnp_ir.append(ir_counts[i])
        else:
            rest_ir.append(ir_counts[i])

    if tnp_ir and rest_ir:
        _, p = stats.mannwhitneyu(tnp_ir, rest_ir, alternative='greater')
        print(f"    Near transposases: {np.mean(tnp_ir):.1f} vs rest: {np.mean(rest_ir):.1f} (p={p:.2e})")

    # Chromosome baseline
    ir_pos_chr, ir_counts_chr = inverted_repeat_density(
        chr_seq[:500_000], min_len=8, max_gap=500, window=5000, step=2000)
    print(f"  Chromosome IR density (first 500kb): mean={ir_counts_chr.mean():.1f}, std={ir_counts_chr.std():.1f}")

    # ===================================================================
    # S6. CUB Anomaly Profile (L3)
    # ===================================================================
    print("\n" + "=" * 70)
    print("[S6] CUB Anomaly Profile (L3 — HGT Detection)")
    print("=" * 70)

    # Use K. pneumoniae (plasmid) CDS as its own host, compare to E. coli
    cub_pos, cub_scores = cub_anomaly_profile(plas_seq, host_rscu, window=3000, step=500)
    print(f"  CUB anomaly (distance from E. coli host, 3kb windows):")
    print(f"    mean = {cub_scores.mean():.3f}, std = {cub_scores.std():.3f}")

    # Anomaly threshold via Otsu
    vals = sorted(cub_scores)
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
    cub_threshold = best_t
    anomalous = np.sum(cub_scores > cub_threshold)
    print(f"    Otsu threshold: {cub_threshold:.3f}")
    print(f"    Anomalous windows: {anomalous}/{len(cub_scores)} ({anomalous/len(cub_scores)*100:.1f}%)")

    # Check overlap with transposases
    tnp_cub = []
    rest_cub = []
    for i, pos in enumerate(cub_pos):
        near = any(abs(pos - ts) < 3000 or abs(pos - te) < 3000 for ts, te, _ in tnp_regions)
        if near:
            tnp_cub.append(cub_scores[i])
        else:
            rest_cub.append(cub_scores[i])

    if tnp_cub and rest_cub:
        _, p = stats.mannwhitneyu(tnp_cub, rest_cub, alternative='greater')
        print(f"    Near transposases: {np.mean(tnp_cub):.3f} vs rest: {np.mean(rest_cub):.3f} (p={p:.2e})")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print("\n" + "=" * 70)
    print("SUMMARY — Genomic-Specific Signals")
    print("=" * 70)

    checks = []

    # S1: GC-skew finds oriC
    orc_pct = orc_err / len(chr_seq) * 100
    ok = orc_pct < 5.0
    checks.append(("S1 GC-skew: oriC within 5%", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] GC-skew oriC error: {orc_pct:.1f}%")

    # S2: Chargaff-2 deviation differs chr vs plas
    _, p_c2 = stats.mannwhitneyu(dev_chr, dev_plas, alternative='two-sided')
    ok = p_c2 < 0.01
    checks.append(("S2 Chargaff-2: chr ≠ plas", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Chargaff-2: chr={dev_chr.mean():.5f} vs plas={dev_plas.mean():.5f} (p={p_c2:.2e})")

    # S3: ρ near transposases > rest
    if tnp_rho_dists and rest_rho_dists:
        ok = np.mean(tnp_rho_dists) > np.mean(rest_rho_dists)
        checks.append(("S3 ρ-distance: transposases > rest", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] ρ near TNP: {np.mean(tnp_rho_dists):.4f} vs rest: {np.mean(rest_rho_dists):.4f}")

    # S4: CUB separates gene types
    if tnp_dists and tra_dists:
        ok = np.mean(tnp_dists) > np.mean(other_dists) or np.mean(tra_dists) > np.mean(other_dists)
        checks.append(("S4 CUB: TNP or TRA differ from other", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] CUB: TNP={np.mean(tnp_dists):.3f}, TRA={np.mean(tra_dists):.3f}, other={np.mean(other_dists):.3f}")

    # S5: IR density higher near transposases
    if tnp_ir and rest_ir:
        ok = np.mean(tnp_ir) > np.mean(rest_ir)
        checks.append(("S5 IR density: higher near TNP", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] IR near TNP: {np.mean(tnp_ir):.1f} vs rest: {np.mean(rest_ir):.1f}")

    # S6: CUB anomaly detects near transposases
    if tnp_cub and rest_cub:
        ok = np.mean(tnp_cub) > np.mean(rest_cub)
        checks.append(("S6 CUB anomaly: higher near TNP", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] CUB anomaly near TNP: {np.mean(tnp_cub):.3f} vs rest: {np.mean(rest_cub):.3f}")

    passed = sum(1 for _, ok in checks if ok)
    print(f"\n  Result: {passed}/{len(checks)} checks passed")
    print(f"  {'SUCCESS' if passed >= 4 else 'PARTIAL — review failures'}")


if __name__ == "__main__":
    main()
