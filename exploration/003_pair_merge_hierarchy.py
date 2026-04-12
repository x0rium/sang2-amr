#!/usr/bin/env python3
"""Exploration 003: PPMI Pair Merge on nucleotide sequences.

Algorithm B1 from SANG2: iterative merging of most-associated pairs.
Like BPE, but scored by PPMI × log2(1 + count) instead of pure frequency.

Expected hierarchy for genomic data:
  Level 0: A, C, G, T                    (4 symbols)
  Level 1: Dinucleotides (CG, AT, ...)   (high PPMI pairs)
  Level 2: Codons (ATG, TAA, ...)        (functional units)
  Level 3: Motifs (TATAAT, AGGAGG, ...)   (regulatory elements)

Success criteria:
  - Lossless reversibility (reconstruct = original)
  - Level 2 contains ≥50% of start codons (ATG)
  - SDR on score distribution gives meaningful threshold

Input: E. coli K-12, first 100 kb.
"""

from __future__ import annotations

import numpy as np
from collections import Counter
from pathlib import Path
import time

from Bio import SeqIO

DATA_DIR = Path(__file__).parent / "data"
NT_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
NT_CHARS = ['A', 'C', 'G', 'T']


def load_genome(limit: int = 100_000) -> list[int]:
    """Load first `limit` nt of E. coli K-12 as list of ints."""
    rec = SeqIO.read(str(DATA_DIR / "ecoli_k12.fasta"), "fasta")
    seq = str(rec.seq).upper()[:limit]
    return [NT_MAP.get(c, 0) for c in seq]


# =========================================================================
# PPMI Pair Merge (SANG2 B1)
# =========================================================================

def count_pairs(sequence: list[int]) -> Counter:
    """Count all adjacent pairs in sequence."""
    return Counter(zip(sequence[:-1], sequence[1:]))


def count_symbols(sequence: list[int]) -> Counter:
    """Count all symbols in sequence."""
    return Counter(sequence)


def compute_ppmi(pair_counts: Counter, sym_counts: Counter, total_pairs: int) -> dict:
    """Compute PPMI for all pairs.

    PMI(a,b) = log2(P(a,b) / (P(a) * P(b)))
    PPMI = max(PMI, 0)
    """
    total_syms = sum(sym_counts.values())
    ppmi = {}
    for (a, b), count in pair_counts.items():
        p_ab = count / total_pairs
        p_a = sym_counts[a] / total_syms
        p_b = sym_counts[b] / total_syms
        if p_a > 0 and p_b > 0 and p_ab > 0:
            pmi = np.log2(p_ab / (p_a * p_b))
            ppmi[(a, b)] = max(pmi, 0.0)
        else:
            ppmi[(a, b)] = 0.0
    return ppmi


def score_pairs(pair_counts: Counter, ppmi: dict) -> dict:
    """Score pairs: PPMI × log2(1 + count)."""
    scores = {}
    for pair, count in pair_counts.items():
        scores[pair] = ppmi.get(pair, 0.0) * np.log2(1 + count)
    return scores


def merge_pair(sequence: list[int], pair: tuple[int, int], new_symbol: int) -> list[int]:
    """Replace all occurrences of `pair` in sequence with `new_symbol`.

    Left-to-right greedy: if pair = (A,B) and sequence has A,B,B,
    merges the first A,B → new, leaving B.
    """
    result = []
    i = 0
    while i < len(sequence):
        if i < len(sequence) - 1 and sequence[i] == pair[0] and sequence[i + 1] == pair[1]:
            result.append(new_symbol)
            i += 2
        else:
            result.append(sequence[i])
            i += 1
    return result


def sdr_threshold(scores: dict) -> float:
    """SDR: find θ where symbols with score > θ carry ~50% of total mass."""
    if not scores:
        return 0.0
    values = list(scores.values())
    total = sum(values)
    if total == 0:
        return 0.0
    lo, hi = 0.0, max(values)
    for _ in range(50):
        mid = (lo + hi) / 2
        mass = sum(v for v in values if v > mid)
        if mass / total > 0.50:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


class MergeHistory:
    """Track the merge hierarchy for lossless reversal."""

    def __init__(self, base_alphabet: list[str]):
        self.symbol_names: dict[int, str] = {i: name for i, name in enumerate(base_alphabet)}
        self.merges: list[tuple[int, int, int, float]] = []  # (a, b, new, score)
        self.next_id = len(base_alphabet)

    def record_merge(self, a: int, b: int, score: float) -> int:
        new_id = self.next_id
        self.next_id += 1
        name_a = self.symbol_names[a]
        name_b = self.symbol_names[b]
        self.symbol_names[new_id] = name_a + name_b
        self.merges.append((a, b, new_id, score))
        return new_id

    def decode(self, sequence: list[int]) -> str:
        """Decode merged sequence back to original nucleotides."""
        return ''.join(self.symbol_names[s] for s in sequence)

    def symbol_length(self, sym_id: int) -> int:
        """Length of symbol in base nucleotides."""
        return len(self.symbol_names.get(sym_id, ''))

    def level_of(self, sym_id: int) -> int:
        """Merge level: 0 for base, 1 for first merge, etc."""
        length = self.symbol_length(sym_id)
        if length <= 1:
            return 0
        return length - 1  # Approximate: each merge adds 1 nt


def ppmi_pair_merge(sequence: list[int], max_merges: int = 200,
                    min_score: float = 0.01, verbose: bool = True):
    """Run PPMI pair merge iteratively.

    Returns (final_sequence, history, level_stats).
    """
    history = MergeHistory(NT_CHARS)
    seq = list(sequence)

    level_stats = []  # (merge_idx, pair_name, score, seq_len, vocab_size)

    for step in range(max_merges):
        pair_counts = count_pairs(seq)
        if not pair_counts:
            break

        sym_counts = count_symbols(seq)
        total_pairs = sum(pair_counts.values())
        ppmi = compute_ppmi(pair_counts, sym_counts, total_pairs)
        scores = score_pairs(pair_counts, ppmi)

        if not scores:
            break

        # Find best pair
        best_pair = max(scores, key=scores.get)
        best_score = scores[best_pair]

        if best_score < min_score:
            break

        # Merge
        new_id = history.record_merge(best_pair[0], best_pair[1], best_score)
        seq = merge_pair(seq, best_pair, new_id)

        name = history.symbol_names[new_id]
        vocab = len(set(seq))
        level_stats.append((step, name, best_score, len(seq), vocab))

        if verbose and step < 30:
            count = pair_counts[best_pair]
            print(f"  Step {step:3d}: merge ({history.symbol_names[best_pair[0]]}, "
                  f"{history.symbol_names[best_pair[1]]}) → {name:10s} "
                  f"score={best_score:.3f} count={count:5d} → len={len(seq):,} vocab={vocab}")
        elif verbose and step == 30:
            print(f"  ... (continuing silently)")

    return seq, history, level_stats


# =========================================================================
# Analysis
# =========================================================================

def main():
    print("=" * 70)
    print("EXPLORATION 003: PPMI Pair Merge Hierarchy")
    print("=" * 70)

    # Load data
    print("\n[1] Loading E. coli K-12 (first 100 kb)...")
    seq_int = load_genome(100_000)
    original_str = ''.join(NT_CHARS[s] for s in seq_int)
    print(f"  Length: {len(seq_int):,} nt")
    print(f"  Alphabet: {sorted(set(seq_int))}")

    # ===================================================================
    # Run pair merge
    # ===================================================================
    print("\n" + "=" * 70)
    print("[2] Running PPMI Pair Merge")
    print("=" * 70)
    print(f"  Max merges: 200, min_score: 0.01\n")

    t0 = time.time()
    final_seq, history, stats = ppmi_pair_merge(seq_int, max_merges=200, min_score=0.01)
    elapsed = time.time() - t0

    print(f"\n  Completed in {elapsed:.1f}s")
    print(f"  Merges: {len(stats)}")
    print(f"  Final sequence: {len(final_seq):,} tokens (from {len(seq_int):,} nt)")
    print(f"  Compression ratio: {len(seq_int)/len(final_seq):.2f}x")
    print(f"  Vocabulary size: {len(set(final_seq))}")

    # ===================================================================
    # Check 1: Lossless reversibility
    # ===================================================================
    print("\n" + "=" * 70)
    print("[3] Lossless Reversibility Check")
    print("=" * 70)

    decoded = history.decode(final_seq)
    is_lossless = decoded == original_str
    print(f"  Original length: {len(original_str):,}")
    print(f"  Decoded length:  {len(decoded):,}")
    print(f"  Match: {is_lossless}")
    if not is_lossless:
        # Find first mismatch
        for i, (a, b) in enumerate(zip(original_str, decoded)):
            if a != b:
                print(f"  First mismatch at position {i}: original={a}, decoded={b}")
                break

    # ===================================================================
    # Check 2: Hierarchy analysis
    # ===================================================================
    print("\n" + "=" * 70)
    print("[4] Merge Hierarchy Analysis")
    print("=" * 70)

    # Group merges by symbol length (= hierarchy level)
    by_length: dict[int, list] = {}
    for step, name, score, seq_len, vocab in stats:
        length = len(name)
        by_length.setdefault(length, []).append((name, score))

    print(f"\n  Merges by resulting symbol length:")
    for length in sorted(by_length.keys()):
        merges = by_length[length]
        names = [m[0] for m in merges]
        scores = [m[1] for m in merges]
        print(f"\n  Length {length} ({len(merges)} merges, avg score {np.mean(scores):.3f}):")
        # Show top 10 by score
        sorted_merges = sorted(merges, key=lambda x: -x[1])
        for name, score in sorted_merges[:10]:
            # Biological interpretation
            bio = ""
            if name == "ATG":
                bio = " ← START CODON"
            elif name in ("TAA", "TAG", "TGA"):
                bio = " ← STOP CODON"
            elif name == "TATAAT":
                bio = " ← Pribnow box (-10)"
            elif name in ("AGGAGG", "GGAGG", "AGGAG"):
                bio = " ← Shine-Dalgarno (RBS)"
            elif name == "TTGACA":
                bio = " ← -35 box"
            elif len(name) == 3 and all(c in 'ACGT' for c in name):
                bio = " (codon)"
            elif len(name) == 2 and all(c in 'ACGT' for c in name):
                bio = " (dinucleotide)"
            print(f"    {name:15s} score={score:.3f}{bio}")

    # ===================================================================
    # Check 3: Codon detection at Level 2
    # ===================================================================
    print("\n" + "=" * 70)
    print("[5] Codon Detection (Level 2 = length 3)")
    print("=" * 70)

    if 3 in by_length:
        codons_found = [name for name, _ in by_length[3]]
        all_codons = [NT_CHARS[a] + NT_CHARS[b] + NT_CHARS[c]
                      for a in range(4) for b in range(4) for c in range(4)]

        print(f"  Codons merged: {len(codons_found)}/64 possible")

        # Check for ATG (start codon)
        atg_found = 'ATG' in codons_found
        print(f"  ATG (start codon): {'FOUND' if atg_found else 'NOT FOUND'}")

        # Check for stop codons
        for stop in ['TAA', 'TAG', 'TGA']:
            found = stop in codons_found
            print(f"  {stop} (stop codon):  {'FOUND' if found else 'NOT FOUND'}")

        # What fraction of standard codons were found?
        standard_found = sum(1 for c in codons_found if c in all_codons)
        print(f"\n  Standard codons found: {standard_found}/{len(codons_found)} merged trimers")

        # What are the NON-codon trimers (spanning codon boundaries)?
        non_standard = [c for c in codons_found if c not in all_codons]
        if non_standard:
            print(f"  Non-ACGT trimers: {non_standard[:10]}")
    else:
        print("  No length-3 merges found!")

    # ===================================================================
    # Check 4: SDR on score distribution
    # ===================================================================
    print("\n" + "=" * 70)
    print("[6] SDR on Merge Score Distribution")
    print("=" * 70)

    all_scores = {name: score for _, name, score, _, _ in stats}
    theta = sdr_threshold(all_scores)
    structural = {name for name, score in all_scores.items() if score > theta}

    print(f"  Total merges: {len(all_scores)}")
    print(f"  SDR threshold: {theta:.3f}")
    print(f"  Structural (score > θ): {len(structural)}")
    print(f"  Structural merges (carrying ~50% score mass):")
    sorted_structural = sorted(structural, key=lambda n: -all_scores[n])
    for name in sorted_structural[:15]:
        bio = ""
        if name == "ATG":
            bio = "START"
        elif name in ("TAA", "TAG", "TGA"):
            bio = "STOP"
        elif len(name) == 2:
            bio = "dinuc"
        elif len(name) == 3:
            bio = "codon"
        elif len(name) >= 4:
            bio = f"motif({len(name)}nt)"
        print(f"    {name:15s} score={all_scores[name]:.3f}  [{bio}]")

    # ===================================================================
    # Check 5: Compression profile — when does compression slow down?
    # ===================================================================
    print("\n" + "=" * 70)
    print("[7] Compression Profile")
    print("=" * 70)

    if stats:
        lengths = [s[3] for s in stats]
        scores_list = [s[2] for s in stats]
        print(f"  Step   Score    Seq len  Compression  Symbol")
        print("  " + "-" * 55)
        milestones = [0, 1, 2, 3, 5, 10, 15, 20, 30, 50, 75, 100, 150, 199]
        for m in milestones:
            if m < len(stats):
                step, name, score, seq_len, vocab = stats[m]
                ratio = len(seq_int) / seq_len
                print(f"  {step:4d}  {score:7.3f}  {seq_len:7,}  {ratio:6.2f}x       {name}")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    checks = []

    # Check 1: Lossless
    checks.append(("Lossless reversibility", is_lossless))
    print(f"  [{'PASS' if is_lossless else 'FAIL'}] Lossless reversibility")

    # Check 2: ATG in level 2
    atg_in_l2 = 3 in by_length and 'ATG' in [n for n, _ in by_length[3]]
    checks.append(("ATG in level 2 (trimers)", atg_in_l2))
    print(f"  [{'PASS' if atg_in_l2 else 'FAIL'}] ATG in level-2 merges")

    # Check 3: ≥50% of all codons appear as trimers
    if 3 in by_length:
        n_codons = len(by_length[3])
        pct = n_codons / 64 * 100
        ok = pct >= 50
        checks.append(("≥50% of 64 codons merged as trimers", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] Codons as trimers: {n_codons}/64 ({pct:.0f}%)")
    else:
        checks.append(("≥50% of 64 codons merged as trimers", False))
        print(f"  [FAIL] No trimers formed")

    # Check 4: SDR gives meaningful threshold
    ok = 0 < theta < max(all_scores.values()) if all_scores else False
    checks.append(("SDR threshold meaningful", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] SDR threshold = {theta:.3f}")

    # Check 5: Compression > 1.5x
    if stats:
        final_ratio = len(seq_int) / len(final_seq)
        ok = final_ratio > 1.5
        checks.append(("Compression > 1.5x", ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] Compression: {final_ratio:.2f}x")

    passed = sum(1 for _, ok in checks if ok)
    print(f"\n  Result: {passed}/{len(checks)} checks passed")
    print(f"  {'SUCCESS — proceed to 004' if passed >= 3 else 'NEEDS INVESTIGATION'}")


if __name__ == "__main__":
    main()
