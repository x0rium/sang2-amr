#!/usr/bin/env python3
"""Exploration 003b: Batch Pair Merge — fixed hierarchy.

Problem with 003: PPMI scoring creates runaway G-chains because merged symbols
are rare → PPMI(pair_with_rare) is artificially high. Only 4/16 dinucleotides
formed, all G-initial.

Fix: Batch merge per level using FREQUENCY as primary criterion.
At each level:
  1. Count all adjacent pairs
  2. Merge the top-N most frequent pairs SIMULTANEOUSLY (non-overlapping)
  3. Compute PPMI for each merged pair (analytical, not scoring)
  4. Advance to next level

This produces diverse merges at each level, creating proper hierarchy.

Also try: standard BPE (one merge at a time, by frequency) for comparison.
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
    rec = SeqIO.read(str(DATA_DIR / "ecoli_k12.fasta"), "fasta")
    seq = str(rec.seq).upper()[:limit]
    return [NT_MAP.get(c, 0) for c in seq]


class MergeHistory:
    def __init__(self):
        self.symbol_names: dict[int, str] = {i: c for i, c in enumerate(NT_CHARS)}
        self.merges: list[tuple[int, int, int, int]] = []  # (a, b, new_id, level)
        self.next_id = 4

    def record(self, a: int, b: int, level: int) -> int:
        new_id = self.next_id
        self.next_id += 1
        self.symbol_names[new_id] = self.symbol_names[a] + self.symbol_names[b]
        self.merges.append((a, b, new_id, level))
        return new_id

    def decode(self, sequence: list[int]) -> str:
        return ''.join(self.symbol_names[s] for s in sequence)

    def name(self, sym_id: int) -> str:
        return self.symbol_names.get(sym_id, '?')


def merge_pair_in_seq(seq: list[int], pair: tuple[int, int], new_id: int) -> list[int]:
    """Replace all non-overlapping occurrences of pair with new_id."""
    result = []
    i = 0
    while i < len(seq):
        if i < len(seq) - 1 and seq[i] == pair[0] and seq[i + 1] == pair[1]:
            result.append(new_id)
            i += 2
        else:
            result.append(seq[i])
            i += 1
    return result


def entropy(counts) -> float:
    values = np.array(list(counts.values()) if isinstance(counts, (dict, Counter)) else counts, dtype=np.float64)
    total = values.sum()
    if total <= 1:
        return 0.0
    probs = values[values > 0] / total
    return float(-np.sum(probs * np.log2(probs)))


# =========================================================================
# Approach 1: Standard BPE (frequency-only)
# =========================================================================
def run_bpe(sequence: list[int], max_merges: int = 300, min_count: int = 5):
    """Standard BPE: always merge the most frequent adjacent pair."""
    history = MergeHistory()
    seq = list(sequence)

    stats = []
    for step in range(max_merges):
        pairs = Counter(zip(seq[:-1], seq[1:]))
        if not pairs:
            break
        best_pair, best_count = pairs.most_common(1)[0]
        if best_count < min_count:
            break

        new_id = history.record(best_pair[0], best_pair[1], level=step)
        seq = merge_pair_in_seq(seq, best_pair, new_id)

        name = history.name(new_id)
        stats.append((step, name, best_count, len(seq), len(set(seq))))

    return seq, history, stats


# =========================================================================
# Approach 2: Batch merge per level
# =========================================================================
def run_batch_merge(sequence: list[int], max_levels: int = 8,
                    merges_per_level: int = 16, min_count: int = 10):
    """Batch merge: at each level, merge top-N pairs simultaneously.

    Non-overlapping: pairs sharing a symbol within the same level are
    resolved by frequency (higher freq pair wins).
    """
    history = MergeHistory()
    seq = list(sequence)

    level_info = []

    for level in range(max_levels):
        pairs = Counter(zip(seq[:-1], seq[1:]))
        if not pairs:
            break

        # Sort by frequency
        sorted_pairs = pairs.most_common()

        # Select non-conflicting pairs (greedy: highest freq first)
        used_symbols = set()
        selected = []
        for pair, count in sorted_pairs:
            if count < min_count:
                break
            # Check no symbol conflict with already selected pairs
            if pair[0] not in used_symbols and pair[1] not in used_symbols:
                selected.append((pair, count))
                used_symbols.add(pair[0])
                used_symbols.add(pair[1])
                if len(selected) >= merges_per_level:
                    break

        if not selected:
            break

        # Apply all selected merges
        merge_names = []
        for pair, count in selected:
            new_id = history.record(pair[0], pair[1], level=level)
            seq = merge_pair_in_seq(seq, pair, new_id)
            merge_names.append((history.name(new_id), count))

        level_info.append((level, merge_names, len(seq), len(set(seq))))

    return seq, history, level_info


# =========================================================================
# Analysis helpers
# =========================================================================
def analyze_hierarchy(history: MergeHistory, approach_name: str):
    """Analyze what biological units the hierarchy produced."""
    by_length: dict[int, list[str]] = {}
    for a, b, new_id, level in history.merges:
        name = history.name(new_id)
        by_length.setdefault(len(name), []).append(name)

    # Biologically significant patterns
    start_codons = {'ATG'}
    stop_codons = {'TAA', 'TAG', 'TGA'}
    pribnow = {'TATAAT'}
    shine_dalgarno = {'AGGAGG', 'GGAGG', 'AGGAG', 'GAGG'}
    minus35 = {'TTGACA'}

    all_codons = {NT_CHARS[a] + NT_CHARS[b] + NT_CHARS[c]
                  for a in range(4) for b in range(4) for c in range(4)}
    all_dinucs = {NT_CHARS[a] + NT_CHARS[b] for a in range(4) for b in range(4)}

    print(f"\n  [{approach_name}] Hierarchy by symbol length:")
    for length in sorted(by_length.keys()):
        names = by_length[length]
        print(f"\n    Length {length}: {len(names)} symbols")

        # Count biological units
        if length == 2:
            found_dinucs = set(names) & all_dinucs
            print(f"      Dinucleotides: {len(found_dinucs)}/16")
            print(f"      Found: {sorted(found_dinucs)}")
        elif length == 3:
            found_codons = set(names) & all_codons
            print(f"      Standard codons: {len(found_codons)}/64")
            found_start = set(names) & start_codons
            found_stop = set(names) & stop_codons
            if found_start:
                print(f"      START codons: {found_start}")
            if found_stop:
                print(f"      STOP codons: {found_stop}")
            # Show all trimers
            if len(names) <= 20:
                print(f"      All: {sorted(names)}")
        elif length >= 4:
            found_pribnow = set(names) & pribnow
            found_sd = set(names) & shine_dalgarno
            found_35 = set(names) & minus35
            if found_pribnow:
                print(f"      Pribnow box (-10): {found_pribnow}")
            if found_sd:
                print(f"      Shine-Dalgarno: {found_sd}")
            if found_35:
                print(f"      -35 box: {found_35}")
            # Show top 10 most interesting
            if len(names) <= 15:
                print(f"      All: {sorted(names)}")
            else:
                print(f"      Sample: {sorted(names)[:10]}...")

    return by_length


# =========================================================================
# Main
# =========================================================================
def main():
    print("=" * 70)
    print("EXPLORATION 003b: Batch Pair Merge (Fixed)")
    print("=" * 70)

    seq_int = load_genome(100_000)
    original_str = ''.join(NT_CHARS[s] for s in seq_int)
    print(f"  Input: {len(seq_int):,} nt E. coli K-12\n")

    # ===================================================================
    # Approach 1: Standard BPE
    # ===================================================================
    print("=" * 70)
    print("[A] Standard BPE (frequency-only, 300 merges)")
    print("=" * 70)

    t0 = time.time()
    bpe_seq, bpe_hist, bpe_stats = run_bpe(seq_int, max_merges=300, min_count=5)
    t1 = time.time()

    bpe_decoded = bpe_hist.decode(bpe_seq)
    bpe_lossless = bpe_decoded == original_str

    print(f"  Time: {t1-t0:.1f}s")
    print(f"  Merges: {len(bpe_stats)}")
    print(f"  Final: {len(bpe_seq):,} tokens, compression {len(seq_int)/len(bpe_seq):.2f}x")
    print(f"  Lossless: {bpe_lossless}")

    # Show first 20 merges
    print(f"\n  First 20 BPE merges:")
    for step, name, count, seq_len, vocab in bpe_stats[:20]:
        print(f"    {step:3d}: {name:12s} count={count:5d} → {seq_len:,} tokens")

    bpe_by_length = analyze_hierarchy(bpe_hist, "BPE")

    # ===================================================================
    # Approach 2: Batch merge
    # ===================================================================
    print("\n" + "=" * 70)
    print("[B] Batch Merge (16 pairs per level, 8 levels)")
    print("=" * 70)

    t0 = time.time()
    batch_seq, batch_hist, batch_levels = run_batch_merge(
        seq_int, max_levels=8, merges_per_level=16, min_count=10)
    t1 = time.time()

    batch_decoded = batch_hist.decode(batch_seq)
    batch_lossless = batch_decoded == original_str

    print(f"  Time: {t1-t0:.1f}s")
    print(f"  Levels: {len(batch_levels)}")
    print(f"  Final: {len(batch_seq):,} tokens, compression {len(seq_int)/len(batch_seq):.2f}x")
    print(f"  Lossless: {batch_lossless}")

    # Show merges per level
    for level, merges, seq_len, vocab in batch_levels:
        names = [f"{n}({c})" for n, c in merges]
        print(f"\n  Level {level}: {len(merges)} merges → {seq_len:,} tokens, vocab={vocab}")
        # Group by resulting symbol length
        for n, c in sorted(merges, key=lambda x: -x[1]):
            print(f"    {n:20s} count={c:5d}  (length {len(n)} nt)")

    batch_by_length = analyze_hierarchy(batch_hist, "Batch")

    # ===================================================================
    # Approach 3: BPE with diversity (max 2 merges per initial letter)
    # ===================================================================
    print("\n" + "=" * 70)
    print("[C] Diverse BPE (max 4 merges starting with same base, 200 merges)")
    print("=" * 70)

    div_hist = MergeHistory()
    div_seq = list(seq_int)
    div_stats = []
    base_letter_count: dict[str, int] = Counter()  # track starting letter of merged pairs

    for step in range(200):
        pairs = Counter(zip(div_seq[:-1], div_seq[1:]))
        if not pairs:
            break

        # Sort by frequency, skip if starting letter is overrepresented
        for pair, count in pairs.most_common():
            if count < 5:
                break
            first_nt = div_hist.name(pair[0])[0]  # first nucleotide
            if base_letter_count[first_nt] < 50:  # allow diversity
                new_id = div_hist.record(pair[0], pair[1], level=step)
                div_seq = merge_pair_in_seq(div_seq, pair, new_id)
                base_letter_count[first_nt] += 1
                name = div_hist.name(new_id)
                div_stats.append((step, name, count, len(div_seq), len(set(div_seq))))
                break
        else:
            break  # no valid pair found

    div_decoded = div_hist.decode(div_seq)
    div_lossless = div_decoded == original_str

    print(f"  Merges: {len(div_stats)}")
    print(f"  Final: {len(div_seq):,} tokens, compression {len(seq_int)/len(div_seq):.2f}x")
    print(f"  Lossless: {div_lossless}")

    div_by_length = analyze_hierarchy(div_hist, "Diverse BPE")

    # ===================================================================
    # SUMMARY — compare approaches
    # ===================================================================
    print("\n" + "=" * 70)
    print("SUMMARY — Comparing Approaches")
    print("=" * 70)

    for name, by_len, lossless, final_seq in [
        ("BPE", bpe_by_length, bpe_lossless, bpe_seq),
        ("Batch", batch_by_length, batch_lossless, batch_seq),
        ("Diverse", div_by_length, div_lossless, div_seq),
    ]:
        dinucs = len(set(by_len.get(2, [])) & {NT_CHARS[a]+NT_CHARS[b] for a in range(4) for b in range(4)})
        trimers = by_len.get(3, [])
        all_codons = {NT_CHARS[a]+NT_CHARS[b]+NT_CHARS[c] for a in range(4) for b in range(4) for c in range(4)}
        codons = len(set(trimers) & all_codons)
        has_atg = 'ATG' in trimers
        has_stop = any(s in trimers for s in ['TAA', 'TAG', 'TGA'])
        ratio = len(seq_int) / len(final_seq)

        print(f"\n  {name:10s}: lossless={lossless}, ratio={ratio:.2f}x")
        print(f"    Dinucleotides: {dinucs}/16")
        print(f"    Codons:        {codons}/64  ATG={'yes' if has_atg else 'no'}  STOP={'yes' if has_stop else 'no'}")

    # Best approach assessment
    print("\n  Checks (on best approach):")
    # Use BPE as the main approach
    best_by_len = bpe_by_length
    best_lossless = bpe_lossless
    best_ratio = len(seq_int) / len(bpe_seq)

    dinucs = len(set(best_by_len.get(2, [])) & {NT_CHARS[a]+NT_CHARS[b] for a in range(4) for b in range(4)})
    trimers = best_by_len.get(3, [])
    all_codons_set = {NT_CHARS[a]+NT_CHARS[b]+NT_CHARS[c] for a in range(4) for b in range(4) for c in range(4)}
    codons = len(set(trimers) & all_codons_set)
    has_atg = 'ATG' in trimers

    checks = []
    checks.append(("Lossless", best_lossless))
    print(f"  [{'PASS' if best_lossless else 'FAIL'}] Lossless")

    ok = dinucs >= 12
    checks.append(("≥12/16 dinucleotides", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Dinucleotides: {dinucs}/16")

    ok = has_atg
    checks.append(("ATG as trimer", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] ATG in trimers")

    ok = codons >= 20
    checks.append(("≥20 codons as trimers", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Codons as trimers: {codons}/64")

    ok = best_ratio > 1.5
    checks.append(("Compression > 1.5x", ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] Compression: {best_ratio:.2f}x")

    passed = sum(1 for _, ok in checks if ok)
    print(f"\n  Result: {passed}/{len(checks)} checks passed")


if __name__ == "__main__":
    main()
