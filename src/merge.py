"""BPE pair merge for nucleotide hierarchy construction.

Frequency-based BPE (NOT PPMI-weighted — see gotcha G17).
Builds lossless hierarchy: nt → dinucleotides → codons → motifs.

Validated in exploration/003b: 3.62x compression, 23/64 codons found,
ATG + TAA + TAG discovered automatically.
"""

from __future__ import annotations

import numpy as np
from collections import Counter


class MergeHistory:
    """Track merge hierarchy for lossless reversal."""

    def __init__(self, base_alphabet: list[str]):
        self.symbol_names: dict[int, str] = {
            i: name for i, name in enumerate(base_alphabet)
        }
        self.merges: list[tuple[int, int, int, int]] = []  # (a, b, new_id, level)
        self.next_id = len(base_alphabet)

    def record(self, a: int, b: int, level: int = 0) -> int:
        new_id = self.next_id
        self.next_id += 1
        self.symbol_names[new_id] = self.symbol_names[a] + self.symbol_names[b]
        self.merges.append((a, b, new_id, level))
        return new_id

    def decode(self, sequence: list[int]) -> str:
        """Decode merged sequence back to original characters."""
        return "".join(self.symbol_names[s] for s in sequence)

    def name(self, sym_id: int) -> str:
        return self.symbol_names.get(sym_id, "?")

    def symbol_length(self, sym_id: int) -> int:
        return len(self.symbol_names.get(sym_id, ""))


def _merge_pair(seq: list[int], pair: tuple[int, int], new_id: int) -> list[int]:
    """Replace all non-overlapping occurrences of pair with new_id (left-to-right greedy)."""
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


def run_bpe(
    sequence: list[int],
    max_merges: int = 300,
    min_count: int = 5,
    base_alphabet: list[str] | None = None,
) -> tuple[list[int], MergeHistory]:
    """Standard BPE: iteratively merge the most frequent adjacent pair.

    Returns (compressed_sequence, history).
    history.decode(compressed) recovers the original.
    """
    if base_alphabet is None:
        base_alphabet = ["A", "C", "G", "T"]
    history = MergeHistory(base_alphabet)
    seq = list(sequence)

    for step in range(max_merges):
        pairs = Counter(zip(seq[:-1], seq[1:]))
        if not pairs:
            break
        best_pair, best_count = pairs.most_common(1)[0]
        if best_count < min_count:
            break
        new_id = history.record(best_pair[0], best_pair[1], level=step)
        seq = _merge_pair(seq, best_pair, new_id)

    return seq, history


def run_batch_merge(
    sequence: list[int],
    max_levels: int = 8,
    merges_per_level: int = 16,
    min_count: int = 10,
    base_alphabet: list[str] | None = None,
) -> tuple[list[int], MergeHistory]:
    """Batch merge: at each level, merge top-N non-conflicting pairs simultaneously.

    Produces more diverse hierarchy per level than sequential BPE.
    Returns (compressed_sequence, history).
    """
    if base_alphabet is None:
        base_alphabet = ["A", "C", "G", "T"]
    history = MergeHistory(base_alphabet)
    seq = list(sequence)

    for level in range(max_levels):
        pairs = Counter(zip(seq[:-1], seq[1:]))
        if not pairs:
            break

        used_symbols: set[int] = set()
        selected: list[tuple[tuple[int, int], int]] = []
        for pair, count in pairs.most_common():
            if count < min_count:
                break
            if pair[0] not in used_symbols and pair[1] not in used_symbols:
                selected.append((pair, count))
                used_symbols.add(pair[0])
                used_symbols.add(pair[1])
                if len(selected) >= merges_per_level:
                    break

        if not selected:
            break

        for pair, _ in selected:
            new_id = history.record(pair[0], pair[1], level=level)
            seq = _merge_pair(seq, pair, new_id)

    return seq, history
