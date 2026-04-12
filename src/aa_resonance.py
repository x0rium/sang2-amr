"""Protein-level resonance scoring for AMR function detection.

Translates ORFs to amino acids, computes AA-pair statistics,
and scores similarity to known AMR protein resonance profiles.

This is the FUNCTION-level detector: composition tells us "foreign DNA",
codon usage tells us "different organism", AA-resonance tells us
"this protein acts like a β-lactamase."

Validated: S-X-N (ΔR=0.229), H-X-X-X-D (ΔR=0.146), S-X-X-K (ΔR=0.085)
all found ab initio on 114 β-lactamases vs 80 kinases.
"""

from __future__ import annotations

import numpy as np
import json
from pathlib import Path

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_MAP = {c: i for i, c in enumerate(AA_ALPHABET)}
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
NT_CHARS = "ACGT"

_PROFILE_CACHE: list | None = None
_PROFILE_PATH = Path(__file__).parent / "amr_profile.json"


def _load_profile() -> list:
    """Load AMR resonance profile (discriminative AA-pairs).

    Format: list of [b1, b2, k, delta_R, R_bl, R_ctrl, P] arrays.
    Converted to list of ((b1, b2, k), {'delta_R': ..., ...}) for compat.
    """
    global _PROFILE_CACHE
    if _PROFILE_CACHE is None:
        with open(_PROFILE_PATH) as f:
            raw = json.load(f)
        _PROFILE_CACHE = [
            ((e[0], e[1], e[2]),
             {"delta_R": e[3], "R_bl": e[4], "R_ctrl": e[5], "P": e[6]})
            for e in raw
        ]
    return _PROFILE_CACHE


def translate_orf(nt_data: np.ndarray, frame: int = 0) -> np.ndarray:
    """Translate nucleotide array to AA index array (0..19).

    Stops at first stop codon. Returns empty array if too short.
    """
    n = len(nt_data)
    start = frame
    aa_list = []
    for i in range(start, n - 2, 3):
        codon = NT_CHARS[nt_data[i]] + NT_CHARS[nt_data[i + 1]] + NT_CHARS[nt_data[i + 2]]
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            break
        if aa in AA_MAP:
            aa_list.append(AA_MAP[aa])
    return np.array(aa_list, dtype=np.uint8)


def find_orfs(nt_data: np.ndarray, min_aa: int = 80) -> list[tuple[int, int, np.ndarray]]:
    """Find all ORFs ≥ min_aa in all 6 reading frames.

    Returns list of (start_nt, end_nt, aa_array).
    """
    orfs = []
    complement = np.array([3, 2, 1, 0], dtype=np.uint8)

    for strand_data, strand_offset in [(nt_data, 0), (complement[nt_data[::-1]], len(nt_data))]:
        for frame in range(3):
            n = len(strand_data)
            i = frame
            while i < n - 2:
                # Look for ATG
                if (strand_data[i] == 0 and strand_data[i + 1] == 3 and strand_data[i + 2] == 2):  # ATG
                    aa = translate_orf(strand_data, frame=i)
                    if len(aa) >= min_aa:
                        if strand_offset == 0:
                            orfs.append((i, i + len(aa) * 3, aa))
                        else:
                            # Reverse strand: convert coordinates
                            real_end = len(nt_data) - i
                            real_start = real_end - len(aa) * 3
                            orfs.append((real_start, real_end, aa))
                        i += len(aa) * 3
                        continue
                i += 3
    return orfs


def score_orf_amr(aa_data: np.ndarray, top_n: int = 50) -> float:
    """Score an ORF's AA sequence against the AMR resonance profile.

    Computes AA-pair frequencies in this ORF, then checks how many
    of the top discriminative AMR pairs are present.

    Returns score in [0, 1]: higher = more AMR-like.
    """
    profile = _load_profile()
    n = len(aa_data)
    if n < 30:
        return 0.0

    # Compute AA-pair counts for this ORF
    pair_counts: dict[tuple[int, int, int], int] = {}
    for k in range(1, 11):
        if k >= n:
            break
        for i in range(n - k):
            key = (int(aa_data[i]), int(aa_data[i + k]), k)
            pair_counts[key] = pair_counts.get(key, 0) + 1

    # Score: weighted match against top discriminative pairs
    total_weight = 0.0
    matched_weight = 0.0

    for (b1, b2, k), info in profile[:top_n]:
        delta_r = info["delta_R"]
        total_weight += delta_r
        if (b1, b2, k) in pair_counts:
            matched_weight += delta_r

    return float(matched_weight / total_weight) if total_weight > 0 else 0.0


def scan_window_aa(nt_data: np.ndarray, min_aa: int = 60, top_n: int = 100) -> float:
    """Score a nucleotide window for AMR-like protein content.

    Finds all ORFs in the window, scores each against AMR profile,
    returns the MAX score (best ORF).
    """
    orfs = find_orfs(nt_data, min_aa=min_aa)
    if not orfs:
        return 0.0
    scores = [score_orf_amr(aa, top_n=top_n) for _, _, aa in orfs]
    return max(scores)
