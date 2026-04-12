"""Internal representation for SANG2-AMR.

Data classes only — zero logic, zero external dependencies.
All pipeline modules reference these types.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal
import json


@dataclass
class Contig:
    """A contiguous DNA sequence (assembled contig or chromosome)."""
    id: str
    sequence: bytes  # A=0, C=1, G=2, T=3 as uint8
    quality: bytes | None = None  # Phred scores if from FASTQ
    source: str = ""  # sample ID


@dataclass
class SignalVector:
    """All SANG2 signals computed for a genomic region.

    Updated from original design based on exploration 001-006:
    - autocorr_period3 replaces mi_halflife (G15: MI decay doesn't work at |Σ|=4)
    - kmer_entropy_3 replaces h_ratio (G16: H/H_max ≈ 0.99 for all regions)
    - Added gc_content, chargaff2_dev, mi_lag1 (validated in 002-005)
    """
    entropy: float = 0.0
    sdr_ratio: float = 0.0
    delta_h: float = 0.0
    autocorr_period3: float = 0.0
    autocorr_peaks: list[int] = field(default_factory=list)
    kmer_entropy_3: float = 0.0
    regime: str = "UNKNOWN"
    hurst: float = 0.5
    gc_content: float = 0.0
    chargaff2_dev: float = 0.0
    mi_lag1: float = 0.0


@dataclass
class GenomicRegion:
    """An annotated region within a contig."""
    contig_id: str
    start: int
    end: int
    strand: Literal["+", "-"] = "+"
    kind: str = "unknown"  # orf, intergenic, repeat, mobile_element, unknown
    signals: SignalVector = field(default_factory=SignalVector)


@dataclass
class WindowScore:
    """Anomaly score for a sliding window along a sequence."""
    position: int  # center of window
    start: int
    end: int
    scores: dict[str, float] = field(default_factory=dict)
    composite: float = 0.0


@dataclass
class AMRCandidate:
    """A candidate antimicrobial resistance region found by SANG2-AMR."""
    region: GenomicRegion
    evidence: list[str] = field(default_factory=list)
    resonance_score: float = 0.0
    mobility_score: float = 0.0
    novelty: str = "novel"  # known, variant, novel
    nearest_known: str | None = None
    confidence: float = 0.0


def to_json(obj) -> str:
    """Serialize any IR dataclass to JSON."""
    def _convert(o):
        if hasattr(o, "__dataclass_fields__"):
            return {k: _convert(v) for k, v in o.__dict__.items()}
        if isinstance(o, bytes):
            return list(o)
        if isinstance(o, list):
            return [_convert(i) for i in o]
        if isinstance(o, dict):
            return {k: _convert(v) for k, v in o.items()}
        return o
    return json.dumps(_convert(obj), indent=2)
