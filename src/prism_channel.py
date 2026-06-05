"""PRISM structural channel — the second voice in two-channel scoring.

Composition (anomaly.py) detects recently-acquired / foreign DNA but goes blind on
HOST-ADAPTED resistance genes (CTX-M, OXA): when a gene's composition matches the
host, GC/codon signals can't separate it. PRISM's richer structural profile (14
metrics: fractal, recurrence, permutation entropy, topology) still can.

Validated in exploration/010-013: scoring a region by the distance between its
PRISM profile and the whole-replicon profile lifts adapted CTX-M from AUC ~0.46 to
~0.77. Fused with composition (soft weight by plasmid heterogeneity) the combined
detector reaches ~0.74 vs ~0.65 for composition alone, without collapsing on
adapted genes.

PRISM (`sang`) lives outside this repo; this module imports it optionally. If it
is unavailable, available() returns False and the pipeline falls back to
composition-only scoring.
"""
from __future__ import annotations

import sys

import numpy as np

# Calibration constants for the heterogeneity -> fusion-weight sigmoid.
# Derived from 76 RefSeq AMR plasmids (exploration/013); see 013_results.md.
# A single genome has no set to rank against, so we map its heterogeneity through
# this fixed reference distribution instead of a transductive rank.
HET_MEDIAN = 0.00844
HET_SCALE = 0.00445

_PRISM_PATH = "/Users/x0rium/CATALYSTO/Prism"
_NT = np.array(list("ACGT"))

try:
    if _PRISM_PATH not in sys.path:
        sys.path.insert(0, _PRISM_PATH)
    import sang  # type: ignore
    _AVAILABLE = True
except Exception:  # pragma: no cover - PRISM optional
    sang = None
    _AVAILABLE = False


def available() -> bool:
    """True if PRISM/sang could be imported."""
    return _AVAILABLE


def _decode(arr: np.ndarray) -> list:
    """uint8 A/C/G/T codes -> char tokens for sang's symbolic metrics."""
    clipped = np.clip(arr, 0, 3)
    return _NT[clipped].tolist()


def profile_of(arr: np.ndarray):
    """PRISM structural profile of a sequence (uint8 A=0..T=3). None if too short."""
    if not _AVAILABLE or len(arr) < 150:
        return None
    return sang.profile(_decode(arr), signal=arr.astype(float), depth="standard")


def region_distance(region: np.ndarray, whole_profile) -> float:
    """Structural distance between a region's profile and the whole-replicon one."""
    if not _AVAILABLE or whole_profile is None or len(region) < 150:
        return 0.0
    try:
        p = sang.profile(_decode(region), signal=region.astype(float), depth="standard")
        return float(sang.compare(p, whole_profile).distance)
    except Exception:  # pragma: no cover
        return 0.0


def heterogeneity(arr: np.ndarray, window: int = 2000, step: int = 500) -> float:
    """Spread of 3-mer entropy across windows — the regime cue (013)."""
    from signals import kmer_entropy
    if len(arr) < window:
        return HET_MEDIAN
    h3 = [kmer_entropy(arr[i:i + window], k=3)
          for i in range(0, len(arr) - window, step)]
    return float(np.std(h3)) if h3 else HET_MEDIAN


def fusion_weight(arr: np.ndarray) -> float:
    """Weight on the PRISM channel in [0,1], from plasmid heterogeneity.

    High heterogeneity (big mosaic replicon, adapted genes) -> trust PRISM more;
    low heterogeneity (compact replicon, recent island) -> trust composition.
    Sign and calibration validated in exploration/013 (soft-fusion 0.743).
    """
    het = heterogeneity(arr)
    z = (het - HET_MEDIAN) / HET_SCALE
    return float(1.0 / (1.0 + np.exp(-z)))
