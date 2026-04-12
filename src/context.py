"""Genomic context analysis — FSM-inspired scoring.

Maps the auto-reverser's FSM concept to genomics: instead of protocol
state transitions (message_type → message_type), we model genomic region
transitions (region_class → region_class) along a sequence.

Key insight: AMR genes are embedded in mobile genetic elements with
characteristic grammar:
  Transposon:    IS → transposase → [cargo/AMR] → transposase → IS
  Integron:      integrase → attI → cassette → attC → cassette → attC
  Genomic island: tRNA → [island content] → tRNA_fragment

A window adjacent to transposase/integrase-like regions gets a "mobility
context boost" to its anomaly score.

Validated on pKpQIL: Tn4401 structure (tnpA → KPC-2 → tnpA).
"""

from __future__ import annotations

import numpy as np

from ir import WindowScore
from signals import (
    gc_content, kmer_entropy, hurst, autocorr_period3,
    chargaff2_dev, otsu_threshold,
)


# ---------------------------------------------------------------------------
# Region classification (window → class label)
# ---------------------------------------------------------------------------
def classify_window(scores: dict[str, float]) -> str:
    """Classify a window into one of 4 classes based on z-scores.

    Classes (genomic "states"):
      NORMAL      — no signal deviates strongly
      ANOMALOUS   — multiple signals deviate (AMR candidate)
      MOBILE_LIKE — high gradient + Chargaff-2 (transposase/IS pattern)
      STRUCTURAL  — high Hurst deviation (genomic island boundary)
    """
    gc_z = scores.get("gc", 0)
    h3_z = scores.get("h3", 0)
    hurst_z = scores.get("hurst", 0)
    c2_z = scores.get("c2", 0)
    gradient_z = scores.get("gradient", 0)
    ac3_z = scores.get("ac3", 0)

    # Mobile-like: high gradient (composition boundary) + high c2 or ac3
    if gradient_z > 2.0 and (c2_z > 1.5 or ac3_z > 1.5):
        return "MOBILE_LIKE"

    # Structural: strong Hurst anomaly (long-range correlation change)
    if hurst_z > 2.5:
        return "STRUCTURAL"

    # Anomalous: multiple moderate deviations
    n_elevated = sum(1 for s in [gc_z, h3_z, hurst_z, c2_z, ac3_z]
                     if s > 1.5)
    if n_elevated >= 2:
        return "ANOMALOUS"

    return "NORMAL"


# ---------------------------------------------------------------------------
# Transition mining along sequence
# ---------------------------------------------------------------------------
def mine_region_transitions(
    windows: list[WindowScore],
) -> list[tuple[str, str, int]]:
    """Mine (class_from → class_to) transitions from ordered windows.

    Returns list of (from_class, to_class, count) sorted by count desc.
    """
    if len(windows) < 2:
        return []

    # Classify each window
    classes = [classify_window(w.scores) for w in windows]

    # Count transitions
    from collections import Counter
    trans = Counter()
    for i in range(len(classes) - 1):
        trans[(classes[i], classes[i + 1])] += 1

    return [(a, b, c) for (a, b), c in trans.most_common()]


# ---------------------------------------------------------------------------
# Context boost scoring
# ---------------------------------------------------------------------------
def compute_context_boost(
    windows: list[WindowScore],
    proximity: int = 3,
) -> list[float]:
    """Compute a mobility context boost for each window.

    A window gets a boost if it's near MOBILE_LIKE or STRUCTURAL windows.
    The boost is proportional to:
      - Number of MOBILE_LIKE neighbors within `proximity` windows
      - Whether the pattern is "flanking" (mobile on both sides)

    Returns boost multiplier per window (1.0 = no boost, >1.0 = boosted).
    """
    n = len(windows)
    classes = [classify_window(w.scores) for w in windows]
    boosts = [1.0] * n

    for i in range(n):
        # Only boost potential cargo regions (ANOMALOUS or NORMAL),
        # NOT mobile elements themselves (they're already highly ranked).
        if classes[i] in ("MOBILE_LIKE", "STRUCTURAL"):
            continue

        # Count neighbors that are anomalous/mobile/structural (= elevated signals)
        elevated_left = 0
        elevated_right = 0
        for j in range(max(0, i - proximity), i):
            if classes[j] != "NORMAL":
                elevated_left += 1
        for j in range(i + 1, min(n, i + proximity + 1)):
            if classes[j] != "NORMAL":
                elevated_right += 1

        total_elevated = elevated_left + elevated_right
        if total_elevated == 0:
            continue

        # Flanking = both sides have elevated regions (classic transposon context).
        # Epsilon boost only: 3-8% per neighbor. Enough to break ties, not to
        # override signal-based ranking.
        if elevated_left > 0 and elevated_right > 0:
            boosts[i] = 1.0 + 0.04 * total_elevated  # flanking: max +24%
        else:
            boosts[i] = 1.0 + 0.02 * total_elevated  # one-sided: max +12%

    return boosts


def apply_context_boost(windows: list[WindowScore]) -> list[WindowScore]:
    """Re-score windows with context boost and re-sort.

    Modifies windows in-place: composite *= context_boost.
    Also adds 'context_class' and 'context_boost' to scores dict.
    """
    # First, need windows in POSITIONAL order for context analysis
    positional = sorted(windows, key=lambda w: w.start)
    boosts = compute_context_boost(positional, proximity=3)
    classes = [classify_window(w.scores) for w in positional]

    for w, boost, cls in zip(positional, boosts, classes):
        w.scores["_context_class"] = {"NORMAL": 0, "ANOMALOUS": 1,
                                       "MOBILE_LIKE": 2, "STRUCTURAL": 3}.get(cls, 0)
        w.scores["_context_boost"] = boost
        w.composite *= boost

    # Re-sort by composite
    windows.sort(key=lambda w: -w.composite)
    return windows


# ---------------------------------------------------------------------------
# FSM summary (diagnostic output)
# ---------------------------------------------------------------------------
def fsm_summary(windows: list[WindowScore]) -> str:
    """Human-readable FSM summary for diagnostic output."""
    positional = sorted(windows, key=lambda w: w.start)
    classes = [classify_window(w.scores) for w in positional]
    transitions = mine_region_transitions(positional)

    lines = [
        "Genomic Context FSM",
        f"  Windows: {len(positional)}",
        f"  Region classes: {dict(zip(*np.unique(classes, return_counts=True)))}",
        "",
        "  Transition map:",
    ]
    for a, b, c in transitions[:10]:
        lines.append(f"    {a:15s} → {b:15s}  n={c}")

    # Show mobile-flanked anomalous regions
    flanked = []
    for i, (w, cls) in enumerate(zip(positional, classes)):
        if cls in ("ANOMALOUS", "NORMAL") and w.scores.get("_context_boost", 1.0) > 1.0:
            flanked.append((w.start, w.end, w.scores.get("_context_boost", 1.0), cls))

    if flanked:
        lines.append(f"\n  Context-boosted regions ({len(flanked)}):")
        for start, end, boost, cls in flanked[:10]:
            lines.append(f"    {start:>6}-{end:>6}  boost={boost:.2f}  class={cls}")

    return "\n".join(lines)
