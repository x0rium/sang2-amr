"""Two-channel anomaly scoring: composition (self) + PRISM structural, soft-fused.

This is the productionized result of exploration/007-013. Composition alone
(anomaly.scan_sequence) caps at ~0.65 AUC and inverts below 0.5 on host-adapted
AMR genes. Adding PRISM's structural profile as a second channel and soft-fusing
the two — weighting PRISM by the replicon's compositional heterogeneity — reaches
~0.74 AUC and no longer collapses on adapted CTX-M/OXA.

Drop-in: same output type as scan_sequence (list[WindowScore], sorted by
composite). If PRISM is unavailable, returns scan_sequence's result unchanged.
"""
from __future__ import annotations

import numpy as np

from ir import WindowScore
from anomaly import scan_sequence
import prism_channel


def _percentile(arr: np.ndarray) -> np.ndarray:
    """Rank to [0,1], higher value -> higher rank."""
    if len(arr) <= 1:
        return np.zeros(len(arr))
    order = np.argsort(np.argsort(arr))
    return order / (len(arr) - 1)


def scan_sequence_two_channel(
    target_seq: np.ndarray,
    baseline: dict,
    window: int = 2000,
    step: int = 500,
    prism_weight: float | None = None,
) -> list[WindowScore]:
    """Composition + PRISM soft-fused window scoring.

    1. Composition candidates via scan_sequence (deduped, ranked).
    2. For each candidate window, PRISM structural distance to the whole replicon.
    3. Soft-fuse the two channels' within-set percentiles, weight = heterogeneity.
    4. Re-rank by the fused score.

    prism_weight overrides the auto heterogeneity weight (mostly for tests).
    Falls back to composition-only if PRISM is unavailable.
    """
    windows = scan_sequence(target_seq, baseline, window=window, step=step)
    if not prism_channel.available() or len(windows) < 2:
        return windows

    whole_profile = prism_channel.profile_of(target_seq)
    if whole_profile is None:
        return windows

    cmp = np.array([
        prism_channel.region_distance(target_seq[w.start:w.end], whole_profile)
        for w in windows
    ])
    comp = np.array([w.composite for w in windows])

    w = prism_weight if prism_weight is not None else prism_channel.fusion_weight(target_seq)
    fused = (1.0 - w) * _percentile(comp) + w * _percentile(cmp)

    for win, f, c in zip(windows, fused, cmp):
        win.scores["composition"] = win.composite
        win.scores["prism_cmp"] = float(c)
        win.scores["_prism_weight"] = float(w)
        # fused is a percentile in [0,1]; the pipeline normalizes confidence as
        # composite/3, so scale by 3 to keep confidence in [0,1].
        win.composite = float(f) * 3.0

    windows.sort(key=lambda x: -x.composite)
    return windows
