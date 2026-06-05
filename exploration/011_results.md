# 011 — Two-channel detector (composition + PRISM profile)

Goal: fuse self (composition) and cmp_bg (PRISM 14-metric profile distance) into
ONE detector that holds across both regimes, without labels or magic thresholds.
18 plasmids, per-genome AUC. Scores cached in data/two_channel_cache.npz.
Run: `python -u exploration/011_two_channel.py`.

## Results (n=18, honest)

| group              | self  | cmp   | mean-fusion | tail-selector | oracle |
|--------------------|------:|------:|------------:|--------------:|-------:|
| ALL                | 0.625 | 0.650 | **0.686**   | 0.625         | 0.786  |
| ADAPTED (ctxm/oxa) | 0.449 | 0.684 | 0.581       | 0.449         | 0.684  |
| RECENT             | 0.714 | 0.632 | 0.739       | 0.714         | 0.837  |

## Findings

1. **MAX rank-fusion fails** (prior run): "one strong signal wins" backfires when
   a channel is pure noise on that regime — it promotes negatives. ADAPTED MAX
   ≈ self, PRISM's advantage lost.

2. **mean-fusion of percentiles is the best honest combiner: ALL 0.686**, beating
   the incumbent self (0.625) by +0.06, and — crucially — it does NOT collapse on
   adapted genes (0.581 vs self's 0.449) or invert (self drops to 0.15/0.27/0.34).
   Averaging damps the noisy channel instead of letting it raise negatives.

3. **The smart tail-contrast selector FAILED (negative result).** It picked `self`
   on all 18 genomes. self is a rank-fusion composite with sharp outliers; cmp is
   a smooth continuous distance — so "top-tail sharpness" is systematically higher
   for self regardless of regime. The cue measures score *shape*, not regime. A
   label-free regime gate needs a feature orthogonal to the score distribution.

4. **oracle 0.786 leaves +0.10 on the table.** The gap mean→oracle is the
   uncaught regime. A working regime selector would convert it.

## Honest standing

- Simplest robust upgrade over the incumbent: cmp_bg alone (0.650, never <0.49,
  no collapse) — or mean-fusion (0.686). Either beats the current engine and kills
  the below-coin failures on CTX-M/OXA.
- Still per-genome mean AUC ~0.69 — useful screening, not diagnostics.
- n=18 is small; the +0.06 and the regime split need a larger replicon set to be
  trustworthy. Caveat stands.

## Open levers (next, not done)

- Find a regime feature orthogonal to score shape (e.g. plasmid-level
  compositional heterogeneity, or self vs cmp top-gene *disagreement*) → recover
  the 0.10 oracle gap.
- Expand to 60–100 replicons + per-gene surrogate significance (PRISM ships it).
- Wire mean-fusion into src/ as the default two-channel score. See [[sang2-honest-auc]].
