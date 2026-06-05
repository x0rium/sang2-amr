# 008 — AUC vs host-donor compositional distance

Hypothesis (from 007's bimodality): SANG2 works exactly when AMR genes are
compositionally FOREIGN relative to the background the baseline is built from.
If true, a cheap distance proxy predicts the per-genome AUC → confidence gate.

Run: `python exploration/008_auc_vs_distance.py` (reuses 007's fetched plasmids).
n = 18 replicons. Distance = AMR-gene vs non-AMR-gene separation in each signal,
background = whole plasmid.

## Correlation with per-genome AUC

| Predictor                          | Spearman r | p     | verdict |
|------------------------------------|-----------:|------:|---------|
| **dinucleotide (rho) dist from bg**| **+0.589** | 0.010 | **significant** |
| GC distance from bg                | +0.494     | 0.037 | significant |
| codon-usage (RSCU) dist from bg    | +0.234     | 0.349 | not sig. |
| \|AMR GC − host genus GC\|          | −0.363     | 0.138 | not sig. |

## Regime gate (median split on rho-distance)

| Regime                          | mean AUC | n |
|---------------------------------|---------:|---|
| AMR more foreign than bg        | **0.757**| 9 |
| AMR not more foreign than bg    | **0.494**| 9 |

## Conclusions (honest)

1. **Hypothesis confirmed in direction, moderate in strength.** The more foreign
   the AMR genes are (dinucleotide / GC) from the plasmid background, the higher
   the AUC. rho gives Spearman +0.59, p=0.01 on n=18. This is a real, free,
   interpretable confidence gate: high distance → trust the hit (AUC ~0.76+);
   low distance → composition is blind, escalate to Mode C / structure.

2. **The gate is noisy — it is not a clean separator.** 2–3 of 18 are
   misclassified (e.g. NZ_CP029714 Pseudomonas-KPC: AUC 0.902 yet rho-distance
   negative — something other than composition is carrying it there). So r=0.59,
   not 0.9. Honest framing: "regime indicator", not "oracle".

3. **Codon-usage distance does NOT predict success (r=0.23, p=0.35).** Yet the
   engine spends two of its tracks (cpd, rscu_d, codon_mi) on codon-level signals.
   Actionable: those tracks may be adding noise, not discrimination — worth an
   ablation (rerun 007 AUC with codon signals removed and see if it moves).

4. **The genus-GC proxy is the wrong reference** (negative, n.s.) — because the
   baseline is the plasmid background, not the host chromosome. Confirms the unit
   of "foreignness" the method actually measures is plasmid-internal.

## Next lever (from here)

- Cheap: ablation — does dropping the codon tracks change the 0.563 pooled AUC?
  If not, simplify the engine and lean on rho+GC + the gate.
- Real: the low-distance regime (AUC 0.49) is unreachable by composition by
  construction → Mode C (cohort/MIC) is the only way in. The gate tells you
  *which genomes* need it.
