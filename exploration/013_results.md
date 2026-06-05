# 013 — Regime selector from plasmid heterogeneity (recover the oracle gap)

012 left an oracle gap: mean-fusion 0.707 vs oracle 0.784. A per-genome regime
selector could recover it IF we can tell, label-free, which channel to trust.
Cue tested: compositional HETEROGENEITY of the plasmid (spread of GC / 3-mer
entropy across windows), computed from sequence, orthogonal to score shape.
76 plasmids, cached scores. Run: `python -u exploration/013_regime_selector.py`.

## Finding 1 — heterogeneity predicts the regime (but opposite to my guess)

| feature      | Spearman r vs (self_auc − cmp_auc) | p       |
|--------------|-----------------------------------:|--------:|
| het_h3_std   | **−0.451**                         | <0.0001 |
| het_gc_std   | −0.365                             | 0.001   |
| het_gc_tail  | +0.336                             | 0.003   |

Sign is OPPOSITE to the "mcr has an island → more heterogeneous → trust self"
hypothesis. Measured: het_gc_std mcr=0.036 vs adapted=0.050 — ADAPTED plasmids
(CTX-M/OXA on big 106–110 kb mosaic replicons, many integrons) are MORE
heterogeneous; mcr sits on compact homogeneous IncX/IncHI. High heterogeneity →
cmp(PRISM) wins. The mechanistic story was wrong; the correlation is real.

## Finding 2 — hard selection fails, SOFT fusion wins

| method        | mean AUC (95% CI)      |
|---------------|------------------------|
| self          | 0.650 [0.598, 0.705]   |
| cmp           | 0.650 [0.604, 0.694]   |
| mean-fusion   | 0.707 [0.664, 0.749]   |
| hard selector | 0.709 [0.664, 0.752]   |  ← = mean-fusion, no gain
| **soft-fusion** | **0.743 [0.702, 0.784]** |  ← recovers ~half the gap
| oracle        | 0.784 [0.752, 0.817]   |

Hard selection only picks the oracle-best channel 63% of the time; at r=0.45 each
wrong switch costs the whole worse channel, so it can't beat blind averaging.
Soft-fusion weights the channels continuously by heterogeneity rank: when the cue
is extreme it gives near-full weight to the right channel, when ambiguous it
relaxes to ~50/50 (= mean-fusion). It never loses everything, and closes the gap
to oracle from 0.077 to 0.041.

## Standing — final detector of this arc

soft-fusion(self, PRISM cmp; weight by plasmid heterogeneity) = **0.743**, up from
the 0.650 incumbent and the fictional 0.667 README claim, and it no longer
collapses below coin on adapted CTX-M/OXA. Honest limits: still ~0.74 = strong
screening, not diagnostics; the remaining 0.04 to oracle needs either a sharper
regime cue or per-region (not per-genome) weighting; weighting uses a transductive
heterogeneity rank over the set (note for productionizing — compute the rank vs a
fixed reference distribution). See [[sang2-honest-auc]].

## The full arc (007→013)

0.667 (claimed, no code) → 0.563 pooled / 0.625 per-genome (honest, 007) →
composition capped, reference swap = net zero (009) → PRISM breaks the blind spot
(010) → two-channel mean-fusion 0.707 stable on 76 (011-012) → soft-fusion 0.743,
near oracle (013).
