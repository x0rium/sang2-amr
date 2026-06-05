# 010 — PRISM/sang vs SANG2 engine (gene-level AUC)

Question: do PRISM's richer, calibrated structural features rank AMR genes above
ordinary genes where SANG2's 12 composition signals go blind (adapted CTX-M/OXA)?

18 plasmids, per-CDS scores, per-genome AUC. PRISM via `/Users/x0rium/CATALYSTO/Prism`.
Scores: self (SANG2 engine), ncd_bg (compression dist gene-vs-plasmid),
cmp_bg (sang.compare(profile(gene), profile(plasmid)).distance over 14 metrics),
anom (max detect_anomalies overlapping gene). Run: `python exploration/010_prism_vs_engine.py`.

## Result: PRISM structural-profile distance reaches the blind spot

| group                     | n  | self  | ncd_bg | **cmp_bg** | anom  |
|---------------------------|----|------:|-------:|-----------:|------:|
| ALL                       | 18 | 0.625 | 0.380  | **0.650**  | 0.497 |
| **ADAPTED (ctxm/oxa)**    | 6  | 0.449 | 0.364  | **0.684**  | 0.483 |
| RECENT (mcr/ndm/kpc/qnr)  | 12 | 0.714 | 0.388  | **0.632**  | 0.504 |

All three CTX-M jumped from ~0.46 to ~0.77 (cmp_bg): 0.467→0.756, 0.462→0.772,
0.451→0.766. OXA CP148085 0.275→0.512. The exact class the README calls invisible
to composition becomes visible to PRISM's 14-metric profile (fractal, recurrence,
permutation entropy, topology) — features GC/codon signals don't carry.

## Two regimes, two channels (the key structure)

- self (composition) is BETTER on recent HGT (0.714 vs 0.632) — a foreign gene
  pops out by raw composition.
- cmp_bg (PRISM profile) is BETTER on adapted genes (0.684 vs 0.449) — when
  composition matches the host, deeper structural metrics still separate it.
- cmp_bg is also ROBUST: never below 0.49; self collapses to 0.15/0.27/0.34 and
  inverts. PRISM doesn't invert.

The channels are anti-correlated by regime → averaging them hurts ((0.625+0.650)/2
≈ 0.64). SELECTING the right channel per genome is what wins:
oracle max(self, cmp_bg) per genome = **mean AUC 0.786** (upper bound, picks best
in hindsight). A real selector exists: the 008 rho-distance gate already predicts
the regime — low distance (adapted) → use cmp_bg, high distance (recent) → use self.

## What does NOT work

- ncd_bg (zlib compression distance): noisy, 0.08–0.91, unreliable.
- detect_anomalies on gene-level: ~0.50, useless — its sparse points don't label
  genes. PRISM's flagship anomaly detector is the wrong tool here; its profile
  *compare* is the right one.

## Caveats

- Still the self-baseline paradigm (distance from plasmid), just with richer
  metrics. Not a new information source — but the richer metrics break the cap
  that GC/codon hit.
- n=6 adapted is small; needs more replicons + surrogate significance per gene.
- ~10–25 s/plasmid (profile per gene). Fine for research, not yet for scale.

## Next lever

Wire the 008 rho-gate to SELECT self vs cmp_bg per region → a single two-channel
detector that holds ~0.7+ across BOTH regimes instead of collapsing on adapted
genes. This is the first method change (007→010) that actually moved the blind
spot. See [[sang2-honest-auc]].
