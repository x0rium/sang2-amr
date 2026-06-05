# 012 — Consolidation on 76 replicons + bootstrap CI

011 (n=18) found mean-fusion(self, PRISM cmp) = 0.686. Real or fluke? Extended to
76 RefSeq plasmids across 8 AMR classes (kpc/ndm/ctxm/oxa/mcr/qnr/sul/van).
Per-plasmid scores cached in data/gene_scores/. Bootstrap CI over genomes (B=3000).
Run: `python -u exploration/012_consolidate.py`.

## Result: mean-fusion advantage is STABLE

| group   | n  | self                 | cmp (PRISM)          | mean-fusion          |
|---------|----|----------------------|----------------------|----------------------|
| ALL     | 76 | 0.650 [0.596, 0.703] | 0.650 [0.602, 0.694] | **0.707 [0.663, 0.748]** |
| ADAPTED | 26 | 0.585 [0.503, 0.670] | 0.661 [0.610, 0.709] | 0.663 [0.586, 0.741] |
| RECENT  | 50 | 0.684 [0.615, 0.752] | 0.645 [0.576, 0.705] | 0.730 [0.679, 0.780] |

PASS: mean-fusion CI lower bound (0.663) > self point estimate (0.650) — the
two-channel gain holds, not an n=18 artifact (0.707 on 76 vs 0.686 on 18).
ADAPTED mean-fusion CI [0.586, 0.741] is clear of 0.5 — no more below-coin
collapse on the genes the README called invisible.

## Per-class breakdown — where fusion helps and where it hurts

| class | regime | self  | cmp   | mean-fusion |
|-------|--------|------:|------:|------------:|
| mcr   | recent | 0.930 | 0.263 | 0.659  (HURT, −0.27) |
| van   | recent | 0.675 | 0.827 | 0.846  (best) |
| ndm   | recent | 0.671 | 0.732 | 0.758 |
| qnr   | recent | 0.637 | 0.717 | 0.757 |
| sul   | adapt  | 0.616 | 0.658 | 0.682 |
| ctxm  | adapt  | 0.551 | 0.685 | 0.676 |
| kpc   | recent | 0.509 | 0.683 | 0.631 |
| oxa   | adapt  | 0.592 | 0.630 | 0.621 |

## The honest weakness (motivates 013)

Blind averaging HURTS classes where one channel clearly dominates. mcr: self 0.930
(composition nails it), PRISM cmp 0.263 (useless), fusion drags it to 0.659. The
mean wins on average only because adapted + comparable-channel classes outvote the
loss. A regime selector that recognizes "mcr → take self whole" recovers this — it
is exactly the +0.10 oracle gap. The 011 selector failed (measured score shape);
013 needs a feature orthogonal to it, e.g. inter-channel DISAGREEMENT on the top
gene, or plasmid compositional heterogeneity.

## Standing

- Two-channel mean-fusion = 0.707 [0.663, 0.748], a stable, honest upgrade over
  the 0.650 incumbent, that no longer dies on CTX-M/OXA. Screening-grade.
- Next: regime selector (013) to stop fusion from taxing single-channel-dominant
  classes like mcr. See [[sang2-honest-auc]].
