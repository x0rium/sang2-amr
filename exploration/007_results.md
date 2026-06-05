# 007 — Honest gene-level AUC (reproducible)

Question: on a plasmid carrying resistance, does SANG2's composite score rank
ACQUIRED AMR genes above ordinary genes on the SAME replicon?

Dataset: 18 RefSeq *complete* plasmids, esearch across KPC/NDM/CTX-M/mcr/OXA/qnr.
No cherry-picking. 1527 CDS scored (34 acquired-AMR / 1493 non-AMR).
Method: self-baseline scan, window=2000 step=500 (prod default). Gene score =
max composite of overlapping windows. Run: `python exploration/007_honest_auc.py`.

## Headline numbers

| Metric                | Value         |
|-----------------------|---------------|
| Mean per-genome AUC   | **0.625** (median 0.584, n=18) |
| Pooled AUC            | **0.563**     |
| README claim          | 0.667 (was never reproduced by code) |
| Recall@1              | 6%  (1/18)    |
| Recall@5              | 22% (4/18)    |
| Recall@10             | 44% (8/18)    |

AUC scale: 0.5 random · 0.7 weak · 0.8 useful · 0.9 strong.
**On average the method sits near noise.** The README's 0.667 was, if anything,
optimistic — and had no script behind it until now.

## The real finding: the result is BIMODAL, and the split is explainable

| AUC band       | Plasmids | AMR class / host |
|----------------|----------|------------------|
| **0.85–0.99 (strong)** | KY689634, KY120365, CP049352, CP137857, CP029714, CP183858, KY689635 | mcr-1, NDM, qnr — recent HGT, compositionally foreign |
| 0.55–0.65 (noise) | KY288024, MF190368, CP149133 | mixed |
| **0.15–0.47 (worse than coin)** | all 3 CTX-M (0.45–0.47), both OXA-48 (0.28–0.39), CP086065 KPC (0.15) | CTX-M endemic in E. coli; OXA on adapted backbones |

Pattern: the method is **not an AMR detector — it is a recent-horizontal-transfer
detector** that catches AMR only insofar as the gene is still compositionally
foreign. When the gene has adapted to its host (all three CTX-M plasmids → AUC
~0.46, exactly the limitation the README already confesses), the signal inverts:
AMR genes look *more* native than the average plasmid gene, so AUC < 0.5.

## Implications for "how to grow this"

1. Stop selling a single composite AUC. The honest pitch is narrower and true:
   "flags compositionally-foreign / recently-acquired resistance."
2. The bimodality is a feature if made explicit — predict *which regime* a hit is
   in (host-donor phylo distance) and report calibrated confidence per regime.
3. Compositionally-adapted genes (CTX-M) are unreachable by composition alone →
   this is precisely what Mode C (cohort / MIC correlation) must solve. It is the
   missing half of the method, not an optional extra.
4. Recall@10 = 44% is a weak-but-real screening signal on raw composition. The
   structural filter (ESMFold→Foldseek) is what could lift precision on the
   top-10 — that pipeline, not the raw score, is the defensible product.
