# 009 — self-baseline vs host-chromosome baseline (method change, not weights)

007/008 showed self-baseline (background = the plasmid) goes blind on adapted
AMR genes. Hypothesis: build the baseline from the HOST CHROMOSOME instead, so
foreignness is measured vs "normal for the species". Same 18 plasmids, same AUC.

One RefSeq chromosome per genus fetched (E. coli, K. pneumoniae, P. aeruginosa,
S. enterica, E. hormaechei, Citrobacter). Run: `python exploration/009_host_chromosome_baseline.py`.

## Result: net zero in the mean, a regime reshuffle in detail

| group                       | self  | host-chrom | delta |
|-----------------------------|------:|-----------:|------:|
| Mean AUC (n=18)             | 0.625 | 0.627      | +0.001|
| Median                      | 0.584 | 0.609      | +0.025|
| **Adapted (CTX-M/OXA, n=6)**| 0.449 | **0.520**  | +0.071|

Better on 9/18, worse on 7. The change HELPS adapted genes (OXA KX636096
0.39→0.65, two CTX-M up to ~0.56) but BREAKS recent-HGT (mcr-1 KY689634
0.99→0.69, NDM CP049352 0.93→0.75) by the same amount.

## Why (mechanism)

With a host-chromosome baseline, the WHOLE plasmid is foreign to the chromosome,
so the entire plasmid background scores anomalous — the AMR gene no longer stands
out *within* it. self-baseline is right for recent HGT (AMR pops out of an
adapted plasmid); host-baseline is right for adapted genes. Same bimodality,
now confirmed from the method side.

## The hard conclusion

Single-genome composition is capped at ~0.6 AUC and **switching the reference
does not break the cap** — one genome lacks the information to separate "AMR"
from "other accessory/mobile material". The most obvious alternative single-genome
method was tried and returned zero net. To go further you must add a NEW
information source, not rearrange the existing one:

1. **Multi-genome contrast (pan-genome):** N genomes of a species → accessory
   genome (absent from core) → rank anomaly within accessory. AMR lives in
   accessory. Information comes from comparison across genomes. Path to Mode C.
2. **Structural channel (already in `structural.py`, ESMFold→Foldseek):** an
   orthogonal axis (3D fold, not composition) to filter the top-k.
