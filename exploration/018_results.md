# 018 — Population structure vs phenotype (is the association causal?)

The main threat to 014-017: a k-mer GWAS with no population-structure correction
can be fooled — if the Resistant isolates are one clone and the Susceptible
another, every clone-specific k-mer associates with R, resistance gene included by
coincidence. Closed here, using cached genomes (no downloads).

Method: MinHash bottom-k (3000) sketch per genome -> pairwise Jaccard distance ->
hierarchical lineages. Then the decisive within-lineage test: in clusters with
BOTH R and S, do the resistance k-mers still separate them?

## Ciprofloxacin (gyrA, n=20+20)

- Mean Jaccard distance: same-phenotype **0.687** vs diff-phenotype **0.711** —
  nearly equal. R and S are genetically INTERMIXED, not separate clones. (A
  confound would show same << diff.)
- 3 phenotype-mixed lineages (k=6): C1 7R/9S, C2 7R/9S, C3 4R/1S.
- **Within-lineage (decisive):** resistance k-mers carried, R vs S:
  - C1: R=2529, S=17 — separates
  - C2: R=2704, S=605 — separates
  - C3: R=3150, S=34 — separates

R isolates from DIFFERENT lineages all carry the resistance k-mers; S from the
same lineages don't. The gyrA association is causal, not a single-clone artifact.

## Gentamicin (AAC(3), n=15+15)

- Mean Jaccard distance: same-phenotype **0.623** vs diff-phenotype **0.660** —
  intermixed, not separate clones.
- 3 phenotype-mixed lineages. Within-lineage resistance k-mers carried, R vs S:
  - C1 (8R/10S): R=1894, S=160 — separates
  - C2 (3R/1S): R=1575, S=242 — separates
  - C5 (3R/2S): R=1799, S=112 — separates

AAC(3) association is causal too. Both claimed results (014-016 acquired gene,
017 SNP) survive the population-structure check.

## Honest scope

- This is a cluster-based within-lineage test, not a full mixed-model GWAS
  (pyseer-style kinship). It directly rules out the dominant failure mode
  (single-clone confound) and is interpretable, but a linear mixed model with a
  kinship matrix would be the rigorous form.
- MinHash Jaccard is a genome-distance proxy, not a phylogeny; lineages are
  approximate. Good enough to show phenotype-mixed clades exist and the signal
  holds within them.
- n=20+20. The within-lineage cell counts are small (e.g. C3 has 1 S); the trend
  is consistent across all three mixed lineages, which is the point.

## Bottom line

The 014-017 reference-free phenotype associations survive a population-structure
check: the resistance signal separates R from S WITHIN lineages, so it is not an
artifact of R and S being different clones. Foundation is solid; safe to build the
phenotype × structural fusion on top.
