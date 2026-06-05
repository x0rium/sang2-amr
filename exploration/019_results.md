# 019 — phenotype × structural fusion: a HONEST NEGATIVE

Hypothesis: a real acquired resistance gene should be BOTH phenotype-associated
(Mode C k-mer GWAS) AND structurally anomalous (PRISM), while a core-gene SNP
(gyrA) is phenotype-associated but structurally native. If true, the structural
channel classifies the mechanism type — a fingerprint pyseer can't produce.

Test: on a gentamicin Resistant genome (562.42768), score a CDS panel (AMR +
housekeeping + core gyr) by phenotype (strong-kmer fraction) and structural (PRISM
profile distance gene vs whole genome). Run: `python -u exploration/019_fusion.py`.

## Result — the hypothesis FAILED

| kind  | phenotype | structural | gene |
|-------|----------:|-----------:|------|
| AMR   | 0.48      | 0.75       | AAC(3) acetyltransferase (DUAL) |
| house | 0.00      | **0.84**   | ribosomal protein L31 (higher than AAC(3)!) |
| house | 0.00      | 0.83       | ribosomal protein L28 |
| core  | 0.00      | 0.55       | gyrA |
| ...   | 0.00      | 0.5–0.8    | most genes |

- **Phenotype axis is perfect:** only AAC(3) scores >0 (0.48); everything else
  exactly 0.00. Specific and causal (consistent with 016/018).
- **Structural axis does NOT discriminate.** Nearly every gene scores 0.5–0.84;
  ribosomal housekeeping proteins score HIGHER than the real resistance gene. The
  acquired-vs-native split the hypothesis predicted does not appear.

## Why (two honest reasons)

1. Scale artifact: a 300–900 bp gene's PRISM profile differs from a 5 Mb genome's
   profile because of SIZE, not foreignness. In 010-013 `cmp` worked because it
   compared a gene to a PLASMID (comparable scale) and used a RELATIVE rank within
   the replicon — not an absolute gene-vs-whole-chromosome distance.
2. Deeper: structural anomaly conflates "compositionally extreme" with "foreign".
   Ribosomal proteins are compositionally extreme (codon-optimized, biased AA
   usage) — native but anomalous-looking. This is the classic composition-based
   genomic-island false-positive, which SANG2's own README acknowledges.

## What this means (re-scoping the project, honestly)

- The dual-evidence "fingerprint" does NOT compose naively. The structural axis at
  gene level is dominated by compositionally-extreme NATIVE genes, not acquired
  ones, so it adds no separating power here.
- The real strength of Mode C is the PHENOTYPE axis alone: specific, causal,
  reference-free, covering acquired genes AND core-gene SNPs (014-018).
- two-channel structural (010-013) remains valid IN ITS OWN setting (ranking AMR
  vs non-AMR WITHIN a plasmid, AUC 0.74) — but it is not a universal
  mechanism-classifier, and bolting it onto phenotype GWAS as a second axis was
  not justified by the data.

A clean negative is worth more than a forced story: it says the project's defensible
contribution is the reference-free phenotype detection (Mode C), not a
phenotype×structure fusion.
