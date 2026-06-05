# 020 — Residual analysis: high residual is NOT novelty (a guardrail result)

Operationalizing "find what databases miss": of the strong phenotype k-mers, how
many sit in the genome's KNOWN resistance genes, and what's the residual?
Run: `python -u exploration/020_residual_novelty.py`.

## Result

| cohort | explained by known AMR genes | residual |
|--------|------------------------------:|---------:|
| gentamicin | 44% (AAC(3) 794, TEM 62, MdtH 22) | 56% |
| ciprofloxacin | 0% | 100% |

## Why the residual is NOT novelty (the honest catch)

Ciprofloxacin shows 0% explained / 100% residual — which looks like "all novel"
but is an artifact:

- **gyrA/parC SNPs are not annotated as "resistance genes".** gyrA is "DNA gyrase
  subunit A", a core gene; the resistance-keyword filter doesn't catch it. So the
  real mechanism (found in 017) lands in the residual, inflating it.
- **Co-selection:** MDR plasmids carry linked bla/sul/tet next to the causal gene;
  their k-mers associate with R without being the cause.
- **Lineage leakage:** the hard carrier-fraction filter doesn't fully remove
  population structure (018 controls the main claim, not every residual k-mer).

So a large residual conflates: unannotated-but-known mechanisms (core-gene SNPs),
co-selected markers, and lineage signal. Residual size alone proves nothing about
novelty.

## What this means for the project

- Where "explained" is high (gentamicin → AAC(3)), the method validly recovers the
  known mechanism reference-free. That part is solid (014-018).
- The novelty claim ("find what databases miss") CANNOT be settled by residual
  mining. It needs a curated genotype-negative panel: isolates phenotypically R
  where NO known mechanism — acquired gene, core-gene SNP, OR co-selection — and no
  lineage confound explains it. The residual k-mers there, localized to an
  uncharacterized ORF and replicated across lineages, would be the real candidate.
  That is a serious study, not a PoC step.

## Honest standing of the whole arc (007-020)

Proven: reference-free, phenotype-anchored detection recovers real resistance
determinants of two mechanism types (acquired AAC(3); core-gene SNP gyrA),
population-structure-controlled (018). NOT proven: discovery of genuinely novel
determinants — the residual is not clean enough, and the structural fusion meant
to disambiguate it does not compose (019). The defensible contribution is
reference-free recovery of known-class mechanisms without a database; true novelty
discovery remains future work with a curated panel.
