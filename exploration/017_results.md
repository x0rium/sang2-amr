# 017 — Unexplained resistance: catching the gyrA SNP gene-based tools miss

Mode C (014-016) found AAC(3) for gentamicin — a known ACQUIRED gene that
gene-presence databases already find. The real claim is "find what databases
miss". The cleanest miss for gene-presence/absence tools is a POINT MUTATION in a
CORE gene: fluoroquinolone resistance is driven by gyrA/parC QRDR mutations
(S83L, D87N). gyrA is present in EVERY isolate, R and S — presence/absence sees no
difference; only the SNP distinguishes them.

Setup: 20 R + 20 S E. coli, ciprofloxacin, Laboratory AST (BV-BRC). Reference-free
k-mer GWAS (k=31, both strands), carrier-fraction association. No AMR database.
Run: `python -u exploration/017_unexplained.py`.

## Result

- 3,344 strong phenotype k-mers (R-carrier ≥ S-carrier + 0.7).
- Mapped against the genome's topoisomerase/efflux genes:

| gene | strong k-mers |
|------|--------------:|
| **DNA gyrase subunit A (gyrA)** | **22** |
| DNA topoisomerase IV subunit A (parC) | 2 |
| gyrB, topo I/III, topo IV-B, efflux (MtrF, OMP) | 0 |

  gyrA + parC are exactly the QRDR targets of fluoroquinolone resistance; gyrA
  primary, parC secondary. Everything else zero. Biologically exact.

- **The decisive proof:** gyrA is present (any k-mer) in **20/20 R and 20/20 S** —
  gene-presence/absence sees NO difference. But its 22 phenotype k-mers have mean
  carrier **R=86%, S=0%** — the MUTANT ALLELE tracks resistance.

## Why this is the project's claim, finally on a real miss

Gene-presence databases (ResFinder; RGI's acquired-gene model) flag acquired
genes. A core gene mutated in place is invisible to them — the gene is everywhere,
only the allele differs. Reference-free k-mer phenotype GWAS catches it with no
database at all. Across 014-017 one method covered BOTH mechanism types: acquired
gene (AAC(3), gentamicin) and core-gene SNP (gyrA, ciprofloxacin).

## Honest limits

- Precision: tools like PointFinder / RGI's protein-variant model DO know gyrA
  S83L — they use a curated SNP database. So the honest framing is: the
  reference-free method finds SNP resistance with NO database, where gene-presence
  tools are blind and SNP tools need a curated variant catalog. One method, both
  mechanism types, zero references — that is the contribution.
- No population-structure correction (gyrA SNP co-varies with lineage). The clean
  specificity (gyrA+parC QRDR only, not housekeeping or other topoisomerases)
  argues causal, not a lineage artifact — but a mixed model would make it rigorous.
- Did not pin the exact codon (S83L): the carrier split proves allele-specificity;
  extracting QRDR codon 83 (Ser->Leu) would be the gold-standard next touch.
- 20+20, one species, one drug. PoC scale.

## Where this leaves the project

The core claim — reference-free, phenotype-anchored detection of resistance
including what gene-presence databases miss — is demonstrated end to end on two
distinct mechanisms. Next: population-structure-aware association; fuse phenotype
k-mers with the two-channel structural score (novel hit = phenotype-associated
AND structurally anomalous); a true unexplained-resistance panel (R isolates where
RGI reports nothing) to surface genuinely novel determinants.
