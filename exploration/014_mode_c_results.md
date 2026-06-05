# 014-016 — Mode C: reference-free, phenotype-anchored AMR detection

The honest test of the project's core claim. 007 showed composition benchmarks
are circular (you need a database to label "AMR", so you can't test "find what
databases miss" that way). The only escape is phenotype: lab-measured R/S and ask
which genetic features track it, no AMR database in the loop.

Scripts: 014_mode_c.py (cohort + GWAS), 015_mode_c_validate.py (localize),
016_mode_c_proof.py (coordinate-free identity). Data in data/mode_c/.

## Setup

- Cohort: 15 Resistant + 15 Susceptible E. coli, gentamicin, **Laboratory Method**
  AST from BV-BRC (not computational predictions — that would be circular again).
- Method: canonical k=31, both strands, presence/absence per genome. Association =
  R-carrier fraction − S-carrier fraction. No AMR database used to detect.

## Result

- 1,032 k-mers carried by ≥80% more Resistant than Susceptible genomes (top tier:
  R=14/15, S=0/15).
- Localized to a single 951 bp / 297 aa ORF.
- **Coordinate-free proof (016):** 794 of the 1,032 phenotype k-mers (77%) fall
  inside one annotated gene — **AAC(3)-II, Aminoglycoside N(3)-acetyltransferase**,
  the enzyme that acetylates gentamicin. It holds 48% of that gene's k-mers.
- **Specificity:** the genome's OTHER aminoglycoside genes — AAC(6')-Ib, APH(3''),
  APH(6), ANT(3'') — show 0% overlap. The signal picked the gentamicin-specific
  determinant, not aminoglycoside genes in general. Not a lineage/co-selection
  artifact: a confounder would light up housekeeping or all co-located markers.

## Verdict

Mode C works end to end: lab phenotype in → resistance gene out, reference-free.
The project's central claim, untestable by composition (007), tests TRUE here.

## Honest limits

- AAC(3) is a KNOWN gene. We found it without a database, which is the necessary
  mechanism check. The full claim ("find UNKNOWN resistance") needs the next step:
  run this on isolates that are phenotypically R but where RGI/CARD finds no known
  gene (unexplained resistance) — the residual phenotype-associated k-mers that
  DON'T map to a known gene are the novel candidates.
- k-mer phenotype GWAS is an established idea (pyseer/DBGWAS). The contribution
  here is a working reference-free Mode C inside this project, not a new algorithm.
- N=15+15, one species, one drug, simple carrier-fraction (no lineage/population-
  structure correction). PoC scale. The clean AAC(3) specificity is reassuring but
  proper GWAS would add a mixed model + more isolates.
- Next: (a) population-structure-aware association; (b) the unexplained-resistance
  run that tests novelty directly; (c) fuse with the two-channel structural score
  so novel hits are both phenotype-associated AND structurally anomalous.
