# SANG2-AMR: Ab Initio Antimicrobial Resistance Detection

Reference-free detection of AMR genes using statistical signal analysis.
No databases. No training. No BLAST. Just composition, symmetry, and protein motifs.

## What it does

SANG2-AMR detects antimicrobial resistance determinants **without any AMR
database** — no CARD, no ResFinder, no BLAST. Two complementary engines:

- **Mode A — structural scan** (single genome): flags compositionally/structurally
  anomalous regions. Good at recently-acquired, foreign DNA; weaker on host-adapted
  genes (the `--two-channel` PRISM channel recovers much of that gap).
- **Mode C — phenotype-anchored** (cohort with R/S labels): reference-free k-mer
  association that recovers the actual resistance determinant from phenotype alone.

The honest, reproducible validation of both is in `exploration/007-020_results.md`.
Headline: Mode C recovers real determinants of two mechanism types (an acquired
gene and a core-gene SNP) reference-free; Mode A is a useful screen (AUC ~0.74 with
two channels) but not a diagnostic. See **Validation** below for the full, candid
picture including what did NOT work.

```
$ python src/cli.py scan --host kp_reference.fasta target_genome.fasta

============================================================
  SANG2-AMR Results
============================================================
  #  1  contig_1:12000-14000 (2000 bp)
        confidence: 0.59  novelty: novel
        evidence: h3: z=8.6, gc: z=4.0, gradient: z=3.2

  #  2  contig_1:7500-9500 (2000 bp)
        confidence: 0.30  novelty: novel
        evidence: ac3: z=2.5, c2: z=1.6, cpd: z=2.0
```

## Validation (honest and reproducible)

Every number below is produced by a script in `exploration/` (007-020), each with
a `*_results.md`. No cherry-picking — datasets are pulled fresh from NCBI/BV-BRC.

### Mode A — structural scan (gene-ranking on a replicon)

Benchmark: rank acquired-AMR genes above ordinary genes on the SAME replicon,
across 76 RefSeq AMR plasmids, 8 antibiotic classes (`007`, `010-013`).

| scorer | mean per-genome AUC (95% CI) |
|--------|------------------------------|
| composition only (the original engine) | 0.65 [0.60, 0.70] |
| **+ PRISM structural channel (`--two-channel`, soft-fused)** | **0.74 [0.70, 0.78]** |

Composition alone caps near AUC 0.6 and **inverts below 0.5 on host-adapted genes**
(all three CTX-M plasmids ~0.46). The PRISM structural channel lifts adapted
CTX-M from ~0.46 to ~0.77 and never collapses below coin. This is a useful
**screen, not a diagnostic**. The README's earlier "AUC=0.667" was never backed by
code; the real number is the one above.

### Mode C — phenotype-anchored, reference-free (the core claim)

Cohorts of lab-AST Resistant + Susceptible E. coli (BV-BRC), reference-free k=31
k-mer association, NO AMR database in detection (`014-018`):

| antibiotic | mechanism found | evidence |
|------------|-----------------|----------|
| gentamicin | **AAC(3)** acetyltransferase | 794/1032 phenotype k-mers in one gene; other aminoglycoside genes 0% |
| ciprofloxacin | **gyrA** QRDR SNP | signal in gyrA (present in 20/20 R *and* S — gene-presence tools are blind; only the mutant allele tracks R: carrier R=86%/S=0%) |

Both associations are **causal, not lineage artifacts**: they separate R from S
*within* phenotype-mixed lineages (`018`, MinHash population-structure check). One
reference-free method covered both an **acquired gene** and a **core-gene SNP** —
the latter is exactly what gene-presence databases (ResFinder) structurally miss.

### What did NOT work (kept honest)

- **phenotype × structural fusion does not compose** (`019`): the structural axis at
  gene level conflates compositionally-extreme *native* genes (ribosomal proteins
  score higher than the real resistance gene) with foreign ones. The defensible
  contribution is the phenotype axis alone.
- **Residual ≠ novelty** (`020`): "unexplained" phenotype k-mers are inflated by
  unannotated-but-known SNP mechanisms, co-selected linked markers, and lineage
  leakage. Discovering genuinely NEW determinants is **not** demonstrated here — it
  needs a curated genotype-negative panel and remains future work.

## How it works

Three levels of detection — each catches what the others miss:

**Level 1: Composition anomaly** — Is this DNA foreign?
12 signals (GC, entropy, Hurst, Chargaff-2, dinucleotide odds, codon usage, ...)
compared to host baseline via z-scores. Foreign DNA stands out.

**Level 2: Codon-level resonance** — Does this region have unusual codon-pair patterns?
Jensen-Shannon divergence of codon-pair frequencies, RSCU deviation.
Catches genes with different codon usage even when GC matches the host.

**Level 3: Protein function** — Does this ORF encode something AMR-like?
Amino acid pair resonance profile built from 114 β-lactamases.
Finds S-X-X-K (serine active site), S-X-N (SDN motif), H-X-X-X-D (zinc binding)
ab initio — without sequence alignment.

Scoring: **rank fusion** across all signals — each signal votes independently,
best rank wins. A gene only needs ONE strong signal to be flagged.

## Key principle: zero magic numbers

Every threshold is derived from data:
- **Otsu:** binary split on any continuous distribution
- **SDR:** find symbols carrying ~50% of mass
- **Z-scores:** deviation from host baseline (mean ± std)
- **Rank fusion:** no weights, no hyperparameters

## Installation

```bash
# Dependencies
pip install numpy scipy biopython scikit-learn

# Run
python src/cli.py scan --host reference.fasta target.fasta
python src/cli.py scan --host reference.fasta target.fasta --format tsv -o results.tsv
```

Requires Python 3.10+, numpy, scipy, biopython, scikit-learn.

## Usage

### Mode A: Single isolate (validated)

Scan a genome or plasmid for anomalous regions:

```bash
python src/cli.py scan target.fasta                     # self-comparison
python src/cli.py scan --host reference.fasta target.fasta  # vs host baseline
python src/cli.py scan --two-channel target.fasta       # + PRISM structural channel
```

**Two-channel mode (`--two-channel`)** adds a structural channel (PRISM/`sang`)
that reaches host-adapted resistance genes composition alone misses. On a
benchmark of 76 RefSeq AMR plasmids it lifts mean per-genome AUC from 0.65
(composition) to 0.74, and stops the below-coin collapse on adapted CTX-M/OXA
(see `exploration/007-013_results.md`). It soft-fuses composition with PRISM,
weighting the structural channel by the replicon's compositional heterogeneity.
Requires the `sang` library; falls back to composition-only if unavailable.

**Important:** Use a same-species reference genome as `--host` for best results.
Cross-species baseline (e.g., E. coli for K. pneumoniae) will flag species
differences, not AMR genes.

### Mode B: Metagenome (validated)

Cluster contigs by signal profile, find AMR per cluster:

```bash
python src/cli.py metagenome assembled_contigs.fasta --clusters 5
```

Validated: ARI=0.875 on 4-species mock metagenome.

### Mode C: Cohort — phenotype-anchored (validated, reference-free)

Given a set of isolates with R/S labels for an antibiotic, find the genetic
determinant by reference-free k-mer association — no AMR database. Validated on
gentamicin (→ AAC(3)) and ciprofloxacin (→ gyrA SNP), population-structure
controlled (`exploration/014-018`). The runnable proof-of-concept lives in
`exploration/014_mode_c.py` … `018_popstructure.py` (pulls cohorts from BV-BRC);
it is not yet wired into `cli.py`.

## Project structure

```
amr/
├── src/                          # Production code (1,800+ lines)
│   ├── ir.py                     # Data classes
│   ├── signals.py                # 12 universal signal functions
│   ├── genome_signals.py         # Genomic-specific signals
│   ├── anomaly.py                # Composite scorer (CORE)
│   ├── resonance.py              # 6-signal pair resonance engine
│   ├── context.py                # FSM context boost
│   ├── merge.py                  # BPE pair merge hierarchy
│   ├── cluster.py                # Metagenomic binning
│   ├── structural.py             # ESMFold + Foldseek structural verification
│   ├── prism_channel.py          # PRISM structural channel (optional sang import)
│   ├── two_channel.py            # composition + PRISM soft-fusion scorer
│   ├── pipeline.py               # Mode A/B/C orchestration
│   ├── report.py                 # Text/JSON/TSV output
│   └── cli.py                    # Command-line interface (--two-channel)
├── exploration/                  # Reproducible experiments + *_results.md
│   ├── 001-006_*.py              # Signal validation, plasmid contrast, metagenome
│   ├── 007_honest_auc.py         # Honest gene-level AUC benchmark (~0.6)
│   ├── 008_auc_vs_distance.py    # AUC tracks host-donor compositional distance
│   ├── 009_host_chromosome_baseline.py  # reference swap = net zero (cap is real)
│   ├── 010_prism_vs_engine.py    # PRISM breaks the adapted-gene blind spot
│   ├── 011_two_channel.py        # composition + PRISM fusion
│   ├── 012_consolidate.py        # 76-plasmid consolidation + bootstrap CI (0.74)
│   ├── 013_regime_selector.py    # soft-fusion by heterogeneity (0.743)
│   ├── 014_mode_c.py             # reference-free phenotype k-mer GWAS
│   ├── 015-016_*.py              # localize + prove: AAC(3) for gentamicin
│   ├── 017_unexplained.py        # gyrA SNP — what gene-presence tools miss
│   ├── 018_popstructure.py       # causal, not lineage (MinHash check)
│   ├── 019_fusion.py             # phenotype×structural — honest negative
│   └── 020_residual_novelty.py   # residual ≠ novelty — guardrail
├── docs/
│   ├── architecture.md           # Pipeline design (L0-L5)
│   ├── algorithms.md             # Signal adaptations from auto-reverser
│   ├── data-sources.md           # Datasets and validation strategy
│   ├── domain-glossary.md        # AMR/genomics terminology
│   ├── commercial.md             # Go-to-market strategy
│   └── gotchas.md                # 18 documented pitfalls
└── prototypes/
```

## What transfers from auto-reverser (and what doesn't)

SANG2-AMR adapts the [auto-reverser](../auto-reverser/) protocol analysis
framework to genomic sequences. Key findings:

| Algorithm | Status | Adaptation |
|-----------|--------|------------|
| Entropy | Works | Use 3-mer level (not nucleotide) |
| SDR | Works | On k-mers |
| Otsu | Works | No changes |
| ΔH (bidirectional) | Works | Finds gene boundaries |
| Autocorrelation | Works (strongest) | 97% CDS have lag-3 peak |
| Hurst | Works | Distinguishes chr/plasmid |
| Resonance (6 signals) | Works at codon level | ATG = magic bytes (R=0.53) |
| **MI half-life** | **Broken** | Replaced by autocorrelation |
| **H/H_max regime** | **Broken** | Replaced by 3-mer entropy |
| **PPMI pair merge** | **Broken** | Replaced by frequency BPE |

**Root cause of failures:** formulas with P(x) in the denominator produce
extreme values when alphabet size ≤ 20 (4 nucleotides vs 256 bytes).

## Discovered signals (not in auto-reverser)

| Signal | What it detects | Validation |
|--------|----------------|------------|
| Chargaff-2 deviation | AMR genes (#1 single detector, KPC rank 1/37) | exploration/005 |
| GC-skew | Replication origin (0.1% accuracy) | exploration/002b |
| Signal gradient | IS-elements (100% recall on 15 transposases) | exploration/002 |
| AA-pair resonance | β-lactamase active sites (S-X-N, H-X-X-X-D) | exploration/validation |

## Limitations

- **Composition-similar genes weak in Mode A (composition-only).** Genes long
  adapted to the host (like CTX-M-15 in Enterobacteriaceae) are near-invisible to
  composition (AUC ~0.46). The `--two-channel` PRISM structural channel recovers
  most of this (adapted AUC ~0.46 → ~0.68); full reach still needs cross-reference
  (Mode C, phenotype-anchored — not yet implemented).
- **Mode A is a screen, not a diagnostic.** Honest composition-only AUC is ~0.6
  (`exploration/007`), ~0.74 with the two-channel structural add-on. Use it to
  shrink a candidate list, not to call resistance.
- **AA-resonance profile covers β-lactamases only.** Other AMR families
  (aminoglycosides, tetracyclines, efflux) need additional profiles.
  Profile generates false positives on ATPases, recombinases, acyltransferases
  (proteins sharing general motifs with β-lactamases).
- **Same-species baseline required.** Cross-species comparison flags species
  differences, not AMR. Demonstrated: E. coli baseline on K. pneumoniae →
  top hits are species-specific genes, not AMR.
- **Novel-determinant discovery is unproven.** Mode C *recovers* known-class
  mechanisms reference-free (validated). Finding genuinely NEW resistance genes —
  in phenotypically-R isolates with no known explanation — is not demonstrated;
  the residual signal is confounded by SNPs, co-selection, and lineage (`020`), and
  needs a curated genotype-negative panel.

## Upstream

Algorithmic core: [auto-reverser](../auto-reverser/) SANG2 framework.
