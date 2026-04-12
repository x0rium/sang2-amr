# SANG2-AMR: Ab Initio Antimicrobial Resistance Detection

Reference-free detection of AMR genes using statistical signal analysis.
No databases. No training. No BLAST. Just composition, symmetry, and protein motifs.

## What it does

SANG2-AMR scans bacterial genomes and flags regions that look like AMR genes
based on 12 statistical signals — without using CARD, ResFinder, or any
reference database. It finds what reference-based tools miss.

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

## Validation Results

Tested on 6 real AMR genes across 5 plasmids, zero prior knowledge:

| Gene | Mechanism | Plasmid | Rank | Found? |
|------|-----------|---------|------|--------|
| blaKPC-2 | Carbapenemase | pKpQIL (20 kb) | 2/16 | Yes |
| blaNDM-1 | Metallo-β-lactamase | pNDM-KN (267 kb) | 1/237 | Yes |
| blaCTX-M-15 | ESBL | pEC_Bactec (3.9 kb) | 1/5 | Yes |
| mcr-1 | Colistin resistance | pHNSHP45 (64 kb) | 2/55 | Yes |
| vanA | Vancomycin resistance | Tn1546 (10.9 kb) | 1/9 | Yes |
| blaCTX-M-15 | ESBL | pNDM-KN (267 kb) | 24/237 | No* |

\* Gene compositionally similar to host — requires cross-reference (Mode B/C) for detection.

**Novel candidate pipeline tested end-to-end:**
On K. pneumoniae HU105 draft genome (5.1 Mb, 36 contigs) with same-species baseline:
SANG2 flagged 5 hypothetical proteins (AMR-motif similarity 0.31–0.57) →
BLAST rejected 4/5 (known non-AMR) → ESMFold confirmed 1/5 as well-folded (pLDDT=76) →
Foldseek identified tegument-like fold (not AMR).
Result: 0 confirmed novel AMR. Honest — the genome had no unexplained resistance.
The pipeline (SANG2 → BLAST → structure) works as designed; needs genomes with
phenotypic resistance but no known genetic explanation.

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
```

**Important:** Use a same-species reference genome as `--host` for best results.
Cross-species baseline (e.g., E. coli for K. pneumoniae) will flag species
differences, not AMR genes.

### Mode B: Metagenome (validated)

Cluster contigs by signal profile, find AMR per cluster:

```bash
python src/cli.py metagenome assembled_contigs.fasta --clusters 5
```

Validated: ARI=0.875 on 4-species mock metagenome.

### Mode C: Cohort (not implemented)

Correlate signals with MIC phenotype across isolates. Stub — requires
validation on PATRIC/BV-BRC data.

## Project structure

```
amr/
├── src/                          # Production code (1,800+ lines)
│   ├── ir.py                     # Data classes
│   ├── signals.py                # 12 universal signal functions
│   ├── genome_signals.py         # Genomic-specific signals
│   ├── anomaly.py                # Composite scorer (CORE)
│   ├── aa_resonance.py           # Protein-level AMR detection
│   ├── resonance.py              # 6-signal pair resonance engine
│   ├── context.py                # FSM context boost
│   ├── merge.py                  # BPE pair merge hierarchy
│   ├── cluster.py                # Metagenomic binning
│   ├── pipeline.py               # Mode A/B/C orchestration
│   ├── report.py                 # Text/JSON/TSV output
│   └── cli.py                    # Command-line interface
├── exploration/                  # Validation experiments (001-006)
│   ├── 001_alphabet_signals.py   # Signal validation on E. coli K-12
│   ├── 002_plasmid_contrast.py   # Chromosome vs plasmid discrimination
│   ├── 002b_genomic_signals.py   # GC-skew, Chargaff-2, ρ, RSCU, IR
│   ├── 003b_batch_pair_merge.py  # BPE hierarchy construction
│   ├── 004_resonance_motifs.py   # Codon-level resonance (ATG = magic bytes)
│   ├── 005_amr_detection.py      # Ab initio KPC detection (rank 8/37)
│   ├── 006_mock_metagenome.py    # 4-species clustering (ARI=0.875)
│   └── 003_tools.md              # Complete experiment log
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

- **Composition-similar genes undetectable in Mode A.** Genes long adapted to
  the host (like CTX-M-15 in Enterobacteriaceae) are invisible by composition.
  Requires cross-reference (Mode B/C) — not yet implemented.
- **False positive rate ~12%** at threshold 0.8 (AUC=0.667 vs self-genome).
  Suitable for screening, not diagnostics.
- **AA-resonance profile covers β-lactamases only.** Other AMR families
  (aminoglycosides, tetracyclines, efflux) need additional profiles.
  Profile generates false positives on ATPases, recombinases, acyltransferases
  (proteins sharing general motifs with β-lactamases).
- **Same-species baseline required.** Cross-species comparison flags species
  differences, not AMR. Demonstrated: E. coli baseline on K. pneumoniae →
  top hits are species-specific genes, not AMR.
- **Not validated on unexplained resistance.** The key experiment — finding
  novel AMR in genomes where RGI fails — requires phenotypically resistant
  isolates with no known genetic explanation (PATRIC/BV-BRC data).

## Upstream

Algorithmic core: [auto-reverser](../auto-reverser/) SANG2 framework.
