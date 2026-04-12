# Verified APIs and Tools

Status: pre-exploration phase. Записи будут добавляться по мере проверки.

---

## Биоинформатические инструменты (для preprocessing)

### fastp
Service: fastp v0.23.4
Req→Resp: `fastp -i reads_R1.fq.gz -I reads_R2.fq.gz -o clean_R1.fq.gz -O clean_R2.fq.gz`
Limits: Single-threaded by default, ~10-50 Mb/s throughput
Sources: https://github.com/OpenGene/fastp
Status: UNRUNNABLE (not installed, to be verified in exploration/001)

### SPAdes (genome assembler)
Service: SPAdes v4.0.0
Req→Resp: `spades.py --isolate -1 clean_R1.fq.gz -2 clean_R2.fq.gz -o assembly/`
Limits: RAM-hungry (8-64 GB for bacterial genomes), hours for metagenomes
Sources: https://github.com/ablab/spades
Status: UNRUNNABLE (not installed, to be verified in exploration/001)

### Prodigal (ORF prediction)
Service: Prodigal v2.6.3 / Pyrodigal v3.x
Req→Resp: `prodigal -i contigs.fasta -o genes.gff -a proteins.faa -p meta`
Limits: No eukaryotic genes; -p meta for metagenomes
Sources: https://github.com/hyattpd/Prodigal
Status: UNRUNNABLE (not installed, to be verified in exploration/001)

### RGI (Resistance Gene Identifier)
Service: RGI v6.0.3 + CARD v3.2.9
Req→Resp: `rgi main -i contigs.fasta -o rgi_output -t contig -a BLAST --clean`
Limits: CARD database must be loaded first (`rgi load`); BLAST can be slow
Sources: https://github.com/arpcard/rgi, https://card.mcmaster.ca
Status: UNRUNNABLE (not installed, to be verified in exploration/001)

---

## Данные

### NCBI datasets CLI
Service: NCBI datasets v16.x
Req→Resp: `datasets download genome accession GCF_000005845.2 --include genome`
Limits: Rate limits apply; large downloads may timeout
Sources: https://www.ncbi.nlm.nih.gov/datasets/
Status: UNRUNNABLE (not installed, to be verified in exploration/001)

### SRA Toolkit
Service: SRA Toolkit v3.x
Req→Resp: `fasterq-dump SRR2726667 --split-files`
Limits: Downloads can be large (10+ GB for metagenomes)
Sources: https://github.com/ncbi/sra-tools
Status: UNRUNNABLE (not installed, to be verified in exploration/001)

---

## Python библиотеки (для ядра)

### BioPython
Service: BioPython v1.87
Req→Resp: `SeqIO.read('file.fasta', 'fasta')` → SeqRecord; `Entrez.efetch(db='nucleotide', id='U00096.3')` → FASTA
Limits: Slow for large files; consider pysam for BAM. Entrez requires email.
Sources: https://biopython.org
Status: verified (exploration/001, 2026-04-12) — FASTA/GenBank parsing + NCBI download

### NumPy (from auto-reverser)
Service: NumPy v1.26+ / v2.x
Req→Resp: Vectorized signal computation (reuse from auto-reverser)
Sources: inherited from auto-reverser
Status: verified in auto-reverser context

---

## Upstream: auto-reverser SANG2 core

Service: auto-reverser/src/ (signals.py, resonance.py, merge.py, cluster.py, xref.py)
Req→Resp: Internal library, imported as Python modules
Limits: Byte-alphabet assumption in some functions (needs parameterization)
Sources: ../auto-reverser/src/
Status: verified (production code in auto-reverser)
