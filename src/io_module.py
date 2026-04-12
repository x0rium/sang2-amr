"""I/O: FASTA/FASTQ/GenBank parsing → IR objects.

Converts BioPython records to SANG2-AMR internal representation.
"""

from __future__ import annotations

import numpy as np
from pathlib import Path

from Bio import SeqIO

from ir import Contig, GenomicRegion
from genome_signals import NT_MAP, seq_to_array


def load_fasta(path: str | Path, source: str = "") -> list[Contig]:
    """Load FASTA file → list of Contig objects."""
    contigs = []
    for rec in SeqIO.parse(str(path), "fasta"):
        arr = seq_to_array(str(rec.seq).upper())
        contigs.append(Contig(
            id=rec.id,
            sequence=arr.tobytes(),
            quality=None,
            source=source or rec.id,
        ))
    return contigs


def load_fastq(path: str | Path, source: str = "") -> list[Contig]:
    """Load FASTQ file → list of Contig objects with quality."""
    contigs = []
    for rec in SeqIO.parse(str(path), "fastq"):
        arr = seq_to_array(str(rec.seq).upper())
        qual = np.array(rec.letter_annotations.get("phred_quality", []),
                        dtype=np.uint8)
        contigs.append(Contig(
            id=rec.id,
            sequence=arr.tobytes(),
            quality=qual.tobytes() if len(qual) > 0 else None,
            source=source or rec.id,
        ))
    return contigs


def load_genbank_regions(path: str | Path) -> list[GenomicRegion]:
    """Load GenBank annotation → list of GenomicRegion (CDS, rRNA, tRNA)."""
    rec = SeqIO.read(str(path), "genbank")
    regions = []
    for feat in rec.features:
        if feat.type not in ("CDS", "rRNA", "tRNA", "mobile_element"):
            continue
        loc = feat.location
        if loc is None:
            continue
        strand = "+" if (loc.strand is None or loc.strand == 1) else "-"
        kind = "orf" if feat.type == "CDS" else feat.type
        regions.append(GenomicRegion(
            contig_id=rec.id,
            start=int(loc.start),
            end=int(loc.end),
            strand=strand,
            kind=kind,
        ))
    return regions


def contig_to_array(contig: Contig) -> np.ndarray:
    """Extract uint8 numpy array from Contig.sequence bytes."""
    return np.frombuffer(contig.sequence, dtype=np.uint8)
