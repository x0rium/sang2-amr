"""Command-line interface for SANG2-AMR.

Usage:
  python -m cli scan --host reference.fasta target.fasta
  python -m cli metagenome contigs.fasta --clusters 5
  python -m cli cohort isolates/ --mic mic_data.tsv
"""

from __future__ import annotations

import argparse
import sys
import os

# Ensure src/ is on path when running as script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from io_module import load_fasta
from pipeline import run_single_isolate, run_metagenome
from report import format_text, format_json, format_tsv


def cmd_scan(args):
    """Mode A: Scan a genome/plasmid for AMR candidates."""
    target_contigs = load_fasta(args.target, source="target")
    host_contigs = load_fasta(args.host, source="host") if args.host else None

    candidates = run_single_isolate(
        target_contigs,
        host_contigs=host_contigs,
        window=args.window,
        step=args.step,
        top_n=args.top,
    )

    output = _format(candidates, args.format)
    if args.output:
        with open(args.output, "w") as f:
            f.write(output)
        print(f"Results written to {args.output}", file=sys.stderr)
    else:
        print(output)


def cmd_metagenome(args):
    """Mode B: Cluster contigs and find AMR per cluster."""
    contigs = load_fasta(args.contigs, source="metagenome")

    result = run_metagenome(
        contigs,
        n_clusters=args.clusters,
        min_contig_len=args.min_len,
    )

    n_clusters = max(result["clusters"]) + 1 if result["clusters"] else 0
    print(f"Contigs: {len(result['clusters'])}, Clusters: {n_clusters}", file=sys.stderr)

    for k, candidates in result["candidates_per_cluster"].items():
        if candidates:
            print(f"\n--- Cluster {k} ---")
            print(format_text(candidates, title=f"Cluster {k}"))


def _format(candidates, fmt: str) -> str:
    if fmt == "json":
        return format_json(candidates)
    if fmt == "tsv":
        return format_tsv(candidates)
    return format_text(candidates)


def main():
    parser = argparse.ArgumentParser(
        prog="sang2-amr",
        description="SANG2-AMR: Ab initio antimicrobial resistance detection",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # scan
    p_scan = sub.add_parser("scan", help="Mode A: Scan genome for AMR candidates")
    p_scan.add_argument("target", help="Target FASTA (genome/plasmid to scan)")
    p_scan.add_argument("--host", help="Host reference FASTA (for baseline)")
    p_scan.add_argument("--window", type=int, default=2000, help="Window size (bp)")
    p_scan.add_argument("--step", type=int, default=500, help="Step size (bp)")
    p_scan.add_argument("--top", type=int, default=20, help="Top N candidates")
    p_scan.add_argument("--format", choices=["text", "json", "tsv"], default="text")
    p_scan.add_argument("-o", "--output", help="Output file")
    p_scan.set_defaults(func=cmd_scan)

    # metagenome
    p_meta = sub.add_parser("metagenome", help="Mode B: Metagenomic binning + AMR")
    p_meta.add_argument("contigs", help="Assembled contigs FASTA")
    p_meta.add_argument("--clusters", type=int, default=None, help="Number of clusters (auto if omitted)")
    p_meta.add_argument("--min-len", type=int, default=2000, help="Minimum contig length")
    p_meta.set_defaults(func=cmd_metagenome)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
