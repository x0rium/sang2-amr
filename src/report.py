"""Output formatting for SANG2-AMR results (L5)."""

from __future__ import annotations

import json

from ir import AMRCandidate, to_json


def format_text(candidates: list[AMRCandidate], title: str = "SANG2-AMR Results") -> str:
    """Human-readable text report."""
    lines = [
        f"{'=' * 60}",
        f"  {title}",
        f"{'=' * 60}",
        f"  Candidates found: {len(candidates)}",
        "",
    ]

    for i, c in enumerate(candidates, 1):
        r = c.region
        lines.append(f"  #{i:3d}  {r.contig_id}:{r.start}-{r.end} ({r.end - r.start} bp)")
        lines.append(f"        confidence: {c.confidence:.2f}  novelty: {c.novelty}")
        if c.nearest_known:
            lines.append(f"        nearest known: {c.nearest_known}")
        if c.evidence:
            lines.append(f"        evidence: {', '.join(c.evidence[:5])}")
        lines.append("")

    return "\n".join(lines)


def format_json(candidates: list[AMRCandidate]) -> str:
    """JSON output."""
    return to_json(candidates)


def format_tsv(candidates: list[AMRCandidate]) -> str:
    """TSV output for downstream tools."""
    header = "rank\tcontig\tstart\tend\tlength\tstrand\tconfidence\tnovelty\tnearest_known\tevidence"
    lines = [header]
    for i, c in enumerate(candidates, 1):
        r = c.region
        evidence_str = "; ".join(c.evidence[:5])
        lines.append(
            f"{i}\t{r.contig_id}\t{r.start}\t{r.end}\t{r.end - r.start}\t"
            f"{r.strand}\t{c.confidence:.3f}\t{c.novelty}\t"
            f"{c.nearest_known or ''}\t{evidence_str}"
        )
    return "\n".join(lines)
