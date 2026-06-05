"""016: Mode C final proof — do the phenotype-found k-mers lie inside a real
aminoglycoside-resistance gene? Coordinate-free, the clean test.

014 found 1,032 k-mers carried by ~all Resistant / no Susceptible (gentamicin),
no AMR database used. 015 localized them to a 951 bp ORF (297 aa). Here we PROVE
identity without coordinates: pull the nucleotide sequences of the genome's
annotated aminoglycoside CDS and ask what fraction of our strong k-mers fall
inside each. Annotation is the answer key, never used to find the signal.

Run: python -u exploration/016_mode_c_proof.py
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "exploration" / "data" / "mode_c"
API = "https://www.bv-brc.org/api"
K = 31
GENOME = "562.42768"

_TAB = np.full(256, 255, np.uint8)
for i, b in enumerate(b"ACGT"):
    _TAB[b] = i
_WEIGHTS = (4 ** np.arange(K - 1, -1, -1, dtype=object)).astype(np.uint64)


def _get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def kmer_set(seq: str) -> np.ndarray:
    codes = _TAB[np.frombuffer(seq.upper().encode(), np.uint8)]
    out = []
    valid = codes != 255
    n = len(codes)
    idx = 0
    while idx < n:
        if not valid[idx]:
            idx += 1
            continue
        j = idx
        while j < n and valid[j]:
            j += 1
        run = codes[idx:j].astype(np.uint64)
        for r in (run, (3 - run)[::-1].astype(np.uint64)):
            if len(r) >= K:
                win = np.lib.stride_tricks.sliding_window_view(r, K)
                out.append((win * _WEIGHTS).sum(axis=1, dtype=np.uint64))
        idx = j
    return np.unique(np.concatenate(out)) if out else np.empty(0, np.uint64)


def main():
    print("=" * 74)
    print("  016: Mode C proof — phenotype k-mers vs annotated AMR genes")
    print("=" * 74)

    g = np.load(DATA / "gwas.npz", allow_pickle=True)
    allk, cntR, cntS = g["allk"], g["cntR"], g["cntS"]
    nR, nS = int(g["nR"]), int(g["nS"])
    score = cntR / nR - cntS / nS
    strong = np.sort(allk[score >= 0.8])
    print(f"\n  strong phenotype k-mers (R-carrier >= S+0.8): {len(strong):,}")

    # Pull annotated aminoglycoside CDS coords (answer key only), then cut their
    # sequence from the contig by accession+coords (BV-BRC coords throughout).
    print(f"\n  Fetching annotated aminoglycoside CDS for {GENOME} (validation key)...")
    q = f"and(eq(genome_id,{GENOME}),eq(feature_type,CDS),keyword(aminoglycoside))"
    url = (f"{API}/genome_feature/?{q}"
           f"&select(product,accession,start,end)&limit(200)&http_accept=application/json")
    meta = _get(url)
    # contig sequences by accession
    contig_cache: dict[str, str] = {}

    def contig_seq(acc: str) -> str:
        if acc not in contig_cache:
            u = (f"{API}/genome_sequence/?and(eq(genome_id,{GENOME}),eq(sequence_id,{acc}))"
                 f"&select(sequence)&limit(1)&http_accept=application/json")
            rows = _get(u)
            contig_cache[acc] = rows[0]["sequence"].upper() if rows else ""
        return contig_cache[acc]

    feats = []
    for m in meta:
        acc = m.get("accession")
        s, e = int(m.get("start", 0)), int(m.get("end", 0))
        cseq = contig_seq(acc)
        lo, hi = min(s, e) - 1, max(s, e)
        sub = cseq[max(0, lo):hi]
        if len(sub) >= K:
            feats.append({"product": m.get("product", "?"), "na_sequence": sub})
    print(f"    {len(feats)} aminoglycoside CDS with sequence")

    # For each, fraction of strong phenotype k-mers it contains.
    print("\n  Strong-kmer overlap per gene (how much of the phenotype signal it holds):")
    rows = []
    for f in feats:
        ks = kmer_set(f["na_sequence"])
        if len(ks) == 0:
            continue
        pos = np.searchsorted(strong, ks)
        pos = np.clip(pos, 0, len(strong) - 1)
        inter = int((strong[pos] == ks).sum())
        rows.append((inter, len(ks), f.get("product", "?")))
    rows.sort(reverse=True)
    total_explained = 0
    for inter, nk, prod in rows[:12]:
        frac_gene = inter / nk
        flag = " <==" if inter >= 20 else ""
        print(f"    strong_in_gene={inter:>4}  gene_kmers={nk:>4}  "
              f"({frac_gene:.0%} of gene)  {prod[:60]}{flag}")
        total_explained += inter
    covered = sum(r[0] for r in rows)
    print(f"\n  Strong k-mers explained by annotated resistance genes: "
          f"{covered}/{len(strong)} = {covered/len(strong):.0%}")

    print("\n" + "=" * 74)
    amino = [r for r in rows if any(a in r[2].lower() for a in
             ("aminoglycoside", "aac", "aph", "gentamicin", "acetyltransferase"))]
    if amino and amino[0][0] >= 20:
        inter, nk, prod = amino[0]
        print("  PROOF: the reference-free phenotype GWAS signal sits inside a real")
        print(f"  aminoglycoside-resistance gene — no AMR database used to find it:")
        print(f"    {prod}")
        print(f"    {inter} of our {len(strong)} phenotype k-mers fall in this gene.")
        print("\n  Mode C works end to end: lab phenotype in -> resistance gene out,")
        print("  reference-free. This is the project's core claim, finally tested true.")
    else:
        print("  Top overlap is not an aminoglycoside gene — inspect the list above.")


if __name__ == "__main__":
    main()
