"""015: Mode C validation — localize phenotype-associated k-mers, identify the gene.

014 found 29,650 k-mers carried by ~all Resistant and NO Susceptible E. coli
(gentamicin), with zero AMR database in the loop. If Mode C works, those k-mers
should localize to a real aminoglycoside-resistance gene.

This script:
  1. Load the GWAS result (gwas.npz) and take the strongly R-associated k-mers.
  2. Localize them in a Resistant genome -> dense clusters = candidate loci
     (contig + coordinates).
  3. VALIDATION ONLY: query BV-BRC genome_feature for CDS overlapping each locus,
     read the product names. Annotation is used to CHECK what we found — never to
     find it. Success = the top locus is an aminoglycoside-modifying gene.

Run: python -u exploration/015_mode_c_validate.py
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

_TAB = np.full(256, 255, np.uint8)
for i, b in enumerate(b"ACGT"):
    _TAB[b] = i
_WEIGHTS = (4 ** np.arange(K - 1, -1, -1, dtype=object)).astype(np.uint64)


def _get(url):
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=90) as r:
            return json.loads(r.read().decode())
    except Exception as e:
        print(f"    GET failed: {e}", flush=True)
        return None


def fwd_hashes(codes_u64, k=K):
    if len(codes_u64) < k:
        return np.empty(0, np.uint64)
    win = np.lib.stride_tricks.sliding_window_view(codes_u64, k)
    return (win * _WEIGHTS).sum(axis=1, dtype=np.uint64)


def localize(fna: Path, strong: np.ndarray):
    """Return list of (contig_id, position, in_strong) hits; positions of strong k-mers."""
    strong_sorted = np.sort(strong)
    hits = {}  # contig -> list of start positions of strong k-mers
    for block in fna.read_text().split(">")[1:]:
        lines = block.splitlines()
        cid = lines[0].split()[0]
        seq = "".join(lines[1:])
        codes = _TAB[np.frombuffer(seq.encode(), np.uint8)]
        valid = codes != 255
        n = len(codes)
        idx = 0
        positions = []
        while idx < n:
            if not valid[idx]:
                idx += 1
                continue
            j = idx
            while j < n and valid[j]:
                j += 1
            run = codes[idx:j].astype(np.uint64)
            if len(run) >= K:
                hs = fwd_hashes(run, K)
                # membership in strong
                pos = np.searchsorted(strong_sorted, hs)
                pos = np.clip(pos, 0, len(strong_sorted) - 1)
                ismem = strong_sorted[pos] == hs
                where = np.where(ismem)[0]
                for w in where:
                    positions.append(idx + int(w))
            idx = j
        if positions:
            hits[cid] = sorted(positions)
    return hits


def cluster_positions(positions, gap=500):
    """Merge nearby strong-kmer positions into loci [start, end, count]."""
    if not positions:
        return []
    loci = []
    s = e = positions[0]
    cnt = 1
    for p in positions[1:]:
        if p - e <= gap:
            e = p
            cnt += 1
        else:
            loci.append([s, e + K, cnt])
            s = e = p
            cnt = 1
    loci.append([s, e + K, cnt])
    return loci


def features_in(genome_id, sequence_id, start, end):
    """BV-BRC annotated CDS overlapping a region (validation only)."""
    q = (f"and(eq(genome_id,{genome_id}),eq(sequence_id,{urllib.parse.quote(sequence_id)}),"
         f"eq(feature_type,CDS),overlaps(location,({start}..{end})))")
    url = f"{API}/genome_feature/?{q}&select(product,start,end,gene)&limit(50)&http_accept=application/json"
    return _get(url) or []


def main():
    print("=" * 72)
    print("  015: Mode C validation — localize + identify")
    print("=" * 72)
    g = np.load(DATA / "gwas.npz", allow_pickle=True)
    allk, cntR, cntS = g["allk"], g["cntR"], g["cntS"]
    nR, nS = int(g["nR"]), int(g["nS"])
    R_ids = [str(x) for x in g["R_ids"]]
    score = cntR / nR - cntS / nS
    strong = allk[score >= 0.8]   # carried by >=80% more R than S
    print(f"\n  strongly R-associated k-mers: {len(strong):,} "
          f"(R-carrier >= S-carrier + 0.8)")

    # localize in the first Resistant genome that has them densely
    print("\n  Localizing in Resistant genomes...")
    best = None
    for gid in R_ids[:5]:
        fna = DATA / f"{gid}.fna"
        if not fna.exists():
            continue
        hits = localize(fna, strong)
        total = sum(len(v) for v in hits.values())
        print(f"    {gid}: {total} strong-kmer hits across {len(hits)} contigs")
        if total and (best is None or total > best[1]):
            best = (gid, total, hits)
    if not best:
        print("  no localization; abort.")
        return

    gid, total, hits = best
    print(f"\n  Using {gid} ({total} hits). Top loci (dense strong-kmer clusters):")
    # rank loci across contigs by count
    all_loci = []
    for cid, positions in hits.items():
        for s, e, cnt in cluster_positions(positions):
            all_loci.append((cnt, cid, s, e))
    all_loci.sort(reverse=True)

    print("\n  Validating top loci against BV-BRC annotation (CHECK only):")
    flagged = []
    AMINO = ("aminoglycoside", "aac", "aph", "aad", "ant", "streptomycin",
             "kanamycin", "gentamicin", "acetyltransferase", "nucleotidyltransferase",
             "phosphotransferase")
    for cnt, cid, s, e in all_loci[:8]:
        feats = features_in(gid, cid, s, e)
        prods = "; ".join(sorted({f.get("product", "?") for f in feats}))[:90]
        hit = any(any(a in (f.get("product", "").lower()) for a in AMINO) for f in feats)
        mark = "  <== AMINOGLYCOSIDE" if hit else ""
        print(f"    {cid}:{s}-{e}  span={e-s}bp  kmers={cnt}  -> {prods}{mark}")
        if hit:
            flagged.append((cid, s, e, prods))

    print("\n" + "=" * 72)
    if flagged:
        print("  RESULT: reference-free phenotype GWAS localized to a real")
        print("  aminoglycoside-resistance gene WITHOUT any AMR database:")
        for cid, s, e, prods in flagged[:5]:
            print(f"    {cid}:{s}-{e}  {prods}")
        print("\n  This is Mode C working end to end: phenotype in -> resistance")
        print("  determinant out, no CARD/ResFinder. The honest core test passes.")
    else:
        print("  RESULT: top loci are not annotated aminoglycoside genes.")
        print("  Either co-selected markers, mobile-element backbone, or a miss.")
        print("  Inspect products above.")


if __name__ == "__main__":
    main()
