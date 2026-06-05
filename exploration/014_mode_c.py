"""014: Mode C — reference-free, phenotype-anchored AMR detection.

THE honest test of the project's core claim. 007 showed every composition
benchmark is circular: to label a gene "AMR" you need a database, so you can't
test "find what databases miss" by composition-vs-known-genes. The only escape is
PHENOTYPE: isolates measured Resistant or Susceptible in the lab, and ask which
genetic features track the phenotype — with no AMR database in the loop.

DESIGN
  1. Cohort: N Resistant + N Susceptible E. coli to gentamicin, LABORATORY AST
     (not computational predictions — that would be circular again). From BV-BRC.
  2. Reference-free k-mer GWAS: canonical k=31 presence/absence per genome,
     association of each k-mer with R/S (Fisher exact). NO AMR database used.
  3. Top-associated k-mers -> localize in a Resistant genome -> the candidate
     resistance determinant.
  4. Validation ONLY: do the top k-mers land on real aminoglycoside genes
     (aac/aph/ant)? Annotation is used to CHECK, never to detect.

This file: stage 1+2 (cohort, download, GWAS). Stage runs incrementally; genomes
cached in data/mode_c/.

Run: python -u exploration/014_mode_c.py
"""
from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "exploration" / "data" / "mode_c"
DATA.mkdir(parents=True, exist_ok=True)

API = "https://www.bv-brc.org/api"
SPECIES = "Escherichia coli"
ANTIBIOTIC = "gentamicin"
N_PER_CLASS = 15
K = 31


def _get(url: str, tries: int = 3):
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if t == tries - 1:
                print(f"    GET failed: {e}", flush=True)
                return None
            time.sleep(1.5)


def cohort_ids():
    """Return {genome_id: 'R'/'S'} for lab-AST gentamicin E. coli, balanced."""
    out = {}
    for ph, tag in [("Resistant", "R"), ("Susceptible", "S")]:
        q = (f"and(eq(genome_name,{urllib.parse.quote(SPECIES)}*),"
             f"eq(antibiotic,{ANTIBIOTIC}),eq(evidence,Laboratory%20Method),"
             f"eq(resistant_phenotype,{ph}))")
        url = f"{API}/genome_amr/?{q}&select(genome_id)&limit(400)&http_accept=application/json"
        rows = _get(url) or []
        # unique genome_ids, take first N that we can download
        ids = []
        seen = set()
        for r in rows:
            gid = r.get("genome_id")
            if gid and gid not in seen:
                seen.add(gid)
                ids.append(gid)
        for gid in ids[: N_PER_CLASS * 3]:   # oversample; some may fail to download
            out[gid] = tag
    return out


def download_genome(gid: str) -> Path | None:
    fna = DATA / f"{gid}.fna"
    if fna.exists() and fna.stat().st_size > 1000:
        return fna
    url = (f"{API}/genome_sequence/?eq(genome_id,{gid})"
           f"&select(sequence_id,sequence)&limit(2000)&http_accept=application/json")
    rows = _get(url)
    if not rows:
        return None
    seqs = [r["sequence"] for r in rows if r.get("sequence")]
    if not seqs or sum(len(s) for s in seqs) < 1_000_000:
        return None
    with open(fna, "w") as f:
        for i, s in enumerate(seqs):
            f.write(f">{gid}_c{i}\n{s.upper()}\n")
    return fna


# ---- k-mer machinery (vectorized, both strands for strand-invariance) ----
_TAB = np.full(256, 255, np.uint8)
for i, b in enumerate(b"ACGT"):
    _TAB[b] = i
# 2-bit positional weights; for k=31 the max packed value is 2^62-1 < 2^64 (exact)
_WEIGHTS = (4 ** np.arange(K - 1, -1, -1, dtype=object)).astype(np.uint64)


def _fwd_hashes(codes_u64: np.ndarray, k: int) -> np.ndarray:
    """Packed 2-bit hashes of all length-k windows of a valid (0..3) run."""
    if len(codes_u64) < k:
        return np.empty(0, np.uint64)
    win = np.lib.stride_tricks.sliding_window_view(codes_u64, k)
    return (win * _WEIGHTS).sum(axis=1, dtype=np.uint64)


def genome_kmer_hashes(fna: Path, k: int = K) -> np.ndarray:
    """Unique k-mer hashes for one genome, both strands (per-contig, no N-span)."""
    parts = []
    for line in fna.read_text().splitlines():
        if line.startswith(">") or not line:
            continue
        codes = _TAB[np.frombuffer(line.encode(), np.uint8)]
        valid = codes != 255
        if not valid.any():
            continue
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
            if len(run) >= k:
                parts.append(_fwd_hashes(run, k))
                rc = (3 - run)[::-1].astype(np.uint64)   # reverse complement
                parts.append(_fwd_hashes(rc, k))
            idx = j
    if not parts:
        return np.empty(0, np.uint64)
    return np.unique(np.concatenate(parts))


def main():
    print("=" * 70)
    print("  014: Mode C — reference-free phenotype GWAS")
    print(f"  {SPECIES} / {ANTIBIOTIC} / lab AST / k={K}")
    print("=" * 70)

    print("\n[1] Cohort from BV-BRC...")
    ids = cohort_ids()
    print(f"    candidate ids: {sum(v=='R' for v in ids.values())} R / "
          f"{sum(v=='S' for v in ids.values())} S")

    print("\n[2] Downloading assemblies (cached)...")
    got = {"R": [], "S": []}
    for gid, tag in ids.items():
        if len(got[tag]) >= N_PER_CLASS:
            continue
        p = download_genome(gid)
        if p:
            got[tag].append((gid, p))
            print(f"    {tag} {gid}  ok ({len(got[tag])}/{N_PER_CLASS})", flush=True)
    R, S = got["R"], got["S"]
    print(f"\n    downloaded {len(R)} R + {len(S)} S genomes")
    if len(R) < 5 or len(S) < 5:
        print("    not enough genomes; aborting.")
        return

    print("\n[3] k-mer presence per genome (canonical k=31)...")
    t0 = time.time()
    R_sets, S_sets = [], []
    for tag, lst, store in [("R", R, R_sets), ("S", S, S_sets)]:
        for gid, p in lst:
            store.append(genome_kmer_hashes(p))
            print(f"    {tag} {gid}: {len(store[-1]):>8} unique kmers "
                  f"({time.time()-t0:.0f}s)", flush=True)

    print("\n[4] Association: count R-genomes and S-genomes carrying each k-mer...")
    nR, nS = len(R_sets), len(S_sets)
    # per-genome unique -> concat -> unique with counts = #genomes carrying
    allR = np.concatenate(R_sets)
    allS = np.concatenate(S_sets)
    kmR, cR = np.unique(allR, return_counts=True)
    kmS, cS = np.unique(allS, return_counts=True)
    # union key space
    allk = np.union1d(kmR, kmS)
    cntR = np.zeros(len(allk), np.int32)
    cntS = np.zeros(len(allk), np.int32)
    cntR[np.searchsorted(allk, kmR)] = cR
    cntS[np.searchsorted(allk, kmS)] = cS
    # keep informative k-mers: present in >=3 genomes total, varying
    tot = cntR + cntS
    keep = (tot >= 3) & (tot <= nR + nS - 1)
    allk, cntR, cntS = allk[keep], cntR[keep], cntS[keep]
    print(f"    informative k-mers: {len(allk):,} (of union)")

    # association score: difference in carrier fraction (R - S), Fisher-ish
    fracR = cntR / nR
    fracS = cntS / nS
    score = fracR - fracS
    np.savez(DATA / "gwas.npz", allk=allk, cntR=cntR, cntS=cntS,
             nR=nR, nS=nS,
             R_ids=np.array([g for g, _ in R]), S_ids=np.array([g for g, _ in S]))

    print("\n[5] Top R-associated k-mers (carrier fraction R vs S):")
    top = np.argsort(-score)[:15]
    for i in top:
        print(f"    hash={int(allk[i]):>20}  R={cntR[i]}/{nR}  S={cntS[i]}/{nS}  "
              f"dScore={score[i]:+.2f}")
    n_strong = int((score >= 0.6).sum())
    print(f"\n    k-mers with R-carrier fraction >= S+0.6: {n_strong}")
    print("    Saved gwas.npz. Stage 015 localizes these to genes for validation.")


if __name__ == "__main__":
    main()
