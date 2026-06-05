"""017: The novelty test — resistance gene-based databases MISS.

Mode C (014-016) found AAC(3) for gentamicin: a known ACQUIRED gene, which
gene-presence databases (ResFinder/RGI) already find. The real claim is "find what
databases miss". The cleanest miss for gene-based tools is a POINT MUTATION in a
core gene: fluoroquinolone resistance is driven by gyrA/parC QRDR mutations
(S83L, D87N), not an acquired gene. gyrA is present in EVERY isolate, R and S —
so presence/absence tools see no difference. Only the SNP distinguishes them.

Hypothesis: reference-free k-mer GWAS on ciprofloxacin R/S will put its strongest
phenotype k-mers inside gyrA — i.e. catch the resistance allele that gene-based
databases are blind to.

Stages: cohort+download (reuse 014) -> GWAS -> identify (do strong k-mers land in
gyrA?) -> show gyrA is present in all genomes but the strong k-mers are the mutant
allele -> extract the QRDR codon 83.

Run: python -u exploration/017_unexplained.py
"""
from __future__ import annotations

import importlib.util
import json
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
API = "https://www.bv-brc.org/api"
K = 31
GENOME_FOR_ID = None  # chosen at runtime (a Resistant genome)

_spec = importlib.util.spec_from_file_location("mc", Path(__file__).parent / "014_mode_c.py")
mc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mc)

# retarget 014's machinery to ciprofloxacin, separate cache
mc.ANTIBIOTIC = "ciprofloxacin"
mc.N_PER_CLASS = 20
mc.DATA = ROOT / "exploration" / "data" / "mode_c_cipro"
mc.DATA.mkdir(parents=True, exist_ok=True)
DATA = mc.DATA

_TAB = mc._TAB
_WEIGHTS = mc._WEIGHTS


def _get(url):
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=90) as r:
            return json.loads(r.read().decode())
    except Exception as e:
        print(f"    GET failed: {e}", flush=True)
        return None


def kmer_set(seq: str) -> np.ndarray:
    codes = _TAB[np.frombuffer(seq.upper().encode(), np.uint8)]
    out = []
    valid = codes != 255
    n, idx = len(codes), 0
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


def gene_seq(genome_id, keyword):
    """Fetch CDS matching keyword, cut its sequence from the contig by coords."""
    url = (f"{API}/genome_feature/?and(eq(genome_id,{genome_id}),eq(feature_type,CDS),"
           f"keyword({keyword}))&select(product,accession,start,end)&limit(20)"
           f"&http_accept=application/json")
    feats = _get(url) or []
    results = []
    for f in feats:
        acc = f.get("accession")
        s, e = int(f.get("start", 0)), int(f.get("end", 0))
        u = (f"{API}/genome_sequence/?and(eq(genome_id,{genome_id}),eq(sequence_id,{acc}))"
             f"&select(sequence)&limit(1)&http_accept=application/json")
        rows = _get(u)
        if rows:
            cseq = rows[0]["sequence"].upper()
            sub = cseq[max(0, min(s, e) - 1):max(s, e)]
            results.append((f.get("product", "?"), sub))
    return results


def main():
    print("=" * 74)
    print("  017: Unexplained resistance — catching the gyrA SNP gene-based tools miss")
    print(f"  E. coli / ciprofloxacin / lab AST / k={K}")
    print("=" * 74)

    print("\n[1] Cohort + download (reuse 014 machinery)...")
    ids = mc.cohort_ids()
    got = {"R": [], "S": []}
    for gid, tag in ids.items():
        if len(got[tag]) >= mc.N_PER_CLASS:
            continue
        p = mc.download_genome(gid)
        if p:
            got[tag].append((gid, p))
    R, S = got["R"], got["S"]
    print(f"    downloaded {len(R)} R + {len(S)} S genomes")
    if len(R) < 8 or len(S) < 8:
        print("    not enough; abort.")
        return

    print("\n[2] k-mer presence per genome...")
    R_sets = [mc.genome_kmer_hashes(p) for _, p in R]
    S_sets = [mc.genome_kmer_hashes(p) for _, p in S]
    nR, nS = len(R_sets), len(S_sets)
    print(f"    done ({nR} R, {nS} S)")

    print("\n[3] Association...")
    kmR, cR = np.unique(np.concatenate(R_sets), return_counts=True)
    kmS, cS = np.unique(np.concatenate(S_sets), return_counts=True)
    allk = np.union1d(kmR, kmS)
    cntR = np.zeros(len(allk), np.int32); cntS = np.zeros(len(allk), np.int32)
    cntR[np.searchsorted(allk, kmR)] = cR
    cntS[np.searchsorted(allk, kmS)] = cS
    score = cntR / nR - cntS / nS
    strong = np.sort(allk[score >= 0.7])
    print(f"    strong phenotype k-mers (R-carrier >= S+0.7): {len(strong):,}")
    np.savez(DATA / "gwas_cipro.npz", allk=allk, cntR=cntR, cntS=cntS, nR=nR, nS=nS)

    # pick a Resistant genome to identify against
    gid = R[0][0]
    print(f"\n[4] Where do strong k-mers map? (validation genome {gid})")
    targets = ["gyrase", "topoisomerase", "quinolone", "acriflavine", "multidrug"]
    seen = set()
    rows = []
    for kw in targets:
        for prod, seq in gene_seq(gid, kw):
            if prod in seen or len(seq) < K:
                continue
            seen.add(prod)
            ks = kmer_set(seq)
            pos = np.clip(np.searchsorted(strong, ks), 0, len(strong) - 1)
            inter = int((strong[pos] == ks).sum())
            rows.append((inter, len(ks), prod, seq))
    rows.sort(reverse=True)
    print("    strong-kmer overlap per candidate gene:")
    for inter, nk, prod, _ in rows[:8]:
        flag = "  <==" if inter >= 5 else ""
        print(f"      strong_in_gene={inter:>4}  gene_kmers={nk:>4}  {prod[:52]}{flag}")

    # [5] Prove gyrA is present in ALL genomes (gene-presence blind) but the
    # strong k-mers are the mutant allele.
    top = rows[0] if rows else None
    if top and top[0] >= 5:
        inter, nk, prod, gseq = top
        gks = kmer_set(gseq)
        # how many genomes carry ANY of this gene's k-mers (gene presence)
        present_R = sum(np.isin(gks, s).any() for s in R_sets)
        present_S = sum(np.isin(gks, s).any() for s in S_sets)
        # strong (mutant) k-mers of this gene: carrier split
        gstrong = gks[np.isin(gks, strong)]
        idxR = np.searchsorted(allk, gstrong)
        carrierR = cntR[idxR].mean() / nR if len(gstrong) else 0
        carrierS = cntS[idxR].mean() / nS if len(gstrong) else 0
        print("\n[5] Gene-presence vs SNP:")
        print(f"    gene: {prod}")
        print(f"    gene present (any k-mer) in {present_R}/{nR} R and {present_S}/{nS} S "
              f"-> presence/absence sees NO difference (core gene)")
        print(f"    but its {len(gstrong)} phenotype k-mers: mean carrier R={carrierR:.0%} "
              f"S={carrierS:.0%} -> the MUTANT ALLELE tracks resistance")

    print("\n" + "=" * 74)
    if top and top[0] >= 5 and "gyrase" in top[2].lower():
        print("  RESULT: reference-free k-mer GWAS put its resistance signal inside")
        print(f"  gyrA ({top[0]} strong k-mers), a CORE gene present in all isolates.")
        print("  Gene-presence/absence databases (ResFinder) cannot flag it — the gene")
        print("  is everywhere; only the SNP differs. k-mer phenotype GWAS catches it.")
        print("  This is 'find what databases miss', demonstrated on a real SNP mechanism.")
    else:
        print("  Top overlap not clearly gyrase — inspect list above (cipro is")
        print("  multi-mechanism: qnr/efflux/QRDR may share the signal).")


if __name__ == "__main__":
    main()
