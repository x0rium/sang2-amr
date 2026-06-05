"""019: Fuse phenotype association with structural anomaly — the SANG2 fingerprint.

Mode C (k-mer GWAS) gives PHENOTYPE evidence: which regions track R/S.
two_channel/prism gives STRUCTURAL evidence: which regions look foreign.
Neither alone is the project's contribution — pyseer does phenotype, IslandViewer
does structure. The contribution is reading them TOGETHER, per region:

  acquired gene (AAC(3))   : phenotype HIGH + structural HIGH  -> foreign + causal
  core-gene SNP (gyrA)     : phenotype HIGH + structural LOW   -> native, mutational
  housekeeping             : phenotype LOW  + structural LOW

If structural separates acquired determinants from core-gene-SNP determinants,
the two-channel score is doing something phenotype-only GWAS cannot: classifying
the MECHANISM TYPE of a resistance hit. That is the unique fingerprint.

This script, on a gentamicin Resistant genome, scores every CDS by:
  phenotype = fraction of its k-mers that are strong gentamicin R-associated
  structural = PRISM profile distance of the gene vs the whole genome
and shows the 2D separation: AAC(3) high-high, housekeeping low-low, and (for
contrast) where the gyrA core gene lands.

Run: python -u exploration/019_fusion.py
"""
from __future__ import annotations

import importlib.util
import json
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
import sys
sys.path.insert(0, "/Users/x0rium/CATALYSTO/Prism")
sys.path.insert(0, str(ROOT / "src"))
import prism_channel  # noqa: E402
from genome_signals import seq_to_array  # noqa: E402

API = "https://www.bv-brc.org/api"
GENOME = "562.42768"   # gentamicin Resistant, carries AAC(3)
K = 31

_spec = importlib.util.spec_from_file_location("mc", Path(__file__).parent / "014_mode_c.py")
mc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mc)
DATA = ROOT / "exploration" / "data" / "mode_c"

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


def kmer_set(seq):
    codes = _TAB[np.frombuffer(seq.upper().encode(), np.uint8)]
    out = []
    valid = codes != 255
    n, idx = len(codes), 0
    while idx < n:
        if not valid[idx]:
            idx += 1; continue
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
    print("  019: phenotype x structural fusion — mechanism-type fingerprint")
    print(f"  genome {GENOME} (gentamicin Resistant)  PRISM={prism_channel.available()}")
    print("=" * 74)

    g = np.load(DATA / "gwas.npz", allow_pickle=True)
    allk, cntR, cntS = g["allk"], g["cntR"], g["cntS"]
    nR, nS = int(g["nR"]), int(g["nS"])
    strong = np.sort(allk[(cntR / nR - cntS / nS) >= 0.8])
    print(f"\n  gentamicin strong phenotype k-mers: {len(strong):,}")

    # Pull a panel of CDS: AMR-annotated + a sample of housekeeping, with coords.
    print("\n  Fetching CDS panel (AMR + housekeeping controls)...")
    panel = []
    for kw, kind in [("aminoglycoside", "AMR"), ("beta-lactamase", "AMR"),
                     ("sulfonamide", "AMR"), ("gyrase", "core"),
                     ("ribosomal+protein", "house"), ("dehydrogenase", "house"),
                     ("synthase", "house")]:
        url = (f"{API}/genome_feature/?and(eq(genome_id,{GENOME}),eq(feature_type,CDS),"
               f"keyword({kw}))&select(product,accession,start,end)&limit(6)"
               f"&http_accept=application/json")
        for f in (_get(url) or []):
            panel.append((kind, f.get("product", "?"), f.get("accession"),
                          int(f.get("start", 0)), int(f.get("end", 0))))

    # whole-genome PRISM profile (structural reference) from the cached assembly
    print("  Building whole-genome PRISM profile...")
    fna = DATA / f"{GENOME}.fna"
    whole = "".join(l for l in fna.read_text().splitlines() if not l.startswith(">"))
    whole_arr = seq_to_array(whole)
    whole_prof = prism_channel.profile_of(whole_arr)

    contig_cache = {}
    def contig(acc):
        if acc not in contig_cache:
            u = (f"{API}/genome_sequence/?and(eq(genome_id,{GENOME}),eq(sequence_id,{acc}))"
                 f"&select(sequence)&limit(1)&http_accept=application/json")
            rows = _get(u)
            contig_cache[acc] = rows[0]["sequence"].upper() if rows else ""
        return contig_cache[acc]

    print(f"\n  {'kind':<6}{'phenotype':>10}{'structural':>11}  product")
    rows = []
    seen = set()
    for kind, prod, acc, s, e in panel:
        if prod in seen:
            continue
        seen.add(prod)
        cseq = contig(acc)
        sub = cseq[max(0, min(s, e) - 1):max(s, e)]
        if len(sub) < K:
            continue
        ks = kmer_set(sub)
        pos = np.clip(np.searchsorted(strong, ks), 0, len(strong) - 1)
        pheno = (strong[pos] == ks).sum() / max(len(ks), 1)
        struct = prism_channel.region_distance(seq_to_array(sub), whole_prof)
        rows.append((kind, pheno, struct, prod))

    rows.sort(key=lambda r: -(r[1] + r[2]))
    for kind, pheno, struct, prod in rows:
        tag = ""
        if pheno > 0.1 and struct > 0.3:
            tag = "  <== DUAL (acquired+causal)"
        elif pheno > 0.1:
            tag = "  <- phenotype-only (SNP-like / native)"
        print(f"  {kind:<6}{pheno:>10.2f}{struct:>11.2f}  {prod[:42]}{tag}")

    print("\n" + "=" * 74)
    dual = [r for r in rows if r[1] > 0.1 and r[2] > 0.3]
    phen_only = [r for r in rows if r[1] > 0.1 and r[2] <= 0.3]
    print(f"  DUAL-evidence hits (phenotype+structural): {len(dual)}")
    for _, p, s, prod in dual:
        print(f"    pheno={p:.2f} struct={s:.2f}  {prod[:50]}")
    print(f"  phenotype-only hits: {len(phen_only)}")
    print("\n  Reading it: the structural channel splits phenotype hits into")
    print("  acquired/foreign (dual evidence) vs native (phenotype-only). pyseer")
    print("  sees only the phenotype axis; this fingerprint adds the mechanism axis.")


if __name__ == "__main__":
    main()
