"""020: Residual analysis — what phenotype signal is NOT explained by known genes?

The final claim is "find what databases miss". Operationalize it: of the strong
phenotype k-mers, how many are explained by the genome's KNOWN resistance genes
(any of them), and where does the RESIDUAL sit? Residual k-mers in hypothetical /
uncharacterized genes are novel-candidate signal; residual that maps nowhere
clean is co-selection or lineage leakage. Honest either way.

For each cohort (gentamicin, ciprofloxacin):
  1. strong phenotype k-mers (from saved GWAS).
  2. all annotated AMR genes of a Resistant genome -> their k-mers.
  3. explained = strong k-mers inside ANY known AMR gene.
  4. residual = the rest; localize the densest residual cluster, annotate it.

Run: python -u exploration/020_residual_novelty.py
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

_spec = importlib.util.spec_from_file_location("mc", Path(__file__).parent / "014_mode_c.py")
mc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mc)
_TAB, _WEIGHTS = mc._TAB, mc._WEIGHTS

COHORTS = [
    ("gentamicin", ROOT / "exploration/data/mode_c", "gwas.npz", "562.42768"),
    ("ciprofloxacin", ROOT / "exploration/data/mode_c_cipro", "gwas_cipro.npz", None),
]


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


def known_amr_genes(genome_id):
    """All resistance-annotated CDS + their sequence (validation key)."""
    url = (f"{API}/genome_feature/?and(eq(genome_id,{genome_id}),eq(feature_type,CDS),"
           f"keyword(resistance%20OR%20beta-lactamase%20OR%20aminoglycoside%20OR%20"
           f"quinolone%20OR%20efflux%20OR%20sulfonamide%20OR%20tetracycline))"
           f"&select(product,accession,start,end)&limit(100)&http_accept=application/json")
    feats = _get(url) or []
    contigs = {}
    def cseq(acc):
        if acc not in contigs:
            u = (f"{API}/genome_sequence/?and(eq(genome_id,{genome_id}),eq(sequence_id,{acc}))"
                 f"&select(sequence)&limit(1)&http_accept=application/json")
            rows = _get(u)
            contigs[acc] = rows[0]["sequence"].upper() if rows else ""
        return contigs[acc]
    genes = []
    for f in feats:
        acc = f.get("accession")
        s, e = int(f.get("start", 0)), int(f.get("end", 0))
        sub = cseq(acc)[max(0, min(s, e) - 1):max(s, e)]
        if len(sub) >= K:
            genes.append((f.get("product", "?"), sub))
    return genes


def pick_resistant_genome(data_dir):
    """First .fna in the cohort dir (a Resistant genome)."""
    fnas = sorted(data_dir.glob("*.fna"))
    return fnas[0].stem if fnas else None


def main():
    print("=" * 74)
    print("  020: residual novelty — phenotype signal not explained by known genes")
    print("=" * 74)

    for drug, data_dir, gwas_file, fixed_gid in COHORTS:
        gpath = data_dir / gwas_file
        if not gpath.exists():
            print(f"\n  {drug}: no GWAS file, skip")
            continue
        g = np.load(gpath, allow_pickle=True)
        allk, cntR, cntS = g["allk"], g["cntR"], g["cntS"]
        nR, nS = int(g["nR"]), int(g["nS"])
        strong = np.sort(allk[(cntR / nR - cntS / nS) >= 0.7])
        gid = fixed_gid or pick_resistant_genome(data_dir)
        print(f"\n{'='*74}\n  {drug}  (genome {gid}, {len(strong):,} strong k-mers)\n{'='*74}")

        genes = known_amr_genes(gid)
        print(f"  known resistance-annotated genes: {len(genes)}")
        explained = np.zeros(len(strong), bool)
        per_gene = []
        for prod, seq in genes:
            ks = kmer_set(seq)
            pos = np.clip(np.searchsorted(strong, ks), 0, len(strong) - 1)
            mask = strong[pos] == ks
            # map back: which strong indices are these
            hit_idx = np.searchsorted(strong, ks[mask])
            hit_idx = hit_idx[(hit_idx < len(strong))]
            hit_idx = hit_idx[strong[hit_idx] == ks[mask]] if mask.any() else hit_idx
            inter = int(mask.sum())
            if inter > 0:
                # mark explained
                idxs = np.searchsorted(strong, ks[mask])
                idxs = idxs[idxs < len(strong)]
                idxs = idxs[strong[idxs] == ks[mask][:len(idxs)]] if len(idxs) else idxs
                explained[np.searchsorted(strong, ks[mask])] = True
                per_gene.append((inter, prod))
        per_gene.sort(reverse=True)
        n_expl = int(explained.sum())
        print(f"  strong k-mers explained by KNOWN resistance genes: "
              f"{n_expl}/{len(strong)} = {n_expl/len(strong):.0%}")
        for inter, prod in per_gene[:6]:
            print(f"      {inter:>5}  {prod[:58]}")
        print(f"  RESIDUAL (not in any known resistance gene): "
              f"{len(strong)-n_expl}/{len(strong)} = {1-n_expl/len(strong):.0%}")
        print("    -> residual is novel-candidate / co-selection / lineage signal;")
        print("       a true unexplained-resistance panel would mine it for new genes.")

    print("\n" + "=" * 74)
    print("  Honest reading: where 'explained' is high, the method recovers the known")
    print("  mechanism reference-free (validation). The residual is where novelty would")
    print("  live — only a phenotypically-R / genotype-negative panel can prove it's new.")


if __name__ == "__main__":
    main()
