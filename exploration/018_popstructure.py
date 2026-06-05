"""018: Is the phenotype association causal, or lineage confounding?

014-017 claimed reference-free k-mer GWAS found real resistance determinants
(AAC(3), gyrA). But a k-mer GWAS with no population-structure correction can be
fooled: if the 20 Resistant isolates are one clone and the 20 Susceptible another,
EVERY clone-specific k-mer associates with R — including the resistance gene by
coincidence. This is the main threat to the 014-017 results. Close it before
building anything on top.

Test (uses cached genomes, no downloads):
  1. MinHash bottom-k sketch per genome from its k-mers.
  2. Pairwise Jaccard -> distance matrix -> hierarchical lineages.
  3. Are R isolates spread ACROSS lineages (causal) or is R one clade / S another
     (confound)? Quantify: silhouette of phenotype vs the tree, and whether
     lineages are phenotype-mixed.
  4. Within-lineage check: in clusters containing BOTH R and S, do the resistance
     k-mers still separate R from S? If yes, association survives structure.

Run: python -u exploration/018_popstructure.py [gentamicin|ciprofloxacin]
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("mc", Path(__file__).parent / "014_mode_c.py")
mc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mc)

DRUG = sys.argv[1] if len(sys.argv) > 1 else "ciprofloxacin"
DATA = ROOT / "exploration" / ("data/mode_c_cipro" if DRUG == "ciprofloxacin"
                               else "data/mode_c")
SKETCH = 3000   # bottom-k MinHash size


def sketch_of(hashes: np.ndarray) -> np.ndarray:
    return np.sort(hashes)[:SKETCH]


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    """Bottom-k MinHash Jaccard estimate."""
    union = np.union1d(a, b)[:SKETCH]
    ina = np.isin(union, a)
    inb = np.isin(union, b)
    both = (ina & inb).sum()
    return both / len(union)


def main():
    print("=" * 72)
    print(f"  018: population structure vs phenotype — {DRUG}")
    print("=" * 72)

    # match genomes to phenotype via the cohort query (cached genomes on disk)
    mc.ANTIBIOTIC = DRUG
    mc.N_PER_CLASS = 20
    mc.DATA = DATA
    ids = mc.cohort_ids()
    genomes = []  # (gid, tag, sketch)
    for gid, tag in ids.items():
        fna = DATA / f"{gid}.fna"
        if not fna.exists():
            continue
        genomes.append((gid, tag, sketch_of(mc.genome_kmer_hashes(fna))))
        if sum(1 for g in genomes if g[1] == "R") >= 20 and \
           sum(1 for g in genomes if g[1] == "S") >= 20:
            pass
    # keep at most 20/20
    R = [g for g in genomes if g[1] == "R"][:20]
    S = [g for g in genomes if g[1] == "S"][:20]
    genomes = R + S
    labels = np.array([g[1] for g in genomes])
    n = len(genomes)
    print(f"\n  {sum(labels=='R')} R + {sum(labels=='S')} S genomes (cached)")

    print("\n[1] Pairwise Jaccard distance (MinHash)...")
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = 1 - jaccard(genomes[i][2], genomes[j][2])
            D[i, j] = D[j, i] = d

    # mean within-phenotype vs across-phenotype distance
    iu = np.triu_indices(n, 1)
    same = labels[iu[0]] == labels[iu[1]]
    print(f"    mean distance: same-phenotype={D[iu][same].mean():.3f}  "
          f"diff-phenotype={D[iu][~same].mean():.3f}")
    print("    (if same << diff, phenotype tracks lineage -> confound risk)")

    print("\n[2] Hierarchical lineages...")
    Z = linkage(squareform(D, checks=False), method="average")
    for k in (4, 6, 8):
        cl = fcluster(Z, k, criterion="maxclust")
        mixed = 0
        rows = []
        for c in range(1, k + 1):
            m = cl == c
            nr, ns = int((labels[m] == "R").sum()), int((labels[m] == "S").sum())
            if nr > 0 and ns > 0:
                mixed += 1
            rows.append(f"C{c}:{nr}R/{ns}S")
        print(f"    k={k}: {'  '.join(rows)}   ({mixed} phenotype-MIXED lineages)")

    print("\n[3] Within-lineage association (the decisive test)...")
    # use k=6 clustering; for mixed lineages, do resistance k-mers still split R/S?
    cl = fcluster(Z, 6, criterion="maxclust")
    # load resistance k-mers: recompute strong from this cohort
    R_sets = [mc.genome_kmer_hashes(DATA / f"{g[0]}.fna") for g in R]
    S_sets = [mc.genome_kmer_hashes(DATA / f"{g[0]}.fna") for g in S]
    nR, nS = len(R_sets), len(S_sets)
    kmR, cR = np.unique(np.concatenate(R_sets), return_counts=True)
    kmS, cS = np.unique(np.concatenate(S_sets), return_counts=True)
    allk = np.union1d(kmR, kmS)
    cntR = np.zeros(len(allk), np.int32); cntS = np.zeros(len(allk), np.int32)
    cntR[np.searchsorted(allk, kmR)] = cR
    cntS[np.searchsorted(allk, kmS)] = cS
    strong = np.sort(allk[(cntR / nR - cntS / nS) >= 0.7])
    print(f"    {len(strong):,} strong phenotype k-mers (whole cohort)")

    # per genome: count strong k-mers carried
    def carried(gid):
        ks = mc.genome_kmer_hashes(DATA / f"{gid}.fna")
        pos = np.clip(np.searchsorted(strong, ks), 0, len(strong) - 1)
        return int((strong[pos] == ks).sum())

    carry = {g[0]: carried(g[0]) for g in genomes}
    mixed_lineages = 0
    survived = 0
    for c in range(1, 7):
        idx = [i for i in range(n) if cl[i] == c]
        labs = labels[idx]
        if (labs == "R").sum() == 0 or (labs == "S").sum() == 0:
            continue
        mixed_lineages += 1
        r_carry = np.mean([carry[genomes[i][0]] for i in idx if labels[i] == "R"])
        s_carry = np.mean([carry[genomes[i][0]] for i in idx if labels[i] == "S"])
        ok = r_carry > s_carry * 1.5 + 1
        survived += ok
        print(f"    lineage C{c} ({(labs=='R').sum()}R/{(labs=='S').sum()}S): "
              f"strong-kmers carried R={r_carry:.0f} S={s_carry:.0f}  "
              f"{'separates' if ok else 'no split'}")

    print("\n" + "=" * 72)
    print(f"  phenotype-mixed lineages: {mixed_lineages}; "
          f"resistance signal separates R/S within {survived} of them")
    if mixed_lineages >= 2 and survived >= max(1, mixed_lineages - 1):
        print("  -> association survives population structure: R isolates from DIFFERENT")
        print("     lineages all carry the resistance k-mers, S from the same lineages")
        print("     don't. Causal, not a single-clone artifact.")
    elif mixed_lineages == 0:
        print("  -> R and S are cleanly separated lineages: CANNOT rule out confound")
        print("     from this cohort. Need phenotype-mixed lineages to test causality.")
    else:
        print("  -> mixed evidence; inspect per-lineage rows above.")


if __name__ == "__main__":
    main()
