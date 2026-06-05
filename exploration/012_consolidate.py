"""012: Consolidate the two-channel result on a larger replicon set + bootstrap CI.

011 (n=18) found mean-fusion(self, PRISM cmp) = 0.686 AUC, beating the incumbent
(0.625) and not collapsing on adapted genes. Is that real or an n=18 artifact?

This run: ~60-80 RefSeq plasmids across more AMR classes, per-genome AUC for
self / cmp / mean-fusion, and BOOTSTRAP CIs over genomes (the honest stability
check on the means — no expensive per-gene surrogate needed for that question).

Per-plasmid scores cached in data/gene_scores/<id>.npz so the heavy PRISM step
runs once and re-runs are incremental.

Run: python -u exploration/012_consolidate.py
"""
from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, "/Users/x0rium/CATALYSTO/Prism")
sys.path.insert(0, str(ROOT / "src"))

from Bio import SeqIO  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "honest_auc", Path(__file__).parent / "007_honest_auc.py")
h = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h)
_spec2 = importlib.util.spec_from_file_location(
    "prism_vs", Path(__file__).parent / "010_prism_vs_engine.py")
pv = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(pv)

CACHE_DIR = ROOT / "exploration" / "data" / "gene_scores"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# broader net: more classes, more per query. adapted = ctxm/oxa/tet/sul/dfr,
# recent-leaning = mcr/ndm/kpc/qnr/aac/van (empirical regime checked later too).
h.QUERIES = {
    "kpc":  'blaKPC[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "ndm":  'blaNDM[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "ctxm": 'blaCTX-M[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "oxa":  'blaOXA-48[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "mcr":  'mcr-1[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "qnr":  'qnrS[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "tet":  'tet(X)[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "sul":  'sul1[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "aac":  'aac(6)[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "van":  'vanA[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "erm":  'erm(B)[All Fields] AND plasmid[Title] AND complete sequence[Title]',
}
ADAPTED = {"ctxm", "oxa", "tet", "sul"}


def percentile(arr):
    order = np.argsort(np.argsort(arr))
    return order / max(len(arr) - 1, 1)


def scores_for(p):
    """Cached (y, self, cmp, cls) for one plasmid."""
    cid = p.stem.split("_", 1)[1] if "_" in p.stem else p.stem
    cls = p.stem.split("_")[0]
    cache = CACHE_DIR / f"{cls}__{cid}.npz"
    if cache.exists():
        d = np.load(cache)
        return d["y"], d["self"], d["cmp"], cls
    contig, genes = h.parse_genes(p)
    n_amr = sum(g[3] for g in genes)
    if n_amr == 0 or len(genes) - n_amr < 3:
        return None
    rec = SeqIO.read(p, "genbank")
    seq = str(rec.seq).upper()
    y = np.array([g[3] for g in genes], dtype=int)
    s = h.score_genes(contig, genes)
    c = pv.prism_gene_scores(seq, genes)["cmp_bg"]
    np.savez(cache, y=y, **{"self": s, "cmp": c})
    return y, s, c, cls


def boot_ci(vals, B=3000):
    vals = np.asarray(vals, dtype=float)
    n = len(vals)
    if n == 0:
        return (float("nan"),) * 3
    idx = np.random.default_rng(0).integers(0, n, size=(B, n))
    means = vals[idx].mean(axis=1)
    return float(vals.mean()), float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def main():
    print("=" * 72)
    print("  012: consolidate two-channel detector + bootstrap CI")
    print("=" * 72)
    print("\n[1] Fetching extended plasmid set...")
    paths = h.fetch_dataset(per_query=10)
    # include everything already on disk too
    disk = sorted(h.DATA.glob("*.gb"))
    seen = {p.name for p in paths}
    paths += [p for p in disk if p.name not in seen]

    print(f"\n[2] Scoring {len(paths)} plasmids (cached per-plasmid)...")
    rows = []
    t0 = time.time()
    for i, p in enumerate(paths):
        try:
            r = scores_for(p)
        except Exception as e:
            print(f"  skip {p.name}: {e}", flush=True)
            continue
        if r is None:
            continue
        y, s, c, cls = r
        ps, pc = percentile(s), percentile(c)
        a_self = h.auc_score(y, s)
        a_cmp = h.auc_score(y, c)
        a_mean = h.auc_score(y, (ps + pc) / 2)
        if None in (a_self, a_cmp, a_mean):
            continue
        rows.append((cls, a_self, a_cmp, a_mean))
        if (i + 1) % 10 == 0:
            print(f"  ...{i+1}/{len(paths)} ({time.time()-t0:.0f}s)", flush=True)

    print(f"\n[3] {len(rows)} replicons scored.  Bootstrap CI over genomes:\n")
    groups = [("ALL", rows),
              ("ADAPTED", [r for r in rows if r[0] in ADAPTED]),
              ("RECENT", [r for r in rows if r[0] not in ADAPTED])]
    print(f"  {'group':<10}{'n':>4}   {'self (95% CI)':>22}{'cmp (95% CI)':>22}"
          f"{'mean-fusion (95% CI)':>24}")
    for name, grp in groups:
        if not grp:
            continue
        def fmt(idx):
            m, lo, hi = boot_ci([r[idx] for r in grp])
            return f"{m:.3f} [{lo:.3f},{hi:.3f}]"
        print(f"  {name:<10}{len(grp):>4}   {fmt(1):>22}{fmt(2):>22}{fmt(3):>24}")

    # per-class breakdown
    print("\n  Per-class mean AUC (self / cmp / mean-fusion):")
    for cls in sorted({r[0] for r in rows}):
        g = [r for r in rows if r[0] == cls]
        reg = "adapt" if cls in ADAPTED else "recent"
        print(f"    {cls:<6} ({reg}) n={len(g):<2}  "
              f"self={np.mean([r[1] for r in g]):.3f}  "
              f"cmp={np.mean([r[2] for r in g]):.3f}  "
              f"mean={np.mean([r[3] for r in g]):.3f}")

    print("\n  Stable if: mean-fusion CI lower bound > self point estimate, and")
    print("  ADAPTED mean-fusion CI stays clear of 0.5. Then 0.686 wasn't an n=18 fluke.")


if __name__ == "__main__":
    main()
