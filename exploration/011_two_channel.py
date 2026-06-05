"""011: Two-channel detector — composition (self) + PRISM profile, rank-fused.

010 showed the two channels are anti-correlated by regime:
  self (composition)    wins on recent HGT      (mean AUC 0.714, blind on adapted 0.449)
  cmp_bg (PRISM profile) wins on adapted genes   (mean AUC 0.684, robust, never <0.49)

Averaging hurts (anti-correlated). The honest combiner is the engine's OWN
philosophy — rank fusion, "a gene needs only ONE strong signal". We add cmp_bg as
an extra independent voter. Per gene: take its within-genome percentile in EACH
channel, combined = max(percentile_self, percentile_cmp). No threshold, no labels,
no gate — best channel wins per gene.

Compare per-genome AUC: self | cmp | combined(max) | combined(mean) | oracle.
oracle = max(self_auc, cmp_auc) per genome (upper bound, picks in hindsight).

Caches per-gene scores to data/two_channel_cache.npz so re-runs are instant.

Run: python -u exploration/011_two_channel.py
"""
from __future__ import annotations

import importlib.util
import sys
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

CACHE = ROOT / "exploration" / "data" / "two_channel_cache.npz"


def percentile(arr: np.ndarray) -> np.ndarray:
    """Within-genome percentile rank in [0,1], higher value -> higher rank."""
    order = np.argsort(np.argsort(arr))
    n = len(arr)
    return order / max(n - 1, 1)


def tail_contrast(arr: np.ndarray) -> float:
    """How sharply the top value stands out from the body (in IQR units).

    Label-free regime cue: a real AMR signal is one/two genes far above the rest.
    The channel with the sharper top tail is the trustworthy one on this genome.
    """
    q25, q50, q75 = np.percentile(arr, [25, 50, 75])
    iqr = q75 - q25
    return float((arr.max() - q50) / (iqr + 1e-9))


def build_cache():
    # read already-downloaded plasmids from disk (robust to NCBI esearch hiccups)
    paths = sorted(h.DATA.glob("*.gb"))
    if not paths:
        paths = h.fetch_dataset()
    order = {"ctxm": 0, "oxa": 1, "mcr": 2, "ndm": 3, "kpc": 4, "qnr": 5}
    paths = sorted(paths, key=lambda p: order.get(p.stem.split("_")[0], 9))
    store = {}
    for p in paths:
        contig, genes = h.parse_genes(p)
        n_amr = sum(g[3] for g in genes)
        if n_amr == 0 or len(genes) - n_amr < 3:
            continue
        cls = p.stem.split("_")[0]
        rec = SeqIO.read(p, "genbank")
        seq = str(rec.seq).upper()
        y = np.array([g[3] for g in genes], dtype=int)
        self_scores = h.score_genes(contig, genes)
        cmp_scores = pv.prism_gene_scores(seq, genes)["cmp_bg"]
        key = contig.id
        store[f"{key}|y"] = y
        store[f"{key}|self"] = self_scores
        store[f"{key}|cmp"] = cmp_scores
        store[f"{key}|cls"] = np.array([cls])
        print(f"  cached {key:<15} {cls:<5} genes={len(genes)} amr={n_amr}", flush=True)
    np.savez(CACHE, **store)
    return store


def load_cache():
    if not CACHE.exists():
        print("[build] computing scores (cached after first run)...", flush=True)
        return build_cache()
    print(f"[load] {CACHE.name}", flush=True)
    d = np.load(CACHE, allow_pickle=True)
    return {k: d[k] for k in d.files}


def main():
    store = load_cache()
    ids = sorted({k.split("|")[0] for k in store})

    print(f"\n  {'replicon':<15}{'cls':>5}{'self':>7}{'cmp':>7}{'mean':>7}"
          f"{'select':>8}{'pick':>6}{'oracle':>8}")
    rows = []
    for gid in ids:
        y = store[f"{gid}|y"]
        s = store[f"{gid}|self"]
        c = store[f"{gid}|cmp"]
        cls = str(store[f"{gid}|cls"][0])
        ps, pc = percentile(s), percentile(c)
        comb_mean = (ps + pc) / 2
        # label-free regime selector: trust the channel with the sharper top tail
        use_self = tail_contrast(s) >= tail_contrast(c)
        sel_scores = s if use_self else c
        a_self = h.auc_score(y, s)
        a_cmp = h.auc_score(y, c)
        a_mean = h.auc_score(y, comb_mean)
        a_sel = h.auc_score(y, sel_scores)
        oracle = max(a_self, a_cmp)
        pick = "self" if use_self else "cmp"
        rows.append((cls, a_self, a_cmp, a_mean, a_sel, oracle, pick))
        print(f"  {gid:<15}{cls:>5}{a_self:>7.3f}{a_cmp:>7.3f}{a_mean:>7.3f}"
              f"{a_sel:>8.3f}{pick:>6}{oracle:>8.3f}", flush=True)

    print("\n" + "=" * 64)
    print("  MEANS")
    print("=" * 64)
    def grpmean(grp, idx):
        return np.mean([r[idx] for r in grp])
    for name, grp in [("ALL", rows),
                      ("ADAPTED (ctxm/oxa)", [r for r in rows if r[0] in ("ctxm", "oxa")]),
                      ("RECENT", [r for r in rows if r[0] not in ("ctxm", "oxa")])]:
        print(f"  {name:<22} n={len(grp):<2} self={grpmean(grp,1):.3f}  "
              f"cmp={grpmean(grp,2):.3f}  mean={grpmean(grp,3):.3f}  "
              f"SELECT={grpmean(grp,4):.3f}  oracle={grpmean(grp,5):.3f}")
    # how often the label-free selector picks the right channel
    correct = sum(1 for r in rows if (r[4] == r[5]))
    print(f"\n  Selector picked the oracle-best channel on {correct}/{len(rows)} genomes.")
    print("  Win: SELECT close to oracle and >0.5 on ADAPTED -> one detector,")
    print("  both regimes, no labels, no magic threshold.")


if __name__ == "__main__":
    main()
