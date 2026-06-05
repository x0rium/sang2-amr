"""013: A label-free regime selector to recover the oracle gap.

012: mean-fusion = 0.707 [0.663,0.748], stable, but blind averaging HURTS
single-channel-dominant classes (mcr: self 0.930 dragged to 0.659). oracle would
recover it if we knew, per genome, which channel to trust.

011's tail-contrast selector failed because it read score SHAPE (self is always a
sharper rank-fusion composite). Need a cue ORTHOGONAL to the score: a property of
the PLASMID itself.

MECHANISTIC HYPOTHESIS: self (composition) is reliable exactly when the plasmid
carries a localized compositional OUTLIER — a recent foreign island (mcr +
transposon). An adapted plasmid is compositionally homogeneous → no island → trust
PRISM. So the regime cue = compositional HETEROGENEITY of the plasmid, computed
from the sequence, independent of any score.

Step 1 (this script): does heterogeneity correlate with which channel wins
(self_auc - cmp_auc)? If yes, build an Otsu-thresholded selector (zero magic
numbers) and measure it vs mean-fusion and oracle, with bootstrap CIs.

Run: python -u exploration/013_regime_selector.py
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, "/Users/x0rium/CATALYSTO/Prism")
sys.path.insert(0, str(ROOT / "src"))

from Bio import SeqIO  # noqa: E402

def _load(name, file):
    s = importlib.util.spec_from_file_location(name, Path(__file__).parent / file)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m

cons = _load("consolidate", "012_consolidate.py")
h = cons.h
from genome_signals import seq_to_array  # noqa: E402
from signals import gc_content, kmer_entropy, otsu_threshold  # noqa: E402

M = {"A": 0, "C": 1, "G": 2, "T": 3}


def heterogeneity(seq: str, window=2000, step=500):
    """Spread of composition across windows — does the plasmid have an island?"""
    arr = seq_to_array(seq)
    gcs, h3s = [], []
    for i in range(0, max(1, len(arr) - window), step):
        seg = arr[i:i + window]
        gcs.append(gc_content(seg))
        h3s.append(kmer_entropy(seg, k=3))
    gcs, h3s = np.array(gcs), np.array(h3s)
    # spread = how far the most-deviant window sits from the body (robust)
    def tail(a):
        q50 = np.median(a)
        mad = np.median(np.abs(a - q50)) + 1e-9
        return float(np.max(np.abs(a - q50)) / mad)
    return {"het_gc_std": float(gcs.std()), "het_h3_std": float(h3s.std()),
            "het_gc_tail": tail(gcs), "het_h3_tail": tail(h3s)}


def percentile(arr):
    order = np.argsort(np.argsort(arr))
    return order / max(len(arr) - 1, 1)


def boot_ci(vals, B=3000):
    vals = np.asarray(vals, float)
    if len(vals) == 0:
        return (float("nan"),) * 3
    idx = np.random.default_rng(0).integers(0, len(vals), size=(B, len(vals)))
    m = vals[idx].mean(axis=1)
    return float(vals.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def main():
    print("=" * 74)
    print("  013: regime selector from plasmid heterogeneity")
    print("=" * 74)
    paths = sorted(h.DATA.glob("*.gb"))
    rows = []
    print(f"\n[1] Loading cached scores + computing heterogeneity for {len(paths)} plasmids...")
    for p in paths:
        r = cons.scores_for(p)
        if r is None:
            continue
        y, s, c, cls = r
        a_self = h.auc_score(y, s)
        a_cmp = h.auc_score(y, c)
        if None in (a_self, a_cmp):
            continue
        rec = SeqIO.read(p, "genbank")
        het = heterogeneity(str(rec.seq).upper())
        rows.append(dict(cls=cls, y=y, s=s, c=c, a_self=a_self, a_cmp=a_cmp, **het))
    print(f"    {len(rows)} replicons.")

    adv = np.array([r["a_self"] - r["a_cmp"] for r in rows])  # >0 => self wins
    print("\n[2] Does heterogeneity predict which channel wins (self_auc - cmp_auc)?")
    best_feat, best_r = None, 0.0
    for feat in ["het_gc_std", "het_h3_std", "het_gc_tail", "het_h3_tail"]:
        x = np.array([r[feat] for r in rows])
        sr, sp = stats.spearmanr(x, adv)
        print(f"    {feat:<14} Spearman r={sr:+.3f} (p={sp:.4f})")
        if abs(sr) > abs(best_r):
            best_feat, best_r = feat, sr

    print(f"\n[3] Selector on best feature: {best_feat} (r={best_r:+.3f})")
    x = np.array([r[feat] for r in rows]) if False else np.array([r[best_feat] for r in rows])
    thr, _ = otsu_threshold(x.tolist())
    # if heterogeneity high -> island present -> trust self; else trust cmp
    hi_self = best_r >= 0   # sign tells us which way the cue points
    sel_auc, fus_auc, oracle = [], [], []
    for i, r in enumerate(rows):
        ps, pc = percentile(r["s"]), percentile(r["c"])
        fus = h.auc_score(r["y"], (ps + pc) / 2)
        use_self = (x[i] >= thr) if hi_self else (x[i] < thr)
        sa = r["a_self"] if use_self else r["a_cmp"]
        sel_auc.append(sa); fus_auc.append(fus)
        oracle.append(max(r["a_self"], r["a_cmp"]))

    def line(name, vals):
        m, lo, hi = boot_ci(vals)
        print(f"    {name:<14} {m:.3f} [{lo:.3f}, {hi:.3f}]")
    print(f"\n  Otsu threshold on {best_feat} = {thr:.4f}  (hi->self={hi_self})")
    print("\n  Mean AUC (bootstrap CI over 76 genomes):")
    line("self", [r["a_self"] for r in rows])
    line("cmp", [r["a_cmp"] for r in rows])
    line("mean-fusion", fus_auc)
    line("SELECTOR", sel_auc)
    line("oracle", oracle)

    # how often selector matches oracle channel
    correct = sum(1 for i, r in enumerate(rows)
                  if ((x[i] >= thr) == (r["a_self"] >= r["a_cmp"])) == hi_self
                  or ((x[i] >= thr) if hi_self else (x[i] < thr)) == (r["a_self"] >= r["a_cmp"]))
    # simpler honest count:
    pick_self = (x >= thr) if hi_self else (x < thr)
    oracle_self = np.array([r["a_self"] >= r["a_cmp"] for r in rows])
    acc = float((pick_self == oracle_self).mean())
    print(f"\n  Selector picks oracle-best channel on {acc:.0%} of genomes.")
    print("  Win: SELECTOR CI lower bound > mean-fusion point estimate, and mcr restored.")
    # --- soft fusion: weight channels by heterogeneity instead of hard switch ---
    # het high -> cmp better (r<0) -> more weight on cmp. Rank-normalize het to [0,1].
    hv = np.array([r[best_feat] for r in rows])
    w_cmp = (np.argsort(np.argsort(hv)) / max(len(hv) - 1, 1))  # 0..1, high het -> ~1
    if hi_self:
        w_cmp = 1 - w_cmp
    soft_auc = []
    for i, r in enumerate(rows):
        ps, pc = percentile(r["s"]), percentile(r["c"])
        w = w_cmp[i]
        soft_auc.append(h.auc_score(r["y"], (1 - w) * ps + w * pc))
    line("soft-fusion", soft_auc)

    # mcr sanity
    mcr = [r for r in rows if r["cls"] == "mcr"]
    if mcr:
        ms = np.array([r["het_gc_std"] for r in mcr]).mean()
        ad = np.array([r["het_gc_std"] for r in rows if r["cls"] in cons.ADAPTED]).mean()
        print(f"\n  het_gc_std: mcr={ms:.4f}  adapted={ad:.4f}  "
              f"(hypothesis: mcr higher = has island)")


if __name__ == "__main__":
    main()
