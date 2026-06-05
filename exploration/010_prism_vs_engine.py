"""010: Does PRISM/sang see AMR where SANG2's 12 signals go blind?

007-009 verdict: single-genome composition caps at ~0.6 AUC; on host-ADAPTED
genes (CTX-M, OXA) it falls below 0.5. PRISM (`sang`) is a richer, calibrated
structural engine — 14 metrics, state machines, forbidden transitions, fused
4-layer anomalies. Question: do its richer features rank AMR genes above ordinary
genes on the SAME replicons, especially the adapted cases where we failed?

For each plasmid, each CDS gets several scores; we compute per-genome AUC for
each and put them next to our engine's self-baseline AUC:

  self      = SANG2 engine composite (the incumbent)
  ncd_bg    = sang.ncd(gene, whole-plasmid)         compression distance
  cmp_bg    = sang.compare(profile(gene), profile(plasmid)).distance
  anom      = max sang.detect_anomalies score overlapping the gene

Focus first on the 6 ADAPTED replicons (ctxm/oxa) — the blind spot.

Run: python exploration/010_prism_vs_engine.py
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

import sang  # noqa: E402
from Bio import SeqIO  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "honest_auc", Path(__file__).parent / "007_honest_auc.py")
h = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h)

M = {"A": 0, "C": 1, "G": 2, "T": 3}


def encode(s: str) -> np.ndarray:
    return np.array([M.get(c, 0) for c in s], dtype=float)


def prism_gene_scores(seq: str, genes):
    """Return dict of per-gene score arrays from PRISM features."""
    bg_data = list(seq)
    bg_sig = encode(seq)
    bg_prof = sang.profile(bg_data, signal=bg_sig, depth="standard")

    # fused anomaly detector over the whole replicon (sparse points)
    anom = sang.detect_anomalies(bg_data, signal=bg_sig)
    apos = np.array([a["position"] for a in anom.anomalies], dtype=float)
    ascore = np.array([a["score"] for a in anom.anomalies], dtype=float)

    ncd_bg, cmp_bg, anom_g = [], [], []
    for (gs, ge, _l, _amr) in genes:
        gseq = seq[gs:ge]
        if len(gseq) < 150:
            ncd_bg.append(0.0); cmp_bg.append(0.0); anom_g.append(0.0); continue
        # compression distance gene vs whole plasmid
        try:
            ncd_bg.append(float(sang.ncd(gseq, seq)))
        except Exception:
            ncd_bg.append(0.0)
        # structural-profile distance gene vs plasmid background
        try:
            gp = sang.profile(list(gseq), signal=encode(gseq), depth="standard")
            cmp_bg.append(float(sang.compare(gp, bg_prof).distance))
        except Exception:
            cmp_bg.append(0.0)
        # max fused-anomaly score landing inside the gene
        if len(apos):
            inside = (apos >= gs) & (apos < ge)
            anom_g.append(float(ascore[inside].max()) if inside.any() else 0.0)
        else:
            anom_g.append(0.0)
    return {"ncd_bg": np.array(ncd_bg), "cmp_bg": np.array(cmp_bg),
            "anom": np.array(anom_g)}


def main():
    print("=" * 92)
    print("  010: PRISM/sang vs SANG2 engine — gene-level AUC (AMR vs ordinary)")
    print("=" * 92)
    paths = h.fetch_dataset()

    # adapted blind spot first, then a couple of recent-HGT controls
    order = {"ctxm": 0, "oxa": 1, "mcr": 2, "ndm": 3, "kpc": 4, "qnr": 5}
    paths = sorted(paths, key=lambda p: order.get(p.stem.split("_")[0], 9))

    print(f"\n  {'replicon':<15}{'cls':>5}{'self':>7}{'ncd_bg':>8}{'cmp_bg':>8}{'anom':>7}  best")
    rows = []
    for p in paths:
        contig, genes = h.parse_genes(p)
        n_amr = sum(g[3] for g in genes)
        if n_amr == 0 or len(genes) - n_amr < 3:
            continue
        y = np.array([g[3] for g in genes], dtype=int)
        cls = p.stem.split("_")[0]
        rec = SeqIO.read(p, "genbank")
        seq = str(rec.seq).upper()

        self_auc = h.auc_score(y, h.score_genes(contig, genes))
        t0 = time.time()
        ps = prism_gene_scores(seq, genes)
        aucs = {k: h.auc_score(y, v) for k, v in ps.items()}
        best_name = max(aucs, key=lambda k: aucs[k] if aucs[k] is not None else 0)
        rows.append((cls, self_auc, aucs, best_name))
        print(f"  {contig.id:<15}{cls:>5}{self_auc:>7.3f}"
              f"{aucs['ncd_bg']:>8.3f}{aucs['cmp_bg']:>8.3f}{aucs['anom']:>7.3f}"
              f"  {best_name} ({time.time()-t0:.1f}s)")

    print("\n" + "=" * 92)
    print("  MEANS")
    print("=" * 92)
    def col(getter):
        v = [getter(r) for r in rows if getter(r) is not None]
        return np.mean(v) if v else float("nan")
    for grp_name, grp in [("ALL", rows),
                          ("ADAPTED (ctxm/oxa)", [r for r in rows if r[0] in ("ctxm", "oxa")]),
                          ("RECENT (mcr/ndm/kpc/qnr)", [r for r in rows if r[0] not in ("ctxm", "oxa")])]:
        if not grp:
            continue
        s = np.mean([r[1] for r in grp])
        ncd = np.mean([r[2]["ncd_bg"] for r in grp])
        cmp = np.mean([r[2]["cmp_bg"] for r in grp])
        an = np.mean([r[2]["anom"] for r in grp])
        print(f"  {grp_name:<26} n={len(grp):<2}  self={s:.3f}  ncd_bg={ncd:.3f}  "
              f"cmp_bg={cmp:.3f}  anom={an:.3f}")
    print("\n  Reading it: if any PRISM column beats self on ADAPTED and stays >0.5,")
    print("  PRISM reaches the blind spot single-genome composition cannot.")


if __name__ == "__main__":
    main()
