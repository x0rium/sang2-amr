"""008: Does per-genome AUC track host-donor compositional distance?

007 showed the AUC is bimodal: strong on recent HGT (mcr/NDM), worse-than-coin
on host-adapted genes (CTX-M/OXA). HYPOTHESIS: the method works exactly when the
AMR genes are compositionally FOREIGN relative to the background the baseline is
built from. If true, a cheap compositional-distance proxy predicts the AUC, and
becomes a per-hit confidence / regime gate that fixes the AUC<0.5 cases.

We can't know the true donor, but the baseline IS the plasmid background, so the
operational "host-donor distance" is: how far do AMR genes sit from the rest of
their own replicon, in the very signals the engine uses.

Predictors per replicon (gene-level, background = whole plasmid):
  sep_gc     = mean|GC_gene - GC_bg| over AMR  minus  same over non-AMR
  sep_rscu   = same, codon-usage (RSCU) distance from plasmid profile
  sep_rho    = same, dinucleotide (rho) distance from plasmid profile
  spcGC      = |GC(AMR genes) - typical chromosomal GC of the host genus|

Then Spearman/Pearson( predictor , per-genome AUC ) across the 18 replicons.

Run: python exploration/008_auc_vs_distance.py
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from io_module import contig_to_array  # noqa: E402
from signals import gc_content, rho_vector, rho_distance  # noqa: E402
from genome_signals import rscu, codon_bias_distance  # noqa: E402

# reuse 007's fetch/parse/scan/auc without duplicating
_spec = importlib.util.spec_from_file_location(
    "honest_auc", Path(__file__).parent / "007_honest_auc.py")
h = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h)

# typical chromosomal GC% by genus (rough, for the species-distance proxy)
GENUS_GC = {
    "Escherichia": 50.8, "Klebsiella": 57.0, "Pseudomonas": 66.0,
    "Salmonella": 52.2, "Enterobacter": 55.0, "Citrobacter": 52.0,
}


def genus_of(desc: str) -> str | None:
    for g in GENUS_GC:
        if g in desc:
            return g
    return None


def gene_distances(contig, genes):
    """Per-gene compositional distance from the plasmid background."""
    seq = contig_to_array(contig)
    bg_gc = gc_content(seq)
    bg_rscu = rscu(seq)
    bg_rho = rho_vector(seq)
    rows = []  # (is_amr, d_gc, d_rscu, d_rho, gene_gc)
    for (gs, ge, _label, amr) in genes:
        sub = seq[gs:ge]
        if len(sub) < 150:
            continue
        rows.append((
            int(amr),
            abs(gc_content(sub) - bg_gc),
            codon_bias_distance(sub, bg_rscu),
            rho_distance(sub, bg_rho),
            gc_content(sub),
        ))
    return np.array(rows, dtype=float)


def main():
    print("=" * 78)
    print("  008: AUC vs host-donor compositional distance")
    print("=" * 78)
    paths = h.fetch_dataset()

    rows = []  # per replicon: dict of metrics
    print(f"\n  {'replicon':<15}{'AUC':>6}{'sep_gc':>8}{'sep_rscu':>9}"
          f"{'sep_rho':>9}{'spcGC':>7}  host")
    for p in paths:
        contig, genes = h.parse_genes(p)
        n_amr = sum(g[3] for g in genes)
        if n_amr == 0 or len(genes) - n_amr < 3:
            continue
        scores = h.score_genes(contig, genes)
        y = np.array([g[3] for g in genes], dtype=int)
        auc = h.auc_score(y, scores)
        if auc is None:
            continue

        D = gene_distances(contig, genes)
        amr_m = D[:, 0] == 1
        sep_gc = D[amr_m, 1].mean() - D[~amr_m, 1].mean()
        sep_rscu = D[amr_m, 2].mean() - D[~amr_m, 2].mean()
        sep_rho = D[amr_m, 3].mean() - D[~amr_m, 3].mean()

        rec = h.parse_genes  # noqa  (keep import warm)
        desc = contig.id
        # read description from genbank for genus
        from Bio import SeqIO
        gbrec = SeqIO.read(p, "genbank")
        genus = genus_of(gbrec.annotations.get("organism", "") + " " + gbrec.description)
        spc_gc = abs(D[amr_m, 4].mean() - GENUS_GC[genus]) if genus else np.nan

        rows.append(dict(id=contig.id, auc=auc, sep_gc=sep_gc, sep_rscu=sep_rscu,
                         sep_rho=sep_rho, spc_gc=spc_gc, genus=genus))
        print(f"  {contig.id:<15}{auc:>6.3f}{sep_gc:>8.3f}{sep_rscu:>9.3f}"
              f"{sep_rho:>9.3f}{spc_gc:>7.1f}  {genus}")

    auc = np.array([r["auc"] for r in rows])
    print("\n" + "=" * 78)
    print("  CORRELATION with per-genome AUC  (n = %d)" % len(rows))
    print("=" * 78)
    for key, name in [("sep_gc", "GC distance from plasmid bg"),
                      ("sep_rscu", "codon-usage distance from bg"),
                      ("sep_rho", "dinucleotide distance from bg"),
                      ("spc_gc", "|AMR GC - host genus GC|")]:
        x = np.array([r[key] for r in rows], dtype=float)
        m = ~np.isnan(x)
        if m.sum() < 4:
            continue
        sr, sp = stats.spearmanr(x[m], auc[m])
        pr, pp = stats.pearsonr(x[m], auc[m])
        print(f"  {name:<34} Spearman r={sr:+.3f} (p={sp:.3f})  "
              f"Pearson r={pr:+.3f} (p={pp:.3f})")

    # best predictor as a regime gate
    print("\n  Regime check (does the strongest predictor separate strong/weak AUC?):")
    best = max(["sep_gc", "sep_rscu", "sep_rho"],
               key=lambda k: abs(stats.spearmanr(
                   [r[k] for r in rows], auc)[0]))
    x = np.array([r[best] for r in rows])
    thr = np.median(x)
    hi = auc[x >= thr]
    lo = auc[x < thr]
    print(f"    predictor = {best}, split at median {thr:+.3f}")
    print(f"    AMR more foreign than bg : mean AUC {hi.mean():.3f}  (n={len(hi)})")
    print(f"    AMR not more foreign     : mean AUC {lo.mean():.3f}  (n={len(lo)})")
    print("\n  If the gap is large, this predictor IS the confidence gate that")
    print("  turns 'AUC near noise on average' into 'AUC ~0.9 on a predictable subset'.")


if __name__ == "__main__":
    main()
