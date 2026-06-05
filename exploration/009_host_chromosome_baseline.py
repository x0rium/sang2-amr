"""009: Change the analysis method, not the weights.

007/008 finding: self-baseline (background = the plasmid itself) goes blind on
host-adapted AMR genes (CTX-M, OXA) because the AMR gene doesn't stand out from
the rest of the mobile plasmid. AUC drops below 0.5.

NEW METHOD: build the baseline from the HOST CHROMOSOME (core genome of the
species), not from the plasmid. Foreignness is then measured relative to "what
is normal for this species". An AMR gene that adapted to the plasmid backbone may
still be foreign to the chromosome — so it should resurface.

Test: same 18 plasmids, same AUC benchmark, two baselines side by side.
Focus on the adapted cases (CTX-M, OXA) where self-baseline failed.

Run: python exploration/009_host_chromosome_baseline.py
"""
from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
from Bio import Entrez, SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from io_module import contig_to_array  # noqa: E402
from genome_signals import seq_to_array  # noqa: E402
from anomaly import compute_host_baseline, scan_sequence  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "honest_auc", Path(__file__).parent / "007_honest_auc.py")
h = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h)

Entrez.email = "tuxera@gmail.com"
CHROM_DIR = h.DATA.parent / "chromosomes"
CHROM_DIR.mkdir(parents=True, exist_ok=True)

# one RefSeq reference chromosome per genus (fetched via esearch, not hardcoded)
GENERA = ["Escherichia coli", "Klebsiella pneumoniae", "Pseudomonas aeruginosa",
          "Salmonella enterica", "Enterobacter hormaechei", "Citrobacter"]


def genus_key(desc: str) -> str | None:
    table = {"Escherichia": "Escherichia coli", "Klebsiella": "Klebsiella pneumoniae",
             "Pseudomonas": "Pseudomonas aeruginosa", "Salmonella": "Salmonella enterica",
             "Enterobacter": "Enterobacter hormaechei", "Citrobacter": "Citrobacter"}
    for g, full in table.items():
        if g in desc:
            return full
    return None


def fetch_chromosome(species: str) -> np.ndarray | None:
    """esearch a RefSeq reference chromosome for the species, cache FASTA, return array."""
    fa = CHROM_DIR / (species.replace(" ", "_") + ".fasta")
    if not fa.exists():
        term = (f'"{species}"[Organism] AND complete genome[Title] '
                f'AND chromosome[Title] AND srcdb_refseq[PROP] AND 3000000:7000000[SLEN]')
        try:
            r = Entrez.read(Entrez.esearch(db="nuccore", term=term, retmax=1, sort="relevance"))
            ids = r["IdList"]
            if not ids:  # fallback without "chromosome" in title
                term2 = (f'"{species}"[Organism] AND complete genome[Title] '
                         f'AND srcdb_refseq[PROP] AND 3000000:7000000[SLEN]')
                ids = Entrez.read(Entrez.esearch(db="nuccore", term=term2, retmax=1))["IdList"]
            if not ids:
                print(f"    no chromosome for {species}")
                return None
            fh = Entrez.efetch(db="nuccore", id=ids[0], rettype="fasta", retmode="text")
            fa.write_text(fh.read())
            fh.close()
            time.sleep(0.4)
        except Exception as e:
            print(f"    chromosome fetch failed {species}: {e}")
            return None
    rec = SeqIO.read(fa, "fasta")
    return seq_to_array(str(rec.seq).upper())


def score_genes_with_baseline(contig, genes, baseline, window=2000, step=500):
    seq = contig_to_array(contig)
    windows = scan_sequence(seq, baseline, window=window, step=step)
    win = [(w.start, w.end, w.composite) for w in windows]
    scores = []
    for (gs, ge, _l, _amr) in genes:
        best = 0.0
        for (ws, we, comp) in win:
            if we > gs and ws < ge and comp > best:
                best = comp
        scores.append(best)
    return np.array(scores)


def main():
    print("=" * 86)
    print("  009: self-baseline  vs  host-chromosome baseline")
    print("=" * 86)
    paths = h.fetch_dataset()

    print("\n[1] Fetching one RefSeq chromosome per genus (cached)...")
    chrom_baselines = {}
    for sp in GENERA:
        arr = fetch_chromosome(sp)
        if arr is not None:
            t0 = time.time()
            chrom_baselines[sp] = compute_host_baseline(arr, window=2000, step=2000)
            print(f"    {sp:<24} {len(arr):>9} bp  baseline built ({time.time()-t0:.1f}s)")

    print(f"\n[2] Scoring plasmids with both baselines...")
    print(f"  {'replicon':<15}{'class':>6}{'self':>8}{'host':>8}{'delta':>8}  host-species")
    rows = []
    for p in paths:
        contig, genes = h.parse_genes(p)
        n_amr = sum(g[3] for g in genes)
        if n_amr == 0 or len(genes) - n_amr < 3:
            continue
        y = np.array([g[3] for g in genes], dtype=int)
        cls = p.stem.split("_")[0]

        # self baseline
        self_scores = h.score_genes(contig, genes)
        auc_self = h.auc_score(y, self_scores)

        # host-chromosome baseline
        gbrec = SeqIO.read(p, "genbank")
        sp = genus_key(gbrec.annotations.get("organism", "") + " " + gbrec.description)
        if sp not in chrom_baselines:
            continue
        host_scores = score_genes_with_baseline(contig, genes, chrom_baselines[sp])
        auc_host = h.auc_score(y, host_scores)

        rows.append((contig.id, cls, auc_self, auc_host, sp))
        d = auc_host - auc_self
        flag = "  <== adapted" if cls in ("ctxm", "oxa") else ""
        print(f"  {contig.id:<15}{cls:>6}{auc_self:>8.3f}{auc_host:>8.3f}{d:>+8.3f}  {sp}{flag}")

    self_a = np.array([r[2] for r in rows])
    host_a = np.array([r[3] for r in rows])
    print("\n" + "=" * 86)
    print("  RESULTS  (n = %d)" % len(rows))
    print("=" * 86)
    print(f"  Mean per-genome AUC   self={self_a.mean():.3f}   host-chrom={host_a.mean():.3f}"
          f"   delta={host_a.mean()-self_a.mean():+.3f}")
    print(f"  Median                self={np.median(self_a):.3f}   host-chrom={np.median(host_a):.3f}")
    won = (host_a > self_a + 0.02).sum()
    lost = (host_a < self_a - 0.02).sum()
    print(f"  host-chrom better on {won}/{len(rows)}, worse on {lost}")

    # the decisive subgroup: adapted genes that self-baseline missed
    adapted = [r for r in rows if r[1] in ("ctxm", "oxa")]
    if adapted:
        sa = np.array([r[2] for r in adapted]); ha = np.array([r[3] for r in adapted])
        print(f"\n  ADAPTED cases (CTX-M / OXA), where self-baseline failed (n={len(adapted)}):")
        print(f"    self mean={sa.mean():.3f}  ->  host-chrom mean={ha.mean():.3f}  "
              f"delta={ha.mean()-sa.mean():+.3f}")
        print("    If host-chrom lifts these above 0.5, the method change works on the")
        print("    exact blind spot. If not, composition can't reach adapted genes at all.")


if __name__ == "__main__":
    main()
