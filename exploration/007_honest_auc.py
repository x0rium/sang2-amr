"""007: Honest gene-level AUC benchmark.

The README claims "AUC=0.667 vs self-genome" but no script ever produced it.
This experiment produces a reproducible number, the right way.

QUESTION
  On a plasmid/genome carrying resistance genes, does SANG2's composite
  anomaly score rank ACQUIRED AMR genes ABOVE ordinary genes on the SAME
  replicon? That is exactly the project's promise: surface AMR candidates
  for an experimentalist without a reference database.

DESIGN (no cherry-picking)
  1. Fetch RefSeq *complete* plasmids via esearch across several AMR classes
     (KPC, NDM, CTX-M, mcr, vanA, OXA-48, qnr, aminoglycoside). We take what
     the query returns, annotation and all.
  2. Parse every CDS + /product. Label each gene:
       AMR    = product matches an ACQUIRED-resistance regex (bla, mcr, van
                ligase, qnr, aminoglycoside-modifying, tet, erm, cat, sul, dfr)
       non-AMR= everything else (rep, mob, tra, par, hypothetical, ...)
     Ambiguous terms (generic "efflux", "acetyltransferase" w/o "aminoglycoside")
     are NOT counted as AMR — if unsure, it is a negative.
  3. Self-baseline scan (the prod default: window=2000, step=500).
     Each CDS gets score = max composite over windows overlapping it.
  4. Metrics:
       - per-genome AUC (AMR vs non-AMR within one replicon), then mean
       - pooled ROC over all genes
       - recall@top-k (does the AMR gene land in the top hits?)

Run:
  python exploration/007_honest_auc.py            # fetch (cached) + benchmark
  python exploration/007_honest_auc.py --fetch    # only download + report dataset
"""
from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from Bio import Entrez, SeqIO  # noqa: E402

from io_module import contig_to_array  # noqa: E402
from genome_signals import seq_to_array  # noqa: E402
from ir import Contig  # noqa: E402
from anomaly import compute_host_baseline, scan_sequence  # noqa: E402

Entrez.email = "tuxera@gmail.com"

DATA = Path(__file__).resolve().parent / "data" / "benchmark"
DATA.mkdir(parents=True, exist_ok=True)

# Each query targets a different AMR class so the dataset spans mechanisms.
# RefSeq + complete + bounded length keeps annotation quality high and scans fast.
QUERIES = {
    "kpc":  'blaKPC[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "ndm":  'blaNDM[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "ctxm": 'blaCTX-M[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "mcr":  'mcr-1[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "oxa":  'blaOXA-48[All Fields] AND plasmid[Title] AND complete sequence[Title]',
    "qnr":  'qnrS[All Fields] AND plasmid[Title] AND complete sequence[Title]',
}
SLEN = "10000:120000[SLEN]"          # avoid megaplasmids (slow) and tiny stubs
REFSEQ = 'srcdb_refseq[PROP]'        # curated annotation

# --- ACQUIRED AMR gene labelling -------------------------------------------
# Deliberately conservative: matches products of horizontally-acquired
# resistance genes, NOT chromosomal housekeeping that happens to share a fold.
AMR_RE = re.compile(
    r"""(
        beta-lactamase | carbapenemase | metallo-beta-lactamase | cephalosporinase |
        \bKPC\b | \bNDM\b | \bOXA-\d | \bCTX-M | \bTEM\b | \bSHV\b | \bVIM\b | \bIMP\b |
        \bGES\b | \bCMY\b | \bDHA\b | extended-spectrum |
        \bmcr- | colistin |
        vancomycin .* ligase | \bvanA\b | \bvanB\b | D-alanine--D-lactate |
        \bqnr | quinolone resistance |
        aminoglycoside .* (transferase|adenylyltransferase|nucleotidyltransferase) |
        \baac\( | \baph\( | \baad[AB] | \bant\( |
        tetracycline (resistance|efflux) | \btet\([A-Z]\) |
        \berm\b | rRNA .* methyltransferase .* macrolide | macrolide .* resistance |
        chloramphenicol acetyltransferase | \bcat[AB]? \b | florfenicol |
        sulfonamide | \bsul[123]\b |
        trimethoprim | \bdfr[AB] |
        rifampin .* (transferase|phosphotransferase) |
        antibiotic resistance | drug resistance protein
    )""",
    re.IGNORECASE | re.VERBOSE,
)
# products that look AMR-ish but are too generic to count as positives
AMBIGUOUS_RE = re.compile(
    r"(efflux|multidrug|transporter|acetyltransferase|hydrolase|^transferase)",
    re.IGNORECASE,
)


def is_amr(product: str) -> bool:
    if not product:
        return False
    if AMR_RE.search(product):
        return True
    return False


def fetch_dataset(per_query: int = 3) -> list[Path]:
    """esearch each AMR class, fetch RefSeq complete plasmids as GenBank. Cached."""
    paths: list[Path] = []
    seen_acc: set[str] = set()
    for tag, q in QUERIES.items():
        term = f"({q}) AND {REFSEQ} AND {SLEN}"
        try:
            h = Entrez.esearch(db="nuccore", term=term, retmax=per_query, sort="relevance")
            ids = Entrez.read(h)["IdList"]
            h.close()
        except Exception as e:
            print(f"  [{tag}] esearch failed: {e}")
            continue
        if not ids:
            print(f"  [{tag}] no hits")
            continue
        for uid in ids:
            gb = DATA / f"{tag}_{uid}.gb"
            if not gb.exists():
                try:
                    fh = Entrez.efetch(db="nuccore", id=uid, rettype="gbwithparts", retmode="text")
                    gb.write_text(fh.read())
                    fh.close()
                    time.sleep(0.4)  # be nice to NCBI
                except Exception as e:
                    print(f"  [{tag}] efetch {uid} failed: {e}")
                    continue
            rec = SeqIO.read(gb, "genbank")
            if rec.id in seen_acc:
                gb.unlink(missing_ok=True)
                continue
            seen_acc.add(rec.id)
            paths.append(gb)
            print(f"  [{tag}] {rec.id}  {len(rec.seq):>7} bp  {rec.description[:55]}")
    return paths


def parse_genes(gb_path: Path):
    """Return (Contig, [(start, end, product, is_amr), ...])."""
    rec = SeqIO.read(gb_path, "genbank")
    seq = str(rec.seq).upper()
    contig = Contig(id=rec.id, sequence=seq_to_array(seq).tobytes())
    genes = []
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        product = " ".join(feat.qualifiers.get("product", [""]))
        gene = " ".join(feat.qualifiers.get("gene", [""]))
        label = f"{gene} {product}".strip()
        start = int(feat.location.start)
        end = int(feat.location.end)
        if end - start < 150:        # skip fragments
            continue
        genes.append((start, end, label, is_amr(label)))
    return contig, genes


def score_genes(contig: Contig, genes, window=2000, step=500):
    """Self-baseline scan; assign each gene the max composite of overlapping windows."""
    seq = contig_to_array(contig)
    baseline = compute_host_baseline(seq, window=window, step=window)
    windows = scan_sequence(seq, baseline, window=window, step=step)
    # windows are deduped & sorted; keep (start, end, composite)
    win = [(w.start, w.end, w.composite) for w in windows]
    scores = []
    for (gs, ge, label, amr) in genes:
        best = 0.0
        for (ws, we, comp) in win:
            if we > gs and ws < ge:          # overlap
                if comp > best:
                    best = comp
        scores.append(best)
    return np.array(scores)


def auc_score(y_true: np.ndarray, y_score: np.ndarray) -> float | None:
    """Rank-based AUC (Mann-Whitney). None if one class absent."""
    pos = y_score[y_true == 1]
    neg = y_score[y_true == 0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    # U statistic via rank-sum
    order = np.argsort(np.concatenate([pos, neg]), kind="mergesort")
    ranks = np.empty(len(order), dtype=float)
    ranks[order] = np.arange(1, len(order) + 1)
    # average ties
    allv = np.concatenate([pos, neg])
    for v in np.unique(allv):
        m = allv == v
        if m.sum() > 1:
            ranks[m] = ranks[m].mean()
    r_pos = ranks[: len(pos)].sum()
    u = r_pos - len(pos) * (len(pos) + 1) / 2
    return float(u / (len(pos) * len(neg)))


def main(fetch_only=False):
    print("=" * 64)
    print("  007: Honest gene-level AUC benchmark")
    print("=" * 64)
    print("\n[1] Fetching RefSeq complete AMR plasmids (cached)...")
    paths = fetch_dataset()
    if not paths:
        print("No data fetched. Check network / NCBI.")
        return
    if fetch_only:
        return

    print(f"\n[2] Scoring {len(paths)} replicons (self-baseline, window=2000 step=500)...")
    all_true: list[int] = []
    all_score: list[float] = []
    per_genome = []
    recall_at = {1: 0, 3: 0, 5: 0, 10: 0}
    n_with_amr = 0

    for p in paths:
        contig, genes = parse_genes(p)
        n_amr = sum(g[3] for g in genes)
        if n_amr == 0 or len(genes) - n_amr < 3:
            print(f"  {contig.id:<14} skipped (amr={n_amr}, total={len(genes)})")
            continue
        t0 = time.time()
        scores = score_genes(contig, genes)
        y = np.array([g[3] for g in genes], dtype=int)
        a = auc_score(y, scores)
        per_genome.append((contig.id, a, n_amr, len(genes)))
        all_true.extend(y.tolist())
        all_score.extend(scores.tolist())
        n_with_amr += 1
        # recall@k: is any AMR gene among the top-k scored genes?
        rank_order = np.argsort(-scores, kind="mergesort")
        amr_positions = [int(np.where(rank_order == i)[0][0]) for i in np.where(y == 1)[0]]
        best_pos = min(amr_positions) + 1 if amr_positions else 999
        for k in recall_at:
            if best_pos <= k:
                recall_at[k] += 1
        print(f"  {contig.id:<14} AUC={a:.3f}  amr={n_amr:<2} genes={len(genes):<3} "
              f"bestAMR_rank={best_pos}  ({time.time()-t0:.1f}s)")

    print("\n" + "=" * 64)
    print("  RESULTS")
    print("=" * 64)
    aucs = [a for (_, a, _, _) in per_genome if a is not None]
    pooled = auc_score(np.array(all_true), np.array(all_score))
    print(f"  Replicons with AMR + >=3 controls : {n_with_amr}")
    print(f"  Total genes scored                : {len(all_true)} "
          f"({sum(all_true)} AMR / {len(all_true)-sum(all_true)} non-AMR)")
    print(f"  Mean per-genome AUC               : {np.mean(aucs):.3f} "
          f"(median {np.median(aucs):.3f}, n={len(aucs)})")
    print(f"  Pooled AUC                        : {pooled:.3f}")
    print(f"  README claim                      : 0.667")
    print(f"\n  Recall@k (AMR gene among top-k anomalous genes):")
    for k in sorted(recall_at):
        print(f"    @{k:<3}: {recall_at[k]}/{n_with_amr} = {recall_at[k]/max(n_with_amr,1):.0%}")
    print("\n  Interpretation: AUC=0.5 random, 0.7 weak, 0.8 useful, 0.9 strong.")


if __name__ == "__main__":
    main(fetch_only="--fetch" in sys.argv)
