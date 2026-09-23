"""The dropout is barcode-driven. Which barcodes, and is it the oligo or the software?

Two mechanisms produce a barcode-specific dropout and they need different fixes:

  (a) the primer oligo never got into the well  -> no amplicon, no reads anywhere.
      Signature: failing indices map to contiguous wells of the primer SOURCE plate
      (a dried/empty region, a missed Echo shot), and their reads are simply absent.

  (b) minibar cannot assign the barcode     -> amplicon exists, reads land in
      'unassigned' or on a near-identical barcode.
      Signature: failing indices are sequence-similar to other indices, and the
      barcode's partner (the other end of the same read) shows no matching deficit.
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
import itertools, collections

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "02_demux_qc" / "outputs"

seq = pd.read_csv(OUT / "d00_sequenced_wells_with_counts.csv")
prim = pd.read_csv(config.SETUP / "primer_layout_20260907.csv")
print("primer_layout columns:", list(prim.columns))
print(prim.head(3).to_string())

def barcode_table(side):
    idx, src = f"{side}_idx", f"{side}_source_well"
    g = seq.groupby([idx, src]).agg(
        n=("dropout", "size"), k=("dropout", "sum"),
        reads=("read_count", "sum"), med=("read_count", "median")).reset_index()
    g["rate"] = g.k / g.n
    p = seq.dropout.mean()
    g["p_binom"] = [stats.binomtest(int(k), int(n), p).pvalue for k, n in zip(g.k, g.n)]
    g = g.sort_values("p_binom")
    g["q"] = np.minimum(1, g.p_binom * len(g) / np.arange(1, len(g) + 1))
    g["src_row"] = g[src].str[0]
    g["src_col"] = g[src].str[1:].astype(int)
    return g.sort_values("rate", ascending=False)

for side in ("fwd", "rev"):
    g = barcode_table(side)
    bad = g[(g.q < 0.05) & (g.rate > seq.dropout.mean())]
    print("\n" + "=" * 72)
    print(f"{side.upper()} barcodes with excess dropout (q<0.05)")
    print("=" * 72)
    print(bad[[f"{side}_idx", f"{side}_source_well", "n", "k", "rate", "reads", "med"]]
          .to_string(index=False))
    print(f"\n  {len(bad)} bad {side} barcodes; they cover {bad.n.sum()} wells "
          f"of which {bad.k.sum()} dropped out")
    # source-plate geometry: are the bad ones clustered on the primer source plate?
    print(f"  source-plate rows of bad {side} barcodes: "
          f"{collections.Counter(bad.src_row).most_common()}")
    print(f"  source-plate cols of bad {side} barcodes: "
          f"{sorted(collections.Counter(bad.src_col).items())}")
    allrows = collections.Counter(g.src_row)
    badrows = collections.Counter(bad.src_row)
    rowtab = pd.DataFrame({"n_bc": pd.Series(allrows), "n_bad": pd.Series(badrows)}).fillna(0)
    rowtab["frac_bad"] = rowtab.n_bad / rowtab.n_bc
    print("\n  bad-barcode fraction by primer source-plate row:")
    print(rowtab.to_string())
    g.to_csv(OUT / f"d03_{side}_barcode_forensics.csv", index=False)

# --- how much of the total dropout do the bad barcodes explain? ---
fg, rg = barcode_table("fwd"), barcode_table("rev")
p = seq.dropout.mean()
badf = set(fg[(fg.q < 0.05) & (fg.rate > p)].fwd_idx)
badr = set(rg[(rg.q < 0.05) & (rg.rate > p)].rev_idx)
seq["hits_bad_bc"] = seq.fwd_idx.isin(badf) | seq.rev_idx.isin(badr)
print("\n" + "=" * 72)
print("ATTRIBUTION")
print("=" * 72)
ct = pd.crosstab(seq.hits_bad_bc, seq.dropout)
print(ct.to_string())
sub = seq[seq.hits_bad_bc]; rest = seq[~seq.hits_bad_bc]
print(f"\nwells touching a flagged barcode : {len(sub):4d}  dropout {sub.dropout.mean():.1%}")
print(f"wells with two clean barcodes    : {len(rest):4d}  dropout {rest.dropout.mean():.1%}")
print(f"\nflagged barcodes are {len(badf)} fwd + {len(badr)} rev out of 96+96.")
print(f"They account for {sub.dropout.sum()} of {seq.dropout.sum()} dropouts "
      f"({sub.dropout.sum()/seq.dropout.sum():.0%}), while covering only "
      f"{len(sub)/len(seq):.0%} of wells.")
print(f"\nIf every flagged barcode had behaved like the clean ones, expected dropout "
      f"would be {rest.dropout.mean():.1%} -> "
      f"{int(round(rest.dropout.mean()*len(seq)))} wells instead of {seq.dropout.sum()}.")
seq.to_csv(OUT / "d00_sequenced_wells_with_counts.csv", index=False)
