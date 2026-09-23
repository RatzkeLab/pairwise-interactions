"""Oligo or software? Compare the same barcode indices across three runs.

A barcode SEQUENCE that minibar cannot resolve is a property of the sequence, so
it fails identically in every run that uses it. A barcode WELL that has run dry,
or that the Echo missed, fails in one run and not the next. The three runs share
one primer plate and one index->sequence mapping, so the comparison is direct.
"""
import re, collections
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats

BASE = Path("/home/rl/scripts/karl/pairwise_interaction_experiments")
OUT = BASE / "20260907" / "02_demux_qc" / "outputs"
FWD_PRIMER, REV_PRIMER = "AGRGTTYGATYMTGGCTCAG", "CGGYTACCTTGTTACGACTT"
ADAPTER = "ATCGCCTACCGTGAC"


def barcodes_from_minibar(path):
    """sample -> (fwd_barcode, rev_barcode), stripping adapter and 16S primer."""
    out = {}
    df = pd.read_csv(path, sep="\t")
    for _, r in df.iterrows():
        f = r.FwIndex[len(ADAPTER):-len(FWD_PRIMER)]
        v = r.RvIndex[len(ADAPTER):-len(REV_PRIMER)]
        out[r.SampleID] = (f, v)
    return out


# --- index -> barcode sequence for this run -------------------------------------------
bc07 = barcodes_from_minibar(BASE / "20260907/01_setup/minibar_primers_20260907.tsv")
lay07 = pd.read_csv(BASE / "20260907/01_setup/strain_layout_20260907.csv")
lay07["sample"] = lay07.apply(lambda r: f"Plate{int(r.dest_plate):02d}_{r.dest_well}", axis=1)
fwd_seq, rev_seq = {}, {}
for _, r in lay07.iterrows():
    f, v = bc07[r["sample"]]
    fwd_seq.setdefault(r.fwd_idx, f); rev_seq.setdefault(r.rev_idx, v)
print(f"barcode length: fwd {set(map(len, fwd_seq.values()))}, rev {set(map(len, rev_seq.values()))}")
print(f"distinct fwd barcodes {len(set(fwd_seq.values()))}/{len(fwd_seq)}, "
      f"rev {len(set(rev_seq.values()))}/{len(rev_seq)}")

# --- per-barcode-SEQUENCE dropout in each run -----------------------------------------
RUNS = {
    "20260630": (BASE / "20260630", "20260710_demultiplex"),
    "20260721": (BASE / "20260721", "20260730_demux"),
    "20260907": (BASE / "20260907", "20260907_demux"),
}
DEMUXROOT = Path("/home/rl/scripts/karl/data_links/interim/demultiplexing")

def run_barcode_rates(minibar_tsv, demux_dir, used_samples=None):
    bc = barcodes_from_minibar(minibar_tsv)
    cnt = pd.read_csv(DEMUXROOT / demux_dir / "summary/demultiplexed_read_counts.tsv", sep="\t")
    cnt = cnt.set_index("sample").read_count
    rows = []
    for s, (f, v) in bc.items():
        if s not in cnt.index: continue
        if used_samples is not None and s not in used_samples: continue
        rows.append((s, f, v, int(cnt[s])))
    d = pd.DataFrame(rows, columns=["sample", "fwd_bc", "rev_bc", "reads"])
    d["dropout"] = d.reads <= 10
    return d

# locate each run's minibar file
mb = {}
for run, (root, _) in RUNS.items():
    hits = sorted(root.rglob("minibar_primers_*.tsv"))
    hits = [h for h in hits if "primary_only" not in h.name]
    mb[run] = hits[0] if hits else None
    print(f"{run}: {mb[run]}")

# --- only wells that were actually PCR'd count as evidence ------------------------------
USED = {}
for run, (root, _) in RUNS.items():
    lay = pd.read_csv(sorted(root.rglob(f"strain_layout_{run}.csv"))[0])
    lay["sample"] = lay.apply(lambda r: f"Plate{int(r.dest_plate):02d}_{r.dest_well}", axis=1)
    if "plate_role" in lay.columns:          # 20260907 marks the frozen extension plates
        lay = lay[lay.plate_role == "primary_sequenced"]
    else:                                    # earlier runs used only the odd-numbered plates
        cnt = pd.read_csv(DEMUXROOT / RUNS[run][1] / "summary/demultiplexed_read_counts.tsv", sep="\t")
        cnt = cnt.set_index("sample").read_count
        lay["reads"] = lay["sample"].map(cnt).fillna(-1)
        byplate = lay.groupby("dest_plate").reads.sum()
        live = byplate[byplate > 0.05 * byplate.max()].index
        print(f"{run}: plates with real signal -> {list(live)}")
        lay = lay[lay.dest_plate.isin(live)]
    USED[run] = set(lay["sample"])
    print(f"{run}: {len(USED[run])} sequenced wells")

per_run = {}
for run in RUNS:
    d = run_barcode_rates(mb[run], RUNS[run][1], USED[run])
    per_run[run] = d
    print(f"{run}: {len(d)} wells, dropout {d.dropout.mean():.1%}")

# --- join on the barcode SEQUENCE (index numbering may differ between runs) --------------
def rate_by(d, col):
    g = d.groupby(col).agg(n=("dropout", "size"), k=("dropout", "sum"))
    g = g[g.n >= 8]
    return (g.k / g.n).rename("rate"), g

print("\n" + "=" * 74)
print("SAME BARCODE SEQUENCE, THREE RUNS -- does a bad barcode stay bad?")
print("=" * 74)
for side in ("fwd_bc", "rev_bc"):
    tab = pd.DataFrame({run: rate_by(per_run[run], side)[0] for run in RUNS}).dropna()
    print(f"\n--- {side}  ({len(tab)} barcodes present in all three runs) ---")
    for a, b in [("20260630", "20260721"), ("20260630", "20260907"), ("20260721", "20260907")]:
        r, p = stats.spearmanr(tab[a], tab[b])
        print(f"  dropout-rate correlation {a} vs {b}: rho={r:+.3f}  p={p:.2g}")
    bad07 = tab[tab["20260907"] > 0.5]
    print(f"\n  {len(bad07)} barcodes failing hard (>50% dropout) in 20260907;"
          f" their rates in the earlier runs:")
    if len(bad07):
        print(bad07.sort_values("20260907", ascending=False)
              .to_string(float_format=lambda v: f"{v:.2f}"))
        print(f"\n  mean rate of those barcodes: 20260630 {bad07['20260630'].mean():.1%}, "
              f"20260721 {bad07['20260721'].mean():.1%}, 20260907 {bad07['20260907'].mean():.1%}")
        print(f"  all other barcodes  20260907: {tab[tab['20260907']<=0.5]['20260907'].mean():.1%}")
    tab.to_csv(OUT / f"d04_crossrun_{side}_dropout.csv")
