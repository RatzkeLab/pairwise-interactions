"""A run-specific mechanism with a free fix: the sample sheet doubled in size.

20260630/20260721 demultiplexed against 3080 barcode pairs. 20260907 was
demultiplexed against all 4928 -- including the 2464 extension wells that were
never PCR'd. minibar matches each read end within edit distance 6; when an
errored read falls within 6 of two different sample rows it is ambiguous and
lost. Adding 2464 unused rows can only add ambiguity, and which barcode suffers
depends on the new design -- exactly the run-specific pattern we observe.

If this is the cause, the fix is free: re-demultiplex against
minibar_primers_20260907_primary_only.tsv.
"""
import numpy as np, pandas as pd, edlib, itertools, collections
from pathlib import Path
from scipy import stats

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "02_demux_qc/outputs"
ADAPTER, FWDP, REVP = "ATCGCCTACCGTGAC", "AGRGTTYGATYMTGGCTCAG", "CGGYTACCTTGTTACGACTT"
ED = 6   # minibar barcode_edit_dist from demultiplex_config.yaml

mb = pd.read_csv(config.SETUP / "minibar_primers_20260907.tsv", sep="\t")
mb["fwd_bc"] = mb.FwIndex.str[len(ADAPTER):-len(FWDP)]
mb["rev_bc"] = mb.RvIndex.str[len(ADAPTER):-len(REVP)]

fwd_bcs = sorted(set(mb.fwd_bc)); rev_bcs = sorted(set(mb.rev_bc))
def dmat(bcs):
    n = len(bcs); M = np.zeros((n, n), int)
    for i in range(n):
        for j in range(i + 1, n):
            M[i, j] = M[j, i] = edlib.align(bcs[i], bcs[j], mode="NW",
                                            task="distance")["editDistance"]
    return M
Mf, Mr = dmat(fwd_bcs), dmat(rev_bcs)
fi = {b: i for i, b in enumerate(fwd_bcs)}; ri = {b: i for i, b in enumerate(rev_bcs)}
off = np.triu_indices(len(fwd_bcs), 1)
print(f"barcode-to-barcode edit distances (24 bp, minibar tolerance {ED}):")
print(f"  fwd: min={Mf[off].min()} 1st pct={np.percentile(Mf[off],1):.0f} "
      f"median={np.median(Mf[off]):.0f}")
print(f"  rev: min={Mr[off].min()} 1st pct={np.percentile(Mr[off],1):.0f} "
      f"median={np.median(Mr[off]):.0f}")
print(f"  pairs within 2*{ED}={2*ED} edits (mutually confusable): "
      f"fwd {(Mf[off] <= 2*ED).sum()}, rev {(Mr[off] <= 2*ED).sum()}")

# --- per well: count confusable sample rows, in the full sheet vs primary-only ----------
lay = pd.read_csv(config.LAYOUT_CSV)
lay["sample"] = lay.apply(lambda r: f"Plate{int(r.dest_plate):02d}_{r.dest_well}", axis=1)
mb = mb.merge(lay[["sample", "plate_role"]], left_on="SampleID", right_on="sample")
F = np.array([fi[b] for b in mb.fwd_bc]); R = np.array([ri[b] for b in mb.rev_bc])
prim = (mb.plate_role == "primary_sequenced").values

def n_confusable(mask_sheet):
    """For each primary well, how many OTHER rows of `mask_sheet` are within ED on
    BOTH ends (i.e. a read could match either)."""
    idx = np.where(prim)[0]
    sheetF, sheetR = F[mask_sheet], R[mask_sheet]
    out = np.zeros(len(idx), int)
    for n, i in enumerate(idx):
        close_f = Mf[F[i], sheetF] <= 2 * ED
        close_r = Mr[R[i], sheetR] <= 2 * ED
        out[n] = (close_f & close_r).sum() - 1
    return out

full = n_confusable(np.ones(len(mb), bool))
only = n_confusable(prim)
cnt = pd.read_csv("/home/rl/scripts/karl/data_links/interim/demultiplexing/"
                  "20260907_demux/summary/demultiplexed_read_counts.tsv", sep="\t"
                  ).set_index("sample").read_count
d = pd.DataFrame({"sample": mb.SampleID[prim].values,
                  "confusable_full": full, "confusable_primary_only": only})
d["reads"] = d["sample"].map(cnt); d["dropout"] = d.reads <= 10
d["extra_from_extension"] = d.confusable_full - d.confusable_primary_only

print(f"\nconfusable sample rows per sequenced well:")
print(f"  against the full 4928-row sheet   : mean {full.mean():.2f}  "
      f"max {full.max()}  wells with >=1: {(full>0).sum()}")
print(f"  against the 2464-row primary sheet: mean {only.mean():.2f}  "
      f"max {only.max()}  wells with >=1: {(only>0).sum()}")
print(f"  extra ambiguity contributed by the frozen extension rows: "
      f"{d.extra_from_extension.sum()} well-collisions across "
      f"{(d.extra_from_extension>0).sum()} wells")

print("\ndropout vs ambiguity introduced by the unused extension rows:")
tab = d.groupby(d.extra_from_extension.clip(0, 3)).agg(
    wells=("dropout", "size"), dropout=("dropout", "mean"))
print(tab.to_string(float_format=lambda v: f"{v:.1%}"))
if d.extra_from_extension.max() > 0:
    a = d[d.extra_from_extension > 0].dropout
    b = d[d.extra_from_extension == 0].dropout
    _, p = stats.fisher_exact([[a.sum(), len(a)-a.sum()], [b.sum(), len(b)-b.sum()]])
    print(f"  Fisher p = {p:.3g}")
d.to_csv(OUT / "d06_barcode_collisions.csv", index=False)
