"""If the failure is physical, it should be visible on the primer SOURCE plate.

Two physical stories, distinguishable:
  depletion  -> failure tracks cumulative draws from that source well across runs,
                and worsens monotonically 20260630 -> 20260721 -> 20260907.
  local      -> failing wells cluster spatially on the source plate (an evaporated
                edge, a mis-seated region, a quadrant the Echo mis-shot).
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
import collections

BASE = Path("/home/rl/scripts/karl/pairwise_interaction_experiments")
OUT = BASE / "20260907/02_demux_qc/outputs"
ROWS = "ABCDEFGHIJKLMNOP"
DEMUXROOT = Path("/home/rl/scripts/karl/data_links/interim/demultiplexing")
RUNS = {"20260630": "20260710_demultiplex", "20260721": "20260730_demux",
        "20260907": "20260907_demux"}


def run_table(run, demux):
    lay = pd.read_csv(sorted((BASE / run).rglob(f"strain_layout_{run}.csv"))[0])
    lay["sample"] = lay.apply(lambda r: f"Plate{int(r.dest_plate):02d}_{r.dest_well}", axis=1)
    cnt = pd.read_csv(DEMUXROOT / demux / "summary/demultiplexed_read_counts.tsv",
                      sep="\t").set_index("sample").read_count
    lay["reads"] = lay["sample"].map(cnt)
    lay = lay[lay.reads.notna()].copy()
    if "plate_role" in lay.columns:
        lay = lay[lay.plate_role == "primary_sequenced"]
    else:
        bp = lay.groupby("dest_plate").reads.sum()
        lay = lay[lay.dest_plate.isin(bp[bp > 0.05 * bp.max()].index)]
    lay["dropout"] = lay.reads <= 10
    lay["run"] = run
    return lay


all_runs = pd.concat([run_table(r, d) for r, d in RUNS.items()], ignore_index=True)

# --- one row per (run, source well), both primer orientations --------------------------
recs = []
for side in ("fwd", "rev"):
    col = f"{side}_source_well"
    g = all_runs.groupby(["run", col]).agg(n=("dropout", "size"),
                                           k=("dropout", "sum")).reset_index()
    g = g.rename(columns={col: "src_well"}); g["side"] = side
    recs.append(g)
src = pd.concat(recs, ignore_index=True)
src["rate"] = src.k / src.n
src["src_row"] = src.src_well.str[0]
src["src_col"] = src.src_well.str[1:].astype(int)
src["row_i"] = src.src_row.map({c: i for i, c in enumerate(ROWS)})

# --- 1. depletion: cumulative draws before this run ------------------------------------
order = ["20260630", "20260721", "20260907"]
src["run_i"] = src.run.map({r: i for i, r in enumerate(order)})
src = src.sort_values(["src_well", "side", "run_i"])
src["cum_draws_before"] = (src.groupby(["src_well", "side"]).n
                           .cumsum().sub(src.n))
print("=" * 74)
print("1. DEPLETION -- does failure track cumulative use of the source well?")
print("=" * 74)
for run in order:
    s = src[src.run == run]
    print(f"  {run}: mean dropout {s.rate.mean():.1%}   "
          f"median cumulative prior draws {s.cum_draws_before.median():.0f}")
s907 = src[src.run == "20260907"]
r, p = stats.spearmanr(s907.cum_draws_before, s907.rate)
print(f"\n  within 20260907: rho(prior draws, dropout rate) = {r:+.3f}  p={p:.2g}  "
      f"(n={len(s907)} source wells)")
print("  -> a source well that has been drawn from more is "
      f"{'more' if r>0 else 'not more'} likely to fail")

# --- 2. spatial clustering on the source plate -----------------------------------------
print("\n" + "=" * 74)
print("2. LOCAL -- do the 20260907 failures cluster on the primer source plate?")
print("=" * 74)
bad = s907[s907.rate > 0.5]
print(f"  {len(bad)} of {len(s907)} source wells fail >50%")

# row / column enrichment
for axis, lab in [("src_row", "row"), ("src_col", "column")]:
    tab = pd.crosstab(s907[axis], s907.rate > 0.5)
    if True in tab.columns and tab.shape[0] > 1:
        chi2, pv, _, _ = stats.chi2_contingency(tab.values)
        print(f"  source-plate {lab:7s}: chi2={chi2:.1f} p={pv:.3g}")
        worst = (tab[True] / tab.sum(axis=1)).sort_values(ascending=False)
        print("     worst: " + ", ".join(f"{k}={v:.0%}" for k, v in worst.head(4).items()))

# nearest-neighbour: are failing wells adjacent to each other more than chance?
def adjacency_enrichment(sub, n_perm=5000, seed=0):
    pts = sub[["row_i", "src_col"]].values
    lab = (sub.rate > 0.5).values
    d = np.abs(pts[:, None, :] - pts[None, :, :]).max(axis=2)
    # same orientation rows are 2 apart, so 'neighbour' = within 2 rows and 2 cols
    nb = (d <= 2) & (d > 0)
    obs = (nb & lab[:, None] & lab[None, :]).sum() / 2
    exp_n = nb.sum() / 2 * lab.mean() ** 2
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        pl = rng.permutation(lab)
        null[i] = (nb & pl[:, None] & pl[None, :]).sum() / 2
    pv = (null >= obs).mean()
    return obs, null.mean(), pv

for side in ("fwd", "rev"):
    sub = s907[s907.side == side]
    obs, exp, pv = adjacency_enrichment(sub)
    print(f"  {side}: failing-failing neighbour pairs obs={obs:.0f} "
          f"exp={exp:.1f} enrichment={obs/max(exp,1e-9):.2f}x  p={pv:.4f}")

src.to_csv(OUT / "d05_source_well_dropout.csv", index=False)

# --- 3. the printable map ---------------------------------------------------------------
print("\n" + "=" * 74)
print("3. PRIMER SOURCE PLATE MAP, 20260907 dropout rate per well")
print("   (odd rows = forward primers, even rows = reverse; . = <20%, "
      "digit = tens%, X = 100%)")
print("=" * 74)
grid = {(r.src_row, r.src_col): r.rate for r in s907.itertuples()}
cols = sorted(s907.src_col.unique())
print("     " + "".join(f"{c:>3}" for c in cols))
for rw in ROWS:
    if not any((rw, c) in grid for c in cols): continue
    line = f"  {rw}  "
    for c in cols:
        v = grid.get((rw, c))
        if v is None: line += "  -"
        elif v >= 1.0: line += "  X"
        elif v < 0.2: line += "  ."
        else: line += f" {int(v*10):>2}"
    print(line)
