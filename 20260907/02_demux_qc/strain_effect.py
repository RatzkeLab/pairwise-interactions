"""After the barcode effect is removed, is there a strain effect left?

The raw strain dispersion (1.49) is inflated because a strain that happens to
draw a dead barcode looks like a dead strain. Restrict to wells whose two
barcodes are both clean, and the question becomes fair.
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats

OUT = Path("/home/rl/scripts/karl/pairwise_interaction_experiments/20260907/02_demux_qc/outputs")
SETUP = Path("/home/rl/scripts/karl/pairwise_interaction_experiments/20260907/01_setup")
NEVER_SEQUENCED = "A4 A8 L14 N16 N2 N5 O16 P9".split()   # from the design README

seq = pd.read_csv(OUT / "d00_sequenced_wells_with_counts.csv")
clean = seq[~seq.hits_bad_bc].copy()
print(f"{len(clean)} wells with two clean barcodes; dropout {clean.dropout.mean():.1%}")

mem = pd.concat([
    clean[["sample", "dropout", "read_count", "well_type", "strain1"]]
        .rename(columns={"strain1": "strain"}),
    clean[["sample", "dropout", "read_count", "well_type", "strain2"]]
        .rename(columns={"strain2": "strain"}),
]).dropna(subset=["strain"])

g = mem.groupby("strain").agg(wells=("dropout", "size"), fails=("dropout", "sum"),
                              median_reads=("read_count", "median"))
g = g[g.wells >= 4]
p = clean.dropout.mean()
g["rate"] = g.fails / g.wells
exp = g.wells * p
chi2 = (((g.fails - exp) ** 2) / (exp * (1 - p))).sum()
dof = len(g) - 1
print(f"\nSTRAIN effect among clean-barcode wells ({len(g)} strains with >=4 wells):")
print(f"  chi2={chi2:.0f} df={dof} dispersion={chi2/dof:.2f} p={stats.chi2.sf(chi2,dof):.3g}")
print(f"  (raw dispersion before removing bad barcodes was 1.49)")

g["p_binom"] = [stats.binomtest(int(k), int(n), p).pvalue for k, n in zip(g.fails, g.wells)]
order = np.argsort(g.p_binom.values); ranks = np.empty(len(g), int)
ranks[order] = np.arange(1, len(g) + 1)
g["q"] = np.minimum(1, g.p_binom.values * len(g) / ranks)
bad = g[(g.q < 0.1) & (g.rate > p)].sort_values("rate", ascending=False)
print(f"\n  strains failing more than chance (q<0.1): {len(bad)}")
if len(bad):
    print(bad[["wells", "fails", "rate", "median_reads", "q"]].to_string())

print(f"\nThe 8 strains never recovered in six previous ONT attempts:")
sub = g.reindex(NEVER_SEQUENCED)
print(sub[["wells", "fails", "rate", "median_reads"]].to_string())
present = sub.dropna()
if len(present):
    print(f"\n  their dropout rate: {present.fails.sum()/present.wells.sum():.1%} "
          f"vs {p:.1%} for everything else")
    k, n = int(present.fails.sum()), int(present.wells.sum())
    print(f"  binomial p = {stats.binomtest(k, n, p).pvalue:.3g}")
    rescued = present[present.rate < 0.5].index.tolist()
    print(f"  strains that finally produced reads this run: {rescued}")

# mono wells get double inoculum -- do they fare better?
print(f"\nby well type (clean-barcode wells only):")
wt = clean.groupby("well_type").agg(wells=("dropout", "size"),
                                    dropout=("dropout", "mean"),
                                    median_reads=("read_count", "median"))
print(wt.to_string(float_format=lambda v: f"{v:.3f}"))
g.to_csv(OUT / "d10_strain_effect_clean_barcodes.csv")
