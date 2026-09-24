"""Original vs corrected demultiplexing of 20260630, through identical code.

    python compare_arms.py        (env karl_seq_analysis; after run_all.py for both arms)

Writes comparison/*.csv and comparison/figures/*.png. Sections:
  c01  demux: reads assigned, even-plate (never sequenced) false-assignment control
  c02  per-well depth and cohort size on the sequenced plates
  c03  per-well agreement: same well, relative abundance under both demuxes
  c04  mapping_validation qc_status
  c05  relative abundance / hierarchy
  c06  genomic ML: dataset, noise ceiling, cross-validation, genome -> strength
  c07  plate-reader tier and read-depth curve
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
import config  # noqa: E402

OUT = BASE / "comparison"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
ARMS = ["original", "corrected"]
DEMUX_RUN = {a: config.DEMUX_DIRS[a].parent for a in ARMS}
pd.set_option("display.width", 200)


def od(step, arm):
    return config.out_dir(step, arm)


def header(t):
    print(f"\n{'=' * 90}\n{t}\n{'=' * 90}")


# --- c01 demux --------------------------------------------------------------------------
header("c01 demultiplexing")
rows = []
for a in ARMS:
    s = DEMUX_RUN[a] / "summary"
    raw = pd.read_csv(s / "raw_read_counts.tsv", sep="\t").iloc[:, 1].sum()
    filt = pd.read_csv(s / "filtered_read_counts.tsv", sep="\t").iloc[:, 1].sum()
    d = pd.read_csv(s / "demux_summary.tsv", sep="\t")
    d = d[d["sample"].str.match(r"Plate\d\d_")]
    d["plate"] = d["sample"].str[5:7].astype(int)
    seq = d.plate.isin(config.SEQUENCED_PLATES)
    n_seq, n_uns = d.loc[seq, "reads"].sum(), d.loc[~seq, "reads"].sum()
    # Barcode combinations are drawn from the same 96x96 pool on every plate, so a random
    # mis-assignment is as likely to land on an unsequenced well as a sequenced one. The
    # false rate on sequenced wells is therefore estimated from the unsequenced-well rate
    # per well, scaled to the number of sequenced wells.
    per_well_false = n_uns / (~seq).sum()
    est_false_on_seq = per_well_false * seq.sum()
    rows.append(dict(arm=a, raw_reads=raw, filtered_reads=filt,
                     assigned_reads=n_seq + n_uns,
                     pct_of_filtered=100 * (n_seq + n_uns) / filt,
                     reads_on_sequenced_plates=n_seq,
                     reads_on_unsequenced_plates=n_uns,
                     est_false_reads_on_sequenced=est_false_on_seq,
                     est_false_rate_pct=100 * est_false_on_seq / n_seq,
                     median_reads_per_sequenced_well=d.loc[seq, "reads"].median(),
                     max_reads_unsequenced_well=d.loc[~seq, "reads"].max()))
c01 = pd.DataFrame(rows).set_index("arm")
c01.loc["fold_change"] = c01.loc["corrected"] / c01.loc["original"]
c01.to_csv(OUT / "c01_demux.csv")
print(c01.T.round(3).to_string())

# --- c02 per-well depth -----------------------------------------------------------------
header("c02 per-well depth and cohort (sequenced plates)")
depth = {}
for a in ARMS:
    s = pd.read_csv(od("04_qc/mapping_validation", a) / "01_samples_gt5reads.csv")
    depth[a] = s.set_index("sample_id")
d = pd.read_csv(DEMUX_RUN["corrected"] / "summary" / "demux_summary.tsv", sep="\t")
d0 = pd.read_csv(DEMUX_RUN["original"] / "summary" / "demux_summary.tsv", sep="\t")
dd = d.merge(d0, on="sample", suffixes=("_corr", "_orig"))
dd = dd[dd["sample"].str.match(r"Plate\d\d_") & dd["sample"].str[5:7].astype(int).isin(config.SEQUENCED_PLATES)]
c02 = pd.DataFrame({
    a: {
        "wells_gt5_reads": len(depth[a]),
        "wells_ge10_reads": int((dd[f"reads_{k}"] >= 10).sum()),
        "wells_ge50_reads": int((dd[f"reads_{k}"] >= 50).sum()),
        "median_depth_in_cohort": depth[a].n_reads.median(),
        "mean_depth_in_cohort": depth[a].n_reads.mean(),
    } for a, k in (("original", "orig"), ("corrected", "corr"))}).T
c02.to_csv(OUT / "c02_depth.csv")
print(c02.to_string())
nz = dd[dd.reads_orig > 0]
ratio = nz.reads_corr / nz.reads_orig
print(f"per-well fold change (wells with >0 original reads): median {ratio.median():.2f}, "
      f"IQR {ratio.quantile(.25):.2f}-{ratio.quantile(.75):.2f}; "
      f"wells that LOST reads: {(nz.reads_corr < nz.reads_orig).sum()}")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
bins = np.logspace(0, np.log10(max(dd.reads_corr.max(), 2)), 40)
ax[0].hist(dd.reads_orig.clip(lower=1), bins=bins, alpha=.6, label="original")
ax[0].hist(dd.reads_corr.clip(lower=1), bins=bins, alpha=.6, label="corrected")
ax[0].set_xscale("log"); ax[0].set_xlabel("reads per well (sequenced plates)"); ax[0].set_ylabel("wells")
ax[0].legend(); ax[0].set_title("Per-well depth")
ax[1].scatter(dd.reads_orig + 1, dd.reads_corr + 1, s=4, alpha=.3)
lim = [1, dd.reads_corr.max() + 1]
ax[1].plot(lim, lim, "k--", lw=.8, label="1x")
ax[1].plot(lim, [4 * x for x in lim], ":", c="grey", lw=.8, label="4x")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].legend()
ax[1].set_xlabel("original reads + 1"); ax[1].set_ylabel("corrected reads + 1"); ax[1].set_title("Same well, both demuxes")
fig.tight_layout(); fig.savefig(FIG / "c02_depth.png", dpi=150); plt.close(fig)

# --- c03 per-well relative abundance agreement ------------------------------------------
header("c03 per-well relative abundance, same well under both demuxes")
w = {a: pd.read_csv(od("05_engineer_relative_abundances/relative_abundance", a) /
                    "r02_well_interaction_scores.csv").set_index("sample_id") for a in ARMS}
common = w["original"].index.intersection(w["corrected"].index)
m = w["original"].loc[common].join(w["corrected"].loc[common], lsuffix="_o", rsuffix="_c")
m = m[(~m.missing_reference_o.astype(bool)) & m.relative_abundance_a_o.notna() & m.relative_abundance_a_c.notna()]
sep = m.ref_pair_bp_dist_o >= 10
diff = (m.relative_abundance_a_c - m.relative_abundance_a_o)
c03 = pd.DataFrame([{
    "subset": name, "n_wells": int(mask.sum()),
    "pearson_r": pearsonr(m.relative_abundance_a_o[mask], m.relative_abundance_a_c[mask])[0],
    "median_abs_diff": diff[mask].abs().median(),
    "p95_abs_diff": diff[mask].abs().quantile(.95),
    "mean_diff_corr_minus_orig": diff[mask].mean(),
    "winner_flips_pct": 100 * ((m.relative_abundance_a_o[mask] - .5) * (m.relative_abundance_a_c[mask] - .5) < 0).mean(),
} for name, mask in (("all", m.index == m.index), (">=10bp refs", sep), ("<10bp refs", ~sep))])
c03.to_csv(OUT / "c03_well_agreement.csv", index=False)
print(c03.round(4).to_string(index=False))

fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(m.relative_abundance_a_o[sep], m.relative_abundance_a_c[sep], s=5, alpha=.35,
           c=np.log10(m.n_reads_o[sep]), cmap="viridis")
ax.plot([0, 1], [0, 1], "k--", lw=.8)
ax.set_xlabel("relative abundance (original demux)"); ax.set_ylabel("relative abundance (corrected demux)")
ax.set_title(f"Same well, both demuxes (refs >=10 bp apart, n={sep.sum()})\ncolour = log10 original depth", fontsize=10)
fig.tight_layout(); fig.savefig(FIG / "c03_well_agreement.png", dpi=150); plt.close(fig)

# --- c04 mapping validation -------------------------------------------------------------
header("c04 mapping_validation qc_status")
qc = {a: pd.read_csv(od("04_qc/mapping_validation", a) / "05_combined_sample_summary.csv") for a in ARMS}
c04 = pd.concat({a: qc[a].groupby(["well_type", "qc_status"]).size() for a in ARMS}, axis=1).fillna(0).astype(int)
c04.to_csv(OUT / "c04_qc_status.csv")
print(c04.to_string())
for a in ARMS:
    c = pd.read_csv(od("04_qc/mapping_validation", a) / "05_contamination_candidates.csv")
    print(f"  {a:9s}: {len(c)} contamination candidates, {int((~c.likely_reference_twin).sum())} not reference-twins")

# --- c05 relative abundance / hierarchy -------------------------------------------------
header("c05 relative abundance and hierarchy")
rows = {}
for a in ARMS:
    r = od("05_engineer_relative_abundances/relative_abundance", a)
    p = pd.read_csv(r / "r03_pair_replicate_stats.csv")
    h = pd.read_csv(r / "r05_hierarchy_summary.csv").set_index("metric").value
    usable = ~p.high_uncertainty_pair & ~p.unstable_replicate
    rows[a] = {
        "pairs_scored": len(p),
        "high_uncertainty_pairs": int(p.high_uncertainty_pair.sum()),
        "unstable_replicate_pairs": int(p.unstable_replicate.sum()),
        "usable_pairs": int(usable.sum()),
        "median_mean_n_reads": p.mean_n_reads.median(),
        "median_replicate_sd_log2 (n_rep>1, usable)": p.loc[usable & (p.n_replicates > 1), "std_log2_ratio_a_over_b"].median(),
        "strains_in_hierarchy": h["n_strains_in_main_component"],
        "bt_pseudo_r2": h["bt_pseudo_r2"], "dci": h["dci"],
        "frac_intransitive_triads": h["frac_intransitive_triads"],
    }
c05 = pd.DataFrame(rows)
c05.to_csv(OUT / "c05_relative_abundance.csv")
print(c05.round(4).to_string())

# strength agreement between arms
bt = {a: pd.read_csv(od("05_engineer_relative_abundances/relative_abundance", a) / "r05_bt_strengths.csv") for a in ARMS}
print("\nBT strength columns:", bt["original"].columns.tolist())
key = bt["original"].columns[0]
val = [c for c in bt["original"].columns if "strength" in c.lower()][0]
b = bt["original"].merge(bt["corrected"], on=key, suffixes=("_o", "_c"))
print(f"BT strengths, same strain under both demuxes: n={len(b)}, "
      f"spearman {spearmanr(b[val + '_o'], b[val + '_c'])[0]:.3f}, pearson {pearsonr(b[val + '_o'], b[val + '_c'])[0]:.3f}")

# --- c06 genomic ML ---------------------------------------------------------------------
header("c06 genomic ML (who wins)")
g = {a: od("06_ml_analysis/genomic_ml", a) for a in ARMS}
ds = pd.concat({a: pd.read_csv(g[a] / "g01_dataset_summary.csv").set_index("step").n_pairs for a in ARMS}, axis=1)
print(ds.to_string())
ceil = pd.concat({a: pd.read_csv(g[a] / "g02_label_noise_ceiling.csv").set_index("estimator")[["n", "ceiling_r2"]] for a in ARMS}, axis=1)
print("\nlabel-noise ceiling\n" + ceil.round(3).to_string())
cv = pd.concat({a: pd.read_csv(g[a] / "g03_cv_summary.csv") for a in ARMS}, names=["arm"]).reset_index(level=0)
keep = ["regime", "model", "n_test", "r2_mean", "r2_sd", "spearman_rho_mean", "spearman_rho_sd", "sign_accuracy_mean"]
c06 = cv[["arm"] + keep].pivot_table(index=["regime", "model"], columns="arm",
                                    values=["r2_mean", "r2_sd", "spearman_rho_mean", "sign_accuracy_mean", "n_test"])
c06 = c06.reindex(columns=pd.MultiIndex.from_product([["n_test", "r2_mean", "r2_sd", "spearman_rho_mean", "sign_accuracy_mean"], ARMS]))
c06.to_csv(OUT / "c06_cv_summary.csv")
print("\n" + c06.round(3).to_string())
g4 = pd.concat({a: pd.read_csv(g[a] / "g04_strength_from_genome_metrics.csv").set_index("feature_source")[["n", "r2", "spearman_rho"]] for a in ARMS}, axis=1)
g4.to_csv(OUT / "c06_strength_from_genome.csv")
print("\ngenome -> per-strain strength\n" + g4.round(3).to_string())
ds.to_csv(OUT / "c06_dataset.csv"); ceil.to_csv(OUT / "c06_noise_ceiling.csv")

# --- c07 plate tier + depth curve -------------------------------------------------------
header("c07 plate-reader tier, read-depth curve")
for a in ARMS:
    f = od("06_ml_analysis/genomic_ml_plate", a) / "p04_paired_vs_baseline.csv"
    if f.exists():
        print(f"--- {a}: p04_paired_vs_baseline"); print(pd.read_csv(f).round(3).to_string(index=False))
for a in ARMS:
    f = g[a] / "g06_read_depth_curve.csv"
    if f.exists():
        print(f"--- {a}: g06_read_depth_curve"); print(pd.read_csv(f).round(3).to_string(index=False))
print(f"\nwrote {OUT}")
