"""Per-well QC: does run b (after the library top-up) describe the same wells as run a?

Needs, for both runs, the corrected demux and the qc+ra steps of run_all.py:
    python ../run_all.py --arm corrected         (pass_a; already done)
    python ../run_all.py --arm pass_b --steps qc,ra

If the top-up was the same pooled library, each well's SHARE of the reads and each pair
well's STRAIN RATIO should agree between runs up to counting noise. Both are tested against
that noise explicitly rather than eyeballed:
    d01  demux: assignment rate, minibar match classes, false rate on unsequenced plates
    d02  well shares: over-dispersion of run-b counts given run-a shares (1 = Poisson only)
    d03  strain ratios: z-score of the a-vs-b difference per pair well under a binomial null
    d04  mapping QC: off-target fraction, qc_status, mono wells' dominant strain

    python demux_ab_qc.py     (env karl_seq_analysis)
"""
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, chi2

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import config  # noqa: E402

OUT = HERE / "outputs"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
ARMS = {"a": "corrected", "b": "pass_b"}          # run -> config arm
COLORS = {"a": "#3b75af", "b": "#ef8636"}


def header(t):
    print(f"\n{'=' * 90}\n{t}\n{'=' * 90}")


def run_dir(run):
    return config.DEMUX_DIRS[ARMS[run]].parent


def od(step, run):
    return config.out_dir(step, ARMS[run])


# --- d01 demux ------------------------------------------------------------------------------
header("d01 demultiplexing")
dem = {}
rows = []
for r in ARMS:
    s = run_dir(r) / "summary"
    filt = pd.read_csv(s / "filtered_read_counts.tsv", sep="\t").iloc[:, 1].sum()
    d = pd.read_csv(s / "demux_summary.tsv", sep="\t")
    d = d[d["sample"].str.match(r"Plate\d\d_")].copy()
    d["plate"] = d["sample"].str[5:7].astype(int)
    d["sequenced"] = d.plate.isin(config.SEQUENCED_PLATES)
    dem[r] = d.set_index("sample")
    n_seq, n_uns = d.loc[d.sequenced, "reads"].sum(), d.loc[~d.sequenced, "reads"].sum()
    false_on_seq = n_uns / (~d.sequenced).sum() * d.sequenced.sum()
    # minibar's match class from each read header: HH = barcode + primer at both ends
    cls = {}
    for f in sorted(config.DEMUX_DIRS[ARMS[r]].glob("Plate*.fastq")):
        if int(f.name[5:7]) not in config.SEQUENCED_PLATES or f.stat().st_size == 0:
            continue
        with open(f) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    m = re.search(r"\|([HhXx?])[-+?x]\([^)]*\),([HhXx?])", line)
                    k = (m.group(1) + m.group(2)) if m else "other"
                    cls[k] = cls.get(k, 0) + 1
    tot = sum(cls.values())
    rows.append({"run": r, "filtered_reads": filt, "assigned": n_seq + n_uns,
                 "pct_of_filtered": 100 * (n_seq + n_uns) / filt,
                 "on_sequenced_plates": n_seq, "on_unsequenced_plates": n_uns,
                 "est_false_rate_pct": 100 * false_on_seq / n_seq,
                 "pct_HH": 100 * cls.get("HH", 0) / tot, "pct_Hh": 100 * cls.get("Hh", 0) / tot,
                 "pct_hh": 100 * cls.get("hh", 0) / tot,
                 "median_reads_per_sequenced_well": d.loc[d.sequenced, "reads"].median()})
d01 = pd.DataFrame(rows).set_index("run")
d01.to_csv(OUT / "d01_demux.csv")
print(d01.round(3).T.to_string())

# --- d02 well shares -------------------------------------------------------------------------
header("d02 do wells get the same SHARE of reads in both runs?")
w = dem["a"][["reads", "plate", "sequenced"]].join(dem["b"][["reads"]], lsuffix="_a", rsuffix="_b")
w = w[w.sequenced]
Na, Nb = w.reads_a.sum(), w.reads_b.sum()
# Two-sample (wells x runs) Pearson chi-square: if b is the same library, each well has one
# share, estimated from both runs. Residuals from BOTH runs count -- summing run b's alone
# deflates the statistic by ~Na/(Na+Nb), which an earlier version of this script did.
share = (w.reads_a + w.reads_b) / (Na + Nb)
exp_a, exp_b = share * Na, share * Nb
ok = share > 0
pearson_chi2 = ((((w.reads_a - exp_a) ** 2) / exp_a) + (((w.reads_b - exp_b) ** 2) / exp_b))[ok].sum()
dof = ok.sum() - 1
disp = pearson_chi2 / dof
lw = np.log2((w.reads_b + 1) / (w.reads_a + 1) * Na / Nb)
by_plate = w.groupby("plate")[["reads_a", "reads_b"]].sum()
by_plate["b_over_a_share"] = (by_plate.reads_b / Nb) / (by_plate.reads_a / Na)
d02 = pd.DataFrame([{
    "wells": len(w), "reads_a": Na, "reads_b": Nb,
    "spearman_well_counts": spearmanr(w.reads_a, w.reads_b)[0],
    "dispersion_vs_poisson": disp,
    "wells_share_changed_2x": int(((lw.abs() > 1) & (w.reads_a + w.reads_b >= 40)).sum()),
    "wells_zero_in_b_but_ge20_in_a": int(((w.reads_b == 0) & (w.reads_a >= 20)).sum()),
    "wells_zero_in_a_but_ge20_in_b": int(((w.reads_a == 0) & (w.reads_b >= 20)).sum()),
}])
d02.to_csv(OUT / "d02_well_shares.csv", index=False)
by_plate.to_csv(OUT / "d02_plate_shares.csv")
print(d02.round(3).T.to_string())
print("\nper-plate share ratio (b/a; 1 = same share)\n" + by_plate.b_over_a_share.round(3).to_string())
print("\nA dispersion near 1 means well-to-well differences are counting noise. Values of a few")
print("are normal for independent library loadings (pore capture is not perfectly random);")
print(">>10, or plate-level ratios far from 1, would mean the pool composition changed.")

# --- d03 strain ratios --------------------------------------------------------------------------
header("d03 do pair wells show the same STRAIN RATIO in both runs?")
ra = {r: pd.read_csv(od("05_engineer_relative_abundances/relative_abundance", r) /
                     "r02_well_interaction_scores.csv").set_index("sample_id") for r in ARMS}
m = ra["a"].join(ra["b"], lsuffix="_a", rsuffix="_b", how="inner")
m = m[~m.missing_reference_a.astype(bool) & (m.ref_pair_bp_dist_a >= 10)]
na = m.n_strain1_a + m.n_strain2_a
nb = m.n_strain1_b + m.n_strain2_b
m = m[(na >= 10) & (nb >= 10)]
na, nb = na[m.index], nb[m.index]
pa, pb = m.n_strain1_a / na, m.n_strain1_b / nb
pp = (m.n_strain1_a + m.n_strain1_b) / (na + nb)
se = np.sqrt(pp * (1 - pp) * (1 / na + 1 / nb))
z = ((pa - pb) / se).where(se > 0)
zz = z.dropna()
d03 = pd.DataFrame([{
    "pair_wells": len(m),
    "pearson_r": pearsonr(pa, pb)[0],
    "median_abs_diff": (pa - pb).abs().median(),
    "mean_diff_a_minus_b": (pa - pb).mean(),
    "z_sd (1 = binomial noise only)": zz.std(),
    "pct_abs_z_gt3 (0.27% expected)": 100 * (zz.abs() > 3).mean(),
    "winner_flips_pct": 100 * ((pa - .5) * (pb - .5) < 0).mean(),
}])
d03.to_csv(OUT / "d03_strain_ratio_agreement.csv", index=False)
m.assign(p_a=pa, p_b=pb, z=z)[["strain1_a", "strain2_a", "p_a", "p_b", "z"]].to_csv(
    OUT / "d03_per_well.csv")
print(d03.round(4).T.to_string())

fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
ax[0].scatter(w.reads_a + 1, w.reads_b + 1, s=4, alpha=.3)
x = np.array([1, max(w.reads_a.max(), 2)])
ax[0].plot(x, x * Nb / Na, "k--", lw=.8, label=f"same share (x{Nb / Na:.2f})")
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].legend()
ax[0].set_xlabel("run a reads + 1"); ax[0].set_ylabel("run b reads + 1"); ax[0].set_title("Reads per well")
ax[1].scatter(pa, pb, s=5, alpha=.35)
ax[1].plot([0, 1], [0, 1], "k--", lw=.8)
ax[1].set_xlabel("strain1 fraction, run a"); ax[1].set_ylabel("strain1 fraction, run b")
ax[1].set_title(f"Pair wells, refs >=10 bp apart (n={len(m)})")
bins = np.linspace(-6, 6, 61)
ax[2].hist(zz.clip(-6, 6), bins=bins, density=True, alpha=.6, label="observed")
g = np.linspace(-6, 6, 200)
ax[2].plot(g, np.exp(-g ** 2 / 2) / np.sqrt(2 * np.pi), "k", lw=1, label="binomial noise only")
ax[2].set_xlabel("z = (p_a - p_b) / SE"); ax[2].legend(); ax[2].set_title("Strain-ratio differences vs counting noise")
fig.tight_layout(); fig.savefig(FIG / "d02_d03_well_agreement.png", dpi=140); plt.close(fig)

# --- d04 mapping QC ------------------------------------------------------------------------------
header("d04 mapping QC")
qc = {r: pd.read_csv(od("04_qc/mapping_validation", r) / "05_combined_sample_summary.csv") for r in ARMS}
d04 = pd.concat({r: qc[r].groupby(["well_type", "qc_status"]).size() for r in ARMS}, axis=1).fillna(0).astype(int)
d04.to_csv(OUT / "d04_qc_status.csv")
print(d04.to_string())
off = pd.DataFrame({r: ra[r].off_target_frac.describe(percentiles=[.5, .9, .99]) for r in ARMS})
print("\noff-target read fraction per pair well\n" + off.round(4).to_string())
off.to_csv(OUT / "d04_off_target.csv")
q = qc["a"].set_index("sample_id").join(qc["b"].set_index("sample_id"), lsuffix="_a", rsuffix="_b", how="inner")
col = [c for c in qc["a"].columns if c.startswith("mm2_") and c.endswith("_top_other_strain")][0]
mono = q[q.well_type_a == "mono"]
same_status = (q.qc_status_a == q.qc_status_b).mean()
print(f"\nwells scored in both runs: {len(q)}; same qc_status: {100 * same_status:.1f}%")
print(f"mono wells: {len(mono)}; same qc_status {100 * (mono.qc_status_a == mono.qc_status_b).mean():.0f}%")
print(f"\nwrote {OUT}")
