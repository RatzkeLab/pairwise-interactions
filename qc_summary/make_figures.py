"""Figures for the QC summary documents. Reads only qc_summary/outputs/ + existing pipeline CSVs."""
import sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
EXPS = HERE.parent
sys.path.insert(0, str(EXPS / "shared_pipelines"))
from genomic_ml import COLOR_BLUE, COLOR_RED, COLOR_GRID, COLOR_TEXT_SECONDARY, COLOR_GOOD, COLOR_CRITICAL
FIG = HERE / "outputs" / "figures"; FIG.mkdir(parents=True, exist_ok=True)
C = {"20260630": COLOR_BLUE, "20260721": COLOR_RED}

# experiment ids are numeric-looking and round-trip from CSV as int64; force str so
# every filter and colour lookup in this file compares like with like
g = pd.read_csv(HERE / "outputs" / "qc02b_two_definitions.csv")
g["experiment"] = g["experiment"].astype(str)

# --- f1: label support, the headline
refs = ["corroborated_db(185)", "corroborated_db_min5(87)", "genome_16S"]
fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
for ax, ref in zip(axes, refs):
    d = g[g.reference == ref]
    labels, vals, cols, ns = [], [], [], []
    for e in ("20260630", "20260721"):
        for c in ("pair", "mono", "merged"):
            r = d[(d.experiment == e) & (d.consensus == c)]
            if not len(r): continue
            labels.append(f"{e}\n{c}"); vals.append(r.pct_consistent.iloc[0])
            cols.append(C[e]); ns.append(int(r.n_testable.iloc[0]))
    b = ax.bar(labels, vals, color=cols)
    ax.bar_label(b, labels=[f"n={n}" for n in ns], fontsize=7.5, padding=2)
    ax.set_ylim(0, 108); ax.set_title(ref, fontsize=10)
    ax.tick_params(labelsize=8)
axes[0].set_ylabel("% of testable labels supported\n(consensus matches its own label at ≥0.99)")
fig.suptitle("Do the 16S consensus sequences support the well labels?  "
             "20260630 (blue) vs 20260721 (red)", fontweight="bold")
fig.tight_layout(); fig.savefig(FIG / "qc_f1_label_support.png", dpi=160); plt.close(fig)

# --- f2: read depth
fig, ax = plt.subplots(figsize=(7.5, 5))
for e in ("20260630", "20260721"):
    d = pd.read_csv(EXPS / e / "analysis/mapping_validation/outputs/01_samples_gt5reads.csv")
    ax.hist(d.n_reads, bins=np.arange(0, 165, 5), histtype="step", lw=2, color=C[e],
            label=f"{e}  n={len(d)}, mean {d.n_reads.mean():.1f}, median {d.n_reads.median():.0f}")
ax.axvline(20, color=COLOR_TEXT_SECONDARY, ls="--", lw=1)
ax.text(21, ax.get_ylim()[1]*.93, "below ~20 reads the\nrelative-abundance label\ndegrades sharply",
        fontsize=8, color=COLOR_TEXT_SECONDARY)
ax.set_xlabel("reads per sample (wells with >5 reads)"); ax.set_ylabel("wells")
ax.legend(frameon=False, fontsize=8.5)
ax.set_title("Sequencing depth is essentially identical between the two runs", fontsize=11, fontweight="bold")
fig.tight_layout(); fig.savefig(FIG / "qc_f2_read_depth.png", dpi=160); plt.close(fig)

# --- f3: the two definitions, and cross-experiment agreement
d = pd.read_csv(HERE / "outputs" / "qc01_consensus_vs_reference_detail.csv")
d["experiment"] = d["experiment"].astype(str)
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5))
ax = axes[0]
sub = g[(g.reference == "corroborated_db(185)")]
x = np.arange(len(sub)); w = 0.38
ax.bar(x - w/2, sub.pct_consistent, w, color=COLOR_GOOD, label="label agrees (≥0.99)")
ax.bar(x + w/2, sub.pct_is_best, w, color=COLOR_TEXT_SECONDARY, label="label is the CLOSEST match")
ax.set_xticks(x); ax.set_xticklabels([f"{r.experiment}\n{r.consensus}" for r in sub.itertuples()], fontsize=8)
ax.set_ylabel("% of testable labels"); ax.legend(frameon=False, fontsize=8.5)
ax.set_title("Two definitions of 'supported'\nthe gap is 16S near-twins, not wrong labels", fontsize=10)
ax = axes[1]
x = d[(d.reference.str.startswith("other_experiment")) & d.self_identity.notna()
      & (d.consensus == "pair")]
for e, s in x.groupby("experiment"):
    ax.hist(s.self_identity, bins=np.linspace(0.80, 1.0, 25), histtype="step", lw=2,
            color=C[e], label=f"{e} (n={len(s)})")
ax.axvline(0.99, color=COLOR_CRITICAL, ls="--", lw=1.2)
ax.text(0.9885, ax.get_ylim()[1]*.6, "same-organism\nthreshold", fontsize=8,
        color=COLOR_CRITICAL, ha="right")
ax.set_xlabel("identity between the two experiments' consensus for the SAME well label")
ax.set_ylabel("labels"); ax.legend(frameon=False, fontsize=8.5)
ax.set_title("The 18 labels both experiments share\nonly 3 are the same organism", fontsize=10)
fig.tight_layout(); fig.savefig(FIG / "qc_f3_definitions_and_crossexp.png", dpi=160); plt.close(fig)
print("wrote:", *[p.name for p in sorted(FIG.glob('*.png'))], sep="\n  ")
