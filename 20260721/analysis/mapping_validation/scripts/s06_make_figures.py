"""06 -- Summary figures for the read-mapping validation.

Color usage follows a fixed, semantic assignment throughout (never re-cycled per plot):
  - well_type   : pair = blue #2a78d6, mono = orange #eb6834   (categorical slots 1 & 2)
  - qc_status   : confirmed = good #0ca30c, single/low-confidence = warning #fab219,
                  neither confirmed = critical #d03b3b          (status palette; semantic, not decorative)
  - magnitude   : single-hue blue ramp (matplotlib 'Blues')     (sequential)

Figures use PRIMARY_DB (merged_consensus_mono_priority) throughout, per the identity
cross-check in s02 / s05: corroborated_db's strain-name labels are only reliable for a
minority of strains in this experiment, so it is not used as the primary QC reference
(figure 00 visualizes that finding directly).
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import OUT_DIR, FIG_DIR
from s02_prepare_references import CROSS_CHECK_RELIABLE_THRESHOLD
from s05_combine_and_flag import PRIMARY_DB

# ---- fixed semantic palette -------------------------------------------------
COLOR_PAIR = "#2a78d6"
COLOR_MONO = "#eb6834"
COLOR_GOOD = "#0ca30c"
COLOR_WARNING = "#fab219"
COLOR_CRITICAL = "#d03b3b"
COLOR_GRID = "#d8d7d2"
COLOR_TEXT_SECONDARY = "#52514e"

QC_COLOR = {
    "mono_confirmed": COLOR_GOOD,
    "pair_both_confirmed": COLOR_GOOD,
    "mono_low_confidence": COLOR_WARNING,
    "pair_single_dominant": COLOR_WARNING,
    "pair_neither_confirmed": COLOR_CRITICAL,
}

plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": COLOR_GRID,
        "axes.grid": True,
        "grid.color": COLOR_GRID,
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
    }
)


def fig00_reference_reliability(cross_check):
    cc = cross_check.sort_values("norm_dist_corrdb_vs_own_consensus").reset_index(drop=True)
    colors = [COLOR_GOOD if r else COLOR_CRITICAL for r in cc["reliable"]]

    fig, ax = plt.subplots(figsize=(11, 4.6))
    ax.bar(range(len(cc)), cc["norm_dist_corrdb_vs_own_consensus"], color=colors, width=0.8)
    ax.axhline(
        CROSS_CHECK_RELIABLE_THRESHOLD, color=COLOR_TEXT_SECONDARY, lw=1, ls="--",
        label=f"reliable cutoff = {CROSS_CHECK_RELIABLE_THRESHOLD}",
    )
    ax.set_xticks(range(len(cc)))
    ax.set_xticklabels(cc["strain"], rotation=90, fontsize=6)
    ax.set_ylabel("norm. edit distance:\ncorroborated_db vs. own consensus")
    n_rel = int(cc["reliable"].sum())
    ax.set_title(
        f"corroborated_db name == same organism as this experiment's? {n_rel}/{len(cc)} agree (green)",
        pad=10,
    )
    fig.text(
        0.5, 0.955,
        "short plate-well-style names (e.g. \"D13\") are reused across unrelated experiments -- most are name collisions",
        ha="center", fontsize=8.5, color=COLOR_TEXT_SECONDARY,
    )
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(FIG_DIR / "00_reference_reliability_check.png", dpi=150)
    plt.close(fig)


def fig01_match_quality(reads_primary, summary):
    thr = float(summary[f"edlib_{PRIMARY_DB}_threshold"].iloc[0])
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    ax.hist(
        reads_primary["best_dist"].dropna(),
        bins=np.linspace(0, 0.4, 81),
        color=COLOR_PAIR,
        edgecolor="white",
        linewidth=0.3,
    )
    ax.axvline(thr, color=COLOR_CRITICAL, ls="--", lw=1.5, label=f"confident-match threshold = {thr}")
    ax.set_yscale("log")
    ax.set_xlabel("normalized edit distance: read → its assigned expected-strain reference")
    ax.set_ylabel("# reads (log scale)")
    ax.set_title(f"Constrained mapping quality ({PRIMARY_DB})", pad=10)
    fig.text(
        0.5, 0.925,
        "left mode = confirms expected strain, right mode = poor match to either expected strain",
        ha="center", fontsize=8.5, color=COLOR_TEXT_SECONDARY,
    )
    ax.legend(frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(FIG_DIR / "01_match_quality_histogram.png", dpi=150)
    plt.close(fig)


def fig02_qc_status_by_well_type(summary):
    order_mono = ["mono_confirmed", "mono_low_confidence"]
    order_pair = ["pair_both_confirmed", "pair_single_dominant", "pair_neither_confirmed"]

    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=False)
    for ax, order, wt, title in zip(
        axes, [order_mono, order_pair], ["mono", "pair"], ["mono wells", "pair wells"]
    ):
        counts = summary.loc[summary.well_type == wt, "qc_status"].value_counts().reindex(order, fill_value=0)
        colors = [QC_COLOR[s] for s in order]
        bars = ax.bar(range(len(order)), counts.values, color=colors, width=0.6)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([s.replace("_", "\n") for s in order], fontsize=8)
        ax.set_title(f"{title} (n={int(counts.sum())})")
        ax.set_ylabel("# wells")
        for b, v in zip(bars, counts.values):
            ax.text(b.get_x() + b.get_width() / 2, v, f"{v}", ha="center", va="bottom", fontsize=9)
    fig.suptitle(
        f"QC call per well ({PRIMARY_DB}, constrained mapping)\n"
        "good = confirmed, warning = only one strain / low confidence, critical = neither expected strain confirmed",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(FIG_DIR / "02_qc_status_by_well_type.png", dpi=150)
    plt.close(fig)


def fig03_purity_histogram(summary):
    summary = summary.copy()
    summary["dominant_frac"] = summary[
        [f"edlib_{PRIMARY_DB}_frac_strain1", f"edlib_{PRIMARY_DB}_frac_strain2"]
    ].max(axis=1)

    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.linspace(0, 1, 41)
    for wt, color, label in [("pair", COLOR_PAIR, "pair"), ("mono", COLOR_MONO, "mono")]:
        vals = summary.loc[summary.well_type == wt, "dominant_frac"].dropna()
        ax.hist(vals, bins=bins, color=color, alpha=0.75, label=f"{label} (n={len(vals)})", edgecolor="white", linewidth=0.3)
    ax.set_xlabel("fraction of confident reads matching the dominant expected strain")
    ax.set_ylabel("# wells")
    ax.set_title(f"Read purity vs. designed well type ({PRIMARY_DB})")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "03_purity_histogram.png", dpi=150)
    plt.close(fig)


def fig04_depth_vs_confidence(summary):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for wt, color, label in [("pair", COLOR_PAIR, "pair"), ("mono", COLOR_MONO, "mono")]:
        sub = summary[summary.well_type == wt]
        ax.scatter(
            sub["n_reads"], sub[f"mm2_{PRIMARY_DB}_frac_expected"],
            s=14, color=color, alpha=0.5, edgecolor="none", label=f"{label} (n={len(sub)})",
        )
    ax.set_xscale("log")
    ax.set_xlabel("reads in well (log scale)")
    ax.set_ylabel("fraction of reads whose best genome-wide hit\nis an expected strain (minimap2)")
    ax.set_title(f"Read depth vs. mapping confidence ({PRIMARY_DB})")
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "04_depth_vs_confidence.png", dpi=150)
    plt.close(fig)


def fig05_cross_db_agreement(summary):
    """Among strains flagged reliable in corroborated_db, does the independent db agree
    with the self-derived primary db on which strain dominates each well?"""
    summary = summary.copy()
    summary["dom_corr"] = summary[
        ["edlib_corroborated_db_frac_strain1", "edlib_corroborated_db_frac_strain2"]
    ].max(axis=1)
    summary["dom_primary"] = summary[
        [f"edlib_{PRIMARY_DB}_frac_strain1", f"edlib_{PRIMARY_DB}_frac_strain2"]
    ].max(axis=1)
    s1_reliable = summary["corrdb_strain1_reliable"] == True  # noqa: E712 (NaN-safe vs fillna+bool downcast)
    s2_reliable = summary["corrdb_strain2_reliable"] == True  # noqa: E712
    both_reliable = s1_reliable & ((summary["well_type"] == "mono") | s2_reliable)
    sub = summary[both_reliable].dropna(subset=["dom_corr", "dom_primary"])

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    for wt, color, label in [("pair", COLOR_PAIR, "pair"), ("mono", COLOR_MONO, "mono")]:
        s = sub[sub.well_type == wt]
        ax.scatter(s["dom_corr"], s["dom_primary"], s=14, color=color, alpha=0.6, edgecolor="none", label=f"{label} (n={len(s)})")
    ax.plot([0, 1], [0, 1], color=COLOR_TEXT_SECONDARY, lw=1, ls="--")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("dominant-strain read fraction (corroborated_db)")
    ax.set_ylabel(f"dominant-strain read fraction ({PRIMARY_DB})")
    ax.set_title("Independent vs. self-derived reference", pad=10)
    fig.text(
        0.55, 0.925,
        f"{sub['sample_id'].nunique()} wells where both expected strains are name-reliable in corroborated_db",
        ha="center", fontsize=8.5, color=COLOR_TEXT_SECONDARY,
    )
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(FIG_DIR / "05_cross_db_agreement.png", dpi=150)
    plt.close(fig)


def fig06_contamination_offenders(contam):
    # reference near-twins (e.g. K6/P17, 1.5% apart) are excluded: 16S genuinely cannot
    # tell those apart, so they aren't evidence of real contamination -- see figure 00-style
    # reasoning, applied to this experiment's own reference set instead of corroborated_db
    genuine = contam[~contam["likely_reference_twin"]]
    if len(genuine) == 0:
        return
    top = genuine[f"mm2_{PRIMARY_DB}_top_other_strain"].value_counts().head(15).sort_values()
    fig, ax = plt.subplots(figsize=(7, max(3, 0.35 * len(top))))
    cmap = plt.get_cmap("Blues")
    norm_vals = top.values / top.values.max()
    colors = [cmap(0.35 + 0.55 * v) for v in norm_vals]
    ax.barh(range(len(top)), top.values, color=colors)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top.index, fontsize=9)
    ax.set_xlabel("# wells where this strain is the top unexpected best-hit")
    ax.set_title("Most frequent unexpected ('contaminant') strains", pad=10)
    fig.text(
        0.5, 0.94,
        f"minimap2 best hit vs. {PRIMARY_DB}; excludes {len(contam) - len(genuine)}/{len(contam)} "
        "reference near-twin false positives",
        ha="center", fontsize=8, color=COLOR_TEXT_SECONDARY,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(FIG_DIR / "06_contamination_offenders.png", dpi=150)
    plt.close(fig)


def make_all_figures():
    summary = pd.read_csv(OUT_DIR / "05_combined_sample_summary.csv")
    contam = pd.read_csv(OUT_DIR / "05_contamination_candidates.csv")
    reads_primary = pd.read_csv(OUT_DIR / f"03_edlib_read_assignments_{PRIMARY_DB}.csv.gz")
    cross_check = pd.read_csv(OUT_DIR / "02_reference_cross_check.csv")

    fig00_reference_reliability(cross_check)
    fig01_match_quality(reads_primary, summary)
    fig02_qc_status_by_well_type(summary)
    fig03_purity_histogram(summary)
    fig04_depth_vs_confidence(summary)
    fig05_cross_db_agreement(summary)
    fig06_contamination_offenders(contam)

    print(f"saved 7 figures -> {FIG_DIR}")


if __name__ == "__main__":
    make_all_figures()
