"""06 -- Summary figures for the relative-abundance / interaction analysis.

Color usage follows a fixed, semantic assignment throughout (never re-cycled per plot):
  - magnitude         : single-hue blue ramp (matplotlib 'Blues')        (sequential)
  - dominance/balance : blue <-> red diverging, gray midpoint at 0.5 / 0  (diverging)
  - categorical splits: blue #2a78d6 / red #e34948                       (categorical slots 1 & 8)
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

from ra_common import OUT_DIR, FIG_DIR

COLOR_BLUE = "#2a78d6"
COLOR_RED = "#e34948"
COLOR_GRID = "#d8d7d2"
COLOR_TEXT_SECONDARY = "#52514e"
COLOR_GOOD = "#0ca30c"
COLOR_CRITICAL = "#d03b3b"

DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "blue_gray_red", [COLOR_BLUE, "#f0efec", COLOR_RED], N=256
)

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


def fig00_read_classification_quality(wells):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(wells["off_target_frac"].dropna(), bins=40, color=COLOR_BLUE, edgecolor="white", linewidth=0.3)
    axes[0].set_xlabel("off-target read fraction")
    axes[0].set_ylabel("# wells")
    axes[0].set_title("Reads matching neither expected strain")

    axes[1].hist(wells["uncertainty_score"].dropna(), bins=40, color=COLOR_BLUE, edgecolor="white", linewidth=0.3)
    axes[1].set_xlabel("uncertainty score (fraction of on-target reads unassignable)")
    axes[1].set_ylabel("# wells")
    axes[1].set_title("Per-well classification uncertainty")
    fig.suptitle("Read-classification quality across all pair wells", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(FIG_DIR / "00_read_classification_quality.png", dpi=150)
    plt.close(fig)


def fig01_abundance_distribution(wells):
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(wells["relative_abundance_a"].dropna(), bins=np.linspace(0, 1, 41), color=COLOR_BLUE, edgecolor="white", linewidth=0.3)
    ax.axvline(0.5, color=COLOR_TEXT_SECONDARY, lw=1, ls="--", label="50/50")
    ax.set_xlabel("relative abundance of strain_a (0 = strain_b wins, 1 = strain_a wins)")
    ax.set_ylabel("# wells")
    ax.set_title("Interaction outcome across all pair wells")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "01_abundance_distribution.png", dpi=150)
    plt.close(fig)


def fig02_replicate_stability(pairs):
    rep = pairs[pairs["n_replicates"] > 1].copy()
    fig, ax = plt.subplots(figsize=(7, 5))
    for flag, color, label in [(True, COLOR_RED, "hard-to-call pair (high uncertainty)"), (False, COLOR_BLUE, "well-resolved pair")]:
        sub = rep[rep["high_uncertainty_pair"] == flag]
        ax.scatter(sub["mean_uncertainty_score"], sub["std_relative_abundance_a"], s=16, color=color, alpha=0.6, edgecolor="none", label=f"{label} (n={len(sub)})")
    ax.axhline(0.15, color=COLOR_TEXT_SECONDARY, lw=1, ls="--", label="unstable cutoff = 0.15")
    ax.set_xlabel("mean uncertainty score across replicates")
    ax.set_ylabel("std of relative_abundance_a across replicates")
    ax.set_title("Replicate stability vs. measurement difficulty")
    fig.text(
        0.5, 0.92,
        "hard-to-call pairs look artificially 'stable' -- forced toward 0.5 every time, not because the biology reproduces",
        ha="center", fontsize=8, color=COLOR_TEXT_SECONDARY,
    )
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(FIG_DIR / "02_replicate_stability.png", dpi=150)
    plt.close(fig)


def fig03_competitiveness_ranking(comp):
    df = comp.dropna(subset=["competitiveness_score"]).sort_values("competitiveness_score")
    vmax = df["competitiveness_score"].abs().max()
    colors = DIVERGING_CMAP((df["competitiveness_score"].values / vmax + 1) / 2)

    fig, ax = plt.subplots(figsize=(7, 14))
    ax.barh(range(len(df)), df["competitiveness_score"], color=colors)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["strain"], fontsize=6.5)
    ax.axvline(0, color=COLOR_TEXT_SECONDARY, lw=1)
    ax.set_xlabel("competitiveness score (mean log2 ratio vs. all tested opponents)")
    ax.set_title("Strain competitiveness ranking")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "03_competitiveness_ranking.png", dpi=150)
    plt.close(fig)


def fig04_mean_vs_dispersion(comp, bt):
    df = comp.merge(bt[["strain", "bt_residual_sd"]], on="strain").dropna(subset=["competitiveness_score", "bt_residual_sd"])
    fig, ax = plt.subplots(figsize=(7.5, 6))
    sc = ax.scatter(
        df["competitiveness_score"], df["bt_residual_sd"],
        s=24, c=df["frac_high_uncertainty_opponents"], cmap="Blues", edgecolor=COLOR_TEXT_SECONDARY, linewidth=0.3,
    )
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("fraction of opponents that are hard-to-call")
    ax.axvline(0, color=COLOR_TEXT_SECONDARY, lw=1, ls="--")
    ax.set_xlabel("competitiveness score (mean log2 ratio)")
    ax.set_ylabel("BT-residual dispersion (opponent-strength-corrected)")
    ax.set_title("Consistent dominance/loss vs. context-dependent outcomes")
    fig.text(
        0.5, 0.93,
        "low residual dispersion = behaves as its BT strength predicts; high = idiosyncratic, opponent-specific outcomes beyond the hierarchy",
        ha="center", fontsize=8, color=COLOR_TEXT_SECONDARY,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(FIG_DIR / "04_mean_vs_dispersion.png", dpi=150)
    plt.close(fig)


def fig05_hierarchy_heatmap(mat, bt):
    order = bt.sort_values("bt_strength", ascending=False)["strain"].tolist()
    order = [s for s in order if s in mat.index]
    m = mat.loc[order, order]

    fig, ax = plt.subplots(figsize=(11, 10))
    im = ax.imshow(m.values, cmap=DIVERGING_CMAP, vmin=0, vmax=1, aspect="equal")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=5)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=5)
    ax.grid(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("relative abundance of row strain vs. column strain")
    ax.set_title("Pairwise outcomes, ordered by Bradley-Terry strength (strongest → weakest)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "05_hierarchy_heatmap.png", dpi=150)
    plt.close(fig)


def fig06_bt_vs_naive(bt, comp):
    merged = bt.merge(comp, on="strain")
    merged = merged.dropna(subset=["bt_strength", "competitiveness_score"])
    corr = merged["bt_strength"].corr(merged["competitiveness_score"])

    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.scatter(merged["competitiveness_score"], merged["bt_strength"], s=20, color=COLOR_BLUE, alpha=0.7, edgecolor="none")
    ax.set_xlabel("naive competitiveness score (mean log2 ratio)")
    ax.set_ylabel("Bradley-Terry strength")
    ax.set_title(f"Two independent competitiveness estimates agree (r = {corr:.2f})")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "06_bt_vs_naive_competitiveness.png", dpi=150)
    plt.close(fig)


def fig07_hierarchy_consistency(summary):
    s = summary.set_index("metric")["value"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))

    consistent, upsets = s["dci_n_consistent"], s["dci_n_inconsistent_upsets"]
    axes[0].bar(["consistent", "upset"], [consistent, upsets], color=[COLOR_GOOD, COLOR_CRITICAL])
    axes[0].set_ylabel("# pairwise outcomes")
    axes[0].set_title(f"DCI = {s['dci']:.2f}")
    for i, v in enumerate([consistent, upsets]):
        axes[0].text(i, v, f"{int(v)}", ha="center", va="bottom")

    transitive = s["n_triads_evaluated"] - s["n_intransitive_triads"]
    cyclic = s["n_intransitive_triads"]
    axes[1].bar(["transitive", "cyclic"], [transitive, cyclic], color=[COLOR_GOOD, COLOR_CRITICAL])
    axes[1].set_ylabel("# strain triads")
    axes[1].set_title(f"cyclic fraction = {s['frac_intransitive_triads']:.3f}")
    for i, v in enumerate([transitive, cyclic]):
        axes[1].text(i, v, f"{int(v)}", ha="center", va="bottom")

    fig.suptitle(f"Hierarchy consistency (Bradley-Terry pseudo-R² = {s['bt_pseudo_r2']:.2f})", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(FIG_DIR / "07_hierarchy_consistency.png", dpi=150)
    plt.close(fig)


def make_all_figures():
    wells = pd.read_csv(OUT_DIR / "r02_well_interaction_scores.csv")
    wells = wells[~wells["missing_reference"]]
    pairs = pd.read_csv(OUT_DIR / "r03_pair_replicate_stats.csv")
    comp = pd.read_csv(OUT_DIR / "r04_strain_competitiveness.csv")
    bt = pd.read_csv(OUT_DIR / "r05_bt_strengths.csv")
    mat = pd.read_csv(OUT_DIR / "r05_pairwise_relative_abundance_matrix.csv", index_col=0)
    summary = pd.read_csv(OUT_DIR / "r05_hierarchy_summary.csv")

    fig00_read_classification_quality(wells)
    fig01_abundance_distribution(wells)
    fig02_replicate_stability(pairs)
    fig03_competitiveness_ranking(comp)
    fig04_mean_vs_dispersion(comp, bt)
    fig05_hierarchy_heatmap(mat, bt)
    fig06_bt_vs_naive(bt, comp)
    fig07_hierarchy_consistency(summary)

    print(f"saved 8 figures -> {FIG_DIR}")


if __name__ == "__main__":
    make_all_figures()
