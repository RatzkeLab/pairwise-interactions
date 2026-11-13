"""Figures for steps 1 and 2. Dumbbells: one row per model, one dot per arm, all arms of a
panel scored on the SAME test pairs fold by fold. Dot = Spearman rho, mean over 100 folds
(20 repeats x 5); bar = SD of the per-repeat mean across repeats.

  f01_leakage.png         published folds vs the same folds with the leak purged from training
                          (and the size-matched control: same data loss, leak kept)
  f02_wrong_genomes.png   leak-free folds; trained on all wells vs without wrong-genome wells vs
                          the same amount removed at random; scored on correct-genome pairs
"""
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import config
from s05_tables import cell, load, sel

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"      # categorical slots 1-3 (validated set)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d8d7d2"
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white", "axes.edgecolor": GRID,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "font.size": 9.5, "axes.titlesize": 10.5, "axes.titleweight": "bold",
    "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
})
FIG = config.OUT / "figures"
FIG.mkdir(exist_ok=True)
for old in FIG.glob("*.png"):
    old.unlink()

ROWS = [("xgboost_raw_ko", "XGBoost, raw KO (best by ρ)"),
        ("two_stage_ridge", "two-stage ridge, KO"),
        ("ridge_pca", "ridge on KO PCs"),
        ("two_stage_ridge_phylo16s_only", "16S only, two-stage"),
        ("ridge_phylo16s_only", "16S only, ridge"),
        ("clone_lookup_no_genomics", "clone lookup, no genomics")]
EXPS = ["20260907", "20260630"]
NOTE = {"20260630": "16S-only rows unstable in 20260630 (PCA whitening; see README)"}


def panel(ax, exp, series, rows):
    y = np.arange(len(rows))[::-1]
    for (model, _), yy in zip(rows, y):
        pts = []
        for i, (frame, _, color) in enumerate(series):
            f = frame(exp, model)
            if not len(f):
                continue
            c = cell(f)
            off = (i - (len(series) - 1) / 2) * 0.2
            ax.errorbar(c["rho"], yy + off, xerr=c["rho_sd"], fmt="o", ms=6.5, color=color, ecolor=color,
                        elinewidth=1.2, markeredgecolor="white", markeredgewidth=1.0, zorder=3)
            pts.append((c["rho"], yy + off))
        if len(pts) > 1:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color=GRID, lw=1.6, zorder=2)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([lab for _, lab in rows])
    ax.set_xlim(-0.15, 0.95)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Spearman ρ, predicted vs observed log2(a/b)")
    ax.set_title(exp, loc="left")
    if exp in NOTE:   # beside the title, clear of the data
        ax.text(1.0, 1.012, NOTE[exp], transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color=INK2)


def figure(series, rows, title, caption, fname):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 0.55 * len(rows) + 2.4), sharey=True)
    for ax, exp in zip(axes, EXPS):
        panel(ax, exp, series, rows)
    axes[1].tick_params(labelleft=False)
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, ms=7, label=lab) for _, lab, c in series]
    fig.suptitle(title, fontweight="bold", x=0.01, ha="left", y=0.99)
    fig.legend(handles=handles, loc="upper left", ncol=len(series), frameon=False, bbox_to_anchor=(0.005, 0.955))
    caption = "\n".join(textwrap.fill(c, 190) for c in caption.split("\n"))
    fig.text(0.01, 0.01, caption, fontsize=8, color=INK2, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.04 + 0.022 * caption.count("\n"), 1, 0.9))
    fig.savefig(FIG / fname, dpi=170)
    plt.close(fig)


def main():
    s02, s03 = load()
    o = sel(s03, folds="by_organism")
    figure([(lambda e, m: sel(s02, exp=e, model=m, arm="by_well"), "published folds (by well)", BLUE),
            (lambda e, m: sel(s02, exp=e, model=m, arm="by_well_random"), "same data loss, leak kept", AQUA),
            (lambda e, m: sel(s02, exp=e, model=m, arm="by_well_purged"), "leak removed", ORANGE)],
           ROWS, "Step 1 — how much of the published score was leakage?",
           "Same test pairs in every arm. 'Leak removed': training pairs involving any organism present in the "
           "test fold are dropped (~1/3 of training). 'Same data loss, leak kept': the same number of training "
           "pairs dropped at random among other organisms.\nLeak = green minus orange. Dot = mean over 100 folds; "
           "bar = SD across 20 repeats.", "f01_leakage.png")
    figure([(lambda e, m: sel(o, exp=e, model=m, arm="all", test_subset="clean"), "trained on all wells", BLUE),
            (lambda e, m: sel(o, exp=e, model=m, arm="drop_wrong_random", test_subset="clean"),
             "same amount removed at random", AQUA),
            (lambda e, m: sel(o, exp=e, model=m, arm="drop_wrong", test_subset="clean"),
             "wrong-genome wells removed", ORANGE)],
           ROWS[:-1], "Step 2 — removing wells whose genome is not the organism in the well  (leak-free CV)",
           "Leak-free folds (organisms held out together). Every arm is scored on the SAME test pairs: pairs whose "
           "two wells both carry a genome that is not known to be wrong.\nDot = mean over 100 folds; bar = SD across "
           "20 repeats.", "f02_wrong_genomes.png")
    print("wrote", sorted(p.name for p in FIG.glob("*.png")))


if __name__ == "__main__":
    main()
