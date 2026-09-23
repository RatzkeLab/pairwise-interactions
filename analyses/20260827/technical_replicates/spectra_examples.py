"""What a correlation of 0.87 (and of 0.08) actually looks like, well by well.

The summary numbers compress each well to one value at one wavelength. These figures show the
full 350-950 nm spectra of the individual wells behind them, so the aggregate statistic can be
checked against what the raw data does.

Three figures:

  s01  random technical replicates within an experiment -- the same strain pair in 2-3 different
       wells on different plates. This is the rho=0.87 / 0.69 baseline.
  s02  the same idea for the cross-experiment shared pairs -- 20260630's wells against
       20260721's wells for the same nominal pair. This is the rho=0.08 case.
  s03  best- and worst-agreeing replicates, since a random draw shows the typical case but not
       the range.

Two things are plotted per panel because they answer different questions:
  - **raw spectra** show whether two wells look like the same culture at all;
  - the annotated **plate z-scored OD600** is the quantity the correlations were computed on,
    so the panels can be tied back to the statistic. A pair can have near-identical raw spectra
    and still disagree in z if both sit near their plates' means -- z measures rank within a
    plate, not absolute similarity.

Selection is seeded, so figures are reproducible; panels are drawn from pairs with >=2 wells and
no missing spectra.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "shared_pipelines"))
from replicate_report import load_od, CONFIG, WAVELENGTHS, OUT, FIG
from genomic_ml import COLOR_BLUE, COLOR_RED, COLOR_GRID, COLOR_TEXT_SECONDARY

REP_COLORS = ["#2a78d6", "#e34948", "#0ca30c"]
OD600 = 25


def zscore_600(lay, M):
    z = np.full(len(lay), np.nan)
    col = M[:, OD600]
    for p in lay.dest_plate.unique():
        m = (lay.dest_plate == p).values
        z[m] = (col[m] - np.nanmean(col[m])) / np.nanstd(col[m])
    return z


def pair_rows(lay, M):
    ok = np.isfinite(M).all(axis=1)
    d = lay[ok].copy()
    d["_row"] = np.where(ok)[0]
    return d.groupby("pair_key")["_row"].apply(list)


def panel(ax, curves, labels, colors, title, sub):
    for c, l, col in zip(curves, labels, colors):
        ax.plot(WAVELENGTHS, c, lw=1.5, color=col, label=l)
    ax.set_title(title, fontsize=8.5, pad=2)
    # bottom-left: the legend occupies the top-right, and these collided
    ax.text(0.03, 0.06, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5,
            color=COLOR_TEXT_SECONDARY,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.75))
    ax.axvline(600, color=COLOR_GRID, lw=0.8, ls="--")
    ax.tick_params(labelsize=6.5)
    ax.legend(fontsize=6, frameon=False, loc="upper right")


def fig_within(name, n=16, seed=1):
    lay, M = load_od(name)
    z = zscore_600(lay, M)
    groups = pair_rows(lay, M)
    groups = groups[groups.map(len) >= 2]
    rng = np.random.default_rng(seed)
    keys = [groups.index[i] for i in rng.choice(len(groups), n, replace=False)]

    fig, axes = plt.subplots(4, 4, figsize=(16, 12), sharex=True)
    for ax, k in zip(axes.ravel(), keys):
        rows = groups[k][:3]
        curves = [M[r] for r in rows]
        labs = [f"P{lay.dest_plate.iloc[r]} {lay.dest_well.iloc[r]}  z={z[r]:+.2f}" for r in rows]
        zz = [z[r] for r in rows]
        panel(ax, curves, labs, REP_COLORS, "|".join(sorted(k)),
              f"Δz(OD600) = {max(zz)-min(zz):+.2f}")
    for ax in axes[-1]:
        ax.set_xlabel("wavelength (nm)")
    for ax in axes[:, 0]:
        ax.set_ylabel("OD")
    fig.suptitle(f"{name} — technical replicates: same strain pair, different plates "
                 f"(overall replicate ρ ≈ {0.869 if name=='20260630' else 0.686:.2f})",
                 fontweight="bold")
    fig.tight_layout()
    p = FIG / f"s01_within_{name}.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def fig_cross(n=16, seed=3):
    la, Ma = load_od("20260630")
    lb, Mb = load_od("20260721")
    za, zb = zscore_600(la, Ma), zscore_600(lb, Mb)
    ga, gb = pair_rows(la, Ma), pair_rows(lb, Mb)
    shared = sorted(set(ga.index) & set(gb.index), key=lambda k: sorted(k))
    rng = np.random.default_rng(seed)
    keys = [shared[i] for i in rng.choice(len(shared), min(n, len(shared)), replace=False)]

    fig, axes = plt.subplots(4, 4, figsize=(16, 12), sharex=True)
    for ax, k in zip(axes.ravel(), keys):
        ra, rb = ga[k][:2], gb[k][:2]
        curves = [Ma[r] for r in ra] + [Mb[r] for r in rb]
        labs = ([f"630  P{la.dest_plate.iloc[r]} {la.dest_well.iloc[r]} z={za[r]:+.2f}" for r in ra]
                + [f"721  P{lb.dest_plate.iloc[r]} {lb.dest_well.iloc[r]} z={zb[r]:+.2f}" for r in rb])
        cols = [COLOR_BLUE, "#7fb0e8"][:len(ra)] + [COLOR_RED, "#f09a99"][:len(rb)]
        panel(ax, curves, labs, cols, "|".join(sorted(k)),
              f"mean z: 630 {np.mean([za[r] for r in ra]):+.2f}  "
              f"721 {np.mean([zb[r] for r in rb]):+.2f}")
    for ax in axes[-1]:
        ax.set_xlabel("wavelength (nm)")
    for ax in axes[:, 0]:
        ax.set_ylabel("OD")
    fig.suptitle("Cross-experiment: the SAME nominal strain pair in each experiment "
                 "(blue = 20260630, red = 20260721; overall ρ = 0.08)", fontweight="bold")
    fig.tight_layout()
    p = FIG / "s02_cross_experiment.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def fig_extremes(name="20260630", n=6):
    """Random draws show the typical case; these show the range the statistic averages over."""
    lay, M = load_od(name)
    z = zscore_600(lay, M)
    groups = pair_rows(lay, M)
    groups = groups[groups.map(len) >= 2]
    spread = groups.map(lambda rs: abs(z[rs[0]] - z[rs[1]]))
    best = spread.nsmallest(n).index
    worst = spread.nlargest(n).index

    fig, axes = plt.subplots(2, n, figsize=(3.1 * n, 6.6), sharex=True)
    for j, k in enumerate(best):
        rows = groups[k][:2]
        panel(axes[0, j], [M[r] for r in rows],
              [f"P{lay.dest_plate.iloc[r]} {lay.dest_well.iloc[r]} z={z[r]:+.2f}" for r in rows],
              REP_COLORS, "|".join(sorted(k)), f"Δz={spread[k]:.3f}")
    for j, k in enumerate(worst):
        rows = groups[k][:2]
        panel(axes[1, j], [M[r] for r in rows],
              [f"P{lay.dest_plate.iloc[r]} {lay.dest_well.iloc[r]} z={z[r]:+.2f}" for r in rows],
              REP_COLORS, "|".join(sorted(k)), f"Δz={spread[k]:.2f}")
    axes[0, 0].set_ylabel("OD — BEST agreeing")
    axes[1, 0].set_ylabel("OD — WORST agreeing")
    for ax in axes[-1]:
        ax.set_xlabel("wavelength (nm)")
    fig.suptitle(f"{name} — the range behind the average: best and worst technical replicates",
                 fontweight="bold")
    fig.tight_layout()
    p = FIG / f"s03_extremes_{name}.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


if __name__ == "__main__":
    outs = [fig_within("20260630"), fig_within("20260721"), fig_cross(),
            fig_extremes("20260630")]
    print("wrote:", *[o.name for o in outs], sep="\n  ")
