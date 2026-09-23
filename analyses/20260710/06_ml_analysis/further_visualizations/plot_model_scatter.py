"""Predicted-vs-observed scatters for models the main report does not plot.

`make_all_figures` deliberately plots only the winner of each regime (f02), which keeps the
report on the top performers. This is the escape hatch: the same panel for ANY model x regime
in g03_cv_predictions.csv -- e.g. xgboost_raw_ko under cv_strain, which ties two_stage_ridge on
Spearman (0.574 vs 0.576) while losing on R² (0.197 vs 0.321), so seeing where its points sit
is the only way to tell whether that is a ranking failure or a calibration failure.

Reads the saved CV predictions; it never refits, so it is instant and can never disagree with
the run that produced the report.

Every score here is POOLED over folds, which is NOT what summ_df reports -- summ_df scores each
fold separately and averages. Both are printed on each panel so the two can never be confused
(see the header of the R² column below for which one to quote in text).

    python plot_model_scatter.py                                  # xgboost_raw_ko + two_stage_ridge, cv_strain
    python plot_model_scatter.py --models xgboost_pca --regime cv_pair
    python plot_model_scatter.py --models all
    python plot_model_scatter.py --format png              # raster, for a quick look
    python plot_model_scatter.py --models xgboost_raw_ko --regime both   # f02 layout, one model
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "shared_scripts"))
import genomic_ml as gm            # palette + rcParams, so these read as one report
import matplotlib.pyplot as plt

SRC = HERE.parents[0] / "genomic_ml" / "outputs"
OUT = HERE / "outputs"
DEFAULT_MODELS = ["xgboost_raw_ko", "two_stage_ridge"]


def pooled_metrics(p):
    """Scored exactly as genomic_ml does: R² is vs. predicting no winner (0), not vs. the mean."""
    y, yhat = p["y_true"].values, p["y_pred"].values
    decisive = np.abs(y) > 0.5
    return {
        "n": len(y),
        "r2": 1 - np.sum((y - yhat) ** 2) / np.sum(y ** 2),
        "spearman_rho": spearmanr(yhat, y)[0],
        "pearson_r": pearsonr(yhat, y)[0],
        "sign_accuracy": np.mean(np.sign(y[decisive]) == np.sign(yhat[decisive])),
    }


def panel(ax, p, model, regime, summ):
    ax.scatter(p["y_pred"], p["y_true"], s=14, alpha=0.45, color=gm.COLOR_BLUE, edgecolor="none")
    lim = max(np.abs(p["y_true"]).max(), np.abs(p["y_pred"]).max()) * 1.05
    ax.plot([-lim, lim], [-lim, lim], color=gm.COLOR_TEXT_SECONDARY, lw=1, ls="--")
    ax.axhline(0, color=gm.COLOR_GRID, lw=1)
    ax.axvline(0, color=gm.COLOR_GRID, lw=1)

    m = pooled_metrics(p)
    ax.set_title(f"{regime} -- {model}\npooled ρ = {m['spearman_rho']:.2f}, "
                 f"pooled R² = {m['r2']:.2f} vs. no winner\nn = {m['n']}")

    # the mean-of-folds numbers, i.e. what summ_df / the report quote, for direct comparison
    row = summ[(summ["regime"] == regime) & (summ["model"] == model)]
    if len(row):
        r = row.iloc[0]
        ax.text(0.03, 0.97,
                f"mean-of-folds (summ_df):\nρ = {r['spearman_rho_mean']:.3f} ± {r['spearman_rho_sd']:.3f}\n"
                f"R² = {r['r2_mean']:.3f} ± {r['r2_sd']:.3f}\nover {int(r['n_folds'])} folds",
                transform=ax.transAxes, va="top", ha="left", fontsize=8,
                color=gm.COLOR_TEXT_SECONDARY,
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=gm.COLOR_GRID, lw=0.8))
    ax.set_xlabel("predicted log2 ratio (a/b)")
    ax.set_ylabel("observed log2 ratio (a/b)")
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default=",".join(DEFAULT_MODELS),
                    help="comma-separated model names, or 'all'")
    ap.add_argument("--regime", default="cv_strain", choices=["cv_strain", "cv_pair", "both"],
                    help="'both' rebuilds f02's layout with ONE model held fixed across the "
                         "two regimes -- f02 itself shows the winner of each, so its panels "
                         "are different models and the drop between them is not attributable")
    ap.add_argument("--src", default=str(SRC), help="directory holding g03_cv_*.csv")
    ap.add_argument("--format", default=gm.FIG_FORMAT, choices=["svg", "png"],
                    help="vector by default, matching the main report")
    args = ap.parse_args()

    src = Path(args.src)
    all_pairs = pd.read_csv(src / "g03_cv_predictions.csv")
    summ = pd.read_csv(src / "g03_cv_summary.csv")

    if args.models == "all":
        models = [m for m in sorted(all_pairs["model"].unique())
                  if m not in ("zero_baseline", "strength_observed_no_genomics")
                  and not m.startswith("SHUFFLED")]
    else:
        models = [m.strip() for m in args.models.split(",") if m.strip()]
    missing = [m for m in models if m not in set(all_pairs["model"])]
    if missing:
        raise SystemExit(f"not in {src.name}/g03_cv_predictions.csv: {missing}\n"
                         f"available: {sorted(all_pairs['model'].unique())}")

    if args.regime == "both":
        if len(models) != 1:
            raise SystemExit("--regime both compares ONE model across the two regimes; "
                             "pass a single --models value")
        specs = [(models[0], r) for r in ("cv_pair", "cv_strain")]
        stem = f"fv02_regime_gap_{models[0]}"
        suptitle = (f"{models[0]}: the same model under both regimes "
                    f"(f02's panels are two different models)")
        csv_tag, prefix = f"regime_gap_{models[0]}", "fv02"
    else:
        specs = [(m, args.regime) for m in models]
        tag = "_".join(models) if len(models) <= 3 else f"{len(models)}models"
        stem = f"fv01_predicted_vs_observed_{args.regime}_{tag}"
        suptitle = "Predicted vs. observed log2 abundance ratio, pooled over CV folds"
        csv_tag, prefix = args.regime, "fv01"

    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, len(specs), figsize=(5.5 * len(specs), 5.4), squeeze=False)
    rows = []
    for ax, (model, regime) in zip(axes[0], specs):
        p = all_pairs[(all_pairs["regime"] == regime) & (all_pairs["model"] == model)]
        rows.append({"regime": regime, "model": model, **panel(ax, p, model, regime, summ)})
    fig.suptitle(suptitle, fontweight="bold")
    fig.tight_layout()
    img = gm._savefig(fig, OUT, stem, args.format)

    tbl = pd.DataFrame(rows)
    csv = OUT / f"{prefix}_pooled_metrics_{csv_tag}.csv"
    tbl.to_csv(csv, index=False)
    print(tbl.round(3).to_string(index=False))
    print(f"\n-> {img}\n-> {csv}")


if __name__ == "__main__":
    main()
