"""The discriminating number: how much of the starting ratio survives to the endpoint?

Each of 50 pairs was set up at five inoculum ratios spanning ~49x
(25:175 ... 175:25), all five on the SAME plate so a dose-response cannot be a
plate effect. Regressing final log2 read ratio on initial log2 inoculum ratio
gives a slope with a direct reading:

    slope ~ 0   the starting ratio is irrelevant; the outcome is set by fitness,
                and the assay measures what it is supposed to measure
    slope ~ 1   the starting ratio passes straight through; the assay is largely
                reporting the inoculum back
    in between  both matter, and the slope says how much

The OD arms already showed total yield is at carrying capacity (1.8% of a
doubling survives), but yield and composition are different questions: a well can
reach the same final density from any starting ratio while its *composition*
still remembers the inoculum. This is the measurement that settles composition.

Pairs were stratified half strong-winner (|log2 ratio| > 3) and half near-neutral
(< 1), because a head start is most likely to flip a strong-winner pair while a
ratio effect is easiest to detect in a neutral one. They are reported separately:
a null in only one stratum would be uninformative.
"""
import sys
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE)); sys.path.insert(0, str(BASE.parent / "shared_scripts"))

import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config
from io_utils import pair_key

OUT = BASE / "06_inoculum_arms" / "outputs"; OUT.mkdir(parents=True, exist_ok=True)
MIN_LEVELS = 4          # a slope from fewer than 4 of the 5 levels is not worth fitting
MAX_UNCERTAINTY = 0.3   # matches relative_abundance's HIGH_UNCERTAINTY_THRESHOLD


def load(cfg):
    wells = pd.read_csv(cfg.relative_abundance_out_dir / "r02_well_interaction_scores.csv")
    lay = pd.read_csv(BASE / "01_setup" / "strain_layout_20260907.csv")
    lay["sample_id"] = lay.apply(
        lambda r: f"Plate{int(r.dest_plate):02d}_{r.dest_well}", axis=1)
    df = wells.merge(lay[["sample_id", "dest_plate", "vol1_nL", "vol2_nL",
                          "log2_inoculum_ratio", "separation_bp"]],
                     on="sample_id", how="left")
    df = df[df.well_type.isin(["ratio", "density"])].copy()
    df = df[~df.missing_reference]

    # layout's log2_inoculum_ratio is strain1/strain2; the score is a/b with a,b sorted
    flip = [pair_key(r.strain1, r.strain2)[0] != r.strain1 for r in df.itertuples()]
    df["log2_inoculum_a_over_b"] = np.where(flip, -df.log2_inoculum_ratio,
                                            df.log2_inoculum_ratio)
    df["pair"] = [pair_key(r.strain1, r.strain2) for r in df.itertuples()]
    return df


def fit_slopes(df):
    rows = []
    for pair, g in df.groupby("pair"):
        g = g[g.uncertainty_score <= MAX_UNCERTAINTY] if "uncertainty_score" in g else g
        g = g.dropna(subset=["log2_ratio_a_over_b", "log2_inoculum_a_over_b"])
        if g.log2_inoculum_a_over_b.nunique() < MIN_LEVELS:
            continue
        x, y = g.log2_inoculum_a_over_b.values, g.log2_ratio_a_over_b.values
        res = stats.linregress(x, y)
        rows.append(dict(strain_a=pair[0], strain_b=pair[1], n_levels=len(g),
                         slope=res.slope, intercept=res.intercept,
                         r2=res.rvalue ** 2, p=res.pvalue,
                         stderr=res.stderr,
                         mean_outcome=float(np.mean(y)),
                         separation_bp=g.separation_bp.iloc[0]))
    return pd.DataFrame(rows)


def report(slopes, df):
    if slopes.empty:
        print("no pair had enough usable levels to fit a slope"); return
    print("=" * 74)
    print("RATIO TITRATION -- does the starting ratio carry through to the endpoint?")
    print("=" * 74)
    print(f"  {len(slopes)} pairs with >= {MIN_LEVELS} usable inoculum levels")
    s = slopes.slope
    print(f"\n  slope of final log2 ratio on initial log2 inoculum ratio:")
    print(f"    median {s.median():+.3f}   mean {s.mean():+.3f}   "
          f"IQR {s.quantile(.25):+.3f} .. {s.quantile(.75):+.3f}")
    t0 = stats.wilcoxon(s) if len(s) > 5 else None
    print(f"    vs 0 (pure fitness): Wilcoxon p = {t0.pvalue:.3g}" if t0 else "")
    t1 = stats.wilcoxon(s - 1) if len(s) > 5 else None
    print(f"    vs 1 (pure pass-through): Wilcoxon p = {t1.pvalue:.3g}" if t1 else "")
    print(f"\n  => roughly {s.median()*100:.0f}% of a change in starting ratio survives "
          f"to the endpoint.")
    if abs(s.median()) < 0.2:
        print("     The assay reports fitness. The inoculum is close to irrelevant, "
              "which is\n     what the design hoped for and what makes the pair wells "
              "interpretable.")
    elif s.median() > 0.7:
        print("     The assay largely reports the inoculum back. Interaction scores "
              "from single\n     wells are confounded by starting density and need "
              "OD-normalised inoculum.")
    else:
        print("     Both matter. Interaction scores carry a real fitness signal but "
              "are shifted by\n     starting density; the slope is the correction factor.")

    # stratify: was the pair a strong winner or near-neutral at 1:1?
    mid = df[np.isclose(df.log2_inoculum_a_over_b, 0)]
    strength = mid.groupby("pair").log2_ratio_a_over_b.mean().abs()
    slopes = slopes.assign(pair=list(zip(slopes.strain_a, slopes.strain_b)))
    slopes["balanced_strength"] = slopes.pair.map(strength)
    strong = slopes[slopes.balanced_strength > 3]
    neutral = slopes[slopes.balanced_strength < 1]
    print(f"\n  by stratum (strength measured at the 1:1 level):")
    for lab, sub in [("strong winner (|log2|>3)", strong), ("near-neutral (<1)", neutral)]:
        if len(sub):
            print(f"    {lab:26s} n={len(sub):3d}  median slope {sub.slope.median():+.3f}")
    print("\n  A null in both strata is the strong result; a null in only one is not.")
    return slopes


def main():
    cfg = config.make_config()
    df = load(cfg)
    slopes = fit_slopes(df)
    slopes = report(slopes, df)
    if slopes is None or slopes.empty:
        return
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for pair, g in df.groupby("pair"):
        g = g.sort_values("log2_inoculum_a_over_b")
        ax[0].plot(g.log2_inoculum_a_over_b, g.log2_ratio_a_over_b,
                   marker="o", ms=3, lw=.7, alpha=.45, color="#2a78d6")
    lim = df.log2_inoculum_a_over_b.abs().max() * 1.1
    ax[0].plot([-lim, lim], [-lim, lim], ls=":", c="crimson", lw=1.2,
               label="slope 1: inoculum passes through")
    ax[0].axhline(0, ls="--", c="grey", lw=1)
    ax[0].set(xlabel="initial log2 inoculum ratio (a/b)",
              ylabel="final log2 read ratio (a/b)",
              title="ratio titration, one line per pair")
    ax[0].legend(fontsize=8)
    ax[1].hist(slopes.slope, bins=20, color="#2a78d6")
    ax[1].axvline(0, c="grey", ls="--", lw=1)
    ax[1].axvline(1, c="crimson", ls=":", lw=1.2)
    ax[1].axvline(slopes.slope.median(), c="#0ca30c", lw=2,
                  label=f"median {slopes.slope.median():+.2f}")
    ax[1].set(xlabel="fitted slope", ylabel="pairs", title="how much of the head start survives")
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "i04_ratio_titration.png", dpi=150)
    slopes.drop(columns=["pair"]).to_csv(OUT / "i05_ratio_slopes.csv", index=False)
    df.to_csv(OUT / "i06_titration_wells.csv", index=False)
    print(f"\nwrote -> {OUT}")


if __name__ == "__main__":
    main()
