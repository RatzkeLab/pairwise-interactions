"""What could a perfect model score? The replicate ceiling, on the ML target's own scale.

`replicate_report.py` answers "how reproducible is a measurement" for two modalities, and does
it in Spearman because OD (plate-z-scored) and relative abundance ([0,1]) have no common scale.
Neither of those is the quantity `genomic_ml` predicts, so neither ceiling can be set against a
cv_strain R². This closes that gap: the same replicate logic, applied to
`mean_log2_ratio_a_over_b` over exactly the pairs in the modeling set, scored in genomic_ml's
own convention (R² against predicting no winner, i.e. 0).

Three numbers, and they are NOT interchangeable:

  rep -> rep        one replicate well used as a literal prediction of another. This is the
                    thing people mean by "what happens when replicates predict each other",
                    but it is NOT the ceiling: both sides carry a full dose of measurement
                    noise, while a model is only ever scored against the *averaged* label.
                    It understates what a perfect model could reach.

  label ceiling     1 - var(label noise)/E[y^2]. The label is the mean of n replicate wells,
                    so its noise is sigma^2/n. This is the ceiling to compare models against.
                    Computed at the ACTUAL replicate mix of the modeling set, which matters a
                    lot here: only 244 of 1465 pairs have a second well, so 83% of labels are
                    single wells carrying twice the variance of a 2-well label. A ceiling
                    quoted from the replicated pairs alone is not the ceiling of this dataset.

  rho ceiling       sqrt(label reliability) -- the classic attenuation bound on how well any
                    predictor of the truth can correlate with a noisy measurement of it.

Read against genomic_ml's independent `binomial_read_depth` estimate (g02), which is derived
from read counts rather than replication and covers all 1465 pairs rather than the 244 with a
second well. Where the two disagree, the binomial one is the more conservative and the more
representative -- see CAVEATS at the bottom of the printout.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

HERE = Path(__file__).resolve().parent
EXPS = HERE.parents[1]
OUT = HERE / "outputs"
sys.path.insert(0, str(EXPS / "shared_scripts"))

EXP = "20260630"
WELLS = EXPS / EXP / "05_engineer_relative_abundances/relative_abundance/outputs/r02_well_interaction_scores.csv"
ML = EXPS / EXP / "06_ml_analysis/genomic_ml/outputs"
TARGET = "mean_log2_ratio_a_over_b"


def r2_vs_zero(y, yhat):
    """genomic_ml's convention: the denominator is the error of predicting 0 (no winner)."""
    return float(1 - np.sum((y - yhat) ** 2) / np.sum(y ** 2))


def main():
    pairs = pd.read_csv(ML / "g01_modeling_pairs.csv")
    wells = pd.read_csv(WELLS)
    for d in (pairs, wells):
        d["key"] = [frozenset((a, b)) for a, b in zip(d.strain_a, d.strain_b)]

    # the well table is already oriented on the same strain_a/strain_b as the modeling set;
    # assert rather than assume, since a silent flip would negate half the replicate deltas
    canon = dict(zip(pairs.key, pairs.strain_a))
    w = wells[wells.key.isin(canon)].copy()
    flipped = w.strain_a.values != np.array([canon[k] for k in w.key])
    assert not flipped.any(), f"{flipped.sum()} wells oriented against the modeling set"

    reps = w.groupby("key")["log2_ratio_a_over_b"].agg(list)
    reps = reps[reps.map(len) >= 2]
    n_obs = reps.map(len)
    assert (n_obs == 2).all(), f"expected 2-well replicates only, saw {sorted(set(n_obs))}"
    x = np.array([v[0] for v in reps])
    y = np.array([v[1] for v in reps])

    y_all = pairs[TARGET].values
    Ey2 = float(np.mean(y_all ** 2))
    var_y = float(np.var(y_all))
    # for n=2, (x-y)^2/2 is an unbiased estimate of the single-well variance
    s2_single = float(np.mean((x - y) ** 2 / 2))
    n_rep = pairs["n_replicates"].values
    s2_label_actual = float(np.mean(s2_single / n_rep))      # at the real 1221x n=1 / 244x n=2 mix

    rows = [
        {"quantity": "rep -> rep (1 well predicts 1 well)", "n": len(x),
         "r2_vs_zero": r2_vs_zero(y, x), "spearman": float(spearmanr(x, y)[0]),
         "pearson": float(pearsonr(x, y)[0])},
        {"quantity": "ceiling: label = mean of 2 wells", "n": len(x),
         "r2_vs_zero": 1 - (s2_single / 2) / Ey2,
         "spearman": float(np.sqrt(max(0.0, 1 - (s2_single / 2) / var_y))), "pearson": np.nan},
        {"quantity": "ceiling: label at ACTUAL replicate mix", "n": len(pairs),
         "r2_vs_zero": 1 - s2_label_actual / Ey2,
         "spearman": float(np.sqrt(max(0.0, 1 - s2_label_actual / var_y))), "pearson": np.nan},
        {"quantity": "ceiling: label = single well", "n": int((n_rep == 1).sum()),
         "r2_vs_zero": 1 - s2_single / Ey2,
         "spearman": float(np.sqrt(max(0.0, 1 - s2_single / var_y))), "pearson": np.nan},
    ]

    ceil = pd.read_csv(ML / "g02_label_noise_ceiling.csv").set_index("estimator")
    rows.append({"quantity": "ceiling: binomial read depth (g02 / all pairs)", "n": len(pairs),
                 "r2_vs_zero": float(ceil.loc["binomial_read_depth", "ceiling_r2"]),
                 "spearman": float(np.sqrt(max(0.0, ceil.loc["binomial_read_depth", "ceiling_r2"]))),
                 "pearson": np.nan})

    summ = pd.read_csv(ML / "g03_cv_summary.csv")
    for regime in ("cv_pair", "cv_strain"):
        d = summ[(summ.regime == regime) & (~summ.model.str.startswith("SHUFFLED"))
                 & (~summ.model.isin(["zero_baseline", "strength_observed_no_genomics"]))]
        best = d.iloc[0]
        rows.append({"quantity": f"MODEL: {regime} best ({best.model})", "n": int(best.n_test),
                     "r2_vs_zero": float(best.r2_mean), "spearman": float(best.spearman_rho_mean),
                     "pearson": float(best.pearson_r_mean)})

    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "r04_ml_target_ceiling.csv", index=False)

    print(f"=== {EXP}: replicate ceiling on {TARGET}, genomic_ml's R² convention ===\n")
    print(out.round(3).to_string(index=False))
    print(f"\nsigma^2 single well {s2_single:.3f} | E[y^2] {Ey2:.3f} | var(y) {var_y:.3f}")
    print(f"replicate mix: {(n_rep == 1).sum()} pairs x 1 well, {(n_rep == 2).sum()} pairs x 2 wells")
    print("""
CAVEATS
  - The 244 replicated pairs are not a random sample: `unstable_replicate` pairs (the ones
    whose wells disagreed most) are dropped upstream, so this estimate is selected for
    agreement and is optimistic. The binomial row, derived from read counts over all 1465
    pairs, is the conservative reading.
  - `rep -> rep` is the honest answer to "what do replicates achieve against each other", but
    it is the wrong bar for a model: it is depressed by noise on BOTH sides.
  - Everything here is for 20260630 only, and R² is vs. predicting no winner -- not comparable
    to genomic_ml_yield, which scores against the mean.""")
    return out


if __name__ == "__main__":
    main()
