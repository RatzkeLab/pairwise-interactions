"""Would more technical replicates make the model better? And how close is it to the ceiling?

Both pipelines train on the replicate MEAN (genomic_ml_yield aggregates od_z with .mean();
genomic_ml uses r03's mean_log2_ratio_a_over_b). So label noise is already suppressed by
averaging, and the question is whether suppressing it further would buy anything.

Answered two ways rather than by argument alone:

  **Attenuation arithmetic.** Single-well reliability r1 is measured 1-vs-1 on pairs with
  exactly two wells. Spearman-Brown then gives the reliability of a mean of k wells,
  r_k = k*r1 / (1 + (k-1)*r1). Observed R^2 relates to the noise-free R^2 by
  R2_obs = R2_true * r_k, so the gain available from k -> infinity is bounded and computable.

  **A replicate learning curve.** Rebuild the target from k = 1, 2, 3 randomly chosen wells per
  pair and re-run cv_strain. Restricted to pairs that HAVE three wells, so the same pairs are
  used at every k and the curve is not confounded with which pairs are included.

And the alternative use of the same effort, for comparison: a **strain learning curve** --
subsample the collection to S genomes and re-run. The binding constraint on every analysis in
this project has been the number of genomes, not the number of measurements per pair, and this
puts a number on that.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import Ridge

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE.parents[3] / "shared_pipelines"))
from config import CFG
import genomic_ml_yield as gy
import genomic_ml as gm

OD = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/Karl_2026/"
          "Karl_20260704_OD_Full")
N_SPLITS, N_REPEATS = 5, 3


def cv_strain_r2(gcfg, pairs, X, summ, seed=0, n_repeats=N_REPEATS):
    """cv_strain R² with ridge_pca -- the best cv_strain model, and fast enough to sweep."""
    out = []
    for rep in range(n_repeats):
        for tr, te in gm.make_folds(pairs, "cv_strain", n_splits=N_SPLITS, seed=seed + 97 * rep):
            if len(te) < 10 or len(tr) < 50:
                continue
            ctx = gy._Ctx(gcfg, pairs, X, summ, tr, te)
            pred = gy.m_ridge(ctx)
            y = ctx.y_test
            ss = np.sum((y - y.mean()) ** 2)
            if ss > 0:
                out.append({"r2": 1 - np.sum((y - pred) ** 2) / ss,
                            "rho": spearmanr(y, pred)[0]})
    d = pd.DataFrame(out)
    return d.r2.mean(), d.rho.mean(), len(d)


def main():
    gcfg = gy.YieldMLConfig(exp_cfg=CFG,
                            layout_csv=CFG.exp_base / "setup" / "strain_layout_20260630.csv",
                            od_dir=OD, out_dir=HERE / "outputs")
    pairs, _, wells = gy.build_dataset(gcfg)
    X, summ = gy.strain_feature_matrix(gcfg, pairs)

    wells = wells.copy()
    wells["pair_key"] = [frozenset((a, b)) for a, b in zip(wells.strain1, wells.strain2)]
    per_pair = wells.groupby("pair_key")["od_z"].apply(list)

    # ---- 1. single-well reliability, measured 1-vs-1 on exactly-2-well pairs
    two = per_pair[per_pair.map(len) == 2]
    r1 = pearsonr([v[0] for v in two], [v[1] for v in two])[0]
    rows = [{"k_replicates": k, "reliability_of_mean": k * r1 / (1 + (k - 1) * r1)}
            for k in [1, 2, 3, 4, 6, 10]]
    rows.append({"k_replicates": np.inf, "reliability_of_mean": 1.0})
    rel = pd.DataFrame(rows)
    rel["max_R2_if_true_R2_matches_observed"] = np.nan
    print(f"single-well reliability r1 = {r1:.3f}  (1-vs-1 on {len(two)} two-well pairs)")
    print(rel.round(4).to_string(index=False))

    # ---- 2. replicate learning curve, on pairs that have 3 wells
    three = per_pair[per_pair.map(len) >= 3]
    keys3 = set(three.index)
    p3 = pairs[[frozenset((a, b)) in keys3 for a, b in zip(pairs.strain_a, pairs.strain_b)]].copy()
    print(f"\nreplicate learning curve on the {len(p3)} pairs with >=3 wells")
    curve = []
    for k in (1, 2, 3):
        r2s, rhos = [], []
        for s in range(3):
            rng = np.random.default_rng(100 + s)
            tgt = {key: float(np.mean(rng.choice(v, size=k, replace=False)))
                   for key, v in three.items()}
            pk = p3.copy()
            pk["target"] = [tgt[frozenset((a, b))] for a, b in zip(pk.strain_a, pk.strain_b)]
            r2, rho, _ = cv_strain_r2(gcfg, pk.reset_index(drop=True), X, summ, seed=s)
            r2s.append(r2)
            rhos.append(rho)
        curve.append({"k_replicates": k, "cv_strain_r2": np.mean(r2s), "r2_sd": np.std(r2s),
                      "cv_strain_rho": np.mean(rhos), "n_pairs": len(p3)})
        print(f"  k={k}: cv_strain R2={np.mean(r2s):+.3f} (sd {np.std(r2s):.3f})  "
              f"rho={np.mean(rhos):.3f}")
    curve = pd.DataFrame(curve)

    # ---- 3. strain learning curve
    strains = sorted(set(pairs.strain_a) | set(pairs.strain_b))
    print(f"\nstrain learning curve (full set = {len(strains)} genomes)")
    scurve = []
    for S in (30, 45, 60, 75, len(strains)):
        r2s, rhos, npairs = [], [], []
        for s in range(3):
            rng = np.random.default_rng(200 + s)
            keep = set(rng.choice(strains, size=min(S, len(strains)), replace=False))
            sub = pairs[pairs.strain_a.isin(keep) & pairs.strain_b.isin(keep)].reset_index(drop=True)
            if len(sub) < 150:
                continue
            r2, rho, _ = cv_strain_r2(gcfg, sub, X, summ, seed=s, n_repeats=2)
            # at small S the ridge extrapolates wildly on held-out strains and R^2 can blow up
            # by many orders of magnitude; rho stays bounded and is the honest readout there
            r2s.append(r2 if np.isfinite(r2) and r2 > -5 else np.nan)
            rhos.append(rho)
            npairs.append(len(sub))
        if r2s:
            scurve.append({"n_strains": S, "cv_strain_r2": np.nanmean(r2s), "r2_sd": np.nanstd(r2s),
                           "r2_unstable": bool(np.isnan(r2s).any()),
                           "cv_strain_rho": np.mean(rhos), "mean_n_pairs": np.mean(npairs)})
            print(f"  S={S:3d} strains ({np.mean(npairs):.0f} pairs): "
                  f"cv_strain R2={np.nanmean(r2s):+.3f}  rho={np.mean(rhos):.3f}"
                  f"{'  [R2 unstable -- read rho]' if np.isnan(r2s).any() else ''}")
    scurve = pd.DataFrame(scurve)

    rel.to_csv(HERE / "outputs" / "y05_reliability_spearman_brown.csv", index=False)
    curve.to_csv(HERE / "outputs" / "y05_replicate_learning_curve.csv", index=False)
    scurve.to_csv(HERE / "outputs" / "y05_strain_learning_curve.csv", index=False)
    return rel, curve, scurve


if __name__ == "__main__":
    main()
