"""The collaborator's panX feature-selection pipeline, re-run INSIDE each CV fold.

His pipeline (`go_feature_selection/pipeline_robust.ipynb`) is ported here as four functions:

    2.1  low-diversity removal     GO whose per-strain profile differs from its mode in fewer
                                   than MIN_SPECIES_DIVERSITY strains carry no usable variation
    2.2  XGBoost gain screen       keep the top GEN_CANDIDATE_GO by gain importance
    2.3  correlation redundancy    drop one of any pair of GO correlated above CORR_THRESHOLD
                                   across strain profiles, tie-broken by gain
    2.4  permutation pruning       iteratively drop GO with non-positive permutation importance
                                   under grouped CV, retrain, repeat until stable

Two adaptations, both forced by the data rather than chosen:

  **Single axis.** His selection ran on two axes (unseen genomes AND unseen nutrients) and kept
  a gene if it helped either. We have only the genome axis, so the union rule collapses to the
  genome axis alone.

  **Selection happens per outer fold, on training strains only.** This is the whole point of
  the exercise. The published 221-GO list was selected once, on all strains, using labels from
  experiments on those same strains -- so scoring it under strain-held-out CV is contaminated:
  the feature set has already seen the held-out strains' phenotypes. Here steps 2.1-2.4 see
  only the training strains of the fold they serve, so the resulting score is honest.

Controls, per fold, matched to whatever count the selection returns that fold:
  `random_matched`  GO drawn at random from the same post-diversity pool, prevalence-matched
  `all_features`    no selection at all

If nested selection lands on top of `random_matched`, the method adds real value and the
published list's apparent edge was leakage. If it ties, the collaborator's own note is right and
these GO are phylogenetic markers.
"""

import gc
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.inspection import permutation_importance
from sklearn.model_selection import GroupKFold

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE.parents[3] / "shared_pipelines"))
from config import CFG
import genomic_ml as gm

OUT = HERE / "outputs"
PANX = "panx_and_strains_table_final_post_filtering_strain_names.csv"

MIN_SPECIES_DIVERSITY = 3
GEN_CANDIDATE_GO = 2000
CORR_THRESHOLD = 0.7
PERM_N_FOLDS = 3
PERM_N_REPEATS = 1
MAX_PERM_ITERATIONS = 4          # his 30 is unaffordable nested inside an outer CV
N_OUTER_SPLITS = 5
XGB_PARAMS = dict(n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8,
                  colsample_bytree=0.3, reg_lambda=5.0, min_child_weight=5,
                  n_jobs=8, random_state=0, verbosity=0)


def pair_design(X, pairs):
    """Antisymmetric + symmetric blocks, as genomic_ml uses for this target."""
    A = X.loc[pairs.strain_a].values
    B = X.loc[pairs.strain_b].values
    return np.hstack([A - B, A + B])


def _go_importance(imp, go_cols):
    """Map importance on the 2*len(go) derived columns back to one score per GO."""
    n = len(go_cols)
    return pd.Series(np.maximum(imp[:n], imp[n:]), index=go_cols)


def select_features(X_tr, pairs_tr, go_cols, log, target_col):
    """Steps 2.1-2.4 on training strains / training pairs only."""
    import xgboost as xgb

    # --- 2.1 low-diversity removal
    prof = X_tr[go_cols]
    keep = []
    for c in go_cols:
        v = prof[c].values
        mode = pd.Series(v).mode().iloc[0]
        if (v != mode).sum() >= MIN_SPECIES_DIVERSITY:
            keep.append(c)
    go = keep
    log(f"    2.1 diversity: {len(go_cols)} -> {len(go)}")
    if len(go) < 20:
        return go

    # --- 2.2 gain screen
    D = pair_design(X_tr[go], pairs_tr)
    # genomic_ml's target column is named by the config (mean_log2_ratio_a_over_b), not "target"
    # -- that is the yield module's convention
    y = pairs_tr[target_col].values
    m = xgb.XGBRegressor(**XGB_PARAMS, importance_type="gain").fit(D, y)
    gain = _go_importance(m.feature_importances_, go)
    go = gain[gain > 0].sort_values(ascending=False).index[:GEN_CANDIDATE_GO].tolist()
    del m, D
    gc.collect()
    log(f"    2.2 gain screen: -> {len(go)}")
    if len(go) < 20:
        return go

    # --- 2.3 correlation redundancy on strain profiles, gain as tie-break
    prof = X_tr[go]
    corr = np.corrcoef(prof.values.T)
    corr = np.nan_to_num(corr)
    np.fill_diagonal(corr, 0)
    names = np.array(go)
    ii, jj = np.where(np.triu(np.abs(corr) > CORR_THRESHOLD, k=1))
    drop = set()
    for a_i, b_i in zip(ii, jj):
        a, b = names[a_i], names[b_i]
        if a in drop or b in drop:
            continue
        drop.add(b if gain.get(a, 0) >= gain.get(b, 0) else a)
    go = [c for c in go if c not in drop]
    del corr
    gc.collect()
    log(f"    2.3 correlation: removed {len(drop)} -> {len(go)}")
    if len(go) < 20:
        return go

    # --- 2.4 iterative permutation pruning, grouped by strain
    groups = pairs_tr.strain_a.values
    for it in range(1, MAX_PERM_ITERATIONS + 1):
        D = pair_design(X_tr[go], pairs_tr)
        gkf = GroupKFold(n_splits=min(PERM_N_FOLDS, len(np.unique(groups))))
        imps = []
        for tr_i, va_i in gkf.split(D, y, groups=groups):
            mm = xgb.XGBRegressor(**XGB_PARAMS).fit(D[tr_i], y[tr_i])
            pi = permutation_importance(mm, D[va_i], y[va_i], n_repeats=PERM_N_REPEATS,
                                        random_state=0, scoring="r2", n_jobs=-1)
            imps.append(pi.importances_mean)
            del mm
        imp = np.mean(imps, axis=0)
        score = _go_importance(imp, go)
        keep = score[score > 0].index.tolist()
        del D
        gc.collect()
        log(f"    2.4 iter {it}: {len(go)} -> {len(keep)}")
        if len(keep) < 20 or len(keep) == len(go):
            go = keep if len(keep) >= 20 else go
            break
        go = keep
    return go


def main():
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    gcfg = gm.GenomicMLConfig(exp_cfg=CFG, out_dir=OUT / "_tmp_nested", feature_table=PANX)
    pairs, _ = gm.build_dataset(gcfg)
    X, summ = gm.strain_feature_matrix(gcfg, pairs)
    go_all = list(X.columns)
    prevalence = X.sum(axis=0)
    log(f"panX: {len(pairs)} pairs, {X.shape[0]} strains, {len(go_all)} prevalence-filtered GO")

    rows, sel_log = [], []
    folds = list(gm.make_folds(pairs, "cv_strain", n_splits=N_OUTER_SPLITS, seed=0))
    for fi, (tr, te) in enumerate(folds):
        tr_strains = sorted(set(pairs.iloc[tr].strain_a) | set(pairs.iloc[tr].strain_b))
        log(f"\n  fold {fi}: {len(tr)} train pairs / {len(tr_strains)} train strains, "
            f"{len(te)} test pairs")
        sel = select_features(X.loc[tr_strains], pairs.iloc[tr].reset_index(drop=True),
                              go_all, log, gcfg.target)
        log(f"    SELECTED {len(sel)} GO  ({time.time()-t0:.0f}s elapsed)")
        sel_log.append({"fold": fi, "n_selected": len(sel), "selected": ";".join(sel[:400])})

        rng = np.random.default_rng(1234 + fi)
        # prevalence-matched random control of the SAME size, from the same pool
        dec = pd.qcut(prevalence.rank(method="first"), 20, labels=False)
        by_dec = {d: dec.index[dec == d].tolist() for d in np.unique(dec)}
        rnd = []
        for d in dec.loc[sel].values:
            pool = [c for c in by_dec[d] if c not in set(sel)]
            if pool:
                rnd.append(pool[int(rng.integers(len(pool)))])
        rnd = list(dict.fromkeys(rnd))

        for name, cols in [("nested_selection", sel), ("random_matched", rnd),
                           ("all_features", go_all)]:
            if len(cols) < 5:
                continue
            ctx = gm._Ctx(gcfg, pairs, X[cols], summ, tr, te)
            for model in ("ridge_pca", "xgboost_pca", "two_stage_ridge"):
                try:
                    pred = gm.MODELS[model](ctx)
                except Exception as e:
                    continue
                y = ctx.y_test
                rows.append({"fold": fi, "feature_set": name, "n_features": len(cols),
                             "model": model,
                             "r2": 1 - np.sum((y - pred) ** 2) / np.sum(y ** 2),
                             "rho": spearmanr(y, pred)[0] if np.std(pred) > 0 else np.nan})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s05_nested_selection_per_fold.csv", index=False)
    pd.DataFrame(sel_log).to_csv(OUT / "s05_selected_features_per_fold.csv", index=False)
    summ_df = (df.groupby(["feature_set", "model"])
                 .agg(n_folds=("fold", "nunique"), n_features=("n_features", "mean"),
                      r2_mean=("r2", "mean"), r2_sd=("r2", "std"),
                      rho_mean=("rho", "mean"), rho_sd=("rho", "std"))
                 .reset_index().sort_values("rho_mean", ascending=False))
    summ_df.to_csv(OUT / "s05_nested_selection_summary.csv", index=False)
    print("\n=== nested selection vs matched random vs all features (cv_strain) ===")
    print(summ_df.round(3).to_string(index=False))
    print(f"\ndone in {time.time()-t0:.0f}s")
    return summ_df


if __name__ == "__main__":
    main()
