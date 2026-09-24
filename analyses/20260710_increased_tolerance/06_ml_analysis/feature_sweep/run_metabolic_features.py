"""Do constructed metabolic pair-features beat the generic KO representation? Both targets.

Compares, under identical folds and models:
    ko_pca_baseline      the incumbent: PCA of KO presence/absence, generic sum/difference
    metabolic_only       ~20 constructed pair features (competition, complementarity, cross-feed)
    ko_plus_metabolic    both
    metabolic_no_crossfeed   as metabolic_only but dropping the block that depends on the
                             unverified BIOSYNTHESIS_MODULES list, so a wrong ID list cannot
                             silently drive the result
plus the usual `taxonomy_16S` and shuffled-label reference rows.

The yield target uses the symmetric features only: an antisymmetric feature cannot inform a
quantity that is unchanged when the two partners are relabelled.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE.parents[3] / "shared_pipelines"))
from config import CFG
import genomic_ml as gm
import genomic_ml_yield as gy
import metabolic_features as mf

OUT = HERE / "outputs"
GEN = gm.GENOMIC_DIR
OD = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/Karl_2026/"
          "Karl_20260704_OD_Full")
N_SPLITS, N_REPEATS = 5, 3


def metabolic_blocks(pairs, strains, w2s, include_crossfeeding=True):
    cazy = mf.load_table(GEN / "CAZy_and_strains_table.csv", strains, w2s)
    mods = mf.load_table(GEN / "KEGG_Module_and_strains_table.csv", strains, w2s)
    ko = mf.load_table(GEN / "KEGG_ko_and_strains_table.csv", strains, w2s)
    return mf.build(pairs, cazy=cazy, modules=mods, ko=ko,
                    include_crossfeeding=include_crossfeeding)


def evaluate(target, pairs, X, summ, gcfg, ctx_cls, models, antisym_ok):
    w2s, _ = gm.load_strain_mapping(gcfg)
    strains = sorted(set(pairs.strain_a) | set(pairs.strain_b))
    F_all, sym, anti = metabolic_blocks(pairs, strains, w2s, True)
    F_nc, sym_nc, anti_nc = metabolic_blocks(pairs, strains, w2s, False)
    cols_all = sym + (anti if antisym_ok else [])
    cols_nc = sym_nc + (anti_nc if antisym_ok else [])
    print(f"  {target}: {len(cols_all)} constructed features "
          f"({len(sym)} symmetric, {len(anti)} antisymmetric)")

    rows = []
    for rep in range(N_REPEATS):
        for regime in ("cv_pair", "cv_strain"):
            for fi, (tr, te) in enumerate(gm.make_folds(pairs, regime, n_splits=N_SPLITS,
                                                        seed=97 * rep)):
                if len(te) < 10:
                    continue
                ctx = ctx_cls(gcfg, pairs, X, summ, tr, te)
                y = ctx.y_test
                ss = np.sum(y ** 2) if target == "relative_abundance" else \
                    np.sum((y - y.mean()) ** 2)

                def score(name, pred):
                    if not np.all(np.isfinite(pred)) or np.std(pred) == 0:
                        return
                    rows.append({"target": target, "regime": regime, "rep": rep, "fold": fi,
                                 "feature_set": name,
                                 "r2": 1 - np.sum((y - pred) ** 2) / ss,
                                 "rho": spearmanr(y, pred)[0]})

                for mname in models:
                    try:
                        score(f"ko_pca_baseline::{mname}", gm.MODELS[mname](ctx)
                              if target == "relative_abundance" else gy.MODELS[mname](ctx))
                    except Exception:
                        pass
                for tag, cols in (("metabolic_only", cols_all),
                                  ("metabolic_no_crossfeed", cols_nc)):
                    Fm = F_all if tag == "metabolic_only" else F_nc
                    sc = StandardScaler().fit(Fm[cols].values[tr])
                    Z = sc.transform(Fm[cols].values)
                    score(f"{tag}::ridge", Ridge(alpha=10).fit(Z[tr], ctx.y_train).predict(Z[te]))
                # KO PCA design + constructed features side by side
                sc = StandardScaler().fit(F_all[cols_all].values[tr])
                Zm = sc.transform(F_all[cols_all].values)
                D = ctx.fwd if target == "relative_abundance" else ctx.D
                Dc = np.hstack([D, Zm])
                score("ko_plus_metabolic::ridge",
                      Ridge(alpha=50).fit(Dc[tr], ctx.y_train).predict(Dc[te]))
    return pd.DataFrame(rows)


def main():
    out = []
    gcfg_ra = gm.GenomicMLConfig(exp_cfg=CFG, out_dir=OUT / "_tmp_m")
    pairs, _ = gm.build_dataset(gcfg_ra)
    X, summ = gm.strain_feature_matrix(gcfg_ra, pairs)
    out.append(evaluate("relative_abundance", pairs, X, summ, gcfg_ra, gm._Ctx,
                        ["ridge_pca", "two_stage_ridge"], antisym_ok=True))

    gcfg_y = gy.YieldMLConfig(exp_cfg=CFG, out_dir=OUT / "_tmp_my",
                              layout_csv=CFG.exp_base / "setup" / "strain_layout_20260630.csv",
                              od_dir=OD)
    py, _, _ = gy.build_dataset(gcfg_y)
    Xy, sy = gy.strain_feature_matrix(gcfg_y, py)
    out.append(evaluate("yield_OD", py, Xy, sy, gcfg_y, gy._Ctx,
                        ["ridge_pca", "two_stage_ridge"], antisym_ok=False))

    df = pd.concat(out, ignore_index=True)
    df.to_csv(OUT / "s06_metabolic_features_per_fold.csv", index=False)
    s = (df.groupby(["target", "regime", "feature_set"])
           .agg(r2_mean=("r2", "mean"), r2_sd=("r2", "std"),
                rho_mean=("rho", "mean"), rho_sd=("rho", "std"), n=("fold", "size"))
           .reset_index().sort_values(["target", "regime", "rho_mean"], ascending=[True, True, False]))
    s.to_csv(OUT / "s06_metabolic_features_summary.csv", index=False)
    for t in s.target.unique():
        for r in ("cv_strain", "cv_pair"):
            d = s[(s.target == t) & (s.regime == r)]
            print(f"\n=== {t} | {r} ===")
            print(d[["feature_set", "r2_mean", "r2_sd", "rho_mean", "rho_sd"]].round(3)
                  .to_string(index=False))
    return s


if __name__ == "__main__":
    main()
