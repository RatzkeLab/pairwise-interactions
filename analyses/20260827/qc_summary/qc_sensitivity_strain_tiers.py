"""Does restricting to sequence-VERIFIED strains change the 20260630 genomic-ML conclusions?

The join succeeding (76/76 well labels have a genomic row) is a lookup, not a verification. Per
strain, only 17 of 83 are confirmed by two independent references; 48 are confirmed by
corroborated_db alone; 2 are contradicted; 33 are untestable because they appear in no reference.

The disagreements are most likely genome_16S's fault -- it is the odd one out 16 times versus
corroborated_db once, and among corroborated-confirmed strains its median identity is 0.966 with
the correct label ranking 20th. But "probably the reference's fault" is not verification, so this
runs the models on progressively stricter strain sets and reports whether the conclusions move.

  all_strains        every strain that joins (the status quo)
  corroborated_ok    corroborated_db confirms the label at >=0.99, and does not contradict it
  both_agree         corroborated_db AND genome_16S both confirm -- the strictest set

If the honest number (cv_strain) survives the strict sets, the untestable strains are not
carrying the result and the status quo is safe. If it collapses, they were.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
EXPS = HERE.parent
sys.path.insert(0, str(EXPS / "20260630" / "analysis"))
sys.path.insert(0, str(EXPS / "shared_pipelines"))
from config import CFG
import genomic_ml as gm
import genomic_ml_yield as gy

OUT = HERE / "outputs"
OD = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/Karl_2026/"
          "Karl_20260704_OD_Full")


def tiers():
    w = pd.read_csv(OUT / "qc06_20260630_per_strain_status.csv")
    corr = "self_identity_corroborated_db(185)"
    gen = "self_identity_genome_16S"
    ok = set(w.loc[w[corr] >= 0.99, "label"])
    bad = set(w.loc[(w[corr].notna()) & (w[corr] < 0.99), "label"])
    both = set(w.loc[(w[corr] >= 0.99) & (w[gen] >= 0.99), "label"])
    return {"all_strains": None, "corroborated_ok": ok - bad, "both_agree": both}


def score(pairs, X, summ, gcfg, ctx_cls, model, target_is_ratio, keep, n_repeats=3):
    if keep is not None:
        pairs = pairs[pairs.strain_a.isin(keep) & pairs.strain_b.isin(keep)].reset_index(drop=True)
    if len(pairs) < 120:
        return {"n_pairs": len(pairs), "n_strains": len(set(pairs.strain_a) | set(pairs.strain_b)),
                "r2": np.nan, "rho": np.nan}
    out = []
    for rep in range(n_repeats):
        for tr, te in gm.make_folds(pairs, "cv_strain", n_splits=5, seed=97 * rep):
            if len(te) < 10 or len(tr) < 50:
                continue
            try:
                ctx = ctx_cls(gcfg, pairs, X, summ, tr, te)
                pred = (gm.MODELS if target_is_ratio else gy.MODELS)[model](ctx)
            except Exception:
                continue
            y = ctx.y_test
            ss = np.sum(y ** 2) if target_is_ratio else np.sum((y - y.mean()) ** 2)
            if ss > 0 and np.std(pred) > 0:
                out.append((1 - np.sum((y - pred) ** 2) / ss, spearmanr(y, pred)[0]))
    if not out:
        return {"n_pairs": len(pairs), "n_strains": len(set(pairs.strain_a) | set(pairs.strain_b)),
                "r2": np.nan, "rho": np.nan}
    a = np.array(out)
    return {"n_pairs": len(pairs), "n_strains": len(set(pairs.strain_a) | set(pairs.strain_b)),
            "r2": a[:, 0].mean(), "rho": a[:, 1].mean()}


def main():
    T = tiers()
    for k, v in T.items():
        print(f"{k}: {'all' if v is None else len(v)} strains")
    rows = []

    g_ra = gm.GenomicMLConfig(exp_cfg=CFG, out_dir=OUT / "_tmp_s")
    p_ra, _ = gm.build_dataset(g_ra)
    X_ra, s_ra = gm.strain_feature_matrix(g_ra, p_ra)
    for name, keep in T.items():
        r = score(p_ra, X_ra, s_ra, g_ra, gm._Ctx, "two_stage_ridge", True, keep)
        rows.append({"target": "relative_abundance", "tier": name, **r})
        print(f"  rel.abund {name:16} pairs={r['n_pairs']:5d} strains={r['n_strains']:3d} "
              f"R2={r['r2']:+.3f} rho={r['rho']:.3f}")

    g_y = gy.YieldMLConfig(exp_cfg=CFG, out_dir=OUT / "_tmp_sy",
                           layout_csv=CFG.exp_base / "setup" / "strain_layout_20260630.csv",
                           od_dir=OD)
    p_y, _, _ = gy.build_dataset(g_y)
    X_y, s_y = gy.strain_feature_matrix(g_y, p_y)
    for name, keep in T.items():
        r = score(p_y, X_y, s_y, g_y, gy._Ctx, "ridge_pca", False, keep)
        rows.append({"target": "yield_OD", "tier": name, **r})
        print(f"  yield     {name:16} pairs={r['n_pairs']:5d} strains={r['n_strains']:3d} "
              f"R2={r['r2']:+.3f} rho={r['rho']:.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "qc07_strain_tier_sensitivity.csv", index=False)
    print("\n=== sensitivity to strain verification tier (cv_strain) ===")
    print(df.round(3).to_string(index=False))
    return df


if __name__ == "__main__":
    main()
