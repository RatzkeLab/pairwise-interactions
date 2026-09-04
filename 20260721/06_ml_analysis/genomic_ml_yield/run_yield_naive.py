"""Yield prediction on 20260721 under the NAIVE strain join, against permuted-mapping nulls.

Deliberately "forgetting" everything established about this experiment: the wells are joined to
the genomic tables on `Well_souce_plate` exactly as if the labels were trusted. The question is
only whether the resulting model beats chance.

This is a better-powered version of `genomic_ml_join_test/` (which used relative abundance):
the yield target does not depend on 16S resolvability, so nothing is discarded as
`high_uncertainty` and the strain count rises. If a naive join carries any signal at all, this
is the analysis most likely to see it.

Two nulls, and the second is the one that matters:

  `SHUFFLED_*`      -- labels permuted. Detects "the model learned nothing at all".
  permuted mapping  -- which GENOME is attached to which well is shuffled, preserving the
                       genomes, the labels, the pair structure and (crucially) the
                       fingerprint-as-ID property. A consistently-WRONG bijection still gives
                       every well a unique, stable KO vector, so under cv_pair a model can score
                       well while knowing no biology. Only this null, read under cv_strain,
                       distinguishes a correct join from a scrambled one.

Layout choice: the physical contents of the destination plates follow `echo_strains_20260721`,
which agrees with `strain_layout_20260721.csv` (unswapped) for all 9120 wells. The
`_plate1_2_swapped` variant is a barcode/sequencing correction, not a statement about what is
physically in a well, so the unswapped layout is the correct one for a plate-reader target.
Both are run to confirm the choice does not change the conclusion.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE.parents[3] / "shared_pipelines"))

from config import CFG
import genomic_ml as gm
import genomic_ml_yield as gy

OD = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/Karl_2026/"
          "Karl_20260723_OD_Full")
LAYOUTS = {"unswapped (matches echo)": CFG.exp_base / "setup" / "strain_layout_20260721.csv",
           "plate1_2_swapped": CFG.exp_base / "setup" / "strain_layout_20260721_plate1_2_swapped.csv"}
N_PERM = 5
MODELS = ["mean_baseline", "additive_yield_no_genomics", "ridge_pca", "xgboost_pca",
          "two_stage_ridge"]


def permute_genomes(X, summ, seed):
    """Reassign genome feature vectors to well labels. Same genomes, same labels, same pair
    structure -- only the well->genome correspondence is destroyed."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(X.index.values)
    Xp, sp = X.copy(), summ.copy()
    Xp.index, sp.index = perm, perm
    return Xp.loc[X.index], sp.loc[summ.index]


def main():
    rows = []
    for lay_name, lay_path in LAYOUTS.items():
        gcfg = gy.YieldMLConfig(exp_cfg=CFG, layout_csv=lay_path, od_dir=OD,
                                out_dir=HERE / "outputs")
        pairs, summary, wells = gy.build_dataset(gcfg)
        X, summ = gy.strain_feature_matrix(gcfg, pairs)
        rel = gy.replicate_reliability(gcfg, wells)
        print(f"\n=== layout: {lay_name} ===")
        print(f"  {len(pairs)} pairs, {len(set(pairs.strain_a)|set(pairs.strain_b))} strains, "
              f"replicate ceiling rho={rel.spearman.iloc[0]:.3f}")

        for variant in ["true_mapping"] + [f"permuted_mapping_{i}" for i in range(N_PERM)]:
            Xv, sv = ((X, summ) if variant == "true_mapping"
                      else permute_genomes(X, summ, 700 + int(variant.split("_")[-1])))
            s, _, _ = gy.cross_validate(gcfg, pairs, Xv, sv, models=MODELS, n_repeats=2,
                                        shuffle_control=(variant == "true_mapping"))
            for _, r in s.iterrows():
                rows.append({"layout": lay_name, "variant": variant, "regime": r.regime,
                             "model": r.model, "r2_mean": r.r2_mean,
                             "rho_mean": r.spearman_rho_mean, "n_test": r.n_test})
            key = s[(s.model == "ridge_pca")]
            print("  " + f"{variant:22} " + "  ".join(
                f"{r.regime}: R2={r.r2_mean:+.3f} rho={r.spearman_rho_mean:.3f}"
                for _, r in key.iterrows()))

    df = pd.DataFrame(rows)
    df.to_csv(HERE / "outputs" / "n01_naive_yield_with_nulls.csv", index=False)

    print("\n=== VERDICT (ridge_pca, the best cv_strain model on 20260630) ===")
    ver = []
    for lay in df.layout.unique():
        for reg in ("cv_pair", "cv_strain"):
            d = df[(df.layout == lay) & (df.regime == reg) & (df.model == "ridge_pca")]
            t = d[d.variant == "true_mapping"].r2_mean
            p = d[d.variant.str.startswith("permuted")].r2_mean.values
            if not len(t) or len(p) < 2:
                continue
            t = float(t.iloc[0])
            ver.append({"layout": lay, "regime": reg, "r2_true": round(t, 4),
                        "r2_permuted_mean": round(float(p.mean()), 4),
                        "r2_permuted_sd": round(float(p.std(ddof=1)), 4),
                        "z": round(float((t - p.mean()) / p.std(ddof=1)), 2),
                        "beats_all_permutations": bool(t > p.max())})
    v = pd.DataFrame(ver)
    v.to_csv(HERE / "outputs" / "n02_verdict.csv", index=False)
    print(v.to_string(index=False))
    return v


if __name__ == "__main__":
    main()
