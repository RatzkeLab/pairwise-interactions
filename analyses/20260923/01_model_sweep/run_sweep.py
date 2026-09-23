"""The tuned models the first pass skipped, under the same two regimes.

The first pass ran one real model (ridge on PCA'd gene content) and got cv_strain
R^2 = 0.372. That was deliberate: with two baselines beside it, the number means
something. This adds the rest of the registry so the question becomes "does anything
beat the simple thing", which is only worth asking once the simple thing is pinned down.

`cv_strain` remains the result. Under `cv_pair` a model can score well by memorising a
strain's competitiveness from its other pairs -- the first pass measured exactly that:
the no-genomics strength model reached 0.641 on held-out pairs and exactly 0.000 on
held-out strains.

16S-phylogeny models are included as a control on the genomic ones. If gene content is
doing real work it should beat "how closely related are these two strains", which is
available for free and carries no metabolic information.
"""
import sys
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE)); sys.path.insert(0, str(BASE.parents[1] / "scripts"))

import numpy as np, pandas as pd
import config, genomic_ml as gml

OUT = BASE / "01_model_sweep" / "outputs"; OUT.mkdir(parents=True, exist_ok=True)

MODELS = [
    "zero_baseline",                    # floor
    "strength_observed_no_genomics",    # the memorisation yardstick
    "ridge_pca",                        # first-pass result, for continuity
    "xgboost_pca",
    "xgboost_raw_ko",
    "two_stage_ridge",
    "two_stage_xgboost",
    "ridge_phylo16s_only",              # 16S-only controls
    "two_stage_ridge_phylo16s_only",
    "two_stage_ridge_ko_plus_phylo",
]
COLS = ["regime", "model", "n_folds", "n_test", "r2_mean", "r2_sd",
        "spearman_rho_mean", "sign_accuracy_mean"]


def main(n_repeats=3):
    cfg = config.make_config()
    gcfg = gml.GenomicMLConfig(exp_cfg=cfg, out_dir=OUT)

    pairs, summary = gml.build_dataset(gcfg)
    print(summary.to_string(index=False))
    X, summ = gml.strain_feature_matrix(gcfg, pairs)
    print(f"\n{X.shape[0]} strains x {X.shape[1]} KO features, {len(pairs)} pairs")

    # 16S distances from this experiment's own consensus -- the same sequences the
    # interaction pipeline mapped against, so the control is exactly comparable
    phylo = gml.phylo_distance_matrix(gcfg, pairs)
    print(f"phylo matrix: {phylo.shape}")

    summ_df, fold_df, pair_df = gml.cross_validate(
        gcfg, pairs, X, summ, models=MODELS, n_repeats=n_repeats,
        phylo=phylo, file_prefix="s01")

    for regime in ("cv_strain", "cv_pair"):
        print("\n" + "=" * 78)
        print(f"{regime}" + ("   <- the result" if regime == "cv_strain"
                             else "   <- inflated by memorisation, for comparison only"))
        print("=" * 78)
        t = summ_df[summ_df.regime == regime].sort_values("r2_mean", ascending=False)
        print(t[COLS].to_string(index=False))

    ceiling = gml.label_noise_ceiling(gcfg, pairs)
    print("\nlabel-noise ceiling:")
    print(ceiling.to_string(index=False))

    # the two-stage model's first half on its own terms: predict a strain's competitiveness
    # from its genome, with no pair structure to inflate anything
    strength, met, ko_df = gml.genome_to_strength(gcfg, pairs, X, summ, phylo=phylo)
    print("\ncan genomics reproduce the competitive hierarchy at all?")
    print(met.to_string(index=False))

    return pairs, X, summ, summ_df, fold_df, pair_df, ceiling, (strength, met, ko_df)


if __name__ == "__main__":
    main()
