"""Driver for the genomic-ML analysis of 20260630. Mirrors the notebook; run either."""
import sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))                       # config.py
sys.path.insert(0, str(HERE.parents[3] / "shared_pipelines"))  # shared library
import pandas as pd
from config import CFG
import genomic_ml as gm

pd.set_option("display.width", 200)
gcfg = gm.GenomicMLConfig(exp_cfg=CFG)

t0 = time.time()
val, ok = gm.validate_strain_join(gcfg)
print("=== g00 strain-join validation ==="); print(val.to_string(index=False))
print("VERDICT:", "PASS" if ok else "FAIL")
assert ok, "strain join failed validation -- refusing to build a modeling set on a bad join"

pairs, summary = gm.build_dataset(gcfg)
print("\n=== g01 dataset ==="); print(summary.to_string(index=False))

X, summ = gm.strain_feature_matrix(gcfg, pairs)
print(f"\nfeature matrix: {X.shape[0]} strains x {X.shape[1]} prevalence-filtered KOs")

phylo = gm.phylo_distance_matrix(gcfg, pairs)
print(f"16S distance matrix: {phylo.shape[0]} strains")

ceiling = gm.label_noise_ceiling(gcfg, pairs)
print("\n=== g02 label-noise ceiling ==="); print(ceiling.to_string(index=False))

_pos = [a for a in sys.argv[1:] if not a.startswith("--")]
models = _pos[0].split(",") if _pos else [m for m in gm.MODELS if "tabpfn" not in m]
summ_df, fold_df, pair_df = gm.cross_validate(gcfg, pairs, X, summ, models=models,
                                              phylo=phylo, n_repeats=3)
print("\n=== g03 cross-validation ===")
print(summ_df[["regime", "model", "n_test", "r2_mean", "r2_sd", "spearman_rho_mean",
               "sign_accuracy_mean", "mae_mean"]].round(3).to_string(index=False))

# TabPFN runs as its own pass: it is the one step that sends data off-machine, and its
# hosted API is slow enough that the repeated folds used above are not worth the round trips.
if "--tabpfn" in sys.argv:
    tf_summ, tf_fold, tf_pair = gm.cross_validate(
        gcfg, pairs, X, summ, models=["tabpfn", "two_stage_tabpfn"],
        phylo=phylo, n_repeats=1, shuffle_control=False, file_prefix="g03tabpfn")
    fold_df = pd.concat([fold_df, tf_fold], ignore_index=True)
    pair_df = pd.concat([pair_df, tf_pair], ignore_index=True)
    summ_df = gm.summarize_folds(fold_df)
    summ_df.to_csv(gcfg.out_dir / "g03_cv_summary.csv", index=False)
    print("\n=== g03 cross-validation, with TabPFN ===")
    print(summ_df[["regime", "model", "n_test", "r2_mean", "r2_sd", "spearman_rho_mean",
                   "sign_accuracy_mean", "mae_mean"]].round(3).to_string(index=False))

sdf, smet, ko_df = gm.genome_to_strength(gcfg, pairs, X, summ, phylo=phylo)
print("\n=== g04 genome -> per-strain competitiveness (leave-strains-out) ===")
print(smet.round(3).to_string(index=False))
print(ko_df.head(10).to_string(index=False))

figs = gm.make_all_figures(gcfg, summ_df, fold_df, pair_df, sdf, ceiling)
print("\nfigures:", figs)
print(f"done in {time.time()-t0:.0f}s -> {gcfg.out_dir}")
