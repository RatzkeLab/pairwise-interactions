"""Driver for the total-yield genomic ML analysis of 20260630."""
import sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))                       # config.py
sys.path.insert(0, str(HERE.parents[3] / "shared_pipelines"))

import pandas as pd
from config import CFG
import genomic_ml as gm
import genomic_ml_yield as gy

pd.set_option("display.width", 200)

PLATE_READER = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/"
                    "Karl_2026/Karl_20260704_OD_Full")

gcfg = gy.YieldMLConfig(
    exp_cfg=CFG,
    layout_csv=CFG.exp_base / "setup" / "strain_layout_20260630.csv",
    od_dir=PLATE_READER,
)

t0 = time.time()
pairs, summary, wells = gy.build_dataset(gcfg)
print("=== y01 dataset ==="); print(summary.to_string(index=False))

X, summ = gy.strain_feature_matrix(gcfg, pairs)
print(f"\nfeatures: {X.shape[0]} strains x {X.shape[1]} prevalence-filtered KOs")

rel = gy.replicate_reliability(gcfg, wells)
print("\n=== y02 replicate reliability (the ceiling) ===")
print(rel.round(4).to_string(index=False))

# 16S distances come from the interaction pipeline's own reference, so "taxonomy" here means
# what it means everywhere else in this project. Strains without a consensus are simply absent
# from the matrix and the phylo-only model is scored on the subset that has one.
phylo = None
try:
    gml = gm.GenomicMLConfig(exp_cfg=CFG, out_dir=gcfg.out_dir)
    phylo = gm.phylo_distance_matrix(gml, pairs, strict=False,
                                     out_name="y00_phylo16s_distance_matrix.csv")
    n_all = pairs.strain_a.nunique() + pairs.strain_b.nunique()
    print(f"16S distance matrix: {phylo.shape[0]} of "
          f"{len(set(pairs.strain_a) | set(pairs.strain_b))} strains have a consensus "
          f"-- the phylo control is scored on that subset only")
except Exception as e:
    print(f"16S matrix unavailable ({type(e).__name__}: {e}); phylo control will be skipped")

summ_df, fold_df, pair_df = gy.cross_validate(gcfg, pairs, X, summ, phylo=phylo, n_repeats=3)
print("\n=== y03 cross-validation ===")
print(summ_df[["regime", "model", "n_test", "r2_mean", "r2_sd", "spearman_rho_mean", "mae_mean"]]
      .round(3).to_string(index=False))

ydf, ymet, ko_df = gy.genome_to_yield(gcfg, pairs, X, summ, phylo=phylo)
print("\n=== y04 genome -> per-strain yield contribution (leave-strains-out) ===")
print(ymet.round(3).to_string(index=False))
print(ko_df.head(10).to_string(index=False))

figs = gy.make_all_figures(gcfg, summ_df, pair_df, ydf, rel)
print("\nfigures:", figs)
print(f"done in {time.time()-t0:.0f}s -> {gcfg.out_dir}")
