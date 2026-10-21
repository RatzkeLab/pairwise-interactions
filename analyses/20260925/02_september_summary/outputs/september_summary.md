| | row | n | R² (vs no winner) | ρ | r | sign acc. |
|---|---|---|---|---|---|---|
| 20260907 ML | best genome model (KEGG KO) (`xgboost_raw_ko`) | 332 | 0.45 ± 0.21 | 0.68 ± 0.14 | 0.68 | 0.81 |
| 20260907 ML | first-pass model (`ridge_pca`) | 332 | 0.37 ± 0.24 | 0.62 ± 0.19 | 0.62 | 0.83 |
| 20260907 ML | genome + preculture OD (not genome-only) (`ridge_pca`) | 332 | 0.50 ± 0.17 | 0.71 ± 0.13 | 0.73 | 0.86 |
| 1 floor | zero baseline (predict no winner) (`zero_baseline`) | 332 | 0.00 ± 0.00 |  |  |  |
| 1 floor | label-shuffled control (`SHUFFLED_ridge_pca`) | 332 | -0.12 ± 0.19 | -0.07 ± 0.30 | -0.06 | 0.48 |
| 2 technical replicates (within 20260907) | one well predicts another well of the same pair | 117 | 0.87 | 0.93 | 0.93 | 0.98 |
| 2 technical replicates (within 20260907) | noise ceiling from replicate spread | 46 | 0.99 |  |  |  |
| 2 technical replicates (within 20260907) | noise ceiling from read depth | 573 | 0.98 |  |  |  |
| 3 biological replicates (20260630 -> 20260907) | same pair, July pair mean predicts September | 53 | 0.89 | 0.87 | 0.96 | 0.97 |
| 3 biological replicates (20260630 -> 20260907) | July per-strain strengths predict September pairs | 64 | 0.89 | 0.91 | 0.95 | 0.95 |
| 3 biological replicates (20260630 -> 20260907) | September strengths predict UNSEEN July pairs | 812 | 0.66 | 0.81 | 0.84 | 0.88 |
| 4 previous experiments | 20260630 (merged demux, current): best genome model (by rho) (`xgboost_raw_ko`) | 855 | 0.12 ± 0.55 | 0.59 ± 0.20 | 0.51 | 0.80 |
| 4 previous experiments | 20260630 (merged demux, current): best genome model (by R^2) (`two_stage_ridge`) | 855 | 0.30 ± 0.25 | 0.58 ± 0.21 | 0.55 | 0.79 |
| 4 previous experiments | 20260630 (merged demux, current): 16S-only (`two_stage_ridge_phylo16s_only`) | 855 | -1.83 ± 4.07 | 0.41 ± 0.36 | 0.39 | 0.69 |
| 4 previous experiments | 20260630 (as published, original demux): best genome model (by rho) (`two_stage_ridge`) | 836 | 0.32 ± 0.23 | 0.58 ± 0.19 | 0.57 | 0.80 |
| 4 previous experiments | 20260630 (as published, original demux): 16S-only (`two_stage_ridge_phylo16s_only`) | 836 | -1.30 ± 3.04 | 0.43 ± 0.36 | 0.38 | 0.72 |
| 4 previous experiments | 20260721 |  |  |  |  |  |
| 5 taxonomy only (same folds) | 16S distance only: two_stage_ridge_phylo16s_only (`two_stage_ridge_phylo16s_only`) | 332 | 0.31 ± 0.13 | 0.61 ± 0.18 | 0.59 | 0.80 |
| 5 taxonomy only (same folds) | 16S distance only: ridge_phylo16s_only (`ridge_phylo16s_only`) | 332 | 0.32 ± 0.20 | 0.60 ± 0.18 | 0.58 | 0.80 |
