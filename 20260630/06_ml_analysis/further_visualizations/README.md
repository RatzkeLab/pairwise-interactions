# further_visualizations

Side plots for models the main `genomic_ml` report does not show. `make_all_figures` plots only
the winner of each regime (f02), which is the right default for a report -- this folder is where
the runners-up get looked at without cluttering it.

Nothing here refits. Everything reads `../genomic_ml/outputs/g03_cv_*.csv`, so these figures can
never disagree with the run that produced the report, and they are instant.

## plot_model_scatter.py

Predicted-vs-observed scatter for any model x regime present in `g03_cv_predictions.csv`.

```bash
python plot_model_scatter.py                                   # xgboost_raw_ko + two_stage_ridge, cv_strain
python plot_model_scatter.py --models xgboost_pca --regime cv_pair
python plot_model_scatter.py --models all
```

Each panel carries BOTH scorings, because they are different numbers and the difference has
already caused confusion once:

- **pooled** (panel title) -- one score over all held-out predictions from all folds at once.
  What the scatter itself shows.
- **mean-of-folds** (inset box) -- scored per fold, then averaged. What `summ_df` and the report
  quote, and the only one carrying a fold SD.

For cv_strain / two_stage_ridge these are ρ 0.596 vs 0.577. Neither is wrong; averaging a bounded
nonlinear statistic is simply not the same as pooling. **Quote the mean-of-folds in text.** The
pooled value is also mildly optimistic in a way the scatter hides: it merges predictions from 15
separately fit models, each with its own stage-1 mean-zero centering, into one ranking.

R² follows `genomic_ml`'s convention throughout: **vs. predicting no winner (0), not vs. the
mean**, so it is comparable to f01 and to the summary table, but *not* to `genomic_ml_yield`,
which scores against the mean.

## Current finding: xgboost_raw_ko vs two_stage_ridge under cv_strain

`fv01_predicted_vs_observed_cv_strain_xgboost_raw_ko_two_stage_ridge.png`

The two are tied on rank correlation (mean-of-folds ρ 0.574 vs 0.576, and paired across folds
p = 0.93) but not on R² (0.197 vs 0.321). The scatter says why, and it is **calibration, not
ranking**:

| | predicted sd | predicted range | slope of obs ~ pred |
|---|---|---|---|
| observed | 3.73 | -7.2 .. 7.2 | -- |
| xgboost_raw_ko | 2.86 | -6.7 .. 6.6 | 0.71 |
| two_stage_ridge | 1.83 | -4.7 .. 4.9 | 1.20 |

xgboost_raw_ko commits to large magnitudes (slope < 1 = it overshoots), two_stage_ridge shrinks
toward zero (slope > 1 = it undershoots). Because R² here is scored **against predicting no
winner**, an overshoot on a pair you ranked correctly-but-not-confidently is expensive while a
timid prediction costs almost nothing -- so ridge's alpha=50 shrinkage is rewarded by this metric
independently of whether it understands the genomes any better. Hence identical ρ and a 0.12 R²
gap that does not survive a paired test (p = 0.42, 7/15 folds).

Read together with the ρ tie, that is a reason to be *more* cautious about the R² ranking in f01,
not less. The claim already in the report still holds: model choice is not the bottleneck, genome
count is.
