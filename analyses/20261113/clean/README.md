# Genomic ML, re-checked: leakage and wrong genomes (2026-11-13)

Experiments 20260907 and 20260630; the published best genome models, rerun with identical
code and data (the rerun reproduces every published fold exactly). Details, controls and
code: `../messy/README.md`.

## Answer

1. **Leakage was real, and it flattered the best model most.** Splitting by well let the same
   organism sit in training and test (70 % of test pairs in 20260907, 58 % in 20260630). With
   the leak removed — same test pairs, and a control for the data the removal costs —
   XGBoost-on-KO loses **0.08 ρ / 0.12 R²** in 20260907 and **0.06 ρ / 0.10 R²** in 20260630.
   The linear models lose ~0.04 (20260907) or nothing (20260630). A lookup that knows no
   biology, only "have I seen this genome", scored ρ 0.41 in the published 20260907 folds.
   The published claim that gene content beats 16S in 20260907 (+0.055 ρ) was the leak.
2. **Removing wells whose genome is wrong helps a lot — more than the leak cost.** Scored on
   the same correct-genome test pairs, leak-free, the best model goes from ρ 0.60 → **0.78**
   (20260907) and 0.61 → **0.81** (20260630), winner called 80 → 90 % and 84 → 93 %, while
   training on a third fewer pairs. Removing the same amount at random makes it slightly
   worse. 16S-only models gain little, as expected: 16S is measured from what is in the well.

## Headline numbers — best model (XGBoost on raw KO), Spearman ρ

| | 20260907 | 20260630 |
|---|---|---|
| A. published (folds by well, all wells) | 0.67 | 0.59 |
| B. same test pairs, leak removed | 0.58 | 0.53 |
| &nbsp;&nbsp;&nbsp;control: same data loss, leak kept | 0.66 | 0.59 |
| C. leak-free, all wells, scored on correct-genome pairs | 0.60 | 0.61 |
| D. leak-free, **wrong-genome wells removed**, same pairs as C | **0.78** | **0.81** |
| D, median R² (vs predicting no winner) / winner called | 0.69 / 90 % | 0.71 / 93 % |

Means over 100 folds (20 repeats × 5). Leak = control − B. Two-stage ridge and the 16S
controls are in `headline_numbers.csv`. A and B/C/D are scored on different pair sets, so
compare A↔B and C↔D, not A↔D.

![step 1](f01_leakage.png)

![step 2](f02_wrong_genomes.png)

## What changes

- **Split by organism from now on.** The per-well organism unit is in
  `../messy/outputs/s01_wells_<exp>.csv` (`cv_group`); `../messy/lib.py:group_folds` does the
  split. 280 wells in 20260907 are only 175 distinct genomes.
- **Drop wrong-genome wells before modelling**, using each experiment's own verdict (the plate
  drifts: 20260630 needed its own check, `../messy/s01_organisms.py`).
- **Gene content vs taxonomy is a near tie once both problems are fixed** (20260907: ρ 0.78 vs
  16S 0.76, +0.02, p 0.005; winner call 90 vs 91 %, n.s.). Don't quote the earlier margins.
- **The 16S control is numerically unstable in 20260630** (PCA whitening; R² < −2 in up to
  half the folds, published run included), so no KO-vs-16S claim for 20260630 stands on this
  pipeline. 20260907 is unaffected. Not fixed — it lives in the shared `genomic_ml.py`.
- **Not rerun:** plate-reader models, yield model, genome→strength. They use the same folds,
  so their published numbers are optimistic too.

## Strains, by strain number

`strain_ids.csv` — well → lab-collection strain number. The 384 sheet numbers freezer
*stocks*: B14 (`M7`) and K19 (`MYb9`) are the same isolate. Here both are **M7**. 384 wells =
372 isolates; agrees with the genome mapping on all 298 genome wells. (Master collection sheet
copied to `/home/rl/data/resources/karl/`.)

`reculture_list.csv` — the 61 wells whose genome is wrong in the current (September) plate,
plus 13 unverified, with strain number, freezer box/position, what the well holds now, and
whether it was already wrong in July (11 were). Commonest invader: C001 (the G01 organism) in
9 wells. Cluster members by strain number: `clone_clusters.csv`.
