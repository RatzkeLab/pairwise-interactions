# feature_sweep — does a richer feature representation help? (20260630)

Three related experiments, all on both targets, all reporting `cv_strain`.

| script | question | verdict |
|---|---|---|
| `run_feature_sweep.py` | do other annotation schemes beat KEGG KO? | **no** — 6 schemes, 57–13154 features, all in a narrow band |
| `nested_selection.py` | does the collaborator's panX selection help, run honestly inside each fold? | **no** — ties prevalence-matched random |
| `run_metabolic_features.py` | do constructed competition/complementarity/cross-feeding features help? | **no**, but 23 features reach ~88% of KO's performance |

| output | contents |
|---|---|
| `s01`–`s02` | full sweep, cv_strain ρ by feature set |
| **`s03_gap_over_taxonomy.csv`** | **the table to read** — each set against **its own** 16S control |
| `s04_cv_strain_r2_by_feature_set.csv` | same in R² |
| `s05_*` | nested selection: per-fold results and the genes each fold chose |
| `s06_*` | constructed metabolic features, with and without the cross-feeding block |

**Compare each feature set to its own taxonomy control, never across rows** — the retained strain
set differs per table and moves the baseline (TRAPS §13). `panX_full` has the best absolute
yield ρ purely for that reason; against its own control it buys less than KO.

The selected genes are near-disjoint across folds (Jaccard 0.08–0.14, **zero** in all five) yet
perform identically — they are interchangeable phylogenetic markers.
