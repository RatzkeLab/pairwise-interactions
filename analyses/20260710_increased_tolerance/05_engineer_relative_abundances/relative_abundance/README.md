# relative_abundance — 20260630

**Question:** when two strains share a well, how do they split it — and is the resulting network
a dominance hierarchy?

Run: `relative_abundance_20260630.ipynb`. Logic in `shared_pipelines/relative_abundance.py`
and `hierarchy_significance.py`.

| output | contents |
|---|---|
| `r01_reference_pair_distances.csv` | 16S bp distance between each tested pair's references |
| `r02_well_interaction_scores.csv`, `r02_read_assignments.csv.gz` | per-well and per-read scoring |
| **`r03_pair_replicate_stats.csv`** | **one row per pair, replicate-averaged — the ML target** |
| `r04_strain_competitiveness.csv`, `r05_bt_strengths.csv` | per-strain competitiveness, Bradley-Terry fit |
| `r05_hierarchy_summary.csv` | pseudo-R² 0.83, DCI, intransitive triads |
| `r07`–`r11` | permutation, bootstrap, goodness-of-fit and sensitivity suites |

**Two columns to respect in `r03`:** `high_uncertainty_pair` marks pairs whose references are too
close to tell apart — the ~50/50 there is the assay's resolution limit, not coexistence, and
training on them teaches an artifact. `n_replicates` is mostly 1.

**Regenerated 2026-08-27** under the `MIN_RESOLVABLE_BP = 10` fix: usable pairs rose from 996 to
1480. Anything computed before that date is stale.
