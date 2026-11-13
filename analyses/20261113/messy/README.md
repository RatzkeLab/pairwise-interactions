# 20261113 (messy) — genomic ML: leakage, and wells whose genome is wrong

Working record. The short version is in `../clean/README.md`.

Asked (`prompt.txt`): (1) rerun the best genomic ML models with and without the train/test
leakage that came from splitting by well instead of by organism; (2) then train and test with
and without the wells whose genome on file is not the organism in the well (overgrowth /
cross-contamination found in `analyses/20261021/03_genome_concordance`). Also: refer to
strains by their lab-collection strain number, not their 384-well position.

## How to run

Env `karl_seq_analysis`. Every script imports `config` first (it caps BLAS at one thread per
process — see "Gotchas").

| script | what | time |
|---|---|---|
| `s00_strain_ids.py` | well → isolate strain number (reconciles three sources) | s |
| `s01_organisms.py` | 20260630's own genome verdicts (BLAST vs the 298 assemblies, 20261021 rules); per-well organism / CV-group tables for both experiments | 20 s |
| `s02_leakage.py 20` | step 1: published folds vs organism folds, and the purge arm | ~15 min* |
| `s02b_purge_control.py 20` | step 1 control: same data loss as the purge, leak kept | ~2 min |
| `s03_wrong_genomes.py 20` | step 2: train with / without wrong-genome wells (organism folds) | ~8 min |
| `s03_wrong_genomes.py 20 --union` | step 2 sensitivity on union folds | ~5 min |
| `s05_tables.py` | all final tables from saved fold metrics (no refitting) | s |
| `s04_figures.py` | f01, f02 | s |
| `s06_reculture_list.py` | wrong / unverified wells by strain number and freezer box | s |

\*with another 60-process job on the machine. s02b/s03 skip refitting if their fold metrics
exist; pass `--rerun` to force it.

`lib.py` is the shared CV driver. It reuses `scripts/genomic_ml.py` **unchanged** for the
dataset, features, per-fold context (in-fold PCA, antisymmetric design), the model zoo and
`_metrics`. Only the fold assignment is new (`group_folds`), plus one diagnostic model.

**Reproduction check (passes).** `group_folds` with identity groups is `genomic_ml.make_folds`;
reps 0–2 of `by_well` reproduce every published per-fold row (20260923 sweep for 20260907,
`20260710_increased_tolerance/.../outputs_merged_ab` for 20260630): 240 fold × model rows,
max |diff| 3e-15 / 1e-12. The modelling sets are byte-identical to the published
`g01_modeling_pairs.csv` (asserted in `lib.load`).

## Models rerun

The best genome models of the previous runs and their controls, all under every scheme:

| model | why |
|---|---|
| `xgboost_raw_ko` | best by ρ in both experiments (20260907 ρ 0.68, 20260630 ρ 0.59) |
| `two_stage_ridge` | best by R² in 20260630 (0.30) |
| `ridge_pca` | first-pass model, continuity |
| `two_stage_ridge_phylo16s_only`, `ridge_phylo16s_only` | 16S-only controls |
| `clone_lookup_no_genomics` | **new.** No gene content: a test well gets the in-fold strength of training wells carrying the *same* genome (≥99.9 % ANI clone cluster), else 0. Everything it scores is leakage by construction |
| `zero_baseline`, `SHUFFLED_*` | floor, label-shuffled control |

Not rerun: the plate-reader models (`genomic_ml_plate`, preculture OD), the yield model and
`genome_to_strength`. The same leak applies to them (same fold code), so treat their
published numbers as optimistic too.

**Metrics.** Spearman ρ is primary (rank-based; the published "best model" was picked by ρ).
R² is vs predicting "no winner" (`genomic_ml._metrics`) and is reported as the **median** over
folds, because a few folds per 100 blow up in 20260630 (see "16S control"). Every number is a
mean over 100 folds (20 repeats × 5 folds, published used 3 × 5); "sd" = SD of the per-repeat
mean across the 20 repeats.

## s00 — strain numbers

`outputs/s00_strain_ids.csv`: `well, strain_no, isolate_name, stock_no, box, box_position, ...`

The 384 sheet's `Strain.Number` is a **freezer stock**, not an isolate: B14 = `M7` and K19 =
`MYb9` are two stocks of MYb9 (master sheet: M7 in WS5/7 and a row numbered `MYb9` in WS4/56).
Rule: isolate = the stock's alternative name in the master sheet if it has one, else the
stock number; `strain_no` = the first M-number of that isolate (Or isolates keep their
integer). 384 wells → **372 isolates**; agrees with `mapping_384_well_plate_collection.csv`'s
box ids on all 298 genome wells. 9 isolates sit in 2–3 wells (MYb9/10/11/14/53/56/57/71/181).

No other well↔strain file exists on disk (searched `karl/`, `/home/rl/data`, the genome
folder). The master sheet (`plate_format/Strain collection including Schulenberg's -
collection.csv`) was copied to `/home/rl/data/resources/karl/` on request.

Side finding: for all 9 duplicated isolates only the K/M-row copy has a genome joined in the
mapping; the B/F/H/J-row copies hold the same isolate and get none. None are in the published
modelling sets.

## s01 — what is in each well

| | 20260907 | 20260630 |
|---|---|---|
| modelled wells | 280 | 74 |
| distinct genomes (clone clusters) | **175** | **52** |
| wells sharing a genome with another modelled well | 148 | 31 |
| genome `ok` / `wrong` / `unverified` | 211 / **57** / 12 | 49 / **13** / 12 |
| verdict source | 20261021 concordance (Sept. plate) | **own**: July references vs assemblies |

20260630 needed its own verdict: where July's and September's 16S agree (≤5 bp) the verdicts
agree perfectly (52 ok/ok, 7 wrong/wrong); where they differ, July has its own calls (4
wrong-in-July that are ok in September, 7 unverified). `outputs/s01_genome_status_20260630.csv`.
Caveat: some July references are coin-flip consensus calls (20260925), so a July `wrong`
may be a bad reference rather than a bad genome — the well's data is unreliable either way.

**CV groups** (`cv_group`, the leak-free unit) = the organism physically in the well:
`ok`/`unverified` wells → their genome's clone cluster; `wrong` wells → the invader's cluster
when one genome unit carries the well's 16S (the 9 C001-overgrown wells join C001), else a
group keyed by the carrier set, else the well alone. A wrong well is *not* tied to its own
genome's cluster (there it shares the KO vector but not the phenotype — that costs the model,
it does not flatter it). Alternatives kept: `cv_group_genome` (clone cluster only, what the
20261021 README prescribed) and `cv_group_union` (both, connected; largest block 37 wells).
20260907: 181 / 175 / 156 groups; 20260630: 54 / 52 / 46.

## Step 1 — leakage (`s02`, `s02b`)

How exposed was the published CV? In the published folds, **70 %** (20260907) and **58 %**
(20260630) of test pairs had at least one strain whose organism was also in training.

Design. The cleanest measure keeps the test pairs fixed: each published fold, with every
training pair involving an organism present in the test fold removed (`by_well_purged`).
That removes about a third of training (366 → 244 pairs; 953 → 609), so `s02b` adds the
control `by_well_random`: the same number of training pairs removed at random among other
organisms (leak kept, same data loss). **Leak = random − purged.**

Spearman ρ (paired over 100 folds):

| exp | model | published | same loss, leak kept | leak removed | **leak Δρ** | leak ΔR² |
|---|---|---|---|---|---|---|
| 20260907 | xgboost_raw_ko | 0.669 | 0.655 | 0.576 | **0.078 ± 0.015** (p 2e-6) | 0.116 ± 0.020 |
| 20260907 | two_stage_ridge | 0.635 | 0.622 | 0.582 | 0.040 ± 0.013 (p 0.015) | 0.033 ± 0.010 |
| 20260907 | ridge_pca | 0.653 | 0.638 | 0.597 | 0.040 ± 0.013 (p 0.012) | 0.053 ± 0.016 |
| 20260907 | 16S two-stage | 0.613 | 0.606 | 0.590 | 0.016 ± 0.010 (n.s.) | 0.013 |
| 20260907 | 16S ridge | 0.609 | 0.596 | 0.594 | 0.002 (n.s.) | 0.013 |
| 20260907 | clone lookup | 0.405 | 0.387 | −0.005 | 0.392 | 0.305 |
| 20260630 | xgboost_raw_ko | 0.594 | 0.594 | 0.532 | **0.062 ± 0.019** (p 0.002) | 0.104 ± 0.047 |
| 20260630 | two_stage_ridge | 0.603 | 0.593 | 0.606 | −0.013 (n.s.) | −0.018 |
| 20260630 | ridge_pca | 0.581 | 0.564 | 0.577 | −0.013 (n.s.) | — |
| 20260630 | clone lookup | 0.184 | 0.190 | −0.008 | 0.190 | 0.103 |

Reading:
- **The leak was real and concentrated in the flexible model.** XGBoost on raw KO can memorise
  a genome as a fingerprint: in 20260907 it lost 0.08 ρ / 0.12 R² when the leak was removed —
  most of its published lead. The linear PCA models lost ~0.04 in 20260907 and nothing in
  20260630; the 16S controls nothing measurable.
- **A model with no genomics at all** (`clone_lookup`) scored ρ 0.41 in 20260907's published
  folds just by recognising genomes it had seen. 20260907 is far more exposed than 20260630:
  ~4 pairs per strain (so a clone-mate's pairs are a large share of what is known) vs ~40.
- **The published "gene content beats 16S" margin in 20260907 was the leak.** Published:
  KO 0.669 vs 16S 0.613 (+0.055, p 4e-4). Leak removed: 0.576 vs 0.590 (−0.014, n.s.).

Full leak-free CV (independent fold partitions, not paired; `t1_leakage.csv`), xgboost_raw_ko ρ:
20260907 by_organism 0.512 / by_genome 0.596 / union 0.596; 20260630 0.482 / 0.596 / 0.609.
`by_organism` is lowest because it can test a wrong-genome well while its genome's genuine
clone-mates are in training — the model then confidently predicts the genome's owner for a
well that holds something else. That is step 2's problem leaking into step 1, which is why
the paired purge is the headline for step 1.

## Step 2 — wrong genomes (`s03`)

Leak-free organism folds. Removing wells changes both training and test sets, so within each
fold the **test pairs are fixed** (both wells' genomes not `wrong`) and only training varies:
all wells / wrong-genome wells removed / the same number of training pairs removed at random
among the other wells. Removing wrong wells drops ~1/3 of training (367 → 241; 951 → 650).

Spearman ρ on correct-genome test pairs:

| exp | model | trained on all | **wrong removed** | random removed | Δρ | median R² all → removed | winner called |
|---|---|---|---|---|---|---|---|
| 20260907 | xgboost_raw_ko | 0.598 | **0.781** | 0.575 | +0.183 ± 0.028 | 0.48 → 0.69 | 80 → 90 % |
| 20260907 | two_stage_ridge | 0.649 | 0.759 | 0.618 | +0.110 ± 0.022 | 0.50 → 0.56 | 85 → 91 % |
| 20260907 | ridge_pca | 0.647 | 0.764 | 0.612 | +0.116 ± 0.022 | 0.61 → 0.69 | 84 → 91 % |
| 20260907 | 16S two-stage | 0.723 | 0.757 | 0.668 | +0.034 ± 0.014 | 0.44 → 0.51 | 89 → 91 % |
| 20260907 | 16S ridge | 0.716 | 0.759 | 0.666 | +0.043 ± 0.014 | 0.60 → 0.68 | 90 → 91 % |
| 20260630 | xgboost_raw_ko | 0.610 | **0.807** | 0.569 | +0.197 ± 0.025 | 0.41 → 0.71 | 84 → 93 % |
| 20260630 | two_stage_ridge | 0.684 | 0.792 | 0.633 | +0.108 ± 0.018 | 0.53 → 0.61 | 88 → 92 % |
| 20260630 | ridge_pca | 0.662 | 0.818 | 0.597 | +0.156 ± 0.021 | 0.51 → 0.67 | 87 → 94 % |

All p < 1e-6 for the KO models. Random removal of the same amount *lowers* ρ by 0.02–0.07, so
the gain is the noise going away, not less data. Stricter `ok_only` (also drop `unverified`):
same picture (20260907 xgboost 0.577 → 0.802; 20260630 0.551 → 0.768).

- **The KO models gain roughly 3–5× more than the 16S models** — the signature of wrong genomes: 16S
  is measured from what is in the well, the KO vector is looked up for what *should* be.
- Pairs involving a wrong well were poorly predicted all along (xgboost ρ 0.30–0.34).
- **Sensitivity, union folds** (a wrong well is held out together with its genome's clone-mates,
  so it can never teach "genome X behaves like the invader" to a test clone of X): gains are
  smaller but every one is still significant — xgboost +0.042 (0.729 → 0.772) in 20260907,
  +0.100 (0.693 → 0.793) in 20260630; two-stage ridge +0.039 / +0.052. The *cleaned* models
  land at the same place in both groupings (ρ 0.77–0.81); only the all-wells baseline moves.
  So most of the damage runs through clone-mates. In real use (predicting a strain whose
  genome matches a contaminated well) that path is open, so organism folds are the realistic
  estimate and union folds the floor.
- **Gene content vs 16S after cleaning (20260907, fold-paired):** xgboost 0.781 vs 16S
  two-stage 0.757, +0.024 ± 0.013 (p 0.005); winner accuracy 90.1 vs 91.1 % (n.s.). Gene
  content is ahead, but barely — not the published margin. Before cleaning, on the same test
  pairs, 16S *beat* KO by 0.125 ρ: the wrong genomes, not taxonomy, were the KO model's
  problem.

## The 16S control in 20260630 is numerically unstable (pre-existing)

In 20260630, `ridge_phylo16s_only` / `two_stage_ridge_phylo16s_only` have R² < −2 in 10–54 of
100 folds under every scheme — including the published one (whose 16S R² is −6.85 ± 13.1).
Cause, checked on the worst published fold: `_Ctx` keeps 30 landmark PCs and whitens each to
unit variance. The last ones carry ~1 bp² of variance (reference noise); after scaling, a
held-out strain sits at up to **76 SD** on such a component (training max 6.4) and the ridge
extrapolates. KO PCA models are hit in ≤5/100 folds; 20260907 in none. Consequences: 20260630's
16S ρ is depressed, and **no KO-vs-16S statement for 20260630 should be taken from this
pipeline** (including the earlier "KO R² 0.41 vs 16S 0.11" from `genome_to_strength`, which
whitens the same way). Fix (not applied — it would change published numbers): keep components
by explained-variance fraction, or do not whiten. Not fixed here; `scripts/genomic_ml.py` is
untouched.

## Re-culture list (`s06`)

`outputs/s06_reculture_list.csv`: all 61 wells whose genome is wrong in the September source
plate (plus 13 unverified), by `strain_no`, isolate, freezer box and position, what the well
holds now (`now_holds`, strain numbers or a clone cluster; clusters resolved in
`outputs/s06_clone_clusters.csv`), and July's verdict where July had a reference. 11 were
already wrong in July. Commonest invaders: C001 (the G01 organism, 16 Schulenberg stocks
M87…M128) in 9 wells, C002 (+ strain 139) in 5, C034 in 4.

## Gotchas hit on the way

- `scripts/paths.py` is stale after the 2026-11-13 move (`PW_ROOT`, `ANALYSES`, `EXPERIMENTS`,
  `CONSENSUS_MERGE` point at the old tree; `paths.check()` → 7 missing). `config.py` builds
  its own paths; `genomic_ml` still works because it needs only `GENOMIC_TABLES`/`MAPPING_DIR`
  and is imported via the package path (`karl/` must be on `sys.path`).
- numpy's BLAS spawns a thread per core in every worker: 14 workers ran ~1000 threads and a
  15-min job stalled for 30+ min. `config` sets OMP/OPENBLAS/MKL threads = 1; xgboost keeps
  `genomic_ml`'s 4 so published folds reproduce exactly.
- `pivot_table` silently drops an arm whose metric is NaN in every fold (constant predictors);
  reindex the columns.
- Sorting `exp` read back from CSV gives ints; cast to str before matching.

## Outputs

`outputs/tables/` — `t3_headline.csv`, `t1_leakage.csv`, `t2_wrong_genomes.csv`,
`t4_ko_vs_16s.csv`, `tables.md` (all as markdown). `outputs/figures/` — `f01_leakage.png`,
`f02_wrong_genomes.png`. Per-fold metrics and per-pair predictions in `outputs/s02_leakage/`,
`outputs/s03_wrong_genomes/`, `outputs/s03_wrong_genomes_union/`.
