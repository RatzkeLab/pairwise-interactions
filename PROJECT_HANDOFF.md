# Handoff: pairwise_interaction_experiments

Written 2026-08-27 to hand off from a read-mapping/statistics session to a new session
focused on a different direction: predicting pairwise strain interactions from genomic data
via machine learning. Read this first; it points to the specific files and gotchas rather
than re-deriving them.

## What this project is

ONT full-length 16S amplicon sequencing of bacterial strains grown pairwise (and as mono
controls) in 384-well plates, across (so far) two experiments: `20260721` and `20260630`.
For each well, we know which two strains (`strain1`, `strain2`) were supposed to be there
(`setup/strain_layout_*.csv`) and have the raw demultiplexed reads. The analysis built so far
answers: does each well contain what it's supposed to (QC), what's the relative abundance of
each strain when both are present, is there a stable per-strain "competitiveness," and is the
resulting network of pairwise outcomes a real dominance hierarchy or a tangled one.

## Code layout

```
pairwise_interaction_experiments/
├── shared_pipelines/                     -- all actual logic, shared across experiments
│   ├── experiment_config.py              ExperimentConfig dataclass
│   ├── io_utils.py                       generic read/reference IO helpers
│   ├── mapping_validation.py             read-mapping QC (does a well contain what it should)
│   ├── relative_abundance.py             per-well interaction scoring, replicate stability,
│   │                                      per-strain competitiveness, Bradley-Terry hierarchy fit
│   └── hierarchy_significance.py         5-test significance suite on the hierarchy claim
├── 20260721/analysis/
│   ├── config.py                         this experiment's ExperimentConfig (paths, reference fasta)
│   ├── mapping_validation/{outputs/, mapping_validation_20260721.ipynb}
│   └── relative_abundance/{outputs/, relative_abundance_20260721.ipynb}
└── 20260630/analysis/                    same shape, different config.py
```

No per-experiment `scripts/` folders -- notebooks import functions from `shared_pipelines`
directly (`sys.path.insert` to `analysis/` for `config.py` and to `shared_pipelines/`). Each
notebook re-runs its whole pipeline top to bottom and is the human-readable report.

## The outputs that matter for the NEW (genomic ML) task

The **prediction target** lives in `relative_abundance/outputs/`, per experiment:

- `r03_pair_replicate_stats.csv` -- **one row per unique tested (strain_a, strain_b) pair**,
  replicate-averaged. Columns of interest: `mean_relative_abundance_a` (0-1, the label),
  `mean_log2_ratio_a_over_b` (symmetric, unbounded -- probably the better ML target than the
  bounded ratio), `mean_uncertainty_score`, `high_uncertainty_pair` (bool -- **near-identical
  16S reference pairs where the label is a forced ~50/50, not a real biological measurement;
  seriously consider excluding these rows or down-weighting them**), `n_replicates`,
  `ref_pair_bp_dist`. This is probably the file to start from.
- `r02_well_interaction_scores.csv` -- the same thing at well (not pair) granularity, if you
  want per-replicate rather than averaged labels.
- `r05_pairwise_relative_abundance_matrix.csv` -- full strain x strain matrix (same info as
  r03, wide format).
- `r05_bt_strengths.csv` -- per-strain Bradley-Terry "strength" -- a plausible auxiliary
  target or sanity-check feature (does genomic similarity predict BT strength too?).
- `r04_strain_competitiveness.csv` -- per-strain competitiveness score, for sanity-checking
  any per-strain genomic feature you engineer (e.g. does genome size/KO count correlate?).

**Only a subset of possible pairs were ever tested** (~58-63% for 20260721, less for
20260630) -- that defines the labeled dataset size. Not every strain has a resolvable
consensus reference either (e.g. `N13` in 20260721 has none at all) -- check
`analysis/consensus.../` coverage before assuming every strain in the layout has both a label
and genomic features available.

## The genomic feature data

`/home/rl/scripts/karl/Link to Karl/final_genomic_tables/`:

- `KEGG_ko_and_strains_table.csv` (298 strains x ~4500 KO columns, counts not just
  presence/absence) -- the file the user wants to start with.
- `KEGG_Module_and_strains_table.csv`, `CAZy_and_strains_table.csv`, `BiGG_and_strains_table.csv`,
  `PFAMs_and_strains_table.csv`, `panx_and_strains_table.csv` -- same strain axis, different
  functional annotation schemes; potential alternative/additional feature sets.
- `mapping_384_well_plate_collection.csv` -- **likely the bridge between genomic-table strain
  IDs and the pairwise-experiment well-coordinate strain labels** (see below) -- not yet fully
  verified, start here.

### ⚠️ The strain-identifier join is the first real problem to solve

The genomic tables key strains as e.g. `BIGb0170`, `CEent1`, `JUb134` -- **not** the short
well-coordinate labels (`D13`, `A11`, `O18`, ...) used as `strain1`/`strain2` in
`setup/strain_layout_*.csv` throughout this project. They are not the same namespace and must
not be assumed interchangeable.

`mapping_384_well_plate_collection.csv` has columns `strain, sequencing_batch,
strain_index_name(as_in_boxes), Well_souce_plate, assembly_name` and contains **three
different `sequencing_batch` values** (`Or`: 206 rows, `Schullenberg`: 63 rows,
`Schullenberg_CeMbio43`: 29 rows) with **different `strain` column conventions per batch**:
for `Schullenberg_CeMbio43`, `strain == assembly_name == "BIGb0170"`-style already; for `Or`,
`strain` is a bare integer (e.g. `172`) that does **not** correspond to the numeric suffix of
a same-numbered `BIGb` strain from a different batch (checked: `strain=172` under `Or` maps
to well `M9`, assembly `P3_A3` -- a completely different, unrelated strain from `BIGb0172`
under `Schullenberg_CeMbio43`, well `O16`). **Do not join on the numeric suffix alone.**

The `Well_souce_plate` column (e.g. `O18`, `M9`, `D13`) uses the same well-coordinate style as
this project's `strain1`/`strain2` labels, which is the promising lead -- but confirm which
`sequencing_batch` (if any single one) corresponds to which pairwise experiment's source
plates before joining, since well-coordinate collisions across unrelated plates/batches are
already a known, confirmed problem in this project (see "corroborated_db naming collision"
note below) -- the exact same failure mode could easily repeat here if the join is done
naively.

## Non-obvious gotchas worth carrying over (also in Claude's memory for this project, but
restated here in case a fresh session doesn't load it)

1. **Strain names are not globally unique identifiers.** They're short plate-well-coordinate
   strings reused across unrelated plates/experiments/years. A name match is not an identity
   match anywhere in this project -- confirmed painfully with `corroborated_db.fasta` (an
   external 16S reference), where only 21/86 names for 20260721 turned out to be the same
   organism the experiment meant, vs. 19/19 for 20260630. Expect the same risk joining
   against the genomic tables' well-coordinate column.
2. **`high_uncertainty_pair`** wells/pairs (near-identical 16S references, forced ~50/50) are
   a measurement ceiling, not a biological coexistence signal -- already flagged in the data,
   but easy to silently poison an ML training set if not filtered.
3. Not every strain has full data everywhere: some lack a 16S consensus reference (excluded
   from relative-abundance labels), and presumably some lack genomic annotation coverage too
   -- check the intersection before building the modeling dataset, don't assume it's complete.
4. The interaction labels have real, quantified **measurement uncertainty**
   (`mean_uncertainty_score`, bootstrap CIs in `hierarchy_significance` outputs) -- worth
   propagating into the ML evaluation (e.g. sample weights, or excluding the least reliable
   labels) rather than treating every label as equally trustworthy ground truth.
5. The interaction network is a **real but imperfect hierarchy** (Bradley-Terry pseudo-R²
   ~0.80-0.83, ~5-6x more residual variation than sampling noise alone) -- so a genomic model
   should expect to explain most, not all, of the variance in pairwise outcomes; residual
   "upsets" may be genuinely opponent-specific (context-dependent) rather than noise.

## Suggested first steps for the new session

1. Resolve the strain-ID join (`mapping_384_well_plate_collection.csv` -> confirm which
   `sequencing_batch` applies, verify a handful of joins against known strain identities
   before trusting it at scale).
2. Build the modeling table: for each labeled pair in `r03_pair_replicate_stats.csv` (probably
   filtering `high_uncertainty_pair == False`), attach both strains' KO feature vectors from
   `KEGG_ko_and_strains_table.csv`.
3. Decide feature engineering for a strain-pair from two per-strain vectors (concatenation,
   difference, presence/absence overlap, a metabolic-complementarity score, etc.) and a
   target (`mean_relative_abundance_a` vs. `mean_log2_ratio_a_over_b` -- the latter is
   symmetric under swapping a/b, which likely matters for model architecture choice).
4. Given the folder-per-analysis convention already established, a new
   `<experiment>/analysis/genomic_ml/` (or a combined cross-experiment one, since the ML
   dataset would benefit from pooling both experiments' labeled pairs) would fit the existing
   structure.
