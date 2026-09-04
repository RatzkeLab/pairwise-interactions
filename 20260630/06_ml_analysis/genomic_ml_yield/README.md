# genomic_ml_yield — 20260630: genome → total yield

**Question:** from two genomes, how much biomass does the pair make together? Target: plate
z-scored OD600 (**symmetric** — a different question from who-wins, not a substitute).

Run: `run_genomic_ml_yield.py`. Logic in `shared_pipelines/genomic_ml_yield.py`.

Why it exists: OD is independent of 16S resolvability, so nothing is discarded as
`high_uncertainty`. **3828 pairs / 88 strains** versus 1465 / 74 for relative abundance.

| output | contents |
|---|---|
| `y01_dataset_summary.csv`, `y01_pairs.csv` | filter cascade, modelling table |
| `y02_replicate_reliability.csv` | split-half ceiling (ρ 0.869) — drawn on the model figure |
| `y03_cv_summary.csv` | model × regime results |
| `y04_yield_from_genome_metrics.csv` | genome → per-strain yield contribution, KO vs 16S |
| `y05_*` (`replicate_value.py`) | replicate and strain learning curves |

**Headline:** `cv_strain` ρ = 0.568, R² = 0.23. Unlike relative abundance, **16S alone nearly
matches KO here** (R² 0.411 vs 0.406) — yield is more phylogenetically determined.

**The strain learning curve is still climbing at 88 genomes** (`y05_strain_learning_curve.csv`):
more genomes is the high-value experiment for this target. More replicates is not (+0.02).
