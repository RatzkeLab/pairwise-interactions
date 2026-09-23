# genomic_ml — 20260630: genome → who wins

**Question:** from two genomes alone, can we predict which strain dominates the well?
Target: `mean_log2_ratio_a_over_b` (antisymmetric).

Run: `run_genomic_ml.py` (headless) or `genomic_ml_20260630.ipynb` (the report).
Logic in `shared_pipelines/genomic_ml.py`.

```bash
/home/rl/mambaforge/envs/karl_seq_analysis/bin/python3 run_genomic_ml.py [--tabpfn]
```

| output | contents |
|---|---|
| **`g00_strain_join_validation.csv`** | **the gate — asserts before any model is fit** |
| `g01_dataset_summary.csv`, `g01_modeling_pairs.csv` | filter cascade and the modelling table |
| `g02_label_noise_ceiling.csv` | how good any model could possibly be |
| `g03_cv_summary.csv`, `g03_cv_predictions.csv` | model × regime results |
| `g04_strength_from_genome*.csv` | genome → per-strain competitiveness, KO vs 16S |
| `g06_*` | replicate / read-depth / strain learning curves (`replicate_value_relabund.py`) |

**Headline:** `cv_strain` ρ = 0.576, R² = 0.32; 16S-only control ρ = 0.458; shuffled ≈ 0.

**Read `cv_strain`, not `cv_pair`** — TRAPS §2. And note `R²` here is referenced to *predicting
no winner*, not to the mean (TRAPS §8).
