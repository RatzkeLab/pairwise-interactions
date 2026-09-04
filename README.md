# pairwise_interaction_experiments

ONT full-length 16S of bacterial strains grown pairwise in 384-well plates, plus plate-reader
optics, plus genomic annotation — asking whether interaction outcomes are predictable from
genomes.

## Start here

| document | what it is |
|---|---|
| **[FINDINGS.md](FINDINGS.md)** | what is established, suggestive, and untested — with status labels |
| **[TRAPS.md](TRAPS.md)** | ways this project produces confident wrong answers. **Read before trusting any number.** |
| [qc_summary/](qc_summary/) | `QC_01` labels · `QC_02` what happened to 20260721 · `QC_03` genomic-feature validity |
| [PROJECT_HANDOFF.md](PROJECT_HANDOFF.md) | the original 2026-08-27 handoff. Largely superseded — kept for history |

## The one thing to know first

**20260630 is the usable experiment. 20260721's strain labels do not name the organisms in its
wells** (FINDINGS §1). Its source plate was the same physical plate (§2) and the experiment was
executed consistently, but the mapping from wells to strain identities is broken somewhere
between the preculture and the Echo source plate (§3). Do not merge the two experiments into one
interaction network.

## Layout

```
shared_pipelines/           all reusable logic; notebooks and runners import from here
  paths.py                  every path pointing outside the repo (import these, don't hardcode)
  experiment_config.py      ExperimentConfig — per-experiment paths and references
  io_utils.py               read/reference IO, edit distances
  mapping_validation.py     does a well contain what it should
  relative_abundance.py     per-well interaction scoring, replicate stability, Bradley-Terry
  hierarchy_significance.py 5-test significance suite on the hierarchy claim
  genomic_ml.py             genome -> WHO WINS (antisymmetric target)
  genomic_ml_yield.py       genome -> TOTAL YIELD from OD (symmetric target)
  genomic_ml_plate.py       genome + plate-reader optics -> who wins
  metabolic_features.py     constructed competition / complementarity / cross-feeding features

<experiment>/analysis/<question>/    one folder per question, each with outputs/ and figures/
qc_summary/                 cross-experiment QC write-ups
plate_reader_comparison/    were the two experiments' source plates the same? (yes)
technical_replicates/       reproducibility baseline for both assays
strain_identity_qc/         sequence-level identity QC (parallel session)
```

Each `analysis/<question>/` folder holds a runner script, numbered output tables and figures.
Output prefixes are stable: `r0*` relative abundance, `g0*` genomic ML, `y0*` yield,
`s0*` feature sweep, `c0*` contamination, `q0*`/`j0*`/`m0*`/`n0*` per-folder.

## Running things

Everything runs in the **`karl_seq_analysis`** conda environment:

```bash
/home/rl/mambaforge/envs/karl_seq_analysis/bin/python3 <script>
```

Runners are standalone and re-runnable; notebooks re-execute top to bottom and are the
human-readable reports. Two conventions worth knowing:

- **Never join on a strain name without an identity check.** `genomic_ml.validate_strain_join()`
  is the gate — it asserts before any model is fit (TRAPS §1).
- **Quote `cv_strain`, not `cv_pair`.** The latter is inflated by strain recognition and a
  provably meaningless mapping scores R²=0.44 on it (TRAPS §2).

## Conventions

- Numeric thresholds live in their module beside the reasoning that justifies them — deliberately
  not centralised (see the note in `paths.py`).
- `R²` means different things in `genomic_ml` and `genomic_ml_yield` (TRAPS §8).
- Figures use one palette, defined in `relative_abundance.py` and imported by the rest.
