# Demultiplexing runs — experiment 20260630

Demultiplexing sits between raw data and analysis: run rarely, referenced often. This
file is the record of every run against this experiment's reads, and which one to use.

**Use `20260710_demultiplex_increased_tolerance`** (2026-09-24). The original run is kept for
comparison; both are analysed side by side in `analyses/20260710_increased_tolerance`.

Output lives under `data_links/interim/demultiplexing/<dir>/`; per-well FASTQs are in
that directory's `unflipped/`.

| run date | output dir | config | `-e` | assigned | status |
|---|---|---|---|---|---|
| 2026-07-11 | `20260710_demultiplex_prelim1` | — | 6 | 119,140 | superseded (partial run) |
| 2026-07-10 | `20260710_demultiplex` | `02_demultiplex_sequencing_data/demultiplex_config_old.yaml` | 6 (59 bp index) | 119,140 — 14.9% of filtered | superseded; under-assigns ~4× |
| 2026-09-24 | `20260710_demultiplex_increased_tolerance` | `02_demultiplex_sequencing_data/demultiplex_config.yaml` | 4 (24 bp index) | 498,260 — **62.2%** of filtered | **CURRENT** |
| 2026-09-24 | `20260710_demultiplex_increased_tolerance_pass_b` | `02_demultiplex_sequencing_data/demultiplex_config_pass_b.yaml` | 4 (24 bp index) | 917,581 — 62.3% of filtered | run b (post top-up) ALONE, for QC; not yet merged with pass_a |

Sheet wells clearing 10 reads: 2,786 (original) → 2,852 (corrected). The gain is depth per
well, not the number of usable wells — most sequenced wells already cleared 10 reads.
Only 3,080 of the 9,120 sheet wells (odd plates 1–19) were sequenced.

## The 2026-09-24 re-run

Same input (`fastq_pass_a` only, 800,746 filtered reads) and same NanoFilt settings; the only
changes are the barcode file (`01_setup/minibar_primers_20260630_corrected.tsv`, 24 bp index
instead of the 59 bp adapter+barcode+primer construct) and `-e 4`, set under the key name the
pipeline actually reads (`barcode_edit_distance`).

**Negative control.** Only plates 1, 3, …, 19 were sequenced, but the sheet covers all 30
(the unsequenced 20 include odd plates 21–29), and every plate draws from the same 96×96
barcode pool. Reads landing on an unsequenced plate are false by construction; scaling the
per-well rate there to the 3,080 sequenced wells estimates the false-assignment rate on real
wells:

| | original | corrected |
|---|---|---|
| assigned | 119,140 (14.9%) | 498,260 (62.2%) — 4.18× |
| on sequenced plates | 117,994 | 492,401 |
| on unsequenced plates | 1,146 | 5,859 |
| est. false-assignment rate | 0.50% | 0.61% |
| median reads / sequenced well | 39 | 166 |

Every sequenced well gained reads (none lost any; median gain 4.18×, IQR 3.8–4.7×). Most of
the false rate was already there at the strict setting (chimeras / index hopping, which no
edit threshold fixes); the corrected matching adds ~0.1 percentage points.

## fastq_pass_b: QC'd, same library, safe to merge (2026-09-24)

`fastq_pass_b` is a second MinKNOW run on the same flow cell after a library top-up (28.4 h vs
12.8 h). It was demultiplexed alone with identical settings and compared to pass_a well by well
(`analyses/20260710_increased_tolerance/02_pass_ab_qc/`):

| | pass_a | pass_b |
|---|---|---|
| assigned | 498,260 (62.2%) | 917,581 (62.3%) |
| est. false-assignment rate | 0.61% | 0.61% |
| minibar `HH` (barcode + primer, both ends) | 99.7% | 99.6% |
| median reads / sequenced well | 166 | 307 |
| full-length share of raw pass reads | 31.8% | 27.6% (top-up material more fragmented) |

- Well shares agree: every plate's share is within 1% between runs; dispersion 1.95× Poisson;
  1 of 3,080 wells changed share by 2×.
- Strain ratios agree to counting noise: 1,856 pair wells, r = 0.993, z-score SD 1.05
  (1 = binomial noise only), 0.29% of wells with |z| > 3 (0.27% expected).
- Mapping QC agrees: same qc_status in 96.5% of wells, same off-target distribution.

So pass_b is the same pool, sequenced deeper. Merging gives ~2.8× pass_a's depth (median
~470 reads/well).

## Reproducing

```
conda activate snakemake_env
cd /home/rl/scripts/ont-demultiplex-pipeline
snakemake -c 48 --use-conda --configfile \
    /home/rl/scripts/karl/pairwise_interaction_experiments/experiments/20260630/02_demultiplex_sequencing_data/demultiplex_config.yaml
```

~1 h on 48 cores: ~12 min just building the DAG, then 9,120 per-sample merge jobs dominate. The corrected barcode file was
made with `analyses/20260921/02_demux_qc/fix_barcode_file.py`.
