# Demultiplexing runs — experiment 20260630

Demultiplexing sits between raw data and analysis: run rarely, referenced often. This
file is the record of every run against this experiment's reads, and which one to use.

**Use `20260710_demultiplex`** — but read the status column before trusting read counts.

Output lives under `data_links/interim/demultiplexing/<dir>/`; per-well FASTQs are in
that directory's `unflipped/`.

| run date | output dir | config | `-e` | assigned | status |
|---|---|---|---|---|---|
| 2026-07-11 | `20260710_demultiplex_prelim1` | — | 6 | 119,140 | superseded (partial run) |
| 2026-07-10 | `20260710_demultiplex` | `02_demultiplex_sequencing_data/demultiplex_config.yaml` | 6 | 119,140 — **14.9%** of filtered | **CURRENT, but under-assigning ~4×** |

2,779 of 9,120 sheet wells cleared 10 reads. Only a subset of the 9,120 barcode
combinations was actually used, so that denominator is not the well count.

## Known under-assignment — not yet re-run

This experiment carries the same malformed barcode file as 20260907 did: the index column
holds the whole 59 bp adapter+barcode+primer construct, so minibar matches a 59-mer at the
barcode edit distance. Measured on 20260907, correcting it took assignment from 13.5% to
64.1% of filtered reads. Expect a similar ~4× gain here.

Nothing already assigned is *wrong* — the reads that were assigned were assigned correctly
— there are simply ~4× more of them waiting in `unk.fastq`. See
[known_issues.md](known_issues.md#minibar-barcode-file) for the fix and the calibration.

## Reproducing (with the fix)

```
python analyses/20260921/02_demux_qc/fix_barcode_file.py \
    experiments/20260630/01_setup/minibar_primers_20260630.tsv \
    experiments/20260630/01_setup/minibar_primers_20260630_corrected.tsv
# then copy the corrected config from 20260907, repoint fastq_input / barcodes_tsv / outdir
conda activate snakemake_env
cd /home/rl/scripts/ont-demultiplex-pipeline
snakemake -c 24 --use-conda --configfile <your config>
```

Use `-e 4`, and mind the config key names — they are `barcode_edit_distance` /
`primer_edit_distance`, not the `_dist` spellings every config in this project uses.
