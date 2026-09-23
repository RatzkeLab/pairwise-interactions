# Demultiplexing runs — experiment 20260721

Demultiplexing sits between raw data and analysis: run rarely, referenced often. This
file is the record of every run against this experiment's reads, and which one to use.

**Use `20260730_demux`** — but read the status column, and read
[known_issues.md](known_issues.md) before trusting anything keyed by well label.

Output lives under `data_links/interim/demultiplexing/<dir>/`; per-well FASTQs are in
that directory's `unflipped/`.

| run date | output dir | config | `-e` | assigned | status |
|---|---|---|---|---|---|
| 2026-07-31 | `20260730_demultiplex_prelim1` | — | 6 | 35,652 | superseded — partial flow-cell output |
| 2026-07-31 | `20260730_demultiplex_prelim2` | — | 6 | 86,574 | superseded — partial |
| 2026-08-01 | `20260730_demultiplex_prelim3` | — | 6 | 100,600 | superseded — partial |
| 2026-08-03 | `20260730_demux` | `02_demultiplexing/demultiplex_config.yaml` | 6 | 114,620 — **15.6%** of filtered | **CURRENT, but under-assigning ~4×** |

The three prelims are progressive snapshots of the same run as the flow cell wrote more
data, not different methods. 2,458 of 9,120 sheet wells cleared 10 reads.

## Known under-assignment — not yet re-run

Same malformed barcode file as 20260630 and 20260907: the index column holds the whole
59 bp adapter+barcode+primer construct. Measured on 20260907, correcting it took
assignment from 13.5% to 64.1% of filtered reads. Expect a similar ~4× gain here.

Nothing already assigned is *wrong* — there are simply ~4× more reads waiting in
`unk.fastq`. See [known_issues.md](known_issues.md#minibar-barcode-file).

**Re-running is worth more here than elsewhere**: the 16S-based recovery of this
experiment's scrambled labels was depth-limited, and that is the approach that worked.

## Reproducing (with the fix)

```
python analyses/20260921/02_demux_qc/fix_barcode_file.py \
    experiments/20260721/01_setup/minibar_primers_20260721.tsv \
    experiments/20260721/01_setup/minibar_primers_20260721_corrected.tsv
# then copy the corrected config from 20260907, repoint fastq_input / barcodes_tsv / outdir
conda activate snakemake_env
cd /home/rl/scripts/ont-demultiplex-pipeline
snakemake -c 24 --use-conda --configfile <your config>
```

Use `-e 4`, and mind the config key names — `barcode_edit_distance` /
`primer_edit_distance`, not the `_dist` spellings.
