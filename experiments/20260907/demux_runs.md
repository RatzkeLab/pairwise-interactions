# Demultiplexing runs — experiment 20260907

Demultiplexing sits between raw data and analysis: run rarely, referenced often. This
file is the record of every run against this experiment's reads, and which one to use.

**Use `20260907_demux_corrected`.** Everything else here is kept for provenance only.

Output lives under `data_links/interim/demultiplexing/<dir>/`; per-well FASTQs are in
that directory's `unflipped/`.

| run date | output dir | config | `-e` | assigned | status |
|---|---|---|---|---|---|
| 2026-09-19 | `20260907_demux` | `01_setup/demultiplex_config.yaml` | 6 | 364,880 — **13.5%** of filtered | **SUPERSEDED** — malformed barcode file, see [known_issues.md](known_issues.md#minibar-barcode-file) |
| 2026-09-21 | `20260907_demux_corrected` | `01_setup/demultiplex_config_corrected.yaml` | 4 | 1,764,185 — **64.1%** of filtered | **CURRENT** |

## What changed between them

The barcode TSV handed minibar the whole 59 bp construct (15 bp adapter + 24 bp barcode
+ 20 bp primer) in the column it reads as the *index*, so it matched a 59-mer at the
barcode edit distance instead of the 24-mer barcode. Correcting the file and dropping
`-e` from 6 to 4 took assignment from 13.5% to 64.1% of filtered reads.

Compare the **rates**, not the absolute counts: the corrected run saw slightly more raw
input (4,237,099 reads against 4,149,781) because more FASTQ files had been written by
the sequencer in between. Median per-well depth went from ~160 to ~750; usable wells
(>10 reads) from 1,946 to 1,995 of the 2,464 sequenced.

## Reading the numbers

- **4,928 wells are barcoded but only 2,464 were amplified.** Plates 9–16 are the frozen
  extension plates. They are deliberately kept in the sample sheet: never having been
  PCR'd, any read landing there is a false assignment by construction, which makes them a
  live mis-assignment monitor. They measured 0.26% in the original run and 0.34% in the
  corrected one.
- **19 of 192 barcodes failed as oligos** and are not recoverable by re-demultiplexing.
  See [known_issues.md](known_issues.md#failed-primer-wells).

## Reproducing

```
conda activate snakemake_env
cd /home/rl/scripts/ont-demultiplex-pipeline
snakemake -c 24 --use-conda --configfile <path to the config above>
```

Roughly 1 h: the demultiplexing itself is ~15 min, the remaining ~9,900 snakemake jobs
are the per-well merge and unflip steps.
