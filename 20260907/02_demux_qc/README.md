# 20260907 demultiplexing QC

Why so many wells look empty, and what to do about it. Three separate things were
going on; only one of them is a real loss.

## Summary

| | |
|---|---|
| wells in the demux output | 4928 |
| wells actually PCR'd and sequenced | **2464** (plates 1–8) |
| the other 2464 | frozen extension plates, never amplified — empty by design |
| dropout among sequenced wells (≤10 reads) | 519 (21.1%), of which 411 at exactly zero |
| dropout attributable to 19 bad barcodes | 394 (76%) |
| dropout among wells with two clean barcodes | 125 (6.3%) |
| dropout wells that nevertheless had cells (OD) | 480 of 519 (92.5%) |
| reads assigned by the original demux | 364,880 = 13.5% of filtered |
| reads assigned after the barcode-file fix | ~1.7M = ~63% of filtered |

## 1. Half the "empty" wells were never meant to have reads

The demux ran against the full 4928-well sample sheet, but only plates 1–8 were
amplified. Keeping the frozen half in the sheet is worth doing — those 2464 wells
carry real barcodes and no DNA, so any read landing there is a false assignment.
They measured 0.26% mis-assignment in the original run, and they are what
calibrated the edit distance for the re-run (`calibrate_minibar.py`).

## 2. A real primer failure: 19 barcodes out of 192

Barcodes are assigned to wells at random in this design, independently of plate
position, which is what makes the causes separable. Dispersion against a binomial
null (1.0 = no structure):

| factor | dispersion |
|---|---|
| reverse barcode index | **9.6** |
| forward barcode index | **6.7** |
| destination plate | 5.9 |
| strain | 1.5 |
| destination row / column | 0.9 / 1.1 |

7 forward + 12 reverse barcodes cover 19% of wells and carry 76% of the dropout.
Wells with two clean barcodes drop out at 6.3%; wells touching a flagged barcode,
82%.

**It is the oligo, not the sequence and not the strain.** Three independent checks:

- those barcodes were fine in 20260630 and 20260721 (20–26% dropout then, 81–85%
  now), so it is not a property of the sequence
- their amplicons are absent from the discarded-read pile too — 1.6 recoverable
  reads per well against 15.2 for clean barcodes, a 10× depletion. The amplicon
  was never made
- **OD says the cells were there**: flagged-barcode wells grew at 95.6% against
  96.7% for the rest (Fisher p = 0.27), median OD 0.317 against 0.325. That
  closes it without using sequencing at all

Not depletion, either — source wells that have been drawn from *more* fail *less*
(rho = −0.22, p = 0.002). The plate was the same one, freshly thawed, so
evaporation or freeze-thaw damage in specific wells fits. Row H of the primer
source plate holds 5 of the 12 bad reverse barcodes; G17/G19/G21, K9/K11 and
H12/H16/H18 are adjacent runs. `outputs/d03_*_barcode_forensics.csv` lists them
with their source wells — worth checking those wells at the bench.

## 3. The big one: the barcode file was malformed, costing 4.6× the data

minibar reads column 2 of the barcode TSV as the forward index. The generated
file puts the whole 59 bp construct there — 15 bp ligation adapter + 24 bp
barcode + 20 bp 16S primer — so minibar matched a 59-mer within the barcode edit
distance instead of the 24-mer barcode. Same tolerance, 2.5× the length.

84% of the discarded reads are clean full-length amplicons whose barcodes match
the sheet at **median edit distance 0**. They were never unreadable; they were
never matched.

Calibration on a representative chunk, false rate measured against the frozen
plates:

| barcode file | `-e` | assigned | est. false rate |
|---|---|---|---|
| as generated | 6 | 13.7% | 0.15% |
| **corrected** | **4** | **63.6%** | **1.23%** |
| corrected | 5 | 67.6% | 1.48% |
| corrected | 6 | 72.7% | 5.15% |

The barcode set's minimum pairwise edit distance is 12, so `-e 6` sits exactly on
the decision boundary and `-e 7`+ assigns noise. `-e 4` is the operating point.

A second, latent bug found in the same place: `demultiplex.smk` reads the config
keys `barcode_edit_distance` / `primer_edit_distance`, while every config in this
project writes `barcode_edit_dist` / `primer_edit_dist`. Those keys have never been
read — the rule has always used its own defaults. They happened to match, so
nothing broke, but a changed setting would have silently done nothing.

**20260630 and 20260721 have the identical barcode file layout** and are
under-demultiplexed the same way. See `TODO_demux_underassignment.txt` in each.

## 4. Two incidental findings

**The flow cell ran in two phases.** Chunks 0–5 pass only 13–17% of reads through
the length filter; chunks 6–22 are empty; chunks 23–39 pass 67% and hold 98.6% of
the data. Looks like a restart or reload partway in.

**30 strains did not grow.** Every strain has two monocultures on *different*
plates, so independent well-level failures would almost never hit both. 30 strains
failed both — a 9× enrichment over the independent expectation of 3.2. Their pair
wells confirm it (median OD 0.259 against 0.349, p = 1.4e-13): a pair carrying a
dead strain grows like a monoculture. Listed in
`outputs/d13_strains_no_growth.csv`. These contribute no interaction data and
should be excluded before modelling.

Conversely, **7 of the 8 strains never recovered in six previous ONT attempts
worked this time** — A4, A8, N16, N2, N5, O16, P9, all with reads and no dropout.

## Files

| script | what it does |
|---|---|
| `dropout_analysis.py` | barcode vs strain vs position, dispersion against a binomial null |
| `barcode_forensics.py` | which barcodes fail, their primer-source wells, attribution |
| `barcode_mechanism.py` | the same barcodes across all three runs — sequence or oligo? |
| `primer_plate_geometry.py` | depletion vs local damage on the primer source plate |
| `collision_test.py` | barcode-ambiguity hypothesis (ruled out: minimum separation is 12) |
| `rescue_probe.py` | are the failing barcodes' amplicons in the discard pile? (no) |
| `fix_barcode_file.py` | rewrites a barcode TSV so the index column holds the index |
| `calibrate_minibar.py` | picks `-e`/`-l` using the frozen plates as a negative control |
| `od_vs_dropout.py` | did the dead wells have cells? (yes) and which strains never grew |
| `demultiplex_config_corrected.yaml` | the config for the re-run |
