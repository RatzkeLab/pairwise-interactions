# 20260907 — analysis

Run `python run_all.py` (env `karl_seq_analysis`, ~5 min) to reproduce everything below
from the corrected demultiplexing. `--original` reruns against the first demux for the
before/after comparison.

## Headline results

| | |
|---|---|
| reads assigned | 364,880 → **1,764,185** (4.83×, after fixing the barcode file) |
| median well depth | ~160 → **~750** |
| usable wells (>10 reads) | 1995 of 2464 sequenced (81.0%) |
| strains with a corroborated consensus | **345 / 384 (89.8%)**, 224 backed by two monocultures |
| two-strain wells with no second sequence | **1.3%** |
| mono wells contaminated ≥10% | **24.3%** |
| label support vs collection (20260630 = 58.7%, 20260721 = 0%) | **51.9%** |
| ratio-titration slope | **+0.31** |

## The four questions

**1. Is the dropout a primer problem or a strain problem? → primer.**
Barcodes are randomised against plate position here, so the causes separate. Dispersion
against a binomial null: reverse barcode **9.86**, forward barcode **6.81**, plate 6.43,
strain 1.41, row/column ~1.0. **19 barcodes of 192 carry 76% of the dropout.** Confirmed
as a failed oligo three independent ways, the cleanest being OD: those wells grew at 95.6%
against 96.7% elsewhere (p=0.27) — the cells were there, the amplicon was not. Details and
the source-plate map in `02_demux_qc/README.md`.

**2. One sequence per mono well, two per pair well? → yes, with a caveat about thresholds.**
Of 1411 two-strain wells, only **1.3% hold no second sequence**; both are present at ≥10%
in 48.5% and at ≥1% in 79.6%. A low minority fraction is competitive exclusion, not a
detection failure. Mono wells: 75.9% single-sequence, 24.3% contaminated at ≥10% — the
per-well contamination readout the 760 monocultures were designed to give.

**3. Reference sequences → 345 strains, anchored on monocultures.**
`shared_scripts/strain_consensus.py` treats the two monocultures as the anchor and the
pair wells as corroboration, rather than taking a cross-well majority as 20260721 had to.
That is what makes "both monos agree" a meaningful status: only **one** strain (I16) has
monocultures that disagree. The set reaches **99% of the resolution ceiling** 16S permits
for this collection.

**4. Is the 20260721 label problem present here? → no.**

| run | full_collection *(primary)* | genome 16S | corroborated *(secondary)* |
|---|---|---|---|
| 20260630 *(positive control)* | 58.7% | 34.1% | 100% |
| 20260721 *(negative control)* | 0.0% | 0.0% | 0.0% |
| **20260907** | **51.9%** | **44.7%** | **89.7%** |

20260907 tracks the positive control on every reference and beats it on the genome link;
the negative control fails at 0% on all three. The numbers barely move between the shallow
and corrected demux, so this is not depth-dependent. Disagreements have no geometry — no
row shift, column shift or rotation explains more than 3.2% — so there is no transform to
recover. The genomic-table link is weak (44.7%), **but weak for everyone**: 20260630, the
run the ML was built on, scores 34.1%. That is a property of the genome-derived reference,
which is assembly-derived and half fragments, not of this run.

## The inoculum arms

Yield and composition gave different answers, and both matter.

| arm | measures | result |
|---|---|---|
| monoculture contrast (200 vs 100 nL, paired within strain, n=346) | is yield inoculum-dependent? | **1.8%** of the doubling survives (p=0.002) |
| density titration (4× span, ratio fixed 1:1, n=30) | does absolute inoculum matter? | **−1.1%** of the span (p=0.54) |
| ratio titration (≈49× span, n=34) | does *composition* remember the inoculum? | **slope +0.31** (p vs 0 = 6e-09; p vs 1 = 5e-09) |

Total yield is at carrying capacity by readout. Composition is not: **about 31% of a
change in starting ratio survives to the endpoint.** Median r² is 0.725, and among the 21
pairs fitting at r²>0.5 the slope is +0.53. The high-ratio levels compress (level medians
−0.93, −0.40, 0.07, 0.34, 0.38), which is saturation and flattens the fit, so +0.31 is
likely an underestimate. Strong-winner pairs are *steeper* (+0.36) than near-neutral ones
(+0.14).

**So a single well's interaction score is partly reporting starting density.** Either
correct with the slope or carry it as a sensitivity analysis — it compounds the earlier
finding that preculture density predicts the winner. OD-normalised inoculum is the real
fix, and only needs the source plate read before the Echo shoot.

## Before modelling

- **Exclude the 30 strains that never grew** (`02_demux_qc/outputs/d13_strains_no_growth.csv`).
  Both monocultures failed, on different plates — 9× the independent expectation — and
  their pair wells grow like monocultures (median OD 0.259 vs 0.349, p=1.4e-13).
- **210 of 1050 tested pairs lack a reference for one strain**, so those wells cannot be scored.
- **Read assignment is clean**: 1,007,148 reads, strain1 469,058 / strain2 480,412 (balanced,
  so assignment is unbiased), 4.5% off-target, median per-well uncertainty 0.000.

## Layout

| directory | contents |
|---|---|
| `01_setup/` | the design (see its own README) |
| `02_demux_qc/` | why wells are empty; the barcode-file fix and its calibration |
| `03_recreate_reference_db/` | per-strain consensus |
| `04_qc/` | sequences per well |
| `05_cross_reference/` | label support against every reference and both prior runs |
| `06_inoculum_arms/` | the three inoculum arms |
| `analysis/` | shared-pipeline outputs (mapping_validation, relative_abundance) |

Shared logic added to `../shared_scripts/`: `well_composition.py`, `strain_consensus.py`,
`strain_identity.py`, `plate_reader.py`.

## Two traps worth knowing

**Compare sequences over their overlap, not end to end.** The corrected demux trims ~20 bp
less at each end (minibar removes whatever it is handed as the "index"), so old and new
reads differ by ~44 bp for reasons unrelated to identity. Half the genome-derived rDNA
entries are assembly fragments as short as 396 bp, with the same effect. Using overlap-aware
distance took the genome-16S scores from 11.5%→34.1% (20260630) and 21.4%→44.0% (20260907).
`strain_identity.flexible_distance` does it.

**Single-linkage clustering chains at depth.** The method inherited from
`20260721/04_qc/check_reads.ipynb` merges the two strains of a pair well once there are
enough error-rich reads to bridge the distance modes — and it gets *worse* with more data
(27.0% → 10.0% of wells resolving both strains, going from 120 to 500 reads, while greedy
seed clustering held at 47.0% → 45.0%). `well_composition.greedy_clusters` replaces it.
**This also means 20260721's own strain-count QC understated its two-strain wells.**
