# Analyses using experiment 20260907

One row per analysis run that touched this experiment. Read this instead of opening the
analysis folders; open them only for method detail or for code worth reusing.

| analysis | date | what it asked | what it found |
|---|---|---|---|
| [`analyses/20260921`](../../analyses/20260921) | 2026-09-21 → 22 | Full QC and first-pass genomic ML on this run | Demux was under-assigning ~5×, fixed. Labels are **not** scrambled (unlike 20260721). 370/384 strains got a corroborated consensus. Genomics predicts held-out strains at R² 0.37. |

## Headline results (from `analyses/20260921`)

| | |
|---|---|
| reads assigned | 364,880 → **1,764,185** after the barcode fix |
| usable wells (>10 reads) | 1,995 of 2,464 sequenced (81.0%) |
| strains with a corroborated consensus | **370 / 384 (96.4%)** — 345 anchored on a monoculture, 25 on pair wells |
| two-strain wells with no second sequence | **1.3%** |
| mono wells contaminated ≥10% | **24.3%** |
| label support vs collection reference | **51.9%** (20260630 positive control 58.7%; 20260721 negative control 0%) |
| well → genome join gate | **PASS**, ρ=+0.387, z=+8.2 (20260630 passed at +0.363/+6.0) |
| genomic ML, held-out strains (`cv_strain`) | R² **0.372**, Spearman 0.621, sign accuracy 83.4% |
| ratio-titration slope | **+0.31** — composition partly remembers the inoculum |

## Where this experiment is used as a reference

`analyses/20260921` also scores **20260630** and **20260721** against the same
references, as positive and negative controls for the label-identity check.
