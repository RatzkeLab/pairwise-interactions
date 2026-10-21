# 20261021 — community experiment, round 1: setup

Higher-order communities, mostly 3–10 members plus bigger ones up to 66, varying in
relatedness, in triplicate, same medium and culture conditions as the pairwise runs,
inoculated on 2026-10-21 from the **same source plate as 20260907**. Goal: data for
predicting community composition from the member genomes.

Built with the shared library `scripts/experiment_setup/` (the structure proposed in
`experiments/README.md`): every choice is in `setup_config.yaml`; code never lives with its
outputs; every output is validated.

| file | what |
|---|---|
| `setup_config.yaml` | every per-experiment choice |
| `run_setup.py` | builds `outputs/` and validates it; refuses to overwrite a dated experiment |
| `build_inspect_notebook.py` → `inspect_setup.ipynb` | plots of the design; writes nothing |
| `outputs/echo_strains/echo_strains_20261021_plate_NN.csv` | **strain picklists, plates 1–5** |
| `outputs/primer_shooting/echo_primers_20261021_plate_NN.csv` | primer picklists, plates 1–8 |
| `outputs/layout.csv` | every well: community, members, replicate, time point, barcodes |
| `outputs/communities.csv` | every community: members, module, relatedness, inoculum |
| `outputs/minibar_sheet_20261021.tsv`, `demultiplex_config_20261021.yaml` | demultiplexing |

## The design

| block | plates | wells | contents |
|---|---|---|---|
| main (24 h, sequenced) | 1–4 | 1224 | 408 communities × 3, each replicate on a different plate |
| passage (daily dilution) | 5 | 306 | 102 communities × 3, all on the plate, replicates never in one row or adjacent; sequenced at 24 h / 48 h / 5 d / 10 d as barcode plates 5 / 6 / 7 / 8 |

2 blank-medium wells per plate; edge wells unused; no monocultures. **One flowcell =
plates 1–8 (2,464 wells).**

**Main plates** (408 communities; share of the 1,224 wells):

| module | communities | wells | |
|---|---|---|---|
| small, 3–10 members | 340 | 83% | ~42 per size, each size spread over three relatedness terciles |
| 12–20 members | 16 | 3.9% | |
| large tests | 9 | 2.2% | 15, 20, 25, 30, 40 members; 51 / 55 = separable ceiling at 30 bp (core / +medium); 61 / 66 = at 15 bp (the +medium ceiling of 77 trimmed to 66 by the 5 µL limit) |
| trios + their 3 pairs each | 5 + 15 | 4.9% | does pairwise predict the trio? 5 of the 15 pairs were also measured in 20260907 |
| pairs re-measured from 20260907 | 8 | 2.0% | identity check on the source plate + cross-experiment replicate |
| leave-one-out series | 9 | 2.2% | an 8-member base and each member removed |
| resolution tests | 6 | 1.5% | one deliberately close pair (10–14 or 16–26 bp), rest ≥30 bp |

**Passage plate** (102 communities, core strains only): **all 7 core large tests (15, 20,
25, 30, 40, 51, 61 members)** — the same communities as on the main plates, so each has a
24 h triplicate there too — plus 95 communities of 3–10 members, ~¼ of them also on the
main plates.

**Who is in the pool.** One well per distinct genome whose genome was verified against the
well (`analyses/20261021/03_genome_concordance`): 75 high-confidence ("core") + 58
medium-confidence ("extended"), 133 in all, including the 7 whose preculture barely grew
in September (B4 D18 F14 I8 J12 K11 M18). Medium-confidence strains share their 16S with a
genomically different relative, so a swap would be invisible; they only appear in the ~⅓ of
communities marked `extended` in `communities.csv`, so every analysis can also be run on core
communities alone.

**Separation.** All members of a community ≥30 bp apart (haplotype-aware 16S), except the
two 15 bp ceiling tests and the 6 resolution tests. Read cross-talk is already at its ~0.7%
floor from ~7 bp, so 30 bp is a wide margin.

**Relatedness.** Mean pairwise 16S distance between members, filled in three terciles per
size — about 40 bp (one clade) to 330 bp (across phyla). It tracks genome content closely
(r = 0.86 with mean pairwise KO-profile Jaccard distance, also in `communities.csv`).

## Inoculation

**Equal volumes, no OD correction this round**: 100 nL per member; 75 nL per member in
communities above 50 members, so no well receives more than 5 µL (10% of 50 µL). Largest
well: 4.95 µL (66 × 75 nL); median 0.7 µL.

**No source refill**: the most-used strain draws 13.3 µL over the run (median 11.2 µL),
under the 25 µL a fill allows. 11,484 transfers / 44,901 droplets (20260907: 9,096 /
38,040). Shoot `outputs/echo_strains/` plate by plate; `echo_strains_20261021_all_plates.csv`
is the same transfers in one file.

Starting cell numbers are therefore NOT equal: source ODs in September spanned ~25-fold
(0.14–3.5), and composition partly remembers the starting ratio (20260907: slope 0.31)
though not the absolute density (0.04 ± 0.07 per doubling, p = 0.6). Record today's
source-plate OD — it is the covariate the analysis needs, and the basis for the OD-corrected
rerun planned for later (`scripts/experiment_setup/echo.normalised_volumes`).

## Barcodes and demultiplexing

- Every well on barcode plates 1–8 has its own forward × reverse combination, plus 308
  decoy combinations (`Plate99_*`) on the sheet that are in no well: reads assigned to them
  measure mis-assignment (the role the frozen plates played in 20260907).
- **The 19 primer oligos that failed in 20260907 are excluded** (fwd 29 30 44 45 46 64 65
  74; rev 7 28 33 38 41 43 44 49 61 64 77). If the primer plate is remade, empty those
  inputs and rebuild.
- The minibar index columns are the **24 bp barcode** — checked by
  `validate.check_minibar`, and the barcode code reproduces the corrected 20260907 sheet for
  4,928/4,928 wells. Barcode set minimum edit distance 12; `-e 4`.
- Fill in `fastq_input` in the demux config after the run.

## Deferred

OD-corrected inoculation (and re-running some of these communities with it); the
leave-two-out series; OD-only frozen extension plates; the large-scale large-community
experiment on its own high-depth flowcell. The passage dilution factor/method is not part of
this setup; the layout is identical at every time point.
