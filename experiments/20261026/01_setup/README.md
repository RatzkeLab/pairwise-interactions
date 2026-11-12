# 20261026 — community experiment, round 2: setup

**20261021 again, well for well, plus a second passage plate.** Plates 1–5 hold exactly
the communities, replicates and wells of 20261021 (same source plate, same equal-volume
inoculation, byte-identical Echo strain picklists). Plate 6 is new: a second
self-contained passage plate built the way plate 5 was and diluted daily in parallel with
it. 31 of its 102 communities are also on plate 5, as technical replicates, to measure how
much two plates diverge under the daily-dilution scheme.

The design rationale (pool, separation, relatedness, inoculum) is in
`../../20261021/01_setup/README.md` and is unchanged; this file covers only what differs.

| file | what |
|---|---|
| `setup_config.yaml` | every per-experiment choice |
| `run_setup.py` | copies plates 1–5 from 20261021's outputs, builds plate 6, writes `outputs/`, validates; refuses to overwrite a dated experiment |
| `build_inspect_notebook.py` → `inspect_setup.ipynb` | plots of the design; writes nothing |
| `outputs/echo_strains/echo_strains_20261026_plate_NN.csv` | **strain picklists, plates 1–6** (1–5 identical in content to 20261021's) |
| `outputs/primer_shooting/echo_primers_20261026_plate_NN.csv` | primer picklists, barcode plates 1–12 (**new barcodes**, see below) |
| `outputs/layout.csv` | every well: community, members, replicate, `culture_plate`, time point, barcodes |
| `outputs/communities.csv` | every community; `passage_plates` = `5 6` marks the between-plate replicates, `plate6_origin` says why a community is on plate 6 |
| `outputs/minibar_sheet_20261026.tsv`, `demultiplex_config_20261026.yaml` | demultiplexing |

## The design

| block | plates | wells | contents |
|---|---|---|---|
| main (24 h, sequenced) | 1–4 | 1224 | 20261021's 408 communities × 3, unchanged |
| passage plate 5 (daily dilution) | 5 | 306 | 20261021's 102 communities × 3, unchanged |
| passage plate 6 (daily dilution) | 6 | 306 | 102 communities × 3, built like plate 5 (below) |

Each passage plate is sequenced at 24 h / 48 h / 5 d / 10 d. **One flowcell = barcode
plates 1–12 (3,696 wells)**, plus 308 decoy combinations (`Plate99_*`) on the sheet.

**Plate 6**, core strains only, same size mix as plate 5 (12 communities of each size 3–9,
11 of size 10, plus the 7 large tests):

| origin | communities | |
|---|---|---|
| core large tests (15–61 members) | 7 | the same 7 as on plate 5, as plate 5 has them |
| shared with plate 5 | 24 | 3 per size 3–10, one from each relatedness tercile of plate 5's communities of that size |
| also on the main plates | 18 | ¼ of the remaining small slots, as on plate 5 |
| new | 53 | built like plate 5's (size × relatedness terciles), never a duplicate of any 20261021 community; IDs K0480 onward |

So **31 communities (93 wells, 30% of each plate) are on both passage plates**, each in
triplicate on each. Their positions on plate 6 are re-randomised under the same spacing
rules (replicates never share a row, never neighbours), not mirrored from plate 5: a
between-plate difference and a within-plate replicate difference then both include well
position, so the excess of between-plate over within-plate variance is the plate effect.
32 of plate 6's communities are also on the main plates (plate 5: 31).

## Barcode plates

Barcode plate N is passage plate N at 24 h; later time points follow:

| barcode plate | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|
| passage plate | 5 | 6 | 5 | 6 | 5 | 6 | 5 | 6 |
| time point | 24 h | 24 h | 48 h | 48 h | 5 d | 5 d | 10 d | 10 d |

**This differs from 20261021, where barcode plates 6–8 were plate 5 at 48 h / 5 d / 10 d.**
`layout.csv` has `culture_plate` and `timepoint` for every barcode plate.

## Barcodes

- **No forward × reverse combination of 20261021 (its decoys included) is reused.** The two
  experiments have identical layouts, so a sample swapped between them would otherwise look
  perfectly plausible; with disjoint combinations every read traces to its experiment, and
  the two could even share a flowcell (sample IDs would need an experiment prefix then).
  Hence the primer picklists are new for every plate, including 1–5.
- Same 19 failed 20260907 oligos excluded; same 24 bp minibar index check; barcode set
  minimum edit distance 12; `-e 4`.
- Primer drawn per oligo source well: max 14.0 µL here; 19.8 µL if 20261021's plates are
  shot from the same primer plate.
- Fill in `fastq_input` in the demux config after the run.

## Inoculation

Unchanged: 100 nL per member (75 nL above 50 members), no OD correction. Six plates draw
at most 17.4 µL from a source well (median 13.8 µL), under the 25 µL a fill allows: **no
refill**. 14,052 transfers / 54,837 droplets. Shoot `outputs/echo_strains/` plate by
plate. As for 20261021, record the day's source-plate OD: it is the covariate the
analysis needs.

## Checks (`outputs/validation.txt`)

Besides 20261021's checks (communities, replicate placement, barcodes, Echo volumes,
minibar sheet): the inputs and 20261021's outputs match 20261021's manifest; plates 1–5
match 20261021 in every well and every strain transfer; every shared community is on both
passage plates in full triplicate with identical per-member volumes; plate 6 is core-only;
no barcode combination overlaps 20261021's.
