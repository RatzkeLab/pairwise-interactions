# 20260907 — pairwise interaction experiment setup

`experiment_setup.ipynb` builds the whole design; `build_notebook.py` regenerates
that notebook (so a parameter change doesn't mean hand-editing JSON), and
`separability.py` / `design.py` hold the logic it imports. Run in
`karl_seq_analysis`; the whole thing executes in about 4 seconds.

Layout mechanics follow `../../20260721/01_setup/experiment_setup.ipynb`. What is
new is *which wells get what*, driven by the reference-database work.

## What the design delivers

**8 plates get PCR'd and sequenced; 8 more are plate-reader-only and frozen.**
Plates 1–8 are self-sufficient — sequencing only those still gives every
monoculture and the complete replicate set, so stopping at 8 costs coverage, not
completeness.

| | sequenced (plates 1–8) | all 16 plates |
|---|---|---|
| wells | 2464 | 4928 |
| monoculture wells | 760 | 760 |
| technical-replicate wells | 250 | 250 |
| ratio-titration wells | 250 | 250 |
| density-titration wells | 90 | 90 |
| distinct interaction pairs | 1114 | 3628 |
| **interactions per strain** | median **6** (377/384 have ≥5) | median **20** |

- **Every one of the 384 strains has a monoculture on the sequenced plates.** 376
  have two, and all 384 have them on *different plates*, so one plate failing
  cannot cost both. That is the contamination check on the fresh master-plate
  copy, and the basis for a properly corroborated reference sequence.
- **50 pairs are replicated 5×, each copy on a different plate** — which is the
  only arrangement that tests plate-to-plate reproducibility rather than
  well-to-well. Everything else is n=1.
- **Zero pair wells are spent on strains that could not be told apart.** Every
  pair is ≥15 bp separated under the conservative combination of all four
  references. Only 4 wells have unknown separability, all involving `N16` and
  `N2`, the two strains with no eligible partner at all.
- **The 8 strains never sequenced in six ONT attempts** (`A4 A8 L14 N16 N2 N5
  O16 P9`) get 1 monoculture and 2 pairs each, not a full allocation.
- **340 wells go to the inoculum arms** (below), which is what pulled the median
  interactions per strain from 8 to 6. 377 of 384 strains still have ≥5.

You asked for 3–5 interactions per strain; 8 plates comfortably supports a median
of 8, so the binding constraint was your PCR effort, not the design. Running the
extension plates later triples that to 21.

## How pairs were chosen

`separability.py` builds a 384×384 matrix from the four independent references
(`ont_v3`, `pool1`, `sanger`, `ngs`), all projected onto one coordinate system by
`external_reference_crosscheck`. Two decisions:

- **Conservative combination** — a pair's separability is the *minimum* over
  sources covering both wells. If any source says they are close, the pair is
  skipped. The errors are asymmetric: an unresolvable pair yields nothing, while
  skipping a resolvable one costs one of ~66,000 alternatives. 65,991 pairs pass,
  5,774 are too close, 1,771 are uncomparable.
- **Random matchings, not random pairs** — each round is a near-perfect matching
  over the strain set, covering every strain exactly once. *k* rounds give every
  strain exactly *k* partners with no balancing step and nobody left short.

Ambiguity codes are handled by bitmask intersection, so pool1's 6% IUPAC content
costs nothing.

## Files

| file | use |
|---|---|
| `strain_layout_20260907.csv` | the design: every well, its contents, and the 16S separation of the pair |
| `echo_strains/echo_strains_20260907_plate_NN.csv` | Echo strain transfers, one file per plate |
| `primer_shooting/echo_primers_20260907_plate_NN.csv` | Echo barcode transfers, one file per plate |
| `echo_strains_{A,B}_20260907.csv` | whole run split in half for a source refill break |
| `minibar_primers_20260907_primary_only.tsv` | demultiplexing for the first sequencing run |
| `minibar_primers_20260907.tsv` | all 4928 wells, for if the extension plates get sequenced later |
| `strain_classes_20260907.csv` | per strain: which references exist, identity confirmed, sequencing history |

Every well across all 16 plates has its own unique forward×reverse barcode
combination (4928 of the 9120 available), so the frozen plates can be sequenced
later without redesigning anything. There is headroom to 24 plates (7392 wells)
if you raise `N_EXTENSION_PLATES`.

## Caveats

- **Monocultures get double the inoculum.** A mono well is 100 nL of the same
  strain twice, matching 20260721 so growth stays comparable with existing data —
  but it means a monoculture starts at twice the density of either member of a
  pair. Worth remembering when comparing mono against pair yields.
- Separability rests on references whose well labels disagree between the two
  plate lineages (see `external_reference_crosscheck`). The conservative
  combination absorbs most of that, but a pair called separable on the strength
  of one lineage could still surprise you.
- `N16` and `N2` have no eligible partner under any reference. Their 4 wells are
  a gamble by construction.

## Contamination technique — what the data actually points at

From `merge_consensus_sequences/cross_contamination_check`, which measured this
rather than guessed:

**Overgrowth is the bigger problem, not splashing.** A few identities spread
across the plate — one is resident in 3 wells but invades 37 — and mixed wells
went 41% → 67% over eight months. The single highest-value change is to **stop
propagating working plates from working plates**. Every new working copy should
come from the freezer master, never from the previous working copy; that is what
converts a one-off contamination into a plate-wide takeover.

**Local transfer is real but modest** (immediate neighbours enriched 1.30×,
p<0.0001). Two specifics worth acting on:

- **Rows are enriched, columns are not** (1.13×, p=0.0002 vs p=0.23, n.s.). That
  asymmetry points at something spanning a row rather than a column — worth
  checking how your Viaflo head is oriented relative to the plate, since a
  12-channel head laid along a row touches 12 wells that then share a
  contamination path. Turning the plate 90° so the head spans a column, or using
  fresh tips per row, would test and likely fix it.
- **96-well-quadrant neighbours (±2 rows/cols) are the most enriched signal at
  1.44×** (p=0.0001). That is the 96→384 stamping geometry, so it implicates the
  stamping step specifically: fresh tips per quadrant, and don't reuse a head
  across quadrants.

The Echo is acoustic and contactless, so it is not your risk — the deepwell
culture, the Viaflo steps, and the stamping are. Other standard measures that fit
what was observed: centrifuge plates before peeling seals (condensate on the seal
underside is a direct well-to-well bridge), peel slowly rather than snapping, and
keep culture times as short as gives usable density, since overgrowth is
time-dependent.

The 760 monoculture wells in this design are also the monitor: after this run you
will have a per-strain contamination readout on the fresh plate, and can compare
it to the working-plate baseline to see whether the fresh copy actually helped.


## The inoculum arms — separating a head start from fitness

A previous analysis found monoculture yield predicts pairwise outcome and
preculture density predicts the winner, but "denser inoculum gives a head start"
and "a fitter strain grows dense everywhere *and* wins" both fit that. The
standard design cannot separate them, so two arms were added. Together they cost
340 wells (14% of the sequenced capacity) and **nothing extra in handling** — the
Echo already takes a per-well volume from the transfer CSV.

**The discriminating quantity is the slope of final log-ratio against initial
log-ratio.** Slope ≈ 0 means the outcome is set by fitness and the starting ratio
is irrelevant; slope ≈ 1 means initial conditions carry straight through and the
assay is largely reporting the inoculum; in between means both. That is a
directly measurable number, not an inference.

| arm | design | what it answers |
|---|---|---|
| **monoculture contrast** | each strain's two monos at 200 nL and 100 nL — **free**, 376 paired observations | is final yield inoculum-dependent *at all*, i.e. are cultures at carrying capacity by readout? |
| **ratio titration** | 50 pairs × `25:175, 50:150, 100:100, 150:50, 175:25` — total held at 200 nL so only the ratio moves, ~49× span | the slope above |
| **density titration** | 30 pairs × `50:50, 100:100, 200:200` — ratio held at 1:1 | does absolute inoculum change the outcome, independent of ratio? |

All levels of one series sit on the **same plate**, so a dose-response is not
confounded by plate effects; different series are spread across all 8 plates so
plate effects average out between series rather than aliasing onto the ratio
axis. Series pairs come from 20260630 (whose well labels join correctly to the
collection, unlike 20260721) and are **stratified half strong-winner
(|log2 ratio| > 3), half near-neutral (< 1)** — a head start is most likely to
flip a strong-winner pair, while a ratio effect is easiest to detect in a neutral
one, and a null in only one stratum would be uninformative.

The monoculture contrast is the one to look at first, and it was free: it also
fixes a pre-existing flaw. A mono has always been 2×100 nL while each member of a
pair gets 100 nL, so the mono was never the matched control for a pair member.
Now one of the two is.

### What these arms do *not* settle

- **They do not normalise the baseline.** Varying the Echo volume varies cell
  number only if source wells are equally dense, and they are not — strains
  differ several-fold in preculture density. The titration is a controlled
  perturbation on top of an uncontrolled between-strain difference. It measures
  *sensitivity to ratio*, which is the causal question, but it does not remove
  the underlying confound.
- **The definitive fix would be OD-normalised inoculum**: read the source plate
  before shooting and compute per-well volumes that deliver equal cell numbers.
  The Echo takes those straight from the CSV, so it is only a scheduling change —
  the OD read has to happen before the transfer. Worth a subset arm next time
  rather than the whole plate, since it breaks comparability with 20260630/20260721.
- **Separating "head start" from "different carrying capacity" needs time-resolved
  growth.** If the plate-reader reads on these plates are kinetic rather than
  endpoint, that comes for free; if they are endpoint only, that specific question
  needs a dedicated run.

### Why not two source plates at different dilutions

It would be strictly worse. The Echo varies volume per destination well already,
so a second plate buys nothing except a manual dilution step (extra variance) and
a second culture in a possibly different growth state. The *only* case where a
diluted source plate helps is pushing below the Echo's 25 nL minimum — i.e.
ratios beyond about 7:1 each way. A 49× span end to end should be plenty to
estimate the slope; if it turns out not to be, that is when a diluted plate earns
its place.
