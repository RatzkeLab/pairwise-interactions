# QC part 2 — what actually happened to 20260721, and everything tried to reconstruct it

Companion to [QC_01_NAIVE_MAPPING.md](QC_01_NAIVE_MAPPING.md), which reports the same QC taking
the labels at face value. This document drops that assumption.

**Bottom line up front.** 20260721's wells hold real, clean, single collection organisms that are
not the ones its labels claim. The frozen stock and the source plate were correct; the pick lists
and barcode files are internally consistent; no plate-handling geometry explains the scramble.
About half the plate can be re-identified from 16S well enough to be predictive, but not to
strain-level certainty. The error sits between the preculture plate and the Echo source plate —
the one step nothing in the repository records.

---

## 1. What was ruled out, and how

| hypothesis | test | verdict |
|---|---|---|
| Bad sequencing / low depth | depth 42.8 vs 41.9 reads/well | **No** — indistinguishable from the good run |
| Two-strain deconvolution error | mono wells scored separately | **No** — mono wells 0% supported, same as pair |
| DNA contamination pre-PCR | excess 16S groups per well vs control | **No** — 18.4% vs 17.8% (pair), and mono wells *cleaner* than control (23% vs 50%) |
| External/reagent contaminant | pairwise identity among low-identity reads | **No** — diffuse (median 0.81, 0.5% of pairs >0.95), not one organism |
| Different frozen stock | plate-reader OD, source plates | **No** — ρ=0.698 at zero offset, ~0 at any one-well shift: *same plate, same layout* |
| Pick list ≠ layout | `echo_strains` vs `strain_layout` | **No** — 9120/9120 wells agree, both experiments |
| Barcode/primer files differ | all primer files diffed | **No** — byte-identical between experiments |
| Plate-handling geometry | 335 transforms (16S) + OD transform search | **No** — see §3 |

## 2. The two decisive positive findings

**The source plates are the same.** Per-well OD, 298 occupied wells, row/column detrended:
20260630's dense NM source plate vs 20260721's dense preculture gives **ρ = 0.698**, against
same-plate positive controls of 0.732 / 0.777 and a within-experiment cross-media comparison of
0.440. Decisively, the correlation is **0.698 at zero offset and ≤0.011 at every one- or two-well
shift** — shared layout, not a shared evaporation gradient.

**The destination wells do not contain the named pairs.** 171 strain pairs appear in both
experiments' layouts. Their OD correlates at **ρ = 0.084** across experiments, against
within-experiment replicate ceilings of **0.816** and **0.622** on those same pairs. Since shared
pairs sit at different plates and wells in each experiment, nothing but strain content links the
two measurements. Attenuation alone would cap a perfect match at ~0.77; observed is 0.08.

Together: right plate, wrong contents downstream of it.

## 3. Why "some plate got rotated" is not the answer

Two independent searches, both null:

- **16S-based** (`strain_identity_qc/qc_layout.py`): 335 candidate transforms — all 24 quadrant
  permutations under both 96→384 recombination conventions, 180° rotations, row/column flips,
  pick-list off-by-N shifts, and compositions. Best for 20260721: `quad_block_2103+rot180_within_96`
  at z=3.01, family-wise **p = 0.234**. The positive control (20260630) recovers `identity` at
  z=8.18. That negative was explicitly underpowered — only 5 usable single-strain calls survive
  16S degeneracy.
- **OD-based** (`plate_reader_comparison/transform_search_od.py`): the same hypothesis space with
  a far better-powered assay (171 pairs reproducing at 0.62–0.82). **No credible hit.** The
  apparent best (`shift_rowmajor_-7`, ρ=0.433 at n=48) is *below* the ρ=0.464 that wrong
  transforms reach by artifact at that sample size in a split-half control where the answer is
  known to be `identity` — and in that control `identity` only ranks 2nd of 142, so the search's
  ranking is biased toward transforms that shrink the overlap. Reported as a genuine negative,
  not a discovery.

## 4. Partial recovery: ~half the plate, validated non-circularly

`20260721/analysis/corroborated_db_mapping/` re-identifies each well by 16S against
`corroborated_db_filtered_min5_edited.fasta`.

- **0 of 86 labels are their own best hit.** 83/86 match *something* at ≥0.99.
- Confidence tiers: 43 high (≥0.99, ≥10 bp margin), 19 medium (3–10 bp), 19 `low_16S_cannot_separate`
  (<3 bp), 2 contradicted by the mono-well consensus, 3 no match.
- **62 of 86 recommended**; 46 of those reach the genomic tables.

**Validated against the interaction labels, which never entered the mapping** (`validate_recovery.py`) —
this avoids the circularity `qc_recovery.py` warns about, where scoring a 16S-derived recovery by
16S-vs-genome correlation answers itself:

| cv_strain (both strains held out) | R² | ρ | winner called |
|---|---|---|---|
| recovered mapping (46 labels, 560 pairs) | **+0.126** | **0.413** | 74.2% |
| 5 permutations of it | −0.465 ± 0.191 | ~0 | ~50% |
| original `Well_souce_plate` mapping | −0.210 | 0.043 | 54.3% |

z = +3.10, beats all five permutations. So the recovered genomes carry real information —
**but read that as "phylogenetically close enough to predict", not "strain-level correct"**: the
recovery picks by 16S similarity, so even a wrong pick lands on a near relative. The 19
degenerate rows are unfixable by any amount of 16S care — six labels all claim `N19`, each
winning by 1 bp, because `D3` and `D13` sit 1 bp from `N19`.

This partially overturns `qc_recovery.py`'s "not recoverable" conclusion — that module measured
against `genome_16S` (the weak source, median margin ~2.7 bp); ONT-vs-ONT against
`corroborated_db` gives 10.5 bp.

### Where the labels moved

![movement](../20260721/analysis/corroborated_db_mapping/outputs/figures/p01_label_movement_plate.png)

Diffuse, with no displacement reused more than twice and a distance distribution sitting on top
of a random-relabelling null — which is what the transform-search negative looks like drawn out.

## 5. The remaining hypothesis

Everything on paper is consistent; the physical contents are not; the discrepancy has no geometry.
The **one step nothing records is preculture → Echo source plate** — the plate reader saw the four
96-well precultures, the files describe the intended mapping, and nothing independently verifies
what physically went onto the Echo deck. A non-geometric error there (a plate in the wrong deck
position, a re-pipetted intermediate) fits every observation.

Confirming it needs evidence outside the sequencing: lab records of which worklist the 20260721
pick used, whether a re-arrayed copy of the collection exists, or shotgun sequencing of the source
wells — the last being the only route that resolves the degenerate clusters, since 16S provably
cannot.

## 6. Is 20260721 worth anything?

Yes, with care. Its interaction data is **internally** strong and self-consistent: a per-strain
strength model using only the well label as an identity explains R²=0.748 / ρ=0.876 / 94.9% winner
accuracy of its own labels, and its cross-plate replicate ceiling (0.812 under the unswapped
layout) is close to 20260630's 0.869. The experiment was executed consistently — only the mapping
from well to organism is broken. So it supports analyses of interaction *structure* that never
name an organism, and nothing that does.

A useful by-product: **the replicate ceiling is a layout diagnostic.** 0.812 under
`strain_layout_20260721.csv` vs 0.687 under `_plate1_2_swapped` — independently confirming the
unswapped layout describes the physical wells (it also matches `echo_strains`), and that the swap
is a barcode correction, not a statement about contents.

## Pointers

| question | where |
|---|---|
| Source plates the same? | `plate_reader_comparison/` — `growth_compare.py`, `g02_offset_control.csv` |
| Destination contents differ? | `plate_reader_comparison/experiment_od_overlap.py` |
| Plate geometry ruled out | `strain_identity_qc/qc_layout.py`; `plate_reader_comparison/transform_search_od.py` |
| Contamination ruled out | `20260721/analysis/contamination_scan/` |
| The recovered mapping | `20260721/analysis/corroborated_db_mapping/outputs/m02_*.csv` |
| Genome join fails | `20260721/analysis/genomic_ml_join_test/` |
| What each plate-reader file is | `Karl_2026/Karl_2026*_OD/FILE_CONTENTS.csv` |
