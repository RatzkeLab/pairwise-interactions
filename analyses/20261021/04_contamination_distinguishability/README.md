# 04 — plate maps: contamination, overgrowth, and what 16S can and cannot tell apart

Six 384-well maps of the 20260907 source plate (the plate the community experiments are
inoculated from), each answering one question. Every well shows its **strain number**
(isolate) on top; the bold label is what the figure is about. Run `python make_figures.py`
(~15 s, env `karl_seq_analysis`) after stages 1–10, or the whole analysis with `../run_all.py`.

| figure | question | headline |
|---|---|---|
| `f01_overgrowth` | which wells hold a different organism than their genome, and what took over? | **61 of 294** wells with a genome. G01's organism (one clone, C001) took 9 wells; G02's took 5 |
| `f02_secondary_contamination` | which source wells hold a second organism next to the strain? | **24** on independent evidence; in **14** the genome confirms the original still dominates |
| `f03a_identical_16S` | which strains share an identical 16S? | **254 of 378**, in 72 groups |
| `f03b_indistinguishable_lt8bp` | which strains can a read not tell apart, at the measured limit? | 323 strains in 65 mutually indistinguishable groups |
| `f03c_indistinguishable_lt15bp` | … at the original 15 bp design rule? | 341 strains in 57 groups |
| `f04a_duplicate_isolates` | the 9 isolates placed in the collection more than once — do the copies still hold them? | 15 of 21 copies do; 4 hold something else; MYb57 is lost from both wells |
| `f04b_genome_clone_clusters` | which genomes on file are identical? | 298 genomes = 181 organisms; 45 clone clusters, none of them deliberate copies |
| `f05_genomic_use` | which wells can go into a genome-predicted community, and why not the rest? | 204 wells = 133 distinct genomes |

PNG (200 dpi) and PDF (vector) for each, in `outputs/`.

## What each map encodes, and the choices behind it

**Colour budget.** A plate map lets any two wells sit side by side, so at most three
categorical hues are used per figure (validated all-pairs: worst colour-blind ΔE 9.2,
normal-vision 24.0). More identities than three are carried by the label in the well, never by
extra hues. Magnitudes (group size) use one blue ramp, validated as an ordinal ramp.

**f01 — overgrowth.** Filled = the well's 16S is not in its genome (stage 9 `wrong`); colour =
the 16S group of what took over (G01 blue, G02 orange, everything else dark with its label).
Ringed = a *home* well: same 16S, and the genome there *is* that organism. Coloured by 16S group
rather than by genome because the 16S is what is known — G02's 16S is carried by two different
genomes (the C002 clone and strain 139), so "which genome took over" is not always answerable.
For G01 it is: one clone, 16 isolates (M87 … M128) at ≥ 99.97% ANI — what the overgrown wells
hold is known, which of the 16 was the donor is not. A white dot = the original strain is still
there as a minority (5 wells).

Spread is at most weakly local: overgrown G01 wells sit 3.2 wells from the nearest home vs 4.9
for random wells (p = 0.06), G02 3.9 vs 4.4 (p = 0.37), all invaders pooled p = 0.17 — in line
with the earlier finding that a few invasive organisms dominate and neighbour transfer is modest.

**f02 — secondary organisms.** A second organism is pinned on the *source* well only when seen
twice independently:
- *both 20260907 monos* (blue, 9) — the two monos sat on different destination plates;
- *20261026* (orange, 10) — a later run from the same source plate: the organism's reads track
  this strain across communities (`community_assembly_experiments/20261026/02_qc`, association test);
- *close haplotypes* (aqua, 5) — reads phase into haplotypes > 10 bp apart. Greedy read
  clustering only separates organisms > ~100 bp apart, so close co-occupants show only this way.

Ring = the genome on file confirms the original strain is the dominant organism (14 — the
"original genome still there, contaminant in the mix" case). Dashed + × = the dominant organism
is not the original either (6, overgrown; see f01). Dots = a second organism in one mono only
(164 wells): it may have entered on the destination plate, so it is shown, not counted.
Labels name the second organism by strain number, by 16S group when several strains share its
16S, or `unk` / `U###` when nothing in the collection matches.

**f03 — indistinguishable groups.** Each group is a *clique*: every pair inside is below the
threshold, so at most one member can share a community. Connected components were not used —
they chain (at < 15 bp one component holds 116 wells whose members are up to 40 bp apart). The
greedy cover is not unique; a dot marks a strain that is also below the threshold to a member of
another group (70 at < 8 bp). f03a uses byte-identical 16S, which is an exact equivalence.

**f04a — deliberate duplicates.** From the strain-ID reconciliation (`analyses/20261113`):
9 isolates in 21 wells; only the K/M-row copy has a genome (ringed). A copy *holds* the isolate if
its 16S is carried by that genome or identical to a verified copy. Copies that hold something else:
F24 (MYb56), M24 (MYb71), and both MYb57 copies (H2, M3 — neither matches the MYb57 genome).
MYb53's genome has no usable 16S, so its two (identical) copies cannot be checked.

**f04b — identical genomes.** Clone clusters at ≥ 99.9% ANI. Not deliberate copies (each of
those has one genome) but separate isolates of one organism; for genome-based models a cluster is
one strain. Dashed + × = a member whose well no longer holds the organism (37).

**f05 — design view.** `genomic_use` from stage 10: dark = usable, and the 16S pins the genome
down; light = usable, but a 16S-identical organism with a different genome exists. Excluded wells
show the first reason that applies (tier D, never grew, mixed, wrong genome, unverified, no genome).

## Tables

| file | contents |
|---|---|
| `outputs/p00_well_table.csv` | one row per well: everything the maps draw |
| `outputs/p00_secondary_organisms.csv` | every second-organism observation, with its evidence |
| `outputs/p00_indistinguishable_groups.csv` | clique groups at 8 and 15 bp, members, max internal distance |
| `outputs/f04a_duplicate_isolates.csv` | the 21 duplicate-isolate wells and what each holds |

## Code

`plate.py` (shared drawing + palette), `p00_well_table.py` (all joins and the three new
computations: invaders and home wells, secondary organisms, clique cover), one `f0*.py` per map,
`make_figures.py` to run them in order.
