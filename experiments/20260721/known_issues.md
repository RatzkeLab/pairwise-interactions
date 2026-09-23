# Known issues — experiment 20260721

Defects in the data or in the methods applied to it, and what they affect.

---

## Well labels do not identify the organisms — CONFIRMED, partially recoverable

**Found** 2026-08 → 09, confirmed repeatedly since. **Affects** everything keyed by well
label, including all genomic-feature joins. **Status** open; ~half the labels are
recoverable by 16S.

This is the headline fact about this experiment. Three independent lines of evidence:

- **The genomic join gate fails.** Two strains with near-identical 16S must have
  near-identical gene content, so 16S divergence should track KO-profile divergence across
  tested pairs. Here ρ=−0.046, z=−0.7 against a permutation null — indistinguishable from
  a random assignment of genomes to wells. 20260630 gives ρ=+0.363, z=+6.0; 20260907
  gives ρ=+0.387, z=+8.2.
- **Label support is zero.** Scored against the collection reference, the genome-derived
  16S and the corroborated DB, this experiment's consensus sequences match their own label
  in **0 of 188 decidable comparisons**. The other two experiments score 34–100%.
- **It does not even share a namespace with 20260630.** Of the 18 well labels the two
  experiments have in common, only 3 are the same organism.

**But the plate is the same stock.** `analyses/20260827` shows per-well OD correlates
across the two experiments' source-plate reads, so both were picked from one physical
384 collection. The mix-up happened downstream of the freezer, in how wells were tracked.

**About half the labels are recoverable** by matching this experiment's own 16S consensus
back to the collection — validated against interaction outcomes rather than circularly.
See `analyses/20260827`.

---

## <a name="minibar-barcode-file"></a>Malformed barcode file — OPEN, not re-run

**Found** 2026-09-21 while working on 20260907. **Affects** `20260730_demux` and every
analysis built on it. **Status** open.

minibar reads column 2 of the barcode TSV as the forward index; the file puts the whole
59 bp adapter+barcode+primer construct there, so minibar matches a 59-mer at the barcode
edit distance. Measured on 20260907: assignment went from 13.5% to **64.1%** of filtered
reads once corrected. This experiment assigns 15.6%.

Nothing already assigned is wrong — there are ~4× more reads available. Fix:
`analyses/20260921/02_demux_qc/fix_barcode_file.py`. Re-running matters more here than
elsewhere because the label recovery is depth-limited.

---

## Single-linkage clustering chains at depth — FIXED in shared code

**Found** 2026-09-22. **Affects** `analyses/20260804/04_qc` (`check_reads.ipynb`).
**Status** method fixed in `scripts/well_composition.py`; the old outputs were not
regenerated.

Single-linkage connected components merge the two strains of a pair well once enough
error-rich reads bridge the same/different distance modes, and it gets *worse* with more
data — on 200 pair wells of 20260907, wells resolving both strains went 27.0% → 10.0%
between 120 and 500 reads, while greedy seed clustering held at 47.0% → 45.0%.

**This experiment's two-strain counts are understated** by that notebook. The replacement
decides membership against a cluster seed rather than a neighbour and cannot chain.

---

## No true monoculture wells — DESIGN LIMITATION

This experiment has no clean single-strain wells, so a strain's consensus can only be
built from the majority argument across pair wells: the sequence that recurs while
partners change. That works, but it cannot distinguish a strain from a contaminant that
travels with it.

20260907 was designed to fix this — every strain has two monocultures on different plates,
and `scripts/strain_consensus.py` anchors on them, falling back to pair-well corroboration
only when the monocultures gave no reads.
