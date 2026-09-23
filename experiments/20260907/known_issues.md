# Known issues — experiment 20260907

Defects in the data or in the methods applied to it, and what they affect. Read this
before trusting any result derived from this experiment.

---

## <a name="minibar-barcode-file"></a>Malformed barcode file cost ~5× the reads — FIXED

**Found** 2026-09-21. **Affects** `20260907_demux` and every analysis built on it.
**Status** fixed for this experiment; the same bug affects 20260630 and 20260721, which
have not been re-run.

minibar reads column 2 of the barcode TSV as the forward index. The generated file put
the whole 59 bp construct there — 15 bp adapter + 24 bp barcode + 20 bp 16S primer — so
minibar matched a 59-mer within the barcode edit distance instead of the 24-mer barcode.
Same tolerance, 2.5× the length, and most reads failed.

84% of the discarded reads were clean full-length amplicons whose barcodes matched the
sample sheet at **median edit distance 0**. Fix: `analyses/20260921/02_demux_qc/fix_barcode_file.py`.
Calibration against the never-amplified frozen plates: `-e 4` gives 63.6% assignment at a
~1.2% false rate; `-e 6` sits exactly on the barcode set's minimum separation of 12 and
pushes the false rate to 5%.

A second, latent bug in the same place: `demultiplex.smk` reads config keys
`barcode_edit_distance`/`primer_edit_distance`, while every config in this project writes
`barcode_edit_dist`/`primer_edit_dist`. Those keys have never been read — the rule always
used its own defaults. They happened to match, so nothing broke, but changing them
silently does nothing.

---

## <a name="failed-primer-wells"></a>19 barcode oligos failed — NOT RECOVERABLE

**Found** 2026-09-21. **Affects** 480 wells (19% of the sequenced plate); 394 of them are
unusable. **Status** permanent data loss for this run; check the source plate before the
next one.

7 forward + 12 reverse barcodes of 192 carry **76% of all dropout**. Wells with two clean
barcodes drop out at 6.3%; wells touching a flagged barcode at 82%.

It is the oligo, not the sequence and not the strain — three independent checks:

- those barcodes were fine in 20260630 and 20260721 (20–26% dropout then, 81–85% now)
- their amplicons are absent from the discarded-read pile too, 10× depleted
- **OD says the cells were there**: flagged-barcode wells grew at 95.6% against 96.7%
  elsewhere (Fisher p=0.27), median OD 0.317 vs 0.325

Not depletion either — source wells drawn from *more* fail *less* (rho −0.22). The plate
was the same one, freshly thawed, so evaporation or freeze–thaw damage in specific wells
fits. Row H of the primer source plate holds 5 of the 12 bad reverse barcodes; G17/G19/G21,
K9/K11 and H12/H16/H18 are adjacent runs.

Affected barcodes with their source wells:
`analyses/20260921/02_demux_qc/outputs/d03_{fwd,rev}_barcode_forensics.csv`.

---

## 30 strains never grew — REAL, exclude or flag

**Found** 2026-09-22. **Affects** any analysis treating their wells as interactions.
**Status** open; kept in the modelling set by choice (see below).

Every strain has two monocultures on *different* plates, so independent well-level
failures would rarely hit both. 30 strains failed both — 9× the independent expectation
of 3.2. They also failed in preculture (median source OD 0.16 vs 1.35, p=3e-12).

They are **present but weak**, not absent: 26 of 30 still produce a corroborated
consensus, their mono wells carry a median 96 reads, and in pair wells they hold a median
5.1% of reads. So their wells are genuine interactions with a poor competitor, not
mislabelled monocultures, and they are currently left in the modelling set.

Two caveats if you use them: their mono wells drop out at 35% against 22%, and **36% of
their pair-well shares fall below 1%**, so the log2 ratio is floor-censored. Prefer a
bounded outcome (relative abundance) over raw log-ratio for those.

List: `analyses/20260921/02_demux_qc/outputs/d13_strains_no_growth.csv`.

---

## Composition partly remembers the inoculum — QUANTIFIED, not corrected

**Found** 2026-09-22. **Affects** every single-well interaction score. **Status** open;
measured, not corrected for.

Total *yield* is at carrying capacity: doubling the monoculture inoculum moves final OD by
1.8% of a doubling, and a 4× density span by −1.1% (p=0.54). But *composition* is a
different question, and the ratio-titration arm gives a **slope of +0.31** on 34 pairs —
about 31% of a change in starting ratio survives to the endpoint (p vs 0 = 6e-09, p vs
1 = 5e-09).

Likely an underestimate: median r² is 0.725 and among the 21 pairs fitting at r²>0.5 the
slope is +0.53, while the extreme levels compress, which flattens the fit. Strong-winner
pairs are *steeper* (+0.36) than near-neutral ones (+0.14).

So a single well's interaction score is partly reporting starting density. The first-pass
model does not correct for it; treat the slope as a sensitivity check.

---

## Method-level issues inherited from the analysis code

- **Single-linkage clustering chains at depth** — fixed 2026-09-22 in
  `scripts/well_composition.py`. Any composition result produced before that date
  understates two-strain wells. See `experiments/20260721/known_issues.md`.
- **`paths.py` pointed at the deleted `Link to Karl`** — fixed 2026-09-22. Genomic-ML
  runs before that fail on the first read rather than producing wrong numbers.
