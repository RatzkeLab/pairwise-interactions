# technical_replicates — the reproducibility baseline for both assays

**Why:** a cross-experiment correlation of 0.08 means one thing if replicates of the same pair
reach 0.9 and something else entirely if they only reach 0.3. Every comparison elsewhere is read
against these numbers.

Run: `replicate_report.py`, then `spectra_examples.py` (figures), `od_target_ml.py` (feasibility
test that led to `genomic_ml_yield`). Full write-up in **`REPORT.md`**.

These are **cross-plate** replicates — 96% of pairs with ≥2 wells have them on different
destination plates — so this measures the whole pipeline, not well-to-well noise.

| | plate reader OD600 | sequencing |
|---|---|---|
| 20260630 | ρ 0.869 | **ρ 0.928** |
| 20260721 | ρ 0.686 | ρ 0.825 |

Sequencing is *more* reproducible than the plate reader, on matched pairs, at ~40 reads/well.
Two independent modalities agree that **20260721 is the noisier run**.

**Wavelength barely matters:** across all 61 channels (350–950 nm) replicate agreement spans
only 0.835–0.878; 410 nm is best but beats OD600 by +0.009. No reason to switch.

`s01`/`s02`/`s03` figures show the individual spectra behind the summary statistics — including
what ρ = 0.08 looks like well-by-well.
