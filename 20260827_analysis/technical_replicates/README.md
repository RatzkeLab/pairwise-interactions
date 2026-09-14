# technical_replicates — the reproducibility baseline for both assays

**Why:** a cross-experiment correlation of 0.08 means one thing if replicates of the same pair
reach 0.9 and something else entirely if they only reach 0.3. Every comparison elsewhere is read
against these numbers.

Run: `replicate_report.py`, then `spectra_examples.py` (figures), `od_target_ml.py` (feasibility
test that led to `genomic_ml_yield`), `ml_target_ceiling.py` (the ceiling on the genomic-ML
target, for comparing against cv_strain scores). Full write-up in **`REPORT.md`**.

These are **cross-plate** replicates — 96% of pairs with ≥2 wells have them on different
destination plates — so this measures the whole pipeline, not well-to-well noise.

| | plate reader OD600 | sequencing |
|---|---|---|
| 20260630 | ρ 0.869 (R² 0.759) | **ρ 0.928 (R² 0.875)** |
| 20260721 | ρ 0.686 (R² 0.302) | ρ 0.825 (R² 0.675) |

Sequencing is *more* reproducible than the plate reader, on matched pairs, at ~40 reads/well.
Two independent modalities agree that **20260721 is the noisier run**.

**Wavelength barely matters — in ρ.** Across all 61 channels (350–950 nm) rank agreement spans
only 0.835–0.878; 410 nm is best but beats OD600 by +0.009. In R² the same sweep spans
0.588–0.786 and the blue end beats OD600 by +0.078, so stay on OD600 for ranking but consider
350–410 nm if the OD value is used quantitatively. See REPORT.md §3.

**Ceiling on the ML target (§6):** cv_pair's R² 0.876 is already at the assay limit (0.897
conservative / 0.966 optimistic); cv_strain's 0.321 is nowhere near it, so what limits
generalization to unseen strains is not measurement noise.

`s01`/`s02`/`s03` figures show the individual spectra behind the summary statistics — including
what ρ = 0.08 looks like well-by-well.
