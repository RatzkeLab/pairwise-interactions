# genomic_ml_plate — plate-reader features on top of the genomic model (20260630)

Logic: `shared_pipelines/genomic_ml_plate.py`. Driver: `run_genomic_ml_plate.py` (`--quick`
for a 1-repeat smoke run). Baseline: `shared_pipelines/genomic_ml.py` /
`../genomic_ml/`, whose design is reproduced verbatim as the `T0_genomic` arm.

Question: does the plate reader add anything to KEGG-KO gene content when predicting
`mean_log2_ratio_a_over_b` — who wins a co-culture well?

Same 1465 pairs / 74 strains as the baseline, same labels, same models, **same folds across
all arms**, so tiers are compared fold-paired (`p04`) rather than by reading two independent
summary rows against each other.

## Result

Under `cv_strain` (both strains of a test pair unseen — the honest regime), paired ΔR² vs
`T0_genomic`, consistent across all four model families:

| tier | what it adds | ΔR² (two_stage_ridge) | p |
|---|---|---|---|
| T1_genomic_precult | source-plate inoculum OD | −0.005 ± 0.022 | 0.85 |
| T2a_genomic_mono_od | monoculture OD600 | +0.010 ± 0.014 | 0.64 |
| **T2b_genomic_mono_spec** | monoculture spectrum SHAPE | **+0.119 ± 0.036** | **0.002** |
| **T2_genomic_mono** | both mono readouts | **+0.170 ± 0.049** | **0.003** |
| T3_genomic_mono_cocult | + the co-culture well itself | +0.170 ± 0.049 | 0.003 |

Best arm: `T2_genomic_mono` + `two_stage_ridge`, R² 0.321 → **0.491**, Spearman ρ 0.577 →
**0.710**, winner called correctly 79.6% → **87.1%** (Δ +0.075 ± 0.021, p = 0.006).

Four things this does and does not say:

1. **The gain is the monoculture SPECTRUM, not monoculture growth.** OD600 alone adds nothing
   (T2a). This is why the mono tier is split — a combined tier would have read as "strains
   that grow better win", which is not what the data says.
2. **The spectral signal is essentially one scalar.** `p02c` — the band profile is a clean
   step through 600 nm (positive across the visible, negative across the near-IR) and PC1
   holds 90% of the shape variance. That is a scattering-slope / cell-morphology phenotype,
   not rich pigment structure. PC2 (8%) and PC4 add some independent signal, so it is not
   literally one number, but do not oversell it as spectroscopy.
3. **The co-culture well adds nothing over the monocultures (T3 ≈ T2).** Spectral unmixing of
   the pair well against the two monoculture spectra reaches only R² 0.10 on its own
   (`plate_unmix_only`). *A plate reader cannot replace the sequencing here.*
4. **Plate features do not replace the genome.** `P2`/`P3` (plate only) land at R² 0.27–0.34,
   at or below the genomic baseline. The two are complementary.

## Controls that had to pass

- `PLATESHUFFLED_*` (same features, strain→plate-well assignment permuted): sits at or below
  the T0 baseline in every tier. The gain is not "90 extra numbers help a 74-strain problem".
- `p05` mono-well QC: mapping_validation puts the 74 mono wells at 25 confirmed / 4
  low-confidence / 45 not-assessed (**not-assessed ≠ wrong** — those wells were below the read
  cutoff, so nothing was checked). Dropping the 4 actively-suspect strains keeps the gain
  (cv_strain ΔR² +0.125, p = 0.049); confirmed-only (25 strains, cv_pair only, too small to
  hold strains out) also keeps it (+0.068, p = 0.002).
- `T0_genomic` reproduces `../genomic_ml/outputs/g03_cv_summary.csv` exactly
  (two_stage_ridge cv_strain R² 0.321, ridge_pca 0.225) — the arms really are the same code.

## Side finding: preculture density predicts who wins

`p02` — the single strongest feature in the whole screen is the **source-plate OD difference**
(NM_dense_02Jul r = +0.529, precult_4x96_30Jun r = +0.515), i.e. how dense each strain was in
the plate the Echo picked from. It is layout-specific, not plate geometry: permuting the
strain→source-well assignment gives 0.003 ± 0.110, z = +4.7, p < 0.0005 (`p02b`). The other
three source reads — including the recovery plate read twice at ρ 0.88, so not a reliability
problem — do not predict (z = 0.7–1.8), which points at the two dense reads being the actual
Echo source rather than at a generic strain-fitness effect.

Note the offset/roll columns in `p02b` are **not** a valid null on this collection: it has
~6-column periodicity, so a 6-column roll leaves each strain anti-correlated with itself
(`selfcorr` ≈ −0.26) and produces a large negative "shifted" r that is an artifact of the
shift, not evidence against the signal. Read those columns only where `selfcorr` ≈ 0 — where
it is, r ≈ 0, as required.

**This cannot be resolved into cause from the data here.** "Denser inoculum → head start →
wins" and "fitter strain → grows dense everywhere AND wins" both fit. It matters either way:
if it is inoculum density, part of the measured competitive hierarchy is a procedural artifact
of the pick step, and normalising inoculum OD would change the results. Worth a dedicated
experiment (same pairs, inoculum densities deliberately varied).

Curiously, T1 adds nothing to the *model* despite that univariate strength — the KO features
already carry it (gene content predicts per-strain strength at R² 0.41, rho 0.63).

## Outputs

```
p01_*   modeling table with plate features attached, and mono-well QC status
p02_*   univariate screen: every plate feature against the target, no model
p02b_*  source-plate specificity control (permutation null + offset diagnostics)
p02c_*  spectral band profile and the PCA structure behind it
p03_*   the tier x model x regime sweep, per fold / summarised / per prediction
p04_*   fold-paired deltas vs T0_genomic  <- the result
p05_*   mono-well QC sensitivity arms
figures/fp01..fp05
```
