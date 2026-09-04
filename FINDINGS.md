# FINDINGS — what is established, what is suggestive, what is untested

Scope note: this consolidates the **modelling** side. The QC/label side is documented in
`qc_summary/QC_01_NAIVE_MAPPING.md`, `QC_02_RECONSTRUCTION.md` and `QC_03_GENOMIC_FEATURES.md`,
written by parallel sessions; findings 1–3 below summarise their conclusions and point there
rather than restating the evidence. Read `TRAPS.md` alongside this.

Status labels: **[established]** multiple independent lines agree · **[suggestive]** one line,
or effect within noise · **[untested]** not attempted.

---

## The gate: only 20260630 is usable

**1. 20260721's strain labels do not name the organisms in its wells. [established]**

Three independent lines: 16S-vs-KO divergence correlates for 20260630 (ρ=+0.363, z=+6.0 vs a
permuted null) and is flat for 20260721 (ρ=−0.046, z=−0.7); the two experiments agree on only
3/18 shared labels by their own consensus; and a naive genomic join predicts 20260721's
interactions no better than a randomly permuted genome mapping (z=−0.06 under `cv_strain`).
→ `qc_summary/QC_01`, `20260721/analysis/genomic_ml_join_test/`

**2. The frozen stock and source plate were the same for both experiments. [established]**

Per-well OD of the source plates correlates across experiments at **ρ=0.698** — at the
same-plate ceiling (0.73–0.78) and above a within-experiment cross-media comparison (0.44).
Decisive control: the correlation is ≤0.011 at *any* one- or two-well offset, so it is shared
layout, not plate geometry.
→ `plate_reader_comparison/`

**3. The error is downstream of the stock, and has no geometric structure. [established]**

Destination wells disagree: the same nominal strain pair reaches uncorrelated OD across
experiments (**ρ=0.084**) against a within-experiment replicate ceiling of 0.62–0.82. Yet
`echo_strains` matches `strain_layout` for all 9120 wells in both experiments, contamination is
not elevated (18.4% vs 17.8% excess-organism rate), the poorly-matching reads are diffuse rather
than one foreign organism, and no plate-handling transform explains it (335 tested on 16S, and a
better-powered OD version whose apparent hit failed its own control — TRAPS §5).

The one step nothing covers is **preculture → Echo source plate**. A non-geometric error there
(wrong plate in a deck position, a re-pipetted intermediate) fits every observation.
→ `qc_summary/QC_02`, `20260721/analysis/contamination_scan/`

**4. ~46 of 86 of 20260721's labels are recoverable by 16S. [suggestive]**

Validated non-circularly — the mapping is built from 16S alone and scored against interaction
labels that never entered its construction: `cv_strain` ρ=0.413 vs permuted 0.
Caveat: recovery assigns *phylogenetically close* genomes, so this shows the recovered genomes
carry real information, **not** that the strains are correctly identified. 19 labels sit in 16S
clusters that cannot be split at all (six labels claim `N19` within 1 bp).
→ `20260721/analysis/corroborated_db_mapping/`

---

## Modelling on 20260630

**5. Gene content predicts interaction outcome on unseen strains. [established]**

| target | best `cv_strain` | 16S-only control | shuffled |
|---|---|---|---|
| relative abundance (who wins) | ρ **0.576**, R² 0.32 | ρ 0.458 | ≈0 |
| total yield (OD) | ρ **0.568**, R² 0.23 | ρ 0.345 | ≈0 |

Both clear their taxonomy control, so this is about gene content, not just the taxon label.
Winner called correctly ~80% of the time.
→ `20260630/analysis/genomic_ml/`, `genomic_ml_yield/`

**6. `cv_pair` numbers are inflated and must not be quoted alone. [established]** See TRAPS §2.

**7. Yield and competitive outcome are different questions with different bottlenecks. [established]**

Yield is symmetric, relative abundance antisymmetric; the modelling machinery does not transfer
(no antisymmetrisation, `|x_a−x_b|` instead of `x_a−x_b`, additive yield instead of
Bradley-Terry strength). Yield covers **3828 pairs / 88 strains** vs 1465 / 74, because OD is
independent of 16S resolvability — every pair dropped as `high_uncertainty` is measurable.

Notably, for **yield** 16S alone nearly matches KO (R² 0.411 vs 0.406); for **relative
abundance** KO crushes it (0.44 vs 0.10). Yield is more phylogenetically determined.

**8. Six levers do not move the relative-abundance model. [established]**

| lever | result |
|---|---|
| more replicates | +0.006 available; labels already at r₁=0.965 |
| more read depth | flat from 10 → 39 reads/well (downsampling test) |
| more genomes | flat past ~45 strains |
| annotation scheme | 6 schemes, 57–13154 features, all in a narrow band; KO on top |
| supervised feature selection | ties prevalence-matched random |
| constructed metabolic features | ρ 0.505 vs KO's 0.576 |

ρ ≈ 0.58 looks like a property of the genotype→phenotype relationship at this scale.
→ `20260630/analysis/genomic_ml/replicate_value_relabund.py`, `feature_sweep/`

**9. For yield, more genomes IS the high-value experiment. [established]**

The strain-learning curve is still climbing steeply at 88 genomes (R² 0.049 → 0.111 → 0.265 for
60 → 75 → 88). Unlike relative abundance, yield is not saturated.
→ `20260630/analysis/genomic_ml_yield/replicate_value.py`

**10. The collaborator's 221 panX GOs are phylogenetic markers, not mechanism. [established]**

Re-running his full four-stage pipeline *inside* each fold: nested selection ρ **0.529** vs
prevalence-matched random ρ **0.529**. Across five folds the selections share Jaccard 0.08–0.14
and **zero genes are selected in all five** (674 distinct genes, 35 overlapping his published
list) — yet all perform identically. Confirms his own `pipeline_notes.txt`.
→ `20260630/analysis/feature_sweep/nested_selection.py`

**11. Constructed metabolic features are compact but not better. [suggestive]**

23 interpretable pair features (CAZy degradative overlap = resource competition, KEGG module
complementarity, biosynthesis-module cross-feeding) reach ρ 0.505 vs KO's 0.576 — ~88% of the
performance from 0.5% of the inputs. Better for *explaining* a result, not for predicting one.
The cross-feeding block contributes +0.034/+0.039 after the module list was corrected.
→ `shared_pipelines/metabolic_features.py`

**12. Plate-reader optics add to the genome. [see that module's own README]**
→ `20260630/analysis/genomic_ml_plate/README.md` (parallel session; not independently verified here)

---

## Untested

- **Flux-based metabolic features** — draft GEMs (CarveMe/gapseq) + SMETANA MRO/MIP. The BiGG
  table cannot supply these (TRAPS §12), so the hypothesis is genuinely open, not refuted.
- **20260721 with a recovered pick list** — if the true source→destination mapping is found,
  `validate_strain_join()` re-tests it in seconds and ~2100 labelled pairs become available.
- **More genomes** — the one change finding 9 predicts will help.
