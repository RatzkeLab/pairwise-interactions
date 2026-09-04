# QC summary — both experiments, taking the well labels at face value

**Scope.** This document reports what the data says when the plate-well→strain labels are
trusted and joined naively, exactly as a first-pass analysis would. It is deliberately the
*naive* view. What actually happened to 20260721, and everything done to reconstruct it, is in
[QC_02_RECONSTRUCTION.md](QC_02_RECONSTRUCTION.md).

Regenerate everything here with `python build_qc_matrix.py && python qc_reference_diagnostics.py
&& python make_figures.py` (~10 min; needs the `karl_seq_analysis` conda env). Tables in `outputs/`, figures in
`outputs/figures/`.

---

## 1. Experiment design and plate selection

Two ONT full-length 16S amplicon runs of the same 384-well strain collection, grown pairwise in
384-well plates, 30 destination plates each.

| | 20260630 (exp 1) | 20260721 (exp 2) |
|---|---|---|
| destination plates | 30 | 30 |
| layout wells | 9120 | 9120 |
| wells with >5 reads | 2812 | 2644 |
| distinct strain labels | 84 (pair consensus) | 81 (pair consensus) |

**Primer/barcode plates: odds (1–19) for exp 1, evens (2–20) for exp 2.**

> ⚠️ **This is not recorded anywhere in the repository.** `primer_layout_*.csv`,
> `minibar_primers_*.tsv` and all 30 `echo_primers_*_part_*.csv` files are **byte-identical
> between the two experiments** — they specify source *wells* (e.g. `E15`), never which physical
> primer plate those wells sat on. The odds/evens split is a bench fact from the experimenter's
> notes only. One known exception, also from notes: plate-1 primers were accidentally shot onto
> what were intended to be plate-2 strains, which is what `strain_layout_20260721_plate1_2_swapped.csv`
> patches.

## 2. Sequencing depth

Depth is essentially identical between runs, so nothing downstream can be blamed on one run
being more deeply sequenced.

| | wells | mean reads | median | min | max |
|---|---|---|---|---|---|
| 20260630 | 2812 | **41.9** | 40 | 6 | 133 |
| 20260721 | 2644 | **42.8** | 41 | 6 | 157 |

![read depth](outputs/figures/qc_f2_read_depth.png)

Depth matters for the *relative-abundance* label: below ~20 reads a replicate pair only
reproduces at ρ≈0.72, versus ρ≈0.91 above 45 reads. It does **not** limit the ML models — see
`technical_replicates/REPORT.md`.

## 3. Do the 16S consensus sequences support the labels?

For every (experiment × consensus build × reference) the same question: for well label *L*, does
the reference agree that *L* is that organism?

Three consensus builds, because they fail differently:

- **pair** — built from two-strain wells. Self-consistent by construction, so it is the one
  source that can be *confidently* wrong.
- **mono** — built from single-strain wells only. Immune to two-strain deconvolution error, and
  built here by the **same code for both experiments**, so it is the like-for-like comparison.
- **merged** — pair + mono combined (20260721 only; 20260630 has no merged build).

![label support](outputs/figures/qc_f1_label_support.png)

### % of testable labels supported (consensus matches its own label at ≥0.99 identity)

| experiment | build | corroborated_db (185) | corroborated_db_min5 (87) | genome_16S | other experiment |
|---|---|---|---|---|---|
| **20260630** | pair | **96.0%** (48/50) | 100% (19/19) | 35.4% (29/82) | 16.7% (3/18) |
| 20260630 | mono | **93.3%** (14/15) | 100% (7/7) | 21.7% (5/23) | 14.3% (1/7) |
| **20260721** | pair | **3.7%** (3/81) | 3.7% (3/81) | 2.0% (1/51) | 16.7% (3/18) |
| 20260721 | mono | **0%** (0/18) | 0% (0/18) | 0% (0/9) | 0% (0/3) |
| 20260721 | merged | **3.5%** (3/86) | 3.5% (3/86) | 1.9% (1/54) | 15.8% (3/19) |

**The headline: 20260630's labels are right, 20260721's are not.** 96% vs 3.7% against the same
external reference, on the same platform, at the same read depth. Mono wells agree with pair
wells within each experiment, so this is not a two-strain deconvolution artifact — 20260721's
single-strain wells are just as wrong as its pair wells.

### Cross-experiment: the two runs disagree with each other

The experiments share 18 well labels. **Only 3 of 18 are the same organism** (identity ≥0.99);
the median identity between the two runs' consensus for the *same label* is 0.904.

![definitions and cross-experiment](outputs/figures/qc_f3_definitions_and_crossexp.png)

So the two experiments do not even share a labelling namespace with each other — a fact that
holds regardless of which external reference you trust.

---

## 4. Two caveats that change how these numbers read

### 4a. "Supported" has two defensible definitions, and they differ a lot

- **label agrees** — the consensus matches its own label at ≥0.99 (this is what the table above
  and `strain_identity_qc` report)
- **label is closest** — the label is the single best match in the reference

| | label agrees | label is closest |
|---|---|---|
| 20260630 pair vs corroborated_db | **96.0%** | **46.0%** |
| 20260630 mono vs corroborated_db | 93.3% | 40.0% |

The 50-point gap is **not** wrong labels. It is 16S degeneracy: for ~25 of 20260630's labels a
near-twin in the reference is marginally closer than the label itself. Only 35 of 294 collection
strains have a unique 16S. Use "label agrees" for QC; do not use "label is closest" as an
identity test.

### 4b. There are two different corroborated databases, used inconsistently

| file | entries | used by |
|---|---|---|
| `corroborated_db.fasta` | 185 | `strain_identity_qc/` |
| `corroborated_db_filtered_min5.fasta` | 87 | `20260721/analysis/config.py`, the mapping work |
| `..._min4 / _min10 / _min15 / _filtered / _edited` | 94 / 76 / 65 / 87 / 87 | various |

`corroborated_db_filtered_min5_edited.fasta` is **content-identical** to the un-`_edited` version
(same 87 ids, zero sequence differences). The 87-entry subset overlaps only 19 of 20260630's 84
labels, which is why its "100%" is on a much smaller denominator than the 185-entry DB's 96%.
Always state which DB a number came from.

---

## 4b. How usable is each reference, really? (`qc_reference_diagnostics.py`)

The 0.99 cutoff is unfair to `genome_16S`: those sequences come from NGS assemblies while the
queries are ONT amplicon consensus, so the whole identity distribution is shifted down by
platform regardless of whether the reference names the right organism. Two threshold-free views
separate "offset" from "uninformative".

![reference diagnostics](outputs/figures/qc_f4_reference_diagnostics.png)

### Threshold sweep — % of labels agreeing, pair consensus

| reference | 0.95 | 0.97 | 0.98 | 0.99 | 1.00 |
|---|---|---|---|---|---|
| 20260630 corroborated_db(185) | 98.0 | 98.0 | 98.0 | **96.0** | 58.0 |
| 20260630 corroborated_db_min5(87) | 100 | 100 | 100 | **100** | 57.9 |
| 20260630 genome_16S | 68.3 | 42.7 | 40.2 | **35.4** | 28.0 |
| 20260721 corroborated_db(185) | 27.2 | 14.8 | 11.1 | **3.7** | 0.0 |
| 20260721 genome_16S | 21.6 | 7.8 | 2.0 | **2.0** | 0.0 |

Relaxing to 0.95 does rescue `genome_16S` a lot (35% → 68%), confirming a real platform offset:
its median self-identity is **0.961** versus **1.000** for the corroborated DBs. It does *not*
rescue 20260721 (3.7% → 27.2% at a threshold so loose it is no longer an identity claim).

### Rank retrieval — is the correct label the *closest* entry? (threshold-free)

| experiment / build / reference | entries | median self-identity | **top-1** | top-5 | median rank |
|---|---|---|---|---|---|
| 20260630 pair → corroborated_db(185) | 185 | 1.000 | **50.0%** | 78.0% | 1.5 |
| 20260630 pair → corroborated_db_min5(87) | 87 | 1.000 | **100%** | 100% | 1.0 |
| 20260630 pair → genome_16S | 233 | 0.961 | **23.2%** | 41.5% | **22** |
| 20260721 pair → corroborated_db(185) | 185 | 0.823 | 0.0% | 3.7% | 106 |
| 20260721 pair → genome_16S | 233 | 0.824 | 0.0% | 3.9% | 116 |

**`genome_16S` is bad on both axes, not just offset.** Even ignoring thresholds entirely, the
correct label ranks **22nd of 233** on average. So its low support figure is not merely a cutoff
artifact — the reference genuinely cannot identify these strains, and should not be used as an
identity reference. (It is also multi-copy, 587 records for 294 strains; the best-matching copy
per strain is taken, which if anything flatters it.)

**`corroborated_db_min5` is the best identification reference despite being the smallest.**
100% top-1 for 20260630. Filtering to ≥5 support removes precisely the near-twin entries that
create ambiguity — smaller, but unambiguous. The trade is coverage: it can only speak to 19 of
20260630's 84 labels.

### The ceiling: 16S can only do so well in this collection

Mapping each experiment's **mono** consensus onto its **own pair** consensus set — same
experiment, same platform, same pipeline, different well type, so no cross-reference confound at
all:

| | n | top-1 | median rank | agree at ≥0.99 |
|---|---|---|---|---|
| 20260630 mono → own pair | 21 | **52.4%** | 1.0 | 81.0% |
| 20260721 mono → own pair | 15 | **66.7%** | 1.0 | 80.0% |

Two things follow.

**A practical ceiling of ~52–67% top-1.** Even querying an experiment against *itself*, 16S puts
the right label first only about half to two-thirds of the time — because only 35 of 294
collection strains have a unique 16S. So `corroborated_db`'s 50% top-1 for 20260630 is *at* that
ceiling: the reference is as good as 16S allows, and the shortfall is the marker, not the DB.

**20260721 is internally consistent.** Its mono and pair wells agree with each other about what
is in each well (66.7% top-1, 80% at ≥0.99) — as well as 20260630 does. The experiment agrees
with itself and disagrees only with the outside world, which is exactly the signature of a
correct experiment on a mislabelled plate.

---

## 5. What this justifies

**20260630 is usable, but "76/76 strains join" is a lookup, not a verification.** Every well
label has a row in `mapping_384_well_plate_collection.csv`; that says nothing about whether the
genome in that row is the organism in that well. Per strain (`qc06_20260630_per_strain_status.csv`):

| status | n of 83 |
|---|---|
| corroborated_db confirms the label (≥0.99) | **48** |
| ...of which genome_16S *also* confirms | 17 |
| corroborated_db **contradicts** the label | **2** |
| untestable — absent from every reference | **33** |

The genome_16S disagreements are most likely genome_16S's fault, not the join's: it is the odd
one out 16 times against corroborated_db's once (`s04_three_way_concordance`), and §4b shows it
cannot identify strains at all. But "probably the reference's fault" is not verification.

**Restricting to the 46 corroborated-confirmed strains makes the relative-abundance model
better, not worse** (`qc07_strain_tier_sensitivity.csv`):

| target | tier | strains | pairs | R² | ρ |
|---|---|---|---|---|---|
| relative abundance | all | 74 | 1465 | 0.321 | 0.576 |
| relative abundance | **corroborated-confirmed** | 46 | 568 | **0.364** | **0.721** |
| yield | all | 88 | 3828 | 0.233 | 0.568 |
| yield | corroborated-confirmed | 46 | 1035 | −0.132 | 0.422 |

For relative abundance this is a real gain, not a sample-size effect: that target's strain
learning curve is flat past ~45 strains, so losing 74→46 costs nothing by itself. The unverified
strains were adding noise. **Default the relative-abundance analysis to this subset.**

For yield the comparison is **confounded and uninterpretable**: yield's strain curve is steep
(R² −0.357 at 45 strains, +0.265 at 88), so halving the strain count is expected to hurt
regardless of verification. It says nothing about whether yield's unverified strains are bad.

The strictest tier (17 strains confirmed by two references) is too small to fit at all, so
"confirmed twice" is not an achievable standard with this collection. The **2 contradicted
strains should be dropped outright**; the **33 untestable ones are neither confirmed nor
refuted** and no reanalysis of existing data can settle them.

The aggregate join validation (16S divergence vs KO-profile divergence, ρ=+0.363, z=+6.0 against
a permuted-assignment null) still holds — but it is an *aggregate* test: the mapping is right on
average, which is not the same as any particular strain being right.

Its mono and pair wells agree, and its labels are supported by an independent external reference
at 96% where testable.

**20260721 is not usable for anything that depends on strain identity.** Its wells contain real,
clean, single organisms — they are simply not the ones the labels claim. Interaction *outcomes*
within 20260721 remain internally valid (a per-strain strength model explains R²=0.75 of its own
labels), so the run is not worthless; only the mapping from well to organism is broken.

**Do not merge the two experiments into one interaction network.**

## Files

| file | contents |
|---|---|
| `outputs/qc01_consensus_vs_reference_detail.csv` | one row per (experiment, build, reference, label) with best match and identities |
| `outputs/qc02_consensus_vs_reference_summary.csv` | verdict counts per cell |
| `outputs/qc02b_two_definitions.csv` | the same cells under both definitions of "supported" |
| `outputs/qc03_reference_ranks.csv` | per label: self-identity, best match, and the RANK of the correct label |
| `outputs/qc04_rank_retrieval.csv` | top-1/3/5 and median rank per reference |
| `outputs/qc05_threshold_sweep.csv` | % agreeing at 0.95 ... 1.00 |
| `outputs/qc06_20260630_per_strain_status.csv` | per-strain verification status for 20260630 |
| `outputs/qc07_strain_tier_sensitivity.csv` | do the ML conclusions survive stricter strain sets |
| `build_qc_matrix.py`, `qc_reference_diagnostics.py`, `qc_sensitivity_strain_tiers.py`, `make_figures.py` | regenerate the above |
