# QC summary — both experiments, taking the well labels at face value

**Scope.** This document reports what the data says when the plate-well→strain labels are
trusted and joined naively, exactly as a first-pass analysis would. It is deliberately the
*naive* view. What actually happened to 20260721, and everything done to reconstruct it, is in
[QC_02_RECONSTRUCTION.md](QC_02_RECONSTRUCTION.md).

Regenerate everything here with `python build_qc_matrix.py && python make_figures.py`
(~4 min; needs the `karl_seq_analysis` conda env). Tables in `outputs/`, figures in
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

## 5. What this justifies

**20260630 is usable.** Its labels are supported by an independent external reference at 96%,
its mono and pair wells agree, and it joins cleanly to the genomic tables (76/76 strains).

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
| `build_qc_matrix.py`, `make_figures.py` | regenerate the above |
