# TRAPS — ways this project produces confident wrong answers

Every entry below actually happened here and produced a plausible-looking result before being
caught. Read this before trusting any number, your own included.

The recurring shape: **a result that looks like a finding but is a property of how the question
was asked.** The fix is almost always a control that must fail — not a better model.

---

## 1. Strain names are not identifiers

Labels like `D13`, `N19`, `A11` are **plate-well coordinates, reused across unrelated plates,
experiments and years**. A name match is not an identity match anywhere in this project.

- `corroborated_db.fasta`: only 21/86 names for 20260721 were the same organism the experiment
  meant, versus 19/19 for 20260630.
- The two experiments share 18 well labels; comparing each one's own 16S consensus, **only 3 of
  18 are the same organism**.

**Rule:** never join on a strain name without an identity check that could fail.

## 2. `cv_pair` cannot detect a wrong genome mapping — only `cv_strain` can

Under `cv_pair` both strains of a test pair appear in training via their *other* pairs. A KO
vector then works as a unique per-strain fingerprint, and the model scores well by recognising
the strain, not by understanding the genome.

A **consistently wrong** mapping is still a bijection, so it supplies equally good fingerprints.
20260721's provably meaningless join scored **R² = 0.44, ρ = 0.65, 79% winner accuracy** under
`cv_pair`. Under `cv_strain` it scored ρ = 0.04 — chance.

**Rule:** quote `cv_strain`. Report `cv_pair` only alongside the no-genomics baseline, which
usually beats it.

## 3. A control that cannot fail is not a control

Two cases here, both nearly reported as findings:

- **Vacuous set membership.** "100% of 20260721's observed organisms are within the expected
  set" — true *by construction*: the reference database contained exactly those 86 labels plus
  one. There was nowhere else to match. The real test needed a wider reference (294 strains),
  which gave 23% not 100%.
- **Self-comparison.** A transform search "positive control" that correlated 20260630 against
  itself returned ρ = 1.0 for `identity` and for two other transforms. It confirmed nothing.

**Rule:** before running a control, state what result would falsify the hypothesis. If nothing
would, it is not a control.

## 4. Random controls must be matched on everything except the thing being tested

A uniform random draw of 221 panX GOs kept only 46–58 columns after the prevalence filter, while
the selected 221 kept 139. The comparison measured *feature count*, not information, and made a
selected list look far better than it was (ρ 0.510 vs 0.335–0.430).

Drawn **prevalence-matched**, the same comparison gave 0.552 vs 0.581 — the selection lost.

**Rule:** a random control must match the real thing in size, prevalence, and anything else that
affects the pipeline.

## 5. Model-selection searches are biased toward whatever shrinks the test set

The 335-transform plate-layout search scored each candidate on however many pairs survived that
transform. Transforms that shrank the overlap were scored on smaller, easier subsets.

In a split-half control where the answer was known to be `identity`, **identity ranked 2nd of
142**. Among transforms reaching n ≤ 60 — all wrong by construction — the best spurious ρ was
0.464, *above* the 0.433 "hit" that had looked significant against a naive family-wise threshold.

**Rule:** when candidates are scored on different sample sizes, the family-wise threshold must
account for that. Verify by running the search where you know the answer.

## 6. Plate-reader filenames are SAVE times, not run times

The internal `Date:`/`Time:` header is the run **start**; the filename timestamp is when the
file was written, up to **44 hours later**. Runs chain back-to-back, so the sequence is
recoverable — but matching a lab-notebook time to a filename is wrong.

Concretely: the note "13:30 read the NM OD" matches `Karl_20260702_134335` (started 13:31:50),
**not** `Karl_20260702_133125` (started 11:17). They are different plates, r = 0.22.

**Rule:** parse the internal header. See `FILE_CONTENTS.csv` in each reader folder.

## 7. Hand-written biological ID lists fail silently

A 47-entry KEGG biosynthesis-module list written from knowledge had **5 wrong entries**, found
by checking `https://rest.kegg.jp/list/module`:

| ID | actually is |
|---|---|
| M00027 | GABA shunt — catabolic |
| M00140 / M00141 | C1-unit interconversion (M00141 eukaryotic) |
| M00128 | ubiquinone biosynthesis, **eukaryotes** |
| M00868 | heme biosynthesis, **animals and fungi** |

The eukaryote-only entries are the dangerous ones: asked of bacteria they produce near-constant
columns that add plausible noise rather than an obvious error. Removing all five cut the
cross-feeding block's measured contribution roughly in half (+0.087 → +0.039 on yield).

**Rule:** verify identifier lists against the source database. Report results with and without
any block that depends on a hand-assembled list.

## 8. R² is not defined the same way in the two ML modules

- `genomic_ml.py` (relative abundance, **antisymmetric** target): `1 − SSE/Σy²`, referenced to
  *predicting no winner*.
- `genomic_ml_yield.py` (yield, **symmetric** target): `1 − SSE/Σ(y−ȳ)²`, conventional.

Here they happen to differ by <1% (the relative-abundance target has mean −0.23 against sd 3.69),
but the definitions are not interchangeable in general.

## 9. R² explodes where ρ stays sane

Ridge extrapolating onto held-out strains from few features or few genomes produces R² values
like **−2×10²⁸**. A single such fold destroys any mean, max or ranking it enters.

Also: the **16S control is rank-informative but scale-miscalibrated** — ρ ≈ 0.45 while
R² ≈ −2.5 to −4.6. A "gap over 16S" computed in R² measures the control blowing up, not the
features' merit.

**Rule:** clip or NaN implausible R², and use ρ for any comparison involving the 16S control.

## 10. Most wells contain two organisms by design

A first contamination scan reported "70% of wells show a minor component" — because **pair wells
contain two strains on purpose**. Split by `well_type` (mono expects 1 group, pair expects 2) and
the real excess rate is 18.4% vs 17.8% in the control experiment: no elevation at all.

## 11. Read identifiers differ between files

`04_reads_combined.fastq` headers are `<sample_id>__<read_id>`; the CSVs carry the bare
`read_id`. Matching on the bare id silently returns **nothing**, and a downstream analysis
reports an empty result rather than an error.

## 12. `BiGG_and_strains_table.csv` is not a metabolic network

9442 columns of `<bigg_model_id>.<gene_locus>` — 73 published models (mostly *E. coli*,
*Shigella*, *Klebsiella*) crossed with their gene loci. **Gene-orthology calls: no reactions, no
metabolites, no exchange fluxes, no stoichiometry.** Uptake overlap, byproduct exchange and
joint FBA cannot be computed from it. It is also redundant (one locus appears under up to 4
models) and E. coli-centric for a *Pseudomonas*/*Sphingomonas* collection — which is why it was
among the worst feature tables tested.

## 13. Different feature tables keep different strains

Each annotation table drops different genomes as unannotated, so the retained strain set — and
therefore the 16S baseline — moves between tables. `genomic_ml._Ctx` indexed the 16S matrix
blindly and raised `KeyError` as soon as the feature table changed; `genomic_ml_yield._Ctx`
intersected defensively and survived. Both now intersect.

**Rule:** compare each feature set against **its own** taxonomy control, never across rows.
`panX_full` posted the best absolute yield ρ (0.591 vs KO's 0.549) purely because its baseline
was higher; against its own control it buys *less* than KO.

---

## Controls that earned their keep

Reusable, and each one caught something here:

| control | catches |
|---|---|
| **permuted mapping** (shuffle which genome attaches to which well) | a consistent-but-wrong join; the sharp null that plain label-shuffling misses |
| **prevalence-matched random features** | feature-count handicaps masquerading as information |
| **positive control that must recover a known answer** | biased search procedures (§5) |
| **spatial offset scan** (shift a plate by 1–2 wells) | smooth plate geometry masquerading as shared layout |
| **split by well type** | design structure masquerading as contamination |
| **with/without a hand-assembled block** | wrong ID lists driving a conclusion |
| **replicate ceiling** | whether a low score is model failure or measurement noise — and, unexpectedly, whether a layout describes physical well contents (0.812 unswapped vs 0.687 swapped) |
