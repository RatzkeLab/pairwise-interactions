# Strain-identity QC — 20260721 vs 20260630

Why this exists: the genomic-ML work found that 20260721's strain labels do not correspond to
the genomes the mapping file assigns them, while 20260630's do. This folder tests that directly,
at sequence level, with the platform confound removed.

Run: `python run_qc.py` (~1 min) or `python run_qc.py --figs`.

```
qc_config.py     paths, per-experiment settings, both 20260721 layout variants
qc_sources.py    load/build the four 16S sources into one long table
qc_compare.py    agreement (do two sources agree about strain S?) and
                 attribution (if not S, then what?) — strand-, length- and N-aware
qc_recovery.py   the 16S-based rescue of 20260721, and why it does not work
qc_readassign.py calibrates the read-assignment resolution limit against mono-well ground
                 truth — sets MIN_RESOLVABLE_BP in shared_pipelines/relative_abundance.py
qc_layout.py     plate-handling hypotheses: 96->384 quadrant mix-ups, backwards plate,
                 pick-list off-by-N — 335 candidate transforms, family-wise corrected
qc_figures.py    four figures
run_qc.py        driver
```

## The four sources

| source | origin | failure mode |
|---|---|---|
| `genome_16S` | 16S pulled from the NGS assemblies (449 records / 294 strains, multi-copy) | NGS assembly of 16S is unreliable — this is the weak source, quantified below |
| `corroborated_db` | external ONT db from prior corroborated experiments | independent of both experiments |
| `mono_consensus` | this experiment's mono wells (one strain per well) | immune to two-strain deconvolution; built here for both experiments by the same code |
| `pair_consensus` | this experiment's pair wells, as the interaction pipeline uses | self-consistent by construction — the one source that can be confidently wrong |

## Findings

**1. 20260721's labels are wrong; 20260630's are right.** Against `corroborated_db` (ONT, same
platform — so the NGS-vs-ONT caveat cannot explain it):

| | labels supported | labels contradicted |
|---|---|---|
| 20260630 | **48/50 (96%)**, median identity 1.0000 | 2 |
| 20260721 | **3/81 (3.7%)**, median identity 0.8228 | **75** |

**2. It is a mix-up, not contamination or bad sequencing.** 20260721's wells contain real,
clean collection organisms — 78/81 match *something* in `corroborated_db` at ≥0.99 — they are
just not what the labels claim. Its mono and pair wells agree with each other (80%), and where
both exist for a label they name the **same** organism 11/14 times. So each source well
consistently contained one wrong-but-real strain: the error happened upstream of the wells, at
source-plate identity, not per-well during the run.

**3. The three-way concordance names the culprit.** On strains present in all three sources:

| | all agree | experiment is odd | genome_16S is odd |
|---|---|---|---|
| 20260721 (n=62) | 1 | **34** | 1 |
| 20260630 (n=49) | 17 | 0 | **16** |

Two independent references agree with each other and disagree with 20260721. For 20260630 the
odd source is `genome_16S`.

**4. Your collaborator's NGS-vs-ONT warning is real, and now quantified.** `corroborated_db` vs
`genome_16S` — two references, no experiment involved — agree for only **65%** of 142 strains.
So genome-derived 16S is a poor identity reference, and 20260630's mediocre 35% agreement with
it is a property of *that source*, not of the experiment. **Never use `genome_16S` as the
identity gate.** Use `corroborated_db`.

**5. The mix-up is not a geometric transform.** No dominant row or column shift (top joint
(Δrow, Δcol) accounts for 2 of 78), not bijective (78 labels → 63 targets), zero self-maps. So
it is not a plate rotation, row/column offset, or transpose. Note 20260721 already ships a
patched layout (`..._plate1_2_swapped.csv`, 616 wells changed on plates 1–2) — a swap was found
before; this is a different and larger problem.

**5b. No plate-handling mistake explains it either.** 335 candidate transforms tested at 16S-
group level, with 20260630 as a positive control:

| | best transform | z | family-wise p |
|---|---|---|---|
| 20260630 (control) | **`identity`** — 88.6% explained | 8.18 | **0.0** |
| 20260721 | `quad_block_2103+rot180_within_96` — 60.7% | 3.01 | **0.234** (n.s.) |

The control recovering `identity` decisively is what makes the negative trustworthy. Tested and
rejected: all 24 permutations of the four 96-well quadrants under **both** recombination
conventions (interleaved checkerboard and spatial block), 180° rotation of the 384 plate
(backwards in the Echo), 180° rotation within the 96-well source, row-only and column-only
flips, pick-list off-by-N shifts (±25, row- and column-major), and every composition of a
quadrant permutation with an orientation error. Nothing survives correction for 335 tests.

Note this had to be asked at **16S-group** level: only 35 of 294 strains have no 16S near-twin,
which leaves 5 usable single-strain calls for 20260721 — far too few for 335 hypotheses.

**6. 20260721 cannot be rescued from 16S — and the number that suggests otherwise is circular.**
Re-identifying wells by observed 16S and re-joining takes the genome-vs-16S validation from
ρ = −0.05 to ρ = +0.91. That is an artifact: the rescue assigns genomes *by* 16S, then the test
asks whether 16S tracks genome content. The honest measure is resolving power:

- 75 labels collapse onto **39 distinct genomes**; one genome is claimed by 10 labels
- 46 of 75 labels sit in a collision; 18 calls are exact ties
- median call beats its runner-up by **2.7 bp in 1420**

Full-length 16S simply cannot separate this collection: at ≥0.99 identity the 294 strains form
just **53 groups, one of which holds 181 of them** (single-linkage, so that giant group is
partly chaining — but the practical consequence stands). Recovery needs a higher-resolution
marker (shotgun of the source wells) or physically identifying the plate that was used.

## 7. A separate bug this QC uncovered: the assay was throwing away resolvable pairs

`relative_abundance.classify_well_reads` called a read ambiguous when the **normalized**
distance margin fell below `AMBIGUOUS_MARGIN_THRESHOLD = 0.02` — about **28 bp** on a 1420 bp
read. Two references 10 bp apart cannot produce a margin that large, so every one of their reads
was filed "ambiguous" *by construction* and the pair was flagged unresolvable no matter how
clean the data. The tell: `mean_uncertainty_score` was exactly 1.000 in every bin below 20 bp.

Only the positions where the two references differ carry information — errors elsewhere add
equally to both distances and cancel — so the fix is to assign in raw bp and call a read
ambiguous only on a genuine tie, with a hard floor below which the pair is declared unresolvable
up front.

**Where the floor belongs, from mono-well ground truth** (`qc_readassign.py`, 68k read-level
comparisons; a substitution model predicts 99.97% at 5 bp and is wrong, because ONT error is
indel-dominated and concentrated in the homopolymers where near-identical 16S sequences differ):

| references apart | 1–2 | 3–5 | 6–10 | 11–20 | 21–40 | >40 |
|---|---|---|---|---|---|---|
| % reads assigned correctly | 39% | 59% | 73% | **92%** | 96% | 96–97% |

Ceiling is ~96–97%, not 100%, because 20260630's mono wells are not pure — by design they were
shot against wells believed to be no-growers, some of which grew. So `MIN_RESOLVABLE_BP = 10`.

**Effect on 20260630**, with no new sequencing:

| | before | after |
|---|---|---|
| usable pairs | 1090 | **1574** (+44%) |
| usable in the 10–20 bp band | 2/228 | **228/228** |
| pairs below the resolution limit | — | 103 (excluded by design) |
| Bradley-Terry pseudo-R² | 0.831 | **0.866** |
| directional consistency (DCI) | 0.670 | **0.823** |

The hierarchy fit *improves* with the recovered pairs, which is the evidence that they are real
signal rather than noise let in by a looser threshold. Previous outputs are preserved in
`20260630/analysis/relative_abundance/outputs_backup_pre_margin_fix/`.

## Bottom line

20260630 is sound and stays usable — and the layout search finds `identity` for it at
family-wise p = 0, so its labels are positively confirmed rather than merely un-refuted. 20260721's interaction measurements are real but their
strain identities are unknown and not reconstructable from these data. For future designs the
strain-selection cutoff should be **≥10 bp** minimum pairwise reference separation (the existing
`corroborated_db_filtered_min10.fasta`, 76 strains), not the 5 bp used for 20260721 — it should not be joined
to genomic tables, and its per-strain conclusions should not be pooled with 20260630's.
