# corroborated_db_mapping — which organisms are actually in 20260721's wells?

**Question:** the labels are wrong; can the wells be re-identified by 16S against
`corroborated_db_filtered_min5_edited.fasta`? (That file is byte-identical in content to the
un-`_edited` one: same 87 ids, zero sequence differences.)

Run: `build_mapping.py`, then `validate_recovery.py`, then `plot_movement.py`.

**Use the `recommended_target` column** of `m02_20260721_to_corroborated_db_mapping.csv` —
populated for 62 of 86 labels and deliberately blank elsewhere. **Zero labels are their own best
hit.**

| confidence | n | meaning |
|---|---|---|
| `high` | 43 | ≥0.99 identity, beats runner-up by ≥10 bp |
| `medium` | 19 | margin 3–10 bp |
| `low_16S_cannot_separate` | 19 | margin <3 bp — six labels all claim `N19` by 1 bp |
| `contradicted_by_mono_well` | 2 | decisive, but the cleaner mono consensus disagrees |
| `no_confident_match` | 3 | nothing at ≥0.99 |

**Validated non-circularly** (`validate_recovery.py`): built from 16S alone, scored against
*interaction labels* that never entered its construction — `cv_strain` ρ = 0.413, z = +3.1 vs
permuted. But recovery assigns *phylogenetically close* genomes, so this shows the recovered
genomes carry real information, **not** that strains are correctly identified.

A one-to-one Hungarian assignment **must** be constrained to ≥0.99-or-unassigned; unconstrained
it maximises the global total by pushing individual labels onto organisms they match at 0.73.
