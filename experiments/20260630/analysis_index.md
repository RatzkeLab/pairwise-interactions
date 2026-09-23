# Analyses using experiment 20260630

One row per analysis run that touched this experiment. Read this instead of opening the
analysis folders; open them only for method detail or for code worth reusing.

| analysis | date | what it asked | what it found |
|---|---|---|---|
| [`analyses/20260710`](../../analyses/20260710) | 2026-07 → 08 | Reference DB, read QC, relative abundance, first genomic ML | Built the 84-strain consensus and the interaction/hierarchy pipeline. Genomic ML: `cv_pair` is inflated by memorised per-strain competitiveness; `cv_strain` is the real number. |
| [`analyses/20260827`](../../analyses/20260827) | 2026-08 → 09 | Cross-experiment QC: is 20260721 the same plate as this one? | OD proves both experiments were picked from one source-plate layout, so the mix-up is downstream of the stock. Plate-reader spectra beat the genome baseline as ML features. |
| [`analyses/20260921`](../../analyses/20260921) | 2026-09 | Used as the **positive control** for 20260907's identity check | 58.7% label support against the collection reference, against 0% for 20260721. |

## This experiment's standing

**It is the reference experiment.** As of 2026-09-22 it is one of only two whose well
labels join correctly to the genomic tables — the join gate gives ρ=+0.363, z=+6.0 against
a permutation null. 20260721 fails that gate decisively; 20260907 now passes it at
ρ=+0.387, z=+8.2.

Best `cv_strain` result from `analyses/20260710`: R² 0.321 (two-stage ridge) over
74 genomes and 1,465 modelling pairs.

## What would change if the demux were re-run

Everything depth-limited: more wells clearing the read threshold, more strains getting a
corroborated consensus, and lower variance on the relative-abundance scores. The labels
and the join would not change. See [demux_runs.md](demux_runs.md).
