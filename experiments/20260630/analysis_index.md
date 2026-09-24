# Analyses using experiment 20260630

One row per analysis run that touched this experiment. Read this instead of opening the
analysis folders; open them only for method detail or for code worth reusing.

| analysis | date | what it asked | what it found |
|---|---|---|---|
| [`analyses/20260710`](../../analyses/20260710) | 2026-07 → 08 | Reference DB, read QC, relative abundance, first genomic ML | Built the 84-strain consensus and the interaction/hierarchy pipeline. Genomic ML: `cv_pair` is inflated by memorised per-strain competitiveness; `cv_strain` is the real number. |
| [`analyses/20260827`](../../analyses/20260827) | 2026-08 → 09 | Cross-experiment QC: is 20260721 the same plate as this one? | OD proves both experiments were picked from one source-plate layout, so the mix-up is downstream of the stock. Plate-reader spectra beat the genome baseline as ML features. |
| [`analyses/20260921`](../../analyses/20260921) | 2026-09 | Used as the **positive control** for 20260907's identity check | 58.7% label support against the collection reference, against 0% for 20260721. |
| [`analyses/20260710_increased_tolerance`](../../analyses/20260710_increased_tolerance) | 2026-09-24 | Re-run 04–06 on the corrected demux (4.18× reads), side by side with the old demux through identical code | Depth 39 → 166 reads/well; per-well relative abundance r=0.988, BT strengths ρ=0.998. **ML headline unchanged**: cv_strain ρ 0.576 → 0.573 (paired, n.s.); R² 0.321 → 0.306 because the log2 label's tails stretch with depth. Label noise was never the bottleneck. |

## This experiment's standing

**It is the reference experiment.** As of 2026-09-22 it is one of only two whose well
labels join correctly to the genomic tables — the join gate gives ρ=+0.363, z=+6.0 against
a permutation null. 20260721 fails that gate decisively; 20260907 now passes it at
ρ=+0.387, z=+8.2.

Best `cv_strain` result from `analyses/20260710`: R² 0.321 (two-stage ridge) over
74 genomes and 1,465 modelling pairs.

## The demux has been re-run (2026-09-24)

See `analyses/20260710_increased_tolerance/README.md`. As predicted: more depth, same labels,
same join. Wells with >5 reads went 2,812 → 2,851; usable pairs 1,557 → 1,580. Nothing
depth-limited in the ML moved, since the ML was never limited by read depth.
