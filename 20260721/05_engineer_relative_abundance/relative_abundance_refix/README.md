# relative_abundance_refix — 20260721 re-run under the MIN_RESOLVABLE_BP fix

**Why:** 20260630's relative abundance was regenerated when the ambiguous-margin bug was fixed
(read assignment now decided in raw bp, ambiguous only on a genuine tie; pairs closer than
`MIN_RESOLVABLE_BP` flagged `below_resolution_limit`). 20260721's outputs still predated it.

Run: `rerun_and_compare.py` — writes here and **leaves `../relative_abundance/outputs/`
untouched**, so old and new can be diffed rather than one silently replacing the other.
Stages `r01`–`r03` only; `r04`/`r05` (competitiveness, Bradley-Terry) were not re-run.

| | pre-fix | post-fix |
|---|---|---|
| high-uncertainty pairs | 494 | **31** |
| usable pairs | 1625 | **2088** |
| labels changed >0.01 | — | **1012 of 2119** (median 0.47, max 7.78) |

`check_upstream_impact.py` answers the companion question — **how far upstream does the fix
reach? Nowhere.** The 16S consensus is built de novo from within-well reads and never assigns a
read to a reference, so the set of 16S sequences cannot change. `mapping_validation` never had
the bug either; its only exposure is normalized-vs-bp argmax, which flips 10/83105 reads
(20260630) and 42/110282 (20260721).
