# plate_reader_comparison — were both experiments inoculated from the same frozen stock?

**Answer: yes.** Independent of 16S, so not limited by the degeneracy that caps every
sequence-based answer.

Run in order: `label_files.py` → `reconstruct.py` → `growth_compare.py` →
`experiment_od_overlap.py` → `transform_search_od.py`.

| finding | evidence |
|---|---|
| **Source plates are the same physical layout** | per-well max OD ρ = **0.698**, at the same-plate ceiling (0.73–0.78) and above a within-experiment cross-media comparison (0.44) |
| ...and that is shared layout, not plate geometry | ρ ≤ **0.011** at *any* one- or two-well offset (`g02_offset_control.csv`) |
| **Destination wells nevertheless disagree** | the same nominal pair reaches ρ = **0.084** across experiments, against a replicate ceiling of 0.62–0.82 |
| No plate-handling transform explains it | best candidate fails its own split-half control (TRAPS §5) |

**`label_files.py` writes `FILE_CONTENTS.csv` into each reader folder** — a reconstruction of
what every plate-reader file contains. Essential, because **filenames are save times and the
internal `Date:` header is the run start**, differing by up to 44 h (TRAPS §6).

Only `max_OD` is a usable cross-plate phenotype here; `delta_OD` and growth rate fail their own
positive control because most runs start at or near stationary.
