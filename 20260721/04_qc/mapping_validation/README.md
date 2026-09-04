# mapping_validation — 20260721

Same pipeline as 20260630's, same outputs. **The verdict is the opposite.**

Only **3.7%** of testable labels are supported by this experiment's own 16S, versus 96% for
20260630. The wells contain real, clean, single collection organisms — they are simply not the
ones the labels name.

Read `qc_summary/QC_01_NAIVE_MAPPING.md` and `QC_02_RECONSTRUCTION.md` before using anything
here, and `../../../FINDINGS.md` §1–3 for the current understanding.

`04_minimap2_read_besthit_corroborated_db.csv` is still useful — it feeds
`../contamination_scan/`, which is what established that the wells are clean rather than mixed.
