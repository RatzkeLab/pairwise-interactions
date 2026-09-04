# mapping_validation — 20260630

**Question:** does each well actually contain the two strains the layout says it should?

Run: `mapping_validation_20260630.ipynb` (re-executes top to bottom; it is the report).
Logic lives in `shared_pipelines/mapping_validation.py`.

Reads are mapped against this experiment's own consensus (`consensus2/`) as the primary
reference, with `corroborated_db` as an external cross-check. Two references are used because a
self-built consensus is self-consistent by construction and so cannot detect a systematic
labelling error on its own.

| output | contents |
|---|---|
| `01_samples_gt5reads.csv` | wells passing the read-count floor, with layout metadata |
| `02_reference_coverage.csv`, `02_reference_cross_check.csv` | which strains have a usable reference; identity check between the two databases |
| `03_edlib_read_assignments_*.csv.gz` | per read, distance to each expected strain |
| `04_minimap2_read_besthit_*.csv` | per read, best hit across the **whole** database — the basis of the contamination scan |
| `05_combined_sample_summary.csv`, `05_contamination_candidates.csv` | per-well verdict |

**Verdict for this experiment: labels are supported** (96% of testable labels). See
`qc_summary/QC_01_NAIVE_MAPPING.md`.

**Caveat:** `03_*` stores *normalized* edit distances; `norm = ed / max(read_len, ref_len)`, so
raw bp is recoverable given `read_len`. See TRAPS §6 and §11 for the two identifier gotchas.
