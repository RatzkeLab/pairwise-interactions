# experiments/

One folder per physical experiment: the plate design, the barcode sheet, and pointers to
the raw data. Analyses live in `../analyses/`; see `../scripts/paths.py` for the layout.

## Setting up the next experiment: do it differently

Every setup so far (20260630, 20260721, 20260907) was a Jupyter notebook that is both the
code and the record of what was done. That made one bug cost three experiments:
`experiment_setup.ipynb` wrote the full ordered oligo (15 bp adapter + 24 bp barcode + 20 bp
primer) into minibar's index column. Nothing checked the output, re-running the notebook
risks regenerating the layout, and the fix could only be applied after the fact
(`analyses/20260921/02_demux_qc/fix_barcode_file.py`). The 2025 barcode file
(`merge_consensus_sequences/250916_resequence/barcodes250807.tsv`) was correct, so this was
introduced when setup moved into the notebook. See `20260630/known_issues.md`.

**Proposed structure:**

```
scripts/experiment_setup/        shared, importable, tested once
    layout.py                    strain pick list -> plate layout (pure functions, seeded)
    barcodes.py                  layout + PrimerPlateSpecs -> minibar TSV
    echo.py                      layout -> Echo transfer CSVs
    validate.py                  checks run on every output (below)
experiments/<date>/01_setup/
    setup_config.yaml            the ONLY per-experiment input: strains, plates, seed,
                                 volumes, which plates get sequenced, paths
    run_setup.py                 reads the yaml, calls the library, writes outputs,
                                 runs validate.py, writes a manifest (inputs, git hash,
                                 checksums); refuses to overwrite existing outputs
    outputs/                     written once, then treated as read-only
    inspect_setup.ipynb          optional: plots and eyeballing, reads outputs/, writes nothing
```

The rule that matters: **code that produces an artifact never lives in the same file as
the record of the artifact.** Notebooks are fine for looking at things, not for making them.

A lighter option if a package feels heavy: keep one notebook, but move every cell that
writes a file into functions in `scripts/`, and have the notebook call them and then
`validate()`. That is most of the benefit for a fraction of the work.

**Checks `validate.py` should run** (each would have caught a real problem here):

- index columns are the barcode alone: all the same length (24), none ending in the
  primer, no shared prefix across all rows (a shared prefix means an adapter got in)
- minimum pairwise edit distance of the barcode set, printed next to the planned `-e`;
  fail if `2 * e >= dmin`
- barcode combinations used / possible (sheet density; see the note on `-e` below)
- a `plate_role` column in the layout (`primary_sequenced` / `not_sequenced`) rather than a
  hand-copied folder of the sequenced plates (20260630's `relevant_fastqs/`)
- strains in any pair at least 10 bp apart in 16S (`MIN_RESOLVABLE_BP`)
- true monoculture wells (strain + blank medium, blank confirmed as a non-grower)
- a demultiplexing config written from the same yaml, using the key names the pipeline
  actually reads (`barcode_edit_distance`, not `barcode_edit_dist`)

**And one check after sequencing:** minibar tags every read with how it matched (`HH` =
barcode + primer found at both ends; `hh` = barcode-only fallback). A healthy run is ~99%
`HH`. 20260630's original demux was 99.99% `hh`. That one-line count would have flagged
the bug on day one.

## Choosing minibar's `-e`

Two things set the safe barcode edit distance: how far apart the barcodes are, and how
many barcode combinations the sheet uses.

- This barcode set's minimum pairwise edit distance is 12. At `-e 6` a read can sit
  exactly between two barcodes, so misreads start being assigned.
- A misread only turns into a wrong sample if the resulting fwd/rev combination is in the
  sheet. In 2025 the sheet used 384 of 6,972 combinations (5.5%), so most misreads became
  "no such sample" and `-e 6` was safe. The 2026 layouts use 9,120 of 9,216 (99%), so
  nearly every misread lands in a real well.

Measured (20260907, against never-amplified wells): `-e 4` 1.2% false, `-e 5` 1.5%,
`-e 6` 5.2%. **Use `-e 4` for dense combinatorial sheets**; `-e 6` is fine for sparse ones.
Keep some unsequenced barcode combinations on the sheet as a standing negative control.
