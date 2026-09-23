20260721 -- pairwise incubation experiment on karls 87 selected strains
(minbp 5, from the previously created corroborated_db)

Start here
----------
  known_issues.md     READ THIS FIRST -- the well labels do not identify the organisms
  analysis_index.md   which analyses used this experiment, and what they concluded
  demux_runs.md       every demultiplexing run, and which one is current
  01_setup/           the plate design and barcode assignment
  02_demultiplexing/  the demultiplexing config and notes

The sequencing is fine; the label-to-organism mapping is not. This experiment scores 0 of
188 decidable comparisons against every reference database, and fails the genomic join
gate at z=-0.7 where 20260630 passes at z=+6.0. The source plate is the same physical
stock as 20260630, so the mix-up happened downstream of the freezer.

About half the labels are recoverable by 16S -- see analyses/20260827.

Where the data lives
--------------------
OD (destination plates)   data_links/raw/plate_reader_csvs/data_ascii/Karl_2026/Karl_20260723_OD_Full
OD (preculture / source)  data_links/raw/plate_reader_csvs/data_ascii/Karl_2026/Karl_20260722_OD
demultiplexed reads       data_links/interim/demultiplexing/20260730_demux/unflipped
                          (see demux_runs.md -- this run under-assigns ~4x, fix available)

Results
-------
See analysis_index.md. Usable for method development and as a negative control; not usable
for genomic ML until the labels are recovered.
