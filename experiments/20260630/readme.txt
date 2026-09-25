20260630 -- pairwise interaction experiment on ameenas strains
(see the pick list in 01_setup/experiment_setup.ipynb)

Start here
----------
  analysis_index.md   which analyses used this experiment, and what they concluded
  known_issues.md     what is wrong with this data and what it affects -- READ BEFORE USE
  demux_runs.md       every demultiplexing run, and which one is current
  01_setup/                          the plate design and barcode assignment
  02_demultiplex_sequencing_data/    the demultiplexing config and notes

This is the reference experiment: as of 2026-09-22 it is one of only two whose well labels
join correctly to the genomic tables (the other is 20260907). It is used as the positive
control whenever a new experiment's strain identity is checked.

Where the data lives
--------------------
OD (destination plates)   data_links/raw/plate_reader_csvs/data_ascii/Karl_2026/Karl_20260704_OD_Full
OD (preculture / source)  data_links/raw/plate_reader_csvs/data_ascii/Karl_2026/Karl_20260623_OD
demultiplexed reads       data_links/interim/demultiplexing/20260710_demultiplex_increased_tolerance_merged/unflipped
                          (see demux_runs.md -- corrected barcode file, pass_a + pass_b, 11.9x the original reads)

Results
-------
See analysis_index.md. Best held-out-strain genomic ML result from analyses/20260710:
R2 0.321, over 74 genomes and 1465 modelling pairs.
