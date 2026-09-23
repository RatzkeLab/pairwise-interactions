20260907 -- whole-collection pairwise interaction experiment
(all 384 strains; see 01_setup/README.md for the design and why each choice was made)

Start here
----------
  analysis_index.md   which analyses used this experiment, and what they concluded
  known_issues.md     what is wrong with this data and what it affects -- READ BEFORE USE
  demux_runs.md       every demultiplexing run, and which one is current
  01_setup/           the plate design, barcode assignment, and both demux configs

What makes this run different
-----------------------------
8 plates PCR'd and sequenced, 8 more plate-reader-only and frozen. Every one of the 384
strains has a monoculture on the sequenced plates, 376 have two on DIFFERENT plates -- so
a strain's reference sequence is anchored on its own monocultures rather than on a
majority vote across pair wells, and a contaminated well is visible instead of averaged
over. Pairs were filtered to >=15 bp of 16S separation, so no wells are spent on strains
that cannot be told apart.

340 wells go to two inoculum arms (ratio and density titration) that measure how much of
the outcome is set by the starting ratio rather than by fitness.

The frozen plates are barcoded but never amplified, which makes them a permanent
false-assignment monitor for any demultiplexing run against this data.

Where the data lives
--------------------
OD (destination plates)   data_links/raw/plate_reader_csvs/data_ascii/Karl_2026/Karl_20260910_ODFull
                          plate number is in the ID1: header field
OD (preculture / source)  data_links/raw/plate_reader_csvs/data_ascii/Karl_2026/Karl_20260908_OD
                          mixed in with an unrelated lactate assay -- match on Testname
raw reads                 /var/lib/minknow/data/kr_pairwise_03/no_sample_id/
                              20260917_1657_MN37816_FBH02602_5ecfeeb4/fastq_pass
demultiplexed reads       data_links/interim/demultiplexing/20260907_demux_corrected/unflipped
                          USE THIS ONE -- see demux_runs.md

Results
-------
See analysis_index.md. Labels are NOT scrambled (unlike 20260721); this experiment passes
the genomic join gate at rho=+0.387, z=+8.2, better than 20260630.
