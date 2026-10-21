20261021 -- community experiment, round 1
(communities of 3-66 members in triplicate; see 01_setup/README.md for the design and why)

Not to be confused with analyses/20261021, which built the 16S reference DB, the
distinguishability table and the genome checks that this design draws on.

Start here
----------
  01_setup/README.md          the design, the inoculation, the barcodes
  01_setup/outputs/           Echo picklists (strains: plates 1-5; primers: plates 1-8),
                              layout.csv, communities.csv, minibar sheet, demux config
  01_setup/inspect_setup.ipynb  plots of the design

Plates
------
  1-4   main, 24 h endpoint, sequenced
  5     passage plate (daily dilution); its 24 h / 48 h / 5 d / 10 d snapshots are
        barcode plates 5 / 6 / 7 / 8 (same layout)
  99    decoy barcodes on the sample sheet only (no physical plate)

Record today's source-plate OD: inoculation used equal volumes, so preculture density is
the covariate the analysis will need.
