20261026 -- community experiment, round 2
(20261021's plates 1-5 well for well, plus a second daily-dilution plate 6 that shares
31 of its 102 communities with plate 5; see 01_setup/README.md)

Start here
----------
  01_setup/README.md          what differs from 20261021, the barcode-plate numbering
  01_setup/outputs/           Echo picklists (strains: plates 1-6; primers: plates 1-12),
                              layout.csv, communities.csv, minibar sheet, demux config
  01_setup/inspect_setup.ipynb  plots of the design

Plates
------
  1-4   main, 24 h endpoint, sequenced (same as 20261021)
  5     passage plate (daily dilution), same as 20261021
  6     second passage plate, diluted in parallel with 5
  barcode plates: 5 / 6 = plates 5 / 6 at 24 h; 7 / 8 at 48 h; 9 / 10 at 5 d; 11 / 12 at 10 d
        (NOT the 20261021 numbering, where 6-8 were plate 5's later time points)
  99    decoy barcodes on the sample sheet only (no physical plate)

Barcodes are new (none shared with 20261021), so use these primer picklists, not 20261021's.
Record the day's source-plate OD: inoculation is equal volumes, as in 20261021.
