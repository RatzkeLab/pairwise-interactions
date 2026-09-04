the setup folder primarily contains experiment_setup.ipynb and its products

two of the most important products are 
 - minibar_primers.tsv, which is used in demultiplexing
  - strain_layout_20260630.csv, which is needed to map data back to the experimental setup, later


the notebook references PrimerPlateSpecs.csv, an artifact of the primer ordering process containing sequences for the primers used in this experiment
the file has now been moved to /home/rl/scripts/karl/Link to Karl/PrimerPlateSpecs.csv and the reference will need to be changed if you want to re-run the notebook

notes on issues and errors with the experiment setup:
1) strains and interactions were not selected carefully enough for all to be distinguishable in coculture with 16s sequences. these are later identified and removed from relative abundance analysis
2) strain_check, which is based on OD readings of the strains growing in minimal media, incorrectly identifies several strains as "no growers" and attempts to use them as blank cells for the monoculture wells. this results in there being no purely reliable monoculture samples for this experiment, which makes finding a decisive consnsensus later more difficult, though we managed to corroboarte between pairs in the end. nonetheless, clean/pure mono wells are recommended for furture experiments
3) although 30 experimental plates were setup and read in the plate reader only the odd numbered plates from 1 to 19 were sequenced since 10 plates is about the maximum that can be sequenced in one run with our current sequencing protocol