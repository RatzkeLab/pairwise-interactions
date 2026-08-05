using check_reads as reference, let's do some complete analysis on the 20260721 data ... here's what i'm thinking to start:

do a constrained mapping on all samples with more than 5 reads i.e. use minimap2 or edlib or similar to map those reads back to what we expect to be there (i.e. the corroborated_db.fasta and/or the newly created merged_consensus_mono_priority_20260721.fasta and/or merged_consensus_20260721.fasta) if they match up with strain1 and strain2 in the pairwise_interaction_experiments/20260721/setup/strain_layout_20260721_plate1_2_swapped.csv )

please do this analysis in a subfolder of the /home/rl/scripts/karl/pairwise_interaction_experiments/20260721/analysis folder, creating additional sub-sub folders as needed

save results as you best see fit, but i reccomend creating numbered python scripts for everything you do to keep a clear record and maybe having a separate outputs folder. alternatively, if a jupyter notebook is the best tool for the job, then that would also work great. in general, perform the analysis so that it is documented, structured, reproducible, and human readable at the end. use the conda environment karl_seq_analysis
