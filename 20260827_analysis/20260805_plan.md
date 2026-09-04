now let's go for a relative abundance analysis. 

1) using 20260721/analysis/consensus/merged_consensus_20260721.fasta as a reference, and 20260721/setup/strain_layout_20260721_plate1_2_swapped.csv as our guide, assume that only the two strains in question are present (within reasonable limits) and give the sample an interaction score. if there are sequences in a sample that really dont align to either of the two strains, then this should be quantified and noted. strains that cant be told apart can be put down as 50 50 but this should also be noted. (perhaps assigned a quality score or an uncertainty score or similar) but even if they are only a few bp (<10) distant from eachother, we should still at least try to see if the mapping algorithm can distinguish two groups and assign relative abundances to each. 

2) then once all samples (all samples with >5 reads, that is) have been analyzed in this way, do a quality check on all samples that have replicates (i.e. it's the same two strains expected to be in that sample) -- is the interaction stable? or does it change from replicate to replicate

3) regardless of replicate reliability we should give every strain a "competitiveness score" by averaging across all pairwise interactions for that strain (or similar), then saving this to a table or similar this should help us answer scientific questions later. (i.e. do some strains just seem to have a hard time growing? / always lose? do some strains tend to co-exist? do some always dominate?)

4) lastly, i would like to do a heirarchy analysis based on this data. is this strain library highly heirarchichal, or is it a more random/tangled network of competitive interactions

once again, make this code as reusable as possible, and do the bulk of the analysis in python functions that can be imported into an ipynb notebook where possible, (ideally also create and run such a notebook, once again in the karl_seq_analysis conda environment)