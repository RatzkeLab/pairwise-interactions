# 20260925: does 20260630 reproduce in 20260907?

First pass at biological replication across experiments: how well does one experiment
predict the other? Inputs are each experiment's existing per-well scores
(20260630: `analyses/20260710_increased_tolerance`, merged demux; 20260907:
`analyses/20260921`). See `config.py` for exactly which wells count.

    python 01_cross_experiment/cross_experiment_replicates.py     (env karl_seq_analysis, ~1 min)

## Answer

| predicting | n | R² (no fit, vs identity) | r | ρ | winner agrees |
|---|---|---|---|---|---|
| 20260907 pairs from 20260630 pairs | 53 pairs | **0.88** [0.80, 0.93] | 0.96 | 0.87 | 81% |
| 20260630 pairs from 20260907 pairs | 53 pairs | **0.91** [0.84, 0.95] | 0.96 | 0.87 | 81% |
| *yardstick:* one well from another, within 20260630 | 275 well pairs | 0.89 | 0.94 | 0.96 | 95% |
| *yardstick:* one well from another, within 20260907 | 234 well pairs | 0.83 | 0.92 | 0.91 | 87% |
| *same footing:* one well from one well, across experiments | 62 well pairs | 0.86 | 0.96 | 0.90 | 82% |

Target log2(a/b); relative abundance gives the same picture (R² 0.85, r 0.93). Brackets are
95% bootstrap CIs over pairs.

**Across-experiment agreement is about as good as replicate wells within one experiment.**
Where the same pair was tested in July and September, the September value is predicted with
R² ≈ 0.88 by the July value with no fitting at all: a slope close to 1 and no offset.

Using everything rather than the overlap, **per-strain competitiveness** transfers too:
Bradley-Terry strengths of the 55 strains fitted in both correlate r = 0.89 / ρ = 0.75.
20260907's strengths, fitted on only 70 pairs among the shared strains, predict **812 20260630
pairs that 20260907 never tested** at R² 0.66, r 0.84, ρ 0.81. That is the honest
"predict an unseen pair from another experiment" number, and it is limited by how sparse
20260907's design is among these strains.

## Caveats

- **Only 53 directly comparable pairs.** 107 pairs were tested in both, then pairs without a
  reference, below the 10 bp resolution limit, with a 20260907 non-grower, or with a
  mismatched strain identity (next point) were removed. The CIs are wide accordingly.
- **14 of the 82 strain names with a consensus in both experiments are not the same
  organism** (`x01_identity.csv`); they were excluded. 7 of them are 20260630 calls from the
  big shared 16S type T01 with coverage ~0.52–0.55, i.e. near coin-flip consensus calls in an
  experiment with no clean mono wells. Their 20260630 reference is probably the pair
  partner's sequence. The other ~6 (A12, E4, J9, E13, O21, H5) are confident in both
  experiments and still disagree. **Needs a follow-up**, since these references also feed
  20260630's own relative-abundance scores.
- **Pearson is flattered by lopsided pairs.** The log2 distribution has a mass near 0 and a
  few strong winners; ρ (0.87) is the more conservative summary.
- **Different depths** (~470 vs ~825 reads/well) stretch the log2 tails differently (see
  `../20260710_increased_tolerance/README.md`). Relative abundance, which is bounded, agrees
  to the same degree, so that is not what drives the result.
- 20260907's ratio-titration and density arms are included only at their 1:1, 100+100 nL
  levels.

## Outputs (`01_cross_experiment/outputs/`)

`x01_identity.csv` · `x02_direct_pairs*.csv` · `x03_well_level_yardstick.csv` ·
`x04_strength_transfer.csv`, `x04_strength_correlation.csv`, `x04_predictions_*.csv` ·
`figures/x_cross_experiment.png`
