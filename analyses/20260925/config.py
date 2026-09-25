"""Cross-experiment comparison: 20260630 (merged demux) vs 20260907.

Reads the per-well scores each experiment's own analysis already produced; nothing upstream
is recomputed here.

    20260630  analyses/20260710_increased_tolerance, arm `merged` (pass_a + pass_b, corrected
              barcodes; median ~486 reads/well). Every pair well is a 1:1, 100+100 nL shot.
    20260907  analyses/20260921 (corrected demux; median ~850 reads/well). Only the wells
              shot like 20260630's count as comparable: 1:1 at 100+100 nL, i.e. well types
              pair / techrep plus the 1:1 levels of the ratio and density arms. The other
              ratio levels are excluded because ~31% of a starting-ratio change survives to
              the endpoint (20260921 README), so they would be different conditions.
"""
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE.parents[1] / "scripts"))

import paths  # noqa: E402

A = "20260630"
B = "20260907"

RA_OUT = {
    A: paths.ANALYSES / "20260710_increased_tolerance" / "05_engineer_relative_abundances"
       / "relative_abundance" / "outputs_merged_ab",
    B: paths.ANALYSES / "20260921" / "analysis" / "relative_abundance" / "outputs",
}
CONSENSUS = {
    A: paths.ANALYSES / "20260710_increased_tolerance" / "03_create_reference_db" / "consensus2"
       / "strain_consensus_20260630.fasta",
    B: paths.ANALYSES / "20260921" / "03_recreate_reference_db" / "outputs"
       / "c02_strain_consensus_corrected.fasta",
}
LAYOUT = {A: paths.layout_csv(A), B: paths.layout_csv(B)}
# strains whose monocultures never grew in 20260907 (20260921 README: exclude before modelling)
NO_GROWTH_B = paths.ANALYSES / "20260921" / "02_demux_qc" / "outputs" / "d13_strains_no_growth.csv"

MIN_RESOLVABLE_BP = 10      # below this the two strains' reads cannot be told apart
MIN_WELL_READS = 20         # on-target reads for a well to count
