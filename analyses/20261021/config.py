"""Paths for the 20261021 analysis: a base-exact 16S reference DB built from experiment
20260907, and the strain-by-strain distinguishability table that the higher-order
community experiment will be designed from.

20260907 is the source because it is the only run where every strain was put into two
monocultures on DIFFERENT plates and into ~8 pair wells with changing partners, so a
strain's sequence can be confirmed by independent wells rather than voted on.
"""
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE.parents[1] / "scripts"))

import paths  # noqa: E402

EXPERIMENT = "20260907"
LAYOUT_CSV = paths.layout_csv(EXPERIMENT)
DEMUX_DIR = Path("/home/rl/data/interim/karl/demultiplexing/20260907_demux_corrected/unflipped")

PRIOR = paths.ANALYSES / "20260921"
PRIOR_CONSENSUS_FASTA = PRIOR / "03_recreate_reference_db" / "outputs" / "c02_strain_consensus_corrected.fasta"
PRIOR_CONSENSUS_SUMMARY = PRIOR / "03_recreate_reference_db" / "outputs" / "c01_strain_consensus_summary_corrected.csv"
NO_GROWTH_CSV = PRIOR / "02_demux_qc" / "outputs" / "d13_strains_no_growth.csv"

FULL_COLLECTION_FASTA = paths.REFERENCE_DBS_16S / "full_collection_db.fasta"
CONTAM_INFORMED_FASTA = (paths.CONSENSUS_MERGE / "contamination_informed_db" / "results"
                         / "contamination_informed_db.fasta")

REFDB_OUT = BASE / "01_reference_db" / "outputs"
DIST_OUT = BASE / "02_distinguishability" / "outputs"

# Any bacterial 16S orients any other: >=75% identity forward vs ~50% reverse-complemented.
ORIENT_ANCHOR_STRAIN = "A1"
