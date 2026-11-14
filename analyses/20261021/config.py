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

# Repo paths are anchored on THIS folder, not on paths.py: after the 2026-11 move under
# sequencing_experiments/, paths.PW_ROOT / CONSENSUS_MERGE still point at the old locations.
PW = BASE.parents[1]                       # .../sequencing_experiments/pairwise_interaction_experiments
SEQ = PW.parent                            # .../sequencing_experiments

EXPERIMENT = "20260907"
LAYOUT_CSV = PW / "experiments" / EXPERIMENT / "01_setup" / f"strain_layout_{EXPERIMENT}.csv"
DEMUX_DIR = Path("/home/rl/data/interim/karl/demultiplexing/20260907_demux_corrected/unflipped")

PRIOR = PW / "analyses" / "20260921"
PRIOR_CONSENSUS_FASTA = PRIOR / "03_recreate_reference_db" / "outputs" / "c02_strain_consensus_corrected.fasta"
PRIOR_CONSENSUS_SUMMARY = PRIOR / "03_recreate_reference_db" / "outputs" / "c01_strain_consensus_summary_corrected.csv"
NO_GROWTH_CSV = PRIOR / "02_demux_qc" / "outputs" / "d13_strains_no_growth.csv"

FULL_COLLECTION_FASTA = paths.REFERENCE_DBS_16S / "full_collection_db.fasta"
CONSENSUS_MERGE = SEQ / "initial_database_creation" / "merge_consensus_sequences"
CONTAM_INFORMED_FASTA = (CONSENSUS_MERGE / "contamination_informed_db" / "results"
                         / "contamination_informed_db.fasta")

REFDB_OUT = BASE / "01_reference_db" / "outputs"
DIST_OUT = BASE / "02_distinguishability" / "outputs"

# Any bacterial 16S orients any other: >=75% identity forward vs ~50% reverse-complemented.
ORIENT_ANCHOR_STRAIN = "A1"

# --- genomes (03_genome_concordance) ------------------------------------------------------
# The genomes behind final_genomic_tables: panX's input GenBank files, one per collection
# strain. Wells join through mapping_384_well_plate_collection.csv (Well_souce_plate ->
# assembly_name; Schulenberg GenBank names carry an extra `_ASM..._genomic` suffix).
import os  # noqa: E402
GENOME_GBK_DIR = Path("/home/rl/Documents/Workstation_Data/Or/KEIO_collection_panX_mapping/"
                      "STRAIN_LIBRARY_SEQ_PROCCSSED_DATA/all_gbk_files_including_CeMbio43/input_GenBank")
GENOME_MAPPING = paths.RESOURCES / "mapping_384_well_plate_collection.csv"
# FASTA extracted from the GenBank files (1.4 GB, ~4 s to rebuild) -- a cache, not an output
GENOME_FNA_CACHE = Path(os.environ.get("TMPDIR", "/tmp")) / "genome_fna_cache_20261021"
SKANI = "/home/rl/mambaforge/envs/nasal-strain/bin/skani"
BLAST_BIN = Path("/home/rl/mambaforge/envs/blast/bin")
GENOME_OUT = BASE / "03_genome_concordance" / "outputs"
BARRNAP_BIN = Path("/home/rl/mambaforge/envs/NanoSeq/bin")   # barrnap 0.9 + nhmmer
