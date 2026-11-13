"""Paths for 20261113: does the genomic ML survive (1) organism-level CV and (2) dropping
wells whose genome on file is not the organism in the well?

Everything is read-only from earlier analyses. Paths are built relative to this file rather
than through scripts/paths.py, because on 2026-11-13 the repository moved under
`sequencing_experiments/` and paths.PW_ROOT / ANALYSES / EXPERIMENTS / CONSENSUS_MERGE still
point at the old location (paths.check() reports 7 missing). genomic_ml itself only needs
GENOMIC_TABLES and MAPPING_DIR, which still resolve.

Run everything with env karl_seq_analysis.
"""
import os
import sys
from pathlib import Path

# One BLAS thread per process: the CV runs one process per fold, and numpy's default of a
# thread per core turned 14 workers into ~1000 threads on a shared 64-core machine. Must be
# set before numpy is first imported, so every script imports config first.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = Path(__file__).resolve().parent                 # .../analyses/20261113/messy
PW = HERE.parents[2]                                   # .../pairwise_interaction_experiments
ANALYSES = PW / "analyses"
SCRIPTS = PW / "scripts"
KARL = PW.parents[1]                                   # /home/rl/scripts/karl
for p in (str(SCRIPTS), str(KARL)):                    # genomic_ml imports via the package path
    if p not in sys.path:
        sys.path.insert(0, p)

CLEAN = HERE.parent / "clean"

# --- the genome-concordance results from 20261021 (describe the 20260907 source plate) -------
CONC = ANALYSES / "20261021" / "03_genome_concordance"
CONC_OUT = CONC / "outputs"
REFDB_0907 = ANALYSES / "20261021" / "01_reference_db" / "outputs" / "reference_db_20260907.fasta"
REFDB_0907_HAP = ANALYSES / "20261021" / "01_reference_db" / "outputs" / "reference_db_20260907_haplotypes.fasta"

COLLECTION_SHEET = (PW.parent / "initial_database_creation" / "merge_consensus_sequences" /
                    "plate_format" / "Full_384_strain_collection_no_seq_info.csv")
STRAIN_MAPPING = KARL / "data_links" / "resources" / "mapping_384_well_plate_collection.csv"

# genomes behind the KO tables, and tools (same as 20261021/config.py)
GENOME_GBK_DIR = Path("/home/rl/Documents/Workstation_Data/Or/KEIO_collection_panX_mapping/"
                      "STRAIN_LIBRARY_SEQ_PROCCSSED_DATA/all_gbk_files_including_CeMbio43/input_GenBank")
BLAST_BIN = Path("/home/rl/mambaforge/envs/blast/bin")

# --- the two experiments whose genomic ML is being re-examined --------------------------------
# r03     the interaction labels the published ML read
# ref     the 16S references those labels were scored against (also the 16S-only control)
# pub     the published CV outputs, used to check this rerun reproduces them fold for fold
EXPERIMENTS = {
    "20260907": dict(
        r03=ANALYSES / "20260921" / "analysis" / "relative_abundance" / "outputs",
        ref=ANALYSES / "20260921" / "03_recreate_reference_db" / "outputs" / "c02_strain_consensus_corrected.fasta",
        pub=ANALYSES / "20260923" / "01_model_sweep" / "outputs",
        pub_prefix="s01",
    ),
    "20260630": dict(
        r03=(ANALYSES / "20260710_increased_tolerance" / "05_engineer_relative_abundances" /
             "relative_abundance" / "outputs_merged_ab"),
        ref=(ANALYSES / "20260710_increased_tolerance" / "03_create_reference_db" / "consensus2" /
             "strain_consensus_20260630.fasta"),
        pub=ANALYSES / "20260710_increased_tolerance" / "06_ml_analysis" / "genomic_ml" / "outputs_merged_ab",
        pub_prefix="g03",
    ),
}

OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)


def gml_config(exp, out_dir):
    """GenomicMLConfig pointing at an experiment's published labels, writing to out_dir.

    exp_cfg is a stand-in carrying only what genomic_ml reads (relative_abundance_out_dir,
    ra_reference_fasta, exp_base): the real ExperimentConfig goes through paths.layout_csv,
    which is broken by the move, and nothing here touches layouts or reads.
    """
    from types import SimpleNamespace
    import genomic_ml as gml
    e = EXPERIMENTS[exp]
    exp_cfg = SimpleNamespace(name=exp, exp_base=out_dir, relative_abundance_out_dir=e["r03"],
                              ra_reference_fasta=e["ref"])
    out_dir.mkdir(parents=True, exist_ok=True)
    return gml.GenomicMLConfig(exp_cfg=exp_cfg, out_dir=out_dir)
