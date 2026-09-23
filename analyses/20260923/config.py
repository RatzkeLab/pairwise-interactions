"""Second-pass genomic ML on experiment 20260907.

Inputs come from `analyses/20260921`, which did the QC, built the per-strain consensus
and scored every well. Nothing here recomputes those: this analysis reads
`r03_pair_replicate_stats.csv` and the consensus FASTA and leaves the first pass
untouched, so its published numbers stay reproducible.

What this pass adds, in order:
  01_model_sweep   the tuned models the first pass deliberately skipped
  02_sensitivity   the caveats -- inoculum density, non-growing strains, plate-reader
                   features -- each as a refit compared against the same baseline
"""
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE.parents[1] / "scripts"))

import paths                                        # noqa: E402
from experiment_config import ExperimentConfig      # noqa: E402

EXPERIMENT = "20260907"
PREV = paths.ANALYSES / "20260921"                  # the first pass, read-only

SETUP = paths.setup_dir(EXPERIMENT)
LAYOUT_CSV = paths.layout_csv(EXPERIMENT)
CONSENSUS_FASTA = PREV / "03_recreate_reference_db/outputs/c02_strain_consensus_corrected.fasta"
STRAINS_NO_GROWTH = PREV / "02_demux_qc/outputs/d13_strains_no_growth.csv"

DEMUX_DIR = paths.DATA_LINKS / "interim/demultiplexing/20260907_demux_corrected/unflipped"
OD_FULL = paths.PLATE_READER / "Karl_20260910_ODFull"     # destination plates
OD_PREP = paths.PLATE_READER / "Karl_20260908_OD"         # preculture / source plate


def make_config():
    """ExperimentConfig whose relative-abundance outputs point at the FIRST pass.

    `exp_base` is this analysis, so anything written lands here, but the labels are read
    from 20260921 -- recomputing them would take ~15 min and produce identical files.
    """
    cfg = ExperimentConfig(
        name=EXPERIMENT,
        exp_base=BASE,
        layout_csv=LAYOUT_CSV,
        demux_dir=DEMUX_DIR,
        reference_dbs={
            "corroborated": paths.REFERENCE_DBS_16S / "corroborated_db_filtered_min5.fasta",
            "full_collection": paths.REFERENCE_DBS_16S / "full_collection_db.fasta",
            "rDNA_all_strains": paths.REFERENCE_DBS_16S / "rDNA_16S_db_all_strains.fasta",
        },
        primary_db="full_collection",
        external_cross_check_dbs=["corroborated", "rDNA_all_strains"],
        ra_reference_fasta=CONSENSUS_FASTA,
        min_reads=10,
    )
    cfg.relative_abundance_out_dir = PREV / "analysis/relative_abundance/outputs"
    cfg.mapping_validation_out_dir = PREV / "analysis/mapping_validation/outputs"
    return cfg
