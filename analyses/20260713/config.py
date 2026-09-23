"""Experiment-specific configuration for 20260721's mapping_validation and
relative_abundance pipelines (both share this one CFG -- see shared_pipelines/).
"""

import sys
from pathlib import Path

_SHARED = Path(__file__).resolve().parents[2] / "shared_pipelines"
sys.path.insert(0, str(_SHARED))

from experiment_config import ExperimentConfig  # noqa: E402

EXP_BASE = Path(__file__).resolve().parents[1]

CFG = ExperimentConfig(
    name="20260721",
    exp_base=EXP_BASE,
    layout_csv=EXP_BASE / "setup" / "strain_layout_20260721_plate1_2_swapped.csv",
    demux_dir=EXP_BASE / "data" / "demux" / "unflipped",  # final, complete demux (9120/9120 wells present)
    reference_dbs={
        # independent, external reference built from prior corroborated experiments
        "corroborated_db": Path(
            "/home/rl/scripts/karl/merge_consensus_sequences/collapse_naive_updated3_15diff/"
            "corroborated_db_filtered_min5.fasta"
        ),
        # this-experiment consensus built by analysis/per_strain_consensus.ipynb
        "merged_consensus": EXP_BASE / "analysis" / "consensus" / "merged_consensus_20260721.fasta",
        "merged_consensus_mono_priority": EXP_BASE / "analysis" / "consensus" / "merged_consensus_mono_priority_20260721.fasta",
    },
    # self-consistent by construction (built from this experiment's own reads) -- see
    # shared_pipelines/mapping_validation.py's identity cross-check for why corroborated_db
    # is not used as the primary reference
    primary_db="merged_consensus_mono_priority",
    external_cross_check_dbs=["corroborated_db"],
    # relative_abundance uses merged_consensus (not the mono-priority variant) per the
    # original analysis request for this experiment
    ra_reference_fasta=EXP_BASE / "analysis" / "consensus" / "merged_consensus_20260721.fasta",
    min_reads=5,
)
