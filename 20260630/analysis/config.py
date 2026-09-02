"""Experiment-specific configuration for 20260630's mapping_validation and
relative_abundance pipelines (both share this one CFG -- see shared_pipelines/).

Differs from 20260721 in three structural ways (see project memory /
shared_pipelines design discussion for why these needed a real decision, not a guess):
  - reads: the user curated and moved the wells to analyze into data/relevant_fastqs
  - reference: analysis/consensus2/strain_consensus_20260630.fasta (84 strains) -- the more
    recent/refined of two independent per-strain consensus builds, per the user. There is no
    mono-priority variant here because (per the user) 20260630 has no *true* mono-culture
    wells, unlike 20260721.
  - no strain_layout "_plate1_2_swapped" variant exists for this experiment -- the plain
    strain_layout_20260630.csv is used directly.
"""

import sys
from pathlib import Path

_SHARED = Path(__file__).resolve().parents[2] / "shared_pipelines"
sys.path.insert(0, str(_SHARED))

from experiment_config import ExperimentConfig  # noqa: E402

EXP_BASE = Path(__file__).resolve().parents[1]

CFG = ExperimentConfig(
    name="20260630",
    exp_base=EXP_BASE,
    layout_csv=EXP_BASE / "setup" / "strain_layout_20260630.csv",
    demux_dir=EXP_BASE / "data" / "relevant_fastqs",
    reference_dbs={
        "corroborated_db": Path(
            "/home/rl/scripts/karl/merge_consensus_sequences/collapse_naive_updated3_15diff/"
            "corroborated_db_filtered_min5.fasta"
        ),
        "consensus2": EXP_BASE / "analysis" / "consensus2" / "strain_consensus_20260630.fasta",
    },
    primary_db="consensus2",
    external_cross_check_dbs=["corroborated_db"],
    ra_reference_fasta=EXP_BASE / "analysis" / "consensus2" / "strain_consensus_20260630.fasta",
    min_reads=5,
)
