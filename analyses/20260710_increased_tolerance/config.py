"""ExperimentConfig for 20260630, re-analysed on the corrected demultiplexing.

This folder started as a copy of analyses/20260710 (2026-09-24). What changed:

  - reads come from 20260710_demultiplex_increased_tolerance, the re-demux with the
    corrected barcode file and -e 4 (see experiments/20260630/demux_runs.md). The
    original demux is kept as a second ARM so both go through the SAME current code:
    the shared library has changed since 20260710 first ran (margin fix, plate_role,
    ...), so comparing against 20260710's stored outputs would mix code drift into
    the demux comparison.

        DEMUX_ARM=corrected   (default) -> outputs/
        DEMUX_ARM=original              -> outputs_original_demux/

  - only odd plates 1-19 were sequenced. The old analysis enforced that by reading a
    hand-copied folder (relevant_fastqs) holding only those plates. Here it is a
    `plate_role` column on the layout, which mapping_validation.gather_samples already
    honours, so the demux output is read directly. Even-plate wells are the negative
    control for mis-assignment.

  - the reference is UNCHANGED: 03_create_reference_db/consensus2 (84 strains), so the
    two arms differ only in which reads each well received. The 16 strains without a
    consensus failed on `no_majority` (wells disagree), not depth, so a rebuild is not
    expected to rescue many of them.

  - paths go through scripts/paths.py (post-2026-09-23 layout, data_links/).
"""
import os
import sys
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent            # this analysis
sys.path.insert(0, str(BASE.parents[1] / "scripts"))

import paths                                       # noqa: E402
from experiment_config import ExperimentConfig     # noqa: E402

EXPERIMENT = "20260630"
SEQUENCED_PLATES = list(range(1, 20, 2))           # odd plates 1-19; 10 plates per run max

DEMUX_ROOT = paths.DATA_LINKS / "interim" / "demultiplexing"
DEMUX_DIRS = {
    "corrected": DEMUX_ROOT / "20260710_demultiplex_increased_tolerance" / "unflipped",
    "original": DEMUX_ROOT / "20260710_demultiplex" / "unflipped",
}
OUT_SUFFIX = {"corrected": "", "original": "_original_demux"}

CONSENSUS_FASTA = BASE / "03_create_reference_db" / "consensus2" / "strain_consensus_20260630.fasta"
LAYOUT_CSV_FULL = paths.layout_csv(EXPERIMENT)     # all 30 plates -- what the OD scripts need
LAYOUT_CSV = BASE / "inputs" / f"strain_layout_{EXPERIMENT}_plate_role.csv"

OD_FULL = paths.OD_FULL[EXPERIMENT]
OD_PREP = paths.OD_PREP[EXPERIMENT]

ARM = os.environ.get("DEMUX_ARM", "corrected")
if ARM not in DEMUX_DIRS:
    raise ValueError(f"DEMUX_ARM must be one of {list(DEMUX_DIRS)}, got {ARM!r}")


def _write_layout_with_plate_role():
    lay = pd.read_csv(LAYOUT_CSV_FULL)
    lay["plate_role"] = lay["dest_plate"].map(
        lambda p: "primary_sequenced" if p in SEQUENCED_PLATES else "not_sequenced")
    LAYOUT_CSV.parent.mkdir(exist_ok=True)
    lay.to_csv(LAYOUT_CSV, index=False)


if not LAYOUT_CSV.exists():
    _write_layout_with_plate_role()


def out_dir(step_dir, arm=ARM):
    """`<step_dir>/outputs` for the corrected arm, `outputs_original_demux` for the original."""
    return BASE / step_dir / f"outputs{OUT_SUFFIX[arm]}"


def make_config(arm=ARM):
    return ExperimentConfig(
        name=EXPERIMENT,
        exp_base=BASE,
        layout_csv=LAYOUT_CSV,
        demux_dir=DEMUX_DIRS[arm],
        reference_dbs={
            "corroborated_db": paths.CORROBORATED_DB,
            "consensus2": CONSENSUS_FASTA,
        },
        primary_db="consensus2",
        external_cross_check_dbs=["corroborated_db"],
        ra_reference_fasta=CONSENSUS_FASTA,
        min_reads=5,
        mapping_validation_out_dir=out_dir("04_qc/mapping_validation", arm),
        relative_abundance_out_dir=out_dir("05_engineer_relative_abundances/relative_abundance", arm),
    )


CFG = make_config()
