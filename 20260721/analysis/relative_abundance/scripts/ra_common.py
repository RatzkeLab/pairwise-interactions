"""Shared paths, config, and helper functions for the 20260721 relative-abundance /
competitive-interaction analysis.

Reuses the already-validated read/reference IO helpers from the sibling
analysis/mapping_validation/scripts/common.py (loaded by file path, not via sys.path, so
there is no risk of a bare `import common` resolving to the wrong module if both analyses
are ever imported in the same process).
"""

import importlib.util
from pathlib import Path

import edlib

# ---------------------------------------------------------------------------
# pull in the proven IO helpers from analysis/mapping_validation/scripts/common.py
# ---------------------------------------------------------------------------
_MV_COMMON_PATH = (
    Path(__file__).resolve().parents[2] / "mapping_validation" / "scripts" / "common.py"
)
_spec = importlib.util.spec_from_file_location("mapping_validation_common", _MV_COMMON_PATH)
_mv_common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mv_common)

EXP_BASE = _mv_common.EXP_BASE
LAYOUT_CSV = _mv_common.LAYOUT_CSV
DEMUX_DIR = _mv_common.DEMUX_DIR

load_layout = _mv_common.load_layout
sample_fastq_path = _mv_common.sample_fastq_path
count_reads_fast = _mv_common.count_reads_fast
load_reads = _mv_common.load_reads
load_reference_db = _mv_common.load_reference_db
norm_edit_distance = _mv_common.norm_edit_distance

# ---------------------------------------------------------------------------
# this analysis's own paths / config
# ---------------------------------------------------------------------------
RA_BASE = EXP_BASE / "analysis" / "relative_abundance"
OUT_DIR = RA_BASE / "outputs"
FIG_DIR = OUT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# per user instruction: this analysis's reference is merged_consensus_20260721.fasta
# specifically (not the mono-priority variant used as PRIMARY_DB in mapping_validation)
REFERENCE_FASTA = EXP_BASE / "analysis" / "consensus" / "merged_consensus_20260721.fasta"

# reuse the exact ">5 reads" cohort already computed by mapping_validation (same definition:
# "all samples with more than 5 reads") instead of recomputing it
SAMPLES_GT5READS_CSV = (
    EXP_BASE / "analysis" / "mapping_validation" / "outputs" / "01_samples_gt5reads.csv"
)

# a read's best (minimum) distance to either candidate must be below this to count as
# "on target" at all -- otherwise it doesn't look like either expected strain. Reused
# verbatim from mapping_validation's data-driven confident-match threshold (same
# instrument, same 16S amplicon, same ONT error profile -- see that project's s05).
OFF_TARGET_DIST_THRESHOLD = 0.128

# a read is "ambiguous" (can't tell which of strain1/strain2 it is) if the *margin*
# between its distances to the two candidates is smaller than this. Below this margin the
# two distances are within each other's ONT-error noise band, so calling a winner would be
# reading noise as signal.
AMBIGUOUS_MARGIN_THRESHOLD = 0.02


def pair_key(strain_a, strain_b):
    """Order-independent identifier for an unordered strain pair."""
    return tuple(sorted((strain_a, strain_b)))


def edit_distance_bp(a, b):
    return edlib.align(a, b, mode="NW", task="distance")["editDistance"]
