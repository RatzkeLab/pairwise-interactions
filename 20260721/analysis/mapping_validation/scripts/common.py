"""Shared paths, config, and helper functions for the 20260721 read-mapping validation.

Every numbered script in this folder imports from here so that paths and parameters
stay consistent across the pipeline (and so the summary notebook can reuse the exact
same functions instead of re-implementing them).
"""

from pathlib import Path
import gzip

import edlib
import pandas as pd
from Bio import SeqIO

# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------
EXP_BASE = Path("/home/rl/scripts/karl/pairwise_interaction_experiments/20260721")
LAYOUT_CSV = EXP_BASE / "setup" / "strain_layout_20260721_plate1_2_swapped.csv"
DEMUX_DIR = EXP_BASE / "demux" / "unflipped"  # final, complete demux (9120/9120 wells present)

VAL_BASE = EXP_BASE / "analysis" / "mapping_validation"
OUT_DIR = VAL_BASE / "outputs"
FIG_DIR = OUT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# reference databases to map reads against ("and/or" from the analysis plan)
REFERENCE_DBS = {
    # independent, external reference built from prior corroborated experiments;
    # complete (87/87 strains in the 20260721 layout have an entry)
    "corroborated_db": Path(
        "/home/rl/scripts/karl/merge_consensus_sequences/collapse_naive_updated3_15diff/"
        "corroborated_db_filtered_min5.fasta"
    ),
    # this-experiment consensus built by analysis/per_strain_consensus.ipynb
    # (86/87 strains -- N13 has no consensus, too few informative wells)
    "merged_consensus": OUT_DIR.parent.parent / "consensus" / "merged_consensus_20260721.fasta",
    "merged_consensus_mono_priority": OUT_DIR.parent.parent
    / "consensus"
    / "merged_consensus_mono_priority_20260721.fasta",
}

MIN_READS = 5  # "samples with more than 5 reads" per the analysis plan

# ---------------------------------------------------------------------------
# sample / fastq helpers
# ---------------------------------------------------------------------------


def sample_fastq_path(plate, well):
    """demux fastqs are named Plate{plate:02d}_{well}.fastq (or .fastq.gz)."""
    stem = f"Plate{int(plate):02d}_{well}"
    for ext in (".fastq", ".fastq.gz"):
        p = DEMUX_DIR / f"{stem}{ext}"
        if p.exists():
            return p
    return None


def count_reads_fast(path):
    """Read count via line counting (4 lines/record) -- avoids full FASTQ parsing."""
    if path is None or not path.exists() or path.stat().st_size == 0:
        return 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        return sum(1 for _ in fh) // 4


def load_reads(path):
    """Return list of (read_id, sequence) tuples from a fastq(.gz)."""
    if path is None:
        return []
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    opener = gzip.open if str(path).endswith(".gz") else open
    out = []
    with opener(path, "rt") as fh:
        for rec in SeqIO.parse(fh, "fastq"):
            out.append((rec.id, str(rec.seq).upper()))
    return out


def load_layout():
    df = pd.read_csv(LAYOUT_CSV)
    df["sample_id"] = df.apply(
        lambda r: f"Plate{int(r['dest_plate']):02d}_{r['dest_well']}", axis=1
    )
    return df


# ---------------------------------------------------------------------------
# reference database helpers
# ---------------------------------------------------------------------------


def load_reference_db(path):
    """Load a reference FASTA -> {strain_name: sequence}.

    Header conventions differ by db:
      - corroborated_db_filtered_min5.fasta   : '>STRAIN_consensus ...'
      - merged_consensus*_20260721.fasta      : '>STRAIN ...' or '>STRAIN|mono'
    """
    seqs = {}
    for rec in SeqIO.parse(path, "fasta"):
        header = rec.id
        if header.endswith("_consensus"):
            name = header[: -len("_consensus")]
        else:
            name = header.split("|")[0]
        seqs[name] = str(rec.seq).upper()
    return seqs


# ---------------------------------------------------------------------------
# edlib helpers
# ---------------------------------------------------------------------------


def norm_edit_distance(a, b):
    """Normalized edit distance: edits / longer(len(a), len(b))."""
    ed = edlib.align(a, b, mode="NW", task="distance")["editDistance"]
    return ed / max(len(a), len(b))
