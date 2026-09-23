"""Generic read/reference IO helpers shared by mapping_validation.py and
relative_abundance.py. Nothing in this module is experiment-specific -- every function
takes the path(s) it needs explicitly.
"""

import gzip
from pathlib import Path

import edlib
import pandas as pd
from Bio import SeqIO


def load_layout(layout_csv):
    df = pd.read_csv(layout_csv)
    df["sample_id"] = df.apply(
        lambda r: f"Plate{int(r['dest_plate']):02d}_{r['dest_well']}", axis=1
    )
    return df


def sample_fastq_path(demux_dir, plate, well):
    """Per-well fastqs are named Plate{plate:02d}_{well}.fastq (or .fastq.gz)."""
    stem = f"Plate{int(plate):02d}_{well}"
    for ext in (".fastq", ".fastq.gz"):
        p = Path(demux_dir) / f"{stem}{ext}"
        if p.exists():
            return p
    return None


def count_reads_fast(path):
    """Read count via line counting (4 lines/record) -- avoids full FASTQ parsing."""
    if path is None or not Path(path).exists() or Path(path).stat().st_size == 0:
        return 0
    path = Path(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        return sum(1 for _ in fh) // 4


def load_reads(path):
    """Return list of (read_id, sequence) tuples from a fastq(.gz))."""
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


def load_reference_db(path):
    """Load a reference FASTA -> {strain_name: sequence}.

    Bio.SeqIO's `.id` is already just the first whitespace-delimited token of the header,
    so this one rule handles every header convention seen so far without per-db
    configuration:
      - '>STRAIN_consensus winner-takes-all ...'        (corroborated_db)   -> strip '_consensus'
      - '>STRAIN|mono' / '>STRAIN'                       (merged_consensus)  -> split on '|'
      - '>STRAIN informative_wells=N support=N/N ...'    (per_strain_consensus builds) -> .id is
        already just 'STRAIN' (metadata is past the first whitespace, so both rules below leave
        it unchanged)
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


def norm_edit_distance(a, b):
    """Normalized edit distance: edits / longer(len(a), len(b))."""
    ed = edlib.align(a, b, mode="NW", task="distance")["editDistance"]
    return ed / max(len(a), len(b))


def edit_distance_bp(a, b):
    return edlib.align(a, b, mode="NW", task="distance")["editDistance"]


def pair_key(strain_a, strain_b):
    """Order-independent identifier for an unordered strain pair."""
    return tuple(sorted((strain_a, strain_b)))
