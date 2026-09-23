"""Paths and per-experiment settings for the strain-identity QC.

Deliberately a separate top-level folder rather than per-experiment: the whole point is to
compare the SAME strain label across experiments and against external references, so nothing
here belongs to one experiment's analysis tree.
"""

from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KARL = Path("/home/rl/scripts/karl")
GENOMIC = KARL / "Link to Karl" / "final_genomic_tables"

MAPPING_CSV = GENOMIC / "mapping_384_well_plate_collection.csv"
GENOME_16S_FASTA = KARL / "Link to Karl" / "rDNA_16S_db_all_strains.fasta"
CORROBORATED_FASTA = (KARL / "merge_consensus_sequences" / "collapse_naive_updated3_15diff"
                      / "corroborated_db.fasta")

OUT = Path(__file__).resolve().parent / "outputs"
FIG = OUT / "figures"


@dataclass
class Exp:
    name: str
    layout_csv: Path
    demux_dir: Path
    pair_consensus: Path            # consensus built from this experiment's own wells
    mono_consensus: Path = None     # consensus built from mono wells only (may not exist yet)
    combined_fastq: Path = None     # all reads, one file, from mapping_validation step 04

    @property
    def out(self):
        d = OUT / self.name
        d.mkdir(parents=True, exist_ok=True)
        return d


EXPERIMENTS = {
    "20260721": Exp(
        name="20260721",
        layout_csv=ROOT / "20260721/setup/strain_layout_20260721_plate1_2_swapped.csv",
        demux_dir=ROOT / "20260721/demux/unflipped",
        pair_consensus=ROOT / "20260721/analysis/consensus/strain_consensus_20260721.fasta",
        mono_consensus=None,       # rebuilt here by the same code as 20260630 so the two are comparable;
        #                            the project's own 26-strain build is kept for cross-check

        combined_fastq=ROOT / "20260721/analysis/mapping_validation/outputs/04_reads_combined.fastq",
    ),
    "20260630": Exp(
        name="20260630",
        layout_csv=ROOT / "20260630/setup/strain_layout_20260630.csv",
        demux_dir=ROOT / "20260630/data/relevant_fastqs",
        pair_consensus=ROOT / "20260630/analysis/consensus2/strain_consensus_20260630.fasta",
        mono_consensus=None,        # never built for this experiment -- qc_sources builds it
        combined_fastq=ROOT / "20260630/analysis/mapping_validation/outputs/04_reads_combined.fastq",
    ),
}

# 20260721 shipped two layouts: someone already found and patched a plate-1/2 swap (616 wells
# differ, all on plates 1 and 2). Both are kept here because "was the patch right?" is one of
# the questions this QC exists to answer.
LAYOUT_VARIANTS_20260721 = {
    "as_patched_plate1_2_swapped": ROOT / "20260721/setup/strain_layout_20260721_plate1_2_swapped.csv",
    "original_unpatched": ROOT / "20260721/setup/strain_layout_20260721.csv",
}

MIN_READS = 5          # a well needs more than this to be called at all
