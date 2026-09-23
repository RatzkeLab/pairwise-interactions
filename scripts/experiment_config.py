"""Per-experiment configuration for the shared mapping_validation / relative_abundance
pipelines.

Each experiment (20260721, 20260630, ...) gets one ExperimentConfig, built in that
experiment's own analysis/mapping_validation/config.py (or shared across both analyses if
convenient). Everything genuinely specific to an experiment lives here: paths, which fastq
directory holds the reads to analyze, which reference fasta(s) to use, and which one is
"primary" (self-consistent, trusted as ground truth) vs. external (needs an identity
cross-check before being trusted -- see mapping_validation.identity_cross_check).

Thresholds that are about the sequencing chemistry/statistics rather than the experiment
(ONT error rate, what counts as a confident match, etc.) are NOT part of this config -- they
live as defaults in mapping_validation.py / relative_abundance.py so they don't need to be
restated per experiment, but every function that uses one accepts it as an optional override.
"""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ExperimentConfig:
    name: str                          # e.g. "20260721"
    exp_base: Path                     # .../pairwise_interaction_experiments/<name>
    layout_csv: Path                   # strain layout CSV (dest_plate, dest_well, well_type, strain1, strain2, ...)
    demux_dir: Path                    # directory of Plate{NN}_{well}.fastq(.gz) files to analyze

    # mapping_validation: one or more reference dbs to check reads against
    reference_dbs: dict                # name -> Path
    primary_db: str                    # key into reference_dbs treated as ground truth (self-consistent)
    external_cross_check_dbs: list = field(default_factory=list)  # keys to identity-check against primary_db

    # relative_abundance: single reference used for interaction scoring
    ra_reference_fasta: Path = None

    min_reads: int = 5                 # "more than N reads" cohort threshold

    mapping_validation_out_dir: Path = None
    relative_abundance_out_dir: Path = None

    def __post_init__(self):
        if self.mapping_validation_out_dir is None:
            self.mapping_validation_out_dir = self.exp_base / "analysis" / "mapping_validation" / "outputs"
        if self.relative_abundance_out_dir is None:
            self.relative_abundance_out_dir = self.exp_base / "analysis" / "relative_abundance" / "outputs"
        self.mapping_validation_fig_dir = self.mapping_validation_out_dir / "figures"
        self.relative_abundance_fig_dir = self.relative_abundance_out_dir / "figures"
        for d in (
            self.mapping_validation_out_dir,
            self.mapping_validation_fig_dir,
            self.relative_abundance_out_dir,
            self.relative_abundance_fig_dir,
        ):
            d.mkdir(parents=True, exist_ok=True)
