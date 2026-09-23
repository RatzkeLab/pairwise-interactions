"""ExperimentConfig for 20260907, plus the paths specific to this run.

Two things differ from 20260630/20260721 and both matter downstream:

  - the layout covers 16 plates but only plates 1-8 were amplified and sequenced.
    `plate_role` marks which; every shared module filters on it. The 2464 frozen
    wells carry real barcodes and no DNA, which makes them a standing negative
    control for demultiplexing mis-assignment.

  - every one of the 384 strains has two monoculture wells, on different plates.
    strain_consensus.py uses those as the anchor rather than a cross-well majority.

DEMUX_DIR points at the re-demultiplexed output. The original run at
.../20260907_demux assigned 13.7% of filtered reads because the barcode file gave
minibar a 59 bp adapter+barcode+primer string as the index; see
02_demux_qc/fix_barcode_file.py. Point DEMUX_DIR_ORIGINAL at the old output only
to reproduce the QC that found it.
"""
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE.parent / "shared_scripts"))

from experiment_config import ExperimentConfig   # noqa: E402

DEMUX_ROOT = Path("/home/rl/data/interim/karl/demultiplexing")
DEMUX_DIR = DEMUX_ROOT / "20260907_demux_corrected" / "unflipped"
DEMUX_DIR_ORIGINAL = DEMUX_ROOT / "20260907_demux" / "unflipped"

REFERENCE_DBS_DIR = Path("/home/rl/scripts/karl/data_links/resources/reference_database_16s")
GENOMIC_TABLES = Path("/home/rl/scripts/karl/data_links/resources/final_genomic_tables")
PLATE_READER = Path("/home/rl/scripts/karl/data_links/raw/plate_reader_csvs/data_ascii/Karl_2026")

OD_FULL = PLATE_READER / "Karl_20260910_ODFull"   # destination plates, 61-wavelength spectra
OD_PREP = PLATE_READER / "Karl_20260908_OD"       # preculture / source-plate reads

CONSENSUS_DIR = BASE / "03_recreate_reference_db" / "outputs"
QC_DIR = BASE / "04_qc" / "outputs"


def make_config(demux_dir=None):
    return ExperimentConfig(
        name="20260907",
        exp_base=BASE,
        layout_csv=BASE / "01_setup" / "strain_layout_20260907.csv",
        demux_dir=demux_dir or DEMUX_DIR,
        reference_dbs={
            "corroborated": REFERENCE_DBS_DIR / "corroborated_db_filtered_min5.fasta",
            "full_collection": REFERENCE_DBS_DIR / "full_collection_db.fasta",
            "rDNA_all_strains": REFERENCE_DBS_DIR / "rDNA_16S_db_all_strains.fasta",
        },
        primary_db="full_collection",
        external_cross_check_dbs=["corroborated", "rDNA_all_strains"],
        min_reads=10,
    )
