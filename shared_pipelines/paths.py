"""Every path that points outside this repository, in one place.

Eleven scripts previously hardcoded `/home/rl/scripts/karl/...` strings. That made the analyses
unmovable and meant a renamed input directory would fail in nine places with nine different
error messages. Import from here instead.

Deliberately NOT moved here: numeric thresholds like `MIN_RESOLVABLE_BP`,
`HIGH_UNCERTAINTY_THRESHOLD`, `CORR_THRESHOLD`. Each of those sits in its module beside the
comment or calibration table that justifies its value -- `MIN_RESOLVABLE_BP = 10` is next to the
mono-well read-accuracy table that produced the number. Centralising the values would separate
them from their reasoning, which is the opposite of what this repository needs.
"""

from pathlib import Path

# --- external data roots -----------------------------------------------------------------
KARL_ROOT = Path("/home/rl/scripts/karl")
LINK_TO_KARL = KARL_ROOT / "Link to Karl"

GENOMIC_TABLES = LINK_TO_KARL / "final_genomic_tables"
PLATE_READER = LINK_TO_KARL / "plate_reader_csvs" / "data_ascii" / "Karl_2026"
CONSENSUS_MERGE = KARL_ROOT / "merge_consensus_sequences"

# --- specific files referenced by name ------------------------------------------------------
STRAIN_MAPPING = GENOMIC_TABLES / "mapping_384_well_plate_collection.csv"
KEGG_KO = GENOMIC_TABLES / "KEGG_ko_and_strains_table.csv"
KEGG_MODULE = GENOMIC_TABLES / "KEGG_Module_and_strains_table.csv"
CAZY = GENOMIC_TABLES / "CAZy_and_strains_table.csv"
BIGG = GENOMIC_TABLES / "BiGG_and_strains_table.csv"          # NOT a metabolic network -- see TRAPS.md
PFAMS = GENOMIC_TABLES / "PFAMs_and_strains_table.csv"
PANX = GENOMIC_TABLES / "panx_and_strains_table_final_post_filtering_strain_names.csv"

CORROBORATED_DB = (CONSENSUS_MERGE / "collapse_naive_updated3_15diff" /
                   "corroborated_db_filtered_min5.fasta")
CORROBORATED_DB_EDITED = (CONSENSUS_MERGE / "collapse_naive_updated3_15diff" /
                          "corroborated_db_filtered_min5_edited.fasta")   # identical content

# --- plate-reader folders, per experiment ---------------------------------------------------
# NOTE: filenames are SAVE times; the internal `Date:` header is the RUN START. See TRAPS.md.
OD_FULL = {
    "20260630": PLATE_READER / "Karl_20260704_OD_Full",     # 30 destination plates, 61-wl spectra
    "20260721": PLATE_READER / "Karl_20260723_OD_Full",     # 30 plates + one 't1' test plate
}
OD_PREP = {
    "20260630": PLATE_READER / "Karl_20260623_OD",          # preculture / source-plate reads
    "20260721": PLATE_READER / "Karl_20260722_OD",
}


def check(verbose=True):
    """Report which of these exist. Cheap first call for any script that uses them."""
    items = {k: v for k, v in globals().items()
             if isinstance(v, Path) and not k.startswith("_")}
    missing = {k: v for k, v in items.items() if not v.exists()}
    for group in (OD_FULL, OD_PREP):
        for k, v in group.items():
            if not v.exists():
                missing[f"OD[{k}]"] = v
    if verbose:
        print(f"paths: {len(items) + 4 - len(missing)} present, {len(missing)} missing")
        for k, v in missing.items():
            print(f"  MISSING {k}: {v}")
    return missing
