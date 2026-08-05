"""01 -- Gather all wells from the strain layout, resolve their demux fastq, count reads,
and keep the subset with more than MIN_READS reads (the "complete analysis" cohort).

Output: outputs/01_samples_gt5reads.csv
    one row per well, with sample_id, plate/well, well_type, strain1, strain2,
    fastq path, and n_reads.
"""

from common import OUT_DIR, MIN_READS
from common import load_layout, sample_fastq_path, count_reads_fast


def gather_samples():
    layout = load_layout()
    layout["path"] = layout.apply(
        lambda r: sample_fastq_path(r["dest_plate"], r["dest_well"]), axis=1
    )
    n_missing = layout["path"].isna().sum()
    layout["n_reads"] = layout["path"].apply(count_reads_fast)

    gt = layout[layout["n_reads"] > MIN_READS].copy()
    gt["path"] = gt["path"].astype(str)

    cols = [
        "sample_id",
        "dest_plate",
        "dest_row",
        "dest_col",
        "dest_well",
        "well_type",
        "strain1",
        "strain2",
        "n_reads",
        "path",
    ]
    gt = gt[cols].sort_values("sample_id").reset_index(drop=True)

    out_path = OUT_DIR / "01_samples_gt5reads.csv"
    gt.to_csv(out_path, index=False)

    print(f"layout wells               : {len(layout)}")
    print(f"  missing fastq files      : {n_missing}")
    print(f"wells with > {MIN_READS} reads      : {len(gt)}")
    print(f"  well_type breakdown      : {gt['well_type'].value_counts().to_dict()}")
    print(f"  total reads in cohort    : {gt['n_reads'].sum()}")
    print(f"saved -> {out_path}")
    return gt


if __name__ == "__main__":
    gather_samples()
