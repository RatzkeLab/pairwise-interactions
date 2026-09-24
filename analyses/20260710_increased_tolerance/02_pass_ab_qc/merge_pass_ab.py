"""Merge the separately demultiplexed pass_a and pass_b into one demux directory.

minibar assigns every read independently, so concatenating each sample's two unflipped FASTQs
is identical to demultiplexing both input folders together, without the ~2 h pipeline run.
Only done after demux_ab_qc.py showed the two runs are the same pool (see
experiments/20260630/demux_runs.md).

Writes <DEMUX_ROOT>/20260710_demultiplex_increased_tolerance_merged/
    unflipped/<sample>.fastq      pass_a reads, then pass_b reads (empty files kept)
    summary/*.tsv                 raw/filtered counts concatenated; demux_summary summed
    MERGE_INFO.txt                provenance

    python merge_pass_ab.py       (refuses to overwrite an existing merge)
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402

A = config.DEMUX_DIRS["corrected"].parent
B = config.DEMUX_DIRS["pass_b"].parent
OUT = config.DEMUX_ROOT / "20260710_demultiplex_increased_tolerance_merged"

if OUT.exists():
    sys.exit(f"{OUT} exists; delete it first if you really want to rebuild")
(OUT / "unflipped").mkdir(parents=True)
(OUT / "summary").mkdir()

samples = sorted(p.name for p in (A / "unflipped").glob("*.fastq"))
assert samples == sorted(p.name for p in (B / "unflipped").glob("*.fastq")), "sample sets differ"
for name in samples:
    with open(OUT / "unflipped" / name, "wb") as out:
        for src in (A, B):
            out.write((src / "unflipped" / name).read_bytes())

for f in ("raw_read_counts.tsv", "filtered_read_counts.tsv"):
    pd.concat([pd.read_csv(s / "summary" / f, sep="\t") for s in (A, B)]).to_csv(
        OUT / "summary" / f, sep="\t", index=False)
d = pd.concat([pd.read_csv(s / "summary" / "demux_summary.tsv", sep="\t") for s in (A, B)])
d = d.groupby("sample", as_index=False)[["reads", "bases"]].sum()
d["mean_len"] = (d.bases / d.reads.where(d.reads > 0)).round(2)
d.to_csv(OUT / "summary" / "demux_summary.tsv", sep="\t", index=False)

(OUT / "MERGE_INFO.txt").write_text(
    f"per-sample concatenation of\n  {A}\n  {B}\n"
    f"made by {Path(__file__).resolve()}\n"
    "both demultiplexed with identical settings (corrected barcode file, -e 4);\n"
    "merge justified by analyses/20260710_increased_tolerance/02_pass_ab_qc/demux_ab_qc.py\n")

# check: every sample's merged read count is the sum of its parts
chk = d.set_index("sample").reads
for s in (A, B):
    chk = chk - pd.read_csv(s / "summary" / "demux_summary.tsv", sep="\t").set_index("sample").reads
assert (chk == 0).all()
n_lines = sum(1 for _ in open(OUT / "unflipped" / "Plate01_B10.fastq"))
assert n_lines // 4 == d.set_index("sample").reads["Plate01_B10"]
print(f"merged {len(samples)} samples, {d.reads.sum():,} reads -> {OUT}")
