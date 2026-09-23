"""Rewrite a minibar barcode TSV so the index columns hold the index.

minibar's `-cols` default for a 5-column file is sample, fwd index, fwd primer,
rev index, rev primer. The files this project generates put the whole
59 bp construct -- ligation adapter (15) + barcode (24) + 16S primer (20) -- in
the index column, and the primer again in the primer column. minibar therefore
matches a 59-mer within edit distance `-e`, instead of the 24-mer barcode it was
meant to match. The tolerance is the same number of edits spread over 2.5x the
length, so most reads fail.

Measured on one 2763-read chunk of 20260907, -e 6 -E 6 -l 150:
    as generated : 236 reads assigned   (8.5%)
    corrected    : 1482 reads assigned  (53.6%)

All three runs (20260630, 20260721, 20260907) have the same file layout.
"""
import sys
import pandas as pd
from pathlib import Path

ADAPTER_LEN = 15          # ATCGCCTACCGTGAC


def fix(path_in, path_out):
    df = pd.read_csv(path_in, sep="\t")
    for ix_col, pr_col in (("FwIndex", "FwPrimer"), ("RvIndex", "RwPrimer")):
        primer_len = int(df[pr_col].str.len().mode()[0])
        full_len = int(df[ix_col].str.len().mode()[0])
        bc_len = full_len - ADAPTER_LEN - primer_len
        if bc_len <= 0:
            raise ValueError(f"{ix_col}: nothing left after stripping adapter+primer; "
                             f"file may already be corrected")
        # verify the construct really is adapter + barcode + primer before cutting
        assert (df[ix_col].str[-primer_len:] == df[pr_col]).all(), \
            f"{ix_col} does not end in {pr_col}; layout is not what this script assumes"
        df[ix_col] = df[ix_col].str[ADAPTER_LEN:-primer_len]
        print(f"  {ix_col}: {full_len} -> {bc_len} bp "
              f"(stripped {ADAPTER_LEN} bp adapter + {primer_len} bp primer)")
    assert df.FwIndex.str.len().nunique() == 1 and df.RvIndex.str.len().nunique() == 1
    df.to_csv(path_out, sep="\t", index=False)
    print(f"  wrote {len(df)} rows -> {path_out}")
    return df


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: fix_barcode_file.py <in.tsv> <out.tsv>")
    fix(Path(sys.argv[1]), Path(sys.argv[2]))
