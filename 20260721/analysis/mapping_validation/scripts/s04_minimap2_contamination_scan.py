"""04 -- Unconstrained cross-check: map every read against *all* 87 strains in a
reference db (not just the 1-2 expected ones), using minimap2. This is the
complement to the constrained edlib mapping in s03: it can tell us not just "does
this read match what's expected" but, when it doesn't, "what does it actually look
like" -- the signature of cross-well contamination or a mislabeled well.

Output (per db):
  outputs/04_reference_normalized_<db>.fasta   strain-name-only headers, for minimap2
  outputs/04_minimap2_<db>.paf                 raw PAF alignments
  outputs/04_minimap2_read_besthit_<db>.csv    one row per read: best genome-wide hit
  outputs/04_reads_combined.fastq              all cohort reads, read name = sample__readid
"""

import subprocess

import pandas as pd

from common import OUT_DIR, REFERENCE_DBS, load_reference_db, load_reads


def write_combined_reads_fastq(samples_df, fastq_out, universe_out):
    """Write every cohort read to one FASTQ (for minimap2) and cache the full
    (sample_id, read_id, read_len) universe as a CSV so downstream steps don't
    need to re-parse every per-well FASTQ for every reference db."""
    if fastq_out.exists() and universe_out.exists():
        return fastq_out, pd.read_csv(universe_out)

    rows = []
    n = 0
    with open(fastq_out, "w") as fh:
        for _, r in samples_df.iterrows():
            for read_id, seq in load_reads(r["path"]):
                qname = f"{r['sample_id']}__{read_id}"
                fh.write(f"@{qname}\n{seq}\n+\n{'I' * len(seq)}\n")
                rows.append({"sample_id": r["sample_id"], "read_id": read_id, "read_len": len(seq)})
                n += 1
    universe = pd.DataFrame(rows)
    universe.to_csv(universe_out, index=False)
    print(f"wrote {n} reads -> {fastq_out} (read universe cached -> {universe_out})")
    return fastq_out, universe


def write_normalized_reference(ref_path, out_path):
    seqs = load_reference_db(ref_path)
    with open(out_path, "w") as fh:
        for strain, seq in seqs.items():
            fh.write(f">{strain}\n{seq}\n")
    return out_path


def run_minimap2(ref_fasta, reads_fastq, paf_out, threads=8):
    cmd = [
        "minimap2",
        "-x", "map-ont",
        "--secondary=no",
        "-t", str(threads),
        str(ref_fasta),
        str(reads_fastq),
    ]
    with open(paf_out, "w") as out_fh:
        proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"minimap2 failed:\n{proc.stderr}")
    return paf_out


PAF_COLS = [
    "qname", "qlen", "qstart", "qend", "strand",
    "tname", "tlen", "tstart", "tend",
    "nmatch", "alnlen", "mapq",
]


def parse_paf(paf_path):
    if paf_path.stat().st_size == 0:
        return pd.DataFrame(columns=PAF_COLS + ["identity"])
    df = pd.read_csv(
        paf_path, sep="\t", header=None, usecols=range(12), names=PAF_COLS, engine="c"
    )
    df["identity"] = df["nmatch"] / df["alnlen"]
    return df


def best_hit_table(samples_df, read_universe, paf_df):
    paf_df = paf_df.copy()
    paf_df[["sample_id", "read_id"]] = paf_df["qname"].str.split("__", n=1, expand=True)
    paf_df = paf_df.rename(columns={"tname": "best_hit_strain"})

    hit = paf_df[
        ["sample_id", "read_id", "best_hit_strain", "nmatch", "alnlen", "identity", "mapq"]
    ].copy()

    # outer-join against the full read universe so unmapped reads (absent from the PAF)
    # get explicit rows instead of being silently dropped
    hit = read_universe.merge(hit, on=["sample_id", "read_id"], how="left")
    hit["best_hit_strain"] = hit["best_hit_strain"].where(hit["best_hit_strain"].notna(), None)
    for c in ("nmatch", "alnlen", "mapq"):
        hit[c] = hit[c].fillna(0).astype(int)
    hit["identity"] = hit["identity"].fillna(0.0)

    expect = samples_df.set_index("sample_id")[["strain1", "strain2"]]
    hit = hit.join(expect, on="sample_id")
    hit["is_expected"] = (hit["best_hit_strain"] == hit["strain1"]) | (
        hit["best_hit_strain"] == hit["strain2"]
    )
    return hit


def run_contamination_scan(samples_df, read_universe, db_name, ref_path):
    norm_ref = write_normalized_reference(ref_path, OUT_DIR / f"04_reference_normalized_{db_name}.fasta")
    reads_fastq = OUT_DIR / "04_reads_combined.fastq"
    paf_out = OUT_DIR / f"04_minimap2_{db_name}.paf"
    run_minimap2(norm_ref, reads_fastq, paf_out)

    paf_df = parse_paf(paf_out)
    hit = best_hit_table(samples_df, read_universe, paf_df)

    out_path = OUT_DIR / f"04_minimap2_read_besthit_{db_name}.csv"
    hit.to_csv(out_path, index=False)

    n_reads = len(hit)
    n_mapped = (hit["best_hit_strain"].notna()).sum()
    n_expected = hit["is_expected"].sum()
    print(
        f"[{db_name}] {n_mapped}/{n_reads} reads mapped ({n_mapped/n_reads:.1%}); "
        f"{n_expected}/{n_reads} best-hit an expected strain ({n_expected/n_reads:.1%}) "
        f"-> {out_path}"
    )
    return hit


def run_all(samples_df):
    reads_fastq, read_universe = write_combined_reads_fastq(
        samples_df, OUT_DIR / "04_reads_combined.fastq", OUT_DIR / "04_read_universe.csv"
    )
    results = {}
    for db_name, ref_path in REFERENCE_DBS.items():
        results[db_name] = run_contamination_scan(samples_df, read_universe, db_name, ref_path)
    return results


if __name__ == "__main__":
    samples_df = pd.read_csv(OUT_DIR / "01_samples_gt5reads.csv")
    run_all(samples_df)
