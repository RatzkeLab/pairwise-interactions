"""03 -- Constrained mapping: for every read in every sample (well) with > 5 reads,
compute the normalized edit distance (edlib) to *only* the reference sequence(s) of the
strain(s) the layout says should be there (strain1 / strain2), for each reference db.

This directly answers "do the reads match what we expect to be there" -- it does not
consider the other ~85 strains, by design (see s04 for the unconstrained check).

Output: outputs/03_edlib_read_assignments_<db>.csv.gz  (one per reference db)
    one row per (sample, read): distance to strain1 ref, distance to strain2 ref
    (NaN where that strain has no reference in this db, or the well is mono),
    the winning strain, its distance, and the margin to the runner-up.
"""

import pandas as pd

from common import OUT_DIR, REFERENCE_DBS, load_reference_db, load_reads, norm_edit_distance


def map_sample_constrained(read_seqs, candidates, ref_seqs):
    """candidates: list of distinct expected strain names present in ref_seqs (1 or 2).
    Returns list of dicts, one per read, with distances to each candidate + winner."""
    rows = []
    for read_id, seq in read_seqs:
        dists = {c: norm_edit_distance(seq, ref_seqs[c]) for c in candidates}
        best_strain = min(dists, key=dists.get)
        best_dist = dists[best_strain]
        sorted_d = sorted(dists.values())
        margin = sorted_d[1] - sorted_d[0] if len(sorted_d) > 1 else float("nan")
        row = {
            "read_id": read_id,
            "read_len": len(seq),
            "best_match_strain": best_strain,
            "best_dist": best_dist,
            "margin_to_runner_up": margin,
        }
        for i, c in enumerate(candidates, start=1):
            row[f"candidate{i}_strain"] = c
            row[f"candidate{i}_dist"] = dists[c]
        rows.append(row)
    return rows


def run_constrained_mapping(samples_df, db_name, ref_path):
    ref_seqs = load_reference_db(ref_path)
    out_rows = []
    n_no_ref = 0
    for _, r in samples_df.iterrows():
        candidates = list(dict.fromkeys([r["strain1"], r["strain2"]]))  # dedupe, keep order
        candidates = [c for c in candidates if c in ref_seqs]
        if not candidates:
            n_no_ref += 1
            continue
        reads = load_reads(r["path"])
        for row in map_sample_constrained(reads, candidates, ref_seqs):
            row["sample_id"] = r["sample_id"]
            row["well_type"] = r["well_type"]
            row["strain1"] = r["strain1"]
            row["strain2"] = r["strain2"]
            row["n_candidates"] = len(candidates)
            out_rows.append(row)

    df = pd.DataFrame(out_rows)
    out_path = OUT_DIR / f"03_edlib_read_assignments_{db_name}.csv.gz"
    df.to_csv(out_path, index=False)
    print(
        f"[{db_name}] {len(samples_df) - n_no_ref}/{len(samples_df)} wells had a reference "
        f"for >=1 expected strain ({n_no_ref} skipped, no reference available); "
        f"{len(df)} reads classified -> {out_path}"
    )
    return df


def run_all(samples_df):
    results = {}
    for db_name, ref_path in REFERENCE_DBS.items():
        results[db_name] = run_constrained_mapping(samples_df, db_name, ref_path)
    return results


if __name__ == "__main__":
    samples_df = pd.read_csv(OUT_DIR / "01_samples_gt5reads.csv")
    run_all(samples_df)
