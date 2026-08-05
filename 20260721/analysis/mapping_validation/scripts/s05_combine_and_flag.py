"""05 -- Combine the constrained edlib mapping (s03) and the unconstrained minimap2
scan (s04) into one master per-well summary, and assign a QC call to every well.

Logic for "is strain X present in this well":
  A data-driven "confident match" threshold is found from the pooled distribution of
  best_dist (constrained edlib, PRIMARY_DB only) the same way check_reads.ipynb finds its
  same/different-strain threshold: look for the empty valley between the low mode (reads
  that really are the reference strain, differing only by ONT error) and the high mode
  (reads that are something else, forced onto the wrong reference because it's the only
  expected candidate). This threshold is learned once, from PRIMARY_DB, and then applied
  to all three dbs -- fitting it separately per db would let corroborated_db's ~65
  name-collision strains (see below) drag its own valley down and silently suppress every
  confident call there, even for the strains whose name genuinely does match.
  A strain then "counts" as present if it has >= max(MIN_STRAIN_READS, MIN_STRAIN_FRAC *
  n_reads) confident reads -- the same rule used in check_reads.ipynb / per_strain_consensus.ipynb.

Which reference is "primary" -- an important finding from s02:
  Strain names here are short plate-well-style labels ("D13", "G5", ...) that get reused
  across many unrelated plate layouts over time. s02's cross-check shows that for only
  21/86 strains does corroborated_db's sequence for a given name actually agree with this
  experiment's own consensus for that name (< 10% normalized edit distance) -- the other
  65 are almost certainly a different organism that happens to share the name. corroborated_db
  is therefore NOT a safe ground truth here. merged_consensus_mono_priority (built entirely
  from this experiment's own reads, so its labels are self-consistent by construction) is
  used as the PRIMARY reference for qc_status; corroborated_db is retained as a secondary,
  informational cross-check, and reported restricted to strains s02 flagged reliable.

Outputs:
  outputs/05_combined_sample_summary.csv   one row per well: fractions from every db,
                                            confident-match calls, and a qc_status label
  outputs/05_contamination_candidates.csv  wells where minimap2's genome-wide best hit
                                            (primary db) is dominated by a strain that
                                            isn't strain1/strain2 -- candidates for
                                            cross-well contamination or a mislabeled well.
                                            Each is also annotated with how similar the
                                            "contaminant" strain's own reference is to the
                                            well's expected strain(s): some of the top
                                            offenders (e.g. K6/P17, 1.5% apart) are reference
                                            near-twins that 16S genuinely cannot discriminate,
                                            not evidence of real contamination.
"""

from math import ceil

import numpy as np
import pandas as pd

from common import OUT_DIR, REFERENCE_DBS, load_reference_db, norm_edit_distance

MIN_STRAIN_READS = 3
MIN_STRAIN_FRAC = 0.10

# self-consistent by construction (built from this experiment's own reads) -- see module
# docstring for why corroborated_db is not used as the primary reference here
PRIMARY_DB = "merged_consensus_mono_priority"
CONTAM_MIN_FRAC = 0.20
CONTAM_MIN_READS = 3

# below this normalized edit distance, two strains' own references are so similar that
# 16S alone can't tell their reads apart -- a "contamination" call against such a strain
# is a reference-ambiguity artifact, not evidence of a real off-target strain
REFERENCE_TWIN_THRESHOLD = 0.05


def find_valley_threshold(values, lo=0.04, hi=0.20, bins=81, empty_frac=0.02):
    """Data-driven same/different threshold: centre of the emptiest band in [lo, hi]
    of the pooled best-match distance histogram (same approach as check_reads.ipynb)."""
    values = np.asarray(values)
    values = values[~np.isnan(values)]
    hist, edges = np.histogram(values, bins=np.linspace(0, max(0.4, values.max()), bins))
    centers = 0.5 * (edges[:-1] + edges[1:])
    band = (centers >= lo) & (centers <= hi)
    empty = band & (hist <= max(1, empty_frac * hist.max()))
    if empty.any():
        return round(float(np.median(centers[empty])), 3)
    return round(float(np.median(values)), 3)  # fallback, shouldn't normally trigger


def summarize_edlib_db(db_name, threshold):
    path = OUT_DIR / f"03_edlib_read_assignments_{db_name}.csv.gz"
    reads = pd.read_csv(path)
    reads["confident"] = reads["best_dist"] < threshold

    rows = []
    for sid, g in reads.groupby("sample_id"):
        strain1, strain2 = g["strain1"].iat[0], g["strain2"].iat[0]
        n = len(g)
        candidates_present = set(g["candidate1_strain"]) | (
            set(g["candidate2_strain"].dropna()) if "candidate2_strain" in g else set()
        )
        gc = g[g["confident"]]

        def frac_for(strain):
            if strain not in candidates_present:
                return np.nan
            return (gc["best_match_strain"] == strain).sum() / n

        def n_for(strain):
            if strain not in candidates_present:
                return np.nan
            return int((gc["best_match_strain"] == strain).sum())

        rows.append(
            {
                "sample_id": sid,
                f"edlib_{db_name}_threshold": threshold,
                f"edlib_{db_name}_n_conf_strain1": n_for(strain1),
                f"edlib_{db_name}_n_conf_strain2": n_for(strain2),
                f"edlib_{db_name}_frac_strain1": frac_for(strain1),
                f"edlib_{db_name}_frac_strain2": frac_for(strain2),
                f"edlib_{db_name}_median_best_dist": g["best_dist"].median(),
            }
        )
    return pd.DataFrame(rows)


def summarize_minimap2_db(db_name):
    path = OUT_DIR / f"04_minimap2_read_besthit_{db_name}.csv"
    hits = pd.read_csv(path)

    rows = []
    for sid, g in hits.groupby("sample_id"):
        n = len(g)
        frac_expected = g["is_expected"].sum() / n
        frac_unmapped = g["best_hit_strain"].isna().sum() / n
        other = g[(~g["is_expected"]) & g["best_hit_strain"].notna()]
        if len(other):
            vc = other["best_hit_strain"].value_counts()
            top_other_strain, top_other_count = vc.index[0], int(vc.iat[0])
        else:
            top_other_strain, top_other_count = None, 0
        rows.append(
            {
                "sample_id": sid,
                f"mm2_{db_name}_frac_expected": frac_expected,
                f"mm2_{db_name}_frac_unmapped": frac_unmapped,
                f"mm2_{db_name}_top_other_strain": top_other_strain,
                f"mm2_{db_name}_top_other_count": top_other_count,
                f"mm2_{db_name}_top_other_frac": top_other_count / n,
            }
        )
    return pd.DataFrame(rows)


def assign_qc_status(row):
    n = row["n_reads"]
    cutoff = max(MIN_STRAIN_READS, ceil(MIN_STRAIN_FRAC * n))
    n1 = row[f"edlib_{PRIMARY_DB}_n_conf_strain1"]
    n2 = row[f"edlib_{PRIMARY_DB}_n_conf_strain2"]
    present1 = (not pd.isna(n1)) and n1 >= cutoff
    present2 = (not pd.isna(n2)) and n2 >= cutoff

    if row["well_type"] == "mono":
        status = "mono_confirmed" if present1 else "mono_low_confidence"
    else:
        if present1 and present2:
            status = "pair_both_confirmed"
        elif present1 or present2:
            status = "pair_single_dominant"
        else:
            status = "pair_neither_confirmed"

    contaminated = (
        row.get(f"mm2_{PRIMARY_DB}_top_other_frac", 0) >= CONTAM_MIN_FRAC
        and row.get(f"mm2_{PRIMARY_DB}_top_other_count", 0) >= CONTAM_MIN_READS
    )
    return pd.Series({"qc_cutoff_reads": cutoff, "qc_status": status, "qc_contaminated": contaminated})


def combine_and_flag():
    samples = pd.read_csv(OUT_DIR / "01_samples_gt5reads.csv")
    summary = samples.copy()

    # Learn the confident-match threshold ONCE, from the primary (self-consistent) db's
    # pooled distances, and apply it to every db. Fitting it separately per db is wrong:
    # corroborated_db's pool is dominated by ~65 name-collision strains (see s02), which
    # drags its own valley down to ~0.05 and silently zeroes out almost every "confident"
    # call -- even for the strains where the name genuinely does match.
    primary_reads = pd.read_csv(OUT_DIR / f"03_edlib_read_assignments_{PRIMARY_DB}.csv.gz")
    shared_threshold = find_valley_threshold(primary_reads["best_dist"].values)
    print(f"confident-match threshold (learned from {PRIMARY_DB}, applied to all dbs): {shared_threshold}")

    for db_name in REFERENCE_DBS:
        summary = summary.merge(summarize_edlib_db(db_name, shared_threshold), on="sample_id", how="left")
        summary = summary.merge(summarize_minimap2_db(db_name), on="sample_id", how="left")

    summary = pd.concat([summary, summary.apply(assign_qc_status, axis=1)], axis=1)

    # annotate whether corroborated_db's entry for strain1/strain2 is trustworthy for this
    # experiment (see s02's identity cross-check) -- corroborated_db columns should only be
    # read at face value where these are True
    cross_check = pd.read_csv(OUT_DIR / "02_reference_cross_check.csv").set_index("strain")["reliable"]
    summary["corrdb_strain1_reliable"] = summary["strain1"].map(cross_check)
    summary["corrdb_strain2_reliable"] = summary["strain2"].map(cross_check)

    out_path = OUT_DIR / "05_combined_sample_summary.csv"
    summary.to_csv(out_path, index=False)

    contam = summary[summary["qc_contaminated"]].sort_values(
        f"mm2_{PRIMARY_DB}_top_other_frac", ascending=False
    ).copy()

    # is the "contaminant" strain actually just a reference near-twin of an expected strain?
    ref_seqs = load_reference_db(REFERENCE_DBS[PRIMARY_DB])
    twin_dist_cache = {}

    def twin_dist(row):
        other = row[f"mm2_{PRIMARY_DB}_top_other_strain"]
        if other is None or pd.isna(other) or other not in ref_seqs:
            return np.nan
        dists = []
        for expected in (row["strain1"], row["strain2"]):
            if expected in ref_seqs:
                key = tuple(sorted((expected, other)))
                if key not in twin_dist_cache:
                    twin_dist_cache[key] = norm_edit_distance(ref_seqs[expected], ref_seqs[other])
                dists.append(twin_dist_cache[key])
        return min(dists) if dists else np.nan

    contam["contaminant_ref_dist_to_expected"] = contam.apply(twin_dist, axis=1)
    contam["likely_reference_twin"] = contam["contaminant_ref_dist_to_expected"] < REFERENCE_TWIN_THRESHOLD

    contam_cols = [
        "sample_id", "well_type", "strain1", "strain2", "n_reads",
        f"mm2_{PRIMARY_DB}_top_other_strain", f"mm2_{PRIMARY_DB}_top_other_count",
        f"mm2_{PRIMARY_DB}_top_other_frac", "contaminant_ref_dist_to_expected",
        "likely_reference_twin", "qc_status",
    ]
    contam_path = OUT_DIR / "05_contamination_candidates.csv"
    contam[contam_cols].to_csv(contam_path, index=False)

    n_twin = int(contam["likely_reference_twin"].sum())
    print(f"combined summary: {len(summary)} wells -> {out_path}")
    print("\nqc_status counts:")
    print(summary["qc_status"].value_counts().to_string())
    print(
        f"\ncontamination-flagged wells: {len(contam)} -> {contam_path}\n"
        f"  of which {n_twin} are reference near-twins (<{REFERENCE_TWIN_THRESHOLD} dist to expected strain's own "
        f"reference -- 16S can't discriminate these, likely not real contamination)\n"
        f"  {len(contam) - n_twin} involve a genuinely divergent unexpected strain -- higher-priority follow-up"
    )
    return summary


if __name__ == "__main__":
    combine_and_flag()
