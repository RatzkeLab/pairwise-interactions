"""02 -- Per pair-well relative-abundance / interaction scoring.

For every read in every pair well (>5 reads), compute the normalized edit distance to both
candidate strain references (merged_consensus_20260721.fasta). Every read is classified as:

  - "off_target"  if its best (smaller) distance is still worse than OFF_TARGET_DIST_THRESHOLD
                  -- it doesn't look like *either* expected strain (contamination / chimera /
                  junk). This directly answers "sequences that really don't align to either
                  of the two strains", quantified per well.
  - "ambiguous"   if it's on-target (matches one of the two well) but the *margin* between its
                  two distances is too small to call a winner -- the two candidates are
                  equally good explanations for this read. When the two reference sequences
                  are (near-)identical, EVERY read ends up here by construction, which is
                  exactly how "can't be told apart -> 50/50" falls out of the same general
                  method rather than needing a special case: with zero separating signal,
                  the classifier can't do better than a coin flip, so it correctly reports
                  "don't know" for every read. We still always attempt the classification
                  first (even for pairs only a few bp apart) rather than assuming indistin-
                  guishability up front from reference distance alone -- see r01 for the
                  reference distance, used here only as context, not as a gate.
  - "strain1"/"strain2" otherwise, whichever reference is closer.

Per well, this gives:
  - off_target_frac                      quantifies unexplained reads
  - uncertainty_score                     fraction of on-target reads that couldn't be called
  - relative_abundance_strain1            (n_strain1 + 0.5*n_ambiguous) / n_on_target -- the
                                           headline interaction score, ambiguous reads split
                                           evenly since there's no directional information in
                                           them (this is what naturally yields exactly 0.5 for
                                           an indistinguishable pair, since n_strain1=n_strain2=0)
  - relative_abundance_strain1_confident   same, but excluding ambiguous reads entirely
  - log2_ratio_strain1_over_strain2        log2((n_strain1+0.5)/(n_strain2+0.5)), a symmetric
                                           score well-suited to averaging across many pairwise
                                           opponents (see r04)

Output:
  outputs/r02_read_assignments.csv.gz     one row per read: distances, margin, class
  outputs/r02_well_interaction_scores.csv one row per well: full scoring above
"""

import numpy as np
import pandas as pd

from ra_common import (
    OUT_DIR, REFERENCE_FASTA, SAMPLES_GT5READS_CSV,
    OFF_TARGET_DIST_THRESHOLD, AMBIGUOUS_MARGIN_THRESHOLD,
    load_reference_db, load_reads, norm_edit_distance, pair_key,
)


def classify_well_reads(reads, ref1, ref2):
    rows = []
    for read_id, seq in reads:
        d1 = norm_edit_distance(seq, ref1)
        d2 = norm_edit_distance(seq, ref2)
        best = min(d1, d2)
        margin = abs(d1 - d2)
        if best > OFF_TARGET_DIST_THRESHOLD:
            cls = "off_target"
        elif margin < AMBIGUOUS_MARGIN_THRESHOLD:
            cls = "ambiguous"
        else:
            cls = "strain1" if d1 < d2 else "strain2"
        rows.append(
            {"read_id": read_id, "dist_strain1": d1, "dist_strain2": d2, "margin": margin, "read_class": cls}
        )
    return rows


def score_well(counts, n_reads):
    n1 = counts.get("strain1", 0)
    n2 = counts.get("strain2", 0)
    n_amb = counts.get("ambiguous", 0)
    n_off = counts.get("off_target", 0)
    n_on_target = n1 + n2 + n_amb

    off_target_frac = n_off / n_reads
    uncertainty_score = (n_amb / n_on_target) if n_on_target > 0 else np.nan
    rel_ab = ((n1 + 0.5 * n_amb) / n_on_target) if n_on_target > 0 else np.nan
    rel_ab_confident = (n1 / (n1 + n2)) if (n1 + n2) > 0 else np.nan
    log2_ratio = np.log2((n1 + 0.5) / (n2 + 0.5))

    return {
        "n_strain1": n1, "n_strain2": n2, "n_ambiguous": n_amb, "n_off_target": n_off,
        "n_on_target": n_on_target,
        "off_target_frac": off_target_frac,
        "uncertainty_score": uncertainty_score,
        "relative_abundance_strain1": rel_ab,
        "relative_abundance_strain1_confident": rel_ab_confident,
        "log2_ratio_strain1_over_strain2": log2_ratio,
    }


def compute_interaction_scores():
    samples = pd.read_csv(SAMPLES_GT5READS_CSV)
    pair_wells = samples[samples["well_type"] == "pair"].copy()
    ref_seqs = load_reference_db(REFERENCE_FASTA)

    ref_dist = pd.read_csv(OUT_DIR / "r01_reference_pair_distances.csv")
    ref_dist_map = {
        pair_key(r.strain_a, r.strain_b): (r.bp_dist, r.norm_dist)
        for r in ref_dist.itertuples()
    }

    read_rows = []
    well_rows = []
    n_missing_ref = 0

    for r in pair_wells.itertuples():
        strain1, strain2 = r.strain1, r.strain2
        ref1, ref2 = ref_seqs.get(strain1), ref_seqs.get(strain2)

        bp_dist, norm_dist = ref_dist_map.get(pair_key(strain1, strain2), (None, None))
        base_row = {
            "sample_id": r.sample_id, "strain1": strain1, "strain2": strain2, "n_reads": r.n_reads,
            "ref_pair_bp_dist": bp_dist, "ref_pair_norm_dist": norm_dist,
        }

        if ref1 is None or ref2 is None:
            n_missing_ref += 1
            well_rows.append({**base_row, "missing_reference": True})
            continue

        reads = load_reads(r.path)
        classified = classify_well_reads(reads, ref1, ref2)
        for row in classified:
            row["sample_id"] = r.sample_id
        read_rows.extend(classified)

        counts = pd.Series([row["read_class"] for row in classified]).value_counts().to_dict()
        well_rows.append({**base_row, "missing_reference": False, **score_well(counts, r.n_reads)})

    reads_df = pd.DataFrame(read_rows)
    wells_df = pd.DataFrame(well_rows)

    # canonical (order-independent) orientation, so replicate wells of the same pair --
    # which may list strain1/strain2 in either order -- can be compared/averaged directly
    def _canonical(row):
        a, b = pair_key(row["strain1"], row["strain2"])
        flip = row["strain1"] != a
        return pd.Series(
            {
                "strain_a": a,
                "strain_b": b,
                "relative_abundance_a": (1 - row["relative_abundance_strain1"]) if flip else row["relative_abundance_strain1"],
                "log2_ratio_a_over_b": (-row["log2_ratio_strain1_over_strain2"]) if flip else row["log2_ratio_strain1_over_strain2"],
            }
        )

    wells_df = pd.concat([wells_df, wells_df.apply(_canonical, axis=1)], axis=1)

    reads_out = OUT_DIR / "r02_read_assignments.csv.gz"
    wells_out = OUT_DIR / "r02_well_interaction_scores.csv"
    reads_df.to_csv(reads_out, index=False)
    wells_df.to_csv(wells_out, index=False)

    resolved = wells_df[~wells_df["missing_reference"]]
    print(f"{len(pair_wells)} pair wells ({n_missing_ref} skipped, missing a reference for one strain)")
    print(f"{len(reads_df)} reads classified -> {reads_out}")
    print(f"read_class breakdown:\n{reads_df['read_class'].value_counts().to_string()}")
    print(f"\nmedian off_target_frac: {resolved['off_target_frac'].median():.3f}")
    print(f"median uncertainty_score: {resolved['uncertainty_score'].median():.3f}")
    print(f"saved -> {wells_out}")
    return wells_df, reads_df


if __name__ == "__main__":
    compute_interaction_scores()
