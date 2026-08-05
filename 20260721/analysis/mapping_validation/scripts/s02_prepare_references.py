"""02 -- Load the three candidate reference databases, report how completely each one
covers the 87 strains used in this experiment, AND check whether corroborated_db's
strain names actually refer to the same organism this experiment calls by that name.

Why the cross-check matters: strain names in this project are short plate-well-style
labels (e.g. "D13", "G5") that get reused across many unrelated plate layouts over the
years. corroborated_db is pooled from many past experiments, so a name match there is
not automatically an identity match -- it could be a different isolate that happened to
land on a well with the same coordinate-style name in some other experiment. We check
this directly: for every strain, compare corroborated_db's sequence to
merged_consensus_mono_priority's sequence for that same name (the latter is built
entirely from *this* experiment's own reads, so its strain labels are self-consistent
by construction). A large edit distance between the two means the corroborated_db entry
almost certainly is NOT the same organism, despite the matching name.

Output:
  outputs/02_reference_coverage.csv     strain x db, sequence length (or missing)
  outputs/02_reference_cross_check.csv  strain, corroborated_db vs. merged_consensus_mono_priority
                                         normalized edit distance, and a `reliable` flag
"""

import pandas as pd

from common import REFERENCE_DBS, OUT_DIR, load_layout, load_reference_db, norm_edit_distance

CROSS_CHECK_RELIABLE_THRESHOLD = 0.10  # normalized edit distance; below this = same organism


def prepare_references():
    layout = load_layout()
    strains = sorted(set(layout["strain1"]) | set(layout["strain2"]))

    dbs = {name: load_reference_db(path) for name, path in REFERENCE_DBS.items()}

    rows = []
    for strain in strains:
        row = {"strain": strain}
        for db_name, seqs in dbs.items():
            row[f"{db_name}_len"] = len(seqs[strain]) if strain in seqs else None
        rows.append(row)
    cov = pd.DataFrame(rows)

    out_path = OUT_DIR / "02_reference_coverage.csv"
    cov.to_csv(out_path, index=False)

    print(f"strains in layout: {len(strains)}")
    for db_name, seqs in dbs.items():
        present = sum(1 for s in strains if s in seqs)
        extra = len(seqs) - present
        print(
            f"  {db_name:32s} {REFERENCE_DBS[db_name]}\n"
            f"    covers {present}/{len(strains)} strains"
            + (f", +{extra} extra entries not in this layout" if extra else "")
        )
    print(f"saved -> {out_path}")

    # ---- cross-check: does corroborated_db's name match this experiment's identity? ----
    corr = dbs["corroborated_db"]
    mcmp = dbs["merged_consensus_mono_priority"]
    cc_rows = []
    for strain in strains:
        if strain in corr and strain in mcmp:
            d = norm_edit_distance(corr[strain], mcmp[strain])
            cc_rows.append({"strain": strain, "norm_dist_corrdb_vs_own_consensus": d})
    cc = pd.DataFrame(cc_rows).sort_values("norm_dist_corrdb_vs_own_consensus").reset_index(drop=True)
    cc["reliable"] = cc["norm_dist_corrdb_vs_own_consensus"] < CROSS_CHECK_RELIABLE_THRESHOLD

    cc_path = OUT_DIR / "02_reference_cross_check.csv"
    cc.to_csv(cc_path, index=False)

    n_reliable = int(cc["reliable"].sum())
    print(
        f"\ncorroborated_db identity cross-check ({len(cc)} strains comparable):\n"
        f"  {n_reliable}/{len(cc)} strain names agree with this experiment's own consensus "
        f"(< {CROSS_CHECK_RELIABLE_THRESHOLD} normalized edit distance)\n"
        f"  {len(cc) - n_reliable}/{len(cc)} strain names in corroborated_db do NOT match this "
        f"experiment's sequence for that name -- likely coincidental reuse of a plate-well-style "
        f"label from an unrelated experiment, not the same organism.\n"
        f"  -> treating merged_consensus_mono_priority as the PRIMARY reference for QC calls; "
        f"corroborated_db is reported only as a secondary check, restricted to reliable strains."
    )
    print(f"saved -> {cc_path}")
    return dbs, cov, cc


if __name__ == "__main__":
    prepare_references()
