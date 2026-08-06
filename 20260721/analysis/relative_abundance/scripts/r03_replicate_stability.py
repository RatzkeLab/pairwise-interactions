"""03 -- Replicate stability: for strain pairs that were plated more than once, is the
interaction score consistent across replicate wells, or does it change?

Every well's relative_abundance_a / log2_ratio_a_over_b is already in the same canonical
(strain_a, strain_b) orientation (see r02), so replicate wells of the same pair -- which may
list strain1/strain2 in either order -- are directly comparable without re-orienting here.

We report replicate spread (std, range) alongside each pair's *own* measurement uncertainty
(mean uncertainty_score across its replicates) so that "unstable" can be told apart into two
different explanations: genuine well-to-well biological variability (low uncertainty, but the
score still moves a lot) vs. simply a hard-to-call pair (near-identical references, so of
course individual replicates bounce around -- the instability is in the measurement, not the
biology).

Output: outputs/r03_pair_replicate_stats.csv
    one row per unique tested strain pair (including singletons, for completeness) with
    replicate count, mean/std/range of the interaction score, and stability flags.
"""

import numpy as np
import pandas as pd

from ra_common import OUT_DIR

HIGH_UNCERTAINTY_THRESHOLD = 0.3   # mean uncertainty_score above this -> "hard to call" pair
UNSTABLE_STD_THRESHOLD = 0.15      # std of relative_abundance_a above this -> "unstable" pair


def compute_replicate_stability():
    wells = pd.read_csv(OUT_DIR / "r02_well_interaction_scores.csv")
    wells = wells[~wells["missing_reference"]].copy()

    rows = []
    for (a, b), g in wells.groupby(["strain_a", "strain_b"]):
        n = len(g)
        row = {
            "strain_a": a, "strain_b": b,
            "n_replicates": n,
            "ref_pair_bp_dist": g["ref_pair_bp_dist"].iloc[0],
            "mean_n_reads": g["n_reads"].mean(),
            "mean_relative_abundance_a": g["relative_abundance_a"].mean(),
            "std_relative_abundance_a": g["relative_abundance_a"].std() if n > 1 else np.nan,
            "range_relative_abundance_a": (g["relative_abundance_a"].max() - g["relative_abundance_a"].min()) if n > 1 else np.nan,
            "mean_log2_ratio_a_over_b": g["log2_ratio_a_over_b"].mean(),
            "std_log2_ratio_a_over_b": g["log2_ratio_a_over_b"].std() if n > 1 else np.nan,
            "mean_uncertainty_score": g["uncertainty_score"].mean(),
            "mean_off_target_frac": g["off_target_frac"].mean(),
        }
        row["high_uncertainty_pair"] = row["mean_uncertainty_score"] > HIGH_UNCERTAINTY_THRESHOLD
        row["unstable_replicate"] = (
            (n > 1) and (row["std_relative_abundance_a"] > UNSTABLE_STD_THRESHOLD)
        )
        rows.append(row)

    df = pd.DataFrame(rows).sort_values("n_replicates", ascending=False)
    out_path = OUT_DIR / "r03_pair_replicate_stats.csv"
    df.to_csv(out_path, index=False)

    replicated = df[df["n_replicates"] > 1]
    print(f"{len(df)} unique pairs total; {len(replicated)} have >=2 replicate wells")
    if len(replicated):
        clean = replicated[~replicated["high_uncertainty_pair"]]
        print(f"  of which {len(clean)} are NOT high-uncertainty (reliable measurement)")
        print(f"  median std(relative_abundance_a) across replicated pairs: {replicated['std_relative_abundance_a'].median():.3f}")
        print(f"  median std, excluding high-uncertainty pairs: {clean['std_relative_abundance_a'].median():.3f}")
        print(f"  flagged unstable (std > {UNSTABLE_STD_THRESHOLD}): {int(replicated['unstable_replicate'].sum())}")
        print(f"  flagged unstable AND not high-uncertainty (genuine biological instability): "
              f"{int((replicated['unstable_replicate'] & ~replicated['high_uncertainty_pair']).sum())}")
    print(f"saved -> {out_path}")
    return df


if __name__ == "__main__":
    compute_replicate_stability()
