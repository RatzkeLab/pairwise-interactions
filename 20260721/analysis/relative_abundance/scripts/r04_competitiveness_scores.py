"""04 -- Per-strain competitiveness score: average each strain's interaction across every
distinct opponent it was tested against.

Built from r03 (one row per unique pair), not the raw per-well table -- so a pair with 3
replicate wells doesn't get 3x the weight of a pair only tested once. For each strain we
report:

  - competitiveness_score        mean log2(self/opponent) ratio across all its opponents.
                                  Positive = tends to dominate; negative = tends to lose;
                                  near zero = roughly even on average. This is the headline
                                  "competitiveness score" requested.
  - sd_log2_ratio                 dispersion of that same per-opponent ratio. This is a
                                  DIFFERENT axis from the mean, and answers a different
                                  question: low sd clustered near 0 = a true "coexister"
                                  (consistently ties everyone); low sd far from 0 = a
                                  consistent dominator/loser; high sd = context-dependent
                                  (crushes some opponents, loses to others) -- a strain can
                                  have the same mean as another strain for entirely
                                  different reasons, which is why both are reported.
  - win_fraction                  fraction of opponents against which this strain has
                                  mean_relative_abundance > 0.5
  - frac_high_uncertainty_opponents  fraction of its opponent-pairs where the pair itself is
                                  hard to call (near-identical references) -- a strain that
                                  happens to sit in a cluster of near-twin strains will have
                                  a competitiveness score built from fewer truly-informative
                                  comparisons than n_opponents suggests.

Output: outputs/r04_strain_competitiveness.csv
"""

import pandas as pd

from ra_common import OUT_DIR, load_layout


def orient_to_strain(df, strain):
    is_a = df["strain_a"] == strain
    out = df.copy()
    out["opponent"] = out["strain_b"].where(is_a, out["strain_a"])
    out["log2_ratio_self_over_opponent"] = out["mean_log2_ratio_a_over_b"].where(is_a, -out["mean_log2_ratio_a_over_b"])
    out["relative_abundance_self"] = out["mean_relative_abundance_a"].where(is_a, 1 - out["mean_relative_abundance_a"])
    return out


def compute_competitiveness():
    pairs = pd.read_csv(OUT_DIR / "r03_pair_replicate_stats.csv")
    layout = load_layout()
    all_strains = sorted(set(layout["strain1"]) | set(layout["strain2"]))

    rows = []
    scored_strains = set(pairs["strain_a"]) | set(pairs["strain_b"])
    for strain in all_strains:
        mask = (pairs["strain_a"] == strain) | (pairs["strain_b"] == strain)
        g = orient_to_strain(pairs[mask], strain)
        if len(g) == 0:
            rows.append({"strain": strain, "n_opponents": 0, "n_wells": 0})
            continue
        rows.append(
            {
                "strain": strain,
                "n_opponents": len(g),
                "n_wells": int(g["n_replicates"].sum()),
                "competitiveness_score": g["log2_ratio_self_over_opponent"].mean(),
                "sd_log2_ratio": g["log2_ratio_self_over_opponent"].std(),
                "mean_relative_abundance": g["relative_abundance_self"].mean(),
                "win_fraction": (g["relative_abundance_self"] > 0.5).mean(),
                "frac_high_uncertainty_opponents": g["high_uncertainty_pair"].mean(),
            }
        )

    df = pd.DataFrame(rows).sort_values("competitiveness_score", ascending=False)
    out_path = OUT_DIR / "r04_strain_competitiveness.csv"
    df.to_csv(out_path, index=False)

    unscored = [s for s in all_strains if s not in scored_strains]
    print(f"{len(df) - len(unscored)}/{len(all_strains)} strains scored")
    if unscored:
        print(f"  NOT scored (no resolvable reference / no tested pairs survived filtering): {unscored}")
    scored = df.dropna(subset=["competitiveness_score"])
    print(f"\ntop 5 most competitive:\n{scored.head(5)[['strain','competitiveness_score','win_fraction','n_opponents']].to_string(index=False)}")
    print(f"\nbottom 5 least competitive:\n{scored.tail(5)[['strain','competitiveness_score','win_fraction','n_opponents']].to_string(index=False)}")
    print(f"\nsaved -> {out_path}")
    return df


if __name__ == "__main__":
    compute_competitiveness()
