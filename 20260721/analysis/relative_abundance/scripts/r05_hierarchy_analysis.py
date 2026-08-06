"""05 -- Is this strain library organized as a dominance hierarchy, or a tangled / cyclic
competitive network?

Three independent, complementary tests, all standard tools from behavioral-ecology
dominance-hierarchy analysis, applied here to microbial pairwise competition:

  1. Bradley-Terry model -- fit a single latent "strength" per strain such that
     P(a beats b) = logistic(strength_a - strength_b), via weighted logistic regression
     (statsmodels GLM, binomial family) on pooled confident-read counts per pair. The
     model's pseudo-R^2 (deviance explained) is a direct measure of "how well does a single
     1-D ranking explain who wins" -- high = hierarchical, low = not.
  2. Directional Consistency Index (DCI; de Vries 1995) -- rank strains by total wins, then
     ask what fraction of pairwise outcomes agree with that ranking. Built for exactly this
     situation: an incomplete tournament (not every strain pair was tested), unlike the
     classic Landau's h which assumes a full round robin.
  3. Intransitive-triad fraction -- among all strain triples where all three pairwise
     outcomes are known, what fraction form a cycle (A beats B beats C beats A) rather than
     a strict order? Cycles are the signature of "rock-paper-scissors" competition that no
     1-D hierarchy can represent, however good its fit.

Also computes, per strain, a BT-residual dispersion: how far a strain's actual outcome
against each opponent deviates from what a single strength value predicts (observed
log2-ratio minus (strength_self - strength_opponent)). This matters because the *raw*
spread of a strain's log2-ratios across opponents is dominated by how strong or weak those
opponents happen to be, not by the focal strain's own consistency -- a strain that obeys the
hierarchy perfectly will still show a wide raw spread if it happens to face both very strong
and very weak opponents. The BT residual factors that out, isolating genuine
context-dependence / idiosyncrasy from simple opponent-strength variation.

Output:
  outputs/r05_pairwise_relative_abundance_matrix.csv  strain x strain, mean_relative_abundance_a
  outputs/r05_bt_strengths.csv                        strain, bt_strength, bt_se, rank,
                                                       bt_residual_mean, bt_residual_sd
  outputs/r05_hierarchy_summary.csv                    the headline numbers from all 3 tests
"""

import itertools

import numpy as np
import pandas as pd
import statsmodels.api as sm

from ra_common import OUT_DIR


def union_find_components(strains, pairs):
    parent = {s: s for s in strains}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for a, b in pairs:
        union(a, b)

    comps = {}
    for s in strains:
        comps.setdefault(find(s), []).append(s)
    return sorted(comps.values(), key=len, reverse=True)


def fit_bradley_terry(pair_reads, strains):
    """pair_reads: rows with strain_a, strain_b, n_reads_a, n_reads_b (pooled confident
    read counts across replicates, in canonical a/b orientation)."""
    strains = sorted(strains)
    baseline, others = strains[0], strains[1:]
    idx = {s: i for i, s in enumerate(others)}

    X = np.zeros((len(pair_reads), len(others)))
    for row_i, r in enumerate(pair_reads.itertuples()):
        if r.strain_a in idx:
            X[row_i, idx[r.strain_a]] = 1
        if r.strain_b in idx:
            X[row_i, idx[r.strain_b]] = -1

    y = pair_reads[["n_reads_a", "n_reads_b"]].to_numpy()
    model = sm.GLM(y, X, family=sm.families.Binomial())
    fit = model.fit()

    strength = pd.Series({baseline: 0.0}, dtype=float)
    strength = pd.concat([strength, pd.Series(fit.params, index=others)])
    se = pd.Series({baseline: np.nan})
    se = pd.concat([se, pd.Series(fit.bse, index=others)])

    pseudo_r2 = 1 - fit.deviance / fit.null_deviance
    return strength, se, pseudo_r2, fit


def compute_directional_consistency(win_matrix, rank_order):
    """DCI (de Vries 1995): rank strains by total wins, then measure what fraction of
    realized outcomes agree with that ranking."""
    rank = {s: i for i, s in enumerate(rank_order)}  # lower index = higher rank (more wins)
    consistent, inconsistent = 0, 0
    for a, b in itertools.combinations(win_matrix.index, 2):
        w = win_matrix.loc[a, b]
        if pd.isna(w) or w == 0.5:
            continue
        winner = a if w == 1 else b
        loser = b if w == 1 else a
        if rank[winner] < rank[loser]:
            consistent += 1
        else:
            inconsistent += 1
    total = consistent + inconsistent
    dci = (consistent - inconsistent) / total if total else np.nan
    return dci, consistent, inconsistent


def count_intransitive_triads(win_matrix):
    strains = list(win_matrix.index)
    n_evaluated = 0
    n_cyclic = 0
    for a, b, c in itertools.combinations(strains, 3):
        wab, wbc, wac = win_matrix.loc[a, b], win_matrix.loc[b, c], win_matrix.loc[a, c]
        if pd.isna(wab) or pd.isna(wbc) or pd.isna(wac):
            continue
        if wab == 0.5 or wbc == 0.5 or wac == 0.5:
            continue
        n_evaluated += 1
        # transitive iff the "beats" relation among {a,b,c} has a total order;
        # a 3-cycle is the only non-transitive possibility for 3 items
        a_beats_b = wab == 1
        b_beats_c = wbc == 1
        a_beats_c = wac == 1
        # cyclic patterns: a>b>c>a  or  a<b<c<a
        cyclic = (a_beats_b and b_beats_c and not a_beats_c) or (
            (not a_beats_b) and (not b_beats_c) and a_beats_c
        )
        if cyclic:
            n_cyclic += 1
    frac = n_cyclic / n_evaluated if n_evaluated else np.nan
    return frac, n_cyclic, n_evaluated


def hierarchy_analysis():
    wells = pd.read_csv(OUT_DIR / "r02_well_interaction_scores.csv")
    wells = wells[~wells["missing_reference"]].copy()

    # pool confident read counts per unique pair, in canonical a/b orientation
    wells["n_reads_a"] = np.where(wells["strain1"] == wells["strain_a"], wells["n_strain1"], wells["n_strain2"])
    wells["n_reads_b"] = np.where(wells["strain1"] == wells["strain_a"], wells["n_strain2"], wells["n_strain1"])
    pooled = (
        wells.groupby(["strain_a", "strain_b"])
        .agg(n_reads_a=("n_reads_a", "sum"), n_reads_b=("n_reads_b", "sum"))
        .reset_index()
    )
    pooled = pooled[(pooled["n_reads_a"] + pooled["n_reads_b"]) > 0]

    strains = sorted(set(pooled["strain_a"]) | set(pooled["strain_b"]))
    components = union_find_components(strains, list(zip(pooled["strain_a"], pooled["strain_b"])))
    main_component = set(components[0])
    print(f"{len(strains)} strains with resolvable pairwise data")
    print(f"connected components: {[len(c) for c in components]}")
    if len(components) > 1:
        excluded = [s for c in components[1:] for s in c]
        print(f"  fitting hierarchy only on the largest component; excluded (too few connections): {excluded}")

    pooled_main = pooled[pooled["strain_a"].isin(main_component) & pooled["strain_b"].isin(main_component)]

    # ---- 1. Bradley-Terry ----
    strength, se, pseudo_r2, fit = fit_bradley_terry(pooled_main, main_component)
    bt_df = pd.DataFrame({"strain": strength.index, "bt_strength": strength.values, "bt_se": se.reindex(strength.index).values})
    bt_df = bt_df.sort_values("bt_strength", ascending=False).reset_index(drop=True)
    bt_df["bt_rank"] = np.arange(1, len(bt_df) + 1)
    print(f"\nBradley-Terry pseudo-R^2 (deviance explained by a single 1-D strength ranking): {pseudo_r2:.3f}")

    # ---- residual dispersion: how much does a strain's outcome deviate from what a single
    # BT strength value predicts, beyond just facing opponents of different strength? Plain
    # sd(log2_ratio) across opponents is dominated by *opponent* strength spread, not the
    # focal strain's own consistency -- this residual isolates the focal strain's part.
    rel_ab_pre = pd.read_csv(OUT_DIR / "r03_pair_replicate_stats.csv")
    rel_ab_main = rel_ab_pre[rel_ab_pre["strain_a"].isin(main_component) & rel_ab_pre["strain_b"].isin(main_component)].copy()
    rel_ab_main["predicted_log2_ratio_a_over_b"] = rel_ab_main["strain_a"].map(strength) - rel_ab_main["strain_b"].map(strength)
    rel_ab_main["bt_residual_a_over_b"] = rel_ab_main["mean_log2_ratio_a_over_b"] - rel_ab_main["predicted_log2_ratio_a_over_b"]

    residual_rows = []
    for s in main_component:
        is_a = rel_ab_main["strain_a"] == s
        is_b = rel_ab_main["strain_b"] == s
        resid_self = pd.concat([rel_ab_main.loc[is_a, "bt_residual_a_over_b"], -rel_ab_main.loc[is_b, "bt_residual_a_over_b"]])
        residual_rows.append({"strain": s, "bt_residual_mean": resid_self.mean(), "bt_residual_sd": resid_self.std()})
    bt_df = bt_df.merge(pd.DataFrame(residual_rows), on="strain", how="left")
    bt_df.to_csv(OUT_DIR / "r05_bt_strengths.csv", index=False)
    print(
        f"median BT-residual sd across strains: {bt_df['bt_residual_sd'].median():.2f} log2 units "
        "(the genuine, opponent-strength-corrected unpredictability of a strain's outcomes)"
    )

    # ---- pairwise relative-abundance matrix (for figures + the discrete win matrix) ----
    rel_ab = rel_ab_pre
    strains_sorted = sorted(main_component)
    mat = pd.DataFrame(np.nan, index=strains_sorted, columns=strains_sorted)
    win = pd.DataFrame(np.nan, index=strains_sorted, columns=strains_sorted)
    for r in rel_ab.itertuples():
        if r.strain_a not in main_component or r.strain_b not in main_component:
            continue
        mat.loc[r.strain_a, r.strain_b] = r.mean_relative_abundance_a
        mat.loc[r.strain_b, r.strain_a] = 1 - r.mean_relative_abundance_a
        w = 1.0 if r.mean_relative_abundance_a > 0.5 else (0.0 if r.mean_relative_abundance_a < 0.5 else 0.5)
        win.loc[r.strain_a, r.strain_b] = w
        win.loc[r.strain_b, r.strain_a] = 1 - w if w != 0.5 else 0.5
    mat.to_csv(OUT_DIR / "r05_pairwise_relative_abundance_matrix.csv")

    # ---- 2. Directional Consistency Index ----
    total_wins = (win == 1).sum(axis=1)
    rank_order = total_wins.sort_values(ascending=False).index.tolist()
    dci, n_consistent, n_inconsistent = compute_directional_consistency(win, rank_order)
    print(f"Directional Consistency Index (rank = total wins): {dci:.3f}  "
          f"({n_consistent} consistent, {n_inconsistent} upsets)")

    # ---- 3. Intransitive triads ----
    frac_cyclic, n_cyclic, n_triads = count_intransitive_triads(win)
    print(f"Intransitive (cyclic) triads: {n_cyclic}/{n_triads} evaluated ({frac_cyclic:.3f})")

    summary = pd.DataFrame(
        [
            {"metric": "n_strains_in_main_component", "value": len(main_component)},
            {"metric": "n_components", "value": len(components)},
            {"metric": "bt_pseudo_r2", "value": pseudo_r2},
            {"metric": "dci", "value": dci},
            {"metric": "dci_n_consistent", "value": n_consistent},
            {"metric": "dci_n_inconsistent_upsets", "value": n_inconsistent},
            {"metric": "frac_intransitive_triads", "value": frac_cyclic},
            {"metric": "n_intransitive_triads", "value": n_cyclic},
            {"metric": "n_triads_evaluated", "value": n_triads},
        ]
    )
    summary.to_csv(OUT_DIR / "r05_hierarchy_summary.csv", index=False)
    print(f"\nsaved -> r05_pairwise_relative_abundance_matrix.csv, r05_bt_strengths.csv, r05_hierarchy_summary.csv")
    return bt_df, mat, win, summary


if __name__ == "__main__":
    hierarchy_analysis()
