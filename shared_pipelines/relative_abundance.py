"""Shared, per-experiment-parameterized relative-abundance / competitive-hierarchy pipeline.

For every pair well, classifies each read against only that well's two expected strains
(edlib normalized edit distance) as strain1/strain2/ambiguous/off_target, turns that into a
per-well interaction score and uncertainty, checks replicate stability, aggregates a
per-strain competitiveness score, and fits a Bradley-Terry hierarchy model (+ Directional
Consistency Index + intransitive-triad fraction) to ask whether the library is organized as
a dominance hierarchy or a tangled/cyclic network. Originally built for
pairwise_interaction_experiments/20260721 as numbered scripts r01-r06; ported here so every
experiment can call the same, already-validated logic through an ExperimentConfig.

Reuses mapping_validation's find_valley_threshold() to learn the "off-target" distance
threshold fresh from THIS experiment's own reference (rather than a value copied from
another experiment's fit) -- see determine_off_target_threshold().

Call order: reference_distances -> well_interaction_scores -> replicate_stability ->
competitiveness_scores -> hierarchy_analysis -> make_all_figures.
"""

import itertools

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import statsmodels.api as sm

from io_utils import load_layout, load_reference_db, load_reads, norm_edit_distance, edit_distance_bp, pair_key
from mapping_validation import find_valley_threshold

# ---- thresholds: about ONT/16S sequencing statistics, not the experiment ----
AMBIGUOUS_MARGIN_THRESHOLD = 0.02   # margin between dist-to-strain1/2 below which a read is unassignable
HIGH_UNCERTAINTY_THRESHOLD = 0.3    # mean uncertainty_score above this -> "hard to call" pair
UNSTABLE_STD_THRESHOLD = 0.15       # std of relative_abundance_a above this -> "unstable" pair

# ---- fixed semantic figure palette ----
COLOR_BLUE = "#2a78d6"
COLOR_RED = "#e34948"
COLOR_GRID = "#d8d7d2"
COLOR_TEXT_SECONDARY = "#52514e"
COLOR_GOOD = "#0ca30c"
COLOR_CRITICAL = "#d03b3b"
DIVERGING_CMAP = LinearSegmentedColormap.from_list("blue_gray_red", [COLOR_BLUE, "#f0efec", COLOR_RED], N=256)

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white", "axes.edgecolor": COLOR_GRID,
    "axes.grid": True, "grid.color": COLOR_GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
})


def _samples_gt5reads(cfg):
    return pd.read_csv(cfg.mapping_validation_out_dir / "01_samples_gt5reads.csv")


# ===========================================================================
# 01 -- reference-pair distances for every unique tested strain pair
# ===========================================================================


def reference_distances(cfg):
    samples = _samples_gt5reads(cfg)
    pair_wells = samples[samples["well_type"] == "pair"].copy()
    pairs = sorted({pair_key(r.strain1, r.strain2) for r in pair_wells.itertuples()})
    ref_seqs = load_reference_db(cfg.ra_reference_fasta)

    rows = []
    for a, b in pairs:
        ref_a, ref_b = ref_seqs.get(a), ref_seqs.get(b)
        if ref_a is None or ref_b is None:
            rows.append({"strain_a": a, "strain_b": b, "bp_dist": None, "norm_dist": None,
                         "len_a": len(ref_a) if ref_a else None, "len_b": len(ref_b) if ref_b else None, "missing_reference": True})
            continue
        bp = edit_distance_bp(ref_a, ref_b)
        rows.append({"strain_a": a, "strain_b": b, "bp_dist": bp, "norm_dist": bp / max(len(ref_a), len(ref_b)),
                     "len_a": len(ref_a), "len_b": len(ref_b), "missing_reference": False})

    df = pd.DataFrame(rows)
    out_path = cfg.relative_abundance_out_dir / "r01_reference_pair_distances.csv"
    df.to_csv(out_path, index=False)

    n_missing = int(df["missing_reference"].sum())
    resolvable = df[~df["missing_reference"]]
    print(f"{len(df)} unique tested strain pairs ({n_missing} missing a reference for one strain)")
    print(f"reference-pair bp distance among resolvable pairs: min={resolvable.bp_dist.min():.0f}, "
          f"median={resolvable.bp_dist.median():.0f}, max={resolvable.bp_dist.max():.0f}")
    print(f"pairs with identical (0 bp) references: {(resolvable.bp_dist == 0).sum()}")
    print(f"pairs with <10 bp difference: {(resolvable.bp_dist < 10).sum()}")
    print(f"saved -> {out_path}")
    return df


# ===========================================================================
# 02 -- per-well interaction scoring
# ===========================================================================


def classify_well_reads(dists, off_target_threshold):
    """dists: list of (read_id, d1, d2) already-computed distances."""
    rows = []
    for read_id, d1, d2 in dists:
        best = min(d1, d2)
        margin = abs(d1 - d2)
        if best > off_target_threshold:
            cls = "off_target"
        elif margin < AMBIGUOUS_MARGIN_THRESHOLD:
            cls = "ambiguous"
        else:
            cls = "strain1" if d1 < d2 else "strain2"
        rows.append({"read_id": read_id, "dist_strain1": d1, "dist_strain2": d2, "margin": margin, "read_class": cls})
    return rows


def score_well(counts, n_reads):
    n1, n2 = counts.get("strain1", 0), counts.get("strain2", 0)
    n_amb, n_off = counts.get("ambiguous", 0), counts.get("off_target", 0)
    n_on_target = n1 + n2 + n_amb

    return {
        "n_strain1": n1, "n_strain2": n2, "n_ambiguous": n_amb, "n_off_target": n_off, "n_on_target": n_on_target,
        "off_target_frac": n_off / n_reads,
        "uncertainty_score": (n_amb / n_on_target) if n_on_target > 0 else np.nan,
        "relative_abundance_strain1": ((n1 + 0.5 * n_amb) / n_on_target) if n_on_target > 0 else np.nan,
        "relative_abundance_strain1_confident": (n1 / (n1 + n2)) if (n1 + n2) > 0 else np.nan,
        "log2_ratio_strain1_over_strain2": np.log2((n1 + 0.5) / (n2 + 0.5)),
    }


def determine_off_target_threshold(cfg, pair_wells, ref_seqs):
    """Data-driven off-target threshold, learned fresh from THIS experiment's own read-to-
    reference distances (reuses mapping_validation's valley-finding method) rather than a
    value carried over from another experiment's fit."""
    all_best = []
    sample_n = 0
    for r in pair_wells.itertuples():
        ref1, ref2 = ref_seqs.get(r.strain1), ref_seqs.get(r.strain2)
        if ref1 is None or ref2 is None:
            continue
        for _, seq in load_reads(r.path):
            d1, d2 = norm_edit_distance(seq, ref1), norm_edit_distance(seq, ref2)
            all_best.append(min(d1, d2))
        sample_n += 1
    threshold = find_valley_threshold(np.array(all_best))
    print(f"data-driven off-target threshold (learned from {sample_n} wells' read distances): {threshold}")
    return threshold


def compute_interaction_scores(cfg):
    samples = _samples_gt5reads(cfg)
    pair_wells = samples[samples["well_type"] == "pair"].copy()
    ref_seqs = load_reference_db(cfg.ra_reference_fasta)

    ref_dist = pd.read_csv(cfg.relative_abundance_out_dir / "r01_reference_pair_distances.csv")
    ref_dist_map = {pair_key(r.strain_a, r.strain_b): (r.bp_dist, r.norm_dist) for r in ref_dist.itertuples()}

    off_target_threshold = determine_off_target_threshold(cfg, pair_wells, ref_seqs)

    read_rows, well_rows = [], []
    n_missing_ref = 0

    for r in pair_wells.itertuples():
        ref1, ref2 = ref_seqs.get(r.strain1), ref_seqs.get(r.strain2)
        bp_dist, norm_dist = ref_dist_map.get(pair_key(r.strain1, r.strain2), (None, None))
        base_row = {"sample_id": r.sample_id, "strain1": r.strain1, "strain2": r.strain2, "n_reads": r.n_reads,
                    "ref_pair_bp_dist": bp_dist, "ref_pair_norm_dist": norm_dist}

        if ref1 is None or ref2 is None:
            n_missing_ref += 1
            well_rows.append({**base_row, "missing_reference": True})
            continue

        dists = [(read_id, norm_edit_distance(seq, ref1), norm_edit_distance(seq, ref2)) for read_id, seq in load_reads(r.path)]
        classified = classify_well_reads(dists, off_target_threshold)
        for row in classified:
            row["sample_id"] = r.sample_id
        read_rows.extend(classified)

        counts = pd.Series([row["read_class"] for row in classified]).value_counts().to_dict()
        well_rows.append({**base_row, "missing_reference": False, **score_well(counts, r.n_reads)})

    reads_df = pd.DataFrame(read_rows)
    wells_df = pd.DataFrame(well_rows)

    def _canonical(row):
        a, b = pair_key(row["strain1"], row["strain2"])
        flip = row["strain1"] != a
        return pd.Series({
            "strain_a": a, "strain_b": b,
            "relative_abundance_a": (1 - row["relative_abundance_strain1"]) if flip else row["relative_abundance_strain1"],
            "log2_ratio_a_over_b": (-row["log2_ratio_strain1_over_strain2"]) if flip else row["log2_ratio_strain1_over_strain2"],
        })

    wells_df = pd.concat([wells_df, wells_df.apply(_canonical, axis=1)], axis=1)

    reads_out = cfg.relative_abundance_out_dir / "r02_read_assignments.csv.gz"
    wells_out = cfg.relative_abundance_out_dir / "r02_well_interaction_scores.csv"
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


# ===========================================================================
# 03 -- replicate stability
# ===========================================================================


def replicate_stability(cfg):
    wells = pd.read_csv(cfg.relative_abundance_out_dir / "r02_well_interaction_scores.csv")
    wells = wells[~wells["missing_reference"]].copy()

    rows = []
    for (a, b), g in wells.groupby(["strain_a", "strain_b"]):
        n = len(g)
        row = {
            "strain_a": a, "strain_b": b, "n_replicates": n, "ref_pair_bp_dist": g["ref_pair_bp_dist"].iloc[0],
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
        row["unstable_replicate"] = (n > 1) and (row["std_relative_abundance_a"] > UNSTABLE_STD_THRESHOLD)
        rows.append(row)

    df = pd.DataFrame(rows).sort_values("n_replicates", ascending=False)
    out_path = cfg.relative_abundance_out_dir / "r03_pair_replicate_stats.csv"
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


# ===========================================================================
# 04 -- per-strain competitiveness score
# ===========================================================================


def orient_to_strain(df, strain):
    is_a = df["strain_a"] == strain
    out = df.copy()
    out["opponent"] = out["strain_b"].where(is_a, out["strain_a"])
    out["log2_ratio_self_over_opponent"] = out["mean_log2_ratio_a_over_b"].where(is_a, -out["mean_log2_ratio_a_over_b"])
    out["relative_abundance_self"] = out["mean_relative_abundance_a"].where(is_a, 1 - out["mean_relative_abundance_a"])
    return out


def competitiveness_scores(cfg):
    pairs = pd.read_csv(cfg.relative_abundance_out_dir / "r03_pair_replicate_stats.csv")
    layout = load_layout(cfg.layout_csv)
    all_strains = sorted(set(layout["strain1"]) | set(layout["strain2"]))

    rows = []
    scored_strains = set(pairs["strain_a"]) | set(pairs["strain_b"])
    for strain in all_strains:
        mask = (pairs["strain_a"] == strain) | (pairs["strain_b"] == strain)
        g = orient_to_strain(pairs[mask], strain)
        if len(g) == 0:
            rows.append({"strain": strain, "n_opponents": 0, "n_wells": 0})
            continue
        rows.append({
            "strain": strain, "n_opponents": len(g), "n_wells": int(g["n_replicates"].sum()),
            "competitiveness_score": g["log2_ratio_self_over_opponent"].mean(),
            "sd_log2_ratio": g["log2_ratio_self_over_opponent"].std(),
            "mean_relative_abundance": g["relative_abundance_self"].mean(),
            "win_fraction": (g["relative_abundance_self"] > 0.5).mean(),
            "frac_high_uncertainty_opponents": g["high_uncertainty_pair"].mean(),
        })

    df = pd.DataFrame(rows).sort_values("competitiveness_score", ascending=False)
    out_path = cfg.relative_abundance_out_dir / "r04_strain_competitiveness.csv"
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


# ===========================================================================
# 05 -- hierarchy analysis
# ===========================================================================


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
    rank = {s: i for i, s in enumerate(rank_order)}
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
    n_evaluated, n_cyclic = 0, 0
    for a, b, c in itertools.combinations(strains, 3):
        wab, wbc, wac = win_matrix.loc[a, b], win_matrix.loc[b, c], win_matrix.loc[a, c]
        if pd.isna(wab) or pd.isna(wbc) or pd.isna(wac):
            continue
        if wab == 0.5 or wbc == 0.5 or wac == 0.5:
            continue
        n_evaluated += 1
        a_beats_b, b_beats_c, a_beats_c = wab == 1, wbc == 1, wac == 1
        cyclic = (a_beats_b and b_beats_c and not a_beats_c) or ((not a_beats_b) and (not b_beats_c) and a_beats_c)
        if cyclic:
            n_cyclic += 1
    frac = n_cyclic / n_evaluated if n_evaluated else np.nan
    return frac, n_cyclic, n_evaluated


def hierarchy_analysis(cfg):
    wells = pd.read_csv(cfg.relative_abundance_out_dir / "r02_well_interaction_scores.csv")
    wells = wells[~wells["missing_reference"]].copy()

    wells["n_reads_a"] = np.where(wells["strain1"] == wells["strain_a"], wells["n_strain1"], wells["n_strain2"])
    wells["n_reads_b"] = np.where(wells["strain1"] == wells["strain_a"], wells["n_strain2"], wells["n_strain1"])
    pooled = wells.groupby(["strain_a", "strain_b"]).agg(n_reads_a=("n_reads_a", "sum"), n_reads_b=("n_reads_b", "sum")).reset_index()
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

    strength, se, pseudo_r2, fit = fit_bradley_terry(pooled_main, main_component)
    bt_df = pd.DataFrame({"strain": strength.index, "bt_strength": strength.values, "bt_se": se.reindex(strength.index).values})
    bt_df = bt_df.sort_values("bt_strength", ascending=False).reset_index(drop=True)
    bt_df["bt_rank"] = np.arange(1, len(bt_df) + 1)
    print(f"\nBradley-Terry pseudo-R^2 (deviance explained by a single 1-D strength ranking): {pseudo_r2:.3f}")

    rel_ab_pre = pd.read_csv(cfg.relative_abundance_out_dir / "r03_pair_replicate_stats.csv")
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
    bt_df.to_csv(cfg.relative_abundance_out_dir / "r05_bt_strengths.csv", index=False)
    print(f"median BT-residual sd across strains: {bt_df['bt_residual_sd'].median():.2f} log2 units "
          "(the genuine, opponent-strength-corrected unpredictability of a strain's outcomes)")

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
    mat.to_csv(cfg.relative_abundance_out_dir / "r05_pairwise_relative_abundance_matrix.csv")

    total_wins = (win == 1).sum(axis=1)
    rank_order = total_wins.sort_values(ascending=False).index.tolist()
    dci, n_consistent, n_inconsistent = compute_directional_consistency(win, rank_order)
    print(f"Directional Consistency Index (rank = total wins): {dci:.3f}  ({n_consistent} consistent, {n_inconsistent} upsets)")

    frac_cyclic, n_cyclic, n_triads = count_intransitive_triads(win)
    print(f"Intransitive (cyclic) triads: {n_cyclic}/{n_triads} evaluated ({frac_cyclic:.3f})")

    summary = pd.DataFrame([
        {"metric": "n_strains_in_main_component", "value": len(main_component)},
        {"metric": "n_components", "value": len(components)},
        {"metric": "bt_pseudo_r2", "value": pseudo_r2},
        {"metric": "dci", "value": dci},
        {"metric": "dci_n_consistent", "value": n_consistent},
        {"metric": "dci_n_inconsistent_upsets", "value": n_inconsistent},
        {"metric": "frac_intransitive_triads", "value": frac_cyclic},
        {"metric": "n_intransitive_triads", "value": n_cyclic},
        {"metric": "n_triads_evaluated", "value": n_triads},
    ])
    summary.to_csv(cfg.relative_abundance_out_dir / "r05_hierarchy_summary.csv", index=False)
    print("\nsaved -> r05_pairwise_relative_abundance_matrix.csv, r05_bt_strengths.csv, r05_hierarchy_summary.csv")
    return bt_df, mat, win, summary


# ===========================================================================
# 06 -- figures
# ===========================================================================


def _fig00_read_classification_quality(wells, fig_dir):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(wells["off_target_frac"].dropna(), bins=40, color=COLOR_BLUE, edgecolor="white", linewidth=0.3)
    axes[0].set_xlabel("off-target read fraction")
    axes[0].set_ylabel("# wells")
    axes[0].set_title("Reads matching neither expected strain")
    axes[1].hist(wells["uncertainty_score"].dropna(), bins=40, color=COLOR_BLUE, edgecolor="white", linewidth=0.3)
    axes[1].set_xlabel("uncertainty score (fraction of on-target reads unassignable)")
    axes[1].set_ylabel("# wells")
    axes[1].set_title("Per-well classification uncertainty")
    fig.suptitle("Read-classification quality across all pair wells", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(fig_dir / "00_read_classification_quality.png", dpi=150)
    plt.close(fig)


def _fig01_abundance_distribution(wells, fig_dir):
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(wells["relative_abundance_a"].dropna(), bins=np.linspace(0, 1, 41), color=COLOR_BLUE, edgecolor="white", linewidth=0.3)
    ax.axvline(0.5, color=COLOR_TEXT_SECONDARY, lw=1, ls="--", label="50/50")
    ax.set_xlabel("relative abundance of strain_a (0 = strain_b wins, 1 = strain_a wins)")
    ax.set_ylabel("# wells")
    ax.set_title("Interaction outcome across all pair wells")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "01_abundance_distribution.png", dpi=150)
    plt.close(fig)


def _fig02_replicate_stability(pairs, fig_dir):
    rep = pairs[pairs["n_replicates"] > 1].copy()
    fig, ax = plt.subplots(figsize=(7, 5))
    for flag, color, label in [(True, COLOR_RED, "hard-to-call pair (high uncertainty)"), (False, COLOR_BLUE, "well-resolved pair")]:
        sub = rep[rep["high_uncertainty_pair"] == flag]
        ax.scatter(sub["mean_uncertainty_score"], sub["std_relative_abundance_a"], s=16, color=color, alpha=0.6, edgecolor="none", label=f"{label} (n={len(sub)})")
    ax.axhline(0.15, color=COLOR_TEXT_SECONDARY, lw=1, ls="--", label="unstable cutoff = 0.15")
    ax.set_xlabel("mean uncertainty score across replicates")
    ax.set_ylabel("std of relative_abundance_a across replicates")
    ax.set_title("Replicate stability vs. measurement difficulty")
    fig.text(0.5, 0.92, "hard-to-call pairs look artificially 'stable' -- forced toward 0.5 every time, not because the biology reproduces",
              ha="center", fontsize=8, color=COLOR_TEXT_SECONDARY)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(fig_dir / "02_replicate_stability.png", dpi=150)
    plt.close(fig)


def _fig03_competitiveness_ranking(comp, fig_dir):
    df = comp.dropna(subset=["competitiveness_score"]).sort_values("competitiveness_score")
    vmax = df["competitiveness_score"].abs().max()
    colors = DIVERGING_CMAP((df["competitiveness_score"].values / vmax + 1) / 2)

    fig, ax = plt.subplots(figsize=(7, max(4, 0.16 * len(df))))
    ax.barh(range(len(df)), df["competitiveness_score"], color=colors)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["strain"], fontsize=6.5)
    ax.axvline(0, color=COLOR_TEXT_SECONDARY, lw=1)
    ax.set_xlabel("competitiveness score (mean log2 ratio vs. all tested opponents)")
    ax.set_title("Strain competitiveness ranking")
    fig.tight_layout()
    fig.savefig(fig_dir / "03_competitiveness_ranking.png", dpi=150)
    plt.close(fig)


def _fig04_mean_vs_dispersion(comp, bt, fig_dir):
    df = comp.merge(bt[["strain", "bt_residual_sd"]], on="strain").dropna(subset=["competitiveness_score", "bt_residual_sd"])
    fig, ax = plt.subplots(figsize=(7.5, 6))
    sc = ax.scatter(df["competitiveness_score"], df["bt_residual_sd"], s=24, c=df["frac_high_uncertainty_opponents"], cmap="Blues", edgecolor=COLOR_TEXT_SECONDARY, linewidth=0.3)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("fraction of opponents that are hard-to-call")
    ax.axvline(0, color=COLOR_TEXT_SECONDARY, lw=1, ls="--")
    ax.set_xlabel("competitiveness score (mean log2 ratio)")
    ax.set_ylabel("BT-residual dispersion (opponent-strength-corrected)")
    ax.set_title("Consistent dominance/loss vs. context-dependent outcomes")
    fig.text(0.5, 0.93, "low residual dispersion = behaves as its BT strength predicts; high = idiosyncratic, opponent-specific outcomes beyond the hierarchy",
              ha="center", fontsize=8, color=COLOR_TEXT_SECONDARY)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(fig_dir / "04_mean_vs_dispersion.png", dpi=150)
    plt.close(fig)


def _fig05_hierarchy_heatmap(mat, bt, fig_dir):
    order = bt.sort_values("bt_strength", ascending=False)["strain"].tolist()
    order = [s for s in order if s in mat.index]
    m = mat.loc[order, order]

    fig, ax = plt.subplots(figsize=(11, 10))
    im = ax.imshow(m.values, cmap=DIVERGING_CMAP, vmin=0, vmax=1, aspect="equal")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=5)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=5)
    ax.grid(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("relative abundance of row strain vs. column strain")
    ax.set_title("Pairwise outcomes, ordered by Bradley-Terry strength (strongest → weakest)")
    fig.tight_layout()
    fig.savefig(fig_dir / "05_hierarchy_heatmap.png", dpi=150)
    plt.close(fig)


def _fig06_bt_vs_naive(bt, comp, fig_dir):
    merged = bt.merge(comp, on="strain").dropna(subset=["bt_strength", "competitiveness_score"])
    corr = merged["bt_strength"].corr(merged["competitiveness_score"])

    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.scatter(merged["competitiveness_score"], merged["bt_strength"], s=20, color=COLOR_BLUE, alpha=0.7, edgecolor="none")
    ax.set_xlabel("naive competitiveness score (mean log2 ratio)")
    ax.set_ylabel("Bradley-Terry strength")
    ax.set_title(f"Two independent competitiveness estimates agree (r = {corr:.2f})")
    fig.tight_layout()
    fig.savefig(fig_dir / "06_bt_vs_naive_competitiveness.png", dpi=150)
    plt.close(fig)


def _fig07_hierarchy_consistency(summary, fig_dir):
    s = summary.set_index("metric")["value"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))

    consistent, upsets = s["dci_n_consistent"], s["dci_n_inconsistent_upsets"]
    axes[0].bar(["consistent", "upset"], [consistent, upsets], color=[COLOR_GOOD, COLOR_CRITICAL])
    axes[0].set_ylabel("# pairwise outcomes")
    axes[0].set_title(f"DCI = {s['dci']:.2f}")
    for i, v in enumerate([consistent, upsets]):
        axes[0].text(i, v, f"{int(v)}", ha="center", va="bottom")

    transitive = s["n_triads_evaluated"] - s["n_intransitive_triads"]
    cyclic = s["n_intransitive_triads"]
    axes[1].bar(["transitive", "cyclic"], [transitive, cyclic], color=[COLOR_GOOD, COLOR_CRITICAL])
    axes[1].set_ylabel("# strain triads")
    axes[1].set_title(f"cyclic fraction = {s['frac_intransitive_triads']:.3f}")
    for i, v in enumerate([transitive, cyclic]):
        axes[1].text(i, v, f"{int(v)}", ha="center", va="bottom")

    fig.suptitle(f"Hierarchy consistency (Bradley-Terry pseudo-R² = {s['bt_pseudo_r2']:.2f})", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(fig_dir / "07_hierarchy_consistency.png", dpi=150)
    plt.close(fig)


def make_all_figures(cfg):
    out_dir, fig_dir = cfg.relative_abundance_out_dir, cfg.relative_abundance_fig_dir
    wells = pd.read_csv(out_dir / "r02_well_interaction_scores.csv")
    wells = wells[~wells["missing_reference"]]
    pairs = pd.read_csv(out_dir / "r03_pair_replicate_stats.csv")
    comp = pd.read_csv(out_dir / "r04_strain_competitiveness.csv")
    bt = pd.read_csv(out_dir / "r05_bt_strengths.csv")
    mat = pd.read_csv(out_dir / "r05_pairwise_relative_abundance_matrix.csv", index_col=0)
    summary = pd.read_csv(out_dir / "r05_hierarchy_summary.csv")

    _fig00_read_classification_quality(wells, fig_dir)
    _fig01_abundance_distribution(wells, fig_dir)
    _fig02_replicate_stability(pairs, fig_dir)
    _fig03_competitiveness_ranking(comp, fig_dir)
    _fig04_mean_vs_dispersion(comp, bt, fig_dir)
    _fig05_hierarchy_heatmap(mat, bt, fig_dir)
    _fig06_bt_vs_naive(bt, comp, fig_dir)
    _fig07_hierarchy_consistency(summary, fig_dir)

    print(f"saved 8 figures -> {fig_dir}")
