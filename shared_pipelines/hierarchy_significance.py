"""Statistical significance testing for the hierarchy analysis in relative_abundance.py.

relative_abundance.hierarchy_analysis() reports three *descriptive* statistics (Bradley-Terry
pseudo-R^2, Directional Consistency Index, intransitive-triad fraction) with no p-values or
confidence intervals attached -- they say "this looks hierarchical" but not "this is
statistically distinguishable from a non-hierarchical null." This module adds two rigorous
tests for that:

1. permutation_test() -- a randomization test against the null "these strains have no real
   consistent ranking at all." For every tested strain pair with a clear (non-tie) winner,
   the winner is redrawn by a fair coin flip, DCI and the cyclic-triad fraction are
   recomputed, and this is repeated thousands of times to build a null distribution. The
   observed DCI/cyclic-fraction are compared against that null: how extreme are they if there
   were truly no hierarchy? This is the standard randomization-test approach used in the
   behavioral-ecology dominance-hierarchy literature (e.g. de Vries 1998), adapted here for
   microbial pairwise competition data.

   Because this runs thousands of iterations, the win/loss matrix arithmetic is done with
   vectorized numpy (fast_dci / fast_cyclic_fraction) rather than the itertools-based
   reference implementation in relative_abundance.py (compute_directional_consistency /
   count_intransitive_triads), which is fine for a single one-off computation but far too
   slow to call thousands of times. fast_dci/fast_cyclic_fraction are validated in
   validate_fast_implementations() to reproduce the reference implementation exactly on the
   observed data before being trusted inside the permutation loop.

2. bt_goodness_of_fit() -- two likelihood-ratio chi-squared tests, both essentially free
   given the already-fit Bradley-Terry GLM:
     a. "do strains differ at all?" -- null_deviance vs. deviance, df = n_strains - 1.
        Should be overwhelmingly significant given the pseudo-R^2 already reported; this
        just makes that formal.
     b. "is a single 1-D strength scale a SUFFICIENT description?" -- the model's residual
        deviance against a chi-squared(df_resid) reference distribution (this is exactly a
        deviance goodness-of-fit test for a grouped-binomial GLM: deviance is defined
        relative to the saturated one-parameter-per-pair model, so this one test already
        answers "does the saturated model fit significantly better than Bradley-Terry" without
        needing to fit it separately). A significant result here is direct, quantified
        evidence of real structure (intransitivity / opponent-specific context-dependence)
        beyond what any single hierarchy ranking can capture.
     Caveat: the chi-squared approximation for GLM deviance goodness-of-fit assumes cell
     counts aren't too small; most pairs here pool tens-to-hundreds of reads, but pairs with
     very few reads make the approximation less reliable -- reported as-is, with that caveat.

Output (per experiment, in cfg.relative_abundance_out_dir):
  r07_permutation_test.csv               observed/null summary stats + p-values
  r07_permutation_null_distributions.csv raw null draws (for the histogram figure)
  r08_bt_goodness_of_fit.csv             both LR test results
"""

import numpy as np
import pandas as pd
from scipy.stats import chi2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from relative_abundance import (
    union_find_components, fit_bradley_terry,
    compute_directional_consistency, count_intransitive_triads,
)

N_PERMUTATIONS_DEFAULT = 5000

COLOR_BLUE = "#2a78d6"
COLOR_CRITICAL = "#d03b3b"
COLOR_TEXT_SECONDARY = "#52514e"


# ===========================================================================
# win-matrix reconstruction (from the already-saved relative-abundance matrix)
# ===========================================================================


def win_matrix_from_relative_abundance(mat):
    """Reconstruct the discrete win/loss/tie matrix from the continuous relative-abundance
    matrix, exactly as relative_abundance.hierarchy_analysis assigns it internally
    (w = 1 if rel_ab > 0.5, 0 if < 0.5, 0.5 if exactly tied)."""
    values = mat.to_numpy(dtype=float)
    win_values = np.where(np.isnan(values), np.nan, np.where(values > 0.5, 1.0, np.where(values < 0.5, 0.0, 0.5)))
    return pd.DataFrame(win_values, index=mat.index, columns=mat.columns)


def _win_to_array(win):
    strains = list(win.index)
    return strains, win.reindex(index=strains, columns=strains).to_numpy(dtype=float)


def _valid_mask(A):
    return ~np.isnan(A) & (A != 0.5)


# ===========================================================================
# fast (numpy-vectorized) DCI / cyclic-triad-fraction, for use inside the
# permutation loop -- validated against the reference implementation below
# ===========================================================================


def fast_dci(A, mask):
    n = A.shape[0]
    wins = mask & (A == 1)
    total_wins = wins.sum(axis=1)
    order = np.argsort(-total_wins, kind="stable")
    rank = np.empty(n, dtype=int)
    rank[order] = np.arange(n)

    upper = np.triu(np.ones((n, n), dtype=bool), k=1)
    valid_pair = mask & upper
    winner_is_row = A == 1
    ri = np.broadcast_to(rank[:, None], (n, n))
    rj = np.broadcast_to(rank[None, :], (n, n))
    consistent = (winner_is_row & (ri < rj)) | (~winner_is_row & (rj < ri))

    n_consistent = int(np.sum(consistent & valid_pair))
    n_total = int(np.sum(valid_pair))
    n_inconsistent = n_total - n_consistent
    dci = (n_consistent - n_inconsistent) / n_total if n_total else np.nan
    return dci, n_consistent, n_inconsistent


def _precompute_triad_indices(n):
    idx = np.arange(n)
    I, J, K = np.meshgrid(idx, idx, idx, indexing="ij")
    valid_triple = (I < J) & (J < K)
    return I, J, K, valid_triple


def _cyclic_fraction_fast(A, I, J, K, mask_ijk):
    a_beats_b = A[I, J] == 1
    b_beats_c = A[J, K] == 1
    a_beats_c = A[I, K] == 1
    cyclic = (a_beats_b & b_beats_c & ~a_beats_c) | (~a_beats_b & ~b_beats_c & a_beats_c)
    n_cyclic = int(np.sum(cyclic & mask_ijk))
    n_evaluated = int(np.sum(mask_ijk))
    frac = n_cyclic / n_evaluated if n_evaluated else np.nan
    return frac, n_cyclic, n_evaluated


def fast_cyclic_fraction(A, mask):
    n = A.shape[0]
    I, J, K, valid_triple = _precompute_triad_indices(n)
    mask_ijk = valid_triple & mask[I, J] & mask[J, K] & mask[I, K]
    return _cyclic_fraction_fast(A, I, J, K, mask_ijk)


def validate_fast_implementations(mat):
    """Sanity check: fast_dci/fast_cyclic_fraction must reproduce the reference
    (itertools/pandas-based) implementation exactly on the observed data before being
    trusted inside a 1000s-of-iterations permutation loop."""
    win = win_matrix_from_relative_abundance(mat)
    strains, A = _win_to_array(win)
    mask = _valid_mask(A)

    dci_fast, nc_fast, ni_fast = fast_dci(A, mask)
    cyc_fast, ncyc_fast, ntri_fast = fast_cyclic_fraction(A, mask)

    total_wins = (win == 1).sum(axis=1)
    rank_order = total_wins.sort_values(ascending=False, kind="stable").index.tolist()
    dci_ref, nc_ref, ni_ref = compute_directional_consistency(win, rank_order)
    cyc_ref, ncyc_ref, ntri_ref = count_intransitive_triads(win)

    assert (nc_fast, ni_fast) == (nc_ref, ni_ref), f"DCI mismatch: fast={nc_fast, ni_fast} ref={nc_ref, ni_ref}"
    assert abs(dci_fast - dci_ref) < 1e-9, f"DCI value mismatch: {dci_fast} vs {dci_ref}"
    assert (ncyc_fast, ntri_fast) == (ncyc_ref, ntri_ref), f"cyclic mismatch: fast={ncyc_fast, ntri_fast} ref={ncyc_ref, ntri_ref}"
    assert abs(cyc_fast - cyc_ref) < 1e-9, f"cyclic fraction mismatch: {cyc_fast} vs {cyc_ref}"
    return {"dci": dci_ref, "n_consistent": nc_ref, "n_inconsistent": ni_ref,
            "cyclic_frac": cyc_ref, "n_cyclic": ncyc_ref, "n_triads": ntri_ref}


# ===========================================================================
# 1. permutation test
# ===========================================================================


def permutation_test(cfg, n_permutations=N_PERMUTATIONS_DEFAULT, seed=0):
    mat = pd.read_csv(cfg.relative_abundance_out_dir / "r05_pairwise_relative_abundance_matrix.csv", index_col=0)
    observed = validate_fast_implementations(mat)

    win = win_matrix_from_relative_abundance(mat)
    strains, A_obs = _win_to_array(win)
    mask = _valid_mask(A_obs)
    n = len(strains)

    upper_mask = mask & np.triu(np.ones((n, n), dtype=bool), k=1)
    dyad_i, dyad_j = np.where(upper_mask)

    I, J, K, valid_triple = _precompute_triad_indices(n)
    mask_ijk = valid_triple & mask[I, J] & mask[J, K] & mask[I, K]

    rng = np.random.default_rng(seed)
    null_dci = np.empty(n_permutations)
    null_cyclic = np.empty(n_permutations)

    for p in range(n_permutations):
        draws = (rng.random(len(dyad_i)) < 0.5).astype(float)
        A = np.full((n, n), np.nan)
        A[dyad_i, dyad_j] = draws
        A[dyad_j, dyad_i] = 1 - draws

        null_dci[p], _, _ = fast_dci(A, mask)
        null_cyclic[p], _, _ = _cyclic_fraction_fast(A, I, J, K, mask_ijk)

    # one-sided: hierarchy predicts DCI unusually HIGH, cyclic fraction unusually LOW
    p_dci = (np.sum(null_dci >= observed["dci"]) + 1) / (n_permutations + 1)
    p_cyclic = (np.sum(null_cyclic <= observed["cyclic_frac"]) + 1) / (n_permutations + 1)

    summary = pd.DataFrame([
        {"metric": "dci_observed", "value": observed["dci"]},
        {"metric": "dci_null_mean", "value": float(np.mean(null_dci))},
        {"metric": "dci_null_sd", "value": float(np.std(null_dci))},
        {"metric": "dci_p_value", "value": p_dci},
        {"metric": "cyclic_frac_observed", "value": observed["cyclic_frac"]},
        {"metric": "cyclic_frac_null_mean", "value": float(np.mean(null_cyclic))},
        {"metric": "cyclic_frac_null_sd", "value": float(np.std(null_cyclic))},
        {"metric": "cyclic_frac_p_value", "value": p_cyclic},
        {"metric": "n_permutations", "value": n_permutations},
        {"metric": "n_dyads_tested", "value": len(dyad_i)},
        {"metric": "n_triads_evaluated", "value": observed["n_triads"]},
    ])
    out_path = cfg.relative_abundance_out_dir / "r07_permutation_test.csv"
    summary.to_csv(out_path, index=False)

    null_path = cfg.relative_abundance_out_dir / "r07_permutation_null_distributions.csv"
    pd.DataFrame({"dci": null_dci, "cyclic_frac": null_cyclic}).to_csv(null_path, index=False)

    print(f"validated fast DCI/cyclic-fraction implementations against the reference (exact match)")
    print(f"observed DCI = {observed['dci']:.3f}  |  null (no-hierarchy) DCI = {np.mean(null_dci):.3f} +/- {np.std(null_dci):.3f}  |  p = {p_dci:.4g}")
    print(f"observed cyclic fraction = {observed['cyclic_frac']:.3f}  |  null cyclic fraction = {np.mean(null_cyclic):.3f} +/- {np.std(null_cyclic):.3f}  |  p = {p_cyclic:.4g}")
    print(f"saved -> {out_path}, {null_path}")
    return summary, null_dci, null_cyclic


# ===========================================================================
# 2. Bradley-Terry goodness-of-fit (likelihood-ratio chi-squared tests)
# ===========================================================================


def bt_goodness_of_fit(cfg):
    wells = pd.read_csv(cfg.relative_abundance_out_dir / "r02_well_interaction_scores.csv")
    wells = wells[~wells["missing_reference"]].copy()
    wells["n_reads_a"] = np.where(wells["strain1"] == wells["strain_a"], wells["n_strain1"], wells["n_strain2"])
    wells["n_reads_b"] = np.where(wells["strain1"] == wells["strain_a"], wells["n_strain2"], wells["n_strain1"])
    pooled = wells.groupby(["strain_a", "strain_b"]).agg(n_reads_a=("n_reads_a", "sum"), n_reads_b=("n_reads_b", "sum")).reset_index()
    pooled = pooled[(pooled["n_reads_a"] + pooled["n_reads_b"]) > 0]

    strains = sorted(set(pooled["strain_a"]) | set(pooled["strain_b"]))
    components = union_find_components(strains, list(zip(pooled["strain_a"], pooled["strain_b"])))
    main_component = set(components[0])
    pooled_main = pooled[pooled["strain_a"].isin(main_component) & pooled["strain_b"].isin(main_component)]

    _, _, pseudo_r2, fit = fit_bradley_terry(pooled_main, main_component)

    n_params = len(main_component) - 1
    lr_stat = fit.null_deviance - fit.deviance
    p_strains_differ = chi2.sf(lr_stat, n_params)

    p_gof = chi2.sf(fit.deviance, fit.df_resid)
    # With thousands of reads, the GoF p-value alone is close to uninformative: it will
    # reject a merely-imperfect-but-good model just as readily as a truly bad one, since
    # power grows with N. The deviance/df ratio is the companion effect-size number -- how
    # many times more residual variation is there than a perfectly-fitting BT model would
    # produce from binomial sampling noise alone (1.0 = perfect fit; bigger = more excess,
    # unmodeled structure, roughly independent of sample size).
    dispersion_ratio = fit.deviance / fit.df_resid

    # smallest per-pair read totals, for the chi-squared-approximation caveat
    pair_totals = pooled_main["n_reads_a"] + pooled_main["n_reads_b"]

    summary = pd.DataFrame([
        {"metric": "n_strains", "value": len(main_component)},
        {"metric": "n_pairs", "value": len(pooled_main)},
        {"metric": "n_params", "value": n_params},
        {"metric": "pseudo_r2", "value": pseudo_r2},
        {"metric": "null_deviance", "value": fit.null_deviance},
        {"metric": "residual_deviance", "value": fit.deviance},
        {"metric": "df_resid", "value": fit.df_resid},
        {"metric": "lr_stat_strains_differ", "value": lr_stat},
        {"metric": "df_strains_differ", "value": n_params},
        {"metric": "p_value_strains_differ", "value": p_strains_differ},
        {"metric": "gof_chi2_stat", "value": fit.deviance},
        {"metric": "df_gof", "value": fit.df_resid},
        {"metric": "p_value_gof_1d_hierarchy_sufficient", "value": p_gof},
        {"metric": "dispersion_ratio", "value": dispersion_ratio},
        {"metric": "min_pair_read_total", "value": int(pair_totals.min())},
        {"metric": "frac_pairs_read_total_lt_10", "value": float((pair_totals < 10).mean())},
    ])
    out_path = cfg.relative_abundance_out_dir / "r08_bt_goodness_of_fit.csv"
    summary.to_csv(out_path, index=False)

    print(f"LR test (do strains differ at all?): stat={lr_stat:.1f}, df={n_params}, p={p_strains_differ:.4g}")
    print(f"deviance GoF test (is a single 1-D hierarchy a SUFFICIENT description?): "
          f"deviance={fit.deviance:.1f}, df_resid={fit.df_resid}, p={p_gof:.4g}, "
          f"dispersion ratio={dispersion_ratio:.2f}x")
    print(f"  (with this many reads, the p-value alone will reject any imperfect model; the "
          f"dispersion ratio is the sample-size-independent effect size -- {dispersion_ratio:.1f}x "
          f"more residual variation than a perfect-fitting hierarchy would show from sampling noise alone)")
    if p_gof < 0.05:
        print("  -> significant: there is real structure beyond a single hierarchy ranking "
              "(consistent with the nonzero upset rate / cyclic-triad fraction already observed)")
    else:
        print("  -> not significant: a single hierarchy ranking adequately explains the data "
              "within sampling noise")
    print(f"saved -> {out_path}")
    return summary


# ===========================================================================
# figures
# ===========================================================================


def make_significance_figures(cfg):
    summary = pd.read_csv(cfg.relative_abundance_out_dir / "r07_permutation_test.csv").set_index("metric")["value"]
    null_dist = pd.read_csv(cfg.relative_abundance_out_dir / "r07_permutation_null_distributions.csv")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    axes[0].hist(null_dist["dci"], bins=40, color=COLOR_BLUE, edgecolor="white", linewidth=0.3, label="null (no hierarchy)")
    axes[0].axvline(summary["dci_observed"], color=COLOR_CRITICAL, lw=2, label=f"observed = {summary['dci_observed']:.2f}")
    axes[0].set_xlabel("Directional Consistency Index")
    axes[0].set_ylabel("# permutations")
    axes[0].set_title(f"DCI vs. random-outcome null (p = {summary['dci_p_value']:.4g})")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].hist(null_dist["cyclic_frac"], bins=40, color=COLOR_BLUE, edgecolor="white", linewidth=0.3, label="null (no hierarchy)")
    axes[1].axvline(summary["cyclic_frac_observed"], color=COLOR_CRITICAL, lw=2, label=f"observed = {summary['cyclic_frac_observed']:.3f}")
    axes[1].set_xlabel("intransitive (cyclic) triad fraction")
    axes[1].set_ylabel("# permutations")
    axes[1].set_title(f"Cyclic fraction vs. random-outcome null (p = {summary['cyclic_frac_p_value']:.4g})")
    axes[1].legend(frameon=False, fontsize=8)

    fig.suptitle(f"Observed hierarchy statistics vs. a {int(summary['n_permutations'])}-permutation no-hierarchy null", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(cfg.relative_abundance_fig_dir / "08_permutation_test.png", dpi=150)
    plt.close(fig)
    print(f"saved -> {cfg.relative_abundance_fig_dir / '08_permutation_test.png'}")
