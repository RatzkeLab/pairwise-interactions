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


def _pool_reads_and_fit_bt(wells):
    """wells: r02-style dataframe already filtered to resolvable (non-missing-reference) rows.
    Pools confident read counts per unique pair (canonical a/b orientation), restricts to the
    largest connected component, and fits Bradley-Terry. Shared by bt_goodness_of_fit,
    parametric_bootstrap_intransitivity, bootstrap_confidence_intervals, and
    sensitivity_reliable_pairs_only -- every test that needs "pool reads -> fit BT" starts here."""
    wells = wells.copy()
    wells["n_reads_a"] = np.where(wells["strain1"] == wells["strain_a"], wells["n_strain1"], wells["n_strain2"])
    wells["n_reads_b"] = np.where(wells["strain1"] == wells["strain_a"], wells["n_strain2"], wells["n_strain1"])
    pooled = wells.groupby(["strain_a", "strain_b"]).agg(n_reads_a=("n_reads_a", "sum"), n_reads_b=("n_reads_b", "sum")).reset_index()
    pooled = pooled[(pooled["n_reads_a"] + pooled["n_reads_b"]) > 0]

    strains = sorted(set(pooled["strain_a"]) | set(pooled["strain_b"]))
    components = union_find_components(strains, list(zip(pooled["strain_a"], pooled["strain_b"])))
    main_component = set(components[0])
    pooled_main = pooled[pooled["strain_a"].isin(main_component) & pooled["strain_b"].isin(main_component)]

    strength, se, pseudo_r2, fit = fit_bradley_terry(pooled_main, main_component)
    return strength, se, pseudo_r2, fit, main_component, pooled_main


def _win_matrix_from_pooled(pooled_main, main_component):
    """Build a discrete win/loss/tie numpy matrix directly from pooled read counts (used
    where we don't have (or don't want to reuse) the saved r05 relative-abundance matrix --
    e.g. inside a bootstrap/simulation loop operating on a resampled or simulated dataset)."""
    strains_sorted = sorted(main_component)
    idx = {s: i for i, s in enumerate(strains_sorted)}
    n = len(strains_sorted)
    A = np.full((n, n), np.nan)
    for r in pooled_main.itertuples():
        i, j = idx[r.strain_a], idx[r.strain_b]
        rel_a = r.n_reads_a / (r.n_reads_a + r.n_reads_b)
        w = 1.0 if rel_a > 0.5 else (0.0 if rel_a < 0.5 else 0.5)
        A[i, j] = w
        A[j, i] = (1 - w) if w != 0.5 else 0.5
    return strains_sorted, A


def bt_goodness_of_fit(cfg):
    wells = pd.read_csv(cfg.relative_abundance_out_dir / "r02_well_interaction_scores.csv")
    wells = wells[~wells["missing_reference"]].copy()
    _, _, pseudo_r2, fit, main_component, pooled_main = _pool_reads_and_fit_bt(wells)

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
# 3. parametric bootstrap: is the observed intransitivity more than a perfect
#    hierarchy would produce from sampling noise alone?
# ===========================================================================


def parametric_bootstrap_intransitivity(cfg, n_simulations=2000, seed=0):
    """The permutation test (test 1) asks "is there a hierarchy at all, vs. pure chance."
    This asks a sharper question: *assuming the fitted Bradley-Terry model is exactly true*,
    would its own sampling noise (at the actual read depth each pair happened to get)
    produce this much apparent intransitivity on its own? If the observed cyclic-triad
    fraction sits comfortably inside this noise band, the small amount of intransitivity we
    see could just be measurement noise around a perfect hierarchy. If it's well above the
    band, that's evidence of real, extra structure (genuine context-dependence) beyond
    binomial noise -- a more targeted companion to test 2's aggregate deviance GoF test.

    Method: take the fitted per-strain strengths and, for every tested dyad, simulate a new
    win count from Binomial(n, p) where n is that dyad's *actual* pooled read total and p is
    the BT-predicted win probability -- same noise level as the real data, but zero
    unexplained structure by construction. Rebuild the win matrix, recompute the cyclic
    fraction and DCI, repeat thousands of times. Uses the exact same dyad/triad mask as the
    observed analysis (from r05's saved matrix) so the two are directly comparable.
    """
    wells = pd.read_csv(cfg.relative_abundance_out_dir / "r02_well_interaction_scores.csv")
    wells = wells[~wells["missing_reference"]].copy()
    strength, _, _, _, main_component, pooled_main = _pool_reads_and_fit_bt(wells)

    mat = pd.read_csv(cfg.relative_abundance_out_dir / "r05_pairwise_relative_abundance_matrix.csv", index_col=0)
    observed = validate_fast_implementations(mat)
    win_obs = win_matrix_from_relative_abundance(mat)
    strains, A_obs = _win_to_array(win_obs)
    mask = _valid_mask(A_obs)
    n = len(strains)
    idx_of = {s: i for i, s in enumerate(strains)}

    upper = np.triu(np.ones((n, n), dtype=bool), k=1)
    valid_pair = mask & upper
    dyad_i, dyad_j = np.where(valid_pair)

    pooled_lookup = {
        (r.strain_a, r.strain_b): (r.n_reads_a + r.n_reads_b, strength[r.strain_a] - strength[r.strain_b])
        for r in pooled_main.itertuples()
    }
    n_totals = np.empty(len(dyad_i))
    p_win = np.empty(len(dyad_i))
    for k, (di, dj) in enumerate(zip(dyad_i, dyad_j)):
        a, b = strains[di], strains[dj]
        ntot, delta = pooled_lookup[(a, b)]
        n_totals[k] = ntot
        p_win[k] = 1 / (1 + np.exp(-delta))
    n_totals_int = n_totals.astype(int)

    I, J, K, valid_triple = _precompute_triad_indices(n)
    mask_ijk = valid_triple & mask[I, J] & mask[J, K] & mask[I, K]

    rng = np.random.default_rng(seed)
    null_cyclic = np.empty(n_simulations)
    null_dci = np.empty(n_simulations)

    for s in range(n_simulations):
        wins_i = rng.binomial(n_totals_int, p_win)
        frac_i = wins_i / n_totals
        w = np.where(frac_i > 0.5, 1.0, np.where(frac_i < 0.5, 0.0, 0.5))
        A = np.full((n, n), np.nan)
        A[dyad_i, dyad_j] = w
        A[dyad_j, dyad_i] = np.where(w == 0.5, 0.5, 1 - w)

        null_dci[s], _, _ = fast_dci(A, mask)
        null_cyclic[s], _, _ = _cyclic_fraction_fast(A, I, J, K, mask_ijk)

    # one-sided: is the observed cyclic fraction unusually HIGH relative to pure BT + noise?
    p_excess_cyclic = (np.sum(null_cyclic >= observed["cyclic_frac"]) + 1) / (n_simulations + 1)

    summary = pd.DataFrame([
        {"metric": "cyclic_frac_observed", "value": observed["cyclic_frac"]},
        {"metric": "cyclic_frac_bt_noise_mean", "value": float(np.mean(null_cyclic))},
        {"metric": "cyclic_frac_bt_noise_sd", "value": float(np.std(null_cyclic))},
        {"metric": "p_value_excess_intransitivity", "value": p_excess_cyclic},
        {"metric": "dci_observed", "value": observed["dci"]},
        {"metric": "dci_bt_noise_mean", "value": float(np.mean(null_dci))},
        {"metric": "dci_bt_noise_sd", "value": float(np.std(null_dci))},
        {"metric": "n_simulations", "value": n_simulations},
        {"metric": "n_dyads_simulated", "value": len(dyad_i)},
    ])
    out_path = cfg.relative_abundance_out_dir / "r09_parametric_bootstrap_intransitivity.csv"
    summary.to_csv(out_path, index=False)

    null_path = cfg.relative_abundance_out_dir / "r09_parametric_bootstrap_null_distributions.csv"
    pd.DataFrame({"cyclic_frac": null_cyclic, "dci": null_dci}).to_csv(null_path, index=False)

    print(f"observed cyclic fraction = {observed['cyclic_frac']:.3f}  |  "
          f"BT-model + sampling-noise-alone cyclic fraction = {np.mean(null_cyclic):.3f} +/- {np.std(null_cyclic):.3f}  |  "
          f"p = {p_excess_cyclic:.4g}")
    if p_excess_cyclic < 0.05:
        print("  -> significant: there is MORE intransitivity than a perfect hierarchy's own sampling "
              "noise would produce -- real excess non-hierarchical structure, not just noise")
    else:
        print("  -> not significant: the observed intransitivity is consistent with pure sampling "
              "noise around a perfect hierarchy")
    print(f"(for comparison: observed DCI = {observed['dci']:.3f} vs. BT+noise DCI = "
          f"{np.mean(null_dci):.3f} +/- {np.std(null_dci):.3f} -- a calibration check, not a formal test)")
    print(f"saved -> {out_path}, {null_path}")
    return summary, null_cyclic, null_dci


# ===========================================================================
# 4. bootstrap confidence intervals on the headline hierarchy statistics
# ===========================================================================


def bootstrap_confidence_intervals(cfg, n_bootstrap=500, seed=0):
    """Resample wells with replacement (the natural unit of the actual data-generating
    process -- reads are nested in wells, wells are nested in pairs) and recompute
    pseudo-R^2, DCI, and cyclic fraction each time, to attach a 95% CI to all three headline
    descriptive statistics rather than reporting them as bare point estimates."""
    wells = pd.read_csv(cfg.relative_abundance_out_dir / "r02_well_interaction_scores.csv")
    wells = wells[~wells["missing_reference"]].reset_index(drop=True)
    n_wells = len(wells)

    rng = np.random.default_rng(seed)
    boot_pseudo_r2 = np.full(n_bootstrap, np.nan)
    boot_dci = np.full(n_bootstrap, np.nan)
    boot_cyclic = np.full(n_bootstrap, np.nan)

    for b in range(n_bootstrap):
        idx = rng.integers(0, n_wells, size=n_wells)
        sample = wells.iloc[idx]
        try:
            _, _, pseudo_r2, _, main_component, pooled_main = _pool_reads_and_fit_bt(sample)
        except Exception:
            continue
        boot_pseudo_r2[b] = pseudo_r2

        _, A = _win_matrix_from_pooled(pooled_main, main_component)
        mask = _valid_mask(A)
        boot_dci[b], _, _ = fast_dci(A, mask)
        boot_cyclic[b], _, _ = fast_cyclic_fraction(A, mask)

    def _ci(x):
        x = x[~np.isnan(x)]
        return float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5)), int(len(x))

    observed = pd.read_csv(cfg.relative_abundance_out_dir / "r05_hierarchy_summary.csv").set_index("metric")["value"]

    rows = []
    for name, boot_arr, obs_val in [
        ("pseudo_r2", boot_pseudo_r2, observed["bt_pseudo_r2"]),
        ("dci", boot_dci, observed["dci"]),
        ("cyclic_frac", boot_cyclic, observed["frac_intransitive_triads"]),
    ]:
        lo, hi, n_valid = _ci(boot_arr)
        rows.append({
            "metric": name, "observed": obs_val, "ci_lower_95": lo, "ci_upper_95": hi,
            "n_valid_bootstraps": n_valid, "n_bootstrap": n_bootstrap,
        })
    summary = pd.DataFrame(rows)
    out_path = cfg.relative_abundance_out_dir / "r10_bootstrap_confidence_intervals.csv"
    summary.to_csv(out_path, index=False)

    boot_path = cfg.relative_abundance_out_dir / "r10_bootstrap_distributions.csv"
    pd.DataFrame({"pseudo_r2": boot_pseudo_r2, "dci": boot_dci, "cyclic_frac": boot_cyclic}).to_csv(boot_path, index=False)

    for _, row in summary.iterrows():
        print(f"{row['metric']:14s} observed={row['observed']:.3f}  95% CI=[{row['ci_lower_95']:.3f}, {row['ci_upper_95']:.3f}]  "
              f"({row['n_valid_bootstraps']}/{n_bootstrap} bootstraps valid)")
    print(f"saved -> {out_path}, {boot_path}")
    return summary


# ===========================================================================
# 5. sensitivity check: how much does restricting to reliable (not hard-to-call)
#    pairs change the hierarchy statistics?
# ===========================================================================


def sensitivity_reliable_pairs_only(cfg):
    """Most tested pairs are singleton-replicate, and a meaningful minority are
    near-identical-reference pairs whose "ties" are a measurement ceiling, not biology (see
    relative_abundance.py's docstring). Recomputes pseudo-R^2/DCI/cyclic-fraction restricted
    to pairs r03 did NOT flag `high_uncertainty_pair`, to see how much of the current
    "15% upset rate"-style tangledness is attributable to those hard-to-call pairs rather
    than genuine biological intransitivity.

    Note: the "all_pairs" column here is recomputed from pooled per-pair read counts (ratio
    of sums, matching how the BT model itself is fit), not read from r05 directly -- r05's
    win matrix instead uses the mean of each replicate well's own ratio (see
    relative_abundance.hierarchy_analysis). The two conventions agree for singleton-replicate
    pairs (most of them) and differ only slightly (a read-count-weighted vs. an
    equally-weighted-per-well average) for the ~445/2159 pairs with multiple replicates --
    e.g. DCI 0.701 here vs. 0.696 in r05 for 20260721. Both are defensible; what matters for
    this test is that "all_pairs" and "reliable_pairs_only" use the identical convention, so
    the comparison between them is fair even though "all_pairs" won't match r05 to the third
    decimal.
    """
    wells = pd.read_csv(cfg.relative_abundance_out_dir / "r02_well_interaction_scores.csv")
    wells = wells[~wells["missing_reference"]].copy()
    pair_stats = pd.read_csv(cfg.relative_abundance_out_dir / "r03_pair_replicate_stats.csv")

    reliable_keys = set(zip(
        pair_stats.loc[~pair_stats["high_uncertainty_pair"], "strain_a"],
        pair_stats.loc[~pair_stats["high_uncertainty_pair"], "strain_b"],
    ))
    wells["pair_key"] = list(zip(wells["strain_a"], wells["strain_b"]))
    wells_reliable = wells[wells["pair_key"].isin(reliable_keys)]

    def _stats_for(w):
        _, _, pseudo_r2, _, main_component, pooled_main = _pool_reads_and_fit_bt(w)
        _, A = _win_matrix_from_pooled(pooled_main, main_component)
        mask = _valid_mask(A)
        dci, n_cons, n_incons = fast_dci(A, mask)
        cyc, n_cyc, n_tri = fast_cyclic_fraction(A, mask)
        return {
            "n_strains": len(main_component), "n_pairs": len(pooled_main), "pseudo_r2": pseudo_r2,
            "dci": dci, "n_consistent": n_cons, "n_inconsistent": n_incons,
            "cyclic_frac": cyc, "n_cyclic": n_cyc, "n_triads": n_tri,
        }

    all_stats = _stats_for(wells)
    reliable_stats = _stats_for(wells_reliable)

    summary = pd.DataFrame(
        [{"metric": k, "all_pairs": all_stats[k], "reliable_pairs_only": reliable_stats[k]} for k in all_stats]
    )
    out_path = cfg.relative_abundance_out_dir / "r11_sensitivity_reliable_pairs.csv"
    summary.to_csv(out_path, index=False)

    print(f"all pairs      : n_strains={all_stats['n_strains']}, n_pairs={all_stats['n_pairs']}, "
          f"pseudo_r2={all_stats['pseudo_r2']:.3f}, DCI={all_stats['dci']:.3f}, cyclic_frac={all_stats['cyclic_frac']:.3f}")
    print(f"reliable only  : n_strains={reliable_stats['n_strains']}, n_pairs={reliable_stats['n_pairs']}, "
          f"pseudo_r2={reliable_stats['pseudo_r2']:.3f}, DCI={reliable_stats['dci']:.3f}, cyclic_frac={reliable_stats['cyclic_frac']:.3f}")
    print(f"saved -> {out_path}")
    return summary


# ===========================================================================
# figures
# ===========================================================================


def _fig08_permutation_test(cfg):
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


def _fig09_parametric_bootstrap(cfg):
    summary = pd.read_csv(cfg.relative_abundance_out_dir / "r09_parametric_bootstrap_intransitivity.csv").set_index("metric")["value"]
    null_dist = pd.read_csv(cfg.relative_abundance_out_dir / "r09_parametric_bootstrap_null_distributions.csv")

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.hist(null_dist["cyclic_frac"], bins=40, color=COLOR_BLUE, edgecolor="white", linewidth=0.3,
            label="BT model + sampling noise only")
    ax.axvline(summary["cyclic_frac_observed"], color=COLOR_CRITICAL, lw=2, label=f"observed = {summary['cyclic_frac_observed']:.3f}")
    ax.set_xlabel("intransitive (cyclic) triad fraction")
    ax.set_ylabel("# simulations")
    ax.set_title(f"Observed intransitivity vs. a perfect hierarchy's own sampling noise\n"
                 f"(p = {summary['p_value_excess_intransitivity']:.4g})")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(cfg.relative_abundance_fig_dir / "09_parametric_bootstrap_intransitivity.png", dpi=150)
    plt.close(fig)


def _fig10_bootstrap_ci(cfg):
    summary = pd.read_csv(cfg.relative_abundance_out_dir / "r10_bootstrap_confidence_intervals.csv").set_index("metric")
    dist = pd.read_csv(cfg.relative_abundance_out_dir / "r10_bootstrap_distributions.csv")

    panels = [("pseudo_r2", "Bradley-Terry pseudo-R²"), ("dci", "Directional Consistency Index"), ("cyclic_frac", "cyclic-triad fraction")]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, (key, label) in zip(axes, panels):
        vals = dist[key].dropna()
        row = summary.loc[key]
        ax.hist(vals, bins=30, color=COLOR_BLUE, edgecolor="white", linewidth=0.3)
        ax.axvline(row["observed"], color=COLOR_CRITICAL, lw=2, label=f"observed = {row['observed']:.3f}")
        ax.axvspan(row["ci_lower_95"], row["ci_upper_95"], color=COLOR_CRITICAL, alpha=0.12,
                   label=f"95% CI [{row['ci_lower_95']:.3f}, {row['ci_upper_95']:.3f}]")
        ax.set_xlabel(label)
        ax.set_ylabel("# bootstrap resamples")
        ax.legend(frameon=False, fontsize=7.5)
    fig.suptitle("Bootstrap (resample wells with replacement) confidence intervals", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(cfg.relative_abundance_fig_dir / "10_bootstrap_confidence_intervals.png", dpi=150)
    plt.close(fig)


def _fig11_sensitivity_comparison(cfg):
    summary = pd.read_csv(cfg.relative_abundance_out_dir / "r11_sensitivity_reliable_pairs.csv").set_index("metric")

    panels = ["pseudo_r2", "dci", "cyclic_frac"]
    labels = ["Bradley-Terry\npseudo-R²", "Directional\nConsistency Index", "cyclic-triad\nfraction"]
    all_vals = [summary.loc[p, "all_pairs"] for p in panels]
    rel_vals = [summary.loc[p, "reliable_pairs_only"] for p in panels]

    x = np.arange(len(panels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    ax.bar(x - width / 2, all_vals, width, color=COLOR_BLUE, label="all tested pairs")
    ax.bar(x + width / 2, rel_vals, width, color=COLOR_CRITICAL, label="excluding hard-to-call pairs")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title("Sensitivity to excluding hard-to-call pairs", pad=10)
    fig.text(0.5, 0.925, "near-identical-reference pairs, whose \"ties\" are a measurement ceiling, not biology",
              ha="center", fontsize=8.5, color=COLOR_TEXT_SECONDARY)
    for i, (a, r) in enumerate(zip(all_vals, rel_vals)):
        ax.text(i - width / 2, a, f"{a:.3f}", ha="center", va="bottom", fontsize=8)
        ax.text(i + width / 2, r, f"{r:.3f}", ha="center", va="bottom", fontsize=8)
    ax.legend(frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(cfg.relative_abundance_fig_dir / "11_sensitivity_reliable_pairs.png", dpi=150)
    plt.close(fig)


def make_significance_figures(cfg):
    _fig08_permutation_test(cfg)
    _fig09_parametric_bootstrap(cfg)
    _fig10_bootstrap_ci(cfg)
    _fig11_sensitivity_comparison(cfg)
    print(f"saved 4 figures -> {cfg.relative_abundance_fig_dir}")
