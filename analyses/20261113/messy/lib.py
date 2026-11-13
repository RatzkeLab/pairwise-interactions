"""CV driver shared by s02 (leakage) and s03 (wrong genomes).

Reuses scripts/genomic_ml.py unchanged for everything that defines a model: the dataset
filters, the feature matrix, the per-fold context (_Ctx: in-fold PCA, antisymmetric design),
the model zoo and _metrics. Only the fold assignment is new:

    group_folds(pairs, group_of)   cv_strain over GROUPS of wells. With group_of = identity it
                                   is genomic_ml.make_folds(kind="cv_strain") exactly (same
                                   sort, same permutation) -- s02 checks this by reproducing
                                   the published per-fold numbers.

and one diagnostic model, clone_lookup_no_genomics, that knows nothing about gene content:
for a test well it returns the in-fold strength of training wells carrying the SAME genome
(same >= 99.9% ANI clone cluster), else 0. Whatever it scores is, by construction, leakage.

Work is farmed out per (fold x training arm); each job returns per-pair predictions, and
metrics are computed afterwards on whichever subset of test pairs a question needs.
"""
import shutil
import warnings
from multiprocessing import get_context

import numpy as np
import pandas as pd

import config
import genomic_ml as gml

MODELS = ["zero_baseline", "clone_lookup_no_genomics", "ridge_pca", "two_stage_ridge",
          "xgboost_raw_ko", "ridge_phylo16s_only", "two_stage_ridge_phylo16s_only",
          "SHUFFLED_ridge_pca", "SHUFFLED_two_stage_ridge"]
GENOME_MODELS = ["ridge_pca", "two_stage_ridge", "xgboost_raw_ko"]
TAXONOMY_MODELS = ["ridge_phylo16s_only", "two_stage_ridge_phylo16s_only"]

_D = {}          # per-process data, filled by load() before the pool forks


def load(exp, out_dir):
    """Rebuild the published modelling set, check it is the published one, return the pieces."""
    gcfg = config.gml_config(exp, out_dir)
    pub = config.EXPERIMENTS[exp]["pub"]
    pairs, _ = gml.build_dataset(gcfg)
    ref = pd.read_csv(pub / "g01_modeling_pairs.csv")
    assert len(pairs) == len(ref) and (pairs.pair_id.values == ref.pair_id.values).all() and \
        np.allclose(pairs[gcfg.target], ref[gcfg.target]), f"{exp}: dataset differs from the published one"
    X, summ = gml.strain_feature_matrix(gcfg, pairs)
    cache = out_dir / "g04_phylo16s_distance_matrix.csv"
    if not cache.exists():
        shutil.copy(pub / "g04_phylo16s_distance_matrix.csv", cache)
    phylo = gml.phylo_distance_matrix(gcfg, pairs)
    wells = pd.read_csv(config.OUT / f"s01_wells_{exp}.csv").set_index("well")
    assert set(wells.index) == set(X.index)
    _D[exp] = dict(gcfg=gcfg, pairs=pairs, X=X, summ=summ, phylo=phylo, wells=wells)
    return _D[exp]


def group_folds(pairs, group_of, n_splits=5, seed=0):
    """(k, train_idx, test_idx): test = both strains' groups in fold k, train = neither."""
    rng = np.random.default_rng(seed)
    strains = sorted(set(pairs["strain_a"]) | set(pairs["strain_b"]))
    groups = sorted({group_of[s] for s in strains})
    fold_of_group = dict(zip(groups, rng.permutation(len(groups)) % n_splits))
    fa = pairs["strain_a"].map(group_of).map(fold_of_group).values
    fb = pairs["strain_b"].map(group_of).map(fold_of_group).values
    for k in range(n_splits):
        test = np.where((fa == k) & (fb == k))[0]
        train = np.where((fa != k) & (fb != k))[0]
        if len(test) >= 5:
            yield k, train, test


def strains_of(pairs, idx):
    p = pairs.iloc[idx]
    return set(p.strain_a) | set(p.strain_b)


def random_matched(pairs, tr, keep_ok, n_remove, rng):
    """Drop whole random wells among those allowed (keep_ok[s] True), until >= n_remove
    training pairs are gone. The "is it just less data?" control for any removal arm."""
    p = pairs.iloc[tr]
    cand = sorted(s for s in strains_of(pairs, tr) if keep_ok[s])
    rng.shuffle(cand)
    gone, removed = set(), np.zeros(len(tr), bool)
    for s in cand:
        if removed.sum() >= n_remove:
            break
        gone.add(s)
        removed = (p.strain_a.isin(gone) | p.strain_b.isin(gone)).values
    return tr[~removed]


def m_clone_lookup(ctx, genome_unit):
    s = gml._fit_additive_strength(ctx.pairs, ctx.train_idx, ctx.train_strains, target=ctx.target)
    unit_mean = s.groupby(s.index.map(genome_unit)).mean()
    te = ctx.pairs.iloc[ctx.test_idx]
    # a test well is never in training under cv_strain, so its unit's mean is its mates' mean
    est = lambda x: unit_mean.get(genome_unit[x], 0.0)
    return np.array([est(a) - est(b) for a, b in zip(te.strain_a, te.strain_b)])


def run_job(job):
    """One fold x one training arm -> per-pair predictions of every requested model."""
    d = _D[job["exp"]]
    pairs, X, summ, phylo, gcfg = d["pairs"], d["X"], d["summ"], d["phylo"], d["gcfg"]
    tr, te = np.asarray(job["train"]), np.asarray(job["test"])
    genome_unit = d["wells"]["genome_unit"].to_dict()
    out = []
    ctx = gml._Ctx(gcfg, pairs, X, summ, tr, te, phylo=phylo)
    ctx_s = None
    for name in job["models"]:
        try:
            if name == "clone_lookup_no_genomics":
                yhat = m_clone_lookup(ctx, genome_unit)
            elif name.startswith("SHUFFLED_"):
                if ctx_s is None:     # genomic_ml.cross_validate's control, same seed rule
                    rng = np.random.default_rng(1000 + job["fold_index"])
                    sh = pairs.copy()
                    perm = rng.permutation(tr)
                    sh.iloc[tr, sh.columns.get_loc(gcfg.target)] = pairs.iloc[perm][gcfg.target].values
                    ctx_s = gml._Ctx(gcfg, sh, X, summ, tr, te, phylo=phylo)
                yhat = gml.MODELS[name.replace("SHUFFLED_", "")](ctx_s)
            else:
                yhat = gml.MODELS[name](ctx)
        except Exception as e:
            warnings.warn(f"{name} failed on {job['tag']}: {type(e).__name__}: {e}")
            continue
        out.append(pd.DataFrame({"model": name, "pair_idx": te, "y_true": ctx.y_test, "y_pred": yhat,
                                 "n_train": len(tr)}))
    res = pd.concat(out, ignore_index=True)
    for k, v in job["tag"].items():
        res[k] = v
    return res


def run_jobs(jobs, processes=12):
    """xgboost keeps genomic_ml's 4 threads per fit (so published folds reproduce exactly);
    everything else is single-threaded (config sets the BLAS thread count)."""
    import time
    t0, parts = time.time(), []
    with get_context("fork").Pool(processes) as pool:
        for i, part in enumerate(pool.imap_unordered(run_job, jobs, chunksize=1), 1):
            parts.append(part)
            if i % 50 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} jobs, {time.time() - t0:.0f}s", flush=True)
    return pd.concat(parts, ignore_index=True)


def fold_metrics(pred, by, subset=None):
    """genomic_ml._metrics per fold (the published unit), on an optional subset of pairs."""
    p = pred if subset is None else pred[subset]
    rows = []
    for key, g in p.groupby(by, sort=False):
        rows.append({**dict(zip(by, key)), **gml._metrics(g.y_true, g.y_pred)})
    return pd.DataFrame(rows)


def summarise(fm, by, rep_col="rep"):
    """Mean over folds (as published) + the spread of that mean across CV repeats.

    Folds within a repeat share data and are not independent; the repeat-to-repeat SD of the
    per-repeat mean is the honest error bar on 'what this CV scheme scores'.
    """
    num = ["r2", "spearman_rho", "pearson_r", "sign_accuracy"]
    per_rep = fm.groupby(by + [rep_col])[num + ["n"]].agg(
        {**{c: "mean" for c in num}, "n": "sum"}).reset_index()
    agg = fm.groupby(by).agg(n_folds=("r2", "size"), n_test=("n", "sum"),
                             **{f"{c}_mean": (c, "mean") for c in num},
                             **{f"{c}_fold_sd": (c, "std") for c in num}).reset_index()
    rep_sd = per_rep.groupby(by)[num].std().add_suffix("_rep_sd").reset_index()
    return agg.merge(rep_sd, on=by)
