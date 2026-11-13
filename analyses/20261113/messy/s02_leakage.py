"""Step 1: how much of the published cv_strain performance was leakage?

Published cv_strain held out WELLS. Wells carrying the same genome (>= 99.9% ANI clone
cluster) -- or the same overgrowing organism -- could sit on both sides of a split, so a test
well's exact genome and phenotype were often already in training. Same data, same models,
same code; only the fold assignment changes:

  by_well            published scheme (reps 0-2 reproduce the published per-fold numbers)
  by_organism        hold out the organism in the well together (s01 cv_group)   <- primary
  by_genome          hold out genome clone clusters together (cv_group_genome)
  by_organism_union  genome clusters + invader links, connected (cv_group_union)

and one PAIRED arm, which takes the test set as given and removes the leak from training:

  by_well_purged     each by_well fold, with every training pair that involves an organism
                     present among the fold's test wells dropped. Identical test pairs to
                     by_well, so the fold-by-fold difference is the leak and nothing else.

N_REP repeats of 5-fold CV per scheme (published: 3). The unit of the published numbers is
"mean over folds"; the error bar here is the SD of that mean across repeats.

Outputs: outputs/s02_leakage/
  s02_predictions.csv.gz   every test pair, every model, every fold/scheme
  s02_fold_metrics.csv     genomic_ml._metrics per fold
  s02_summary.csv          per experiment x scheme x model
  s02_purge_paired.csv     by_well vs by_well_purged, fold-paired
  s02_exposure.csv         by_well test pairs split by whether their organism was in training
"""
import sys

import config  # first: caps BLAS threads before numpy loads
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import lib
import genomic_ml as gml

N_REP = int(sys.argv[1]) if len(sys.argv) > 1 else 20
OUT = config.OUT / "s02_leakage"
SCHEMES = {"by_well": None, "by_organism": "cv_group", "by_genome": "cv_group_genome",
           "by_organism_union": "cv_group_union"}


def jobs_for(exp):
    d = lib.load(exp, OUT / exp)
    pairs, wells = d["pairs"], d["wells"]
    org = wells["cv_group"].to_dict()
    jobs = []
    for scheme, col in SCHEMES.items():
        group_of = {s: s for s in wells.index} if col is None else wells[col].to_dict()
        fi = 0
        for rep in range(N_REP):
            for k, tr, te in lib.group_folds(pairs, group_of, seed=97 * rep):
                tag = dict(exp=exp, scheme=scheme, rep=rep, k=k, fold_index=fi)
                jobs.append(dict(exp=exp, train=tr, test=te, models=lib.MODELS, fold_index=fi, tag=tag))
                if scheme == "by_well":
                    test_orgs = {org[s] for s in lib.strains_of(pairs, te)}
                    p = pairs.iloc[tr]
                    keep = ~(p.strain_a.map(org).isin(test_orgs) | p.strain_b.map(org).isin(test_orgs)).values
                    jobs.append(dict(exp=exp, train=tr[keep], test=te, models=lib.MODELS, fold_index=fi,
                                     tag=dict(tag, scheme="by_well_purged")))
                fi += 1
    return jobs


def exposure(pred, exp):
    """For by_well folds: was the test pair's organism (either strain) present in training?"""
    d = lib._D[exp]
    pairs, org = d["pairs"], d["wells"]["cv_group"].to_dict()
    rows = []
    for (rep, k), g in pred[(pred.exp == exp) & (pred.scheme == "by_well")].groupby(["rep", "k"]):
        te = np.unique(g.pair_idx)
        # rebuild this fold's training strains from the fold assignment
        folds = {kk: (tr, tt) for kk, tr, tt in lib.group_folds(pairs, {s: s for s in d["wells"].index}, seed=97 * rep)}
        tr_orgs = {org[s] for s in lib.strains_of(pairs, folds[k][0])}
        p = pairs.iloc[te]
        n_exp = p.strain_a.map(org).isin(tr_orgs).astype(int) + p.strain_b.map(org).isin(tr_orgs).astype(int)
        rows.append(pd.DataFrame({"exp": exp, "rep": rep, "k": k, "pair_idx": te, "n_strains_exposed": n_exp.values}))
    return pd.concat(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = jobs_for("20260907") + jobs_for("20260630")
    print(f"{len(jobs)} fold-jobs, {N_REP} repeats per scheme")
    pred = lib.run_jobs(jobs)
    pred.to_csv(OUT / "s02_predictions.csv.gz", index=False)

    fm = lib.fold_metrics(pred, ["exp", "scheme", "rep", "k", "model"])
    fm.to_csv(OUT / "s02_fold_metrics.csv", index=False)
    summ = lib.summarise(fm, ["exp", "scheme", "model"])
    summ.to_csv(OUT / "s02_summary.csv", index=False)

    # --- reproduction check: reps 0-2 of by_well are the published per-fold numbers ---------
    for exp in ("20260907", "20260630"):
        e = config.EXPERIMENTS[exp]
        pub = pd.read_csv(e["pub"] / f"{e['pub_prefix']}_cv_per_fold.csv")
        pub = pub[pub.regime == "cv_strain"].rename(columns={"fold": "fold_index"})
        mine = fm[(fm.exp == exp) & (fm.scheme == "by_well") & (fm.rep < 3)].copy()
        mine["fold_index"] = mine.rep * 5 + mine.k
        j = mine.merge(pub[["fold_index", "model", "r2", "spearman_rho"]], on=["fold_index", "model"],
                       suffixes=("", "_pub"))
        dev = np.nanmax(np.abs(j[["r2", "spearman_rho"]].values - j[["r2_pub", "spearman_rho_pub"]].values))
        print(f"{exp}: reproduced {len(j)} published fold x model rows, max |diff| = {dev:.2e}")

    # --- paired: same test pairs, leak purged from training -----------------------------------
    rows = []
    for (exp, model), g in fm[fm.scheme.isin(["by_well", "by_well_purged"])].groupby(["exp", "model"]):
        if model == "zero_baseline":
            continue
        w = g.pivot_table(index=["rep", "k"], columns="scheme", values=["r2", "spearman_rho"])
        for met in ("r2", "spearman_rho"):
            dlt = (w[(met, "by_well")] - w[(met, "by_well_purged")]).dropna()
            rows.append(dict(exp=exp, model=model, metric=met, by_well=w[(met, "by_well")].mean(),
                             by_well_purged=w[(met, "by_well_purged")].mean(), leak=dlt.mean(),
                             leak_se=dlt.std() / np.sqrt(len(dlt)), folds_leaky_better=int((dlt > 0).sum()),
                             n_folds=len(dlt), wilcoxon_p=wilcoxon(dlt).pvalue if (dlt != 0).any() else np.nan))
    paired = pd.DataFrame(rows)
    paired.to_csv(OUT / "s02_purge_paired.csv", index=False)

    # --- where the leaky score came from: exposed vs unexposed test pairs ----------------------
    ex = pd.concat([exposure(pred, e) for e in ("20260907", "20260630")])
    bw = pred[pred.scheme == "by_well"].merge(ex, on=["exp", "rep", "k", "pair_idx"])
    bw["exposed"] = np.where(bw.n_strains_exposed > 0, "organism seen in training", "not seen")
    erows = []
    for (exp, model, e), g in bw.groupby(["exp", "model", "exposed"]):
        erows.append(dict(exp=exp, model=model, exposed=e, n_pair_folds=len(g),
                          **{k: v for k, v in gml._metrics(g.y_true, g.y_pred).items() if k != "n"}))
    expo = pd.DataFrame(erows)
    expo.to_csv(OUT / "s02_exposure.csv", index=False)

    pd.set_option("display.width", 220)
    show = ["exp", "scheme", "model", "n_folds", "n_test", "r2_mean", "r2_rep_sd",
            "spearman_rho_mean", "spearman_rho_rep_sd", "sign_accuracy_mean"]
    print(summ[show].sort_values(["exp", "model", "scheme"]).round(3).to_string(index=False))
    print("\npaired (same test pairs):")
    print(paired.round(4).to_string(index=False))
    print("\nexposure (by_well folds, pooled):")
    print(expo.round(3).to_string(index=False))
    share = bw.drop_duplicates(["exp", "rep", "k", "pair_idx"]).groupby("exp").exposed.value_counts(normalize=True)
    print("\nshare of by_well test pairs whose organism was in training:\n", share.round(3).to_string())


if __name__ == "__main__":
    main()
