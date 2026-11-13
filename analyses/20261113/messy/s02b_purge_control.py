"""Step 1b: the purge's data-loss control.

s02's `by_well_purged` arm drops every training pair involving an organism present in the test
fold -- about a third of training (366 -> 244 pairs in 20260907, 953 -> 609 in 20260630). So
its drop mixes two things: the leak going away, and simply having less data.

`by_well_random` removes the same number of training pairs per fold, by dropping random whole
wells whose organism is NOT in the test fold. It keeps the leak and loses as much data, so

    leak            = by_well_random - by_well_purged       (same test pairs, same training size)
    data-loss part  = by_well - by_well_random

Outputs: outputs/s02_leakage/
  s02b_predictions.csv.gz, s02b_fold_metrics.csv
  s02b_leak_decomposition.csv
"""
import sys

import config  # first: caps BLAS threads before numpy loads
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import lib

N_REP = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 20
OUT = config.OUT / "s02_leakage"


def jobs_for(exp):
    d = lib.load(exp, OUT / exp)
    pairs, wells = d["pairs"], d["wells"]
    org = wells["cv_group"].to_dict()
    ident = {s: s for s in wells.index}
    jobs, fi = [], 0
    for rep in range(N_REP):
        for k, tr, te in lib.group_folds(pairs, ident, seed=97 * rep):
            test_orgs = {org[s] for s in lib.strains_of(pairs, te)}
            p = pairs.iloc[tr]
            n_purged = int((p.strain_a.map(org).isin(test_orgs) | p.strain_b.map(org).isin(test_orgs)).sum())
            allowed = {s: org[s] not in test_orgs for s in wells.index}
            tr_r = lib.random_matched(pairs, tr, allowed, n_purged, np.random.default_rng(20_000 + fi))
            tag = dict(exp=exp, scheme="by_well_random", rep=rep, k=k, fold_index=fi)
            jobs.append(dict(exp=exp, train=tr_r, test=te, models=lib.MODELS, fold_index=fi, tag=tag))
            fi += 1
    return jobs


def main():
    if (OUT / "s02b_fold_metrics.csv").exists() and "--rerun" not in sys.argv:
        pred = pd.read_csv(OUT / "s02b_predictions.csv.gz")    # folds already run: only summarise
        fm = pd.read_csv(OUT / "s02b_fold_metrics.csv")
    else:
        jobs = jobs_for("20260907") + jobs_for("20260630")
        print(f"{len(jobs)} fold-jobs")
        pred = lib.run_jobs(jobs)
        pred.to_csv(OUT / "s02b_predictions.csv.gz", index=False)
        fm = lib.fold_metrics(pred, ["exp", "scheme", "rep", "k", "model"])
        fm.to_csv(OUT / "s02b_fold_metrics.csv", index=False)
    fm = pd.concat([pd.read_csv(OUT / "s02_fold_metrics.csv"), fm], ignore_index=True)
    fm = fm[fm.scheme.isin(["by_well", "by_well_random", "by_well_purged"])]

    ntr = pred.drop_duplicates(["exp", "rep", "k"]).groupby("exp").n_train.mean()
    print("mean training pairs, random arm:", ntr.round(0).to_dict())

    rows = []
    for (exp, model), g in fm.groupby(["exp", "model"]):
        if model == "zero_baseline":
            continue
        for met in ("r2", "spearman_rho", "sign_accuracy"):
            # pivot_table drops a scheme whose metric is NaN in every fold (a constant predictor)
            w = g.pivot_table(index=["rep", "k"], columns="scheme", values=met, dropna=False).reindex(
                columns=["by_well", "by_well_random", "by_well_purged"])
            if w.isna().all().any():
                continue
            leak = (w.by_well_random - w.by_well_purged).dropna()
            loss = (w.by_well - w.by_well_random).dropna()
            rows.append(dict(exp=exp, model=model, metric=met, by_well=w.by_well.mean(),
                             by_well_random=w.by_well_random.mean(), by_well_purged=w.by_well_purged.mean(),
                             leak=leak.mean(), leak_se=leak.std() / np.sqrt(len(leak)),
                             leak_wilcoxon_p=wilcoxon(leak).pvalue if (leak != 0).any() else np.nan,
                             data_loss=loss.mean(), data_loss_se=loss.std() / np.sqrt(len(loss)),
                             n_folds=len(leak)))
    dec = pd.DataFrame(rows)
    dec.to_csv(OUT / "s02b_leak_decomposition.csv", index=False)
    pd.set_option("display.width", 220)
    print(dec[dec.metric != "sign_accuracy"].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
