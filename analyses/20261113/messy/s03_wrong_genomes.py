"""Step 2: does removing wells whose genome is not the organism in the well help the model?

Leak-free CV throughout (by_organism folds from s02). The genome verdict per well is s01's:
20260907 from 20261021's concordance, 20260630 from its own references (s01).

Removing wells changes BOTH the training data and the test set, and a smaller or easier test
set can move a score without the model getting any better. So within each fold the test set
is fixed and only training varies:

  training arm          drops from training ...
    all                 nothing (published data)
    drop_wrong          every pair involving a `wrong` well
    drop_wrong_random   random NON-wrong wells, until at least as many training pairs are gone
                        as drop_wrong removed -- the "is it just less data?" control
    ok_only             every pair involving a well that is not verified `ok` (wrong + unverified)
    ok_only_random      matched random control for ok_only

  test subset (same pairs for every arm)
    clean               neither well wrong          <- drop_wrong is judged here
    ok                  both wells verified ok      <- ok_only is judged here
    wrong               at least one well wrong     (how well were these ever predicted?)
    all                 every test pair

Outputs: outputs/s03_wrong_genomes/
  s03_predictions.csv.gz, s03_fold_metrics.csv
  s03_summary.csv        experiment x arm x test subset x model
  s03_paired.csv         arm vs `all` on the arm's own test subset, fold-paired
"""
import sys

import config  # first: caps BLAS threads before numpy loads
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import lib

N_REP = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 20
# --union: sensitivity run on cv_group_union folds, where a wrong well is held out together with
# BOTH its genome's clone-mates and its invader, so no wrong well can ever sit in training while
# a clone-mate of its genome is being tested. Isolates the "general noise" part of the gain.
UNION = "--union" in sys.argv
GROUP_COL = "cv_group_union" if UNION else "cv_group"
OUT = config.OUT / ("s03_wrong_genomes_union" if UNION else "s03_wrong_genomes")
ARMS = {"drop_wrong": lambda st: st != "wrong"}
if not UNION:
    ARMS["ok_only"] = lambda st: st == "ok"
JUDGED_ON = {"all": "all", "drop_wrong": "clean", "drop_wrong_random": "clean",
             "ok_only": "ok", "ok_only_random": "ok"}
if "--union" in sys.argv:
    JUDGED_ON = {k: v for k, v in JUDGED_ON.items() if not k.startswith("ok_only")}


def jobs_for(exp):
    d = lib.load(exp, OUT / exp)
    pairs, wells = d["pairs"], d["wells"]
    st = wells["genome_status"].to_dict()
    group_of = wells[GROUP_COL].to_dict()
    jobs = []
    fi = 0
    for rep in range(N_REP):
        for k, tr, te in lib.group_folds(pairs, group_of, seed=97 * rep):
            tag = dict(exp=exp, rep=rep, k=k, fold_index=fi)
            jobs.append(dict(exp=exp, train=tr, test=te, models=lib.MODELS, fold_index=fi,
                             tag=dict(tag, arm="all")))
            p = pairs.iloc[tr]
            for arm, rule in ARMS.items():
                keep_ok = {s: rule(v) for s, v in st.items()}
                keep = (p.strain_a.map(keep_ok) & p.strain_b.map(keep_ok)).values
                jobs.append(dict(exp=exp, train=tr[keep], test=te, models=lib.MODELS, fold_index=fi,
                                 tag=dict(tag, arm=arm, n_train_removed=int((~keep).sum()))))
                rng = np.random.default_rng(10_000 + fi)
                tr_r = lib.random_matched(pairs, tr, keep_ok, int((~keep).sum()), rng)
                jobs.append(dict(exp=exp, train=tr_r, test=te, models=lib.MODELS, fold_index=fi,
                                 tag=dict(tag, arm=f"{arm}_random", n_train_removed=len(tr) - len(tr_r))))
            fi += 1
    return jobs


def label_subsets(pred):
    out = []
    for exp, g in pred.groupby("exp"):
        d = lib._D[exp]
        st = d["wells"]["genome_status"]
        p = d["pairs"]
        sa = p.strain_a.map(st).values[g.pair_idx.values]
        sb = p.strain_b.map(st).values[g.pair_idx.values]
        g = g.copy()
        g["any_wrong"] = (sa == "wrong") | (sb == "wrong")
        g["both_ok"] = (sa == "ok") & (sb == "ok")
        out.append(g)
    return pd.concat(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "s03_fold_metrics.csv").exists() and "--rerun" not in sys.argv:
        for exp in ("20260907", "20260630"):              # folds already run: only summarise
            lib.load(exp, OUT / exp)
        pred = pd.read_csv(OUT / "s03_predictions.csv.gz")
        fm = pd.read_csv(OUT / "s03_fold_metrics.csv")
    else:
        jobs = jobs_for("20260907") + jobs_for("20260630")
        print(f"{len(jobs)} fold x arm jobs, {N_REP} repeats")
        for exp in ("20260907", "20260630"):
            w = lib._D[exp]["wells"]
            print(f"  {exp}: wells {len(w)}, genome status {w.genome_status.value_counts().to_dict()}")
        pred = label_subsets(lib.run_jobs(jobs))
        pred.to_csv(OUT / "s03_predictions.csv.gz", index=False)

        subsets = {"all": None, "clean": ~pred.any_wrong, "ok": pred.both_ok, "wrong": pred.any_wrong}
        fms = []
        for name, mask in subsets.items():
            fm = lib.fold_metrics(pred, ["exp", "arm", "rep", "k", "model"], subset=mask)
            fm["test_subset"] = name
            fms.append(fm)
        fm = pd.concat(fms, ignore_index=True)
        fm.to_csv(OUT / "s03_fold_metrics.csv", index=False)
    summ = lib.summarise(fm, ["exp", "arm", "test_subset", "model"])
    summ.to_csv(OUT / "s03_summary.csv", index=False)

    rows = []
    for arm, sub in JUDGED_ON.items():
        if arm == "all":
            continue
        f = fm[(fm.test_subset == sub) & fm.arm.isin(["all", arm])]
        for (exp, model), g in f.groupby(["exp", "model"]):
            if model == "zero_baseline":
                continue
            for met in ("r2", "spearman_rho", "sign_accuracy"):
                # pivot_table drops an arm whose metric is NaN in every fold (a constant predictor)
                w = g.pivot_table(index=["rep", "k"], columns="arm", values=met, dropna=False).reindex(
                    columns=["all", arm])
                if w.isna().all().any():
                    continue
                dlt = (w[arm] - w["all"]).dropna()
                rows.append(dict(exp=exp, model=model, arm=arm, test_subset=sub, metric=met,
                                 trained_on_all=w["all"].mean(), trained_on_arm=w[arm].mean(),
                                 delta=dlt.mean(), delta_se=dlt.std() / np.sqrt(len(dlt)),
                                 folds_arm_better=int((dlt > 0).sum()), n_folds=len(dlt),
                                 wilcoxon_p=wilcoxon(dlt).pvalue if (dlt != 0).any() else np.nan))
    paired = pd.DataFrame(rows)
    paired.to_csv(OUT / "s03_paired.csv", index=False)

    pd.set_option("display.width", 220)
    show = ["exp", "arm", "test_subset", "model", "n_test", "r2_mean", "r2_rep_sd",
            "spearman_rho_mean", "spearman_rho_rep_sd", "sign_accuracy_mean"]
    s = summ[[a in JUDGED_ON and (t == JUDGED_ON[a] or t in ("all", "wrong")) for a, t in zip(summ.arm, summ.test_subset)]]
    print(s[show].sort_values(["exp", "model", "test_subset", "arm"]).round(3).to_string(index=False))
    print("\npaired, arm minus all, on the arm's test subset:")
    print(paired[paired.metric != "sign_accuracy"].round(4).to_string(index=False))
    rem = pred.drop_duplicates(["exp", "rep", "k", "arm"]).groupby(["exp", "arm"]).n_train_removed.mean()
    print("\nmean training pairs removed per fold:\n", rem.round(1).to_string())


if __name__ == "__main__":
    main()
