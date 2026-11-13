"""Final tables, from the saved fold metrics of s02 / s02b / s03 (nothing is refitted).

Metrics. Spearman rho is primary: rank-based, and it is what the published "best model" was
chosen by. R^2 (vs predicting no winner, genomic_ml._metrics) is reported as the MEDIAN over
folds, not the mean: in 20260630 a few folds per 100 blow up to R^2 ~ -1e19 (see README,
"PCA whitening"), which makes a fold mean meaningless. Every rho is a mean over 100 folds
(20 repeats x 5); `sd` is the SD of the per-repeat mean across the 20 repeats.

Outputs (outputs/tables/):
  t1_leakage.csv            step 1 per experiment x model
  t2_wrong_genomes.csv      step 2 per experiment x model (organism folds; union-fold sensitivity)
  t3_headline.csv           the three numbers per experiment for the best models
  t4_ko_vs_16s.csv          fold-paired gene content vs 16S, in each setting
  tables.md                 all of the above as markdown
"""
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import config

OUT = config.OUT / "tables"
OUT.mkdir(exist_ok=True)
EXPS = ["20260907", "20260630"]
MODELS = {"xgboost_raw_ko": "XGBoost on raw KO", "two_stage_ridge": "two-stage ridge (KO)",
          "ridge_pca": "ridge on KO PCs", "two_stage_ridge_phylo16s_only": "16S only, two-stage",
          "ridge_phylo16s_only": "16S only, ridge", "clone_lookup_no_genomics": "clone lookup (no genomics)"}
KO_BEST = "xgboost_raw_ko"
TAX_BEST = "two_stage_ridge_phylo16s_only"


def load():
    s02 = pd.concat([pd.read_csv(config.OUT / "s02_leakage" / "s02_fold_metrics.csv"),
                     pd.read_csv(config.OUT / "s02_leakage" / "s02b_fold_metrics.csv")])
    s02 = s02.rename(columns={"scheme": "arm"}).assign(test_subset="all", folds="by_well/various")
    s03 = pd.read_csv(config.OUT / "s03_wrong_genomes" / "s03_fold_metrics.csv").assign(folds="by_organism")
    u = config.OUT / "s03_wrong_genomes_union" / "s03_fold_metrics.csv"
    s03u = pd.read_csv(u).assign(folds="by_organism_union") if u.exists() else pd.DataFrame()
    s03 = pd.concat([s03, s03u], ignore_index=True)
    for d in (s02, s03):
        d["exp"] = d["exp"].astype(str)          # read back from CSV as int
    return s02, s03


def cell(f, met="spearman_rho"):
    """mean over folds, SD of the per-repeat mean across repeats, median R^2, sign accuracy."""
    rep = f.groupby("rep")[met].mean()
    return dict(rho=f[met].mean(), rho_sd=rep.std(), r2_median=f.r2.median(),
                r2_mean=f.r2.mean(), sign_acc=f.sign_accuracy.mean(), n_folds=len(f))


def paired(fa, fb, met="spearman_rho"):
    """fa - fb on identical folds (rep, k)."""
    a = fa.set_index(["rep", "k"])[met]
    b = fb.set_index(["rep", "k"])[met]
    d = (a - b).dropna()
    return dict(delta=d.mean(), se=d.std() / np.sqrt(len(d)), p=wilcoxon(d).pvalue if (d != 0).any() else np.nan,
                frac_better=(d > 0).mean())


def sel(df, **kw):
    m = np.ones(len(df), bool)
    for k, v in kw.items():
        m &= (df[k] == v).values
    return df[m]


def t1(s02):
    rows = []
    for exp in EXPS:
        for model, lab in MODELS.items():
            g = sel(s02, exp=exp, model=model)
            r = dict(exp=exp, model=model, label=lab)
            for arm in ("by_well", "by_well_random", "by_well_purged", "by_organism", "by_genome", "by_organism_union"):
                c = cell(sel(g, arm=arm))
                r[f"{arm}_rho"], r[f"{arm}_rho_sd"], r[f"{arm}_r2_med"] = c["rho"], c["rho_sd"], c["r2_median"]
            r["published_r2_mean"] = cell(sel(g, arm="by_well"))["r2_mean"]
            for met, tag in (("spearman_rho", "rho"), ("r2", "r2")):
                lk = paired(sel(g, arm="by_well_random"), sel(g, arm="by_well_purged"), met)
                r[f"leak_{tag}"], r[f"leak_{tag}_se"], r[f"leak_{tag}_p"] = lk["delta"], lk["se"], lk["p"]
            rows.append(r)
    return pd.DataFrame(rows)


def t2(s03):
    rows = []
    for folds in s03.folds.unique():
        for exp in EXPS:
            for model, lab in MODELS.items():
                if model == "clone_lookup_no_genomics":
                    continue
                g = sel(s03, exp=exp, model=model, folds=folds)
                for arm, sub in (("drop_wrong", "clean"), ("ok_only", "ok")):
                    if not len(sel(g, arm=arm)):
                        continue
                    base, new, rnd = (sel(g, arm=a, test_subset=sub) for a in ("all", arm, f"{arm}_random"))
                    c0, c1, cr = cell(base), cell(new), cell(rnd)
                    d, dr = paired(new, base), paired(rnd, base)
                    rows.append(dict(folds=folds, exp=exp, model=model, label=lab, removed=arm, test_subset=sub,
                                     trained_all_rho=c0["rho"], trained_all_rho_sd=c0["rho_sd"],
                                     removed_rho=c1["rho"], removed_rho_sd=c1["rho_sd"],
                                     random_rho=cr["rho"], delta_rho=d["delta"], delta_rho_se=d["se"],
                                     delta_rho_p=d["p"], random_delta_rho=dr["delta"], random_delta_rho_se=dr["se"],
                                     trained_all_r2_med=c0["r2_median"], removed_r2_med=c1["r2_median"],
                                     random_r2_med=cr["r2_median"],
                                     trained_all_sign=c0["sign_acc"], removed_sign=c1["sign_acc"],
                                     random_sign=cr["sign_acc"],
                                     # how well were the wrong wells' pairs ever predicted?
                                     wrong_pairs_rho=cell(sel(g, arm="all", test_subset="wrong"))["rho"]))
    return pd.DataFrame(rows)


def t3(s02, s03):
    """Per experiment, best KO model and the 16S control: published -> leak removed -> wrong
    genomes removed. Each step is fold-paired against the one it changes."""
    rows = []
    o = sel(s03, folds="by_organism")
    for exp in EXPS:
        for model in (KO_BEST, "two_stage_ridge", TAX_BEST):
            g2, g3 = sel(s02, exp=exp, model=model), sel(o, exp=exp, model=model)
            steps = [
                ("A  published: folds by well, all wells, all test pairs", sel(g2, arm="by_well")),
                ("B  leak removed, same test pairs (purged; size-matched control = A')", sel(g2, arm="by_well_purged")),
                ("A' size-matched control for B (leak kept, same data loss)", sel(g2, arm="by_well_random")),
                ("C  leak-free folds by organism, all wells, all test pairs", sel(g3, arm="all", test_subset="all")),
                ("D  leak-free, all wells, test on correct-genome pairs only", sel(g3, arm="all", test_subset="clean")),
                ("E  leak-free, wrong-genome wells removed, same test pairs as D", sel(g3, arm="drop_wrong", test_subset="clean")),
            ]
            for name, f in steps:
                c = cell(f)
                rows.append(dict(exp=exp, model=model, step=name, rho=c["rho"], rho_sd=c["rho_sd"],
                                 r2_median=c["r2_median"], sign_acc=c["sign_acc"], n_folds=c["n_folds"],
                                 n_test_pair_folds=int(f.n.sum())))
    return pd.DataFrame(rows)


def t4(s02, s03):
    rows = []
    o = sel(s03, folds="by_organism")
    settings = [("published (by_well)", s02, dict(arm="by_well")),
                ("leak removed (purged)", s02, dict(arm="by_well_purged")),
                ("leak-free, all wells, correct-genome test pairs", o, dict(arm="all", test_subset="clean")),
                ("leak-free, wrong genomes removed, same test pairs", o, dict(arm="drop_wrong", test_subset="clean"))]
    for exp in EXPS:
        for name, df, kw in settings:
            ko, tx = sel(df, exp=exp, model=KO_BEST, **kw), sel(df, exp=exp, model=TAX_BEST, **kw)
            for met in ("spearman_rho", "sign_accuracy"):
                p = paired(ko, tx, met)
                rows.append(dict(exp=exp, setting=name, metric=met, ko=ko[met].mean(), s16=tx[met].mean(),
                                 ko_minus_16s=p["delta"], se=p["se"], p=p["p"], frac_folds_ko_better=p["frac_better"]))
    return pd.DataFrame(rows)


def md(df, cols, fmt):
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(fmt(c, r[c]) for c in cols) + " |")
    return "\n".join(lines)


def f3(c, v):
    if isinstance(v, (float, np.floating)):
        if c.endswith("_p") or c == "p":
            return f"{v:.1e}" if v < 1e-3 else f"{v:.3f}"
        return "" if pd.isna(v) else f"{v:.3f}"
    return str(v)


def main():
    s02, s03 = load()
    a, b, c, d = t1(s02), t2(s03), t3(s02, s03), t4(s02, s03)
    for name, df in (("t1_leakage", a), ("t2_wrong_genomes", b), ("t3_headline", c), ("t4_ko_vs_16s", d)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    txt = ["# 20261113 tables (generated by s05_tables.py)\n",
           "## t3 headline\n", md(c, ["exp", "model", "step", "rho", "rho_sd", "r2_median", "sign_acc"], f3),
           "\n## t1 leakage\n", md(a, ["exp", "label", "by_well_rho", "by_well_random_rho", "by_well_purged_rho",
                                       "leak_rho", "leak_rho_se", "leak_rho_p", "leak_r2", "leak_r2_se",
                                       "by_organism_rho", "by_genome_rho", "by_organism_union_rho"], f3),
           "\n## t2 wrong genomes\n", md(b, ["folds", "exp", "label", "removed", "trained_all_rho", "removed_rho",
                                             "random_rho", "delta_rho", "delta_rho_se", "delta_rho_p",
                                             "trained_all_r2_med", "removed_r2_med", "trained_all_sign",
                                             "removed_sign", "wrong_pairs_rho"], f3),
           "\n## t4 KO vs 16S (fold-paired)\n", md(d, ["exp", "setting", "metric", "ko", "s16", "ko_minus_16s",
                                                       "se", "p", "frac_folds_ko_better"], f3)]
    (OUT / "tables.md").write_text("\n".join(txt) + "\n")
    print("\n".join(txt))


if __name__ == "__main__":
    main()
