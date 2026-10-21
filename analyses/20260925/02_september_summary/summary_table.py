"""One table: 20260907 genomic-ML performance against every reference point.

Nothing is refitted. ML rows are read from existing cross-validation outputs; the replicate
rows are computed here, but with genomic_ml._metrics itself, so every row uses the same R^2
and sign definitions as the ML:

    R^2         1 - SSE / sum(y^2), i.e. against predicting "no winner" (log2 ratio 0),
                NOT against the mean. The zero baseline scores exactly 0 by construction.
    sign acc.   winner called correctly, over pairs with |log2 ratio| > 1 only.
    target      log2(a/b) per pair, mean over 1:1 standard wells.

ML rows are the mean over 15 cv_strain folds (5 folds x 3 repeats, both strains of every
test pair unseen in training); replicate rows are pooled over all pairs. Sources:

    20260907 ML               analyses/20260923/01_model_sweep/outputs/s01_*   (sweep; same folds for all models)
    20260907 + inoculum OD    analyses/20260923/02_sensitivity/outputs/s02_C_density_as_feature_*
    20260907 wells            analyses/20260921/analysis/relative_abundance/outputs/r02_*
    cross-experiment          analyses/20260925/01_cross_experiment/outputs/
    20260630 ML               analyses/20260710_increased_tolerance/06_ml_analysis/genomic_ml/outputs_merged_ab/
                              analyses/20260710/06_ml_analysis/genomic_ml/outputs/  (as published)

    python summary_table.py      (env karl_seq_analysis, seconds)
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import config  # noqa: E402  (puts scripts/ on sys.path)
import paths   # noqa: E402
from genomic_ml import _metrics  # noqa: E402

OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)
AN = paths.ANALYSES
SWEEP = AN / "20260923" / "01_model_sweep" / "outputs"
SENS = AN / "20260923" / "02_sensitivity" / "outputs"
XEXP = AN / "20260925" / "01_cross_experiment" / "outputs"
ML_0630 = {
    "20260630 (merged demux, current)": AN / "20260710_increased_tolerance" / "06_ml_analysis" / "genomic_ml" / "outputs_merged_ab",
    "20260630 (as published, original demux)": AN / "20260710" / "06_ml_analysis" / "genomic_ml" / "outputs",
}
COLS = ["r2", "pearson_r", "spearman_rho", "sign_accuracy"]
rows = []


def cv_row(group, label, summ, model, n_pairs=None, note=""):
    r = summ[(summ.regime == "cv_strain") & (summ.model == model)].iloc[0]
    rows.append({"group": group, "row": label, "model": model, "n_test": int(r.n_test),
                 "r2": r.r2_mean, "r2_sd": r.r2_sd, "pearson_r": r.pearson_r_mean,
                 "spearman_rho": r.spearman_rho_mean, "rho_sd": r.spearman_rho_sd,
                 "sign_accuracy": r.sign_accuracy_mean, "note": note})


def pooled_row(group, label, y, yhat, note=""):
    m = _metrics(y, yhat)
    rows.append({"group": group, "row": label, "model": "", "n_test": m["n"],
                 **{k: m[k] for k in COLS}, "note": note})


# --- the September ML ---------------------------------------------------------------------------
s01 = pd.read_csv(SWEEP / "s01_cv_summary.csv")
genome_models = ["xgboost_raw_ko", "two_stage_xgboost", "ridge_pca", "two_stage_ridge",
                 "xgboost_pca", "two_stage_ridge_ko_plus_phylo"]
cs = s01[(s01.regime == "cv_strain") & s01.model.isin(genome_models)]
BEST = cs.sort_values("spearman_rho_mean", ascending=False).model.iloc[0]
cv_row("20260907 ML", "best genome model (KEGG KO)", s01, BEST,
       note="best genome-only model in the 20260923 sweep, by rho")
cv_row("20260907 ML", "first-pass model", s01, "ridge_pca", note="20260921 first pass; same folds")
s02c = pd.read_csv(SENS / "s02_C_density_as_feature_cv_summary.csv")
cv_row("20260907 ML", "genome + preculture OD (not genome-only)", s02c, "ridge_pca",
       note="source-plate OD added as a feature; 20260923 sensitivity arm C")

# --- 1: floor ------------------------------------------------------------------------------------
cv_row("1 floor", "zero baseline (predict no winner)", s01, "zero_baseline", note="R^2 = 0 by definition")
cv_row("1 floor", "label-shuffled control", s01, "SHUFFLED_ridge_pca",
       note="one permutation per fold; its spread across seeds is ~+-0.1 rho")

# --- 2: technical replicates within 20260907 -------------------------------------------------------
w = pd.read_csv(config.RA_OUT[config.B] / "r02_well_interaction_scores.csv")
lay = pd.read_csv(config.LAYOUT[config.B])
lay["sample_id"] = [f"Plate{int(p):02d}_{x}" for p, x in zip(lay.dest_plate, lay.dest_well)]
w = w.drop(columns="well_type").merge(
    lay[["sample_id", "well_type", "vol1_nL", "vol2_nL", "log2_inoculum_ratio"]], on="sample_id")
w = w[(w.vol1_nL == 100) & (w.vol2_nL == 100) & (w.log2_inoculum_ratio == 0) & ~w.missing_reference.astype(bool)]
model_pairs = pd.read_csv(SWEEP / "g01_modeling_pairs.csv")
keys = set(zip(model_pairs.strain_a, model_pairs.strain_b))
w = w[[k in keys for k in zip(w.strain_a, w.strain_b)]]
rng = np.random.default_rng(0)
xs, ys = [], []
for _, g in w.groupby(["strain_a", "strain_b"]):
    v = g.log2_ratio_a_over_b.to_numpy()
    for i, j in itertools.combinations(range(len(v)), 2):
        a, b = (i, j) if rng.random() < .5 else (j, i)
        xs.append(v[a]); ys.append(v[b])
n_rep_pairs = w.groupby(["strain_a", "strain_b"]).size().ge(2).sum()
pooled_row("2 technical replicates (within 20260907)", "one well predicts another well of the same pair",
           ys, xs, note=f"{len(xs)} well pairs from {n_rep_pairs} modelled pairs with >=2 wells")
g02 = pd.read_csv(SWEEP / "g02_label_noise_ceiling.csv").set_index("estimator")
for est, lab in (("replicate_spread", "noise ceiling from replicate spread"),
                 ("binomial_read_depth", "noise ceiling from read depth")):
    rows.append({"group": "2 technical replicates (within 20260907)", "row": lab, "model": "",
                 "n_test": int(g02.loc[est, "n"]), "r2": g02.loc[est, "ceiling_r2"],
                 "note": "best R^2 any model could reach given label noise (mean-referenced)"})

# --- 3: biological replicates (cross-experiment) -----------------------------------------------------
d = pd.read_csv(XEXP / "x02_direct_pairs_values.csv")
pooled_row("3 biological replicates (20260630 -> 20260907)", "same pair, July pair mean predicts September",
           d.log2_ab_B, d.log2_ab_A, note="53 pairs tested in both; 14 mismatched strain identities excluded")
t = pd.read_csv(XEXP / "x04_predictions_bt_20260630_to_20260907.csv")
pooled_row("3 biological replicates (20260630 -> 20260907)", "July per-strain strengths predict September pairs",
           t.log2_ab, t.pred_log2, note="Bradley-Terry fit on July; mostly pairs July also tested")
t = pd.read_csv(XEXP / "x04_predictions_bt_20260907_to_20260630.csv")
t = t[~t.in_fit_exp]
pooled_row("3 biological replicates (20260630 -> 20260907)", "September strengths predict UNSEEN July pairs",
           t.log2_ab, t.pred_log2, note="pairs September never tested; limited by its sparse design (70 fitted pairs)")

# --- 4: previous experiments --------------------------------------------------------------------------
for lab, p in ML_0630.items():
    s = pd.read_csv(p / "g03_cv_summary.csv")
    c = s[(s.regime == "cv_strain") & s.model.isin(genome_models)]
    best = c.sort_values("spearman_rho_mean", ascending=False).model.iloc[0]
    best_r2 = c.sort_values("r2_mean", ascending=False).model.iloc[0]
    cv_row("4 previous experiments", f"{lab}: best genome model (by rho)", s, best)
    if best_r2 != best:   # the two criteria pick different models in 20260630
        cv_row("4 previous experiments", f"{lab}: best genome model (by R^2)", s, best_r2)
    cv_row("4 previous experiments", f"{lab}: 16S-only", s, "two_stage_ridge_phylo16s_only")
rows.append({"group": "4 previous experiments", "row": "20260721", "model": "",
             "note": "not modellable: its well labels fail the genome-join gate (source-plate mix-up)"})

# --- 5: taxonomy only, same folds ------------------------------------------------------------------------
for m in ("two_stage_ridge_phylo16s_only", "ridge_phylo16s_only"):
    cv_row("5 taxonomy only (same folds)", f"16S distance only: {m}", s01, m,
           note="same 15 folds as the genome models")

# paired, fold-by-fold: best genome model vs the better taxonomy model
f = pd.read_csv(SWEEP / "s01_cv_per_fold.csv")
f = f[f.regime == "cv_strain"]
tax = cs = s01[(s01.regime == "cv_strain") & s01.model.str.contains("phylo16s_only")].sort_values(
    "spearman_rho_mean", ascending=False).model.iloc[0]
paired = []
for met in ("spearman_rho", "r2"):
    a = f[f.model == BEST].set_index("fold")[met]
    b = f[f.model == tax].set_index("fold")[met]
    dlt = a - b
    paired.append({"comparison": f"{BEST} minus {tax}", "metric": met, "mean_delta": dlt.mean(),
                   "se": dlt.std() / np.sqrt(len(dlt)), "folds_better": int((dlt > 0).sum()),
                   "n_folds": len(dlt), "wilcoxon_p": wilcoxon(dlt).pvalue})
paired = pd.DataFrame(paired)

# --- write ------------------------------------------------------------------------------------------------
tab = pd.DataFrame(rows)
tab = tab[["group", "row", "model", "n_test", "r2", "r2_sd", "pearson_r", "spearman_rho", "rho_sd",
           "sign_accuracy", "note"]]
tab.to_csv(OUT / "september_summary.csv", index=False)
paired.to_csv(OUT / "september_best_vs_taxonomy_paired.csv", index=False)


def fmt(v, sd=None):
    if pd.isna(v):
        return ""
    return f"{v:.2f}" + (f" ± {sd:.2f}" if sd is not None and pd.notna(sd) else "")


lines = ["| | row | n | R² (vs no winner) | ρ | r | sign acc. |", "|---|---|---|---|---|---|---|"]
for r in tab.itertuples():
    lines.append(f"| {r.group} | {r.row}{' (`' + r.model + '`)' if r.model else ''} | "
                 f"{'' if pd.isna(r.n_test) else int(r.n_test)} | {fmt(r.r2, r.r2_sd)} | "
                 f"{fmt(r.spearman_rho, r.rho_sd)} | {fmt(r.pearson_r)} | {fmt(r.sign_accuracy)} |")
(OUT / "september_summary.md").write_text("\n".join(lines) + "\n")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(tab.drop(columns="note").round(3).to_string(index=False))
print("\npaired per fold, best genome model vs best taxonomy model\n" + paired.round(4).to_string(index=False))
print(f"\nwrote {OUT}")
