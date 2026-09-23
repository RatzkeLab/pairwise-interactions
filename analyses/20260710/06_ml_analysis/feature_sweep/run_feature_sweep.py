"""Does a richer annotation scheme beat KEGG KO -- on either target?

Every analysis so far used KEGG KO. This sweeps the annotation schemes that share the 298-strain
axis, on BOTH targets, with the same models and the same controls, so the results sit in one
comparable table instead of a series of one-off runs.

    KEGG_ko        7,164   the incumbent
    BiGG           9,442   metabolic reconstruction -- the only table that could support genuine
                           complementarity/competition rather than generic feature overlap, and
                           so the most likely to beat the taxonomy control on something other
                           than phylogeny
    CAZy              95   carbohydrate-active enzymes. Tiny, so almost no overfitting risk with
                           74-88 strains, and directly on point if carbon competition drives
                           these interactions
    KEGG_Module      622   curated pathway modules -- coarser than KO
    PFAMs          7,389   domain-level
    panX_full     58,424   pangenome gene clusters; finer-grained than KO
    panX_221_*       221   see below

**The two panX-221 rows are the point of including panX at all.** A collaborator's supervised
pipeline (gain screen -> redundancy removal -> dual-axis permutation selection -> Optuna -> SHAP)
selected 221 GOs on Or's carbon-source data. Two reasons not to treat that list as a clean
result here:

  1. It was selected using labels from experiments on **the same 298 strains**. Our honest
     metric holds strains out, so a feature set chosen because it predicts these particular
     strains' growth has already seen the held-out strains' phenotypes in a related assay.
     Growth and competitive ability are not independent, so that leaks into cv_strain.
  2. The collaborator's own note says it plainly: *"a random set of 221 genes does also a quite
     ok job, because these genes are basically markers for large scale patterns distributed
     across the whole genome."*

So `panX_221_random` (three independent draws from the same 58,424-column pool) is run beside
it. If the two match, the selection is capturing phylogeny, not mechanism -- which is what both
his note and our own 16S controls predict.

Read the table against two reference rows, never against zero:
  `ridge_phylo16s_only`  taxonomy alone -- what any feature set must beat to be about gene
                         content rather than the taxon label. It does not depend on the feature
                         table, so it is the same number in every block.
  `SHUFFLED_*`           the empirical zero.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE.parents[3] / "shared_pipelines"))

from config import CFG
import genomic_ml as gm
import genomic_ml_yield as gy

OUT = HERE / "outputs"
GEN = gm.GENOMIC_DIR
PANX = "panx_and_strains_table_final_post_filtering_strain_names.csv"
OD = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/Karl_2026/"
          "Karl_20260704_OD_Full")
SEL = json.loads((HERE.parents[0] / "go_feature_selection" /
                  "robust_feature_list.json").read_text())["go_cols"]

MODELS_RA = ["ridge_pca", "xgboost_pca", "two_stage_ridge", "ridge_phylo16s_only"]
MODELS_Y = ["ridge_pca", "xgboost_pca", "two_stage_ridge", "ridge_phylo16s_only"]
N_REPEATS = 2


def feature_sets():
    panx_cols = pd.read_csv(GEN / PANX, index_col=0, nrows=1).columns
    fs = {
        "KEGG_ko": ("KEGG_ko_and_strains_table.csv", None),
        "BiGG": ("BiGG_and_strains_table.csv", None),
        "CAZy": ("CAZy_and_strains_table.csv", None),
        "KEGG_Module": ("KEGG_Module_and_strains_table.csv", None),
        "PFAMs": ("PFAMs_and_strains_table.csv", None),
        "panX_full": (PANX, None),
        "panX_221_collaborator": (PANX, list(SEL)),
    }
    # A uniformly random draw is NOT a fair control: the selected 221 survive our prevalence
    # filter far better (139 vs ~52 columns), so an unmatched random set is handicapped on
    # feature count rather than on information. Draw prevalence-matched columns instead --
    # for each selected GO, a random GO from the same prevalence decile.
    panx = pd.read_csv(GEN / PANX, index_col=0)
    prev = (panx > 0).sum(axis=0)
    dec = pd.qcut(prev.rank(method="first"), 20, labels=False)
    by_dec = {d: prev.index[dec == d].tolist() for d in np.unique(dec)}
    sel_dec = dec.loc[[c for c in SEL if c in dec.index]]
    for i in range(3):
        rng = np.random.default_rng(900 + i)
        pick = []
        for d in sel_dec.values:
            pool = [c for c in by_dec[d] if c not in set(SEL)]
            if pool:
                pick.append(pool[int(rng.integers(len(pool)))])
        fs[f"panX_221_random_{i}"] = (PANX, list(dict.fromkeys(pick)))
        fs[f"panX_221_uniformrandom_{i}"] = (PANX, list(rng.choice(panx_cols, size=len(SEL),
                                                                   replace=False)))
    return fs


def run_relative_abundance(name, table, cols, phylo_cache):
    gcfg = gm.GenomicMLConfig(exp_cfg=CFG, out_dir=OUT / "_tmp_ra",
                              feature_table=table, feature_columns=cols)
    pairs, _ = gm.build_dataset(gcfg)
    X, summ = gm.strain_feature_matrix(gcfg, pairs)
    key = ("ra", tuple(sorted(set(pairs.strain_a) | set(pairs.strain_b))))
    phylo = phylo_cache.get(key)
    if phylo is None:
        phylo = gm.phylo_distance_matrix(gm.GenomicMLConfig(exp_cfg=CFG, out_dir=OUT),
                                         pairs, strict=False, out_name="s00_phylo_ra.csv")
        phylo_cache[key] = phylo
    s, _, _ = gm.cross_validate(gcfg, pairs, X, summ, models=MODELS_RA, phylo=phylo,
                                n_repeats=N_REPEATS, file_prefix="_tmp_ra")
    s["target"] = "relative_abundance"
    s["feature_set"] = name
    s["n_features"] = X.shape[1]
    s["n_pairs"] = len(pairs)
    s["n_strains"] = len(set(pairs.strain_a) | set(pairs.strain_b))
    return s


def run_yield(name, table, cols, phylo_cache):
    gcfg = gy.YieldMLConfig(exp_cfg=CFG, out_dir=OUT / "_tmp_y",
                            layout_csv=CFG.exp_base / "setup" / "strain_layout_20260630.csv",
                            od_dir=OD, feature_table=table, feature_columns=cols)
    pairs, _, _ = gy.build_dataset(gcfg)
    X, summ = gy.strain_feature_matrix(gcfg, pairs)
    key = ("y", tuple(sorted(set(pairs.strain_a) | set(pairs.strain_b))))
    phylo = phylo_cache.get(key)
    if phylo is None:
        phylo = gm.phylo_distance_matrix(gm.GenomicMLConfig(exp_cfg=CFG, out_dir=OUT),
                                         pairs, strict=False, out_name="s00_phylo_yield.csv")
        phylo_cache[key] = phylo
    s, _, _ = gy.cross_validate(gcfg, pairs, X, summ, models=MODELS_Y, phylo=phylo,
                                n_repeats=N_REPEATS)
    s["target"] = "yield_OD"
    s["feature_set"] = name
    s["n_features"] = X.shape[1]
    s["n_pairs"] = len(pairs)
    s["n_strains"] = len(set(pairs.strain_a) | set(pairs.strain_b))
    return s


def main():
    rows, cache = [], {}
    for name, (table, cols) in feature_sets().items():
        for label, fn in (("relative_abundance", run_relative_abundance), ("yield_OD", run_yield)):
            try:
                s = fn(name, table, cols, cache)
                rows.append(s)
                best = s[(s.regime == "cv_strain") & (s.model == "ridge_pca")]
                if len(best):
                    b = best.iloc[0]
                    print(f"  {name:24} {label:18} {int(s.n_features.iloc[0]):6d} feats  "
                          f"cv_strain ridge R2={b.r2_mean:+.3f} rho={b.spearman_rho_mean:.3f}")
            except Exception as e:
                print(f"  {name:24} {label:18} FAILED: {type(e).__name__}: {e}")
    df = pd.concat(rows, ignore_index=True)
    df.to_csv(OUT / "s01_feature_sweep_full.csv", index=False)

    key = df[(df.regime == "cv_strain") & (~df.model.str.startswith("SHUFFLED"))]
    piv = (key.pivot_table(index=["target", "feature_set", "n_features"], columns="model",
                           values="spearman_rho_mean")
              .round(3).reset_index())
    piv.to_csv(OUT / "s02_cv_strain_rho_by_feature_set.csv", index=False)
    print("\n=== cv_strain Spearman rho by feature set ===")
    print(piv.to_string(index=False))
    for f in OUT.glob("_tmp*"):
        import shutil
        shutil.rmtree(f) if f.is_dir() else f.unlink()
    return df


if __name__ == "__main__":
    main()
