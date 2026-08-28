"""Does the Well_souce_plate -> genome join carry ANY signal for 20260721?

Direct version of the question: forget structural arguments about 16S vs. gene content -- just
train XGBoost on 20260721's interaction labels with KO features attached through
`mapping_384_well_plate_collection.csv`, and see whether it beats chance.

The trap this script is built around
-----------------------------------
If the mapping is *consistently* wrong -- every well points at some genome, just the wrong one --
it is still a bijection. The KO vector then works as a unique per-strain fingerprint, and a model
that has seen a strain's other pairs during training can memorise "this fingerprint wins a lot"
and post a high score while knowing nothing about biology. So **`cv_pair` cannot tell a correct
mapping from a scrambled one**, and a good `cv_pair` number here would prove nothing.

Two things fix that:

1.  Score under `cv_strain` (both strains of a test pair never seen in training). A fingerprint
    is useless for a strain the model has never met; only a real genotype -> phenotype relation
    survives.
2.  Compare against a **permuted mapping**: shuffle which genome vector is attached to which
    well. This preserves the genomes, the labels, the pair structure and the
    fingerprint-as-ID property, and destroys *only* the well->genome correspondence. It is the
    exact null for "does this join mean anything", and is a much sharper control than shuffling
    the labels.

Both experiments are run: 20260630 is the positive control. Without it, a null result on
20260721 is unreadable -- it could just mean the whole approach is too weak to detect anything.

Usage:  python run_join_test.py [--quick]
"""

import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPS = HERE.parents[2]
sys.path.insert(0, str(EXPS / "shared_pipelines"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import genomic_ml as gm
from genomic_ml import COLOR_BLUE, COLOR_CRITICAL, COLOR_GRID, COLOR_GOOD, COLOR_TEXT_SECONDARY

OUT = HERE / "outputs"
FIG = OUT / "figures"
for d in (OUT, FIG):
    d.mkdir(parents=True, exist_ok=True)

N_PERMUTED_MAPPINGS = 5
MODEL = "xgboost_pca"


def load_experiment(exp_name):
    """Build the modeling table for one experiment via the Well_souce_plate join."""
    cfg_dir = EXPS / exp_name / "analysis"
    sys.path.insert(0, str(cfg_dir))
    for m in ("config",):
        sys.modules.pop(m, None)
    import config
    importlib_cfg = config.CFG
    sys.path.remove(str(cfg_dir))
    sys.modules.pop("config", None)

    gcfg = gm.GenomicMLConfig(exp_cfg=importlib_cfg, out_dir=OUT / f"_{exp_name}_dataset")
    pairs, summary = gm.build_dataset(gcfg)
    X, summ = gm.strain_feature_matrix(gcfg, pairs)
    return gcfg, pairs, X, summ, summary


def permute_mapping(X, summ, seed):
    """Reassign genome feature vectors to well labels at random.

    Same genomes, same labels, same pair structure -- only the correspondence is destroyed.
    """
    rng = np.random.default_rng(seed)
    perm = rng.permutation(X.index.values)
    Xp, sp = X.copy(), summ.copy()
    Xp.index = perm
    sp.index = perm
    return Xp.loc[X.index], sp.loc[summ.index]


def run(exp_name, n_repeats=3, quick=False):
    gcfg, pairs, X, summ, summary = load_experiment(exp_name)
    print(f"\n=== {exp_name}: {len(pairs)} pairs, {X.shape[0]} strains, {X.shape[1]} KOs")

    variants = [("true_mapping", X, summ, False)]
    for s in range(N_PERMUTED_MAPPINGS if not quick else 2):
        Xp, sp = permute_mapping(X, summ, seed=1234 + s)
        variants.append((f"permuted_mapping_{s}", Xp, sp, False))
    variants.append(("shuffled_labels", X, summ, True))

    rows = []
    for name, Xv, sv, shuffle_y in variants:
        pv = pairs
        if shuffle_y:
            rng = np.random.default_rng(7)
            pv = pairs.copy()
            pv[gcfg.target] = rng.permutation(pairs[gcfg.target].values)
        summ_df, fold_df, _ = gm.cross_validate(
            gcfg, pv, Xv, sv, models=[MODEL], n_repeats=n_repeats,
            shuffle_control=False, file_prefix=f"_tmp_{exp_name}_{name}")
        for _, r in summ_df.iterrows():
            rows.append({"experiment": exp_name, "variant": name, "regime": r["regime"],
                         "r2_mean": r["r2_mean"], "r2_sd": r["r2_sd"],
                         "spearman_rho_mean": r["spearman_rho_mean"],
                         "sign_accuracy_mean": r["sign_accuracy_mean"],
                         "n_test": r["n_test"]})
        print(f"  {name:22} " + "  ".join(
            f"{r['regime']}: R2={r['r2_mean']:+.3f} rho={r['spearman_rho_mean']:.3f} "
            f"sign={r['sign_accuracy_mean']:.3f}" for _, r in summ_df.iterrows()))
    for f in OUT.glob("_tmp_*"):
        f.unlink()
    return pd.DataFrame(rows), summary


def verdict(df, exp_name, regime="cv_strain"):
    """Is the true mapping distinguishable from a random one, under this regime?"""
    d = df[(df["experiment"] == exp_name) & (df["regime"] == regime)]
    true = d[d["variant"] == "true_mapping"]["r2_mean"].iloc[0]
    perm = d[d["variant"].str.startswith("permuted_mapping")]["r2_mean"].values
    z = (true - perm.mean()) / perm.std(ddof=1) if perm.std(ddof=1) > 0 else np.nan
    return {"experiment": exp_name, "regime": regime, "r2_true_mapping": true,
            "r2_permuted_mean": perm.mean(), "r2_permuted_sd": perm.std(ddof=1),
            "n_permutations": len(perm), "z_true_vs_permuted": z,
            "beats_random_mapping": bool(true > perm.max())}


def make_figure(df, verdicts):
    exps = list(dict.fromkeys(df["experiment"]))
    fig, axes = plt.subplots(2, len(exps), figsize=(6.2 * len(exps), 8.4), squeeze=False)
    for j, exp in enumerate(exps):
        for i, regime in enumerate(["cv_pair", "cv_strain"]):
            ax = axes[i][j]
            d = df[(df["experiment"] == exp) & (df["regime"] == regime)]
            perm = d[d["variant"].str.startswith("permuted_mapping")]
            true = d[d["variant"] == "true_mapping"]
            shuf = d[d["variant"] == "shuffled_labels"]
            ax.bar(["true\nmapping"], true["r2_mean"], yerr=true["r2_sd"], color=COLOR_BLUE,
                   error_kw=dict(ecolor=COLOR_TEXT_SECONDARY, lw=1, capsize=4), width=0.55)
            ax.bar([f"permuted\nmapping\n(n={len(perm)})"], [perm["r2_mean"].mean()],
                   yerr=[perm["r2_mean"].std(ddof=1)], color=COLOR_TEXT_SECONDARY,
                   error_kw=dict(ecolor="black", lw=1, capsize=4), width=0.55)
            ax.scatter([1] * len(perm), perm["r2_mean"], color="white", edgecolor="black",
                       zorder=5, s=26, linewidth=0.8)
            ax.bar(["shuffled\nlabels"], shuf["r2_mean"], yerr=shuf["r2_sd"], color=COLOR_CRITICAL,
                   error_kw=dict(ecolor=COLOR_TEXT_SECONDARY, lw=1, capsize=4), width=0.55)
            ax.axhline(0, color="black", lw=1)
            v = [x for x in verdicts if x["experiment"] == exp and x["regime"] == regime]
            tag = ""
            if v:
                tag = ("  ✓ beats random mapping" if v[0]["beats_random_mapping"]
                       else "  ✗ NOT better than random mapping")
            ax.set_title(f"{exp} — {regime}{tag}", fontsize=10,
                         color=COLOR_GOOD if "✓" in tag else (COLOR_CRITICAL if "✗" in tag else "black"))
            ax.set_ylabel("R² vs. predicting no winner")
    fig.suptitle("Does the Well_souce_plate → genome join carry real signal?\n"
                 "cv_pair cannot tell a correct mapping from a consistently scrambled one — "
                 "read the bottom row", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "j01_true_vs_permuted_mapping.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    quick = "--quick" in sys.argv
    t0 = time.time()
    frames, summaries = [], []
    for exp in ("20260721", "20260630"):
        f, s = run(exp, n_repeats=1 if quick else 3, quick=quick)
        frames.append(f)
        summaries.append(s.assign(experiment=exp))
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(OUT / "j01_mapping_permutation_results.csv", index=False)
    pd.concat(summaries, ignore_index=True).to_csv(OUT / "j00_dataset_summaries.csv", index=False)

    verdicts = [verdict(df, e, r) for e in ("20260721", "20260630")
                for r in ("cv_pair", "cv_strain")]
    vdf = pd.DataFrame(verdicts)
    vdf.to_csv(OUT / "j02_verdicts.csv", index=False)
    print("\n=== VERDICTS ===")
    print(vdf.round(3).to_string(index=False))
    make_figure(df, verdicts)
    print(f"\ndone in {time.time()-t0:.0f}s -> {OUT}")
