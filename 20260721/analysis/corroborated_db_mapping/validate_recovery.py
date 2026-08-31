"""Is the recovered mapping actually right? The one non-circular test available.

`strain_identity_qc/qc_recovery.py` warns that validating a 16S-derived recovery by asking
whether 16S similarity tracks genome similarity is circular -- the recovery assigned genomes
*by* 16S, so that test answers itself. Its rho jumps from -0.05 to +0.91 and means nothing.

This test avoids that: the recovered mapping was built from **16S sequence only**, and is scored
against the **interaction labels**, which played no part in building it. If the recovery names
the right organisms, KO content should now predict how those organisms compete -- and it should
beat mappings that are permuted the same way as in genomic_ml_join_test/.

Controls are identical to that test: 5 permutations of the recovered label->genome assignment,
scored under cv_pair (which a consistent-but-wrong bijection can still pass) and cv_strain
(which it cannot).
"""

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXPS = HERE.parents[2]
sys.path.insert(0, str(EXPS / "shared_pipelines"))
sys.path.insert(0, str(EXPS / "20260721/analysis"))

from dataclasses import replace
import genomic_ml as gm
from config import CFG

OUT = HERE / "outputs"
REFIX = EXPS / "20260721/analysis/relative_abundance_refix/outputs"
MODEL = "xgboost_pca"
N_PERM = 5


def recovered_mapping():
    d = pd.read_csv(OUT / "m02_20260721_to_corroborated_db_mapping.csv")
    d = d.dropna(subset=["recommended_genome"])
    return dict(zip(d["label_20260721"], d["recommended_genome"].astype(str)))


def run():
    rec = recovered_mapping()
    print(f"recovered mapping: {len(rec)} labels -> {len(set(rec.values()))} genomes")

    real_loader = gm.load_strain_mapping
    dummy = pd.DataFrame({"Well_souce_plate": list(rec), "strain": list(rec.values())})

    cfg = replace(CFG, relative_abundance_out_dir=REFIX)
    rows = []
    for variant in ["recovered_mapping"] + [f"permuted_recovered_{i}" for i in range(N_PERM)]:
        if variant == "recovered_mapping":
            mp = dict(rec)
        else:
            rng = np.random.default_rng(500 + int(variant.split("_")[-1]))
            mp = dict(zip(rec.keys(), rng.permutation(list(rec.values()))))
        gm.load_strain_mapping = lambda g, _m=mp, _d=dummy: (_m, _d)
        try:
            gcfg = gm.GenomicMLConfig(exp_cfg=cfg, out_dir=OUT / "_recovery_tmp")
            pairs, _ = gm.build_dataset(gcfg)
            X, summ = gm.strain_feature_matrix(gcfg, pairs)
            s, _, _ = gm.cross_validate(gcfg, pairs, X, summ, models=[MODEL], n_repeats=3,
                                        shuffle_control=False, file_prefix="_rec")
        finally:
            gm.load_strain_mapping = real_loader
        for _, r in s.iterrows():
            rows.append({"variant": variant, "regime": r["regime"], "n_pairs": len(pairs),
                         "r2_mean": r["r2_mean"], "r2_sd": r["r2_sd"],
                         "spearman_rho_mean": r["spearman_rho_mean"],
                         "sign_accuracy_mean": r["sign_accuracy_mean"]})
        print(f"  {variant:22} " + "  ".join(
            f"{r['regime']}: R2={r['r2_mean']:+.3f} rho={r['spearman_rho_mean']:.3f} "
            f"sign={r['sign_accuracy_mean']:.3f}" for _, r in s.iterrows()))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "m04_recovery_validation.csv", index=False)
    ver = []
    for reg in ("cv_pair", "cv_strain"):
        d = df[df.regime == reg]
        t = d[d.variant == "recovered_mapping"]["r2_mean"].iloc[0]
        p = d[d.variant.str.startswith("permuted")]["r2_mean"].values
        ver.append({"regime": reg, "r2_recovered": t, "r2_permuted_mean": p.mean(),
                    "r2_permuted_sd": p.std(ddof=1),
                    "z": (t - p.mean()) / p.std(ddof=1) if p.std(ddof=1) > 0 else np.nan,
                    "beats_all_permutations": bool(t > p.max())})
    v = pd.DataFrame(ver)
    v.to_csv(OUT / "m05_recovery_verdict.csv", index=False)
    print("\n=== VERDICT ===")
    print(v.round(3).to_string(index=False))
    for f in OUT.glob("_rec*"):
        shutil.rmtree(f) if f.is_dir() else f.unlink()
    return v


if __name__ == "__main__":
    run()
