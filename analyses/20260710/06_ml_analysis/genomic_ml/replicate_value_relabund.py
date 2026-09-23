"""More replicates, or more reads? The same question as the yield analysis, for relative abundance.

The two targets have different bottlenecks, so the answer is not transferable:

  yield (OD)          every pair is measured, 2-3 replicate wells each, reliability 0.85 single
                      / 0.94 averaged. Replication is already saturated.
  relative abundance  only 413 pairs have TWO sequenced wells and NONE has three -- most wells
                      never clear the read threshold. So a replicate curve barely exists here,
                      and the interesting axis is READ DEPTH, which the yield target does not
                      have at all.

Four analyses:

  1. Single-well reliability + Spearman-Brown, as for yield. Bounded at k=2 by the data.
  2. Replicate learning curve k=1 vs k=2, on the 2-well subset. Small (256 usable pairs) and
     reported with that caveat rather than as a headline.
  3. **Read-depth curve by downsampling.** Reads are re-drawn from r02_read_assignments at
     25/50/75/100% and the label recomputed from the resulting class counts. This is the
     honest way to ask whether depth is the limit: it measures the actual shape of the curve
     rather than assuming a binomial model. Downsampling only probes BELOW current depth, so
     what it establishes is whether we are on the steep part or the plateau.
  4. Strain learning curve, for comparability with the yield analysis.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE.parents[3] / "shared_pipelines"))
from config import CFG
import genomic_ml as gm

OUT = HERE / "outputs"
N_SPLITS, N_REPEATS = 5, 3
MODEL = "two_stage_ridge"      # best cv_strain model on this target for 20260630


def cv_strain_score(gcfg, pairs, X, summ, seed=0, n_repeats=N_REPEATS):
    out = []
    for rep in range(n_repeats):
        for tr, te in gm.make_folds(pairs, "cv_strain", n_splits=N_SPLITS, seed=seed + 97 * rep):
            if len(te) < 10 or len(tr) < 50:
                continue
            try:
                ctx = gm._Ctx(gcfg, pairs, X, summ, tr, te)
                pred = gm.MODELS[MODEL](ctx)
            except Exception:
                continue
            y = ctx.y_test
            if np.sum(y ** 2) > 0 and np.std(pred) > 0:
                out.append({"r2": 1 - np.sum((y - pred) ** 2) / np.sum(y ** 2),
                            "rho": spearmanr(y, pred)[0]})
    if not out:
        return np.nan, np.nan
    d = pd.DataFrame(out)
    return d.r2.mean(), d.rho.mean()


def _pairs_from_labels(base, label_by_pair, gcfg):
    p = base.copy()
    key = [frozenset((a, b)) for a, b in zip(p.strain_a, p.strain_b)]
    p[gcfg.target] = [label_by_pair.get(k, np.nan) for k in key]
    return p.dropna(subset=[gcfg.target]).reset_index(drop=True)


def main():
    gcfg = gm.GenomicMLConfig(exp_cfg=CFG, out_dir=OUT)
    pairs, _ = gm.build_dataset(gcfg)
    X, summ = gm.strain_feature_matrix(gcfg, pairs)
    print(f"baseline: {len(pairs)} pairs, {len(set(pairs.strain_a)|set(pairs.strain_b))} strains")
    base_r2, base_rho = cv_strain_score(gcfg, pairs, X, summ)
    print(f"  cv_strain ({MODEL}): R2={base_r2:+.3f} rho={base_rho:.3f}")

    w = pd.read_csv(CFG.relative_abundance_out_dir / "r02_well_interaction_scores.csv")
    w["pk"] = [frozenset((a, b)) for a, b in zip(w.strain_a, w.strain_b)]
    keep = set(frozenset((a, b)) for a, b in zip(pairs.strain_a, pairs.strain_b))
    w = w[w.pk.isin(keep)]

    # ---- 1. single-well reliability, 1-vs-1 on exactly-two-well pairs
    per = w.groupby("pk")["mean_log2_ratio" if "mean_log2_ratio" in w.columns
                          else "log2_ratio_strain1_over_strain2"].apply(list)
    two = per[per.map(len) == 2]
    r1 = pearsonr([v[0] for v in two], [v[1] for v in two])[0]
    rel = pd.DataFrame([{"k_replicates": k, "reliability_of_mean": k * r1 / (1 + (k - 1) * r1)}
                        for k in [1, 2, 3, 4, 6, 10]] + [{"k_replicates": np.inf,
                                                          "reliability_of_mean": 1.0}])
    print(f"\nsingle-well reliability r1 = {r1:.3f} (1-vs-1 on {len(two)} two-well pairs)")
    print(rel.round(4).to_string(index=False))
    print(f"  attenuation-corrected ceiling: with PERFECT labels cv_strain R2 would be "
          f"~{base_r2 / (2*r1/(1+r1)):.3f} (currently {base_r2:+.3f})")

    # ---- 2. replicate curve, k=1 vs k=2, on the 2-well subset only
    print(f"\nreplicate curve on the {len(two)} pairs with two sequenced wells "
          f"(small -- treat as indicative)")
    rc = []
    for k in (1, 2):
        r2s, rhos = [], []
        for s in range(3):
            rng = np.random.default_rng(300 + s)
            lab = {key: float(np.mean(rng.choice(v, size=k, replace=False)))
                   for key, v in two.items()}
            pk = _pairs_from_labels(pairs, lab, gcfg)
            if len(pk) < 100:
                continue
            r2, rho = cv_strain_score(gcfg, pk, X, summ, seed=s, n_repeats=2)
            r2s.append(r2); rhos.append(rho)
        if r2s:
            rc.append({"k_replicates": k, "cv_strain_r2": np.nanmean(r2s),
                       "cv_strain_rho": np.nanmean(rhos), "n_pairs": len(two)})
            print(f"  k={k}: R2={np.nanmean(r2s):+.3f} rho={np.nanmean(rhos):.3f}")
    rc = pd.DataFrame(rc)

    # ---- 3. read-depth curve by downsampling reads
    print("\nread-depth curve (downsampling reads, label recomputed from class counts)")
    ra = pd.read_csv(CFG.relative_abundance_out_dir / "r02_read_assignments.csv.gz")
    ra = ra.merge(w[["sample_id", "strain_a", "strain_b", "strain1", "pk"]], on="sample_id")
    dc = []
    for frac in (0.25, 0.5, 0.75, 1.0):
        r2s, rhos, depth = [], [], []
        for s in range(2):
            rng = np.random.default_rng(400 + s)
            sub = (ra.sample(frac=frac, random_state=int(rng.integers(1e6)))
                   if frac < 1.0 else ra)
            cnt = (sub[sub.read_class.isin(["strain1", "strain2"])]
                   .groupby(["sample_id", "read_class"]).size().unstack(fill_value=0))
            if "strain1" not in cnt or "strain2" not in cnt:
                continue
            tot = cnt.sum(axis=1)
            cnt = cnt[tot >= 5]
            p1 = cnt["strain1"] / (cnt["strain1"] + cnt["strain2"])
            meta = w.set_index("sample_id").loc[cnt.index]
            # orient to strain_a, then log2 ratio, matching the pipeline's target
            p_a = np.where(meta.strain1.values == meta.strain_a.values, p1.values, 1 - p1.values)
            p_a = np.clip(p_a, 1e-3, 1 - 1e-3)
            lab = pd.Series(np.log2(p_a / (1 - p_a)), index=meta.pk.values).groupby(level=0).mean()
            pk = _pairs_from_labels(pairs, lab.to_dict(), gcfg)
            r2, rho = cv_strain_score(gcfg, pk, X, summ, seed=s, n_repeats=2)
            r2s.append(r2); rhos.append(rho); depth.append(float(tot.median()))
        if r2s:
            dc.append({"read_fraction": frac, "median_reads_per_well": np.mean(depth),
                       "cv_strain_r2": np.nanmean(r2s), "cv_strain_rho": np.nanmean(rhos)})
            print(f"  {int(frac*100):3d}% of reads (median {np.mean(depth):.0f}/well): "
                  f"R2={np.nanmean(r2s):+.3f} rho={np.nanmean(rhos):.3f}")
    dc = pd.DataFrame(dc)

    # ---- 4. strain learning curve
    print("\nstrain learning curve")
    strains = sorted(set(pairs.strain_a) | set(pairs.strain_b))
    sc = []
    for S in (30, 45, 60, len(strains)):
        r2s, rhos, npair = [], [], []
        for s in range(3):
            rng = np.random.default_rng(500 + s)
            kp = set(rng.choice(strains, size=min(S, len(strains)), replace=False))
            sub = pairs[pairs.strain_a.isin(kp) & pairs.strain_b.isin(kp)].reset_index(drop=True)
            if len(sub) < 120:
                continue
            r2, rho = cv_strain_score(gcfg, sub, X, summ, seed=s, n_repeats=2)
            r2s.append(r2 if np.isfinite(r2) and r2 > -5 else np.nan)
            rhos.append(rho); npair.append(len(sub))
        if rhos:
            sc.append({"n_strains": S, "cv_strain_r2": np.nanmean(r2s),
                       "cv_strain_rho": np.nanmean(rhos), "mean_n_pairs": np.mean(npair)})
            print(f"  S={S:3d} ({np.mean(npair):.0f} pairs): R2={np.nanmean(r2s):+.3f} "
                  f"rho={np.nanmean(rhos):.3f}")
    sc = pd.DataFrame(sc)

    rel.to_csv(OUT / "g06_reliability_spearman_brown.csv", index=False)
    rc.to_csv(OUT / "g06_replicate_curve.csv", index=False)
    dc.to_csv(OUT / "g06_read_depth_curve.csv", index=False)
    sc.to_csv(OUT / "g06_strain_curve.csv", index=False)
    return rel, rc, dc, sc


if __name__ == "__main__":
    main()
