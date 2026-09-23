"""Do the known caveats change the answer? One refit per caveat, same folds throughout.

The first pass reported cv_strain R^2 = 0.372 with none of these applied. Each arm below
changes exactly one thing and is compared against that same baseline, so a difference is
attributable. A caveat that moves nothing is a caveat that can stop being mentioned.

The arms:

  A  baseline                 no corrections -- the first-pass configuration
  B  density-corrected y      see below; this is NOT the correction I first proposed
  C  + preculture density     the same information as a FEATURE instead of a correction
  D  drop non-growing strains the 30 strains that failed both monocultures
  E  + mono OD                per-strain monoculture yield, which is known to carry signal

**On B, and why the obvious version of it is a no-op.** The ratio titration measured a
slope of +0.31 of final log-ratio on *initial inoculum log-ratio*. The natural correction
is to subtract slope x initial_ratio from the outcome -- but every pair in the modelling
set comes from a well at 1:1 by design, where initial_ratio is exactly 0, so that
subtraction does nothing at all.

The confound that actually survives into the modelling set is different: equal VOLUMES are
not equal CELL NUMBERS. Source wells differ several-fold in preculture density, so a
nominally 1:1 well starts skewed by however much the two strains' precultures differed.
That is measurable -- it is the source-plate OD - so the effective initial log-ratio is
log2(OD_source_a / OD_source_b), and B subtracts slope x that. Arm C then asks whether the
model can use the same quantity better as a feature than as a fixed correction.
"""
import sys
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE)); sys.path.insert(0, str(BASE.parents[1] / "scripts"))

import numpy as np, pandas as pd
import config, genomic_ml as gml, plate_reader as pr

OUT = BASE / "02_sensitivity" / "outputs"; OUT.mkdir(parents=True, exist_ok=True)
ROWS16 = "ABCDEFGHIJKLMNOP"
TITRATION_SLOPE = 0.31          # from analyses/20260921/06_inoculum_arms
MODELS = ["ridge_pca", "two_stage_ridge"]
N_REPEATS = 3


def source_od():
    """Per-strain preculture density, median over the source-plate reads.

    The collection plate is a 384 whose well coordinate IS the strain name, so the grid
    maps straight onto strain labels. Several reads exist at different times; the median
    is used because individual wells saturate at the reader's ceiling of 3.5.
    """
    files = ["Karl_20260908_115807_ODFull_10.csv", "Karl_20260908_123137_ODFull_10.csv",
             "Karl_20260909_121752_ODFull_10.csv", "Karl_20260909_123000_ODFull_10.csv"]
    grids = [pr.read_plate(config.OD_PREP / f, 600)[0] for f in files]
    grids = [g for g in grids if g is not None and g.shape == (16, 24)]
    med = np.median(np.dstack(grids), axis=2)
    return pd.Series({f"{ROWS16[r]}{c+1}": float(med[r, c])
                      for r in range(16) for c in range(24)}, name="source_od")


def mono_od():
    """Per-strain monoculture yield on the destination plates (assay medium, endpoint)."""
    od, _ = pr.load_folder(config.OD_FULL, testname="ODFull")
    lay = pd.read_csv(config.LAYOUT_CSV)
    lay = lay[(lay.plate_role == "primary_sequenced") & (lay.well_type == "mono")]
    m = lay.merge(od[["dest_plate", "dest_well", "od"]], on=["dest_plate", "dest_well"])
    return m.groupby("strain1").od.median().rename("mono_od")


def run(gcfg, pairs, X, summ, label, phylo=None):
    s, f, _ = gml.cross_validate(gcfg, pairs, X, summ, models=MODELS,
                                 n_repeats=N_REPEATS, phylo=phylo,
                                 shuffle_control=False, file_prefix=f"s02_{label}")
    s = s[s.regime == "cv_strain"].copy(); s["arm"] = label
    return s


def main():
    cfg = config.make_config()
    gcfg = gml.GenomicMLConfig(exp_cfg=cfg, out_dir=OUT)
    pairs, _ = gml.build_dataset(gcfg)
    X, summ = gml.strain_feature_matrix(gcfg, pairs)
    phylo = gml.phylo_distance_matrix(gcfg, pairs)
    print(f"{len(pairs)} pairs, {X.shape[0]} strains\n")

    src, mono = source_od(), mono_od()
    strains = list(X.index)
    src_s = src.reindex(strains); mono_s = mono.reindex(strains)
    print(f"preculture OD available for {src_s.notna().sum()}/{len(strains)} strains "
          f"(median {src_s.median():.2f})")
    print(f"mono OD available for       {mono_s.notna().sum()}/{len(strains)} strains "
          f"(median {mono_s.median():.2f})")

    # effective starting imbalance of a nominally 1:1 well
    a, b = pairs.strain_a.map(src).values, pairs.strain_b.map(src).values
    eff = np.log2(np.clip(a, 1e-3, None) / np.clip(b, 1e-3, None))
    eff = np.nan_to_num(eff, nan=0.0)
    print(f"\neffective initial log2 ratio from preculture density:")
    print(f"  median |imbalance| {np.median(np.abs(eff)):.2f}  "
          f"90th pct {np.percentile(np.abs(eff), 90):.2f}  "
          f"(a nominally 1:1 well is not 1:1 in cells)")
    print(f"  correlation with the observed outcome: "
          f"r={np.corrcoef(eff, pairs[gcfg.target])[0,1]:+.3f}")

    arms = []
    arms.append(run(gcfg, pairs, X, summ, "A_baseline", phylo))

    p_b = pairs.copy()
    p_b[gcfg.target] = p_b[gcfg.target] - TITRATION_SLOPE * eff
    arms.append(run(gcfg, p_b, X, summ, "B_density_corrected_y", phylo))

    summ_c = summ.copy(); summ_c["source_od"] = np.log2(src_s.clip(lower=1e-3)).fillna(0).values
    arms.append(run(gcfg, pairs, X, summ_c, "C_density_as_feature", phylo))

    dead = set(pd.read_csv(config.STRAINS_NO_GROWTH).strain)
    keep = ~(pairs.strain_a.isin(dead) | pairs.strain_b.isin(dead))
    p_d = pairs[keep].reset_index(drop=True)
    print(f"\nD: dropping non-growing strains removes {(~keep).sum()} of {len(pairs)} pairs")
    arms.append(run(gcfg, p_d, X, summ, "D_drop_nongrowing", phylo))

    summ_e = summ.copy(); summ_e["mono_od"] = mono_s.fillna(mono_s.median()).values
    arms.append(run(gcfg, pairs, X, summ_e, "E_mono_od_feature", phylo))

    res = pd.concat(arms, ignore_index=True)
    piv = res.pivot(index="arm", columns="model", values="r2_mean").round(3)
    sd = res.pivot(index="arm", columns="model", values="r2_sd").round(3)
    print("\n" + "=" * 78)
    print("cv_strain R^2 by arm (mean over folds)")
    print("=" * 78)
    print(piv.to_string())
    print("\nfold-to-fold sd:")
    print(sd.to_string())
    base = piv.loc["A_baseline"]
    print("\nchange vs baseline:")
    print((piv - base).round(3).drop(index="A_baseline").to_string())
    res.to_csv(OUT / "s02_sensitivity_summary.csv", index=False)
    print(f"\nwrote -> {OUT}")
    return res


if __name__ == "__main__":
    main()
