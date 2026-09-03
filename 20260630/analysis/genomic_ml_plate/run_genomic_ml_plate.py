"""Driver: does plate-reader optics add anything to the genome for 20260630's who-wins model?

The baseline arm (`T0_genomic`) is genomic_ml's design recomputed on these folds, so the
comparison is paired rather than a read-off against the published g03 table -- see
shared_pipelines/genomic_ml_plate.py for what each tier is allowed to see.

Usage:
    python run_genomic_ml_plate.py                 # full sweep (~ minutes)
    python run_genomic_ml_plate.py --quick         # 1 repeat, no shuffle control, no QC sweep
"""
import sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))                       # config.py
sys.path.insert(0, str(HERE.parents[3] / "shared_pipelines"))

import pandas as pd
from config import CFG
import genomic_ml as gm
import genomic_ml_plate as gp

pd.set_option("display.width", 220)
QUICK = "--quick" in sys.argv

PR = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/Karl_2026")

gml = gm.GenomicMLConfig(exp_cfg=CFG)
cfg = gp.PlateMLConfig(
    gml=gml,
    layout_csv=CFG.exp_base / "setup" / "strain_layout_20260630.csv",
    od_full_dir=PR / "Karl_20260704_OD_Full",     # 61-channel scans of the 30 destination plates
    source_od_dir=PR / "Karl_20260623_OD",        # the pre-experiment 384-collection reads
)

t0 = time.time()

# ---- the same labels the baseline uses, and the same strain-join gate ----------------
val, ok = gm.validate_strain_join(gml)
assert ok, "strain join failed validation -- the baseline itself would be invalid"
pairs0, dsum = gm.build_dataset(gml)
print("=== labels (from genomic_ml.build_dataset) ===")
print(dsum.to_string(index=False))

# ---- plate-reader features ----------------------------------------------------------
strains0 = sorted(set(pairs0.strain_a) | set(pairs0.strain_b))
pf = gp.build_plate_features(cfg, strains0)
print(f"\nplate reader: {pf['n_plates']} destination plates x {len(pf['wavelengths'])} channels "
      f"({pf['wavelengths'][0]:.0f}-{pf['wavelengths'][-1]:.0f} nm); "
      f"{len(pf['mono_scalar'])}/{len(strains0)} modeled strains have a mono well; "
      f"{pf['precult'].shape[1]} source-plate reads")
if pf["missing_mono"]:
    print("  strains with no mono well:", pf["missing_mono"])

pairs, psum = gp.attach_to_pairs(cfg, pairs0, pf)
print("\n=== p01 dataset ==="); print(psum.to_string(index=False))
qc = pf["qc"].reindex(sorted(set(pairs.strain_a) | set(pairs.strain_b)))
print("mono-well QC (mapping_validation):",
      ", ".join(f"{k}={v}" for k, v in qc.mono_qc_status.value_counts().items()),
      "-- only mono_low_confidence is evidence AGAINST a well; not_assessed just means too few reads")

X, summ = gm.strain_feature_matrix(gml, pairs)
print(f"KO features: {X.shape[0]} strains x {X.shape[1]} prevalence-filtered KOs")

# ---- what each feature says on its own, before any model ----------------------------
screen = gp.univariate_feature_screen(cfg, pairs)
print("\n=== p02 univariate screen ===")
print(screen.round(3).to_string(index=False))

# ---- the control that decides whether the source-plate feature may be believed ------
off = gp.source_plate_specificity_control(cfg, pairs)
print("\n=== p02b source-plate specificity control ===")
print(off[["source_read", "r_true", "perm_null_mean", "perm_null_sd", "z_vs_permutation",
           "p_perm", "layout_specific"]].round(3).to_string(index=False))
print("offset rows (read only where |grid_selfcorr| ~ 0 -- this collection has 6-column periodicity):")
cols = [c for c in off.columns if c.startswith(("r_offset", "selfcorr"))]
print(off[["source_read"] + cols].round(3).to_string(index=False))

band = gp.spectral_band_profile(cfg, pairs)
print("\n=== p02c spectral band profile (strongest bands) ===")
print(band.reindex(band.pearson_r.abs().sort_values(ascending=False).index)
      .head(8).round(3).to_string(index=False))
print("600 nm check (must be 0.0):",
      float(band.loc[band.wavelength_nm == 600, "pearson_r"].iloc[0]))
print(pd.read_csv(cfg.out_dir / "p02c_spectral_pc_structure.csv").round(4).to_string(index=False))

# ---- the tier sweep -----------------------------------------------------------------
summ_df, fold_df, pair_df = gp.cross_validate_tiers(
    cfg, pairs, X, summ, n_repeats=1 if QUICK else 3,
    plate_shuffle_control=not QUICK)
print("\n=== p03 cross-validation, all tiers, identical folds ===")
for regime in ("cv_strain", "cv_pair"):
    d = summ_df[summ_df.regime == regime]
    print(f"\n-- {regime} --")
    print(d[["tier", "model", "n_test", "r2_mean", "r2_sd", "spearman_rho_mean",
             "sign_accuracy_mean", "mae_mean"]].round(3).to_string(index=False))

# ---- the actual result: fold-paired deltas against the genomic baseline --------------
paired = gp.paired_tier_comparison(cfg, fold_df)
print("\n=== p04 paired vs T0_genomic (R², cv_strain) ===")
d = paired[(paired.metric == "r2") & (paired.regime == "cv_strain")]
print(d[["tier", "model", "baseline_mean", "tier_mean", "delta_mean", "delta_se",
         "wilcoxon_p"]].round(4).to_string(index=False))
print("\n=== p04 paired vs T0_genomic (sign accuracy, cv_strain) ===")
d = paired[(paired.metric == "sign_accuracy") & (paired.regime == "cv_strain")]
print(d[["tier", "model", "baseline_mean", "tier_mean", "delta_mean", "delta_se",
         "wilcoxon_p"]].round(4).to_string(index=False))

# ---- sensitivity: only strains whose mono well was confirmed by mapping_validation ---
if not QUICK:
    info, cmp_conf = gp.mono_qc_sensitivity(cfg, pairs, X, summ)
    print("\n=== p05 mono-well QC sensitivity ===")
    print(info.to_string(index=False))
    if len(cmp_conf):
        d = cmp_conf[cmp_conf.metric.isin(["r2", "sign_accuracy"])]
        print(d[["arm", "regime", "metric", "tier", "model", "baseline_mean", "tier_mean",
                 "delta_mean", "delta_se", "wilcoxon_p"]].round(4).to_string(index=False))

figs = gp.make_all_figures(cfg, summ_df, fold_df, pair_df, screen, paired, band=band)
print("\nfigures:", figs)
print(f"done in {time.time()-t0:.0f}s -> {cfg.out_dir}")
