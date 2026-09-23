"""First-pass genomic ML on 20260907: predict pairwise outcome from gene content.

Deliberately bare-bones. No slope correction, no plate-reader features, no dropping of
the strains that grew poorly -- those are the next pass, and the point of doing this one
first is to have something to compare them against.

Two things are not optional and are done here:

  the strain-ID join gate   `validate_strain_join` asks whether joining this experiment's
                            well labels to the genomic tables gives the RIGHT genomes,
                            using a fact the labels cannot fake: two strains with
                            near-identical 16S must have near-identical gene content. Only
                            20260630 has ever passed it; 20260721 failed decisively. If
                            20260907 fails, nothing below is worth reading and the script
                            stops.

  cv_strain, not cv_pair    the pairs are combinations of a few hundred genomes, not
                            independent samples, so a model can score well on held-out
                            PAIRS by memorising a strain's competitiveness from its other
                            pairs. Held-out STRAINS is the honest number; cv_pair is
                            reported alongside only to show the size of that gap.

Labels come from wells at the standard 1:1 inoculum only -- the titration wells are a
different condition, not replicates (see relative_abundance.replicate_stability).
"""
import sys
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE)); sys.path.insert(0, str(BASE.parent / "shared_scripts"))

import numpy as np, pandas as pd
import config, relative_abundance as ra, genomic_ml as gml

OUT = BASE / "07_genomic_ml" / "outputs"; OUT.mkdir(parents=True, exist_ok=True)


def main():
    cfg = config.make_config()
    cfg.ra_reference_fasta = (BASE / "03_recreate_reference_db" / "outputs" /
                              "c02_strain_consensus_corrected.fasta")

    print("=" * 74); print("labels: aggregate each pair's 1:1 wells"); print("=" * 74)
    ra.replicate_stability(cfg, balanced_inoculum_only=True)

    gcfg = gml.GenomicMLConfig(exp_cfg=cfg, out_dir=OUT)

    print("\n" + "=" * 74)
    print("GATE: does the well -> genome join give the right genomes?")
    print("=" * 74)
    val, passed = gml.validate_strain_join(gcfg)
    print(val.to_string(index=False))
    print(f"\n  verdict: {'PASS' if passed else 'FAIL'}")
    if not passed:
        print("  20260907 cannot be joined to the genomic tables on these labels. "
              "Stopping:\n  every number below would be a model of a random "
              "genome-to-well assignment.")
        return None
    print("  (20260630 passed at rho=+0.363, z=+6.0; 20260721 failed at rho=-0.046, z=-0.7)")

    print("\n" + "=" * 74); print("dataset"); print("=" * 74)
    pairs, summary = gml.build_dataset(gcfg)
    print(summary.to_string(index=False))

    X, summ = gml.strain_feature_matrix(gcfg, pairs)
    print(f"\nfeature matrix: {X.shape[0]} strains x {X.shape[1]} features")

    print("\n" + "=" * 74)
    print("cross-validation -- cv_strain is the result, cv_pair is the comparison")
    print("=" * 74)
    # bare-bones: the two baselines that say what "good" even means, plus one real model.
    #   zero_baseline                  predicts 0 for everything -- the floor
    #   strength_observed_no_genomics  fits a per-strain strength from the TRAINING pairs and
    #                                  uses no genome at all. Under cv_strain a held-out strain
    #                                  has no fitted strength, so this collapses -- which is
    #                                  exactly the point: it shows how much of cv_pair is
    #                                  memorised competitiveness rather than genomics.
    #   ridge_pca                      gene content -> outcome, the actual question
    models = ["zero_baseline", "strength_observed_no_genomics", "ridge_pca"]
    per_model, folds, preds = gml.cross_validate(gcfg, pairs, X, summ, models=models,
                                                 n_repeats=3)
    cols = ["regime", "model", "n_folds", "n_test", "r2_mean", "r2_sd",
            "spearman_rho_mean", "sign_accuracy_mean"]
    print(per_model[cols].to_string(index=False))

    ceiling = gml.label_noise_ceiling(gcfg, pairs)
    print("\nlabel-noise ceiling (how well anything could do given replicate disagreement):")
    print(ceiling.to_string(index=False))

    return pairs, X, summ, folds, per_model


if __name__ == "__main__":
    main()
