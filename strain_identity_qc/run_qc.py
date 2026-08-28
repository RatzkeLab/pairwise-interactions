"""Run the whole strain-identity QC end to end.

    python run_qc.py            # everything (~20 min, dominated by all-vs-all 16S alignment)
    python run_qc.py --figs     # figures only, from existing outputs

Writes to outputs/<experiment>/ and outputs/figures/.
"""

import sys
import time

import pandas as pd

import qc_config as C
import qc_sources as S
import qc_compare as Q
import qc_recovery as R
import qc_layout as L
import qc_readassign as RA
import qc_figures as F

pd.set_option("display.width", 240)


def main(figs_only=False):
    t0 = time.time()
    if not figs_only:
        att_all = []
        for name, exp in C.EXPERIMENTS.items():
            print(f"\n{'='*70}\n{name}\n{'='*70}")
            df = S.load_all_sources(exp)
            df.to_csv(exp.out / f"s00_all_16S_sources_{name}.csv", index=False)
            print(df.groupby("source").agg(n_seqs=("seq", "size"),
                                           n_strains=("strain_label", "nunique"),
                                           median_len=("length", "median")).to_string())

            ag = Q.agreement_table(df)
            ag["experiment"] = name
            ag.to_csv(exp.out / f"s01_source_agreement_{name}.csv", index=False)

            for tgt in ("corroborated_db", "genome_16S"):
                att = Q.attribution_table(df, "pair_consensus", tgt)
                att.to_csv(exp.out / f"s02_attribution_pair_vs_{tgt}_{name}.csv", index=False)
                att = att.copy(); att["query_source"] = f"pair@{name}"
                att_all.append(att)

            # three-way: with two independent references, blame can be assigned
            reps = {s: Q._rep_seqs(df, s) for s in ("pair_consensus", "corroborated_db", "genome_16S")}
            common = sorted(set(reps["pair_consensus"]) & set(reps["corroborated_db"]) & set(reps["genome_16S"]))
            rows = [{"strain_label": st,
                     "pair_vs_corrob": Q.best_identity(reps["pair_consensus"][st], reps["corroborated_db"][st])[0],
                     "pair_vs_genome": Q.best_identity(reps["pair_consensus"][st], reps["genome_16S"][st])[0],
                     "corrob_vs_genome": Q.best_identity(reps["corroborated_db"][st], reps["genome_16S"][st])[0]}
                    for st in common]
            tw = pd.DataFrame(rows)
            tw["odd_one_out"] = tw.apply(_odd_one_out, axis=1)
            tw.to_csv(exp.out / f"s04_three_way_concordance_{name}.csv", index=False)
            print(f"\nthree-way concordance on {len(tw)} strains present in all three sources:")
            print(" ", tw["odd_one_out"].value_counts().to_dict())

        summ = Q.summarize_attribution(pd.concat(att_all, ignore_index=True))
        summ.to_csv(C.OUT / "s03_attribution_summary.csv", index=False)
        print(f"\n{'='*70}\nattribution summary\n{'='*70}")
        print(summ.to_string(index=False))

        # the rescue that does not work, kept on record
        exp = C.EXPERIMENTS["20260721"]
        mp, conf = R.recovered_map(exp)
        print("\n16S-based re-identification of 20260721 -- resolving power:")
        print(R.resolution_report(exp, conf).to_string(index=False))
        w2s = dict(zip(*R._mapping_cols()))
        circ = pd.DataFrame([R.circularity_demo(exp, w2s, "original_labels"),
                             R.circularity_demo(exp, mp, "16S_recovered")])
        circ.to_csv(exp.out / f"s06_recovery_circularity_{exp.name}.csv", index=False)
        print("\ncircularity demonstration (the recovered rho is an artifact, not a result):")
        print(circ.to_string(index=False))

        # can a named plate-handling mistake explain 20260721? (20260630 is the positive control)
        print(f"\n{'='*70}\nplate-layout transform search\n{'='*70}")
        res = L.resolvable_strains()
        print(f"resolvable strains (no 16S near-twin): {int(res['resolvable'].sum())}/{len(res)}")
        group_of, sizes = L.sixteen_s_groups()
        print(f"16S groups: {len(sizes)}; largest holds {sizes.max()} of {len(group_of)} strains")
        for name, exp in C.EXPERIMENTS.items():
            r = L.group_level_transform_search(exp, group_of)
            r.to_csv(exp.out / f"s09_layout_transform_search_{name}.csv", index=False)
            b = r.iloc[0]
            print(f"\n{name}: {int(b['n_transforms_tested'])} transforms tested, "
                  f"family-wise z threshold (p=0.05) = {b['familywise_z_threshold_p05']}")
            print(f"  best = {b['transform']}  z = {b['z_vs_null']}  "
                  f"family-wise p = {b['familywise_p_for_best']}")

        # calibrate the read-assignment resolution limit against mono-well ground truth
        print(f"\n{'='*70}\nread-assignment calibration (mono-well ground truth)\n{'='*70}")
        e630 = C.EXPERIMENTS["20260630"]
        pp, pr = RA.mono_ground_truth(
            e630, C.ROOT / "20260630/analysis/consensus2/strain_consensus_20260630.fasta")
        print(RA.read_level_curve(pr).to_string(index=False))

    print("\nfigures:", F.make_all())
    print(f"done in {time.time()-t0:.0f}s -> {C.OUT}")


def _odd_one_out(r):
    T = Q.SAME_ORGANISM_IDENTITY
    pc, pg, cg = r.pair_vs_corrob >= T, r.pair_vs_genome >= T, r.corrob_vs_genome >= T
    if pc and pg and cg:
        return "all_agree"
    if pc and not pg and not cg:
        return "genome_16S_is_odd"
    if pg and not pc and not cg:
        return "corroborated_is_odd"
    if cg and not pc and not pg:
        return "experiment_is_odd"
    return "no_two_agree"


if __name__ == "__main__":
    main(figs_only="--figs" in sys.argv)
