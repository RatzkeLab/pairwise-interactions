"""Can 20260721's strains be re-identified from what its wells actually contain?

The QC establishes that 20260721's wells hold real, clean collection organisms that are simply
not the ones its labels claim. That invites an obvious rescue: ignore the labels, identify each
well by its 16S, and re-join to the genomic tables on the observed identity instead.

**This module exists to show that the rescue does not work, and to keep the reason on record.**

Run naively it looks like a triumph -- re-joining this way takes the genome-vs-16S validation
correlation from rho=-0.05 to rho=+0.91. That number is an artifact. The recovery assigns
genomes *by* 16S similarity, and the validation then asks whether 16S similarity tracks genome
similarity, so it is asking a question it has already answered. `resolution_report()` measures
the thing that actually matters instead: whether 16S has the resolving power to make the call
at all. It does not -- 75 labels collapse onto 39 genomes, one genome is claimed by 10
different labels, and the median call beats its runner-up by ~2.7 bp in 1420.

So the honest conclusion is that 20260721 is not recoverable from 16S. It would take a marker
with more resolution than 16S (shotgun of the source wells), or physically identifying the
plate that was actually used.
"""

import collections

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import qc_config as C
import qc_sources as S


def recovered_map(exp, target_source="genome_16S", min_identity=0.99):
    """{experiment strain label -> genomic strain name} implied by observed 16S, ignoring labels."""
    att = pd.read_csv(exp.out / f"s02_attribution_pair_vs_{target_source}_{exp.name}.csv")
    conf = att[att["best_identity"] >= min_identity].copy()
    w2s = dict(zip(*_mapping_cols()))
    conf["genome"] = conf["best_match"].map(w2s)
    conf = conf.dropna(subset=["genome"])
    return dict(zip(conf["strain_label"], conf["genome"])), conf


def _mapping_cols():
    m = S.load_mapping()
    return m["Well_souce_plate"], m["strain"]


def resolution_report(exp, conf):
    """Does 16S actually have the power to name these strains? (The answer decides everything.)"""
    counts = collections.Counter(conf["genome"])
    multi = {g: n for g, n in counts.items() if n > 1}
    rows = [
        {"metric": "n_labels_with_confident_call", "value": len(conf)},
        {"metric": "n_distinct_genomes_assigned", "value": conf["genome"].nunique()},
        {"metric": "n_genomes_claimed_by_multiple_labels", "value": len(multi)},
        {"metric": "n_labels_in_a_collision", "value": int(sum(multi.values()))},
        {"metric": "max_labels_on_one_genome", "value": int(max(counts.values()))},
        {"metric": "median_margin_over_runner_up", "value": round(float(conf["margin_over_runner_up"].median()), 5)},
        {"metric": "median_margin_in_bp_per_1420", "value": round(float(conf["margin_over_runner_up"].median()) * 1420, 1)},
        {"metric": "n_calls_exactly_tied", "value": int((conf["margin_over_runner_up"] <= 0.0005).sum())},
    ]
    out = pd.DataFrame(rows)
    out.to_csv(exp.out / f"s05_recovery_resolution_{exp.name}.csv", index=False)
    return out


def circularity_demo(exp, mapping, label, n_perm=200, seed=0):
    """Re-run the genome-vs-16S join test under a given well->genome map.

    Reported for both the original labels and the 16S-recovered map purely to show how large
    and how misleading the recovered number is. Do not use it as evidence either way.
    """
    ko = pd.read_csv(C.GENOMIC / "KEGG_ko_and_strains_table.csv", index_col=0)
    ko.index = ko.index.astype(str)
    pa = (ko > 0).astype(float)
    pa = pa[pa.sum(axis=1) > 0]
    r = pd.read_csv(C.ROOT / exp.name / "analysis/relative_abundance/outputs/r03_pair_replicate_stats.csv")
    d = r.assign(ga=r["strain_a"].map(mapping), gb=r["strain_b"].map(mapping)) \
         .dropna(subset=["ga", "gb", "ref_pair_bp_dist"])
    d = d[d["ga"].isin(pa.index) & d["gb"].isin(pa.index)]
    if len(d) < 20:
        return {"map": label, "n_pairs": len(d), "spearman": np.nan, "z": np.nan}
    A, B = pa.loc[d["ga"]].values, pa.loc[d["gb"]].values
    jac = 1 - (A * B).sum(1) / ((A + B) > 0).sum(1)
    rho = spearmanr(d["ref_pair_bp_dist"], jac)[0]
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_perm):
        perm = dict(zip(mapping.keys(), rng.permutation(list(mapping.values()))))
        ga, gb = d["strain_a"].map(perm), d["strain_b"].map(perm)
        k = ga.isin(pa.index) & gb.isin(pa.index)
        A2, B2 = pa.loc[ga[k]].values, pa.loc[gb[k]].values
        null.append(spearmanr(d.loc[k, "ref_pair_bp_dist"],
                              1 - (A2 * B2).sum(1) / ((A2 + B2) > 0).sum(1))[0])
    null = np.array(null)
    return {"map": label, "n_pairs": len(d), "spearman": round(float(rho), 3),
            "z": round(float((rho - null.mean()) / null.std()), 1),
            "median_ko_jaccard_at_16S_le_5bp": round(float(np.median(jac[d["ref_pair_bp_dist"].values <= 5])), 3)
            if (d["ref_pair_bp_dist"].values <= 5).any() else np.nan}
