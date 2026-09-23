"""How far upstream does the MIN_RESOLVABLE_BP / raw-bp read-assignment fix actually reach?

Answers "if we re-ran the earlier pipeline steps under the fix, what would change?" without
re-running them, by exploiting the fact that mapping_validation stores enough to reconstruct
the alternative answer exactly.

Two upstream steps, two different answers:

1.  **Per-strain 16S consensus** -- unaffected *by construction*, nothing to run. It is built de
    novo from the reads within each well (NanoFilt -> MAFFT -> plurality consensus -> HDBSCAN
    cluster split, plus greedy clustering of the resulting consensus sequences against each
    other). It never assigns a read to a reference, so no ambiguity margin and no resolution
    limit enter into it. The set of 16S sequences found cannot move.

2.  **mapping_validation** -- never had the margin bug either: `map_sample_constrained` assigns
    every read to its nearest candidate with no ambiguity cutoff at all. It does choose the
    winner on *normalized* distance (`ed / max(len(read), len(ref))`), and since the two
    candidate references differ in length the divisor differs between them, so the winner can
    in principle differ from the raw-bp winner. This script measures how often it actually does.

    `norm = ed / max(read_len, ref_len)` is invertible given the stored `read_len` and the
    reference lengths, so the exact bp distances come back without re-aligning anything.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
EXPS = HERE.parents[2]
OUT = HERE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

JOBS = [
    ("20260630",
     EXPS / "20260630/analysis/consensus2/strain_consensus_20260630.fasta",
     EXPS / "20260630/analysis/mapping_validation/outputs/03_edlib_read_assignments_consensus2.csv.gz"),
    ("20260721",
     EXPS / "20260721/analysis/consensus/merged_consensus_mono_priority_20260721.fasta",
     EXPS / "20260721/analysis/mapping_validation/outputs/03_edlib_read_assignments_merged_consensus_mono_priority.csv.gz"),
]


def ref_lengths(path):
    return {r.id.split("|")[0].removesuffix("_consensus"): len(r.seq)
            for r in SeqIO.parse(path, "fasta")}


def analyse(exp, fasta, reads_csv):
    ref = ref_lengths(fasta)
    d = pd.read_csv(reads_csv)
    d = d[d.n_candidates == 2].copy()
    for i in (1, 2):
        L = d[f"candidate{i}_strain"].map(ref)
        d[f"ed{i}"] = (d[f"candidate{i}_dist"] * np.maximum(d["read_len"], L)).round()
    d = d.dropna(subset=["ed1", "ed2"])

    d["norm_win"] = np.where(d.candidate1_dist < d.candidate2_dist, 1, 2)
    d["bp_win"] = np.where(d.ed1 < d.ed2, 1, np.where(d.ed2 < d.ed1, 2, 0))
    d["bp_sep"] = (d.ed1 - d.ed2).abs()
    flip = (d.bp_win != 0) & (d.norm_win != d.bp_win)

    # well-level: restricted to non-tie reads, so this is the effect of genuine disagreements
    # rather than of bp declining to call a tie that normalized distance would have called
    nt = d[d.bp_win != 0].copy()
    nt["f_norm"] = (nt.norm_win == 1).astype(float)
    nt["f_bp"] = (nt.bp_win == 1).astype(float)
    w = nt.groupby("sample_id").agg(n=("f_norm", "size"), f_norm=("f_norm", "mean"),
                                    f_bp=("f_bp", "mean"))
    w["abs_shift"] = (w.f_norm - w.f_bp).abs()

    rows = [
        {"metric": "two_candidate_reads", "value": len(d)},
        {"metric": "assignments_flipped_norm_vs_bp", "value": int(flip.sum())},
        {"metric": "pct_flipped", "value": round(100 * flip.mean(), 4)},
        {"metric": "pct_bp_ties", "value": round(100 * (d.bp_win == 0).mean(), 3)},
        {"metric": "wells_total", "value": int(nt.sample_id.nunique())},
        {"metric": "wells_with_any_changed_read", "value": int((w.abs_shift > 0).sum())},
        {"metric": "max_abs_shift_in_strain1_fraction", "value": round(float(w.abs_shift.max()), 4)},
        {"metric": "wells_shifting_gt_0.05", "value": int((w.abs_shift > 0.05).sum())},
    ]
    by_sep = (d.assign(flip=flip,
                       bin=pd.cut(d.bp_sep, [-.1, 0, 1, 2, 5, 10, 20, 50, 1e9],
                                  labels=["0", "1", "2", "3-5", "6-10", "11-20", "21-50", ">50"]))
              .groupby("bin", observed=True)
              .agg(n_reads=("flip", "size"), n_flipped=("flip", "sum")))
    by_sep["pct_flipped"] = (100 * by_sep.n_flipped / by_sep.n_reads).round(3)
    by_sep = by_sep.reset_index().assign(experiment=exp)

    w[w.abs_shift > 0].reset_index().assign(experiment=exp).to_csv(
        OUT / f"up02_{exp}_affected_wells.csv", index=False)
    return pd.DataFrame(rows).assign(experiment=exp), by_sep


if __name__ == "__main__":
    summaries, bins = [], []
    for exp, fa, rd in JOBS:
        if not Path(rd).exists():
            print(f"{exp}: per-read file missing, skipped -> {rd}")
            continue
        s, b = analyse(exp, fa, rd)
        summaries.append(s)
        bins.append(b)
        print(f"\n=== {exp} ===")
        print(s.to_string(index=False))
        print(b.to_string(index=False))
    pd.concat(summaries, ignore_index=True).to_csv(OUT / "up01_norm_vs_bp_summary.csv", index=False)
    pd.concat(bins, ignore_index=True).to_csv(OUT / "up03_flip_rate_by_bp_separation.csv", index=False)
    print(f"\nwritten -> {OUT}")
