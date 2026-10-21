"""Stage 6: how large can a community be and still be read out? A preview for design.

The largest set of strains in which every pair is at least T bp apart (haplotype-aware)
is a maximum independent set on the "too close" graph; solved exactly as an integer
program. It is the ceiling on community size at threshold T, before any other design
constraint, and it sets how much the choice of T costs.

Eligibility (the minimum that any design would apply):
  core       tier A or B, grew in 20260907, not mixed_source, haplotypes within 10 bp of
             each other (more than that is two organisms, not operon variation)
  extended   core + tier C

Output: d03_max_separable.csv, and one example set per (eligibility, T) -- an example,
not a recommendation: many sets reach the same size.
"""
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import config  # noqa: E402

import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

OUT = config.DIST_OUT
THRESHOLDS = (5, 7, 8, 10, 15)
MAX_HAP_DIVERGENCE = 10


def eligible(ref, tiers):
    r = ref[ref.tier.isin(tiers) & ~ref.no_growth.astype(bool)
            & ~ref.mixed_source.astype(bool)
            & (ref.haplotype_divergence_bp <= MAX_HAP_DIVERGENCE)]
    return sorted(r.index)


def max_separable(names, pairs, thr):
    idx = {n: i for i, n in enumerate(names)}
    e = pairs[(pairs.bp < thr) & pairs.strain_a.isin(idx) & pairs.strain_b.isin(idx)]
    cons = []
    if len(e):
        A = lil_matrix((len(e), len(names)))
        for k, (a, b) in enumerate(zip(e.strain_a, e.strain_b)):
            A[k, idx[a]] = 1; A[k, idx[b]] = 1
        cons = [LinearConstraint(A.tocsr(), -np.inf, 1)]
    res = milp(c=-np.ones(len(names)), constraints=cons, integrality=np.ones(len(names)),
               bounds=Bounds(0, 1), options={"time_limit": 300})
    chosen = [n for n, v in zip(names, res.x) if v > 0.5]
    return chosen, res.status == 0


def main():
    ref = pd.read_csv(config.REFDB_OUT / "s01_strain_reference.csv").set_index("strain")
    pairs = pd.read_csv(OUT / "d01_pairs.csv")
    rows = []
    for label, tiers in (("core", ["A", "B"]), ("extended", ["A", "B", "C"])):
        names = eligible(ref, tiers)
        for thr in THRESHOLDS:
            chosen, optimal = max_separable(names, pairs, thr)
            rows.append(dict(eligibility=label, n_eligible=len(names), min_separation_bp=thr,
                             max_community=len(chosen), proven_optimal=optimal,
                             example_set=" ".join(chosen)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "d03_max_separable.csv", index=False)
    print(df.pivot(index="min_separation_bp", columns="eligibility", values="max_community")
          .rename(columns=lambda c: f"{c} (n={df[df.eligibility == c].n_eligible.iloc[0]})"))
    return df


if __name__ == "__main__":
    main()
