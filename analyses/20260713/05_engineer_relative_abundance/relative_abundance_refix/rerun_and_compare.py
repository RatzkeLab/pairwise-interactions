"""Re-run 20260721's relative_abundance under the MIN_RESOLVABLE_BP fix, and diff it.

20260630's relative_abundance was re-run when the ambiguous-margin bug was fixed on 2026-08-27
(read assignment now decided in raw bp, ambiguous only on a genuine tie; pairs closer than
MIN_RESOLVABLE_BP flagged `below_resolution_limit` up front). 20260721's outputs still predate
that fix -- its r03 has no `below_resolution_limit` column.

This writes the re-run into its OWN directory and leaves
`analysis/relative_abundance/outputs/` untouched, so old and new can be compared rather than
one silently replacing the other.

Scope note: the fix cannot touch anything upstream of here.
  - The per-strain 16S consensus is built de novo from within-well reads (NanoFilt -> MAFFT ->
    plurality consensus -> HDBSCAN cluster split). It never assigns reads to a reference, so
    the set of 16S sequences found is unaffected by construction.
  - mapping_validation assigns each read to its nearest reference with no ambiguity cutoff, so
    it never had the bug. It does pick the winner on *normalized* distance, and because the two
    references differ in length the divisor differs between them -- but that flips only
    10/83105 reads (20260630) and 42/110282 (20260721), all at 1-2 bp reference separation.
"""

import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
sys.path.insert(0, str(ANALYSIS))
sys.path.insert(0, str(ANALYSIS.parents[1] / "shared_pipelines"))

import pandas as pd
from dataclasses import replace

from config import CFG
import relative_abundance as ra

OUT = HERE / "outputs"
OLD = ANALYSIS / "relative_abundance" / "outputs"


def main(stages=("r01", "r02", "r03")):
    cfg = replace(CFG, relative_abundance_out_dir=OUT)
    t0 = time.time()
    if "r01" in stages:
        ra.reference_distances(cfg)
        print(f"  r01 done {time.time()-t0:.0f}s")
    if "r02" in stages:
        ra.compute_interaction_scores(cfg)
        print(f"  r02 done {time.time()-t0:.0f}s")
    if "r03" in stages:
        ra.replicate_stability(cfg)
        print(f"  r03 done {time.time()-t0:.0f}s")

    new = pd.read_csv(OUT / "r03_pair_replicate_stats.csv")
    old = pd.read_csv(OLD / "r03_pair_replicate_stats.csv")
    key = ["strain_a", "strain_b"]
    m = old.merge(new, on=key, suffixes=("_old", "_new"))

    rows = [
        {"metric": "n_pairs", "old": len(old), "new": len(new)},
        {"metric": "n_high_uncertainty_pair", "old": int(old.high_uncertainty_pair.sum()),
         "new": int(new.high_uncertainty_pair.sum())},
        {"metric": "n_usable (not high_uncertainty)", "old": int((~old.high_uncertainty_pair).sum()),
         "new": int((~new.high_uncertainty_pair).sum())},
        {"metric": "n_unstable_replicate", "old": int(old.unstable_replicate.sum()),
         "new": int(new.unstable_replicate.sum())},
    ]
    if "below_resolution_limit" in new.columns:
        rows.append({"metric": "n_below_resolution_limit (new only)", "old": float("nan"),
                     "new": int(new.below_resolution_limit.sum())})
    d = (m.mean_log2_ratio_a_over_b_new - m.mean_log2_ratio_a_over_b_old).abs()
    rows += [
        {"metric": "pairs in both", "old": len(m), "new": len(m)},
        {"metric": "labels changed >0.01 (log2 ratio)", "old": float("nan"), "new": int((d > 0.01).sum())},
        {"metric": "median |change| among changed", "old": float("nan"),
         "new": round(float(d[d > 0.01].median()), 3) if (d > 0.01).any() else 0.0},
        {"metric": "max |change|", "old": float("nan"), "new": round(float(d.max()), 3)},
    ]
    cmp = pd.DataFrame(rows)
    cmp.to_csv(OUT / "cmp01_old_vs_new_r03.csv", index=False)
    m.to_csv(OUT / "cmp02_pair_level_join.csv", index=False)
    print("\n=== 20260721 r03: pre-fix vs post-fix ===")
    print(cmp.to_string(index=False))
    return cmp


if __name__ == "__main__":
    main()
