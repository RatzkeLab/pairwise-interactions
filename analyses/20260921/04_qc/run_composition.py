"""QC: how many distinct 16S sequences are in each well?

A mono well should hold one, a pair well two. Deviations are the interesting part:
one sequence in a pair well means a strain was lost or the two cannot be told
apart; three or more anywhere means contamination.
"""
import sys, pickle
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE)); sys.path.insert(0, str(BASE.parent / "shared_scripts"))

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config, well_composition

OUT = BASE / "04_qc" / "outputs"; OUT.mkdir(parents=True, exist_ok=True)


def main(demux_dir=None, tag="corrected", processes=20):
    cfg = config.make_config(demux_dir)
    print(f"demux dir: {cfg.demux_dir}")
    df, reps, thr, pooled = well_composition.analyze(cfg, processes=processes)
    well_composition.summarize(df)

    df.to_csv(OUT / f"q01_well_composition_{tag}.csv", index=False)
    with open(OUT / f"q02_well_representatives_{tag}.pkl", "wb") as fh:
        pickle.dump({"reps": reps, "threshold": thr}, fh)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].hist(pooled, bins=np.linspace(0, 0.4, 81), color="#4C72B0")
    ax[0].axvline(thr, color="crimson", ls="--", lw=2, label=f"threshold = {thr}")
    ax[0].set(yscale="log", xlabel="normalized edit distance (all read pairs, pooled)",
              ylabel="read pairs", title="left mode = same strain, right = different")
    ax[0].legend()
    order = ["mono", "pair", "techrep", "ratio", "density"]
    present = [w for w in order if w in set(df.well_type)]
    for i, wt in enumerate(present):
        sub = df[df.well_type == wt]
        vals, counts = np.unique(sub.n_strains.clip(upper=4), return_counts=True)
        ax[1].bar(vals + (i - len(present) / 2) * 0.15, counts / counts.sum(),
                  width=0.15, label=f"{wt} (n={len(sub)})")
    ax[1].set(xlabel="distinct sequences detected", ylabel="fraction of wells",
              xticks=[1, 2, 3, 4], title="sequences per well by well type")
    ax[1].set_xticklabels(["1", "2", "3", "4+"])
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / f"q03_composition_{tag}.png", dpi=150)
    print(f"\nwrote -> {OUT}")
    return df, reps


if __name__ == "__main__":
    use_original = "--original" in sys.argv
    main(config.DEMUX_DIR_ORIGINAL if use_original else None,
         tag="original" if use_original else "corrected")
