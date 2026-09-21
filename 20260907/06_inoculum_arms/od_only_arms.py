"""Does the inoculum set the outcome, or does fitness? The half you can answer from OD.

The design added two arms to separate "a denser inoculum gives a head start" from
"a fitter strain grows dense everywhere and also wins". Two of the three questions
need only the endpoint OD and can be answered before any sequencing:

  monoculture contrast   each strain has two monos, 200 nL and 100 nL of itself.
                         If yield depends on inoculum at all, these differ. If the
                         cultures are at carrying capacity by readout, they do not.
                         Paired within strain, so between-strain density differences
                         cancel -- this is the cleanest measurement in the plate.

  density titration      30 pairs at 50:50, 100:100 and 200:200. Ratio is held at
                         1:1, so anything that moves is absolute inoculum, not
                         composition.

The third question -- the slope of final log-ratio on initial log-ratio, from the
ratio titration -- needs read counts and waits for the sequencing.

Reading this the right way round matters. A flat monoculture contrast means the
assay reports fitness rather than starting conditions, which is what you want. A
steep one means the readout is partly just reporting the inoculum.
"""
import sys
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE)); sys.path.insert(0, str(BASE.parent / "shared_scripts"))

import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config, plate_reader as pr

OUT = BASE / "06_inoculum_arms" / "outputs"; OUT.mkdir(parents=True, exist_ok=True)
QC = BASE / "02_demux_qc" / "outputs"


def load():
    od, _ = pr.load_folder(config.OD_FULL, testname="ODFull")
    lay = pd.read_csv(BASE / "01_setup" / "strain_layout_20260907.csv")
    lay = lay[lay.plate_role == "primary_sequenced"]
    df = lay.merge(od[["dest_plate", "dest_well", "od"]],
                   on=["dest_plate", "dest_well"], how="left", validate="1:1")
    blank = od[(od.dest_plate <= 8)].od.quantile(0.01)
    df["od_net"] = df.od - blank
    # strains that never grew carry no information about inoculum response
    dead = set(pd.read_csv(QC / "d13_strains_no_growth.csv").strain)
    df["has_dead"] = df.strain1.isin(dead) | df.strain2.isin(dead)
    print(f"blank OD = {blank:.3f}; {len(dead)} non-growing strains excluded from the arms")
    return df


def mono_contrast(df):
    print("\n" + "=" * 72)
    print("MONOCULTURE CONTRAST -- is final yield inoculum-dependent at all?")
    print("=" * 72)
    m = df[(df.well_type == "mono") & (~df.has_dead)]
    w = m.pivot_table(index="strain1", columns="total_inoculum_nL", values="od_net")
    w = w.dropna()
    if not {100, 200}.issubset(w.columns):
        print("  both inoculum levels not present"); return None
    lo, hi = w[100], w[200]
    print(f"  {len(w)} strains with both a 100 nL and a 200 nL monoculture")
    print(f"  median OD  100 nL {lo.median():.3f}   200 nL {hi.median():.3f}   "
          f"ratio {hi.median()/lo.median():.3f}")
    d = hi - lo
    t = stats.wilcoxon(hi, lo)
    print(f"  paired difference (200 - 100): median {d.median():+.4f}, "
          f"Wilcoxon p = {t.pvalue:.3g}")
    # effect size relative to spread: what fraction of a doubling shows up?
    lr = np.log2(hi.clip(lower=1e-3) / lo.clip(lower=1e-3))
    print(f"  log2(OD200/OD100): median {lr.median():+.3f}  "
          f"(a pass-through of the 2x inoculum would be +1.000, "
          f"carrying capacity would be 0.000)")
    print(f"  -> {lr.median()/1.0:.1%} of the doubling survives to the endpoint")
    if abs(lr.median()) < 0.15:
        print("  => cultures are essentially at carrying capacity by readout: the "
              "endpoint\n     reports fitness, not the head start. That is the good case.")
    else:
        print("  => a measurable part of the inoculum carries through to the endpoint; "
              "the\n     assay is partly reporting starting density.")
    return w


def density_titration(df):
    print("\n" + "=" * 72)
    print("DENSITY TITRATION -- does absolute inoculum change the outcome at 1:1?")
    print("=" * 72)
    d = df[(df.well_type == "density") & (~df.has_dead)].copy()
    d["pair"] = d.apply(lambda r: tuple(sorted((r.strain1, r.strain2))), axis=1)
    print(f"  {d.pair.nunique()} pairs x {sorted(d.total_inoculum_nL.unique())} nL")
    g = d.groupby("total_inoculum_nL").od_net.agg(["size", "median"])
    print(g.round(4).to_string())
    w = d.pivot_table(index="pair", columns="total_inoculum_nL", values="od_net").dropna()
    if w.shape[1] >= 2:
        lo, hi = w[w.columns.min()], w[w.columns.max()]
        lr = np.log2(hi.clip(lower=1e-3) / lo.clip(lower=1e-3))
        fold = w.columns.max() / w.columns.min()
        span = np.log2(fold)
        print(f"\n  {len(w)} pairs with all levels; log2(OD at {w.columns.max()} nL / "
              f"OD at {w.columns.min()} nL) median {lr.median():+.3f}")
        print(f"  a pass-through of the {fold:.0f}-fold inoculum span would give "
              f"{span:+.3f}; observed is {lr.median()/span:.1%} of it")
        print(f"  Wilcoxon p = {stats.wilcoxon(hi, lo).pvalue:.3g}")
    return w


def main():
    df = load()
    w_mono = mono_contrast(df)
    w_dens = density_titration(df)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    if w_mono is not None:
        ax[0].scatter(w_mono[100], w_mono[200], s=14, alpha=.6, color="#2a78d6")
        lim = [0, max(w_mono.max()) * 1.05]
        ax[0].plot(lim, lim, ls="--", c="grey", lw=1, label="no inoculum effect")
        ax[0].plot(lim, [2 * x for x in lim], ls=":", c="crimson", lw=1,
                   label="inoculum passes straight through")
        ax[0].set(xlim=lim, ylim=lim, xlabel="OD, 100 nL monoculture",
                  ylabel="OD, 200 nL monoculture",
                  title="monoculture contrast (paired within strain)")
        ax[0].legend(fontsize=8)
    if w_dens is not None and w_dens.shape[1] >= 2:
        for i, c in enumerate(w_dens.columns):
            ax[1].scatter(np.full(len(w_dens), i) + np.random.uniform(-.08, .08, len(w_dens)),
                          w_dens[c], s=14, alpha=.6, color="#2a78d6")
        for _, row in w_dens.iterrows():
            ax[1].plot(range(len(w_dens.columns)), row.values, color="grey", lw=.4, alpha=.4)
        ax[1].set(xticks=range(len(w_dens.columns)),
                  xticklabels=[f"{int(c)} nL" for c in w_dens.columns],
                  ylabel="OD", title="density titration, ratio held at 1:1")
    fig.tight_layout(); fig.savefig(OUT / "i01_inoculum_arms_od.png", dpi=150)
    if w_mono is not None: w_mono.to_csv(OUT / "i02_mono_contrast.csv")
    if w_dens is not None: w_dens.to_csv(OUT / "i03_density_titration.csv")
    print(f"\nwrote -> {OUT}")


if __name__ == "__main__":
    main()
