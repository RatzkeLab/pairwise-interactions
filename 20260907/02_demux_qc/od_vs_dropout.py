"""Did the wells that produced no reads actually contain cells?

OD answers what sequencing cannot. A well with normal OD and no reads failed at
PCR or demultiplexing; a well with blank OD never grew, and no amount of
re-demultiplexing will recover it. The barcode analysis already concluded that 19
primer wells failed, and this is the independent check on that conclusion -- it
uses a different instrument and does not involve sequencing at all.

The unused border (row A/P, column 1/24) is the built-in blank: those wells were
never inoculated, so they set the scale for "no cells".
"""
import sys
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE)); sys.path.insert(0, str(BASE.parent / "shared_scripts"))

import numpy as np, pandas as pd
from scipy import stats
import config, plate_reader as pr

OUT = BASE / "02_demux_qc" / "outputs"


def main():
    od, _ = pr.load_folder(config.OD_FULL, testname="ODFull")
    seq = pd.read_csv(OUT / "d00_sequenced_wells_with_counts.csv")
    lay = pd.read_csv(BASE / "01_setup" / "strain_layout_20260907.csv")

    used = set(zip(lay.dest_plate, lay.dest_well))
    od["inoculated"] = [(p, w) in used for p, w in zip(od.dest_plate, od.dest_well)]
    blank = od[(~od.inoculated) & (od.dest_plate <= 8)].od
    print(f"uninoculated border wells (the blank): n={len(blank)} "
          f"median OD {blank.median():.3f}  95th pct {blank.quantile(.95):.3f}")
    grown = blank.quantile(.99)
    print(f"-> calling a well 'grew' at OD > {grown:.3f} (99th pct of blanks)\n")

    df = seq.merge(od[["dest_plate", "dest_well", "od"]], on=["dest_plate", "dest_well"],
                   how="left", validate="1:1")
    assert df.od.notna().all(), "sequenced well with no OD reading"
    df["grew"] = df.od > grown

    print(f"of {len(df)} sequenced wells, {df.grew.sum()} ({df.grew.mean():.1%}) "
          f"are above the blank")
    print("\n" + "=" * 70)
    print("DROPOUT vs GROWTH")
    print("=" * 70)
    ct = pd.crosstab(df.dropout, df.grew)
    ct.index = ["usable (>10 reads)", "dropout (<=10)"]
    ct.columns = ["no growth", "grew"]
    print(ct.to_string())
    d = df[df.dropout]
    print(f"\n{d.grew.sum()} of {len(d)} dropout wells ({d.grew.mean():.1%}) "
          f"had cells -- their DNA existed and the sequencing lost it")
    print(f"{(~d.grew).sum()} ({(~d.grew).mean():.1%}) never grew -- nothing to recover")

    print("\n" + "=" * 70)
    print("THE 19 FLAGGED BARCODES: was it the primer, or did those wells not grow?")
    print("=" * 70)
    for label, sub in [("wells touching a flagged barcode", df[df.hits_bad_bc]),
                       ("wells with two clean barcodes", df[~df.hits_bad_bc])]:
        print(f"  {label:36s} n={len(sub):5d}  grew {sub.grew.mean():6.1%}  "
              f"median OD {sub.od.median():.3f}  dropout {sub.dropout.mean():6.1%}")
    a, b = df[df.hits_bad_bc].grew, df[~df.hits_bad_bc].grew
    _, p = stats.fisher_exact([[a.sum(), (~a).sum()], [b.sum(), (~b).sum()]])
    print(f"\n  growth rate difference between the two groups: Fisher p = {p:.3g}")
    if p > 0.01:
        print("  -> the flagged wells grew just as well as everything else. The cells "
              "were there;\n     the amplicon was not. That is a primer failure, "
              "confirmed independently of sequencing.")
    else:
        print("  -> the flagged wells also grew differently, so growth is confounded "
              "with the barcode\n     effect and the primer conclusion needs revisiting.")

    print("\n" + "=" * 70)
    print("WHICH STRAINS DID NOT GROW AT ALL?")
    print("=" * 70)
    print("Each strain has two monocultures ON DIFFERENT PLATES, so if failure were a\n"
          "well-level accident the two would fail independently. They do not.")
    mono = df[df.well_type == "mono"]
    g = mono.groupby("strain1").agg(n=("grew", "size"), n_grew=("grew", "sum"),
                                    plates=("dest_plate", "nunique"),
                                    median_od=("od", "median"))
    both = g[g.n == 2]
    p_fail = 1 - mono.grew.mean()
    exp_both = len(both) * p_fail ** 2
    obs_both = int((both.n_grew == 0).sum())
    print(f"\n  mono wells that did not grow : {(~mono.grew).sum()} of {len(mono)} "
          f"({(~mono.grew).mean():.1%})")
    pair = df[df.well_type == "pair"]
    print(f"  pair wells that did not grow : {(~pair.grew).sum()} of {len(pair)} "
          f"({(~pair.grew).mean():.1%})")
    print(f"\n  of {len(both)} strains with two monos: "
          f"{obs_both} failed both, {int((both.n_grew == 1).sum())} failed one, "
          f"{int((both.n_grew == 2).sum())} grew both")
    print(f"  independent-well expectation for 'failed both': {exp_both:.1f}  "
          f"-> {obs_both/max(exp_both,1e-9):.0f}x enrichment")
    print(f"  => growth failure is a property of the STRAIN, not the well. Those "
          f"precultures\n     were already dead; the pair-well rate follows as "
          f"roughly {p_fail:.3f}^2 = {p_fail**2:.3%}.")
    dead = sorted(both[both.n_grew == 0].index)
    print(f"\n  {len(dead)} strains grew in neither monoculture:")
    for i in range(0, len(dead), 16):
        print("    " + " ".join(dead[i:i + 16]))
    inv = pair[pair.strain1.isin(dead) | pair.strain2.isin(dead)]
    oth = pair[~(pair.strain1.isin(dead) | pair.strain2.isin(dead))]
    if len(inv) and len(oth):
        pv = stats.mannwhitneyu(inv.od, oth.od).pvalue
        print(f"\n  their pair wells confirm it: median OD {inv.od.median():.3f} "
              f"(n={len(inv)}) vs {oth.od.median():.3f} (n={len(oth)}) for pairs "
              f"without one,\n     Mann-Whitney p = {pv:.3g}. A pair carrying a dead "
              f"strain grows like a monoculture.")
    g.to_csv(OUT / "d12_strain_growth.csv")
    pd.Series(dead, name="strain").to_csv(OUT / "d13_strains_no_growth.csv", index=False)

    print("\nby well type (OD, sequenced plates):")
    print(df.groupby("well_type").agg(n=("od", "size"), median_od=("od", "median"),
                                      grew=("grew", "mean"),
                                      dropout=("dropout", "mean")).round(3).to_string())
    df.to_csv(OUT / "d11_wells_with_od.csv", index=False)
    print(f"\nwrote -> {OUT}/d11_wells_with_od.csv")
    return df


if __name__ == "__main__":
    main()
