"""Why are so many 20260907 wells empty after demultiplexing?

Three candidate explanations, and the design lets us separate them because
barcodes were assigned to wells at random (fwd_idx/rev_idx are uncorrelated
with plate position), unlike a conventional row=forward / column=reverse layout:

  primer/PCR   -> dropout concentrates on particular barcode indices
  strain       -> dropout concentrates on particular strains, across plates
  handling     -> dropout concentrates on plate positions (row/col/edge/plate)

Each factor gets the same treatment: observed between-group variance in dropout
rate against the binomial null, plus per-group FDR-corrected binomial tests.
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats

BASE = Path("/home/rl/scripts/karl/pairwise_interaction_experiments/20260907")
DEMUX_ROOT = Path("/home/rl/data/interim/karl/demultiplexing")
DEMUX_ORIGINAL = DEMUX_ROOT / "20260907_demux"
DEMUX_CORRECTED = DEMUX_ROOT / "20260907_demux_corrected"
OUT = BASE / "02_demux_qc" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

DROPOUT_MAX_READS = 10      # a well with <=10 reads cannot be called; that is the failure
                            # (read_counts elsewhere in this project use >10 as usable)


def load(demux=None):
    demux = demux or DEMUX_CORRECTED
    lay = pd.read_csv(BASE / "01_setup" / "strain_layout_20260907.csv")
    cnt = pd.read_csv(demux / "summary" / "demultiplexed_read_counts.tsv", sep="\t")
    lay["sample"] = lay.apply(
        lambda r: f"Plate{int(r.dest_plate):02d}_{r.dest_well}", axis=1)
    df = lay.merge(cnt, on="sample", how="left", validate="1:1")
    assert df.read_count.notna().all(), "layout well with no demux entry"
    df["read_count"] = df.read_count.astype(int)
    df["row_idx"] = df.dest_row.map({c: i for i, c in enumerate("ABCDEFGHIJKLMNOP")})
    df["dropout"] = df.read_count <= DROPOUT_MAX_READS
    return df


def overdispersion(df, key, label):
    """Between-group variance in dropout rate vs the binomial null.

    Under 'this factor does not matter', each group's failure count is
    Binomial(n_g, p) with p the global rate. Pearson chi2 over groups tests that.
    The dispersion ratio (chi2/df) is the effect size: 1.0 = no structure.
    """
    g = df.groupby(key).agg(n=("dropout", "size"), k=("dropout", "sum"))
    g = g[g.n >= 5]
    p = df.dropout.mean()
    exp = g.n * p
    chi2 = (((g.k - exp) ** 2) / (exp * (1 - p))).sum()
    dof = len(g) - 1
    pval = stats.chi2.sf(chi2, dof)
    g["rate"] = g.k / g.n
    # per-group two-sided binomial test, BH-corrected
    g["p_binom"] = [stats.binomtest(int(k), int(n), p).pvalue
                    for k, n in zip(g.k, g.n)]
    order = np.argsort(g.p_binom.values)
    ranks = np.empty(len(g), int); ranks[order] = np.arange(1, len(g) + 1)
    g["q"] = np.minimum(1, g.p_binom.values * len(g) / ranks)
    g = g.sort_values("rate", ascending=False)
    print(f"\n{label}  ({len(g)} groups, global dropout {p:.1%})")
    print(f"  chi2={chi2:.0f}  df={dof}  dispersion={chi2/dof:.2f}  p={pval:.3g}")
    print(f"  groups significant at q<0.05: {(g.q < 0.05).sum()}")
    print(f"  rate range: {g.rate.min():.1%} .. {g.rate.max():.1%}"
          f"   (binomial null sd of a rate ~ {np.sqrt(p*(1-p)/g.n.median()):.1%})")
    return g, dict(factor=label, n_groups=len(g), chi2=chi2, df=dof,
                   dispersion=chi2 / dof, p=pval, n_sig=int((g.q < 0.05).sum()))


def main(demux=None):
    demux = demux or DEMUX_CORRECTED
    print(f"demux: {demux.name}")
    df = load(demux)

    print("=" * 72)
    print("STRUCTURAL ACCOUNTING")
    print("=" * 72)
    for role, sub in df.groupby("plate_role"):
        z = (sub.read_count == 0).sum()
        d = sub.dropout.sum()
        print(f"{role:20s} wells={len(sub):5d}  reads={sub.read_count.sum():>8,}"
              f"  zero={z:5d} ({z/len(sub):5.1%})  <=10 reads={d:5d} ({d/len(sub):5.1%})")
    ext = df[df.plate_role == "extension_frozen"]
    print(f"\nExtension plates were never PCR'd or sequenced, yet carry "
          f"{ext.read_count.sum():,} reads.")
    print(f"  -> that is the barcode mis-assignment floor: "
          f"{ext.read_count.sum()/df.read_count.sum():.2%} of all assigned reads.")

    seq = df[df.plate_role == "primary_sequenced"].copy()
    print(f"\nFrom here on: the {len(seq)} sequenced wells only.")
    print(f"  usable (>{DROPOUT_MAX_READS} reads): {(~seq.dropout).sum()} "
          f"({(~seq.dropout).mean():.1%})")
    print(f"  dropout  (<={DROPOUT_MAX_READS} reads): {seq.dropout.sum()} "
          f"({seq.dropout.mean():.1%})   of which exactly 0: {(seq.read_count==0).sum()}")

    print("\n" + "=" * 72)
    print("IS IT THE PRIMERS, THE STRAINS, OR THE PLATE POSITION?")
    print("=" * 72)

    stats_rows, tables = [], {}
    for key, label in [("fwd_idx", "FORWARD barcode index"),
                       ("rev_idx", "REVERSE barcode index"),
                       ("dest_plate", "destination PLATE"),
                       ("dest_row", "destination ROW"),
                       ("dest_col", "destination COLUMN")]:
        g, s = overdispersion(seq, key, label)
        stats_rows.append(s); tables[key] = g

    # strain: a well fails for strain reasons if either member is a bad strain,
    # so score every (well, strain) membership
    mem = pd.concat([
        seq[["sample", "dropout", "strain1"]].rename(columns={"strain1": "strain"}),
        seq[["sample", "dropout", "strain2"]].rename(columns={"strain2": "strain"}),
    ]).dropna(subset=["strain"])
    g, s = overdispersion(mem, "strain", "STRAIN (per membership)")
    stats_rows.append(s); tables["strain"] = g

    summary = pd.DataFrame(stats_rows).sort_values("dispersion", ascending=False)
    print("\n" + "=" * 72)
    print("RANKED BY EFFECT SIZE (dispersion = chi2/df; 1.0 means no structure)")
    print("=" * 72)
    print(summary.to_string(index=False,
          float_format=lambda v: f"{v:.3g}"))

    summary.to_csv(OUT / "d01_factor_dispersion.csv", index=False)
    for k, g in tables.items():
        g.to_csv(OUT / f"d02_dropout_by_{k}.csv")
    seq.to_csv(OUT / "d00_sequenced_wells_with_counts.csv", index=False)
    print(f"\nwrote -> {OUT}")
    return df, seq, tables, summary


if __name__ == "__main__":
    import sys
    main(DEMUX_ORIGINAL if "--original" in sys.argv else DEMUX_CORRECTED)
