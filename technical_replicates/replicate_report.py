"""How reproducible is a technical replicate -- by plate reader, and by sequencing?

This is the baseline every cross-experiment comparison has to be read against. A
cross-experiment correlation of 0.08 means one thing if replicates of the SAME pair in the SAME
experiment reach 0.9, and something else entirely if they only reach 0.3.

Design points:

  **These are cross-plate replicates.** 96% of pairs with >=2 wells have them on different
  destination plates, so the number measures the whole pipeline -- pipetting, plate, position,
  reader session -- not just well-to-well noise on one plate. That is the right baseline for
  comparing across experiments.

  **All 61 wavelengths, not just 600 nm.** The spectra run 350-950 nm in 10 nm steps. OD600 is
  convention, not necessarily the most reproducible channel for these cultures: scattering rises
  toward the blue while absorbance from media components and pigments falls off in the red, so
  the best signal-to-noise may sit elsewhere. Reproducibility is computed per wavelength so the
  choice is made from data.

  **Sequencing replicates are scored on the same footing** -- per-well relative abundance for
  the same pair in different wells -- so the two modalities are directly comparable.

  20260721's sequencing uses the post-MIN_RESOLVABLE_BP re-run (relative_abundance_refix),
  since its original outputs predate that fix.
"""

import glob
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

HERE = Path(__file__).resolve().parent
EXPS = HERE.parent
OUT = HERE / "outputs"
FIG = OUT / "figures"
for d in (OUT, FIG):
    d.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(EXPS / "shared_pipelines"))
from genomic_ml import COLOR_BLUE, COLOR_RED, COLOR_GRID, COLOR_TEXT_SECONDARY, COLOR_GOOD

PR = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/Karl_2026")
WAVELENGTHS = np.arange(350, 350 + 61 * 10, 10)

CONFIG = {
    "20260630": {"od": PR / "Karl_20260704_OD_Full",
                 "layout": EXPS / "20260630/setup/strain_layout_20260630.csv",
                 "seq": EXPS / "20260630/analysis/relative_abundance/outputs"},
    "20260721": {"od": PR / "Karl_20260723_OD_Full",
                 "layout": EXPS / "20260721/setup/strain_layout_20260721_plate1_2_swapped.csv",
                 "seq": EXPS / "20260721/analysis/relative_abundance_refix/outputs"},
}


def read_spectrum(path):
    txt = Path(path).read_text(encoding="latin-1").splitlines()
    plate = next((l.split("ID1:")[1].split("ID2")[0].strip() for l in txt[:6] if "ID1:" in l), None)

    def num(x):
        try:
            return float(x.strip())
        except ValueError:
            return np.nan
    rows = [[num(v) for v in l.split(",")] for l in txt if re.match(r"^\s*[\d.\-\s]+,", l)]
    rows = [r for r in rows if len(r) == 24]
    a = np.array(rows)
    n_wl = a.shape[0] // 16
    return plate, a.reshape(n_wl, 16, 24)


def load_od(name):
    """(n_layout_wells, 61) OD matrix aligned to the layout, plus the layout."""
    cfg = CONFIG[name]
    cube = {}
    for f in sorted(glob.glob(str(cfg["od"] / "*.csv"))):
        p, c = read_spectrum(f)
        if p and p.isdigit() and c.shape[0] == 61:
            cube.setdefault(int(p), []).append(c)
    cube = {k: np.nanmean(v, axis=0) for k, v in cube.items()}
    lay = pd.read_csv(cfg["layout"])
    M = np.full((len(lay), 61), np.nan)
    for i, r in enumerate(lay.itertuples()):
        c = cube.get(int(r.dest_plate))
        if c is None:
            continue
        row = ord(str(r.dest_row).upper()[0]) - 65
        col = int(r.dest_col) - 1
        if 0 <= row < 16 and 0 <= col < 24:
            M[i] = c[:, row, col]
    lay["pair_key"] = [frozenset((a, b)) for a, b in zip(lay.strain1, lay.strain2)]
    print(f"  {name}: {len(cube)} plates, {np.isfinite(M[:, 25]).sum()}/{len(lay)} wells with OD600")
    return lay, M


def split_half(values, groups, plates=None):
    """First replicate vs mean of the rest, per group. Returns rho, r, n, median |diff|."""
    df = pd.DataFrame({"g": groups, "v": values})
    df = df.dropna()
    g = df.groupby("g")["v"].agg(list)
    g = g[g.map(len) >= 2]
    if len(g) < 20:
        return dict(spearman=np.nan, pearson=np.nan, n_pairs=0, median_abs_diff=np.nan)
    x = np.array([v[0] for v in g])
    y = np.array([np.mean(v[1:]) for v in g])
    return dict(spearman=float(spearmanr(x, y)[0]), pearson=float(pearsonr(x, y)[0]),
                n_pairs=int(len(g)), median_abs_diff=float(np.median(np.abs(x - y))))


def od_replicates():
    rows, spectra = [], {}
    for name in CONFIG:
        lay, M = load_od(name)
        # z-score within plate at each wavelength: plate-to-plate offsets are batch, not signal
        pl = lay["dest_plate"].values
        Z = np.full_like(M, np.nan)
        for p in np.unique(pl):
            m = pl == p
            sub = M[m]
            Z[m] = (sub - np.nanmean(sub, axis=0)) / np.nanstd(sub, axis=0)
        for j, wl in enumerate(WAVELENGTHS):
            r = split_half(Z[:, j], lay["pair_key"].values)
            rows.append({"experiment": name, "wavelength_nm": int(wl), **r,
                         "median_raw_OD": float(np.nanmedian(M[:, j]))})
        spectra[name] = (lay, M, Z)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "r01_od_replicate_by_wavelength.csv", index=False)
    return df, spectra


def seq_replicates():
    rows, detail = [], {}
    for name, cfg in CONFIG.items():
        f = Path(cfg["seq"]) / "r02_well_interaction_scores.csv"
        if not f.exists():
            print(f"  {name}: no r02 at {f}")
            continue
        w = pd.read_csv(f)
        w["pair_key"] = [frozenset((a, b)) for a, b in zip(w.strain1, w.strain2)]
        # orient consistently: abundance of the alphabetically-first strain
        first = [sorted(k)[0] for k in w.pair_key]
        w["ra_oriented"] = np.where(w.strain1.values == np.array(first),
                                    w.relative_abundance_a, 1 - w.relative_abundance_a)
        base = split_half(w.ra_oriented.values, w.pair_key.values)
        rows.append({"experiment": name, "subset": "all sequenced wells", **base,
                     "median_reads": float(w.n_reads.median())})
        for lo, hi, lab in [(0, 20, "<20 reads"), (20, 50, "20-50 reads"), (50, 1e9, ">50 reads")]:
            s = w[(w.n_reads >= lo) & (w.n_reads < hi)]
            rows.append({"experiment": name, "subset": lab,
                         **split_half(s.ra_oriented.values, s.pair_key.values),
                         "median_reads": float(s.n_reads.median()) if len(s) else np.nan})
        detail[name] = w
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "r02_sequencing_replicates.csv", index=False)
    return df, detail


def make_figures(od_df, seq_df, spectra, matched=None):
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    colors = {"20260630": COLOR_BLUE, "20260721": COLOR_RED}

    ax = axes[0]
    for name, d in od_df.groupby("experiment"):
        ax.plot(d.wavelength_nm, d.spearman, "-", lw=2, color=colors[name], label=name)
        best = d.loc[d.spearman.idxmax()]
        ax.scatter([best.wavelength_nm], [best.spearman], s=60, color=colors[name], zorder=5)
        ax.annotate(f"{int(best.wavelength_nm)} nm\nρ={best.spearman:.3f}",
                    (best.wavelength_nm, best.spearman), textcoords="offset points",
                    xytext=(6, -20), fontsize=8, color=colors[name])
    ax.axvline(600, color=COLOR_TEXT_SECONDARY, ls="--", lw=1)
    ax.text(605, ax.get_ylim()[0] + .02, "OD600\n(convention)", fontsize=8,
            color=COLOR_TEXT_SECONDARY)
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("technical-replicate Spearman")
    ax.set_title("Plate reader: replicate agreement by wavelength", fontsize=10)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    for name, d in od_df.groupby("experiment"):
        ax.plot(d.wavelength_nm, d.median_raw_OD, "-", lw=2, color=colors[name], label=name)
    ax.axvline(600, color=COLOR_TEXT_SECONDARY, ls="--", lw=1)
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("median raw OD")
    ax.set_title("Signal level by wavelength\n(reproducibility is not just signal size)",
                 fontsize=10)
    ax.legend(frameon=False, fontsize=9)

    # depth stratified on each PAIR's minimum well depth -- stratifying per well leaves the
    # low-depth bin empty, since few pairs have both replicates under 20 reads
    ax = axes[2]
    subs = ["pair min <25 reads", "pair min 25-45", "pair min >=45"]
    w = 0.35
    m = matched[matched.subset.isin(subs)]
    for i, (name, d) in enumerate(m.groupby("experiment")):
        d = d.set_index("subset").reindex(subs)
        b = ax.bar(np.arange(len(subs)) + i * w, d.spearman, w, color=colors[name], label=name)
        ax.bar_label(b, labels=[f"n={int(v)}" if np.isfinite(v) else "" for v in d.n_pairs],
                     fontsize=7, padding=2)
    for name, d in matched[matched.modality == "plate reader OD600"].groupby("experiment"):
        ax.axhline(d.spearman.iloc[0], color=colors[name], ls=":", lw=1.5)
    ax.set_xticks(np.arange(len(subs)) + w / 2)
    ax.set_xticklabels(["<25", "25-45", ">=45"])
    ax.set_xlabel("pair's minimum well depth (reads)")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("technical-replicate Spearman")
    ax.set_title("Sequencing vs read depth\n(dotted = same experiment's OD600 ceiling)", fontsize=10)
    ax.legend(frameon=False, fontsize=9, loc="lower right")

    fig.suptitle("Technical-replicate reproducibility — the baseline for any cross-experiment "
                 "comparison", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "t01_replicate_reproducibility.png", dpi=160)
    plt.close(fig)


def matched_comparison():
    """OD vs sequencing on the SAME pairs, and depth stratified by the PAIR's minimum depth.

    Comparing OD reproducibility over 4000 pairs against sequencing over 300 is not a
    comparison. Restricting OD to the pairs that actually have two sequenced wells puts the two
    modalities on one footing. Depth is stratified on each pair's *minimum* well depth, because
    a replicate comparison is only as good as its shallower member -- stratifying per well left
    the low-depth bin empty, since few pairs have both wells under 20 reads.
    """
    rows = []
    for name in CONFIG:
        lay, M = load_od(name)
        w = pd.read_csv(Path(CONFIG[name]["seq"]) / "r02_well_interaction_scores.csv")
        w["pair_key"] = [frozenset((a, b)) for a, b in zip(w.strain1, w.strain2)]
        first = [sorted(k)[0] for k in w.pair_key]
        w["ra"] = np.where(w.strain1.values == np.array(first),
                           w.relative_abundance_a, 1 - w.relative_abundance_a)
        seq_pairs = w.groupby("pair_key")["ra"].agg(list)
        seq_pairs = set(seq_pairs[seq_pairs.map(len) >= 2].index)

        pl = lay.dest_plate.values
        col = M[:, 25]
        Z = np.full(len(lay), np.nan)
        for pp in np.unique(pl):
            m = pl == pp
            Z[m] = (col[m] - np.nanmean(col[m])) / np.nanstd(col[m])
        sub = lay.pair_key.isin(seq_pairs).values
        od_r = split_half(Z[sub], lay.pair_key.values[sub])
        seq_r = split_half(w.ra.values, w.pair_key.values)
        rows.append({"experiment": name, "modality": "plate reader OD600", "subset": "matched pairs",
                     **od_r})
        rows.append({"experiment": name, "modality": "sequencing rel. abundance",
                     "subset": "matched pairs", **seq_r})
        mn = w.groupby("pair_key")["n_reads"].min()
        for lo, hi, lab in [(0, 25, "pair min <25 reads"), (25, 45, "pair min 25-45"),
                            (45, 1e9, "pair min >=45")]:
            keep = set(mn[(mn >= lo) & (mn < hi)].index)
            s2 = w[w.pair_key.isin(keep)]
            rows.append({"experiment": name, "modality": "sequencing rel. abundance",
                         "subset": lab, **split_half(s2.ra.values, s2.pair_key.values)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "r03_matched_modality_comparison.csv", index=False)
    return df


def main():
    print("plate reader ...")
    od_df, spectra = od_replicates()
    print("sequencing ...")
    seq_df, _ = seq_replicates()
    matched = matched_comparison()
    make_figures(od_df, seq_df, spectra, matched)
    print("\n=== MATCHED-PAIR COMPARISON ===")
    print(matched.round(3).to_string(index=False))

    print("\n=== PLATE READER: best and conventional wavelengths ===")
    for name, d in od_df.groupby("experiment"):
        b = d.loc[d.spearman.idxmax()]
        c = d[d.wavelength_nm == 600].iloc[0]
        print(f"  {name}: best {int(b.wavelength_nm)} nm rho={b.spearman:.3f} (n={b.n_pairs}) | "
              f"600 nm rho={c.spearman:.3f} | gain {b.spearman - c.spearman:+.3f}")
    print("\n=== SEQUENCING ===")
    print(seq_df.round(3).to_string(index=False))
    return od_df, seq_df


if __name__ == "__main__":
    main()
