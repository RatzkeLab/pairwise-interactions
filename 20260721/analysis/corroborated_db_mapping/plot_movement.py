"""Where did 20260721's labels actually come from? The scramble drawn on the plate.

Both namespaces are 384-well coordinates -- 20260721's own labels and corroborated_db's entry
names -- so every recovered call is an arrow from "where the label says it is" to "where the
organism it actually contains lives in the corroborated collection".

Drawn twice, because the two mappings answer different questions:
  trustworthy -- only the calls that stand on their own evidence (high/medium confidence).
                 Sparser, but every arrow is meant.
  forced      -- the unconstrained Hungarian bijection over all 86 labels, arrows shaded by how
                 bad the forced match is. Complete, and mostly not to be believed.

What to look for, and what has already been ruled out: `strain_identity_qc/qc_layout.py` tested
335 plate-handling hypotheses (all 24 quadrant permutations under both recombination
conventions, 180 deg rotations, row/column flips, pick-list off-by-N shifts, and compositions)
and none survives family-wise correction for 20260721 (best z=3.01, p=0.234), while 20260630 as
positive control recovers `identity` at z=8.18. So the expectation here is a diffuse field with
no coherent direction -- these figures are what that negative result looks like, not a search
for a pattern that was already rejected. The displacement panel is the direct check: a shift or
rotation would pile the arrows onto one displacement, and it does not.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "shared_pipelines"))
from genomic_ml import COLOR_BLUE, COLOR_RED, COLOR_GRID, COLOR_TEXT_SECONDARY, COLOR_GOOD

OUT = HERE / "outputs"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

N_ROWS, N_COLS = 16, 24     # A-P x 1-24


def well_xy(w):
    """'D13' -> (col=13, row_index=3). Returns (nan, nan) for anything unparseable."""
    if not isinstance(w, str) or len(w) < 2:
        return np.nan, np.nan
    r, c = w[0].upper(), w[1:]
    if not ("A" <= r <= "P") or not c.isdigit():
        return np.nan, np.nan
    col = int(c)
    if not 1 <= col <= N_COLS:
        return np.nan, np.nan
    return col, ord(r) - ord("A")


def _plate_base(ax, occupied=None):
    xs, ys = np.meshgrid(np.arange(1, N_COLS + 1), np.arange(N_ROWS))
    ax.scatter(xs, ys, s=13, facecolor="none", edgecolor=COLOR_GRID, linewidth=0.7, zorder=1)
    if occupied is not None and len(occupied):
        ox, oy = zip(*occupied)
        ax.scatter(ox, oy, s=15, color=COLOR_GRID, zorder=2)
    ax.set_xlim(0.2, N_COLS + 0.8)
    ax.set_ylim(N_ROWS - 0.5, -0.5)                     # row A at the top, as the plate is read
    ax.set_xticks([1, 4, 8, 12, 16, 20, 24])
    ax.set_yticks(range(N_ROWS))
    ax.set_yticklabels([chr(ord("A") + i) for i in range(N_ROWS)], fontsize=7)
    ax.tick_params(labelsize=7)
    ax.set_aspect("equal")
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_color(COLOR_GRID)


def _arrows(ax, src, dst, colors, lw=1.1, alpha=0.85):
    """Arrows source->target. Drawn as segments plus a head, so 86 of them stay readable."""
    segs = [[(a[0], a[1]), (b[0], b[1])] for a, b in zip(src, dst)]
    ax.add_collection(LineCollection(segs, colors=colors, linewidths=lw, alpha=alpha, zorder=3))
    for (x0, y0), (x1, y1), c in zip(src, dst, colors):
        dx, dy = x1 - x0, y1 - y0
        n = np.hypot(dx, dy)
        if n < 1e-9:
            ax.scatter([x0], [y0], s=42, facecolor="none", edgecolor=c, linewidth=1.4, zorder=5)
            continue
        ax.annotate("", xy=(x1, y1), xytext=(x1 - dx / n * 0.45, y1 - dy / n * 0.45),
                    arrowprops=dict(arrowstyle="-|>", color=c, lw=lw, shrinkA=0, shrinkB=0),
                    zorder=4)
    ax.scatter([p[0] for p in src], [p[1] for p in src], s=9, color="white",
               edgecolor=COLOR_TEXT_SECONDARY, linewidth=0.5, zorder=6)


def load():
    d = pd.read_csv(OUT / "m02_20260721_to_corroborated_db_mapping.csv")
    for col, pre in [("label_20260721", "src"), ("recommended_target", "rec"),
                     ("forced_target", "frc")]:
        xy = d[col].map(well_xy)
        d[f"{pre}_x"] = [p[0] for p in xy]
        d[f"{pre}_y"] = [p[1] for p in xy]
    return d


def figure_plate(d):
    occ = list(zip(d["src_x"], d["src_y"]))
    fig, axes = plt.subplots(1, 2, figsize=(17.5, 6.4))

    t = d.dropna(subset=["rec_x"])
    ax = axes[0]
    _plate_base(ax, occ)
    cols = [COLOR_BLUE if c == "high" else COLOR_GOOD for c in t["confidence"]]
    _arrows(ax, list(zip(t.src_x, t.src_y)), list(zip(t.rec_x, t.rec_y)), cols)
    n_gen = int(t["recommended_genome"].notna().sum())
    ax.set_title(f"Trustworthy calls only — {len(t)} of 86 labels\n"
                 f"blue = high confidence ({(t.confidence=='high').sum()}), "
                 f"green = medium ({(t.confidence=='medium').sum()}); "
                 f"{n_gen} reach the genomic tables", fontsize=10)

    ax = axes[1]
    _plate_base(ax, occ)
    f = d.dropna(subset=["frc_x"])
    norm = Normalize(vmin=0.72, vmax=1.0)
    cmap = plt.get_cmap("viridis")
    cols = [cmap(norm(v)) for v in f["forced_identity"]]
    _arrows(ax, list(zip(f.src_x, f.src_y)), list(zip(f.frc_x, f.frc_y)), cols, lw=1.0, alpha=0.8)
    ax.set_title(f"Forced bijection — all {len(f)} labels, however poor the match\n"
                 f"{(f.forced_identity < 0.90).sum()} arrows below 0.90 identity "
                 f"(min {f.forced_identity.min():.3f}) are noise, not signal", fontsize=10)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("16S identity of the forced match", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    fig.suptitle("Where 20260721's labels actually came from — arrows point to the corroborated_db "
                 "well holding the organism the label really contains", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "p01_label_movement_plate.png", dpi=160)
    plt.close(fig)


def figure_displacement(d):
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4))
    sets = [("Trustworthy", d.dropna(subset=["rec_x"]), "rec", COLOR_BLUE),
            ("Forced (all 86)", d.dropna(subset=["frc_x"]), "frc", COLOR_RED)]

    for ax, (name, sub, pre, color) in zip(axes[:2], sets):
        dx = sub[f"{pre}_x"] - sub["src_x"]
        dy = sub[f"{pre}_y"] - sub["src_y"]
        ax.axhline(0, color=COLOR_GRID, lw=1)
        ax.axvline(0, color=COLOR_GRID, lw=1)
        ax.scatter(dx + np.random.default_rng(0).normal(0, .12, len(dx)),
                   dy + np.random.default_rng(1).normal(0, .12, len(dy)),
                   s=34, color=color, alpha=.65, edgecolor="white", linewidth=.6)
        ax.scatter([0], [0], marker="x", s=90, color="black", zorder=5)
        ax.set_xlabel("column displacement"); ax.set_ylabel("row displacement")
        reuse = int(pd.Series(list(zip(dx, dy))).value_counts().max())
        ax.set_title(f"{name}  (n={len(sub)})\n"
                     f"most-reused displacement: {reuse}×\n"
                     f"a shift or rotation would pile up on one point", fontsize=9)
        ax.set_aspect("equal")

    ax = axes[2]
    t = d.dropna(subset=["rec_x"])
    dist = np.hypot(t.rec_x - t.src_x, t.rec_y - t.src_y)
    rng = np.random.default_rng(3)
    null = [np.hypot(*(np.array([rng.integers(1, N_COLS + 1) - rng.integers(1, N_COLS + 1),
                                 rng.integers(0, N_ROWS) - rng.integers(0, N_ROWS)])))
            for _ in range(20000)]
    bins = np.linspace(0, 28, 29)
    ax.hist(null, bins=bins, density=True, color=COLOR_GRID, label="random relabelling")
    ax.hist(dist, bins=bins, density=True, histtype="step", lw=2, color=COLOR_BLUE,
            label="observed (trustworthy)")
    ax.set_xlabel("distance moved on the plate (wells)")
    ax.set_ylabel("density")
    ax.set_title("Distance moved vs. a random relabelling\n"
                 "overlapping = the mix-up carries no spatial structure", fontsize=9)
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle("Does the mix-up have plate geometry? (335 layout transforms already rejected "
                 "at family-wise p=0.23)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "p02_displacement.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    d = load()
    figure_plate(d)
    figure_displacement(d)
    print("wrote:", *(p.name for p in sorted(FIG.glob("*.png"))))
