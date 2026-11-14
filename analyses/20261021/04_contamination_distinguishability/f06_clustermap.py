"""f06: every strain against every other, clustered -- with side bars that show where the
groups break at each level of granularity.

Built from the distances themselves. (The older clustermap in
initial_database_creation/.../combine_consensuses.ipynb used seaborn's default, which
clusters on Euclidean distances BETWEEN ROWS of the distance matrix, not on the distances.)

COMPLETE linkage, on purpose: cutting that tree at T gives groups in which EVERY pair is
below T -- mutually indistinguishable groups -- and cuts at different T are nested, so the
side bars line up with the dendrogram and each finer group sits inside a coarser one. These
groups are not identical to f03's greedy clique cover (deterministic and nested instead of
greedy), but close: < 8 bp gives 66 groups over 324 strains here vs 65 over 323 there.

Distance = haplotype-aware bp (02_distinguishability/d01_distance_matrix_bp.csv), with one
adjustment: a pair at 0 bp only through a shared operon haplotype (majority 16S differs) is
set to 0.5 bp. Unadjusted, complete linkage split 9 of the 72 byte-identical groups of f03a --
a strain matching PART of a group via one operon variant got absorbed first. With the
adjustment, byte-identical pairs are the only exact zeros, so the finest cut ("identical 16S")
reproduces f03a's groups exactly, and no coarser cut can change (0.5 crosses no threshold).
Conversely, byte-identical pairs are set to 0 even where phasing put them 1 bp apart (one
pair, I21/L19: L19's haplotypes sit 1 and 3 bp from its own consensus).

Outputs: f06_clustermap.png / .pdf (vector -- zoom for the strain labels), and
f06_cluster_order.csv (leaf order + group id at every level + genome facts).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage, optimal_leaf_ordering
from scipy.spatial.distance import squareform

import plate as P

OUT = HERE / "outputs"
BASE = HERE.parent
LEVELS = [1, 5, 8, 15, 30]                       # bars: groups with every pair < T bp;
                                                 # T = 1 is cut at 0.25 -> byte-identical only
BIN_EDGES = [0, 1, 5, 8, 15, 30, 10 ** 4]        # heatmap shades, same thresholds
SHADES = [P.RAMP[4], P.RAMP[3], P.RAMP[2], P.RAMP[1], P.RAMP[0], P.NEUTRAL]
BLOCK = [P.INK_2, P.MUTED]                       # alternating shades: a break = a change
STATUS = {"ok": ("#0ca30c", "genome on file is the well's organism"),
          "wrong": ("#d03b3b", "genome on file is another organism"),
          "unverified": ("#fab219", "genome cannot be checked"),
          "no_genome": (P.NEUTRAL, "no genome on file")}
LABEL_MIN = 6                                    # name < 30 bp... no: identical blocks this big


def block_colors(labels):
    """Clusters of >= 2 alternate between two greys in leaf order; singletons stay blank."""
    sizes = pd.Series(labels).value_counts()
    out, k, prev = [], -1, None
    for lab in labels:
        if sizes[lab] < 2:
            out.append(P.SURFACE); prev = None; continue
        if lab != prev:
            k += 1; prev = lab
        out.append(BLOCK[k % 2])
    return out


def main():
    D = pd.read_csv(BASE / "02_distinguishability/outputs/d01_distance_matrix_bp.csv", index_col=0)
    plain = pd.read_csv(BASE / "02_distinguishability/outputs/d01_distance_matrix_plain_bp.csv",
                        index_col=0).loc[D.index, D.columns]
    D = D.astype(float).mask((D == 0) & (plain > 0), 0.5).mask(plain == 0, 0.0)
    t = pd.read_csv(OUT / "p00_well_table.csv").set_index("well")
    names = list(D.index)
    cond = squareform(D.values, checks=False)
    Z = optimal_leaf_ordering(linkage(cond, "complete"), cond)
    order = dendrogram(Z, no_plot=True)["leaves"]
    leaves = [names[i] for i in order]
    M = D.loc[leaves, leaves].values

    tab = pd.DataFrame({"leaf": range(len(leaves)), "well": leaves,
                        "strain_no": t.strain_no.reindex(leaves).values})
    for T in LEVELS:
        lab = fcluster(Z, 0.25 if T == 1 else T - 0.5, criterion="distance")
        tab[f"group_lt{T}bp"] = pd.Series(lab, index=names).reindex(leaves).values
    tab["identity_group"] = t.identity_group.reindex(leaves).values
    tab["clone_cluster"] = t.clone_cluster.reindex(leaves).values
    tab["genome_status"] = t.genome_status.reindex(leaves).fillna("no_genome").values
    # renumber groups 1..n in leaf order so ids read top to bottom; singletons keep an id too
    for T in LEVELS:
        c = f"group_lt{T}bp"
        tab[c] = tab[c].map({g: k + 1 for k, g in enumerate(dict.fromkeys(tab[c]))})
    tab.to_csv(OUT / "f06_cluster_order.csv", index=False)

    bars = [(f"< {T} bp" if T > 1 else "identical 16S", block_colors(tab[f"group_lt{T}bp"].values))
            for T in LEVELS]
    clone = tab.clone_cluster.where(tab.clone_cluster.notna(), None)
    clone_lab = [c if isinstance(c, str) else f"_u{k}" for k, c in enumerate(clone)]
    bars_g = [("genome clone", block_colors(clone_lab)),
              ("genome status", [STATUS[s][0] for s in tab.genome_status])]

    # ---- layout (figure coordinates; the heatmap is square in inches) -----------------------
    W, H = 18.0, 19.2
    fig = P.plt.figure(figsize=(W, H), facecolor=P.SURFACE)
    hm_l, hm_b, hm_w = 0.215, 0.115, 0.70
    hm_h = hm_w * W / H
    ax = fig.add_axes([hm_l, hm_b, hm_w, hm_h])
    axd = fig.add_axes([hm_l, hm_b + hm_h + 0.004, hm_w, 0.135])
    bw, gap = 0.0105, 0.0022
    xs = [hm_l - gap - (len(bars) + len(bars_g) - i) * (bw + gap) - 0.008 * (i < len(bars))
          for i in range(len(bars) + len(bars_g))]

    cmap = ListedColormap(SHADES)
    norm = BoundaryNorm(BIN_EDGES, cmap.N)
    ax.imshow(M, cmap=cmap, norm=norm, interpolation="nearest", aspect="auto")
    ax.set_xticks([]); ax.yaxis.tick_right()
    ax.set_yticks(range(len(leaves)),
                  [f"{s} · {w}" for s, w in zip(tab.strain_no, tab.well)], fontsize=2.6, color=P.INK_2)
    ax.tick_params(length=0, pad=1)
    for s in ax.spines.values():
        s.set_visible(False)

    dendrogram(Z, ax=axd, no_labels=True, link_color_func=lambda k: P.INK_2,
               above_threshold_color=P.INK_2)
    for col in axd.collections:
        col.set_linewidth(0.35)
    # linear up to 30 bp, where every threshold lives; log above it, where only the shape matters
    axd.set_yscale("symlog", linthresh=30, linscale=2.5)
    axd.set_ylim(0, Z[:, 2].max() * 1.05)
    axd.set_yticks([0, 5, 8, 15, 30, 100, 300], ["0", "5", "8", "15", "30", "100", "300"])
    axd.yaxis.tick_right(); axd.yaxis.set_label_position("right")
    for T in (5, 8, 15, 30):
        axd.axhline(T - 0.5, color=P.SLOT[1], lw=0.7, ls="--", alpha=0.9)
    axd.set_ylabel("cut height, bp", fontsize=8, color=P.INK_2, rotation=270, labelpad=11)
    axd.tick_params(axis="y", labelsize=7, colors=P.MUTED, length=2)
    axd.set_xticks([])
    for s in ("top", "left", "bottom"):
        axd.spines[s].set_visible(False)
    axd.spines["right"].set_color(P.HAIRLINE)

    for x, (title, colors) in zip(xs, bars + bars_g):
        a = fig.add_axes([x, hm_b, bw, hm_h])
        a.imshow(np.array([P.matplotlib.colors.to_rgb(c) for c in colors])[:, None, :],
                 aspect="auto", interpolation="nearest")
        a.set_xticks([]); a.set_yticks([])
        for s in a.spines.values():
            s.set_visible(False)
        fig.text(x + bw / 2, hm_b + hm_h + 0.004, title, rotation=90, ha="center", va="bottom",
                 fontsize=8.5, color=P.INK)

    # name the big identical blocks (by their 16S identity group) left of the bars
    ident = tab.group_lt1bp.values
    for g in dict.fromkeys(ident):
        rows = np.where(ident == g)[0]
        if len(rows) < LABEL_MIN:
            continue
        ig = tab.identity_group.iloc[rows].dropna()
        cl = tab.clone_cluster.iloc[rows].dropna()
        lab = (ig.mode().iloc[0] if len(ig) else "") + (f" · {cl.mode().iloc[0]}" if len(cl) else "")
        y = hm_b + hm_h * (1 - (rows.mean() + 0.5) / len(leaves))
        fig.text(xs[0] - 0.004, y, f"{lab}  ({len(rows)})", ha="right", va="center", fontsize=7.6,
                 color=P.INK)

    n_at = {T: int((tab[f"group_lt{T}bp"].map(tab[f"group_lt{T}bp"].value_counts()) > 1).sum())
            for T in LEVELS}
    g_at = {T: int((tab[f"group_lt{T}bp"].value_counts() > 1).sum()) for T in LEVELS}
    fig.text(0.035, 0.975, f"{len(leaves)} strains clustered by 16S distance — where the groups break at "
                           "each threshold", fontsize=16, fontweight="bold", color=P.INK, va="top")
    fig.text(0.035, 0.957,
             "Complete linkage on haplotype-aware bp: cutting the tree at T leaves groups in which EVERY pair is "
             "< T bp (mutually indistinguishable), and the cuts are nested — each bar's\nblocks split only "
             "inside the blocks of the bar to its right. Grey blocks alternate so a change of shade is a "
             "break; blank = a strain with no partner below T. Dashed lines on the dendrogram\nare the cuts; 'identical 16S' = byte-identical majority sequence, exactly f03a's groups. "
             + "  ".join(f"< {T} bp: {g_at[T]} groups / {n_at[T]} strains" for T in LEVELS[1:])
             + f"  (identical: {g_at[1]} / {n_at[1]}).",
             fontsize=9.4, color=P.INK_2, va="top", linespacing=1.45)

    hm_items = [Patch(facecolor=c, edgecolor=P.HAIRLINE, label=l) for c, l in
                zip(SHADES, ["0 bp", "1–4 bp", "5–7 bp", "8–14 bp", "15–29 bp", "≥ 30 bp"])]
    lg1 = fig.legend(handles=hm_items, title="pair distance (heatmap)", loc="upper left",
                     bbox_to_anchor=(0.035, hm_b - 0.012), ncol=6, frameon=False, fontsize=8.6,
                     title_fontsize=8.8, handlelength=1.6, columnspacing=1.6)
    lg1._legend_box.align = "left"
    side = [Patch(facecolor=BLOCK[0], label="one group (shades alternate between neighbouring groups)"),
            Patch(facecolor=P.SURFACE, edgecolor=P.HAIRLINE, label="no partner below the threshold")]
    side += [Patch(facecolor=c, edgecolor=P.HAIRLINE, label=l) for c, l in STATUS.values()]
    lg2 = fig.legend(handles=side, title="side bars", loc="upper left",
                     bbox_to_anchor=(0.035, hm_b - 0.045), ncol=3, frameon=False, fontsize=8.6,
                     title_fontsize=8.8, handlelength=1.6, columnspacing=1.6)
    lg2._legend_box.align = "left"
    P.footnote(fig, "Row labels (right; zoom the PDF) = strain number · well, in leaf order. Left of the bars: blocks of "
                    f">= {LABEL_MIN} identical strains, named by 16S identity group (f03a) · genome clone cluster "
                    "(f04b).\nDistances: 02_distinguishability/d01_distance_matrix_bp.csv. Leaf order and every group id: "
                    "outputs/f06_cluster_order.csv.", y=0.006)
    P.save(fig, OUT, "f06_clustermap")
    return tab


if __name__ == "__main__":
    main()
