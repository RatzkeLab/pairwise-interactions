"""Write inspect_setup.ipynb: plots and tables over outputs/. The notebook writes nothing.

    python build_inspect_notebook.py && jupyter nbconvert --execute --inplace inspect_setup.ipynb
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
md = lambda s: nb.cells.append(nbf.v4.new_markdown_cell(s))
code = lambda s: nb.cells.append(nbf.v4.new_code_cell(s))

md("""# 20261026 setup -- inspection

Plates 1-5 are 20261021's, well for well; plate 6 is a second passage plate built the same
way as plate 5, with 31 of its 102 communities shared with plate 5 (between-plate replicates
under daily dilution). Reads `outputs/` (written by `run_setup.py`) and plots it. Nothing
here writes a file.""")

code(r'''from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
OUT = Path("outputs")
comm = pd.read_csv(OUT / "communities.csv")
lay = pd.read_csv(OUT / "layout.csv")
pool = pd.read_csv(OUT / "pool.csv")
tr = pd.read_csv(OUT / "strain_transfers.csv")
# validated categorical slots 1-2 (blue, orange); recessive ink and grid
C = {"core": "#2a78d6", "extended": "#eb6834"}
SH = {"this plate only": "#2a78d6", "on both passage plates": "#eb6834"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "font.size": 9})
print(open(OUT / "summary.txt").read())
print(open(OUT / "validation.txt").read())''')

md("""## Size x relatedness coverage

Each dot is one community of up to 10 members (larger ones are listed below).
Relatedness is the mean pairwise 16S distance between members: low = one clade, high =
across phyla. Each size is filled in three relatedness terciles, so every size spans the
range it can reach.""")

code(r'''rng = np.random.default_rng(0)
fig, ax = plt.subplots(figsize=(4.4, 3.6))
c = comm[comm.blocks.str.contains("main") & (comm["size"] <= 10)]
for p in ["core", "extended"]:
    s = c[c.pool == p]
    ax.scatter(s["size"] + rng.uniform(-0.28, 0.28, len(s)), s.mpd_bp, s=12, color=C[p],
               alpha=0.75, linewidths=0, label=f"{p} (n={len(s)})")
ax.set(title=f"main plates: {c.shape[0]} communities", xlabel="members", xticks=range(2, 11),
       ylabel="mean pairwise 16S distance (bp)")
ax.legend(frameon=False, loc="lower right", fontsize=8)
fig.tight_layout()''')

md("""## The two passage plates

Both hold core strains only, the same 7 large tests and the same size mix (12 communities
of each size 3-9, 11 of size 10). Orange: the 24 small communities on both plates (3 per
size, one from each relatedness tercile of plate 5's communities of that size). Their
positions on plate 6 are re-randomised under the same spacing rules, so a between-plate
difference and a within-plate replicate difference both include well position and can be
compared directly.""")

code(r'''fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=True)
for ax, p in zip(axes, [5, 6]):
    ids = set(lay[(lay.plate == p) & (lay.well_type == "community")].community_id)
    c = comm[comm.community_id.isin(ids) & (comm["size"] <= 10)]
    both = c.passage_plates == "5 6"
    for lab, s in [("this plate only", c[~both]), ("on both passage plates", c[both])]:
        ax.scatter(s["size"] + rng.uniform(-0.28, 0.28, len(s)), s.mpd_bp, s=14, color=SH[lab],
                   alpha=0.85, linewidths=0, label=f"{lab} (n={len(s)})")
    ax.set(title=f"passage plate {p}: {len(ids)} communities (7 large not shown)",
           xlabel="members", xticks=range(3, 11))
    ax.legend(frameon=False, loc="lower right", fontsize=8)
axes[0].set_ylabel("mean pairwise 16S distance (bp)")
fig.tight_layout()''')

code(r'''six = comm[comm.plate6_origin != ""]
print("plate 6 by origin:")
display(six.groupby("plate6_origin")["size"].value_counts().unstack(fill_value=0))
sh = comm[comm.passage_plates == "5 6"]
w = lay[lay.plate.isin([5, 6]) & lay.community_id.isin(sh.community_id)]
wells = w.groupby(["community_id", "plate"]).well.apply(lambda s: " ".join(s)).unstack()
wells.columns = [f"wells on plate {p}" for p in wells.columns]
print("communities on both passage plates (between-plate replicates):")
display(sh.set_index("community_id")[["size", "module", "mpd_bp", "blocks"]].round(0).join(wells))
print("barcode plates (same layout at every time point):")
display(lay.groupby("plate")[["culture_plate", "block", "timepoint", "inoculated"]].first())''')

code(r'''cols = ["community_id", "size", "pool", "blocks", "target_min_sep", "min_sep_bp", "mpd_bp", "member_nL", "inoculum_uL"]
print("large tests:")
display(comm[comm.module == "large_test"][cols].round(2))
print("12-20 member communities:")
display(comm[comm.module == "mid"][["community_id", "size", "pool", "mpd_bp"]].round(0))
print("resolution tests (one deliberately close pair):")
display(comm[comm.module == "resolution_test"][["community_id", "members", "close_pair", "close_pair_bp", "min_sep_bp"]])
print("trios and their pairs:")
display(comm[comm.module.isin(["trio", "trio_pair"])][["community_id", "module", "members", "parent", "pairs_measured_20260907", "mpd_bp"]].round(0))
print("sanity pairs (re-measured from 20260907):")
display(comm[comm.module == "sanity_pair"][["community_id", "members", "min_sep_bp"]])
print("leave-out series:")
display(comm[comm.module == "leave_out"][["community_id", "size", "removed"]])''')

md("""## Strain usage

Wells containing each strain on the main (sequenced 24 h) plates. Extended strains appear
less often by construction: they are confined to the extended third of the communities.""")

code(r'''m = lay[(lay.block == "main") & (lay.well_type == "community")]
inc = pd.Series([x for r in m.members for x in r.split()]).value_counts()
pp = pool.set_index("strain")["pool"]
fig, ax = plt.subplots(figsize=(11, 3))
for p in ["core", "extended"]:
    s = inc[inc.index.map(pp) == p]
    x = [list(inc.index).index(k) for k in s.index]
    ax.scatter(x, s.values, s=14, color=C[p], linewidths=0, label=f"{p} ({len(s)} strains)")
ax.set(xlabel="strains, most used first", ylabel="main-plate wells", xticks=[])
ax.legend(frameon=False)
fig.tight_layout()''')

md("""## Plate maps

Community size per well on every inoculated plate (blank-medium wells marked B).
Replicates of a community sit on different plates, except on the passage plates (5, 6),
where all three are on the plate but never in the same row or next to each other. Orange
rings on plates 5 and 6: communities on both passage plates.""")

code(r'''plates_ = sorted(lay[lay.inoculated].plate.unique())
ROWS = "ABCDEFGHIJKLMNOP"
cmap = ListedColormap(["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
shared = set(comm.community_id[comm.passage_plates == "5 6"])
fig, axes = plt.subplots(2, 3, figsize=(13, 6))
for ax, p in zip(axes.flat, plates_):
    g = lay[lay.plate == p]
    M = np.full((16, 24), np.nan)
    for r in g.itertuples():
        M[ROWS.index(r.well[0]), int(r.well[1:]) - 1] = min(r.size, 12) if r.well_type == "community" else np.nan
    im = ax.imshow(M, cmap=cmap, vmin=2, vmax=12, aspect="auto")
    for r in g[g.well_type == "blank"].itertuples():
        ax.text(int(r.well[1:]) - 1, ROWS.index(r.well[0]), "B", ha="center", va="center", fontsize=7, color=INK)
    big = g[g["size"] > 12]
    for r in big.itertuples():
        ax.text(int(r.well[1:]) - 1, ROWS.index(r.well[0]), str(r.size), ha="center", va="center", fontsize=6, color="#fcfcfb")
    if p in (5, 6):
        s = g[g.community_id.isin(shared)]
        ax.scatter(s.col - 1, s.well.str[0].map(ROWS.index), s=20, facecolors="none",
                   edgecolors=SH["on both passage plates"], linewidths=1.1)
    blk = g.block.iloc[0]
    ax.set(title=f"plate {p} ({blk})", xticks=range(0, 24, 4), xticklabels=range(1, 25, 4),
           yticks=range(16), yticklabels=list(ROWS))
    ax.grid(False)
fig.subplots_adjust(hspace=0.45, wspace=0.25, right=0.86)
cb = fig.colorbar(im, ax=list(axes.flat[:len(plates_)]), shrink=0.6, label="members (above 12: size printed)")
plt.show()''')

md("""## Inoculum

Equal volume per member (100 nL; 75 nL in communities above 50 members, so no well gets
more than 5 uL = 10% of the culture). Left: preculture carried into each well. Right:
volume drawn from each source well over all six plates -- all under the 25 uL that a fill
allows, so no refill.""")

code(r'''fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
w = tr.groupby(["dest_plate", "dest_well"]).volume_nL.sum() / 1000
ax[0].hist(w, bins=np.linspace(0, 5.2, 53), color=C["core"], rwidth=0.85)
ax[0].set(xlabel="inoculum per well (uL)", ylabel="wells", yscale="log")
d = (tr.groupby("strain").volume_nL.sum() / 1000).sort_values(ascending=False)
ax[1].scatter(range(len(d)), d.values, s=12, color=C["core"], linewidths=0)
ax[1].axhline(25, color=MUTED, lw=1, ls="--"); ax[1].text(len(d) - 1, 25.5, "fill limit", color=MUTED, fontsize=8, ha="right")
ax[1].set(xlabel="source wells, most drawn first", ylabel="drawn (uL)", ylim=(0, 28), xticks=[])
fig.tight_layout()''')

nbf.write(nb, "inspect_setup.ipynb")
print("wrote inspect_setup.ipynb")
