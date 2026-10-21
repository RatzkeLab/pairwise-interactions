"""Write inspect_setup.ipynb: plots and tables over outputs/. The notebook writes nothing.

    python build_inspect_notebook.py && jupyter nbconvert --execute --inplace inspect_setup.ipynb
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
md = lambda s: nb.cells.append(nbf.v4.new_markdown_cell(s))
code = lambda s: nb.cells.append(nbf.v4.new_code_cell(s))

md("""# Community experiment setup -- inspection

Reads `outputs/` (written by `run_setup.py`) and plots it. Nothing here writes a file:
the design is changed in `setup_config.yaml` and rebuilt with `python run_setup.py --force`.""")

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

code(r'''fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
rng = np.random.default_rng(0)
axes[2].remove()
for ax, blk in zip(axes, ["main", "passage"]):
    c = comm[comm.blocks.str.contains(blk) & (comm["size"] <= 10)]
    for p in ["core", "extended"]:
        s = c[c.pool == p]
        ax.scatter(s["size"] + rng.uniform(-0.28, 0.28, len(s)), s.mpd_bp, s=12, color=C[p],
                   alpha=0.75, linewidths=0, label=f"{p} (n={len(s)})")
    ax.set(title=f"{blk}: {c.shape[0]} communities", xlabel="members", xticks=range(2, 11))
    ax.legend(frameon=False, loc="lower right", fontsize=8)
axes[0].set_ylabel("mean pairwise 16S distance (bp)")
fig.tight_layout()''')

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
Replicates of a community sit on different plates, except on the passage plate (5), where
all three are on the plate but never in the same row or next to each other.""")

code(r'''plates_ = sorted(lay[lay.inoculated].plate.unique())
ROWS = "ABCDEFGHIJKLMNOP"
cmap = ListedColormap(["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
fig, axes = plt.subplots(2, 3, figsize=(13, 6))
axes.flat[-1].remove()
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
volume drawn from each source well over the whole run -- all under the 25 uL that a fill
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
