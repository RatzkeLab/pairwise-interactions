"""f04a: isolates placed in the collection more than once ON PURPOSE -- do the copies still
hold the isolate?   f04b: genomes on file that are identical (clone clusters) -- mostly
separate isolates of one organism, not deliberate copies.

f04a. 9 isolates sit in 21 wells (strain-ID reconciliation, analyses/20261113 s00). Only
one copy of each (the K/M-row well) has a genome joined. A copy HOLDS the isolate if its
16S is carried by the isolate's genome (blast, stage 9) or is identical to a copy whose
genome verifies. Connectors join the copies of one isolate.

f04b. 298 assemblies collapse into 181 organisms; 162 sit in 45 clone clusters at
>= 99.9% ANI. None of these clusters is a deliberate copy (the duplicated isolates have one
genome each), so each is several isolates of the same organism. A dashed well with x is a
member whose WELL no longer holds that organism (stage 9 `wrong`).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pandas as pd

import plate as P

OUT = HERE / "outputs"
BASE = HERE.parent
STRAIN_IDS = BASE.parents[0] / "20261113/messy/outputs/s00_strain_ids.csv"
CLONE_BINS = [(2, 2), (3, 4), (5, 9), (10, 19)]


def duplicates():
    ids = pd.read_csv(STRAIN_IDS)
    t = pd.read_csv(OUT / "p00_well_table.csv").set_index("well")
    gs = pd.read_csv(BASE / "03_genome_concordance/outputs/g02_strain_genome_status.csv").set_index("strain")
    M = pd.read_csv(BASE / "02_distinguishability/outputs/d01_distance_matrix_plain_bp.csv", index_col=0)
    rows = []
    for iso, g in ids[ids.n_wells_same_isolate > 1].groupby("isolate_name"):
        gw = g[g.genome_table_id.notna()].well.iloc[0]
        gst = t.at[gw, "genome_status"]
        for r in g.itertuples():
            w = r.well
            carriers = str(gs.genomes_carrying_16S.get(w, "")).split()
            same_as_genome_copy = w in M.index and gw in M.index and M.at[w, gw] == 0
            if w == gw:
                state = {"ok": "holds", "wrong": "other"}.get(gst, "unverified")
            elif gw in carriers or (gst == "ok" and same_as_genome_copy):
                state = "holds"
            elif gst == "unverified" and same_as_genome_copy:
                state = "unverified"
            else:
                state = "other"
            rows.append(dict(isolate=iso, well=w, stock=r.stock_no, strain_no=r.strain_no,
                             genome_copy=w == gw, state=state,
                             bp_to_genome_copy=int(M.at[w, gw]) if w in M.index and gw in M.index else None))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "f04a_duplicate_isolates.csv", index=False)
    return d, t


def fig_a():
    d, t = duplicates()
    fill = {"holds": P.SLOT[0], "other": P.SLOT[1], "unverified": P.SLOT[2]}
    cells = {}
    for w, r in t.iterrows():
        cells[w] = dict(top=r.strain_no, fill=P.NEUTRAL if r.has_reference else P.SURFACE,
                        ink=P.MUTED, hatch=None if r.has_reference else "////")
    for r in d.itertuples():
        c = dict(top=r.stock, mid=r.isolate, mid_size=5.9, fill=fill[r.state], edge=fill[r.state])
        if r.genome_copy:
            c.update(edge=P.INK, edge_w=2.2)
        cells[r.well] = c
    n = d.state.value_counts()
    lost = sorted(d[d.state == "other"].isolate.unique())
    fig, ax = P.new_plate_figure()
    P.draw_plate(ax, cells)
    for iso, g in d.groupby("isolate"):
        P.connect(ax, list(g.well))
    both_lost = [i for i in lost if (d[d.isolate == i].state == "other").all()]
    P.header(fig, f"9 isolates are in the collection on purpose more than once (21 wells) — "
                  f"{n.get('holds', 0)} copies still hold the isolate",
             "Connected wells are copies of one isolate (strain-ID reconciliation, analyses/20261113). Only one copy per "
             "isolate has a genome on file (ringed). A copy holds the isolate when its 20260907 16S is carried by that "
             "genome.\n"
             f"Copies that hold something else: {', '.join(f'{r.well} ({r.isolate})' for r in d[d.state == 'other'].itertuples())}"
             + (f" — for {', '.join(both_lost)} NO copy matches the genome." if both_lost else ".")
             + "\nTop label = the stock number of that copy (the 384 sheet's Strain.Number); bold = the isolate.")
    P.legend(fig, [
        ("patch", f"copy holds the isolate — its 16S is in the isolate's genome  ({n.get('holds', 0)})",
         dict(fill=P.SLOT[0], edge=P.SLOT[0])),
        ("patch", f"copy holds another organism  ({n.get('other', 0)})", dict(fill=P.SLOT[1], edge=P.SLOT[1])),
        ("patch", f"copies agree with each other, but the genome has no usable 16S to check  "
                  f"({n.get('unverified', 0)})", dict(fill=P.SLOT[2], edge=P.SLOT[2])),
        ("patch", "ring: the copy whose genome is on file", dict(fill=P.SURFACE, edge=P.INK, edge_w=2.2)),
        ("patch", "not a duplicated isolate (top = strain number)", dict(fill=P.NEUTRAL)),
    ], ncol=2)
    P.footnote(fig, "Table: outputs/f04a_duplicate_isolates.csv. Copy identity from "
                    "analyses/20261113/messy/outputs/s00_strain_ids.csv; 16S-in-genome from "
                    "03_genome_concordance (g02_strain_genome_status.csv).")
    P.save(fig, OUT, "f04a_duplicate_isolates")


def fig_b():
    t = pd.read_csv(OUT / "p00_well_table.csv")
    cl = pd.read_csv(BASE / "03_genome_concordance/outputs/g02_genome_clusters.csv").set_index("well")
    cells, per_bin, n_wrong = {}, [0] * len(CLONE_BINS), 0
    sizes = cl.dropna(subset=["clone_cluster"]).groupby("clone_cluster").size()
    for k, (lo, hi) in enumerate(CLONE_BINS):
        per_bin[k] = int(((sizes >= lo) & (sizes <= hi)).sum())
    for r in t.itertuples():
        d = dict(top=r.strain_no)
        if r.well not in cl.index:
            d.update(fill=P.SURFACE, dashed=True, edge=P.MUTED)
        else:
            c = cl.at[r.well, "clone_cluster"]
            if isinstance(c, str):
                k = [i for i, (lo, hi) in enumerate(CLONE_BINS) if lo <= sizes[c] <= hi][0]
                d.update(fill=P.RAMP[k + 1], mid=c, mid_size=6.0)
                if r.genome_status == "wrong":
                    d.update(edge=P.INK, edge_w=1.6, dashed=True, mark="x",
                             mark_color=P.text_ink(d["fill"]), mark_edge=P.text_ink(d["fill"]))
                    n_wrong += 1
            else:
                d.update(fill=P.NEUTRAL)
        cells[r.well] = d
    n_cl = int(cl.clone_cluster.notna().sum())
    fig, ax = P.new_plate_figure()
    P.draw_plate(ax, cells)
    P.header(fig, f"The 298 genomes on file are {len(cl) - n_cl + cl.clone_cluster.nunique()} distinct organisms: "
                  f"{n_cl} assemblies sit in {cl.clone_cluster.nunique()} clone clusters",
             "Clone cluster = genomes ≥ 99.9% ANI (skani; min 99.97% within a cluster, alignment ≥ 94%). These are "
             "separate isolates of one organism, not the deliberate copies of f04a (each of those has only one "
             "genome).\nFor genome-based models a cluster is ONE strain: cross-validation must hold its members out "
             f"together. Dashed + × ({n_wrong}) = a member whose well now holds a different organism (f01) — its "
             "genome is in the cluster, its contents are not.\nC001 is G01's organism (16 isolates).")
    items = [("patch", f"clone cluster of {lo if lo == hi else f'{lo}–{hi}'} genomes  ({per_bin[k]} clusters)",
              dict(fill=P.RAMP[k + 1], edge=P.RAMP[k + 1])) for k, (lo, hi) in enumerate(CLONE_BINS)]
    items += [("patch", f"genome on file, no clone  ({len(cl) - n_cl})", dict(fill=P.NEUTRAL)),
              ("patch", "no genome on file", dict(fill=P.SURFACE, edge=P.MUTED, dashed=True)),
              ("patch", f"member whose well holds another organism  ({n_wrong})",
               dict(fill=P.SURFACE, edge=P.INK, edge_w=1.6, dashed=True))]
    P.legend(fig, items, ncol=3)
    P.footnote(fig, "Top number = strain number (isolate). Clusters: 03_genome_concordance/outputs/"
                    "g02_genome_clusters.csv (skani ANI over the panX input GenBank assemblies).")
    P.save(fig, OUT, "f04b_genome_clone_clusters")


if __name__ == "__main__":
    fig_a()
    fig_b()
