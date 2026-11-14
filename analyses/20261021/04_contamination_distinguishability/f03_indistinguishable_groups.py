"""f03a-c: which strains a nanopore 16S readout cannot tell apart, at three levels.

  f03a  identical 16S (majority sequence base for base) -- no readout separates them
  f03b  < 8 bp, haplotype-aware -- the measured limit: read cross-talk reaches its ~0.7%
        floor from 7 bp (02_distinguishability/d02)
  f03c  < 15 bp -- the original design rule

Each group is a CLIQUE (every pair inside is below the threshold), so at most one member
of a group can share a community. Connected components would chain (< 15 bp: one component
of 116 wells, members up to 40 bp apart). The cover is greedy and not unique: a dot marks a
strain that is also below the threshold to a member of ANOTHER group.

Fill = group size (one-hue ordinal ramp); bold label = group id; top = strain number.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pandas as pd

import plate as P

OUT = HERE / "outputs"
SIZE_BINS = [(2, 2), (3, 4), (5, 9), (10, 19), (20, 10 ** 6)]


def size_fill(n):
    for k, (lo, hi) in enumerate(SIZE_BINS):
        if lo <= n <= hi:
            return P.RAMP[k]


def bin_label(k):
    lo, hi = SIZE_BINS[k]
    return f"{lo}" if lo == hi else (f"{lo}+" if hi > 10 ** 5 else f"{lo}–{hi}")


def figure(t, col, outside_col, stem, title, subtitle, prefix, none_label="no other strain within reach"):
    h = t[t.has_reference]
    sizes = h[col].value_counts()
    cells = {}
    for r in t.itertuples():
        d = dict(top=r.strain_no)
        if not r.has_reference:
            d.update(fill=P.SURFACE, hatch="////", mid="no ref", mid_size=5.4, ink=P.MUTED)
        else:
            g = getattr(r, col)
            if isinstance(g, (str, float)) and pd.notna(g):
                gs = str(g) if prefix == "" else f"{prefix}{int(float(g))}"
                d.update(fill=size_fill(int(sizes[g])), mid=gs)
            else:
                d.update(fill=P.NEUTRAL)
            if outside_col and getattr(r, outside_col, 0) and getattr(r, outside_col) > 0:
                d.update(mark="o", mark_color=P.text_ink(d["fill"]), mark_edge=P.text_ink(d["fill"]))
        cells[r.well] = d
    n_grouped = int(h[col].notna().sum())
    per_bin = [0] * len(SIZE_BINS)
    for g, n in sizes.items():
        per_bin[[k for k, (lo, hi) in enumerate(SIZE_BINS) if lo <= n <= hi][0]] += 1
    fig, ax = P.new_plate_figure()
    P.draw_plate(ax, cells)
    P.header(fig, title.format(n_grouped=n_grouped, n_groups=len(sizes), total=len(h)), subtitle)
    items = [("patch", f"group of {bin_label(k)} strains  ({per_bin[k]} group{'' if per_bin[k] == 1 else 's'})",
              dict(fill=P.RAMP[k], edge=P.RAMP[k])) for k in range(len(SIZE_BINS))]
    items += [("patch", f"{none_label}  ({len(h) - n_grouped})", dict(fill=P.NEUTRAL))]
    if outside_col:
        items += [("mark", f"also within reach of a strain in another group  "
                           f"({int((h[outside_col] > 0).sum())})",
                   dict(mark="o", mark_color=P.INK_2, mark_edge=P.INK_2))]
    items += [("patch", "no 16S reference from 20260907", dict(fill=P.SURFACE, hatch="////"))]
    P.legend(fig, items, ncol=3)
    P.footnote(fig, "Top number = strain number (isolate); bold = group id. Distances: haplotype-aware bp between the "
                    "20260907 references (02_distinguishability/d01_pairs.csv).\nGroups and members: "
                    "outputs/p00_indistinguishable_groups.csv (clique groups) and "
                    "02_distinguishability/outputs/d01_identity_groups.csv (identical 16S).")
    P.save(fig, OUT, stem)


def main():
    t = pd.read_csv(OUT / "p00_well_table.csv")
    figure(t, "identity_group", None, "f03a_identical_16S",
           "{n_grouped} of {total} strains share an identical 16S with at least one other — {n_groups} groups",
           "Byte-identical majority 16S (20260907 references). No nanopore readout can tell members of a group apart: "
           "a community can hold at most one of them.\nIdentical 16S is not the same genome — see f04b for which "
           "groups are one clone on the genome side. Fill = group size; label = group id (as in "
           "d01_identity_groups.csv).", "", none_label="no other strain shares its 16S")
    figure(t, "group_lt8bp", "confusable_outside_group_lt8bp", "f03b_indistinguishable_lt8bp",
           "{n_grouped} of {total} strains fall in {n_groups} mutually indistinguishable groups at < 8 bp",
           "< 8 bp (haplotype-aware) is the MEASURED limit: on this run's monocultures, read cross-talk falls from "
           "16% (p90, 5–6 bp) to its ~0.7% floor at 7 bp and stays there.\nEach group is a clique — every pair "
           "inside is < 8 bp — so at most one member per community. The grouping is greedy and not unique; a dot "
           "marks a strain also < 8 bp from a member of another group.", "#")
    figure(t, "group_lt15bp", "confusable_outside_group_lt15bp", "f03c_indistinguishable_lt15bp",
           "{n_grouped} of {total} strains fall in {n_groups} mutually indistinguishable groups at < 15 bp",
           "< 15 bp (haplotype-aware) is the original design rule — conservative: 8–14 bp pairs separate as cleanly "
           "as 15–40 bp pairs in this run's monocultures.\nEach group is a clique — every pair inside is < 15 bp — "
           "so at most one member per community. The grouping is greedy and not unique; a dot marks a strain also "
           "< 15 bp from a member of another group.", "#")


if __name__ == "__main__":
    main()
