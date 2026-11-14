"""f02: source wells that hold a SECOND organism next to the strain.

Colour = the kind of evidence that pins the second organism on the source well (three
kinds, three hues). A ring = the genome on file confirms that the dominant organism is still
the original strain -- the case "original genome still there, contaminant in the mix".
Dashed outline + x = the dominant organism is not the original either (overgrown; f01).

A dot = a second organism in ONE of the strain's two 20260907 monos. The two monos sat on
different destination plates, so a one-off may have entered there rather than in the source
well; it is shown, but not counted as source contamination.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pandas as pd

import plate as P

OUT = HERE / "outputs"
KIND = {"both_monos": 0, "run_20261026": 1, "close_haplo": 2}


def main():
    t = pd.read_csv(OUT / "p00_well_table.csv")
    sec = pd.read_csv(OUT / "p00_secondary_organisms.csv")
    src = sec[sec.evidence.isin(KIND)]
    one = set(sec[sec.evidence == "one_mono"].strain)
    cells, n = {}, {"src": 0, "ring": 0, "hatch": 0, "dot": 0}
    per_kind = {k: 0 for k in KIND}
    for r in t.itertuples():
        d = dict(top=r.strain_no)
        if not r.has_reference:
            d.update(fill=P.SURFACE, hatch="////", mid="no ref", mid_size=5.4, ink=P.MUTED)
            cells[r.well] = d; continue
        s = src[src.strain == r.well]
        if len(s):
            best = s.sort_values("evidence", key=lambda e: e.map(KIND)).iloc[0]
            k = KIND[best.evidence]
            per_kind[best.evidence] += 1
            labs = list(dict.fromkeys(s[s.evidence == best.evidence].label))
            lab = labs[0] if best.evidence == "close_haplo" else "+" + labs[0]
            if best.evidence == "both_monos":
                lab += f" {100 * s[s.evidence == 'both_monos'].frac.max():.0f}%"
            if len(labs) > 1:
                lab += f"+{len(labs) - 1}"
            d.update(fill=P.SLOT[k], edge=P.SLOT[k], mid=lab, mid_size=5.7)
            n["src"] += 1
            if r.genome_status == "ok":
                d.update(edge=P.INK, edge_w=2.0); n["ring"] += 1
            elif r.genome_status == "wrong":
                d.update(edge=P.INK, edge_w=1.6, dashed=True, mark="x", mark_color=P.SURFACE,
                         mark_edge=P.SURFACE); n["hatch"] += 1
        else:
            d.update(fill=P.NEUTRAL)
            if r.well in one:
                d.update(mark="o", mark_color=P.INK_2, mark_edge=P.INK_2); n["dot"] += 1
        cells[r.well] = d

    fig, ax = P.new_plate_figure()
    P.draw_plate(ax, cells)
    P.header(fig, f"{n['src']} source wells hold a second organism — in {n['ring']} the original strain is "
                  f"still the dominant one",
             "A second organism counts as being in the SOURCE well only if it is seen twice independently: in both "
             "20260907 monocultures (which sat on different plates), or tracking the strain across communities in a "
             "later run\nfrom the same source plate (20261026), or as a second haplotype > 10 bp from the strain's "
             "own. Ringed = the genome on file confirms the dominant organism is the original strain. Label = the "
             "second organism\n(strain number; its 16S group if several strains share that 16S; unk / U### = matches "
             "nothing in the collection; % = its share of the mono reads). Dots = seen in one mono only — "
             "may have entered on the destination plate.")
    P.legend(fig, [
        ("patch", f"same second organism in BOTH 20260907 monos  ({per_kind['both_monos']})",
         dict(fill=P.SLOT[0], edge=P.SLOT[0])),
        ("patch", f"travels with this strain across 20261026 communities (same source plate)  "
                  f"({per_kind['run_20261026']})", dict(fill=P.SLOT[1], edge=P.SLOT[1])),
        ("patch", f"two close organisms: haplotypes > 10 bp apart, label = how far  ({per_kind['close_haplo']})",
         dict(fill=P.SLOT[2], edge=P.SLOT[2])),
        ("patch", f"ring: genome on file confirms the original strain dominates  ({n['ring']})",
         dict(fill=P.SURFACE, edge=P.INK, edge_w=2.0)),
        ("patch", f"dashed + ×: the dominant organism is not the original either — overgrown, see f01  "
                  f"({n['hatch']})", dict(fill=P.SURFACE, edge=P.INK, edge_w=1.6, dashed=True)),
        ("mark", f"second organism in one mono only  ({n['dot']} wells without stronger evidence)",
         dict(mark="o", mark_color=P.INK_2, mark_edge=P.INK_2)),
        ("patch", "nothing found", dict(fill=P.NEUTRAL)),
        ("patch", "no 16S reference from 20260907", dict(fill=P.SURFACE, hatch="////")),
    ], ncol=2)
    P.footnote(fig, "Top number in each well = strain number (isolate). Mono populations: >= 8 reads and >= 3% of a "
                    "well; greedy read clustering only separates organisms more than ~100 bp apart in 16S, so closer "
                    "co-occupants appear only as haplotypes (aqua).\n"
                    "20261026 evidence: the organism's reads track one community member (association test, "
                    "community_assembly_experiments/20261026/02_qc q04_foreign_association.csv). "
                    "Full list: outputs/p00_secondary_organisms.csv.")
    P.save(fig, OUT, "f02_secondary_contamination")


if __name__ == "__main__":
    main()
