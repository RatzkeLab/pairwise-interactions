"""f05: the design view -- which wells can go into a community whose composition is to be
predicted from genomes (stage 10, g04_design_eligibility.csv), and why the rest cannot.

high   = genome verified, and no genomically distinct organism on file shares its 16S
medium = genome verified, but a 16S-identical organism with another genome exists, so a swap
         would be invisible to the readout
Excluded wells carry a code for the FIRST reason that applies, in this order: tier D,
never grew, mixed well, wrong genome, unverified genome, no genome.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pandas as pd

import plate as P

OUT = HERE / "outputs"
CODE = {"tier_D": ("D", "tier D: weak or conflicting 16S evidence"),
        "no_growth": ("ng", "never grew in 20260907"),
        "mixed_well": ("mix", "mixed well: two organisms"),
        "genome_wrong": ("×", "genome on file is another organism (f01)"),
        "genome_unverified": ("?", "genome cannot be checked"),
        "genome_no_genome": ("–", "no genome on file")}


def main():
    t = pd.read_csv(OUT / "p00_well_table.csv")
    cells, n = {}, {}
    for r in t.itertuples():
        d = dict(top=r.strain_no)
        if not r.has_reference:
            d.update(fill=P.SURFACE, hatch="////", mid="no ref", mid_size=5.4, ink=P.MUTED); key = "no_ref"
        elif r.genomic_use == "high":
            d.update(fill=P.RAMP[3], mid=r.clone_cluster if isinstance(r.clone_cluster, str) else None,
                     mid_size=5.6); key = "high"
        elif r.genomic_use == "medium":
            d.update(fill=P.RAMP[1], mid=r.clone_cluster if isinstance(r.clone_cluster, str) else None,
                     mid_size=5.6); key = "medium"
        else:
            code = CODE.get(r.exclude_reason, ("", r.exclude_reason))[0]
            d.update(fill=P.NEUTRAL, mid=code, ink=P.INK_2); key = r.exclude_reason
        n[key] = n.get(key, 0) + 1
        cells[r.well] = d
    use = t[t.genomic_use.isin(["high", "medium"])]
    units = use.clone_cluster.fillna(use.well).nunique()
    fig, ax = P.new_plate_figure()
    P.draw_plate(ax, cells)
    P.header(fig, f"{len(use)} wells — {units} distinct genomes — can go into genome-predicted communities",
             "Usable = tier A–C 16S reference, grew, a single organism, and the genome on file verified as that "
             "organism. Dark = the 16S also pins the genome down; light = a 16S-identical organism with a\ndifferent "
             "genome exists, so a swap would go unseen. A bold Cxxx = the well's genome is one of several identical "
             "genomes (f04b): those wells are ONE strain for modelling, and 16S-identical,\nso at most one per "
             "community. Excluded wells show the first reason that applies.")
    items = [("patch", f"usable, high confidence  ({n.get('high', 0)})", dict(fill=P.RAMP[3], edge=P.RAMP[3])),
             ("patch", f"usable, medium: 16S shared with a different genome  ({n.get('medium', 0)})",
              dict(fill=P.RAMP[1], edge=P.RAMP[1]))]
    for k, (code, text) in CODE.items():
        items.append(("patch", f"{code}  {text}  ({n.get(k, 0)})", dict(fill=P.NEUTRAL)))
    items.append(("patch", f"no 16S reference from 20260907  ({n.get('no_ref', 0)})",
                  dict(fill=P.SURFACE, hatch="////")))
    P.legend(fig, items, ncol=3)
    P.footnote(fig, "Top number = strain number (isolate). Source: 03_genome_concordance/outputs/"
                    "g04_design_eligibility.csv (genomic_use, exclude_reason, genomic_unit).")
    P.save(fig, OUT, "f05_genomic_use")


if __name__ == "__main__":
    main()
