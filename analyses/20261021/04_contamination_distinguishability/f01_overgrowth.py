"""f01: which wells hold a different organism from the genome on file, and what took over.

Colour carries only the two big invaders (G01's organism = the C001 clone, G02's); every
other overtaken well is folded into one dark neutral and named by its label. A ring marks a
HOME well: same 16S, and the genome on file there is that organism -- so for G01 the
filled wells are where its organism spread, the rings are where it belongs. Which of the 16
home isolates was the donor cannot be said: they are one clone (>= 99.97% ANI).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd

import plate as P

OUT = HERE / "outputs"


def locality(t):
    """Mean distance (wells) from an overtaken well to the nearest home well of the
    organism it now holds, vs the same number of random wells. Returns {group: (obs, null, p)}."""
    pos = {w: (r, c) for w, r, c in zip(t.well, t.row_idx, t.col)}
    h = t[t.has_reference]
    rng = np.random.default_rng(0)

    def near(w, homes):
        return min(np.hypot(pos[w][0] - pos[x][0], pos[w][1] - pos[x][1]) for x in homes)
    res = {}
    for g in ("G01", "G02"):
        homes = list(h[(h.identity_group == g) & (h.genome_status == "ok")].well)
        lost = list(h[(h.identity_group == g) & (h.genome_status == "wrong")].well)
        pool = [w for w in h.well if w not in homes]
        obs = np.mean([near(w, homes) for w in lost])
        null = np.array([np.mean([near(w, homes) for w in rng.choice(pool, len(lost), replace=False)])
                         for _ in range(4000)])
        res[g] = (obs, null.mean(), (null <= obs).mean())
    return res


def short(label, limit=7):
    return label if len(label) <= limit else label[:limit - 1] + "…"


def main():
    t = pd.read_csv(OUT / "p00_well_table.csv")
    cells, n = {}, {}
    for r in t.itertuples():
        d = dict(top=r.strain_no)
        st = r.genome_status if r.has_reference else "no_reference"
        g = r.identity_group if isinstance(r.identity_group, str) else None
        if st == "no_reference":
            d.update(fill=P.SURFACE, hatch="////", edge=P.HAIRLINE, mid="no ref", mid_size=5.4, ink=P.MUTED)
            key = "no_ref"
        elif st == "wrong":
            if g == "G01":
                d.update(fill=P.SLOT[0], mid="G01"); key = "wrong_G01"
            elif g == "G02":
                d.update(fill=P.SLOT[1], mid="G02"); key = "wrong_G02"
            else:
                toks = r.invader_label.split(" / ")
                more = f"+{len(toks) - 1}" if len(toks) > 1 else ""
                lab = g if g else ("?" if r.invader_label == "none on file"
                                   else "→" + short(toks[0], 6) + more)
                d.update(fill=P.OTHER, mid=lab, mid_size=6.0); key = "wrong_other"
            if r.original_still_minor:
                d.update(mark="o", mark_color=P.SURFACE, mark_edge=P.SURFACE)
        elif st == "ok":
            d.update(fill=P.NEUTRAL); key = "ok"
            if g in ("G01", "G02"):
                d.update(edge=P.SLOT[0] if g == "G01" else P.SLOT[1], edge_w=2.4, mid=g)
                key = f"home_{g}"
        elif st == "unverified":
            d.update(fill=P.NEUTRAL, mid="?", ink=P.INK_2); key = "unverified"
        else:
            d.update(fill=P.SURFACE, dashed=True, edge=P.MUTED, edge_w=0.6); key = "no_genome"
        n[key] = n.get(key, 0) + 1
        cells[r.well] = d

    loc = locality(t)
    with_genome = n["ok"] + n.get("home_G01", 0) + n.get("home_G02", 0) + n["unverified"] + \
        n["wrong_G01"] + n["wrong_G02"] + n["wrong_other"]
    wrong = n["wrong_G01"] + n["wrong_G02"] + n["wrong_other"]
    fig, ax = P.new_plate_figure()
    P.draw_plate(ax, cells)
    P.header(fig, f"{wrong} of {with_genome} wells with a genome on file now hold a different organism",
             "Source plate as sequenced in 20260907. Filled = the well's 16S is not its genome's; colour = what took over. "
             "Ringed = HOME wells: same 16S, and the genome there IS that organism.\n"
             f"G01's organism is one clone (C001, 16 isolates, ≥ 99.97% ANI) — so it is known what the 9 overgrown G01 wells hold, "
             f"but not which of the 16 isolates was the donor.\n"
             f"Spread is at most weakly local: overgrown G01 wells sit {loc['G01'][0]:.1f} wells from the nearest home vs "
             f"{loc['G01'][1]:.1f} by chance (p = {loc['G01'][2]:.2f}); G02 {loc['G02'][0]:.1f} vs {loc['G02'][1]:.1f} "
             f"(p = {loc['G02'][2]:.2f}).")
    P.legend(fig, [
        ("patch", f"holds G01's organism (C001 clone) — genome on file is another organism  ({n['wrong_G01']})",
         dict(fill=P.SLOT[0], edge=P.SLOT[0])),
        ("patch", f"G01 home: genome on file is the C001 clone  ({n.get('home_G01', 0)})",
         dict(fill=P.NEUTRAL, edge=P.SLOT[0], edge_w=2.4)),
        ("patch", f"holds G02's organism (C002 clone or strain 139 — same 16S)  ({n['wrong_G02']})",
         dict(fill=P.SLOT[1], edge=P.SLOT[1])),
        ("patch", f"G02 home  ({n.get('home_G02', 0)})", dict(fill=P.NEUTRAL, edge=P.SLOT[1], edge_w=2.4)),
        ("patch", f"holds another organism: label = its 16S group, or → the genome(s) carrying its 16S "
                  f"(+n = n more), ? = none on file  ({n['wrong_other']})", dict(fill=P.OTHER, edge=P.OTHER)),
        ("mark", "the original strain is still there, as a minority (≥ 3% of a mono)",
         dict(mark="o", mark_color=P.SURFACE, mark_edge=P.INK)),
        ("patch", f"genome on file matches the well  ({n['ok']} + {n.get('home_G01', 0) + n.get('home_G02', 0)} homes)",
         dict(fill=P.NEUTRAL)),
        ("mark", f"unverified: no usable 16S in the assembly or its clones  ({n['unverified']})",
         dict(mark="$?$", size=8, mark_color=P.INK_2, mark_edge=P.INK_2)),
        ("patch", f"no genome on file  ({n['no_genome']})", dict(fill=P.SURFACE, edge=P.MUTED, dashed=True)),
        ("patch", f"no 16S reference from 20260907  ({n['no_ref']})", dict(fill=P.SURFACE, hatch="////")),
    ], ncol=2)
    P.footnote(fig, "Top number in each well = strain number (isolate). A genome is 'wrong' when its assembly carries a "
                    "clearly different 16S (≥ 1% vs ≤ 0.2% for a match), only a contaminant operon matches, a clone-mate's "
                    "16S differs, or it is < 80% ANI to every genome that does carry the well's 16S.\n"
                    "Source: analyses/20261021/03_genome_concordance (g02_strain_genome_status.csv) on the 20260907 run; "
                    "table: 04_contamination_distinguishability/outputs/p00_well_table.csv.")
    P.save(fig, OUT, "f01_overgrowth")


if __name__ == "__main__":
    main()
