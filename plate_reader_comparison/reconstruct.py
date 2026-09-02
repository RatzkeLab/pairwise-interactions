"""Were 20260630's and 20260721's source plates the same physical plate?

Independent of 16S -- this asks the question with growth phenotype instead, so it cannot be
defeated by the 16S degeneracy that limits every sequence-based answer (only 35 of 294 strains
have a unique 16S).

The logic: if both experiments were picked from the same frozen 384 collection, the same well
holds the same organism in both, so per-well OD should correlate across the two prep dates.
If they are different plates, it should not.

The 30/06 prep was read as FOUR 96-well plates 15 min apart (11:16, 11:32, 11:43, 11:53),
mutually uncorrelated (r -0.22..+0.12) at near-identical mean OD -- the signature of four
quadrants of the 384. They are reassembled here in the order the experimenter recalls using:
A1, A2, B1, B2.

Two controls decide whether any of this is readable:

  **Positive control (must pass first).** Correlate the reconstruction against 20260701's 384
  read -- same experiment, same plate, one day later. If the quadrant order is right this must
  be clearly positive. If it is not, the reconstruction is wrong and the cross-experiment
  comparison below means nothing, so that is checked before anything else is claimed.

  **Empty wells.** The collection occupies 298 of 384 wells. Empty wells read low in *every*
  plate, so including them manufactures correlation out of nothing. Every statistic here is
  reported both over all wells and over occupied wells only; the occupied-only number is the
  real one.

Spearman throughout: media, timepoint and growth state differ between reads, so only the
ordering of wells is comparable, not the OD scale.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

PR = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/Karl_2026")
D630, D721 = PR / "Karl_20260623_OD", PR / "Karl_20260722_OD"
MAPPING = Path("/home/rl/scripts/karl/Link to Karl/final_genomic_tables/"
               "mapping_384_well_plate_collection.csv")

# the four 96-well quadrant reads of 30/06, chronological, assigned A1 A2 B1 B2 as recalled
QUADS_630 = ["Karl_20260630_111648_OD600_6m.csv", "Karl_20260630_113201_OD600_6m.csv",
             "Karl_20260630_114305_OD600_6m.csv", "Karl_20260630_115357_OD600_6m.csv"]
QUAD_NAMES = ["A1", "A2", "B1", "B2"]


def read_plate(path, cycle=0):
    """Numeric grid for one cycle. Multi-cycle files stack cycles vertically."""
    txt = Path(path).read_text(encoding="latin-1").splitlines()
    rows = [[float(x) for x in l.split(",")]
            for l in txt if re.match(r"^\s*[\d.-]+\s*,", l)]
    g = np.array(rows)
    h = 8 if g.shape[1] == 12 else 16
    n_cycles = g.shape[0] // h
    return g[cycle * h:(cycle + 1) * h], n_cycles


def assemble_384(quad_grids, order=("A1", "A2", "B1", "B2")):
    """Interleaved 96->384: 384 well (R,C) comes from quadrant (R%2, C%2) at (R//2, C//2).

    A1 = (even row, even col), A2 = (even, odd), B1 = (odd, even), B2 = (odd, odd) -- the
    standard checkerboard a 4-channel head produces when stamping four 96s into one 384.
    """
    slot = {"A1": (0, 0), "A2": (0, 1), "B1": (1, 0), "B2": (1, 1)}
    out = np.full((16, 24), np.nan)
    for name, g in zip(order, quad_grids):
        ro, co = slot[name]
        out[ro::2, co::2] = g
    return out


def well_names():
    return np.array([[f"{chr(65+r)}{c+1}" for c in range(24)] for r in range(16)])


def occupied_mask():
    m = pd.read_csv(MAPPING)
    occ = set(m["Well_souce_plate"])
    W = well_names()
    return np.isin(W, list(occ))


def compare(a, b, mask, label_a, label_b):
    fa, fb = a.ravel(), b.ravel()
    ok = np.isfinite(fa) & np.isfinite(fb)
    rows = []
    for scope, sel in [("all_wells", ok), ("occupied_only", ok & mask.ravel())]:
        if sel.sum() < 10:
            continue
        rows.append({"comparison": f"{label_a} vs {label_b}", "scope": scope, "n": int(sel.sum()),
                     "spearman": round(float(spearmanr(fa[sel], fb[sel])[0]), 4),
                     "spearman_p": float(spearmanr(fa[sel], fb[sel])[1]),
                     "pearson": round(float(pearsonr(fa[sel], fb[sel])[0]), 4)})
    return rows


def main():
    mask = occupied_mask()
    print(f"collection occupies {mask.sum()} of 384 wells\n")

    quads = [read_plate(D630 / f)[0] for f in QUADS_630]
    recon = assemble_384(quads, order=tuple(QUAD_NAMES))
    np.savetxt(OUT / "q01_reconstructed_384_20260630.csv", recon, delimiter=",", fmt="%.4f")
    print("reconstructed 20260630 source plate from 4x96 in order "
          f"{'/'.join(QUAD_NAMES)}: mean OD {np.nanmean(recon):.3f}\n")

    rows = []
    # --- positive control: same plate, next day, 384 read
    for f in ["Karl_20260701_132952_OD 384.csv", "Karl_20260701_171954_OD600_48h_384_shake.csv"]:
        g, nc = read_plate(D630 / f)
        rows += compare(recon, g, mask, "recon_20260630", f"{f[5:18]} (384)")

    # --- the actual question: the other experiment's prep
    for f in ["Karl_20260722_172409_OD 384.csv", "Karl_20260722_173208_OD 384.csv"]:
        g, nc = read_plate(D721 / f)
        rows += compare(recon, g, mask, "recon_20260630", f"{f[5:18]} (384, 20260721 prep)")

    # --- and the two 20260721 reads against each other, as a reproducibility floor
    g1 = read_plate(D721 / "Karl_20260722_172409_OD 384.csv")[0]
    g2 = read_plate(D721 / "Karl_20260722_173208_OD 384.csv")[0]
    rows += compare(g1, g2, mask, "20260722_172409", "20260722_173208 (same plate, 8 min)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q02_correlations.csv", index=False)
    print(df.to_string(index=False))
    return df


if __name__ == "__main__":
    main()
