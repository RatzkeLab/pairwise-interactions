"""Same frozen stock? Per-well growth phenotype, compared across the two experiments.

If both experiments were inoculated from the same frozen 384 collection, the same well holds
the same organism in both, so per-well growth should agree. This is independent of 16S, so it
is not limited by the degeneracy that caps every sequence-based answer.

Four things decide whether a correlation here means anything, and all four are applied:

1.  **Occupied wells only.** The collection fills 298 of 384 wells. Empty wells read low in
    every plate, so including them manufactures agreement out of nothing.

2.  **Plate geometry removed.** Edge evaporation and reader gradients make ANY two plates
    correlate. Row and column means are subtracted from each plate (additive two-way
    detrending) before comparing, so what remains is well-specific.

3.  **A permutation floor.** Residuals are shuffled among occupied wells to get the
    no-shared-layout null for this exact data.

4.  **A positive-control ceiling.** Two reads known to be the same physical plate, compared the
    same way. Without it a correlation of 0.5 is uninterpretable -- media, shaking and timing
    all differ between the two experiments, so even a perfect layout match cannot reach 1.0.

Growth parameters beat single endpoints: max OD, total change, and maximum growth rate are
strain properties, whereas one timepoint mostly reports how far along the plate happened to be.

Caveat carried into the output: which physical plate each read is of, is itself reconstructed
(see FILE_CONTENTS.csv). The preculture/source-plate reads are the confident ones; whether
experiment 2's plate-reader plate holds monocultures or the pairwise experimental layout is
NOT established, so pairings involving it are reported but flagged.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
PR = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/Karl_2026")
D1, D2 = PR / "Karl_20260623_OD", PR / "Karl_20260722_OD"
MAPPING = Path("/home/rl/scripts/karl/Link to Karl/final_genomic_tables/"
               "mapping_384_well_plate_collection.csv")
N_PERM = 500


def read_cycles(path):
    txt = Path(path).read_text(encoding="latin-1").splitlines()

    def num(x):
        try:
            return float(x.strip())
        except ValueError:
            return np.nan
    rows = [[num(x) for x in l.split(",")] for l in txt if re.match(r"^\s*[\d.\-\s]+,", l)]
    rows = [r for r in rows if len(r) in (12, 24)]
    a = np.array(rows)
    h = 8 if a.shape[1] == 12 else 16
    return np.stack([a[i * h:(i + 1) * h] for i in range(a.shape[0] // h)])


def growth_params(stack):
    """max OD, change from start, and max growth rate, per well."""
    s = stack.astype(float)
    t0 = np.nanmedian(s[:max(1, len(s) // 20)], axis=0)
    mx = np.nanmax(s, axis=0)
    out = {"max_OD": mx, "delta_OD": mx - t0}
    if len(s) >= 8:
        ls = np.log(np.clip(s, 1e-3, None))
        w = max(2, len(s) // 10)
        out["max_growth_rate"] = np.nanmax(ls[w:] - ls[:-w], axis=0) / w
    return out


def detrend(x):
    x = x.astype(float)
    return x - np.nanmean(x, 1, keepdims=True) - np.nanmean(x, 0, keepdims=True) + np.nanmean(x)


def occupied():
    occ = set(pd.read_csv(MAPPING)["Well_souce_plate"])
    W = np.array([[f"{chr(65+r)}{c+1}" for c in range(24)] for r in range(16)])
    return np.isin(W, list(occ))


def assemble_384(quads):
    out = np.full((16, 24), np.nan)
    for g, (ro, co) in zip(quads, [(0, 0), (0, 1), (1, 0), (1, 1)]):
        out[ro::2, co::2] = g
    return out


def compare(a, b, mask, rng):
    da, db = detrend(a).ravel(), detrend(b).ravel()
    sel = mask.ravel() & np.isfinite(da) & np.isfinite(db)
    if sel.sum() < 30:
        return None
    x, y = da[sel], db[sel]
    rho = spearmanr(x, y)[0]
    null = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(N_PERM)])
    return {"n_wells": int(sel.sum()), "spearman_detrended": round(float(rho), 4),
            "null_mean": round(float(null.mean()), 4), "null_sd": round(float(null.std()), 4),
            "z_vs_permutation": round(float((rho - null.mean()) / null.std()), 2),
            "p_perm": float((np.abs(null) >= abs(rho)).mean())}


def main():
    rng = np.random.default_rng(0)
    mask = occupied()

    quads = [read_cycles(D1 / f)[0] for f in
             ["Karl_20260630_111648_OD600_6m.csv", "Karl_20260630_113201_OD600_6m.csv",
              "Karl_20260630_114305_OD600_6m.csv", "Karl_20260630_115357_OD600_6m.csv"]]
    plates = {
        "E1_4x96_NM_preculture_30Jun": {"max_OD": assemble_384(quads)},
        "E1_scrapped_0.01gluc_01Jul": growth_params(read_cycles(D1 / "Karl_20260701_171215_OD600_48h_384_shake.csv")),
        "E1_NM_recovery_overnight_01Jul": growth_params(read_cycles(D1 / "Karl_20260702_111117_OD600_48h_384_shake.csv")),
        "E1_NM_dense_02Jul": growth_params(read_cycles(D1 / "Karl_20260702_182146_OD600_48h_384_shake.csv")),
        "E1_NM_recovery_overnight2_02Jul": growth_params(read_cycles(D1 / "Karl_20260703_112451_OD600_48h_384_shake.csv")),
        "E2_preculture_dense_22Jul": {"max_OD": read_cycles(D2 / "Karl_20260722_172409_OD 384.csv")[0]},
        "E2_experiment_plate_22Jul": {"max_OD": read_cycles(D2 / "Karl_20260722_173208_OD 384.csv")[0]},
        "E2_growth_overnight_22Jul": growth_params(read_cycles(D2 / "Karl_20260723_132335_OD600_48h_384_noshake.csv")),
        "E2_growth_overnight_23Jul": growth_params(read_cycles(D2 / "Karl_20260724_105211_OD600_48h_384_noshake.csv")),
    }

    pairs = [
        # positive controls -- reads established as the same physical plate
        ("POS CTRL same plate", "E1_NM_recovery_overnight_01Jul", "E1_NM_recovery_overnight2_02Jul"),
        ("POS CTRL same plate", "E1_4x96_NM_preculture_30Jun", "E1_scrapped_0.01gluc_01Jul"),
        # cross-experiment: the question
        ("CROSS-EXP", "E1_4x96_NM_preculture_30Jun", "E2_preculture_dense_22Jul"),
        ("CROSS-EXP", "E1_NM_dense_02Jul", "E2_preculture_dense_22Jul"),
        ("CROSS-EXP", "E1_NM_recovery_overnight_01Jul", "E2_growth_overnight_22Jul"),
        ("CROSS-EXP", "E1_NM_recovery_overnight2_02Jul", "E2_growth_overnight_23Jul"),
        ("CROSS-EXP (flagged)", "E1_4x96_NM_preculture_30Jun", "E2_experiment_plate_22Jul"),
        # within-experiment-2 reproducibility
        ("E2 internal", "E2_growth_overnight_22Jul", "E2_growth_overnight_23Jul"),
        ("E2 internal", "E2_preculture_dense_22Jul", "E2_experiment_plate_22Jul"),
    ]

    rows = []
    for kind, a, b in pairs:
        for param in ("max_OD", "delta_OD", "max_growth_rate"):
            if param not in plates[a] or param not in plates[b]:
                continue
            r = compare(plates[a][param], plates[b][param], mask, rng)
            if r:
                rows.append({"kind": kind, "plate_a": a, "plate_b": b, "parameter": param, **r})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "g01_growth_correlations.csv", index=False)
    pd.set_option("display.width", 220)
    print(df.to_string(index=False))
    return df


if __name__ == "__main__":
    main()


# ===========================================================================
# spatial-offset control + figure
# ===========================================================================

def offset_control(P, mask, offsets=((0, 0), (0, 1), (0, 2), (1, 0), (2, 0), (1, 1),
                                     (0, -1), (-1, 0), (0, 6), (4, 0))):
    """Does the agreement survive sliding one plate by a well or two?

    Additive row/column detrending removes gradients but not smooth 2-D structure, so a
    correlation between plates that share no layout is still possible. Shared LAYOUT lives at
    exactly one alignment: shift either plate by a single well and it must vanish. Smooth
    geometry does not care about a one-well shift and survives. This is the control that
    decides whether a cross-experiment correlation means "same strains in the same wells".
    """
    def shifted(a, b, dr, dc):
        A, B = detrend(a), detrend(b)
        B = np.roll(np.roll(B, dr, axis=0), dc, axis=1)
        m = mask & np.roll(np.roll(mask, dr, 0), dc, 1) & np.isfinite(A) & np.isfinite(B)
        return spearmanr(A[m], B[m])[0] if m.sum() > 50 else np.nan

    rows = []
    for (a, b) in [("E1_NM_dense", "E2_preculture_dense"),
                   ("E1_4x96_preculture", "E2_preculture_dense"),
                   ("E1_recov_overnight", "E1_NM_dense"),
                   ("E1_4x96_preculture", "E2_experiment_plate"),
                   ("E1_recov_overnight", "E2_growth_overnight")]:
        r = {"pair": f"{a} x {b}"}
        for o in offsets:
            r[f"offset_{o[0]}_{o[1]}"] = round(float(shifted(P[a], P[b], *o)), 4)
        rows.append(r)
    return pd.DataFrame(rows)
