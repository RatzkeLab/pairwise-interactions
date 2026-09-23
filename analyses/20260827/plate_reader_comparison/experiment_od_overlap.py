"""Do the strain pairs shared by both experiments grow the same way in each?

The source plates are the same (see growth_compare.py), the strains were shot to protocol and
the barcodes were shot to protocol. So if 20260721's wells really hold what the layout says,
then a strain pair appearing in BOTH experiments should reach a similar OD in both.

Why this test is clean:

  The shared pairs sit at **different plate positions in the two experiments** -- different
  destination plate, different well. So a correlation here cannot be edge effects, evaporation
  or reader gradients, the confounds that dogged the source-plate comparison. Position is
  matched only by strain content, and nothing else links the two measurements.

  It is also independent of 16S, so the degeneracy that caps every sequence-based answer does
  not apply.

Reading the result:
  correlation ~ the within-experiment replicate ceiling -> the pairs really are the same
                organisms; 20260721's wells contain what the layout claims
  correlation ~ 0                                       -> they do not

The within-experiment ceiling is essential: OD of a two-strain co-culture is noisy, so the
honest question is not "is it 1.0" but "does it reach what a replicate of the SAME pair in the
SAME experiment reaches".

OD600 is wavelength index 25 of the 61-point spectrum (350 nm + 10 nm steps). Plate identity
comes from the `ID1` header field, not file order -- exp1's files are out of order at the end
and plate 2 is missing entirely (plate 1 was read twice instead).
"""

import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
EXPS = HERE.parent
PR = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/Karl_2026")

OD600_INDEX = 25          # (600 - 350) / 10
N_PERM = 2000

CONFIG = {
    "20260630": (PR / "Karl_20260704_OD_Full", EXPS / "20260630/setup/strain_layout_20260630.csv"),
    "20260721": (PR / "Karl_20260723_OD_Full",
                 EXPS / "20260721/setup/strain_layout_20260721_plate1_2_swapped.csv"),
}


def read_od600(path):
    txt = Path(path).read_text(encoding="latin-1").splitlines()
    plate = next((l.split("ID1:")[1].split("ID2")[0].strip() for l in txt[:6] if "ID1:" in l), None)

    def num(x):
        try:
            return float(x.strip())
        except ValueError:
            return np.nan
    rows = [[num(v) for v in l.split(",")] for l in txt if re.match(r"^\s*[\d.\-\s]+,", l)]
    rows = [r for r in rows if len(r) == 24]
    a = np.array(rows)
    if a.shape[0] < (OD600_INDEX + 1) * 16:
        return plate, None
    return plate, a[OD600_INDEX * 16:(OD600_INDEX + 1) * 16]


def load_experiment(name):
    folder, layout_csv = CONFIG[name]
    grids = {}
    for f in sorted(glob.glob(str(folder / "*.csv"))):
        p, g = read_od600(f)
        if g is None or not p or not p.isdigit():
            continue
        grids.setdefault(int(p), []).append(g)          # plate 1 read twice in exp1 -> average
    grids = {k: np.nanmean(v, axis=0) for k, v in grids.items()}

    lay = pd.read_csv(layout_csv)
    od = []
    for r in lay.itertuples():
        g = grids.get(int(r.dest_plate))
        if g is None:
            od.append(np.nan)
            continue
        row = ord(str(r.dest_row).upper()[0]) - 65
        od.append(g[row, int(r.dest_col) - 1] if 0 <= row < 16 and 1 <= int(r.dest_col) <= 24
                  else np.nan)
    lay["od600"] = od
    # plate-level batch differences are not biology -- z-score within plate before comparing
    lay["od_z"] = lay.groupby("dest_plate")["od600"].transform(
        lambda s: (s - s.mean()) / s.std(ddof=0))
    lay["pair_key"] = [frozenset((a, b)) for a, b in zip(lay.strain1, lay.strain2)]
    lay["experiment"] = name
    print(f"  {name}: {len(grids)} plates with OD, "
          f"{lay.od600.notna().sum()}/{len(lay)} layout wells matched to a reading")
    missing = sorted(set(lay.dest_plate) - set(grids))
    if missing:
        print(f"    NO OD for plate(s): {missing}")
    return lay


def main():
    print("loading OD600 (wavelength index 25 of 61) ...")
    L = {k: load_experiment(k) for k in CONFIG}

    # ---- within-experiment replicate ceiling: same pair, different wells, same experiment
    ceil = []
    for name, lay in L.items():
        d = lay.dropna(subset=["od_z"])
        d = d[d.well_type == "pair"]
        g = d.groupby("pair_key")["od_z"].agg(list)
        g = g[g.map(len) >= 2]
        x = [v[0] for v in g]
        y = [np.mean(v[1:]) for v in g]
        rho = spearmanr(x, y)[0]
        ceil.append({"scope": f"{name} replicate ceiling", "n_pairs": len(g),
                     "spearman": round(float(rho), 4)})
        print(f"  {name}: replicate ceiling rho={rho:.3f} over {len(g)} pairs with >=2 wells")

    # ---- cross-experiment on the shared pairs
    a, b = L["20260630"], L["20260721"]
    sa = a.dropna(subset=["od_z"]).groupby("pair_key")["od_z"].mean()
    sb = b.dropna(subset=["od_z"]).groupby("pair_key")["od_z"].mean()
    shared = sorted(set(sa.index) & set(sb.index), key=lambda k: sorted(k))
    x = sa.loc[shared].values
    y = sb.loc[shared].values
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    rho = spearmanr(x, y)[0]
    rng = np.random.default_rng(0)
    null = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(N_PERM)])
    z = (rho - null.mean()) / null.std()
    p = float((np.abs(null) >= abs(rho)).mean())

    res = pd.DataFrame(ceil + [{
        "scope": "CROSS-EXPERIMENT shared pairs", "n_pairs": int(ok.sum()),
        "spearman": round(float(rho), 4), "null_mean": round(float(null.mean()), 4),
        "null_sd": round(float(null.std()), 4), "z": round(float(z), 2), "p_perm": p}])
    res.to_csv(OUT / "e01_overlap_summary.csv", index=False)

    tbl = pd.DataFrame({
        "pair": ["|".join(sorted(k)) for k, o in zip(shared, ok) if o],
        "od_z_20260630": x, "od_z_20260721": y})
    tbl.to_csv(OUT / "e02_shared_pair_od.csv", index=False)

    print()
    print(res.to_string(index=False))
    return res


if __name__ == "__main__":
    main()
