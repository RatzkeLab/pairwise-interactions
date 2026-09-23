"""Read the plate-reader ASCII exports into a tidy per-well table.

The exports are a text header followed by one comma-separated grid per channel.
`ODFull_*` files are wavelength scans: their repeated grids are 61 WAVELENGTHS, not
timepoints, which is the single easiest thing to get wrong here -- treating grid 0
as "the OD" silently reads 350 nm.

Plate identity comes from the `ID1:` header field, which 20260907 fills in with the
plate number. Earlier runs left it blank and had to have plate identity reconstructed
from timing and correlation (see 20260827_analysis/plate_reader_comparison), so do not
assume this works for them.

Filenames are SAVE times; the internal `Date:`/`Time:` header is the RUN START. For a
44 h growth curve those differ by 44 h. Always use the header.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROWS16 = "ABCDEFGHIJKLMNOP"
OD600_NM = 600


def parse_header(path):
    txt = Path(path).read_text(encoding="latin-1")
    head = txt[:2000]
    def grab(pat, cast=str, default=None):
        m = re.search(pat, head)
        return cast(m.group(1).strip()) if m else default
    return dict(
        path=Path(path),
        testname=grab(r"Testname:\s*([^\n,\"]+)"),
        date=grab(r"Date:\s*([\d/]+)"),
        time=grab(r"Time:\s*([\d:]+)"),
        # ID1/ID2/ID3 share one line, so stop at the next field rather than
        # at the newline -- otherwise ID1 swallows "ID2:   ID3:"
        id1=grab(r"ID1:\s*(.*?)\s*(?:ID2:|[\n,\"])"),
        n_channels=grab(r"No\. of Channels / Multichromatics:\s*(\d+)", int, 1),
        n_intervals=grab(r"No\. of Intervals:\s*(\d+)", int, 1),
    )


def parse_grids(path):
    """Return (stack, wavelengths). stack is (n_grids, rows, cols)."""
    txt = Path(path).read_text(encoding="latin-1").splitlines()
    rows, waves = [], []
    for ln in txt:
        w = re.match(r"^\s*Wavelength:\s*(\d+)\s*nm", ln)
        if w:
            waves.append(int(w.group(1)))
        elif re.match(r"^\s*[\d.-]+\s*,", ln):
            rows.append([float(x) for x in ln.split(",") if x.strip() != ""])
    g = np.array(rows, dtype=float)
    if g.size == 0:
        return np.empty((0, 0, 0)), []
    height = 16 if g.shape[1] == 24 else 8
    n = g.shape[0] // height
    stack = g[:n * height].reshape(n, height, g.shape[1])
    return stack, waves


def read_plate(path, wavelength=OD600_NM):
    """One (16, 24) grid at the requested wavelength, or grid 0 if not a scan."""
    stack, waves = parse_grids(path)
    if stack.shape[0] == 0:
        return None, None
    if waves and wavelength in waves:
        idx = waves.index(wavelength)
    elif waves:
        idx = int(np.argmin([abs(w - wavelength) for w in waves]))
    else:
        idx = 0
    return stack[idx], (waves[idx] if waves else None)


def to_long(grid, plate, value_name="od"):
    h, w = grid.shape
    recs = []
    for r in range(h):
        for c in range(w):
            recs.append((plate, f"{ROWS16[r]}{c+1}", float(grid[r, c])))
    return pd.DataFrame(recs, columns=["dest_plate", "dest_well", value_name])


def load_folder(folder, wavelength=OD600_NM, plate_from_id1=True, testname=None):
    """Every file in `folder` -> one long table, plus the header index.

    Plates whose ID1 is a bare integer are taken as that plate number; anything else
    (blank, 'plate 1 duplicate', 'plate reader') is kept in the index with a null
    plate so a repeat read is visible rather than silently overwriting the first.
    """
    heads = [parse_header(p) for p in sorted(Path(folder).glob("*.csv"))]
    idx = pd.DataFrame(heads)
    if testname:
        idx = idx[idx.testname.astype(str).str.contains(testname, na=False)]
    idx["plate"] = pd.to_numeric(idx.id1, errors="coerce") if plate_from_id1 else np.nan
    frames = []
    for r in idx.itertuples():
        if pd.isna(r.plate):
            continue
        grid, wl = read_plate(r.path, wavelength)
        if grid is None or grid.shape != (16, 24):
            continue
        t = to_long(grid, int(r.plate))
        t["wavelength_nm"] = wl
        t["source_file"] = r.path.name
        frames.append(t)
    long = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    dup = long.duplicated(["dest_plate", "dest_well"]).sum() if len(long) else 0
    if dup:
        raise ValueError(f"{dup} duplicate plate/well rows -- two files claim the same "
                         f"plate number; resolve in the header index before using")
    return long, idx


def spectrum_table(folder, plates=None):
    """Full 61-wavelength spectrum per well: (plate, well) x wavelength."""
    heads = [parse_header(p) for p in sorted(Path(folder).glob("*.csv"))]
    out = []
    for h in heads:
        plate = pd.to_numeric(h["id1"], errors="coerce")
        if pd.isna(plate) or (plates is not None and int(plate) not in plates):
            continue
        stack, waves = parse_grids(h["path"])
        if stack.shape[0] != len(waves) or stack.shape[1:] != (16, 24):
            continue
        for wi, wl in enumerate(waves):
            t = to_long(stack[wi], int(plate))
            t["wavelength_nm"] = wl
            out.append(t)
    if not out:
        return pd.DataFrame()
    long = pd.concat(out, ignore_index=True)
    return long.pivot_table(index=["dest_plate", "dest_well"],
                            columns="wavelength_nm", values="od")
