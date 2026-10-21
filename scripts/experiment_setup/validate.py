"""Checks run on every setup output. Each returns a list of failure messages (empty = pass).

Each check exists because the thing it checks went wrong before, or would silently cost a
whole experiment if it did: a minibar index carrying adapter + primer (three runs, ~5x
reads lost), barcodes reused, communities whose members cannot be told apart, replicates
sharing a plate, a source well drawn dry.
"""
import itertools
import os

import numpy as np
import pandas as pd

from . import plates
from .barcodes import BARCODE_LEN, FWD_PRIMER, REV_PRIMER, min_edit_distance


def check_communities(comm, pool, D, default_sep):
    errs = []
    p = pool.set_index("strain")
    for r in comm.itertuples():
        m = r.members.split()
        if len(m) != r.size or len(set(m)) != len(m):
            errs.append(f"{r.community_id}: size/member mismatch")
        bad = [x for x in m if x not in p.index]
        if bad:
            errs.append(f"{r.community_id}: strains not in the pool {bad}")
            continue
        if r.pool == "core" and (p.loc[m, "pool"] == "extended").any():
            errs.append(f"{r.community_id}: core community holds an extended strain")
        need = r.target_min_sep if not pd.isna(getattr(r, "target_min_sep", np.nan)) else default_sep
        allowed = set()
        if isinstance(getattr(r, "close_pair", None), str):
            allowed = {frozenset(r.close_pair.split("-"))}
        if r.module in ("sanity_pair", "trio_pair"):
            need = 15                      # pairs: resolvable is enough, as in 20260907
        for a, b in itertools.combinations(m, 2):
            if D.at[a, b] < need and frozenset((a, b)) not in allowed:
                errs.append(f"{r.community_id}: {a}-{b} only {D.at[a, b]:.0f} bp apart (< {need})")
    if comm.members.duplicated().any():
        errs.append(f"{comm.members.duplicated().sum()} duplicate communities")
    return errs


def check_layout(layout, n_reps, within_plate_blocks, blanks_per_plate):
    errs = []
    inner = set(plates.inner_wells())
    if not set(layout.well) <= inner:
        errs.append("edge wells used")
    if layout.duplicated(["plate", "well"]).any():
        errs.append("a well is used twice")
    com = layout[layout.well_type == "community"]
    for (blk, cid), g in com.groupby(["block", "community_id"]):
        if len(g) != n_reps:
            errs.append(f"{blk}/{cid}: {len(g)} replicates, expected {n_reps}")
        if blk in within_plate_blocks:
            rows = list(g.well.str[0])
            if len(set(rows)) != len(rows):
                errs.append(f"{blk}/{cid}: replicates share a row")
            for a, b in itertools.combinations(g.well, 2):
                if b in plates.neighbours(a):
                    errs.append(f"{blk}/{cid}: replicates {a},{b} are neighbours")
        elif g.plate.nunique() != len(g):
            errs.append(f"{blk}/{cid}: replicates share a plate")
    nb = layout[layout.well_type == "blank"].groupby("plate").size()
    for p in layout.plate.unique():
        if nb.get(p, 0) != blanks_per_plate:
            errs.append(f"plate {p}: {nb.get(p, 0)} blanks, expected {blanks_per_plate}")
    return errs


def check_barcodes(assigned, failed):
    errs = []
    combos = [(i, j) for _, _, i, j in assigned]
    if len(set(combos)) != len(combos):
        errs.append("a barcode combination is used twice")
    if any(i == j for i, j in combos):
        errs.append("equal forward/reverse index used")
    if any(i in failed["fwd"] or j in failed["rev"] for i, j in combos):
        errs.append("a barcode that failed in 20260907 is used")
    if len({(p, w) for p, w, _, _ in assigned}) != len(assigned):
        errs.append("a sample (plate, well) has two barcode combinations")
    return errs


def check_minibar(sheet, planned_e):
    """The 20260630-20260907 bug, and the -e / barcode-distance relationship."""
    errs, info = [], {}
    for col, primer in (("FwIndex", FWD_PRIMER), ("RvIndex", REV_PRIMER)):
        lens = sheet[col].str.len().unique()
        if list(lens) != [BARCODE_LEN]:
            errs.append(f"{col}: lengths {sorted(lens)}, expected only {BARCODE_LEN}")
        if sheet[col].str.endswith(primer[-8:]).any():
            errs.append(f"{col}: index ends in the 16S primer -- construct, not barcode")
        pre = os.path.commonprefix(list(sheet[col].unique()))
        if len(pre) >= 4:
            errs.append(f"{col}: every index starts with '{pre}' -- an adapter got in")
        info[f"{col}_min_edit_distance"] = min_edit_distance(sheet[col].unique())
    dmin = min(info.values())
    if 2 * planned_e >= dmin:
        errs.append(f"-e {planned_e} too permissive for barcode set with min distance {dmin}")
    if sheet.SampleID.duplicated().any():
        errs.append("duplicate SampleID")
    info["planned_e"] = planned_e
    return errs, info


def check_echo(plan, inoc, min_nL):
    errs = []
    if (plan.volume_nL % 25).any():
        errs.append("a transfer volume is not a multiple of 25 nL")
    d = plan[plan.stage == "direct"]
    if (d.volume_nL < min_nL).any():
        errs.append(f"a direct shot is below the {min_nL} nL reproducibility floor")
    per_well = plan[plan.stage.isin(["direct", "premix_2"])].groupby(
        ["dest_plate", "dest_well"]).volume_nL.sum() / 1000
    if (per_well > inoc["max_well_inoculum_uL"]).any():
        errs.append(f"{(per_well > inoc['max_well_inoculum_uL']).sum()} wells receive more than "
                    f"{inoc['max_well_inoculum_uL']} uL of inoculum")
    draw = plan[plan.source_plate == "SOURCE"].groupby(["fill", "strain"]).volume_nL.sum() / 1000
    if (draw > inoc["fill_cap_uL"]).any():
        errs.append("a strain is drawn past the per-fill cap")
    mix = plan[plan.stage == "premix_1"].groupby("dest_well").volume_nL.sum() / 1000
    if (mix > 60).any():
        errs.append("a MIX well would exceed 60 uL")
    return errs
