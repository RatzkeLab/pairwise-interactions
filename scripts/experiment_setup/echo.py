"""Transfer plan -> Echo CSVs, with OD-normalised volumes and source-refill batches.

The design fixes WHAT goes WHERE; volumes are only computed on the day, from that day's
source-plate OD reading, because preculture density differs run to run.

Normalisation. Composition remembers the starting ratio (20260907 ratio titration: slope
0.31) but not the absolute density (density titration: 0.04 +- 0.07 per doubling, p=0.6),
and source wells differ up to 25-fold in OD. So each member's shot aims at equal cells:

    shot_i = clip(ceil25(floor_nL * OD_ref / OD_i), floor_nL, cap_nL)

OD_i is blank-subtracted; OD_ref is a quantile of the strains used, so strains denser than
it are slightly over-dosed at the floor instead of forcing everyone else's volume up. The
floor (100 nL = 4 droplets) is the reproducibility minimum; the cap bounds the faintest
strains and the carry-over into a 50 uL well. Measured on the September reads: 90th-pct
pairwise starting imbalance 10x at equal volumes, ~4x with this scheme at cap 300 nL.

Premix (large communities). A 50-member community shot member by member is >= 5 uL of
preculture in a 50 uL well. Instead: stage 1 shoots every member, normalised, into one
well of an Echo-qualified MIX plate (scaled up so the well holds >= mix_target_uL); stage 2
shoots `shot_nL` of that mix into each replicate well. Carry-over ~2%, one transfer per
well, and the replicates receive literally the same mix.

Refills. Wells are shot plate by plate; a new source fill starts whenever continuing would
take any strain past `fill_cap_uL` within the current fill.
"""
import numpy as np
import pandas as pd

STEP = 25          # Echo 525 droplet, nL


def ceil_step(v):
    return (np.ceil(np.asarray(v, float) / STEP) * STEP).astype(int)


def transfer_plan(layout, communities, mix_wells):
    """One row per transfer, volumes left blank."""
    meta = communities.set_index("community_id")
    rows = []
    for cid, mw in mix_wells.items():
        for s in meta.at[cid, "members"].split():
            rows.append(dict(stage="premix_1", source_plate="SOURCE", source_well=s, strain=s,
                             dest_plate="MIX", dest_well=mw, community_id=cid))
    inoc = layout[(layout.well_type == "community") & layout.inoculated]
    for r in inoc.itertuples():
        if r.community_id in mix_wells:
            rows.append(dict(stage="premix_2", source_plate="MIX",
                             source_well=mix_wells[r.community_id], strain="",
                             dest_plate=str(r.plate), dest_well=r.well, community_id=r.community_id))
        else:
            for s in meta.at[r.community_id, "members"].split():
                rows.append(dict(stage="direct", source_plate="SOURCE", source_well=s, strain=s,
                                 dest_plate=str(r.plate), dest_well=r.well,
                                 community_id=r.community_id))
    return pd.DataFrame(rows)


def equal_volumes(plan, communities, member_nL, reduced_nL, max_well_uL):
    """Same volume for every member: `member_nL`, or `reduced_nL` where that would put
    more than `max_well_uL` of preculture into the well (large communities)."""
    size = communities.set_index("community_id")["size"]
    per = {}
    for cid, k in size.items():
        if k * member_nL <= max_well_uL * 1000:
            per[cid] = member_nL
        elif k * reduced_nL <= max_well_uL * 1000:
            per[cid] = reduced_nL
        else:
            raise ValueError(f"{cid}: {k} members x {reduced_nL} nL exceeds {max_well_uL} uL")
    plan = plan.copy()
    plan["volume_nL"] = plan.community_id.map(per).astype(int)
    return plan


def normalised_volumes(od, strains, blank, ref_quantile, floor_nL, cap_nL):
    x = (od.reindex(strains) - blank).clip(lower=0.01)
    ref = x.quantile(ref_quantile)
    v = pd.Series(np.clip(ceil_step(floor_nL * ref / x), floor_nL, cap_nL), index=x.index)
    return v, x, ref


def volumes(plan, od, cfg):
    """Fill in Transfer Volume for every row of the plan."""
    inoc = cfg["inoculum"]
    strains = sorted(set(plan.strain) - {""})
    v, x, ref = normalised_volumes(od, strains, inoc["blank_od"], inoc["ref_quantile"],
                                   inoc["floor_nL"], inoc["cap_nL"])
    plan = plan.copy()
    plan["volume_nL"] = 0
    d = plan.stage == "direct"
    plan.loc[d, "volume_nL"] = plan.loc[d, "strain"].map(v).values
    p1 = plan.stage == "premix_1"
    target = inoc["premix"]["mix_target_uL"] * 1000
    for cid, g in plan[p1].groupby("community_id"):
        base = g.strain.map(v).astype(float)
        scale = max(1.0, target / base.sum())
        plan.loc[g.index, "volume_nL"] = ceil_step(base * scale)
    plan.loc[plan.stage == "premix_2", "volume_nL"] = inoc["premix"]["shot_nL"]
    plan["volume_nL"] = plan.volume_nL.astype(int)
    info = pd.DataFrame({"od": od.reindex(strains), "od_minus_blank": x, "shot_nL": v})
    info["dose_vs_ref"] = (info.od_minus_blank * info.shot_nL) / (ref * inoc["floor_nL"])
    return plan, info, ref


def fills(plan, plate_order, fill_cap_uL):
    """Assign each SOURCE-drawing transfer to a source fill (refill between fills)."""
    plan = plan.copy()
    plan["fill"] = 0
    cap = fill_cap_uL * 1000
    cur = pd.Series(dtype=float)
    k = 1
    for p in plate_order:
        sel = (plan.dest_plate == p) & (plan.source_plate == "SOURCE")
        if not sel.any():
            continue
        need = plan[sel].groupby("strain").volume_nL.sum()
        if (need > cap).any():
            raise ValueError(f"plate {p} alone draws more than {fill_cap_uL} uL of "
                             f"{list(need[need > cap].index)}; split the plate or raise the cap")
        if len(cur) and (cur.add(need, fill_value=0) > cap).any():
            k += 1
            cur = pd.Series(dtype=float)
        cur = cur.add(need, fill_value=0)
        plan.loc[sel, "fill"] = k
    return plan


def write_echo(plan, outdir, prefix):
    """One CSV per (fill, destination plate) for SOURCE transfers, one per destination plate
    for MIX transfers; columns as the Echo 525 picklists used so far."""
    outdir.mkdir(parents=True, exist_ok=True)
    cols = {"source_well": "Source Well", "dest_plate": "Destination Plate Name",
            "dest_well": "Destination Well", "volume_nL": "Transfer Volume"}
    files = []
    for (src, fill, dest), g in plan.groupby(["source_plate", "fill", "dest_plate"], sort=False):
        tag = f"fill{fill}_source_to_{dest}" if src == "SOURCE" else f"mix_to_plate_{dest}"
        f = outdir / f"{prefix}_{tag}.csv"
        g[list(cols)].rename(columns=cols).to_csv(f, index=False)
        files.append((f.name, src, fill, dest, len(g), int(g.volume_nL.sum())))
    return pd.DataFrame(files, columns=["file", "source_plate", "fill", "dest_plate",
                                        "transfers", "total_nL"])


def primer_transfers(assigned, volume_nL, fwd_wells, rev_wells):
    rows = []
    for p, w, i, j in assigned:
        for src in (fwd_wells[i], rev_wells[j]):
            rows.append({"Source Well": src, "Destination Plate Name": p,
                         "Destination Well": w, "Transfer Volume": volume_nL})
    return pd.DataFrame(rows)
