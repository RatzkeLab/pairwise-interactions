"""Communities x replicates -> plates and wells."""
import numpy as np
import pandas as pd

from . import plates


def across_plates(communities, plate_ids, n_reps, rng, blanks_per_plate, block):
    """Each community's replicates on DIFFERENT plates -- the only arrangement that
    measures plate-to-plate rather than well-to-well reproducibility. Larger communities
    are placed first so every plate gets a similar size mix; positions are random."""
    wells = plates.inner_wells()
    cap = {p: len(wells) - blanks_per_plate for p in plate_ids}
    if sum(cap.values()) < n_reps * len(communities):
        raise ValueError(f"{block}: {len(communities)} communities x {n_reps} do not fit "
                         f"{sum(cap.values())} wells")
    order = sorted(communities, key=lambda c: (-c["size"], rng.random()))
    content = {p: [] for p in plate_ids}
    for c in order:
        avail = sorted((p for p in plate_ids if cap[p] > 0), key=lambda p: (-cap[p], rng.random()))
        if len(avail) < n_reps:
            raise ValueError(f"{block}: no {n_reps} distinct plates left for {c['community_id']}")
        for rep, p in enumerate(avail[:n_reps], start=1):
            content[p].append((c["community_id"], rep)); cap[p] -= 1
    rows = []
    for p in plate_ids:
        items = content[p] + [("BLANK", 0)] * blanks_per_plate
        pos = rng.permutation(wells)[:len(items)]
        for (cid, rep), w in zip(items, pos):
            rows.append(dict(plate=p, well=w, block=block, community_id=cid, replicate=rep))
    return rows


def within_plate(communities, plate_id, n_reps, rng, blanks, block):
    """All replicates on ONE plate (the passage plate, which is diluted daily as a unit).
    Replicates of a community never share a row and are never neighbours, so a splash or
    a row-wise contamination path (rows were enriched 1.13x in the collection, columns
    were not) cannot hit two replicates of one community."""
    wells = plates.inner_wells()
    need = n_reps * len(communities) + blanks
    if need > len(wells):
        raise ValueError(f"{block}: {need} wells needed, {len(wells)} available")
    order = sorted(communities, key=lambda c: (-c["size"], rng.random()))
    for _ in range(500):
        free = list(rng.permutation(wells))
        rows, ok = [], True
        for c in order:
            picked = []
            for w in free:
                if all(w[0] != x[0] and w not in plates.neighbours(x) for x in picked):
                    picked.append(w)
                    if len(picked) == n_reps:
                        break
            if len(picked) < n_reps:
                ok = False
                break
            for rep, w in enumerate(picked, start=1):
                rows.append(dict(plate=plate_id, well=w, block=block,
                                 community_id=c["community_id"], replicate=rep))
                free.remove(w)
        if ok:
            for w in free[:blanks]:
                rows.append(dict(plate=plate_id, well=w, block=block,
                                 community_id="BLANK", replicate=0))
            return rows
    raise RuntimeError(f"{block}: could not satisfy the replicate-spacing constraint")


def to_frame(rows, communities):
    df = pd.DataFrame(rows)
    meta = communities.set_index("community_id")
    df["well_type"] = np.where(df.community_id == "BLANK", "blank", "community")
    for col in ("size", "members", "module", "pool", "mpd_bp", "min_sep_bp", "premix"):
        if col in meta:
            df[col] = df.community_id.map(meta[col])
    df.loc[df.well_type == "blank", "size"] = 0
    df["size"] = df["size"].astype(int)
    r = df.well.str[0]; c = df.well.str[1:].astype(int)
    df["row"], df["col"] = r, c
    return df.sort_values(["plate", "row", "col"]).reset_index(drop=True)
