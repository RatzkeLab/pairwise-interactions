"""Does 20260721's mix-up correspond to a plate-handling mistake we can name?

The general permutation search already failed (no dominant row/column shift), but that search
assumed the error lives in 384-well coordinates. The source collection was cultured in 96-well
format and recombined into 384 for the Echo, so the natural unit of error is not a row or a
column -- it is a **96-well quadrant**, and quadrants interleave in a checkerboard, not in
spatial blocks. A quadrant swap therefore looks like noise to a shift-based search while being
a single, simple handling mistake.

Two conventions exist for the 96 -> 384 recombination and both are tested:

  interleaved  384 (R, C) belongs to quadrant (R mod 2, C mod 2); its 96-well address is
               (R div 2, C div 2). This is what a 4-tip head or a Mosquito produces, and is the
               "four-corner checkerboard" pattern.
  block        each 96 plate occupies one spatial 8x12 corner of the 384.

On top of those, plate-orientation errors: 180 deg rotation (plate loaded backwards), and row-
or column-only flips.

Everything here is restricted to **resolvable** strains -- ones whose 16S has no near-twin
anywhere in the collection. That is the only way the question can be asked honestly: most of
this collection is too similar at 16S to attribute at all (see qc_recovery), and a positional
pattern computed from ambiguous calls would be a pattern in the ambiguity.
"""

import itertools
import re

import numpy as np
import pandas as pd

import qc_config as C

N_ROW, N_COL = 16, 24          # 384-well
RESOLVABLE_MAX_NEIGHBOUR = 0.99  # a strain with any neighbour at or above this is NOT resolvable


def well_to_rc(w):
    m = re.match(r"^([A-P])(\d+)$", str(w))
    if not m:
        return None
    r, c = ord(m.group(1)) - 65, int(m.group(2)) - 1
    return (r, c) if 0 <= r < N_ROW and 0 <= c < N_COL else None


def rc_to_well(rc):
    return f"{chr(65 + rc[0])}{rc[1] + 1}"


# ---- candidate transforms -------------------------------------------------

def _quadrant_perm(perm, mode):
    """Permute which 96-well source plate landed in which quadrant of the 384."""
    quads = [(0, 0), (0, 1), (1, 0), (1, 1)]
    mapping = dict(zip(quads, [quads[i] for i in perm]))

    def f(rc, mapping=mapping, mode=mode):
        r, c = rc
        if mode == "interleaved":
            q, inner = (r % 2, c % 2), (r // 2, c // 2)
            nq = mapping[q]
            return (inner[0] * 2 + nq[0], inner[1] * 2 + nq[1])
        q, inner = (r // 8, c // 12), (r % 8, c % 12)
        nq = mapping[q]
        return (nq[0] * 8 + inner[0], nq[1] * 12 + inner[1])
    return f


def _rot180(rc):
    return (N_ROW - 1 - rc[0], N_COL - 1 - rc[1])


def _flip_rows(rc):
    return (N_ROW - 1 - rc[0], rc[1])


def _flip_cols(rc):
    return (rc[0], N_COL - 1 - rc[1])


def _rot180_96(rc):
    """The 96-well source plate loaded backwards: rotate within quadrant, keep the quadrant."""
    r, c = rc
    q, inner = (r % 2, c % 2), (r // 2, c // 2)
    ni = (7 - inner[0], 11 - inner[1])
    return (ni[0] * 2 + q[0], ni[1] * 2 + q[1])


def _list_shift(k, order):
    """An off-by-N in a pick list: everything slides by k positions in some LINEAR ordering.

    Geometric transforms cannot express this, but it is one of the most common ways a plate map
    goes wrong -- a row inserted or deleted in a spreadsheet shifts every subsequent entry. Two
    orderings are plausible: the plate read row-major or column-major.
    """
    cells = ([(r, c) for r in range(N_ROW) for c in range(N_COL)] if order == "rowmajor"
             else [(r, c) for c in range(N_COL) for r in range(N_ROW)])
    idx = {rc: i for i, rc in enumerate(cells)}
    n = len(cells)

    def f(rc, idx=idx, cells=cells, k=k, n=n):
        return cells[(idx[rc] + k) % n]
    return f


def candidate_transforms():
    """-> {name: rc -> rc}. Orientation errors, quadrant mix-ups, list shifts, compositions."""
    base = {"identity": lambda rc: rc, "rot180_384": _rot180,
            "flip_rows": _flip_rows, "flip_cols": _flip_cols, "rot180_within_96": _rot180_96}
    out = dict(base)
    for order in ("rowmajor", "colmajor"):
        for k in list(range(-25, 0)) + list(range(1, 26)):
            out[f"shift_{order}_{k:+d}"] = _list_shift(k, order)
    for mode in ("interleaved", "block"):
        for perm in itertools.permutations(range(4)):
            if perm == (0, 1, 2, 3):
                continue
            name = f"quad_{mode}_{''.join(map(str, perm))}"
            f = _quadrant_perm(perm, mode)
            out[name] = f
            for bname, bf in base.items():
                if bname != "identity":
                    out[f"{name}+{bname}"] = (lambda rc, f=f, bf=bf: bf(f(rc)))
    return out


# ---- resolvable strains ---------------------------------------------------

def resolvable_strains(max_neighbour=RESOLVABLE_MAX_NEIGHBOUR):
    """Strains with no 16S near-twin anywhere in the collection, and their nearest neighbour."""
    M = pd.read_csv(C.OUT / "s07_collection_16S_identity_matrix.csv", index_col=0)
    A = M.values.copy()
    np.fill_diagonal(A, -1)
    nn = A.max(axis=1)
    df = pd.DataFrame({"strain_label": M.index, "nearest_neighbour_identity": nn,
                       "nearest_neighbour": M.columns[A.argmax(axis=1)]})
    df["resolvable"] = df["nearest_neighbour_identity"] < max_neighbour
    df.to_csv(C.OUT / "s08_strain_resolvability.csv", index=False)
    return df


def positional_calls(exp, resolvability, min_identity=0.99, min_margin=0.005):
    """label position -> observed position, for calls that are actually trustworthy.

    Three filters, all necessary: the match must be at organism level, it must beat its runner-up
    by more than ONT consensus error, and the strain it matches must have no near-twin that could
    have produced the same match.
    """
    att = pd.read_csv(exp.out / f"s02_attribution_pair_vs_genome_16S_{exp.name}.csv")
    res = set(resolvability.loc[resolvability["resolvable"], "strain_label"])
    d = att[(att["best_identity"] >= min_identity)
            & (att["margin_over_runner_up"] >= min_margin)
            & (att["best_match"].isin(res))].copy()
    d["from_rc"] = d["strain_label"].map(well_to_rc)
    d["to_rc"] = d["best_match"].map(well_to_rc)
    d = d.dropna(subset=["from_rc", "to_rc"])
    return d.reset_index(drop=True)


def transform_search(calls):
    """Score every candidate: what fraction of trustworthy calls does it explain?"""
    if not len(calls):
        return pd.DataFrame()
    rows = []
    for name, f in candidate_transforms().items():
        hit = sum(f(a) == b for a, b in zip(calls["from_rc"], calls["to_rc"]))
        rows.append({"transform": name, "n_calls": len(calls), "n_explained": hit,
                     "pct_explained": round(100 * hit / len(calls), 1)})
    out = pd.DataFrame(rows).sort_values("n_explained", ascending=False)
    # chance level: a random relabelling of 294 collection positions
    out["pct_expected_by_chance"] = round(100 / 294, 2)
    return out


# ---- group-level search ---------------------------------------------------
# Only 35 of 294 collection strains have no 16S near-twin, which leaves ~5 usable positional
# calls for 20260721 -- far too few to test ~200 candidate transforms without fitting noise.
# So the question is asked one level coarser: not "which strain is in this well" but "which
# 16S GROUP", where a group is a set of strains 16S cannot tell apart. A quadrant swap still
# has to show up at group level, because a group occupies specific plate positions.

def sixteen_s_groups(threshold=RESOLVABLE_MAX_NEIGHBOUR):
    """Connected components of the collection at >=threshold identity: what 16S can distinguish."""
    M = pd.read_csv(C.OUT / "s07_collection_16S_identity_matrix.csv", index_col=0)
    strains = list(M.index)
    A = M.values >= threshold
    seen, group_of = {}, {}
    for i, s in enumerate(strains):
        if s in group_of:
            continue
        stack, comp = [i], []
        seen[i] = True
        while stack:
            k = stack.pop()
            comp.append(k)
            for j in np.where(A[k])[0]:
                if j not in seen:
                    seen[j] = True
                    stack.append(j)
        gid = len(set(group_of.values()))
        for k in comp:
            group_of[strains[k]] = gid
    sizes = pd.Series(group_of).value_counts()
    return group_of, sizes


def group_level_transform_search(exp, group_of, min_identity=0.99, n_perm=500, seed=0):
    """Which plate-handling mistake, if any, explains where strains actually turned up?

    A transform T explains label L if the collection strain sitting at T(position(L)) is in the
    same 16S group as what L's well actually contained.

    Roughly 200 candidate transforms are tested, so a per-transform z-score is not
    interpretable on its own -- the best of 200 draws from a null is large by construction. The
    permutations are therefore SHARED across transforms, and the family-wise null is the
    distribution of the *maximum* z over all transforms. `p_familywise` is the fraction of
    permutations whose best transform beat the observed best, which is the number to read.
    """
    att = pd.read_csv(exp.out / f"s02_attribution_pair_vs_genome_16S_{exp.name}.csv")
    d = att[att["best_identity"] >= min_identity].copy()
    d["from_rc"] = d["strain_label"].map(well_to_rc)
    d["obs_group"] = d["best_match"].map(group_of)
    d = d.dropna(subset=["from_rc", "obs_group"])
    at_rc = {well_to_rc(s): g for s, g in group_of.items() if well_to_rc(s)}
    obs = d["obs_group"].to_numpy(float)

    names, T = [], []
    for name, f in candidate_transforms().items():
        names.append(name)
        T.append([at_rc.get(f(rc), np.nan) for rc in d["from_rc"]])
    T = np.array(T, dtype=float)
    covered = ~np.isnan(T)
    hits = np.nansum(T == obs[None, :], axis=1)

    rng = np.random.default_rng(seed)
    null = np.empty((n_perm, len(names)))
    for i in range(n_perm):
        null[i] = np.nansum(T == rng.permutation(obs)[None, :], axis=1)
    mu, sd = null.mean(0), null.std(0)
    z = np.where(sd > 0, (hits - mu) / np.where(sd > 0, sd, 1), np.nan)
    null_z = (null - mu) / np.where(sd > 0, sd, 1)
    fam = float((np.nanmax(null_z, axis=1) >= np.nanmax(z)).mean())

    out = pd.DataFrame({
        "transform": names, "n_calls": covered.sum(1), "n_explained": hits,
        "pct_explained": (100 * hits / np.maximum(covered.sum(1), 1)).round(1),
        "null_mean": mu.round(1), "z_vs_null": z.round(2),
    }).sort_values("z_vs_null", ascending=False)
    out["n_transforms_tested"] = len(names)
    out["familywise_p_for_best"] = round(fam, 4)
    out["familywise_z_threshold_p05"] = round(float(np.nanpercentile(np.nanmax(null_z, axis=1), 95)), 2)
    return out
