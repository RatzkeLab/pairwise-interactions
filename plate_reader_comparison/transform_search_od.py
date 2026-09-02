"""Is 20260721's source-well assignment a plate-handling transform of the truth? Asked with OD.

`strain_identity_qc/qc_layout.py` already searched 335 plate-handling hypotheses using 16S and
found nothing -- but it said so itself: only 5 usable single-strain calls survive 16S
degeneracy, "far too few for 335 hypotheses". That negative was underpowered, not conclusive.

OD is a far better instrument for this particular question: 171 strain pairs are shared between
the experiments, and the same pair reproduces within an experiment at Spearman 0.62-0.82. So a
transform search that was hopeless on 16S has real power here.

The hypothesis: 20260721's wells hold the right collection, mis-addressed -- the 4x96
precultures consolidated into the Echo source plate in the wrong quadrant order, the plate
loaded rotated, or a pick-list offset. Each is a permutation T of source-well labels. Under the
true T, relabelling 20260721's strains should make its shared pairs' ODs line up with
20260630's.

Scoring: for each T, relabel 20260721's pairs, intersect with 20260630's, and correlate mean
per-pair OD (plate z-scored).

**RESULT: no transform is identified, and the search itself is biased -- read this before
trusting any ranking it produces.** Different transforms leave different numbers of shared
pairs (n from 40 to 4000+), and a transform that shrinks the overlap is scored on an easier,
noisier subset. The split-half positive control (`t03_positive_control_split_half.csv`) makes
this concrete: relabelling one half of 20260630 and searching against the other, where the
answer is known to be `identity`, ranks **identity only 2nd of 142** -- beaten by
`quad_interleaved_0132` on a third of the data. Among transforms reaching n<=60 in that control,
where every one of them is wrong by construction, the best spurious rho is **0.464**.

The best 20260721 candidate, `shift_rowmajor_-7`, scores rho=0.433 at n=48 -- BELOW the 0.464
that wrong transforms reach by artifact at that sample size, though above the naive family-wise
threshold of 0.299 computed here (that threshold permutes labels within each transform but does
not capture the subset-selection bias, so it is too lenient). It is not a finding.

Large-n transforms are where this search does have power -- `identity` reaches rho=0.84 in the
control at n=2157 -- and no 20260721 transform with substantial overlap gets near that
(best: shift_colmajor_-22, rho=0.288 at n=232; untransformed, rho=0.084 at n=171). So the
conclusion is a genuine negative: **no plate-handling transform relabels 20260721 into
agreement with 20260630.** This agrees with qc_layout.py's 16S-based negative, but unlike that
one it is not limited by 16S degeneracy.
"""

import sys
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from experiment_od_overlap import CONFIG, load_experiment, EXPS

OUT = HERE / "outputs"
MIN_SHARED = 40


def wells():
    return [f"{chr(65+r)}{c+1}" for r in range(16) for c in range(24)]


def rc(w):
    return ord(w[0]) - 65, int(w[1:]) - 1


def build_transforms():
    """Quadrant permutations under both 96->384 conventions, x plate orientations, x offsets."""
    W = wells()
    idx = {w: rc(w) for w in W}
    grid = np.array(W).reshape(16, 24)
    T = {}

    def add(name, mapping):
        if len(set(mapping.values())) == len(mapping):
            T[name] = mapping

    T["identity"] = {w: w for w in W}
    # orientations
    for name, g in [("rot180", grid[::-1, ::-1]), ("flip_rows", grid[::-1, :]),
                    ("flip_cols", grid[:, ::-1])]:
        add(name, {W[i]: g.ravel()[i] for i in range(len(W))})
    # pick-list offsets, row- and column-major
    flat_r = grid.ravel()
    flat_c = grid.T.ravel()
    for n in list(range(-25, 26)):
        if n == 0:
            continue
        add(f"shift_rowmajor_{n:+d}", {flat_r[i]: flat_r[(i + n) % 384] for i in range(384)})
        add(f"shift_colmajor_{n:+d}", {flat_c[i]: flat_c[(i + n) % 384] for i in range(384)})
    # quadrant permutations, interleaved (checkerboard) and block conventions
    inter = [(0, 0), (0, 1), (1, 0), (1, 1)]
    quad_wells_inter = [[grid[r::2, c::2].ravel() for (r, c) in inter]][0]
    blocks = [grid[:8, :12].ravel(), grid[:8, 12:].ravel(),
              grid[8:, :12].ravel(), grid[8:, 12:].ravel()]
    for perm in permutations(range(4)):
        if perm == (0, 1, 2, 3):
            continue
        m = {}
        for src, dst in enumerate(perm):
            for a, b in zip(quad_wells_inter[src], quad_wells_inter[dst]):
                m[a] = b
        add(f"quad_interleaved_{''.join(map(str,perm))}", m)
        m = {}
        for src, dst in enumerate(perm):
            for a, b in zip(blocks[src], blocks[dst]):
                m[a] = b
        add(f"quad_block_{''.join(map(str,perm))}", m)
    return T


def pair_od(lay):
    d = lay.dropna(subset=["od_z"])
    return d.groupby("pair_key")["od_z"].mean()


def score(sa, lay_b, tmap):
    b = lay_b.dropna(subset=["od_z"]).copy()
    s1 = b.strain1.map(tmap)
    s2 = b.strain2.map(tmap)
    keep = s1.notna() & s2.notna()
    b = b[keep]
    b = b.assign(pk=[frozenset((x, y)) for x, y in zip(s1[keep], s2[keep])])
    sb = b.groupby("pk")["od_z"].mean()
    shared = sorted(set(sa.index) & set(sb.index), key=lambda k: sorted(k))
    if len(shared) < MIN_SHARED:
        return np.nan, len(shared)
    x, y = sa.loc[shared].values, sb.loc[shared].values
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < MIN_SHARED:
        return np.nan, int(ok.sum())
    return float(spearmanr(x[ok], y[ok])[0]), int(ok.sum())


def cache_transform(lay_b, tmap):
    """Per transform, the relabelled pair->OD series. Computed once so the permutation null
    only has to re-index, not re-group -- the difference between minutes and hours."""
    b = lay_b.dropna(subset=["od_z"]).copy()
    s1, s2 = b.strain1.map(tmap), b.strain2.map(tmap)
    keep = s1.notna() & s2.notna()
    b = b[keep].assign(pk=[frozenset((x, y)) for x, y in zip(s1[keep], s2[keep])])
    return b.groupby("pk")["od_z"].mean()


def main():
    L = {k: load_experiment(k) for k in CONFIG}
    sa = pair_od(L["20260630"])
    T = build_transforms()
    print(f"\ntesting {len(T)} transforms\n")

    rows = []
    for name, tmap in T.items():
        rho, n = score(sa, L["20260721"], tmap)
        rows.append({"transform": name, "n_shared_pairs": n, "spearman": rho})
    df = pd.DataFrame(rows).dropna(subset=["spearman"]).sort_values("spearman", ascending=False)

    # positive control: the same search against 20260630 itself must recover `identity`
    ctrl = []
    for name, tmap in T.items():
        rho, n = score(sa, L["20260630"], tmap)
        ctrl.append({"transform": name, "n_shared_pairs": n, "spearman": rho})
    cdf = pd.DataFrame(ctrl).dropna(subset=["spearman"]).sort_values("spearman", ascending=False)

    # family-wise threshold: best-of-N under label permutation
    rng = np.random.default_rng(0)
    cached = []
    for name, tmap in T.items():
        sb = cache_transform(L["20260721"], tmap)
        shared = sorted(set(sa.index) & set(sb.index), key=lambda k: sorted(k))
        if len(shared) >= MIN_SHARED:
            cached.append((np.array([sa[k] for k in shared]),
                           np.array([sb[k] for k in shared])))
    best_null = []
    for _ in range(300):
        best = -1.0
        for x, y in cached:
            yp = rng.permutation(y)
            r = spearmanr(x, yp)[0]
            if r > best:
                best = r
        best_null.append(best)
    thr = float(np.quantile(best_null, 0.95))
    print(f"  (family-wise null built over {len(cached)} transforms x 300 permutations)")

    df.to_csv(OUT / "t01_transform_search_20260721.csv", index=False)
    cdf.to_csv(OUT / "t02_transform_search_control_20260630.csv", index=False)
    print("POSITIVE CONTROL (20260630 vs itself) -- top 3:")
    print(cdf.head(3).to_string(index=False))
    print("\n20260721 -- top 8 transforms:")
    print(df.head(8).to_string(index=False))
    print(f"\nfamily-wise 95% threshold from best-of-search under permuted labels: rho = {thr:.3f}")
    best = df.iloc[0]
    print(f"best real transform: {best.transform} rho={best.spearman:.3f} "
          f"-> {'EXCEEDS' if best.spearman > thr else 'does NOT exceed'} the threshold")
    return df, cdf, thr


if __name__ == "__main__":
    main()
