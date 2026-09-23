"""Pairwise 16S separability for the 384-well master plate.

Builds, for every pair of master-plate wells, how far apart their 16S sequences
are according to each independent reference, and caches the result. The
experiment design uses this to avoid spending wells on pairs whose members
cannot be told apart in the sequencing readout.

Sources come from ../../../merge_consensus_sequences/external_reference_crosscheck,
already projected onto one coordinate system:

  ont_v3, pool1   ONT runs on the WORKING plates (carry the contamination)
  sanger, ngs     original-isolate references (predate the working plates)

Two decisions worth stating:

* **Conservative combination.** A pair's separability is the MINIMUM across the
  sources that cover both wells. If any source says the two are close, the pair
  is treated as risky. The cost of being wrong is asymmetric -- an unresolvable
  pair wastes a well and yields nothing, while skipping a resolvable pair costs
  only one of many thousands of alternatives -- and there are far more eligible
  pairs than the design can use.

* **Absolute counts, not percentages.** What the ONT readout resolves is the
  number of discriminating positions, so a raw count over the compared columns
  is the right unit. A minimum overlap is enforced so a short reference cannot
  manufacture a small distance.

Ambiguity codes are handled exactly: each base is a 4-bit mask and two positions
agree when their masks intersect, so pool1's 6% IUPAC content costs nothing.
"""
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
XREF = os.path.abspath(os.path.join(
    HERE, "..", "..", "..", "merge_consensus_sequences", "external_reference_crosscheck"))
CACHE = os.path.join(HERE, "separability_cache.npz")

MASK = {"A": 1, "C": 2, "G": 4, "T": 8, "-": 0, "R": 5, "Y": 10, "S": 6, "W": 9,
        "K": 12, "M": 3, "B": 14, "D": 13, "H": 11, "V": 7, "N": 15}
SOURCES = ["ont_v3", "pool1", "sanger", "ngs"]
WORKING, ORIGINAL = ["ont_v3", "pool1"], ["sanger", "ngs"]
MIN_OVERLAP = 900


def _load():
    M = np.load(os.path.join(XREF, "work", "02_matrix.npy"))
    ids = open(os.path.join(XREF, "work", "02_ids.txt")).read().split("\n")
    enc = np.zeros(M.shape, dtype=np.uint8)
    for ch, m in MASK.items():
        enc[M == ch] = m
    by = defaultdict(list)
    for i, s in enumerate(ids):
        w, src, c = s.split("|")
        by[(w, src)].append(i)
    return enc, by


def build(wells, force=False):
    """Return (dist, cov) dicts keyed by source: (n,n) int arrays.

    dist[src][i, j] is the mismatch count between wells i and j under that
    source, -1 where the source cannot compare them. For ngs, which carries
    several rRNA copies per strain, the minimum over copy pairs is taken: two
    strains sharing any one operon copy are confusable in a mixed readout, so
    the minimum is the honest answer.
    """
    wells = list(wells)
    if os.path.exists(CACHE) and not force:
        z = np.load(CACHE, allow_pickle=True)
        if list(z["wells"]) == wells:
            return {s: z[f"d_{s}"] for s in SOURCES}, list(z["wells"])

    enc, by = _load()
    n = len(wells)
    out = {}
    for src in SOURCES:
        D = np.full((n, n), -1, dtype=np.int32)
        rows = [by.get((w, src), []) for w in wells]
        for i in range(n):
            if not rows[i]:
                continue
            for j in range(i + 1, n):
                if not rows[j]:
                    continue
                best = None
                for a in rows[i]:
                    ea = enc[a]
                    for b in rows[j]:
                        eb = enc[b]
                        both = (ea != 0) & (eb != 0)
                        m = int(both.sum())
                        if m < MIN_OVERLAP:
                            continue
                        d = int(((ea & eb) == 0)[both].sum())
                        if best is None or d < best:
                            best = d
                if best is not None:
                    D[i, j] = D[j, i] = best
        out[src] = D
        print(f"  {src:8s} comparable pairs: {int((D >= 0).sum() // 2)}")
    np.savez_compressed(CACHE, wells=np.array(wells),
                        **{f"d_{s}": out[s] for s in SOURCES})
    return out, wells


def combine(dist, wells, mode="conservative"):
    """Collapse the per-source matrices into one separability matrix.

    -1 marks a pair no source can compare. 'conservative' takes the minimum over
    covering sources; 'original' uses only the original-isolate lineage, which
    is the better predictor for a fresh copy taken from the freezer master plate.
    """
    srcs = {"conservative": SOURCES, "original": ORIGINAL, "working": WORKING}[mode]
    n = len(wells)
    stack = np.stack([dist[s] for s in srcs])
    has = stack >= 0
    big = np.where(has, stack, np.iinfo(np.int32).max)
    out = big.min(axis=0).astype(np.int32)
    out[~has.any(axis=0)] = -1
    np.fill_diagonal(out, 0)
    return out, has.sum(axis=0)
