"""How many distinct 16S sequences are actually in each well, and in what proportion?

Shared version of 20260721/04_qc/check_reads.ipynb, with the per-well representative
extraction from per_strain_consensus.ipynb folded into the same pass -- both need the
same clustering, and doing it twice was the slow part of the notebook.

**Clustering method changed, and it matters.** The notebook used single-linkage
connected components, on the argument that the same-strain and different-strain
distance modes are cleanly separated so nothing can chain across the gap. That holds
at ~100 reads per well and fails at ~500: error-rich and chimeric reads become
numerous enough to bridge the two modes, the two strains merge into one component,
and the well is reported as containing a single strain. Measured on 200 pair wells
of 20260907:

    reads used   single-linkage   greedy seeds
       120           27.0%           47.0%       <- wells resolving both strains
       500           10.0%           45.0%

Single-linkage gets *worse* with more data, which is the signature of chaining.
Greedy seed clustering decides membership against a cluster seed rather than against
a neighbour, so it cannot chain, and it is stable in depth. It is also cheaper --
O(n*k) rather than O(n^2) -- which is what makes a larger read cap affordable.

The pooled distance threshold is still learned from the data, from a sample of wells,
rather than assumed.
"""
import collections
import numpy as np
import pandas as pd
import edlib
from multiprocessing import Pool

from io_utils import load_layout, sample_fastq_path, load_reads

# ---- thresholds about ONT/16S statistics, not the experiment ----
LEN_LO, LEN_HI = 1300, 1600     # full-length 16S window, after barcode/primer trimming
READ_CAP = 400                  # reads clustered per well
THRESHOLD_CAP = 120             # reads per well in the all-pairwise threshold estimate
THRESHOLD_WELLS = 300           # wells sampled to learn the threshold
MIN_STRAIN_READS = 3            # a cluster is a strain if it has >= this many reads ...
MIN_STRAIN_FRAC = 0.10          # ... and >= this fraction of the well
MEDOID_CAP = 25
FALLBACK_THRESHOLD = 0.08


def _norm(a, b):
    return edlib.align(a, b, mode="NW", task="distance")["editDistance"] / max(len(a), len(b))


def _pairwise(seqs):
    """Condensed (i<j) normalized edit distances. Only used to learn the threshold."""
    n = len(seqs)
    lens = [len(s) for s in seqs]
    d = np.empty(n * (n - 1) // 2, dtype=np.float32)
    k = 0
    for i in range(n):
        si, li = seqs[i], lens[i]
        for j in range(i + 1, n):
            ed = edlib.align(si, seqs[j], mode="NW", task="distance")["editDistance"]
            d[k] = ed / max(li, lens[j]); k += 1
    return d


def greedy_clusters(seqs, thr):
    """Each read joins its NEAREST seed if within `thr`, else becomes a new seed.

    Membership is decided against a seed, never against an arbitrary neighbour, so a
    chain of intermediate reads cannot merge two genuinely distinct strains.
    Returns a list of member-index lists, largest first.
    """
    seeds, members = [], []
    for i, s in enumerate(seqs):
        best, bj = 1e9, -1
        for j, cs in enumerate(seeds):
            d = _norm(s, cs)
            if d < best:
                best, bj = d, j
        if best <= thr:
            members[bj].append(i)
        else:
            seeds.append(s); members.append([i])
    return sorted(members, key=len, reverse=True)


def _medoid(seqs):
    if len(seqs) == 1:
        return seqs[0]
    sub = seqs[:MEDOID_CAP]
    best, best_tot = sub[0], None
    for a in sub:
        tot = sum(edlib.align(a, b, mode="NW", task="distance")["editDistance"] for b in sub)
        if best_tot is None or tot < best_tot:
            best_tot, best = tot, a
    return best


def _load_well(args):
    sample_id, path, seed, cap = args
    seqs = [s for _, s in load_reads(path) if LEN_LO <= len(s) <= LEN_HI]
    if len(seqs) > cap:
        rng = np.random.default_rng(seed)
        seqs = [seqs[i] for i in rng.choice(len(seqs), cap, replace=False)]
    return sample_id, seqs


def _threshold_job(seqs):
    return _pairwise(seqs[:THRESHOLD_CAP]) if len(seqs) > 1 else np.empty(0, np.float32)


def _cluster_job(args):
    sample_id, seqs, thr = args
    n = len(seqs)
    comps = greedy_clusters(seqs, thr)
    cutoff = max(MIN_STRAIN_READS, MIN_STRAIN_FRAC * n)
    strain_comps = [c for c in comps if len(c) >= cutoff]
    reps = [_medoid([seqs[i] for i in c]) for c in strain_comps]
    sizes = [len(c) for c in comps]
    return sample_id, dict(
        n_reads=n, n_strains=len(strain_comps), n_raw_clusters=len(comps),
        cluster_sizes=sizes[:6],
        top_frac=sizes[0] / n if sizes else 0.0,
        minority_frac=sizes[1] / n if len(sizes) > 1 else 0.0,
    ), reps


def pooled_threshold(dist_vectors, lo=0.03, hi=0.15):
    """Place the same/different cut in the empty valley between the two modes."""
    pooled = np.concatenate([d for d in dist_vectors if d.size])
    hist, edges = np.histogram(pooled, bins=np.linspace(0, 0.4, 81))
    centers = 0.5 * (edges[:-1] + edges[1:])
    band = (centers >= lo) & (centers <= hi)
    empty = band & (hist <= max(1, 0.001 * hist.max()))
    thr = round(float(np.median(centers[empty])), 3) if empty.any() else FALLBACK_THRESHOLD
    return thr, pooled


def analyze(cfg, min_reads=None, processes=12, seed=0, read_cap=READ_CAP):
    """Returns (per_well_df, well_reps, threshold, pooled_distances)."""
    min_reads = cfg.min_reads if min_reads is None else min_reads
    lay = load_layout(cfg.layout_csv)
    if "plate_role" in lay.columns:
        lay = lay[lay.plate_role == "primary_sequenced"]
    jobs = []
    for i, r in enumerate(lay.itertuples()):
        p = sample_fastq_path(cfg.demux_dir, r.dest_plate, r.dest_well)
        if p is not None:
            jobs.append((r.sample_id, p, seed + i, read_cap))

    with Pool(processes) as pool:
        loaded = dict(pool.map(_load_well, jobs, chunksize=32))
    usable = {k: v for k, v in loaded.items() if len(v) > min_reads}
    print(f"{len(jobs)} wells with a fastq, {len(usable)} with more than {min_reads} "
          f"length-filtered reads (median {np.median([len(v) for v in usable.values()]):.0f} "
          f"reads used per well, cap {read_cap})")

    rng = np.random.default_rng(seed)
    keys = list(usable)
    pick = [keys[i] for i in rng.choice(len(keys), min(THRESHOLD_WELLS, len(keys)),
                                        replace=False)]
    with Pool(processes) as pool:
        dvecs = pool.map(_threshold_job, [usable[k] for k in pick], chunksize=8)
    thr, pooled = pooled_threshold(dvecs)
    print(f"data-driven same/different-strain threshold = {thr} (normalized edit "
          f"distance, from {pooled.size:,} read pairs over {len(pick)} wells)")

    with Pool(processes) as pool:
        out = pool.map(_cluster_job, [(k, v, thr) for k, v in usable.items()], chunksize=8)

    meta = lay.set_index("sample_id")
    rows, well_reps = [], {}
    for sid, stats, reps in out:
        well_reps[sid] = reps
        m = meta.loc[sid]
        rows.append(dict(sample_id=sid, well_type=m.well_type, strain1=m.strain1,
                         strain2=m.strain2,
                         separation_bp=m.get("separation_bp", np.nan), **stats))
    return pd.DataFrame(rows), well_reps, thr, pooled


def summarize(df):
    """The QC answer: do mono wells hold one sequence and pair wells two?"""
    print("\nestimated distinct 16S sequences per well:")
    tab = pd.crosstab(df.well_type, df.n_strains.clip(upper=4), normalize="index")
    tab.columns = [f"{c}{'+' if c == 4 else ''} seq" for c in tab.columns]
    counts = df.well_type.value_counts().rename("wells")
    print(pd.concat([counts, (tab * 100).round(1)], axis=1).to_string())

    mono = df[df.well_type == "mono"]
    two = df[df.well_type.isin(["pair", "techrep", "ratio", "density"])]
    if len(mono):
        print(f"\nmono wells ({len(mono)}): {(mono.n_strains == 1).mean():.1%} hold "
              f"exactly one sequence.")
        print("  A monoculture has one organism by construction, so anything else in it "
              "is\n  contamination -- this is the per-well readout the 760 mono wells "
              "were for:")
        print(f"    minority fraction: median {mono.minority_frac.median():.3f}, "
              f"75th {mono.minority_frac.quantile(.75):.3f}, "
              f"90th {mono.minority_frac.quantile(.9):.3f}")
        for t in (0.20, 0.10, 0.05, 0.02):
            print(f"    contaminated at >={t:>4.0%}: {(mono.minority_frac >= t).mean():6.1%}")
        print(f"    completely clean (no second sequence at all): "
              f"{(mono.minority_frac == 0).mean():.1%}")
    if len(two):
        print(f"\ntwo-strain wells: the binary count depends entirely on where you put "
              f"the\nminority cutoff, so here is the distribution instead "
              f"(n={len(two)}):")
        print(f"  minority cluster fraction: median {two.minority_frac.median():.3f}, "
              f"75th {two.minority_frac.quantile(.75):.3f}, "
              f"90th {two.minority_frac.quantile(.9):.3f}")
        for t in (0.20, 0.10, 0.05, 0.02, 0.01):
            print(f"    both sequences present at >={t:>4.0%}: "
                  f"{(two.minority_frac >= t).mean():6.1%}")
        print(f"    second sequence entirely absent: "
              f"{(two.minority_frac == 0).mean():.1%}")
        print("  A low minority fraction is a real result -- competitive exclusion -- "
              "not a\n  failure to detect, provided the two references are separable; "
              "check\n  separation_bp before reading any individual well that way.")
    return tab
