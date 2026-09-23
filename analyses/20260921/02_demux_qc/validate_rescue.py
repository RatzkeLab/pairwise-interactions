"""Are the 1.4M newly-assigned reads landing in the right wells?

The earlier attempt to validate by re-parsing barcodes failed: minibar trims the
barcode and primer off every read it writes, so the output carries no barcode to
check. The biological check works instead, and is stronger.

For each well present in both demultiplexings, take the dominant sequence from the
ORIGINAL reads (assigned under the strict 59-mer match, and therefore trustworthy)
and the dominant sequence from the CORRECTED reads (4.8x more of them). If the
extra reads belong to that well, the two agree to within ONT consensus error. If
the relaxed matching were pulling in reads from elsewhere, the new dominant
sequence would drift away from the old one, and it would drift most in the wells
that gained the most reads -- so the gain is checked as a covariate, not just the
overall agreement.
"""
import sys, random
from pathlib import Path
from multiprocessing import Pool

import numpy as np, pandas as pd, edlib
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE)); sys.path.insert(0, str(BASE.parents[1] / "scripts"))
import config
from io_utils import load_reads
from strain_identity import flexible_distance

OUT = BASE / "02_demux_qc" / "outputs"
LEN_LO, LEN_HI = 1300, 1600
CAP = 40                 # reads per well used for the medoid
AGREE_BP = 15            # two medoids this close are the same organism

# The two demultiplexings trim differently: minibar removes the "index" it was given,
# which was the 59 bp construct before the fix and the 24 bp barcode after it, so the
# corrected reads keep ~20 bp more at each end. Comparing end-to-end would score every
# well as a 44 bp mismatch for reasons that have nothing to do with which organism is
# in it, so the medoids are compared over their overlap.


def medoid(path, cap=CAP, seed=0):
    seqs = [s for _, s in load_reads(path) if LEN_LO <= len(s) <= LEN_HI]
    if len(seqs) < 5:
        return None, len(seqs)
    n_all = len(seqs)
    if len(seqs) > cap:
        seqs = random.Random(seed).sample(seqs, cap)
    best, best_tot = None, None
    for a in seqs:
        tot = sum(edlib.align(a, b, mode="NW", task="distance")["editDistance"] for b in seqs)
        if best_tot is None or tot < best_tot:
            best_tot, best = tot, a
    return best, n_all


def one(args):
    sample, old_p, new_p = args
    m_old, n_old = medoid(old_p)
    m_new, n_new = medoid(new_p)
    if m_old is None or m_new is None:
        return None
    dist, overlap = flexible_distance(m_old, m_new)
    return dict(sample=sample, n_old=n_old, n_new=n_new,
                medoid_dist_bp=round(dist * overlap), overlap_bp=overlap,
                len_old=len(m_old), len_new=len(m_new),
                gain=n_new / max(n_old, 1))


def main(max_wells=400):
    old_dir, new_dir = config.DEMUX_DIR_ORIGINAL, config.DEMUX_DIR
    jobs = []
    for p in sorted(new_dir.glob("*.fastq")):
        o = old_dir / p.name
        if o.exists() and o.stat().st_size > 0 and p.stat().st_size > 0:
            jobs.append((p.stem, o, p))
    print(f"{len(jobs)} wells present and non-empty in both demultiplexings")
    if len(jobs) > max_wells:
        jobs = random.Random(0).sample(jobs, max_wells)
    with Pool(16) as pool:
        res = [r for r in pool.map(one, jobs) if r]
    df = pd.DataFrame(res)
    if df.empty:
        print("no well had enough reads in both runs yet"); return df

    agree = df.medoid_dist_bp <= AGREE_BP
    print(f"\n{len(df)} wells compared (>=5 usable reads in both)")
    print(f"  median reads: original {df.n_old.median():.0f} -> corrected "
          f"{df.n_new.median():.0f}  (median gain {df.gain.median():.1f}x)")
    print(f"  median read length: original {df.len_old.median():.0f} bp -> corrected "
          f"{df.len_new.median():.0f} bp (less aggressive trimming, so more 16S kept)")
    print(f"\n  dominant sequence agrees within {AGREE_BP} bp: "
          f"{agree.sum()}/{len(df)} ({agree.mean():.1%})")
    print(f"  medoid distance: median {df.medoid_dist_bp.median():.0f} bp, "
          f"90th pct {df.medoid_dist_bp.quantile(.9):.0f} bp")

    # if relaxed matching were leaking reads, disagreement would grow with the gain
    hi = df[df.gain > df.gain.median()]
    lo = df[df.gain <= df.gain.median()]
    print(f"\n  wells gaining MORE reads than median: {(hi.medoid_dist_bp<=AGREE_BP).mean():.1%} agree")
    print(f"  wells gaining less                  : {(lo.medoid_dist_bp<=AGREE_BP).mean():.1%} agree")
    from scipy import stats
    rho, p = stats.spearmanr(df.gain, df.medoid_dist_bp)
    print(f"  rho(read gain, medoid drift) = {rho:+.3f}  p={p:.3g}")
    if p > 0.05 or rho <= 0:
        print("  -> disagreement does not grow with the number of extra reads. The "
              "newly\n     assigned reads match what was already in the well.")
    else:
        print("  -> wells that gained most also drifted most; the relaxed matching may "
              "be leaking.")
    # A monoculture has one organism, so its dominant sequence cannot legitimately
    # change; that is the clean measure of assignment correctness. A two-strain well's
    # dominant sequence CAN flip when the read count quintuples, because the medoid
    # simply follows whichever cluster is larger -- so those are reported separately
    # rather than counted as errors.
    lay = pd.read_csv(config.LAYOUT_CSV)
    lay["sample"] = lay.apply(lambda r: f"Plate{int(r.dest_plate):02d}_{r.dest_well}", axis=1)
    df = df.merge(lay[["sample", "well_type"]], on="sample", how="left")
    df["agree"] = df.medoid_dist_bp <= AGREE_BP
    by = df.groupby(df.well_type.eq("mono").map({True: "mono", False: "two-strain"})).agree
    print("\n  split by well type:")
    for k, v in by.agg(["size", "mean"]).iterrows():
        print(f"    {k:11s} n={int(v['size']):4d}  agree {v['mean']:.1%}")
    if "mono" in by.groups:
        print(f"  -> the monoculture figure is the one that measures assignment: one "
              f"organism,\n     so its dominant sequence has no legitimate reason to "
              f"move.")

    df.to_csv(OUT / "d14_rescue_validation.csv", index=False)
    print(f"\nwrote -> {OUT}/d14_rescue_validation.csv")
    return df


if __name__ == "__main__":
    main()
