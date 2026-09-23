"""Do the 19 'failing' barcodes reappear among the reads minibar discarded?

If a barcode's amplicons are sitting in unk.fastq, the PCR worked and the
barcode's dropout is a demultiplexing artefact, not a primer failure. If they
are genuinely absent from the unk pile too, the oligo really did fail.
"""
import random, collections, sys
from multiprocessing import Pool
import numpy as np, pandas as pd, edlib
from pathlib import Path

D = Path("/home/rl/scripts/karl/data_links/interim/demultiplexing/20260907_demux/demuxed")
SETUP = config.SETUP
OUT = Path("/home/rl/scripts/karl/pairwise_interaction_experiments/20260907/02_demux_qc/outputs")
ADAPTER, FWDP, REVP = "ATCGCCTACCGTGAC", "AGRGTTYGATYMTGGCTCAG", "CGGYTACCTTGTTACGACTT"
EQ = [(k, v) for k, vs in {"R":"AG","Y":"CT","M":"AC"}.items() for v in vs]
BC = 24
N_PER_FILE = int(sys.argv[1]) if len(sys.argv) > 1 else 4000

def rc(s):
    c = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
    return "".join(c.get(x, "N") for x in reversed(s))

mb = pd.read_csv(SETUP / "minibar_primers_20260907.tsv", sep="\t")
mb["fwd_bc"] = mb.FwIndex.str[len(ADAPTER):-len(FWDP)]
mb["rev_bc"] = mb.RvIndex.str[len(ADAPTER):-len(REVP)]
lay = pd.read_csv(SETUP / "strain_layout_20260907.csv")
lay["sample"] = lay.apply(lambda r: f"Plate{int(r.dest_plate):02d}_{r.dest_well}", axis=1)
mb = mb.merge(lay[["sample", "plate_role", "fwd_idx", "rev_idx"]],
              left_on="SampleID", right_on="sample")
FWD, REV = sorted(set(mb.fwd_bc)), sorted(set(mb.rev_bc))
FIDX = dict(zip(mb.fwd_bc, mb.fwd_idx)); RIDX = dict(zip(mb.rev_bc, mb.rev_idx))
PAIR = {(f, r): (s, pr) for f, r, s, pr in
        zip(mb.fwd_bc, mb.rev_bc, mb.SampleID, mb.plate_role)}


def loc(pat, seq):
    r = edlib.align(pat, seq, mode="HW", task="locations", additionalEqualities=EQ)
    return r["editDistance"], (r["locations"][0] if r["locations"] else None)

def nearest(obs, pool):
    bd, bb = 99, None
    for b in pool:
        d = edlib.align(obs, b, mode="NW", task="distance", k=8)["editDistance"]
        if 0 <= d < bd: bd, bb = d, b
    return bd, bb

def call(s):
    for seq in (s, rc(s)):
        df, lf = loc(FWDP, seq[:300])
        if df > 5: continue
        if lf[0] - BC < 0: return None
        obs_f = seq[lf[0] - BC:lf[0]]
        tail = seq[-300:]
        dr, lr = loc(rc(REVP), tail)
        if dr > 5 or lr[1] + 1 + BC > len(tail): return None
        obs_r = rc(tail[lr[1] + 1:lr[1] + 1 + BC])
        d1, b1 = nearest(obs_f, FWD)
        d2, b2 = nearest(obs_r, REV)
        if b1 is None or b2 is None: return None
        return b1, b2
    return None


def process(f):
    rng = random.Random(hash(f.name) & 0xffff)
    with open(f) as fh:
        buf = [l.strip() for i, l in enumerate(fh) if i % 4 == 1]
    if len(buf) > N_PER_FILE: buf = rng.sample(buf, N_PER_FILE)
    fw, rv, sam = collections.Counter(), collections.Counter(), collections.Counter()
    n_ok = 0
    for s in buf:
        res = call(s)
        if res is None: continue
        b1, b2 = res
        fw[FIDX[b1]] += 1; rv[RIDX[b2]] += 1
        hit = PAIR.get((b1, b2))
        if hit: sam[hit[0]] += 1
        n_ok += 1
    return len(buf), n_ok, fw, rv, sam


if __name__ == "__main__":
    files = sorted(D.glob("*/unk.fastq"))
    print(f"{len(files)} unk files, sampling up to {N_PER_FILE} reads each")
    with Pool(12) as p:
        res = p.map(process, files)
    n_tot = sum(r[0] for r in res); n_ok = sum(r[1] for r in res)
    fw, rv, sam = collections.Counter(), collections.Counter(), collections.Counter()
    for _, _, a, b, c in res: fw += a; rv += b; sam += c
    print(f"parsed {n_ok:,}/{n_tot:,} sampled discarded reads ({n_ok/n_tot:.1%})")

    seq = pd.read_csv(OUT / "d00_sequenced_wells_with_counts.csv")
    bad_f = pd.read_csv(OUT / "d03_fwd_barcode_forensics.csv")
    bad_r = pd.read_csv(OUT / "d03_rev_barcode_forensics.csv")
    p0 = seq.dropout.mean()
    for side, tab, cnt, idxcol in [("FWD", bad_f, fw, "fwd_idx"),
                                   ("REV", bad_r, rv, "rev_idx")]:
        tab = tab.copy()
        tab["rescued_reads"] = tab[idxcol].map(cnt).fillna(0).astype(int)
        flagged = tab[(tab.q < 0.05) & (tab.rate > p0)]
        clean = tab[~tab.index.isin(flagged.index)]
        print(f"\n=== {side} ===")
        print(f"  flagged barcodes ({len(flagged)}): "
              f"{flagged.reads.sum():>7,} reads assigned by minibar | "
              f"{flagged.rescued_reads.sum():>7,} recovered from the discard pile")
        print(f"  clean barcodes   ({len(clean)}): "
              f"{clean.reads.sum():>7,} assigned | {clean.rescued_reads.sum():>7,} recovered")
        pw_f = flagged.rescued_reads.sum() / max(flagged.n.sum(), 1)
        pw_c = clean.rescued_reads.sum() / max(clean.n.sum(), 1)
        print(f"  recovered reads PER WELL: flagged {pw_f:.1f}  vs  clean {pw_c:.1f}"
              f"   ({pw_f/pw_c:.2f}x)")
        print(f"\n  per flagged barcode:")
        print(flagged[[idxcol, f"{side.lower()}_source_well", "n", "rate",
                       "reads", "rescued_reads"]]
              .assign(rescued_per_well=lambda d: (d.rescued_reads / d.n).round(1))
              .to_string(index=False))
        tab.to_csv(OUT / f"d07_{side.lower()}_rescue.csv", index=False)

    pd.Series(sam, name="rescued_reads").rename_axis("sample").to_csv(
        OUT / "d08_rescued_reads_per_sample.csv")
    print(f"\nwrote per-sample rescue counts -> {OUT}/d08_rescued_reads_per_sample.csv")
