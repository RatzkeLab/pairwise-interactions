"""Before trusting 1.3M rescued reads: does the primer-anchored parser agree with
minibar on the reads minibar *did* assign?

Those reads have a known ground-truth label (the file they were written to). If
the parser reproduces that label it is safe to apply to the discard pile.
"""
import random, collections
from multiprocessing import Pool
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from rescue_probe import call, PAIR, D          # reuse the exact same caller

DEMUX = Path("/home/rl/scripts/karl/data_links/interim/demultiplexing/20260907_demux")
UNFLIP = DEMUX / "unflipped"


def check(path):
    truth = path.stem                      # e.g. Plate01_B10
    with open(path) as fh:
        seqs = [l.strip() for i, l in enumerate(fh) if i % 4 == 1]
    if len(seqs) > 40:
        seqs = random.Random(0).sample(seqs, 40)
    ok = wrong = unparsed = 0
    wrong_ex = []
    for s in seqs:
        res = call(s)
        if res is None: unparsed += 1; continue
        hit = PAIR.get(res)
        if hit and hit[0] == truth: ok += 1
        else:
            wrong += 1
            wrong_ex.append((truth, hit[0] if hit else "no-such-well"))
    return ok, wrong, unparsed, wrong_ex


if __name__ == "__main__":
    counts = pd.read_csv(DEMUX / "summary/demultiplexed_read_counts.tsv", sep="\t")
    live = counts[counts.read_count >= 20]["sample"]
    files = [UNFLIP / f"{s}.fastq" for s in live]
    files = [f for f in files if f.exists()]
    files = random.Random(0).sample(files, min(400, len(files)))
    print(f"validating on {len(files)} wells that minibar populated (<=40 reads each)")
    with Pool(12) as p:
        res = p.map(check, files)
    ok = sum(r[0] for r in res); wrong = sum(r[1] for r in res)
    unp = sum(r[2] for r in res)
    n = ok + wrong + unp
    print(f"\n  reads checked        {n:,}")
    print(f"  parser agrees        {ok:,}  ({ok/n:.1%} of all, {ok/max(ok+wrong,1):.2%} of parsed)")
    print(f"  parser disagrees     {wrong:,}  ({wrong/max(ok+wrong,1):.2%} of parsed)")
    print(f"  parser cannot parse  {unp:,}  ({unp/n:.1%})  <- these stay unassigned, no harm")
    ex = [e for r in res for e in r[3]][:10]
    if ex:
        print("\n  examples of disagreement (minibar -> parser):")
        for a, b in ex: print(f"    {a}  ->  {b}")
    print(f"\n  false-assignment rate is the number that matters: "
          f"{wrong/max(ok+wrong,1):.2%} of rescued reads would go to the wrong well.")
