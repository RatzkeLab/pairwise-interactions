"""Stage 2: an exact consensus for each strain in each of its pair wells.

Each pair well is a second, independent look at the same source well, with a different
partner and on a different destination plate. If a strain's monocultures say sequence S,
then S should reappear -- exactly -- in its pair wells.

Reads cannot simply be clustered here: pairs were chosen >=15 bp apart, and greedy read
clustering only separates organisms more than ~60 bp apart (a Q24 read sits a median
~12 bp from its own consensus). So reads are assigned to the nearest SEED -- every
candidate sequence either strain could be -- and each group then gets its own consensus,
which is compared back to the seeds.

Seeds per strain: the dominant population of each of its mono wells, any population
holding >=10% of one of its mono wells (a strain overgrown in one mono still shows up),
and the 20260921 consensus (the only candidate for strains whose monos gave no reads).

A group is credited to a strain only if its consensus is within CREDIT_MAX_BP of one of
that strain's seeds AND closer to them than to the partner's. The first half drops groups
that drifted 18-51 bp from every seed -- a minority strain's group at a few % of the well
is mostly chimeras and stray reads, and its consensus is not an observation of anything.
The second half matters under competitive exclusion: at a 97:3 split
with the partners 20 bp apart, the minority's group can be mostly misassigned reads of
the winner, and its consensus is then the winner's sequence, not evidence about the
minority.

Outputs:
  p01_pair_observations.csv   one row per (well, group)
  p02_pair_observations.pkl   group consensus, IUPAC version, het sites, top reads
"""
import collections
import pickle
import sys
import time
from multiprocessing import Pool
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import config  # noqa: E402

import numpy as np
import pandas as pd
from Bio import SeqIO

import deep_consensus as dc
from io_utils import load_layout, sample_fastq_path

OUT = config.REFDB_OUT
PAIR_TYPES = ("pair", "techrep", "ratio", "density")
SEED_MIN_FRAC = 0.10
SEED_DEDUP_BP = 10       # seeds are only for read assignment, so near-duplicates collapse
SHARED_BP = 5            # a seed this close to a partner seed cannot split reads
CREDIT_MAX_BP = 15
KEEP_READS = 30


def build_seeds():
    """strain -> [(seed_id, seq)]"""
    with open(OUT / "m02_mono_populations.pkl", "rb") as fh:
        mono = pickle.load(fh)
    wells = pd.read_csv(OUT / "m01_mono_wells.csv")
    prior = {r.id: str(r.seq) for r in SeqIO.parse(config.PRIOR_CONSENSUS_FASTA, "fasta")}
    seeds = collections.defaultdict(list)
    for r in wells.itertuples():
        for k, p in enumerate(mono.get(r.sample_id, [])):
            if k == 0 or p["frac"] >= SEED_MIN_FRAC:
                seeds[r.strain].append((f"{r.sample_id}#{k}", p["consensus"]))
    for s, seq in prior.items():
        seeds[s].append(("prior", seq))
    out = {}
    for s, lst in seeds.items():
        keep = []
        for sid, seq in lst:
            if all(dc.core_bp(seq, k[1]) > SEED_DEDUP_BP for k in keep):
                keep.append((sid, seq))
        out[s] = keep
    return out


def _job(args):
    sid, plate, wtype, a, b, seeds_a, seeds_b, path, anchor, seed = args
    reads = dc.prepare_reads(path, anchor, seed=seed)
    rows, obs = [], []
    if not reads:
        return rows, obs
    labelled = [(a, i, s) for i, (_, s) in enumerate(seeds_a)] + \
               [(b, i, s) for i, (_, s) in enumerate(seeds_b)]
    groups = collections.defaultdict(list)
    n_off = 0
    for seq, q in reads:
        d = [dc._ed(seq, s, dc.OUTLIER_BP) for _, _, s in labelled]
        d = [x if x >= 0 else 10 ** 9 for x in d]
        j = int(np.argmin(d))
        if d[j] > dc.OUTLIER_BP:
            n_off += 1
        else:
            groups[j].append((seq, q))
    n = len(reads)
    for j, mem in groups.items():
        strain, k, seedseq = labelled[j]
        partner = b if strain == a else a
        if len(mem) < max(dc.MIN_POP_READS, dc.MIN_POP_FRAC * n):
            continue
        cons = dc.top_q_consensus(mem)
        het_n, het = dc.het_sites(cons, [s for s, _ in sorted(mem, key=lambda x: -x[1])])
        own = min(dc.core_bp(cons, s) for st, _, s in labelled if st == strain)
        other = min(dc.core_bp(cons, s) for st, _, s in labelled if st == partner)
        seed_gap = min(dc.core_bp(seedseq, s) for st, _, s in labelled if st == partner)
        credited = own <= CREDIT_MAX_BP and own < other and seed_gap > SHARED_BP
        rows.append(dict(sample_id=sid, plate=plate, well_type=wtype, strain=strain,
                         partner=partner, seed=(seeds_a if strain == a else seeds_b)[k][0],
                         n_reads=len(mem), frac=round(len(mem) / n, 4), n_well_reads=n,
                         n_off_target=n_off, to_own_seed_bp=own, to_partner_seed_bp=other,
                         seed_to_partner_bp=seed_gap, het_sites=het_n, credited=credited))
        obs.append(dict(sample_id=sid, strain=strain, consensus=cons,
                        iupac=dc.iupac(cons, het), het=het, n_reads=len(mem),
                        frac=len(mem) / n, plate=plate, credited=credited, to_own_seed_bp=own,
                        reads=sorted(mem, key=lambda x: -x[1])[:KEEP_READS]))
    return rows, obs


def main(processes=48):
    lay = load_layout(config.LAYOUT_CSV)
    pw = lay[(lay.plate_role == "primary_sequenced") & lay.well_type.isin(PAIR_TYPES)]
    seeds = build_seeds()
    anchor = {r.id: str(r.seq) for r in SeqIO.parse(config.PRIOR_CONSENSUS_FASTA, "fasta")}[
        config.ORIENT_ANCHOR_STRAIN]
    print(f"seeds per strain: {pd.Series({k: len(v) for k, v in seeds.items()}).value_counts().sort_index().to_dict()}")
    jobs = []
    for i, r in enumerate(pw.itertuples()):
        p = sample_fastq_path(config.DEMUX_DIR, r.dest_plate, r.dest_well)
        if p is None or r.strain1 not in seeds or r.strain2 not in seeds:
            continue
        jobs.append((r.sample_id, r.dest_plate, r.well_type, r.strain1, r.strain2,
                     seeds[r.strain1], seeds[r.strain2], str(p), anchor, i))
    print(f"{len(pw)} pair-type wells, {len(jobs)} with a fastq and seeds for both strains")
    t = time.time()
    with Pool(processes) as pool:
        res = pool.map(_job, jobs, chunksize=2)
    print(f"  done in {time.time() - t:.0f}s")
    rows = [r for rr, _ in res for r in rr]
    obs = [o for _, oo in res for o in oo]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "p01_pair_observations.csv", index=False)
    with open(OUT / "p02_pair_observations.pkl", "wb") as fh:
        pickle.dump(obs, fh)
    print(f"  {len(df)} strain groups in {df.sample_id.nunique()} wells; "
          f"{df.credited.sum()} credited to their strain")
    c = df[df.credited]
    print(f"  credited group consensus vs its own seed: "
          f"{pd.cut(c.to_own_seed_bp, [-1, 0, 1, 2, 5, 15, 2000]).value_counts().sort_index().to_dict()}")
    return df, obs


if __name__ == "__main__":
    main()
