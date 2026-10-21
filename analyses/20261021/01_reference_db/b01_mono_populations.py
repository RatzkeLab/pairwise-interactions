"""Stage 1: every population in every monoculture well, with a base-exact consensus.

A mono well should hold one population. When it holds more, the extra ones are
contamination of THAT well (if the strain's other mono, on another plate, is clean) or
of the source well (if both show it). Stage 3 decides which; this stage only records
what each well holds, so the decision can be audited.

Outputs (01_reference_db/outputs/):
  m01_mono_wells.csv          one row per mono well: depth, populations, dominant share,
                              het sites in the dominant population
  m01_mono_populations.csv    one row per population (consensus length, reads, share)
  m02_mono_populations.pkl    consensus sequences + each population's top-qscore reads,
                              which stage 3 pools across a strain's wells
"""
import pickle
import sys
import time
from multiprocessing import Pool
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import config  # noqa: E402  (also puts scripts/ on sys.path)

import pandas as pd
from Bio import SeqIO

import deep_consensus as dc
from io_utils import load_layout, sample_fastq_path

OUT = config.REFDB_OUT
KEEP_READS = 60          # top-qscore reads kept per population for pooling in stage 3


def _anchor():
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(config.PRIOR_CONSENSUS_FASTA, "fasta")}
    return seqs[config.ORIENT_ANCHOR_STRAIN]


def _job(args):
    sid, strain, plate, vol, path, anchor, seed = args
    reads = dc.prepare_reads(path, anchor, seed=seed)
    pops, n_out = dc.well_populations(reads)
    for p in pops:
        p["het_n"], p["het"] = dc.het_sites(p["consensus"], [s for s, _ in p["reads"]])
        p["iupac"] = dc.iupac(p["consensus"], p["het"])
    row = dict(sample_id=sid, strain=strain, plate=plate, vol_nL=vol, n_reads=len(reads),
               n_outlier_reads=n_out, n_pops=len(pops))
    if pops:
        top = pops[0]
        n_het, sites = top["het_n"], top["het"]
        row.update(top_frac=round(top["frac"], 4), top_reads=top["n_reads"],
                   second_frac=round(pops[1]["frac"], 4) if len(pops) > 1 else 0.0,
                   second_vs_top_bp=dc.core_bp(pops[1]["consensus"], top["consensus"])
                   if len(pops) > 1 else None,
                   top_len=len(top["consensus"]), het_sites=n_het,
                   het_detail=";".join(f"{p}:{f}:{a}" for p, f, a in sites) if sites else "")
    slim = [dict(consensus=p["consensus"], iupac=p["iupac"], het=p["het"],
                 n_reads=p["n_reads"], frac=p["frac"],
                 reads=sorted(p["reads"], key=lambda x: -x[1])[:KEEP_READS]) for p in pops]
    return row, sid, slim


def main(processes=48):
    OUT.mkdir(parents=True, exist_ok=True)
    lay = load_layout(config.LAYOUT_CSV)
    mono = lay[(lay.plate_role == "primary_sequenced") & (lay.well_type == "mono")]
    anchor = _anchor()
    jobs = []
    for i, r in enumerate(mono.itertuples()):
        p = sample_fastq_path(config.DEMUX_DIR, r.dest_plate, r.dest_well)
        if p is not None:
            jobs.append((r.sample_id, r.strain1, r.dest_plate, r.vol1_nL, str(p), anchor, i))
    print(f"{len(mono)} mono wells, {len(jobs)} with a fastq")
    t = time.time()
    with Pool(processes) as pool:
        res = pool.map(_job, jobs, chunksize=2)
    print(f"  populations built in {time.time() - t:.0f}s")

    rows = [r for r, _, _ in res]
    pops = {sid: p for _, sid, p in res}
    wells = pd.DataFrame(rows)
    wells.to_csv(OUT / "m01_mono_wells.csv", index=False)
    prow = []
    for sid, plist in pops.items():
        for k, p in enumerate(plist):
            prow.append(dict(sample_id=sid, pop=k, n_reads=p["n_reads"],
                             frac=round(p["frac"], 4), length=len(p["consensus"])))
    pd.DataFrame(prow).to_csv(OUT / "m01_mono_populations.csv", index=False)
    with open(OUT / "m02_mono_populations.pkl", "wb") as fh:
        pickle.dump(pops, fh)

    has = wells[wells.n_pops > 0]
    print(f"  {len(has)} wells yielded at least one population "
          f"(>= {dc.MIN_POP_READS} reads, >= {dc.MIN_POP_FRAC:.0%})")
    print(f"  populations per well: {has.n_pops.value_counts().sort_index().to_dict()}")
    print(f"  dominant share: median {has.top_frac.median():.3f}, "
          f"<0.9 in {(has.top_frac < 0.9).mean():.1%}, <0.5 in {(has.top_frac < 0.5).mean():.1%}")
    print(f"  het sites in dominant population (MAF >= {dc.HET_MAF}): "
          f"{has.het_sites.value_counts().sort_index().to_dict()}")
    return wells, pops


if __name__ == "__main__":
    main()
