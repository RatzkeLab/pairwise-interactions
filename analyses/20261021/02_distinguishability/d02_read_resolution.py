"""Stage 5: at what separation can a nanopore read tell two strains apart? Measured.

The 15 bp design rule came from 20260630's monocultures, which were not pure (they were
shot against believed-dead wells, some of which grew), and its 10 bp MIN_RESOLVABLE_BP
from the same table. 20260907 has true monocultures, two per strain on different plates,
so the question can be asked cleanly:

  for strain X with two agreeing monos (tier A), take reads from mono well 1 and race
  them against X's reference built from mono well 2 -- an independent well, so the reads
  never helped build the sequence they are scored against -- and against every other
  strain Y whose reference lies within MAX_BP of X's.

A read goes to whichever strain it is closer to (edit distance of the read's interior as
an infix, so read-end trimming cannot bias it); ties split. The share of X's reads going
to Y is the cross-talk a community containing both would suffer. Two races are run:

  misassigned             haplotype-aware: distance to a strain = to its nearest operon
                          haplotype (reference_db_20260907_haplotypes.fasta). This is how
                          reads should be assigned in the community experiment.
  misassigned_single_ref  X's single consensus from the OTHER mono well vs Y's consensus --
                          fully independent of the reads being scored, but blind to
                          linked operon variants.

Result (clean strains: <=10 het sites, not mixed_source, both sides). Cross-talk falls
off a cliff between 5-6 bp (p90 16%) and 7-8 bp (p90 0%), and from there on sits at a
flat ~0.7% p90 floor that does NOT improve with distance (it is the same at 19-40 bp):
chimeras and stray reads, not resolution. So with base-exact references, 8 bp separates
as well as 15 does. The earlier 10/15 bp rules came from references built over a handful
of reads and from impure monos.

The tail that remains is strain-specific, not distance-specific: E23 (23 het sites) sends
~45% of its reads to strains 11 bp away, reproducibly in both of its monos. Strains like
that are most likely stable mixtures of two close organisms, and are excluded from the
clean table rather than averaged in.

Outputs:
  d02_pair_crosstalk.csv           per (X, Y): reads, correct / tie / wrong, misassigned
  d02_resolution_by_bp.csv         binned by separation, all tier-A strains
  d02_resolution_by_bp_clean.csv   binned, both strains clean -- THE calibration table
  d02_resolution.png
"""
import pickle
import sys
from multiprocessing import Pool
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import config  # noqa: E402

import edlib
import numpy as np
import pandas as pd
from Bio import SeqIO

import deep_consensus as dc
from io_utils import load_layout, sample_fastq_path

OUT = config.DIST_OUT
MAX_BP = 40
READS_PER_WELL = 300
READ_TRIM = 25
BINS = [-1, 0, 2, 4, 6, 8, 10, 14, 19, 24, 29, 40]
CLEAN_MAX_HET = 10


def _ed(read, ref):
    return edlib.align(read, ref, mode="HW", task="distance")["editDistance"]


def _job(args):
    x, path, ref_x, haps_x, neighbours, anchor, seed = args
    reads = dc.prepare_reads(path, anchor, seed=seed, cap=READS_PER_WELL)
    core = [s[READ_TRIM:-READ_TRIM] for s, _ in reads]
    dx = np.array([_ed(r, ref_x) for r in core])
    keep = dx <= dc.OUTLIER_BP            # far contaminants are not X's reads
    core = [c for c, k in zip(core, keep) if k]; dx = dx[keep]
    dxh = np.array([min(_ed(r, h) for h in haps_x) for r in core])
    out = []
    n = len(core)
    if n == 0:
        return out
    for y, ref_y, haps_y, bp, plain in neighbours:
        dy = np.array([_ed(r, ref_y) for r in core])
        dyh = np.array([min(_ed(r, h) for h in haps_y) for r in core])
        c, t, w = int((dxh < dyh).sum()), int((dxh == dyh).sum()), int((dxh > dyh).sum())
        cs, ts, ws = int((dx < dy).sum()), int((dx == dy).sum()), int((dx > dy).sum())
        out.append(dict(strain_x=x, strain_y=y, bp=bp, plain_bp=plain, n_reads=n,
                        correct=c, tie=t, wrong=w, misassigned=(w + 0.5 * t) / n,
                        misassigned_single_ref=(ws + 0.5 * ts) / n))
    return out


def main(processes=48):
    ref = pd.read_csv(config.REFDB_OUT / "s01_strain_reference.csv").set_index("strain")
    pairs = pd.read_csv(OUT / "d01_pairs.csv")
    plain = {r.id: str(r.seq) for r in SeqIO.parse(
        config.REFDB_OUT / "reference_db_20260907.fasta", "fasta")}
    haps = {}
    for r in SeqIO.parse(config.REFDB_OUT / "reference_db_20260907_haplotypes.fasta", "fasta"):
        haps.setdefault(r.id.rsplit("_h", 1)[0], []).append(str(r.seq))
    with open(config.REFDB_OUT / "m02_mono_populations.pkl", "rb") as fh:
        mono = pickle.load(fh)
    obs = pd.read_csv(config.REFDB_OUT / "s02_strain_observations.csv")
    lay = load_layout(config.LAYOUT_CSV).set_index("sample_id")
    anchor = plain[config.ORIENT_ANCHOR_STRAIN]

    near = pairs[pairs.bp <= MAX_BP]
    nb = {}
    for r in near.itertuples():
        nb.setdefault(r.strain_a, []).append((r.strain_b, r.bp, r.plain_bp))
        nb.setdefault(r.strain_b, []).append((r.strain_a, r.bp, r.plain_bp))

    jobs = []
    for x in ref.index[ref.tier == "A"]:
        if x not in nb:
            continue
        wells = obs[(obs.strain == x) & (obs.cluster == 0) & (obs.kind == "mono_dominant")]
        wells = sorted(wells.sample_id.unique())
        if len(wells) < 2:
            continue
        # each mono's reads are scored against the OTHER mono's consensus
        for w_reads, w_ref in ((wells[0], wells[1]), (wells[1], wells[0])):
            ref_x = mono[w_ref][0]["consensus"]
            path = sample_fastq_path(config.DEMUX_DIR, lay.at[w_reads, "dest_plate"],
                                     lay.at[w_reads, "dest_well"])
            ngh = [(y, plain[y], haps[y], bp, pb) for y, bp, pb in nb[x]]
            jobs.append((x, str(path), ref_x, haps[x], ngh, anchor, len(jobs)))
    print(f"{len(jobs)} mono wells from {len(jobs) // 2} tier-A strains, "
          f"{sum(len(j[4]) for j in jobs):,} (well, neighbour) races")
    with Pool(processes) as pool:
        res = pool.map(_job, jobs, chunksize=1)
    df = pd.DataFrame([r for rr in res for r in rr])
    comp = pd.concat([pairs.rename(columns={"strain_a": "strain_x", "strain_b": "strain_y"}),
                      pairs.rename(columns={"strain_b": "strain_x", "strain_a": "strain_y"})])
    df = df.merge(comp[["strain_x", "strain_y", "subs", "indel_bp", "hp_indel_bp"]],
                  on=["strain_x", "strain_y"], how="left")
    # the two directions of the same well pair are independent reads; keep both

    clean = set(ref.index[(ref.het_sites <= CLEAN_MAX_HET) & ~ref.mixed_source.astype(bool)])
    df["clean_pair"] = df.strain_x.isin(clean) & df.strain_y.isin(clean)
    df.to_csv(OUT / "d02_pair_crosstalk.csv", index=False)
    tab = binned(df)
    tab.to_csv(OUT / "d02_resolution_by_bp.csv", index=False)
    tab_c = binned(df[df.clean_pair])
    tab_c.to_csv(OUT / "d02_resolution_by_bp_clean.csv", index=False)
    print(f"clean pairs only ({df[df.clean_pair].strain_x.nunique()} X strains):")
    print(tab_c[tab_c.distance == "bp"].drop(columns="distance").to_string(index=False))
    plot(df)
    return df, tab_c


def binned(df):
    rows = []
    for col in ("bp", "plain_bp"):
        for b, g in df.groupby(pd.cut(df[col], BINS), observed=True):
            rows.append(dict(distance=col, bin=str(b), n_races=len(g),
                             reads=int(g.n_reads.sum()),
                             median_misassigned=round(g.misassigned.median(), 4),
                             p90_misassigned=round(g.misassigned.quantile(0.9), 4),
                             p95_misassigned=round(g.misassigned.quantile(0.95), 4),
                             p90_single_ref=round(g.misassigned_single_ref.quantile(0.9), 4),
                             pooled_correct=round(g.correct.sum() / g.n_reads.sum(), 4),
                             pooled_tie=round(g.tie.sum() / g.n_reads.sum(), 4)))
    return pd.DataFrame(rows)


def plot(df):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    c = df[df.clean_pair & (df.bp >= 1)].groupby("bp").misassigned
    a = df[df.bp >= 1].groupby("bp").misassigned
    ax[0].plot(a.quantile(0.9).index, a.quantile(0.9).clip(lower=1e-3), color="0.6",
               lw=1, ls="--", label="p90, all tier-A strains")
    ax[0].plot(c.quantile(0.9).index, c.quantile(0.9).clip(lower=1e-3), color="#C44E52",
               lw=2, label="p90, clean strains")
    ax[0].plot(c.median().index, c.median().clip(lower=1e-3), color="#4C72B0", lw=2,
               label="median, clean strains")
    for v, ls in ((5, ":"), (15, "--")):
        ax[0].axvline(v, color="0.3", ls=ls, lw=1)
    ax[0].set(yscale="log", xlabel="separation, het-aware bp", ylim=(1e-3, 0.6),
              ylabel="share of X's reads assigned to Y (floored 0.001)",
              title="cross-talk vs separation (tier-A monos, independent reference)")
    ax[0].legend()
    sub = df[df.clean_pair & (df.bp >= 3) & (df.bp < 25)].copy()
    sub["sub_share"] = sub.subs / sub[["subs", "indel_bp"]].sum(axis=1).clip(lower=1)
    sc = ax[1].scatter(sub.bp + np.random.default_rng(0).uniform(-.3, .3, len(sub)),
                       sub.misassigned.clip(lower=1e-3), s=6, alpha=0.5,
                       c=sub.sub_share, cmap="viridis", vmin=0, vmax=1)
    fig.colorbar(sc, ax=ax[1], label="share of differences that are substitutions")
    ax[1].set(yscale="log", ylim=(9e-4, 0.6), xlabel="separation, het-aware bp",
              ylabel="share misassigned (floored at 0.001)", title="each race, clean strains, 3-24 bp")
    fig.tight_layout(); fig.savefig(OUT / "d02_resolution.png", dpi=150)


if __name__ == "__main__":
    main()
