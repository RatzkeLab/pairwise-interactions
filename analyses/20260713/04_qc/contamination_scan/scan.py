"""Is there low-level contamination across 20260721's wells, and could it explain the mix-up?

Hypothesis under test: DNA carried over before PCR, so each well amplifies mostly the right
thing plus some of something else.

Design decisions that matter, each of which reverses the reading if skipped:

1.  **Split by well type.** Most wells are PAIR wells and contain two organisms *by design*.
    A "second organism detected" rate over all wells measures the experiment, not
    contamination. Mono wells (26 here, 34 in the control) expect exactly 1 group; pair wells
    expect exactly 2; anything above that is the contamination signal.

2.  **Work on 16S GROUPS, not strain names.** Only 35 of 294 collection strains have a unique
    16S -- db entries `D3` and `D13` sit 1 bp from `N19`. A read best-hitting its 1 bp twin is
    the resolution limit, not contamination. Single-linkage grouping at MIN_RESOLVABLE_BP
    collapses those (87 entries -> 73 groups).

3.  **20260630 is the control, run identically.** A contamination rate means nothing without
    the rate in an experiment whose labels are known good.

4.  **Low-identity reads are tested for identity, not just counted.** Essentially every read
    maps to *something* (16S is conserved, so even a foreign organism hits a collection member
    at ~0.8), so a foreign contaminant shows up as reads with poor identity -- but so do
    chimeras and truncated reads. The two are separated by asking whether the low-identity
    reads resemble *each other*: a real foreign organism is one tight cluster, junk is diffuse.

5.  Depth is ~41 reads/well, so a 5% contaminant is ~2 reads. Per-well calls are weak; the
    power is in the aggregate over ~2600 wells.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
EXPS = HERE.parents[2]
sys.path.insert(0, str(EXPS / "strain_identity_qc"))
sys.path.insert(0, str(EXPS / "shared_pipelines"))
from qc_compare import identity
from genomic_ml import COLOR_BLUE, COLOR_RED, COLOR_GRID, COLOR_TEXT_SECONDARY

OUT = HERE / "outputs"
FIG = OUT / "figures"
for d in (OUT, FIG):
    d.mkdir(parents=True, exist_ok=True)

DB = Path("/home/rl/scripts/karl/merge_consensus_sequences/collapse_naive_updated3_15diff/"
          "corroborated_db_filtered_min5_edited.fasta")
MIN_RESOLVABLE_BP = 10
LOW_IDENTITY = 0.88
MIN_COMPONENT_READS, MIN_COMPONENT_FRAC = 2, 0.05
EXPECTED_GROUPS = {"mono": 1, "pair": 2}

MV = "analysis/mapping_validation/outputs"
EXPERIMENTS = ["20260721", "20260630"]


def group_16s():
    seqs = {r.id.split("|")[0].removesuffix("_consensus"): str(r.seq).upper()
            for r in SeqIO.parse(DB, "fasta")}
    names = sorted(seqs)
    parent = {n: n for n in names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]] ; x = parent[x]
        return x
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if identity(seqs[a], seqs[b])[1] < MIN_RESOLVABLE_BP:
                parent[find(a)] = find(b)
    g = {n: find(n) for n in names}
    print(f"  {len(names)} db entries -> {len(set(g.values()))} 16S groups")
    return g


def composition(exp, groups):
    base = EXPS / exp / MV
    d = pd.read_csv(base / "04_minimap2_read_besthit_corroborated_db.csv")
    d = d.drop_duplicates("read_id")
    meta = pd.read_csv(base / "01_samples_gt5reads.csv")[["sample_id", "well_type"]]
    d = d.merge(meta, on="sample_id", how="left")
    d["group"] = d["best_hit_strain"].map(groups)
    d["low_id"] = d["identity"] < LOW_IDENTITY

    rows = []
    for sid, s in d.groupby("sample_id"):
        good = s[~s.low_id & s.group.notna()]
        if good.empty:
            continue
        vc = good.group.value_counts()
        keep = vc[(vc >= MIN_COMPONENT_READS) & (vc / len(good) >= MIN_COMPONENT_FRAC)]
        rows.append({"sample_id": sid, "well_type": s.well_type.iloc[0], "n_reads": len(s),
                     "frac_low_identity": s.low_id.mean(),
                     "n_groups": int(len(keep)), "dominant_frac": vc.iloc[0] / len(good),
                     "groups": ",".join(keep.index[:4])})
    w = pd.DataFrame(rows)
    w["expected_groups"] = w.well_type.map(EXPECTED_GROUPS)
    w["excess_groups"] = w.n_groups - w.expected_groups
    w["experiment"] = exp
    return w, d


def low_identity_character(exp, d, n_sample=250, seed=0):
    """Are the low-identity reads one foreign organism, or just bad reads?

    A genuine foreign contaminant is a single sequence: sampled low-identity reads would sit
    close to each other while sitting far from the collection. Chimeras and truncations are
    unrelated to each other, so they stay far from everything including themselves.
    """
    rng = np.random.default_rng(seed)
    low = d[d.low_id]
    if low.empty:
        return {}
    # fastq headers are "<sample_id>__<read_id>", not the bare read_id -- matching on the bare
    # id silently finds nothing and the whole test comes back empty
    key = low.sample_id + "__" + low.read_id
    ids = set(rng.choice(key.values, size=min(n_sample, len(low)), replace=False))
    seqs = []
    fq = EXPS / exp / MV / "04_reads_combined.fastq"
    with open(fq) as fh:
        for rec in SeqIO.parse(fh, "fastq"):
            if rec.id in ids:
                seqs.append(str(rec.seq).upper())
                if len(seqs) == len(ids):
                    break
    if len(seqs) < 20:
        return {}
    idx = rng.choice(len(seqs), size=min(120, len(seqs)), replace=False)
    sub = [seqs[i] for i in idx]
    dists = []
    for i in range(len(sub)):
        for j in range(i + 1, len(sub)):
            dists.append(identity(sub[i], sub[j])[0])
    dists = np.array(dists)
    return {"n_low_identity_reads": int(len(low)),
            "n_sampled": len(sub),
            "median_pairwise_identity_among_low_reads": round(float(np.median(dists)), 4),
            "frac_pairs_above_0.95": round(float((dists > 0.95).mean()), 4),
            "median_read_len_low": float(low.read_len.median()),
            "median_read_len_all": float(d.read_len.median())}


def main():
    groups = group_16s()
    wells, chars = [], []
    for exp in EXPERIMENTS:
        print(f"\n=== {exp} ===")
        w, d = composition(exp, groups)
        wells.append(w)
        for wt in ("mono", "pair"):
            s = w[w.well_type == wt]
            if s.empty:
                continue
            print(f"  {wt:5} wells n={len(s):5d}  expected {EXPECTED_GROUPS[wt]} group(s)  "
                  f"observed median {s.n_groups.median():.0f}  "
                  f"| EXCESS groups: {100*(s.excess_groups > 0).mean():5.1f}% of wells  "
                  f"| median low-identity read frac {s.frac_low_identity.median():.3f}")
        c = low_identity_character(exp, d)
        c["experiment"] = exp
        chars.append(c)
        print("  low-identity read character:",
              {k: v for k, v in c.items() if k != "experiment"})
    W = pd.concat(wells, ignore_index=True)
    W.to_csv(OUT / "c01_well_composition.csv", index=False)
    C = pd.DataFrame(chars)
    C.to_csv(OUT / "c03_low_identity_read_character.csv", index=False)
    summ = (W.groupby(["experiment", "well_type"])
              .agg(n_wells=("sample_id", "size"), median_groups=("n_groups", "median"),
                   pct_excess=("excess_groups", lambda s: round(100 * (s > 0).mean(), 2)),
                   median_dominant_frac=("dominant_frac", "median"),
                   median_frac_low_identity=("frac_low_identity", "median"))
              .reset_index())
    summ.to_csv(OUT / "c00_summary.csv", index=False)
    print("\n=== SUMMARY ===")
    print(summ.to_string(index=False))
    make_figure(W)
    print("\nwritten ->", OUT)
    return summ


def make_figure(W):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5))
    colors = {"20260721": COLOR_RED, "20260630": COLOR_BLUE}

    ax = axes[0]
    p = W[W.well_type == "pair"]
    for e, s in p.groupby("experiment"):
        vc = s.n_groups.value_counts(normalize=True).sort_index()
        ax.plot(vc.index, 100 * vc.values, "o-", color=colors[e], label=f"{e} (n={len(s)})")
    ax.axvline(2, color=COLOR_TEXT_SECONDARY, ls="--", lw=1)
    ax.text(2.05, ax.get_ylim()[1] * .9, "expected", fontsize=8, color=COLOR_TEXT_SECONDARY)
    ax.set_xlabel("distinct 16S groups detected in a pair well")
    ax.set_ylabel("% of wells"); ax.legend(frameon=False, fontsize=8)
    ax.set_title("Pair wells: a 3rd organism would sit right of the line", fontsize=9)

    ax = axes[1]
    m = W[W.well_type == "mono"]
    for e, s in m.groupby("experiment"):
        vc = s.n_groups.value_counts(normalize=True).sort_index()
        ax.plot(vc.index, 100 * vc.values, "o-", color=colors[e], label=f"{e} (n={len(s)})")
    ax.axvline(1, color=COLOR_TEXT_SECONDARY, ls="--", lw=1)
    ax.set_xlabel("distinct 16S groups detected in a mono well")
    ax.set_ylabel("% of wells"); ax.legend(frameon=False, fontsize=8)
    ax.set_title("Mono wells: the clean test\n(only 26 / 34 wells — weak but direct)", fontsize=9)

    ax = axes[2]
    for e, s in W.groupby("experiment"):
        ax.hist(s.frac_low_identity, bins=np.linspace(0, .6, 41), density=True,
                histtype="step", lw=2, color=colors[e], label=e)
    ax.set_xlabel("fraction of a well's reads below 0.88 identity to any reference")
    ax.set_ylabel("density"); ax.legend(frameon=False, fontsize=8)
    ax.set_title("Poorly-matching reads\nforeign DNA, chimeras, or truncations", fontsize=9)

    fig.suptitle("Contamination scan — 20260721 vs 20260630 control, split by well type",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "c01_contamination_scan.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
