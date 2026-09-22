"""Build one 16S consensus per strain from the wells that strain was put into.

Shared version of 20260721/03_recreate_reference_db/per_strain_consensus.ipynb.
The logic changes in one important way, because 20260907's design changed: every
strain now has monoculture wells, two of them, on different plates.

20260721 had no true monocultures, so the only way to find a strain's sequence was
that it recurs across wells while its partners change -- take the sequence spanning
the most wells. That works but it is a majority argument, and it cannot tell a
strain from a contaminant that travelled with it.

Here the monocultures are the anchor and the pair wells are the corroboration:

    the strain is what its monocultures say it is;
    the pair wells then confirm it, and how often they do is the quality score.

That ordering matters for the question this run was designed to answer. A strain
whose two monos agree with each other but whose pair wells show something else is
a labelling problem. A strain whose two monos disagree is a contaminated well.
Both are invisible to a pure majority vote, which would simply report whichever
sequence occurred most.
"""
import collections
import numpy as np
import pandas as pd
import edlib

from io_utils import load_layout

CROSS_THR = 0.92          # identity to merge per-well representatives across wells
MIN_CORROBORATING = 2     # wells that must carry a sequence before it is a strain

# A strain whose monocultures produced no usable reads is not lost: its sequence still
# recurs across every pair well it was put into, while its partners change. That is the
# only evidence 20260721 ever had, and it worked. It is accepted here too, but held to a
# stricter bar than the monoculture route, because a pair well shows two organisms and the
# recurring one could in principle be a travelling contaminant rather than the strain --
# requiring several DIFFERENT partners across several PLATES is what makes that unlikely.
MIN_PAIR_ONLY_WELLS = 3
MIN_PAIR_ONLY_PLATES = 2


def identity(a, b):
    r = edlib.align(a, b, mode="NW", task="distance")
    return 1 - r["editDistance"] / max(len(a), len(b))


def majority_consensus(seqs):
    """Column-wise majority over an alignment anchored on the medoid.

    Anchoring on one sequence rather than a full MSA keeps this O(n) alignments.
    Insertions are kept only where a majority of sequences carry one, which is what
    removes ONT's homopolymer-length noise.
    """
    if len(seqs) == 1:
        return seqs[0]
    anchor = min(seqs, key=lambda a: sum(
        edlib.align(a, b, mode="NW", task="distance")["editDistance"] for b in seqs))
    cols = [collections.Counter() for _ in range(len(anchor) + 1)]
    ins = [collections.Counter() for _ in range(len(anchor) + 1)]
    for s in seqs:
        r = edlib.align(s, anchor, mode="NW", task="path")
        na = edlib.getNiceAlignment(r, s, anchor)
        q, t = na["query_aligned"], na["target_aligned"]
        pos, pending = 0, []
        for qc, tc in zip(q, t):
            if tc == "-":
                pending.append(qc)
            else:
                if pending:
                    ins[pos]["".join(pending)] += 1
                    pending = []
                cols[pos][qc] += 1
                pos += 1
        if pending:
            ins[pos]["".join(pending)] += 1
    n = len(seqs)
    out = []
    for i in range(len(anchor) + 1):
        if ins[i]:
            block, cnt = ins[i].most_common(1)[0]
            if cnt > n / 2:
                out.append(block.replace("-", ""))
        if i < len(anchor) and cols[i]:
            base, cnt = cols[i].most_common(1)[0]
            if base != "-" and cnt > 0:
                out.append(base)
    return "".join(out)


def _cross_cluster(reps):
    """reps: list of (seq, sample_id, well_type, plate). Greedy identity clustering."""
    clusters = []
    for seq, sid, wt, pl in reps:
        best, bj = -1.0, -1
        for j, c in enumerate(clusters):
            v = identity(seq, c["seed"])
            if v > best: best, bj = v, j
        if best >= CROSS_THR:
            clusters[bj]["mem"].append((seq, sid, wt, pl))
        else:
            clusters.append({"seed": seq, "mem": [(seq, sid, wt, pl)]})
    for c in clusters:
        c["wells"] = {m[1] for m in c["mem"]}
        c["monos"] = {m[1] for m in c["mem"] if m[2] == "mono"}
        c["plates"] = {m[3] for m in c["mem"]}
    return clusters


def build(cfg, well_reps, well_df):
    """Returns (summary DataFrame, {strain: consensus sequence})."""
    lay = load_layout(cfg.layout_csv)
    if "plate_role" in lay.columns:
        lay = lay[lay.plate_role == "primary_sequenced"]
    meta = lay.set_index("sample_id")[["well_type", "dest_plate"]]

    # mono wells carry strain1 == strain2, so a plain loop over both columns counts every
    # monoculture twice -- which inflates `informative` and understates `coverage`
    strain_wells = collections.defaultdict(list)
    for r in lay.itertuples():
        for s in {r.strain1, r.strain2}:
            if isinstance(s, str) and s:
                strain_wells[s].append(r.sample_id)

    rows, consensus = [], {}
    for strain in sorted(strain_wells):
        wells = strain_wells[strain]
        reps, informative, mono_total = [], 0, 0
        for sid in wells:
            wt = meta.at[sid, "well_type"] if sid in meta.index else "?"
            pl = meta.at[sid, "dest_plate"] if sid in meta.index else -1
            if wt == "mono":
                mono_total += 1
            for seq in well_reps.get(sid, []):
                reps.append((seq, sid, wt, pl))
            if well_reps.get(sid):
                informative += 1

        if not reps:
            rows.append(dict(strain=strain, assigned_wells=len(wells), informative_wells=0,
                             mono_wells=mono_total, status="no_data"))
            continue

        clusters = _cross_cluster(reps)
        # monocultures decide; pair wells only break ties among clusters with no mono
        clusters.sort(key=lambda c: (len(c["monos"]), len(c["wells"])), reverse=True)
        top = clusters[0]
        second = clusters[1] if len(clusters) > 1 else None

        n_mono_seen = len({m[1] for m in reps if m[2] == "mono"})
        mono_agree = len(top["monos"])
        pair_wells = top["wells"] - top["monos"]
        coverage = len(top["wells"]) / informative if informative else 0.0

        n_plates = len(top["plates"])
        if mono_agree >= 2:
            status = "ok_monos_agree"
        elif mono_agree == 1 and len(pair_wells) >= MIN_CORROBORATING:
            status = "ok_mono_plus_pairs"
        elif mono_agree == 1 and n_mono_seen >= 2:
            status = "monos_disagree"          # both monos gave data but different sequences
        elif mono_agree == 1:
            status = "mono_only"
        elif (n_mono_seen == 0 and len(top["wells"]) >= MIN_PAIR_ONLY_WELLS
              and n_plates >= MIN_PAIR_ONLY_PLATES):
            status = "ok_pairs_corroborated"   # no mono data, but many partners agree
        elif len(top["wells"]) >= MIN_CORROBORATING:
            status = "pairs_only_weak"
        else:
            status = "single_well"

        consensus[strain] = majority_consensus([m[0] for m in top["mem"]])
        rows.append(dict(
            strain=strain, assigned_wells=len(wells), informative_wells=informative,
            mono_wells=mono_total, mono_wells_with_data=n_mono_seen,
            monos_supporting=mono_agree, pair_wells_supporting=len(pair_wells),
            plates_supporting=len(top["plates"]), coverage=round(coverage, 3),
            n_clusters=len(clusters),
            second_cluster_wells=len(second["wells"]) if second else 0,
            second_cluster_monos=len(second["monos"]) if second else 0,
            consensus_len=len(consensus[strain]), status=status))

    summary = pd.DataFrame(rows).sort_values("strain").reset_index(drop=True)
    return summary, consensus


# Statuses trusted enough to write out. The first two are anchored on a monoculture; the
# third is the pair-corroborated fallback, kept separate so it can be dropped or audited
# on its own (`report` prints its reference agreement next to the mono-anchored strains).
MONO_ANCHORED = ("ok_monos_agree", "ok_mono_plus_pairs")
GOOD = MONO_ANCHORED + ("ok_pairs_corroborated",)


def write_fasta(path, summary, consensus, statuses=GOOD):
    n = 0
    with open(path, "w") as fh:
        for r in summary.itertuples():
            if r.status in statuses and r.strain in consensus:
                fh.write(f">{r.strain} monos={r.monos_supporting}/{r.mono_wells_with_data} "
                         f"pairs={r.pair_wells_supporting} plates={r.plates_supporting} "
                         f"coverage={r.coverage} status={r.status}\n{consensus[r.strain]}\n")
                n += 1
    print(f"wrote {n} consensus sequences -> {path}")
    return n


def report(summary):
    print("\nper-strain consensus status:")
    for k, v in summary.status.value_counts().items():
        print(f"  {k:24s} {v:4d}")
    ok = summary[summary.status.isin(GOOD)]
    anchored = summary[summary.status.isin(MONO_ANCHORED)]
    print(f"\n{len(ok)} of {len(summary)} strains have a corroborated consensus "
          f"({len(ok)/len(summary):.1%})")
    print(f"  {len(anchored)} anchored on a monoculture, "
          f"{len(ok) - len(anchored)} corroborated by pair wells only "
          f"(>={MIN_PAIR_ONLY_WELLS} wells on >={MIN_PAIR_ONLY_PLATES} plates)")
    bad = summary[summary.status == "monos_disagree"]
    if len(bad):
        print(f"\n{len(bad)} strains whose two monocultures disagree -- a contaminated or "
              f"mislabelled well, worth looking at directly:")
        print(bad[["strain", "mono_wells_with_data", "monos_supporting",
                   "pair_wells_supporting", "second_cluster_wells"]].to_string(index=False))
    return summary.status.value_counts()
