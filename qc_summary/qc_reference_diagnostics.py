"""How usable is each reference for IDENTIFYING a strain -- separate from its absolute identity?

`build_qc_matrix.py` scores support at a fixed 0.99 identity. That is the project's
same-organism convention, but it is unfair to `genome_16S`: those sequences come from NGS
assemblies while the queries are ONT amplicon consensus, so the whole identity distribution is
shifted down by platform, independent of whether the reference names the right organism. A
reference can be systematically offset and still perfectly discriminative.

Three threshold-free or threshold-swept views:

  A  **Threshold sweep.** % of labels supported at 0.95 ... 1.00. Shows how much of a verdict is
     the cutoff rather than the data.

  B  **Rank retrieval (threshold-free).** Rank every reference entry by identity to the query and
     ask where the CORRECT label lands. top-1 / top-3 / top-5 and median rank. This is the
     honest measure of "can this reference identify this strain": a reference whose identities
     are all low but whose correct answer still ranks first is perfectly usable.

  C  **Internal consistency.** Map each experiment's MONO consensus against its own PAIR
     consensus set. Same experiment, same platform, same pipeline, different well type -- so it
     removes every cross-platform and cross-reference confound and asks only whether the
     experiment agrees with itself.

Identity: strain_identity_qc.qc_compare.identity (infix, both strands, N as wildcard).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
EXPS = HERE.parent
KARL = EXPS.parent
sys.path.insert(0, str(EXPS / "strain_identity_qc"))
from qc_compare import identity

OUT = HERE / "outputs"
CORR = KARL / "merge_consensus_sequences" / "collapse_naive_updated3_15diff"
THRESHOLDS = [0.95, 0.96, 0.97, 0.98, 0.985, 0.99, 0.995, 1.0]


def load(path, multi=False):
    out = {}
    for r in SeqIO.parse(path, "fasta"):
        n = r.id.split("|")[0]
        if n.endswith("_consensus"):
            n = n[: -len("_consensus")]
        s = str(r.seq).upper()
        out.setdefault(n, []).append(s) if multi else out.setdefault(n, s)
    return out


def rank_table(queries, refs, exp, qname, rname):
    """Per label: identity to its own entry, the best identity, and the RANK of its own entry."""
    rows = []
    names = list(refs)
    for lab, q in queries.items():
        qs = q[0] if isinstance(q, list) else q
        ids = []
        for n in names:
            rs = refs[n]
            ids.append(max(identity(qs, s)[0] for s in (rs if isinstance(rs, list) else [rs])))
        ids = np.array(ids)
        order = np.argsort(-ids)
        ranked = [names[i] for i in order]
        if lab not in refs:
            rows.append({"experiment": exp, "consensus": qname, "reference": rname,
                         "label": lab, "in_reference": False, "self_identity": np.nan,
                         "self_rank": np.nan, "best_identity": float(ids.max()),
                         "best_match": ranked[0], "n_reference": len(names)})
            continue
        rows.append({"experiment": exp, "consensus": qname, "reference": rname, "label": lab,
                     "in_reference": True,
                     "self_identity": float(ids[names.index(lab)]),
                     "self_rank": int(ranked.index(lab)) + 1,
                     "best_identity": float(ids.max()), "best_match": ranked[0],
                     "n_reference": len(names)})
    return rows


def main():
    corr185 = load(CORR / "corroborated_db.fasta")
    corr87 = load(CORR / "corroborated_db_filtered_min5.fasta")
    gen = load(KARL / "Link to Karl" / "rDNA_16S_db_all_strains.fasta", multi=True)
    mp = pd.read_csv(KARL / "Link to Karl" / "final_genomic_tables" /
                     "mapping_384_well_plate_collection.csv")
    a2w = dict(zip(mp["assembly_name"].astype(str), mp["Well_souce_plate"]))
    gen_w = {}
    for k, v in gen.items():
        w = a2w.get(k)
        if w:
            gen_w.setdefault(w, []).extend(v)

    builds = {
        "20260630": {"pair": EXPS / "20260630/analysis/consensus2/strain_consensus_20260630.fasta",
                     "mono": EXPS / "strain_identity_qc/outputs/20260630/mono_consensus_20260630.fasta"},
        "20260721": {"pair": EXPS / "20260721/analysis/consensus/strain_consensus_20260721.fasta",
                     "mono": EXPS / "strain_identity_qc/outputs/20260721/mono_consensus_20260721.fasta"},
    }
    L = {e: {k: load(p) for k, p in d.items()} for e, d in builds.items()}

    rows = []
    for exp, d in L.items():
        refs = {"corroborated_db(185)": corr185, "corroborated_db_min5(87)": corr87,
                "genome_16S": gen_w,
                "own_pair_consensus": d["pair"]}          # test C
        for qname, q in d.items():
            for rname, r in refs.items():
                if qname == "pair" and rname == "own_pair_consensus":
                    continue                              # self vs self is trivially rank 1
                rows.append(pd.DataFrame(rank_table(q, r, exp, qname, rname)))
                print(f"  {exp} {qname:5} -> {rname}")
    det = pd.concat(rows, ignore_index=True)
    det.to_csv(OUT / "qc03_reference_ranks.csv", index=False)

    t = det[det.in_reference].copy()
    summ = (t.groupby(["experiment", "consensus", "reference"])
             .agg(n_testable=("label", "size"),
                  n_reference=("n_reference", "max"),
                  median_self_identity=("self_identity", "median"),
                  median_best_identity=("best_identity", "median"),
                  top1=("self_rank", lambda s: 100 * (s == 1).mean()),
                  top3=("self_rank", lambda s: 100 * (s <= 3).mean()),
                  top5=("self_rank", lambda s: 100 * (s <= 5).mean()),
                  median_rank=("self_rank", "median"))
             .round(3).reset_index())
    summ.to_csv(OUT / "qc04_rank_retrieval.csv", index=False)

    sweep = []
    for (e, c, r), s in t.groupby(["experiment", "consensus", "reference"]):
        for th in THRESHOLDS:
            sweep.append({"experiment": e, "consensus": c, "reference": r, "threshold": th,
                          "pct_label_agrees": round(100 * (s.self_identity >= th).mean(), 1),
                          "n_testable": len(s)})
    sw = pd.DataFrame(sweep)
    sw.to_csv(OUT / "qc05_threshold_sweep.csv", index=False)

    print("\n=== B: rank retrieval -- does the reference put the CORRECT label first? ===")
    print(summ.to_string(index=False))
    print("\n=== A: threshold sweep (% label agrees) ===")
    piv = sw.pivot_table(index=["experiment", "consensus", "reference"],
                         columns="threshold", values="pct_label_agrees")
    print(piv.to_string())
    return summ, sw


if __name__ == "__main__":
    main()
