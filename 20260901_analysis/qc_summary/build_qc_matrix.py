"""One table behind both QC documents: does each experiment's 16S consensus support its labels?

For every (experiment x consensus build x reference) cell, asks the same question: taking the
consensus sequence the experiment produced for well label L, does the reference agree that L is
that organism?

  label_supported      L exists in the reference AND is the best match at >= 0.99 identity
  label_contradicted   something OTHER than L is the best match at >= 0.99
  no_confident_match   nothing in the reference reaches 0.99
  label_absent         L is not in the reference at all, so the question cannot be asked

Consensus builds compared, because they fail differently:
  pair    built from two-strain wells -- self-consistent by construction, so it is the one
          source that can be confidently wrong
  mono    built from single-strain wells only -- immune to two-strain deconvolution error, and
          built here by the SAME code for both experiments so the comparison is like-for-like
  merged  pair + mono combined (20260721 only; 20260630 has no merged build)

References:
  corroborated_db            185 entries, external, from prior corroborated experiments
  corroborated_db_min5       87 entries, the >=5-support subset of the same
  other_experiment           the other experiment's own pair consensus -- do the two runs agree
  genome_16S                 16S extracted from the NGS assemblies (multi-copy; the weak source,
                             ~65% concordant with corroborated_db even before any experiment)

Identity uses strain_identity_qc.qc_compare.identity: infix alignment of the shorter sequence
into the longer, both strands, N treated as a wildcard.
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
from qc_compare import identity, SAME_ORGANISM_IDENTITY

OUT = HERE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
CORR = KARL / "merge_consensus_sequences" / "collapse_naive_updated3_15diff"


def load(path, multi=False):
    """{label: seq} or {label: [seqs]} for multi-copy sources."""
    out = {}
    for r in SeqIO.parse(path, "fasta"):
        name = r.id.split("|")[0]
        if name.endswith("_consensus"):
            name = name[: -len("_consensus")]
        s = str(r.seq).upper()
        if multi:
            out.setdefault(name, []).append(s)
        else:
            out.setdefault(name, s)
    return out


def best_hit(q, refs):
    best, bid = None, -1.0
    for name, seqs in refs.items():
        for s in (seqs if isinstance(seqs, list) else [seqs]):
            i, _ = identity(q, s)
            if i > bid:
                bid, best = i, name
    return best, bid


def compare(queries, refs, exp, qname, rname):
    rows = []
    for lab, q in queries.items():
        qs = q[0] if isinstance(q, list) else q
        in_ref = lab in refs
        b, bid = best_hit(qs, refs)
        self_id = np.nan
        if in_ref:
            rs = refs[lab]
            self_id = max(identity(qs, s)[0] for s in (rs if isinstance(rs, list) else [rs]))
        if not in_ref:
            verdict = "label_absent"
        elif bid < SAME_ORGANISM_IDENTITY:
            verdict = "no_confident_match"
        elif b == lab or (self_id >= SAME_ORGANISM_IDENTITY and self_id >= bid - 1e-9):
            verdict = "label_supported"
        else:
            verdict = "label_contradicted"
        rows.append({"experiment": exp, "consensus": qname, "reference": rname,
                     "label": lab, "best_match": b, "best_identity": round(bid, 5),
                     "self_identity": round(self_id, 5) if self_id == self_id else np.nan,
                     "verdict": verdict})
    return rows


def main():
    corr185 = load(CORR / "corroborated_db.fasta")
    corr87 = load(CORR / "corroborated_db_filtered_min5.fasta")
    gen = load(KARL / "Link to Karl" / "rDNA_16S_db_all_strains.fasta", multi=True)
    mp = pd.read_csv(KARL / "Link to Karl" / "final_genomic_tables" /
                     "mapping_384_well_plate_collection.csv")
    # genome 16S is keyed by assembly name; re-key to well label so labels are comparable
    a2w = dict(zip(mp["assembly_name"].astype(str), mp["Well_souce_plate"]))
    gen_w = {}
    for k, v in gen.items():
        w = a2w.get(k)
        if w:
            gen_w.setdefault(w, []).extend(v)

    builds = {
        "20260630": {
            "pair": EXPS / "20260630/analysis/consensus2/strain_consensus_20260630.fasta",
            "mono": EXPS / "strain_identity_qc/outputs/20260630/mono_consensus_20260630.fasta",
        },
        "20260721": {
            "pair": EXPS / "20260721/analysis/consensus/strain_consensus_20260721.fasta",
            "mono": EXPS / "strain_identity_qc/outputs/20260721/mono_consensus_20260721.fasta",
            "merged": EXPS / "20260721/analysis/consensus/merged_consensus_20260721.fasta",
        },
    }
    loaded = {e: {k: load(p) for k, p in d.items() if Path(p).exists()}
              for e, d in builds.items()}
    for e, d in loaded.items():
        print(f"{e}: " + ", ".join(f"{k}={len(v)}" for k, v in d.items()))

    rows = []
    for exp, d in loaded.items():
        other = "20260721" if exp == "20260630" else "20260630"
        refs = {"corroborated_db(185)": corr185,
                "corroborated_db_min5(87)": corr87,
                "genome_16S": gen_w,
                f"other_experiment_pair({other})": loaded[other]["pair"]}
        for qname, q in d.items():
            for rname, r in refs.items():
                rows.append(pd.DataFrame(compare(q, r, exp, qname, rname)))
                print(f"  {exp} {qname:7} vs {rname}")
    det = pd.concat(rows, ignore_index=True)
    det.to_csv(OUT / "qc01_consensus_vs_reference_detail.csv", index=False)

    piv = (det.groupby(["experiment", "consensus", "reference", "verdict"]).size()
              .unstack("verdict", fill_value=0).reset_index())
    for c in ("label_supported", "label_contradicted", "no_confident_match", "label_absent"):
        if c not in piv:
            piv[c] = 0
    piv["n_labels"] = piv[["label_supported", "label_contradicted",
                           "no_confident_match", "label_absent"]].sum(axis=1)
    testable = piv.n_labels - piv.label_absent
    piv["pct_supported_of_testable"] = (100 * piv.label_supported / testable.replace(0, np.nan)).round(1)
    piv.to_csv(OUT / "qc02_consensus_vs_reference_summary.csv", index=False)
    print("\n=== label support by consensus build and reference ===")
    print(piv[["experiment", "consensus", "reference", "n_labels", "label_supported",
               "label_contradicted", "no_confident_match", "label_absent",
               "pct_supported_of_testable"]].to_string(index=False))
    return piv


if __name__ == "__main__":
    main()
