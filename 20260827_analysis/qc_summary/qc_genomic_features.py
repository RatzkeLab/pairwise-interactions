"""QC centred on the GENOMIC side: is the genome attached to each well the right organism?

Everything here compares sequences that are keyed to a well label via
`mapping_384_well_plate_collection.csv`, so a disagreement means the mapping, the genome
assembly, or the reference is wrong -- not the experiment.

Sections, each writing its own CSV and recording a run timestamp in `gq00_provenance.csv`:

  1  name overlap        how many corroborated_db entries, how many annotated genomes, and how
                         many labels the two could possibly share given the plate mapping
  2  genome_16S vs corroborated_db   the two REFERENCES against each other, no experiment
                         involved. Identity for same-named strains, plus rank retrieval: does a
                         strain's own genome rank first among all genomes? Worst offenders named.
  3  genome_16S vs each experiment's consensus   the same statistics per experiment

Rank retrieval matters more than identity here. Identity is confounded by platform (NGS assembly
vs ONT amplicon), but rank is not: if a strain's own genome ranks 40th of 233, the mapping cannot
be verified for that strain no matter how the identity threshold is set.
"""

import sys
from datetime import datetime
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
OUT.mkdir(parents=True, exist_ok=True)
CORR = KARL / "merge_consensus_sequences" / "collapse_naive_updated3_15diff"
GEN = KARL / "Link to Karl" / "final_genomic_tables"
STAMP = datetime.now().strftime("%Y-%m-%d %H:%M")
PROV = []


def prov(section, func, out_files, note=""):
    PROV.append({"section": section, "generated_by": f"qc_genomic_features.py::{func}",
                 "run_at": STAMP, "outputs": "; ".join(out_files), "note": note})


def load(path, multi=False):
    out = {}
    for r in SeqIO.parse(path, "fasta"):
        n = r.id.split("|")[0]
        if n.endswith("_consensus"):
            n = n[: -len("_consensus")]
        s = str(r.seq).upper()
        out.setdefault(n, []).append(s) if multi else out.setdefault(n, s)
    return out


def genome_16s_by_well():
    gen = load(KARL / "Link to Karl" / "rDNA_16S_db_all_strains.fasta", multi=True)
    mp = pd.read_csv(GEN / "mapping_384_well_plate_collection.csv")
    a2w = dict(zip(mp["assembly_name"].astype(str), mp["Well_souce_plate"]))
    out = {}
    for k, v in gen.items():
        w = a2w.get(k)
        if w:
            out.setdefault(w, []).extend(v)
    return out, mp, gen


def section1(corr185, corr87, gen_w, mp, gen_raw):
    ko = pd.read_csv(GEN / "KEGG_ko_and_strains_table.csv", index_col=0)
    ko.index = ko.index.astype(str)
    annotated = set(ko.index[ko.sum(axis=1) > 0])
    mp_wells = set(mp["Well_souce_plate"])
    ann_wells = {w for w, s in zip(mp["Well_souce_plate"], mp["strain"].astype(str))
                 if s in annotated}
    rows = [
        ("corroborated_db entries", len(corr185)),
        ("corroborated_db_min5 entries", len(corr87)),
        ("wells in the plate mapping", len(mp_wells)),
        ("wells whose genome has non-empty KO annotation", len(ann_wells)),
        ("genome 16S records (multi-copy)", sum(len(v) for v in gen_raw.values())),
        ("assemblies with a 16S record", len(gen_raw)),
        ("wells with a genome 16S after mapping", len(gen_w)),
        ("IDEAL overlap: corroborated_db AND genome 16S", len(set(corr185) & set(gen_w))),
        ("IDEAL overlap: corroborated_db_min5 AND genome 16S", len(set(corr87) & set(gen_w))),
        ("corroborated_db labels NOT in the plate mapping", len(set(corr185) - mp_wells)),
    ]
    df = pd.DataFrame(rows, columns=["quantity", "n"])
    df.to_csv(OUT / "gq01_name_overlap.csv", index=False)
    prov("1 name overlap", "section1", ["gq01_name_overlap.csv"])
    print(df.to_string(index=False))
    return df


def rank_compare(queries, refs, qname, rname):
    names = list(refs)
    rows = []
    for lab, q in queries.items():
        qs = q[0] if isinstance(q, list) else q
        ids = np.array([max(identity(qs, s)[0]
                            for s in (refs[n] if isinstance(refs[n], list) else [refs[n]]))
                        for n in names])
        order = np.argsort(-ids)
        ranked = [names[i] for i in order]
        r = {"query_set": qname, "reference_set": rname, "label": lab,
             "best_match": ranked[0], "best_identity": round(float(ids.max()), 5),
             "n_reference": len(names)}
        if lab in refs:
            r["self_identity"] = round(float(ids[names.index(lab)]), 5)
            r["self_rank"] = int(ranked.index(lab)) + 1
        else:
            r["self_identity"] = np.nan
            r["self_rank"] = np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def summarise(df):
    t = df[df.self_rank.notna()]
    if not len(t):
        return {}
    return {"n_shared_labels": len(t),
            "n_reference": int(df.n_reference.max()),
            "median_identity": round(float(t.self_identity.median()), 4),
            "pct_identity_ge_0.99": round(100 * (t.self_identity >= 0.99).mean(), 1),
            "pct_identity_ge_0.97": round(100 * (t.self_identity >= 0.97).mean(), 1),
            "pct_identity_ge_0.95": round(100 * (t.self_identity >= 0.95).mean(), 1),
            "pct_identity_ge_0.90": round(100 * (t.self_identity >= 0.90).mean(), 1),
            "top1": round(100 * (t.self_rank == 1).mean(), 1),
            "top5": round(100 * (t.self_rank <= 5).mean(), 1),
            "top10": round(100 * (t.self_rank <= 10).mean(), 1),
            "median_rank": float(t.self_rank.median())}


def main():
    corr185 = load(CORR / "corroborated_db.fasta")
    corr87 = load(CORR / "corroborated_db_filtered_min5.fasta")
    gen_w, mp, gen_raw = genome_16s_by_well()
    print("=== 1. name overlap ===")
    section1(corr185, corr87, gen_w, mp, gen_raw)

    print("\n=== 2/3. rank retrieval against the genomes ===")
    queries = {
        "corroborated_db(185)": corr185,
        "20260630_pair_consensus": load(EXPS / "20260630/analysis/consensus2/strain_consensus_20260630.fasta"),
        "20260721_pair_consensus": load(EXPS / "20260721/analysis/consensus/strain_consensus_20260721.fasta"),
    }
    details, summ = [], []
    for qname, q in queries.items():
        d = rank_compare(q, gen_w, qname, "genome_16S_by_well")
        details.append(d)
        s = summarise(d)
        s.update({"query_set": qname, "reference_set": "genome_16S_by_well"})
        summ.append(s)
        print(f"  {qname:26} shared={s.get('n_shared_labels')} "
              f"median_id={s.get('median_identity')} top1={s.get('top1')}% "
              f"top10={s.get('top10')}% median_rank={s.get('median_rank')}")
    det = pd.concat(details, ignore_index=True)
    det.to_csv(OUT / "gq02_vs_genome_16S_detail.csv", index=False)
    sdf = pd.DataFrame(summ)[["query_set", "reference_set", "n_shared_labels", "n_reference",
                              "median_identity", "pct_identity_ge_0.99", "pct_identity_ge_0.97",
                              "pct_identity_ge_0.95", "pct_identity_ge_0.90",
                              "top1", "top5", "top10", "median_rank"]]
    sdf.to_csv(OUT / "gq03_vs_genome_16S_summary.csv", index=False)
    prov("2/3 genome comparisons", "rank_compare/summarise",
         ["gq02_vs_genome_16S_detail.csv", "gq03_vs_genome_16S_summary.csv"])
    print()
    print(sdf.to_string(index=False))

    worst = (det[det.self_rank.notna()]
             .sort_values(["query_set", "self_rank"], ascending=[True, False])
             .groupby("query_set").head(8)
             [["query_set", "label", "self_identity", "self_rank", "n_reference",
               "best_match", "best_identity"]])
    worst.to_csv(OUT / "gq04_worst_strains.csv", index=False)
    prov("2/3 worst strains", "main", ["gq04_worst_strains.csv"],
         "8 worst-ranked labels per query set")
    print("\n=== worst-ranked labels (their own genome is far down the list) ===")
    print(worst.to_string(index=False))

    pd.DataFrame(PROV).to_csv(OUT / "gq00_provenance.csv", index=False)
    print(f"\nprovenance written; all stats generated {STAMP}")
    return sdf


if __name__ == "__main__":
    main()
