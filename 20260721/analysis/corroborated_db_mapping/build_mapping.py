"""A best-effort mapping from 20260721's (scrambled) well labels to corroborated_db entries.

20260721's labels do not name the organisms its wells actually contain -- established three
independent ways (see strain_identity_qc/README.md and genomic_ml_join_test/). The wells do
hold real, clean collection organisms, so the obvious question is: which ones?

This builds that mapping by 16S, against
`collapse_naive_updated3_15diff/corroborated_db_filtered_min5_edited.fasta` (87 entries; note
this file is byte-for-byte identical in content to the un-`_edited` version -- same 87 ids,
zero sequence differences).

**Read the confidence column before using any row.** 16S frequently cannot separate members of
this collection: `strain_identity_qc/qc_recovery.py` already showed that independent best-hit
recovery collapses 75 labels onto 39 targets, one target claimed by 10 labels, with a median
winning margin of ~2.7 bp in 1420. A mapping is offered here anyway because it is the useful
form of the question, not because 16S has become adequate.

Three candidate columns, in increasing order of how much you should trust them:

  `best_hit_*`   -- independent argmax per label. Physically impossible as a plate map (it lets
                    six wells claim one organism), but it is the raw evidence.
  `assigned_*`   -- ONE-TO-ONE assignment over the identity matrix (Hungarian), **constrained**
                    so a label may only be assigned to a target it matches at
                    >= SAME_ORGANISM_IDENTITY, and left unassigned otherwise. The bijection is
                    the right structural prior -- a source plate holds distinct strains in
                    distinct wells -- but it must be constrained: unconstrained, the solver
                    maximises the total and happily pays for it by pushing individual labels
                    onto organisms they match at 0.73, which is worse than no answer.
  `recommended_target` -- the column to actually use. Populated only where the call is
                    defensible; deliberately blank for ~half the plate.

Where `best_hit` and `assigned_one_to_one` disagree, the disagreement itself is the finding:
those labels sit in a cluster of near-identical db entries that 16S cannot split.

Evidence columns let you audit any row: identity and bp distance to the winner, the runner-up
and the margin over it, whether the label's own name exists in the db and how it ranks, and
whether an independent mono-well consensus (26 strains) agrees with the pair-well consensus.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.optimize import linear_sum_assignment

HERE = Path(__file__).resolve().parent
EXPS = HERE.parents[2]
sys.path.insert(0, str(EXPS / "strain_identity_qc"))

from qc_compare import identity, SAME_ORGANISM_IDENTITY   # HW infix, both strands, N-wildcard

TARGET_FASTA = Path("/home/rl/scripts/karl/merge_consensus_sequences/"
                    "collapse_naive_updated3_15diff/corroborated_db_filtered_min5_edited.fasta")
PAIR_FASTA = EXPS / "20260721/analysis/consensus/merged_consensus_20260721.fasta"
MONO_FASTA = EXPS / "20260721/analysis/consensus/mono_strain_consensus_20260721.fasta"
MAPPING_CSV = Path("/home/rl/scripts/karl/Link to Karl/final_genomic_tables/"
                   "mapping_384_well_plate_collection.csv")
OUT = HERE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

# a win by less than this many bp is not a call -- calibrated against mono-well ground truth in
# strain_identity_qc/qc_readassign.py, same constant as relative_abundance.MIN_RESOLVABLE_BP
MIN_RESOLVABLE_BP = 10
WEAK_MARGIN_BP = 3


def load(path, strip=True):
    out = {}
    for r in SeqIO.parse(path, "fasta"):
        name = r.id.split("|")[0]
        if strip and name.endswith("_consensus"):
            name = name[: -len("_consensus")]
        out.setdefault(name, str(r.seq).upper())
    return out


def identity_matrix(queries, targets):
    qn, tn = sorted(queries), sorted(targets)
    I = np.zeros((len(qn), len(tn)))
    D = np.full((len(qn), len(tn)), np.nan)
    for i, q in enumerate(qn):
        for j, t in enumerate(tn):
            I[i, j], D[i, j] = identity(queries[q], targets[t])
    return pd.DataFrame(I, index=qn, columns=tn), pd.DataFrame(D, index=qn, columns=tn)


def build():
    tgt = load(TARGET_FASTA)
    pair = load(PAIR_FASTA)
    mono = load(MONO_FASTA)
    print(f"targets {len(tgt)} | 20260721 pair-consensus {len(pair)} | mono-consensus {len(mono)}")

    I, D = identity_matrix(pair, tgt)
    I.to_csv(OUT / "m01_identity_matrix_pair_vs_corroborated.csv")
    Im, _ = identity_matrix(mono, tgt) if mono else (None, None)

    order = np.argsort(-I.values, axis=1)
    rows = []
    for i, lab in enumerate(I.index):
        b, r = order[i, 0], order[i, 1]
        best, runner = I.columns[b], I.columns[r]
        # margin in bp, the unit the resolution limit is calibrated in
        margin_bp = float(D.iloc[i, r] - D.iloc[i, b])
        self_id = float(I.loc[lab, lab]) if lab in I.columns else np.nan
        self_rank = (int((I.iloc[i] > I.loc[lab, lab]).sum()) + 1) if lab in I.columns else np.nan
        mono_best = np.nan
        if Im is not None and lab in Im.index:
            mono_best = Im.columns[int(np.argmax(Im.loc[lab].values))]
        rows.append({
            "label_20260721": lab,
            "best_hit": best,
            "best_hit_identity": round(float(I.iloc[i, b]), 5),
            "best_hit_bp_dist": float(D.iloc[i, b]),
            "runner_up": runner,
            "runner_up_identity": round(float(I.iloc[i, r]), 5),
            "margin_over_runner_up_bp": margin_bp,
            "label_name_in_db": lab in I.columns,
            "self_identity": round(self_id, 5) if self_id == self_id else np.nan,
            "self_rank": self_rank,
            "mono_best_hit": mono_best,
            "mono_agrees_with_pair": (mono_best == best) if isinstance(mono_best, str) else np.nan,
        })
    df = pd.DataFrame(rows)

    # one-to-one assignment: a source plate cannot put one organism in six wells.
    # CONSTRAINED -- forbidden pairs get a cost worse than any real pairing, and are stripped
    # after solving, so a label with no acceptable partner comes back unassigned instead of
    # being forced onto an unrelated organism to improve the global total.
    allowed = I.values >= SAME_ORGANISM_IDENTITY
    cost = np.where(allowed, -I.values, 1e6)
    ri, ci = linear_sum_assignment(cost)
    assign = {I.index[a]: I.columns[b] for a, b in zip(ri, ci)
              if I.values[a, b] >= SAME_ORGANISM_IDENTITY}
    df["assigned_one_to_one"] = df["label_20260721"].map(assign)
    df["assigned_identity"] = [
        round(float(I.loc[l, a]), 5) if isinstance(a, str) else np.nan
        for l, a in zip(df["label_20260721"], df["assigned_one_to_one"])]
    df["assignment_equals_best_hit"] = df["assigned_one_to_one"] == df["best_hit"]

    # how contested is each target under independent best-hit?
    claims = df["best_hit"].value_counts()
    df["n_labels_claiming_best_hit"] = df["best_hit"].map(claims)

    def tier(r):
        if r["best_hit_identity"] < SAME_ORGANISM_IDENTITY:
            return "no_confident_match"
        if r["margin_over_runner_up_bp"] < WEAK_MARGIN_BP:
            return "low_16S_cannot_separate"
        # an independent mono-well consensus that names a different organism outranks a wide
        # margin: the margin only says the pair-well consensus is decisive, not that it is right
        if r["mono_agrees_with_pair"] is False:
            return "contradicted_by_mono_well"
        if r["margin_over_runner_up_bp"] >= MIN_RESOLVABLE_BP:
            return "high"
        return "medium"
    df["confidence"] = df.apply(tier, axis=1)

    # the usable answer: only where the call stands on its own evidence
    df["recommended_target"] = np.where(
        df["confidence"].isin(["high", "medium"]), df["best_hit"], None)

    # does the recovered name reach the genomic tables?
    m = pd.read_csv(MAPPING_CSV)
    w2s = dict(zip(m["Well_souce_plate"], m["strain"].astype(str)))
    df["best_hit_genome"] = df["best_hit"].map(w2s)
    df["assigned_genome"] = df["assigned_one_to_one"].map(w2s)
    df["recommended_genome"] = df["recommended_target"].map(w2s)

    df = df.sort_values(["confidence", "margin_over_runner_up_bp"],
                        ascending=[True, False]).reset_index(drop=True)
    df.to_csv(OUT / "m02_20260721_to_corroborated_db_mapping.csv", index=False)

    summ = pd.DataFrame([
        {"metric": "labels", "value": len(df)},
        {"metric": "db_entries", "value": len(tgt)},
        {"metric": "label name is its own best hit (i.e. label was right)",
         "value": int((df["best_hit"] == df["label_20260721"]).sum())},
        {"metric": "best hit identity >= 0.99", "value": int((df["best_hit_identity"] >= 0.99).sum())},
        {"metric": "confidence=high", "value": int((df.confidence == "high").sum())},
        {"metric": "confidence=medium", "value": int((df.confidence == "medium").sum())},
        {"metric": "confidence=low_16S_cannot_separate",
         "value": int((df.confidence == "low_16S_cannot_separate").sum())},
        {"metric": "confidence=contradicted_by_mono_well",
         "value": int((df.confidence == "contradicted_by_mono_well").sum())},
        {"metric": "confidence=no_confident_match",
         "value": int((df.confidence == "no_confident_match").sum())},
        {"metric": "distinct targets used by best-hit", "value": int(df["best_hit"].nunique())},
        {"metric": "max labels claiming one target", "value": int(claims.max())},
        {"metric": "one-to-one assignment == best hit",
         "value": int(df["assignment_equals_best_hit"].sum())},
        {"metric": "median margin over runner-up (bp)",
         "value": round(float(df["margin_over_runner_up_bp"].median()), 2)},
        {"metric": "mono consensus available", "value": int(df["mono_best_hit"].notna().sum())},
        {"metric": "mono agrees with pair best-hit",
         "value": int((df["mono_agrees_with_pair"] == True).sum())},
        {"metric": "RECOMMENDED (usable) rows", "value": int(df["recommended_target"].notna().sum())},
        {"metric": "one-to-one assigned (constrained)", "value": int(df["assigned_one_to_one"].notna().sum())},
        {"metric": "recommended name reaches genomic tables",
         "value": int(df["recommended_genome"].notna().sum())},
    ])
    summ.to_csv(OUT / "m03_mapping_summary.csv", index=False)
    print(summ.to_string(index=False))
    print(f"\nwritten -> {OUT}")
    return df, summ


if __name__ == "__main__":
    build()
