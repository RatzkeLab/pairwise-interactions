"""Compare 16S sequences across sources, and ask what each sequence actually matches.

Two different questions, and keeping them apart is the whole design:

  AGREEMENT  -- for strain S, does source X's sequence look like source Y's sequence?
                Answers "is this strain's identity consistent?"
  ATTRIBUTION -- for strain S in source X, which strain in source Y does it match BEST, out of
                the entire collection? Answers "if it isn't S, what is it?" -- the question
                that separates a label mix-up (it is confidently some OTHER specific strain)
                from contamination or a failed consensus (it is nothing in particular).

Three things every comparison here must handle, each of which silently corrupts the result if
skipped:
  - **Strand.** Genome-extracted 16S comes off either strand; amplicon consensus is oriented.
    Every comparison tries both and keeps the better.
  - **Length.** Genome 16S records run 396-1563 bp (fragmented contig ends) vs ~1420 bp
    amplicons, so raw edit distance is not comparable across pairs. The shorter sequence is
    aligned into the longer (edlib 'HW', free end gaps) and the metric is per-base identity.
  - **Ns.** Low-coverage consensus carries internal Ns; they are treated as wildcards rather
    than mismatches, so a thin well is not mistaken for a different organism.
"""

import numpy as np
import pandas as pd
import edlib

# N matches anything, both directions, so an ambiguous base is never counted as a mismatch
_WILDCARD = ([(a, "N") for a in "ACGT"] + [("N", a) for a in "ACGT"]
             + [(a, "?") for a in "ACGTN"] + [("?", a) for a in "ACGTN"])
_COMP = str.maketrans("ACGTN", "TGCAN")

SAME_ORGANISM_IDENTITY = 0.99   # ~14 bp over a 1420 bp amplicon; ONT consensus error sits well
                                # below this, real congeners sit well above


def revcomp(s):
    return s.translate(_COMP)[::-1]


def identity(a, b):
    """Best per-base identity of the shorter sequence within the longer, over both strands."""
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    best = None
    for cand in (q, revcomp(q)):
        d = edlib.align(cand, t, mode="HW", task="distance",
                        additionalEqualities=_WILDCARD)["editDistance"]
        if d >= 0 and (best is None or d < best):
            best = d
    if best is None:
        return 0.0, np.nan
    return 1.0 - best / len(q), best


def _rep_seqs(df, source):
    """{strain: [seqs]} for one source. Lists, because genome 16S is multi-copy."""
    sub = df[df["source"] == source]
    out = {}
    for r in sub.itertuples():
        out.setdefault(r.strain_label, []).append(r.seq)
    return out


def best_identity(seqs_a, seqs_b):
    """Best identity over every copy-vs-copy combination (multi-copy 16S)."""
    best, bd = -1.0, np.nan
    for x in seqs_a:
        for y in seqs_b:
            i, d = identity(x, y)
            if i > best:
                best, bd = i, d
    return best, bd


def agreement_table(df, sources=None):
    """Per strain, per ordered source pair: do the two sources agree on this strain?"""
    sources = sources or sorted(df["source"].unique())
    reps = {s: _rep_seqs(df, s) for s in sources}
    rows = []
    for i, sa in enumerate(sources):
        for sb in sources[i + 1:]:
            shared = sorted(set(reps[sa]) & set(reps[sb]))
            for st in shared:
                ident, dist = best_identity(reps[sa][st], reps[sb][st])
                rows.append({"strain_label": st, "source_a": sa, "source_b": sb,
                             "identity": ident, "edit_distance": dist,
                             "same_organism": ident >= SAME_ORGANISM_IDENTITY})
    return pd.DataFrame(rows)


def attribution_table(df, query_source, target_source, restrict_to=None):
    """For each strain in query_source, which target_source strain does it match best?

    Searches the ENTIRE target collection, not just the same label -- that is what turns
    "these disagree" into "this well actually contains that other strain".
    """
    q, t = _rep_seqs(df, query_source), _rep_seqs(df, target_source)
    if restrict_to is not None:
        q = {k: v for k, v in q.items() if k in set(restrict_to)}
    targets = sorted(t)
    rows = []
    for st, qs in sorted(q.items()):
        scores = np.array([best_identity(qs, t[tt])[0] for tt in targets])
        order = np.argsort(-scores)
        best_i = order[0]
        self_i = targets.index(st) if st in t else None
        second = scores[order[1]] if len(order) > 1 else np.nan
        rows.append({
            "strain_label": st,
            "best_match": targets[best_i],
            "best_identity": scores[best_i],
            "self_identity": scores[self_i] if self_i is not None else np.nan,
            "self_rank": int(np.where(order == self_i)[0][0]) + 1 if self_i is not None else np.nan,
            "runner_up_identity": second,
            "margin_over_runner_up": scores[best_i] - second,
            "self_is_best": (targets[best_i] == st),
            "query_source": query_source, "target_source": target_source,
        })
    return pd.DataFrame(rows)


def summarize_attribution(att):
    """One line per source pair, separating three genuinely different outcomes.

    `self_is_best` alone is a trap: this collection contains near-identical strains (the same
    ones that force ~50/50 in the interaction assay), so a correct label routinely loses "best
    match" to a congener by a base or two. That is an ambiguity, not an error. Hence:

      label_consistent  -- the label's own sequence matches at organism level. The label is
                           supported, whether or not it also wins outright.
      label_unique      -- ...and no other strain matches as well. Fully resolved.
      label_wrong       -- the label's own sequence does NOT match, but some other strain
                           does. This is the one that means a mix-up.
      unresolved        -- nothing matches at organism level; bad consensus or a strain
                           missing from the reference.
    """
    rows = []
    for (qs, ts), g in att.groupby(["query_source", "target_source"]):
        # a label the reference has never heard of is NOT evidence of a mix-up -- it is simply
        # untestable, and every rate below is computed over the testable labels only
        testable = g[g["self_identity"].notna()]
        consistent = testable["self_identity"] >= SAME_ORGANISM_IDENTITY
        unique = consistent & testable["self_is_best"]
        wrong = (~consistent) & (testable["best_identity"] >= SAME_ORGANISM_IDENTITY)
        unresolved = (~consistent) & (~wrong)
        n = max(len(testable), 1)
        rows.append({
            "query_source": qs, "target_source": ts, "n_labels": len(g),
            "n_testable": len(testable),
            "n_label_not_in_reference": int(len(g) - len(testable)),
            "n_label_consistent": int(consistent.sum()),
            "n_label_unique": int(unique.sum()),
            "n_label_wrong": int(wrong.sum()),
            "n_unresolved": int(unresolved.sum()),
            "pct_label_consistent": round(100 * consistent.sum() / n, 1),
            "pct_label_wrong": round(100 * wrong.sum() / n, 1),
            "median_self_identity": round(float(testable["self_identity"].median()), 4),
        })
    return pd.DataFrame(rows)
