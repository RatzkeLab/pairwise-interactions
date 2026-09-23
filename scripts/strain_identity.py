"""Does the strain in well X actually match what every other source says is in well X?

This is the question 20260721 failed: its consensus sequences were self-consistent
but only 3.7% of well labels were supported by the collection references, against
96% for 20260630 -- a source-plate mix-up, proven twice.

The test has to be run carefully, because a naive "best hit has the wrong name"
over-reports badly. 234 of the 371 strains in the full collection share an identical
16S with at least one other strain, so for those a wrong-named best hit is not
evidence of anything. Every strain is therefore placed in one of four classes:

  supported     best hit is the strain's own label
  ambiguous     best hit is a different label, but that label's reference sequence
                is within the ONT noise floor of this one -- 16S cannot separate them
  contradicted  best hit is a different label that IS separable. Real disagreement.
  absent        the strain has no entry in this reference

Only `contradicted` counts against a run, and only the supported/contradicted split
is comparable between runs and references.

When a run is contradicted wholesale, `permutation_structure` asks the follow-up:
is the mismatch random, or does it map well X onto well f(X) systematically? A
row shift, a column shift, a quadrant swap or a plate rotation all leave a
signature, and knowing which one it is turns an unusable plate into a fixable one.
"""
import collections
import numpy as np
import pandas as pd
import edlib

from io_utils import load_reference_db

NOISE_FLOOR = 0.02        # normalized edit distance; below this two references are
                          # indistinguishable given ONT consensus accuracy. Used for the
                          # twin sets in `compare`, where the question is whether a
                          # best hit could have gone either way.
RESOLVABLE_BP = 10        # absolute edit distance, matching relative_abundance's
                          # MIN_RESOLVABLE_BP: two 16S sequences this far apart CAN be
                          # told apart by read assignment. 0.02 normalized is ~28 bp and
                          # is far too coarse to use as a resolution limit -- the
                          # experiment was designed around >=15 bp separation.
ROWS = "ABCDEFGHIJKLMNOP"


def load_multicopy_db(path):
    """name -> [sequences]. The genome-derived rDNA db carries one entry per rRNA
    operon, so a genome appears several times under the same assembly name; keying
    by name and overwriting would silently drop copies that differ."""
    from Bio import SeqIO
    out = collections.defaultdict(list)
    for rec in SeqIO.parse(path, "fasta"):
        out[rec.id].append(str(rec.seq).upper())
    return dict(out)


def collapse_multicopy(db, query=None):
    """Pick one sequence per name: the copy closest to `query` if given (a genome
    matches if ANY of its operons matches), else the longest."""
    if query is None:
        return {n: max(v, key=len) for n, v in db.items()}
    out = {}
    for n, v in db.items():
        out[n] = min(v, key=lambda s: edlib.align(query, s, mode="NW",
                                                  task="distance")["editDistance"])
    return out


def flexible_distance(a, b):
    """Divergence over the overlapping region, not over the longer sequence.

    Half the entries in the genome-derived rDNA db are partial 16S -- assembly
    fragments as short as 396 bp against our ~1420 bp consensus. Global (NW)
    alignment scores those at ~0.65 normalized distance no matter how well the
    shared region actually matches, which makes every fragment look like a
    different organism. Aligning the shorter sequence as an infix of the longer
    and normalizing by the shorter length asks the right question: over the part
    both sequences cover, do they agree?

    Returns (normalized_distance, overlap_bp). A short overlap is weak evidence
    either way, so the caller keeps the length and can require a minimum.
    """
    if len(a) <= len(b):
        short, long_ = a, b
    else:
        short, long_ = b, a
    ed = edlib.align(short, long_, mode="HW", task="distance")["editDistance"]
    return ed / len(short), len(short)


def _best_hits(query_seq, ref, k=3, min_overlap=0):
    d = []
    for name, seq in ref.items():
        dist, ov = flexible_distance(query_seq, seq)
        if ov < min_overlap:
            continue
        d.append((dist, name, ov))
    d.sort()
    return d[:k]


def _twin_sets(ref, floor=NOISE_FLOOR):
    """name -> set of reference names indistinguishable from it."""
    names = list(ref)
    twins = {n: {n} for n in names}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            dist, _ = flexible_distance(ref[a], ref[b])
            if dist <= floor:
                twins[a].add(b); twins[b].add(a)
    return twins


def compare(consensus, reference_fasta, ref_name, name_map=None, floor=NOISE_FLOOR,
            twins=None, min_overlap=0):
    """consensus: {strain: seq} from this run. name_map: our label -> reference label."""
    ref = load_reference_db(reference_fasta) if not isinstance(reference_fasta, dict) \
        else reference_fasta
    if twins is None:
        twins = _twin_sets(ref, floor)
    rows = []
    for strain, seq in consensus.items():
        expect = name_map.get(strain) if name_map else strain
        if expect not in ref:
            rows.append(dict(strain=strain, reference=ref_name, expected=expect,
                             verdict="absent"))
            continue
        hits = _best_hits(seq, ref, min_overlap=min_overlap)
        if not hits:
            rows.append(dict(strain=strain, reference=ref_name, expected=expect,
                             verdict="absent"))
            continue
        best_d, best_n, best_ov = hits[0]
        own_d, own_ov = flexible_distance(seq, ref[expect])
        if best_n == expect:
            verdict = "supported"
        elif best_n in twins.get(expect, {expect}):
            verdict = "ambiguous"
        else:
            verdict = "contradicted"
        rows.append(dict(strain=strain, reference=ref_name, expected=expect,
                         best_hit=best_n, best_dist=round(best_d, 4),
                         best_overlap_bp=best_ov,
                         own_dist=round(own_d, 4), own_overlap_bp=own_ov,
                         margin=round(own_d - best_d, 4),
                         n_twins=len(twins.get(expect, {expect})) - 1,
                         verdict=verdict))
    return pd.DataFrame(rows)


def summarize(df, label=""):
    """The comparable number: supported / (supported + contradicted)."""
    v = df.verdict.value_counts()
    s, c = v.get("supported", 0), v.get("contradicted", 0)
    amb, ab = v.get("ambiguous", 0), v.get("absent", 0)
    rate = s / (s + c) if (s + c) else float("nan")
    print(f"  {label:22s} supported {s:4d}  contradicted {c:4d}  "
          f"ambiguous {amb:4d}  absent {ab:4d}   ->  label support {rate:6.1%}")
    return dict(reference=label, supported=s, contradicted=c, ambiguous=amb,
                absent=ab, label_support=rate, n_decidable=s + c)


def permutation_structure(df):
    """If labels are wrong, are they wrong in a pattern?

    Takes the contradicted rows and looks at the offset between the well we think a
    strain is in and the well its sequence actually matches. A consistent (drow, dcol)
    means the plate was shifted; a consistent row-only or column-only offset narrows
    it further; a flat distribution means the labels are genuinely scrambled and no
    single transform recovers them.
    """
    bad = df[df.verdict == "contradicted"].dropna(subset=["best_hit"])
    if bad.empty:
        print("  no contradicted strains -- nothing to explain")
        return None
    offs = collections.Counter()
    for r in bad.itertuples():
        try:
            er, ec = ROWS.index(r.expected[0]), int(r.expected[1:])
            br, bc = ROWS.index(r.best_hit[0]), int(r.best_hit[1:])
        except (ValueError, IndexError):
            continue
        offs[(br - er, bc - ec)] += 1
    if not offs:
        print("  best hits are not well-coordinate names; no geometry to test")
        return None
    n = sum(offs.values())
    top = offs.most_common(5)
    print(f"  offset (row, col) between claimed well and matching well, "
          f"{n} contradicted strains:")
    for (dr, dc), c in top:
        print(f"     {dr:+3d},{dc:+3d}   {c:4d}  ({c/n:5.1%})")
    frac = top[0][1] / n
    if frac > 0.3:
        print(f"  -> {frac:.0%} share one offset: this is a systematic shift, "
              f"not a scramble. Applying it would recover the labels.")
    else:
        print(f"  -> the most common offset covers only {frac:.0%}; no single shift "
              f"explains it. The labels are scrambled well-by-well.")
    return pd.DataFrame([dict(d_row=dr, d_col=dc, n=c, frac=c / n)
                         for (dr, dc), c in offs.most_common()])


def self_consistency(consensus, reference=None, floor_bp=RESOLVABLE_BP):
    """Is the run internally coherent, independent of any external reference?

    Reports how many consensus sequences are distinguishable from each other. The
    number is only interpretable against a baseline, because 16S genuinely cannot
    separate much of this collection: `reference` (a {name: seq} dict, normally
    full_collection_db) is scored the same way so the two can be compared. A run
    resolving about as many strains as the reference does is doing as well as 16S
    allows -- it is not a defect of the run.

    A run can be perfectly self-consistent and still have every label wrong, which
    is exactly the 20260721 situation, so this never merges with `compare`.
    """
    def score(seqs):
        names = list(seqs)
        twinned = set()
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                dist, ov = flexible_distance(seqs[a], seqs[b])
                if dist * ov < floor_bp:
                    twinned.add(a); twinned.add(b)
        return len(names) - len(twinned), len(names)

    uniq, n = score(consensus)
    print(f"  this run : {uniq}/{n} consensus sequences are separable at "
          f"{floor_bp} bp ({uniq/n:.1%})")
    if reference:
        runiq, rn = score(reference)
        print(f"  reference: {runiq}/{rn} ({runiq/rn:.1%}) -- the ceiling 16S allows "
              f"for this collection")
        ratio = (uniq / n) / (runiq / rn)
        print(f"  -> the run reaches {ratio:.0%} of that ceiling. Most of what it "
              f"cannot separate, nothing can: 16S is simply degenerate across this "
              f"collection. Only a shortfall well below ~85% would point at "
              f"consensus quality rather than at 16S.")
    return uniq, n
