"""01 -- For every unique strain pair actually tested (>5-read pair wells), compute the
distance between the two strains' own reference sequences in merged_consensus_20260721.fasta.

This distance is the key piece of context for everything downstream: if two strains'
references are (near-)identical, no read-mapping method can tell their reads apart, and any
per-well "interaction score" for that pair should be reported with high uncertainty -- not
because the measurement failed, but because the question ("which strain is this read?") has
no answer. See r02 for how this plays out per read.

Output: outputs/r01_reference_pair_distances.csv
    one row per unique tested (strain_a, strain_b) pair: bp and normalized edit distance
    between their reference sequences, or missing_reference=True if either lacks one.
"""

import pandas as pd

from ra_common import OUT_DIR, REFERENCE_FASTA, SAMPLES_GT5READS_CSV
from ra_common import load_reference_db, edit_distance_bp, pair_key


def compute_reference_distances():
    samples = pd.read_csv(SAMPLES_GT5READS_CSV)
    pair_wells = samples[samples["well_type"] == "pair"].copy()

    pairs = sorted({pair_key(r.strain1, r.strain2) for r in pair_wells.itertuples()})

    ref_seqs = load_reference_db(REFERENCE_FASTA)

    rows = []
    for a, b in pairs:
        ref_a, ref_b = ref_seqs.get(a), ref_seqs.get(b)
        if ref_a is None or ref_b is None:
            rows.append(
                {
                    "strain_a": a, "strain_b": b,
                    "bp_dist": None, "norm_dist": None,
                    "len_a": len(ref_a) if ref_a else None,
                    "len_b": len(ref_b) if ref_b else None,
                    "missing_reference": True,
                }
            )
            continue
        bp = edit_distance_bp(ref_a, ref_b)
        rows.append(
            {
                "strain_a": a, "strain_b": b,
                "bp_dist": bp, "norm_dist": bp / max(len(ref_a), len(ref_b)),
                "len_a": len(ref_a), "len_b": len(ref_b),
                "missing_reference": False,
            }
        )

    df = pd.DataFrame(rows)
    out_path = OUT_DIR / "r01_reference_pair_distances.csv"
    df.to_csv(out_path, index=False)

    n_missing = int(df["missing_reference"].sum())
    resolvable = df[~df["missing_reference"]]
    print(f"{len(df)} unique tested strain pairs ({n_missing} missing a reference for one strain)")
    print(
        f"reference-pair bp distance among resolvable pairs: "
        f"min={resolvable.bp_dist.min():.0f}, median={resolvable.bp_dist.median():.0f}, "
        f"max={resolvable.bp_dist.max():.0f}"
    )
    print(f"pairs with identical (0 bp) references: {(resolvable.bp_dist == 0).sum()}")
    print(f"pairs with <10 bp difference: {(resolvable.bp_dist < 10).sum()}")
    print(f"saved -> {out_path}")
    return df


if __name__ == "__main__":
    compute_reference_distances()
