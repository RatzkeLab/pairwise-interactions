"""Stage 10: which strains can go into a community whose composition is to be predicted
from genomes -- and how large can such a community be?

A strain is usable for genome-based modelling only if the genome on file is the organism
in its well. Everything upstream feeds one column, `genomic_use`:

  high      genome verified (own or clone 16S), and no genomically distinct organism on
            file carries the same 16S -- so the readout itself pins down the genome
  medium    genome verified, but other genomes on file carry the same 16S: consistent,
            yet a swap with one of them would be invisible to a 16S readout
  exclude   wrong genome, unverified, no genome, a mixed well (mixed_source or haplotypes
            > 10 bp apart), never grew, or tier D

`genomic_unit` is the clone cluster (or the well itself): wells sharing one are one genome,
so a genome-based model must treat them as one strain, and a community should hold at
most one of them (they are 16S-identical anyway).

Outputs:
  g04_design_eligibility.csv     per strain: everything a design step needs, one row
  g04_max_separable_genomic.csv  largest mutually separable set, genome-usable strains only
"""
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
sys.path.insert(0, str(BASE / "02_distinguishability"))
import config  # noqa: E402

import pandas as pd

from d03_separable_sets import THRESHOLDS, max_separable  # noqa: E402

OUT = config.GENOME_OUT
MAX_HAP_DIVERGENCE = 10


def main():
    ref = pd.read_csv(config.REFDB_OUT / "s01_strain_reference.csv")
    dist = pd.read_csv(config.DIST_OUT / "d01_strain_distinguishability.csv")
    gs = pd.read_csv(OUT / "g02_strain_genome_status.csv")
    pairs = pd.read_csv(config.DIST_OUT / "d01_pairs.csv")

    df = ref[ref.tier != "none"][["strain", "tier", "status", "no_growth", "mixed_source",
                                   "many_het", "haplotype_divergence_bp", "n_haplotypes",
                                   "competing"]]
    df = df.merge(dist[["strain", "identity_group", "nearest", "nearest_bp", "n_identical",
                        "n_within_5", "n_within_15"]], on="strain", how="left")
    df = df.merge(gs[["strain", "genome_status", "evidence", "clone_cluster", "clone_size",
                      "distinct_genomes_with_this_16S", "genome_is_of", "assembly_organisms"]],
                  on="strain", how="left")
    df["genomic_unit"] = df.clone_cluster.fillna(df.strain)
    mixed = df.mixed_source.astype(bool) | (df.haplotype_divergence_bp > MAX_HAP_DIVERGENCE)
    reason = []
    for r, mx in zip(df.itertuples(), mixed):
        if r.tier == "D":
            reason.append("tier_D")
        elif r.no_growth:
            reason.append("no_growth")
        elif mx:
            reason.append("mixed_well")
        elif r.genome_status != "ok":
            reason.append(f"genome_{r.genome_status}")
        else:
            reason.append(None)
    df["exclude_reason"] = reason
    df["genomic_use"] = [
        "exclude" if rs else ("high" if d <= 1 else "medium")
        for rs, d in zip(df.exclude_reason, df.distinct_genomes_with_this_16S.fillna(0))]
    df = df.sort_values(["genomic_use", "identity_group", "strain"])
    df.to_csv(OUT / "g04_design_eligibility.csv", index=False)

    rows = []
    for label, keep in (("high", ["high"]), ("high+medium", ["high", "medium"])):
        names = sorted(df[df.genomic_use.isin(keep)].strain)
        for thr in THRESHOLDS:
            chosen, optimal = max_separable(names, pairs, thr)
            units = df.set_index("strain").genomic_unit.reindex(chosen)
            rows.append(dict(eligible=label, n_eligible=len(names),
                             distinct_genomes_eligible=df[df.strain.isin(names)].genomic_unit.nunique(),
                             min_separation_bp=thr, max_community=len(chosen),
                             distinct_genomes_in_set=units.nunique(), proven_optimal=optimal,
                             example_set=" ".join(chosen)))
    ms = pd.DataFrame(rows)
    ms.to_csv(OUT / "g04_max_separable_genomic.csv", index=False)
    report(df, ms)
    return df, ms


def report(df, ms):
    print("genomic use:", df.genomic_use.value_counts().to_dict())
    print("why excluded:", df.exclude_reason.value_counts().to_dict())
    u = df[df.genomic_use != "exclude"]
    print(f"usable wells {len(u)} = {u.genomic_unit.nunique()} distinct genomes")
    print("\nlargest separable community, genome-usable strains only:")
    print(ms.pivot(index="min_separation_bp", columns="eligible", values="max_community"))


if __name__ == "__main__":
    main()
