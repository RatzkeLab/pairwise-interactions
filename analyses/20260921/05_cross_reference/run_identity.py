"""Are 20260907's strains who the labels say they are?

The comparison that matters is against the genome-derived 16S, because that is the
only reference carrying through to the genomic feature tables. The chain is

    our well label  ->  Well_souce_plate  ->  assembly_name  ->  rDNA_16S_db
                                          ->  strain id      ->  KEGG/CAZy/panX/...

If the first link fails, every genomic feature is joined to the wrong organism.
That is what happened to 20260721.

A low score against the genome db is NOT by itself evidence against this run,
because the genome db may simply be a poor reference -- it is assembly-derived,
multi-copy, and built years earlier. So 20260630 (whose labels are known good, 96%
supported by the collection) and 20260721 (known scrambled, 3.7%) are scored
against exactly the same reference in the same pass. Those two are the ruler:

    if 20260630 also scores low, the genome db is the problem, not the run
    if 20260630 scores high and 20260907 does not, the run has a labelling problem

Everything is converted into well-coordinate space first, so a failure can be
tested for geometry (a plate shift or rotation) rather than just counted.

The three references are NOT equally trustworthy, and the project already settled
the ranking (see mapping_validation.py):

  full_collection  PRIMARY. One entry per collection well, same namespace as the
                   layout. The number to quote.
  corroborated     SECONDARY. Only ~87 entries, and its strain names are pooled
                   plate-well labels from unrelated past experiments, so a name
                   match is not by itself proof of the same organism. It scores
                   high partly because it only contains strains that were already
                   corroborated. Good corroboration, bad headline.
  genome_16S       The one that matters for the genomic tables, and the weakest
                   as a sequence reference: assembly-derived, half of it partial.
"""
import sys, pickle
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE)); sys.path.insert(0, str(BASE.parent / "shared_scripts"))

import pandas as pd
import config, strain_identity as si, strain_consensus
from io_utils import load_reference_db

OUT = BASE / "05_cross_reference" / "outputs"; OUT.mkdir(parents=True, exist_ok=True)

# half the genome-derived rDNA entries are assembly fragments; below this the shared
# region is too short for a match or a mismatch to mean anything either way
MIN_OVERLAP_BP = 800
MAPPING = Path("/home/rl/data/resources/karl/mapping_384_well_plate_collection.csv")
PRIOR = BASE.parent
# Each prior run's main cross-well consensus, so the two controls are comparable in
# size (84 and 81 strains). 20260721's mono-only file has just 26 entries and makes a
# much weaker negative control.
PRIOR_CONSENSUS = {
    "20260630": PRIOR / "20260630/03_create_reference_db/consensus2/strain_consensus_20260630.fasta",
    "20260721": PRIOR / "20260721/03_recreate_reference_db/consensus/strain_consensus_20260721.fasta",
}


def genome_db_in_well_space(cfg):
    """rDNA db keyed by assembly_name -> keyed by source-plate well coordinate.

    Also returns the wells with no genomic entry at all: those can never be scored,
    and counting them as failures would understate every run equally.
    """
    m = pd.read_csv(MAPPING)
    asm2well = dict(zip(m["assembly_name"], m["Well_souce_plate"]))
    multi = si.load_multicopy_db(cfg.reference_dbs["rDNA_all_strains"])
    flat = si.collapse_multicopy(multi)
    well_db = {asm2well[a]: s for a, s in flat.items() if a in asm2well}
    print(f"genome 16S db: {len(multi)} assemblies -> {len(well_db)} of 384 wells "
          f"have a genome-derived 16S")
    return well_db


def score_run(name, consensus, references, twins_cache):
    out = []
    for ref_name, ref in references.items():
        if ref_name not in twins_cache:
            twins_cache[ref_name] = si._twin_sets(ref)
        df = si.compare(consensus, ref, ref_name, twins=twins_cache[ref_name],
                        min_overlap=MIN_OVERLAP_BP)
        df["run"] = name
        s = si.summarize(df, f"{name} vs {ref_name}")
        s["run"], s["reference"] = name, ref_name
        out.append((s, df))
    return out


def main(tag="corrected"):
    cfg = config.make_config()
    cdir = BASE / "03_recreate_reference_db" / "outputs"
    consensus = pickle.load(open(cdir / f"c03_consensus_seqs_{tag}.pkl", "rb"))
    summary = pd.read_csv(cdir / f"c01_strain_consensus_summary_{tag}.csv")
    good = set(summary[summary.status.isin(strain_consensus.GOOD)].strain)
    runs = {"20260907": {s: q for s, q in consensus.items() if s in good}}
    for k, p in PRIOR_CONSENSUS.items():
        if p.exists():
            runs[k] = load_reference_db(p)
    print(f"consensus sets: " + ", ".join(f"{k} n={len(v)}" for k, v in runs.items()))

    full = load_reference_db(cfg.reference_dbs["full_collection"])
    print("\nself-consistency of 20260907 (does it resolve its own strains?)")
    si.self_consistency(runs["20260907"], full)

    references = {                       # order = how much weight to give them
        "full_collection": full,                        # primary
        "genome_16S": genome_db_in_well_space(cfg),     # the genomic-table link
        "corroborated": load_reference_db(cfg.reference_dbs["corroborated"]),  # secondary
    }

    print("\n" + "=" * 78)
    print("LABEL SUPPORT: does each run's well X match reference well X?")
    print("  20260630 is the positive control, 20260721 the negative control.")
    print("=" * 78)
    stats, frames, twins = [], [], {}
    for run_name, cons in runs.items():
        print(f"\n{run_name}:")
        for s, df in score_run(run_name, cons, references, twins):
            stats.append(s); frames.append(df)

    stat_df = pd.DataFrame(stats)
    piv = stat_df.pivot(index="run", columns="reference", values="label_support")
    print("\n" + "=" * 78)
    print("label support (supported / (supported + contradicted))")
    print("=" * 78)
    print((piv * 100).round(1).to_string())

    print("\ngeometry of the disagreement against the genome 16S:")
    for run_name in runs:
        d = pd.concat(frames)
        d = d[(d.run == run_name) & (d.reference == "genome_16S")]
        print(f"\n  {run_name}:")
        si.permutation_structure(d)

    stat_df.to_csv(OUT / f"x01_label_support_{tag}.csv", index=False)
    pd.concat(frames).to_csv(OUT / f"x02_per_strain_verdicts_{tag}.csv", index=False)
    print(f"\nwrote -> {OUT}")
    return stat_df, frames


if __name__ == "__main__":
    main("original" if "--original" in sys.argv else "corrected")
