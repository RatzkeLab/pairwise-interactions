"""Run the whole 20260907 analysis on the corrected demultiplexing, in order.

    python run_all.py              # the corrected demux
    python run_all.py --original   # the original, for the before/after comparison

Order matters: well composition produces the per-well representative sequences
that the consensus builder consumes, the consensus becomes the reference that
relative abundance scores reads against, and the identity check needs the
consensus to exist. Each step writes CSVs, so any one can be re-run alone once
the steps above it have run once.
"""
import runpy
import sys, time
from pathlib import Path
BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE)); sys.path.insert(0, str(BASE.parents[1] / "scripts"))

import config


def banner(n, title):
    print("\n" + "#" * 78)
    print(f"# {n}. {title}")
    print("#" * 78)


def main(tag="corrected"):
    demux = config.DEMUX_DIR if tag == "corrected" else config.DEMUX_DIR_ORIGINAL
    if not demux.exists():
        sys.exit(f"demux directory not found: {demux}\n"
                 f"(the corrected re-run writes here; see 02_demux_qc/README.md)")
    t0 = time.time()

    banner(1, "demultiplexing QC: why wells are empty, and did they have cells")
    sys.path.insert(0, str(BASE / "02_demux_qc"))
    import dropout_analysis, od_vs_dropout
    dropout_analysis.main(dropout_analysis.DEMUX_ORIGINAL if tag == "original"
                          else dropout_analysis.DEMUX_CORRECTED)
    # barcode_forensics is a script, and it is what adds the `hits_bad_bc` column that
    # od_vs_dropout needs, so it has to run between the two
    runpy.run_path(str(BASE / "02_demux_qc" / "barcode_forensics.py"), run_name="__main__")
    od_vs_dropout.main()

    banner(2, "well composition: one sequence per mono well, two per pair well?")
    sys.path.insert(0, str(BASE / "04_qc"))
    import run_composition
    run_composition.main(demux, tag=tag)

    banner(3, "per-strain consensus, anchored on the monocultures")
    sys.path.insert(0, str(BASE / "03_recreate_reference_db"))
    import run_consensus
    run_consensus.main(tag)

    banner(4, "strain identity: is this run's plate labelled correctly?")
    sys.path.insert(0, str(BASE / "05_cross_reference"))
    import run_identity
    run_identity.main(tag)

    banner(5, "relative abundance and the inoculum arms")
    import mapping_validation as mv, relative_abundance as ra
    cfg = config.make_config(demux)
    cfg.ra_reference_fasta = (BASE / "03_recreate_reference_db" / "outputs" /
                              f"c02_strain_consensus_{tag}.fasta")
    mv.gather_samples(cfg)
    mv.prepare_references(cfg)
    two_strain = ("pair", "techrep", "ratio", "density")
    ra.reference_distances(cfg, well_types=two_strain)
    ra.compute_interaction_scores(cfg, well_types=two_strain)

    sys.path.insert(0, str(BASE / "06_inoculum_arms"))
    import od_only_arms, ratio_titration
    od_only_arms.main()
    ratio_titration.main()

    print(f"\nall steps finished in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main("original" if "--original" in sys.argv else "corrected")
