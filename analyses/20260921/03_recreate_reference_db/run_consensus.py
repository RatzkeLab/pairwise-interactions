"""Build a 16S consensus per strain from this run's wells.

Depends on 04_qc/run_composition.py having run first -- it produces the per-well
representative sequences this consumes, and doing the clustering once for both is
what makes the pass affordable.
"""
import sys, pickle
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE)); sys.path.insert(0, str(BASE.parents[1] / "scripts"))

import pandas as pd
import config, strain_consensus, strain_identity

OUT = BASE / "03_recreate_reference_db" / "outputs"; OUT.mkdir(parents=True, exist_ok=True)
QC = BASE / "04_qc" / "outputs"


def main(tag="corrected"):
    cfg = config.make_config(config.DEMUX_DIR if tag == "corrected"
                             else config.DEMUX_DIR_ORIGINAL)
    with open(QC / f"q02_well_representatives_{tag}.pkl", "rb") as fh:
        blob = pickle.load(fh)
    well_df = pd.read_csv(QC / f"q01_well_composition_{tag}.csv")
    print(f"{len(blob['reps'])} wells carry representatives "
          f"(clustered at {blob['threshold']})")

    summary, consensus = strain_consensus.build(cfg, blob["reps"], well_df)
    strain_consensus.report(summary)
    summary.to_csv(OUT / f"c01_strain_consensus_summary_{tag}.csv", index=False)
    strain_consensus.write_fasta(OUT / f"c02_strain_consensus_{tag}.fasta",
                                 summary, consensus)
    with open(OUT / f"c03_consensus_seqs_{tag}.pkl", "wb") as fh:
        pickle.dump(consensus, fh)

    print("\ninternal coherence (independent of any external reference):")
    good = {r.strain: consensus[r.strain] for r in summary.itertuples()
            if r.status in strain_consensus.GOOD and r.strain in consensus}
    from io_utils import load_reference_db
    ref = load_reference_db(cfg.reference_dbs["full_collection"])
    strain_identity.self_consistency(good, ref)
    return summary, consensus


if __name__ == "__main__":
    main("original" if "--original" in sys.argv else "corrected")
