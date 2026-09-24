"""Re-run 04 -> 05 -> 06 for one demultiplexing arm, headless.

    python run_all.py                       # corrected demux  -> outputs/
    python run_all.py --arm original        # original demux   -> outputs_original_demux/
    python run_all.py --steps ra,ml         # subset: qc, ra, ml, plate, depth

Env: karl_seq_analysis. The notebooks in 04/05 are the same calls with figures shown
inline; this driver exists so both arms run through identical code. Not re-run here:
03 (the reference is held fixed, see config.py), 04 check_reads (superseded by
scripts/well_composition.py), genomic_ml_yield (OD-only, independent of the demux),
feature_sweep, and the hierarchy_significance suite (slow; the hierarchy summary itself
is in r05_hierarchy_summary.csv).
"""
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--arm", default="corrected", choices=["corrected", "original", "pass_b"])
ap.add_argument("--steps", default="qc,ra,ml,plate,depth")
args = ap.parse_args()
steps = args.steps.split(",")

os.environ["DEMUX_ARM"] = args.arm          # read by config.py at import
BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
from config import CFG  # noqa: E402  (also puts scripts/ on sys.path)

print(f"arm={args.arm}  demux_dir={CFG.demux_dir}")
t0 = time.time()


def stamp(msg):
    print(f"\n[{time.time() - t0:7.0f}s] ===== {msg} =====", flush=True)


if "qc" in steps:
    import mapping_validation as mv
    stamp("04 mapping_validation")
    samples = mv.gather_samples(CFG)
    mv.prepare_references(CFG)
    mv.constrained_edlib_mapping(CFG, samples)
    mv.minimap2_contamination_scan(CFG, samples)
    summary = mv.combine_and_flag(CFG)
    print(summary.groupby(["well_type", "qc_status"]).size().to_string())
    mv.make_all_figures(CFG)

if "ra" in steps:
    import relative_abundance as ra
    stamp("05 relative_abundance")
    ra.reference_distances(CFG)
    ra.compute_interaction_scores(CFG)
    ra.replicate_stability(CFG)
    ra.competitiveness_scores(CFG)
    _, _, _, hsummary = ra.hierarchy_analysis(CFG)
    print(hsummary.to_string())
    ra.make_all_figures(CFG)   # after r03-r05: the figures read them (the notebook
                               # only got away with its order because old outputs existed)

scripts = {
    "ml": ["06_ml_analysis/genomic_ml/run_genomic_ml.py"],
    "plate": ["06_ml_analysis/genomic_ml_plate/run_genomic_ml_plate.py"],
    "depth": ["06_ml_analysis/genomic_ml/replicate_value_relabund.py"],
}
for step, cmd in scripts.items():
    if step in steps:
        stamp(f"06 {cmd[0]}")
        subprocess.run([sys.executable, str(BASE / cmd[0]), *cmd[1:]], check=True,
                       cwd=BASE / Path(cmd[0]).parent, env=os.environ)

stamp("done")
