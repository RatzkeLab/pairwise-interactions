"""Pick -e and -l using the frozen plates as a negative control.

Plates 9-16 were never PCR'd. Any read minibar assigns to one of their 2464
wells is a false assignment, and they are half the sample sheet, so doubling
that count estimates the total false rate. That makes the edit-distance choice
a measurement rather than a guess.

Result on a representative chunk (10,000 filtered reads from chunk 30), with the
corrected barcode file. False rate is measured, not assumed: plates 9-16 were
never PCR'd, so reads landing there are false, and they are half the sheet.

      barcode file           -e   -l   assigned   est. false rate
      as generated            6  150      13.7%            0.15%
      corrected               4  150      63.6%            1.23%
      corrected               5  150      67.6%            1.48%
      corrected               6  150      72.7%            5.15%
      corrected               6  250      78.8%           14.65%

The barcode set has a minimum pairwise edit distance of 12, so -e 6 sits exactly
on the decision boundary and -e 7 or more is meaningless. -e 4 -l 150 is the
operating point: 4.6x the reads the run currently has, at a false rate well
below the well-to-well contamination the plates already carry.
"""
import subprocess, sys, tempfile
from pathlib import Path
import pandas as pd

MINIBAR = "/home/rl/scripts/ont-demultiplex-pipeline/workflow/scripts/minibar_parallel.py"
PY = "/home/rl/mambaforge/envs/minibar_env/bin/python"
SETUP = Path("/home/rl/scripts/karl/pairwise_interaction_experiments/20260907/01_setup")


def run(barcodes, fastq, e, E, l):
    out = subprocess.run([PY, MINIBAR, str(barcodes), str(fastq),
                          "-e", str(e), "-E", str(E), "-l", str(l), "-w", "-S"],
                         capture_output=True, text=True)
    ids = [ln.split("|")[-1].strip()
           for i, ln in enumerate(out.stdout.splitlines()) if i % 4 == 0]
    return [x for x in ids if x], out.stdout.count("\n") // 4


def calibrate(barcodes_all, fastq, grid=((4, 150), (4, 250), (5, 150), (5, 250),
                                         (6, 150), (6, 250), (7, 250))):
    lay = pd.read_csv(SETUP / "strain_layout_20260907.csv")
    lay["s"] = lay.apply(lambda r: f"Plate{int(r.dest_plate):02d}_{r.dest_well}", axis=1)
    primary = set(lay[lay.plate_role == "primary_sequenced"].s)
    rows = []
    for e, l in grid:
        ids, n = run(barcodes_all, fastq, e, 6, l)
        called = [i for i in ids if i not in ("unk", "Multiple_Matches")]
        false = sum(1 for i in called if i not in primary)
        rows.append(dict(e=e, l=l, reads=n, assigned=len(called),
                         assigned_frac=len(called) / n,
                         frozen_hits=false,
                         est_false_rate=2 * false / max(len(called), 1)))
        print(f"  -e {e} -l {l}: {len(called)}/{n} assigned "
              f"({rows[-1]['assigned_frac']:.1%}), est false {rows[-1]['est_false_rate']:.2%}")
    return pd.DataFrame(rows)


if __name__ == "__main__":
    df = calibrate(Path(sys.argv[1]), Path(sys.argv[2]))
    df.to_csv(Path(__file__).parent / "outputs" / "d09_minibar_calibration.csv", index=False)
