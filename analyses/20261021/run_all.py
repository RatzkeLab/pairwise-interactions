"""Build the 20260907 reference DB and the distinguishability tables, in order (~3 min).

    TMPDIR=<scratch> python run_all.py

Order matters: stage 2 seeds its read assignment from stage 1's monoculture
populations, stage 3 pools both, stage 4 needs the haplotype FASTA from stage 3, and
stage 5 needs stage 4's pair list. Each step writes its outputs, so any one can be
re-run alone once the steps above it have run once.
"""
import sys
import time
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
sys.path.insert(0, str(BASE / "01_reference_db"))
sys.path.insert(0, str(BASE / "02_distinguishability"))

import config  # noqa: E402,F401
import b01_mono_populations, b02_pair_observations, b03_strain_reference  # noqa: E402
import d01_distances, d02_read_resolution, d03_separable_sets  # noqa: E402


def banner(n, title):
    print("\n" + "#" * 78 + f"\n# {n}. {title}\n" + "#" * 78)


if __name__ == "__main__":
    t0 = time.time()
    banner(1, "monoculture populations, base-exact consensus")
    b01_mono_populations.main()
    banner(2, "pair wells: each strain's own consensus, seeded read assignment")
    b02_pair_observations.main()
    banner(3, "one reference per strain, by independent-well agreement; haplotypes")
    b03_strain_reference.main()
    banner(4, "all-against-all distances, identity groups, operon evidence")
    d01_distances.main()
    banner(5, "read-level resolution, calibrated on this run's monocultures")
    d02_read_resolution.main()
    banner(6, "largest mutually separable set at each threshold")
    d03_separable_sets.main()
    print(f"\ndone in {time.time() - t0:.0f}s")
