"""Build the 20260907 reference DB, the distinguishability tables and the genome checks,
in order (~9 min).

    TMPDIR=<scratch> python run_all.py

Order matters: stage 2 seeds its read assignment from stage 1's monoculture
populations, stage 3 pools both, stage 4 needs the haplotype FASTA from stage 3, and
stage 5 needs stage 4's pair list; stages 7-10 need the references and the genomes
(GenBank -> FASTA cache under $TMPDIR, rebuilt in seconds). Each step writes its outputs, so any one can be
re-run alone once the steps above it have run once.
"""
import sys
import time
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
sys.path.insert(0, str(BASE / "01_reference_db"))
sys.path.insert(0, str(BASE / "02_distinguishability"))
sys.path.insert(0, str(BASE / "03_genome_concordance"))

import config  # noqa: E402,F401
import b01_mono_populations, b02_pair_observations, b03_strain_reference  # noqa: E402
import d01_distances, d02_read_resolution, d03_separable_sets  # noqa: E402
import g01_genomes, g03_assembly_purity, g02_verdicts, g04_design_eligibility  # noqa: E402


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
    banner(7, "genomes: extract, all-vs-all ANI, 16S of every well against every assembly")
    g01_genomes.main()
    banner(8, "assembly purity: one organism per assembly?")
    g03_assembly_purity.main()
    banner(9, "does each well's genome describe what is in the well?")
    g02_verdicts.main()
    banner(10, "design eligibility for genome-based communities")
    g04_design_eligibility.main()
    print(f"\ndone in {time.time() - t0:.0f}s")
