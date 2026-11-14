"""Stage 11: build the per-well table, then every plate map (~15 s).

    python make_figures.py

p00 must run first; each f* script only reads outputs/p00_*.csv and earlier stages' outputs.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p00_well_table, f01_overgrowth, f02_secondary_contamination  # noqa: E402
import f03_indistinguishable_groups, f04_identical_on_purpose, f05_genomic_use  # noqa: E402
import f06_clustermap  # noqa: E402


def main():
    p00_well_table.main()
    f01_overgrowth.main()
    f02_secondary_contamination.main()
    f03_indistinguishable_groups.main()
    f04_identical_on_purpose.fig_a()
    f04_identical_on_purpose.fig_b()
    f05_genomic_use.main()
    f06_clustermap.main()


if __name__ == "__main__":
    main()
