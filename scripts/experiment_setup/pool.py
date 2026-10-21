"""Which strains a community may contain.

Built on analyses/20261021 (the 20260907 reference DB and the genome checks):

  core      `genomic_use == high`: genome verified against the well, and no genomically
            distinct organism on file shares its 16S -- the readout pins down the genome
  extended  `genomic_use == medium`: genome verified, but a 16S-identical relative with a
            different genome exists, and 78 of 97 such strains share their 16S with an
            organism that has taken over at least one other well

Wells that hold the same genome (a clone cluster) are one strain here: they are
16S-identical, so a community could only ever hold one of them, and a genome-based model
must count them once. The representative is the best-evidenced well (tier A > B > C, no
competing sequence, then alphabetical).
"""
import numpy as np
import pandas as pd

from . import plates

TIER_RANK = {"A": 0, "B": 1, "C": 2}


def source_plate_od(files, wavelength=600):
    """Per-well OD of the 384 source plate, median over repeat reads (single reads
    saturate at the reader ceiling, 3.5). The source plate's well coordinate IS the
    strain name."""
    import plate_reader as pr
    grids = [pr.read_plate(f, wavelength)[0] for f in files]
    grids = [g for g in grids if g is not None and g.shape == (plates.N_ROWS, plates.N_COLS)]
    if not grids:
        raise ValueError(f"no 16x24 OD grid in {files}")
    med = np.median(np.dstack(grids), axis=2)
    return pd.Series({f"{plates.ROWS[r]}{c + 1}": float(med[r, c])
                      for r in range(plates.N_ROWS) for c in range(plates.N_COLS)},
                     name="source_od")


def load_pool(eligibility_csv, source_od=None, min_source_od=None, include_extended=True):
    """One row per usable genome, with `pool` = core / extended."""
    e = pd.read_csv(eligibility_csv)
    keep = ["high", "medium"] if include_extended else ["high"]
    e = e[e.genomic_use.isin(keep)].copy()
    e["tier_rank"] = e.tier.map(TIER_RANK)
    e["competing_rank"] = e.competing.astype(int)
    reps = (e.sort_values(["genomic_unit", "tier_rank", "competing_rank", "strain"])
              .groupby("genomic_unit", as_index=False).first())
    reps["pool"] = np.where(reps.genomic_use == "high", "core", "extended")
    reps["clone_mates_not_used"] = [
        " ".join(sorted(set(e.strain[e.genomic_unit == u]) - {s}))
        for u, s in zip(reps.genomic_unit, reps.strain)]
    if source_od is not None:
        reps["source_od"] = reps.strain.map(source_od)
        if min_source_od is not None:
            reps["low_preculture_density"] = reps.source_od < min_source_od
    cols = ["strain", "pool", "genomic_unit", "tier", "status", "identity_group",
            "n_haplotypes", "haplotype_divergence_bp", "clone_mates_not_used",
            "source_od", "low_preculture_density"]
    cols = [c for c in cols if c in reps]
    return reps[cols].sort_values(["pool", "strain"]).reset_index(drop=True)


def distance_matrix(path, names):
    D = pd.read_csv(path, index_col=0)
    return D.loc[names, names].astype(float)
