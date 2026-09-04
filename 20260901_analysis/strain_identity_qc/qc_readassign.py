"""Calibrate the two-reference read-assignment rule against mono-well ground truth.

The interaction pipeline decides, per read, which of a well's two strains it came from. Its
original rule called a read ambiguous when the *normalized* distance margin fell below 0.02 --
about 28 bp on a 1420 bp read. Two references 10 bp apart can never produce a margin that
large, so every read in a close pair was filed "ambiguous" by construction, and the pair was
flagged unresolvable no matter how clean the data were. That is a property of the threshold,
not of the measurement.

This module replaces the guess with a measurement. Mono wells contain one intended strain, so
the correct answer is known: score their reads against the true strain and against every other
strain, and see how far apart two references must be before "assign to the nearer one" is
reliable.

Two corrections matter for the ground truth to mean anything:

  * 20260630's mono wells are not pure -- by design they were shot against wells believed to be
    no-growers, some of which grew. So a read is only counted if its best match across the WHOLE
    reference set is the intended strain; contaminant reads are excluded rather than being
    forced into a binary choice they don't belong in.
  * Reads far from every reference are dropped on the same off-target threshold the pipeline
    learns for itself.

The theory says 5 bp should be plenty (at 4.7% per-base error, per-read accuracy at d=5 is
99.97% under a substitution model). The data disagree, because ONT error is indel-dominated and
concentrated in homopolymers -- which is exactly where two nearly identical 16S sequences
differ. Trust the mono wells, not the binomial.
"""

import numpy as np
import pandas as pd
import edlib

import qc_config as C
import qc_sources as S


def _reads(path, lo=1000, hi=1800):
    lines = open(path).read().splitlines()
    return [lines[i + 1] for i in range(0, len(lines) - 3, 4) if lo <= len(lines[i + 1]) <= hi]


def mono_ground_truth(exp, ref_fasta, off_target_norm=0.20, min_reads=20):
    """-> per (mono strain, other strain) accuracy of 'assign to the nearer reference'.

    Also returns the per-read table so a margin rule can be calibrated, not just a distance cut.
    """
    ref = {n.split("|")[0]: s for n, s in S.read_fasta(ref_fasta)}
    lay = pd.read_csv(exp.layout_csv)
    mono = lay[lay["well_type"] == "mono"]
    names = list(ref)

    rows, read_rows = [], []
    for r in mono.itertuples():
        X = str(r.strain1 if pd.notna(r.strain1) else r.strain2)
        if X not in ref:
            continue
        path = exp.demux_dir / ("Plate%02d_%s.fastq" % (r.dest_plate, r.dest_well))
        if not path.exists():
            continue
        reads = _reads(path)
        if len(reads) < min_reads:
            continue

        # distance of every read to every reference, once
        D = np.array([[edlib.align(q, ref[n], mode="NW", task="distance")["editDistance"]
                       for n in names] for q in reads], dtype=float)
        L = np.array([max(len(q), 1400) for q in reads], dtype=float)[:, None]
        norm = D / L
        xi = names.index(X)

        # Ground truth by ABSOLUTE closeness to the intended strain, never by argmin: selecting
        # reads whose global best match is X would force d(X) < d(Y) for every Y and make the
        # accuracy 1.0 by construction. An absolute cut still drops the no-grower contaminants
        # (a read from a third organism is far from X) without deciding the comparison in advance.
        keep = norm[:, xi] <= off_target_norm
        if keep.sum() < min_reads:
            continue
        for yi, Y in enumerate(names):
            if Y == X:
                continue
            bp = edlib.align(ref[X], ref[Y], mode="NW", task="distance")["editDistance"]
            if bp == 0 or bp > 250:
                continue
            dx, dy = D[keep, xi], D[keep, yi]
            rows.append({"mono_strain": X, "other": Y, "ref_bp_dist": bp,
                         "n_reads_true": int(keep.sum()),
                         "frac_correct": float((dx < dy).mean()),
                         "frac_tied": float((dx == dy).mean())})
            read_rows.append(pd.DataFrame({"ref_bp_dist": bp, "margin_bp": np.abs(dx - dy),
                                           "correct": (dx < dy).astype(int),
                                           "tied": (dx == dy).astype(int)}))

    per_pair = pd.DataFrame(rows)
    per_read = pd.concat(read_rows, ignore_index=True) if read_rows else pd.DataFrame()
    per_pair.to_csv(C.OUT / f"s11_mono_discrimination_{exp.name}.csv", index=False)
    per_read.to_csv(C.OUT / f"s12_mono_discrimination_reads_{exp.name}.csv.gz", index=False)
    return per_pair, per_read


def read_level_curve(per_read, bins=((1, 2), (3, 5), (6, 10), (11, 20), (21, 40), (41, 80), (81, 250))):
    """Read-level accuracy binned by reference separation -- far more power than the per-pair
    view at the small distances that actually decide the experimental design."""
    rows = []
    for lo, hi in bins:
        m = (per_read["ref_bp_dist"] >= lo) & (per_read["ref_bp_dist"] <= hi)
        if not m.sum():
            continue
        rows.append({"ref_bp_range": f"{lo}-{hi}", "n_read_tests": int(m.sum()),
                     "pct_reads_correct": round(100 * per_read.loc[m, "correct"].mean(), 2),
                     "pct_tied": round(100 * per_read.loc[m, "tied"].mean(), 2)})
    return pd.DataFrame(rows)


def resolvability_curve(per_pair, bins=((1, 2), (3, 5), (6, 10), (11, 20), (21, 40), (41, 80), (81, 250))):
    rows = []
    for lo, hi in bins:
        m = (per_pair["ref_bp_dist"] >= lo) & (per_pair["ref_bp_dist"] <= hi)
        if not m.sum():
            continue
        rows.append({"ref_bp_range": f"{lo}-{hi}", "n_tests": int(m.sum()),
                     "median_pct_correct": round(100 * per_pair.loc[m, "frac_correct"].median(), 1),
                     "pct_tests_above_95": round(100 * (per_pair.loc[m, "frac_correct"] > 0.95).mean(), 1),
                     "pct_tests_above_90": round(100 * (per_pair.loc[m, "frac_correct"] > 0.90).mean(), 1)})
    return pd.DataFrame(rows)


def margin_rule_scan(per_read, min_bp_options=(1, 3, 6, 8, 10, 12, 15, 20)):
    """Where to put the hard floor: accuracy vs. how many pairs survive it."""
    rows = []
    for mb in min_bp_options:
        m = per_read["ref_bp_dist"] >= mb
        if not m.sum():
            continue
        rows.append({"min_ref_bp_dist": mb,
                     "pct_of_read_tests_retained": round(100 * m.mean(), 1),
                     "read_accuracy": round(100 * per_read.loc[m, "correct"].mean(), 2),
                     "pct_tied": round(100 * per_read.loc[m, "tied"].mean(), 2)})
    return pd.DataFrame(rows)
