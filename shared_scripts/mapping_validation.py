"""Shared, per-experiment-parameterized read-mapping validation pipeline.

For every well with more than `cfg.min_reads` reads, checks whether the reads actually match
the strain(s) the layout says should be there: constrained edlib mapping (does a read match
*only* the well's expected strain(s)) plus an unconstrained minimap2 scan (what does an
off-target read actually look like -- the contamination signature). Originally built for
pairwise_interaction_experiments/20260721 as numbered scripts s01-s06; ported here so every
experiment can call the same, already-validated logic through an ExperimentConfig instead of
copy-pasting it.

Assumption (true for every experiment so far): exactly one reference db in
`cfg.reference_dbs` needs an identity cross-check before being trusted --
`cfg.external_cross_check_dbs[0]` (e.g. corroborated_db, whose strain-name labels are pooled
across many unrelated past experiments and are only sometimes the same organism as this
experiment's own use of that name -- see identity_cross_check()).

Call order: gather_samples -> prepare_references -> constrained_edlib_mapping ->
minimap2_contamination_scan -> combine_and_flag -> make_all_figures. Each step reads/writes
CSVs under cfg.mapping_validation_out_dir so they can be re-run independently once upstream
steps have completed once.
"""

import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

from io_utils import (
    load_layout, sample_fastq_path, count_reads_fast, load_reads, load_reference_db,
    norm_edit_distance,
)

# ---- thresholds: about ONT/16S sequencing statistics, not the experiment, so these are
# shared defaults rather than per-experiment config ----
CROSS_CHECK_RELIABLE_THRESHOLD = 0.10   # normalized edit distance; below this = same organism
MIN_STRAIN_READS = 3
MIN_STRAIN_FRAC = 0.10
CONTAM_MIN_FRAC = 0.20
CONTAM_MIN_READS = 3
REFERENCE_TWIN_THRESHOLD = 0.05         # below this, "contaminant" is a reference near-twin, not real

# ---- fixed semantic figure palette --------------------------------------------------
COLOR_PAIR = "#2a78d6"
COLOR_MONO = "#eb6834"
COLOR_GOOD = "#0ca30c"
COLOR_WARNING = "#fab219"
COLOR_CRITICAL = "#d03b3b"
COLOR_GRID = "#d8d7d2"
COLOR_TEXT_SECONDARY = "#52514e"

QC_COLOR = {
    "mono_confirmed": COLOR_GOOD,
    "pair_both_confirmed": COLOR_GOOD,
    "mono_low_confidence": COLOR_WARNING,
    "pair_single_dominant": COLOR_WARNING,
    "pair_neither_confirmed": COLOR_CRITICAL,
}

plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": COLOR_GRID,
        "axes.grid": True,
        "grid.color": COLOR_GRID,
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
    }
)


# ===========================================================================
# 01 -- gather the >min_reads cohort
# ===========================================================================


def gather_samples(cfg):
    layout = load_layout(cfg.layout_csv)
    layout["path"] = layout.apply(
        lambda r: sample_fastq_path(cfg.demux_dir, r["dest_plate"], r["dest_well"]), axis=1
    )
    n_missing = layout["path"].isna().sum()
    layout["n_reads"] = layout["path"].apply(count_reads_fast)

    gt = layout[layout["n_reads"] > cfg.min_reads].copy()
    gt["path"] = gt["path"].astype(str)

    cols = ["sample_id", "dest_plate", "dest_row", "dest_col", "dest_well", "well_type", "strain1", "strain2", "n_reads", "path"]
    gt = gt[cols].sort_values("sample_id").reset_index(drop=True)

    out_path = cfg.mapping_validation_out_dir / "01_samples_gt5reads.csv"
    gt.to_csv(out_path, index=False)

    print(f"layout wells               : {len(layout)}")
    print(f"  missing fastq files      : {n_missing}")
    print(f"wells with > {cfg.min_reads} reads      : {len(gt)}")
    print(f"  well_type breakdown      : {gt['well_type'].value_counts().to_dict()}")
    print(f"  total reads in cohort    : {gt['n_reads'].sum()}")
    print(f"saved -> {out_path}")
    return gt


# ===========================================================================
# 02 -- reference coverage + identity cross-check
# ===========================================================================


def prepare_references(cfg):
    layout = load_layout(cfg.layout_csv)
    strains = sorted(set(layout["strain1"]) | set(layout["strain2"]))

    dbs = {name: load_reference_db(path) for name, path in cfg.reference_dbs.items()}

    rows = []
    for strain in strains:
        row = {"strain": strain}
        for db_name, seqs in dbs.items():
            row[f"{db_name}_len"] = len(seqs[strain]) if strain in seqs else None
        rows.append(row)
    cov = pd.DataFrame(rows)

    out_path = cfg.mapping_validation_out_dir / "02_reference_coverage.csv"
    cov.to_csv(out_path, index=False)

    print(f"strains in layout: {len(strains)}")
    for db_name, seqs in dbs.items():
        present = sum(1 for s in strains if s in seqs)
        extra = len(seqs) - present
        print(
            f"  {db_name:32s} {cfg.reference_dbs[db_name]}\n"
            f"    covers {present}/{len(strains)} strains"
            + (f", +{extra} extra entries not in this layout" if extra else "")
        )
    print(f"saved -> {out_path}")

    # ---- identity cross-check: does the external db's name match this experiment's own? ----
    external_db = cfg.external_cross_check_dbs[0]
    ext = dbs[external_db]
    primary = dbs[cfg.primary_db]
    cc_rows = []
    for strain in strains:
        if strain in ext and strain in primary:
            d = norm_edit_distance(ext[strain], primary[strain])
            cc_rows.append({"strain": strain, "norm_dist_corrdb_vs_own_consensus": d})
    cc = pd.DataFrame(cc_rows).sort_values("norm_dist_corrdb_vs_own_consensus").reset_index(drop=True)
    cc["reliable"] = cc["norm_dist_corrdb_vs_own_consensus"] < CROSS_CHECK_RELIABLE_THRESHOLD

    cc_path = cfg.mapping_validation_out_dir / "02_reference_cross_check.csv"
    cc.to_csv(cc_path, index=False)

    n_reliable = int(cc["reliable"].sum())
    print(
        f"\n{external_db} identity cross-check ({len(cc)} strains comparable):\n"
        f"  {n_reliable}/{len(cc)} strain names agree with this experiment's own consensus "
        f"(< {CROSS_CHECK_RELIABLE_THRESHOLD} normalized edit distance)\n"
        f"  {len(cc) - n_reliable}/{len(cc)} strain names in {external_db} do NOT match this "
        f"experiment's sequence for that name -- likely coincidental reuse of a plate-well-style "
        f"label from an unrelated experiment, not the same organism.\n"
        f"  -> treating {cfg.primary_db} as the PRIMARY reference for QC calls; "
        f"{external_db} is reported only as a secondary check, restricted to reliable strains."
    )
    print(f"saved -> {cc_path}")
    return dbs, cov, cc


# ===========================================================================
# 03 -- constrained edlib mapping (read vs. only the well's expected strain(s))
# ===========================================================================


def map_sample_constrained(read_seqs, candidates, ref_seqs):
    """candidates: list of distinct expected strain names present in ref_seqs (1 or 2).
    Returns list of dicts, one per read, with distances to each candidate + winner."""
    rows = []
    for read_id, seq in read_seqs:
        dists = {c: norm_edit_distance(seq, ref_seqs[c]) for c in candidates}
        best_strain = min(dists, key=dists.get)
        best_dist = dists[best_strain]
        sorted_d = sorted(dists.values())
        margin = sorted_d[1] - sorted_d[0] if len(sorted_d) > 1 else float("nan")
        row = {
            "read_id": read_id, "read_len": len(seq),
            "best_match_strain": best_strain, "best_dist": best_dist, "margin_to_runner_up": margin,
        }
        for i, c in enumerate(candidates, start=1):
            row[f"candidate{i}_strain"] = c
            row[f"candidate{i}_dist"] = dists[c]
        rows.append(row)
    return rows


def _run_constrained_mapping(cfg, samples_df, db_name, ref_path):
    ref_seqs = load_reference_db(ref_path)
    out_rows = []
    n_no_ref = 0
    for _, r in samples_df.iterrows():
        candidates = list(dict.fromkeys([r["strain1"], r["strain2"]]))
        candidates = [c for c in candidates if c in ref_seqs]
        if not candidates:
            n_no_ref += 1
            continue
        reads = load_reads(r["path"])
        for row in map_sample_constrained(reads, candidates, ref_seqs):
            row["sample_id"] = r["sample_id"]
            row["well_type"] = r["well_type"]
            row["strain1"] = r["strain1"]
            row["strain2"] = r["strain2"]
            row["n_candidates"] = len(candidates)
            out_rows.append(row)

    df = pd.DataFrame(out_rows)
    out_path = cfg.mapping_validation_out_dir / f"03_edlib_read_assignments_{db_name}.csv.gz"
    df.to_csv(out_path, index=False)
    print(
        f"[{db_name}] {len(samples_df) - n_no_ref}/{len(samples_df)} wells had a reference "
        f"for >=1 expected strain ({n_no_ref} skipped, no reference available); "
        f"{len(df)} reads classified -> {out_path}"
    )
    return df


def constrained_edlib_mapping(cfg, samples_df):
    return {db_name: _run_constrained_mapping(cfg, samples_df, db_name, ref_path) for db_name, ref_path in cfg.reference_dbs.items()}


# ===========================================================================
# 04 -- unconstrained minimap2 scan (read vs. every strain in a db)
# ===========================================================================

PAF_COLS = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alnlen", "mapq"]


def _write_combined_reads_fastq(samples_df, fastq_out, universe_out):
    if fastq_out.exists() and universe_out.exists():
        return fastq_out, pd.read_csv(universe_out)

    rows = []
    n = 0
    with open(fastq_out, "w") as fh:
        for _, r in samples_df.iterrows():
            for read_id, seq in load_reads(r["path"]):
                qname = f"{r['sample_id']}__{read_id}"
                fh.write(f"@{qname}\n{seq}\n+\n{'I' * len(seq)}\n")
                rows.append({"sample_id": r["sample_id"], "read_id": read_id, "read_len": len(seq)})
                n += 1
    universe = pd.DataFrame(rows)
    universe.to_csv(universe_out, index=False)
    print(f"wrote {n} reads -> {fastq_out} (read universe cached -> {universe_out})")
    return fastq_out, universe


def _write_normalized_reference(ref_path, out_path):
    seqs = load_reference_db(ref_path)
    with open(out_path, "w") as fh:
        for strain, seq in seqs.items():
            fh.write(f">{strain}\n{seq}\n")
    return out_path


def _run_minimap2(ref_fasta, reads_fastq, paf_out, threads=8):
    cmd = ["minimap2", "-x", "map-ont", "--secondary=no", "-t", str(threads), str(ref_fasta), str(reads_fastq)]
    with open(paf_out, "w") as out_fh:
        proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"minimap2 failed:\n{proc.stderr}")
    return paf_out


def _parse_paf(paf_path):
    if paf_path.stat().st_size == 0:
        return pd.DataFrame(columns=PAF_COLS + ["identity"])
    df = pd.read_csv(paf_path, sep="\t", header=None, usecols=range(12), names=PAF_COLS, engine="c")
    df["identity"] = df["nmatch"] / df["alnlen"]
    return df


def _best_hit_table(samples_df, read_universe, paf_df):
    paf_df = paf_df.copy()
    paf_df[["sample_id", "read_id"]] = paf_df["qname"].str.split("__", n=1, expand=True)
    paf_df = paf_df.rename(columns={"tname": "best_hit_strain"})

    hit = paf_df[["sample_id", "read_id", "best_hit_strain", "nmatch", "alnlen", "identity", "mapq"]].copy()
    hit = read_universe.merge(hit, on=["sample_id", "read_id"], how="left")
    hit["best_hit_strain"] = hit["best_hit_strain"].where(hit["best_hit_strain"].notna(), None)
    for c in ("nmatch", "alnlen", "mapq"):
        hit[c] = hit[c].fillna(0).astype(int)
    hit["identity"] = hit["identity"].fillna(0.0)

    expect = samples_df.set_index("sample_id")[["strain1", "strain2"]]
    hit = hit.join(expect, on="sample_id")
    hit["is_expected"] = (hit["best_hit_strain"] == hit["strain1"]) | (hit["best_hit_strain"] == hit["strain2"])
    return hit


def _run_contamination_scan(cfg, samples_df, read_universe, db_name, ref_path):
    norm_ref = _write_normalized_reference(ref_path, cfg.mapping_validation_out_dir / f"04_reference_normalized_{db_name}.fasta")
    reads_fastq = cfg.mapping_validation_out_dir / "04_reads_combined.fastq"
    paf_out = cfg.mapping_validation_out_dir / f"04_minimap2_{db_name}.paf"
    _run_minimap2(norm_ref, reads_fastq, paf_out)

    paf_df = _parse_paf(paf_out)
    hit = _best_hit_table(samples_df, read_universe, paf_df)

    out_path = cfg.mapping_validation_out_dir / f"04_minimap2_read_besthit_{db_name}.csv"
    hit.to_csv(out_path, index=False)

    n_reads = len(hit)
    n_mapped = (hit["best_hit_strain"].notna()).sum()
    n_expected = hit["is_expected"].sum()
    print(
        f"[{db_name}] {n_mapped}/{n_reads} reads mapped ({n_mapped/n_reads:.1%}); "
        f"{n_expected}/{n_reads} best-hit an expected strain ({n_expected/n_reads:.1%}) -> {out_path}"
    )
    return hit


def minimap2_contamination_scan(cfg, samples_df):
    reads_fastq, read_universe = _write_combined_reads_fastq(
        samples_df, cfg.mapping_validation_out_dir / "04_reads_combined.fastq", cfg.mapping_validation_out_dir / "04_read_universe.csv"
    )
    return {db_name: _run_contamination_scan(cfg, samples_df, read_universe, db_name, ref_path) for db_name, ref_path in cfg.reference_dbs.items()}


# ===========================================================================
# 05 -- combine + QC-flag every well
# ===========================================================================


def find_valley_threshold(values, lo=0.04, hi=0.20, bins=81, empty_frac=0.02):
    """Data-driven same/different threshold: centre of the emptiest band in [lo, hi]
    of the pooled best-match distance histogram."""
    values = np.asarray(values)
    values = values[~np.isnan(values)]
    hist, edges = np.histogram(values, bins=np.linspace(0, max(0.4, values.max()), bins))
    centers = 0.5 * (edges[:-1] + edges[1:])
    band = (centers >= lo) & (centers <= hi)
    empty = band & (hist <= max(1, empty_frac * hist.max()))
    if empty.any():
        return round(float(np.median(centers[empty])), 3)
    return round(float(np.median(values)), 3)


def _summarize_edlib_db(cfg, db_name, threshold):
    reads = pd.read_csv(cfg.mapping_validation_out_dir / f"03_edlib_read_assignments_{db_name}.csv.gz")
    reads["confident"] = reads["best_dist"] < threshold

    rows = []
    for sid, g in reads.groupby("sample_id"):
        strain1, strain2 = g["strain1"].iat[0], g["strain2"].iat[0]
        n = len(g)
        candidates_present = set(g["candidate1_strain"]) | (set(g["candidate2_strain"].dropna()) if "candidate2_strain" in g else set())
        gc = g[g["confident"]]

        def frac_for(strain):
            return (gc["best_match_strain"] == strain).sum() / n if strain in candidates_present else np.nan

        def n_for(strain):
            return int((gc["best_match_strain"] == strain).sum()) if strain in candidates_present else np.nan

        rows.append({
            "sample_id": sid, f"edlib_{db_name}_threshold": threshold,
            f"edlib_{db_name}_n_conf_strain1": n_for(strain1), f"edlib_{db_name}_n_conf_strain2": n_for(strain2),
            f"edlib_{db_name}_frac_strain1": frac_for(strain1), f"edlib_{db_name}_frac_strain2": frac_for(strain2),
            f"edlib_{db_name}_median_best_dist": g["best_dist"].median(),
        })
    return pd.DataFrame(rows)


def _summarize_minimap2_db(cfg, db_name):
    hits = pd.read_csv(cfg.mapping_validation_out_dir / f"04_minimap2_read_besthit_{db_name}.csv")
    rows = []
    for sid, g in hits.groupby("sample_id"):
        n = len(g)
        frac_expected = g["is_expected"].sum() / n
        frac_unmapped = g["best_hit_strain"].isna().sum() / n
        other = g[(~g["is_expected"]) & g["best_hit_strain"].notna()]
        if len(other):
            vc = other["best_hit_strain"].value_counts()
            top_other_strain, top_other_count = vc.index[0], int(vc.iat[0])
        else:
            top_other_strain, top_other_count = None, 0
        rows.append({
            "sample_id": sid, f"mm2_{db_name}_frac_expected": frac_expected, f"mm2_{db_name}_frac_unmapped": frac_unmapped,
            f"mm2_{db_name}_top_other_strain": top_other_strain, f"mm2_{db_name}_top_other_count": top_other_count,
            f"mm2_{db_name}_top_other_frac": top_other_count / n,
        })
    return pd.DataFrame(rows)


def _assign_qc_status(cfg, row):
    from math import ceil
    n = row["n_reads"]
    cutoff = max(MIN_STRAIN_READS, ceil(MIN_STRAIN_FRAC * n))
    n1 = row[f"edlib_{cfg.primary_db}_n_conf_strain1"]
    n2 = row[f"edlib_{cfg.primary_db}_n_conf_strain2"]
    present1 = (not pd.isna(n1)) and n1 >= cutoff
    present2 = (not pd.isna(n2)) and n2 >= cutoff

    if row["well_type"] == "mono":
        status = "mono_confirmed" if present1 else "mono_low_confidence"
    else:
        if present1 and present2:
            status = "pair_both_confirmed"
        elif present1 or present2:
            status = "pair_single_dominant"
        else:
            status = "pair_neither_confirmed"

    contaminated = (
        row.get(f"mm2_{cfg.primary_db}_top_other_frac", 0) >= CONTAM_MIN_FRAC
        and row.get(f"mm2_{cfg.primary_db}_top_other_count", 0) >= CONTAM_MIN_READS
    )
    return pd.Series({"qc_cutoff_reads": cutoff, "qc_status": status, "qc_contaminated": contaminated})


def combine_and_flag(cfg):
    samples = pd.read_csv(cfg.mapping_validation_out_dir / "01_samples_gt5reads.csv")
    summary = samples.copy()

    primary_reads = pd.read_csv(cfg.mapping_validation_out_dir / f"03_edlib_read_assignments_{cfg.primary_db}.csv.gz")
    shared_threshold = find_valley_threshold(primary_reads["best_dist"].values)
    print(f"confident-match threshold (learned from {cfg.primary_db}, applied to all dbs): {shared_threshold}")

    for db_name in cfg.reference_dbs:
        summary = summary.merge(_summarize_edlib_db(cfg, db_name, shared_threshold), on="sample_id", how="left")
        summary = summary.merge(_summarize_minimap2_db(cfg, db_name), on="sample_id", how="left")

    summary = pd.concat([summary, summary.apply(lambda row: _assign_qc_status(cfg, row), axis=1)], axis=1)

    external_db = cfg.external_cross_check_dbs[0]
    cross_check = pd.read_csv(cfg.mapping_validation_out_dir / "02_reference_cross_check.csv").set_index("strain")["reliable"]
    summary["corrdb_strain1_reliable"] = summary["strain1"].map(cross_check)
    summary["corrdb_strain2_reliable"] = summary["strain2"].map(cross_check)

    out_path = cfg.mapping_validation_out_dir / "05_combined_sample_summary.csv"
    summary.to_csv(out_path, index=False)

    contam = summary[summary["qc_contaminated"]].sort_values(f"mm2_{cfg.primary_db}_top_other_frac", ascending=False).copy()

    ref_seqs = load_reference_db(cfg.reference_dbs[cfg.primary_db])
    twin_dist_cache = {}

    def twin_dist(row):
        other = row[f"mm2_{cfg.primary_db}_top_other_strain"]
        if other is None or pd.isna(other) or other not in ref_seqs:
            return np.nan
        dists = []
        for expected in (row["strain1"], row["strain2"]):
            if expected in ref_seqs:
                key = tuple(sorted((expected, other)))
                if key not in twin_dist_cache:
                    twin_dist_cache[key] = norm_edit_distance(ref_seqs[expected], ref_seqs[other])
                dists.append(twin_dist_cache[key])
        return min(dists) if dists else np.nan

    contam["contaminant_ref_dist_to_expected"] = contam.apply(twin_dist, axis=1)
    contam["likely_reference_twin"] = contam["contaminant_ref_dist_to_expected"] < REFERENCE_TWIN_THRESHOLD

    contam_cols = [
        "sample_id", "well_type", "strain1", "strain2", "n_reads",
        f"mm2_{cfg.primary_db}_top_other_strain", f"mm2_{cfg.primary_db}_top_other_count",
        f"mm2_{cfg.primary_db}_top_other_frac", "contaminant_ref_dist_to_expected", "likely_reference_twin", "qc_status",
    ]
    contam_path = cfg.mapping_validation_out_dir / "05_contamination_candidates.csv"
    contam[contam_cols].to_csv(contam_path, index=False)

    n_twin = int(contam["likely_reference_twin"].sum())
    print(f"combined summary: {len(summary)} wells -> {out_path}")
    print("\nqc_status counts:")
    print(summary["qc_status"].value_counts().to_string())
    print(
        f"\ncontamination-flagged wells: {len(contam)} -> {contam_path}\n"
        f"  of which {n_twin} are reference near-twins (<{REFERENCE_TWIN_THRESHOLD} dist to expected strain's own "
        f"reference -- 16S can't discriminate these, likely not real contamination)\n"
        f"  {len(contam) - n_twin} involve a genuinely divergent unexpected strain -- higher-priority follow-up"
    )
    return summary


# ===========================================================================
# 06 -- figures
# ===========================================================================

DIVERGING_CMAP = LinearSegmentedColormap.from_list("blue_gray_red", [COLOR_PAIR, "#f0efec", COLOR_CRITICAL], N=256)


def _fig00_reference_reliability(cfg, cross_check, fig_dir):
    cc = cross_check.sort_values("norm_dist_corrdb_vs_own_consensus").reset_index(drop=True)
    colors = [COLOR_GOOD if r else COLOR_CRITICAL for r in cc["reliable"]]

    fig, ax = plt.subplots(figsize=(11, 4.6))
    ax.bar(range(len(cc)), cc["norm_dist_corrdb_vs_own_consensus"], color=colors, width=0.8)
    ax.axhline(CROSS_CHECK_RELIABLE_THRESHOLD, color=COLOR_TEXT_SECONDARY, lw=1, ls="--", label=f"reliable cutoff = {CROSS_CHECK_RELIABLE_THRESHOLD}")
    ax.set_xticks(range(len(cc)))
    ax.set_xticklabels(cc["strain"], rotation=90, fontsize=6)
    ax.set_ylabel("norm. edit distance:\nexternal db vs. own consensus")
    n_rel = int(cc["reliable"].sum())
    ax.set_title(f"{cfg.external_cross_check_dbs[0]} name == same organism as this experiment's? {n_rel}/{len(cc)} agree (green)", pad=10)
    fig.text(0.5, 0.955, "short plate-well-style names are reused across unrelated experiments -- most are name collisions",
              ha="center", fontsize=8.5, color=COLOR_TEXT_SECONDARY)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(fig_dir / "00_reference_reliability_check.png", dpi=150)
    plt.close(fig)


def _fig01_match_quality(cfg, reads_primary, summary, fig_dir):
    thr = float(summary[f"edlib_{cfg.primary_db}_threshold"].iloc[0])
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    ax.hist(reads_primary["best_dist"].dropna(), bins=np.linspace(0, 0.4, 81), color=COLOR_PAIR, edgecolor="white", linewidth=0.3)
    ax.axvline(thr, color=COLOR_CRITICAL, ls="--", lw=1.5, label=f"confident-match threshold = {thr}")
    ax.set_yscale("log")
    ax.set_xlabel("normalized edit distance: read → its assigned expected-strain reference")
    ax.set_ylabel("# reads (log scale)")
    ax.set_title(f"Constrained mapping quality ({cfg.primary_db})", pad=10)
    fig.text(0.5, 0.925, "left mode = confirms expected strain, right mode = poor match to either expected strain",
              ha="center", fontsize=8.5, color=COLOR_TEXT_SECONDARY)
    ax.legend(frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(fig_dir / "01_match_quality_histogram.png", dpi=150)
    plt.close(fig)


def _fig02_qc_status_by_well_type(cfg, summary, fig_dir):
    order_mono = ["mono_confirmed", "mono_low_confidence"]
    order_pair = ["pair_both_confirmed", "pair_single_dominant", "pair_neither_confirmed"]

    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=False)
    for ax, order, wt, title in zip(axes, [order_mono, order_pair], ["mono", "pair"], ["mono wells", "pair wells"]):
        counts = summary.loc[summary.well_type == wt, "qc_status"].value_counts().reindex(order, fill_value=0)
        colors = [QC_COLOR[s] for s in order]
        bars = ax.bar(range(len(order)), counts.values, color=colors, width=0.6)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([s.replace("_", "\n") for s in order], fontsize=8)
        ax.set_title(f"{title} (n={int(counts.sum())})")
        ax.set_ylabel("# wells")
        for b, v in zip(bars, counts.values):
            ax.text(b.get_x() + b.get_width() / 2, v, f"{v}", ha="center", va="bottom", fontsize=9)
    fig.suptitle(f"QC call per well ({cfg.primary_db}, constrained mapping)\n"
                 "good = confirmed, warning = only one strain / low confidence, critical = neither expected strain confirmed", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(fig_dir / "02_qc_status_by_well_type.png", dpi=150)
    plt.close(fig)


def _fig03_purity_histogram(cfg, summary, fig_dir):
    summary = summary.copy()
    summary["dominant_frac"] = summary[[f"edlib_{cfg.primary_db}_frac_strain1", f"edlib_{cfg.primary_db}_frac_strain2"]].max(axis=1)

    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.linspace(0, 1, 41)
    for wt, color, label in [("pair", COLOR_PAIR, "pair"), ("mono", COLOR_MONO, "mono")]:
        vals = summary.loc[summary.well_type == wt, "dominant_frac"].dropna()
        ax.hist(vals, bins=bins, color=color, alpha=0.75, label=f"{label} (n={len(vals)})", edgecolor="white", linewidth=0.3)
    ax.set_xlabel("fraction of confident reads matching the dominant expected strain")
    ax.set_ylabel("# wells")
    ax.set_title(f"Read purity vs. designed well type ({cfg.primary_db})")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "03_purity_histogram.png", dpi=150)
    plt.close(fig)


def _fig04_depth_vs_confidence(cfg, summary, fig_dir):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for wt, color, label in [("pair", COLOR_PAIR, "pair"), ("mono", COLOR_MONO, "mono")]:
        sub = summary[summary.well_type == wt]
        ax.scatter(sub["n_reads"], sub[f"mm2_{cfg.primary_db}_frac_expected"], s=14, color=color, alpha=0.5, edgecolor="none", label=f"{label} (n={len(sub)})")
    ax.set_xscale("log")
    ax.set_xlabel("reads in well (log scale)")
    ax.set_ylabel("fraction of reads whose best genome-wide hit\nis an expected strain (minimap2)")
    ax.set_title(f"Read depth vs. mapping confidence ({cfg.primary_db})")
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(fig_dir / "04_depth_vs_confidence.png", dpi=150)
    plt.close(fig)


def _fig05_cross_db_agreement(cfg, summary, fig_dir):
    external_db = cfg.external_cross_check_dbs[0]
    summary = summary.copy()
    summary["dom_corr"] = summary[[f"edlib_{external_db}_frac_strain1", f"edlib_{external_db}_frac_strain2"]].max(axis=1)
    summary["dom_primary"] = summary[[f"edlib_{cfg.primary_db}_frac_strain1", f"edlib_{cfg.primary_db}_frac_strain2"]].max(axis=1)
    s1_reliable = summary["corrdb_strain1_reliable"] == True  # noqa: E712
    s2_reliable = summary["corrdb_strain2_reliable"] == True  # noqa: E712
    both_reliable = s1_reliable & ((summary["well_type"] == "mono") | s2_reliable)
    sub = summary[both_reliable].dropna(subset=["dom_corr", "dom_primary"])

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    for wt, color, label in [("pair", COLOR_PAIR, "pair"), ("mono", COLOR_MONO, "mono")]:
        s = sub[sub.well_type == wt]
        ax.scatter(s["dom_corr"], s["dom_primary"], s=14, color=color, alpha=0.6, edgecolor="none", label=f"{label} (n={len(s)})")
    ax.plot([0, 1], [0, 1], color=COLOR_TEXT_SECONDARY, lw=1, ls="--")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(f"dominant-strain read fraction ({external_db})")
    ax.set_ylabel(f"dominant-strain read fraction ({cfg.primary_db})")
    ax.set_title("Independent vs. self-derived reference", pad=10)
    fig.text(0.55, 0.925, f"{sub['sample_id'].nunique()} wells where both expected strains are name-reliable in {external_db}",
              ha="center", fontsize=8.5, color=COLOR_TEXT_SECONDARY)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(fig_dir / "05_cross_db_agreement.png", dpi=150)
    plt.close(fig)


def _fig06_contamination_offenders(cfg, contam, fig_dir):
    genuine = contam[~contam["likely_reference_twin"]]
    if len(genuine) == 0:
        return
    top = genuine[f"mm2_{cfg.primary_db}_top_other_strain"].value_counts().head(15).sort_values()
    fig, ax = plt.subplots(figsize=(7, max(3, 0.35 * len(top))))
    cmap = plt.get_cmap("Blues")
    norm_vals = top.values / top.values.max()
    colors = [cmap(0.35 + 0.55 * v) for v in norm_vals]
    ax.barh(range(len(top)), top.values, color=colors)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top.index, fontsize=9)
    ax.set_xlabel("# wells where this strain is the top unexpected best-hit")
    ax.set_title("Most frequent unexpected ('contaminant') strains", pad=10)
    fig.text(0.5, 0.94, f"minimap2 best hit vs. {cfg.primary_db}; excludes {len(contam) - len(genuine)}/{len(contam)} reference near-twin false positives",
              ha="center", fontsize=8, color=COLOR_TEXT_SECONDARY)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(fig_dir / "06_contamination_offenders.png", dpi=150)
    plt.close(fig)


def make_all_figures(cfg):
    out_dir, fig_dir = cfg.mapping_validation_out_dir, cfg.mapping_validation_fig_dir
    summary = pd.read_csv(out_dir / "05_combined_sample_summary.csv")
    contam = pd.read_csv(out_dir / "05_contamination_candidates.csv")
    reads_primary = pd.read_csv(out_dir / f"03_edlib_read_assignments_{cfg.primary_db}.csv.gz")
    cross_check = pd.read_csv(out_dir / "02_reference_cross_check.csv")

    _fig00_reference_reliability(cfg, cross_check, fig_dir)
    _fig01_match_quality(cfg, reads_primary, summary, fig_dir)
    _fig02_qc_status_by_well_type(cfg, summary, fig_dir)
    _fig03_purity_histogram(cfg, summary, fig_dir)
    _fig04_depth_vs_confidence(cfg, summary, fig_dir)
    _fig05_cross_db_agreement(cfg, summary, fig_dir)
    _fig06_contamination_offenders(cfg, contam, fig_dir)

    print(f"saved 7 figures -> {fig_dir}")
