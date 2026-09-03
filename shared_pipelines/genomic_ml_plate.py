"""Does cheap plate-reader optics add anything to the genome when predicting who wins a well?

Companion to `genomic_ml.py`, which is the baseline this module is written to be compared
against: same experiment, same labels (`mean_log2_ratio_a_over_b` from
`r03_pair_replicate_stats.csv`), same two CV regimes, same models, same folds. The only thing
that changes between arms is *which features the model may see*. Everything else is held fixed
so that a per-fold difference is attributable to the features and nothing else.

The plate-reader data is a 61-wavelength absorbance scan of every destination well
(`ODFull_10`, 350-950 nm, 30 plates), plus the 384-well source-collection reads that preceded
the experiment. Three things can be pulled out of it, and they are NOT interchangeable --
they differ in *what you must already have measured* to use them, which is the whole point of
splitting them into tiers:

  precult   per-strain OD of the source/preculture plate wells. Free: it was read before the
            experiment existed. (Empirically worthless here -- see below. Kept because a null
            result on a free feature is worth having on the record.)
  mono      per-strain readout of that strain's own monoculture well in the assay condition:
            OD600 and the shape of its 61-point spectrum (pigment + scattering fingerprint).
            Needs one well per strain, no co-culture and no sequencing.
  cocult    the co-culture well's own OD and spectrum, and what they say when compared against
            the two monoculture spectra: a non-negative least-squares unmixing of the pair
            spectrum into its two monoculture components, and which monoculture the pair's OD
            sits closer to. Needs the pair to have been grown -- but not sequenced.

So the tiers answer different questions, in increasing order of what they cost to answer:

  T0_genomic              "I have two genomes."                        <- the baseline
  T1_genomic_precult      "...and the inoculum ODs I already had."
  T2a_genomic_mono_od     "...and how dense each strain grows alone."
  T2b_genomic_mono_spec   "...and the SHAPE of each strain's spectrum alone."
  T2_genomic_mono         "...and both."
  T3_genomic_mono_cocult  "...and I grew the pair; can I skip sequencing it?"
  T4_genomic_all_plate    everything the plate reader has.
  P2_mono_only            plate only, no genome -- how much of any gain is even genomic?
  P3_mono_cocult_only     plate only, including the co-culture well.

T2a/T2b exist because the answer for 20260630 turns out to be counter-intuitive: the
monoculture *spectrum* carries the gain (paired dR^2 +0.12 to +0.19 under cv_strain) and the
monoculture OD carries essentially none (+0.01 to +0.04). A single combined mono tier could
not have shown that, and would have been read as "growth predicts competition".

Two traps this module is built around:

1.  **T3 is not a prediction in the same sense as T0-T2.** Its features are measured on the
    very well whose sequencing produced the label. That is not label leakage (absorbance and
    16S read counts are independent measurement channels), and "replace sequencing with a
    plate reader" is a genuinely useful result -- but it is a *cheaper assay*, not a
    prediction from first principles, and it must never be quoted as "genomes predict
    interactions". T0-T2 are the rows for that claim.

2.  **Adding 90 numbers to a 74-strain problem can help by accident.** Every tier therefore
    also runs with the strain -> plate-well assignment permuted (`PLATESHUFFLED_*`): same
    feature count, same marginals, correspondence destroyed. A gain that survives that is
    real. Tiers are additionally compared *fold-paired* (identical folds and seeds across
    arms) so the delta has a per-fold standard error rather than being two independent
    numbers with overlapping error bars.

3.  **A plate feature is only as good as the well it was read from.** The mono wells are the
    soft spot: 20260630 has one mono well per strain and no replicate, and the user's position
    is that this experiment has no *verified* true monoculture wells. mapping_validation sorts
    them into confirmed / low_confidence / not-assessed (25 / 4 / 45 of the 74 modeled strains
    here) -- and only `low_confidence` is evidence AGAINST a well; `not_assessed` just means
    the well had too few reads to check. `mono_qc_sensitivity()` refits without the suspect
    strains, which is the arm that matters, and separately on the confirmed-only subset, which
    is too small to hold strains out and is therefore scored under cv_pair only.

The source-plate feature gets its own control (`source_plate_specificity_control`) because it
turns out to be the strongest single feature in the univariate screen and the obvious
objection -- that it is plate geometry, not strain identity -- has to be closed before it is
believed. Note the offsets there are NOT a valid null on this collection; read that function.

Call order: build_plate_features -> attach_to_pairs -> univariate_feature_screen ->
source_plate_specificity_control -> spectral_band_profile -> cross_validate_tiers ->
paired_tier_comparison -> mono_qc_sensitivity -> make_all_figures.
"""

import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import spearmanr, pearsonr, wilcoxon
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import genomic_ml as gm
from genomic_ml import (COLOR_BLUE, COLOR_RED, COLOR_GRID, COLOR_TEXT_SECONDARY, COLOR_GOOD,
                        COLOR_CRITICAL, _scaled, _metrics, make_folds, summarize_folds)

OD600_INDEX = 25          # 61 channels from 350 nm in 10 nm steps -> index 25 is 600 nm
SPECTRUM_LO = 5           # channels below 400 nm are dominated by medium absorbance, not cells


# ===========================================================================
# 00 -- configuration
# ===========================================================================

# The 384-well source-collection reads that precede the experiment. Only max_OD is used:
# delta_OD and max growth rate fail even a same-plate positive control here (rho -0.17 / -0.01,
# see plate_reader_comparison/outputs/g01_growth_correlations.csv) because most of these runs
# start at or near stationary.
SOURCE_READS_20260630 = {
    "precult_4x96_30Jun": ("assemble4", ["Karl_20260630_111648_OD600_6m.csv",
                                         "Karl_20260630_113201_OD600_6m.csv",
                                         "Karl_20260630_114305_OD600_6m.csv",
                                         "Karl_20260630_115357_OD600_6m.csv"]),
    "scrapped_gluc_01Jul": ("maxod", ["Karl_20260701_171215_OD600_48h_384_shake.csv"]),
    "NM_recovery_01Jul": ("maxod", ["Karl_20260702_111117_OD600_48h_384_shake.csv"]),
    "NM_dense_02Jul": ("maxod", ["Karl_20260702_182146_OD600_48h_384_shake.csv"]),
    "NM_recovery2_02Jul": ("maxod", ["Karl_20260703_112451_OD600_48h_384_shake.csv"]),
}


@dataclass
class PlateMLConfig:
    """Paths and knobs for the plate-reader-feature analysis of one experiment.

    `gml` is the experiment's GenomicMLConfig -- the baseline arm is run through exactly that
    object, so the T0 rows produced here are the published baseline recomputed on these folds,
    not a reimplementation of it.
    """
    gml: object                       # genomic_ml.GenomicMLConfig -- supplies labels + KO features
    layout_csv: Path                  # strain layout (dest_plate/dest_row/dest_col -> strain1/2)
    od_full_dir: Path                 # folder of 61-channel ODFull csvs, one per destination plate
    source_od_dir: Path = None        # folder of 384-well source-collection reads (optional)
    source_reads: dict = field(default_factory=lambda: dict(SOURCE_READS_20260630))
    mapping_validation_csv: Path = None   # 05_combined_sample_summary.csv, for mono-well QC

    wavelength_index: int = OD600_INDEX
    n_spec_pca: int = 6               # spectral shape is smooth; 6 PCs already saturate it
    detrend_plate: bool = True        # remove additive row/column gradients before z-scoring
    out_dir: Path = None

    def __post_init__(self):
        if self.out_dir is None:
            self.out_dir = self.gml.exp_cfg.exp_base / "analysis" / "genomic_ml_plate" / "outputs"
        self.out_dir = Path(self.out_dir)
        self.fig_dir = self.out_dir / "figures"
        for d in (self.out_dir, self.fig_dir):
            d.mkdir(parents=True, exist_ok=True)
        if self.mapping_validation_csv is None:
            p = self.gml.exp_cfg.mapping_validation_out_dir / "05_combined_sample_summary.csv"
            self.mapping_validation_csv = p if p.exists() else None


# ===========================================================================
# 01 -- reading the plate reader
# ===========================================================================

def read_spectrum_plate(path):
    """One ODFull csv -> (plate id, wavelengths, cube[channel, row, col]).

    The reader writes each channel as its own 16x24 block under a `Wavelength:` header, so the
    'cycles' of these files are WAVELENGTHS, not timepoints -- reading them as a time series is
    the standard way to get nonsense out of this format.
    """
    txt = Path(path).read_text(encoding="latin-1").splitlines()
    plate = next((l.split("ID1:")[1].split("ID2")[0].strip() for l in txt[:8] if "ID1:" in l), None)
    wl = [float(l.split(":")[1].replace("nm", "").strip()) for l in txt if l.startswith("Wavelength:")]
    rows = [[float(x) for x in l.split(",")] for l in txt if re.match(r"^\s*[\d.\-]+\s*,", l)]
    rows = [r for r in rows if len(r) == 24]
    a = np.array(rows, dtype=float)
    n = a.shape[0] // 16
    return plate, np.array(wl), a[: n * 16].reshape(n, 16, 24)


def _read_od_cycles(path):
    """A plain OD csv -> stack[cycle, row, col] (8x12 or 16x24 depending on plate format)."""
    txt = Path(path).read_text(encoding="latin-1").splitlines()
    rows = [[float(x) for x in l.split(",")] for l in txt if re.match(r"^\s*[\d.\-]+\s*,", l)]
    g = np.array(rows, dtype=float)
    h = 8 if g.shape[1] == 12 else 16
    return np.stack([g[i * h:(i + 1) * h] for i in range(g.shape[0] // h)])


def _assemble_384(quads):
    """Four 96-well quadrant reads -> one 384 grid (A1/A2/B1/B2 checkerboard interleave)."""
    out = np.full((16, 24), np.nan)
    for g, (ro, co) in zip(quads, [(0, 0), (0, 1), (1, 0), (1, 1)]):
        out[ro::2, co::2] = g
    return out


def _detrend(grid):
    """Remove additive row and column gradients (evaporation, reader/edge effects)."""
    with warnings.catch_warnings():          # an all-NaN row is legitimate for a partial read
        warnings.simplefilter("ignore", RuntimeWarning)
        return (grid - np.nanmean(grid, 1, keepdims=True)
                - np.nanmean(grid, 0, keepdims=True) + np.nanmean(grid))


def load_well_table(cfg):
    """Every layout well with its blank-corrected OD600 and its 61-channel spectrum.

    Two corrections, both per destination plate:
      - **blank**: the median spectrum of the plate's *unoccupied* wells (76 of 384 per plate
        here). Subtracting it removes medium absorbance and the reader's own offset, which are
        wavelength-dependent and would otherwise dominate the low end of every spectrum.
      - **row/column detrend + z-score** for the scalar OD600, so plate identity and plate
        geometry cannot enter the features. Spectra are left on the blank-corrected scale
        because their *shape* is what is used, and shape is already scale-free.
    """
    cubes, wls = {}, None
    for f in sorted(Path(cfg.od_full_dir).glob("*.csv")):
        pid, wl, cube = read_spectrum_plate(f)
        if pid is None or not str(pid).strip().isdigit():
            continue
        cubes.setdefault(int(pid), []).append(cube)
        wls = wl
    cubes = {k: np.nanmean(v, axis=0) for k, v in cubes.items()}

    lay = pd.read_csv(cfg.layout_csv)
    lay["row_i"] = [ord(str(r).upper()[0]) - 65 for r in lay["dest_row"]]
    lay["col_i"] = lay["dest_col"].astype(int) - 1

    occupied = {(int(p), r, c) for p, r, c in zip(lay.dest_plate, lay.row_i, lay.col_i)}
    blanks, spectra, od = {}, [], []
    for p, cube in cubes.items():
        m = np.ones((16, 24), bool)
        for (pp, r, c) in occupied:
            if pp == p:
                m[r, c] = False
        blanks[p] = np.nanmedian(cube[:, m], axis=1) if m.any() else np.zeros(cube.shape[0])

    for t in lay.itertuples():
        cube = cubes.get(int(t.dest_plate))
        if cube is None or not (0 <= t.row_i < 16 and 0 <= t.col_i < 24):
            spectra.append(np.full(len(wls), np.nan)); od.append(np.nan); continue
        s = cube[:, t.row_i, t.col_i] - blanks[int(t.dest_plate)]
        spectra.append(s)
        od.append(s[cfg.wavelength_index])
    lay["od"] = od
    spectra = np.array(spectra)

    if cfg.detrend_plate:
        dt = np.full(len(lay), np.nan)
        for p, sub in lay.groupby("dest_plate"):
            g = np.full((16, 24), np.nan)
            g[sub.row_i.values, sub.col_i.values] = sub.od.values
            g = _detrend(g)
            dt[sub.index.values] = g[sub.row_i.values, sub.col_i.values]
        lay["od_dt"] = dt
    else:
        lay["od_dt"] = lay["od"]
    lay["od_z"] = lay.groupby("dest_plate")["od_dt"].transform(lambda s: (s - s.mean()) / s.std(ddof=0))
    return lay, spectra, wls, len(cubes)


def _source_grids(cfg):
    """name -> detrended, z-scored 16x24 grid, one per pre-experiment 384-collection read."""
    if cfg.source_od_dir is None:
        return {}
    D = Path(cfg.source_od_dir)
    out = {}
    for name, (how, files) in cfg.source_reads.items():
        paths = [D / f for f in files]
        if not all(p.exists() for p in paths):
            continue
        if how == "assemble4":
            grid = _assemble_384([_read_od_cycles(p)[0] for p in paths])
        else:
            grid = np.nanmax(_read_od_cycles(paths[0]), axis=0)
        g = _detrend(grid)
        out[name] = (g - np.nanmean(g)) / np.nanstd(g)
    return out


def load_source_plate_od(cfg):
    """Per-source-well max OD for each pre-experiment 384-collection read.

    Indexed by source-plate well coordinate, which for 20260630 *is* this experiment's strain
    label namespace (see genomic_ml.validate_strain_join). Returned detrended and z-scored so
    it is on the same footing as the destination features.
    """
    out = _source_grids(cfg)
    if not out:
        return pd.DataFrame()
    wells = [f"{chr(65 + r)}{c + 1}" for r in range(16) for c in range(24)]
    return pd.DataFrame({k: v.ravel() for k, v in out.items()}, index=wells)


# ===========================================================================
# 02 -- per-strain and per-pair plate features
# ===========================================================================

def build_plate_features(cfg, strains):
    """All plate-derived features for one experiment, split by the tier that may use them.

    Returns a dict with:
      mono_scalar  strains x {od_z, spectral brightness/slope summaries}
      mono_spec    strains x 61, each strain's monoculture spectrum normalised to its own OD600
                   (so it carries SHAPE -- pigmentation and scattering -- not biomass, which is
                   already in mono_scalar and would otherwise be counted twice)
      precult      strains x source-plate reads
      wells/spectra/... the raw tables, for the pair-level features and for QC
    """
    lay, spectra, wls, n_plates = load_well_table(cfg)
    lay = lay.reset_index(drop=True)

    # `lay` was reset above, so a row's positional index is its row in `spectra`.
    mono_rows = lay[lay.well_type == "mono"].drop_duplicates("strain1")
    spec_of = {s: spectra[i] for s, i in zip(mono_rows["strain1"], mono_rows.index)}
    mono = mono_rows.set_index("strain1")

    have = [s for s in strains if s in mono.index]
    sl = slice(SPECTRUM_LO, len(wls))

    def _shape(s):
        v = spec_of[s]
        return v / max(v[cfg.wavelength_index], 1e-3)

    mono_spec = pd.DataFrame({s: _shape(s)[sl] for s in have}).T
    mono_spec.columns = [f"wl{int(w)}" for w in wls[sl]]

    # A few interpretable scalars alongside the spectral shape: how much biomass, how strongly
    # pigmented (visible-band absorbance above the 600 nm scattering baseline), and how steep
    # the scattering slope is (cell size / refractive index proxy).
    vis = (wls >= 400) & (wls <= 550)
    nir = wls >= 800
    mono_scalar = pd.DataFrame({
        "mono_od_z": mono.loc[have, "od_z"].values,
        "mono_od_raw": mono.loc[have, "od"].values,
        "mono_pigment": [np.nanmean(_shape(s)[vis]) for s in have],
        "mono_nir_slope": [np.nanmean(_shape(s)[nir]) for s in have],
    }, index=have)

    precult_all = load_source_plate_od(cfg)
    precult = (precult_all.reindex(have) if len(precult_all)
               else pd.DataFrame(index=have))

    qc = pd.DataFrame(index=have)
    qc["has_mono_well"] = True
    if cfg.mapping_validation_csv is not None and Path(cfg.mapping_validation_csv).exists():
        mv = pd.read_csv(cfg.mapping_validation_csv)
        mv = mv[mv.well_type == "mono"].drop_duplicates("strain1").set_index("strain1")
        # Three states, and conflating them would misstate the risk: a strain absent from the
        # mapping_validation summary was never ASSESSED (its mono well fell below the read
        # cutoff), which is not the same as its mono well being shown to hold something else.
        # Only `mono_low_confidence` is positive evidence against the well.
        qc["mono_qc_status"] = mv.reindex(have)["qc_status"].fillna("not_assessed")
        qc["mono_confirmed"] = qc["mono_qc_status"].eq("mono_confirmed")
        qc["mono_suspect"] = qc["mono_qc_status"].eq("mono_low_confidence")
    else:
        qc["mono_qc_status"] = "not_assessed"
        qc["mono_confirmed"] = np.nan
        qc["mono_suspect"] = False

    pf = {"wells": lay, "spectra": spectra, "wavelengths": wls, "n_plates": n_plates,
          "mono_spec": mono_spec, "mono_scalar": mono_scalar, "precult": precult,
          "spec_of": spec_of, "qc": qc, "missing_mono": sorted(set(strains) - set(have))}
    cfg._pf = pf          # every downstream helper reads the bundle off the config
    return pf


def pair_plate_features(cfg, pf, pairs):
    """Per-pair co-culture readouts, split into symmetric and antisymmetric columns.

    Symmetric (unchanged when a and b swap):
      `cocult_od_z`         the pair well's own biomass
      `cocult_yield_dev`    that biomass minus the mean of the two monocultures -- the classic
                            over/under-yielding score
      `cocult_unmix_resid`  how badly the pair spectrum fails to be a mixture of the two
                            monoculture spectra; large means something else is in the well, or
                            one strain looks different in co-culture

    Antisymmetric (negates when a and b swap -- these are the ones that can say WHO won):
      `cocult_unmix_log2`   log2(w_a / w_b) from a non-negative least-squares fit of the pair
                            spectrum onto the two monoculture spectra. This is an optical
                            estimate of the same quantity the sequencing measures, from a
                            completely independent physical channel.
      `cocult_od_closer_a`  |OD_pair - OD_mono_b| - |OD_pair - OD_mono_a|: a well that ends up
                            looking like a's monoculture is a well a took over.
      `cocult_spec_closer_a` the same argument on spectral shape (cosine distance), which does
                            not care about biomass at all.

    Averaged over the pair's replicate wells first; those sit on different destination plates,
    so the average already integrates out plate position.
    """
    lay = pf["wells"]
    spectra = pf["spectra"]
    sl = slice(SPECTRUM_LO, spectra.shape[1])
    wi = cfg.wavelength_index
    spec_of = pf["spec_of"]

    pw = lay[lay.well_type == "pair"].copy()
    pw["key"] = [frozenset((a, b)) for a, b in zip(pw.strain1, pw.strain2)]
    pw = pw[pw.key.map(len) == 2]

    od_z, od_raw, n_wells, mean_spec = {}, {}, {}, {}
    for k, sub in pw.groupby("key"):
        od_z[k] = float(np.nanmean(sub.od_z.values))
        od_raw[k] = float(np.nanmean(sub.od.values))
        n_wells[k] = len(sub)
        mean_spec[k] = np.nanmean(spectra[sub.index.values], axis=0)

    monoz = pf["mono_scalar"]["mono_od_z"].to_dict()
    mono_raw = lay[lay.well_type == "mono"].drop_duplicates("strain1").set_index("strain1")["od"].to_dict()

    def _cos(u, v):
        nu, nv = np.linalg.norm(u), np.linalg.norm(v)
        return 1.0 - float(u @ v) / (nu * nv) if nu > 0 and nv > 0 else np.nan

    rows = []
    for t in pairs.itertuples():
        a, b, k = t.strain_a, t.strain_b, frozenset((t.strain_a, t.strain_b))
        r = {"pair_id": t.pair_id, "cocult_n_wells": n_wells.get(k, 0)}
        ok = k in mean_spec and a in spec_of and b in spec_of
        if ok:
            v = mean_spec[k]
            A, B = spec_of[a], spec_of[b]
            w, _ = nnls(np.column_stack([A[sl], B[sl]]), v[sl])
            resid = float(np.linalg.norm(v[sl] - np.column_stack([A[sl], B[sl]]) @ w)
                          / max(np.linalg.norm(v[sl]), 1e-6))
            sh_v, sh_a, sh_b = (v / max(v[wi], 1e-3), A / max(A[wi], 1e-3), B / max(B[wi], 1e-3))
            r.update({
                "cocult_od_z": od_z[k],
                "cocult_yield_dev": od_z[k] - 0.5 * (monoz.get(a, np.nan) + monoz.get(b, np.nan)),
                "cocult_unmix_resid": resid,
                "cocult_unmix_log2": float(np.log2(max(w[0], 1e-3) / max(w[1], 1e-3))),
                "cocult_unmix_frac_a": float(w[0] / max(w.sum(), 1e-9)) - 0.5,
                "cocult_od_closer_a": float(abs(od_raw[k] - mono_raw.get(b, np.nan))
                                            - abs(od_raw[k] - mono_raw.get(a, np.nan))),
                "cocult_spec_closer_a": float(_cos(sh_v[sl], sh_b[sl]) - _cos(sh_v[sl], sh_a[sl])),
            })
        rows.append(r)
    out = pd.DataFrame(rows).set_index("pair_id")
    return out


SYM_PAIR_COLS = ["cocult_od_z", "cocult_yield_dev", "cocult_unmix_resid"]
ANTI_PAIR_COLS = ["cocult_unmix_log2", "cocult_unmix_frac_a", "cocult_od_closer_a",
                  "cocult_spec_closer_a"]


def attach_to_pairs(cfg, pairs, pf):
    """Join the per-pair plate features onto the modeling table and report what was lost.

    Pairs with no plate-reader coverage (no co-culture well read, or a strain with no mono
    well) are dropped from ALL arms including the baseline -- otherwise the tiers would be
    scored on different pairs and the comparison would be meaningless.
    """
    ppf = pair_plate_features(cfg, pf, pairs)
    p = pairs.set_index("pair_id").join(ppf).reset_index()
    steps = [("modeling_pairs_in", len(pairs))]
    need = SYM_PAIR_COLS + ANTI_PAIR_COLS
    keep = p[need].notna().all(axis=1) & p.strain_a.isin(pf["mono_scalar"].index) \
        & p.strain_b.isin(pf["mono_scalar"].index)
    dropped = p[~keep]
    p = p[keep].reset_index(drop=True)
    steps += [("dropped_no_plate_coverage", len(dropped)), ("pairs_final", len(p))]
    strains = sorted(set(p.strain_a) | set(p.strain_b))
    steps += [("n_strains_final", len(strains)),
              ("n_strains_missing_mono_well", len(pf["missing_mono"])),
              ("n_destination_plates_read", pf["n_plates"]),
              ("n_mono_wells_confirmed", int(pf["qc"].reindex(strains)["mono_confirmed"].sum(skipna=True))),
              ("n_precult_reads", pf["precult"].shape[1])]
    summary = pd.DataFrame(steps, columns=["step", "value"])
    summary.to_csv(cfg.out_dir / "p01_dataset_summary.csv", index=False)
    p.to_csv(cfg.out_dir / "p01_modeling_pairs_with_plate.csv", index=False)
    pf["qc"].reindex(strains).to_csv(cfg.out_dir / "p01_mono_well_qc.csv")
    return p, summary


# ===========================================================================
# 03 -- feature tiers and the fold context that builds them
# ===========================================================================

# ko: KEGG KO PCs + genome summary scalars + 16S bp distance (exactly genomic_ml's design)
# precult / mono / cocult: as described in the module docstring
# `mono` takes "od" (the two OD600 scalars only), "spectrum" (the 61-channel shape PCs only)
# or True (both) -- the split exists because the gain turns out NOT to come from growth alone,
# and a single combined tier could not have shown that.
TIERS = {
    "T0_genomic":             dict(ko=True,  precult=False, mono=False,      cocult=False),
    "T1_genomic_precult":     dict(ko=True,  precult=True,  mono=False,      cocult=False),
    "T2a_genomic_mono_od":    dict(ko=True,  precult=False, mono="od",       cocult=False),
    "T2b_genomic_mono_spec":  dict(ko=True,  precult=False, mono="spectrum", cocult=False),
    "T2_genomic_mono":        dict(ko=True,  precult=False, mono=True,       cocult=False),
    "T3_genomic_mono_cocult": dict(ko=True,  precult=False, mono=True,       cocult=True),
    "T4_genomic_all_plate":   dict(ko=True,  precult=True,  mono=True,       cocult=True),
    "P2_mono_only":           dict(ko=False, precult=False, mono=True,       cocult=False),
    "P3_mono_cocult_only":    dict(ko=False, precult=False, mono=True,       cocult=True),
}


def _scale_anti(M, train_rows):
    """Scale-only standardisation for antisymmetric columns.

    Centring an antisymmetric feature would destroy the antisymmetry it exists to carry (a
    constant offset does not flip sign when a and b swap), so these get divided by their
    training standard deviation and nothing else.
    """
    sd = np.nanstd(M[train_rows], axis=0)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    return M / sd


class TierCtx(gm._Ctx):
    """genomic_ml's fold context, with a tier-selected feature set.

    Subclassed rather than reimplemented so that every model in `genomic_ml.MODELS` keeps
    working unchanged, and so that the `T0_genomic` arm is byte-for-byte the baseline design
    (the plate blocks are *appended* after genomic_ml's own block layout, never interleaved).
    `self.Z` and `self.S` are rewritten in place for the tier, which is what makes the
    two-stage models -- whose regression is per strain, not per pair -- pick the plate
    features up automatically.
    """

    def __init__(self, cfg, gcfg, pairs, X, summ, train_idx, test_idx, tier, pf,
                 phylo=None, plate_perm=None):
        self.cfg, self.tier, self.spec_ = cfg, tier, TIERS[tier]
        self.pf, self.plate_perm = pf, plate_perm
        super().__init__(gcfg, pairs, X, summ, train_idx, test_idx, phylo=phylo)

    # -- per-strain plate blocks, standardised on the training strains only ----------
    def _strain_plate(self):
        strains = list(self.X.index)
        perm = self.plate_perm or {}
        src = [perm.get(s, s) for s in strains]

        hi, sc = [], []
        if self.spec_["mono"] in ("spectrum", True):
            spec = self.pf["mono_spec"].reindex(src)
            spec.index = strains
            k = min(self.cfg.n_spec_pca, len(self.train_strains) - 1, spec.shape[1])
            p = PCA(n_components=k, random_state=0).fit(spec.loc[self.train_strains].values)
            hi.append(_scaled(pd.DataFrame(p.transform(spec.values), index=strains,
                                           columns=[f"spec_pc{i+1}" for i in range(k)]),
                              self.train_strains))
        if self.spec_["mono"] in ("od", True):
            cols = ["mono_od_z", "mono_od_raw"] if self.spec_["mono"] == "od" else None
            s = self.pf["mono_scalar"].reindex(src)
            s.index = strains
            sc.append(_scaled(s[cols] if cols else s, self.train_strains))
        if self.spec_["precult"] and self.pf["precult"].shape[1]:
            s = self.pf["precult"].reindex(src)
            s.index = strains
            sc.append(_scaled(s.fillna(s.loc[self.train_strains].mean()), self.train_strains))
        return hi, sc

    # -- per-pair plate blocks ------------------------------------------------------
    def _pair_plate(self):
        if not self.spec_["cocult"]:
            n = len(self.pairs)
            return np.zeros((n, 0)), np.zeros((n, 0)), [], []
        sym = self.pairs[SYM_PAIR_COLS].values.astype(float)
        anti = self.pairs[ANTI_PAIR_COLS].values.astype(float)
        sc = StandardScaler().fit(sym[self.train_idx])
        return sc.transform(sym), _scale_anti(anti, self.train_idx), SYM_PAIR_COLS, ANTI_PAIR_COLS

    def design(self, Z):
        """[genomic_ml's blocks if ko] + [plate strain blocks] + [plate pair blocks]."""
        hi, sc = self._strain_plate()
        Zp = pd.concat(hi, axis=1) if hi else pd.DataFrame(index=Z.index)
        Sp = pd.concat(sc, axis=1) if sc else pd.DataFrame(index=Z.index)

        a, b = self.pairs["strain_a"].values, self.pairs["strain_b"].values
        blocks_f, blocks_r, names = [], [], []

        if self.spec_["ko"]:
            Za, Zb = Z.loc[a].values, Z.loc[b].values
            Sa, Sb = self.S.loc[a].values, self.S.loc[b].values
            symsc = np.column_stack([np.asarray(v) for v in gm.pair_feature_blocks(self.pairs, Z, self.S)])
            blocks_f += [Za - Zb, Za + Zb, symsc, Sa - Sb]
            blocks_r += [Zb - Za, Za + Zb, symsc, Sb - Sa]
            names += ([f"anti_{c}" for c in Z.columns] + [f"sym_{c}" for c in Z.columns]
                      + [f"sym_x{i}" for i in range(symsc.shape[1])]
                      + [f"anti_{c}" for c in self.S.columns])

        if Zp.shape[1] or Sp.shape[1]:
            Pa = np.hstack([Zp.loc[a].values, Sp.loc[a].values])
            Pb = np.hstack([Zp.loc[b].values, Sp.loc[b].values])
            cols = list(Zp.columns) + list(Sp.columns)
            blocks_f += [Pa - Pb, Pa + Pb]
            blocks_r += [Pb - Pa, Pa + Pb]
            names += [f"anti_{c}" for c in cols] + [f"sym_{c}" for c in cols]

        psym, panti, sname, aname = self._pair_plate()
        if psym.shape[1] or panti.shape[1]:
            blocks_f += [psym, panti]
            blocks_r += [psym, -panti]
            names += [f"sym_{c}" for c in sname] + [f"anti_{c}" for c in aname]

        # what the two-stage (per-strain) models see: the tier's per-strain representation only
        self.Z = pd.concat([Z, Zp], axis=1) if self.spec_["ko"] else (
            Zp if Zp.shape[1] else pd.DataFrame(np.zeros((len(Z), 1)), index=Z.index, columns=["const"]))
        self.S = pd.concat([self.S, Sp], axis=1) if self.spec_["ko"] else Sp

        fwd = np.hstack(blocks_f) if blocks_f else np.zeros((len(self.pairs), 1))
        rev = np.hstack(blocks_r) if blocks_r else np.zeros((len(self.pairs), 1))
        return fwd, rev, names

    @property
    def raw(self):
        raise RuntimeError("raw-KO models are not part of the tier sweep -- the ~9.5k-column "
                           "design would swamp the handful of plate columns being tested")


# -- one extra model: the plate reader with no machine learning at all ----------------
def m_plate_unmix_only(ctx):
    """Rescale the spectral unmixing log-ratio to the target with a single in-fold slope.

    Not a model so much as a reading of the instrument: it asks what the plate reader alone
    says about who won, with one training-set-fitted number (the slope) and nothing else. If a
    full model on the same features cannot beat this, the machine learning is decoration.
    """
    if not TIERS[ctx.tier]["cocult"]:
        return np.zeros(len(ctx.test_idx))
    x = ctx.pairs["cocult_unmix_log2"].values.astype(float)
    xt, yt = x[ctx.train_idx], ctx.y_train
    ok = np.isfinite(xt) & np.isfinite(yt)
    slope = float(np.sum(xt[ok] * yt[ok]) / max(np.sum(xt[ok] ** 2), 1e-9))   # through the origin
    return slope * np.nan_to_num(x[ctx.test_idx])


TIER_MODELS = ["zero_baseline", "strength_observed_no_genomics", "ridge_pca", "xgboost_pca",
               "two_stage_ridge", "two_stage_xgboost"]

# Models whose prediction cannot change with the tier (they never look at features), so running
# them once under T0 is enough and repeating them per tier would just pad the table.
TIER_INVARIANT = {"zero_baseline", "strength_observed_no_genomics"}


# ===========================================================================
# 04 -- the tier sweep
# ===========================================================================

def cross_validate_tiers(cfg, pairs, X, summ, tiers=None, models=None,
                         regimes=("cv_pair", "cv_strain"), n_splits=5, seed=0, n_repeats=3,
                         plate_shuffle_control=True, file_prefix="p03"):
    """Every tier under every regime, on IDENTICAL folds.

    The folds are generated once per (regime, repeat) and reused across all tiers, so a tier
    difference is a paired difference -- the noise from which strains happened to land in which
    fold is common to both arms and cancels in `paired_tier_comparison`. With ~40 test pairs per
    cv_strain fold that pairing is not a nicety; unpaired, the fold-to-fold spread is larger
    than the effect being measured.
    """
    tiers = tiers or list(TIERS)
    models = models or TIER_MODELS
    gcfg = cfg.gml
    pf = cfg._pf
    per_fold, per_pair = [], []

    rng = np.random.default_rng(4242)
    strains_all = list(X.index)

    for regime in regimes:
        folds = [(rep, tr, te) for rep in range(n_repeats)
                 for tr, te in make_folds(pairs, regime, n_splits=n_splits, seed=seed + 97 * rep)]
        for fi, (rep, tr, te) in enumerate(folds):
            for tier in tiers:
                ctx = TierCtx(cfg, gcfg, pairs, X, summ, tr, te, tier, pf)
                run = [m for m in models if tier == tiers[0] or m not in TIER_INVARIANT]
                if TIERS[tier]["cocult"]:
                    run = run + ["plate_unmix_only"]
                for name in run:
                    fn = m_plate_unmix_only if name == "plate_unmix_only" else gm.MODELS[name]
                    try:
                        yhat = fn(ctx)
                    except Exception as e:
                        print(f"  ! {tier}/{name} failed on {regime} fold {fi}: "
                              f"{type(e).__name__}: {e}")
                        continue
                    per_fold.append({"regime": regime, "fold": fi, "tier": tier, "model": name,
                                     "n_train": len(tr), **_metrics(ctx.y_test, yhat)})
                    per_pair.append(pd.DataFrame({
                        "regime": regime, "fold": fi, "tier": tier, "model": name,
                        "pair_id": pairs.iloc[te]["pair_id"].values,
                        "y_true": ctx.y_test, "y_pred": yhat}))

                # negative control: same features, wrong strains. Only for tiers that actually
                # add plate columns -- under T0 it would be a no-op.
                if plate_shuffle_control and tier != "T0_genomic":
                    perm = dict(zip(strains_all, rng.permutation(strains_all)))
                    ctx_s = TierCtx(cfg, gcfg, pairs, X, summ, tr, te, tier, pf, plate_perm=perm)
                    for name in ("ridge_pca", "two_stage_ridge"):
                        if name not in models:
                            continue
                        try:
                            yhat = gm.MODELS[name](ctx_s)
                        except Exception:
                            continue
                        per_fold.append({"regime": regime, "fold": fi, "tier": tier,
                                         "model": f"PLATESHUFFLED_{name}", "n_train": len(tr),
                                         **_metrics(ctx.y_test, yhat)})

    fold_df = pd.DataFrame(per_fold)
    pair_df = pd.concat(per_pair, ignore_index=True) if per_pair else pd.DataFrame()
    # A two-stage model regresses per STRAIN, so per-PAIR co-culture columns are invisible to
    # it: T2/T3 and P2/P3 must come out identical for those models. Flagged rather than hidden,
    # since the equality is itself the evidence that the pair features are wired where intended.
    fold_df["model_can_use_pair_features"] = ~fold_df.model.str.contains("two_stage")
    num = ["pearson_r", "spearman_rho", "r2", "mae", "sign_accuracy"]
    summ_df = (fold_df.groupby(["regime", "tier", "model"])
               .agg(n_folds=("fold", "nunique"), n_test=("n", "sum"),
                    can_use_pair_features=("model_can_use_pair_features", "first"),
                    **{f"{c}_mean": (c, "mean") for c in num},
                    **{f"{c}_sd": (c, "std") for c in num})
               .reset_index().sort_values(["regime", "r2_mean"], ascending=[True, False]))

    fold_df.to_csv(cfg.out_dir / f"{file_prefix}_cv_per_fold.csv", index=False)
    summ_df.to_csv(cfg.out_dir / f"{file_prefix}_cv_summary.csv", index=False)
    pair_df.to_csv(cfg.out_dir / f"{file_prefix}_cv_predictions.csv", index=False)
    return summ_df, fold_df, pair_df


def paired_tier_comparison(cfg, fold_df, baseline="T0_genomic", file_prefix="p04"):
    """Fold-paired delta of every tier against the baseline tier, per regime and model.

    Reported as mean delta with the standard error OF THE DELTA (not of either arm), plus a
    Wilcoxon signed-rank test over folds. The folds are shared across tiers by construction, so
    this is the comparison the sweep was designed to support; reading two rows of
    `p03_cv_summary.csv` against each other instead would overstate the uncertainty several
    fold.
    """
    rows = []
    d = fold_df[~fold_df.model.str.startswith(("PLATESHUFFLED", "SHUFFLED"))]
    for (regime, model), sub in d.groupby(["regime", "model"]):
        base = sub[sub.tier == baseline].set_index("fold")
        if not len(base):
            continue
        for tier, t in sub.groupby("tier"):
            if tier == baseline:
                continue
            t = t.set_index("fold")
            common = base.index.intersection(t.index)
            if len(common) < 3:
                continue
            for metric in ("r2", "spearman_rho", "sign_accuracy", "mae"):
                delta = (t.loc[common, metric] - base.loc[common, metric]).dropna()
                if len(delta) < 3:
                    continue
                try:
                    p = float(wilcoxon(delta)[1])
                except ValueError:            # all-zero differences
                    p = 1.0
                rows.append({"regime": regime, "model": model, "tier": tier, "metric": metric,
                             "n_folds": len(delta), "baseline_mean": float(base.loc[common, metric].mean()),
                             "tier_mean": float(t.loc[common, metric].mean()),
                             "delta_mean": float(delta.mean()),
                             "delta_se": float(delta.std(ddof=1) / np.sqrt(len(delta))),
                             "wilcoxon_p": p})
    out = pd.DataFrame(rows).sort_values(["regime", "metric", "delta_mean"],
                                         ascending=[True, True, False])
    out.to_csv(cfg.out_dir / f"{file_prefix}_paired_vs_baseline.csv", index=False)
    return out


def univariate_feature_screen(cfg, pairs, file_prefix="p02"):
    """Each plate feature against the target on its own, before any model touches it.

    Antisymmetric features are correlated with the signed target; symmetric ones with |target|
    (a symmetric feature cannot say who won, only how lopsided the well was). Cheap, and it
    makes the model results interpretable rather than a black box that got better.
    """
    y = pairs[cfg.gml.target].values
    rows = []
    for c in ANTI_PAIR_COLS:
        v = pairs[c].values.astype(float)
        ok = np.isfinite(v) & np.isfinite(y)
        dec = ok & (np.abs(y) > 1.0)
        rows.append({"feature": c, "kind": "antisymmetric", "vs": "log2_ratio", "n": int(ok.sum()),
                     "pearson_r": float(pearsonr(y[ok], v[ok])[0]),
                     "spearman_rho": float(spearmanr(y[ok], v[ok])[0]),
                     "sign_accuracy": float(np.mean(np.sign(y[dec]) == np.sign(v[dec])))})
    for c in SYM_PAIR_COLS:
        v = pairs[c].values.astype(float)
        ok = np.isfinite(v) & np.isfinite(y)
        rows.append({"feature": c, "kind": "symmetric", "vs": "|log2_ratio|", "n": int(ok.sum()),
                     "pearson_r": float(pearsonr(np.abs(y[ok]), v[ok])[0]),
                     "spearman_rho": float(spearmanr(np.abs(y[ok]), v[ok])[0]),
                     "sign_accuracy": np.nan})
    pf = cfg._pf
    for c in list(pf["mono_scalar"].columns) + list(pf["precult"].columns):
        src = pf["mono_scalar"] if c in pf["mono_scalar"].columns else pf["precult"]
        d = pairs["strain_a"].map(src[c]) - pairs["strain_b"].map(src[c])
        ok = np.isfinite(d.values) & np.isfinite(y)
        dec = ok & (np.abs(y) > 1.0)
        rows.append({"feature": f"{c}__a_minus_b", "kind": "antisymmetric", "vs": "log2_ratio",
                     "n": int(ok.sum()), "pearson_r": float(pearsonr(y[ok], d.values[ok])[0]),
                     "spearman_rho": float(spearmanr(y[ok], d.values[ok])[0]),
                     "sign_accuracy": float(np.mean(np.sign(y[dec]) == np.sign(d.values[dec])))})
    out = pd.DataFrame(rows).sort_values("pearson_r", key=np.abs, ascending=False)
    out.to_csv(cfg.out_dir / f"{file_prefix}_univariate_screen.csv", index=False)
    return out


def source_plate_specificity_control(cfg, pairs, n_perm=2000, seed=0,
                                     offsets=((0, 1), (1, 0), (0, -1), (-1, 0), (1, 1),
                                              (0, 6), (4, 0)), file_prefix="p02b"):
    """Is the source-plate OD -> outcome link about the LAYOUT, or about plate geometry?

    The per-strain source-plate OD is the single strongest plate feature here, and that is a
    claim worth attacking before it is believed: a plate read carries smooth spatial structure
    (edge drying, thermal gradients) that additive row/column detrending does not remove, and
    which would manufacture a correlation between any two reads of any two plates.

    Two nulls, and the difference between them is the point:

    `permutation` (the verdict) reassigns the 74 strains to source wells at random and
        recomputes r. It preserves both marginals -- the OD distribution and the pair structure
        -- and destroys only the strain-to-well correspondence, which is exactly the hypothesis.
        Its spread is wide (the pairs are 74 genomes, not n independent rows), and that width is
        the honest error bar.

    `offset` slides the grid by a well or two. Classically the sharper control -- shared layout
        lives at exactly one alignment -- but it is NOT trustworthy on this plate and the output
        says so per row: `grid_selfcorr` reports how much of the true per-strain OD ordering
        survives each shift. This collection has ~6-column periodicity from the way it was
        consolidated, so a 6-column roll leaves the strains ANTI-correlated with themselves
        (selfcorr ~ -0.26) rather than uncorrelated, and the "shifted" r comes back large and
        negative -- an artifact of the shift being a bad null, not evidence against the signal.
        Read an offset row only where |grid_selfcorr| is near zero.
    """
    y = pairs[cfg.gml.target].values
    strains = sorted(set(pairs["strain_a"]) | set(pairs["strain_b"]))
    pos = {s: i for i, s in enumerate(strains)}
    A = pairs["strain_a"].map(pos).values
    B = pairs["strain_b"].map(pos).values
    rc = [(ord(s[0]) - 65, int(s[1:]) - 1) for s in strains]
    rng = np.random.default_rng(seed)
    rows = []
    for name, g in _source_grids(cfg).items():
        v = np.array([g[i, j] for i, j in rc])
        obs = float(pearsonr(y, v[A] - v[B])[0])
        null = np.array([pearsonr(y, (p := rng.permutation(v))[A] - p[B])[0] for _ in range(n_perm)])
        z = (obs - null.mean()) / null.std()
        r = {"source_read": name, "r_true": obs,
             "perm_null_mean": float(null.mean()), "perm_null_sd": float(null.std()),
             "z_vs_permutation": float(z),
             "p_perm": float((np.abs(null) >= abs(obs)).mean()),
             "layout_specific": bool(z > 3)}
        for dr, dc in offsets:
            gs = np.roll(np.roll(g, dr, 0), dc, 1)
            vs = np.array([gs[i, j] for i, j in rc])
            r[f"r_offset_{dr}_{dc}"] = float(pearsonr(y, vs[A] - vs[B])[0])
            r[f"selfcorr_offset_{dr}_{dc}"] = float(pearsonr(v, vs)[0])
        rows.append(r)
    out = pd.DataFrame(rows).sort_values("z_vs_permutation", ascending=False)
    out.to_csv(cfg.out_dir / f"{file_prefix}_source_plate_specificity.csv", index=False)
    return out


def spectral_band_profile(cfg, pairs, file_prefix="p02c"):
    """Which wavelengths carry the who-wins signal, band by band, with no model fitted.

    The monoculture spectrum turns out to beat the monoculture OD as a predictor, which only
    means something if the useful part is somewhere other than the 600 nm scattering channel
    the OD already is. This resolves that: for each of the 56 usable channels, correlate
    (shape_a - shape_b) with log2(a/b), where shape is the spectrum divided by its own OD600 --
    so channel 600 is 1.0 for every strain by construction and must come out at exactly zero.
    That built-in zero is the internal check that the normalisation did what it claims.

    Also writes the PCA structure of the shape matrix, because the band profile alone invites
    an over-reading. For 20260630 the profile is a clean step through 600 nm -- positive across
    the visible, negative across the near-IR -- and PC1 carries 90% of the shape variance. That
    is the signature of ONE latent scalar, the steepness of the absorbance spectrum (a
    scattering exponent, i.e. cell size / refractive index), not of rich pigment structure. PC2
    (8%) and PC4 do add independent signal, so it is not literally one number, but the honest
    summary is "a scattering-slope phenotype plus a little", not "56 channels of information".
    """
    pf = cfg._pf
    spec = pf["mono_spec"]
    y = pairs[cfg.gml.target].values
    A = spec.reindex(pairs["strain_a"]).values
    B = spec.reindex(pairs["strain_b"]).values
    d = A - B
    rows = []
    for k, col in enumerate(spec.columns):
        v = d[:, k]
        ok = np.isfinite(v) & np.isfinite(y)
        if ok.sum() < 10 or np.std(v[ok]) == 0:
            rows.append({"channel": col, "wavelength_nm": int(col[2:]), "pearson_r": 0.0,
                         "spearman_rho": 0.0})
            continue
        rows.append({"channel": col, "wavelength_nm": int(col[2:]),
                     "pearson_r": float(pearsonr(y[ok], v[ok])[0]),
                     "spearman_rho": float(spearmanr(y[ok], v[ok])[0])})
    out = pd.DataFrame(rows)
    out.to_csv(cfg.out_dir / f"{file_prefix}_spectral_band_profile.csv", index=False)

    k = min(6, spec.shape[0] - 1, spec.shape[1])
    pca = PCA(n_components=k, random_state=0).fit(spec.values)
    Z = pd.DataFrame(pca.transform(spec.values), index=spec.index,
                     columns=[f"pc{i+1}" for i in range(k)])
    pc_rows = []
    for i, c in enumerate(Z.columns):
        dz = (pairs["strain_a"].map(Z[c]) - pairs["strain_b"].map(Z[c])).values
        ok = np.isfinite(dz) & np.isfinite(y)
        pc_rows.append({"component": c,
                        "explained_variance_ratio": float(pca.explained_variance_ratio_[i]),
                        "pearson_r_with_target": float(pearsonr(y[ok], dz[ok])[0])})
    pd.DataFrame(pc_rows).to_csv(cfg.out_dir / f"{file_prefix}_spectral_pc_structure.csv",
                                 index=False)
    return out


def mono_qc_sensitivity(cfg, pairs, X, summ, tiers=("T0_genomic", "T2_genomic_mono"),
                        model="two_stage_ridge", n_repeats=2, file_prefix="p05"):
    """Does the monoculture gain survive dropping the mono wells the sequencing distrusts?

    mapping_validation puts each mono well in one of three states, and the distinction matters
    because only one of them is evidence against the feature:

      mono_confirmed       the reads say the well holds the strain the layout claims
      mono_low_confidence  the reads say it holds something else -- a WRONG monoculture readout
      not_assessed         the well fell below the read-count cutoff, so nothing was checked

    Two arms, because they trade off differently:

      `drop_suspect`    removes only the `mono_low_confidence` strains. Keeps almost every
                        strain, so leave-strains-out stays powered. This is the real test.
      `confirmed_only`  keeps only positively confirmed strains. Far stricter, but it leaves
                        too few strains to hold any out -- five per fold, and a ridge
                        extrapolating off five genomes produces R² in the 1e24 range, which is
                        noise wearing a number. So this arm is scored under `cv_pair` ONLY, and
                        must be read as "does the signal survive at all", not as a
                        generalization estimate.
    """
    qc = cfg._pf["qc"]
    suspect = set(qc.index[qc.get("mono_suspect", pd.Series(False, index=qc.index)).fillna(False)])
    confirmed = set(qc.index[qc["mono_confirmed"].fillna(False)])

    arms = {
        "drop_suspect": (pairs[~pairs.strain_a.isin(suspect) & ~pairs.strain_b.isin(suspect)],
                         ("cv_pair", "cv_strain")),
        "confirmed_only": (pairs[pairs.strain_a.isin(confirmed) & pairs.strain_b.isin(confirmed)],
                           ("cv_pair",)),
    }
    info, cmps = [], []
    counts = qc["mono_qc_status"].value_counts().to_dict()
    for k, v in counts.items():
        info.append({"arm": f"qc_status::{k}", "n_pairs": np.nan, "n_strains": int(v),
                     "regimes": ""})
    info.append({"arm": "all", "n_pairs": len(pairs),
                 "n_strains": len(set(pairs.strain_a) | set(pairs.strain_b)), "regimes": "both"})

    for arm, (sub, regimes) in arms.items():
        sub = sub.reset_index(drop=True)
        strains = sorted(set(sub.strain_a) | set(sub.strain_b))
        info.append({"arm": arm, "n_pairs": len(sub), "n_strains": len(strains),
                     "regimes": "+".join(regimes)})
        if len(sub) < 120 or len(strains) < 12:
            continue
        _, f, _ = cross_validate_tiers(cfg, sub, X.loc[strains], summ.loc[strains],
                                       tiers=list(tiers), models=[model], regimes=regimes,
                                       n_repeats=n_repeats, plate_shuffle_control=False,
                                       file_prefix=f"{file_prefix}_{arm}")
        c = paired_tier_comparison(cfg, f, file_prefix=f"{file_prefix}_{arm}_paired")
        cmps.append(c.assign(arm=arm))
    out = pd.DataFrame(info)
    out.to_csv(cfg.out_dir / f"{file_prefix}_mono_qc_sensitivity.csv", index=False)
    cmp_df = pd.concat(cmps, ignore_index=True) if cmps else pd.DataFrame()
    if len(cmp_df):
        cmp_df.to_csv(cfg.out_dir / f"{file_prefix}_mono_qc_paired.csv", index=False)
    return out, cmp_df


# ===========================================================================
# 05 -- figures
# ===========================================================================

TIER_ORDER = ["T0_genomic", "T1_genomic_precult", "T2a_genomic_mono_od", "T2b_genomic_mono_spec",
              "T2_genomic_mono", "T3_genomic_mono_cocult", "T4_genomic_all_plate",
              "P2_mono_only", "P3_mono_cocult_only"]


def make_all_figures(cfg, summ_df, fold_df, pair_df, screen, paired, band=None):
    fd = cfg.fig_dir

    # -- f1: tier x model, one panel per regime ---------------------------------------
    # Error bars are the SE of the fold mean, not the fold-to-fold SD: under cv_strain the SD
    # is larger than every effect on the plot (~40 test pairs a fold), so drawing it here would
    # hide the result rather than qualify it. The honest between-tier comparison is fp03, which
    # is paired; this panel is for reading levels, not differences.
    for metric, lab, fn in [("r2", "R² (vs. predicting 0)", "fp01_tier_r2.png"),
                            ("spearman_rho", "Spearman ρ", "fp02_tier_rho.png")]:
        fig, axes = plt.subplots(1, 2, figsize=(15, 6.4))
        for ax, regime, title in zip(axes, ["cv_pair", "cv_strain"],
                                     ["cv_pair: unseen PAIR, strains seen",
                                      "cv_strain: unseen STRAINS — the honest test"]):
            d = summ_df[(summ_df.regime == regime) & (~summ_df.model.str.startswith("PLATESHUFFLED"))
                        & (~summ_df.model.isin(TIER_INVARIANT))]
            if not len(d):
                continue
            tiers = [t for t in TIER_ORDER if t in set(d.tier)]
            mods = sorted(set(d.model))
            w = 0.8 / len(mods)
            for i, m in enumerate(mods):
                v, e = [], []
                for t in tiers:
                    row = d[(d.tier == t) & (d.model == m)]
                    v.append(row[f"{metric}_mean"].mean() if len(row) else np.nan)
                    n = row["n_folds"].mean() if len(row) else np.nan
                    e.append(row[f"{metric}_sd"].mean() / np.sqrt(max(n, 1)) if len(row) else np.nan)
                ax.bar(np.arange(len(tiers)) + i * w, v, w, yerr=e, label=m,
                       error_kw=dict(lw=.8, capsize=2, ecolor=COLOR_TEXT_SECONDARY))
            base = d[d.tier == "T0_genomic"][f"{metric}_mean"].max()
            ax.axhline(base, color=COLOR_RED, lw=1.2, ls="--")
            ax.text(len(tiers) - .45, base, "best genomic baseline ", color=COLOR_RED, fontsize=7.5,
                    va="bottom", ha="right")
            ax.axhline(0, color="black", lw=.8)
            ax.set_xticks(np.arange(len(tiers)) + (0.8 - w) / 2)
            ax.set_xticklabels(tiers, fontsize=8, rotation=25, ha="right")
            ax.set_ylabel(lab); ax.set_title(title, fontsize=10)
            ax.grid(axis="y", color=COLOR_GRID, lw=.5)
            ax.set_axisbelow(True)
        axes[0].legend(fontsize=7.5, frameon=False, loc="upper left")
        fig.suptitle("Plate-reader features added to the genomic baseline, same folds throughout",
                     fontweight="bold")
        fig.tight_layout(); fig.savefig(fd / fn, dpi=160); plt.close(fig)

    # -- f3: the paired deltas, which is the actual result ---------------------------
    d = paired[(paired.metric == "r2") & (paired.regime == "cv_strain")]
    if len(d):
        fig, ax = plt.subplots(figsize=(8.5, 0.42 * len(d) + 2))
        d = d.sort_values("delta_mean")
        lab = [f"{r.tier}  ·  {r.model}" for r in d.itertuples()]
        cols = [COLOR_GOOD if r.delta_mean > 0 and r.wilcoxon_p < .05 else
                (COLOR_BLUE if r.delta_mean > 0 else COLOR_CRITICAL) for r in d.itertuples()]
        ax.barh(lab, d.delta_mean, xerr=d.delta_se, color=cols,
                error_kw=dict(ecolor=COLOR_TEXT_SECONDARY, lw=1, capsize=3))
        ax.axvline(0, color="black", lw=1)
        for i, r in enumerate(d.itertuples()):
            ax.text(r.delta_mean + np.sign(r.delta_mean) * (r.delta_se + .004), i,
                    f"p={r.wilcoxon_p:.3f}", va="center", fontsize=7,
                    ha="left" if r.delta_mean > 0 else "right", color=COLOR_TEXT_SECONDARY)
        ax.set_xlabel("Δ R² vs. T0_genomic, paired over folds (mean ± SE of the difference)")
        ax.set_title("cv_strain: what each plate-reader tier buys over the genomic baseline",
                     fontweight="bold", fontsize=11)
        fig.tight_layout(); fig.savefig(fd / "fp03_paired_delta_cv_strain.png", dpi=160)
        plt.close(fig)

    # -- f4: the univariate screen ----------------------------------------------------
    if len(screen):
        s = screen[screen.kind == "antisymmetric"].sort_values("pearson_r")
        fig, ax = plt.subplots(figsize=(8, 0.4 * len(s) + 1.8))
        ax.barh(s.feature, s.pearson_r, color=[COLOR_BLUE if v > 0 else COLOR_RED for v in s.pearson_r])
        ax.axvline(0, color="black", lw=1)
        ax.set_xlabel("Pearson r with log2(a/b), no model fitted")
        ax.set_title("Single plate-reader features against the target", fontweight="bold", fontsize=11)
        fig.tight_layout(); fig.savefig(fd / "fp04_univariate_screen.png", dpi=160); plt.close(fig)

    # -- f4b: where in the spectrum the signal lives ---------------------------------
    if band is not None and len(band):
        fig, ax = plt.subplots(figsize=(8.6, 4.2))
        ax.plot(band.wavelength_nm, band.pearson_r, color=COLOR_BLUE, lw=1.8)
        ax.fill_between(band.wavelength_nm, 0, band.pearson_r, color=COLOR_BLUE, alpha=.15)
        ax.axhline(0, color="black", lw=.8)
        ax.axvline(600, color=COLOR_RED, lw=1, ls="--")
        ax.text(603, ax.get_ylim()[1] * .9, "600 nm — pinned to 0 by the\nnormalisation (internal check)",
                fontsize=7.5, color=COLOR_RED, va="top")
        ax.set_xlabel("wavelength (nm)")
        ax.set_ylabel("Pearson r of (shape$_a$ − shape$_b$) with log2(a/b)")
        ax.set_title("Where in the monoculture spectrum the who-wins signal sits",
                     fontweight="bold", fontsize=11)
        ax.grid(color=COLOR_GRID, lw=.5); ax.set_axisbelow(True)
        fig.tight_layout(); fig.savefig(fd / "fp04b_spectral_band_profile.png", dpi=160)
        plt.close(fig)

    # -- f5: predicted vs observed, baseline vs best plate tier -----------------------
    if len(pair_df):
        d = summ_df[(summ_df.regime == "cv_strain") & (~summ_df.model.str.startswith("PLATESHUFFLED"))
                    & (~summ_df.model.isin(TIER_INVARIANT))]
        picks = []
        for tier in ["T0_genomic", "T2_genomic_mono", "T4_genomic_all_plate"]:
            dd = d[d.tier == tier]
            if len(dd):
                picks.append((tier, dd.sort_values("r2_mean").iloc[-1].model))
        if picks:
            fig, axes = plt.subplots(1, len(picks), figsize=(4.7 * len(picks), 4.8), squeeze=False)
            for ax, (tier, model) in zip(axes[0], picks):
                p = pair_df[(pair_df.regime == "cv_strain") & (pair_df.tier == tier)
                            & (pair_df.model == model)]
                ax.scatter(p.y_pred, p.y_true, s=9, alpha=.3, color=COLOR_BLUE, edgecolor="none")
                lim = np.nanpercentile(np.abs(p.y_true), 99) * 1.1
                ax.plot([-lim, lim], [-lim, lim], color=COLOR_TEXT_SECONDARY, lw=1, ls="--")
                ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
                ax.axhline(0, color=COLOR_GRID, lw=.8); ax.axvline(0, color=COLOR_GRID, lw=.8)
                ax.set_title(f"{tier}\n{model} — ρ = {spearmanr(p.y_pred, p.y_true)[0]:.2f}",
                             fontsize=9)
                ax.set_xlabel("predicted log2(a/b)"); ax.set_ylabel("observed log2(a/b)")
            fig.suptitle("cv_strain (both strains unseen): baseline vs. plate-reader tiers",
                         fontweight="bold")
            fig.tight_layout(); fig.savefig(fd / "fp05_predicted_vs_observed.png", dpi=160)
            plt.close(fig)
    return sorted(p.name for p in fd.glob("*.png"))
