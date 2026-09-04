"""Predict pairwise strain interactions from genomic (KEGG KO) content.

Question: given the gene content of two strains, can we predict how they split a well when
grown together -- i.e. predict `mean_log2_ratio_a_over_b` from `relative_abundance`'s
r03_pair_replicate_stats.csv using only the two genomes?

Two things make this harder than the row count suggests, and both shape every design choice
below:

1.  **The strain-ID join is not free.** The genomic tables key strains by names from a
    384-well *source collection* (`mapping_384_well_plate_collection.csv`), while this
    project's pairwise experiments key strains by *their own* plate well-coordinates. The
    namespaces collide (same-looking labels, different organisms) -- see
    `validate_strain_join()`, which is not optional bookkeeping but the gate that decides
    whether an experiment may enter the modeling set at all. As of 2026-08-27 only 20260630
    passes; 20260721 fails decisively (see that function's docstring).

2.  **n is 74 genomes, not ~1000 pairs.** The pairs are not independent samples: they are all
    the pairwise combinations of 74 genomes. A model can score well on held-out *pairs*
    purely by memorizing per-strain competitiveness from the strain's other pairs. So every
    model is evaluated under two regimes (`cv_pair` and `cv_strain`) and the gap between them
    is the actual result -- see `cross_validate()`.

The target is antisymmetric (swapping a/b negates it). Rather than hoping a model learns
that, we build it in: features are split into an antisymmetric block (x_a - x_b) and a
symmetric block (x_a + x_b, overlap scores), training rows are augmented with both
orientations, and predictions are antisymmetrized as (f(a,b) - f(b,a))/2 -- exact, free, and
strictly better than relying on augmentation alone.

Call order: validate_strain_join -> build_dataset -> strain_feature_matrix ->
cross_validate -> make_all_figures.
"""

import os
import warnings
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

# ---- palette, matched to relative_abundance.py so figures read as one report ----
COLOR_BLUE = "#2a78d6"
COLOR_RED = "#e34948"
COLOR_GRID = "#d8d7d2"
COLOR_TEXT_SECONDARY = "#52514e"
COLOR_GOOD = "#0ca30c"
COLOR_CRITICAL = "#d03b3b"

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white", "axes.edgecolor": COLOR_GRID,
    "axes.grid": True, "grid.color": COLOR_GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
})

from paths import GENOMIC_TABLES
GENOMIC_DIR = GENOMIC_TABLES     # kept as an alias: several modules and runners import this name


@dataclass
class GenomicMLConfig:
    """Paths and knobs for the genomic-ML analysis of one experiment.

    Separate from ExperimentConfig because the genomic tables are shared across all
    experiments and are not part of any single experiment's sequencing run.
    """
    exp_cfg: object                                   # ExperimentConfig of the experiment supplying labels
    genomic_dir: Path = GENOMIC_DIR
    feature_table: str = "KEGG_ko_and_strains_table.csv"
    feature_columns: list = None   # optional subset of feature columns to use
    mapping_csv: str = "mapping_384_well_plate_collection.csv"

    # label filtering
    drop_high_uncertainty: bool = True   # near-identical 16S refs -> label is a forced ~50/50 artifact
    drop_unstable: bool = True           # replicates disagree beyond UNSTABLE_STD_THRESHOLD
    target: str = "mean_log2_ratio_a_over_b"

    # features
    min_prevalence: int = 4              # KO must be present in >=this many and <=n-this many strains
    n_pca: int = 30

    out_dir: Path = None

    def __post_init__(self):
        if self.out_dir is None:
            self.out_dir = self.exp_cfg.exp_base / "analysis" / "genomic_ml" / "outputs"
        self.fig_dir = self.out_dir / "figures"
        for d in (self.out_dir, self.fig_dir):
            d.mkdir(parents=True, exist_ok=True)

    @property
    def feature_path(self):
        return self.genomic_dir / self.feature_table

    @property
    def mapping_path(self):
        return self.genomic_dir / self.mapping_csv


# ===========================================================================
# 00 -- the strain-ID join, and proving it is real
# ===========================================================================

def load_strain_mapping(gcfg):
    """well-coordinate (source plate) -> genomic-table strain id.

    `Well_souce_plate` (sic) is unique across all 298 rows even though the `strain` column
    switches naming convention per `sequencing_batch` (bare integers for `Or`, Schulenberg
    names otherwise). Verified consistent with
    merge_consensus_sequences/plate_format/Full_384_strain_collection_no_seq_info.csv: all
    298 wells match, so there is exactly one source-collection namespace here.
    """
    m = pd.read_csv(gcfg.mapping_path)
    m["strain"] = m["strain"].astype(str)
    assert not m["Well_souce_plate"].duplicated().any(), "source-plate wells are not unique"
    return dict(zip(m["Well_souce_plate"], m["strain"])), m


def load_ko_table(gcfg):
    """Per-strain feature count matrix, with all-zero (annotation-failed) genomes dropped.

    Named for KEGG KO, but used for every annotation scheme sharing the 298-strain axis (BiGG,
    CAZy, KEGG Modules, PFAMs, panX). An optional `feature_columns` on the config restricts to a
    named subset -- used to compare a pre-selected feature list against a random list of the
    same size drawn from the same pool.
    """
    ko = pd.read_csv(gcfg.feature_path, index_col=0)
    ko.index = ko.index.astype(str)
    cols = getattr(gcfg, "feature_columns", None)
    if cols is not None:
        keep = [c for c in cols if c in ko.columns]
        if len(keep) < len(cols):
            warnings.warn(f"{len(cols) - len(keep)} requested feature columns absent from "
                          f"{gcfg.feature_path.name}")
        ko = ko[keep]
    empty = ko.index[ko.sum(axis=1) == 0]
    return ko.loc[ko.sum(axis=1) > 0], list(empty)


def validate_strain_join(gcfg, n_perm=200, seed=0):
    """Does joining this experiment's well labels to the genomic tables give the *right* genomes?

    A name match is not an identity match anywhere in this project (well-coordinate labels are
    reused across unrelated plates), so the join has to be earned. The test exploits a fact the
    labels cannot fake: two strains with near-identical 16S are near-identical organisms and
    must therefore have near-identical gene content. So across tested pairs, 16S divergence
    (`ref_pair_bp_dist`, measured from *this* experiment's own reads) should track KO-profile
    divergence (Jaccard, measured from the *joined* genomes). Under a wrong join the two are
    unrelated by construction.

    Null: permute the well -> genome assignment (preserves both marginals, destroys only the
    correspondence), recompute the correlation.

    Result as of 2026-08-27:
      20260630  rho=+0.363  z=+6.0 vs null   16S<=5bp -> KO Jaccard 0.024, >100bp -> 0.569  PASS
      20260721  rho=-0.046  z=-0.7 vs null   16S<=5bp -> KO Jaccard 0.613, >100bp -> 0.589  FAIL
    20260721's labels come from a different physical collection: flat across every 16S bin,
    indistinguishable from a random assignment of genomes to wells. Independently confirmed by
    comparing the two experiments' own 16S consensus sequences for the 18 well labels they
    share -- only 3/18 are the same organism (<=15bp), so the two experiments do not even
    share a namespace with each other. There is no mapping available that would rescue
    20260721, so it is excluded from modeling rather than joined on a name match.
    """
    w2s, mapping = load_strain_mapping(gcfg)
    ko, _ = load_ko_table(gcfg)
    pa = (ko > 0).astype(float)
    rng = np.random.default_rng(seed)

    r = pd.read_csv(gcfg.exp_cfg.relative_abundance_out_dir / "r03_pair_replicate_stats.csv")
    r = r.assign(genome_a=r["strain_a"].map(w2s), genome_b=r["strain_b"].map(w2s))
    d = r.dropna(subset=["genome_a", "genome_b", "ref_pair_bp_dist"])
    d = d[d["genome_a"].isin(pa.index) & d["genome_b"].isin(pa.index)].reset_index(drop=True)

    def _jaccard(ga, gb):
        A, B = pa.loc[ga].values, pa.loc[gb].values
        return 1 - (A * B).sum(1) / ((A + B) > 0).sum(1)

    jac = _jaccard(d["genome_a"], d["genome_b"])
    rho, p = spearmanr(d["ref_pair_bp_dist"], jac)

    null = []
    for _ in range(n_perm):
        perm = dict(zip(mapping["Well_souce_plate"], rng.permutation(mapping["strain"].values)))
        ga, gb = d["strain_a"].map(perm), d["strain_b"].map(perm)
        keep = ga.isin(pa.index) & gb.isin(pa.index)
        null.append(spearmanr(d.loc[keep, "ref_pair_bp_dist"], _jaccard(ga[keep], gb[keep]))[0])
    null = np.array(null)
    z = (rho - null.mean()) / null.std()

    bp = d["ref_pair_bp_dist"].values
    bins = {"16S<=5bp": bp <= 5, "16S 6-50bp": (bp > 5) & (bp <= 50),
            "16S 51-100bp": (bp > 100) & (bp <= 100), "16S>100bp": bp > 100}
    rows = [{"metric": f"median_ko_jaccard__{k}", "value": float(np.median(jac[v])) if v.sum() else np.nan,
             "n": int(v.sum())} for k, v in bins.items()]
    rows += [
        {"metric": "n_pairs_joinable", "value": float(len(d)), "n": len(r)},
        {"metric": "spearman_bpdist_vs_kojaccard", "value": float(rho), "n": len(d)},
        {"metric": "permutation_null_mean", "value": float(null.mean()), "n": n_perm},
        {"metric": "permutation_null_sd", "value": float(null.std()), "n": n_perm},
        {"metric": "z_vs_permutation_null", "value": float(z), "n": n_perm},
        {"metric": "join_verdict_pass", "value": float(z > 3 and rho > 0.15), "n": np.nan},
    ]
    out = pd.DataFrame(rows)
    out.to_csv(gcfg.out_dir / "g00_strain_join_validation.csv", index=False)
    return out, bool(z > 3 and rho > 0.15)


# ===========================================================================
# 01 -- the modeling dataset
# ===========================================================================

def build_dataset(gcfg):
    """One row per unique tested pair that survives QC, with both strains' genome ids attached.

    Dropped, and why:
      - pairs where either strain has no genome in the KO table (or an all-zero row)
      - `high_uncertainty_pair`: the two 16S references are too close to tell the strains
        apart, so the reported ~50/50 is a measurement ceiling, not a coexistence result.
        Training on these teaches "similar genomes -> 50/50", which is an artifact of the
        assay, not biology. This costs ~37% of the rows and is the single most important
        filter here.
      - `unstable_replicate`: replicates disagree more than the pipeline's stability threshold.
    """
    w2s, _ = load_strain_mapping(gcfg)
    ko, empty = load_ko_table(gcfg)

    r = pd.read_csv(gcfg.exp_cfg.relative_abundance_out_dir / "r03_pair_replicate_stats.csv")
    steps = [("r03_pairs_total", len(r))]

    r = r.assign(genome_a=r["strain_a"].map(w2s), genome_b=r["strain_b"].map(w2s))
    r = r.dropna(subset=["genome_a", "genome_b"])
    steps.append(("after_well_to_genome_join", len(r)))

    r = r[r["genome_a"].isin(ko.index) & r["genome_b"].isin(ko.index)]
    steps.append(("after_dropping_genomes_without_annotation", len(r)))

    if gcfg.drop_high_uncertainty:
        r = r[~r["high_uncertainty_pair"]]
        steps.append(("after_dropping_high_uncertainty_pairs", len(r)))
    if gcfg.drop_unstable:
        r = r[~r["unstable_replicate"]]
        steps.append(("after_dropping_unstable_replicates", len(r)))

    r = r.reset_index(drop=True)
    r["pair_id"] = [f"{min(a, b)}|{max(a, b)}" for a, b in zip(r["strain_a"], r["strain_b"])]
    # label reliability, for sample weighting / sensitivity analysis
    r["weight"] = r["n_replicates"] * np.log1p(r["mean_n_reads"]) * (1 - r["mean_uncertainty_score"])

    strains = sorted(set(r["strain_a"]) | set(r["strain_b"]))
    summary = pd.DataFrame(steps, columns=["step", "n_pairs"])
    summary.loc[len(summary)] = ["n_strains_final", len(strains)]
    summary.loc[len(summary)] = ["n_possible_pairs", len(strains) * (len(strains) - 1) // 2]
    summary.loc[len(summary)] = ["pct_of_possible_pairs_tested",
                                 round(100 * len(r) / (len(strains) * (len(strains) - 1) / 2), 1)]
    summary.loc[len(summary)] = ["n_genomes_dropped_all_zero_ko", len(empty)]
    summary.to_csv(gcfg.out_dir / "g01_dataset_summary.csv", index=False)
    r.to_csv(gcfg.out_dir / "g01_modeling_pairs.csv", index=False)
    return r, summary


def strain_feature_matrix(gcfg, pairs):
    """Per-strain genomic feature matrix, indexed by *well label* (the experiment's namespace).

    Presence/absence rather than raw counts: the counts are gene-copy numbers whose dynamic
    range (1-90) is dominated by a handful of repeat families, and with 74 genomes the extra
    variance is not affordable. Copy number survives only as the two summary columns.
    """
    w2s, _ = load_strain_mapping(gcfg)
    ko, _ = load_ko_table(gcfg)
    strains = sorted(set(pairs["strain_a"]) | set(pairs["strain_b"]))
    genomes = [w2s[s] for s in strains]

    sub = ko.loc[genomes]
    pa = (sub > 0).astype(float)
    n = len(strains)
    prev = pa.sum(axis=0)
    keep = (prev >= gcfg.min_prevalence) & (prev <= n - gcfg.min_prevalence)
    X = pa.loc[:, keep]
    X.index = strains
    X.columns = [c.replace("ko:", "") for c in X.columns]

    summary = pd.DataFrame({
        "n_ko_present": pa.sum(axis=1).values,
        "total_gene_copies": sub.sum(axis=1).values,
    }, index=strains)
    return X, summary


# ===========================================================================
# 02 -- pair features: an antisymmetric block and a symmetric block
# ===========================================================================

def _pair_design(Z, summ, pairs, extra_sym=None):
    """Stack per-pair features for BOTH orientations of every pair.

    Columns are [Z_a - Z_b | Z_a + Z_b | symmetric scalars | antisymmetric scalars]. The
    first block carries "who has which genes"; the second carries "what kind of pair is this
    at all" (a strong Serratia/Pseudomonas pair behaves differently from a strong/weak one).
    """
    a, b = pairs["strain_a"].values, pairs["strain_b"].values
    Za, Zb = Z.loc[a].values, Z.loc[b].values
    sa, sb = summ.loc[a].values, summ.loc[b].values

    def _block(Za, Zb, sa, sb, sym_scalars):
        anti = Za - Zb
        sym = Za + Zb
        anti_s = sa - sb
        return np.hstack([anti, sym, sym_scalars, anti_s])

    sym_scalars = np.column_stack([np.asarray(v) for v in extra_sym]) if extra_sym else np.zeros((len(a), 0))
    fwd = _block(Za, Zb, sa, sb, sym_scalars)
    rev = _block(Zb, Za, sb, sa, sym_scalars)
    names = ([f"anti_{c}" for c in Z.columns] + [f"sym_{c}" for c in Z.columns]
             + [f"sym_x{i}" for i in range(sym_scalars.shape[1])]
             + [f"anti_{c}" for c in summ.columns])
    return fwd, rev, names


def pair_feature_blocks(pairs, Z, summ):
    """Symmetric scalar covariates available for every pair, orientation-independent."""
    return [pairs["ref_pair_bp_dist"].values]


# ===========================================================================
# 03 -- evaluation: two regimes, because only one of them is a real test
# ===========================================================================

def make_folds(pairs, kind, n_splits=5, seed=0):
    """Yield (train_pair_idx, test_pair_idx) under one of two very different questions.

    `cv_pair` -- hold out *pairs*, keep all strains visible. "Fill in the ~64% of pairs that
        were never tested, in a collection we have already characterized." Both orientations
        of a pair always land in the same fold (otherwise the model sees the test label with
        the inputs swapped -- total leakage). Note this regime still lets a model learn a
        strain's competitiveness from that strain's *other* pairs, so a good score here does
        NOT show that genome content predicts interactions; it mostly shows the outcome is
        strain-driven and the model found the strain.

    `cv_strain` -- hold out *strains*, and test only on pairs where **both** strains were
        unseen. "Sequence a new isolate; predict how it will do against another new isolate."
        Pairs straddling the train/test strain split are dropped entirely rather than used,
        since half their information is already in the training set. This is the honest test
        of genome -> interaction generalization, and the only one where a per-strain lookup
        baseline is unavailable by construction.
    """
    rng = np.random.default_rng(seed)
    if kind == "cv_pair":
        groups = pairs["pair_id"].values
        uniq = np.array(sorted(set(groups)))
        fold_of = dict(zip(uniq, rng.permutation(len(uniq)) % n_splits))
        f = np.array([fold_of[g] for g in groups])
        for k in range(n_splits):
            yield np.where(f != k)[0], np.where(f == k)[0]
    elif kind == "cv_strain":
        strains = np.array(sorted(set(pairs["strain_a"]) | set(pairs["strain_b"])))
        fold_of = dict(zip(strains, rng.permutation(len(strains)) % n_splits))
        fa = pairs["strain_a"].map(fold_of).values
        fb = pairs["strain_b"].map(fold_of).values
        for k in range(n_splits):
            test = np.where((fa == k) & (fb == k))[0]
            train = np.where((fa != k) & (fb != k))[0]
            if len(test) >= 5:
                yield train, test
    else:
        raise ValueError(kind)


def _metrics(y, yhat, sign_threshold=1.0):
    y, yhat = np.asarray(y, float), np.asarray(yhat, float)
    ok = np.isfinite(y) & np.isfinite(yhat)
    y, yhat = y[ok], yhat[ok]
    decisive = np.abs(y) > sign_threshold
    # r2 and mae stay defined for a constant predictor -- the all-zero baseline scores exactly
    # r2 = 0 and is the reference every other row is read against, so it must not come back NaN.
    out = {
        "n": len(y),
        "r2": float(1 - np.sum((y - yhat) ** 2) / np.sum(y ** 2)) if len(y) else np.nan,
        "mae": float(np.mean(np.abs(y - yhat))) if len(y) else np.nan,
        "n_sign_evaluated": int(decisive.sum()),
    }
    if len(y) < 3 or np.std(yhat) == 0:   # correlations/sign are undefined without spread
        return {**out, "pearson_r": np.nan, "spearman_rho": np.nan, "sign_accuracy": np.nan}
    return {
        **out,
        "pearson_r": float(pearsonr(y, yhat)[0]),
        "spearman_rho": float(spearmanr(y, yhat)[0]),
        "sign_accuracy": float(np.mean(np.sign(y[decisive]) == np.sign(yhat[decisive]))) if decisive.sum() else np.nan,
    }


def _fit_additive_strength(pairs, idx, strains, ridge_lambda=1.0, target="mean_log2_ratio_a_over_b"):
    """Least-squares per-strain strength s with log2ratio(a,b) ~ s_a - s_b, fit IN-FOLD.

    Same shape as the Bradley-Terry fit in relative_abundance.py, but refit from training
    pairs only -- the published r05_bt_strengths.csv was fit on every pair including the ones
    we are about to predict, so using it anywhere in a CV loop would leak the answer.
    """
    pos = {s: i for i, s in enumerate(strains)}
    sub = pairs.iloc[idx]
    D = np.zeros((len(sub), len(strains)))
    D[np.arange(len(sub)), sub["strain_a"].map(pos).values] = 1
    D[np.arange(len(sub)), sub["strain_b"].map(pos).values] = -1
    y = sub[target].values
    # + a mean-zero constraint row, since s is only identified up to an additive constant
    A = np.vstack([D, np.ones((1, len(strains)))])
    b = np.concatenate([y, [0.0]])
    s = np.linalg.lstsq(A.T @ A + ridge_lambda * np.eye(len(strains)), A.T @ b, rcond=None)[0]
    return pd.Series(s, index=strains)


# ---- the model zoo -------------------------------------------------------
# Every model exposes fit_predict(ctx) -> antisymmetrized predictions for ctx.test pairs.

def _scaled(Z, train_strains):
    """Standardize a per-strain representation on the training strains.

    Ridge penalizes coefficients on the raw feature scale, so a representation whose columns
    happen to be large (16S PCs are in units of base pairs, ~10^2; KO PCs are ~10^0) is
    effectively penalized far harder. Without this, the taxonomy control loses on units rather
    than on information -- which would make the KO-vs-16S comparison meaningless.
    """
    sc = StandardScaler().fit(Z.loc[train_strains].values)
    return pd.DataFrame(sc.transform(Z.values), index=Z.index, columns=Z.columns)


class _Ctx:
    """Everything a model needs for one fold, with PCA already fit on training strains only."""
    def __init__(self, gcfg, pairs, X, summ, train_idx, test_idx, phylo=None):
        self.gcfg, self.pairs, self.X, self.summ = gcfg, pairs, X, summ
        self.train_idx, self.test_idx = train_idx, test_idx
        self.target = gcfg.target
        self.train_strains = sorted(set(pairs.iloc[train_idx]["strain_a"]) | set(pairs.iloc[train_idx]["strain_b"]))
        self.y_train = pairs.iloc[train_idx][self.target].values
        self.y_test = pairs.iloc[test_idx][self.target].values

        # PCA fit on training strains ONLY -- unsupervised, but keeping it in-fold costs
        # nothing and removes any argument about the held-out genomes shaping the basis
        n_comp = min(gcfg.n_pca, len(self.train_strains) - 1, X.shape[1])
        pca = PCA(n_components=n_comp, random_state=0).fit(X.loc[self.train_strains].values)
        Z = _scaled(pd.DataFrame(pca.transform(X.values), index=X.index,
                                 columns=[f"pc{i+1}" for i in range(n_comp)]), self.train_strains)
        sc = StandardScaler().fit(summ.loc[self.train_strains].values)
        S = pd.DataFrame(sc.transform(summ.values), index=summ.index, columns=summ.columns)
        self.Z, self.S, self.pca = Z, S, pca

        # 16S landmark representation: each strain described by its distances to the
        # TRAINING strains, then PCA'd -- an in-fold stand-in for "where does this sit on the
        # tree", so a genomic model can be asked whether it beats plain taxonomy
        self.Zp = None
        if phylo is not None:
            # the 16S matrix is built for one strain set; a different feature table keeps a
            # different set of genomes, so intersect rather than index blindly (genomic_ml_yield
            # already did this -- which is why the yield half of the sweep survived and this
            # half raised KeyError)
            tr_p = [t for t in self.train_strains if t in phylo.columns and t in phylo.index]
            rows_p = [r for r in X.index if r in phylo.index]
            if len(tr_p) < 3 or not rows_p:
                self.Zp = None
                self.fwd, self.rev, self.feature_names = self.design(Z)
                self._raw = None
                return
            L = phylo.loc[rows_p, tr_p]
            kp = min(gcfg.n_pca, len(tr_p) - 1)
            pca_p = PCA(n_components=kp, random_state=0).fit(L.loc[tr_p].values)
            self.Zp = _scaled(pd.DataFrame(pca_p.transform(L.values), index=L.index,
                                           columns=[f"ph{i+1}" for i in range(kp)]), tr_p)

        self.fwd, self.rev, self.feature_names = self.design(Z)
        self._raw = None

    def design(self, Z):
        sym = pair_feature_blocks(self.pairs, Z, self.S)
        return _pair_design(Z, self.S, self.pairs, extra_sym=sym)

    @property
    def raw(self):
        """Un-reduced KO design (~9.5k columns), built on first use only.

        Materialising it costs a few hundred MB per fold, and only xgboost_raw_ko wants it --
        so folds that never touch it should not pay for it.
        """
        if self._raw is None:
            self._raw = self.design(self.X)
        return self._raw

    def xy(self, fwd, rev):
        """Training rows in both orientations + the matching test matrices, for any design."""
        i = self.train_idx
        Xtr = np.vstack([fwd[i], rev[i]])
        ytr = np.concatenate([self.y_train, -self.y_train])
        return Xtr, ytr, fwd[self.test_idx], rev[self.test_idx]

    def train_xy(self, raw=False):
        """Training rows in BOTH orientations, labels negated for the reversed copy."""
        F, R = self.raw[:2] if raw else (self.fwd, self.rev)
        i = self.train_idx
        return np.vstack([F[i], R[i]]), np.concatenate([self.y_train, -self.y_train])

    def test_x(self, raw=False):
        F, R = self.raw[:2] if raw else (self.fwd, self.rev)
        return F[self.test_idx], R[self.test_idx]


def _antisym(model, Xf, Xr):
    """(f(a,b) - f(b,a))/2 -- makes the prediction exactly antisymmetric regardless of model."""
    return (model.predict(Xf) - model.predict(Xr)) / 2.0


def m_zero(ctx):
    return np.zeros(len(ctx.test_idx))


def m_strength_observed(ctx):
    """Upper-bound reference: per-strain strength fit from the training pairs, no genomics.

    Under cv_pair this is the model to beat -- it is what you can already do with the
    interaction data alone. Under cv_strain the test strains have no training pairs at all,
    so it necessarily falls back to 0; that collapse is the point.
    """
    s = _fit_additive_strength(ctx.pairs, ctx.train_idx, ctx.train_strains, target=ctx.target)
    te = ctx.pairs.iloc[ctx.test_idx]
    return (te["strain_a"].map(s).fillna(0.0) - te["strain_b"].map(s).fillna(0.0)).values


def m_ridge(ctx, alpha=100.0):
    Xtr, ytr = ctx.train_xy()
    Xf, Xr = ctx.test_x()
    mdl = Ridge(alpha=alpha).fit(Xtr, ytr)
    return _antisym(mdl, Xf, Xr)


def m_xgb(ctx, **kw):
    from xgboost import XGBRegressor
    Xtr, ytr = ctx.train_xy()
    Xf, Xr = ctx.test_x()
    params = dict(n_estimators=400, max_depth=3, learning_rate=0.05, subsample=0.8,
                  colsample_bytree=0.5, reg_lambda=5.0, min_child_weight=5,
                  n_jobs=4, random_state=0, verbosity=0)
    params.update(kw)
    mdl = XGBRegressor(**params).fit(Xtr, ytr)
    return _antisym(mdl, Xf, Xr)


def m_xgb_raw_ko(ctx, **kw):
    """XGBoost straight on the ~4400 KO presence/absence differences, no PCA."""
    from xgboost import XGBRegressor
    Xtr, ytr = ctx.train_xy(raw=True)
    Xf, Xr = ctx.test_x(raw=True)
    params = dict(n_estimators=400, max_depth=4, learning_rate=0.05, subsample=0.8,
                  colsample_bytree=0.2, reg_lambda=5.0, min_child_weight=5,
                  n_jobs=4, random_state=0, verbosity=0)
    params.update(kw)
    mdl = XGBRegressor(**params).fit(Xtr, ytr)
    return _antisym(mdl, Xf, Xr)


def _tabpfn_login(gcfg=None):
    """Authenticate the hosted TabPFN client from an env var or a local, unversioned token file.

    The token is a credential, so it is never written into this module. Order: $TABPFN_TOKEN,
    then <genomic_ml>/.tabpfn_token (mode 600).
    """
    import tabpfn_client
    tok = os.environ.get("TABPFN_TOKEN")
    if not tok and gcfg is not None:
        f = gcfg.out_dir.parent / ".tabpfn_token"
        if f.exists():
            tok = f.read_text().strip()
    if not tok:
        raise RuntimeError("no TabPFN token: set $TABPFN_TOKEN or write <genomic_ml>/.tabpfn_token")
    for fn in ("set_access_token", "init"):
        if hasattr(tabpfn_client, fn):
            try:
                getattr(tabpfn_client, fn)(tok) if fn == "set_access_token" else getattr(tabpfn_client, fn)(access_token=tok)
                return
            except TypeError:
                continue
    raise RuntimeError("tabpfn_client exposes no usable login entry point")


def m_tabpfn(ctx):
    """TabPFN via the hosted API.

    NOTE: this sends the fold's feature and label rows to TabPFN's servers -- the only step
    here that moves this project's data off the machine. It is a good fit on paper (TabPFN is
    built for exactly this regime: <10k rows, <500 features, small-n tabular regression), which
    is why it is worth the trip.
    """
    from tabpfn_client import TabPFNRegressor
    _tabpfn_login(ctx.gcfg)
    Xtr, ytr = ctx.train_xy()
    Xf, Xr = ctx.test_x()
    mdl = TabPFNRegressor().fit(Xtr, ytr)
    return _antisym(mdl, Xf, Xr)


def _two_stage(ctx, learner, Z=None):
    """Genome -> per-strain strength -> pairwise outcome.

    Far better posed than predicting pairs directly: the regression has one row per *strain*
    (74, matching the true independent-sample count) and one target, and the resulting
    prediction s_a - s_b is exactly antisymmetric by construction. It also encodes the thing
    relative_abundance already established -- the network is ~83% a dominance hierarchy -- as
    a structural prior rather than something the model must rediscover from pairs.
    """
    Z = ctx.Z if Z is None else Z
    s = _fit_additive_strength(ctx.pairs, ctx.train_idx, ctx.train_strains, target=ctx.target)
    Ztr = Z.loc[ctx.train_strains].values
    Str = ctx.S.loc[ctx.train_strains].values
    Atr = np.hstack([Ztr, Str])
    mdl = learner().fit(Atr, s.loc[ctx.train_strains].values)
    Aall = np.hstack([Z.values, ctx.S.values])
    shat = pd.Series(mdl.predict(Aall), index=Z.index)
    te = ctx.pairs.iloc[ctx.test_idx]
    return (te["strain_a"].map(shat) - te["strain_b"].map(shat)).values


def m_two_stage_ridge(ctx):
    return _two_stage(ctx, lambda: Ridge(alpha=50.0))


def m_two_stage_xgb(ctx):
    from xgboost import XGBRegressor
    return _two_stage(ctx, lambda: XGBRegressor(
        n_estimators=300, max_depth=2, learning_rate=0.05, subsample=0.8,
        colsample_bytree=0.5, reg_lambda=5.0, min_child_weight=3,
        n_jobs=4, random_state=0, verbosity=0))


def m_two_stage_tabpfn(ctx):
    """TabPFN on the 74-strain genome -> competitiveness regression (see _two_stage).

    The best-matched use of TabPFN here: a genuinely small, wide-ish tabular problem, which is
    the regime it was trained for. Same off-machine caveat as m_tabpfn.
    """
    from tabpfn_client import TabPFNRegressor
    _tabpfn_login(ctx.gcfg)
    return _two_stage(ctx, TabPFNRegressor)


def m_ridge_phylo16s(ctx, alpha=100.0):
    """Control: same machinery, but the strain is described only by its 16S position.

    KEGG content and 16S are both largely phylogenetic, so a genomic model beating the
    zero baseline proves little on its own -- it may only have rediscovered "Pseudomonas beats
    Sphingomonas". The comparison that matters is genomic vs. this row: whatever KO content
    buys *over* 16S is the part that is about gene content rather than taxonomy.
    """
    if ctx.Zp is None:
        raise RuntimeError("no 16S distance matrix supplied")
    Xtr, ytr, Xf, Xr = ctx.xy(*ctx.design(ctx.Zp)[:2])
    return _antisym(Ridge(alpha=alpha).fit(Xtr, ytr), Xf, Xr)


def m_two_stage_ridge_phylo16s(ctx):
    if ctx.Zp is None:
        raise RuntimeError("no 16S distance matrix supplied")
    return _two_stage(ctx, lambda: Ridge(alpha=50.0), Z=ctx.Zp)


def m_two_stage_ridge_ko_plus_phylo(ctx):
    if ctx.Zp is None:
        return m_two_stage_ridge(ctx)
    Z = pd.concat([ctx.Z, ctx.Zp], axis=1)
    return _two_stage(ctx, lambda: Ridge(alpha=50.0), Z=Z)


MODELS = {
    "zero_baseline": m_zero,
    "strength_observed_no_genomics": m_strength_observed,
    "ridge_pca": m_ridge,
    "xgboost_pca": m_xgb,
    "xgboost_raw_ko": m_xgb_raw_ko,
    "ridge_phylo16s_only": m_ridge_phylo16s,
    "two_stage_ridge": m_two_stage_ridge,
    "two_stage_ridge_phylo16s_only": m_two_stage_ridge_phylo16s,
    "two_stage_ridge_ko_plus_phylo": m_two_stage_ridge_ko_plus_phylo,
    "two_stage_xgboost": m_two_stage_xgb,
    "tabpfn": m_tabpfn,
    "two_stage_tabpfn": m_two_stage_tabpfn,
}


def summarize_folds(fold_df):
    """Collapse per-fold metrics into the model x regime table, ranked by R²."""
    num = ["pearson_r", "spearman_rho", "r2", "mae", "sign_accuracy"]
    return (fold_df.groupby(["regime", "model"])
            .agg(n_folds=("fold", "nunique"), n_test=("n", "sum"),
                 **{f"{c}_mean": (c, "mean") for c in num},
                 **{f"{c}_sd": (c, "std") for c in num})
            .reset_index().sort_values(["regime", "r2_mean"], ascending=[True, False]))


def label_noise_ceiling(gcfg, pairs):
    """How much of the target is even predictable, given how noisy the labels are?

    An R^2 of 0.3 means something completely different if the labels themselves only reproduce
    to R^2 = 0.4. Two independent estimates of the measurement variance:

      replicate -- from the 2-replicate pairs' within-pair spread of log2 ratio. Optimistic:
                   only ~15% of pairs have 2 replicates, and `unstable_replicate` pairs (the
                   ones whose replicates disagreed most) were already filtered out upstream,
                   so this subset is selected for agreement. Treat as a loose upper bound.
      binomial  -- from read depth alone: a well read ~40x cannot resolve a 55/45 split, and
                   the delta-method variance of log2(p/(1-p)) at the observed p and n blows up
                   near the extremes, which is exactly where this target lives.

    Ceiling R^2 = 1 - sigma^2_noise / var(y): the best any model could score against these
    labels.
    """
    y = pairs[gcfg.target].values
    var_y = float(np.var(y))
    rows = []

    rep = pairs[pairs["n_replicates"] >= 2]
    if len(rep):
        # std over n replicates -> variance of the replicate MEAN that forms the label
        sigma2_rep = float(np.mean(rep["std_log2_ratio_a_over_b"].values ** 2 / rep["n_replicates"].values))
        rows.append({"estimator": "replicate_spread", "n": len(rep), "sigma2_noise": sigma2_rep,
                     "var_label": var_y, "ceiling_r2": 1 - sigma2_rep / var_y})

    p = np.clip(pairs["mean_relative_abundance_a"].values, 1e-6, 1 - 1e-6)
    n_reads = np.maximum(pairs["mean_n_reads"].values * pairs["n_replicates"].values, 1)
    pe = np.clip(p, 0.5 / n_reads, 1 - 0.5 / n_reads)   # unresolvable beyond +-1 read
    # nanmean: a handful of pairs carry a log2 ratio but no mean_relative_abundance_a
    sigma2_bin = float(np.nanmean((1 / np.log(2) ** 2) * (1 / (n_reads * pe * (1 - pe)))))
    rows.append({"estimator": "binomial_read_depth", "n": len(pairs), "sigma2_noise": sigma2_bin,
                 "var_label": var_y, "ceiling_r2": 1 - sigma2_bin / var_y})

    out = pd.DataFrame(rows)
    out.to_csv(gcfg.out_dir / "g02_label_noise_ceiling.csv", index=False)
    return out


def cross_validate(gcfg, pairs, X, summ, models=None, regimes=("cv_pair", "cv_strain"),
                   n_splits=5, seed=0, shuffle_control=True, phylo=None, n_repeats=1,
                   file_prefix="g03"):
    """Run every model under every regime; also run a label-shuffled control per regime.

    The shuffled control (permute y across training pairs, keep everything else) is the
    empirical zero: with 74 genomes and thousands of features, a model can post a positive
    R^2 on noise, and the control shows how much.
    """
    models = models or list(MODELS)
    per_fold, per_pair = [], []

    for regime in regimes:
        # cv_strain leaves only ~40 pairs per fold, so a single 5-fold pass is a noisy readout
        # of a noisy quantity; repeats with different strain partitions stabilise the mean
        folds = [(rep, tr, te) for rep in range(n_repeats)
                 for tr, te in make_folds(pairs, regime, n_splits=n_splits, seed=seed + 97 * rep)]
        for fi, (rep, tr, te) in enumerate(folds):
            ctx = _Ctx(gcfg, pairs, X, summ, tr, te, phylo=phylo)
            for name in models:
                try:
                    yhat = MODELS[name](ctx)
                except Exception as e:            # a missing optional dep must not kill the sweep
                    warnings.warn(f"{name} failed on {regime} fold {fi}: {type(e).__name__}: {e}")
                    continue
                per_fold.append({"regime": regime, "fold": fi, "model": name,
                                 "n_train": len(tr), **_metrics(ctx.y_test, yhat)})
                per_pair.append(pd.DataFrame({
                    "regime": regime, "fold": fi, "model": name,
                    "pair_id": pairs.iloc[te]["pair_id"].values,
                    "strain_a": pairs.iloc[te]["strain_a"].values,
                    "strain_b": pairs.iloc[te]["strain_b"].values,
                    "y_true": ctx.y_test, "y_pred": yhat}))

            if shuffle_control:
                rng = np.random.default_rng(1000 + fi)
                sh = pairs.copy()
                perm = rng.permutation(tr)
                sh.iloc[tr, sh.columns.get_loc(gcfg.target)] = pairs.iloc[perm][gcfg.target].values
                ctx_s = _Ctx(gcfg, sh, X, summ, tr, te, phylo=phylo)
                for name in ("ridge_pca", "xgboost_pca", "two_stage_ridge"):
                    if name not in models:
                        continue
                    try:
                        yhat = MODELS[name](ctx_s)
                    except Exception:
                        continue
                    per_fold.append({"regime": regime, "fold": fi, "model": f"SHUFFLED_{name}",
                                     "n_train": len(tr), **_metrics(ctx.y_test, yhat)})

    fold_df = pd.DataFrame(per_fold)
    pair_df = pd.concat(per_pair, ignore_index=True) if per_pair else pd.DataFrame()

    summ_df = summarize_folds(fold_df)

    fold_df.to_csv(gcfg.out_dir / f"{file_prefix}_cv_per_fold.csv", index=False)
    summ_df.to_csv(gcfg.out_dir / f"{file_prefix}_cv_summary.csv", index=False)
    pair_df.to_csv(gcfg.out_dir / f"{file_prefix}_cv_predictions.csv", index=False)
    return summ_df, fold_df, pair_df


# ===========================================================================
# 04 -- which genes, and can genomics reproduce the hierarchy at all
# ===========================================================================

def phylo_distance_matrix(gcfg, pairs, cache=True, strict=True, out_name=None):
    """16S edit-distance matrix over the modeled strains, from this experiment's own consensus.

    Deliberately the same sequences the interaction pipeline used as its mapping reference, so
    "taxonomy" here means exactly what it means everywhere else in this project rather than an
    outside tree with its own naming problems.
    """
    import edlib
    out = gcfg.out_dir / (out_name or "g04_phylo16s_distance_matrix.csv")
    strains = sorted(set(pairs["strain_a"]) | set(pairs["strain_b"]))
    if cache and out.exists():
        d = pd.read_csv(out, index_col=0)
        if list(d.index) == strains:
            return d

    seqs, name = {}, None
    for line in open(gcfg.exp_cfg.ra_reference_fasta):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line
    missing = [x for x in strains if x not in seqs]
    if missing and strict:
        raise KeyError(f"no 16S consensus for {missing[:5]}")
    if missing:
        # strict=False: the OD-based analysis covers strains the sequencing never resolved, so
        # dropping them here keeps the taxonomy control available on the subset that has a
        # consensus rather than discarding the control entirely
        strains = [x for x in strains if x in seqs]

    D = np.zeros((len(strains), len(strains)))
    for i, a in enumerate(strains):
        for j in range(i + 1, len(strains)):
            d = edlib.align(seqs[a], seqs[strains[j]], mode="NW", task="distance")["editDistance"]
            D[i, j] = D[j, i] = d
    df = pd.DataFrame(D, index=strains, columns=strains)
    df.to_csv(out)
    return df

def genome_to_strength(gcfg, pairs, X, summ, phylo=None, n_splits=10, seed=0, n_repeats=3):
    """Leave-strains-out: predict a strain's competitiveness from its genome alone.

    This is the two-stage model's first half examined on its own terms, and the cleanest
    statement of the underlying question -- 74 rows, one target, no pair structure to inflate
    anything. Run for each available strain representation (KO content, 16S position, both),
    because the interesting number is not "can we predict it" but "does gene content say more
    than the taxon label".

    Also reports the KOs whose presence most separates strong from weak strains, refit on all
    strains: descriptive only, not cross-validated, and -- given that KO content tracks
    phylogeny -- these mostly mark clades rather than mechanisms.
    """
    strains = sorted(set(pairs["strain_a"]) | set(pairs["strain_b"]))
    s_full = _fit_additive_strength(pairs, np.arange(len(pairs)), strains, target=gcfg.target)

    sources = {"ko_content": "ko"}
    if phylo is not None:
        sources.update({"phylo16s_only": "phylo", "ko_plus_phylo16s": "both"})

    rows, metrics = [], []
    for src_name, src in sources.items():
        acc = {}
        for rep in range(n_repeats):
            rng = np.random.default_rng(seed + 97 * rep)
            fold_of = dict(zip(strains, rng.permutation(len(strains)) % n_splits))
            for k in range(n_splits):
                tr_s = [x for x in strains if fold_of[x] != k]
                te_s = [x for x in strains if fold_of[x] == k]
                if not te_s:
                    continue
                keep = pairs["strain_a"].isin(tr_s) & pairs["strain_b"].isin(tr_s)
                s_tr = _fit_additive_strength(pairs, np.where(keep)[0], tr_s, target=gcfg.target)
                sc = StandardScaler().fit(summ.loc[tr_s].values)

                blocks = []
                if src in ("ko", "both"):
                    pk = PCA(n_components=min(gcfg.n_pca, len(tr_s) - 1, X.shape[1]),
                             random_state=0).fit(X.loc[tr_s].values)
                    sk = StandardScaler().fit(pk.transform(X.loc[tr_s].values))
                    blocks.append(lambda ss, pk=pk, sk=sk: sk.transform(pk.transform(X.loc[ss].values)))
                if src in ("phylo", "both"):
                    L = phylo.loc[:, tr_s]
                    pp = PCA(n_components=min(gcfg.n_pca, len(tr_s) - 1),
                             random_state=0).fit(L.loc[tr_s].values)
                    sp = StandardScaler().fit(pp.transform(L.loc[tr_s].values))
                    blocks.append(lambda ss, pp=pp, sp=sp, L=L: sp.transform(pp.transform(L.loc[ss].values)))
                A = lambda ss: np.hstack([b(ss) for b in blocks] + [sc.transform(summ.loc[ss].values)])

                mdl = Ridge(alpha=50.0).fit(A(tr_s), s_tr.loc[tr_s].values)
                for s_, v in zip(te_s, mdl.predict(A(te_s))):
                    acc.setdefault(s_, []).append(v)

        pred = {k: float(np.mean(v)) for k, v in acc.items()}
        df = pd.DataFrame({"strain": list(pred), "feature_source": src_name,
                           "strength_observed": [s_full[x] for x in pred],
                           "strength_predicted_from_genome": [pred[x] for x in pred]})
        rows.append(df)
        metrics.append({"feature_source": src_name,
                        **_metrics(df["strength_observed"], df["strength_predicted_from_genome"],
                                   sign_threshold=0.0)})

    out = pd.concat(rows, ignore_index=True)
    out["n_pairs"] = out["strain"].map(pd.concat([pairs["strain_a"], pairs["strain_b"]]).value_counts())
    met = pd.DataFrame(metrics)

    sv = s_full.loc[X.index].values
    corr = np.array([spearmanr(X.iloc[:, j].values, sv)[0] if X.iloc[:, j].std() > 0 else np.nan
                     for j in range(X.shape[1])])
    ko_df = (pd.DataFrame({"ko": X.columns, "spearman_vs_strength": corr,
                           "prevalence": X.sum(axis=0).values})
             .dropna().sort_values("spearman_vs_strength", key=np.abs, ascending=False))

    out.to_csv(gcfg.out_dir / "g04_strength_from_genome.csv", index=False)
    met.to_csv(gcfg.out_dir / "g04_strength_from_genome_metrics.csv", index=False)
    ko_df.to_csv(gcfg.out_dir / "g04_ko_strength_correlations.csv", index=False)
    return out, met, ko_df


# ===========================================================================
# 05 -- figures
# ===========================================================================

def make_all_figures(gcfg, summ_df, fold_df, pair_df, strength_df, ceiling_df):
    fd = gcfg.fig_dir

    # f01 -- model comparison across the two regimes
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8))
    XLO = -0.8   # floor: the 16S-only pair models miss by R^2 ~ -7 on unseen strains and would
                 # otherwise compress every informative bar into a sliver
    for ax, regime, title in zip(axes, ["cv_pair", "cv_strain"],
                                 ["cv_pair: unseen PAIR, strains already seen\n(a per-strain lookup is available -- inflated)",
                                  "cv_strain: unseen STRAINS, both held out\n(no lookup possible -- the honest test)"]):
        d = summ_df[summ_df["regime"] == regime].sort_values("r2_mean")
        if not len(d):
            continue
        cols = [COLOR_CRITICAL if m.startswith("SHUFFLED") else
                (COLOR_TEXT_SECONDARY if m in ("zero_baseline", "strength_observed_no_genomics")
                 else COLOR_BLUE) for m in d["model"]]
        vals = d["r2_mean"].values
        clipped = vals < XLO
        ax.barh(d["model"], np.where(clipped, XLO, vals), color=cols,
                xerr=np.where(clipped, 0, d["r2_sd"].values),
                error_kw=dict(ecolor=COLOR_TEXT_SECONDARY, lw=1, capsize=3))
        for i, (v, c) in enumerate(zip(vals, clipped)):
            if c:
                ax.text(XLO + 0.02, i, f"R² = {v:.1f}  ◄ off scale", va="center", ha="left",
                        fontsize=8, color="white", fontweight="bold")
        ax.axvline(0, color="black", lw=1)
        # the conservative (read-depth) ceiling only -- the replicate estimate is optimistic,
        # and two near-identical dashed lines just crowd the panel
        cr = float(ceiling_df.set_index("estimator").loc["binomial_read_depth", "ceiling_r2"])
        ax.axvline(cr, color=COLOR_GOOD, lw=1.2, ls="--")
        ax.text(cr - 0.015, len(d) - 0.5, f"label-noise ceiling {cr:.2f} ", color=COLOR_GOOD,
                fontsize=8, ha="right", va="center")
        ax.set_xlim(XLO, 1.02)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("R² vs. predicting no winner  (mean ± sd over folds)")
    fig.suptitle("Predicting log2 abundance ratio from KEGG KO content", fontweight="bold")
    fig.tight_layout()
    fig.savefig(fd / "f01_model_comparison.png", dpi=160)
    plt.close(fig)

    # f02 -- predicted vs observed for the best genomic model in each regime
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    for ax, regime in zip(axes, ["cv_pair", "cv_strain"]):
        d = summ_df[(summ_df["regime"] == regime) & (~summ_df["model"].str.startswith("SHUFFLED"))
                    & (~summ_df["model"].isin(["zero_baseline", "strength_observed_no_genomics"]))]
        if not len(d) or not len(pair_df):
            continue
        best = d.iloc[0]["model"]
        p = pair_df[(pair_df["regime"] == regime) & (pair_df["model"] == best)]
        ax.scatter(p["y_pred"], p["y_true"], s=14, alpha=0.45, color=COLOR_BLUE, edgecolor="none")
        lim = max(np.abs(p["y_true"]).max(), np.abs(p["y_pred"]).max()) * 1.05
        ax.plot([-lim, lim], [-lim, lim], color=COLOR_TEXT_SECONDARY, lw=1, ls="--")
        ax.axhline(0, color=COLOR_GRID, lw=1); ax.axvline(0, color=COLOR_GRID, lw=1)
        rho = spearmanr(p["y_pred"], p["y_true"])[0]
        ax.set_title(f"{regime} -- {best}\nSpearman ρ = {rho:.2f}, n = {len(p)}")
        ax.set_xlabel("predicted log2 ratio (a/b)"); ax.set_ylabel("observed log2 ratio (a/b)")
    fig.tight_layout()
    fig.savefig(fd / "f02_predicted_vs_observed.png", dpi=160)
    plt.close(fig)

    # f03 -- per-strain competitiveness predicted from genome alone (leave-strains-out)
    if len(strength_df):
        srcs = list(dict.fromkeys(strength_df["feature_source"]))
        fig, axes = plt.subplots(1, len(srcs), figsize=(5.2 * len(srcs), 5.4), squeeze=False)
        for ax, src in zip(axes[0], srcs):
            d = strength_df[strength_df["feature_source"] == src]
            ax.scatter(d["strength_predicted_from_genome"], d["strength_observed"],
                       s=np.clip(d["n_pairs"], 8, 60), color=COLOR_BLUE, alpha=0.7,
                       edgecolor="white", linewidth=0.6)
            rho = spearmanr(d["strength_predicted_from_genome"], d["strength_observed"])[0]
            ax.axhline(0, color=COLOR_GRID, lw=1); ax.axvline(0, color=COLOR_GRID, lw=1)
            ax.set_xlabel("predicted (held-out strain)")
            ax.set_ylabel("observed competitiveness (in-fold fit)")
            ax.set_title(f"{src}\nSpearman ρ = {rho:.2f}, n = {len(d)}")
        fig.suptitle("Can a genome alone say how competitive a strain is?  (point size = pairs tested)",
                     fontweight="bold")
        fig.tight_layout()
        fig.savefig(fd / "f03_strength_from_genome.png", dpi=160)
        plt.close(fig)

    # f04 -- sign accuracy: does the model at least call the winner?
    fig, ax = plt.subplots(figsize=(9, 5))
    d = summ_df[~summ_df["model"].str.startswith("SHUFFLED")]
    w = 0.38
    mods = list(dict.fromkeys(d["model"]))
    for i, (regime, color) in enumerate([("cv_pair", COLOR_BLUE), ("cv_strain", COLOR_RED)]):
        dd = d[d["regime"] == regime].set_index("model").reindex(mods)
        ax.bar(np.arange(len(mods)) + i * w, dd["sign_accuracy_mean"], w,
               yerr=dd["sign_accuracy_sd"], color=color, label=regime,
               error_kw=dict(ecolor=COLOR_TEXT_SECONDARY, lw=1, capsize=2))
    ax.axhline(0.5, color=COLOR_CRITICAL, lw=1.2, ls="--")
    ax.text(0.995, 0.508, "coin flip ", color=COLOR_CRITICAL, fontsize=8, ha="right",
            transform=ax.get_yaxis_transform())
    ax.set_xticks(np.arange(len(mods)) + w / 2)
    ax.set_xticklabels(mods, rotation=30, ha="right")
    ax.set_ylabel("accuracy calling the winner\n(pairs with |log2 ratio| > 1)")
    ax.set_ylim(0, 1); ax.legend(frameon=False)
    ax.set_title("Which strain wins -- the question a bench scientist actually asks")
    fig.tight_layout()
    fig.savefig(fd / "f04_sign_accuracy.png", dpi=160)
    plt.close(fig)

    return sorted(p.name for p in fd.glob("*.png"))
