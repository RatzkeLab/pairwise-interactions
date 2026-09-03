"""Predict a pair's TOTAL YIELD (plate-reader OD) from the two genomes.

Companion to genomic_ml.py, which predicts the *competitive outcome* (who wins, from 16S
relative abundance). These are orthogonal halves of a pairwise result and neither replaces the
other:

    genomic_ml.py         log2(a/b)   ANTISYMMETRIC   who wins
    this module           OD of well  SYMMETRIC       how much biomass the pair makes together

Because the target is symmetric, essentially none of genomic_ml's antisymmetric apparatus
carries over. The `(f(a,b) - f(b,a))/2` antisymmetrisation is wrong here and is not used; the
antisymmetric feature block `x_a - x_b` is replaced by `|x_a - x_b|`; and the per-strain
Bradley-Terry "strength" decomposition becomes an additive per-strain YIELD contribution. What
IS shared -- the strain-ID join, the KO feature matrix, the two CV regimes, the 16S control --
is imported rather than reimplemented.

Why bother when relative abundance already works:

  **OD is independent of 16S resolvability.** Every pair genomic_ml must discard as
  `high_uncertainty` (references too close to tell the strains apart, label forced to ~50/50)
  is measured perfectly well by a plate reader. That recovers ~2400 pairs for 20260630 and
  stops degenerate-16S strains being second-class. This, not the raw pair count, is the reason
  to run it -- the binding constraint on the other analysis is the number of GENOMES (74 -> 88
  here, +19%), not the number of pairs.

Reliability, measured from cross-plate technical replicates (see technical_replicates/):
OD600 rho = 0.869 for 20260630, against 0.928 for sequencing relative abundance. Slightly
noisier, far more plentiful.

The same evaluation trap applies and is handled the same way: under `cv_pair` a model can score
by recognising a strain rather than understanding a genome, so `cv_strain` (both strains of a
test pair unseen) is the honest number. For this target the no-genomics additive baseline is
especially strong under cv_pair -- expect it to WIN there. That is the point of reporting it.

Call order: build_dataset -> strain_feature_matrix -> replicate_reliability ->
cross_validate -> genome_to_yield -> make_all_figures.
"""

import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from genomic_ml import (COLOR_BLUE, COLOR_RED, COLOR_GRID, COLOR_TEXT_SECONDARY, COLOR_GOOD,
                        COLOR_CRITICAL, GENOMIC_DIR, load_strain_mapping, load_ko_table,
                        make_folds, _scaled)

OD600_INDEX = 25          # 61-point spectra, 350 nm + 10 nm steps


@dataclass
class YieldMLConfig:
    """Everything specific to running the yield analysis for one experiment."""
    exp_cfg: object                  # ExperimentConfig -- supplies exp_base
    layout_csv: Path                 # strain layout (dest_plate/dest_well -> strain1/strain2)
    od_dir: Path                     # folder of full-spectrum plate-reader csvs, one per plate
    genomic_dir: Path = GENOMIC_DIR
    feature_table: str = "KEGG_ko_and_strains_table.csv"
    mapping_csv: str = "mapping_384_well_plate_collection.csv"
    min_prevalence: int = 4
    n_pca: int = 30
    wavelength_index: int = OD600_INDEX
    out_dir: Path = None

    def __post_init__(self):
        if self.out_dir is None:
            self.out_dir = self.exp_cfg.exp_base / "analysis" / "genomic_ml_yield" / "outputs"
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
# 01 -- dataset: plate-reader wells -> one yield value per strain pair
# ===========================================================================

def _read_plate(path, wl_index):
    txt = Path(path).read_text(encoding="latin-1").splitlines()
    plate = next((l.split("ID1:")[1].split("ID2")[0].strip() for l in txt[:6] if "ID1:" in l), None)

    def num(x):
        try:
            return float(x.strip())
        except ValueError:
            return np.nan
    rows = [[num(v) for v in l.split(",")] for l in txt if re.match(r"^\s*[\d.\-\s]+,", l)]
    rows = [r for r in rows if len(r) == 24]
    a = np.array(rows)
    if a.shape[0] < (wl_index + 1) * 16:
        return plate, None
    return plate, a[wl_index * 16:(wl_index + 1) * 16]


def build_dataset(gcfg):
    """One row per unique tested strain pair, target = mean plate-z-scored OD.

    Two normalisations that matter:
      - **z-scored within destination plate.** Plates differ in overall level for reasons that
        are not biology (inoculum age, evaporation, reader drift); without this, plate identity
        leaks into the target.
      - **mean over replicate wells**, which sit on different plates 96% of the time, so the
        per-pair value already averages out plate position.

    Mono wells are dropped (they are the per-strain reference, not a pair), as are pairs where
    either strain lacks an annotated genome.
    """
    grids = {}
    for f in sorted(Path(gcfg.od_dir).glob("*.csv")):
        p, g = _read_plate(f, gcfg.wavelength_index)
        if g is not None and p and p.isdigit():
            grids.setdefault(int(p), []).append(g)
    grids = {k: np.nanmean(v, axis=0) for k, v in grids.items()}

    lay = pd.read_csv(gcfg.layout_csv)
    od = []
    for r in lay.itertuples():
        g = grids.get(int(r.dest_plate))
        row = ord(str(r.dest_row).upper()[0]) - 65
        col = int(r.dest_col) - 1
        od.append(g[row, col] if (g is not None and 0 <= row < 16 and 0 <= col < 24) else np.nan)
    lay["od"] = od
    lay["od_z"] = lay.groupby("dest_plate")["od"].transform(lambda s: (s - s.mean()) / s.std(ddof=0))

    steps = [("layout_wells", len(lay)), ("wells_with_od", int(lay.od.notna().sum()))]
    w = lay[(lay.well_type == "pair") & lay.od_z.notna()].copy()
    steps.append(("pair_wells_with_od", len(w)))
    w["pair_key"] = [frozenset((a, b)) for a, b in zip(w.strain1, w.strain2)]
    w = w[w.pair_key.map(len) == 2]

    g = (w.groupby("pair_key")
           .agg(target=("od_z", "mean"), sd_target=("od_z", "std"), n_wells=("od_z", "size"),
                mean_raw_od=("od", "mean"))
           .reset_index())
    steps.append(("distinct_pairs", len(g)))

    w2s, _ = load_strain_mapping(gcfg)
    ko, empty = load_ko_table(gcfg)
    g["strain_a"] = [sorted(k)[0] for k in g.pair_key]
    g["strain_b"] = [sorted(k)[1] for k in g.pair_key]
    g["genome_a"] = g.strain_a.map(w2s)
    g["genome_b"] = g.strain_b.map(w2s)
    g = g[g.genome_a.isin(ko.index) & g.genome_b.isin(ko.index)].reset_index(drop=True)
    g["pair_id"] = [f"{a}|{b}" for a, b in zip(g.strain_a, g.strain_b)]
    steps.append(("pairs_with_both_genomes", len(g)))

    strains = sorted(set(g.strain_a) | set(g.strain_b))
    steps += [("n_strains", len(strains)), ("n_plates_read", len(grids)),
              ("n_genomes_dropped_all_zero_ko", len(empty))]
    summary = pd.DataFrame(steps, columns=["step", "value"])
    summary.to_csv(gcfg.out_dir / "y01_dataset_summary.csv", index=False)
    g.drop(columns=["pair_key"]).to_csv(gcfg.out_dir / "y01_pairs.csv", index=False)
    return g, summary, w


def strain_feature_matrix(gcfg, pairs):
    """KO presence/absence per strain, prevalence-filtered, indexed by well label."""
    w2s, _ = load_strain_mapping(gcfg)
    ko, _ = load_ko_table(gcfg)
    strains = sorted(set(pairs.strain_a) | set(pairs.strain_b))
    sub = ko.loc[[w2s[s] for s in strains]]
    pa = (sub > 0).astype(float)
    n = len(strains)
    prev = pa.sum(axis=0)
    keep = (prev >= gcfg.min_prevalence) & (prev <= n - gcfg.min_prevalence)
    X = pa.loc[:, keep]
    X.index = strains
    X.columns = [c.replace("ko:", "") for c in X.columns]
    summ = pd.DataFrame({"n_ko_present": pa.sum(axis=1).values,
                         "total_gene_copies": sub.sum(axis=1).values}, index=strains)
    return X, summ


def replicate_reliability(gcfg, wells):
    """Split-half agreement over replicate wells -- the ceiling any model is measured against.

    Replicates sit on different destination plates 96% of the time, so this is whole-pipeline
    reproducibility, not well-to-well noise. R^2 cannot exceed roughly this value.
    """
    wells = wells.copy()
    wells["pair_key"] = [frozenset((a, b)) for a, b in zip(wells.strain1, wells.strain2)]
    g = wells.groupby("pair_key")["od_z"].agg(list)
    g = g[g.map(len) >= 2]
    x = np.array([v[0] for v in g])
    y = np.array([np.mean(v[1:]) for v in g])
    out = pd.DataFrame([{"n_pairs_with_replicates": len(g),
                         "spearman": float(spearmanr(x, y)[0]),
                         "pearson": float(pearsonr(x, y)[0]),
                         "median_abs_diff": float(np.median(np.abs(x - y))),
                         "ceiling_r2_estimate": float(pearsonr(x, y)[0])}])
    out.to_csv(gcfg.out_dir / "y02_replicate_reliability.csv", index=False)
    return out


# ===========================================================================
# 02 -- features and models for a SYMMETRIC target
# ===========================================================================

def _design(Z, S, pairs, extra=None):
    """Pair features that do not change when the two strains are swapped.

    `Z_a + Z_b` carries "what is in this well between them"; `|Z_a - Z_b|` carries "how
    different are they", which is where complementarity vs redundancy would live. A signed
    difference would be wrong: the target is unchanged by swapping a and b, so a feature that
    flips sign can only add noise.
    """
    A, B = Z.loc[pairs.strain_a].values, Z.loc[pairs.strain_b].values
    sa, sb = S.loc[pairs.strain_a].values, S.loc[pairs.strain_b].values
    blocks = [A + B, np.abs(A - B), sa + sb, np.abs(sa - sb)]
    if extra:
        blocks.append(np.column_stack([np.asarray(e) for e in extra]))
    names = ([f"sum_{c}" for c in Z.columns] + [f"absdiff_{c}" for c in Z.columns]
             + [f"sum_{c}" for c in S.columns] + [f"absdiff_{c}" for c in S.columns])
    return np.hstack(blocks), names


def _metrics(y, yhat):
    y, yhat = np.asarray(y, float), np.asarray(yhat, float)
    ok = np.isfinite(y) & np.isfinite(yhat)
    y, yhat = y[ok], yhat[ok]
    # conventional R^2 against the mean -- NOTE this differs from genomic_ml.py, where the
    # antisymmetric target makes "predict zero" the right reference
    ss = np.sum((y - y.mean()) ** 2)
    out = {"n": len(y), "r2": float(1 - np.sum((y - yhat) ** 2) / ss) if ss > 0 else np.nan,
           "mae": float(np.mean(np.abs(y - yhat)))}
    if len(y) < 3 or np.std(yhat) == 0:
        return {**out, "pearson_r": np.nan, "spearman_rho": np.nan}
    return {**out, "pearson_r": float(pearsonr(y, yhat)[0]),
            "spearman_rho": float(spearmanr(y, yhat)[0])}


def _fit_additive_yield(pairs, idx, strains, ridge_lambda=1.0):
    """yield(a,b) ~ c_a + c_b, fit in-fold. The no-genomics baseline.

    Under cv_pair this is expected to WIN: knowing each strain's typical contribution is simply
    a better predictor of a co-culture's biomass than any genome model, once you have seen that
    strain. It is reported precisely so that a genomic model is never credited for it.
    """
    pos = {s: i for i, s in enumerate(strains)}
    sub = pairs.iloc[idx]
    D = np.zeros((len(sub), len(strains)))
    for r, (a, b) in enumerate(zip(sub.strain_a, sub.strain_b)):
        D[r, pos[a]] += 1
        D[r, pos[b]] += 1
    y = sub.target.values
    c = np.linalg.lstsq(D.T @ D + ridge_lambda * np.eye(len(strains)), D.T @ y, rcond=None)[0]
    return pd.Series(c, index=strains)


class _Ctx:
    """One fold: PCA and scalers fit on training strains only."""
    def __init__(self, gcfg, pairs, X, summ, train_idx, test_idx, phylo=None):
        self.gcfg, self.pairs, self.X = gcfg, pairs, X
        self.train_idx, self.test_idx = train_idx, test_idx
        self.train_strains = sorted(set(pairs.iloc[train_idx].strain_a)
                                    | set(pairs.iloc[train_idx].strain_b))
        self.y_train = pairs.target.values[train_idx]
        self.y_test = pairs.target.values[test_idx]

        k = min(gcfg.n_pca, len(self.train_strains) - 1, X.shape[1])
        pca = PCA(n_components=k, random_state=0).fit(X.loc[self.train_strains].values)
        self.Z = _scaled(pd.DataFrame(pca.transform(X.values), index=X.index,
                                      columns=[f"pc{i+1}" for i in range(k)]), self.train_strains)
        sc = StandardScaler().fit(summ.loc[self.train_strains].values)
        self.S = pd.DataFrame(sc.transform(summ.values), index=summ.index, columns=summ.columns)

        self.Zp = None
        if phylo is not None:
            L = phylo.loc[:, [s for s in self.train_strains if s in phylo.columns]]
            L = L.loc[[s for s in X.index if s in L.index]]
            if L.shape[1] > 2:
                kp = min(gcfg.n_pca, L.shape[1] - 1)
                pp = PCA(n_components=kp, random_state=0).fit(
                    L.loc[[s for s in self.train_strains if s in L.index]].values)
                Zp = pd.DataFrame(pp.transform(L.values), index=L.index,
                                  columns=[f"ph{i+1}" for i in range(kp)])
                self.Zp = _scaled(Zp, [s for s in self.train_strains if s in Zp.index])

        self.D, self.names = _design(self.Z, self.S, pairs)
        self._raw = None

    @property
    def raw(self):
        """Un-reduced KO design, built on first use only (it is ~10k columns)."""
        if self._raw is None:
            self._raw = _design(self.X, self.S, self.pairs)[0]
        return self._raw


def m_mean(ctx):
    return np.full(len(ctx.test_idx), ctx.y_train.mean())


def m_additive_no_genomics(ctx):
    c = _fit_additive_yield(ctx.pairs, ctx.train_idx, ctx.train_strains)
    te = ctx.pairs.iloc[ctx.test_idx]
    m = c.mean()
    return (te.strain_a.map(c).fillna(m) + te.strain_b.map(c).fillna(m)).values


def m_ridge(ctx, alpha=50.0):
    return Ridge(alpha=alpha).fit(ctx.D[ctx.train_idx], ctx.y_train).predict(ctx.D[ctx.test_idx])


def m_xgb(ctx, raw=False, **kw):
    from xgboost import XGBRegressor
    D = ctx.raw if raw else ctx.D
    p = dict(n_estimators=400, max_depth=3, learning_rate=0.05, subsample=0.8,
             colsample_bytree=0.2 if raw else 0.5, reg_lambda=5.0, min_child_weight=5,
             n_jobs=4, random_state=0, verbosity=0)
    p.update(kw)
    return XGBRegressor(**p).fit(D[ctx.train_idx], ctx.y_train).predict(D[ctx.test_idx])


def m_xgb_raw(ctx):
    return m_xgb(ctx, raw=True)


def m_ridge_phylo16s(ctx, alpha=50.0):
    """Taxonomy-only control. KO content and 16S are both largely phylogenetic, so what KO buys
    OVER this row is the part that is about gene content rather than the taxon label."""
    if ctx.Zp is None:
        raise RuntimeError("no 16S distance matrix supplied")
    common = [s for s in ctx.Zp.index]
    sub = ctx.pairs[ctx.pairs.strain_a.isin(common) & ctx.pairs.strain_b.isin(common)]
    if len(sub) < 50:
        raise RuntimeError("too few pairs with 16S")
    Dp, _ = _design(ctx.Zp, ctx.S.loc[common], sub)
    idx = {v: i for i, v in enumerate(sub.index)}
    tr = [idx[i] for i in ctx.pairs.index[ctx.train_idx] if i in idx]
    te = [idx[i] for i in ctx.pairs.index[ctx.test_idx] if i in idx]
    if len(tr) < 30 or len(te) < 5:
        raise RuntimeError("fold too small after 16S restriction")
    mdl = Ridge(alpha=alpha).fit(Dp[tr], sub.target.values[tr])
    pred = np.full(len(ctx.test_idx), np.nan)
    hit = [j for j, i in enumerate(ctx.pairs.index[ctx.test_idx]) if i in idx]
    pred[hit] = mdl.predict(Dp[te])
    return pred


def _two_stage(ctx, learner, Z=None):
    """genome -> per-strain yield contribution -> pair yield as the sum.

    Better posed than predicting pairs directly: one row per STRAIN (the real independent-sample
    count), one target, and the resulting c_a + c_b is symmetric by construction.
    """
    Z = ctx.Z if Z is None else Z
    c = _fit_additive_yield(ctx.pairs, ctx.train_idx, ctx.train_strains)
    A = np.hstack([Z.loc[ctx.train_strains].values, ctx.S.loc[ctx.train_strains].values])
    mdl = learner().fit(A, c.loc[ctx.train_strains].values)
    allA = np.hstack([Z.values, ctx.S.loc[Z.index].values])
    chat = pd.Series(mdl.predict(allA), index=Z.index)
    te = ctx.pairs.iloc[ctx.test_idx]
    return (te.strain_a.map(chat) + te.strain_b.map(chat)).values


def m_two_stage_ridge(ctx):
    return _two_stage(ctx, lambda: Ridge(alpha=50.0))


def m_two_stage_xgb(ctx):
    from xgboost import XGBRegressor
    return _two_stage(ctx, lambda: XGBRegressor(
        n_estimators=300, max_depth=2, learning_rate=0.05, subsample=0.8, colsample_bytree=0.5,
        reg_lambda=5.0, min_child_weight=3, n_jobs=4, random_state=0, verbosity=0))


MODELS = {
    "mean_baseline": m_mean,
    "additive_yield_no_genomics": m_additive_no_genomics,
    "ridge_pca": m_ridge,
    "xgboost_pca": m_xgb,
    "xgboost_raw_ko": m_xgb_raw,
    "ridge_phylo16s_only": m_ridge_phylo16s,
    "two_stage_ridge": m_two_stage_ridge,
    "two_stage_xgboost": m_two_stage_xgb,
}


# ===========================================================================
# 03 -- evaluation
# ===========================================================================

def summarize_folds(fold_df):
    num = ["pearson_r", "spearman_rho", "r2", "mae"]
    return (fold_df.groupby(["regime", "model"])
            .agg(n_folds=("fold", "nunique"), n_test=("n", "sum"),
                 **{f"{c}_mean": (c, "mean") for c in num},
                 **{f"{c}_sd": (c, "std") for c in num})
            .reset_index().sort_values(["regime", "r2_mean"], ascending=[True, False]))


def cross_validate(gcfg, pairs, X, summ, models=None, regimes=("cv_pair", "cv_strain"),
                   n_splits=5, seed=0, n_repeats=3, phylo=None, shuffle_control=True):
    """Both regimes, because only one of them is a real test.

    cv_pair  -- unseen PAIR, both strains seen elsewhere. Answers "fill in the pairs we never
                tested". A model can score here by recognising the strain, so a good number is
                NOT evidence that gene content predicts yield. The additive no-genomics
                baseline is the thing to beat, and usually wins.
    cv_strain -- unseen STRAINS, both held out, pairs straddling the split dropped entirely.
                Answers "sequence a new isolate; predict its co-culture yield". The additive
                baseline is undefined here by construction and collapses; that collapse is the
                point.
    """
    models = models or [m for m in MODELS]
    per_fold, per_pair = [], []
    for regime in regimes:
        folds = [(rep, tr, te) for rep in range(n_repeats)
                 for tr, te in make_folds(pairs, regime, n_splits=n_splits, seed=seed + 97 * rep)]
        for fi, (rep, tr, te) in enumerate(folds):
            ctx = _Ctx(gcfg, pairs, X, summ, tr, te, phylo=phylo)
            for name in models:
                try:
                    yhat = MODELS[name](ctx)
                except Exception as e:
                    continue
                per_fold.append({"regime": regime, "fold": fi, "model": name,
                                 "n_train": len(tr), **_metrics(ctx.y_test, yhat)})
                per_pair.append(pd.DataFrame({
                    "regime": regime, "fold": fi, "model": name,
                    "pair_id": pairs.pair_id.values[te],
                    "y_true": ctx.y_test, "y_pred": yhat}))
            if shuffle_control:
                rng = np.random.default_rng(1000 + fi)
                sh = pairs.copy()
                sh.iloc[tr, sh.columns.get_loc("target")] = rng.permutation(pairs.target.values[tr])
                ctx_s = _Ctx(gcfg, sh, X, summ, tr, te, phylo=phylo)
                for name in ("ridge_pca", "xgboost_pca", "two_stage_ridge"):
                    if name not in models:
                        continue
                    try:
                        per_fold.append({"regime": regime, "fold": fi,
                                         "model": f"SHUFFLED_{name}", "n_train": len(tr),
                                         **_metrics(ctx.y_test, MODELS[name](ctx_s))})
                    except Exception:
                        pass
    fold_df = pd.DataFrame(per_fold)
    pair_df = pd.concat(per_pair, ignore_index=True) if per_pair else pd.DataFrame()
    summ_df = summarize_folds(fold_df)
    fold_df.to_csv(gcfg.out_dir / "y03_cv_per_fold.csv", index=False)
    summ_df.to_csv(gcfg.out_dir / "y03_cv_summary.csv", index=False)
    pair_df.to_csv(gcfg.out_dir / "y03_cv_predictions.csv", index=False)
    return summ_df, fold_df, pair_df


def genome_to_yield(gcfg, pairs, X, summ, phylo=None, n_splits=10, seed=0, n_repeats=3):
    """Leave-strains-out: predict a strain's yield contribution from its genome alone.

    The cleanest form of the question -- one row per strain, one target, no pair structure to
    inflate anything -- and the fair place to compare KO content against 16S taxonomy, since
    both face an identical regression.
    """
    all_strains = sorted(set(pairs.strain_a) | set(pairs.strain_b))
    c_full = _fit_additive_yield(pairs, np.arange(len(pairs)), all_strains)
    sources = {"ko_content": "ko"}
    if phylo is not None:
        sources.update({"phylo16s_only": "phylo", "ko_plus_phylo16s": "both"})

    rows, metrics = [], []
    for src_name, src in sources.items():
        # only 74 of 88 strains have a 16S consensus, so the phylo sources are run on that
        # subset rather than skipping folds -- otherwise the KO-vs-taxonomy comparison, the
        # single most important control here, silently disappears from the output
        strains = ([s for s in all_strains if s in phylo.index and s in phylo.columns]
                   if src in ("phylo", "both") else all_strains)
        acc = {}
        for rep in range(n_repeats):
            rng = np.random.default_rng(seed + 97 * rep)
            fold_of = dict(zip(strains, rng.permutation(len(strains)) % n_splits))
            for k in range(n_splits):
                tr_s = [s for s in strains if fold_of[s] != k]
                te_s = [s for s in strains if fold_of[s] == k]
                if not te_s:
                    continue
                keep = pairs.strain_a.isin(tr_s) & pairs.strain_b.isin(tr_s)
                if keep.sum() < 30:
                    continue
                c_tr = _fit_additive_yield(pairs, np.where(keep)[0], tr_s)
                sc = StandardScaler().fit(summ.loc[tr_s].values)
                blocks = []
                if src in ("ko", "both"):
                    pk = PCA(n_components=min(gcfg.n_pca, len(tr_s) - 1, X.shape[1]),
                             random_state=0).fit(X.loc[tr_s].values)
                    sk = StandardScaler().fit(pk.transform(X.loc[tr_s].values))
                    blocks.append(lambda ss, pk=pk, sk=sk: sk.transform(pk.transform(X.loc[ss].values)))
                if src in ("phylo", "both"):
                    cols = [s for s in tr_s if s in phylo.columns]
                    if len(cols) < 3:
                        continue
                    L = phylo.loc[strains, cols]
                    pp = PCA(n_components=min(gcfg.n_pca, len(cols) - 1),
                             random_state=0).fit(L.loc[[s for s in tr_s if s in L.index]].values)
                    sp = StandardScaler().fit(pp.transform(L.loc[[s for s in tr_s if s in L.index]].values))
                    blocks.append(lambda ss, pp=pp, sp=sp, L=L: sp.transform(pp.transform(L.loc[ss].values)))
                if not blocks:
                    continue
                A = lambda ss: np.hstack([b(ss) for b in blocks] + [sc.transform(summ.loc[ss].values)])
                try:
                    mdl = Ridge(alpha=50.0).fit(A(tr_s), c_tr.loc[tr_s].values)
                    for s_, v in zip(te_s, mdl.predict(A(te_s))):
                        acc.setdefault(s_, []).append(v)
                except Exception:
                    continue
        if not acc:
            continue
        pred = {k: float(np.mean(v)) for k, v in acc.items()}
        df = pd.DataFrame({"strain": list(pred), "feature_source": src_name,
                           "yield_observed": [c_full[s] for s in pred],
                           "yield_predicted_from_genome": [pred[s] for s in pred]})
        rows.append(df)
        metrics.append({"feature_source": src_name,
                        **_metrics(df.yield_observed, df.yield_predicted_from_genome)})

    out = pd.concat(rows, ignore_index=True)
    met = pd.DataFrame(metrics)
    sv = c_full.loc[X.index].values
    corr = np.array([spearmanr(X.iloc[:, j].values, sv)[0] if X.iloc[:, j].std() > 0 else np.nan
                     for j in range(X.shape[1])])
    ko_df = (pd.DataFrame({"ko": X.columns, "spearman_vs_yield": corr,
                           "prevalence": X.sum(axis=0).values})
             .dropna().sort_values("spearman_vs_yield", key=np.abs, ascending=False))
    out.to_csv(gcfg.out_dir / "y04_yield_from_genome.csv", index=False)
    met.to_csv(gcfg.out_dir / "y04_yield_from_genome_metrics.csv", index=False)
    ko_df.to_csv(gcfg.out_dir / "y04_ko_yield_correlations.csv", index=False)
    return out, met, ko_df


# ===========================================================================
# 04 -- figures
# ===========================================================================

def make_all_figures(gcfg, summ_df, pair_df, yield_df, reliability):
    fd = gcfg.fig_dir
    ceiling = float(reliability.ceiling_r2_estimate.iloc[0]) if len(reliability) else np.nan

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8))
    for ax, regime, title in zip(axes, ["cv_pair", "cv_strain"],
                                 ["cv_pair: unseen PAIR, strains seen\n(a per-strain lookup is available)",
                                  "cv_strain: unseen STRAINS\n(no lookup possible -- the honest test)"]):
        d = summ_df[summ_df.regime == regime].sort_values("r2_mean")
        if not len(d):
            continue
        cols = [COLOR_CRITICAL if m.startswith("SHUFFLED") else
                (COLOR_TEXT_SECONDARY if m in ("mean_baseline", "additive_yield_no_genomics")
                 else COLOR_BLUE) for m in d.model]
        lo = -1.0
        v = d.r2_mean.values
        clip = v < lo
        ax.barh(d.model, np.where(clip, lo, v), color=cols,
                xerr=np.where(clip, 0, d.r2_sd.values),
                error_kw=dict(ecolor=COLOR_TEXT_SECONDARY, lw=1, capsize=3))
        for i, (val, c) in enumerate(zip(v, clip)):
            if c:
                ax.text(lo + 0.02, i, f"R² = {val:.1f} ◄", va="center", fontsize=8,
                        color="white", fontweight="bold")
        ax.axvline(0, color="black", lw=1)
        if np.isfinite(ceiling):
            ax.axvline(ceiling, color=COLOR_GOOD, lw=1.2, ls="--")
            ax.text(ceiling, len(d) - 0.4, " replicate ceiling", color=COLOR_GOOD, fontsize=7.5,
                    rotation=90, va="top")
        ax.set_xlim(lo, 1.02)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("R² (vs. predicting the mean), mean ± sd over folds")
    fig.suptitle("Predicting a pair's total yield (OD) from KEGG KO content", fontweight="bold")
    fig.tight_layout()
    fig.savefig(fd / "fy01_model_comparison.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    for ax, regime in zip(axes, ["cv_pair", "cv_strain"]):
        d = summ_df[(summ_df.regime == regime) & (~summ_df.model.str.startswith("SHUFFLED"))
                    & (~summ_df.model.isin(["mean_baseline", "additive_yield_no_genomics"]))]
        if not len(d) or not len(pair_df):
            continue
        best = d.iloc[0].model
        p = pair_df[(pair_df.regime == regime) & (pair_df.model == best)]
        ax.scatter(p.y_pred, p.y_true, s=10, alpha=.3, color=COLOR_BLUE, edgecolor="none")
        lim = np.nanpercentile(np.abs(p.y_true), 99) * 1.1
        ax.plot([-lim, lim], [-lim, lim], color=COLOR_TEXT_SECONDARY, lw=1, ls="--")
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_title(f"{regime} — {best}\nSpearman ρ = {spearmanr(p.y_pred, p.y_true)[0]:.2f}, "
                     f"n = {len(p)}", fontsize=10)
        ax.set_xlabel("predicted yield (plate z)"); ax.set_ylabel("observed yield (plate z)")
    fig.tight_layout()
    fig.savefig(fd / "fy02_predicted_vs_observed.png", dpi=160)
    plt.close(fig)

    if len(yield_df):
        srcs = list(dict.fromkeys(yield_df.feature_source))
        fig, axes = plt.subplots(1, len(srcs), figsize=(5.2 * len(srcs), 5.2), squeeze=False)
        for ax, src in zip(axes[0], srcs):
            d = yield_df[yield_df.feature_source == src]
            ax.scatter(d.yield_predicted_from_genome, d.yield_observed, s=26, color=COLOR_BLUE,
                       alpha=.7, edgecolor="white", linewidth=.6)
            ax.axhline(0, color=COLOR_GRID, lw=1); ax.axvline(0, color=COLOR_GRID, lw=1)
            ax.set_xlabel("predicted (held-out strain)")
            ax.set_ylabel("observed yield contribution")
            ax.set_title(f"{src}\nSpearman ρ = "
                         f"{spearmanr(d.yield_predicted_from_genome, d.yield_observed)[0]:.2f}, "
                         f"n = {len(d)}", fontsize=10)
        fig.suptitle("Can a genome alone say how much a strain contributes to co-culture yield?",
                     fontweight="bold")
        fig.tight_layout()
        fig.savefig(fd / "fy03_yield_from_genome.png", dpi=160)
        plt.close(fig)
    return sorted(p.name for p in fd.glob("*.png"))
