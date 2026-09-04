"""Would an OD-based target be worth redoing the genomic ML on? Feasibility test.

Three candidate targets, and they are not interchangeable:

  relative abundance (current)  -- who wins. ANTISYMMETRIC. 1465 pairs / 74 strains.
  total yield (OD)              -- how much biomass the pair makes together. SYMMETRIC, a
                                   different biological question, not a substitute. 3916 pairs
                                   / 98 strains.
  per-strain absolute (OD x rel.abund.) -- the literal reading of "absolute abundance", but it
                                   multiplies two measurements, so it is limited to the
                                   sequenced subset (no gain in n) and compounds their noise.

Because total yield is symmetric, the whole antisymmetric apparatus in genomic_ml.py does not
apply: no (f(a,b)-f(b,a))/2, and the antisymmetric feature block is replaced by |x_a - x_b|.
The per-strain "strength" decomposition becomes an additive per-strain YIELD model, refit
in-fold, which is the baseline to beat exactly as before.

Evaluated under the same two regimes, because the same trap applies: under cv_pair a model can
score by recognising the strain, not by understanding the genome.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
EXPS = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(EXPS / "shared_pipelines"))
from replicate_report import load_od, CONFIG
import genomic_ml as gm

OUT = HERE / "outputs"
N_PCA, N_SPLITS, N_REPEATS = 30, 5, 3


def build_pairs():
    lay, M = load_od("20260630")
    z = np.full(len(lay), np.nan)
    col = M[:, 25]
    for p in lay.dest_plate.unique():
        m = (lay.dest_plate == p).values
        z[m] = (col[m] - np.nanmean(col[m])) / np.nanstd(col[m])
    lay["od_z"] = z
    lay = lay[lay.well_type == "pair"].dropna(subset=["od_z"])
    lay["pair_key"] = [frozenset((a, b)) for a, b in zip(lay.strain1, lay.strain2)]
    g = lay.groupby("pair_key").agg(target=("od_z", "mean"), n_wells=("od_z", "size")).reset_index()
    g = g[g.pair_key.map(len) == 2]
    g["strain_a"] = [sorted(k)[0] for k in g.pair_key]
    g["strain_b"] = [sorted(k)[1] for k in g.pair_key]
    g["pair_id"] = [f"{a}|{b}" for a, b in zip(g.strain_a, g.strain_b)]

    w2s, _ = gm.load_strain_mapping(_cfg())
    ko, _ = gm.load_ko_table(_cfg())
    g = g[g.strain_a.map(w2s).isin(ko.index) & g.strain_b.map(w2s).isin(ko.index)].reset_index(drop=True)
    return g, ko, w2s


class _Cfg:
    genomic_dir = gm.GENOMIC_DIR
    feature_table = "KEGG_ko_and_strains_table.csv"
    mapping_csv = "mapping_384_well_plate_collection.csv"
    min_prevalence = 4

    @property
    def feature_path(self):
        return self.genomic_dir / self.feature_table

    @property
    def mapping_path(self):
        return self.genomic_dir / self.mapping_csv


def _cfg():
    return _Cfg()


def features(pairs, ko, w2s):
    strains = sorted(set(pairs.strain_a) | set(pairs.strain_b))
    sub = ko.loc[[w2s[s] for s in strains]]
    pa = (sub > 0).astype(float)
    prev = pa.sum(0)
    n = len(strains)
    keep = (prev >= 4) & (prev <= n - 4)
    X = pa.loc[:, keep]
    X.index = strains
    return X


def additive_yield(pairs, idx, strains, lam=1.0):
    """Pair yield ~ s_a + s_b, refit in-fold. The no-genomics baseline for a symmetric target."""
    pos = {s: i for i, s in enumerate(strains)}
    sub = pairs.iloc[idx]
    D = np.zeros((len(sub), len(strains)))
    for r, (a, b) in enumerate(zip(sub.strain_a, sub.strain_b)):
        D[r, pos[a]] += 1
        D[r, pos[b]] += 1
    y = sub.target.values
    s = np.linalg.lstsq(D.T @ D + lam * np.eye(len(strains)), D.T @ y, rcond=None)[0]
    return pd.Series(s, index=strains)


def design(Z, pairs):
    """SYMMETRIC pair features: sum and absolute difference. No antisymmetric block -- the
    target does not change sign when the two strains are swapped."""
    A = Z.loc[pairs.strain_a].values
    B = Z.loc[pairs.strain_b].values
    return np.hstack([A + B, np.abs(A - B)])


def run():
    pairs, ko, w2s = build_pairs()
    X = features(pairs, ko, w2s)
    strains = sorted(set(pairs.strain_a) | set(pairs.strain_b))
    print(f"OD target: {len(pairs)} pairs, {len(strains)} strains, {X.shape[1]} KOs\n")

    rows = []
    for regime in ("cv_pair", "cv_strain"):
        for rep in range(N_REPEATS):
            for tr, te in gm.make_folds(pairs, regime, n_splits=N_SPLITS, seed=rep * 97):
                tr_s = sorted(set(pairs.iloc[tr].strain_a) | set(pairs.iloc[tr].strain_b))
                k = min(N_PCA, len(tr_s) - 1, X.shape[1])
                pca = PCA(n_components=k, random_state=0).fit(X.loc[tr_s].values)
                Zv = pca.transform(X.values)
                sc = StandardScaler().fit(pca.transform(X.loc[tr_s].values))
                Z = pd.DataFrame(sc.transform(Zv), index=X.index)
                D = design(Z, pairs)
                ytr, yte = pairs.target.values[tr], pairs.target.values[te]

                preds = {}
                preds["ridge"] = Ridge(alpha=50).fit(D[tr], ytr).predict(D[te])
                try:
                    from xgboost import XGBRegressor
                    preds["xgboost"] = XGBRegressor(
                        n_estimators=400, max_depth=3, learning_rate=0.05, subsample=0.8,
                        colsample_bytree=0.5, reg_lambda=5.0, min_child_weight=5,
                        n_jobs=4, random_state=0, verbosity=0).fit(D[tr], ytr).predict(D[te])
                except Exception:
                    pass
                s = additive_yield(pairs, tr, tr_s)
                te_df = pairs.iloc[te]
                preds["additive_yield_no_genomics"] = (te_df.strain_a.map(s).fillna(s.mean())
                                                       + te_df.strain_b.map(s).fillna(s.mean())).values
                rng = np.random.default_rng(1000 + rep)
                preds["SHUFFLED_ridge"] = Ridge(alpha=50).fit(D[tr], rng.permutation(ytr)).predict(D[te])

                for m, p in preds.items():
                    ss = np.sum((yte - yte.mean()) ** 2)
                    rows.append({"regime": regime, "model": m,
                                 "r2": 1 - np.sum((yte - p) ** 2) / ss,
                                 "spearman": spearmanr(yte, p)[0] if np.std(p) > 0 else np.nan,
                                 "n": len(te)})
    df = pd.DataFrame(rows)
    summ = (df.groupby(["regime", "model"])
              .agg(r2_mean=("r2", "mean"), r2_sd=("r2", "std"),
                   rho_mean=("spearman", "mean"), n_test=("n", "sum"))
              .reset_index().sort_values(["regime", "r2_mean"], ascending=[True, False]))
    summ.to_csv(OUT / "od_target_ml_summary.csv", index=False)
    print(summ.round(3).to_string(index=False))
    return summ


if __name__ == "__main__":
    run()
