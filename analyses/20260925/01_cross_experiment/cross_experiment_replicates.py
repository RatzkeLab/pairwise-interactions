"""How well does one experiment predict the other? 20260630 (merged) vs 20260907.

Two ways in, because only ~100 pairs were tested in both experiments:

  x01  identity gate: are the shared strain names the same organisms? Each experiment's own
       16S consensus, compared over the overlap (strain_identity.flexible_distance).
  x02  DIRECT replicates: pairs tested in both. One experiment's pair mean predicts the
       other's, with no fitting (R^2 against the identity line) and with a linear fit (r^2).
  x03  the yardstick: the same statistics for replicate wells WITHIN each experiment, and
       for single wells ACROSS experiments, so across-experiment agreement can be read
       against what re-running the same experiment gives.
  x04  STRENGTH TRANSFER: fit per-strain Bradley-Terry strengths on one experiment, predict
       every scoreable pair of the other (predicted log2 ratio = (s_a - s_b) / ln 2). Uses all
       the data, not just the overlap; pairs shared by both are also reported separately.

Targets: log2(a/b) with pseudocount 0.5 (the ML target) and relative abundance of a. The
log2 tails grow with depth (one-sided wells cap near log2(2*depth)), and 20260907 is ~1.8x
deeper, so relative abundance is the depth-robust check.

    python cross_experiment_replicates.py      (env karl_seq_analysis, ~1 min)
"""
import itertools
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from config import (A, B, RA_OUT, CONSENSUS, LAYOUT, NO_GROWTH_B,  # noqa: E402
                    MIN_RESOLVABLE_BP, MIN_WELL_READS)
from io_utils import load_reference_db                               # noqa: E402
from strain_identity import flexible_distance, NOISE_FLOOR           # noqa: E402
from relative_abundance import fit_bradley_terry                     # noqa: E402

OUT = HERE / "outputs"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(0)
N_BOOT = 2000
LN2 = np.log(2)
pd.set_option("display.width", 200)


def header(t):
    print(f"\n{'=' * 90}\n{t}\n{'=' * 90}")


# --- x01 identity gate --------------------------------------------------------------------------
header("x01 are the shared strain names the same organisms?")
ref = {e: load_reference_db(CONSENSUS[e]) for e in (A, B)}
shared = sorted(set(ref[A]) & set(ref[B]))
rows = []
for s in shared:
    d, ov = flexible_distance(ref[A][s], ref[B][s])
    # best match of A's sequence among ALL of B's references, to tell "different organism"
    # from "16S cannot tell" (a twin elsewhere in the collection)
    best = min(((flexible_distance(ref[A][s], ref[B][t])[0], t) for t in ref[B]), key=lambda x: x[0])
    rows.append({"strain": s, "dist_same_name": d, "overlap_bp": ov,
                 "best_match_in_B": best[1], "best_dist": best[0]})
x01 = pd.DataFrame(rows)
x01["same_organism"] = x01.dist_same_name <= NOISE_FLOOR
x01.to_csv(OUT / "x01_identity.csv", index=False)
print(f"strains with a consensus in both: {len(x01)} (of {len(set(ref[A]))} in {A}, {len(set(ref[B]))} in {B})")
print(f"same-name sequences within the noise floor ({NOISE_FLOOR}): {x01.same_organism.sum()}/{len(x01)}")
print(f"median distance same name: {x01.dist_same_name.median():.4f}")
bad = x01[~x01.same_organism]
if len(bad):
    print("\nnames whose sequences disagree (excluded below):")
    print(bad.round(4).to_string(index=False))
IDENTICAL = set(x01.loc[x01.same_organism, "strain"])


# --- wells ------------------------------------------------------------------------------------
def load_wells(exp):
    w = pd.read_csv(RA_OUT[exp] / "r02_well_interaction_scores.csv")
    lay = pd.read_csv(LAYOUT[exp])
    lay["sample_id"] = [f"Plate{int(p):02d}_{wl}" for p, wl in zip(lay.dest_plate, lay.dest_well)]
    keep = ["sample_id", "well_type"] + [c for c in ("vol1_nL", "vol2_nL", "log2_inoculum_ratio")
                                         if c in lay.columns]
    w = w.drop(columns=[c for c in ("well_type",) if c in w.columns]).merge(lay[keep], on="sample_id")
    n0 = len(w)
    if exp == B:
        std = (w.vol1_nL == 100) & (w.vol2_nL == 100) & (w.log2_inoculum_ratio == 0)
        w = w[std & w.well_type.isin(["pair", "techrep", "ratio", "density"])]
    else:
        w = w[w.well_type == "pair"]
    w = w[~w.missing_reference.astype(bool)]
    w = w[w.ref_pair_bp_dist >= MIN_RESOLVABLE_BP]
    w["n_a"] = np.where(w.strain1 == w.strain_a, w.n_strain1, w.n_strain2)
    w["n_b"] = np.where(w.strain1 == w.strain_a, w.n_strain2, w.n_strain1)
    w = w[(w.n_a + w.n_b) >= MIN_WELL_READS].copy()
    w["log2_ab"] = np.log2((w.n_a + .5) / (w.n_b + .5))
    w["ra_a"] = w.n_a / (w.n_a + w.n_b)
    w["exp"] = exp
    print(f"{exp}: {n0} scored wells -> {len(w)} standard, resolvable, >= {MIN_WELL_READS} reads "
          f"({w.groupby(['strain_a', 'strain_b']).ngroups} pairs, median {int((w.n_a + w.n_b).median())} reads)")
    return w[["sample_id", "strain_a", "strain_b", "n_a", "n_b", "log2_ab", "ra_a", "ref_pair_bp_dist", "exp"]]


header("wells used")
no_growth = set(pd.read_csv(NO_GROWTH_B).strain)
W = {}
for e in (A, B):
    w = load_wells(e)
    before = len(w)
    w = w[~w.strain_a.isin(no_growth) & ~w.strain_b.isin(no_growth)]
    print(f"   dropping pairs with a {B} non-grower: {before} -> {len(w)} wells")
    W[e] = w


def pair_means(w):
    g = w.groupby(["strain_a", "strain_b"])
    return g.agg(n_wells=("sample_id", "size"), log2_ab=("log2_ab", "mean"), ra_a=("ra_a", "mean"),
                 n_a=("n_a", "sum"), n_b=("n_b", "sum")).reset_index()


P = {e: pair_means(W[e]) for e in (A, B)}


def metrics(x, y):
    """y predicted by x: no-fit R^2 against identity, linear-fit r^2, correlations, winner."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    sst = ((y - y.mean()) ** 2).sum()
    decided = (np.abs(x) > 0) & (np.abs(y) > 0)
    return {"n": len(x),
            "r2_identity": 1 - ((y - x) ** 2).sum() / sst,
            "r2_linear": pearsonr(x, y)[0] ** 2,
            "pearson_r": pearsonr(x, y)[0],
            "spearman_rho": spearmanr(x, y)[0],
            "slope_y_on_x": np.polyfit(x, y, 1)[0],
            "winner_agreement": np.mean(np.sign(x[decided]) == np.sign(y[decided]))}


def boot_ci(x, y, key):
    x, y = np.asarray(x), np.asarray(y)
    v = []
    for _ in range(N_BOOT):
        i = RNG.integers(0, len(x), len(x))
        try:
            v.append(metrics(x[i], y[i])[key])
        except Exception:
            pass
    return np.percentile(v, [2.5, 97.5])


def metric_rows(label, x, y, target, ci=True):
    m = metrics(x, y)
    row = {"comparison": label, "target": target, **m}
    if ci:
        for k in ("r2_identity", "pearson_r", "spearman_rho"):
            lo, hi = boot_ci(x, y, k)
            row[f"{k}_ci"] = f"[{lo:.2f}, {hi:.2f}]"
    return row


# --- x02 direct replicates -------------------------------------------------------------------------
header("x02 DIRECT: pairs tested in both experiments (pair means)")
D = P[A].merge(P[B], on=["strain_a", "strain_b"], suffixes=("_A", "_B"))
D_all = D.copy()
D = D[D.strain_a.isin(IDENTICAL) & D.strain_b.isin(IDENTICAL)]
print(f"pairs scoreable in both: {len(D_all)}; after identity gate: {len(D)}")
# relative abundance on the logit-free scale is bounded, so compare it as is; for "predict the
# other experiment" use RA centred at 0.5 so winner agreement is defined the same way
rows = []
for tgt, xa, xb in (("log2_ab", "log2_ab_A", "log2_ab_B"), ("ra_a - 0.5", "ra_a_A", "ra_a_B")):
    xA = D[xa] - (0.5 if tgt != "log2_ab" else 0)
    xB = D[xb] - (0.5 if tgt != "log2_ab" else 0)
    rows.append(metric_rows(f"{A} -> {B}", xA, xB, tgt))
    rows.append(metric_rows(f"{B} -> {A}", xB, xA, tgt))
x02 = pd.DataFrame(rows)
x02.to_csv(OUT / "x02_direct_pairs.csv", index=False)
D.to_csv(OUT / "x02_direct_pairs_values.csv", index=False)
print(x02.round(3).to_string(index=False))


# --- x03 yardstick: single wells within and across experiments ------------------------------------
header("x03 yardstick: one WELL predicting another well of the same pair")


def within_pairs(w, restrict=None):
    """All ordered-at-random well pairs of the same strain pair within one experiment."""
    xs, ys = [], []
    for key, g in w.groupby(["strain_a", "strain_b"]):
        if restrict is not None and key not in restrict:
            continue
        if len(g) < 2:
            continue
        for i, j in itertools.combinations(range(len(g)), 2):
            a, b = (i, j) if RNG.random() < .5 else (j, i)
            xs.append(g.iloc[a]); ys.append(g.iloc[b])
    return pd.DataFrame(xs).reset_index(drop=True), pd.DataFrame(ys).reset_index(drop=True)


def across_pairs(wa, wb, restrict):
    xs, ys = [], []
    ga, gb = dict(list(wa.groupby(["strain_a", "strain_b"]))), dict(list(wb.groupby(["strain_a", "strain_b"])))
    for key in restrict:
        for _, ra_ in ga[key].iterrows():
            for _, rb_ in gb[key].iterrows():
                xs.append(ra_); ys.append(rb_)
    return pd.DataFrame(xs).reset_index(drop=True), pd.DataFrame(ys).reset_index(drop=True)


shared_keys = set(zip(D.strain_a, D.strain_b))
rows = []
for label, (xw, yw) in {
    f"within {A} (all pairs)": within_pairs(W[A]),
    f"within {B} (all pairs)": within_pairs(W[B]),
    f"within {A} (shared pairs)": within_pairs(W[A], shared_keys),
    f"within {B} (shared pairs)": within_pairs(W[B], shared_keys),
    f"across, well vs well (shared pairs)": across_pairs(W[A], W[B], shared_keys),
}.items():
    if len(xw) < 5:
        print(f"{label}: only {len(xw)} well pairs, skipped")
        continue
    for tgt in ("log2_ab", "ra_a"):
        x = xw[tgt] - (0.5 if tgt == "ra_a" else 0)
        y = yw[tgt] - (0.5 if tgt == "ra_a" else 0)
        rows.append(metric_rows(label, x, y, tgt if tgt == "log2_ab" else "ra_a - 0.5", ci=False))
x03 = pd.DataFrame(rows)
x03.to_csv(OUT / "x03_well_level_yardstick.csv", index=False)
print(x03.round(3).to_string(index=False))


# --- x04 strength transfer ---------------------------------------------------------------------------
header("x04 STRENGTH TRANSFER: one experiment's Bradley-Terry strengths predict the other's pairs")


def fit_strengths(w, strains_ok):
    pooled = w.groupby(["strain_a", "strain_b"]).agg(n_reads_a=("n_a", "sum"), n_reads_b=("n_b", "sum")).reset_index()
    pooled = pooled[pooled.strain_a.isin(strains_ok) & pooled.strain_b.isin(strains_ok)]
    strains = sorted(set(pooled.strain_a) | set(pooled.strain_b))
    # keep the largest connected component, as hierarchy_analysis does
    from relative_abundance import union_find_components
    comp = set(union_find_components(strains, list(zip(pooled.strain_a, pooled.strain_b)))[0])
    pooled = pooled[pooled.strain_a.isin(comp) & pooled.strain_b.isin(comp)]
    s, se, pr2, _ = fit_bradley_terry(pooled, comp)
    return s, pr2, len(comp)


rows = []
S = {}
for fit_on, pred in ((A, B), (B, A)):
    s, pr2, n = fit_strengths(W[fit_on], IDENTICAL)
    S[fit_on] = s
    t = P[pred][P[pred].strain_a.isin(s.index) & P[pred].strain_b.isin(s.index)].copy()
    t["pred_log2"] = (t.strain_a.map(s) - t.strain_b.map(s)) / LN2        # BT is in natural log-odds
    t["pred_ra"] = 1 / (1 + np.exp(-(t.strain_a.map(s) - t.strain_b.map(s))))
    t["in_fit_exp"] = [k in set(zip(P[fit_on].strain_a, P[fit_on].strain_b)) for k in zip(t.strain_a, t.strain_b)]
    print(f"fit on {fit_on}: {n} strains, pseudo-R2 {pr2:.3f}; predicting {len(t)} {pred} pairs "
          f"({(~t.in_fit_exp).sum()} of them never tested in {fit_on})")
    for subset, m in (("all pairs", t.index == t.index), (f"pairs NOT tested in {fit_on}", ~t.in_fit_exp)):
        rows.append(metric_rows(f"BT({fit_on}) -> {pred}, {subset}", t.pred_log2[m], t.log2_ab[m], "log2_ab"))
        rows.append(metric_rows(f"BT({fit_on}) -> {pred}, {subset}", t.pred_ra[m] - .5, t.ra_a[m] - .5, "ra_a - 0.5"))
    # in-experiment reference: the same model scored on the data it was fitted to
    own = P[fit_on][P[fit_on].strain_a.isin(s.index) & P[fit_on].strain_b.isin(s.index)]
    rows.append(metric_rows(f"BT({fit_on}) -> {fit_on} (in-sample reference)",
                            (own.strain_a.map(s) - own.strain_b.map(s)) / LN2, own.log2_ab, "log2_ab", ci=False))
    t.to_csv(OUT / f"x04_predictions_bt_{fit_on}_to_{pred}.csv", index=False)
x04 = pd.DataFrame(rows)
x04.to_csv(OUT / "x04_strength_transfer.csv", index=False)
print(x04.round(3).to_string(index=False))

common = S[A].index.intersection(S[B].index)
sa, sb = S[A][common] - S[A][common].mean(), S[B][common] - S[B][common].mean()
x04s = {"n_strains": len(common), "pearson_r": pearsonr(sa, sb)[0], "spearman_rho": spearmanr(sa, sb)[0],
        "slope_B_on_A": np.polyfit(sa, sb, 1)[0]}
pd.DataFrame([x04s]).to_csv(OUT / "x04_strength_correlation.csv", index=False)
print(f"\nper-strain BT strength, {A} vs {B}: n={len(common)}, r={x04s['pearson_r']:.3f}, "
      f"rho={x04s['spearman_rho']:.3f}, slope={x04s['slope_B_on_A']:.2f}")

# --- figure ----------------------------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))
ax[0].scatter(D.log2_ab_A, D.log2_ab_B, s=14, alpha=.6)
lim = [min(D.log2_ab_A.min(), D.log2_ab_B.min()) - .5, max(D.log2_ab_A.max(), D.log2_ab_B.max()) + .5]
ax[0].plot(lim, lim, "k--", lw=.8)
ax[0].set_xlabel(f"log2(a/b), {A}"); ax[0].set_ylabel(f"log2(a/b), {B}")
m = x02[(x02.target == "log2_ab")].iloc[0]
ax[0].set_title(f"Pairs tested in both (n={len(D)})\nr={m.pearson_r:.2f}, rho={m.spearman_rho:.2f}, R2 vs identity={m.r2_identity:.2f}", fontsize=10)
ax[1].scatter(D.ra_a_A, D.ra_a_B, s=14, alpha=.6)
ax[1].plot([0, 1], [0, 1], "k--", lw=.8)
ax[1].set_xlabel(f"relative abundance of a, {A}"); ax[1].set_ylabel(f"relative abundance of a, {B}")
ax[1].set_title("Same pairs, relative abundance", fontsize=10)
ax[2].scatter(sa, sb, s=14, alpha=.6)
ax[2].set_xlabel(f"BT strength, {A} (centred)"); ax[2].set_ylabel(f"BT strength, {B} (centred)")
ax[2].set_title(f"Per-strain competitiveness (n={len(common)})\nr={x04s['pearson_r']:.2f}, rho={x04s['spearman_rho']:.2f}", fontsize=10)
fig.tight_layout(); fig.savefig(FIG / "x_cross_experiment.png", dpi=140); plt.close(fig)
print(f"\nwrote {OUT}")
