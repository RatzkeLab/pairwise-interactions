"""Four figures, one per question the QC asks. Palette matched to the rest of the project."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import qc_config as C
from qc_compare import SAME_ORGANISM_IDENTITY as THR

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

PAIR_LABEL = {
    ("corroborated_db", "pair_consensus"): "experiment (pair wells)\nvs corroborated_db",
    ("corroborated_db", "mono_consensus"): "experiment (mono wells)\nvs corroborated_db",
    ("genome_16S", "pair_consensus"): "experiment (pair wells)\nvs genome 16S",
    ("genome_16S", "mono_consensus"): "experiment (mono wells)\nvs genome 16S",
    ("mono_consensus", "pair_consensus"): "mono wells\nvs pair wells (internal)",
    ("corroborated_db", "genome_16S"): "corroborated_db\nvs genome 16S (both external)",
}


def _agree(exp):
    return pd.read_csv(exp.out / f"s01_source_agreement_{exp.name}.csv")


def f01_source_agreement():
    """The headline: which experiment agrees with the outside world, and which does not."""
    fig, ax = plt.subplots(figsize=(12, 5.6))
    exps = list(C.EXPERIMENTS.values())
    keys = list(PAIR_LABEL)
    w = 0.38
    for i, (exp, color) in enumerate(zip(exps, [COLOR_BLUE, COLOR_RED])):
        d = _agree(exp)
        vals, ns = [], []
        for k in keys:
            g = d[(d.source_a == k[0]) & (d.source_b == k[1])]
            vals.append(100 * g["same_organism"].mean() if len(g) else np.nan)
            ns.append(len(g))
        pos = np.arange(len(keys)) + i * w
        ax.bar(pos, vals, w, color=color, label=exp.name)
        for p, v, n in zip(pos, vals, ns):
            if np.isfinite(v):
                ax.text(p, v + 1.5, f"n={n}", ha="center", fontsize=7, color=COLOR_TEXT_SECONDARY)
    ax.set_xticks(np.arange(len(keys)) + w / 2)
    ax.set_xticklabels([PAIR_LABEL[k] for k in keys], fontsize=8.5)
    ax.set_ylabel(f"% of strains matching at ≥{THR:.2f} identity")
    ax.set_ylim(0, 108)
    ax.legend(frameon=False, title="experiment")
    ax.set_title("Does each 16S source agree about which strain is which?")
    fig.tight_layout()
    fig.savefig(C.FIG / "f01_source_agreement.png", dpi=160)
    plt.close(fig)


def f02_identity_distributions():
    """Not just how often sources agree, but how badly they disagree when they do."""
    keys = [("corroborated_db", "pair_consensus"), ("genome_16S", "pair_consensus"),
            ("mono_consensus", "pair_consensus"), ("corroborated_db", "genome_16S")]
    fig, axes = plt.subplots(1, len(keys), figsize=(4.0 * len(keys), 4.4), sharey=True)
    for ax, k in zip(axes, keys):
        for exp, color in zip(C.EXPERIMENTS.values(), [COLOR_BLUE, COLOR_RED]):
            d = _agree(exp)
            g = d[(d.source_a == k[0]) & (d.source_b == k[1])]
            if not len(g):
                continue
            x = np.random.default_rng(0).normal(0, 0.055, len(g)) + (0 if exp.name == "20260721" else 1)
            ax.scatter(x, g["identity"], s=16, alpha=0.55, color=color, edgecolor="none")
        ax.axhline(THR, color=COLOR_GOOD, ls="--", lw=1.1)
        ax.set_xticks([0, 1]); ax.set_xticklabels(list(C.EXPERIMENTS), fontsize=9)
        ax.set_title(PAIR_LABEL[k], fontsize=9)
        ax.set_xlim(-0.5, 1.5)
    axes[0].set_ylabel("16S identity (best over strands & copies)")
    axes[0].text(-0.45, THR + 0.004, "same organism", color=COLOR_GOOD, fontsize=7.5)
    fig.suptitle("How far apart are the sources when they disagree?", fontweight="bold")
    fig.tight_layout()
    fig.savefig(C.FIG / "f02_identity_distributions.png", dpi=160)
    plt.close(fig)


def f03_attribution():
    """Self vs best match: a label that is merely noisy sits on the diagonal; a mixed-up one
    sits far below it, because something ELSE matches nearly perfectly."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharex=True, sharey=True)
    for ax, exp in zip(axes, C.EXPERIMENTS.values()):
        f = exp.out / f"s02_attribution_pair_vs_corroborated_db_{exp.name}.csv"
        d = pd.read_csv(f).dropna(subset=["self_identity"])
        ok = d["self_identity"] >= THR
        ax.scatter(d.loc[ok, "best_identity"], d.loc[ok, "self_identity"], s=26,
                   color=COLOR_GOOD, alpha=0.75, edgecolor="white", lw=0.5, label="label supported")
        ax.scatter(d.loc[~ok, "best_identity"], d.loc[~ok, "self_identity"], s=26,
                   color=COLOR_CRITICAL, alpha=0.8, edgecolor="white", lw=0.5,
                   label="label contradicted")
        ax.plot([0.6, 1.01], [0.6, 1.01], color=COLOR_TEXT_SECONDARY, lw=1, ls="--")
        ax.axhline(THR, color=COLOR_GOOD, lw=0.9, ls=":")
        ax.set_title(f"{exp.name} — {int((~ok).sum())}/{len(d)} labels contradicted")
        ax.set_xlabel("identity to its BEST match in corroborated_db")
        ax.legend(frameon=False, fontsize=8, loc="lower left")
    axes[0].set_ylabel("identity to the strain the LABEL claims")
    fig.suptitle("Is the well what its label says, or is it confidently something else?",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(C.FIG / "f03_attribution.png", dpi=160)
    plt.close(fig)


def f04_resolution_limit():
    """Why re-identifying 20260721 from 16S cannot work, in two numbers."""
    exp = C.EXPERIMENTS["20260721"]
    att = pd.read_csv(exp.out / f"s02_attribution_pair_vs_genome_16S_{exp.name}.csv")
    conf = att[att["best_identity"] >= THR]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

    counts = conf["best_match"].value_counts()
    axes[0].hist(counts.values, bins=np.arange(0.5, counts.max() + 1.5), color=COLOR_RED)
    axes[0].set_xlabel("experiment labels assigned to the SAME genome")
    axes[0].set_ylabel("number of genomes")
    axes[0].set_title(f"{len(conf)} labels collapse onto {conf['best_match'].nunique()} genomes")

    axes[1].hist(conf["margin_over_runner_up"] * 1420, bins=25, color=COLOR_RED)
    axes[1].axvline(0, color=COLOR_CRITICAL, lw=1.2)
    axes[1].set_xlabel("how much the best call beats the runner-up (bp per 1420)")
    axes[1].set_ylabel("number of labels")
    axes[1].set_title(f"median margin {conf['margin_over_runner_up'].median()*1420:.1f} bp "
                      f"— 16S cannot separate them")
    fig.suptitle("Why 20260721 cannot be rescued by re-identifying wells from 16S",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(C.FIG / "f04_resolution_limit.png", dpi=160)
    plt.close(fig)


def make_all():
    C.FIG.mkdir(parents=True, exist_ok=True)
    f01_source_agreement(); f02_identity_distributions(); f03_attribution(); f04_resolution_limit()
    f05_source_plate_map()
    return sorted(p.name for p in C.FIG.glob("*.png"))


def f05_source_plate_map():
    """Where the distinguishable strains actually live, and where 20260721 found them.

    Left: the 384 source plate coloured by 16S group. One group holds 181 of 294 strains, which
    is the single fact behind every resolution failure in this QC -- most of this collection is
    one 16S-indistinguishable blob, and no amount of care with the amplicon changes that.

    Right: for the strains that ARE distinguishable, an arrow from where the label says the
    strain should be to where its sequence was actually found. If a quadrant swap or a
    backwards plate were responsible, these arrows would be parallel.
    """
    import qc_layout as L
    group_of, sizes = L.sixteen_s_groups()
    res = L.resolvable_strains()
    resolvable = set(res.loc[res["resolvable"], "strain_label"])

    fig, axes = plt.subplots(1, 2, figsize=(17, 6.4))

    ax = axes[0]
    big = sizes.idxmax()
    for s, g in group_of.items():
        rc = L.well_to_rc(s)
        if rc is None:
            continue
        if g == big:
            col, sz = "#c9c8c4", 26            # the 181-strain indistinguishable blob
        elif s in resolvable:
            col, sz = COLOR_GOOD, 54
        else:
            col, sz = COLOR_BLUE, 34
        ax.scatter(rc[1], rc[0], s=sz, color=col, edgecolor="white", lw=0.5)
    ax.set_title(f"Source collection: 16S groups\n{len(resolvable)} of {len(group_of)} strains are "
                 f"individually distinguishable (green); {sizes.max()} share one group (grey)")
    for a in (ax,):
        a.set_xlim(-1, 24); a.set_ylim(16, -1)
        a.set_xticks(range(0, 24, 2)); a.set_xticklabels(range(1, 25, 2), fontsize=7)
        a.set_yticks(range(16)); a.set_yticklabels([chr(65 + i) for i in range(16)], fontsize=7)

    ax = axes[1]
    exp = C.EXPERIMENTS["20260721"]
    calls = L.positional_calls(exp, res)
    att = pd.read_csv(exp.out / f"s02_attribution_pair_vs_genome_16S_{exp.name}.csv")
    loose = att[att["best_identity"] >= 0.99].copy()
    loose["from_rc"] = loose["strain_label"].map(L.well_to_rc)
    loose["to_rc"] = loose["best_match"].map(L.well_to_rc)
    loose = loose.dropna(subset=["from_rc", "to_rc"])
    for r in loose.itertuples():
        (r0, c0), (r1, c1) = r.from_rc, r.to_rc
        ax.annotate("", xy=(c1, r1), xytext=(c0, r0),
                    arrowprops=dict(arrowstyle="->", color="#c9c8c4", lw=0.7, alpha=0.85))
    for r in calls.itertuples():
        (r0, c0), (r1, c1) = r.from_rc, r.to_rc
        ax.annotate("", xy=(c1, r1), xytext=(c0, r0),
                    arrowprops=dict(arrowstyle="->", color=COLOR_CRITICAL, lw=2.0))
        ax.scatter([c0], [r0], s=46, color=COLOR_GOOD, zorder=5, edgecolor="white", lw=0.6)
        ax.text(c1 + 0.4, r1, f"{r.strain_label}→{r.best_match}", fontsize=7.5, color=COLOR_CRITICAL,
                va="center")
    ax.set_xlim(-1, 24); ax.set_ylim(16, -1)
    ax.set_xticks(range(0, 24, 2)); ax.set_xticklabels(range(1, 25, 2), fontsize=7)
    ax.set_yticks(range(16)); ax.set_yticklabels([chr(65 + i) for i in range(16)], fontsize=7)
    ax.set_title(f"20260721: label position → where the strain was actually found\n"
                 f"grey = all {len(loose)} confident calls (16S groups); "
                 f"red = the {len(calls)} calls on distinguishable strains")
    fig.suptitle("If the plate were rotated or its quadrants swapped, these arrows would be parallel",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(C.FIG / "f05_source_plate_map.png", dpi=160)
    plt.close(fig)
