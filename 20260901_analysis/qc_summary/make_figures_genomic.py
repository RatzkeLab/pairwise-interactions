"""Figures for QC_03_GENOMIC_FEATURES.md. Appends its own provenance row."""
import sys
from datetime import datetime
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "shared_pipelines"))
from genomic_ml import COLOR_BLUE, COLOR_RED, COLOR_GRID, COLOR_TEXT_SECONDARY, COLOR_GOOD, COLOR_CRITICAL
OUT, FIG = HERE / "outputs", HERE / "outputs" / "figures"
STAMP = datetime.now().strftime("%Y-%m-%d %H:%M")

d = pd.read_csv(OUT / "gq02_vs_genome_16S_detail.csv")
t = d[d.self_rank.notna()]
COL = {"corroborated_db(185)": COLOR_GOOD, "20260630_pair_consensus": COLOR_BLUE,
       "20260721_pair_consensus": COLOR_RED}

fig, axes = plt.subplots(1, 3, figsize=(17, 5))
ax = axes[0]
for q, s in t.groupby("query_set"):
    ax.hist(s.self_identity, bins=np.linspace(0.70, 1.0, 31), histtype="step", lw=2,
            color=COL[q], label=f"{q}  (n={len(s)}, median {s.self_identity.median():.3f})")
ax.axvline(0.99, color=COLOR_CRITICAL, ls="--", lw=1)
ax.set_xlabel("identity to its OWN genome's 16S"); ax.set_ylabel("labels")
ax.legend(fontsize=7, frameon=False)
ax.set_title("A. Identity to the mapped genome\n(confounded by NGS-vs-ONT platform)", fontsize=10)

ax = axes[1]
for q, s in t.groupby("query_set"):
    r = np.sort(s.self_rank.values)
    ax.plot(r, 100 * np.arange(1, len(r) + 1) / len(r), lw=2, color=COL[q], label=q)
ax.axvline(10, color=COLOR_TEXT_SECONDARY, ls="--", lw=1)
ax.text(11, 20, "top-10", fontsize=8, color=COLOR_TEXT_SECONDARY)
ax.set_xscale("log"); ax.set_xlabel("rank of the correct genome (of 233)")
ax.set_ylabel("% of labels at or better than this rank")
ax.legend(fontsize=7.5, frameon=False)
ax.set_title("B. Rank retrieval (threshold-free)\nthe honest measure", fontsize=10)

ax = axes[2]
s = pd.read_csv(OUT / "gq03_vs_genome_16S_summary.csv")
x = np.arange(len(s)); w = 0.27
for i, (c, lab) in enumerate([("top1", "top-1"), ("top5", "top-5"), ("top10", "top-10")]):
    ax.bar(x + i * w, s[c], w, label=lab,
           color=[COLOR_TEXT_SECONDARY, COLOR_GRID, COLOR_BLUE][i])
ax.set_xticks(x + w); ax.set_xticklabels([q.replace("_pair_consensus", "\npair")
                                          .replace("corroborated_db(185)", "corroborated_db\n(reference vs reference)")
                                          for q in s.query_set], fontsize=7.5)
ax.set_ylabel("% of shared labels"); ax.legend(frameon=False, fontsize=8)
ax.set_title("C. Even reference-vs-reference tops out\nat 27% top-1 / 69% top-10", fontsize=10)
fig.suptitle("Does the genome attached to each well match that well's 16S?", fontweight="bold")
fig.tight_layout(); fig.savefig(FIG / "gq_f1_genome_identity_and_rank.png", dpi=160); plt.close(fig)

p = pd.read_csv(OUT / "gq00_provenance.csv")
p = pd.concat([p, pd.DataFrame([{"section": "figures", "generated_by": "make_figures_genomic.py",
                                 "run_at": STAMP, "outputs": "gq_f1_genome_identity_and_rank.png",
                                 "note": ""}])], ignore_index=True)
p.to_csv(OUT / "gq00_provenance.csv", index=False)
print("wrote gq_f1_genome_identity_and_rank.png at", STAMP)
