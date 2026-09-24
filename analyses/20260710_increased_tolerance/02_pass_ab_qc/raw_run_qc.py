"""Raw-level QC: MinKNOW run a (fastq_pass_a) vs run b (fastq_pass_b, after the library top-up).

Same flow cell (FBF36772), two MinKNOW runs:
    a  eb6d246d  2026-07-10 00:37 -> 15:00   (14.4 h)
    b  ad553134  2026-07-10 15:14 -> 07-11 19:51  (28.6 h), after loading more library

Uses MinKNOW's own per-read sequencing_summary (every read, pass and fail) and pore-activity
logs, so nothing here depends on demultiplexing. Sections:
    q01  per-run totals: reads, pass rate, length, quality, fraction full-length
    q02  end_reason, overall and by length class
    q03  trends over clock time (1 h bins): throughput, full-length fraction, quality,
         active channels, pore occupancy
    q04  fragment structure: where the 16S primers sit in short reads (read start vs read
         end), joined to end_reason -- distinguishes molecules that were already broken in
         the tube from reads that stopped early in the pore

    python raw_run_qc.py        (env karl_seq_analysis, ~3 min)
"""
import gzip
from pathlib import Path

import edlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

MINKNOW = Path("/var/lib/minknow/data/kr_pairwise_01/no_sample_id")
RUNS = {
    "a": MINKNOW / "20260710_0032_MN37816_FBF36772_eb6d246d",
    "b": MINKNOW / "20260710_1510_MN37816_FBF36772_ad553134",
}
FASTQ_PASS = {r: Path(f"/home/rl/scripts/karl/data_links/raw/sequencing_fastqs/20260710_fastq_pass_{r}")
              for r in RUNS}
LEN_LO, LEN_HI, MIN_Q = 1400, 1800, 10       # NanoFilt settings used by the demux
LEN_BINS = [0, 500, 1000, 1400, 1800, 3000, np.inf]
LEN_LABELS = ["<500", "500-1000", "1000-1400", "1400-1800", "1800-3000", ">3000"]
COLORS = {"a": "#3b75af", "b": "#ef8636"}

F_PRIMER = "AGRGTTYGATYMTGGCTCAG"   # 27F
R_PRIMER = "CGGYTACCTTGTTACGACTT"   # 1492R
IUPAC = [("Y", "C"), ("Y", "T"), ("R", "A"), ("R", "G"), ("M", "A"), ("M", "C")]
_COMP = str.maketrans("ACGTYRMK", "TGCARYKM")


def header(t):
    print(f"\n{'=' * 90}\n{t}\n{'=' * 90}")


def load_summary(run):
    f = next(RUNS[run].glob("sequencing_summary_*.txt"))
    cols = ["read_id", "channel", "start_time", "duration", "passes_filtering",
            "sequence_length_template", "mean_qscore_template", "end_reason"]
    d = pd.read_csv(f, sep="\t", usecols=cols)
    d = d.rename(columns={"sequence_length_template": "length", "mean_qscore_template": "q"})
    d["run"] = run
    return d


def run_start(run):
    s = next(RUNS[run].glob("final_summary_*.txt")).read_text()
    started = [l.split("=", 1)[1] for l in s.splitlines() if l.startswith("started=")][0]
    return pd.Timestamp(started)


# --- load ------------------------------------------------------------------------------------
summ = {r: load_summary(r) for r in RUNS}
t0 = {r: run_start(r) for r in RUNS}
origin = t0["a"]
for r, d in summ.items():
    # clock time since run a started, so both runs share one axis
    d["hours"] = (d["start_time"] + (t0[r] - origin).total_seconds()) / 3600
    d["len_class"] = pd.cut(d["length"], LEN_BINS, labels=LEN_LABELS, right=False)
    d["full_length"] = d["length"].between(LEN_LO, LEN_HI)
all_reads = pd.concat(summ.values(), ignore_index=True)

# --- q01 totals ----------------------------------------------------------------------------
header("q01 per-run totals (all reads in sequencing_summary; pass = MinKNOW qscore filter)")
rows = []
for r, d in summ.items():
    p = d[d.passes_filtering]
    rows.append({
        "run": r,
        "duration_h": round(d.hours.max() - d.hours.min(), 1),
        "reads_total": len(d), "reads_pass": len(p),
        "pass_pct": 100 * len(p) / len(d),
        "median_len_pass": p.length.median(),
        "median_q_pass": p.q.median(),
        "full_length_pct_of_pass": 100 * p.full_length.mean(),
        "nanofilt_survivors": int((p.full_length & (p.q >= MIN_Q)).sum()),
        "reads_per_hour_pass": len(p) / (d.hours.max() - d.hours.min()),
        "active_channels": d.channel.nunique(),
        "fastq_pass_files": len(list(FASTQ_PASS[r].glob("*.fastq.gz"))),
    })
q01 = pd.DataFrame(rows).set_index("run")
q01.to_csv(OUT / "q01_run_totals.csv")
print(q01.round(2).T.to_string())

lc = pd.crosstab(all_reads.loc[all_reads.passes_filtering, "run"],
                 all_reads.loc[all_reads.passes_filtering, "len_class"], normalize="index") * 100
lc.to_csv(OUT / "q01_length_classes_pct.csv")
print("\nlength class, % of pass reads\n" + lc.round(1).to_string())

# --- q02 end reasons -----------------------------------------------------------------------
header("q02 end_reason (pass reads)")
p = all_reads[all_reads.passes_filtering]
er = pd.crosstab(p.run, p.end_reason, normalize="index") * 100
er.to_csv(OUT / "q02_end_reason_pct.csv")
print(er.round(2).to_string())
er_len = (pd.crosstab([p.run, p.len_class], p.end_reason, normalize="index") * 100)
er_len.to_csv(OUT / "q02_end_reason_by_length_pct.csv")
print("\nby length class\n" + er_len.round(1).to_string())

# --- q03 time trends -------------------------------------------------------------------------
header("q03 trends over clock time (1 h bins, pass reads)")
p = p.assign(hour=np.floor(p.hours).astype(int))
tr = p.groupby(["run", "hour"]).agg(
    reads=("read_id", "size"),
    full_length_pct=("full_length", lambda x: 100 * x.mean()),
    median_q=("q", "median"),
    median_len=("length", "median"),
    active_channels=("channel", "nunique"),
).reset_index()


def pore_occupancy(run):
    f = next(RUNS[run].glob("pore_activity_*.csv"))
    pa = pd.read_csv(f)
    pa.columns = ["state", "minute", "samples"]
    w = pa.pivot_table(index="minute", columns="state", values="samples", aggfunc="sum").fillna(0)
    w["hours"] = (w.index / 60) + (t0[run] - origin).total_seconds() / 3600
    # occupancy: of the pores that are usable, how much time is spent reading DNA
    usable = w["strand"] + w["adapter"] + w["pore"]
    w["occupancy_pct"] = 100 * (w["strand"] + w["adapter"]) / usable.replace(0, np.nan)
    w["usable_pore_time"] = usable
    w["run"] = run
    return w.reset_index()[["run", "minute", "hours", "occupancy_pct", "usable_pore_time"]]


occ = pd.concat([pore_occupancy(r) for r in RUNS], ignore_index=True)
occ["hour"] = np.floor(occ.hours).astype(int)
tr = tr.merge(occ.groupby(["run", "hour"]).occupancy_pct.mean().rename("occupancy_pct").reset_index(),
              on=["run", "hour"], how="left")
tr.to_csv(OUT / "q03_hourly_trends.csv", index=False)
with pd.option_context("display.max_rows", 100):
    print(tr.round(1).to_string(index=False))

fig, ax = plt.subplots(5, 1, figsize=(10, 12), sharex=True)
for r in RUNS:
    t = tr[tr.run == r]
    ax[0].plot(t.hour + .5, t.reads, "o-", ms=3, c=COLORS[r], label=f"run {r}")
    ax[1].plot(t.hour + .5, t.full_length_pct, "o-", ms=3, c=COLORS[r])
    ax[2].plot(t.hour + .5, t.median_q, "o-", ms=3, c=COLORS[r])
    ax[3].plot(t.hour + .5, t.active_channels, "o-", ms=3, c=COLORS[r])
    o = occ[occ.run == r]
    ax[4].plot(o.hours, o.occupancy_pct.rolling(15, min_periods=1).mean(), c=COLORS[r], lw=1)
ax[0].set_ylabel("pass reads / h"); ax[0].legend()
ax[1].set_ylabel(f"% {LEN_LO}-{LEN_HI} bp")
ax[2].set_ylabel("median Q")
ax[3].set_ylabel("channels producing reads")
ax[4].set_ylabel("pore occupancy %\n(strand+adapter)/usable")
ax[4].set_xlabel("hours since run a started")
for a in ax:
    a.axvline((t0["b"] - origin).total_seconds() / 3600, c="grey", ls=":", lw=1)
ax[0].set_title("Run a vs run b over clock time (dotted line = run b start, after library top-up)")
fig.tight_layout(); fig.savefig(FIG / "q03_time_trends.png", dpi=140); plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4))
bins = np.logspace(np.log10(100), np.log10(10000), 80)
for r in RUNS:
    x = summ[r].loc[summ[r].passes_filtering, "length"]
    ax.hist(x, bins=bins, histtype="step", density=True, color=COLORS[r], lw=1.5, label=f"run {r}")
ax.axvspan(LEN_LO, LEN_HI, color="grey", alpha=.15, label="length filter")
ax.set_xscale("log"); ax.set_xlabel("read length (bp)"); ax.set_ylabel("density"); ax.legend()
ax.set_title("Read length distribution, pass reads")
fig.tight_layout(); fig.savefig(FIG / "q01_length_distribution.png", dpi=140); plt.close(fig)

# --- q04 fragment structure ------------------------------------------------------------------
header("q04 fragment structure: where are the primers in short reads?")


def rc(s):
    return s.translate(_COMP)[::-1]


def has_primer(seg):
    for pr in (F_PRIMER, R_PRIMER):
        if edlib.align(pr, seg, "HW", "distance", 5, additionalEqualities=IUPAC)["editDistance"] > -1:
            return True
    return False


def sample_reads(run, n_target=30000):
    """Every 4th fastq_pass file, first reads of each, so the sample spans the run."""
    files = sorted(FASTQ_PASS[run].glob("*.fastq.gz"), key=lambda f: int(f.stem.split("_")[-1].split(".")[0]))
    files = files[::4]
    per_file = n_target // len(files)
    rows = []
    for f in files:
        with gzip.open(f, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    rid = line[1:].split()[0]
                elif i % 4 == 1:
                    s = line.strip()
                    # primer at the read START means the adapter-ligated end the pore grabbed
                    # was an amplicon end; at the read END means the broken end went in first
                    rows.append({"read_id": rid, "length": len(s),
                                 "primer_at_start": has_primer(s[:150]),
                                 "primer_at_end": has_primer(rc(s[-150:]))})
                    if len(rows) % per_file == 0:
                        break
    return pd.DataFrame(rows)


frag = []
for r in RUNS:
    s = sample_reads(r).merge(summ[r][["read_id", "end_reason"]], on="read_id", how="left")
    s["run"] = r
    frag.append(s)
frag = pd.concat(frag, ignore_index=True)
frag["len_class"] = pd.cut(frag.length, LEN_BINS, labels=LEN_LABELS, right=False)
frag["structure"] = np.select(
    [frag.primer_at_start & frag.primer_at_end, frag.primer_at_start, frag.primer_at_end],
    ["both ends", "start only", "end only"], "neither")
q04 = (pd.crosstab([frag.run, frag.len_class], frag.structure, normalize="index") * 100)
q04["n"] = frag.groupby(["run", "len_class"], observed=True).size()
q04.to_csv(OUT / "q04_fragment_structure_pct.csv")
print(q04.round(1).to_string())
print("\nIf molecules broke in the tube, each fragment carrying one amplicon end is equally likely")
print("to enter the pore from either end: 'start only' ~= 'end only'. If reads stop early in the")
print("pore, the primer sits at the start: 'start only' >> 'end only'.")

one = frag[frag.structure.isin(["start only", "end only"]) & (frag.length < LEN_LO)]
q04b = pd.crosstab([one.run, one.structure], one.end_reason, normalize="index") * 100
q04b.to_csv(OUT / "q04_one_end_reads_end_reason_pct.csv")
print("\nend_reason of short one-end reads\n" + q04b.round(1).to_string())
frag.to_csv(OUT / "q04_fragment_sample.csv.gz", index=False)
print(f"\nwrote {OUT}")
