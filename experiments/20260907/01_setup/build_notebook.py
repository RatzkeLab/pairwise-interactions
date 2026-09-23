"""Regenerates experiment_setup.ipynb from the cell definitions below.

Kept as a script so the notebook can be rebuilt after changing a parameter
without hand-editing JSON.
"""
import json

CELLS = []
def md(s): CELLS.append(("markdown", s))
def code(s): CELLS.append(("code", s))

md(r'''# 20260907 pairwise interaction experiment — setup

Follows `../../20260721/01_setup/experiment_setup.ipynb` for the layout mechanics
(barcode grid, Echo transfer files, minibar table). What changes is **which wells
get which content**, driven by three findings from the reference-database work:

1. **The master plate cannot be assumed clean.** `merge_consensus_sequences/cross_contamination_check`
   found a few identities spreading across the collection, and wells splitting into
   more than one read cluster rose from 41% to 67% between Dec 2024 and Aug 2025.
   So *every* strain gets a monoculture well, and most get two.
2. **Pairs must be resolvable.** Only ~90% of strain pairs are >15 bp apart in 16S;
   the rest cannot be deconvolved from a mixed readout at all, so a well spent on
   them yields nothing. `separability.py` filters these out.
3. **8 strains have never yielded an ONT consensus in six attempts.** They get one
   mono and two pairs rather than a full allocation.

Plates **1-8 are self-sufficient**: sequencing only those still gives every
monoculture and the complete technical-replicate set. Plates 9+ add pair coverage
and can be read on the plate reader and frozen for later sequencing.''')

code(r'''import os
import string
import itertools
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

import separability as sep
import design

randseed_date = 20260907
rng = np.random.default_rng(randseed_date)

# --- capacity -------------------------------------------------------------
N_ROWS, N_COLS = 16, 24
N_PRIMARY_PLATES   = 8    # PCR'd + sequenced. Capped here because 10 was already a slog.
N_EXTENSION_PLATES = 8    # plate reader + freeze; sequence later. Can be raised toward 16.
WELLS_PER_PLATE    = (N_ROWS - 2) * (N_COLS - 2)   # 308 non-edge wells

# --- content targets ------------------------------------------------------
MONO_PER_STRAIN    = 2    # two independent monocultures -> a corroborated reference
MONO_LIMITED       = 1    # for strains that have never sequenced
N_TECHREP_PAIRS    = 50   # pairs replicated to spot-check reproducibility
TECHREP_REPLICATES = 5    # each on a different plate
MIN_SEPARATION_BP  = 15   # a pair must be at least this far apart to be worth a well
SEPARABILITY_MODE  = "conservative"   # min across all sources; see separability.py

# --- inoculum arms --------------------------------------------------------
# A monoculture has always been 2x100 nL of the same strain while each member of
# a pair gets 100 nL, so the mono was never the matched control for the pair.
# Giving one of a strain's two monos 100 nL fixes that AND buys a within-strain
# 2x density contrast across 376 strains for free.
MONO_VOLUMES = (200, 100)   # the first keeps the 20260721 volume for comparability

# Ratio titration: total inoculum held at 200 nL, only the split varies, so the
# 100:100 level IS the standard condition. Multiples of 25 nL (the Echo minimum
# for 384PP); 25:175 is a 7x ratio each way, ~49x span end to end.
RATIO_LEVELS   = [(25, 175), (50, 150), (100, 100), (150, 50), (175, 25)]
N_RATIO_SERIES = 50

# Density titration: ratio held at 1:1, total inoculum varied 4x.
DENSITY_LEVELS   = [(50, 50), (100, 100), (200, 200)]
N_DENSITY_SERIES = 30

N_PLATES = N_PRIMARY_PLATES + N_EXTENSION_PLATES
PRIMARY_CAPACITY = N_PRIMARY_PLATES * WELLS_PER_PLATE
TOTAL_CAPACITY   = N_PLATES * WELLS_PER_PLATE
print(f"{N_PLATES} plates x {WELLS_PER_PLATE} wells = {TOTAL_CAPACITY} "
      f"({PRIMARY_CAPACITY} of them sequenced)")

wells384 = np.array([r + str(c) for r in string.ascii_uppercase[:N_ROWS]
                     for c in range(1, N_COLS + 1)])
wells384_2d = wells384.reshape(N_ROWS, N_COLS)''')

md(r'''## Step 1 - strains, and what we know about each

The source plate is the 384-well master collection. Strains are split by how much
we can trust a 16S reference for them, which sets how many wells each deserves.''')

code(r'''PLATE_CSV = "/home/rl/scripts/karl/merge_consensus_sequences/plate_format/Full_384_strain_collection_no_seq_info.csv"
XREF = "/home/rl/scripts/karl/merge_consensus_sequences/external_reference_crosscheck"

df_plate = pd.read_csv(PLATE_CSV)
strains = list(df_plate["Well.position"])
strain_number = dict(zip(df_plate["Well.position"], df_plate["Strain.Number"]))
print(f"{len(strains)} wells on the master plate")

status = pd.read_csv(os.path.join(XREF, "results", "well_identity_status.csv")).set_index("well")

# Which wells has ONT ever produced a consensus for (our five runs, or the
# collaborator's pool1)? Those that never have are the hard-to-lyse/sequence
# ones -- six attempts is enough evidence to stop spending wells on them.
src_of = defaultdict(set)
for line in open(os.path.join(XREF, "work", "02_ids.txt")):
    w, s, _ = line.strip().split("|")
    src_of[w].add(s)

limited = [w for w in strains if not (src_of[w] & {"ont_v3", "pool1"})]
regular = [w for w in strains if w not in limited]

df_strains = pd.DataFrame({
    "strain": strains,
    "strain_number": [strain_number[w] for w in strains],
    "n_sources": [len(src_of[w]) for w in strains],
    "sources": [";".join(sorted(src_of[w])) for w in strains],
    "identity_confirmed": [str(status.loc[w, "identity_confirmed"]) == "True"
                           if w in status.index else False for w in strains],
    "ont_sequenced_before": [w in regular for w in strains],
})
df_strains["class"] = np.where(~df_strains.ont_sequenced_before, "never_sequenced",
                       np.where(df_strains.identity_confirmed, "confirmed", "uncertain"))
print(df_strains["class"].value_counts().to_string())
print(f"\nnever sequenced despite six ONT attempts: {sorted(limited)}")
display(df_strains.head())''')

md(r'''## Step 2 - which pairs are worth a well

A pair only earns a well if we expect to tell its members apart in the readout.
`separability.py` combines the four independent references conservatively: the
distance is the **minimum** over sources covering both wells, so if any source
says they are close, the pair is skipped. That direction is right because the
errors are asymmetric - an unresolvable pair yields nothing, while skipping a
resolvable one costs one of ~66,000 alternatives.''')

code(r'''dist_by_source, sep_wells = sep.build(strains)
S, n_src = sep.combine(dist_by_source, sep_wells, SEPARABILITY_MODE)
idx = {w: i for i, w in enumerate(sep_wells)}
off = ~np.eye(len(strains), dtype=bool)

print(f"pairs >= {MIN_SEPARATION_BP} bp apart (eligible): {int(((S >= MIN_SEPARATION_BP) & off).sum() // 2)}")
print(f"pairs closer than that (skipped):       {int(((S >= 0) & (S < MIN_SEPARATION_BP) & off).sum() // 2)}")
print(f"pairs no source can compare (skipped):  {int(((S < 0) & off).sum() // 2)}")

def eligible(a, b):
    return S[idx[a], idx[b]] >= MIN_SEPARATION_BP

def unknown(a, b):
    """No reference can compare these two -- distinct from 'known to be close'."""
    return S[idx[a], idx[b]] < 0

partners = ((S >= MIN_SEPARATION_BP) & off).sum(axis=1)
print(f"\neligible partners per strain: median {np.median(partners):.0f}, "
      f"min {partners.min()}, strains with none: {(partners == 0).sum()}")

fig, ax = plt.subplots(1, 2, figsize=(12, 4))
ax[0].hist(partners, bins=40, color="steelblue")
ax[0].set_xlabel("eligible partners"); ax[0].set_ylabel("strains")
ax[0].set_title("Partners each strain can be told apart from")
d = S[off & (S >= 0)]
ax[1].hist(np.clip(d, 0, 400), bins=80, color="firebrick")
ax[1].axvline(MIN_SEPARATION_BP, ls="--", c="k")
ax[1].set_xlabel("16S separation (bp, clipped at 400)"); ax[1].set_ylabel("pairs")
ax[1].set_title("Pairwise separation across the collection")
plt.tight_layout(); plt.show()''')

md(r'''## Step 2b - pairs for the inoculum arms

A previous analysis found that monoculture yield predicts pairwise outcome, and
that preculture density predicts the winner - but "denser inoculum gives a head
start" and "a fitter strain grows dense everywhere AND wins" both fit that, and
the standard design cannot separate them.

The discriminating quantity is the **slope of final log-ratio against initial
log-ratio**. Slope near 0 means the outcome is set by fitness and the starting
ratio is irrelevant; slope near 1 means initial conditions simply carry through
and the assay is largely reporting the inoculum; in between means both. A
titration over a ~49x span of starting ratios measures that slope directly.

Pairs come from 20260630, whose well labels join correctly to the collection
(20260721's do not - confirmed source-plate mix-up), and are **stratified by
prior outcome**: strong-winner pairs are where a head start could plausibly flip
the result, near-neutral pairs are where a ratio effect is easiest to detect.
Both strata are needed, since a null result in only one would be uninformative.''')

code(r'''PRIOR = ("/home/rl/scripts/karl/pairwise_interaction_experiments/20260630/"
         "05_engineer_relative_abundances/relative_abundance/outputs/"
         "r03_pair_replicate_stats.csv")
if not os.path.exists(PRIOR):
    PRIOR = PRIOR.replace("/outputs/", "/outputs_backup_pre_margin_fix/")

prior = pd.read_csv(PRIOR)
prior = prior[~prior.high_uncertainty_pair]
prior = prior[prior.strain_a.isin(strains) & prior.strain_b.isin(strains)]
# and still resolvable under the current, much better references
prior = prior[[eligible(a, b) for a, b in zip(prior.strain_a, prior.strain_b)]]
prior["abs_log2"] = prior.mean_log2_ratio_a_over_b.abs()
print(f"{len(prior)} prior pairs usable: both strains in the collection, "
      f">= {MIN_SEPARATION_BP} bp apart, not high-uncertainty")

strong  = prior[prior.abs_log2 > 3]
neutral = prior[prior.abs_log2 < 1]
print(f"  strong winner (|log2 ratio| > 3): {len(strong)}")
print(f"  near neutral  (|log2 ratio| < 1): {len(neutral)}")

def stratified(n):
    """Half strong-winner, half near-neutral, sampled without replacement."""
    k = n // 2
    a = strong.sample(min(k, len(strong)), random_state=int(rng.integers(1e6)))
    b = neutral.sample(min(n - len(a), len(neutral)), random_state=int(rng.integers(1e6)))
    return [(r.strain_a, r.strain_b) for r in pd.concat([a, b]).itertuples()]

ratio_series   = stratified(N_RATIO_SERIES)
density_series = stratified(N_DENSITY_SERIES)
print(f"\nratio titration:   {len(ratio_series)} pairs x {len(RATIO_LEVELS)} levels "
      f"= {len(ratio_series) * len(RATIO_LEVELS)} wells")
print(f"density titration: {len(density_series)} pairs x {len(DENSITY_LEVELS)} levels "
      f"= {len(density_series) * len(DENSITY_LEVELS)} wells")
print(f"strains involved in the inoculum arms: "
      f"{len({s for p in ratio_series + density_series for s in p})}")''')

md(r'''## Step 3 - allocate content

Budget, in priority order:

1. **Monocultures** - 2 per sequenceable strain, 1 for the never-sequenced ones.
   The two go on *different plates*, so one plate failing cannot cost both.
2. **Technical replicates** - 50 pairs at 5 copies each, every copy on a
   *different plate*. That is the only arrangement that tests plate-to-plate
   reproducibility rather than well-to-well.
3. **Distinct pairs** - successive random matchings over the strain set. One
   matching covers every strain exactly once, so *k* rounds give every strain
   exactly *k* partners with no balancing step and nobody left short.

All of 1-3 land in plates 1-8. Plates 9+ carry additional distinct pairs only.''')

code(r'''n_mono = MONO_PER_STRAIN * len(regular) + MONO_LIMITED * len(limited)
n_techrep_wells = N_TECHREP_PAIRS * TECHREP_REPLICATES
n_ratio_wells   = len(ratio_series) * len(RATIO_LEVELS)
n_density_wells = len(density_series) * len(DENSITY_LEVELS)
n_primary_pairs = (PRIMARY_CAPACITY - n_mono - n_techrep_wells
                   - n_ratio_wells - n_density_wells)
n_extension_pairs = N_EXTENSION_PLATES * WELLS_PER_PLATE

print(f"primary capacity       {PRIMARY_CAPACITY}")
print(f"  monocultures         {n_mono}   ({MONO_PER_STRAIN}x{len(regular)} + {MONO_LIMITED}x{len(limited)})")
print(f"  technical replicates {n_techrep_wells}   ({N_TECHREP_PAIRS} pairs x {TECHREP_REPLICATES})")
print(f"  ratio titration      {n_ratio_wells}   ({len(ratio_series)} pairs x {len(RATIO_LEVELS)})")
print(f"  density titration    {n_density_wells}   ({len(density_series)} pairs x {len(DENSITY_LEVELS)})")
print(f"  distinct pairs       {n_primary_pairs}")
print(f"extension capacity     {n_extension_pairs}   (distinct pairs; plate reader + freeze)")
assert n_primary_pairs > 0, "monos + replicates already exceed the primary plates"

# N_TECHREP_PAIRS extra, because the replicated pairs are pulled out of this
# list rather than also being placed as ordinary single wells
pairs, rounds = design.build_pairs(regular, limited,
                                   n_primary_pairs + n_extension_pairs + N_TECHREP_PAIRS,
                                   eligible, rng, unknown=unknown, limited_rounds=2)
print(f"\ngenerated {len(pairs)} distinct pairs over {max(rounds) + 1} matching rounds")

appear = Counter()
for a, b in pairs:
    appear[a] += 1; appear[b] += 1
v = np.array([appear[s] for s in strains])
print(f"pair appearances per strain: min {v.min()}, median {np.median(v):.0f}, max {v.max()}")
print(f"  never_sequenced strains: {[appear[s] for s in limited]}")
bad = [(a, b) for a, b in pairs if not eligible(a, b)]
print(f"pairs used that are NOT >= {MIN_SEPARATION_BP} bp apart: {len(bad)} "
      f"(only allowed for never-sequenced strains) -> {bad[:6]}")''')

md(r'''## Step 4 - barcode layout

Same scheme as 20260721: 96 forward and 96 reverse barcodes in a checkerboard on
the source plate, every forward x reverse combination valid except equal indices,
giving 9120 unique combinations. Each destination well gets its own combination,
so every well is individually demultiplexable - including the extension plates,
which lets them be sequenced later without redesigning anything.''')

code(r'''fprimers_full = wells384_2d[::2, ::2].flatten()     # 96 forward source wells
rprimers_full = wells384_2d[1::2, 1::2].flatten()   # 96 reverse source wells
n_fwd, n_rev = len(fprimers_full), len(rprimers_full)

combos = np.array([(i, j) for i in range(n_fwd) for j in range(n_rev) if i != j])
print(f"valid barcode combinations: {len(combos)}; wells to fill: {TOTAL_CAPACITY}")
assert len(combos) >= TOTAL_CAPACITY, "not enough unique barcode combinations"
rng.shuffle(combos)

n_cols_nonedge = N_COLS - 2
records = []
for k in range(TOTAL_CAPACITY):
    fi, rj = combos[k]
    plate  = k // WELLS_PER_PLATE + 1
    within = k %  WELLS_PER_PLATE
    dest_row = string.ascii_uppercase[(within // n_cols_nonedge) + 1]   # skip row A
    dest_col = (within % n_cols_nonedge) + 2                            # skip col 1
    records.append({"dest_plate": plate, "dest_row": dest_row, "dest_col": dest_col,
                    "dest_well": f"{dest_row}{dest_col}",
                    "fwd_idx": int(fi), "rev_idx": int(rj),
                    "fwd_source_well": fprimers_full[fi],
                    "rev_source_well": rprimers_full[rj]})
df_layout = pd.DataFrame(records)
df_layout["plate_role"] = np.where(df_layout.dest_plate <= N_PRIMARY_PLATES,
                                   "primary_sequenced", "extension_frozen")

assert (df_layout.fwd_idx == df_layout.rev_idx).sum() == 0
assert df_layout[["fwd_idx", "rev_idx"]].drop_duplicates().shape[0] == len(df_layout)
assert set(df_layout.fwd_source_well).isdisjoint(set(df_layout.rev_source_well))
print(f"layout rows: {len(df_layout)}")
display(df_layout.head())''')

md(r'''## Step 5 - place content into wells

Monocultures and technical replicates are scattered across plates 1-8 with the
per-group spreading described above; distinct pairs fill whatever is left, the
primary plates first so that plates 1-8 stand alone.''')

code(r'''primary_rows   = df_layout.index[df_layout.dest_plate <= N_PRIMARY_PLATES].to_numpy()
extension_rows = df_layout.index[df_layout.dest_plate >  N_PRIMARY_PLATES].to_numpy()
plate_of = df_layout.dest_plate.to_numpy()
content = [None] * len(df_layout)

free_by_plate = defaultdict(list)
for i in primary_rows:
    free_by_plate[plate_of[i]].append(i)
for p in free_by_plate:
    rng.shuffle(free_by_plate[p])

def take(plate):
    """Claim a free well on the given primary plate, else on the emptiest one."""
    if free_by_plate[plate]:
        return free_by_plate[plate].pop()
    alt = max(free_by_plate, key=lambda p: len(free_by_plate[p]))
    return free_by_plate[alt].pop()

primary_plates = list(range(1, N_PRIMARY_PLATES + 1))

# --- monocultures: the replicates of a strain go on different plates ---------
# The two monos of a strain get different total inoculum (200 and 100 nL). The
# 200 nL one matches 20260721 so the series stays comparable; the 100 nL one is
# the properly matched control for a pair member, and the two together give a
# within-strain density contrast at no extra well cost.
mono_groups = [[s] * (MONO_LIMITED if s in limited else MONO_PER_STRAIN) for s in strains]
for placed in design.spread_across_plates(mono_groups, primary_plates, rng):
    for k, (s, p) in enumerate(placed):
        content[take(p)] = ("mono", s, s, MONO_VOLUMES[k % len(MONO_VOLUMES)], 0)

# --- technical replicates: each copy on a different plate --------------------
techrep_idx = set(rng.choice(len(pairs), size=N_TECHREP_PAIRS, replace=False).tolist())
techrep_pairs = [pairs[i] for i in sorted(techrep_idx)]
single_pairs  = [p for i, p in enumerate(pairs) if i not in techrep_idx]
placements = design.spread_across_plates(
    [[p] * TECHREP_REPLICATES for p in techrep_pairs], primary_plates, rng)
for (a, b), placed in zip(techrep_pairs, placements):
    for _, p in placed:
        content[take(p)] = ("techrep", a, b, 100, 100)

# --- inoculum titrations -----------------------------------------------------
# All levels of one series go on the SAME plate, so the dose-response within a
# series is not confounded by plate effects; different series are spread across
# plates so plate effects average out between series rather than aliasing onto
# the ratio axis.
series_plates = [primary_plates[i % len(primary_plates)]
                 for i in range(len(ratio_series) + len(density_series))]
rng.shuffle(series_plates)
si = 0
for (a, b) in ratio_series:
    p = series_plates[si]; si += 1
    for v1, v2 in RATIO_LEVELS:
        content[take(p)] = ("ratio", a, b, v1, v2)
for (a, b) in density_series:
    p = series_plates[si]; si += 1
    for v1, v2 in DENSITY_LEVELS:
        content[take(p)] = ("density", a, b, v1, v2)

# --- distinct pairs: primary plates first, then the extension plates ---------
remaining = [i for p in free_by_plate for i in free_by_plate[p]]
rng.shuffle(remaining)
slots = list(remaining) + list(rng.permutation(extension_rows))
for (a, b), i in zip(single_pairs, slots):
    content[i] = ("pair", a, b, 100, 100)

assert all(c is not None for c in content), f"{sum(c is None for c in content)} wells unassigned"
df_layout["well_type"] = [c[0] for c in content]
df_layout["strain1"]   = [c[1] for c in content]
df_layout["strain2"]   = [c[2] for c in content]
df_layout["vol1_nL"]   = [c[3] for c in content]
df_layout["vol2_nL"]   = [c[4] for c in content]
df_layout["total_inoculum_nL"] = df_layout.vol1_nL + df_layout.vol2_nL
df_layout["log2_inoculum_ratio"] = np.where(
    df_layout.well_type == "mono", np.nan,
    np.log2(df_layout.vol1_nL / df_layout.vol2_nL.replace(0, np.nan)))
df_layout["separation_bp"] = [int(S[idx[c[1]], idx[c[2]]]) if c[0] != "mono" else -1
                              for c in content]
print(df_layout.groupby(["plate_role", "well_type"]).size().to_string())''')

md(r'''## Step 6 - checks that the design actually delivers

These assert the promises made above rather than merely describing the output.''')

code(r'''prim = df_layout[df_layout.dest_plate <= N_PRIMARY_PLATES]

mono_counts = Counter(prim.loc[prim.well_type == "mono", "strain1"])
assert set(mono_counts) == set(strains), "a strain has no monoculture on the primary plates"
print(f"monocultures: every one of {len(strains)} strains present on plates 1-{N_PRIMARY_PLATES}")
print(f"  {sum(1 for s in strains if mono_counts[s] == 2)} strains with 2, "
      f"{sum(1 for s in strains if mono_counts[s] == 1)} with 1")

mono_plates = defaultdict(list)
for r in prim[prim.well_type == "mono"].itertuples():
    mono_plates[r.strain1].append(r.dest_plate)
split_ok = sum(1 for s, ps in mono_plates.items() if len(ps) == len(set(ps)))
print(f"  {split_ok}/{len(mono_plates)} strains have their monocultures on distinct plates")

tr = prim[prim.well_type == "techrep"].copy()
tr["pair"] = [frozenset((a, b)) for a, b in zip(tr.strain1, tr.strain2)]
g = tr.groupby("pair").agg(n=("dest_plate", "size"), plates=("dest_plate", "nunique"))
assert (g.n == TECHREP_REPLICATES).all(), "a technical replicate set is incomplete"
print(f"technical replicates: {len(g)} pairs x {TECHREP_REPLICATES}, each spanning "
      f"{g.plates.min()}-{g.plates.max()} distinct plates")

# monoculture inoculum contrast
mv = prim[prim.well_type == "mono"].groupby("strain1").total_inoculum_nL.apply(
    lambda x: tuple(sorted(x)))
n_contrast = (mv == tuple(sorted(MONO_VOLUMES))).sum()
print(f"  {n_contrast} strains have both mono volumes {sorted(MONO_VOLUMES)} nL "
      f"(a within-strain density contrast)")

# titration series intact, and each series on one plate
for name, levels in [("ratio", RATIO_LEVELS), ("density", DENSITY_LEVELS)]:
    t = prim[prim.well_type == name].copy()
    t["pair"] = [frozenset((a, b)) for a, b in zip(t.strain1, t.strain2)]
    gg = t.groupby("pair").agg(n=("dest_plate", "size"), plates=("dest_plate", "nunique"),
                               levels=("total_inoculum_nL", "size"))
    assert (gg.n == len(levels)).all(), f"{name} series incomplete"
    assert (gg.plates == 1).all(), f"{name} series split across plates"
    print(f"{name} titration: {len(gg)} series x {len(levels)} levels, "
          f"each series on a single plate, spanning "
          f"{t.dest_plate.nunique()} plates overall")

spent = df_layout[df_layout.well_type != "mono"]
viol = spent[(spent.separation_bp >= 0) & (spent.separation_bp < MIN_SEPARATION_BP)]
unk = spent[spent.separation_bp < 0]
print(f"pair wells below {MIN_SEPARATION_BP} bp separation: {len(viol)}")
print(f"pair wells of unknown separability: {len(unk)} "
      f"(strains: {sorted(set(unk.strain1) | set(unk.strain2))})")

INTERACTION_TYPES = ["pair", "techrep", "ratio", "density"]
ap = Counter()
for r in prim[prim.well_type.isin(INTERACTION_TYPES)].itertuples():
    ap[r.strain1] += 1; ap[r.strain2] += 1
v = np.array([ap[s] for s in strains])
print(f"\ninteractions per strain on the sequenced plates: "
      f"min {v.min()}, median {np.median(v):.0f}, max {v.max()}")
print(f"  strains below 3 interactions: {(v < 3).sum()}  "
      f"(expected: only the {len(limited)} never-sequenced ones)")

apall = Counter()
for r in df_layout[df_layout.well_type.isin(INTERACTION_TYPES)].itertuples():
    apall[r.strain1] += 1; apall[r.strain2] += 1
va = np.array([apall[s] for s in strains])
print(f"interactions per strain across all {N_PLATES} plates: "
      f"min {va.min()}, median {np.median(va):.0f}, max {va.max()}")

fig, ax = plt.subplots(1, 2, figsize=(12, 4))
ax[0].hist(v, bins=range(0, int(v.max()) + 2), color="steelblue", align="left")
ax[0].set_title(f"Interactions per strain, plates 1-{N_PRIMARY_PLATES} (sequenced)")
ax[0].set_xlabel("interaction wells"); ax[0].set_ylabel("strains")
ax[1].hist(va, bins=range(0, int(va.max()) + 2), color="seagreen", align="left")
ax[1].set_title(f"Interactions per strain, all {N_PLATES} plates")
ax[1].set_xlabel("interaction wells")
plt.tight_layout(); plt.show()''')

md(r'''## Step 7 - Echo transfer files

Two 100 nL transfers per destination well. Monocultures get 100 nL of the same
strain twice, matching 20260721 so growth stays comparable with the existing data
- note this means a monoculture receives double the inoculum of either member of
a pair, which matters when comparing mono against pair yields.''')

code(r'''PRIMER_VOL = 250   # nL per primer transfer

strain_rows = []
for r in df_layout.itertuples():
    for src, vol in ((r.strain1, r.vol1_nL), (r.strain2, r.vol2_nL)):
        if vol <= 0:          # a 100 nL monoculture is a single transfer
            continue
        strain_rows.append({"Source Well": src, "Destination Plate Name": r.dest_plate,
                            "Destination Well": r.dest_well, "Transfer Volume": int(vol)})
df_strain_echo = pd.DataFrame(strain_rows)
assert (df_strain_echo["Transfer Volume"] % 25 == 0).all(), "volume below Echo 25 nL granularity"
print(f"strain transfers: {len(df_strain_echo)}  "
      f"(volumes used: {sorted(df_strain_echo['Transfer Volume'].unique())} nL)")

vol_per_source = df_strain_echo.groupby("Source Well")["Transfer Volume"].sum() / 1000
print(f"volume drawn per source well (uL): min {vol_per_source.min():.1f}, "
      f"median {vol_per_source.median():.1f}, max {vol_per_source.max():.1f}")
display(vol_per_source.sort_values(ascending=False).head())

primer_rows = []
for r in df_layout.itertuples():
    for src in (r.fwd_source_well, r.rev_source_well):
        primer_rows.append({"Source Well": src, "Destination Plate Name": r.dest_plate,
                            "Destination Well": r.dest_well, "Transfer Volume": PRIMER_VOL})
df_primer_echo = pd.DataFrame(primer_rows)
print(f"\nprimer transfers: {len(df_primer_echo)}")''')

code(r'''df_layout.to_csv(f"strain_layout_{randseed_date}.csv", index=False)
df_strains.to_csv(f"strain_classes_{randseed_date}.csv", index=False)
df_layout.to_csv(f"primer_layout_{randseed_date}.csv", index=False)

os.makedirs("echo_strains", exist_ok=True)
os.makedirs("primer_shooting", exist_ok=True)
for p in range(1, N_PLATES + 1):
    df_strain_echo[df_strain_echo["Destination Plate Name"] == p].to_csv(
        f"echo_strains/echo_strains_{randseed_date}_plate_{p:02d}.csv", index=False)
    df_primer_echo[df_primer_echo["Destination Plate Name"] == p].to_csv(
        f"primer_shooting/echo_primers_{randseed_date}_plate_{p:02d}.csv", index=False)

df_strain_echo.to_csv(f"echo_strains_{randseed_date}.csv", index=False)
half = len(df_strain_echo) // 2
df_strain_echo.iloc[:half].to_csv(f"echo_strains_A_{randseed_date}.csv", index=False)
df_strain_echo.iloc[half:].to_csv(f"echo_strains_B_{randseed_date}.csv", index=False)
print(f"wrote per-plate Echo files for {N_PLATES} plates")''')

md(r'''## Step 8 - minibar demultiplexing table

Every well across all plates gets a row, so the extension plates can be sequenced
later without regenerating anything. A primary-plates-only table is written too,
for the first sequencing run.''')

code(r'''# 20260721 read this from ../../all_experiments/, which no longer exists;
# this is the same file, kept where it actually lives.
primer_specs_path = "/home/rl/scripts/karl/drug_experiments/260130_8commtest/01_pcr_setup/PrimerPlateSpecs.csv"

df_barcodes = pd.read_csv(primer_specs_path,
                          usecols=["Well Position", "Sequence Name", "Sequence"]).copy()
df_barcodes["Sequence"] = df_barcodes["Sequence"].str.replace(" ", "", regex=False)
df_barcodes["Storage"] = np.concatenate((rprimers_full, fprimers_full))[:len(df_barcodes)]

def lookup(col, name):
    return (df_layout[["dest_plate", "dest_well", col]]
            .merge(df_barcodes[["Storage", "Sequence"]], left_on=col,
                   right_on="Storage", how="left", validate="m:1")
            .rename(columns={"Sequence": name}).drop(columns=["Storage", col]))

df_minibar = lookup("fwd_source_well", "FwIndex").merge(
    lookup("rev_source_well", "RvIndex"), on=["dest_plate", "dest_well"], how="inner")
df_minibar["SampleID"] = [f"Plate{int(p):02d}_{w}" for p, w in
                          zip(df_minibar.dest_plate, df_minibar.dest_well)]
df_minibar["FwPrimer"] = "AGRGTTYGATYMTGGCTCAG"
df_minibar["RwPrimer"] = "CGGYTACCTTGTTACGACTT"
df_minibar = df_minibar[["SampleID", "FwIndex", "FwPrimer", "RvIndex", "RwPrimer"]]

assert df_minibar.FwIndex.notna().all() and df_minibar.RvIndex.notna().all()
assert df_minibar.SampleID.nunique() == len(df_minibar)
df_minibar.to_csv(f"minibar_primers_{randseed_date}.tsv", sep="\t", index=False)

primary_ids = {f"Plate{int(p):02d}_{w}" for p, w in zip(prim.dest_plate, prim.dest_well)}
df_minibar[df_minibar.SampleID.isin(primary_ids)].to_csv(
    f"minibar_primers_{randseed_date}_primary_only.tsv", sep="\t", index=False)
print(f"minibar rows: {len(df_minibar)} total, "
      f"{df_minibar.SampleID.isin(primary_ids).sum()} on the sequenced plates")
display(df_minibar.head())''')

md('## Step 9 - summary')

code(r'''summary = pd.DataFrame([
    {"quantity": "plates sequenced (PCR)",           "value": N_PRIMARY_PLATES},
    {"quantity": "plates plate-reader + frozen",     "value": N_EXTENSION_PLATES},
    {"quantity": "wells sequenced",                  "value": len(prim)},
    {"quantity": "wells total",                      "value": len(df_layout)},
    {"quantity": "strains",                          "value": len(strains)},
    {"quantity": "monoculture wells",                "value": int((prim.well_type == "mono").sum())},
    {"quantity": "strains with 2 monocultures",      "value": sum(1 for s in strains if mono_counts[s] == 2)},
    {"quantity": "technical replicate wells",        "value": int((prim.well_type == "techrep").sum())},
    {"quantity": "ratio-titration wells",            "value": int((prim.well_type == "ratio").sum())},
    {"quantity": "density-titration wells",          "value": int((prim.well_type == "density").sum())},
    {"quantity": "strains with a mono density contrast", "value": int(n_contrast)},
    {"quantity": "distinct interaction pairs",       "value": len(pairs)},
    {"quantity": "  of which replicated 5x",         "value": N_TECHREP_PAIRS},
    {"quantity": "median interactions/strain (seq)", "value": float(np.median(v))},
    {"quantity": "median interactions/strain (all)", "value": float(np.median(va))},
])
display(summary)
summary.to_csv(f"design_summary_{randseed_date}.csv", index=False)''')

def to_source(s):
    lines = s.split("\n")
    return [l + "\n" for l in lines[:-1]] + [lines[-1]]

cells = []
for kind, text in CELLS:
    c = {"cell_type": kind, "metadata": {}, "source": to_source(text)}
    if kind == "code":
        c["execution_count"] = None
        c["outputs"] = []
    cells.append(c)

nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                  "name": "python3"},
                   "language_info": {"name": "python", "version": "3.11"}},
      "nbformat": 4, "nbformat_minor": 5}
json.dump(nb, open("experiment_setup.ipynb", "w"), indent=1)
print(f"wrote experiment_setup.ipynb with {len(cells)} cells")
