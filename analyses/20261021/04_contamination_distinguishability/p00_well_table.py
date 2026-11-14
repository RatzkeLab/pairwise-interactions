"""Stage 11a: one row per source-plate well, with everything the plate maps draw.

Joins what stages 1-10 found (reference tier, genome verdict, identity groups) with
three things computed here:

  invader          for a `wrong` well, which organism the well now holds: a genome clone
                   cluster (C###), a single genome (its strain number), or `none on file`
                   when no assembly carries the well's 16S. Its HOME wells are the wells
                   where that organism's genome is on file and still matches the well.
  secondary        a second organism in the source well. Evidence, strongest first:
                     both_monos    the same organism (<= 2 bp) in both 20260907 monos,
                                   which sat on different plates -> it came from the source
                     run_20261026  a later run from the SAME source plate: the organism
                                   tracks this strain across communities (association test,
                                   `community_assembly_experiments/20261026/02_qc`)
                     close_haplo   the strain's reads phase into haplotypes > 10 bp apart --
                                   two close organisms (greedy read clustering only separates
                                   organisms >~100 bp apart, so these are invisible above)
                     one_mono      in one 20260907 mono only: may have entered on the
                                   destination plate, not the source
  indistinguishable groups
                   a CLIQUE cover at < T bp (haplotype-aware): every pair inside a group is
                   below T. Connected components are not used because they chain -- at
                   < 15 bp the largest component holds 116 wells up to 40 bp apart, which
                   is not "mutually indistinguishable". The cover is greedy and not unique;
                   `confusable_outside_group` counts partners < T bp left in other groups.

Outputs (04_contamination_distinguishability/outputs/):
  p00_well_table.csv            384 rows
  p00_secondary_organisms.csv   every secondary-organism observation with its evidence
  p00_indistinguishable_groups.csv  clique groups at each threshold
"""
import collections
import re
import pickle
import sys
from multiprocessing import Pool
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import config  # noqa: E402

import numpy as np
import pandas as pd
from Bio import SeqIO

import deep_consensus as dc

OUT = BASE / "04_contamination_distinguishability" / "outputs"
SEQ = config.SEQ
STRAIN_IDS = SEQ / "pairwise_interaction_experiments/analyses/20261113/messy/outputs/s00_strain_ids.csv"
ASSOC_20261026 = SEQ / "community_assembly_experiments/20261026/02_qc/outputs/q04_foreign_association.csv"
ROWS = "ABCDEFGHIJKLMNOP"
CLIQUE_THRESHOLDS = (8, 15)
SAME_ORGANISM_BP = 2       # same as stage 3's observation clustering
KNOWN_BP = 2               # a secondary population this close to a reference IS that organism
CLOSE_HAPLO_BP = 10


def strain_label(well, ids):
    return str(ids.get(well, well))


def organism_label(nearest, ties, bp, gid, ids):
    """How to name a secondary organism: its strain number if one reference matches,
    its identity group if several do (16S cannot say which), `unk` if none is close."""
    if bp > KNOWN_BP:
        return "unk"
    if len(ties) > 1 and gid.get(nearest):
        return gid[nearest]
    return strain_label(nearest, ids)


def _nearest(job):
    key, seq, refs = job
    d = sorted((dc.core_bp(seq, s), n) for n, s in refs.items())
    return key, d[0][0], d[0][1], [n for b, n in d if b == d[0][0]]


def secondary_from_monos(ids, gid, processes=48):
    with open(config.REFDB_OUT / "m02_mono_populations.pkl", "rb") as fh:
        pops = pickle.load(fh)
    wells = pd.read_csv(config.REFDB_OUT / "m01_mono_wells.csv").set_index("sample_id")
    refs = {r.id: str(r.seq) for r in SeqIO.parse(config.REFDB_OUT / "reference_db_20260907.fasta", "fasta")}
    jobs = [((sid, k), p["consensus"], refs) for sid, pl in pops.items() for k, p in enumerate(pl) if k > 0]
    with Pool(processes) as pool:
        hits = {k: (bp, n, t) for k, bp, n, t in pool.map(_nearest, jobs, chunksize=4)}
    rows = []
    for (sid, k), (bp, nearest, ties) in hits.items():
        p = pops[sid][k]
        rows.append(dict(strain=wells.at[sid, "strain"], sample_id=sid, pop=k, frac=p["frac"],
                         n_reads=p["n_reads"], nearest=nearest, nearest_bp=bp,
                         ties=" ".join(ties[:8]),
                         label=organism_label(nearest, ties, bp, gid, ids)))
    m = pd.DataFrame(rows)
    # same organism in both monos of a strain -> came from the source well
    m["both_monos"] = False
    for s, g in m.groupby("strain"):
        sids = sorted(g.sample_id.unique())
        if len(sids) < 2:
            continue
        for i in g.index[g.sample_id == sids[0]]:
            for j in g.index[g.sample_id == sids[1]]:
                a, b = m.loc[i], m.loc[j]
                if dc.core_bp(pops[a.sample_id][a["pop"]]["consensus"],
                              pops[b.sample_id][b["pop"]]["consensus"]) <= SAME_ORGANISM_BP:
                    m.loc[[i, j], "both_monos"] = True
    return m


def clique_cover(names, pairs, thr):
    """Greedy: seed with the unassigned strain that has most unassigned partners, then
    add partners that are < thr to EVERY current member."""
    e = pairs[pairs.bp < thr]
    nb = collections.defaultdict(set)
    for a, b in zip(e.strain_a, e.strain_b):
        nb[a].add(b); nb[b].add(a)
    left = {n for n in names if nb[n]}
    groups = []
    while left:
        v = max(sorted(left), key=lambda n: len(nb[n] & left))
        g = [v]
        cand = sorted(nb[v] & left, key=lambda c: (-len(nb[c] & nb[v] & left), c))
        for c in cand:
            if all(c in nb[x] for x in g):
                g.append(c)
        groups.append(sorted(g))
        left -= set(g)
    groups = [g for g in groups if len(g) > 1]
    groups.sort(key=lambda g: (-len(g), g[0]))
    gid = {s: f"{k + 1}" for k, g in enumerate(groups) for s in g}
    outside = {n: len([x for x in nb[n] if gid.get(x) != gid.get(n) or gid.get(n) is None])
               for n in names}
    return gid, outside, groups, nb


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ids_df = pd.read_csv(STRAIN_IDS)
    ids = dict(zip(ids_df.well, ids_df.strain_no.astype(str)))
    ref = pd.read_csv(config.REFDB_OUT / "s01_strain_reference.csv").set_index("strain")
    dist = pd.read_csv(config.DIST_OUT / "d01_strain_distinguishability.csv").set_index("strain")
    gid = dist.identity_group.dropna().to_dict()
    gs = pd.read_csv(config.GENOME_OUT / "g02_strain_genome_status.csv").set_index("strain")
    clus = pd.read_csv(config.GENOME_OUT / "g02_genome_clusters.csv").set_index("well")
    elig = pd.read_csv(config.GENOME_OUT / "g04_design_eligibility.csv").set_index("strain")
    pairs = pd.read_csv(config.DIST_OUT / "d01_pairs.csv")

    # --- invaders and their home wells ---------------------------------------------------
    def unit_label(token):
        return token if re.fullmatch(r"C\d{3}", token) else strain_label(token, ids)
    invader, invader_label = {}, {}
    for s, r in gs[gs.genome_status == "wrong"].iterrows():
        toks = str(r.genome_is_of).split() if isinstance(r.genome_is_of, str) else []
        if not toks or r.genome_is_of == "no genome on file":
            invader[s], invader_label[s] = "none on file", "none on file"
        else:
            invader[s] = " ".join(toks)
            invader_label[s] = " / ".join(unit_label(t) for t in toks)
    ok = set(gs.index[gs.genome_status == "ok"])
    home = collections.defaultdict(list)       # invader unit -> wells that hold it with its genome
    for w in ok:
        u = clus.clone_cluster.get(w)
        home[u if isinstance(u, str) else w].append(w)

    # --- secondary organisms -------------------------------------------------------------
    m = secondary_from_monos(ids, gid)
    sec = []
    for r in m.itertuples():
        if r.nearest == r.strain:          # the strain's OWN sequence as a minority: original
            continue                       # persisting under an overgrower (drawn on f01)
        sec.append(dict(strain=r.strain, evidence="both_monos" if r.both_monos else "one_mono",
                        label=r.label, organism=r.nearest, organism_bp=r.nearest_bp,
                        frac=round(r.frac, 4), sample_id=r.sample_id))
    a = pd.read_csv(ASSOC_20261026)
    a = a[a.verdict.str.startswith("tracks")]
    for (ent, mem), g in a.groupby(["entity", "best_member"]):
        lab = ent if ent.startswith("U") else (gid.get(ent) if gid.get(ent) and dist.n_identical.get(ent, 0) > 0
                                               else strain_label(ent, ids))
        sec.append(dict(strain=mem, evidence="run_20261026", label=lab, organism=ent,
                        organism_bp=None, frac=round(float(g.coverage.max()), 4), sample_id=None))
    for s, r in ref.iterrows():
        if r.tier != "none" and r.haplotype_divergence_bp > CLOSE_HAPLO_BP:
            sec.append(dict(strain=s, evidence="close_haplo", label=f"+{int(r.haplotype_divergence_bp)}bp",
                            organism=None, organism_bp=r.haplotype_divergence_bp, frac=None, sample_id=None))
    sec = pd.DataFrame(sec)
    rank = {"both_monos": 0, "run_20261026": 1, "close_haplo": 2, "one_mono": 3}
    sec["rank"] = sec.evidence.map(rank)
    sec = sec.sort_values(["strain", "rank", "frac"], ascending=[True, True, False])
    sec.drop(columns="rank").to_csv(OUT / "p00_secondary_organisms.csv", index=False)
    original_minor = set(m[m.nearest == m.strain].strain)

    # --- indistinguishable groups ----------------------------------------------------------
    names = sorted(ref.index[ref.tier != "none"])
    cliques, crow = {}, []
    for thr in CLIQUE_THRESHOLDS:
        g, outside, groups, _ = clique_cover(names, pairs, thr)
        cliques[thr] = (g, outside)
        for k, grp in enumerate(groups):
            sub = pairs[pairs.strain_a.isin(grp) & pairs.strain_b.isin(grp)]
            crow.append(dict(threshold_bp=thr, group=k + 1, size=len(grp), members=" ".join(grp),
                             strain_nos=" ".join(strain_label(x, ids) for x in grp),
                             max_internal_bp=int(sub.bp.max()) if len(sub) else 0))
    pd.DataFrame(crow).to_csv(OUT / "p00_indistinguishable_groups.csv", index=False)

    # --- one row per well -------------------------------------------------------------------
    rows = []
    for i, row in enumerate(ROWS):
        for col in range(1, 25):
            w = f"{row}{col}"
            has = w in ref.index and ref.at[w, "tier"] != "none"
            s_rows = sec[sec.strain == w] if has else sec.iloc[0:0]
            top = s_rows.iloc[0] if len(s_rows) else None
            r = dict(well=w, row=row, col=col, row_idx=i, strain_no=strain_label(w, ids),
                     n_wells_same_isolate=int(ids_df.set_index("well").n_wells_same_isolate.get(w, 1)),
                     isolate_name=ids_df.set_index("well").isolate_name.get(w),
                     has_reference=has)
            if has:
                r.update(tier=ref.at[w, "tier"], no_growth=bool(ref.at[w, "no_growth"]),
                         identity_group=gid.get(w), n_identical=int(dist.n_identical.get(w, 0)),
                         genome_status=gs.genome_status.get(w, "no_genome"),
                         clone_cluster=clus.clone_cluster.get(w) if w in clus.index else None,
                         clone_size=int(clus.clone_size.get(w)) if w in clus.index else None,
                         invader=invader.get(w), invader_label=invader_label.get(w),
                         original_still_minor=w in original_minor,
                         genomic_use=elig.genomic_use.get(w), exclude_reason=elig.exclude_reason.get(w),
                         secondary_evidence=top.evidence if top is not None else None,
                         secondary_label=top.label if top is not None else None,
                         secondary_frac=top.frac if top is not None else None,
                         secondary_all=" ; ".join(f"{x.label}({x.evidence})" for x in s_rows.itertuples()))
                for thr, (g, outside) in cliques.items():
                    r[f"group_lt{thr}bp"] = g.get(w)
                    r[f"confusable_outside_group_lt{thr}bp"] = outside.get(w, 0)
            rows.append(r)
    t = pd.DataFrame(rows)
    # home wells: wells whose own (verified) genome is an invader somewhere else
    t["home_of"] = None
    for u, ws in home.items():
        hits = [x for x in invader.values() if u in x.split()]
        if hits:
            for w in ws:
                t.loc[t.well == w, "home_of"] = u
    t.to_csv(OUT / "p00_well_table.csv", index=False)
    report(t, sec)
    return t


def report(t, sec):
    h = t[t.has_reference]
    print(f"{len(t)} wells, {len(h)} with a reference; genome status "
          f"{h.genome_status.value_counts().to_dict()}")
    w = h[h.genome_status == "wrong"]
    print(f"wrong-genome wells by invader: {w.invader_label.value_counts().head(8).to_dict()}")
    print(f"  ... of which the original still shows as a minority: {w.original_still_minor.sum()}")
    print(f"secondary-organism evidence (best per well): {h.secondary_evidence.value_counts().to_dict()}")
    for thr in CLIQUE_THRESHOLDS:
        g = h[f"group_lt{thr}bp"].dropna()
        print(f"< {thr} bp clique groups: {g.nunique()} groups covering {len(g)} wells, "
              f"largest {g.value_counts().head(5).tolist()}")


if __name__ == "__main__":
    main()
