"""Stage 8: for which wells does the genome on file describe what is actually in the well?

Per well with a 20260907 reference and a genome:

  own 16S call    the well's sequence (consensus or any operon haplotype, whichever fits
                  best) against its OWN genome. HSPs are tiled along the query, because
                  short-read assemblies often break the 16S where operons collapse.
                    match      <= MATCH_DIV divergence over >= MIN_TILED_BP
                    mismatch   >= MISMATCH_DIV divergence over >= MIN_TILED_BP
                    ambiguous  in between
                    no_16S     the assembly carries no usable 16S
  mixed assembly  9% of assemblies carry a second organism's rRNA operon (g03). A match
                  counts only if it is to the assembly's MAIN organism; matching only the
                  contaminant's operon (E22: own 16S at 9.3x bulk coverage, the well's at
                  0.21x) means the bulk genome is something else -> wrong.
  clone vouching  genomes >= CLONE_ANI apart are one organism (45 clusters, min within-
                  cluster ANI 99.97%, alignment fraction >= 94%), so a clone-mate's 16S
                  stands in when the own assembly has none.
  ANI exclusion   still undecided, but the well's 16S IS carried by other genomes, and the
                  own genome is < FOREIGN_ANI to every one of them (skani reports nothing
                  below ~80%): an organism that distant cannot carry an identical 16S, so the
                  genome is not the well's organism. Between FOREIGN_ANI and CLONE_ANI it
                  stays unverified -- identical 16S across close species does happen.
                  Validated on wells where the 16S answer is known: flags 15/16 mismatches
                  and 2/2 contaminant-operon cases, and 2/158 matches (G22, a carrier at
                  87.8% ANI, same genus; L6, a mixed well). All 35 calls it makes are
                  < 80% ANI -- different genera -- so the exact cut does not matter.

  genome_status   ok          the genome (or its clone) carries the well's 16S
                  wrong       it carries a clearly different 16S: the well holds another
                              organism -- overgrowth, contamination or a label error. The KO
                              vector of this well describes something that is not there.
                  unverified  no usable 16S in the genome or its clones, or ambiguous
                  no_genome

Two further facts per well, both about what 16S CANNOT see:
  clone_cluster          wells whose genomes are one organism. For genome-based ML these are
                         one data point, not several, and cv_strain must hold them out
                         together or it leaks.
  distinct_genomes_with_this_16S
                         how many genomically distinct organisms (clone clusters) carry this
                         well's 16S. >1 means a swap between those wells would be invisible
                         to a 16S readout; status `ok` then means "consistent", not "proven".

Outputs:
  g02_genome_clusters.csv        well -> clone cluster, size, members
  g02_strain_genome_status.csv   per strain: the calls above + the genomes that DO carry its 16S
  g02_identity_group_genomes.csv per 16S identity group: how many genomes it really holds
"""
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import config  # noqa: E402

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

OUT = config.GENOME_OUT
CLONE_ANI = 99.9
MATCH_DIV = 0.002          # observed: 222 of 261 own-genome hits are 0 diffs, then a gap
MISMATCH_DIV = 0.01        #   from ~0.5% to >1%; >1% is ~15 bp, a different organism
MIN_TILED_BP = 400
FOREIGN_ANI = 85.0
MAX_OVERLAP = 20
CARRIER_MIN_BP = 1000      # a genome "carries" a 16S only over near-full length: conserved
                           # 400-800 bp stretches match across distant taxa (A14/A16, K19/E6)


def clone_clusters(idx, ani):
    W = list(idx.well); wi = {w: i for i, w in enumerate(W)}
    hi = ani[ani.ani >= CLONE_ANI]
    A = coo_matrix((np.ones(len(hi)), ([wi[a] for a in hi.well_a], [wi[b] for b in hi.well_b])),
                   shape=(len(W), len(W)))
    _, lab = connected_components(A, directed=False)
    sizes = pd.Series(lab).value_counts()
    # name clusters by size, then by first member, so ids are stable across reruns
    groups = {}
    for w, l in zip(W, lab):
        groups.setdefault(l, []).append(w)
    order = sorted(groups.values(), key=lambda g: (-len(g), sorted(g)[0]))
    rows = []
    for k, g in enumerate(order):
        cid = f"C{k + 1:03d}" if len(g) > 1 else None
        for w in g:
            rows.append(dict(well=w, clone_cluster=cid, clone_size=len(g),
                             clone_members=" ".join(sorted(g)) if len(g) > 1 else ""))
    return pd.DataFrame(rows)


def tile(h):
    """Greedy, best-first, non-overlapping HSPs along the query: (tiled_bp, diffs)."""
    taken, L, D = [], 0, 0
    for r in h.sort_values("bitscore", ascending=False).itertuples():
        lo, hi = min(r.qstart, r.qend), max(r.qstart, r.qend)
        if all(min(hi, b) - max(lo, a) < MAX_OVERLAP for a, b in taken):
            taken.append((lo, hi)); L += r.length; D += r.diffs
    return L, D


def per_strain_genome(blast, purity, organism_of):
    """Best fit of each strain's 16S (consensus or any haplotype) to each genome, kept
    separately for the genome's MAIN organism and for a contaminant organism in a mixed
    assembly (g03). A match to a contaminant's operon says the well's organism was in the
    sequenced culture, not that the genome -- and its KO vector -- is that organism."""
    blast = blast.copy()
    blast["strain"] = blast.qseqid.str.replace(r"_h\d+$", "", regex=True)
    blast["contig"] = blast.sseqid.str.split("|").str[1]
    org = [organism_of.get((g, c), -1) for g, c in zip(blast.genome, blast.contig)]

    def role(g, o):
        p = purity.get(g)
        if p is None or p["n_organisms"] <= 1:
            return "main"
        if p["main_called_by"] != "coverage" or o == -1:
            return "unresolved"
        return "main" if o == p["main_organism"] else "other"
    blast["role"] = [role(g, o) for g, o in zip(blast.genome, org)]
    rows = []
    for (q, s, g, rl), h in blast.groupby(["qseqid", "strain", "genome", "role"], sort=False):
        L, D = tile(h)
        rows.append((s, g, rl, q, L, D))
    t = pd.DataFrame(rows, columns=["strain", "genome", "role", "query", "tiled_bp", "diffs"])
    t["div"] = t.diffs / t.tiled_bp
    t["usable"] = t.tiled_bp >= MIN_TILED_BP
    # usable tiles first, then lowest divergence, then longest
    t = t.sort_values(["usable", "div", "tiled_bp"], ascending=[False, True, False])
    return t.groupby(["strain", "genome", "role"], as_index=False).first()


def call(row):
    if row is None or not row.usable:
        return "no_16S"
    if row["div"] <= MATCH_DIV:
        return "match"
    if row["div"] >= MISMATCH_DIV:
        return "mismatch"
    return "ambiguous"


def main():
    idx = pd.read_csv(OUT / "g01_genome_index.csv")
    ani = pd.read_csv(OUT / "g01_ani_pairs.csv")
    blast = pd.read_csv(OUT / "g01_blast_16S_vs_genomes.csv.gz")
    ref = pd.read_csv(config.REFDB_OUT / "s01_strain_reference.csv")
    dist = pd.read_csv(config.DIST_OUT / "d01_strain_distinguishability.csv").set_index("strain")

    cl = clone_clusters(idx, ani)
    cl.merge(idx[["well", "genome_table_id", "assembly_name"]], on="well").to_csv(
        OUT / "g02_genome_clusters.csv", index=False)
    clone_of = dict(zip(cl.well, cl.clone_cluster))
    ani_of = {(a, b): v for a, b, v in zip(ani.well_a, ani.well_b, ani.ani)}
    ani_of.update({(b, a): v for (a, b), v in ani_of.items()})
    mates = {w: [m for m in cl.clone_members[cl.well == w].iloc[0].split() if m != w]
             for w in cl.well}

    pur = pd.read_csv(OUT / "g03_assembly_purity.csv")
    purity = {r.well: dict(n_organisms=r.n_organisms, main_organism=r.main_organism,
                           main_called_by=r.main_called_by) for r in pur.itertuples()}
    genes = pd.read_csv(OUT / "g03_assembly_16S_genes.csv")
    organism_of = {(w, c): o for w, c, o in zip(genes.well, genes.contig, genes.organism)}
    sg = per_strain_genome(blast, purity, organism_of)
    sg["call"] = [call(r) for _, r in sg.iterrows()]
    fit = sg.set_index(["strain", "genome", "role"])

    def fit_of(s, g, role="main"):
        return fit.loc[(s, g, role)] if (s, g, role) in fit.index else None

    rows = []
    has_genome = set(idx.well)
    for r in ref[ref.tier != "none"].itertuples():
        s = r.strain
        carriers = sg[(sg.strain == s) & (sg.role == "main") & (sg.call == "match")
                      & (sg.tiled_bp >= CARRIER_MIN_BP)].genome.tolist()
        carrier_units = {clone_of.get(g) or g for g in carriers}
        row = dict(strain=s, tier=r.tier, status=r.status, no_growth=r.no_growth,
                   identity_group=dist.identity_group.get(s),
                   genomes_carrying_16S=" ".join(sorted(carriers)),
                   distinct_genomes_with_this_16S=len(carrier_units))
        if s not in has_genome:
            rows.append(dict(row, genome_status="no_genome")); continue
        own = fit_of(s, s, "main")
        own_call = call(own)
        if own_call != "match":
            if call(fit_of(s, s, "other")) == "match":
                own_call = "match_contaminant_operon"
            elif call(fit_of(s, s, "unresolved")) == "match":
                own_call = "match_unresolved_assembly"
        p = purity.get(s, {})
        row.update(clone_cluster=clone_of.get(s), clone_size=len(mates[s]) + 1,
                   assembly_organisms=p.get("n_organisms"),
                   assembly_main_called_by=p.get("main_called_by"),
                   own_call=own_call,
                   own_div_pct=round(100 * own["div"], 3) if own is not None else None,
                   own_tiled_bp=int(own.tiled_bp) if own is not None else 0,
                   own_best_query=own["query"] if own is not None else None)
        via = None
        if own_call in ("no_16S", "ambiguous") and mates[s]:
            mc = [call(fit_of(s, m)) for m in mates[s]]
            if "match" in mc:
                via = "match"
            elif "mismatch" in mc:
                via = "mismatch"
        row["clone_call"] = via
        others = [g for g in carriers if g != s]
        row["max_ani_to_carriers"] = max((ani_of.get((s, g), 0.0) for g in others), default=None)
        if via is None and own_call in ("no_16S", "ambiguous") and others \
                and row["max_ani_to_carriers"] < FOREIGN_ANI:
            via = "foreign_by_ani"
        row["evidence"] = ("own_16S" if own_call in ("match", "mismatch", "match_contaminant_operon")
                           else "clone_16S" if via in ("match", "mismatch")
                           else "ani" if via == "foreign_by_ani" else "none")
        eff = ("mismatch" if via == "foreign_by_ani" or own_call == "match_contaminant_operon"
               else (via or own_call))
        row["genome_status"] = {"match": "ok", "mismatch": "wrong"}.get(eff, "unverified")
        if row["genome_status"] == "wrong":
            row["genome_is_of"] = " ".join(
                sorted({clone_of.get(g) or g for g in carriers})) or "no genome on file"
        rows.append(row)
    st = pd.DataFrame(rows)
    st.to_csv(OUT / "g02_strain_genome_status.csv", index=False)

    groups = pd.read_csv(config.DIST_OUT / "d01_identity_groups.csv")
    grow = []
    for g in groups.itertuples():
        mem = g.members.split()
        sub = st[st.strain.isin(mem)]
        ok = sub[sub.genome_status == "ok"]
        units = {u if isinstance(u, str) else w for u, w in zip(ok.clone_cluster, ok.strain)}
        grow.append(dict(group=g.group, size=g.size, members=g.members,
                         with_genome=int(sub.genome_status.ne("no_genome").sum()),
                         ok=int((sub.genome_status == "ok").sum()),
                         wrong=int((sub.genome_status == "wrong").sum()),
                         unverified=int((sub.genome_status == "unverified").sum()),
                         distinct_genomes_among_ok=len(units),
                         ok_members=" ".join(ok.strain),
                         wrong_members=" ".join(sub[sub.genome_status == "wrong"].strain)))
    gdf = pd.DataFrame(grow)
    gdf.to_csv(OUT / "g02_identity_group_genomes.csv", index=False)
    report(st, gdf, cl)
    return st, gdf


def report(st, gdf, cl):
    print(f"{len(cl)} genomes on file; {cl.clone_cluster.nunique()} clone clusters cover "
          f"{cl.clone_cluster.notna().sum()} of them -> "
          f"{cl.clone_cluster.nunique() + cl.clone_cluster.isna().sum()} distinct organisms")
    print("\ngenome status by tier:")
    print(pd.crosstab(st.tier, st.genome_status, margins=True).to_string())
    g = st[st.genome_status != "no_genome"]
    print(f"\nown-genome 16S call: {g.own_call.value_counts().to_dict()}")
    print(f"evidence behind each status: "
          f"{g.groupby(['genome_status', 'evidence']).size().to_dict()}")
    ok = st[st.genome_status == "ok"]
    print(f"\nok wells: {len(ok)}, but as distinct genomes: "
          f"{ok.clone_cluster.nunique() + ok.clone_cluster.isna().sum()}")
    print(f"ok wells whose 16S is also carried by a genomically distinct organism: "
          f"{(ok.distinct_genomes_with_this_16S > 1).sum()}")
    print("\nlargest identity groups:")
    print(gdf.head(8)[["group", "size", "with_genome", "ok", "wrong", "unverified",
                       "distinct_genomes_among_ok"]].to_string(index=False))


if __name__ == "__main__":
    main()
