"""Per experiment and per well: which organism is the genome on file, which organism is in the
well, and does the genome describe the well?

The 20261021 genome concordance answered this for the 20260907 source plate. 20260630 used the
plate as it was in July, and the two disagree for some wells (analyses/20260925 found 14 of 82
shared names to be different organisms, E4 among them). So 20260630 gets its OWN verdict here:
its consensus2 16S references -- the sequences its interaction labels were scored against --
blasted against all 298 assemblies, then judged with 20261021's g02 rules unchanged
(tiled HSPs, match <= 0.2% / mismatch >= 1% divergence over >= 400 bp, mixed-assembly roles,
clone vouching, ANI exclusion).

Columns that matter downstream:
  genome_unit     the genome on file, as an organism: its >= 99.9% ANI clone cluster
                  (C001...) or the well itself. Wells sharing one carry the SAME KO vector.
  genome_status   ok / wrong / unverified -- does the genome carry the well's 16S
  genome_is_of    for `wrong` wells, the genome units that DO carry the well's 16S
  cv_group        the unit a leak-free CV must hold out together (see below)

cv_group = organism_unit, the organism physically in the well. Leakage needs a test strain
whose feature vector AND phenotype are both present in training:
  ok / unverified wells   the genome on file is (or may be) the organism: genome_unit. Same
                          unit -> identical KO vector and the same organism -> the leak.
  wrong wells             the organism is the invader: its genome unit when ONE unit carries
                          the well's 16S (the 9 wells overgrown by C001 join C001), else a
                          group keyed by the set of carriers (wells overgrown by the same
                          unidentified 16S type stay together), else the well alone. A wrong
                          well is NOT tied to its own genome_unit: there it shares the KO
                          vector but not the phenotype, which costs the model rather than
                          flattering it -- and tying it in would chain clusters (A10 carries
                          C002's genome and holds C001, which would merge 37 wells into one
                          block).
Two alternatives are kept for sensitivity:
  cv_group_genome   genome_unit alone (what the 20261021 README prescribed)
  cv_group_union    connected components of genome_unit + invader links (most conservative)

Outputs (outputs/):
  s01_blast_20260630_vs_genomes.csv.gz
  s01_genome_status_20260630.csv     20260630's own verdicts, with the 20260907 verdict beside
  s01_wells_<exp>.csv                per modelled well: strain_no, genome_unit, status, cv_group
"""
import importlib.util
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

import edlib
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

import config

SCRATCH = Path("/tmp/claude-1000/-home-rl-scripts-karl/2c0eda65-b5b0-4305-b7bc-1ea751f6759d/scratchpad")
CACHE = SCRATCH / "genome_fna_cache"
OUT = config.OUT


def load_g02():
    """Import 20261021's g02_verdicts with ITS config, then put ours back."""
    mine = sys.modules.pop("config")
    sys.path.insert(0, str(config.ANALYSES / "20261021"))
    try:
        spec = importlib.util.spec_from_file_location("g02_verdicts", config.CONC / "g02_verdicts.py")
        g02 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(g02)
    finally:
        sys.path.remove(str(config.ANALYSES / "20261021"))
        sys.modules["config"] = mine
    return g02


def _extract(args):
    well, gbk = args
    lines = []
    for rec in SeqIO.parse(config.GENOME_GBK_DIR / f"{gbk}.gbk", "genbank"):
        s = str(rec.seq)
        if s and set(s) != {"N"}:
            lines.append(f">{well}|{rec.id}\n{s}\n")
    return "".join(lines)


def blast_db():
    db = CACHE / "all_genomes"
    if (CACHE / "all_genomes.nsq").exists() or (CACHE / "all_genomes.00.nsq").exists():
        return db
    CACHE.mkdir(parents=True, exist_ok=True)
    idx = pd.read_csv(config.CONC_OUT / "g01_genome_index.csv")
    with Pool(32) as pool:
        chunks = pool.map(_extract, list(zip(idx.well, idx.gbk)))
    (CACHE / "all_genomes.fna").write_text("".join(chunks))
    subprocess.run([str(config.BLAST_BIN / "makeblastdb"), "-in", str(CACHE / "all_genomes.fna"),
                    "-dbtype", "nucl", "-out", str(db)], check=True, capture_output=True)
    return db


def blast_0630():
    out = OUT / "s01_blast_20260630_vs_genomes.csv.gz"
    if out.exists():
        return pd.read_csv(out)
    db = blast_db()
    q = CACHE / "q_20260630.fa"
    with open(q, "w") as fh:
        for r in SeqIO.parse(config.EXPERIMENTS["20260630"]["ref"], "fasta"):
            fh.write(f">{r.id}\n{r.seq}\n")
    cols = "qseqid sseqid pident length qlen mismatch gaps qstart qend bitscore"
    raw = CACHE / "blast_0630.tsv"
    # identical settings to 20261021/g01
    subprocess.run([str(config.BLAST_BIN / "blastn"), "-query", str(q), "-db", str(db),
                    "-outfmt", f"6 {cols}", "-perc_identity", "85", "-max_target_seqs", "100000",
                    "-max_hsps", "20", "-num_threads", "48", "-out", str(raw)], check=True)
    b = pd.read_csv(raw, sep="\t", names=cols.split())
    b["genome"] = b.sseqid.str.split("|").str[0]
    b["diffs"] = b.mismatch + b.gaps
    b = b[["qseqid", "genome", "sseqid", "pident", "length", "qlen", "diffs", "qstart", "qend", "bitscore"]]
    b.to_csv(out, index=False)
    return b


def verdicts_0630(g02, blast):
    """g02.main()'s per-well loop, for 20260630's references (no tiers, no haplotypes)."""
    idx = pd.read_csv(config.CONC_OUT / "g01_genome_index.csv")
    ani = pd.read_csv(config.CONC_OUT / "g01_ani_pairs.csv")
    cl = pd.read_csv(config.CONC_OUT / "g02_genome_clusters.csv")
    clone_of = dict(zip(cl.well, cl.clone_cluster.where(cl.clone_cluster.notna(), None)))
    mates = {w: [m for m in str(mm).split() if m != w] if isinstance(mm, str) else []
             for w, mm in zip(cl.well, cl.clone_members)}
    ani_of = {(a, b): v for a, b, v in zip(ani.well_a, ani.well_b, ani.ani)}
    ani_of.update({(b, a): v for (a, b), v in ani_of.items()})
    pur = pd.read_csv(config.CONC_OUT / "g03_assembly_purity.csv")
    purity = {r.well: dict(n_organisms=r.n_organisms, main_organism=r.main_organism,
                           main_called_by=r.main_called_by) for r in pur.itertuples()}
    genes = pd.read_csv(config.CONC_OUT / "g03_assembly_16S_genes.csv")
    organism_of = {(w, c): o for w, c, o in zip(genes.well, genes.contig, genes.organism)}

    sg = g02.per_strain_genome(blast, purity, organism_of)
    sg["call"] = [g02.call(r) for _, r in sg.iterrows()]
    fit = sg.set_index(["strain", "genome", "role"])

    def fit_of(s, g, role="main"):
        return fit.loc[(s, g, role)] if (s, g, role) in fit.index else None

    has_genome = set(idx.well)
    rows = []
    for s in sorted(set(blast.qseqid) | {r.id for r in SeqIO.parse(config.EXPERIMENTS["20260630"]["ref"], "fasta")}):
        carriers = sg[(sg.strain == s) & (sg.role == "main") & (sg.call == "match")
                      & (sg.tiled_bp >= g02.CARRIER_MIN_BP)].genome.tolist()
        row = dict(strain=s, genomes_carrying_16S=" ".join(sorted(carriers)),
                   distinct_genomes_with_this_16S=len({clone_of.get(g) or g for g in carriers}))
        if s not in has_genome:
            rows.append(dict(row, genome_status="no_genome")); continue
        own = fit_of(s, s, "main")
        own_call = g02.call(own)
        if own_call != "match":
            if g02.call(fit_of(s, s, "other")) == "match":
                own_call = "match_contaminant_operon"
            elif g02.call(fit_of(s, s, "unresolved")) == "match":
                own_call = "match_unresolved_assembly"
        row.update(own_call=own_call,
                   own_div_pct=round(100 * own["div"], 3) if own is not None else None,
                   own_tiled_bp=int(own.tiled_bp) if own is not None else 0)
        via = None
        if own_call in ("no_16S", "ambiguous") and mates.get(s):
            mc = [g02.call(fit_of(s, m)) for m in mates[s]]
            via = "match" if "match" in mc else ("mismatch" if "mismatch" in mc else None)
        others = [g for g in carriers if g != s]
        row["max_ani_to_carriers"] = max((ani_of.get((s, g), 0.0) for g in others), default=None)
        if via is None and own_call in ("no_16S", "ambiguous") and others \
                and row["max_ani_to_carriers"] < g02.FOREIGN_ANI:
            via = "foreign_by_ani"
        row["evidence"] = ("own_16S" if own_call in ("match", "mismatch", "match_contaminant_operon")
                           else "clone_16S" if via in ("match", "mismatch")
                           else "ani" if via == "foreign_by_ani" else "none")
        eff = ("mismatch" if via == "foreign_by_ani" or own_call == "match_contaminant_operon"
               else (via or own_call))
        row["genome_status"] = {"match": "ok", "mismatch": "wrong"}.get(eff, "unverified")
        if row["genome_status"] == "wrong":
            row["genome_is_of"] = " ".join(sorted({clone_of.get(g) or g for g in carriers})) or "no genome on file"
        rows.append(row)
    return pd.DataFrame(rows)


def _seqs(path):
    return {r.id.split()[0]: str(r.seq).upper() for r in SeqIO.parse(path, "fasta")}


def bp_distance(a, b, trim=20):
    """20261021/d01's distance: ends trimmed, smaller of the two infix directions."""
    a, b = a[trim:-trim], b[trim:-trim]
    return min(edlib.align(a, b, mode="HW", task="distance")["editDistance"],
               edlib.align(b, a, mode="HW", task="distance")["editDistance"])


def cv_groups(wells, genome_unit, is_of):
    """Connected components over wells: same genome_unit, or a wrong well -> its one invader."""
    nodes = sorted(set(wells) | set(genome_unit.values()) | {v for v in is_of.values() if v})
    ix = {n: i for i, n in enumerate(nodes)}
    e = [(ix[w], ix[genome_unit[w]]) for w in wells]
    e += [(ix[w], ix[is_of[w]]) for w in wells if is_of.get(w)]
    A = coo_matrix((np.ones(len(e)), tuple(zip(*e))), shape=(len(nodes), len(nodes)))
    _, lab = connected_components(A, directed=False)
    comp = {w: lab[ix[w]] for w in wells}
    # name each group by its largest member unit, so ids read as organisms
    name = {}
    for w in wells:
        name.setdefault(comp[w], []).append(genome_unit[w])
    name = {k: pd.Series(v).value_counts().index[0] for k, v in name.items()}
    return {w: f"G:{name[comp[w]]}" for w in wells}


def well_table(exp, status):
    """Per modelled well of one experiment."""
    pub = pd.read_csv(config.EXPERIMENTS[exp]["pub"] / "g01_modeling_pairs.csv")
    wells = sorted(set(pub.strain_a) | set(pub.strain_b))
    cl = pd.read_csv(config.CONC_OUT / "g02_genome_clusters.csv")
    ids = pd.read_csv(config.OUT / "s00_strain_ids.csv", dtype=str)
    t = pd.DataFrame({"well": wells})
    t = t.merge(ids[["well", "strain_no", "isolate_name", "genome_table_id"]], on="well", how="left")
    t = t.merge(cl[["well", "clone_cluster", "clone_size"]], on="well", how="left")
    t["genome_unit"] = t.clone_cluster.fillna(t.well)
    st = status.set_index("strain")
    for c in ("genome_status", "genome_is_of", "evidence", "distinct_genomes_with_this_16S"):
        t[c] = t.well.map(st[c]) if c in st else None
    t["genome_status"] = t.genome_status.fillna("not_assessed")

    # an invader named as a WELL may itself sit in a clone cluster: resolve to its unit
    unit_of = dict(zip(cl.well, cl.clone_cluster.fillna(cl.well)))

    def invader(v):
        if not isinstance(v, str) or v == "no genome on file":
            return None, None
        units = sorted({unit_of.get(u, u) for u in v.split()})
        return (units[0], units[0]) if len(units) == 1 else (None, "16S:" + "+".join(units))
    inv = [invader(v) if s == "wrong" else (None, None) for s, v in zip(t.genome_status, t.genome_is_of)]
    t["invader_unit"] = [a for a, _ in inv]
    t["organism_unit"] = [
        gu if s != "wrong" else (a or b or f"well:{w}")
        for w, s, gu, (a, b) in zip(t.well, t.genome_status, t.genome_unit, inv)]
    t["cv_group"] = "O:" + t.organism_unit
    t["cv_group_genome"] = "G:" + t.genome_unit
    gu = dict(zip(t.well, t.genome_unit))
    t["cv_group_union"] = t.well.map(cv_groups(list(t.well), gu, dict(zip(t.well, t.invader_unit))))
    return t


def main():
    g02 = load_g02()
    blast = blast_0630()
    v30 = verdicts_0630(g02, blast)

    # how does each 20260630 reference compare with the 20260907 one for the same well?
    s30, s07 = _seqs(config.EXPERIMENTS["20260630"]["ref"]), _seqs(config.REFDB_0907)
    v30["bp_to_20260907_ref"] = [bp_distance(s30[w], s07[w]) if w in s30 and w in s07 else np.nan
                                 for w in v30.strain]
    g07 = pd.read_csv(config.CONC_OUT / "g02_strain_genome_status.csv").set_index("strain")
    v30["genome_status_20260907"] = v30.strain.map(g07.genome_status)
    v30.to_csv(OUT / "s01_genome_status_20260630.csv", index=False)

    print("20260630 own verdicts:", v30.genome_status.value_counts().to_dict())
    print("\nagreement with the 20260907 verdict, by whether the two references match (<=5 bp):")
    v30["same_16S_as_0907"] = np.where(v30.bp_to_20260907_ref.isna(), "no 0907 ref",
                                       np.where(v30.bp_to_20260907_ref <= 5, "same 16S", "different 16S"))
    print(pd.crosstab([v30.same_16S_as_0907, v30.genome_status], v30.genome_status_20260907.fillna("-")).to_string())

    g07r = g07.reset_index()
    for exp, status in (("20260907", g07r), ("20260630", v30)):
        t = well_table(exp, status)
        t.to_csv(OUT / f"s01_wells_{exp}.csv", index=False)
        n_units = t.genome_unit.nunique()
        print(f"\n=== {exp}: {len(t)} modelled wells = {n_units} distinct genomes")
        print("genome status:", t.genome_status.value_counts().to_dict())
        for c in ("cv_group", "cv_group_genome", "cv_group_union"):
            sizes = t[c].value_counts()
            print(f"  {c:16s} {sizes.size:4d} groups, largest {sizes.head(6).to_dict()}")
        print("wrong wells assigned to an identified invader:", t.invader_unit.notna().sum(),
              "| to an unidentified 16S type:", t.organism_unit.str.startswith("16S:").sum())


if __name__ == "__main__":
    main()
