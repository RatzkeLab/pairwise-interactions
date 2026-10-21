"""Stage 9: does each assembly hold ONE organism?

Found while checking G01: E21's and E22's assemblies share an identical full-length 16S,
yet skani puts them below 80% ANI and they share only ~4.8 kb at >=98% identity -- one
rRNA operon and nothing else. That is the signature of a contaminant in the culture that
was sequenced: rRNA operons are multi-copy, so a contaminant at a few % still has enough
operon coverage to assemble while the rest of its genome does not. Such an assembly can
make a well look `ok` (its 16S is "in the genome") when the bulk genome -- and its KO
vector -- is another organism.

Per assembly:
  1. every 16S gene, by barrnap;
  2. a coordinate frame: the longest 16S if >= 1300 bp, else the nearest full-length
     20260907 reference -- used ONLY to place genes, never to judge them;
  3. two genes are compared only where they actually overlap in that frame (>= MIN_OVERLAP
     bp), base for base, and are different organisms if they differ by > ORGANISM_DIV
     there. Short-read assemblies break collapsed operons into fragments covering
     different regions; comparing fragments against each other, or each against a
     reference, invented 50-87 "mixed" assemblies (P15/F5/F7/P17: three fragments of one
     16S, 447-520 bp). Intragenomic operon variation is <1% (259/294 assemblies at 0);
     E22 shows two organisms only 2.0% apart that skani puts below 80% ANI, so the cut
     sits at 1.5%;
  4. for Or's own assemblies (P*), contig coverage from SPAdes: GenBank contigs map 1:1,
     in order, onto pilon.fasta's NODE_<n>_length_<L>_cov_<C> -- checked by length.
     A 16S on contigs far below the bulk coverage is the contaminant's.

Outputs:
  g03_assembly_16S_genes.csv   every 16S gene: contig, length, coverage ratio, organism
  g03_assembly_purity.csv      per assembly: n organisms, which is main, how sure
"""
import glob
import os
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import config  # noqa: E402

import edlib
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

OUT = config.GENOME_OUT
FNA = config.GENOME_FNA_CACHE / "fna"
BARR = config.GENOME_FNA_CACHE / "barrnap"
ASSEMBLIES = config.GENOME_GBK_DIR.parents[1] / "assmebled_and_proccessed_genomes"
MIN_GENE_BP = 400
ORGANISM_DIV = 0.015
MIN_OVERLAP = 300
PARTIAL_TRIM = 15        # barrnap's partial hits can run into the junk at a contig end
BULK_MIN_CONTIG = 10_000


def _barrnap(well):
    gff = BARR / f"{well}.gff"
    if not gff.exists() or gff.stat().st_size == 0:
        env = dict(os.environ, PATH=f"{config.BARRNAP_BIN}:{os.environ['PATH']}")
        with open(gff, "w") as fh:
            subprocess.run([str(config.BARRNAP_BIN / "barrnap"), "--kingdom", "bac",
                            "--threads", "1", "--quiet", str(FNA / f"{well}.fna")],
                           stdout=fh, stderr=subprocess.DEVNULL, env=env, check=True)
    return well


def contig_coverage(assembly, gbk):
    """GenBank contig id -> SPAdes coverage, for assemblies that carry it."""
    pil = ASSEMBLIES / assembly / "spades_assembled" / "pilon.fasta"
    if not pil.exists():
        return {}
    p = [(r.id, len(r.seq)) for r in SeqIO.parse(pil, "fasta")]
    g = [(r.id, len(r.seq)) for r in SeqIO.parse(config.GENOME_GBK_DIR / f"{gbk}.gbk", "genbank")]
    if len(p) != len(g) or any(a[1] != b[1] for a, b in zip(p, g)):
        return {}
    return {gid: (float(pid.split("_cov_")[1].split("_")[0]), L) for (pid, L), (gid, _) in zip(p, g)}


def _div_into(gene, frame):
    r = edlib.align(gene, frame, mode="HW", task="distance")
    return r["editDistance"] / len(gene)


def genes_of(well):
    seqs = {r.id.split("|", 1)[1]: str(r.seq) for r in SeqIO.parse(FNA / f"{well}.fna", "fasta")}
    out = []
    for line in open(BARR / f"{well}.gff"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or "16S_rRNA" not in f[8]:
            continue
        contig = f[0].split("|", 1)[1]
        s = seqs[contig][int(f[3]) - 1:int(f[4])]
        if f[6] == "-":
            s = str(Seq(s).reverse_complement())
        if len(s) >= MIN_GENE_BP:
            out.append(dict(well=well, contig=contig, length=len(s), seq=s))
    return out


def _place(gene, frame):
    """gene -> (frame_start, frame_end, {frame_pos: gene_pos}) via an infix alignment."""
    r = edlib.align(gene, frame, mode="HW", task="path")
    import re
    fpos = r["locations"][0][0]
    gpos, m = 0, {}
    for n, op in re.findall(r"(\d+)([=XID])", r["cigar"]):
        n = int(n)
        if op in "=X":
            for k in range(n):
                m[fpos + k] = gpos + k
            fpos += n; gpos += n
        elif op == "I":
            gpos += n
        else:
            fpos += n
    return min(m), max(m), m


def _overlap_div(pa, pb, sa, sb):
    lo, hi = max(pa[0], pb[0]), min(pa[1], pb[1])
    if hi - lo < MIN_OVERLAP:
        return None
    ka = [pa[2][f] for f in range(lo, hi + 1) if f in pa[2]]
    kb = [pb[2][f] for f in range(lo, hi + 1) if f in pb[2]]
    if not ka or not kb:
        return None
    x, y = sa[min(ka):max(ka) + 1], sb[min(kb):max(kb) + 1]
    return edlib.align(x, y, mode="NW", task="distance")["editDistance"] / max(len(x), len(y))


def assign_organisms(genes, refs):
    """Organism 0 = the longest gene's; a gene joins an organism it overlaps within
    ORGANISM_DIV, founds a new one if it overlaps an existing organism and differs from
    all of them, and is left with organism 0 (`overlap=False`) if it overlaps nothing --
    a fragment that cannot be compared is not evidence of a second organism."""
    for g in genes:
        if g["length"] < 1300:
            g["seq"] = g["seq"][PARTIAL_TRIM:-PARTIAL_TRIM]
    order = sorted(range(len(genes)), key=lambda i: -genes[i]["length"])
    seed = genes[order[0]]["seq"]
    frame = seed if len(seed) >= 1300 else min(
        refs.values(), key=lambda r: edlib.align(seed, r, mode="HW", task="distance")["editDistance"])
    place = {i: _place(genes[i]["seq"], frame) for i in order}
    orgs = []                       # list of member index lists
    for i in order:
        best, compared = None, False
        for k, mem in enumerate(orgs):
            ds = [d for j in mem if (d := _overlap_div(place[i], place[j], genes[i]["seq"],
                                                       genes[j]["seq"])) is not None]
            if ds:
                compared = True
                if min(ds) <= ORGANISM_DIV and (best is None or min(ds) < best[1]):
                    best = (k, min(ds))
        if best is not None:
            orgs[best[0]].append(i); genes[i]["overlap"] = True
            genes[i]["div_to_organism"] = round(best[1], 4)
        elif compared or not orgs:
            orgs.append([i]); genes[i]["overlap"] = bool(compared)
            genes[i]["div_to_organism"] = 0.0
        else:
            orgs[0].append(i); genes[i]["overlap"] = False
            genes[i]["div_to_organism"] = None
    for k, mem in enumerate(orgs):
        for i in mem:
            genes[i]["organism"] = k
    return genes


def main(processes=48):
    BARR.mkdir(parents=True, exist_ok=True)
    idx = pd.read_csv(OUT / "g01_genome_index.csv")
    with Pool(processes) as pool:
        pool.map(_barrnap, list(idx.well))
    refs = {r.id: str(r.seq) for r in SeqIO.parse(config.REFDB_OUT / "reference_db_20260907.fasta", "fasta")}

    rows, arow = [], []
    for r in idx.itertuples():
        genes = genes_of(r.well)
        cov = contig_coverage(r.assembly_name, r.gbk)
        bulk = None
        if cov:
            c = np.array([v[0] for v in cov.values() if v[1] >= BULK_MIN_CONTIG])
            L = np.array([v[1] for v in cov.values() if v[1] >= BULK_MIN_CONTIG])
            if len(c):
                o = np.argsort(c); cw = np.cumsum(L[o]) / L.sum()
                bulk = float(c[o][np.searchsorted(cw, 0.5)])     # length-weighted median
        genes = assign_organisms(genes, refs) if genes else []
        for g in genes:
            cv = cov.get(g["contig"])
            g["contig_bp"] = cv[1] if cv else None
            g["cov_ratio"] = round(cv[0] / bulk, 3) if (cv and bulk) else None
            rows.append({k: v for k, v in g.items() if k != "seq"} | dict(seq=g["seq"]))
        orgs = sorted({g["organism"] for g in genes})
        a = dict(well=r.well, n_16S_genes=len(genes), n_organisms=len(orgs), bulk_cov=bulk)
        if len(orgs) > 1:
            # the main organism's rRNA sits at or above bulk coverage (multi-copy, collapsed);
            # a contaminant's sits well below it
            score = {o: max((g["cov_ratio"] or 0) for g in genes if g["organism"] == o) for o in orgs}
            ranked = sorted(orgs, key=lambda o: -score[o])
            a.update(main_organism=ranked[0], main_cov_ratio=score[ranked[0]],
                     other_cov_ratio=score[ranked[1]],
                     main_called_by="coverage" if bulk else "undetermined")
            if bulk and not (score[ranked[0]] >= 0.5 and score[ranked[1]] < 0.5 * score[ranked[0]]):
                a["main_called_by"] = "coverage_ambiguous"
        else:
            a.update(main_organism=0 if genes else None, main_called_by="single_organism" if genes else "no_16S")
        arow.append(a)
    g = pd.DataFrame(rows)
    g.drop(columns="seq").to_csv(OUT / "g03_assembly_16S_genes.csv", index=False)
    g.to_pickle(config.GENOME_FNA_CACHE / "g03_assembly_16S_genes.pkl")
    a = pd.DataFrame(arow)
    a.to_csv(OUT / "g03_assembly_purity.csv", index=False)
    print(f"{(a.n_16S_genes > 0).sum()} / {len(a)} assemblies carry a 16S gene >= {MIN_GENE_BP} bp")
    print(f"organisms per assembly: {a.n_organisms.value_counts().sort_index().to_dict()}")
    mixed = a[a.n_organisms > 1]
    print(f"mixed assemblies: {len(mixed)}; main organism called by "
          f"{mixed.main_called_by.value_counts().to_dict()}")
    return g, a


if __name__ == "__main__":
    main()
