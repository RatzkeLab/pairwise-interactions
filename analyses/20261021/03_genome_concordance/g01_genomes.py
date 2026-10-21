"""Stage 7: the genomes on file, compared with each other and with the wells' 16S.

The genomic tables (KO etc.) were built from assemblies of the ORIGINAL isolates. The
community experiment will be inoculated from the 20260907 working plate, which has
drifted -- a few invasive organisms have overgrown many wells. A well is only usable for
genome-based prediction if the organism in it is the one whose genome we have. Two
independent measurements answer that:

  genome vs genome   skani ANI over the 298 assemblies. >=99.9% is one organism (clonal
                     duplicates in the collection); <95% different species.
  genome vs well     every 20260907 reference (consensus + operon haplotypes) blasted
                     against every assembly: does a well's OWN genome carry the 16S the
                     working plate showed, and which genomes do?

Outputs (03_genome_concordance/outputs/):
  g01_genome_index.csv           well -> genome-table strain id, assembly, size
  g01_ani_pairs.csv              every genome pair skani reports (ANI >~80%)
  g01_blast_16S_vs_genomes.csv.gz  every HSP (>=85% identity) of every query in every genome
"""
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import config  # noqa: E402

import pandas as pd
from Bio import SeqIO

OUT = config.GENOME_OUT
CACHE = config.GENOME_FNA_CACHE
FNA = CACHE / "fna"          # one file per genome -- skani globs this, so nothing else goes here
BLAST = CACHE / "blast"


def genome_index():
    m = pd.read_csv(config.GENOME_MAPPING)
    m["strain"] = m["strain"].astype(str)
    gbk = [f[:-4] for f in __import__("os").listdir(config.GENOME_GBK_DIR) if f.endswith(".gbk")]

    def resolve(a):
        hit = [g for g in gbk if g == a or g.startswith(a + "_")]
        return hit[0] if len(hit) == 1 else None
    m["gbk"] = m.assembly_name.map(resolve)
    assert m.gbk.notna().all() and m.gbk.is_unique, "every mapped assembly must resolve to one GenBank file"
    return m.rename(columns={"Well_souce_plate": "well", "strain": "genome_table_id"})[
        ["well", "genome_table_id", "sequencing_batch", "assembly_name", "gbk"]]


def _extract(args):
    well, gbk = args
    out = FNA / f"{well}.fna"
    n = L = 0
    with open(out, "w") as fh:
        for rec in SeqIO.parse(config.GENOME_GBK_DIR / f"{gbk}.gbk", "genbank"):
            s = str(rec.seq)
            if s and set(s) != {"N"}:
                fh.write(f">{well}|{rec.id}\n{s}\n"); n += 1; L += len(s)
    return well, n, L


def main(processes=48):
    for d in (OUT, FNA, BLAST):
        d.mkdir(parents=True, exist_ok=True)
    idx = genome_index()
    with Pool(32) as pool:
        sizes = pool.map(_extract, list(zip(idx.well, idx.gbk)))
    idx = idx.merge(pd.DataFrame(sizes, columns=["well", "contigs", "genome_bp"]), on="well")
    idx.to_csv(OUT / "g01_genome_index.csv", index=False)
    print(f"{len(idx)} genomes extracted, median {idx.genome_bp.median() / 1e6:.2f} Mb, "
          f"{idx.contigs.median():.0f} contigs")

    sk = BLAST / "skani_sparse.tsv"
    subprocess.run([config.SKANI, "triangle", "-t", str(processes), "--sparse",
                    *sorted(str(p) for p in FNA.glob("*.fna")), "-o", str(sk)],
                   check=True, capture_output=True)
    k = pd.read_csv(sk, sep="\t")
    k["well_a"] = k.Ref_file.map(lambda p: Path(p).stem)
    k["well_b"] = k.Query_file.map(lambda p: Path(p).stem)
    k = k.rename(columns={"ANI": "ani", "Align_fraction_ref": "af_a", "Align_fraction_query": "af_b"})
    k[["well_a", "well_b", "ani", "af_a", "af_b"]].to_csv(OUT / "g01_ani_pairs.csv", index=False)
    print(f"skani: {len(k)} genome pairs above the ~80% ANI screen, "
          f"{(k.ani >= 99.9).sum()} at >=99.9%")

    allfa = BLAST / "all_genomes.fna"
    with open(allfa, "w") as out:
        for p in sorted(FNA.glob("*.fna")):
            out.write(p.read_text())
    subprocess.run([str(config.BLAST_BIN / "makeblastdb"), "-in", str(allfa), "-dbtype", "nucl",
                    "-out", str(BLAST / "all_genomes")], check=True, capture_output=True)
    q = BLAST / "queries.fa"
    with open(q, "w") as fh:
        for fa in ("reference_db_20260907.fasta", "reference_db_20260907_haplotypes.fasta"):
            for r in SeqIO.parse(config.REFDB_OUT / fa, "fasta"):
                fh.write(f">{r.id}\n{r.seq}\n")
    cols = "qseqid sseqid pident length qlen mismatch gaps qstart qend bitscore"
    raw = BLAST / "blast.tsv"
    subprocess.run([str(config.BLAST_BIN / "blastn"), "-query", str(q), "-db", str(BLAST / "all_genomes"),
                    "-outfmt", f"6 {cols}", "-perc_identity", "85", "-max_target_seqs", "100000",
                    "-max_hsps", "20", "-num_threads", str(processes), "-out", str(raw)], check=True)
    b = pd.read_csv(raw, sep="\t", names=cols.split())
    b["genome"] = b.sseqid.str.split("|").str[0]
    b["diffs"] = b.mismatch + b.gaps
    b[["qseqid", "genome", "sseqid", "pident", "length", "qlen", "diffs", "qstart", "qend",
       "bitscore"]].to_csv(OUT / "g01_blast_16S_vs_genomes.csv.gz", index=False)
    print(f"blast: {b.qseqid.nunique()} queries x {b.genome.nunique()} genomes, "
          f"{len(b):,} HSPs")

if __name__ == "__main__":
    main()
