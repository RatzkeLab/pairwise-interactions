"""Load every available 16S sequence for every strain into one long table.

Four independent sources, each with a different failure mode -- which is the point: agreement
between two of them is only meaningful because they cannot fail the same way.

  genome_16S       extracted from the assembled genomes (Illumina/NGS assemblies). External
                   ground truth: nothing about the ONT run can corrupt it. Multi-copy -- 587
                   records over 375 assemblies, up to 7 rRNA operons each, and 290 records are
                   <1200 bp (fragmented contig ends), so length is not comparable to an
                   amplicon and every comparison must be identity-based, not raw bp.
  mono_consensus   from this experiment's MONO wells: one strain shot into a well alone. The
                   cleanest in-experiment evidence of what a source well actually contained,
                   and the only source that is immune to the two-strain deconvolution problem.
  pair_consensus   from this experiment's PAIR wells, as used by the interaction pipeline.
                   Self-consistent by construction and therefore the one source that can be
                   confidently wrong.
  corroborated_db  external 16S db built from prior corroborated experiments.

Everything is keyed by the experiment's own well-coordinate strain label. The genome source
reaches that namespace through mapping_384_well_plate_collection.csv, which is exactly the
join under suspicion -- so genome_16S rows carry the assembly they came from, letting a
mismatch be attributed to either the join or the biology.
"""

import re
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

import qc_config as C


def read_fasta(path, split_header=True):
    """-> list of (name, seq). Keeps duplicates: multi-copy 16S is real signal, not a bug."""
    out, name = [], None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:]
            if split_header:
                name = re.split(r"[\t ]", name)[0]
            out.append([name, ""])
        elif name is not None:
            out[-1][1] += line.strip()
    return [(n, s.upper()) for n, s in out if s]


def norm_assembly(x):
    """'GCF_002966835.1_ASM296683v1_genomic' -> 'GCF_002966835.1'.

    The 16S fasta and the mapping table label the same assembly differently; normalising to the
    accession lifts the join from 233/375 to 294/298 mapping rows.
    """
    m = re.match(r"^(GC[AF]_\d+\.\d+)", str(x))
    return m.group(1) if m else str(x)


def load_mapping():
    m = pd.read_csv(C.MAPPING_CSV)
    m["strain"] = m["strain"].astype(str)
    m["asm_norm"] = m["assembly_name"].astype(str).map(norm_assembly)
    return m


def genome_16S_by_well():
    """-> DataFrame(strain_label, source, copy, seq, assembly, genomic_strain_name)."""
    m = load_mapping()
    by_asm = {}
    for name, seq in read_fasta(C.GENOME_16S_FASTA):
        by_asm.setdefault(norm_assembly(name), []).append(seq)
    rows = []
    for r in m.itertuples():
        for i, seq in enumerate(by_asm.get(r.asm_norm, [])):
            rows.append({"strain_label": r.Well_souce_plate, "source": "genome_16S",
                         "copy": i, "seq": seq, "assembly": r.asm_norm,
                         "genomic_strain_name": r.strain})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# mono-well consensus -- built here for whichever experiment lacks one
# ---------------------------------------------------------------------------

def _sample_id(plate, well):
    return f"Plate{int(plate):02d}_{well}"


def mono_well_reads(exp, min_reads=C.MIN_READS):
    """-> {strain: [fastq paths]} for wells the layout says contain exactly one strain."""
    lay = pd.read_csv(exp.layout_csv)
    mono = lay[lay["well_type"] == "mono"]
    out = {}
    for r in mono.itertuples():
        sid = _sample_id(r.dest_plate, r.dest_well)
        hits = list(exp.demux_dir.glob(f"{sid}.fastq*"))
        if not hits:
            continue
        strain = r.strain1 if pd.notna(r.strain1) else r.strain2
        out.setdefault(str(strain), []).append(hits[0])
    return out


def build_mono_consensus(exp, min_reads=C.MIN_READS, force=False):
    """Per-strain consensus from mono wells, via minimap2+racon-free majority pileup.

    Uses the simplest thing that is defensible: map the well's reads to their own longest read
    and take the samtools consensus. Mono wells are single-strain by design, so there is no
    deconvolution to get wrong -- which is the whole reason this source is worth building.
    """
    out_fa = exp.out / f"mono_consensus_{exp.name}.fasta"
    if out_fa.exists() and not force:
        return out_fa
    reads = mono_well_reads(exp, min_reads)
    written = 0
    with open(out_fa, "w") as fh, tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for strain, fqs in sorted(reads.items()):
            recs = []
            for fq in fqs:
                opener = subprocess.run(["zcat", "-f", str(fq)], capture_output=True, text=True).stdout
                lines = opener.splitlines()
                recs += [lines[i + 1] for i in range(0, len(lines) - 3, 4)]
            recs = [r for r in recs if 1000 <= len(r) <= 1800]
            if len(recs) <= min_reads:
                continue
            backbone = sorted(recs, key=len)[len(recs) // 2]
            (td / "b.fa").write_text(f">b\n{backbone}\n")
            (td / "r.fq").write_text("".join(f"@r{i}\n{s}\n+\n{'I'*len(s)}\n" for i, s in enumerate(recs)))
            sam = subprocess.run(["minimap2", "-ax", "map-ont", "-t", "4", str(td / "b.fa"), str(td / "r.fq")],
                                 capture_output=True, text=True).stdout
            (td / "a.sam").write_text(sam)
            subprocess.run(f"samtools sort -o {td}/a.bam {td}/a.sam && samtools index {td}/a.bam",
                           shell=True, capture_output=True)
            cons = subprocess.run(["samtools", "consensus", "-a", "-f", "fasta", str(td / "a.bam")],
                                  capture_output=True, text=True).stdout
            seq = "".join(l.strip() for l in cons.splitlines() if not l.startswith(">"))
            # trim uncovered ends only -- deleting INTERNAL Ns would splice the sequence and
            # manufacture false deletions, which is exactly the artifact this QC hunts for
            seq = seq.strip("N")
            # internal Ns are KEPT and recorded: every comparison in qc_compare treats N as a
            # wildcard, so a low-coverage well stays usable instead of being silently dropped
            n_internal = seq.count("N")
            if len(seq) >= 1200 and n_internal <= 0.10 * len(seq):
                fh.write(f">{strain}|mono|n={len(recs)}|internalN={n_internal}\n{seq}\n")
                written += 1
    print(f"[{exp.name}] mono consensus: {written} strains -> {out_fa.name}")
    return out_fa


def strain_label_from_header(h):
    """Header -> bare well-coordinate label.

    The three fasta sources decorate the label differently ('A6_consensus', 'C8|Plate01_F15|
    reads=52', plain 'C8'); normalising here keeps every downstream join on one namespace.
    """
    lab = str(h).split("|")[0].split()[0]
    return lab[: -len("_consensus")] if lab.endswith("_consensus") else lab


def load_all_sources(exp):
    """Every 16S sequence available for this experiment's strain labels, one row per sequence."""
    frames = [genome_16S_by_well()]

    for src, path in [("pair_consensus", exp.pair_consensus),
                      ("corroborated_db", C.CORROBORATED_FASTA)]:
        frames.append(pd.DataFrame([
            {"strain_label": strain_label_from_header(n), "source": src, "copy": 0, "seq": s}
            for n, s in read_fasta(path, split_header=False)]))

    mono = exp.mono_consensus or build_mono_consensus(exp)
    if Path(mono).exists():
        frames.append(pd.DataFrame([
            {"strain_label": strain_label_from_header(n), "source": "mono_consensus", "copy": 0, "seq": s}
            for n, s in read_fasta(Path(mono), split_header=False)]))

    df = pd.concat(frames, ignore_index=True)
    df["experiment"] = exp.name
    df["length"] = df["seq"].str.len()
    return df
