"""Base-exact 16S consensus per population, and bp-level distances between consensuses.

Why this exists next to strain_consensus.py: that module answers "which organism is in
this well?" and does it well, but its per-well representative is a single medoid READ and
the per-strain consensus is a column vote over ~5-10 such reads. At Q24 a read sits a
median ~12 bp from its own strain's consensus, so a vote over a handful of reads is
accurate to a few bp -- fine for identity at the ~100 bp level, too noisy to classify
strain pairs as identical / <5 bp / <15 bp, which is what community design needs.

Here every population gets a spoa consensus over its highest-quality reads. Checked on
P15's two monocultures (different plates): spoa over the top-20, top-50 and top-100 reads
by qscore gives the SAME sequence in both wells, 0 bp apart. Medaka polishing moved one of
them by 1 bp, i.e. it made the result less reproducible, so it is not used. A random
50-read subset (low-qscore reads included) drifted 6 bp, which is why reads are ranked by
qscore before the consensus.
"""
import collections
import os
import re
import subprocess
import tempfile

import edlib
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

SPOA = "/home/rl/mambaforge/envs/NGSpeciesID/bin/spoa"

LEN_LO, LEN_HI = 1300, 1600   # full-length amplicon, primers already trimmed by minibar
READ_CAP = 800                # reads used per well (random subsample above this)
CLUSTER_THR = 0.08            # greedy seed clustering, same value 20260921 learned from data
MIN_POP_READS = 8             # a population needs this many reads ...
MIN_POP_FRAC = 0.03           # ... and this share of the well to get a consensus
TOP_Q_READS = 50              # spoa input: the highest-qscore reads of the population
OUTLIER_BP = 60               # a read further than this from every population consensus
                              # is chimeric or from an organism below MIN_POP_*
HET_MIN_READS = 30            # site scan needs this depth
HET_SCAN_READS = 200
HET_MAF = 0.20                # minor substitution allele at >= this frequency = het site


def revcomp(s):
    return str(Seq(s).reverse_complement())


def load_reads_q(path):
    """[(seq, qscore)] from a dorado fastq; qscore parsed from the `qs:f:` tag."""
    out = []
    if path is None or not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    for rec in SeqIO.parse(str(path), "fastq"):
        m = re.search(r"qs:f:([0-9.]+)", rec.description)
        q = float(m.group(1)) if m else float(np.mean(rec.letter_annotations["phred_quality"]))
        out.append((str(rec.seq).upper(), q))
    return out


def orient(seq, anchor):
    """Reads arrive in both orientations. Any bacterial 16S aligns to any other at
    >=75% identity forward and ~50% reverse-complemented, so one anchor serves all."""
    f = edlib.align(seq, anchor, mode="NW", task="distance")["editDistance"]
    r = edlib.align(revcomp(seq), anchor, mode="NW", task="distance")["editDistance"]
    return seq if f <= r else revcomp(seq)


def prepare_reads(path, anchor, seed=0, cap=READ_CAP):
    reads = [(s, q) for s, q in load_reads_q(path) if LEN_LO <= len(s) <= LEN_HI]
    if len(reads) > cap:
        rng = np.random.default_rng(seed)
        reads = [reads[i] for i in rng.choice(len(reads), cap, replace=False)]
    return [(orient(s, anchor), q) for s, q in reads]


def spoa(seqs):
    if len(seqs) == 1:
        return seqs[0]
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for i, s in enumerate(seqs):
            fh.write(f">r{i}\n{s}\n")
        fa = fh.name
    try:
        out = subprocess.run([SPOA, fa, "-l", "1", "-r", "0"], capture_output=True,
                             text=True, check=True).stdout
    finally:
        os.unlink(fa)
    return out.strip().split("\n")[-1]


def top_q_consensus(reads, n=TOP_Q_READS):
    """reads: [(seq, q)]. spoa over the n highest-qscore reads."""
    best = sorted(reads, key=lambda x: -x[1])[:n]
    return spoa([s for s, _ in best])


def _ed(a, b, k=-1):
    return edlib.align(a, b, mode="NW", task="distance", k=k)["editDistance"]


def greedy_clusters(seqs, thr=CLUSTER_THR):
    """Each read joins its nearest seed if within thr, else seeds a new cluster. Seeds,
    not neighbours, decide membership, so intermediate reads cannot chain two strains
    together (see well_composition.py for the measurement)."""
    seeds, members = [], []
    for i, s in enumerate(seqs):
        best, bj = None, -1
        for j, cs in enumerate(seeds):
            k = int(thr * max(len(s), len(cs)))
            d = _ed(s, cs, k)
            if d >= 0 and (best is None or d < best):
                best, bj = d, j
        if best is not None:
            members[bj].append(i)
        else:
            seeds.append(s); members.append([i])
    return sorted(members, key=len, reverse=True)


def core_diff(a, b):
    """Differences between two consensus sequences, ignoring end jitter.

    The interior of each sequence (END_TRIM bp dropped at both ends) is aligned as an
    infix of the other, in both directions, and the SMALLER count is kept: an internal
    difference shows up in both directions, while one sequence overhanging the other only
    shows up in one. (Taking the larger reported a 19-23 bp floor against the older DBs,
    whose sequences still carry primer tails at both ends.) Consensus ends
    wobble by 1-2 bp with read-end trimming, and ~2% of wells keep a few bases of 27F
    primer tail (`...TGGCTCAG`); a global alignment interleaves those with matches and
    reports up to 9 bp of "difference" between two monocultures of the same strain.
    The 20 bp dropped are conserved flanks next to the primers, where reads are noisiest.

    IUPAC codes (from `iupac`) match any base they stand for, so a site where a strain's
    own operons disagree is not counted as a difference -- reads from that strain carry
    either base, so the site cannot tell it apart from a strain carrying one of them.

    Every difference is classified, because they are not equally visible to a nanopore:
    a substitution is read correctly ~99.6% of the time at Q24, while a homopolymer-length
    change is exactly the error ONT makes most.

    Returns dict(bp, subs, indel_bp, hp_indel_bp).
    """
    best = None
    for q, t in ((a, b), (b, a)):
        d = _core_diff_one(q, t)
        if best is None or d["bp"] < best["bp"]:
            best = d
    return best


END_TRIM = 20
_IUPAC = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
          "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
_CODE = {frozenset(v): k for k, v in _IUPAC.items()}
IUPAC_EQ = ([(k, b) for k, v in _IUPAC.items() for b in v]
            + [(k1, k2) for k1, v1 in _IUPAC.items() for k2, v2 in _IUPAC.items()
               if k1 != k2 and set(v1) & set(v2)])


def _core_diff_one(a, b):
    core = a[END_TRIM:len(a) - END_TRIM]
    r = edlib.align(core, b, mode="HW", task="path", additionalEqualities=IUPAC_EQ)
    s, e = r["locations"][0]
    return _classify(core, b[s:e + 1], r["cigar"])


def _classify(a, b, cigar):
    ops = [(int(n), op) for n, op in re.findall(r"(\d+)([=XID])", cigar)]
    ia = ib = 0
    subs = indel = hp = 0
    for n, op in ops:
        if op == "=":
            ia += n; ib += n
        elif op == "X":
            subs += n; ia += n; ib += n
        elif op == "I":                  # bases present in a, absent in b
            run = a[ia:ia + n]
            indel += n
            if _is_hp(run, a, ia, n):
                hp += n
            ia += n
        elif op == "D":                  # bases present in b, absent in a
            run = b[ib:ib + n]
            indel += n
            if _is_hp(run, b, ib, n):
                hp += n
            ib += n
    return dict(bp=subs + indel, subs=subs, indel_bp=indel, hp_indel_bp=hp)


def _is_hp(run, seq, pos, n):
    """An indel run is a homopolymer-length change if it is one base repeated and that
    base also flanks it -- i.e. it lengthens or shortens an existing run."""
    if len(set(run)) != 1:
        return False
    base = run[0]
    left = seq[pos - 1] if pos > 0 else ""
    right = seq[pos + n] if pos + n < len(seq) else ""
    return base == left or base == right


def core_bp(a, b):
    return core_diff(a, b)["bp"]


def allele_profile(consensus, seqs, full=False):
    """Per consensus position: (most frequent non-consensus SUBSTITUTION base, its
    frequency, coverage). Indels are excluded because at this depth they are dominated
    by read error; the first and last 15 positions are zeroed (read ends are noisy).

    full=True also returns the L x 4 (A,C,G,T) frequency matrix. Comparing two strains
    needs it: a tri-allelic site reports a different "top minor base" in two wells of
    the same isolate (MYb181, position 180: T at 0.20 vs A at 0.28) though both carry
    both."""
    L = len(consensus)
    counts = np.zeros((L, 4), dtype=np.int32)
    cov = np.zeros(L, dtype=np.int32)
    for s in seqs:
        r = edlib.align(s, consensus, mode="NW", task="path")
        na = edlib.getNiceAlignment(r, s, consensus)
        pos = 0
        for qc, tc in zip(na["query_aligned"], na["target_aligned"]):
            if tc == "-":
                continue
            if qc in _IDX:
                counts[pos, _IDX[qc]] += 1
                cov[pos] += 1
            pos += 1
    cons_idx = np.array([_IDX.get(c, -1) for c in consensus])
    alt = counts.copy()
    alt[np.arange(L)[cons_idx >= 0], cons_idx[cons_idx >= 0]] = 0
    frac = alt.max(axis=1) / np.maximum(cov, 1)
    frac[:15] = 0; frac[-15:] = 0
    out = (np.array(list("ACGT"))[alt.argmax(axis=1)], frac, cov)
    if full:
        return out + ((counts / np.maximum(cov, 1)[:, None]).astype(np.float32),)
    return out


def het_sites(consensus, seqs, maf=HET_MAF):
    """Positions where a minority SUBSTITUTION is carried by >= maf of reads.

    Within-population heterogeneity is either intragenomic (rRNA operons that differ --
    then it recurs at a stable frequency in every well of the strain) or a second
    organism too close to cluster apart (then it varies between wells). Either way it
    is the thing that makes a single consensus an incomplete description of a strain.

    Returns (n_sites, [(position, minor_frac, minor_base)]), or (None, None) if the
    population is too shallow to say.
    """
    if len(seqs) < HET_MIN_READS:
        return None, None
    seqs = seqs[:HET_SCAN_READS]
    alt, frac, cov = allele_profile(consensus, seqs)
    sites = np.where((cov >= 0.8 * len(seqs)) & (frac >= maf))[0]
    return len(sites), [(int(p), round(float(frac[p]), 3), str(alt[p])) for p in sites]


_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def iupac(consensus, sites):
    """The consensus with each het site replaced by the IUPAC code for both alleles."""
    s = list(consensus)
    for pos, _, alt in sites or []:
        if s[pos] in _IDX and alt != s[pos]:
            s[pos] = _CODE[frozenset((s[pos], alt))]
    return "".join(s)


HAP_MIN_FRAC = 0.10


def _bases_at(seq, consensus, positions):
    """The read's base at each consensus position (None where it has a gap)."""
    r = edlib.align(seq, consensus, mode="NW", task="path")
    na = edlib.getNiceAlignment(r, seq, consensus)
    want, out, pos = set(positions), {}, 0
    for qc, tc in zip(na["query_aligned"], na["target_aligned"]):
        if tc == "-":
            continue
        if pos in want:
            out[pos] = qc if qc != "-" else None
        pos += 1
    return [out.get(p) for p in positions]


def phase_haplotypes(consensus, het, seqs, min_frac=HAP_MIN_FRAC):
    """Split a strain's reads into operon haplotypes at its het sites.

    A single consensus is the per-column majority, and it misrepresents a strain whose
    het sites are LINKED: with two operon types at 60/40, the 40% type differs from the
    consensus at every het site at once, and can sit closer to another strain than the
    consensus does. Measured: G15 and O9 (5 het sites, clean monos) send ~26% of their
    reads to strains 10-11 bp from their consensus.

    Each read is read out at the het positions as major (0) / minor (1); patterns
    carried by >= min_frac of fully-called reads become haplotypes, every read is then
    given to its nearest haplotype (Hamming over the sites it covers), and each
    haplotype's sequence is the consensus with its minor alleles substituted.

    Returns [(sequence, frac, pattern)], most frequent first. No het sites -> the
    consensus alone.
    """
    if not het:
        return [(consensus, 1.0, "")]
    positions = [p for p, _, _ in het]
    minor = [m for _, _, m in het]
    major = [consensus[p] for p in positions]
    calls = []
    for s in seqs:
        b = _bases_at(s, consensus, positions)
        calls.append([0 if x == M else 1 if x == m else None
                      for x, M, m in zip(b, major, minor)])
    full = collections.Counter(tuple(c) for c in calls if None not in c)
    n_full = sum(full.values())
    seeds = [p for p, c in full.most_common() if c >= min_frac * n_full]
    if not seeds:
        return [(consensus, 1.0, "0" * len(positions))]
    tally = collections.Counter()
    for c in calls:
        d = [sum(1 for x, y in zip(c, sd) if x is not None and x != y) for sd in seeds]
        if sum(x is not None for x in c) >= len(c) / 2:
            tally[int(np.argmin(d))] += 1
    tot = sum(tally.values()) or 1
    out = []
    for k, sd in enumerate(seeds):
        s = list(consensus)
        for p, bit, m in zip(positions, sd, minor):
            if bit == 1:
                s[p] = m
        out.append(("".join(s), tally[k] / tot, "".join(map(str, sd))))
    return sorted(out, key=lambda h: -h[1])


def pos_map(a, b):
    """position in a -> position in b, over aligned (match/mismatch) columns."""
    r = edlib.align(a, b, mode="NW", task="path")
    m, ia, ib = {}, 0, 0
    for n, op in re.findall(r"(\d+)([=XID])", r["cigar"]):
        n = int(n)
        if op in "=X":
            for k in range(n):
                m[ia + k] = ib + k
            ia += n; ib += n
        elif op == "I":
            ia += n
        else:
            ib += n
    return m


def well_populations(reads):
    """reads: [(seq, q)] oriented. Returns list of population dicts, largest first.

    1. greedy clusters at CLUSTER_THR to find candidate populations,
    2. consensus per cluster over its top-qscore reads,
    3. one reassignment of every read to its nearest consensus (fixes reads the
       greedy pass put with an early noisy seed), then consensus again.
    """
    if not reads:
        return [], 0
    seqs = [s for s, _ in reads]
    comps = greedy_clusters(seqs)
    n = len(seqs)
    cut = max(MIN_POP_READS, MIN_POP_FRAC * n)
    big = [c for c in comps if len(c) >= cut]
    if not big:
        return [], n
    cons = [top_q_consensus([reads[i] for i in c]) for c in big]

    assign = collections.defaultdict(list)
    n_out = 0
    for i, s in enumerate(seqs):
        d = [_ed(s, c, OUTLIER_BP) for c in cons]
        d = [x if x >= 0 else 10 ** 9 for x in d]
        j = int(np.argmin(d))
        if d[j] > OUTLIER_BP:
            n_out += 1
        else:
            assign[j].append(i)
    pops = []
    for j in sorted(assign, key=lambda j: -len(assign[j])):
        mem = assign[j]
        if len(mem) < MIN_POP_READS:
            continue
        sub = [reads[i] for i in mem]
        c = top_q_consensus(sub)
        pops.append(dict(consensus=c, n_reads=len(mem), frac=len(mem) / n,
                         reads=sub))
    # two clusters can converge on the same consensus after reassignment; merge them
    merged = []
    for p in pops:
        for m in merged:
            if core_bp(p["consensus"], m["consensus"]) <= 1:
                m["n_reads"] += p["n_reads"]; m["frac"] += p["frac"]
                m["reads"] += p["reads"]
                break
        else:
            merged.append(p)
    merged.sort(key=lambda p: -p["n_reads"])
    return merged, n_out
