"""Stage 4: every strain against every other -- how far apart, and how they differ.

Three distances per pair:
  bp          HAPLOTYPE-AWARE, the one the classes use: the smallest distance between any
              operon haplotype of one strain and any of the other. A read comes from one
              operon, so this is what read assignment has to resolve. A single consensus
              misses it when het sites are linked: G15's 31% haplotype differs from its
              consensus at all 5 sites at once and leaked ~26% of G15's reads to strains
              10-11 bp from the consensus.
  plain_bp    between the majority-base consensuses.
  iupac_bp    between the IUPAC references (each het site matches either allele, sites
              treated independently). Kept for comparison; it ignores linkage.

Classes follow the thresholds set for the community design:

    identical          plain_bp == 0     same majority 16S, base for base
    close              bp < 5            closely related
    indistinguishable  5 <= bp < 15      functionally inseparable by nanopore read assignment
    distinguishable    bp >= 15

These nest (an identical pair is also <5 and <15), so per-strain counts are cumulative.
`d02_read_resolution.py` checks the 5/15 cut-offs against this run's own monoculture reads.

Identity groups are built on plain_bp == 0, which is an equivalence relation. Building
them on the het-aware distance chains: one strain with 23 IUPAC codes is "compatible" with
sequences that differ from each other, and the largest group came out at 31 strains with
two members 23 bp apart.

"Identical 16S" is not "same strain". rRNA operon variation is genome-specific, so for
each identical pair the minor-allele frequencies are compared position by position
(from the reads behind each reference, s04_allele_profiles.pkl):
  different_operons   some site differs by >= OPERON_DIFF in minor-allele frequency, both
                      sides covered by >= OPERON_MIN_COV reads -> different genomes
  same_operons        both carry >= 1 het site (>= HET_MAF) and none differs
  uninformative       neither carries a het site -- 16S has nothing more to say
The collection sheet's `duplicated` column and Strain.Number are reported alongside.

Outputs (02_distinguishability/outputs/):
  d01_distance_matrix_bp.csv         square, haplotype-aware bp (the one to design with)
  d01_distance_matrix_plain_bp.csv   square, majority-base bp
  d01_distance_matrix_iupac_bp.csv   square, IUPAC bp
  d01_pairs.csv                      every pair: all three distances; subs / indel /
                                     homopolymer-indel of the closest haplotype pair
  d01_close_pairs.csv                pairs with bp < 15, with tiers and operon evidence
  d01_identity_groups.csv            groups of byte-identical majority 16S
  d01_strain_distinguishability.csv  per strain: nearest neighbour, counts within 0/5/15 bp
"""
import itertools
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

OUT = config.DIST_OUT
CLOSE_BP, RESOLVE_BP = 5, 15
OPERON_DIFF = 0.20
OPERON_MIN_COV = 30
COLLECTION_CSV = (config.paths.CONSENSUS_MERGE / "plate_format"
                  / "Full_384_strain_collection_no_seq_info.csv")


def klass(bp, plain_bp):
    if plain_bp == 0:
        return "identical"
    if bp < CLOSE_BP:
        return "close"
    if bp < RESOLVE_BP:
        return "indistinguishable"
    return "distinguishable"


def _load(path):
    return {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}


_IU = _PL = _HAP = None


def _init(iu, pl, hap):
    global _IU, _PL, _HAP
    _IU, _PL, _HAP = iu, pl, hap


def _job(pair):
    a, b = pair
    best = None
    for ha in _HAP[a]:
        for hb in _HAP[b]:
            d = dc.core_diff(ha, hb)
            if best is None or d["bp"] < best["bp"]:
                best = d
    return (a, b, best["bp"], dc.core_bp(_PL[a], _PL[b]), dc.core_bp(_IU[a], _IU[b]),
            best["subs"], best["indel_bp"], best["hp_indel_bp"])


def load_haplotypes():
    hap = {}
    for r in SeqIO.parse(config.REFDB_OUT / "reference_db_20260907_haplotypes.fasta", "fasta"):
        hap.setdefault(r.id.rsplit("_h", 1)[0], []).append(str(r.seq))
    return hap


def operon_compare(pa, pb):
    """(verdict, max_freq_diff, position) for two strains with the same majority 16S.

    Compares the full A/C/G/T frequency vector at every position both cover with
    >= OPERON_MIN_COV reads. Read error contributes ~0.5% per base and cancels; a real
    operon difference shows as one base at different frequency in the two strains."""
    if pa["consensus"] == pb["consensus"]:
        idx = jdx = np.arange(len(pa["consensus"]))
    else:
        m = dc.pos_map(pa["consensus"], pb["consensus"])
        idx = np.array(sorted(m)); jdx = np.array([m[i] for i in idx])
    ok = (pa["cov"][idx] >= OPERON_MIN_COV) & (pb["cov"][jdx] >= OPERON_MIN_COV)
    ok[:15] = False; ok[-15:] = False
    diff = np.abs(pa["freqs"][idx] - pb["freqs"][jdx]).max(axis=1) * ok
    k = int(diff.argmax())
    if diff[k] >= OPERON_DIFF:
        return "different_operons", round(float(diff[k]), 3), int(idx[k])
    if ((pa["frac"][idx] * ok) >= dc.HET_MAF).any() or ((pb["frac"][jdx] * ok) >= dc.HET_MAF).any():
        return "same_operons", round(float(diff[k]), 3), None
    return "uninformative", round(float(diff[k]), 3), None


def main(processes=48):
    OUT.mkdir(parents=True, exist_ok=True)
    iu = _load(config.REFDB_OUT / "reference_db_20260907_iupac.fasta")
    pl = _load(config.REFDB_OUT / "reference_db_20260907.fasta")
    ref = pd.read_csv(config.REFDB_OUT / "s01_strain_reference.csv").set_index("strain")
    with open(config.REFDB_OUT / "s04_allele_profiles.pkl", "rb") as fh:
        prof = pickle.load(fh)
    names = sorted(iu)
    pairs = list(itertools.combinations(names, 2))
    print(f"{len(names)} references, {len(pairs):,} pairs")
    hap = load_haplotypes()
    with Pool(processes, initializer=_init, initargs=(iu, pl, hap)) as pool:
        res = pool.map(_job, pairs, chunksize=200)
    df = pd.DataFrame(res, columns=["strain_a", "strain_b", "bp", "plain_bp", "iupac_bp",
                                    "subs", "indel_bp", "hp_indel_bp"])
    df["class"] = [klass(b, p) for b, p in zip(df.bp, df.plain_bp)]
    df["tier_a"] = df.strain_a.map(ref.tier)
    df["tier_b"] = df.strain_b.map(ref.tier)
    df.to_csv(OUT / "d01_pairs.csv", index=False)

    for col, fn in (("bp", "d01_distance_matrix_bp.csv"),
                    ("plain_bp", "d01_distance_matrix_plain_bp.csv"),
                    ("iupac_bp", "d01_distance_matrix_iupac_bp.csv")):
        M = pd.DataFrame(0, index=names, columns=names, dtype=int)
        for a, b, v in zip(df.strain_a, df.strain_b, df[col]):
            M.at[a, b] = M.at[b, a] = v
        M.to_csv(OUT / fn)

    parent = {n: n for n in names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for r in df[df.plain_bp == 0].itertuples():
        parent[find(r.strain_a)] = find(r.strain_b)
    comps = {}
    for n in names:
        comps.setdefault(find(n), []).append(n)
    groups = sorted([sorted(v) for v in comps.values() if len(v) > 1], key=lambda g: (-len(g), g))
    gid = {s: f"G{k + 1:02d}" for k, g in enumerate(groups) for s in g}

    coll = pd.read_csv(COLLECTION_CSV)
    dup_flag = dict(zip(coll["Well.position"], coll["duplicated"]))
    strain_no = dict(zip(coll["Well.position"], coll["Strain.Number"]))

    close = df[df.bp < RESOLVE_BP].copy()
    close["identity_group"] = [gid.get(a) if p == 0 else None
                               for a, p in zip(close.strain_a, close.plain_bp)]
    ops = [operon_compare(prof[a], prof[b]) if p == 0 else (None, None, None)
           for a, b, p in zip(close.strain_a, close.strain_b, close.plain_bp)]
    close["operon_evidence"] = [o[0] for o in ops]
    close["operon_max_freq_diff"] = [o[1] for o in ops]
    close["operon_diff_position"] = [o[2] for o in ops]
    close["collection_dup_a"] = close.strain_a.map(dup_flag)
    close["collection_dup_b"] = close.strain_b.map(dup_flag)
    close = close.sort_values(["bp", "plain_bp", "strain_a", "strain_b"])
    close.to_csv(OUT / "d01_close_pairs.csv", index=False)

    grow = []
    ident = close[close.plain_bp == 0]
    for k, g in enumerate(groups):
        G = f"G{k + 1:02d}"
        ev = ident[ident.identity_group == G].operon_evidence.value_counts()
        # members whose operon pattern matches every other member's
        grow.append(dict(group=G, size=len(g), members=" ".join(g),
                         tiers="".join(sorted(ref.tier.reindex(g).fillna("-"))),
                         pairs=len(g) * (len(g) - 1) // 2,
                         same_operons=int(ev.get("same_operons", 0)),
                         different_operons=int(ev.get("different_operons", 0)),
                         uninformative=int(ev.get("uninformative", 0)),
                         collection_dup_flagged=" ".join(s for s in g
                                                         if str(dup_flag.get(s, "no")) != "no")))
    gdf = pd.DataFrame(grow)
    gdf.to_csv(OUT / "d01_identity_groups.csv", index=False)

    rows = []
    full = pd.concat([df.rename(columns={"strain_a": "s", "strain_b": "o"}),
                      df.rename(columns={"strain_b": "s", "strain_a": "o"})])
    for s, g in full.groupby("s"):
        nn = g.loc[g.bp.idxmin()]
        rows.append(dict(strain=s, tier=ref.tier.get(s), status=ref.status.get(s),
                         no_growth=ref.no_growth.get(s), identity_group=gid.get(s),
                         nearest=nn.o, nearest_bp=int(nn.bp), nearest_plain_bp=int(nn.plain_bp),
                         nearest_tier=ref.tier.get(nn.o),
                         n_identical=int((g.plain_bp == 0).sum()),
                         n_within_5=int((g.bp < CLOSE_BP).sum()),
                         n_within_15=int((g.bp < RESOLVE_BP).sum()),
                         within_5=" ".join(sorted(g[g.bp < CLOSE_BP].o)),
                         het_sites=ref.het_sites.get(s),
                         n_haplotypes=ref.n_haplotypes.get(s),
                         many_het=ref.many_het.get(s), mixed_source=ref.mixed_source.get(s),
                         competing=ref.competing.get(s),
                         collection_strain_number=strain_no.get(s),
                         collection_duplicated=dup_flag.get(s)))
    st = pd.DataFrame(rows)
    st.to_csv(OUT / "d01_strain_distinguishability.csv", index=False)
    report(df, st, gdf, close)
    return df, st


def report(df, st, groups, close):
    order = ["identical", "close", "indistinguishable", "distinguishable"]
    print("\npairs by class:", df["class"].value_counts().reindex(order).to_dict())
    print("strains with at least one other strain that is ...")
    for col, lab in (("n_identical", "identical 16S"), ("n_within_5", "< 5 bp"),
                     ("n_within_15", "< 15 bp")):
        print(f"  {lab:14s} {(st[col] > 0).sum():4d} / {len(st)}")
    print(f"  nearest neighbour >= 15 bp: {(st.nearest_bp >= RESOLVE_BP).sum()}")
    print(f"\n{len(groups)} identity groups covering {groups['size'].sum()} strains; "
          f"sizes {groups['size'].tolist()[:12]} ...")
    ev = close[close.plain_bp == 0].operon_evidence.value_counts().to_dict()
    print(f"  identical pairs by operon evidence: {ev}")
    c = df[(df.bp > 0) & (df.bp < RESOLVE_BP)]
    print(f"\ndifferences in pairs 1-14 bp (n={len(c)}): substitutions {c.subs.sum()}, "
          f"indel bp {c.indel_bp.sum()} (homopolymer {c.hp_indel_bp.sum()})")


if __name__ == "__main__":
    main()
