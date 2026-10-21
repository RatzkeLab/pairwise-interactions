"""Stage 3: one reference per strain, chosen by how many independent wells agree on it.

Every exact consensus a strain produced -- its monocultures' populations and its credited
pair-well groups -- is an OBSERVATION. Observations are clustered at "the same sequence":

    plain distance <= SAME_BP, or
    distance <= SAME_BP once sites that are het in BOTH observations are IUPAC-coded.

Requiring the het site in both is what keeps a mixed well honest. A well holding two
organisms 20 bp apart has ~20 het sites, and if its own het sites were honoured it would
match either organism at 0 bp (A17's two monos: 21 bp plain, 2 bp het-aware). Intragenomic
operon variation recurs in every well of the strain; a mixture does not.

Calibration for SAME_BP (230 strains with two usable monos, different plates): 191 are
0 bp apart het-aware, 199 within 2 bp, then a gap to 9+ bp. The consensus adds ~no noise.

The winning cluster is the one seen in the most distinct wells; monocultures break ties.
Tiers, all about independent evidence for THIS run's source plate:

  A  monos_agree          both monocultures, on different plates, give this sequence
  B  mono_plus_pairs      one monoculture + >=2 pair wells
  C  pairs_only           no monoculture, but >=3 pair wells on >=2 plates
     single_mono_ext      one monoculture only, but an older, independent run's call
                          for the same well agrees within EXT_AGREE_BP
  D  single_mono          one monoculture only, nothing independent agrees or disagrees
     monos_disagree       the monos give different sequences and nothing breaks the tie
     few_pair_wells       no monoculture, 1-2 pair wells
     weak                 anything else

Why so many strains have no pair-well support: a strong loser is simply not there. E1 is
0% of reads in all six of its pair wells, every partner taking 95-98%. That is a result,
not a failure -- but it leaves the monocultures as the only evidence.

When the two leading clusters tie (typically two monos that disagree and no pair well to
break it), the older full_collection and contamination_informed DBs are consulted: they
come from earlier runs on other plate copies, so they are independent of this run's wells.
Their label support is only ~52% for this run -- 16S degeneracy plus a known lineage
scramble -- so agreement counts FOR a sequence but disagreement is never held against one.

and flags that cut across tiers:
  competing    a second sequence, >COMPETING_MIN_BP from the winner, seen in >=2 wells or
               as a mono's dominant population (closer than that it is the same organism
               split by het-site coin flips)
  mixed_source both monos hold a second population at >=20% -- the source well is a mix
  many_het     >10 consistent het sites: unusual for operon variation, could be a stable
               mixture of two close organisms
  no_growth    never grew in 20260907 (both monos on different plates failed)

Outputs:
  s01_strain_reference.csv          per strain: tier, evidence, flags
  s02_strain_observations.csv       every observation and the cluster it joined (audit)
  s03_strain_het_sites.csv          consistent intragenomic het sites per strain
  s05_haplotypes.csv                per strain: operon haplotypes phased at its het sites
  reference_db_20260907_haplotypes.fasta   one entry per haplotype (>= 10% of reads),
                                    named <strain>_h<k>. Use THIS for read assignment:
                                    a read goes to its nearest haplotype, then to the strain
  s04_allele_profiles.pkl           per strain: minor-allele frequency at every position,
                                    from the reads behind the final sequence -- lets two
                                    strains with the same 16S be compared operon by operon
  reference_db_20260907.fasta       plain-base reference, all tiers (header carries tier)
  reference_db_20260907_iupac.fasta same, het sites IUPAC-coded -- for distances only
"""
import collections
import pickle
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import config  # noqa: E402

import numpy as np
import pandas as pd

import deep_consensus as dc

OUT = config.REFDB_OUT
SAME_BP = 2
COMPETING_MIN_WELLS = 2
COMPETING_MIN_BP = 5
MIXED_SOURCE_FRAC = 0.20
MANY_HET = 10
MONO_MINOR_FRAC = 0.10
EXT_AGREE_BP = 10        # older DBs were built from lower-depth, older-chemistry reads


def agree(x, y):
    """(same?, plain_bp, shared_het_bp)"""
    plain = dc.core_bp(x["consensus"], y["consensus"])
    if plain <= SAME_BP:
        return True, plain, plain
    if not x["het"] or not y["het"]:
        return False, plain, plain
    m = dc.pos_map(x["consensus"], y["consensus"])
    yh = {p: alt for p, _, alt in y["het"]}
    xs, ys = [], []
    for p, f, alt in x["het"]:
        q = m.get(p)
        if q is not None and q in yh:
            xs.append((p, f, alt)); ys.append((q, f, yh[q]))
    if not xs:
        return False, plain, plain
    d = dc.core_bp(dc.iupac(x["consensus"], xs), dc.iupac(y["consensus"], ys))
    return d <= SAME_BP, plain, d


def load_external():
    from io_utils import load_reference_db
    return {"full_collection": load_reference_db(config.FULL_COLLECTION_FASTA),
            "contam_informed": load_reference_db(config.CONTAM_INFORMED_FASTA)}


def external_bp(seq, strain, ext):
    return {k: (dc.core_bp(seq, db[strain]) if strain in db else None) for k, db in ext.items()}


def ext_agrees(d):
    return any(v is not None and v <= EXT_AGREE_BP for v in d.values())


def gather():
    with open(OUT / "m02_mono_populations.pkl", "rb") as fh:
        mono = pickle.load(fh)
    with open(OUT / "p02_pair_observations.pkl", "rb") as fh:
        pair = pickle.load(fh)
    wells = pd.read_csv(OUT / "m01_mono_wells.csv").set_index("sample_id")
    obs = collections.defaultdict(list)
    for sid, plist in mono.items():
        strain, plate = wells.at[sid, "strain"], int(wells.at[sid, "plate"])
        for k, p in enumerate(plist):
            if k == 0 or p["frac"] >= MONO_MINOR_FRAC:
                obs[strain].append(dict(p, sample_id=sid, plate=plate,
                                        kind="mono_dominant" if k == 0 else "mono_minor"))
    for o in pair:
        if o["credited"]:
            obs[o["strain"]].append(dict(o, kind="pair"))
    return obs, wells


def cluster(olist):
    rank = {"mono_dominant": 0, "pair": 1, "mono_minor": 2}
    olist = sorted(olist, key=lambda o: (rank[o["kind"]], -o["n_reads"]))
    clusters = []
    for o in olist:
        for c in clusters:
            ok, plain, het = agree(o, c["rep"])
            if ok:
                c["mem"].append(o); o["_plain"], o["_het"] = plain, het
                break
        else:
            clusters.append(dict(rep=o, mem=[o])); o["_plain"] = o["_het"] = 0
    for c in clusters:
        c["wells"] = {o["sample_id"] for o in c["mem"]}
        c["mono_dom"] = {o["sample_id"] for o in c["mem"] if o["kind"] == "mono_dominant"}
        c["pair_wells"] = {o["sample_id"] for o in c["mem"] if o["kind"] == "pair"}
        c["plates"] = {o["plate"] for o in c["mem"]}
    clusters.sort(key=lambda c: (len(c["wells"]), len(c["mono_dom"]),
                                 sum(o["n_reads"] for o in c["mem"])), reverse=True)
    return clusters


def final_sequence(win):
    """Monoculture reads if the winner has any (cleaner -- no partner to misassign),
    otherwise the pooled pair-well groups."""
    monos = [o for o in win["mem"] if o["kind"] == "mono_dominant"]
    src = monos if monos else win["mem"]
    reads = [r for o in src for r in o["reads"]]
    return dc.top_q_consensus(reads, n=80), reads


def consistent_het(final, members):
    """Het sites seen in >=2 observations and in >=half of those deep enough to call.

    The minor allele is re-expressed against the FINAL consensus. At a ~50/50 site an
    observation's majority base can be the other one, and its "minor" is then the final
    consensus's own base -- carried over unchanged, that recorded G15 position 431 as
    A/A and E23 positions 427-447 likewise, which silently dropped those sites from the
    IUPAC reference and from haplotype phasing."""
    called = [o for o in members if o["het"] is not None]
    if len(called) < 2:
        return []
    votes = collections.defaultdict(list)
    for o in called:
        m = dc.pos_map(o["consensus"], final)
        for p, f, alt in o["het"]:
            q = m.get(p)
            if q is None:
                continue
            alleles = {o["consensus"][p], alt}
            if final[q] not in alleles:
                continue
            other = (alleles - {final[q]}).pop()
            votes[(q, other)].append(f if o["consensus"][p] == final[q] else 1 - f)
    need = max(2, int(np.ceil(0.5 * len(called))))
    return sorted((q, round(float(np.median(fs)), 3), alt, len(fs), len(called))
                  for (q, alt), fs in votes.items() if len(fs) >= need)


def main():
    obs, wells = gather()
    ext = load_external()
    lay_strains = sorted(set(pd.read_csv(config.LAYOUT_CSV).query(
        "plate_role == 'primary_sequenced'")[["strain1", "strain2"]].stack()))
    no_growth = set(pd.read_csv(config.NO_GROWTH_CSV).strain)
    mono_by_strain = wells.reset_index().groupby("strain")

    rows, orow, hrow, fasta, fasta_iu = [], [], [], [], []
    profiles, haprow, fasta_hap = {}, [], []
    for strain in lay_strains:
        olist = obs.get(strain, [])
        mw = mono_by_strain.get_group(strain) if strain in mono_by_strain.groups else None
        n_mono_data = int((mw.n_pops > 0).sum()) if mw is not None else 0
        base = dict(strain=strain, mono_wells_with_data=n_mono_data,
                    no_growth=strain in no_growth)
        if not olist:
            rows.append(dict(base, tier="none")); continue
        cl = cluster(olist)
        tie = len(cl) > 1 and ((len(cl[0]["wells"]), len(cl[0]["mono_dom"]))
                               == (len(cl[1]["wells"]), len(cl[1]["mono_dom"])))
        tiebreak = ""
        if tie:
            e0 = ext_agrees(external_bp(cl[0]["rep"]["consensus"], strain, ext))
            e1 = ext_agrees(external_bp(cl[1]["rep"]["consensus"], strain, ext))
            if e1 and not e0:
                cl[0], cl[1] = cl[1], cl[0]
            tiebreak = "external" if e0 != e1 else "unresolved"
        win = cl[0]
        n_md, n_pw, n_pl = len(win["mono_dom"]), len(win["pair_wells"]), len(win["plates"])
        win_ext = ext_agrees(external_bp(win["rep"]["consensus"], strain, ext))
        if n_md >= 2:
            tier, status = "A", "monos_agree"
        elif n_md == 1 and n_pw >= 2:
            tier, status = "B", "mono_plus_pairs"
        elif n_md == 0 and n_pw >= 3 and n_pl >= 2:
            tier, status = "C", "pairs_only"
        elif tie and tiebreak == "unresolved":
            tier, status = "D", "monos_disagree"
        elif n_md == 1 and win_ext:
            tier, status = "C", "single_mono_ext"
        elif n_md == 1:
            tier, status = "D", "single_mono"
        elif n_md == 0:
            tier, status = "D", "few_pair_wells"
        else:
            tier, status = "D", "weak"
        second = cl[1] if len(cl) > 1 else None
        second_bp = (dc.core_bp(second["rep"]["iupac"], win["rep"]["iupac"])
                     if second else None)
        competing = bool(second and second_bp > COMPETING_MIN_BP
                         and (len(second["wells"]) >= COMPETING_MIN_WELLS
                              or second["mono_dom"]))
        mixed_source = False
        if mw is not None and n_mono_data >= 2:
            mixed_source = bool((mw[mw.n_pops > 0].second_frac >= MIXED_SOURCE_FRAC).all())

        seq, reads = final_sequence(win)
        het = consistent_het(seq, win["mem"])
        prof_reads = [r for r, _ in sorted(reads, key=lambda x: -x[1])[:dc.HET_SCAN_READS]]
        alt, frac, cov, freqs = dc.allele_profile(seq, prof_reads, full=True)
        profiles[strain] = dict(consensus=seq, alt=alt, frac=frac.astype(np.float32),
                                cov=cov, freqs=freqs, n_reads=len(prof_reads),
                                source="mono" if win["mono_dom"] else "pair")
        haps = dc.phase_haplotypes(seq, [(q, f, alt) for q, f, alt, _, _ in het], prof_reads)
        hap_div = max((dc.core_bp(a[0], b[0]) for a in haps for b in haps), default=0)
        for k, (hs, hf, pat) in enumerate(haps):
            haprow.append(dict(strain=strain, haplotype=f"{strain}_h{k}", frac=round(hf, 3),
                               pattern=pat, to_consensus_bp=dc.core_bp(hs, seq)))
            fasta_hap.append(f">{strain}_h{k} strain={strain} frac={hf:.3f} tier={tier}\n{hs}\n")
        iu = dc.iupac(seq, [(q, f, alt) for q, f, alt, _, _ in het])
        drift = min(dc.core_bp(iu, o["iupac"]) for o in win["mem"])
        ext_bp = external_bp(seq, strain, ext)

        rows.append(dict(
            base, tier=tier, status=status,
            wells_supporting=len(win["wells"]), monos_supporting=n_md,
            pair_wells_supporting=n_pw, plates_supporting=n_pl,
            observations=len(olist), n_clusters=len(cl),
            competing=competing,
            second_wells=len(second["wells"]) if second else 0,
            second_monos=len(second["mono_dom"]) if second else 0,
            second_vs_winner_bp=second_bp,
            tiebreak=tiebreak or None,
            mixed_source=mixed_source, het_sites=len(het), many_het=len(het) > MANY_HET,
            n_haplotypes=len(haps), haplotype_divergence_bp=hap_div,
            final_vs_nearest_obs_bp=drift, length=len(seq),
            **{f"{k}_bp": v for k, v in ext_bp.items()},
            external_agrees=ext_agrees(ext_bp)))
        for k, c in enumerate(cl):
            for o in c["mem"]:
                orow.append(dict(strain=strain, cluster=k, kind=o["kind"],
                                 sample_id=o["sample_id"], plate=o["plate"],
                                 n_reads=o["n_reads"], frac=round(o["frac"], 4),
                                 to_rep_plain_bp=o["_plain"], to_rep_shared_het_bp=o["_het"],
                                 het_sites=len(o["het"]) if o["het"] is not None else None))
        for q, f, alt, nv, nc in het:
            hrow.append(dict(strain=strain, position=q, major=seq[q], minor=alt,
                             median_minor_frac=f, wells_het=nv, wells_called=nc))
        hdr = (f"{strain} tier={tier} status={status} wells={len(win['wells'])} "
               f"monos={n_md} pairs={n_pw} plates={n_pl} het={len(het)}"
               + (" competing" if competing else "") + (" mixed_source" if mixed_source else "")
               + (" no_growth" if strain in no_growth else ""))
        fasta.append(f">{hdr}\n{seq}\n"); fasta_iu.append(f">{hdr}\n{iu}\n")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s01_strain_reference.csv", index=False)
    pd.DataFrame(orow).to_csv(OUT / "s02_strain_observations.csv", index=False)
    pd.DataFrame(hrow).to_csv(OUT / "s03_strain_het_sites.csv", index=False)
    pd.DataFrame(haprow).to_csv(OUT / "s05_haplotypes.csv", index=False)
    (OUT / "reference_db_20260907_haplotypes.fasta").write_text("".join(fasta_hap))
    with open(OUT / "s04_allele_profiles.pkl", "wb") as fh:
        pickle.dump(profiles, fh)
    (OUT / "reference_db_20260907.fasta").write_text("".join(fasta))
    (OUT / "reference_db_20260907_iupac.fasta").write_text("".join(fasta_iu))
    report(df)
    return df


def report(df):
    print(f"{len(df)} strains in the 20260907 layout")
    print("tier:", df.tier.value_counts().sort_index().to_dict())
    print(df.groupby(["tier", "status"], dropna=False).size().to_string())
    has = df[df.tier != "none"]
    for flag in ("competing", "mixed_source", "many_het", "no_growth"):
        sub = has[has[flag].astype(bool)]
        print(f"  {flag:13s} {len(sub):4d}   by tier {sub.tier.value_counts().sort_index().to_dict()}")
    print(f"  het sites per strain: {has.het_sites.value_counts().sort_index().to_dict()}")
    print(f"  haplotypes per strain: {has.n_haplotypes.value_counts().sort_index().to_dict()}")
    print(f"  max divergence between a strain's haplotypes: "
          f"{pd.cut(has.haplotype_divergence_bp, [-1, 0, 2, 5, 10, 20, 1000]).value_counts().sort_index().to_dict()}")
    for t, g in has.groupby("tier"):
        print(f"  tier {t}: agrees with an older independent DB (<= {EXT_AGREE_BP} bp) "
              f"{g.external_agrees.mean():.0%}")
    print(f"  final sequence vs nearest supporting observation: "
          f"{has.final_vs_nearest_obs_bp.value_counts().sort_index().to_dict()}")


if __name__ == "__main__":
    main()
