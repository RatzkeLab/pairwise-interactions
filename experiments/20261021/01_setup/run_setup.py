"""Build the community experiment's design from setup_config.yaml.

    python run_setup.py            # refuses to overwrite outputs/
    python run_setup.py --force    # allowed only while `experiment: DRAFT`

Writes outputs/ (layout, Echo picklists, primer picklists, minibar sheet, demux config),
runs every check in scripts/experiment_setup/validate.py, and exits non-zero (after writing
validation.txt) if any fails. Deterministic for a given config.
"""
import argparse
import datetime
import hashlib
import itertools
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
KARL = Path("/home/rl/scripts/karl")
sys.path.insert(0, str(KARL / "pairwise_interaction_experiments" / "scripts"))

import numpy as np
import pandas as pd
import yaml

from experiment_setup import (barcodes as B, communities as C, echo as E, layout as L,
                              plates, pool as P, validate as V)


def kp(rel):
    return KARL / rel


def split_counts(total, sizes, ext_frac):
    per, rem = divmod(total, len(sizes))
    counts = {k: per + (1 if i < rem else 0) for i, k in enumerate(sizes)}
    ext = {k: int(round(n * ext_frac)) for k, n in counts.items()}
    core = {k: counts[k] - ext[k] for k in counts}
    return core, ext


def sanity_pair_candidates(cfg, used):
    """Pairs measured in 20260907 with both strains in today's pool and both clearly present
    (5-95%), so both 16S are visible: a check that the source plate still holds the same
    strains, and a cross-experiment replicate. Technical-replicate pairs (5 wells on
    different plates) first, then the rest, most reproducible / most replicated first.
    Returns (ordered candidate pairs, every pair measured in 20260907)."""
    lay = pd.read_csv(kp(cfg["inputs"]["layout_20260907"]))
    t = lay[lay.well_type == "techrep"]
    tr = {frozenset(p) for p in zip(t.strain1, t.strain2)}
    r = pd.read_csv(kp(cfg["inputs"]["pairs_20260907"]))
    r["key"] = [frozenset(p) for p in zip(r.strain_a, r.strain_b)]
    measured = set(r.key)
    r = r[r.strain_a.isin(used) & r.strain_b.isin(used)
          & r.mean_relative_abundance_a.between(0.05, 0.95)].copy()
    r["techrep"] = r.key.isin(tr)
    r = r.sort_values(["techrep", "n_replicates", "std_relative_abundance_a"],
                      ascending=[False, False, True], na_position="last")
    return list(zip(r.strain_a, r.strain_b)), measured


def ko_distance(cfg, strains):
    idx = pd.read_csv(kp(cfg["inputs"]["genome_index"])).set_index("well").genome_table_id.astype(str)
    ko = pd.read_csv(kp(cfg["inputs"]["ko_table"]), index_col=0)
    ko.index = ko.index.astype(str)
    X = (ko.reindex(idx.reindex(strains).values).fillna(0).values > 0)
    inter = X.astype(int) @ X.T.astype(int)
    size = X.sum(1)
    union = size[:, None] + size[None, :] - inter
    return pd.DataFrame(1 - inter / np.maximum(union, 1), index=strains, columns=strains)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(HERE / "setup_config.yaml"))
    exp = str(cfg["experiment"])
    out = HERE / "outputs"
    if out.exists():
        if not (args.force and exp.startswith("DRAFT")):
            sys.exit("outputs/ exists; --force is only allowed for a DRAFT experiment")
        shutil.rmtree(out)
    out.mkdir()
    rng = np.random.default_rng(cfg["seed"])
    sep, nrep, blanks = cfg["separation_bp"], cfg["replicates"], cfg["blanks_per_plate"]
    per_plate = (len(plates.inner_wells()) - blanks) // nrep * nrep
    inoc = cfg["inoculum"]
    max_members = int(inoc["max_well_inoculum_uL"] * 1000 // inoc["reduced_nL"])

    # ---- pool -------------------------------------------------------------------------
    od_sept = P.source_plate_od([kp(f) for f in cfg["inputs"]["source_od_files"]])
    pool = P.load_pool(kp(cfg["inputs"]["eligibility"]), od_sept, cfg["pool"]["min_source_od"],
                       cfg["pool"]["include_extended"])
    used = pool[~pool.get("low_preculture_density", pd.Series(False, index=pool.index))].copy()
    pool.to_csv(out / "pool.csv", index=False)
    D = P.distance_matrix(kp(cfg["inputs"]["distances"]), list(used.strain))
    bld = C.Builder(used, D, rng)
    blocks = {"main": [], "passage": []}

    def take(recs, block):
        for r in recs:
            if r["community_id"] not in blocks[block]:
                blocks[block].append(r["community_id"])

    # ---- main plates ------------------------------------------------------------------
    m = cfg["modules"]
    n_main = len(cfg["blocks"]["main"]["plates"]) * per_plate // nrep
    pairs, measured = sanity_pair_candidates(cfg, set(used.strain))
    take(bld.known_pairs(pairs, m["sanity_pairs"]), "main")
    take(bld.trios_with_pairs(m["trios"], "core", sep, measured), "main")
    take([r for r in bld.order if r["modules"][0] == "trio_pair"], "main")
    for lo in m["leave_out"]:
        take(bld.leave_out_series(lo["base"], lo["depth"], "core", sep), "main")
    rt = m["resolution_tests"]
    take(bld.resolution_tests([tuple(b) for b in rt["bands"]], rt["n_per_band"], rt["size"],
                              "core", sep), "main")
    large = bld.large(m["large_tests"], cap=max_members)
    take(large, "main")
    mid = m["mid"]
    core_mid, ext_mid = split_counts(mid["n"], mid["sizes"], cfg["small"]["extended_fraction"])
    take(bld.stratified(core_mid, "core", sep, "mid", 1), "main")
    take(bld.stratified(ext_mid, "extended", sep, "mid", 1), "main")
    small = cfg["small"]
    while len(blocks["main"]) < n_main:
        core, ext = split_counts(n_main - len(blocks["main"]), small["sizes"],
                                 small["extended_fraction"])
        before = len(blocks["main"])
        take(bld.stratified(core, "core", sep, "small", small["relatedness_strata"]), "main")
        take(bld.stratified(ext, "extended", sep, "small", small["relatedness_strata"]), "main")
        if len(blocks["main"]) == before:
            raise RuntimeError("could not fill the main plates")
    blocks["main"] = blocks["main"][:n_main]

    # ---- passage plate ---------------------------------------------------------------
    pc = cfg["blocks"]["passage"]
    n_pass = per_plate // nrep
    recs = {r["community_id"]: r for r in bld.order}
    pool_of = used.set_index("strain")["pool"]
    is_core = lambda c: all(pool_of[x] == "core" for x in recs[c]["members"])
    if pc.get("large_from_main"):
        take([r for r in large if is_core(r["community_id"])], "passage")
    main_small_core = [c for c in blocks["main"] if recs[c]["modules"][0] == "small" and is_core(c)]
    n_over = int(round(pc["overlap_with_main"] * (n_pass - len(blocks["passage"]))))
    by_size = {}
    for c in main_small_core:
        by_size.setdefault(recs[c]["size"], []).append(c)
    over = []
    while len(over) < n_over:
        for k in sorted(by_size):
            if by_size[k] and len(over) < n_over:
                over.append(by_size[k].pop(int(rng.integers(len(by_size[k])))))
    blocks["passage"] += over
    while len(blocks["passage"]) < n_pass:
        core, _ = split_counts(n_pass - len(blocks["passage"]), small["sizes"], 0.0)
        take(bld.stratified(core, pc["pool"], sep, "small", small["relatedness_strata"]), "passage")
    blocks["passage"] = blocks["passage"][:n_pass]

    # ---- community table -------------------------------------------------------------
    keep = set(itertools.chain.from_iterable(blocks.values()))
    comm = bld.table()
    comm = comm[comm.community_id.isin(keep)].copy()
    comm["blocks"] = [";".join(b for b in blocks if c in blocks[b]) for c in comm.community_id]
    comm["member_nL"] = [inoc["member_nL"] if k * inoc["member_nL"] <= inoc["max_well_inoculum_uL"] * 1000
                         else inoc["reduced_nL"] for k in comm["size"]]
    comm["inoculum_uL"] = comm["size"] * comm.member_nL / 1000
    Kd = ko_distance(cfg, list(used.strain))
    comm["ko_jaccard_distance"] = [
        float(np.mean([Kd.at[a, b] for a, b in itertools.combinations(r.split(), 2)]))
        for r in comm.members]
    comm["pairs_measured_20260907"] = [
        sum(frozenset(p) in measured for p in itertools.combinations(r.split(), 2))
        for r in comm.members]
    comm.to_csv(out / "communities.csv", index=False)

    # ---- layout ------------------------------------------------------------------------
    cdict = comm.set_index("community_id").to_dict("index")
    cl = lambda ids: [dict(community_id=c, size=cdict[c]["size"]) for c in ids]
    rows = L.across_plates(cl(blocks["main"]), cfg["blocks"]["main"]["plates"], nrep, rng, blanks, "main")
    rows += L.within_plate(cl(blocks["passage"]), pc["plate"], nrep, rng, blanks, "passage")
    lay = L.to_frame(rows, comm.assign(premix=False))
    lay = lay.drop(columns=["premix"])
    lay["inoculated"] = True
    lay["timepoint"] = "24h"
    tp = {int(k): v for k, v in pc["timepoints"].items()}
    lay.loc[lay.plate == pc["plate"], "timepoint"] = tp[pc["plate"]]
    copies = []
    for p, t in tp.items():
        if p != pc["plate"]:
            c = lay[lay.plate == pc["plate"]].copy()
            c["plate"], c["timepoint"], c["inoculated"] = p, t, False
            copies.append(c)
    lay = pd.concat([lay] + copies, ignore_index=True).sort_values(["plate", "row", "col"])
    lay = lay.reset_index(drop=True)
    lay["sample_id"] = [f"Plate{int(p):02d}_{w}" for p, w in zip(lay.plate, lay.well)]

    # ---- barcodes ----------------------------------------------------------------------
    table = B.barcode_table(kp(cfg["inputs"]["primer_specs"]))
    failed = B.failed_barcodes(kp(cfg["inputs"]["failed_barcodes_fwd"]),
                               kp(cfg["inputs"]["failed_barcodes_rev"]))
    decoy = [(cfg["decoy_barcode_plate"], w) for w in plates.inner_wells()]
    assigned = B.assign(list(zip(lay.plate, lay.well)) + decoy, rng, failed)
    real = assigned[:len(lay)]
    fwd_w, rev_w = B.primer_source_wells()
    lay["fwd_idx"] = [a[2] for a in real]; lay["rev_idx"] = [a[3] for a in real]
    lay["fwd_source_well"] = [fwd_w[i] for i in lay.fwd_idx]
    lay["rev_source_well"] = [rev_w[j] for j in lay.rev_idx]
    lay.to_csv(out / "layout.csv", index=False)
    sheet = B.minibar_sheet(assigned, table)
    sheet.to_csv(out / f"minibar_sheet_{exp}.tsv", sep="\t", index=False)
    demux = B.demux_config(out / f"minibar_sheet_{exp}.tsv",
                           f"/home/rl/data/interim/karl/demultiplexing/{exp}_demux/",
                           barcode_edit_distance=cfg["demux"]["barcode_edit_distance"])
    with open(out / f"demultiplex_config_{exp}.yaml", "w") as fh:
        fh.write(f"# Plate{cfg['decoy_barcode_plate']}_* rows are decoys: barcode combinations put in no\n"
                 "# well. Reads assigned to them measure demultiplexing mis-assignment.\n")
        yaml.safe_dump(demux, fh, sort_keys=False)
    pdir = out / "primer_shooting"; pdir.mkdir()
    prim = E.primer_transfers(real, cfg["primers"]["volume_nL"], fwd_w, rev_w)
    for p, g in prim.groupby("Destination Plate Name"):
        g.to_csv(pdir / f"echo_primers_{exp}_plate_{int(p):02d}.csv", index=False)

    # ---- strain transfers (equal volumes) -------------------------------------------------
    plan = E.transfer_plan(lay, comm, {})
    plan = E.equal_volumes(plan, comm, inoc["member_nL"], inoc["reduced_nL"],
                           inoc["max_well_inoculum_uL"])
    order = [str(p) for p in sorted(lay[lay.inoculated].plate.unique())]
    plan = E.fills(plan, order, inoc["fill_cap_uL"])
    plan.to_csv(out / "strain_transfers.csv", index=False)
    edir = out / "echo_strains"; edir.mkdir()
    cols = {"source_well": "Source Well", "dest_plate": "Destination Plate Name",
            "dest_well": "Destination Well", "volume_nL": "Transfer Volume"}
    nfill = int(plan.fill.max())
    for (f, p), g in plan.groupby(["fill", "dest_plate"]):
        tag = f"_fill{f}" if nfill > 1 else ""
        g[list(cols)].rename(columns=cols).to_csv(
            edir / f"echo_strains_{exp}{tag}_plate_{int(p):02d}.csv", index=False)
    plan[list(cols)].rename(columns=cols).to_csv(out / f"echo_strains_{exp}_all_plates.csv", index=False)

    # ---- validation ------------------------------------------------------------------------
    inoc_plates = cfg["blocks"]["main"]["plates"] + [pc["plate"]]
    checks = {
        "communities": V.check_communities(comm, used, D, sep),
        "layout": V.check_layout(lay[lay.plate.isin(inoc_plates)], nrep, {"passage"}, blanks),
        "barcodes": V.check_barcodes(assigned, failed),
        "echo": V.check_echo(plan, inoc, inoc["reduced_nL"]),
    }
    mb_errs, mb_info = V.check_minibar(sheet, cfg["demux"]["barcode_edit_distance"])
    checks["minibar"] = mb_errs
    report = summarize(cfg, pool, used, comm, lay, blocks, failed, mb_info, plan, max_members)
    vlines = [f"{k}: {'PASS' if not v else 'FAIL'}" + "".join(f"\n    {e}" for e in v[:20])
              for k, v in checks.items()]
    (out / "validation.txt").write_text("\n".join(vlines) + "\n")
    (out / "summary.txt").write_text(report)
    manifest = dict(
        experiment=exp, created=datetime.datetime.now().isoformat(timespec="seconds"),
        config=cfg, inputs={k: hashlib.md5(open(kp(v), "rb").read()).hexdigest()
                            for k, v in cfg["inputs"].items() if isinstance(v, str)},
        outputs={str(p.relative_to(out)): hashlib.md5(p.read_bytes()).hexdigest()
                 for p in sorted(out.rglob("*")) if p.is_file()})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
    print(report)
    print("\n".join(vlines))
    if any(checks.values()):
        sys.exit("VALIDATION FAILED -- see outputs/validation.txt")


def summarize(cfg, pool, used, comm, lay, blocks, failed, mb_info, plan, max_members):
    S = [f"experiment {cfg['experiment']}  seed {cfg['seed']}",
         f"pool: {len(used)} genomes ({(used['pool'] == 'core').sum()} core, "
         f"{(used['pool'] == 'extended').sum()} extended)"]
    for b, ids in blocks.items():
        c = comm[comm.community_id.isin(ids)]
        S.append(f"\n{b}: {len(ids)} communities x {cfg['replicates']} = {len(ids) * cfg['replicates']} wells")
        S.append("  by module: " + str(c.module.value_counts().to_dict()))
        S.append("  by size:   " + str(c["size"].value_counts().sort_index().to_dict()))
        S.append(f"  core / extended: {(c.pool == 'core').sum()} / {(c.pool == 'extended').sum()}")
        S.append(f"  MPD (bp) 5/50/95%: {c.mpd_bp.quantile([.05, .5, .95]).round(0).tolist()}")
    mm = comm[comm.blocks.str.contains("main")]
    for mod in ("trio", "trio_pair", "sanity_pair", "leave_out", "mid", "large_test", "resolution_test"):
        n = (mm.module == mod).sum()
        S.append(f"  main-plate wells, {mod:15s} {n * cfg['replicates']:4d} "
                 f"({n * cfg['replicates'] / (len(blocks['main']) * cfg['replicates']):.1%})")
    inc = pd.Series([s for r in lay[(lay.well_type == "community") & lay.inoculated].members
                     for s in r.split()]).value_counts()
    S.append(f"\nwells per strain (all inoculated plates): min {inc.min()}, median {inc.median():.0f}, "
             f"max {inc.max()}; strains used {len(inc)} / {len(used)}")
    draw = plan.groupby("strain").volume_nL.sum() / 1000
    well = plan.groupby(["dest_plate", "dest_well"]).volume_nL.sum() / 1000
    S.append(f"inoculum: {cfg['inoculum']['member_nL']} nL per member, {cfg['inoculum']['reduced_nL']} nL "
             f"for communities over {int(cfg['inoculum']['max_well_inoculum_uL'] * 1000 // cfg['inoculum']['member_nL'])} "
             f"members (cap {max_members} members)")
    S.append(f"  per well: median {well.median():.2f} uL, max {well.max():.2f} uL")
    S.append(f"  drawn per source well: median {draw.median():.1f} uL, max {draw.max():.1f} uL "
             f"({draw.idxmax()}); source fills needed: {int(plan.fill.max())}"
             + ("  -> NO REFILL" if plan.fill.max() == 1 else "  -> refill between fills"))
    S.append(f"  Echo strain transfers: {len(plan):,}, droplets {int(plan.volume_nL.sum() / 25):,} "
             f"(20260907: 9,096 / 38,040)")
    S.append(f"plates: inoculated {sorted(int(p) for p in lay[lay.inoculated].plate.unique())}; "
             f"barcoded {sorted(int(p) for p in lay.plate.unique())} + decoy plate {cfg['decoy_barcode_plate']}")
    S.append(f"barcodes excluded (failed in 20260907): fwd {sorted(failed['fwd'])}, rev {sorted(failed['rev'])}")
    S.append(f"minibar: {mb_info}")
    return "\n".join(S) + "\n"


if __name__ == "__main__":
    main()
