"""Build 20261026 from setup_config.yaml: 20261021's plates 1-5 well for well, plus a second
passage plate (6) built the way plate 5 was, part of it shared with plate 5.

    python run_setup.py            # refuses to overwrite outputs/
    python run_setup.py --force    # allowed only while `experiment: DRAFT`

Plates 1-5 are read from 20261021's outputs, not regenerated, so they cannot drift; the
checks confirm they (and their Echo strain transfers) are identical. Writes outputs/,
runs every check in scripts/experiment_setup/validate.py plus the copy checks, and exits
non-zero (after writing validation.txt) if any fails. Deterministic for a given config.
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
                              plates, validate as V)
from experiment_setup import pool as P

LAYOUT_COLS = ["plate", "well", "block", "community_id", "replicate", "well_type", "size",
               "members", "module", "pool", "mpd_bp", "min_sep_bp", "row", "col",
               "inoculated", "timepoint"]


def kp(rel):
    return KARL / rel


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def measured_pairs(cfg):
    r = pd.read_csv(kp(cfg["inputs"]["pairs_20260907"]))
    return {frozenset(p) for p in zip(r.strain_a, r.strain_b)}


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
    main_plates = cfg["blocks"]["main"]["plates"]
    pc = cfg["blocks"]["passage"]
    new_p = pc["new_plate"]
    bmap = {int(k): (int(v[0]), str(v[1])) for k, v in pc["barcode_plates"].items()}
    for p in pc["plates"]:
        if bmap.get(p, (None,))[0] != p:
            sys.exit(f"barcode plate {p} must be passage plate {p} itself (its first time point)")

    # ---- the experiment being copied ----------------------------------------------------
    cp = cfg["copy_from"]
    ref_exp, ref_out, copied = str(cp["experiment"]), kp(cp["outputs"]), cp["plates"]
    ref_manifest = json.load(open(ref_out / "manifest.json"))
    ref_comm = pd.read_csv(ref_out / "communities.csv")
    ref_lay = pd.read_csv(ref_out / "layout.csv")
    ref_tr = pd.read_csv(ref_out / "strain_transfers.csv")
    ref_sheet = pd.read_csv(ref_out / f"minibar_sheet_{ref_exp}.tsv", sep="\t")
    pool = pd.read_csv(ref_out / "pool.csv")
    shutil.copy(ref_out / "pool.csv", out / "pool.csv")
    used = pool
    if "low_preculture_density" in pool:
        used = pool[~pool.low_preculture_density.fillna(False).astype(bool)].copy()
    D = P.distance_matrix(kp(cfg["inputs"]["distances"]), list(used.strain))
    pool_of = used.set_index("strain")["pool"]

    old = ref_lay[ref_lay.plate.isin(copied) & ref_lay.inoculated][LAYOUT_COLS].copy()
    oc = old[old.well_type == "community"]
    on_plate = oc.groupby("community_id").plate.agg(lambda s: set(s)).to_dict()
    p5_ids = set(oc[oc.plate.isin([p for p in pc["plates"] if p != new_p])].community_id)
    main_ids = set(oc[oc.plate.isin(main_plates)].community_id)

    bld = C.Builder(used, D, rng)
    bld.preload(ref_comm)
    R = {}

    def rec(c):
        if c not in R:
            R.update({r["community_id"]: r for r in bld.order})
        return R[c]

    is_core = lambda c: all(pool_of[x] == "core" for x in rec(c)["members"])

    # ---- passage plate 6 -----------------------------------------------------------------
    n_pass = per_plate // nrep
    blk6, origin = [], {}

    def take(ids, why):
        for c in ids:
            if c not in blk6:
                blk6.append(c)
                origin[c] = why

    if pc.get("large_from_main"):
        take([c for c in ref_comm.community_id[ref_comm.module == "large_test"] if is_core(c)],
             "large test (also on plate 5)")
    p5_small = ref_comm[ref_comm.community_id.isin(p5_ids) & (ref_comm.module == "small")]
    k_sh = pc["shared_with_plate_5"]["small_per_size"]
    shared = []
    for k, g in p5_small.sort_values(["mpd_bp", "community_id"]).groupby("size"):
        for tercile in np.array_split(g.community_id.values, k_sh):   # relatedness within size
            shared.append(str(rng.choice(tercile)))
    take(shared, "shared with plate 5")

    n_over = int(round(pc["overlap_with_main"] * (n_pass - len(blk6))))
    main_small_core = [c for c in ref_comm.community_id
                       if c in main_ids and c not in p5_ids and c not in blk6
                       and rec(c)["modules"][0] == "small" and is_core(c)]
    by_size = {}
    for c in main_small_core:
        by_size.setdefault(rec(c)["size"], []).append(c)
    over = []
    while len(over) < n_over:
        for k in sorted(by_size):
            if by_size[k] and len(over) < n_over:
                over.append(by_size[k].pop(int(rng.integers(len(by_size[k])))))
    take(over, "also on main plates")

    small = cfg["small"]
    target = p5_small["size"].value_counts()          # plate 5's size mix
    while len(blk6) < n_pass:
        have = pd.Series([rec(c)["size"] for c in blk6 if rec(c)["modules"][0] == "small"],
                         dtype=int).value_counts()
        need = {k: int(max(0, target.get(k, 0) - have.get(k, 0))) for k in small["sizes"]}
        if not sum(need.values()):
            raise RuntimeError("plate 6 size quotas met but the plate is not full")
        before = len(blk6)
        take([r["community_id"] for r in bld.stratified(need, pc["pool"], sep, "small",
                                                        small["relatedness_strata"])], "new")
        if len(blk6) == before:
            raise RuntimeError("could not fill plate 6")
    blk6 = blk6[:n_pass]
    for c in blk6:
        on_plate.setdefault(c, set()).add(new_p)

    # ---- community table ---------------------------------------------------------------
    ref_ids = set(ref_comm.community_id)
    new_ids = [c for c in blk6 if c not in ref_ids]
    nt = bld.table()
    nt = nt[nt.community_id.isin(new_ids)].copy()
    nt["member_nL"] = [inoc["member_nL"] if k * inoc["member_nL"] <= inoc["max_well_inoculum_uL"] * 1000
                       else inoc["reduced_nL"] for k in nt["size"]]
    nt["inoculum_uL"] = nt["size"] * nt.member_nL / 1000
    Kd = ko_distance(cfg, list(used.strain))
    measured = measured_pairs(cfg)
    nt["ko_jaccard_distance"] = [
        float(np.mean([Kd.at[a, b] for a, b in itertools.combinations(r.split(), 2)]))
        for r in nt.members]
    nt["pairs_measured_20260907"] = [
        sum(frozenset(p) in measured for p in itertools.combinations(r.split(), 2))
        for r in nt.members]
    comm = pd.concat([ref_comm, nt], ignore_index=True)[list(ref_comm.columns)]
    comm = comm[comm.community_id.isin(on_plate)].copy()
    comm["blocks"] = [";".join(b for b, ps in (("main", main_plates), ("passage", pc["plates"]))
                               if on_plate[c] & set(ps)) for c in comm.community_id]
    comm["passage_plates"] = [" ".join(str(p) for p in sorted(on_plate[c] & set(pc["plates"])))
                              for c in comm.community_id]
    comm[f"in_{ref_exp}"] = comm.community_id.isin(ref_ids)
    comm["plate6_origin"] = comm.community_id.map(origin).fillna("")
    comm.to_csv(out / "communities.csv", index=False)

    # ---- layout --------------------------------------------------------------------------
    sizes = comm.set_index("community_id")["size"]
    rows6 = L.within_plate([dict(community_id=c, size=int(sizes[c])) for c in blk6], new_p,
                           nrep, rng, blanks, "passage")
    new6 = L.to_frame(rows6, comm)
    new6["inoculated"] = True
    new6["timepoint"] = bmap[new_p][1]
    lay = pd.concat([old, new6[LAYOUT_COLS]], ignore_index=True)
    lay.insert(1, "culture_plate", lay.plate)
    copies = []
    for bp, (cplate, t) in bmap.items():
        if bp != cplate:
            c = lay[lay.plate == cplate].copy()
            c["plate"], c["timepoint"], c["inoculated"] = bp, t, False
            copies.append(c)
    lay = pd.concat([lay] + copies, ignore_index=True).sort_values(["plate", "row", "col"])
    lay = lay.reset_index(drop=True)
    lay["sample_id"] = [f"Plate{int(p):02d}_{w}" for p, w in zip(lay.plate, lay.well)]

    # ---- barcodes ------------------------------------------------------------------------
    table = B.barcode_table(kp(cfg["inputs"]["primer_specs"]))
    failed = B.failed_barcodes(kp(cfg["inputs"]["failed_barcodes_fwd"]),
                               kp(cfg["inputs"]["failed_barcodes_rev"]))
    ref_combos = B.combos_in_sheet(ref_sheet, table) if cp.get("disjoint_barcodes") else []
    decoy = [(cfg["decoy_barcode_plate"], w) for w in plates.inner_wells()]
    assigned = B.assign(list(zip(lay.plate, lay.well)) + decoy, rng, failed,
                        exclude_combos=ref_combos)
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
                 "# well. Reads assigned to them measure demultiplexing mis-assignment.\n"
                 f"# No combination on this sheet is on the {ref_exp} sheet.\n")
        yaml.safe_dump(demux, fh, sort_keys=False)
    pdir = out / "primer_shooting"; pdir.mkdir()
    prim = E.primer_transfers(real, cfg["primers"]["volume_nL"], fwd_w, rev_w)
    for p, g in prim.groupby("Destination Plate Name"):
        g.to_csv(pdir / f"echo_primers_{exp}_plate_{int(p):02d}.csv", index=False)

    # ---- strain transfers (equal volumes) ------------------------------------------------
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

    # ---- validation ----------------------------------------------------------------------
    inoc_plates = main_plates + pc["plates"]
    checks = {
        "inputs": check_inputs(cfg, ref_manifest, ref_out, ref_exp),
        "communities": V.check_communities(comm, used, D, sep),
        "layout": V.check_layout(lay[lay.plate.isin(inoc_plates)], nrep, {"passage"}, blanks),
        "copy": check_copied(lay, ref_lay, plan, ref_tr, comm, ref_comm, copied),
        "plate 6": check_plate6(lay, plan, comm, pc, new_p, k_sh, len(small["sizes"])),
        "barcodes": V.check_barcodes(assigned, failed, ref_combos),
        "echo": V.check_echo(plan, inoc, inoc["reduced_nL"]),
    }
    mb_errs, mb_info = V.check_minibar(sheet, cfg["demux"]["barcode_edit_distance"])
    checks["minibar"] = mb_errs
    report = summarize(cfg, used, comm, lay, plan, prim, failed, mb_info, ref_combos, ref_lay,
                       bmap, new_p)
    vlines = [f"{k}: {'PASS' if not v else 'FAIL'}" + "".join(f"\n    {e}" for e in v[:20])
              for k, v in checks.items()]
    (out / "validation.txt").write_text("\n".join(vlines) + "\n")
    (out / "summary.txt").write_text(report)
    manifest = dict(
        experiment=exp, created=datetime.datetime.now().isoformat(timespec="seconds"),
        config=cfg,
        inputs={k: md5(kp(v)) for k, v in cfg["inputs"].items() if isinstance(v, str)},
        copied_from={f: md5(ref_out / f) for f in
                     ("communities.csv", "layout.csv", "pool.csv", "strain_transfers.csv",
                      f"minibar_sheet_{ref_exp}.tsv", "manifest.json")},
        outputs={str(p.relative_to(out)): md5(p) for p in sorted(out.rglob("*")) if p.is_file()})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
    print(report)
    print("\n".join(vlines))
    if any(checks.values()):
        sys.exit("VALIDATION FAILED -- see outputs/validation.txt")


def check_inputs(cfg, ref_manifest, ref_out, ref_exp):
    """Plate 6 is built on the same files as the copied plates, and the copied outputs are
    the ones that experiment wrote."""
    errs = []
    for k, v in cfg["inputs"].items():
        if k in ref_manifest["inputs"] and md5(kp(v)) != ref_manifest["inputs"][k]:
            errs.append(f"input {k} changed since {ref_exp} was built")
    for f in ("communities.csv", "layout.csv", "pool.csv", "strain_transfers.csv",
              f"minibar_sheet_{ref_exp}.tsv"):
        if md5(ref_out / f) != ref_manifest["outputs"][f]:
            errs.append(f"{ref_exp} {f} differs from its manifest")
    return errs


def check_copied(lay, ref_lay, plan, ref_tr, comm, ref_comm, copied):
    """Plates 1-5: same community, replicate and members in every well, and the same Echo
    strain transfers; every copied community keeps its members."""
    cols = ["plate", "well", "block", "community_id", "replicate", "well_type", "size", "members"]
    sel = lambda d: d[d.plate.isin(copied) & d.inoculated]
    errs = V.check_copy(sel(lay), sel(ref_lay), cols, "layout")
    tcols = ["dest_plate", "dest_well", "community_id", "strain", "source_well", "stage", "volume_nL"]
    a = plan[plan.dest_plate.astype(int).isin(copied)][tcols].astype(str)
    errs += V.check_copy(a, ref_tr[tcols].astype(str), tcols, "Echo strain transfers")
    m = comm.set_index("community_id")
    r = ref_comm.set_index("community_id")
    for c in ("members", "size", "module", "member_nL"):
        if not (m.loc[r.index, c] == r[c]).all():
            errs.append(f"communities: `{c}` changed for a copied community")
    return errs


def check_plate6(lay, plan, comm, pc, new_p, k_sh, n_sizes):
    """Shared communities: on both passage plates, in full triplicate, with the same
    per-member volumes; core strains only on the plate."""
    errs = []
    sh = comm[comm.passage_plates == " ".join(str(p) for p in sorted(pc["plates"]))]
    n_large = (sh.module == "large_test").sum()
    want_large = ((comm.module == "large_test") & (comm.pool == "core")).sum()
    want = k_sh * n_sizes
    if (sh.module == "small").sum() != want:
        errs.append(f"{(sh.module == 'small').sum()} small communities shared, expected {want}")
    if pc.get("large_from_main") and n_large != want_large:
        errs.append(f"{n_large} large tests shared, expected the {want_large} core ones")
    six = lay[(lay.plate == new_p) & (lay.well_type == "community")]
    if (six.pool != "core").any():
        errs.append("plate 6 holds an extended-pool community")
    d = plan[plan.dest_plate.astype(int).isin(pc["plates"])]
    vol = d.groupby(["dest_plate", "community_id"]).apply(
        lambda g: tuple(sorted(g.groupby("strain").volume_nL.sum().items())), include_groups=False)
    for c in sh.community_id:
        v = {vol.get((str(p), c)) for p in pc["plates"]}
        if len(v) != 1 or None in v:
            errs.append(f"{c}: inoculum differs between passage plates")
    return errs


def summarize(cfg, used, comm, lay, plan, prim, failed, mb_info, ref_combos, ref_lay, bmap, new_p):
    S = [f"experiment {cfg['experiment']}  seed {cfg['seed']}",
         f"plates {cfg['copy_from']['plates']} copied well for well from {cfg['copy_from']['experiment']}; "
         f"plate {new_p} new",
         f"pool: {len(used)} genomes ({(used['pool'] == 'core').sum()} core, "
         f"{(used['pool'] == 'extended').sum()} extended)"]
    inoc = lay[lay.inoculated & (lay.well_type == "community")]
    groups = {"main (plates " + ", ".join(map(str, cfg["blocks"]["main"]["plates"])) + ")":
              set(inoc[inoc.block == "main"].community_id)}
    for p in cfg["blocks"]["passage"]["plates"]:
        groups[f"passage plate {p}"] = set(inoc[inoc.plate == p].community_id)
    for b, ids in groups.items():
        c = comm[comm.community_id.isin(ids)]
        S.append(f"\n{b}: {len(ids)} communities x {cfg['replicates']} = {len(ids) * cfg['replicates']} wells")
        S.append("  by module: " + str(c.module.value_counts().to_dict()))
        S.append("  by size:   " + str(c["size"].value_counts().sort_index().to_dict()))
        S.append(f"  core / extended: {(c.pool == 'core').sum()} / {(c.pool == 'extended').sum()}")
        S.append(f"  MPD (bp) 5/50/95%: {c.mpd_bp.quantile([.05, .5, .95]).round(0).tolist()}")
    six = comm[comm.plate6_origin != ""]
    S.append(f"\nplate {new_p} by origin: " + str(six.plate6_origin.value_counts().to_dict()))
    sh = comm[comm.passage_plates.str.contains(" ")]
    S.append(f"  on both passage plates: {len(sh)} communities = {len(sh) * cfg['replicates']} wells "
             f"per plate ({len(sh) * cfg['replicates'] / ((len(plates.inner_wells()) - cfg['blanks_per_plate']) // 3 * 3):.0%}); "
             f"sizes {sh['size'].value_counts().sort_index().to_dict()}")
    S.append(f"  plate {new_p} communities also on the main plates: "
             f"{six.blocks.str.contains('main').sum()}")
    inc = pd.Series([s for r in inoc.members for s in r.split()]).value_counts()
    S.append(f"\nwells per strain (all inoculated plates): min {inc.min()}, median {inc.median():.0f}, "
             f"max {inc.max()}; strains used {len(inc)} / {len(used)}")
    draw = plan.groupby("strain").volume_nL.sum() / 1000
    well = plan.groupby(["dest_plate", "dest_well"]).volume_nL.sum() / 1000
    S.append(f"inoculum: {cfg['inoculum']['member_nL']} nL per member, {cfg['inoculum']['reduced_nL']} nL "
             f"above {int(cfg['inoculum']['max_well_inoculum_uL'] * 1000 // cfg['inoculum']['member_nL'])} members")
    S.append(f"  per well: median {well.median():.2f} uL, max {well.max():.2f} uL")
    S.append(f"  drawn per source well: median {draw.median():.1f} uL, max {draw.max():.1f} uL "
             f"({draw.idxmax()}); source fills needed: {int(plan.fill.max())}"
             + ("  -> NO REFILL" if plan.fill.max() == 1 else "  -> refill between fills"))
    S.append(f"  Echo strain transfers: {len(plan):,}, droplets {int(plan.volume_nL.sum() / 25):,} "
             f"(20261021: 11,484 / 44,901)")
    S.append(f"plates: inoculated {sorted(int(p) for p in lay[lay.inoculated].plate.unique())}; "
             f"barcoded {sorted(int(p) for p in lay.plate.unique())} + decoy plate {cfg['decoy_barcode_plate']}")
    S.append("  barcode plate -> passage plate @ time point: "
             + ", ".join(f"{b}: {p}@{t}" for b, (p, t) in sorted(bmap.items())))
    S.append(f"barcodes excluded (failed in 20260907): fwd {sorted(failed['fwd'])}, rev {sorted(failed['rev'])}")
    S.append(f"barcode combinations: {len(lay) + len(plates.inner_wells()):,} used here, "
             f"{len(ref_combos):,} of {cfg['copy_from']['experiment']} excluded (decoys included); none shared")
    pv = prim.groupby("Source Well")["Transfer Volume"].sum() / 1000
    rv = (pd.concat([ref_lay.fwd_source_well, ref_lay.rev_source_well]).value_counts()
          * cfg["primers"]["volume_nL"] / 1000)
    both = pv.add(rv, fill_value=0)
    S.append(f"primer drawn per oligo source well: median {pv.median():.1f} uL, max {pv.max():.1f} uL "
             f"(with {cfg['copy_from']['experiment']} on the same primer plate: max {both.max():.1f} uL)")
    S.append(f"minibar: {mb_info}")
    return "\n".join(S) + "\n"


if __name__ == "__main__":
    main()
