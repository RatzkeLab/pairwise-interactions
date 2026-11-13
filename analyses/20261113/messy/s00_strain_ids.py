"""Well of the 384-well collection plate -> the strain's number in the lab collection.

Every earlier analysis keyed strains by their well (A10, K14, ...). From here on outputs carry
`strain_no`, the ISOLATE's number in the lab collection. Wells stay alongside as `well`.

Three sources, none complete on its own:
  384 sheet      merge_consensus_sequences/plate_format/Full_384_strain_collection_no_seq_info.csv
                 well -> `Strain.Number`, which is a freezer STOCK, not an isolate: B14 = M7 and
                 K19 = "MYb9" are two stocks of the same isolate MYb9, so the column is unique
                 per well even though the strain is not.
  master sheet   plate_format/Strain collection including Schulenberg's - collection.csv
                 stock number -> alternative name, box, position (464 stocks). MYb9 has two
                 stocks (M7 in WS5/7, and a row numbered "MYb9" in WS4/56); MYb10 has three
                 (M8, M135, "MYb10").
  mapping        data_links/resources/mapping_384_well_plate_collection.csv
                 well -> genome, for the 298 wells with a genome; its `strain_index_name
                 (as_in_boxes)` already gives the isolate number (K19 -> M7).

Rule: isolate = the master sheet's alternative name for the stock if it has one (M7 -> MYb9),
else the stock number itself ("MYb9" -> MYb9, 8 -> 8). strain_no = the first M-number that the
master sheet gives that isolate (MYb9 -> M7, MYb10 -> M8 for B16, K2 and M20 alike); isolates
with no M-number keep their stock number (all Or isolates: integers). Checked against the
mapping's box ids below -- they must agree on every one of its 298 wells.

Output: outputs/s00_strain_ids.csv
  well, strain_no, isolate_name, stock_no, box, box_position, n_wells_same_isolate,
  genome_table_id, assembly_name
"""
import re

import pandas as pd

import config

PLATE_FORMAT = config.COLLECTION_SHEET.parent
c = pd.read_csv(config.COLLECTION_SHEET, dtype=str)
master = pd.read_csv(PLATE_FORMAT / "Strain collection including Schulenberg's - collection.csv", dtype=str)
master = master.rename(columns={"Strain number": "stock_no", "Alternative strain name": "alt_name",
                                "Box": "box", "Position in Box": "box_position"})
master = master[master.stock_no.notna()]
master["stock_no"] = master.stock_no.str.strip()
master["alt_name"] = master.alt_name.str.strip()
assert master.stock_no.is_unique
mrow = master.set_index("stock_no")


def m_number(s):
    m = re.fullmatch(r"M(\d+)", s)
    return int(m.group(1)) if m else None


ids = pd.DataFrame({"well": c["Well.position"], "stock_no": c["Strain.Number"].str.strip()})
missing = sorted(set(ids.stock_no) - set(mrow.index))
assert not missing, f"stock numbers absent from the master sheet: {missing}"
ids["isolate_name"] = [mrow.alt_name.get(s) if isinstance(mrow.alt_name.get(s), str) else s
                       for s in ids.stock_no]
first_m = (master[master.alt_name.notna() & master.stock_no.map(m_number).notna()]
           .assign(n=lambda d: d.stock_no.map(m_number)).sort_values("n")
           .drop_duplicates("alt_name").set_index("alt_name").stock_no)
ids["strain_no"] = [first_m.get(i, s) for i, s in zip(ids.isolate_name, ids.stock_no)]
ids["box"] = ids.stock_no.map(mrow.box)
ids["box_position"] = ids.stock_no.map(mrow.box_position)
ids["n_wells_same_isolate"] = ids.groupby("strain_no").well.transform("size")

m = pd.read_csv(config.STRAIN_MAPPING, dtype=str).rename(columns={
    "Well_souce_plate": "well", "strain": "genome_table_id", "strain_index_name(as_in_boxes)": "box_id"})
ids = ids.merge(m[["well", "box_id", "genome_table_id", "assembly_name"]], on="well", how="left")
chk = ids[ids.box_id.notna()]
bad = chk[chk.box_id != chk.strain_no]
assert bad.empty, f"disagrees with the mapping's box ids:\n{bad.to_string()}"
ids = ids.drop(columns="box_id")
assert ids.well.is_unique and len(ids) == 384
cols = ["well", "strain_no", "isolate_name", "stock_no", "box", "box_position", "n_wells_same_isolate",
        "genome_table_id", "assembly_name"]
ids[cols].to_csv(config.OUT / "s00_strain_ids.csv", index=False)

print(f"{len(ids)} wells -> {ids.strain_no.nunique()} isolates "
      f"({ids.stock_no.nunique()} freezer stocks); agrees with the mapping on all {len(chk)} genome wells")
dup = ids[ids.n_wells_same_isolate > 1].sort_values(["strain_no", "well"])
print(f"\n{dup.strain_no.nunique()} isolates in more than one well:")
print(dup.groupby(["strain_no", "isolate_name"]).apply(
    lambda g: ", ".join(f"{w} ({s}{', genome' if isinstance(gt, str) else ''})"
                        for w, s, gt in zip(g.well, g.stock_no, g.genome_table_id)),
    include_groups=False).to_string())


def lookup():
    """well -> strain_no, for relabelling any table produced in this analysis."""
    t = pd.read_csv(config.OUT / "s00_strain_ids.csv", dtype=str)
    return dict(zip(t.well, t.strain_no))
