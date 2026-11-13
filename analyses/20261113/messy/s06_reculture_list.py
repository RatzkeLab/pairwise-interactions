"""Which strains to re-culture from an older stock: every well whose genome on file is NOT
the organism in the 20260907 source plate (20261021 concordance, all 378 wells with a
reference -- not only the ones that were modelled), by strain number and freezer position.

  status_20260907   20261021's verdict for the September source plate (the current plate)
  status_20260630   s01's verdict from July's references, where July had one -- shows whether
                    the well was already wrong in July (then the problem predates both runs)
  now_holds         the genome unit(s) that DO carry the well's 16S, i.e. the likely invader
                    ("C001" = the G01 clone; "no genome on file" = invader not sequenced)

Well names inside `now_holds` are translated to strain numbers; C### ids are genome clone
clusters, whose member strains are listed in s06_clone_clusters.csv.

Outputs: outputs/s06_reculture_list.csv   (sorted by box and position, for pulling from the freezer)
         outputs/s06_clone_clusters.csv   (clone cluster -> member strain numbers)
"""
import re

import pandas as pd

import config

g = pd.read_csv(config.CONC_OUT / "g02_strain_genome_status.csv")
ids = pd.read_csv(config.OUT / "s00_strain_ids.csv", dtype=str)
v30 = pd.read_csv(config.OUT / "s01_genome_status_20260630.csv").set_index("strain")

w = g[g.genome_status.isin(["wrong", "unverified"])].rename(columns={"strain": "well"})
w = w.merge(ids, on="well", how="left")
w["status_20260630"] = w.well.map(v30.genome_status).fillna("no July reference")
w = w.rename(columns={"genome_status": "status_20260907", "genome_is_of": "now_holds"})
to_no = dict(zip(ids.well, ids.strain_no))
w["now_holds"] = w.now_holds.map(lambda v: v if not isinstance(v, str) else " ".join(
    to_no.get(t, t) if re.fullmatch(r"[A-P]\d{1,2}", t) else t for t in v.split())
    if v != "no genome on file" else v)

cl = pd.read_csv(config.CONC_OUT / "g02_genome_clusters.csv").dropna(subset=["clone_cluster"])
cl = (cl.assign(strain_no=cl.well.map(to_no)).groupby("clone_cluster")
      .agg(n_wells=("well", "size"), strain_nos=("strain_no", lambda x: " ".join(sorted(x, key=lambda v: (len(v), v)))),
           wells=("well", " ".join)).reset_index())
cl.to_csv(config.OUT / "s06_clone_clusters.csv", index=False)
w["box_position"] = pd.to_numeric(w.box_position, errors="coerce")
cols = ["strain_no", "isolate_name", "well", "box", "box_position", "status_20260907", "status_20260630",
        "now_holds", "evidence", "tier", "genome_table_id"]
w = w[cols].sort_values(["status_20260907", "box", "box_position"], ascending=[False, True, True])
w.to_csv(config.OUT / "s06_reculture_list.csv", index=False)

print(pd.crosstab(w.status_20260907, w.status_20260630).to_string())
wr = w[w.status_20260907 == "wrong"]
print(f"\n{len(wr)} wrong-genome wells (= {wr.strain_no.nunique()} strains) to re-culture; "
      f"{(w.status_20260907 == 'unverified').sum()} more unverified")
print("most common invaders:", wr.now_holds.value_counts().head(6).to_dict())
print(wr.head(12).to_string(index=False))
