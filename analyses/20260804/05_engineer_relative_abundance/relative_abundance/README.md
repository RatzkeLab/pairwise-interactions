# relative_abundance — 20260721 (PRE-FIX — see `../relative_abundance_refix/`)

Same pipeline as 20260630's. **These outputs predate the `MIN_RESOLVABLE_BP = 10` fix of
2026-08-27** — `r03_pair_replicate_stats.csv` here has no `below_resolution_limit` column.

Corrected outputs are in **`../relative_abundance_refix/outputs/`**: high-uncertainty pairs fall
494 → 31, usable pairs rise 1625 → 2088, and **~48% of the labels move** (median |Δlog2| 0.47).

This folder is kept because earlier analyses consumed it and the diff between the two is itself
documented (`../relative_abundance_refix/outputs/cmp01_old_vs_new_r03.csv`).

Separately: this experiment's **strain labels do not name the organisms in its wells**
(FINDINGS §1). The interaction data is internally excellent — a per-strain strength model
explains it at R² 0.74 — but it cannot be joined to genomes or merged with 20260630.
