"""Reconstruct what each plate-reader file in the two prep folders actually contains.

The lab notes give start times; the files are named by *save* time. Those differ by the whole
run duration (up to 44 h), so matching a note to a filename by its name is wrong. Every file's
internal `Date:`/`Time:` header is the RUN START, and the runs chain -- each starts within
minutes of the previous file being written -- which is what makes the reconstruction possible.

Plate identity is then established from the data rather than assumed: reads of the same
physical plate correlate at Spearman ~0.9-1.0 across the well grid, reads of different plates
at ~0-0.4. That is what separates, for example, the two "OD 384" reads on 22/07 (r=0.38 --
two different plates, not one plate read twice) from the 173208/173918 pair (r=0.97 -- one
plate, OD then spectrum).

Caveats kept in the output rather than smoothed over:
  - `ODFull_*` files are wavelength scans: their "cycles" are wavelengths, not timepoints, so
    correlating their first cycle against an OD600 read compares different wavelengths and
    understates agreement. They are paired to a plate by timing, not by correlation.
  - `TRF_1` / `FI_1` are fluorescence, not OD (values in the thousands); they carry no OD.
  - The 22-26/06 chain predates the earliest note supplied and is labelled by its own
    characteristics only.
  - The note dated "2026 06 09" is 3 weeks before the 30/06 preculture it feeds; a 29/06
    reading of the same activity fits the file record far better, and is flagged as such.

Writes `FILE_CONTENTS.csv` into each source folder.
"""

import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

PR = Path("/home/rl/scripts/karl/Link to Karl/plate_reader_csvs/data_ascii/Karl_2026")
FOLDERS = {"Karl_20260623_OD": "experiment 1 prep (-> 20260630 run)",
           "Karl_20260722_OD": "experiment 2 prep (-> 20260721 run)"}

# key: save-time stem -> (plate group, description). Plate groups from the correlation matrix.
LABELS = {
 # ---------------- Karl_20260623_OD ----------------
 "Karl_20260623_120022": ("pilot-A", "Growth curve, 384, shaking, 16.3 h (started 22/06 19:30). "
   "OD 0.20->0.77, so real growth. PREDATES the supplied notes; first of a 22-26/06 chain that "
   "looks like pilot growth-condition testing, not the documented prep."),
 "Karl_20260624_220526": ("pilot-B", "Growth curve, 384, shaking, 33.9 h (started 23/06 12:00). "
   "OD 0.10->0.46. Continues the pilot chain; only r=0.39 to the previous plate, so a different "
   "or re-inoculated plate."),
 "Karl_20260626_184629": ("pilot-C", "Growth curve, 384, shaking, 44.5 h (started 24/06 22:07). "
   "OD 0.20->0.19 -- essentially NO growth. Last of the pilot chain; likely a failed or "
   "no-carbon condition."),
 "Karl_20260629_191418": ("29/06-fresh", "Quick 2-cycle 384 check (started 29/06 19:03), mean OD "
   "0.08 -- freshly inoculated, no growth yet. Pairs with the 19:14 read (r=0.90, same plate). "
   "Timing fits the 'retrieve stock plates, inoculate NM preculture' note dated 2026-06-09, "
   "which is very likely 2026-06-29 -- 06-09 would leave a 3-week incubation before 30/06."),
 "Karl_20260629_192909": ("29/06-fresh", "Second quick 2-cycle 384 check (started 29/06 19:14), "
   "mean OD 0.14. Same plate as the 19:03 read (r=0.90). Sanity check right after inoculation."),
 "Karl_20260630_111648": ("4x96-Q1", "**96-well** read (started 30/06 11:10). Quadrant 1 of 4 "
   "read 15 min apart -- the NM precultures recovered from frozen stock, in the 4x96 working "
   "format the notes describe. Assigned quadrant A1 by chronological order."),
 "Karl_20260630_113201": ("4x96-Q2", "**96-well** read (started 30/06 11:25). Quadrant 2 of 4. "
   "Assigned A2. Uncorrelated with the other three (r -0.22..+0.12), as different wells should be."),
 "Karl_20260630_114305": ("4x96-Q3", "**96-well** read (started 30/06 11:39). Quadrant 3 of 4. "
   "Assigned B1. NOT a lid/no-lid repeat of the 11:49 read -- they share no per-well pattern "
   "(r=0.07 under every orientation), while cycles within one file repeat at r=0.998."),
 "Karl_20260630_115357": ("4x96-Q4", "**96-well** read (started 30/06 11:49). Quadrant 4 of 4. "
   "Assigned B2. Assembling these four as A1/A2/B1/B2 reproduces the 01/07 384 read at "
   "Spearman 0.56 (occupied wells), which validates both the quadrant order and the "
   "4x96 -> 384 consolidation."),
 "Karl_20260701_132952": ("scrapped-MM", "Single 384 read (started 01/07 13:26), mean OD 0.36. "
   "The FIRST, SCRAPPED experimental plate -- inoculated 30/06 09:30 into 0.01% glucose by "
   "mistake, realised 01/07 10:00. Same plate as the 13:42 run (r=0.98). Consolidated 4x96->384."),
 "Karl_20260701_171215": ("scrapped-MM", "3.4 h, 35-cycle run on the scrapped 0.01%-glucose "
   "plate (started 01/07 13:42). OD 0.35->0.36: flat, consistent with the wrong carbon source."),
 "Karl_20260701_171954": ("NM-recovery-2", "Single 384 read (started 01/07 17:13), mean OD 0.15 "
   "-- freshly inoculated. This is the RE-START: 4x96 stock retrieved from the freezer at 15:50 "
   "and inoculated into NM. Same plate as the 17:20 overnight run (r=0.96)."),
 "Karl_20260702_111117": ("NM-recovery-2", "The overnight sanity-check growth curve the notes "
   "put at 01/07 17:30 -- 178 cycles, 17.7 h, started 01/07 17:20. OD 0.15->0.71: the freezer "
   "recovery working. **This is the file that matches that note**, not one named 0701."),
 "Karl_20260702_133125": ("NM-recovery-2", "2.1 h, 22-cycle continuation of the same recovery "
   "plate (started 02/07 11:17), holding at OD 0.71. NOTE: despite its name this is NOT the "
   "13:30 NM recovery read -- see Karl_20260702_134335."),
 "Karl_20260702_134335": ("NM-dense", "**The 'read the NM OD to check which strains recovered "
   "from stock' read** -- started 02/07 13:31, matching the 13:30 note exactly. Mean OD 1.39, "
   "a dense overnight NM culture. A DIFFERENT plate from the 0.71 one (r=0.22)."),
 "Karl_20260702_182146": ("NM-dense", "4.5 h, 46-cycle run on that same dense NM plate (started "
   "02/07 13:44), OD 1.40->1.58. Same plate (r=0.99)."),
 "Karl_20260703_112451": ("NM-recovery-2", "16.9 h, 170-cycle overnight run (started 02/07 "
   "18:22) on the 0.71 recovery plate, holding 0.71->0.74 -- stationary."),
 "Karl_20260703_154632": ("multispec", "1.3 h, 14-cycle 384 read (started 03/07 14:22), mean OD "
   "0.17. Start of the plate carried into that night's reader-settings tests (r=0.90-0.92 with "
   "the 04/07 reads)."),
 "Karl_20260703_234757": ("multispec", "**Fluorescence (TRF), not OD** -- values in the "
   "thousands. Part of the 'test multispectrum reading settings' session noted at 03/07 23:00."),
 "Karl_20260703_235236": ("multispec", "**Fluorescence intensity (FI), not OD.** Reader-settings "
   "test. Contains no usable optical density."),
 "Karl_20260703_235350": ("multispec", "Full-spectrum wavelength scan (ODFull_1). Its 'cycles' "
   "are WAVELENGTHS, not timepoints. Reader-settings test."),
 "Karl_20260704_001122": ("multispec", "Full-spectrum scan, 50 flashes (ODFull_50). Settings "
   "test -- flash-count comparison."),
 "Karl_20260704_001637": ("multispec", "Full-spectrum scan, 10 flashes (ODFull_10). Settings "
   "test -- the 10-vs-50 flash comparison."),
 "Karl_20260704_002409": ("multispec", "Single OD 384 read (started 04/07 00:18) closing the "
   "settings session, mean OD 0.27."),
 # ---------------- Karl_20260722_OD ----------------
 "Karl_20260722_172409": ("exp2-preculture", "Single 384 read (started 22/07 17:18), mean OD "
   "**1.63** -- a DENSE grown culture. Best read as the MM+glucose preculture inoculated 21/07 "
   "10:30 (~31 h growth), consolidated 4x96->384. This is the closest analogue to experiment "
   "1's dense NM read."),
 "Karl_20260722_172802": ("exp2-preculture", "Full-spectrum scan of that same dense preculture "
   "plate (started 17:24). Cycles are WAVELENGTHS. Paired to the 17:18 plate by timing."),
 "Karl_20260722_173208": ("exp2-experiment", "Single 384 read (started 22/07 17:28), mean OD "
   "**0.46** -- a DIFFERENT plate from the 17:18 one (r=0.38). Diluted/freshly-set experiment "
   "plate: the 'plate for plate reader' prepared 22/07 09:00."),
 "Karl_20260722_173918": ("exp2-experiment", "Full-spectrum scan of that experiment plate "
   "(started 17:35). Same plate as the 17:28 read (r=0.97). Cycles are WAVELENGTHS."),
 "Karl_20260723_132335": ("exp2-experiment", "**The experiment growth curve** -- 197 cycles, "
   "19.6 h, no shaking, started 22/07 17:42 immediately after the plate was characterised. "
   "OD 0.10->0.36."),
 "Karl_20260724_105211": ("exp2-experiment", "Second 18.3 h growth run (started 23/07 16:24), "
   "OD 0.30->0.33 -- stationary. Continuation/repeat monitoring of the experiment plate."),
}


def meta(path):
    txt = Path(path).read_text(encoding="latin-1").splitlines()
    g = {}
    for l in txt[:12]:
        if l.startswith("Testname:"):
            g["test"] = l.split(":", 1)[1].strip()
        m = re.search(r"Date:\s*(\S+)\s+Time:\s*(\S+)", l)
        if m:
            g["run_started"] = f"{m.group(1)} {m.group(2)}"
        if "No. of Cycles" in l:
            g["cycles"] = l.split(":")[1].strip()
    times = [float(m.group(1)) for m in
             (re.match(r"Time \[s\]:\s*([\d.]+)", l) for l in txt) if m]
    g["duration_h"] = round((max(times) - min(times)) / 3600, 2) if times else 0.0

    def num(x):
        try:
            return float(x.strip())
        except ValueError:
            return np.nan
    rows = [[num(x) for x in l.split(",")] for l in txt if re.match(r"^\s*[\d.\-\s]+,", l)]
    rows = [r for r in rows if len(r) in (12, 24)]
    if rows:
        a = np.array(rows)
        g["plate_format"] = {12: "96-well", 24: "384-well"}[a.shape[1]]
        h = 8 if a.shape[1] == 12 else 16
        g["mean_OD_first_cycle"] = round(float(np.nanmean(a[:h])), 3)
        g["mean_OD_last_cycle"] = round(float(np.nanmean(a[-h:])), 3)
    return g


def main():
    for folder, ctx in FOLDERS.items():
        rows = []
        for f in sorted(glob.glob(str(PR / folder / "*.csv"))):
            name = os.path.basename(f)
            stem = name[:20].rstrip("_")
            stem = re.match(r"(Karl_\d{8}_\d{6})", name).group(1)
            grp, desc = LABELS.get(stem, ("?", "UNLABELLED - not matched to any note"))
            rows.append({"filename": name,
                         "what_this_file_most_likely_contains": desc,
                         "plate_group": grp, **meta(f)})
        df = pd.DataFrame(rows)
        cols = ["filename", "what_this_file_most_likely_contains", "plate_group", "test",
                "run_started", "cycles", "duration_h", "plate_format",
                "mean_OD_first_cycle", "mean_OD_last_cycle"]
        df = df[[c for c in cols if c in df.columns]]
        out = PR / folder / "FILE_CONTENTS.csv"
        df.to_csv(out, index=False)
        print(f"{folder}  ({ctx})  ->  {out.name}   [{len(df)} files]")
        print(df[["filename", "plate_group", "run_started", "plate_format",
                  "mean_OD_first_cycle"]].to_string(index=False))
        print()


if __name__ == "__main__":
    main()
