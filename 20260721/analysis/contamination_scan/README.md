# contamination_scan — is 20260721 contaminated?

**Question:** could pre-PCR DNA carryover explain the label mismatch?

Run: `scan.py` (runs 20260630 as the control in the same pass).

**Verdict: no.** Excess-organism rate in pair wells **18.4% vs 17.8%** in the control, and
20260721's mono wells are *cleaner* than the control's (23.1% vs 50.0%).

Two design points that reverse the reading if skipped:

- **Split by well type.** Most wells are pair wells and hold two organisms *by design*; a rate
  over all wells measures the experiment, not contamination (TRAPS §10).
- **Work on 16S groups.** Only 35 of 294 strains have a unique 16S; a read best-hitting its 1 bp
  twin is the resolution limit, not contamination. 87 db entries collapse to 73 groups.

The one real difference — 20260721 has ~2× the low-identity reads — is **not** a foreign
organism: sampled low-identity reads are diffuse (median pairwise identity 0.81, <0.5% of pairs
above 0.95) and full length. It is consistent with wells holding organisms from outside the
87-entry reference, which is exactly what FINDINGS §1 concluded.
