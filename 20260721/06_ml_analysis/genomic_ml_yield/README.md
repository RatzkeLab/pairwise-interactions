# genomic_ml_yield — 20260721, naive join (expected null)

**Question:** the same yield model that works on 20260630, run on 20260721 with the naive
`Well_souce_plate` join — deliberately "forgetting" that the labels are wrong.

Run: `run_yield_naive.py`. Controls: 5 permuted genome mappings, both layout variants.

**Verdict: null.** `cv_strain` R² = −0.47, ρ = 0.108 against permuted nulls at −0.46 ± 0.18
(z = −0.06). The same pipeline gives R² = +0.265, ρ = 0.581 on 20260630.

Two things worth keeping from it:

- Only **2016 pairs / 64 strains** are available, because just 66 of 86 labels join to the
  genomic tables at all — far fewer than the 3828/88 that 20260630 yields.
- **The replicate ceiling works as a layout diagnostic.** 0.812 under the unswapped layout vs
  0.687 under `_plate1_2_swapped`, independently confirming that the unswapped layout (which
  matches `echo_strains`) describes the *physical* well contents, and that the swap is a
  barcode/sequencing correction. At 0.812 versus 20260630's 0.869, **this experiment was executed
  consistently** — only the identity mapping is broken.
