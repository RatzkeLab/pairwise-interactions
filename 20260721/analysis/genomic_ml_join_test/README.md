# genomic_ml_join_test — is 20260721's genomic join real?

**Question:** train on 20260721's interaction labels with genomes attached naively via
`Well_souce_plate`, and ask whether it beats chance.

Run: `run_join_test.py [--quick] [--refix]` (`--refix` uses the corrected labels).

**The trap this is built around:** a *consistently wrong* mapping is still a bijection, so each
KO vector works as a unique per-strain fingerprint. Under `cv_pair` a model then scores well
while knowing no biology — 20260721 posts R² = 0.44 there. Only `cv_strain`, plus a
**permuted-mapping** null (shuffle which genome attaches to which well), can tell a correct join
from a scrambled one.

| | 20260721 | 20260630 (positive control) |
|---|---|---|
| `cv_pair` z vs permuted | −0.03 | **+3.5** |
| `cv_strain` z vs permuted | +0.69 | **+4.1** |

**Verdict: null.** The true mapping sits on the permutation mean. The control experiment passes
decisively through the same code, so this is not the method failing.

Unchanged when re-run on the corrected labels (`*_refix` outputs, 1116 pairs instead of 874).
