"""Precompute the all-vs-all 16S identity matrix for the source collection (cached)."""
import numpy as np, pandas as pd, qc_config as C, qc_sources as S, qc_compare as Q

g = S.genome_16S_by_well()
reps = {}
for r in g.itertuples():
    reps.setdefault(r.strain_label, []).append(r.seq)
strains = sorted(reps)
n = len(strains)
M = np.eye(n)
for i in range(n):
    for j in range(i + 1, n):
        M[i, j] = M[j, i] = Q.best_identity(reps[strains[i]], reps[strains[j]])[0]
pd.DataFrame(M, index=strains, columns=strains).to_csv(C.OUT / "s07_collection_16S_identity_matrix.csv")
print("done", M.shape)
