"""Community generators.

All of them share one `Builder`, which holds the strain pool, the 16S distance matrix,
and how often each strain has been used so far -- so every generator draws toward the
least-used strains and the design stays balanced across modules without a separate
balancing step.

Relatedness is the mean pairwise 16S distance between members (MPD, haplotype-aware bp).
In this collection most strains are far apart (median pair ~300 bp), so random draws all
land at MPD ~220-320 bp; related communities have to be BUILT. `walk` does it: the next
member is drawn with weight exp(-(mean distance to current members)/tau) among strains
still >= min_sep from every member. Small tau stays inside a clade, tau = inf is uniform.
Candidates from a range of tau are binned into per-size MPD terciles, and each tercile
gets an equal quota -- so every size spans the whole relatedness range it can reach.
"""
import itertools

import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

USAGE_ALPHA = 1.0          # draw weight ~ 1 / (1 + times used) ** alpha
TAUS = (5.0, 15.0, 40.0, 100.0, np.inf)


class Builder:
    def __init__(self, pool, D, rng):
        self.pool = pool.set_index("strain")
        self.D = D
        self.names = list(D.index)
        self.idx = {n: i for i, n in enumerate(self.names)}
        self.Dv = D.values
        self.rng = rng
        self.usage = np.zeros(len(self.names))
        self.records = {}          # frozenset(members) -> record
        self.order = []

    # ---- bookkeeping -------------------------------------------------------------------
    def strains(self, pool="core"):
        p = self.pool
        if pool == "core":
            return list(p.index[p["pool"] == "core"])
        if pool == "extended_only":
            return list(p.index[p["pool"] == "extended"])
        return list(p.index)                     # "extended": core + extended

    def add(self, members, module, **kw):
        """Register a community; a set already present gets the new module tag added."""
        key = frozenset(members)
        if key in self.records:
            rec = self.records[key]
            rec["modules"].append(module)
            return rec
        rec = dict(community_id=f"K{len(self.order) + 1:04d}", modules=[module],
                   members=tuple(sorted(members)), size=len(members), **kw)
        self.records[key] = rec
        self.order.append(rec)
        for m in members:
            self.usage[self.idx[m]] += 1
        return rec

    def table(self):
        rows = []
        for r in self.order:
            sub = self.Dv[np.ix_([self.idx[m] for m in r["members"]],
                                 [self.idx[m] for m in r["members"]])]
            iu = np.triu_indices(len(r["members"]), 1)
            d = sub[iu] if len(r["members"]) > 1 else np.array([np.nan])
            pools = self.pool.loc[list(r["members"]), "pool"]
            rows.append({**{k: v for k, v in r.items() if k not in ("members", "modules")},
                         "module": r["modules"][0], "also_in": ";".join(r["modules"][1:]),
                         "members": " ".join(r["members"]),
                         "min_sep_bp": float(np.nanmin(d)), "mpd_bp": float(np.nanmean(d)),
                         "n_extended": int((pools == "extended").sum())})
        df = pd.DataFrame(rows)
        df["pool"] = np.where(df.n_extended > 0, "extended", "core")
        return df

    # ---- generators -------------------------------------------------------------------
    def walk(self, k, candidates, min_sep, tau, start=(), seed_from=None, tries=200):
        """Grow a community of size k; every pair >= min_sep. None if it cannot."""
        cand = np.array([self.idx[c] for c in candidates])
        u = self.usage
        for _ in range(tries):
            members = [self.idx[s] for s in start]
            if not members:
                pool0 = (np.array([self.idx[c] for c in seed_from]) if seed_from else cand)
                w = 1.0 / (1.0 + u[pool0]) ** USAGE_ALPHA
                members = [int(self.rng.choice(pool0, p=w / w.sum()))]
            while len(members) < k:
                ok = np.all(self.Dv[np.ix_(cand, members)] >= min_sep, axis=1)
                avail = cand[ok]
                if not len(avail):
                    break
                score = self.Dv[np.ix_(avail, members)].mean(axis=1)
                wt = np.ones(len(avail)) if np.isinf(tau) else np.exp(-(score - score.min()) / tau)
                wt = wt / (1.0 + u[avail]) ** USAGE_ALPHA
                members.append(int(self.rng.choice(avail, p=wt / wt.sum())))
            if len(members) == k:
                return [self.names[i] for i in members]
        return None

    def mpd(self, members):
        ix = [self.idx[m] for m in members]
        sub = self.Dv[np.ix_(ix, ix)]
        return float(sub[np.triu_indices(len(ix), 1)].mean())

    def stratified(self, counts, pool, min_sep, module, strata=3,
                   n_candidates=150, jaccard_max=0.6, exclude=()):
        """counts: {size: n}. n communities of each size, split evenly across MPD terciles,
        least-used strains first. Picks round-robin over (size, stratum) so usage stays
        balanced across sizes rather than being exhausted by whichever size runs first."""
        cands_all = self.strains("core" if pool == "core" else "extended")
        cands_all = [c for c in cands_all if c not in set(exclude)]
        seed_from = ([c for c in self.strains("extended_only") if c not in set(exclude)]
                     if pool == "extended" else None)
        bins, quota = {}, {}
        sizes = [k for k, n in counts.items() if n > 0]
        for k in sizes:
            pool_k = {}
            for tau in TAUS:
                for _ in range(n_candidates):
                    m = self.walk(k, cands_all, min_sep, tau, seed_from=seed_from)
                    if m:
                        pool_k[frozenset(m)] = m
            cl = [m for key, m in pool_k.items() if key not in self.records]
            if not cl:
                continue
            mp = np.array([self.mpd(m) for m in cl])
            edges = np.quantile(mp, np.linspace(0, 1, strata + 1)[1:-1])
            st = np.digitize(mp, edges)
            base, extra = divmod(counts[k], strata)
            order = self.rng.permutation(strata)
            for j, s in enumerate(order):
                bins[(k, s)] = [m for m, x in zip(cl, st) if x == s]
                quota[(k, s)] = base + (1 if j < extra else 0)
        chosen = {k: [] for k in sizes}
        progress = True
        while progress:
            progress = False
            for key in sorted(bins, key=lambda x: (x[0], x[1])):
                if quota[key] <= 0:
                    continue
                k, s = key
                best, best_score = None, None
                for m in bins[key]:
                    fs = frozenset(m)
                    if fs in self.records:
                        continue
                    if k >= 4 and any(len(fs & c) / len(fs | c) >= jaccard_max for c in chosen[k]):
                        continue
                    score = self.usage[[self.idx[x] for x in m]].sum() + self.rng.random()
                    if best is None or score < best_score:
                        best, best_score = m, score
                if best is None:
                    quota[key] = 0
                    continue
                self.add(best, module, relatedness_stratum=int(s), stratum_of=strata)
                chosen[k].append(frozenset(best))
                quota[key] -= 1
                progress = True
        return [self.records[c] for k in sizes for c in chosen[k]]

    def max_separable(self, names, min_sep, jitter=0.0):
        """Largest set with every pair >= min_sep (exact ILP). `jitter` > 0 breaks ties at
        random, so repeated calls return different maximum sets."""
        idx = {n: i for i, n in enumerate(names)}
        sub = self.D.loc[names, names].values
        e = [(i, j) for i in range(len(names)) for j in range(i + 1, len(names)) if sub[i, j] < min_sep]
        A = lil_matrix((max(len(e), 1), len(names)))
        for r, (i, j) in enumerate(e):
            A[r, i] = 1; A[r, j] = 1
        c = -(1 + jitter * self.rng.random(len(names)))
        res = milp(c=c, constraints=[LinearConstraint(A.tocsr(), -np.inf, 1)] if e else [],
                   integrality=np.ones(len(names)), bounds=Bounds(0, 1),
                   options={"time_limit": 120})
        return [n for n, v in zip(names, res.x) if v > 0.5]

    def large(self, specs, module="large_test", cap=None):
        """specs: [{"size": int | "max", "pool": core/extended, "min_sep": bp}]. A "max"
        set larger than `cap` is randomly trimmed to `cap` (inoculum volume limit)."""
        out = []
        for sp in specs:
            names = self.strains(sp["pool"])
            if sp["size"] == "max":
                m = self.max_separable(names, sp["min_sep"])
                if cap and len(m) > cap:
                    m = list(self.rng.choice(m, cap, replace=False))
            else:
                base = self.max_separable(names, sp["min_sep"], jitter=0.01)
                if sp["pool"] == "extended":
                    ext = [x for x in base if self.pool.at[x, "pool"] == "extended"]
                    keep = list(self.rng.choice(ext, 1)) if ext else []
                    rest = [x for x in base if x not in keep]
                    m = keep + list(self.rng.choice(rest, sp["size"] - len(keep), replace=False))
                else:
                    m = list(self.rng.choice(base, sp["size"], replace=False))
            out.append(self.add(m, module, target_min_sep=sp["min_sep"]))
        return out

    def resolution_tests(self, bands, n_per_band, size, pool, min_sep_rest, module="resolution_test"):
        """Small communities built around ONE deliberately close pair (in `band` bp), all
        other pairs >= min_sep_rest: does community-level read assignment hold up at the
        separations the calibration says should work?"""
        names = self.strains(pool)
        out, used = [], set()
        for lo, hi in bands:
            pairs = [(a, b) for a, b in itertools.combinations(names, 2)
                     if lo <= self.D.at[a, b] < hi]
            self.rng.shuffle(pairs)
            pairs.sort(key=lambda p: self.usage[self.idx[p[0]]] + self.usage[self.idx[p[1]]])
            n = 0
            for a, b in pairs:
                if n >= n_per_band or a in used or b in used:
                    continue
                rest = [c for c in names if c not in (a, b)
                        and self.D.at[a, c] >= min_sep_rest and self.D.at[b, c] >= min_sep_rest]
                m = self._extend((a, b), size, rest, min_sep_rest)
                if m is None:
                    continue
                out.append(self.add(m, module, close_pair=f"{a}-{b}",
                                    close_pair_bp=float(self.D.at[a, b])))
                used.update((a, b)); n += 1
        return out

    def _extend(self, start, k, candidates, min_sep):
        members = list(start)
        cand = [c for c in candidates if c not in members]
        for _ in range(200):
            m = list(members)
            pool = list(cand)
            while len(m) < k:
                ok = [c for c in pool if all(self.D.at[c, x] >= min_sep for x in m if x not in start)
                      and all(self.D.at[c, s] >= min_sep for s in start)]
                if not ok:
                    break
                w = np.array([1.0 / (1.0 + self.usage[self.idx[c]]) for c in ok])
                pick = ok[int(self.rng.choice(len(ok), p=w / w.sum()))]
                m.append(pick); pool.remove(pick)
            if len(m) == k:
                return m
        return None

    def trios_with_pairs(self, n, pool, min_sep, measured=frozenset(), module="trio"):
        """`n` trios, each with all three of its pairs as communities of their own: the
        direct test of whether pairwise outcomes predict a 3-member community. Trios whose
        pairs were already measured (in 20260907) are preferred, so the pairs double as a
        cross-experiment check."""
        names = self.strains(pool)
        cands = {}
        for tau in TAUS:
            for _ in range(200):
                m = self.walk(3, names, min_sep, tau)
                if m:
                    cands[frozenset(m)] = m
        cl = list(cands.values())
        mp = np.array([self.mpd(m) for m in cl])
        edges = np.quantile(mp, np.linspace(0, 1, n + 1))
        out = []
        for j in range(n):                       # one trio per MPD slice: spans relatedness
            sl = [m for m, x in zip(cl, mp) if edges[j] <= x <= edges[j + 1]
                  and frozenset(m) not in self.records]
            if not sl:
                continue
            sl.sort(key=lambda m: (-sum(frozenset(p) in measured for p in itertools.combinations(m, 2)),
                                   self.usage[[self.idx[x] for x in m]].sum() + self.rng.random()))
            trio = self.add(sl[0], module)
            trio_id = trio["community_id"]
            for p in itertools.combinations(sl[0], 2):
                self.add(list(p), "trio_pair", parent=trio_id,
                         measured_in_20260907=frozenset(p) in measured)
            out.append(trio)
        return out

    def leave_out_series(self, base_size, depth, pool, min_sep, tau=40.0, module="leave_out"):
        """A base community plus every version with 1..depth members removed (down to 2)."""
        names = self.strains(pool)
        base = self.walk(base_size, names, min_sep, tau)
        series = f"S{base_size}d{depth}"
        out = [self.add(base, module, series=series, removed="")]
        for d in range(1, depth + 1):
            for rem in itertools.combinations(base, d):
                m = [x for x in base if x not in rem]
                if len(m) >= 2:
                    out.append(self.add(m, module, series=series, removed=" ".join(rem)))
        return out

    def known_pairs(self, pairs, n, module="sanity_pair"):
        """Re-measure pairs from a previous run: a check that the source plate still holds
        the same strains, and a cross-experiment replicate of the interaction."""
        out, used = [], set()
        for a, b in pairs:
            if len(out) >= n:
                break
            if a in used or b in used:
                continue
            out.append(self.add([a, b], module)); used.update((a, b))
        return out
