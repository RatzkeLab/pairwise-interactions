"""Allocation helpers for the 20260907 pairwise-interaction design.

The layout mechanics (barcode grid, Echo transfer files, minibar table) follow
20260721/01_setup/experiment_setup.ipynb. What changes here is WHICH wells get
which content:

* every strain needs a monoculture well, because the master plate itself cannot
  be assumed clean -- see merge_consensus_sequences/cross_contamination_check
* pairs are drawn only from partners whose 16S is far enough apart to be
  resolved in the readout (see separability.py)
* plates 1..N_PRIMARY must be self-sufficient, so that sequencing only those
  still yields every mono and the full replicate set

Pairs are generated as successive random MATCHINGS over the strain set rather
than by sampling pairs independently. One matching covers every strain exactly
once, so k rounds give every strain exactly k partners with no balancing step
and no strain left short.
"""
import numpy as np
from collections import defaultdict


def random_matching(strains, eligible, used, rng, unknown=None, allow_unknown=False):
    """One near-perfect matching over `strains`.

    Returns (pairs, unmatched). Greedy over a random order, which is enough here
    because the eligibility graph is dense (median 350 partners per strain).
    Preference order for each strain: an unused eligible partner, then -- only if
    `allow_unknown` -- an unused partner whose separability NO reference can
    assess, then leave it unmatched. A partner that a reference actively says is
    too close is never selected, and a repeated pair is never created silently.
    """
    order = list(strains)
    rng.shuffle(order)
    free = set(order)
    pairs, unmatched = [], []
    for a in order:
        if a not in free:
            continue
        free.discard(a)
        cands = [b for b in free if eligible(a, b) and frozenset((a, b)) not in used]
        if not cands and allow_unknown and unknown is not None:
            # unknown separability only -- never a pair a reference says is close
            cands = [b for b in free
                     if unknown(a, b) and frozenset((a, b)) not in used]
        if not cands:
            unmatched.append(a)
            continue
        b = cands[rng.integers(len(cands))]
        free.discard(b)
        pairs.append((a, b))
        used.add(frozenset((a, b)))
    return pairs, unmatched


def build_pairs(regular, limited, n_wells, eligible, rng, unknown=None, limited_rounds=2):
    """Fill `n_wells` with distinct pairs, balanced across strains.

    `limited` strains (never successfully ONT-sequenced despite six attempts)
    take part in only the first `limited_rounds` matchings: they are unlikely to
    yield readable sequence, so spending more wells on them costs coverage
    elsewhere for little expected return.
    """
    used, pairs, rounds = set(), [], []
    r = 0
    while len(pairs) < n_wells:
        pool = list(regular) + (list(limited) if r < limited_rounds else [])
        got, un = random_matching(pool, eligible, used, rng, unknown=unknown,
                                  allow_unknown=(r < limited_rounds))
        if not got:
            break
        for p in got:
            if len(pairs) < n_wells:
                pairs.append(p); rounds.append(r)
        r += 1
    return pairs, rounds


def spread_across_plates(items, plates, rng, key=None):
    """Assign items to plates so that members of a group land on distinct plates.

    Used for the technical-replicate sets (5 copies of one pair, which only test
    plate-to-plate reproducibility if they sit on 5 different plates) and for the
    two monocultures of a strain (so one plate failing cannot cost both).
    """
    out = []
    for group in items:
        n = len(group)
        chosen = list(rng.permutation(plates)[:n]) if n <= len(plates) else \
                 list(rng.permutation(plates)) + list(rng.choice(plates, n - len(plates)))
        out.append(list(zip(group, chosen)))
    return out
