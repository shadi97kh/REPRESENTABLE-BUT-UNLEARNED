"""Independent enumerators.

INDEPENDENCE IS THE POINT. Nothing here imports, calls, or reimplements the dynamic
programme. Members are produced from the *definition* of the family — a predicate applied
to candidate sets — and the quantities are then computed by iterating over the explicit
list. If these agree with the recurrence, the recurrence is doing what it claims; if the
grammar were ambiguous, `count` would exceed `len(enumerate_*)` and the tests would fail.

This module deliberately does not memoise over intervals, does not decompose spans, and
does not use a semiring.
"""
from __future__ import annotations

import math

CANONICAL = frozenset({("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"),
                       ("G", "U"), ("U", "G")})


# ------------------------------------------------------------------ RNA, by definition
def rna_candidates(seq, min_loop=3, canonical_only=True):
    seq = seq.upper().replace("T", "U")
    out = []
    for i in range(len(seq)):
        for j in range(i + min_loop + 1, len(seq)):
            if canonical_only and (seq[i], seq[j]) not in CANONICAL:
                continue
            out.append((i, j))
    return out


def rna_is_valid(pairs, seq, min_loop=3, canonical_only=True):
    """The definition, checked directly: bounds, min-loop, canonicity, one partner per
    index, and no crossing. No recurrence anywhere."""
    seq = seq.upper().replace("T", "U")
    n = len(seq)
    ps = sorted(tuple(p) for p in pairs)
    used = set()
    for (i, j) in ps:
        if not (0 <= i < j < n):
            return False
        if j - i <= min_loop:
            return False
        if canonical_only and (seq[i], seq[j]) not in CANONICAL:
            return False
        if i in used or j in used:
            return False
        used.add(i); used.add(j)
    for a in range(len(ps)):
        i, j = ps[a]
        for b in range(a + 1, len(ps)):
            k, l = ps[b]
            if i < k < j < l or k < i < l < j:
                return False
    return True


def _compatible(new, chosen):
    (i, j) = new
    for (k, l) in chosen:
        if i == k or i == l or j == k or j == l:
            return False
        if i < k < j < l or k < i < l < j:
            return False
    return True


def enumerate_rna(seq, min_loop=3, canonical_only=True, forced=(), forbidden=()):
    """Every valid structure satisfying the mask. Backtracking over the candidate list
    with the pairwise-compatibility definition; no interval decomposition."""
    forced = [tuple(e) for e in forced]
    forbidden = {tuple(e) for e in forbidden}
    cands = [c for c in rna_candidates(seq, min_loop, canonical_only)
             if c not in forbidden]
    # forced pairs must themselves be admissible and mutually compatible
    for e in forced:
        if e in forbidden or not rna_is_valid([e], seq, min_loop, canonical_only):
            return []
    base = []
    for e in forced:
        if not _compatible(e, base):
            return []
        base.append(e)
    free = [c for c in cands if c not in set(forced)]

    out = []

    def rec(start, chosen):
        out.append(tuple(sorted(chosen)))
        for idx in range(start, len(free)):
            c = free[idx]
            if _compatible(c, chosen):
                chosen.append(c)
                rec(idx + 1, chosen)
                chosen.pop()

    rec(0, list(base))
    # a final definitional pass, so nothing reaches the tests unvalidated
    return [s for s in out if rna_is_valid(s, seq, min_loop, canonical_only)]


# ------------------------------------------------------- layered DAG paths, by definition
def enumerate_dag_paths(layer_sizes, edges, forced=(), forbidden=()):
    """Every source-to-sink path, built by extending one layer at a time and checking
    that consecutive edges meet. No semiring, no accumulated chart."""
    forced = {tuple(e) for e in forced}
    forbidden = {tuple(e) for e in forbidden}
    if forced & forbidden:
        return []
    T = len(layer_sizes) - 1
    by_layer = {l: [] for l in range(T)}
    for e in edges:
        by_layer[e[0]].append(tuple(e))
    forced_at = {}
    for e in forced:
        if e[0] in forced_at and forced_at[e[0]] != e:
            return []
        forced_at[e[0]] = e

    out = []

    def rec(l, node, acc):
        if l == T:
            out.append(tuple(acc))
            return
        for e in by_layer[l]:
            if e in forbidden:
                continue
            if l in forced_at and e != forced_at[l]:
                continue
            if e[1] != node:
                continue
            acc.append(e)
            rec(l + 1, e[2], acc)
            acc.pop()

    for start in range(layer_sizes[0]):
        rec(0, start, [])
    return out


# ------------------------------------------------------------- quantities, by iteration
def brute_count(members):
    return len(members)


def brute_support(members, coefficients):
    """(max score, witness). Returns (-inf, None) for an empty family."""
    if not members:
        return float("-inf"), None
    best, arg = None, None
    for m in members:
        s = sum(coefficients.get(e, 0) for e in m)
        if best is None or s > best:
            best, arg = s, m
    return best, arg


def brute_log_partition(members, coefficients):
    if not members:
        return float("-inf")
    scores = [sum(coefficients.get(e, 0) for e in m) for m in members]
    hi = max(scores)
    return hi + math.log(sum(math.exp(s - hi) for s in scores))


def brute_marginals(members, coefficients, ground):
    if not members:
        return {e: 0.0 for e in ground}
    scores = [sum(coefficients.get(e, 0) for e in m) for m in members]
    hi = max(scores)
    w = [math.exp(s - hi) for s in scores]
    Z = sum(w)
    out = {}
    for e in ground:
        out[e] = sum(wi for wi, m in zip(w, members) if e in m) / Z
    return out
