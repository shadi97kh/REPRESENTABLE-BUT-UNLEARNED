"""Non-crossing RNA secondary structures as a conditionable support oracle.

FAMILY. Given a sequence of length n, a member is a set of pairs (i, j), i < j, with
  * j - i > min_loop            (hairpins need room)
  * (s_i, s_j) canonical         (optional; disabled gives the combinatorial family)
  * every index in at most one pair
  * no crossing: never i < k < j < l for two pairs (i, j), (k, l)

RECURRENCE — leftmost-anchored, unambiguous. Over the interval [i, j]:

    W[i][j] =  W[i+1][j]                                     ... i is unpaired
             + sum_k  w(i,k) * W[i+1][k-1] * W[k+1][j]       ... i pairs with k

    W[i][j] =  1  whenever i > j                             ... the empty structure

Unambiguity is by case on the *leftmost* index i: in any structure i is either unpaired
or paired to exactly one k, and those cases are disjoint and exhaustive. Each structure
therefore has exactly one derivation, which is what makes `count()` equal to |F| rather
than to a number of derivations. `tests/oracles` checks that against brute force, which
is the only way to be sure the grammar is not quietly ambiguous.

PRODUCTION-TO-EDGE MAP P. `unpaired(i)` introduces no edge; `pair(i,k)` introduces
exactly the edge (i, k). So P is injective on edge-introducing productions and the edge
multiset of a derivation is exactly the structure it derives — both are asserted in the
tests rather than assumed here.

COMPLEXITY. Chart values O(n^2), materialised productions O(n^3). Both are counted at
run time and reported in `stats`, so the claim is measured rather than asserted.
"""
from __future__ import annotations

import math

from .base import (COUNTING, EXACT_INTEGER, FLOAT_DIAGNOSTIC, FLOAT_UNVERIFIED,
                   INFEASIBLE, LOG_SUM_EXP, MAX_PLUS, NEG_INF, PROVED, Feasibility,
                   OracleResult, Production, SupportOracle, family_hash, normalise_mask)

CANONICAL = frozenset({("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"),
                       ("G", "U"), ("U", "G")})


def crosses(a, b) -> bool:
    (i, j), (k, l) = a, b
    return (i < k < j < l) or (k < i < l < j)


class RNANonCrossingOracle(SupportOracle):
    family_id = "rna_noncrossing"

    def __init__(self, seq: str, min_loop: int = 3, canonical_only: bool = True):
        self.seq = seq.upper().replace("T", "U")
        self.n = len(self.seq)
        self.min_loop = int(min_loop)
        self.canonical_only = bool(canonical_only)

        # ground set, and the per-index partner lists the recurrence walks
        self._pairable = [[] for _ in range(self.n)]
        ground = []
        for i in range(self.n):
            for j in range(i + self.min_loop + 1, self.n):
                if self.canonical_only and (self.seq[i], self.seq[j]) not in CANONICAL:
                    continue
                self._pairable[i].append(j)
                ground.append((i, j))
        self._ground = tuple(ground)
        self._ground_set = frozenset(ground)
        self._last_stats: dict = {}

    # ------------------------------------------------------------------ description
    def ground_set(self): return self._ground

    def describe(self):
        return {"family": self.family_id, "seq": self.seq, "n": self.n,
                "min_loop": self.min_loop, "canonical_only": self.canonical_only}

    def assumptions(self):
        return (f"min_loop={self.min_loop}",
                f"canonical_only={self.canonical_only}",
                "non-crossing (nested) structures only",
                "each index in at most one pair")

    def _hash(self, forced, forbidden):
        d = self.describe()
        d["forced"] = sorted(forced)
        d["forbidden"] = sorted(forbidden)
        return family_hash(d)

    # ------------------------------------------------------------------ feasibility
    def feasibility(self, forced=(), forbidden=()) -> Feasibility:
        forced, forbidden = normalise_mask(forced, forbidden)
        h = self._hash(forced, forbidden)
        A = self.assumptions()

        def no(reason):
            return Feasibility(False, reason, INFEASIBLE, h, A)

        both = forced & forbidden
        if both:
            return no(f"pair(s) both forced and forbidden: {sorted(both)}")

        for e in sorted(forced):
            i, j = e
            if not (0 <= i < j < self.n):
                return no(f"forced pair {e} out of range for n={self.n}")
            if j - i <= self.min_loop:
                return no(f"forced pair {e} violates min_loop={self.min_loop}")
            if e not in self._ground_set:
                return no(f"forced pair {e} is not canonical")

        seen = {}
        for (i, j) in sorted(forced):
            for idx in (i, j):
                if idx in seen:
                    return no(f"forced pairs {seen[idx]} and {(i, j)} share index {idx}")
                seen[idx] = (i, j)

        fl = sorted(forced)
        for a in range(len(fl)):
            for b in range(a + 1, len(fl)):
                if crosses(fl[a], fl[b]):
                    return no(f"forced pairs {fl[a]} and {fl[b]} cross")

        # the pre-checks above are sufficient for this family: pairwise non-crossing
        # forced pairs with distinct endpoints already form a valid member. Confirm
        # with the recurrence anyway, so a disagreement surfaces as a bug not a result.
        c = self._run(COUNTING, {}, forced, forbidden)
        if c["value"] == 0:
            return no("no structure satisfies the mask (recurrence returned count 0)")
        return Feasibility(True, f"{c['value']} structure(s) satisfy the mask",
                           PROVED, h, A)

    # ------------------------------------------------------------------ the recurrence
    def _run(self, semiring, coefficients, forced, forbidden, want_backtrace=False):
        n = self.n
        # Two forced pairs sharing an index are contradictory. Building `partner` as a
        # plain dict would let the later assignment overwrite the earlier one, and the
        # recurrence would then happily count structures satisfying only one of them --
        # a silently wrong count. Detect the conflict and return zero.
        partner = {}
        conflict = False
        for (i, j) in forced:
            for x, y in ((i, j), (j, i)):
                if partner.get(x, y) != y:
                    conflict = True
                partner[x] = y
        if conflict:
            self._last_stats = {"chart_values": 0, "materialized_productions": 0,
                                "n": n, "ground_set_size": len(self._ground)}
            return {"value": semiring.zero, "W": None, "choice": None,
                    "stats": dict(self._last_stats)}

        zero, one = semiring.zero, semiring.one
        add, mul = semiring.add, semiring.mul

        # W[i][j+1] holds the interval [i..j]; j = i-1 is the empty interval.
        W = [[zero] * (n + 2) for _ in range(n + 2)]
        choice = [[None] * (n + 2) for _ in range(n + 2)] if want_backtrace else None
        chart_values = 0
        productions = 0

        for i in range(n, -1, -1):
            for j in range(i - 1, n):
                chart_values += 1
                if i > j:                       # empty interval
                    W[i][j + 1] = one
                    continue

                p = partner.get(i)
                val = zero
                best = None

                if p is None:
                    # production: i unpaired
                    productions += 1
                    val = W[i + 1][j + 1]
                    if want_backtrace and val != zero:
                        best = Production("unpaired", (i, j), None)
                elif p < i or p > j:
                    # i is forced to a partner outside this interval: dead cell
                    W[i][j + 1] = zero
                    continue

                for k in self._pairable[i]:
                    if k > j:
                        break
                    if p is not None and k != p:
                        continue
                    if (i, k) in forbidden:
                        continue
                    q = partner.get(k)
                    if q is not None and q != i:
                        continue                 # k is forced to someone else
                    productions += 1
                    inner = W[i + 1][k]          # [i+1 .. k-1]
                    outer = W[k + 1][j + 1]      # [k+1 .. j]
                    if inner == zero or outer == zero:
                        continue
                    w = semiring.weight(coefficients.get((i, k), 0))
                    cand = mul(w, mul(inner, outer))
                    if want_backtrace:
                        if val == zero or cand > val:
                            best = Production("pair", (i, j), (i, k))
                        val = add(val, cand)
                    else:
                        val = add(val, cand)

                W[i][j + 1] = val
                if want_backtrace:
                    choice[i][j + 1] = best

        self._last_stats = {"chart_values": chart_values,
                            "materialized_productions": productions,
                            "n": n, "ground_set_size": len(self._ground)}
        return {"value": W[0][n] if n > 0 else one, "W": W, "choice": choice,
                "stats": dict(self._last_stats)}

    def _backtrace(self, choice, i, j, out):
        if i > j:
            return
        c = choice[i][j + 1]
        if c is None:
            return
        out.append(c)
        if c.kind == "unpaired":
            self._backtrace(choice, i + 1, j, out)
        else:
            k = c.edge[1]
            self._backtrace(choice, i + 1, k - 1, out)
            self._backtrace(choice, k + 1, j, out)

    # ------------------------------------------------------------------ public queries
    def support(self, coefficients, forced=(), forbidden=()) -> OracleResult:
        forced, forbidden = normalise_mask(forced, forbidden)
        h = self._hash(forced, forbidden)
        coefficients = dict(coefficients or {})
        exact = all(isinstance(v, int) for v in coefficients.values())
        status = EXACT_INTEGER if exact else FLOAT_UNVERIFIED

        r = self._run(MAX_PLUS, coefficients, forced, forbidden, want_backtrace=True)
        if r["value"] == NEG_INF:
            return OracleResult(NEG_INF, None, status, INFEASIBLE, h, "none",
                                self.assumptions(), r["stats"])
        prods: list[Production] = []
        if self.n > 0:
            self._backtrace(r["choice"], 0, self.n - 1, prods)
        pairs = tuple(sorted(p.edge for p in prods if p.edge is not None))
        st = dict(r["stats"]); st["productions"] = tuple(prods)
        return OracleResult(r["value"], pairs, status, PROVED, h, "exact",
                            self.assumptions(), st)

    def count(self, forced=(), forbidden=()) -> OracleResult:
        forced, forbidden = normalise_mask(forced, forbidden)
        h = self._hash(forced, forbidden)
        r = self._run(COUNTING, {}, forced, forbidden)
        term = PROVED if r["value"] > 0 else INFEASIBLE
        return OracleResult(r["value"], None, EXACT_INTEGER, term, h, "exact",
                            self.assumptions(), r["stats"])

    def log_partition(self, coefficients, forced=(), forbidden=()) -> OracleResult:
        forced, forbidden = normalise_mask(forced, forbidden)
        h = self._hash(forced, forbidden)
        r = self._run(LOG_SUM_EXP, dict(coefficients or {}), forced, forbidden)
        term = PROVED if r["value"] != NEG_INF else INFEASIBLE
        return OracleResult(r["value"], None, FLOAT_DIAGNOSTIC, term, h, "exact",
                            self.assumptions(), r["stats"])

    def log_inside(self, coefficients, forced=(), forbidden=()):
        """The inside chart under log-sum-exp, for stochastic backtrace.

        Public because sampling needs it: drawing from the Gibbs measure by conditioning
        one edge at a time costs |E| partition evaluations per draw, which is hopeless
        beyond a few hundred candidates. With the chart in hand a draw is a single
        top-down pass.
        """
        forced, forbidden = normalise_mask(forced, forbidden)
        r = self._run(LOG_SUM_EXP, dict(coefficients or {}), forced, forbidden)
        return {"W": r["W"], "value": r["value"], "forced": forced,
                "forbidden": forbidden, "stats": r["stats"]}

    def stats(self): return dict(self._last_stats)
